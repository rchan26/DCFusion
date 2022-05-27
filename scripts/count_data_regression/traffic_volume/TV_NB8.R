library(DCFusion)
library(HMCGLMR)

##### Initialise example #####
seed <- 2022
set.seed(seed)
nsamples_MCF <- 10000
nsamples_GBF <- 10000
warmup <- 10000
time_choice <- 1
phi_rate <- 1
prior_means <- rep(0, 9)
prior_variances <- rep(10, 9)
ESS_threshold <- 0.5
CESS_0_threshold <- 0.5
CESS_j_threshold <- 0.05
diffusion_estimator <- 'NB'
n_cores <- parallel::detectCores()

##### Loading in Data #####

load_tv_data <- function(file, seed = NULL, standardise_variables = TRUE) {
  original_data <- read.csv(file)
  original_data$hour <- format(strptime(original_data$date_time, format = "%Y-%m-%d %H:%M:%S"), format = "%H")
  traffic_volume <- data.frame(tv = original_data$traffic_volume,
                               holiday = as.numeric(original_data$holiday!="None"),
                               weekend = as.numeric(weekdays(as.Date(original_data$date_time)) %in% c("Saturday", "Sunday")),
                               weather = as.numeric(original_data$weather_main %in% c("Clear", "Clouds")),
                               rush_hour = as.numeric(original_data$hour %in% c("07", "08", "09", "16", "17", "18", "19")),
                               temp = original_data$temp,
                               rain = original_data$rain_1h,
                               snow = original_data$snow_1h,
                               clouds = original_data$clouds_all)
  traffic_volume <- traffic_volume[complete.cases(traffic_volume),]
  if (!is.null(seed)) {
    set.seed(seed)
    traffic_volume <- traffic_volume[sample(1:nrow(traffic_volume)),]
  }
  if (standardise_variables) {
    X <- subset(traffic_volume, select = -c(holiday, weekend, weather, rush_hour, tv))
    variable_means <- rep(NA, ncol(X))
    variable_sds <- rep(NA, ncol(X))
    for (col in 1:ncol(X)) {
      variable_means[col] <- mean(X[,col])
      variable_sds[col] <- sd(X[,col])
      X[,col] <- (X[,col]-variable_means[col])/variable_sds[col]
    }
    design_mat <- cbind('intercept' = rep(1, nrow(X)),
                        subset(traffic_volume, select = c(holiday, weekend, weather, rush_hour)),
                        X)
    return(list('data' = cbind(subset(design_mat, select = -c(intercept)),
                               'tv' = traffic_volume$tv),
                'y' = traffic_volume$tv,
                'X' = design_mat,
                'variable_means' = variable_means,
                'variable_sds' = variable_sds))
  } else {
    X <- subset(traffic_volume, select = -c(tv))
    design_mat <- as.matrix(cbind(rep(1, nrow(X)), X))
    colnames(design_mat)[1] <- 'intercept'
    return(list('data' = traffic_volume,
                'y' = traffic_volume$tv,
                'X' = design_mat))
  }
}

traffic_volume <- load_tv_data('scripts/robust_regression/traffic_volume/traffic_volume.csv', seed = seed)

##### Sampling from full posterior #####

full_posterior <- hmc_sample_GLMR(likelihood = 'NB',
                                  y = traffic_volume$y,
                                  X = traffic_volume$X,
                                  C = 1,
                                  phi = phi_rate,
                                  prior_means = prior_means,
                                  prior_variances = prior_variances,
                                  iterations = nsamples_MCF + 10000,
                                  warmup = 10000,
                                  chains = 1,
                                  seed = seed,
                                  output = T)

##### Sampling from sub-posterior C=8 #####

data_split_8 <- split_data(traffic_volume$data,
                           y_col_index = 9,
                           X_col_index = 1:8,
                           C = 8,
                           as_dataframe = F)
sub_posteriors_8 <- hmc_base_sampler_GLMR(likelihood = 'NB',
                                          nsamples = nsamples_MCF,
                                          data_split = data_split_8,
                                          C = 8,
                                          phi = phi_rate,
                                          prior_means = prior_means,
                                          prior_variances = prior_variances,
                                          warmup = 10000,
                                          seed = seed,
                                          output = T)

##### Applying other methodologies #####

print('Applying other methodologies')
consensus_mat_8 <- consensus_scott(S = 8, samples_to_combine = sub_posteriors_8, indep = F)
consensus_sca_8 <- consensus_scott(S = 8, samples_to_combine = sub_posteriors_8, indep = T)
neiswanger_true_8 <- neiswanger(S = 8,
                                samples_to_combine = sub_posteriors_8,
                                anneal = TRUE)
neiswanger_false_8 <- neiswanger(S = 8,
                                 samples_to_combine = sub_posteriors_8,
                                 anneal = FALSE)
weierstrass_importance_8 <- weierstrass(Samples = sub_posteriors_8,
                                        method = 'importance')
weierstrass_rejection_8 <- weierstrass(Samples = sub_posteriors_8,
                                       method = 'reject')

##### NB (Hypercube Centre) #####
print('NB Fusion (hypercube centre)')
NB_hc_8 <- bal_binary_fusion_SMC_BNBR(N_schedule = rep(nsamples_MCF, 3),
                                      m_schedule = rep(2, 3),
                                      time_schedule = rep(time_choice, 3),
                                      base_samples = sub_posteriors_8,
                                      L = 4,
                                      dim = 9,
                                      data_split = data_split_8,
                                      phi = phi_rate,
                                      prior_means = prior_means,
                                      prior_variances = prior_variances,
                                      C = 8,
                                      precondition = TRUE,
                                      resampling_method = 'resid',
                                      ESS_threshold = ESS_threshold,
                                      diffusion_estimator = 'NB',
                                      seed = seed,
                                      n_cores = n_cores)
NB_hc_8$particles <- resample_particle_y_samples(particle_set = NB_hc_8$particles[[1]],
                                                 multivariate = TRUE,
                                                 resampling_method = 'resid',
                                                 seed = seed)
NB_hc_8$proposed_samples <- NB_hc_8$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, NB_hc_8$particles$y_samples))

##### Generalised Bayesian Fusion #####

##### all at once #####
GBF_8 <- list('reg' = bal_binary_GBF_BNBR(N_schedule = nsamples_GBF,
                                          m_schedule = 8,
                                          time_mesh = NULL,
                                          base_samples = sub_posteriors_8,
                                          L = 2,
                                          dim = 9,
                                          data_split = data_split_8,
                                          phi = phi_rate,
                                          prior_means = prior_means,
                                          prior_variances = prior_variances,
                                          C = 8,
                                          precondition = TRUE,
                                          resampling_method = 'resid',
                                          ESS_threshold = ESS_threshold,
                                          adaptive_mesh = FALSE,
                                          mesh_parameters = list('condition' = 'SH',
                                                                 'CESS_0_threshold' = CESS_0_threshold,
                                                                 'CESS_j_threshold' = CESS_j_threshold,
                                                                 'vanilla' = FALSE),
                                          diffusion_estimator = diffusion_estimator,
                                          seed = seed))
GBF_8$adaptive <-  bal_binary_GBF_BNBR(N_schedule = nsamples_GBF,
                                       m_schedule = 8,
                                       time_mesh = NULL,
                                       base_samples = sub_posteriors_8,
                                       L = 2,
                                       dim = 9,
                                       data_split = data_split_8,
                                       phi = phi_rate,
                                       prior_means = prior_means,
                                       prior_variances = prior_variances,
                                       C = 8,
                                       precondition = TRUE,
                                       resampling_method = 'resid',
                                       ESS_threshold = ESS_threshold,
                                       adaptive_mesh = TRUE,
                                       mesh_parameters = list('condition' = 'SH',
                                                              'CESS_0_threshold' = CESS_0_threshold,
                                                              'CESS_j_threshold' = CESS_j_threshold,
                                                              'vanilla' = FALSE),
                                       diffusion_estimator = diffusion_estimator,
                                       seed = seed)

# regular mesh
GBF_8$reg$particles <- resample_particle_y_samples(particle_set = GBF_8$reg$particles[[1]],
                                                   multivariate = TRUE,
                                                   resampling_method = 'resid',
                                                   seed = seed)
print(integrated_abs_distance(full_posterior, GBF_8$reg$particles$y_samples)) 
# adaptive mesh
GBF_8$adaptive$particles <- resample_particle_y_samples(particle_set = GBF_8$adaptive$particles[[1]],
                                                        multivariate = TRUE,
                                                        resampling_method = 'resid',
                                                        seed = seed)
print(integrated_abs_distance(full_posterior, GBF_8$adaptive$particles$y_samples))

##### bal binary combining two sub-posteriors at a time #####
balanced_C8 <- list('reg' = bal_binary_GBF_BNBR(N_schedule = rep(nsamples_GBF, 3),
                                                m_schedule = rep(2, 3),
                                                time_mesh = NULL,
                                                base_samples = sub_posteriors_8,
                                                L = 4,
                                                dim = 9,
                                                data_split = data_split_8,
                                                phi = phi_rate,
                                                prior_means = prior_means,
                                                prior_variances = prior_variances,
                                                C = 8,
                                                precondition = TRUE,
                                                resampling_method = 'resid',
                                                ESS_threshold = ESS_threshold,
                                                adaptive_mesh = FALSE,
                                                mesh_parameters = list('condition' = 'SH',
                                                                       'CESS_0_threshold' = CESS_0_threshold,
                                                                       'CESS_j_threshold' = CESS_j_threshold,
                                                                       'vanilla' = FALSE),
                                                diffusion_estimator = diffusion_estimator,
                                                seed = seed))
balanced_C8$adaptive <- bal_binary_GBF_BNBR(N_schedule = rep(nsamples_GBF, 3),
                                            m_schedule = rep(2, 3),
                                            time_mesh = NULL,
                                            base_samples = sub_posteriors_8,
                                            L = 4,
                                            dim = 9,
                                            data_split = data_split_8,
                                            phi = phi_rate,
                                            prior_means = prior_means,
                                            prior_variances = prior_variances,
                                            C = 8,
                                            precondition = TRUE,
                                            resampling_method = 'resid',
                                            ESS_threshold = ESS_threshold,
                                            adaptive_mesh = TRUE,
                                            mesh_parameters = list('condition' = 'SH',
                                                                   'CESS_0_threshold' = CESS_0_threshold,
                                                                   'CESS_j_threshold' = CESS_j_threshold,
                                                                   'vanilla' = FALSE),
                                            diffusion_estimator = diffusion_estimator,
                                            seed = seed)

# regular mesh
balanced_C8$reg$particles <- resample_particle_y_samples(particle_set = balanced_C8$reg$particles[[1]],
                                                         multivariate = TRUE,
                                                         resampling_method = 'resid',
                                                         seed = seed)
print(integrated_abs_distance(full_posterior, balanced_C8$reg$particles$y_samples))
# adaptive mesh
balanced_C8$adaptive$particles <- resample_particle_y_samples(particle_set = balanced_C8$adaptive$particles[[1]],
                                                              multivariate = TRUE,
                                                              resampling_method = 'resid',
                                                              seed = seed)
print(integrated_abs_distance(full_posterior, balanced_C8$adaptive$particles$y_samples))

##### IAD #####

integrated_abs_distance(full_posterior, GBF_8$reg$particles$y_samples)
integrated_abs_distance(full_posterior, GBF_8$adaptive$particles$y_samples)
integrated_abs_distance(full_posterior, balanced_C8$reg$particles$y_samples)
integrated_abs_distance(full_posterior, balanced_C8$adaptive$particles$y_samples)
integrated_abs_distance(full_posterior, Poisson_hc_8$particles$y_samples)
integrated_abs_distance(full_posterior, NB_hc_8$particles$y_samples)
integrated_abs_distance(full_posterior, consensus_mat_8$samples)
integrated_abs_distance(full_posterior, consensus_sca_8$samples)
integrated_abs_distance(full_posterior, neiswanger_true_8$samples)
integrated_abs_distance(full_posterior, neiswanger_false_8$samples)
integrated_abs_distance(full_posterior, weierstrass_importance_8$samples)
integrated_abs_distance(full_posterior, weierstrass_rejection_8$samples)

save.image('TV8.RData')
