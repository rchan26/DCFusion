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
CESS_j_threshold <- 0.2
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

##### Sampling from sub-posterior C=128 #####

data_split_128 <- split_data(traffic_volume$data,
                             y_col_index = 9,
                             X_col_index = 1:8,
                             C = 128,
                             as_dataframe = F)
sub_posteriors_128 <- hmc_base_sampler_GLMR(likelihood = 'NB',
                                            nsamples = nsamples_MCF,
                                            data_split = data_split_128,
                                            C = 128,
                                            phi = phi_rate,
                                            prior_means = prior_means,
                                            prior_variances = prior_variances,
                                            warmup = 10000,
                                            seed = seed,
                                            output = T)

##### Applying other methodologies #####

print('Applying other methodologies')
consensus_mat_128 <- consensus_scott(S = 128, samples_to_combine = sub_posteriors_128, indep = F)
consensus_sca_128 <- consensus_scott(S = 128, samples_to_combine = sub_posteriors_128, indep = T)
neiswanger_true_128 <- neiswanger(S = 128,
                                  samples_to_combine = sub_posteriors_128,
                                  anneal = TRUE)
neiswanger_false_128 <- neiswanger(S = 128,
                                   samples_to_combine = sub_posteriors_128,
                                   anneal = FALSE)
weierstrass_importance_128 <- weierstrass(Samples = sub_posteriors_128,
                                          method = 'importance')
weierstrass_rejection_128 <- weierstrass(Samples = sub_posteriors_128,
                                         method = 'reject')

##### Poisson (Hypercube Centre) #####
print('Poisson Fusion (hypercube centre)')
Poisson_hc_128 <- bal_binary_fusion_SMC_BNBR(N_schedule = rep(nsamples_MCF, 7),
                                             m_schedule = rep(2, 7),
                                             time_schedule = rep(time_choice, 7),
                                             base_samples = sub_posteriors_128,
                                             L = 8,
                                             dim = 9,
                                             data_split = data_split_128,
                                             phi = phi_rate,
                                             prior_means = prior_means,
                                             prior_variances = prior_variances,
                                             C = 128,
                                             precondition = TRUE,
                                             resampling_method = 'resid',
                                             ESS_threshold = ESS_threshold,
                                             diffusion_estimator = 'Poisson',
                                             seed = seed,
                                             n_cores = n_cores)
Poisson_hc_128$particles <- resample_particle_y_samples(particle_set = Poisson_hc_128$particles[[1]],
                                                        multivariate = TRUE,
                                                        resampling_method = 'resid',
                                                        seed = seed)
Poisson_hc_128$proposed_samples <- Poisson_hc_128$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, Poisson_hc_128$particles$y_samples))

##### NB (Hypercube Centre) #####
print('NB Fusion (hypercube centre)')
NB_hc_128 <- bal_binary_fusion_SMC_BNBR(N_schedule = rep(nsamples_MCF, 7),
                                        m_schedule = rep(2, 7),
                                        time_schedule = rep(time_choice, 7),
                                        base_samples = sub_posteriors_128,
                                        L = 8,
                                        dim = 9,
                                        data_split = data_split_128,
                                        phi = phi_rate,
                                        prior_means = prior_means,
                                        prior_variances = prior_variances,
                                        C = 128,
                                        precondition = TRUE,
                                        resampling_method = 'resid',
                                        ESS_threshold = ESS_threshold,
                                        diffusion_estimator = 'NB',
                                        seed = seed,
                                        n_cores = n_cores)
NB_hc_128$particles <- resample_particle_y_samples(particle_set = NB_hc_128$particles[[1]],
                                                   multivariate = TRUE,
                                                   resampling_method = 'resid',
                                                   seed = seed)
NB_hc_128$proposed_samples <- NB_hc_128$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, NB_hc_128$particles$y_samples))

##### Generalised Bayesian Fusion #####

##### all at once #####
GBF_128 <- list('reg' = bal_binary_GBF_BNBR(N_schedule = nsamples_GBF,
                                            m_schedule = 128,
                                            time_mesh = NULL,
                                            base_samples = sub_posteriors_128,
                                            L = 2,
                                            dim = 9,
                                            data_split = data_split_128,
                                            phi = phi_rate,
                                            prior_means = prior_means,
                                            prior_variances = prior_variances,
                                            C = 128,
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
GBF_128$adaptive <-  bal_binary_GBF_BNBR(N_schedule = nsamples_GBF,
                                         m_schedule = 128,
                                         time_mesh = NULL,
                                         base_samples = sub_posteriors_128,
                                         L = 2,
                                         dim = 9,
                                         data_split = data_split_128,
                                         phi = phi_rate,
                                         prior_means = prior_means,
                                         prior_variances = prior_variances,
                                         C = 128,
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
GBF_128$reg$particles <- resample_particle_y_samples(particle_set = GBF_128$reg$particles[[1]],
                                                     multivariate = TRUE,
                                                     resampling_method = 'resid',
                                                     seed = seed)
print(integrated_abs_distance(full_posterior, GBF_128$reg$particles$y_samples)) 
# adaptive mesh
GBF_128$adaptive$particles <- resample_particle_y_samples(particle_set = GBF_128$adaptive$particles[[1]],
                                                          multivariate = TRUE,
                                                          resampling_method = 'resid',
                                                          seed = seed)
print(integrated_abs_distance(full_posterior, GBF_128$adaptive$particles$y_samples))

##### bal binary combining two sub-posteriors at a time #####
balanced_C128 <- list('reg' = bal_binary_GBF_BNBR(N_schedule = rep(nsamples_GBF, 7),
                                                  m_schedule = rep(2, 7),
                                                  time_mesh = NULL,
                                                  base_samples = sub_posteriors_128,
                                                  L = 8,
                                                  dim = 9,
                                                  data_split = data_split_128,
                                                  phi = phi_rate,
                                                  prior_means = prior_means,
                                                  prior_variances = prior_variances,
                                                  C = 128,
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
balanced_C128$adaptive <- bal_binary_GBF_BNBR(N_schedule = rep(nsamples_GBF, 7),
                                              m_schedule = rep(2, 7),
                                              time_mesh = NULL,
                                              base_samples = sub_posteriors_128,
                                              L = 8,
                                              dim = 9,
                                              data_split = data_split_128,
                                              phi = phi_rate,
                                              prior_means = prior_means,
                                              prior_variances = prior_variances,
                                              C = 128,
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
balanced_C128$reg$particles <- resample_particle_y_samples(particle_set = balanced_C128$reg$particles[[1]],
                                                           multivariate = TRUE,
                                                           resampling_method = 'resid',
                                                           seed = seed)
print(integrated_abs_distance(full_posterior, balanced_C128$reg$particles$y_samples))
# adaptive mesh
balanced_C128$adaptive$particles <- resample_particle_y_samples(particle_set = balanced_C128$adaptive$particles[[1]],
                                                                multivariate = TRUE,
                                                                resampling_method = 'resid',
                                                                seed = seed)
print(integrated_abs_distance(full_posterior, balanced_C128$adaptive$particles$y_samples))

##### IAD #####

integrated_abs_distance(full_posterior, GBF_128$reg$particles$y_samples)
integrated_abs_distance(full_posterior, GBF_128$adaptive$particles$y_samples)
integrated_abs_distance(full_posterior, balanced_C128$reg$particles$y_samples)
integrated_abs_distance(full_posterior, balanced_C128$adaptive$particles$y_samples)
integrated_abs_distance(full_posterior, Poisson_hc_128$particles$y_samples)
integrated_abs_distance(full_posterior, NB_hc_128$particles$y_samples)
integrated_abs_distance(full_posterior, consensus_mat_128$samples)
integrated_abs_distance(full_posterior, consensus_sca_128$samples)
integrated_abs_distance(full_posterior, neiswanger_true_128$samples)
integrated_abs_distance(full_posterior, neiswanger_false_128$samples)
integrated_abs_distance(full_posterior, weierstrass_importance_128$samples)
integrated_abs_distance(full_posterior, weierstrass_rejection_128$samples)

save.image('TV128.RData')
