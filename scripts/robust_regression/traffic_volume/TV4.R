library(DCFusion)
library(HMCBRR)

##### Initialise example #####
seed <- 2022
set.seed(seed)
nsamples_MCF <- 100000
nsamples_GBF <- 10000
warmup <- 10000
time_choice <- 1
nu <- 5
sigma <- 1
prior_means <- rep(0, 6)
prior_variances <- rep(10, 6)
ESS_threshold <- 0.5
CESS_0_threshold <- 0.5
CESS_j_threshold <- 0.2
diffusion_estimator <- 'NB'
n_cores <- parallel::detectCores()

##### Loading in Data #####

load_tv_data <- function(file, standardise_variables = TRUE) {
  original_data <- read.csv('scripts/robust_regression/traffic_volume/traffic_volume.csv')
  traffic_volume <- data.frame(tv = original_data$traffic_volume,
                               temp = original_data$temp,
                               holiday = as.numeric(original_data$holiday!="None"),
                               weather = as.numeric(original_data$weather_main %in% c("Clear", "Clouds")),
                               hour = format(strptime(original_data$date_time, format = "%Y-%m-%d %H:%M:%S"), format = "%H"))
  traffic_volume$weekday <- weekdays(as.Date(original_data$date_time))
  traffic_volume$weekend <- as.numeric(traffic_volume$weekday %in% c("Saturday", "Sunday"))
  traffic_volume$rush_hour <- as.numeric(!(traffic_volume$weekday %in% c("Saturday", "Sunday")) 
                                        & (traffic_volume$hour %in% c("07", "08", "09", "16", "17", "18", "19")))
  if (standardise_variables) {
    X <- subset(traffic_volume, select = c(temp))
    variable_means <- rep(NA, ncol(X))
    variable_sds <- rep(NA, ncol(X))
    for (col in 1:ncol(X)) {
      variable_means[col] <- mean(X[,col])
      variable_sds[col] <- sd(X[,col])
      X[,col] <- (X[,col]-variable_means[col])/variable_sds[col]
    }
    design_mat <- cbind('intercept' = rep(1, nrow(X)),
                        X,
                        subset(traffic_volume, select = c(weather, holiday, weekend, rush_hour)))
    return(list('data' = cbind(subset(design_mat, select = -c(intercept)),
                               'tv' = traffic_volume$tv),
                'y' = traffic_volume$tv,
                'X' = design_mat,
                'variable_means' = variable_means,
                'variable_sds' = variable_sds))
  } else {
    X <- subset(traffic_volume, select = c(temp, weather, holiday, weekend, rush_hour))
    design_mat <- as.matrix(cbind(rep(1, nrow(X)), X))
    colnames(design_mat)[1] <- 'intercept'
    return(list('data' = cbind(subset(design_mat, select = -c(intercept)),
                               'tv' = traffic_volume$tv),
                'y' = traffic_volume$tv,
                'X' = design_mat))
  }
}

traffic_volume <- load_tv_data('scripts/robust_regression/traffic_volume/traffic_volume.csv')

##### Sampling from full posterior #####

full_posterior <- hmc_sample_BRR(noise_error = 'student_t',
                                 y = traffic_volume$y,
                                 X = traffic_volume$X,
                                 C = 1,
                                 nu = nu,
                                 sigma = sigma,
                                 prior_means = prior_means,
                                 prior_variances = prior_variances,
                                 iterations = nsamples_MCF + 10000,
                                 warmup = 10000,
                                 chains = 1,
                                 seed = seed,
                                 output = T)

##### Sampling from sub-posterior C=4 #####

data_split_4 <- split_data(traffic_volume$data,
                           y_col_index = 6,
                           X_col_index = 1:5,
                           C = 4,
                           as_dataframe = F)
sub_posteriors_4 <- hmc_base_sampler_BRR(nsamples = nsamples_MCF,
                                         data_split = data_split_4,
                                         C = 4,
                                         nu = nu,
                                         sigma = sigma,
                                         prior_means = prior_means,
                                         prior_variances = prior_variances,
                                         warmup = 10000,
                                         seed = seed,
                                         output = T)

##### Applying other methodologies #####

print('Applying other methodologies')
consensus_mat_4 <- consensus_scott(S = 4, samples_to_combine = sub_posteriors_4, indep = F)
consensus_sca_4 <- consensus_scott(S = 4, samples_to_combine = sub_posteriors_4, indep = T)
neiswanger_true_4 <- neiswanger(S = 4,
                                samples_to_combine = sub_posteriors_4,
                                anneal = TRUE)
neiswanger_false_4 <- neiswanger(S = 4,
                                 samples_to_combine = sub_posteriors_4,
                                 anneal = FALSE)
weierstrass_importance_4 <- weierstrass(Samples = sub_posteriors_4,
                                        method = 'importance')
weierstrass_rejection_4 <- weierstrass(Samples = sub_posteriors_4,
                                       method = 'reject')

# ##### Poisson (Hypercube Centre) #####
# print('Poisson Fusion (hypercube centre)')
# Poisson_hc_4 <- bal_binary_fusion_SMC_BRR(N_schedule = rep(nsamples_MCF, 2),
#                                           m_schedule = rep(2, 2),
#                                           time_schedule = rep(time_choice, 2),
#                                           base_samples = sub_posteriors_4,
#                                           L = 3,
#                                           dim = 8,
#                                           data_split = data_split_4,
#                                           nu = nu,
#                                           sigma = sigma,
#                                           prior_means = prior_means,
#                                           prior_variances = prior_variances,
#                                           C = 4,
#                                           precondition = TRUE,
#                                           resampling_method = 'resid',
#                                           ESS_threshold = ESS_threshold,
#                                           cv_location = 'hypercube_centre',
#                                           diffusion_estimator = 'Poisson',
#                                           seed = seed,
#                                           n_cores = n_cores)
# Poisson_hc_4$particles <- resample_particle_y_samples(particle_set = Poisson_hc_4$particles[[1]],
#                                                       multivariate = TRUE,
#                                                       resampling_method = 'resid',
#                                                       seed = seed)
# Poisson_hc_4$proposed_samples <- Poisson_hc_4$proposed_samples[[1]]
# print(integrated_abs_distance(full_posterior, Poisson_hc_4$particles$y_samples))

##### NB (Hypercube Centre) #####
print('NB Fusion (hypercube centre)')
NB_hc_4 <- bal_binary_fusion_SMC_BRR(N_schedule = rep(nsamples_MCF, 2),
                                     m_schedule = rep(2, 2),
                                     time_schedule = rep(time_choice, 2),
                                     base_samples = sub_posteriors_4,
                                     L = 3,
                                     dim = 8,
                                     data_split = data_split_4,
                                     nu = nu,
                                     sigma = sigma,
                                     prior_means = prior_means,
                                     prior_variances = prior_variances,
                                     C = 4,
                                     precondition = TRUE,
                                     resampling_method = 'resid',
                                     ESS_threshold = ESS_threshold,
                                     cv_location = 'hypercube_centre',
                                     diffusion_estimator = 'NB',
                                     seed = seed,
                                     n_cores = n_cores)
NB_hc_4$particles <- resample_particle_y_samples(particle_set = NB_hc_4$particles[[1]],
                                                 multivariate = TRUE,
                                                 resampling_method = 'resid',
                                                 seed = seed)
NB_hc_4$proposed_samples <- NB_hc_4$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, NB_hc_4$particles$y_samples))

##### Generalised Bayesian Fusion #####

##### all at once #####
GBF_4 <- list('reg' = bal_binary_GBF_BRR(N_schedule = nsamples_GBF,
                                         m_schedule = 4,
                                         time_mesh = NULL,
                                         base_samples = sub_posteriors_4,
                                         L = 2,
                                         dim = 8,
                                         data_split = data_split_4,
                                         nu = nu,
                                         sigma = sigma,
                                         prior_means = prior_means,
                                         prior_variances = prior_variances,
                                         C = 4,
                                         precondition = TRUE,
                                         resampling_method = 'resid',
                                         ESS_threshold = ESS_threshold,
                                         adaptive_mesh = FALSE,
                                         mesh_parameters = list('condition' = 'SH',
                                                                'CESS_0_threshold' = CESS_0_threshold,
                                                                'CESS_j_threshold' = CESS_j_threshold,
                                                                'vanilla' = FALSE),
                                         diffusion_estimator = diffusion_estimator,
                                         seed = seed),
              'adaptive' = bal_binary_GBF_BRR(N_schedule = nsamples_GBF,
                                              m_schedule = 4,
                                              time_mesh = NULL,
                                              base_samples = sub_posteriors_4,
                                              L = 2,
                                              dim = 8,
                                              data_split = data_split_4,
                                              nu = nu,
                                              sigma = sigma,
                                              prior_means = prior_means,
                                              prior_variances = prior_variances,
                                              C = 4,
                                              precondition = TRUE,
                                              resampling_method = 'resid',
                                              ESS_threshold = ESS_threshold,
                                              adaptive_mesh = TRUE,
                                              mesh_parameters = list('condition' = 'SH',
                                                                     'CESS_0_threshold' = CESS_0_threshold,
                                                                     'CESS_j_threshold' = CESS_j_threshold,
                                                                     'vanilla' = FALSE),
                                              diffusion_estimator = diffusion_estimator,
                                              seed = seed))

# regular mesh
GBF_4$reg$particles <- resample_particle_y_samples(particle_set = GBF_4$reg$particles[[1]],
                                                   multivariate = TRUE,
                                                   resampling_method = 'resid',
                                                   seed = seed)
print(integrated_abs_distance(full_posterior, GBF_4$reg$particles$y_samples)) 
compare_samples_bivariate(posteriors = list(full_posterior,
                                            GBF_4$reg$proposed_samples[[1]],
                                            GBF_4$reg$particles$y_samples),
                          colours = c('black', 'green', 'red'),
                          common_limit = c(-4, 4))
# adaptive mesh
GBF_4$adaptive$particles <- resample_particle_y_samples(particle_set = GBF_4$adaptive$particles[[1]],
                                                        multivariate = TRUE,
                                                        resampling_method = 'resid',
                                                        seed = seed)
print(integrated_abs_distance(full_posterior, GBF_4$adaptive$particles$y_samples))
compare_samples_bivariate(posteriors = list(full_posterior,
                                            GBF_4$adaptive$proposed_samples[[1]],
                                            GBF_4$adaptive$particles$y_samples),
                          colours = c('black', 'green', 'red'),
                          common_limit = c(-4, 4))

##### bal binary combining two sub-posteriors at a time #####
balanced_C4 <- list('reg' = bal_binary_GBF_BRR(N_schedule = rep(nsamples_GBF, 2),
                                               m_schedule = rep(2, 2),
                                               time_mesh = NULL,
                                               base_samples = sub_posteriors_4,
                                               L = 3,
                                               dim = 8,
                                               data_split = data_split_4,
                                               nu = nu,
                                               sigma = sigma,
                                               prior_means = prior_means,
                                               prior_variances = prior_variances,
                                               C = 4,
                                               precondition = TRUE,
                                               resampling_method = 'resid',
                                               ESS_threshold = ESS_threshold,
                                               adaptive_mesh = FALSE,
                                               mesh_parameters = list('condition' = 'SH',
                                                                      'CESS_0_threshold' = CESS_0_threshold,
                                                                      'CESS_j_threshold' = CESS_j_threshold,
                                                                      'vanilla' = FALSE),
                                               diffusion_estimator = diffusion_estimator,
                                               seed = seed),
                    'adaptive' = bal_binary_GBF_BRR(N_schedule = rep(nsamples_GBF, 2),
                                                    m_schedule = rep(2, 2),
                                                    time_mesh = NULL,
                                                    base_samples = sub_posteriors_4,
                                                    L = 3,
                                                    dim = 8,
                                                    data_split = data_split_4,
                                                    nu = nu,
                                                    sigma = sigma,
                                                    prior_means = prior_means,
                                                    prior_variances = prior_variances,
                                                    C = 4,
                                                    precondition = TRUE,
                                                    resampling_method = 'resid',
                                                    ESS_threshold = ESS_threshold,
                                                    adaptive_mesh = TRUE,
                                                    mesh_parameters = list('condition' = 'SH',
                                                                           'CESS_0_threshold' = CESS_0_threshold,
                                                                           'CESS_j_threshold' = CESS_j_threshold,
                                                                           'vanilla' = FALSE),
                                                    diffusion_estimator = diffusion_estimator,
                                                    seed = seed))

# regular mesh
balanced_C4$reg$particles <- resample_particle_y_samples(particle_set = balanced_C4$reg$particles[[1]],
                                                         multivariate = TRUE,
                                                         resampling_method = 'resid',
                                                         seed = seed)
print(integrated_abs_distance(full_posterior, balanced_C4$reg$particles$y_samples))
compare_samples_bivariate(posteriors = list(full_posterior,
                                            balanced_C4$reg$proposed_samples[[1]],
                                            balanced_C4$reg$particles$y_samples),
                          colours = c('black', 'green', 'red'),
                          common_limit = c(-4, 4))
# adaptive mesh
balanced_C4$adaptive$particles <- resample_particle_y_samples(particle_set = balanced_C4$adaptive$particles[[1]],
                                                              multivariate = TRUE,
                                                              resampling_method = 'resid',
                                                              seed = seed)
print(integrated_abs_distance(full_posterior, balanced_C4$adaptive$particles$y_samples))
compare_samples_bivariate(posteriors = list(full_posterior,
                                            balanced_C4$adaptive$proposed_samples[[1]],
                                            balanced_C4$adaptive$particles$y_samples),
                          colours = c('black', 'green', 'red'),
                          common_limit = c(-4, 4))

##### IAD #####

integrated_abs_distance(full_posterior, GBF_4$reg$particles$y_samples)
integrated_abs_distance(full_posterior, GBF_4$adaptive$particles$y_samples)
integrated_abs_distance(full_posterior, balanced_C4$reg$particles$y_samples)
integrated_abs_distance(full_posterior, balanced_C4$adaptive$particles$y_samples)
# integrated_abs_distance(full_posterior, Poisson_hc_4$particles$y_samples)
integrated_abs_distance(full_posterior, NB_hc_4$particles$y_samples)
integrated_abs_distance(full_posterior, consensus_mat_4$samples)
integrated_abs_distance(full_posterior, consensus_sca_4$samples)
integrated_abs_distance(full_posterior, neiswanger_true_4$samples)
integrated_abs_distance(full_posterior, neiswanger_false_4$samples)
integrated_abs_distance(full_posterior, weierstrass_importance_4$samples)
integrated_abs_distance(full_posterior, weierstrass_rejection_4$samples)

save.image('TV4.RData')
