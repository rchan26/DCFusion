library(DCFusion)
library(HMCBLR)

seed <- 2016
nsamples <- 30000
nsamples_GBF <- 5000
time_choice <- 0.5
C <- 8
n_cores <- parallel::detectCores()
ESS_threshold <- 0.5
CESS_0_threshold <- 0.5
CESS_j_threshold <- 0.05
diffusion_estimator <- 'NB'

##### Loading in Data #####

load_nycflights_data <- function() {
  nyc_flights <- subset(nycflights13::flights, select = c("arr_delay", "year", "day", "month", "hour", "distance"))
  nyc_flights <- nyc_flights[complete.cases(nyc_flights),]
  nyc_flights$weekday <- weekdays(as.Date(paste(nyc_flights$year, "-", nyc_flights$month, "-", nyc_flights$day, sep = "")))
  nyc_flights$delayed <- as.numeric(nyc_flights$arr_delay > 15)
  nyc_flights$weekend <- as.numeric(nyc_flights$weekday %in% c("Saturday", "Sunday"))
  nyc_flights$night <- as.numeric(nyc_flights$hour >= 20 | nyc_flights$hour <= 5)
  distance_min <- min(nyc_flights$distance)
  distance_range <- max(nyc_flights$distance)-distance_min
  nyc_flights$distance_standardised <- (nyc_flights$distance-distance_min) / distance_range
  X <- subset(nyc_flights, select = c("weekend", "night", "distance_standardised"))
  design_mat <- as.matrix(cbind(rep(1, nrow(X)), X))
  colnames(design_mat)[1] <- 'intercept'
  return(list('data' = cbind('delayed' = nyc_flights$delayed, X),
              'y' = nyc_flights$delayed,
              'X' = design_mat,
              'distance_min' = distance_min,
              'distance_range' = distance_range))
}

nyc_flights <- load_nycflights_data()

##### Sampling from full posterior #####

full_data_count <- unique_row_count(nyc_flights$y, nyc_flights$X)$full_data_count
full_posterior <- hmc_sample_BLR(full_data_count = full_data_count,
                                 C = 1,
                                 prior_means = rep(0, 4),
                                 prior_variances = rep(1, 4),
                                 iterations = nsamples + 10000,
                                 warmup = 10000,
                                 chains = 1,
                                 seed = seed,
                                 output = T)

##### Sampling from sub-posterior C=4 #####

data_split_8 <- split_data(nyc_flights$data, y_col_index = 1, X_col_index = 2:4, C = C, as_dataframe = F)
sub_posteriors_8 <- hmc_base_sampler_BLR(nsamples = nsamples,
                                         data_split = data_split_8,
                                         C = C,
                                         prior_means = rep(0, 4),
                                         prior_variances = rep(1, 4),
                                         warmup = 10000,
                                         seed = seed,
                                         output = T)

##### Applying other methodologies #####

print('Applying other methodologies')
consensus_mat_8 <- consensus_scott(S = C, samples_to_combine = sub_posteriors_8, indep = F)
consensus_sca_8 <- consensus_scott(S = C, samples_to_combine = sub_posteriors_8, indep = T)
neiswanger_true_8 <- neiswanger(S = C,
                                samples_to_combine = sub_posteriors_8,
                                anneal = TRUE)
neiswanger_false_8 <- neiswanger(S = C,
                                 samples_to_combine = sub_posteriors_8,
                                 anneal = FALSE)
weierstrass_importance_8 <- weierstrass(Samples = sub_posteriors_8,
                                        method = 'importance')
weierstrass_rejection_8 <- weierstrass(Samples = sub_posteriors_8,
                                       method = 'reject')

integrated_abs_distance(full_posterior, consensus_mat_8$samples)
integrated_abs_distance(full_posterior, consensus_sca_8$samples)
integrated_abs_distance(full_posterior, neiswanger_true_8$samples)
integrated_abs_distance(full_posterior, neiswanger_false_8$samples)
integrated_abs_distance(full_posterior, weierstrass_importance_8$samples)
integrated_abs_distance(full_posterior, weierstrass_rejection_8$samples)

##### NB (Hypercube Centre) #####
print('NB Fusion (hypercube centre)')
NB_hc_8 <- bal_binary_fusion_SMC_BLR(N_schedule = rep(nsamples, 3),
                                     m_schedule = rep(2, 3),
                                     time_schedule = rep(time_choice, 3),
                                     base_samples = sub_posteriors_8,
                                     L = 4,
                                     dim = 4,
                                     data_split = data_split_8,
                                     prior_means = rep(0, 4),
                                     prior_variances = rep(1, 4),
                                     C = C,
                                     precondition = TRUE,
                                     resampling_method = 'resid',
                                     ESS_threshold = 0.5,
                                     cv_location = 'hypercube_centre',
                                     diffusion_estimator = 'NB',
                                     seed = seed,
                                     n_cores = n_cores,
                                     print_progress_iters = 500)
NB_hc_8$particles <- resample_particle_y_samples(particle_set = NB_hc_8$particles[[1]],
                                                 multivariate = TRUE,
                                                 resampling_method = 'resid',
                                                 seed = seed)
NB_hc_8$proposed_samples <- NB_hc_8$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, NB_hc_8$particles$y_samples))

##### bal binary combining two sub-posteriors at a time #####
balanced_C4 <- list('reg' = bal_binary_GBF_BLR(N_schedule = rep(nsamples_GBF, 3),
                                               m_schedule = rep(2, 3),
                                               time_mesh = NULL,
                                               base_samples = sub_posteriors_8,
                                               L = 4,
                                               dim = 4,
                                               data_split = data_split_8,
                                               prior_means = rep(0, 4),
                                               prior_variances = rep(1, 4),
                                               C = C,
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
balanced_C4$adaptive <- bal_binary_GBF_BLR(N_schedule = rep(nsamples_GBF, 3),
                                           m_schedule = rep(2, 3),
                                           time_mesh = NULL,
                                           base_samples = sub_posteriors_8,
                                           L = 4,
                                           dim = 4,
                                           data_split = data_split_8,
                                           prior_means = rep(0, 4),
                                           prior_variances = rep(1, 4),
                                           C = C,
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
balanced_C4$reg$particles <- resample_particle_y_samples(particle_set = balanced_C4$reg$particles[[1]],
                                                         multivariate = TRUE,
                                                         resampling_method = 'resid',
                                                         seed = seed)
balanced_C4$reg$proposed_samples <- balanced_C4$reg$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, balanced_C4$reg$particles$y_samples))
# adaptive mesh
balanced_C4$adaptive$particles <- resample_particle_y_samples(particle_set = balanced_C4$adaptive$particles[[1]],
                                                              multivariate = TRUE,
                                                              resampling_method = 'resid',
                                                              seed = seed)
balanced_C4$adaptive$proposed_samples <- balanced_C4$adaptive$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, balanced_C4$adaptive$particles$y_samples))

save.image('NYC4_DCGBF.RData')
