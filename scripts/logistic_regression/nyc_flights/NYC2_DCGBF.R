library(DCFusion)
library(HMCBLR)

seed <- 2016
nsamples <- 30000
nsamples_GBF <- 5000
time_choice <- 0.5
C <- 2
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
                                 prior_variances = rep(10, 4),
                                 iterations = nsamples + 10000,
                                 warmup = 10000,
                                 chains = 1,
                                 seed = seed,
                                 output = T)

##### Sampling from sub-posterior C=2 #####

data_split_2 <- split_data(nyc_flights$data, y_col_index = 1, X_col_index = 2:4, C = C, as_dataframe = F)
sub_posteriors_2 <- hmc_base_sampler_BLR(nsamples = nsamples,
                                         data_split = data_split_2,
                                         C = C,
                                         prior_means = rep(0, 4),
                                         prior_variances = rep(10, 4),
                                         warmup = 10000,
                                         seed = seed,
                                         output = T)

# compare_samples_bivariate(sub_posteriors_2, c('red', 'red', 'blue', 'blue'), c(-4, 4))

##### Applying other methodologies #####

print('Applying other methodologies')
consensus_mat_2 <- consensus_scott(S = C, samples_to_combine = sub_posteriors_2, indep = F)
consensus_sca_2 <- consensus_scott(S = C, samples_to_combine = sub_posteriors_2, indep = T)
neiswanger_true_2 <- neiswanger(S = C,
                                samples_to_combine = sub_posteriors_2,
                                anneal = TRUE)
neiswanger_false_2 <- neiswanger(S = C,
                                 samples_to_combine = sub_posteriors_2,
                                 anneal = FALSE)
weierstrass_importance_2 <- weierstrass(Samples = sub_posteriors_2,
                                        method = 'importance')
weierstrass_rejection_2 <- weierstrass(Samples = sub_posteriors_2,
                                       method = 'reject')

integrated_abs_distance(full_posterior, consensus_mat_2$samples)
integrated_abs_distance(full_posterior, consensus_sca_2$samples)
integrated_abs_distance(full_posterior, neiswanger_true_2$samples)
integrated_abs_distance(full_posterior, neiswanger_false_2$samples)
integrated_abs_distance(full_posterior, weierstrass_importance_2$samples)
integrated_abs_distance(full_posterior, weierstrass_rejection_2$samples)

# ##### NB (Hypercube Centre) #####
# print('NB Fusion (hypercube centre)')
# NB_hc_2 <- bal_binary_fusion_SMC_BLR(N_schedule = nsamples,
#                                      m_schedule = 2,
#                                      time_schedule = time_choice,
#                                      base_samples = sub_posteriors_2,
#                                      L = 2,
#                                      dim = 4,
#                                      data_split = data_split_2,
#                                      prior_means = rep(0, 4),
#                                      prior_variances = rep(10, 4),
#                                      C = C,
#                                      precondition = TRUE,
#                                      resampling_method = 'resid',
#                                      ESS_threshold = 0.5,
#                                      cv_location = 'hypercube_centre',
#                                      diffusion_estimator = 'NB',
#                                      seed = seed,
#                                      n_cores = n_cores,
#                                      print_progress_iters = 500)
# NB_hc_2$particles <- resample_particle_y_samples(particle_set = NB_hc_2$particles[[1]],
#                                                  multivariate = TRUE,
#                                                  resampling_method = 'resid',
#                                                  seed = seed)
# NB_hc_2$proposed_samples <- NB_hc_2$proposed_samples[[1]]
# print(integrated_abs_distance(full_posterior, NB_hc_2$particles$y_samples))

##### bal binary combining two sub-posteriors at a time #####
# balanced_C2 <- list('reg' = bal_binary_GBF_BLR(N_schedule = nsamples_GBF,
#                                                m_schedule = 2,
#                                                time_mesh = NULL,
#                                                base_samples = sub_posteriors_2,
#                                                L = 2,
#                                                dim = 4,
#                                                data_split = data_split_2,
#                                                prior_means = rep(0, 4),
#                                                prior_variances = rep(10, 4),
#                                                C = C,
#                                                precondition = TRUE,
#                                                resampling_method = 'resid',
#                                                ESS_threshold = ESS_threshold,
#                                                adaptive_mesh = FALSE,
#                                                mesh_parameters = list('condition' = 'SH',
#                                                                       'lambda' = 5,
#                                                                       'CESS_0_threshold' = CESS_0_threshold,
#                                                                       'CESS_j_threshold' = CESS_j_threshold,
#                                                                       'vanilla' = FALSE),
#                                                diffusion_estimator = diffusion_estimator,
#                                                seed = seed))
balanced_C2 <- list('adaptive' = bal_binary_GBF_BLR(N_schedule = nsamples_GBF,
                                                    m_schedule = 2,
                                                    time_mesh = NULL,
                                                    base_samples = sub_posteriors_2,
                                                    L = 2,
                                                    dim = 4,
                                                    data_split = data_split_2,
                                                    prior_means = rep(0, 4),
                                                    prior_variances = rep(10, 4),
                                                    C = C,
                                                    precondition = TRUE,
                                                    resampling_method = 'resid',
                                                    ESS_threshold = ESS_threshold,
                                                    adaptive_mesh = TRUE,
                                                    mesh_parameters = list('condition' = 'SSH',
                                                                           'CESS_0_threshold' = CESS_0_threshold,
                                                                           'CESS_j_threshold' = CESS_j_threshold,
                                                                           'vanilla' = FALSE),
                                                    diffusion_estimator = diffusion_estimator,
                                                    local_bounds = FALSE,
                                                    record = TRUE,
                                                    seed = seed))

# # regular mesh
# balanced_C2$reg$particles <- resample_particle_y_samples(particle_set = balanced_C2$reg$particles[[1]],
#                                                          multivariate = TRUE,
#                                                          resampling_method = 'resid',
#                                                          seed = seed)
# balanced_C2$reg$proposed_samples <- balanced_C2$reg$proposed_samples[[1]]
# print(integrated_abs_distance(full_posterior, balanced_C2$reg$particles$y_samples))
# adaptive mesh
balanced_C2$adaptive$particles <- resample_particle_y_samples(particle_set = balanced_C2$adaptive$particles[[1]],
                                                              multivariate = TRUE,
                                                              resampling_method = 'resid',
                                                              seed = seed)
balanced_C2$adaptive$proposed_samples <- balanced_C2$adaptive$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, balanced_C2$adaptive$particles$y_samples))
compare_samples_bivariate(list(full_posterior,
                               balanced_C2$adaptive$particles$y_samples,
                               balanced_C2$adaptive$proposed_samples),
                          c('black', 'red', 'green'),
                          c(-4, 4))

save.image('NYC2_DCGBF.RData')
