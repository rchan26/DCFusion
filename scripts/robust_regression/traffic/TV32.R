library(DCFusion)
library(HMCBRR)

seed <- 2016
nsamples <- 30000
time_choice <- 0.5
n_cores <- parallel::detectCores()

# note that the reproducibility of this script depends on the number of cores available
# various functions will be done in parallel
# can be changed by changing n_cores above

##### Loading in Data #####

traffic_volume <- read.csv('traffic_volume.csv')
# holiday: 1 for yes (i.e. not "None") or 0 for no
traffic_volume_full <- data.frame(holiday = as.numeric(traffic_volume$holiday!="None"))
# temp: continuous variable
temp_mean <- mean(traffic_volume$temp)
temp_sd <- sd(traffic_volume$temp)
traffic_volume_full$temp <- (traffic_volume$temp-temp_mean)/temp_sd
# rain: continuous variable
rain_mean <- mean(traffic_volume$rain_1h)
rain_sd <- sd(traffic_volume$rain_1h)
traffic_volume_full$rain <- (traffic_volume$rain_1h-rain_mean)/rain_sd
# snow: continuous variable
snow_mean <- mean(traffic_volume$snow_1h)
snow_sd <- sd(traffic_volume$snow_1h)
traffic_volume_full$snow <- (traffic_volume$snow_1h-snow_mean)/snow_sd
# clouds: continuous variable
clouds_mean <- mean(traffic_volume$clouds_all)
clouds_sd <- sd(traffic_volume$clouds_all)
traffic_volume_full$clouds <- (traffic_volume$clouds_all-clouds_mean)/clouds_sd
# weather_main: 1 for good (i.e. "Clear" or "Clouds") or 0 for bad
traffic_volume_full$weather <- as.numeric(traffic_volume$weather_main %in% c("Clear", "Clouds"))
# time: 1 for rush hour (i.e. between 07:00-09:00 and 16:00-19:00) or 0 for no
traffic_volume$hour <- format(strptime(traffic_volume$date_time, format = "%Y-%m-%d %H:%M:%S"), format = "%H")
traffic_volume_full$rush_hour <- as.numeric(traffic_volume$hour %in% c("07", "08", "09", "16", "17", "18", "19"))
# y (response variable): traffic volume
traffic_volume_full$y <- traffic_volume$traffic_volume

##### Sampling from full posterior #####

traffic_volume_data <- list()
traffic_volume_data$y <- traffic_volume_full$y
traffic_volume_data$X <- as.matrix(cbind(rep(1, nrow(traffic_volume_full[,1:7])), traffic_volume_full[,1:7]))
colnames(traffic_volume_data$X)[1] <- 'intercept'

full_posterior <- hmc_sample_BRR(data = traffic_volume_data,
                                 C = 1,
                                 nu = 5,
                                 sigma = 1,
                                 prior_means = rep(0, 8),
                                 prior_variances = rep(1, 8),
                                 iterations = nsamples + 10000,
                                 warmup = 10000,
                                 chains = 1,
                                 seed = seed,
                                 output = T)

##### Sampling from sub-posterior C=32 #####

data_split_32 <- split_data(traffic_volume_full,
                            y_col_index = 8,
                            X_col_index = 1:7,
                            C = 32,
                            as_dataframe = F)
sub_posteriors_32 <- hmc_base_sampler_BRR(nsamples = nsamples,
                                          data_split = data_split_32,
                                          C = 32,
                                          nu = 5,
                                          sigma = 1,
                                          prior_means = rep(0, 8),
                                          prior_variances = rep(1, 8),
                                          warmup = 10000,
                                          seed = seed,
                                          output = T)

##### Applying other methodologies #####

print('Applying other methodologies')
consensus_mat_32 <- consensus_scott(S = 32, samples_to_combine = sub_posteriors_32, indep = F)
consensus_sca_32 <- consensus_scott(S = 32, samples_to_combine = sub_posteriors_32, indep = T)
neiswanger_true_32 <- neiswanger(S = 32,
                                 samples_to_combine = sub_posteriors_32,
                                 anneal = TRUE)
neiswanger_false_32 <- neiswanger(S = 32,
                                  samples_to_combine = sub_posteriors_32,
                                  anneal = FALSE)
weierstrass_importance_32 <- weierstrass(Samples = sub_posteriors_32,
                                         method = 'importance')
weierstrass_rejection_32 <- weierstrass(Samples = sub_posteriors_32,
                                        method = 'reject')

##### Applying Fusion #####

##### Poisson (Hypercube Centre) #####
print('Poisson Fusion (hypercube centre)')
Poisson_hc_32 <- bal_binary_fusion_SMC_BLR(N_schedule = rep(nsamples, 5),
                                           m_schedule = rep(2, 5),
                                           time_schedule = rep(time_choice, 5),
                                           base_samples = sub_posteriors_32,
                                           L = 6,
                                           dim = 5,
                                           data_split = data_split_32,
                                           prior_means = rep(0, 8),
                                           prior_variances = rep(1, 8),
                                           C = 32,
                                           precondition = TRUE,
                                           resampling_method = 'resid',
                                           ESS_threshold = 0.5,
                                           cv_location = 'hypercube_centre',
                                           diffusion_estimator = 'Poisson',
                                           seed = seed,
                                           n_cores = n_cores,
                                           print_progress_iters = 500)
Poisson_hc_32$particles <- resample_particle_y_samples(particle_set = Poisson_hc_32$particles[[1]],
                                                       multivariate = TRUE,
                                                       resampling_method = 'resid',
                                                       seed = seed)
Poisson_hc_32$proposed_samples <- Poisson_hc_32$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior,
                              Poisson_hc_32$particles$y_samples))

##### NB (Hypercube Centre) #####
print('NB Fusion (hypercube centre)')
NB_hc_32 <- bal_binaryfusion_SMC_BLR(N_schedule = rep(nsamples, 5),
                                     m_schedule = rep(2, 5),
                                     time_schedule = rep(time_choice, 5),
                                     base_samples = sub_posteriors_32,
                                     L = 6,
                                     dim = 5,
                                     data_split = data_split_32,
                                     prior_means = rep(0, 8),
                                     prior_variances = rep(1, 8),
                                     C = 32,
                                     precondition = TRUE,
                                     resampling_method = 'resid',
                                     ESS_threshold = 0.5,
                                     cv_location = 'hypercube_centre',
                                     diffusion_estimator = 'NB',
                                     seed = seed,
                                     n_cores = n_cores,
                                     print_progress_iters = 500)
NB_hc_32$particles <- resample_particle_y_samples(particle_set = NB_hc_32$particles[[1]],
                                                  multivariate = TRUE,
                                                  resampling_method = 'resid',
                                                  seed = seed)
NB_hc_32$proposed_samples <- NB_hc_32$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior,
                              NB_hc_32$particles$y_samples))

##### IAD #####

integrated_abs_distance(full_posterior, Poisson_hc_32$particles$y_samples)
integrated_abs_distance(full_posterior, NB_hc_32$particles$y_samples)
integrated_abs_distance(full_posterior, consensus_mat_32$samples)
integrated_abs_distance(full_posterior, consensus_sca_32$samples)
integrated_abs_distance(full_posterior, neiswanger_true_32$samples)
integrated_abs_distance(full_posterior, neiswanger_false_32$samples)
integrated_abs_distance(full_posterior, weierstrass_importance_32$samples)
integrated_abs_distance(full_posterior, weierstrass_rejection_32$samples)

##### Save data #####

save.image('TV32.RData')
