library(DCFusion)
library(HMCBRR)

seed <- 2016
nsamples <- 30000
time_choice <- 0.5
nu <- 10
sigma <- 1
n_cores <- parallel::detectCores()

# note that the reproducibility of this script depends on the number of cores available
# various functions will be done in parallel
# can be changed by changing n_cores above

##### Loading in Data #####

original_data <- read.csv('traffic_volume.csv')
# holiday: 1 for yes (i.e. not "None") or 0 for no
traffic_volume <- data.frame(holiday = as.numeric(original_data$holiday!="None"))
# temp: continuous variable
temp_mean <- mean(original_data$temp)
temp_sd <- sd(original_data$temp)
# traffic_volume <- data.frame('temp' = (original_data$temp-temp_mean)/temp_sd)
traffic_volume$temp <- (original_data$temp-temp_mean)/temp_sd
# rain: continuous variable
rain_mean <- mean(original_data$rain_1h)
rain_sd <- sd(original_data$rain_1h)
traffic_volume$rain <- (original_data$rain_1h-rain_mean)/rain_sd
# snow: continuous variable
snow_mean <- mean(original_data$snow_1h)
snow_sd <- sd(original_data$snow_1h)
traffic_volume$snow <- (original_data$snow_1h-snow_mean)/snow_sd
# clouds: continuous variable
clouds_mean <- mean(original_data$clouds_all)
clouds_sd <- sd(original_data$clouds_all)
traffic_volume$clouds <- (original_data$clouds_all-clouds_mean)/clouds_sd
# weather_main: 1 for good (i.e. "Clear" or "Clouds") or 0 for bad
traffic_volume$weather <- as.numeric(original_data$weather_main %in% c("Clear", "Clouds"))
# time: 1 for rush hour (i.e. between 07:00-09:00 and 16:00-19:00) or 0 for no
original_data$hour <- format(strptime(original_data$date_time, format = "%Y-%m-%d %H:%M:%S"), format = "%H")
traffic_volume$rush_hour <- as.numeric(original_data$hour %in% c("07", "08", "09", "16", "17", "18", "19"))
# y (response variable): traffic volume
traffic_volume$y <- original_data$traffic_volume
rm(original_data)

##### Sampling from full posterior #####

traffic_volume_data <- list()
traffic_volume_data$y <- traffic_volume$y
traffic_volume_data$X <- as.matrix(cbind(rep(1, nrow(traffic_volume[,1:7])), traffic_volume[,1:7]))
colnames(traffic_volume_data$X)[1] <- 'intercept'
full_data_count <- unique_row_count(traffic_volume_data$y, traffic_volume_data$X)$full_data_count

# lm(y ~., data = traffic_volume)

full_posterior <- hmc_sample_BRR(full_data_count = full_data_count,
                                 C = 1,
                                 nu = nu,
                                 sigma = sigma,
                                 prior_means = rep(0, 8),
                                 prior_variances = rep(1, 8),
                                 iterations = nsamples + 10000,
                                 warmup = 10000,
                                 chains = 1,
                                 seed = seed,
                                 output = T)

##### Sampling from sub-posterior C=8 #####

data_split_8 <- split_data(traffic_volume,
                           y_col_index = 8,
                           X_col_index = 1:7,
                           C = 8,
                           as_dataframe = F)
sub_posteriors_8 <- hmc_base_sampler_BRR(nsamples = nsamples,
                                         data_split = data_split_8,
                                         C = 8,
                                         nu = nu,
                                         sigma = sigma,
                                         prior_means = rep(0, 8),
                                         prior_variances = rep(1, 8),
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

##### Applying Fusion #####

##### Poisson (Hypercube Centre) #####
print('Poisson Fusion (hypercube centre)')
Poisson_hc_8 <- bal_binary_fusion_SMC_BLR(N_schedule = rep(nsamples, 3),
                                          m_schedule = rep(2, 3),
                                          time_schedule = rep(time_choice, 3),
                                          base_samples = sub_posteriors_8,
                                          L = 6,
                                          dim = 5,
                                          data_split = data_split_8,
                                          prior_means = rep(0, 8),
                                          prior_variances = rep(1, 8),
                                          C = 8,
                                          precondition = TRUE,
                                          resampling_method = 'resid',
                                          ESS_threshold = 0.5,
                                          cv_location = 'hypercube_centre',
                                          diffusion_estimator = 'Poisson',
                                          seed = seed,
                                          n_cores = n_cores,
                                          print_progress_iters = 500)
Poisson_hc_8$particles <- resample_particle_y_samples(particle_set = Poisson_hc_8$particles[[1]],
                                                      multivariate = TRUE,
                                                      resampling_method = 'resid',
                                                      seed = seed)
Poisson_hc_8$proposed_samples <- Poisson_hc_8$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior,
                              Poisson_hc_8$particles$y_samples))

##### NB (Hypercube Centre) #####
print('NB Fusion (hypercube centre)')
NB_hc_8 <- bal_binaryfusion_SMC_BLR(N_schedule = rep(nsamples, 3),
                                    m_schedule = rep(2, 3),
                                    time_schedule = rep(time_choice, 3),
                                    base_samples = sub_posteriors_8,
                                    L = 6,
                                    dim = 5,
                                    data_split = data_split_8,
                                    prior_means = rep(0, 8),
                                    prior_variances = rep(1, 8),
                                    C = 8,
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
print(integrated_abs_distance(full_posterior,
                              NB_hc_8$particles$y_samples))

##### IAD #####

integrated_abs_distance(full_posterior, Poisson_hc_8$particles$y_samples)
integrated_abs_distance(full_posterior, NB_hc_8$particles$y_samples)
integrated_abs_distance(full_posterior, consensus_mat_8$samples)
integrated_abs_distance(full_posterior, consensus_sca_8$samples)
integrated_abs_distance(full_posterior, neiswanger_true_8$samples)
integrated_abs_distance(full_posterior, neiswanger_false_8$samples)
integrated_abs_distance(full_posterior, weierstrass_importance_8$samples)
integrated_abs_distance(full_posterior, weierstrass_rejection_8$samples)

##### Save data #####

save.image('TV8.RData')
