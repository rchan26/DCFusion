library(DCFusion)
library(HMCBLR)

seed <- 2016
nsamples <- 10000
nsamples_GBF <- 5000
time_choice <- 1
C <- 128
prior_means <- rep(0, 11)
prior_variances <- rep(1, 11)
n_cores <- parallel::detectCores()
ESS_threshold <- 0.5
CESS_0_threshold <- 0.25
CESS_j_threshold <- 0.01
diffusion_estimator <- 'NB'

##### Loading in Data #####

load_nycflights_data <- function(seed = NULL) {
  nyc_flights <- subset(nycflights13::flights %>% left_join(nycflights13::weather),
                        select = c("arr_delay", "origin", "dep_delay",
                                   "year", "day", "month", "hour", "distance",
                                   "temp", "humid", "wind_speed", "wind_gust",
                                   "precip", "pressure", "visib"))
  nyc_flights <- nyc_flights[complete.cases(nyc_flights),]
  if (!is.null(seed)) {
    set.seed(seed)
    nyc_flights <- nyc_flights[sample(1:nrow(nyc_flights)),]
  }
  nyc_flights$delayed <- as.numeric(nyc_flights$arr_delay > 0)
  # binary variables
  nyc_flights$dep_delayed <- as.numeric(nyc_flights$dep_delay > 0)
  nyc_flights$weekday <- weekdays(as.Date(paste(nyc_flights$year, "-", nyc_flights$month, "-", nyc_flights$day, sep = "")))
  nyc_flights$weekend <- as.numeric(nyc_flights$weekday %in% c("Saturday", "Sunday"))
  nyc_flights$night <- as.numeric(nyc_flights$hour >= 20 | nyc_flights$hour <= 5)
  nyc_flights$EWR <- as.numeric(nyc_flights$origin=="EWR")
  nyc_flights$JFK <- as.numeric(nyc_flights$origin=="JFK")
  nyc_flights$windy <- as.numeric(nyc_flights$wind_speed > 25 | nyc_flights$wind_gust > 35)
  nyc_flights$rain <- as.numeric(nyc_flights$precip > 0.01)
  nyc_flights$visibility <- as.numeric(nyc_flights$visib < 8)
  # standardise continuous variables
  minimums <- list('distance' = min(nyc_flights$distance),
                   'temp' = min(nyc_flights$temp),
                   'humid' = min(nyc_flights$humid),
                   'wind_speed' = min(nyc_flights$wind_speed),
                   'wind_gust' = min(nyc_flights$wind_gust),
                   'precip' = min(nyc_flights$precip),
                   'pressure' = min(nyc_flights$pressure),
                   'visib' = min(nyc_flights$visib))
  ranges <- list('distance' = max(nyc_flights$distance)-minimums$distance,
                 'temp' = max(nyc_flights$temp)-minimums$temp,
                 'humid' = max(nyc_flights$humid)-minimums$humid,
                 'wind_speed' = max(nyc_flights$wind_speed)-minimums$wind_speed,
                 'wind_gust' = max(nyc_flights$wind_gust)-minimums$wind_gust,
                 'precip' = max(nyc_flights$precip)-minimums$precip,
                 'pressure' = max(nyc_flights$pressure)-minimums$pressure,
                 'visib' = max(nyc_flights$visib)-minimums$visib)
  nyc_flights$distance_standardised <- (nyc_flights$distance-minimums$distance) / ranges$distance
  nyc_flights$temp_standardised <- (nyc_flights$temp-minimums$temp) / ranges$temp
  nyc_flights$humid_standardised <- (nyc_flights$humid-minimums$humid) / ranges$humid
  nyc_flights$wind_speed_standardised <- (nyc_flights$wind_speed-minimums$wind_speed) / ranges$wind_speed
  nyc_flights$wind_gust_standardised <- (nyc_flights$wind_gust-minimums$wind_gust) / ranges$wind_gust
  nyc_flights$precip_standardised <- (nyc_flights$precip-minimums$precip) / ranges$precip
  nyc_flights$pressure_standardised <- (nyc_flights$pressure-minimums$pressure) / ranges$pressure
  nyc_flights$visib_standardised <- (nyc_flights$visib-minimums$visib) / ranges$visib
  X <- subset(nyc_flights, select = c("dep_delayed",
                                      "weekend",
                                      "night",
                                      "EWR",
                                      "JFK",
                                      "windy",
                                      "rain",
                                      "visibility",
                                      "distance_standardised",
                                      "temp_standardised"
                                      # "humid_standardised",
                                      # "wind_speed_standardised",
                                      # "wind_gust_standardised",
                                      # "precip_standardised",
                                      # "pressure_standardised",
                                      # "visib_standardised"
  ))
  design_mat <- as.matrix(cbind(rep(1, nrow(X)), X))
  colnames(design_mat)[1] <- 'intercept'
  return(list('data' = cbind('delayed' = nyc_flights$delayed, X),
              'y' = nyc_flights$delayed,
              'X' = design_mat,
              'minimums' = minimums,
              'ranges' = ranges))
}

nyc_flights <- load_nycflights_data(seed)

##### Sampling from full posterior #####

full_data_count <- unique_row_count(nyc_flights$y, nyc_flights$X)$full_data_count
full_posterior <- hmc_sample_BLR(full_data_count = full_data_count,
                                 C = 1,
                                 prior_means = prior_means,
                                 prior_variances = prior_variances,
                                 iterations = nsamples + 10000,
                                 warmup = 10000,
                                 chains = 1,
                                 seed = seed,
                                 output = T)

##### Sampling from sub-posterior C=128 #####

data_split_128 <- split_data(nyc_flights$data, y_col_index = 1, X_col_index = 2:11, C = C, as_dataframe = F)
sub_posteriors_128 <- hmc_base_sampler_BLR(nsamples = nsamples,
                                           data_split = data_split_128,
                                           C = C,
                                           prior_means = prior_means,
                                           prior_variances = prior_variances,
                                           warmup = 10000,
                                           seed = seed,
                                           output = T)

##### Applying other methodologies #####

print('Applying other methodologies')
consensus_mat_128 <- consensus_scott(S = C, samples_to_combine = sub_posteriors_128, indep = F)
consensus_sca_128 <- consensus_scott(S = C, samples_to_combine = sub_posteriors_128, indep = T)
neiswanger_true_128 <- neiswanger(S = C,
                                  samples_to_combine = sub_posteriors_128,
                                  anneal = TRUE)
neiswanger_false_128 <- neiswanger(S = C,
                                   samples_to_combine = sub_posteriors_128,
                                   anneal = FALSE)
weierstrass_importance_128 <- weierstrass(Samples = sub_posteriors_128,
                                          method = 'importance')
weierstrass_rejection_128 <- weierstrass(Samples = sub_posteriors_128,
                                         method = 'reject')

integrated_abs_distance(full_posterior, consensus_mat_128$samples)
integrated_abs_distance(full_posterior, consensus_sca_128$samples)
integrated_abs_distance(full_posterior, neiswanger_true_128$samples)
integrated_abs_distance(full_posterior, neiswanger_false_128$samples)
integrated_abs_distance(full_posterior, weierstrass_importance_128$samples)
integrated_abs_distance(full_posterior, weierstrass_rejection_128$samples)

##### NB (Hypercube Centre) #####
print('NB Fusion (hypercube centre)')
NB_hc_128 <- bal_binary_fusion_SMC_BLR(N_schedule = rep(nsamples, 7),
                                       m_schedule = rep(2, 7),
                                       time_schedule = rep(time_choice, 7),
                                       base_samples = sub_posteriors_128,
                                       L = 8,
                                       dim = 4,
                                       data_split = data_split_128,
                                       prior_means = prior_means,
                                       prior_variances = prior_variances,
                                       C = C,
                                       precondition = TRUE,
                                       resampling_method = 'resid',
                                       ESS_threshold = 0.5,
                                       cv_location = 'hypercube_centre',
                                       diffusion_estimator = 'NB',
                                       seed = seed,
                                       n_cores = n_cores,
                                       print_progress_iters = 500)
NB_hc_128$particles <- resample_particle_y_samples(particle_set = NB_hc_128$particles[[1]],
                                                   multivariate = TRUE,
                                                   resampling_method = 'resid',
                                                   seed = seed)
NB_hc_128$proposed_samples <- NB_hc_128$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, NB_hc_128$particles$y_samples))

##### bal binary combining two sub-posteriors at a time #####
balanced_C128 <- list('reg' = bal_binary_GBF_BLR(N_schedule = rep(nsamples_GBF, 7),
                                                 m_schedule = rep(2, 7),
                                                 time_mesh = NULL,
                                                 base_samples = sub_posteriors_128,
                                                 L = 8,
                                                 dim = 4,
                                                 data_split = data_split_128,
                                                 prior_means = prior_means,
                                                 prior_variances = prior_variances,
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
balanced_C128$adaptive <- bal_binary_GBF_BLR(N_schedule = rep(nsamples_GBF, 7),
                                             m_schedule = rep(2, 7),
                                             time_mesh = NULL,
                                             base_samples = sub_posteriors_128,
                                             L = 8,
                                             dim = 4,
                                             data_split = data_split_128,
                                             prior_means = prior_means,
                                             prior_variances = prior_variances,
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
balanced_C128$reg$particles <- resample_particle_y_samples(particle_set = balanced_C128$reg$particles[[1]],
                                                           multivariate = TRUE,
                                                           resampling_method = 'resid',
                                                           seed = seed)
balanced_C128$reg$proposed_samples <- balanced_C128$reg$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, balanced_C128$reg$particles$y_samples))
# adaptive mesh
balanced_C128$adaptive$particles <- resample_particle_y_samples(particle_set = balanced_C128$adaptive$particles[[1]],
                                                                multivariate = TRUE,
                                                                resampling_method = 'resid',
                                                                seed = seed)
balanced_C128$adaptive$proposed_samples <- balanced_C128$adaptive$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, balanced_C128$adaptive$particles$y_samples))

save.image('NYC128_DCGBF.RData')
