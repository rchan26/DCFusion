library(DCFusion)
library(HMCBLR)

load('nycflights13_sub_posteriors.RData')
rm(data_split_4, data_split_8, data_split_32, data_split_64, data_split_128, data_split_256)
rm(sub_posteriors_4, sub_posteriors_8, sub_posteriors_32, sub_posteriors_64, sub_posteriors_128)
seed <- 2016
nsamples <- 30000
C <- 16
dim <- 21
prior_means <- rep(0, dim)
prior_variances <- rep(1, dim)
n_cores <- parallel::detectCores()
ESS_threshold <- 0.5
CESS_0_threshold <- 0.2
CESS_j_threshold <- 0.05
diffusion_estimator <- 'NB'

##### Loading in Data #####

load_nycflights_data <- function(seed = NULL) {
  nyc_flights <- subset(nycflights13::flights,
                        select = c("arr_delay", "origin", "carrier", "dep_delay",
                                   "year", "day", "month", "hour", "distance"))
  nyc_flights <- nyc_flights[complete.cases(nyc_flights),]
  if (!is.null(seed)) {
    set.seed(seed)
    nyc_flights <- nyc_flights[sample(1:nrow(nyc_flights)),]
  }
  nyc_flights$delayed <- as.numeric(nyc_flights$arr_delay > 0)
  nyc_flights$dep_delayed <- as.numeric(nyc_flights$dep_delay > 0)
  nyc_flights$weekday <- weekdays(as.Date(paste(nyc_flights$year, "-", nyc_flights$month, "-", nyc_flights$day, sep = "")))
  nyc_flights$weekend <- as.numeric(nyc_flights$weekday %in% c("Saturday", "Sunday"))
  nyc_flights$night <- as.numeric(nyc_flights$hour >= 20 | nyc_flights$hour <= 5)
  carrier_names <- names(table(nyc_flights$carrier))
  for (i in 1:(length(carrier_names)-1)) {
    nyc_flights[carrier_names[i]] <- as.numeric(nyc_flights$carrier==carrier_names[i])
  }
  nyc_flights$EWR <- as.numeric(nyc_flights$origin=="EWR")
  nyc_flights$JFK <- as.numeric(nyc_flights$origin=="JFK")
  X <- subset(nyc_flights, select = c("dep_delayed",
                                      "weekend",
                                      "night",
                                      carrier_names[1:(length(carrier_names)-1)],
                                      "EWR",
                                      "JFK"))
  design_mat <- as.matrix(cbind(rep(1, nrow(X)), X))
  colnames(design_mat)[1] <- 'intercept'
  return(list('data' = cbind('delayed' = nyc_flights$delayed, X),
              'y' = nyc_flights$delayed,
              'X' = design_mat))
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

##### Sampling from sub-posterior C=16 #####

data_split_16 <- split_data(nyc_flights$data, y_col_index = 1, X_col_index = 2:dim, C = C, as_dataframe = F)
sub_posteriors_16 <- hmc_base_sampler_BLR(nsamples = nsamples,
                                          data_split = data_split_16,
                                          C = C,
                                          prior_means = prior_means,
                                          prior_variances = prior_variances,
                                          warmup = 10000,
                                          seed = seed,
                                          output = T)

print(paste('##### C:', C))

##### Applying other methodologies #####

print('Applying other methodologies')
consensus_mat_16 <- consensus_scott(S = C, samples_to_combine = sub_posteriors_16, indep = F)
consensus_sca_16 <- consensus_scott(S = C, samples_to_combine = sub_posteriors_16, indep = T)
neiswanger_true_16 <- neiswanger(S = C,
                                 samples_to_combine = sub_posteriors_16,
                                 anneal = TRUE)
neiswanger_false_16 <- neiswanger(S = C,
                                  samples_to_combine = sub_posteriors_16,
                                  anneal = FALSE)
weierstrass_importance_16 <- weierstrass(Samples = sub_posteriors_16,
                                         method = 'importance')
weierstrass_rejection_16 <- weierstrass(Samples = sub_posteriors_16,
                                        method = 'reject')

integrated_abs_distance(full_posterior, consensus_mat_16$samples)
integrated_abs_distance(full_posterior, consensus_sca_16$samples)
integrated_abs_distance(full_posterior, neiswanger_true_16$samples)
integrated_abs_distance(full_posterior, neiswanger_false_16$samples)
integrated_abs_distance(full_posterior, weierstrass_importance_16$samples)
integrated_abs_distance(full_posterior, weierstrass_rejection_16$samples)

##### bal binary combining two sub-posteriors at a time #####

# regular mesh
balanced_C16 <- list('reg' = bal_binary_GBF_BLR(N_schedule = rep(nsamples, 4),
                                                m_schedule = rep(2, 4),
                                                time_mesh = NULL,
                                                base_samples = sub_posteriors_16,
                                                L = 5,
                                                dim = dim,
                                                data_split = data_split_16,
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
                                                seed = seed,
                                                n_cores = n_cores,
                                                print_progress_iters = 500))
balanced_C16$reg$particles <- resample_particle_y_samples(particle_set = balanced_C16$reg$particles[[1]],
                                                          multivariate = TRUE,
                                                          resampling_method = 'resid',
                                                          seed = seed)
balanced_C16$reg$proposed_samples <- balanced_C16$reg$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, balanced_C16$reg$particles$y_samples))

# adaptive mesh
balanced_C16$adaptive <- bal_binary_GBF_BLR(N_schedule = rep(nsamples, 4),
                                            m_schedule = rep(2, 4),
                                            time_mesh = NULL,
                                            base_samples = sub_posteriors_16,
                                            L = 5,
                                            dim = dim,
                                            data_split = data_split_16,
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
                                            seed = seed,n_cores = n_cores,
                                            print_progress_iters = 500)
balanced_C16$adaptive$particles <- resample_particle_y_samples(particle_set = balanced_C16$adaptive$particles[[1]],
                                                               multivariate = TRUE,
                                                               resampling_method = 'resid',
                                                               seed = seed)
balanced_C16$adaptive$proposed_samples <- balanced_C16$adaptive$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, balanced_C16$adaptive$particles$y_samples))

save.image('NYC16.RData')
