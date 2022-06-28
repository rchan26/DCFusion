library(DCFusion)
library(HMCGLMR)

##### Initialise example #####
seed <- 2022
set.seed(seed)
nsamples <- 10000
warmup <- 10000
phi_rate <- 1
prior_means <- rep(0, 10)
prior_variances <- rep(10, 10)
ESS_threshold <- 0.5
CESS_0_threshold <- 0.2
CESS_j_threshold <- 0.05
diffusion_estimator <- 'NB'
n_cores <- parallel::detectCores()

##### Loading in Data #####

load_bs_data <- function(file, seed = NULL) {
  original_data <- read.csv(file)
  bike_sharing <- data.frame(rented = original_data$cnt,
                             spring = as.numeric(original_data$season == 2),
                             summer = as.numeric(original_data$season == 3),
                             fall = as.numeric(original_data$season == 4),
                             weekend = as.numeric(original_data$weekday %in% c(6, 0)),
                             holiday = original_data$holiday,
                             rush_hour = as.numeric(!(original_data$weekday %in% c(6, 0)) 
                                                    & (original_data$hr %in% c(7, 8, 9, 16, 17, 18, 19))),
                             weather = as.numeric(original_data$weathersit == 1),
                             atemp = original_data$atemp,
                             wind = original_data$windspeed)
  bike_sharing <- bike_sharing[complete.cases(bike_sharing),]
  if (!is.null(seed)) {
    set.seed(seed)
    bike_sharing <- bike_sharing[sample(1:nrow(bike_sharing)),]
  }
  X <- subset(bike_sharing, select = -c(rented))
  design_mat <- as.matrix(cbind(rep(1, nrow(X)), X))
  colnames(design_mat)[1] <- 'intercept'
  return(list('data' = cbind(subset(design_mat, select = -c(intercept)),
                             'rented' = bike_sharing$rented),
              'y' = bike_sharing$rented,
              'X' = design_mat))
}

bike_sharing <- load_bs_data('scripts/count_data_regression/bike_sharing/hour.csv', seed = seed)

# MASS::glm.nb(rented~., data = as.data.frame(bike_sharing$data))

##### Sampling from full posterior #####

full_posterior <- hmc_sample_GLMR(likelihood = 'NB',
                                  y = bike_sharing$y,
                                  X = bike_sharing$X,
                                  C = 1,
                                  phi = phi_rate,
                                  prior_means = prior_means,
                                  prior_variances = prior_variances,
                                  iterations = nsamples + 10000,
                                  warmup = 10000,
                                  chains = 1,
                                  seed = seed,
                                  output = T)

apply(full_posterior, 2, mean)

##### Sampling from sub-posterior C=32 #####

data_split_32 <- split_data(bike_sharing$data,
                            y_col_index = 10,
                            X_col_index = 1:9,
                            C = 32,
                            as_dataframe = F)
sub_posteriors_32 <- hmc_base_sampler_GLMR(likelihood = 'NB',
                                           nsamples = nsamples,
                                           warmup = 10000,
                                           data_split = data_split_32,
                                           C = 32,
                                           phi = phi_rate,
                                           prior_means = prior_means,
                                           prior_variances = prior_variances,
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

##### bal binary combining two sub-posteriors at a time #####
balanced_C32 <- list('reg' = bal_binary_GBF_BNBR(N_schedule = rep(nsamples, 5),
                                                 m_schedule = rep(2, 5),
                                                 time_mesh = NULL,
                                                 base_samples = sub_posteriors_32,
                                                 L = 6,
                                                 dim = 10,
                                                 phi_rate = phi_rate,
                                                 data_split = data_split_32,
                                                 prior_means = prior_means,
                                                 prior_variances = prior_variances,
                                                 C = 32,
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
                                                 print_progress_iters = 100))
balanced_C32$adaptive <- bal_binary_GBF_BNBR(N_schedule = rep(nsamples, 5),
                                             m_schedule = rep(2, 5),
                                             time_mesh = NULL,
                                             base_samples = sub_posteriors_32,
                                             L = 6,
                                             dim = 10,
                                             phi_rate = phi_rate,
                                             data_split = data_split_32,
                                             prior_means = prior_means,
                                             prior_variances = prior_variances,
                                             C = 32,
                                             precondition = TRUE,
                                             resampling_method = 'resid',
                                             ESS_threshold = ESS_threshold,
                                             adaptive_mesh = TRUE,
                                             mesh_parameters = list('condition' = 'SH',
                                                                    'CESS_0_threshold' = CESS_0_threshold,
                                                                    'CESS_j_threshold' = CESS_j_threshold,
                                                                    'vanilla' = FALSE),
                                             diffusion_estimator = diffusion_estimator,
                                             seed = seed,
                                             print_progress_iters = 100)

# regular mesh
balanced_C32$reg$particles <- resample_particle_y_samples(particle_set = balanced_C32$reg$particles[[1]],
                                                          multivariate = TRUE,
                                                          resampling_method = 'resid',
                                                          seed = seed)
balanced_C32$reg$proposed_samples <- balanced_C32$reg$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, balanced_C32$reg$particles$y_samples))
# adaptive mesh
balanced_C32$adaptive$particles <- resample_particle_y_samples(particle_set = balanced_C32$adaptive$particles[[1]],
                                                               multivariate = TRUE,
                                                               resampling_method = 'resid',
                                                               seed = seed)
balanced_C32$reg$proposed_samples <- balanced_C32$adaptive$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, balanced_C32$adaptive$particles$y_samples))

##### IAD #####

integrated_abs_distance(full_posterior, balanced_C32$reg$particles$y_samples)
integrated_abs_distance(full_posterior, balanced_C32$adaptive$particles$y_samples)
integrated_abs_distance(full_posterior, consensus_mat_32$samples)
integrated_abs_distance(full_posterior, consensus_sca_32$samples)
integrated_abs_distance(full_posterior, neiswanger_true_32$samples)
integrated_abs_distance(full_posterior, neiswanger_false_32$samples)
integrated_abs_distance(full_posterior, weierstrass_importance_32$samples)
integrated_abs_distance(full_posterior, weierstrass_rejection_32$samples)

save.image('BS32_NB.RData')
