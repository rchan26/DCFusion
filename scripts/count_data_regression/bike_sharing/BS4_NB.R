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

bike_sharing <- load_bs_data('scripts/count_data_regression//bike_sharing/hour.csv', seed = seed)

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

##### Sampling from sub-posterior C=4 #####

data_split_4 <- split_data(bike_sharing$data,
                           y_col_index = 10,
                           X_col_index = 1:9,
                           C = 4,
                           as_dataframe = F)
sub_posteriors_4 <- hmc_base_sampler_GLMR(likelihood = 'NB',
                                          nsamples = nsamples,
                                          warmup = 10000,
                                          data_split = data_split_4,
                                          C = 4,
                                          phi = phi_rate,
                                          prior_means = prior_means,
                                          prior_variances = prior_variances,
                                          seed = seed,
                                          output = T)

for (i in 1:4) {
  print(apply(sub_posteriors_4[[i]], 2, mean))
}
# compare_samples_bivariate(sub_posteriors_4, c('red', 'red', 'blue', 'blue'), c(-5, 5))
# integrated_abs_distance(sub_posteriors_4[[1]], sub_posteriors_4[[2]])
# integrated_abs_distance(sub_posteriors_4[[3]], sub_posteriors_4[[4]])

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

##### bal binary combining two sub-posteriors at a time #####
balanced_C4 <- list('reg' = bal_binary_GBF_BNBR(N_schedule = rep(nsamples, 2),
                                                m_schedule = rep(2, 2),
                                                time_mesh = NULL,
                                                base_samples = sub_posteriors_4,
                                                L = 3,
                                                dim = 10,
                                                phi_rate = phi_rate,
                                                data_split = data_split_4,
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
                                                seed = seed,
                                                print_progress_iters = 100))
balanced_C4$adaptive <- bal_binary_GBF_BNBR(N_schedule = rep(nsamples, 2),
                                            m_schedule = rep(2, 2),
                                            time_mesh = NULL,
                                            base_samples = sub_posteriors_4,
                                            L = 3,
                                            dim = 10,
                                            phi_rate = phi_rate,
                                            data_split = data_split_4,
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
                                            seed = seed,
                                            print_progress_iters = 100)

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
balanced_C4$reg$proposed_samples <- balanced_C4$adaptive$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, balanced_C4$adaptive$particles$y_samples))

##### IAD #####

integrated_abs_distance(full_posterior, balanced_C4$reg$particles$y_samples)
integrated_abs_distance(full_posterior, balanced_C4$adaptive$particles$y_samples)
integrated_abs_distance(full_posterior, consensus_mat_4$samples)
integrated_abs_distance(full_posterior, consensus_sca_4$samples)
integrated_abs_distance(full_posterior, neiswanger_true_4$samples)
integrated_abs_distance(full_posterior, neiswanger_false_4$samples)
integrated_abs_distance(full_posterior, weierstrass_importance_4$samples)
integrated_abs_distance(full_posterior, weierstrass_rejection_4$samples)

save.image('BS4_NB.RData')
