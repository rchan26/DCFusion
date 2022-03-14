library(DCFusion)
library(HMCBLR)

##### Initialise example #####
seed <- 2021
set.seed(seed)
nsamples <- 10000
ndata <- 1000
C <- 16
n_cores <- parallel::detectCores()
true_beta <- c(-3, 1.2, -0.5, 0.8, 3)
frequencies <- c(0.2, 0.3, 0.5, 0.01)
diffusion_estimator <- 'NB'
ESS_threshold <- 0.5
CESS_0_threshold <- 0.2
CESS_j_threshold <- 0.1
k1 <- NULL
k2 <- NULL
k3 <- -log(CESS_j_threshold)/2
k4 <- -log(CESS_j_threshold)/2

# simulate data set
simulated_data <- simulate_LR_data(N = ndata,
                                   alpha = true_beta[1],
                                   frequencies = frequencies,
                                   coefficients = true_beta[2:length(true_beta)],
                                   seed = seed)

# check activity of the parameters
check_activity(simulated_data)

##### Sampling from full posterior #####t

full_data_count <- unique_row_count(y = simulated_data[,1],
                                    X = cbind('intercept' = rep(1, ndata), simulated_data[,2:ncol(simulated_data)]))$full_data_count
full_posterior <- hmc_sample_BLR(full_data_count = full_data_count,
                                 C = 1,
                                 prior_means = rep(0, 5),
                                 prior_variances = rep(1, 5),
                                 iterations = nsamples + 10000,
                                 warmup = 10000,
                                 chains = 1,
                                 seed = seed,
                                 output = T)

##### Sampling from sub-posterior C=16 #####

data_split_16 <- split_data(simulated_data, y_col_index = 1, X_col_index = 2:ncol(simulated_data), C = 16, as_dataframe = F)
sub_posteriors_16 <- hmc_base_sampler_BLR(nsamples = nsamples,
                                          data_split = data_split_16,
                                          C = 16, 
                                          prior_means = rep(0, 5),
                                          prior_variances = rep(1, 5),
                                          warmup = 10000,
                                          seed = seed,
                                          output = T)

##### Applying other methodologies #####

print('Applying other methodologies')
consensus_mat_16 <- consensus_scott(S = 16, samples_to_combine = sub_posteriors_16, indep = F)
consensus_sca_16 <- consensus_scott(S = 16, samples_to_combine = sub_posteriors_16, indep = T)
neiswanger_true_16 <- neiswanger(S = 16,
                                 samples_to_combine = sub_posteriors_16,
                                 anneal = TRUE)
neiswanger_false_16 <- neiswanger(S = 16,
                                  samples_to_combine = sub_posteriors_16,
                                  anneal = FALSE)
weierstrass_importance_16 <- weierstrass(Samples = sub_posteriors_16,
                                         method = 'importance')
weierstrass_rejection_16 <- weierstrass(Samples = sub_posteriors_16,
                                        method = 'reject')

##### all at once #####
GBF_16 <- list('reg' = bal_binary_GBF_BLR(N_schedule = nsamples,
                                          m_schedule = 16,
                                          time_mesh = NULL,
                                          base_samples = sub_posteriors_16,
                                          L = 2,
                                          dim = 5,
                                          data_split = data_split_16,
                                          prior_means = rep(0, 5),
                                          prior_variances = rep(1, 5),
                                          C = C,
                                          precondition = TRUE,
                                          resampling_method = 'resid',
                                          ESS_threshold = ESS_threshold,
                                          adaptive_mesh = FALSE,
                                          mesh_parameters = list('condition' = 'SH',
                                                                 'CESS_0_threshold' = CESS_0_threshold,
                                                                 'CESS_j_threshold' = CESS_j_threshold,
                                                                 'k1' = k1,
                                                                 'k2' = k2,
                                                                 'k3' = k3,
                                                                 'k4' = k4,
                                                                 'vanilla' = FALSE),
                                          diffusion_estimator = diffusion_estimator,
                                          seed = seed),
               'adaptive' = bal_binary_GBF_BLR(N_schedule = nsamples,
                                               m_schedule = 16,
                                               time_mesh = NULL,
                                               base_samples = sub_posteriors_16,
                                               L = 2,
                                               dim = 5,
                                               data_split = data_split_16,
                                               prior_means = rep(0, 5),
                                               prior_variances = rep(1, 5),
                                               C = C,
                                               precondition = TRUE,
                                               resampling_method = 'resid',
                                               ESS_threshold = ESS_threshold,
                                               adaptive_mesh = TRUE,
                                               mesh_parameters = list('condition' = 'SH',
                                                                      'CESS_0_threshold' = CESS_0_threshold,
                                                                      'CESS_j_threshold' = CESS_j_threshold,
                                                                      'k1' = k1,
                                                                      'k2' = k2,
                                                                      'k3' = k3,
                                                                      'k4' = k4,
                                                                      'vanilla' = FALSE),
                                               diffusion_estimator = diffusion_estimator,
                                               seed = seed))

# regular mesh
GBF_16$reg$particles <- resample_particle_y_samples(particle_set = GBF_16$reg$particles[[1]],
                                                    multivariate = TRUE,
                                                    resampling_method = 'resid',
                                                    seed = seed)
print(integrated_abs_distance(full_posterior, GBF_16$reg$particles$y_samples))
compare_samples_bivariate(posteriors = list(full_posterior,
                                            GBF_16$reg$proposed_samples[[1]],
                                            GBF_16$reg$particles$y_samples),
                          colours = c('black', 'green', 'red'),
                          common_limit = c(-4, 4))
# adaptive mesh
GBF_16$adaptive$particles <- resample_particle_y_samples(particle_set = GBF_16$adaptive$particles[[1]],
                                                         multivariate = TRUE,
                                                         resampling_method = 'resid',
                                                         seed = seed)
print(integrated_abs_distance(full_posterior, GBF_16$adaptive$particles$y_samples))
compare_samples_bivariate(posteriors = list(full_posterior,
                                            GBF_16$adaptive$proposed_samples[[1]],
                                            GBF_16$adaptive$particles$y_samples),
                          colours = c('black', 'green', 'red'),
                          common_limit = c(-4, 4))

##### bal binary combining two sub-posteriors at a time #####
balanced_C16 <- list('reg' = bal_binary_GBF_BLR(N_schedule = rep(nsamples, 4),
                                                m_schedule = rep(2, 4),
                                                time_mesh = NULL,
                                                base_samples = sub_posteriors_16,
                                                L = 5,
                                                dim = 5,
                                                data_split = data_split_16,
                                                prior_means = rep(0, 5),
                                                prior_variances = rep(1, 5),
                                                C = C,
                                                precondition = TRUE,
                                                resampling_method = 'resid',
                                                ESS_threshold = ESS_threshold,
                                                adaptive_mesh = FALSE,
                                                mesh_parameters = list('condition' = 'SH',
                                                                       'CESS_0_threshold' = CESS_0_threshold,
                                                                       'CESS_j_threshold' = CESS_j_threshold,
                                                                       'k1' = k1,
                                                                       'k2' = k2,
                                                                       'k3' = k3,
                                                                       'k4' = k4,
                                                                       'vanilla' = FALSE),
                                                diffusion_estimator = diffusion_estimator,
                                                seed = seed),
                     'adaptive' = bal_binary_GBF_BLR(N_schedule = rep(nsamples, 4),
                                                     m_schedule = rep(2, 4),
                                                     time_mesh = NULL,
                                                     base_samples = sub_posteriors_16,
                                                     L = 5,
                                                     dim = 5,
                                                     data_split = data_split_16,
                                                     prior_means = rep(0, 5),
                                                     prior_variances = rep(1, 5),
                                                     C = C,
                                                     precondition = TRUE,
                                                     resampling_method = 'resid',
                                                     ESS_threshold = ESS_threshold,
                                                     adaptive_mesh = TRUE,
                                                     mesh_parameters = list('condition' = 'SH',
                                                                            'CESS_0_threshold' = CESS_0_threshold,
                                                                            'CESS_j_threshold' = CESS_j_threshold,
                                                                            'k1' = k1,
                                                                            'k2' = k2,
                                                                            'k3' = k3,
                                                                            'k4' = k4,
                                                                            'vanilla' = FALSE),
                                                     diffusion_estimator = diffusion_estimator,
                                                     seed = seed))

# regular mesh
balanced_C16$reg$particles <- resample_particle_y_samples(particle_set = balanced_C16$reg$particles[[1]],
                                                          multivariate = TRUE,
                                                          resampling_method = 'resid',
                                                          seed = seed)
print(integrated_abs_distance(full_posterior, balanced_C16$reg$particles$y_samples))
compare_samples_bivariate(posteriors = list(full_posterior,
                                            balanced_C16$reg$proposed_samples[[1]],
                                            balanced_C16$reg$particles$y_samples),
                          colours = c('black', 'green', 'red'),
                          common_limit = c(-4, 4))
# adaptive mesh
balanced_C16$adaptive$particles <- resample_particle_y_samples(particle_set = balanced_C16$adaptive$particles[[1]],
                                                               multivariate = TRUE,
                                                               resampling_method = 'resid',
                                                               seed = seed)
print(integrated_abs_distance(full_posterior, balanced_C16$adaptive$particles$y_samples))
compare_samples_bivariate(posteriors = list(full_posterior,
                                            balanced_C16$adaptive$proposed_samples[[1]],
                                            balanced_C16$adaptive$particles$y_samples),
                          colours = c('black', 'green', 'red'),
                          common_limit = c(-4, 4))

save.image('SD16_DCGBF.RData')
