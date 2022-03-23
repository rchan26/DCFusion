library(DCFusion)
library(HMCBLR)

##### Initialise example #####
seed <- 2021
set.seed(seed)
nsamples <- 10000
ndata <- 1000
C <- 64
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

##### Sampling from sub-posterior C=64 #####

data_split_64 <- split_data(simulated_data, y_col_index = 1, X_col_index = 2:ncol(simulated_data), C = 64, as_dataframe = F)
sub_posteriors_64 <- hmc_base_sampler_BLR(nsamples = nsamples,
                                          data_split = data_split_64,
                                          C = 64, 
                                          prior_means = rep(0, 5),
                                          prior_variances = rep(1, 5),
                                          warmup = 10000,
                                          seed = seed,
                                          output = T)

##### Applying other methodologies #####

print('Applying other methodologies')
consensus_mat_64 <- consensus_scott(S = 64, samples_to_combine = sub_posteriors_64, indep = F)
consensus_sca_64 <- consensus_scott(S = 64, samples_to_combine = sub_posteriors_64, indep = T)
neiswanger_true_64 <- neiswanger(S = 64,
                                 samples_to_combine = sub_posteriors_64,
                                 anneal = TRUE)
neiswanger_false_64 <- neiswanger(S = 64,
                                  samples_to_combine = sub_posteriors_64,
                                  anneal = FALSE)
weierstrass_importance_64 <- weierstrass(Samples = sub_posteriors_64,
                                         method = 'importance')
weierstrass_rejection_64 <- weierstrass(Samples = sub_posteriors_64,
                                        method = 'reject')

##### all at once #####
GBF_64 <- list('reg' = bal_binary_GBF_BLR(N_schedule = nsamples,
                                          m_schedule = 64,
                                          time_mesh = NULL,
                                          base_samples = sub_posteriors_64,
                                          L = 2,
                                          dim = 5,
                                          data_split = data_split_64,
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
                                               m_schedule = 64,
                                               time_mesh = NULL,
                                               base_samples = sub_posteriors_64,
                                               L = 2,
                                               dim = 5,
                                               data_split = data_split_64,
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
GBF_64$reg$particles <- resample_particle_y_samples(particle_set = GBF_64$reg$particles[[1]],
                                                    multivariate = TRUE,
                                                    resampling_method = 'resid',
                                                    seed = seed)
print(integrated_abs_distance(full_posterior, GBF_64$reg$particles$y_samples))
compare_samples_bivariate(posteriors = list(full_posterior,
                                            GBF_64$reg$proposed_samples[[1]],
                                            GBF_64$reg$particles$y_samples),
                          colours = c('black', 'green', 'red'),
                          common_limit = c(-4, 4))
# adaptive mesh
GBF_64$adaptive$particles <- resample_particle_y_samples(particle_set = GBF_64$adaptive$particles[[1]],
                                                         multivariate = TRUE,
                                                         resampling_method = 'resid',
                                                         seed = seed)
print(integrated_abs_distance(full_posterior, GBF_64$adaptive$particles$y_samples))
compare_samples_bivariate(posteriors = list(full_posterior,
                                            GBF_64$adaptive$proposed_samples[[1]],
                                            GBF_64$adaptive$particles$y_samples),
                          colours = c('black', 'green', 'red'),
                          common_limit = c(-4, 4))

##### bal binary combining two sub-posteriors at a time #####
balanced_C64 <- list('reg' = bal_binary_GBF_BLR(N_schedule = rep(nsamples, 6),
                                                m_schedule = rep(2, 6),
                                                time_mesh = NULL,
                                                base_samples = sub_posteriors_64,
                                                L = 7,
                                                dim = 5,
                                                data_split = data_split_64,
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
                     'adaptive' = bal_binary_GBF_BLR(N_schedule = rep(nsamples, 6),
                                                     m_schedule = rep(2, 6),
                                                     time_mesh = NULL,
                                                     base_samples = sub_posteriors_64,
                                                     L = 7,
                                                     dim = 5,
                                                     data_split = data_split_64,
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
balanced_C64$reg$particles <- resample_particle_y_samples(particle_set = balanced_C64$reg$particles[[1]],
                                                          multivariate = TRUE,
                                                          resampling_method = 'resid',
                                                          seed = seed)
print(integrated_abs_distance(full_posterior, balanced_C64$reg$particles$y_samples))
compare_samples_bivariate(posteriors = list(full_posterior,
                                            balanced_C64$reg$proposed_samples[[1]],
                                            balanced_C64$reg$particles$y_samples),
                          colours = c('black', 'green', 'red'),
                          common_limit = c(-4, 4))
# adaptive mesh
balanced_C64$adaptive$particles <- resample_particle_y_samples(particle_set = balanced_C64$adaptive$particles[[1]],
                                                               multivariate = TRUE,
                                                               resampling_method = 'resid',
                                                               seed = seed)
print(integrated_abs_distance(full_posterior, balanced_C64$adaptive$particles$y_samples))
compare_samples_bivariate(posteriors = list(full_posterior,
                                            balanced_C64$adaptive$proposed_samples[[1]],
                                            balanced_C64$adaptive$particles$y_samples),
                          colours = c('black', 'green', 'red'),
                          common_limit = c(-4, 4))

save.image('SD64_DCGBF.RData')