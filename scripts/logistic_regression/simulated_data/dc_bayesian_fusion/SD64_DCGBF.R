library(DCFusion)
library(HMCBLR)

##### Initialise example #####
seed <- 2022
set.seed(seed)
nsamples <- 10000
ndata <- 1000
time_choice <- 0.5
prior_means <- rep(0, 5)
prior_variances <- rep(1, 5)
C <- 64
n_cores <- parallel::detectCores()
true_beta <- c(-3, 1.2, -0.5, 0.8, 3)
frequencies <- c(0.2, 0.3, 0.5, 0.01)
ESS_threshold <- 0.5
CESS_0_threshold <- 0.5
CESS_j_threshold <- 0.2
diffusion_estimator <- 'NB'

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
                                 prior_means = prior_means,
                                 prior_variances = prior_variances,
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
                                          prior_means = prior_means,
                                          prior_variances = prior_variances,
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

integrated_abs_distance(full_posterior, consensus_mat_64$samples)
integrated_abs_distance(full_posterior, consensus_sca_64$samples)
integrated_abs_distance(full_posterior, neiswanger_true_64$samples)
integrated_abs_distance(full_posterior, neiswanger_false_64$samples)
integrated_abs_distance(full_posterior, weierstrass_importance_64$samples)
integrated_abs_distance(full_posterior, weierstrass_rejection_64$samples)

##### Generalised Bayesian Fusion #####

##### bal binary combining two sub-posteriors at a time #####
balanced_C64 <- list('reg' = bal_binary_GBF_BLR(N_schedule = rep(nsamples, 6),
                                                m_schedule = rep(2, 6),
                                                time_mesh = NULL,
                                                base_samples = sub_posteriors_64,
                                                L = 7,
                                                dim = 5,
                                                data_split = data_split_64,
                                                prior_means = prior_means,
                                                prior_variances = prior_variances,
                                                C = 64,
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
balanced_C64$adaptive <- bal_binary_GBF_BLR(N_schedule = rep(nsamples, 6),
                                            m_schedule = rep(2, 6),
                                            time_mesh = NULL,
                                            base_samples = sub_posteriors_64,
                                            L = 7,
                                            dim = 5,
                                            data_split = data_split_64,
                                            prior_means = prior_means,
                                            prior_variances = prior_variances,
                                            C = 64,
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
balanced_C64$reg$particles <- resample_particle_y_samples(particle_set = balanced_C64$reg$particles[[1]],
                                                          multivariate = TRUE,
                                                          resampling_method = 'resid',
                                                          seed = seed)
balanced_C64$reg$proposed_samples <- balanced_C64$reg$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, balanced_C64$reg$particles$y_samples))
# adaptive mesh
balanced_C64$adaptive$particles <- resample_particle_y_samples(particle_set = balanced_C64$adaptive$particles[[1]],
                                                               multivariate = TRUE,
                                                               resampling_method = 'resid',
                                                               seed = seed)
balanced_C64$adaptive$proposed_samples <- balanced_C64$adaptive$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, balanced_C64$adaptive$particles$y_samples))

##### IAD #####

integrated_abs_distance(full_posterior, balanced_C64$reg$particles$y_samples)
integrated_abs_distance(full_posterior, balanced_C64$adaptive$particles$y_samples)

save.image('SD64_DCGBF.RData')
