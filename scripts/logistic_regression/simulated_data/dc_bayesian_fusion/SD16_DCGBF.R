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
C <- 16
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

##### Sampling from sub-posterior C=16 #####

data_split_16 <- split_data(simulated_data, y_col_index = 1, X_col_index = 2:ncol(simulated_data), C = 16, as_dataframe = F)
sub_posteriors_16 <- hmc_base_sampler_BLR(nsamples = nsamples,
                                          data_split = data_split_16,
                                          C = 16, 
                                          prior_means = prior_means,
                                          prior_variances = prior_variances,
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

integrated_abs_distance(full_posterior, consensus_mat_16$samples)
integrated_abs_distance(full_posterior, consensus_sca_16$samples)
integrated_abs_distance(full_posterior, neiswanger_true_16$samples)
integrated_abs_distance(full_posterior, neiswanger_false_16$samples)
integrated_abs_distance(full_posterior, weierstrass_importance_16$samples)
integrated_abs_distance(full_posterior, weierstrass_rejection_16$samples)

##### NB (Hypercube Centre) #####
print('NB Fusion (hypercube centre)')
NB_hc_16 <- bal_binary_fusion_SMC_BLR(N_schedule = rep(nsamples, 4),
                                      m_schedule = rep(2, 4),
                                      time_schedule = rep(time_choice, 4),
                                      base_samples = sub_posteriors_16,
                                      L = 5,
                                      dim = 5,
                                      data_split = data_split_16,
                                      prior_means = prior_means,
                                      prior_variances = prior_variances,
                                      C = 16,
                                      precondition = TRUE,
                                      resampling_method = 'resid',
                                      ESS_threshold = ESS_threshold,
                                      cv_location = 'hypercube_centre',
                                      diffusion_estimator = 'NB',
                                      seed = seed,
                                      n_cores = n_cores)
NB_hc_16$particles <- resample_particle_y_samples(particle_set = NB_hc_16$particles[[1]],
                                                  multivariate = TRUE,
                                                  resampling_method = 'resid',
                                                  seed = seed)
NB_hc_16$proposed_samples <- NB_hc_16$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, NB_hc_16$particles$y_samples))

##### Generalised Bayesian Fusion #####

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
                                                                 'vanilla' = FALSE),
                                          diffusion_estimator = diffusion_estimator,
                                          seed = seed))
GBF_16$adaptive <- bal_binary_GBF_BLR(N_schedule = nsamples,
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
                                                             'vanilla' = FALSE),
                                      diffusion_estimator = diffusion_estimator,
                                      seed = seed)

# regular mesh
GBF_16$reg$particles <- resample_particle_y_samples(particle_set = GBF_16$reg$particles[[1]],
                                                    multivariate = TRUE,
                                                    resampling_method = 'resid',
                                                    seed = seed)
GBF_16$reg$proposed_samples <- GBF_16$reg$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, GBF_16$reg$particles$y_samples))
# adaptive mesh
GBF_16$adaptive$particles <- resample_particle_y_samples(particle_set = GBF_16$adaptive$particles[[1]],
                                                         multivariate = TRUE,
                                                         resampling_method = 'resid',
                                                         seed = seed)
GBF_16$adaptive$proposed_samples <- GBF_16$adaptive$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, GBF_16$adaptive$particles$y_samples))

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
                                                                       'vanilla' = FALSE),
                                                diffusion_estimator = diffusion_estimator,
                                                seed = seed))
balanced_C16$adaptive <- bal_binary_GBF_BLR(N_schedule = rep(nsamples, 4),
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
                                                                   'vanilla' = FALSE),
                                            diffusion_estimator = diffusion_estimator,
                                            seed = seed)

# regular mesh
balanced_C16$reg$particles <- resample_particle_y_samples(particle_set = balanced_C16$reg$particles[[1]],
                                                          multivariate = TRUE,
                                                          resampling_method = 'resid',
                                                          seed = seed)
balanced_C16$reg$proposed_samples <- balanced_C16$reg$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, balanced_C16$reg$particles$y_samples))
# adaptive mesh
balanced_C16$adaptive$particles <- resample_particle_y_samples(particle_set = balanced_C16$adaptive$particles[[1]],
                                                               multivariate = TRUE,
                                                               resampling_method = 'resid',
                                                               seed = seed)
balanced_C16$adaptive$proposed_samples <- balanced_C16$adaptive$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, balanced_C16$adaptive$particles$y_samples))

##### IAD #####

integrated_abs_distance(full_posterior, GBF_16$reg$particles$y_samples)
integrated_abs_distance(full_posterior, GBF_16$adaptive$particles$y_samples)
integrated_abs_distance(full_posterior, balanced_C16$reg$particles$y_samples)
integrated_abs_distance(full_posterior, balanced_C16$adaptive$particles$y_samples)
integrated_abs_distance(full_posterior, NB_hc_16$particles$y_samples)

save.image('SD16_DCGBF.RData')
