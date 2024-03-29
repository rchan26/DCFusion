library(DCFusion)
library(HMCBLR)

##### Initialise example #####
seed <- 2021
set.seed(seed)
nsamples <- 10000
ndata <- 1000
time_choice <- 0.5
n_cores <- parallel::detectCores()
true_beta <- c(-3, 1.2, -0.5, 0.8, 3)
frequencies <- c(0.2, 0.3, 0.5, 0.01) # must have length = length(true_beta)-1

# simulate data set
simulated_data <- simulate_LR_data(N = ndata,
                                   alpha = true_beta[1],
                                   frequencies = frequencies,
                                   coefficients = true_beta[2:length(true_beta)],
                                   seed = seed) 

# check activity of the parameters
check_activity(simulated_data)

##### Sampling from full posterior #####

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

##### Sampling from sub-posterior C=4 #####

data_split_4 <- split_data(simulated_data, y_col_index = 1, X_col_index = 2:ncol(simulated_data), C = 4, as_dataframe = F)
sub_posteriors_4 <- hmc_base_sampler_BLR(nsamples = nsamples,
                                         data_split = data_split_4,
                                         C = 4, 
                                         prior_means = rep(0, 5),
                                         prior_variances = rep(1, 5),
                                         warmup = 10000,
                                         seed = seed,
                                         output = T)

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

##### Applying Fusion #####

##### Poisson (Hypercube Centre) #####
print('Poisson Fusion (hypercube centre)')
Poisson_hc_4 <- bal_binary_fusion_SMC_BLR(N_schedule = rep(nsamples, 2),
                                          m_schedule = rep(2, 2),
                                          time_schedule = rep(time_choice, 2),
                                          base_samples = sub_posteriors_4,
                                          L = 3,
                                          dim = 5,
                                          data_split = data_split_4,
                                          prior_means = rep(0, 5),
                                          prior_variances = rep(1, 5),
                                          C = 4,
                                          precondition = TRUE,
                                          resampling_method = 'resid',
                                          ESS_threshold = 0.5,
                                          cv_location = 'hypercube_centre',
                                          diffusion_estimator = 'Poisson',
                                          seed = seed,
                                          n_cores = n_cores,
                                          print_progress_iters = 500)
Poisson_hc_4$particles <- resample_particle_y_samples(particle_set = Poisson_hc_4$particles[[1]],
                                                      multivariate = TRUE,
                                                      resampling_method = 'resid',
                                                      seed = seed)
Poisson_hc_4$proposed_samples <- Poisson_hc_4$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, Poisson_hc_4$particles$y_samples))

##### NB (Hypercube Centre) #####
print('NB Fusion (hypercube centre)')
NB_hc_4 <- bal_binary_fusion_SMC_BLR(N_schedule = rep(nsamples, 2),
                                     m_schedule = rep(2, 2),
                                     time_schedule = rep(time_choice, 2),
                                     base_samples = sub_posteriors_4,
                                     L = 3,
                                     dim = 5,
                                     data_split = data_split_4,
                                     prior_means = rep(0, 5),
                                     prior_variances = rep(1, 5),
                                     C = 4,
                                     precondition = TRUE,
                                     resampling_method = 'resid',
                                     ESS_threshold = 0.5,
                                     cv_location = 'hypercube_centre',
                                     diffusion_estimator = 'NB',
                                     seed = seed,
                                     n_cores = n_cores,
                                     print_progress_iters = 500)
NB_hc_4$particles <- resample_particle_y_samples(particle_set = NB_hc_4$particles[[1]],
                                                 multivariate = TRUE,
                                                 resampling_method = 'resid',
                                                 seed = seed)
NB_hc_4$proposed_samples <- NB_hc_4$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, NB_hc_4$particles$y_samples))

##### IAD #####

integrated_abs_distance(full_posterior, Poisson_hc_4$particles$y_samples)
integrated_abs_distance(full_posterior, NB_hc_4$particles$y_samples)
integrated_abs_distance(full_posterior, consensus_mat_4$samples)
integrated_abs_distance(full_posterior, consensus_sca_4$samples)
integrated_abs_distance(full_posterior, neiswanger_true_4$samples)
integrated_abs_distance(full_posterior, neiswanger_false_4$samples)
integrated_abs_distance(full_posterior, weierstrass_importance_4$samples)
integrated_abs_distance(full_posterior, weierstrass_rejection_4$samples)

##### Save data #####

save.image('SD4.RData')
