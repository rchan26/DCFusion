library(DCFusion)
library(HMCBLR)

##### Initialise example #####
seed <- 2021
set.seed(seed)
nsamples <- 10000
ndata <- 5000
time_choice <- 0.5
n_cores <- parallel::detectCores()
ESS_threshold <- 0.5
resampling_method <- 'resid'
dim <- 20
prior_means <- rep(0, dim)
prior_variances <- rep(1, dim)
true_beta <- round(rnorm(dim, 0, 1), 2)
frequencies <- round(runif(dim-1, 0, 1), 2)

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
                                 prior_means = prior_means,
                                 prior_variances = prior_variances,
                                 iterations = nsamples + 10000,
                                 warmup = 10000,
                                 chains = 1,
                                 seed = seed,
                                 output = T)

apply(full_posterior, 2, mean)
true_beta
abs(apply(full_posterior, 2, mean)-true_beta)

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

##### Applying Fusion #####

##### Poisson (Hypercube Centre) #####
print('Poisson Fusion (hypercube centre)')
Poisson_hc_16 <- bal_binary_fusion_SMC_BLR(N_schedule = rep(nsamples, 4),
                                           m_schedule = rep(2, 4),
                                           time_schedule = rep(time_choice, 4),
                                           base_samples = sub_posteriors_16,
                                           L = 5,
                                           dim = dim,
                                           data_split = data_split_16,
                                           prior_means = prior_means,
                                           prior_variances = prior_variances,
                                           C = 16,
                                           precondition = TRUE,
                                           resampling_method = resampling_method,
                                           ESS_threshold = ESS_threshold,
                                           cv_location = 'hypercube_centre',
                                           diffusion_estimator = 'Poisson',
                                           seed = seed,
                                           n_cores = n_cores,
                                           print_progress_iters = 10)
Poisson_hc_16$particles <- resample_particle_y_samples(particle_set = Poisson_hc_16$particles[[1]],
                                                       multivariate = TRUE,
                                                       resampling_method = resampling_method,
                                                       seed = seed)
Poisson_hc_16$proposed_samples <- Poisson_hc_16$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, Poisson_hc_16$particles$y_samples))

##### NB (Hypercube Centre) #####
print('NB Fusion (hypercube centre)')
NB_hc_16 <- bal_binary_fusion_SMC_BLR(N_schedule = rep(nsamples, 4),
                                      m_schedule = rep(2, 4),
                                      time_schedule = rep(time_choice, 4),
                                      base_samples = sub_posteriors_16,
                                      L = 5,
                                      dim = dim,
                                      data_split = data_split_16,
                                      prior_means = prior_means,
                                      prior_variances = prior_variances,
                                      C = 16,
                                      precondition = TRUE,
                                      resampling_method = resampling_method,
                                      ESS_threshold = ESS_threshold,
                                      cv_location = 'hypercube_centre',
                                      diffusion_estimator = 'NB',
                                      seed = seed,
                                      n_cores = n_cores,
                                      print_progress_iters = 10)
NB_hc_16$particles <- resample_particle_y_samples(particle_set = NB_hc_16$particles[[1]],
                                                  multivariate = TRUE,
                                                  resampling_method = resampling_method,
                                                  seed = seed)
NB_hc_16$proposed_samples <- NB_hc_16$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, NB_hc_16$particles$y_samples))

##### IAD #####

integrated_abs_distance(full_posterior, Poisson_hc_16$particles$y_samples)
integrated_abs_distance(full_posterior, NB_hc_16$particles$y_samples)
integrated_abs_distance(full_posterior, consensus_mat_16$samples)
integrated_abs_distance(full_posterior, consensus_sca_16$samples)
integrated_abs_distance(full_posterior, neiswanger_true_16$samples)
integrated_abs_distance(full_posterior, neiswanger_false_16$samples)
integrated_abs_distance(full_posterior, weierstrass_importance_16$samples)
integrated_abs_distance(full_posterior, weierstrass_rejection_16$samples)

##### Save data #####

save.image('SD16_d20.RData')
