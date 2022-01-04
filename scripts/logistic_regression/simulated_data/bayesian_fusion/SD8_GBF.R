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

##### Sampling from sub-posterior C=8 #####

data_split_8 <- split_data(simulated_data, y_col_index = 1, X_col_index = 2:ncol(simulated_data), C = 8, as_dataframe = F)
sub_posteriors_8 <- hmc_base_sampler_BLR(nsamples = nsamples,
                                         data_split = data_split_8,
                                         C = 8, 
                                         prior_means = rep(0, 5),
                                         prior_variances = rep(1, 5),
                                         warmup = 10000,
                                         seed = seed,
                                         output = T)

##### Applying Vanilla Bayesian Fusion (by passing in identity matrices) #####

time_choice <- 2
time_mesh <- seq(0, time_choice, 0.005)
particles_to_fuse <- initialise_particle_sets(samples_to_fuse = sub_posteriors_8,
                                              multivariate = TRUE,
                                              number_of_steps = length(time_mesh))

##### Poisson (Hypercube Centre) #####
print('Poisson Fusion (hypercube centre)')
Poisson_VBF_8 <- parallel_generalised_BF_SMC_BLR(particles_to_fuse = particles_to_fuse,
                                                 N = nsamples,
                                                 m = 8,
                                                 time_mesh = time_mesh,
                                                 dim = 5,
                                                 data_split = data_split_8,
                                                 prior_means = rep(0, 5),
                                                 prior_variances = rep(1, 5),
                                                 C = 8,
                                                 precondition_matrices = rep(list(diag(1, 5)), 8),
                                                 resampling_method = 'resid',
                                                 ESS_threshold = 0.5,
                                                 cv_location = 'hypercube_centre',
                                                 diffusion_estimator = 'Poisson',
                                                 seed = seed,
                                                 n_cores = n_cores,
                                                 print_progress_iters = 500)
Poisson_VBF_8$particles <- resample_particle_y_samples(particle_set = Poisson_VBF_8$particles,
                                                       multivariate = TRUE,
                                                       resampling_method = 'resid',
                                                       seed = seed)
print(integrated_abs_distance(full_posterior, Poisson_VBF_8$particles$y_samples))
compare_samples_bivariate(posteriors = list(full_posterior,
                                            Poisson_VBF_8$particles$y_samples,
                                            Poisson_VBF_8$proposed_samples),
                          colours = c('black', 'red', 'green'),
                          common_limit = c(-4, 4))

##### NB (Hypercube Centre) #####
print('NB Fusion (hypercube centre)')
NB_VBF_8 <- parallel_generalised_BF_SMC_BLR(particles_to_fuse = particles_to_fuse,
                                            N = nsamples,
                                            m = 8,
                                            time_mesh = time_mesh,
                                            dim = 5,
                                            data_split = data_split_8,
                                            prior_means = rep(0, 5),
                                            prior_variances = rep(1, 5),
                                            C = 8,
                                            precondition_matrices = rep(list(diag(1, 5)), 8),
                                            resampling_method = 'resid',
                                            ESS_threshold = 0.5,
                                            cv_location = 'hypercube_centre',
                                            diffusion_estimator = 'NB',
                                            seed = seed,
                                            n_cores = n_cores,
                                            print_progress_iters = 500)
NB_VBF_8$particles <- resample_particle_y_samples(particle_set = NB_VBF_8$particles,
                                                  multivariate = TRUE,
                                                  resampling_method = 'resid',
                                                  seed = seed)
print(integrated_abs_distance(full_posterior, NB_VBF_8$particles$y_samples))
compare_samples_bivariate(posteriors = list(full_posterior,
                                            NB_VBF_8$particles$y_samples,
                                            NB_VBF_8$proposed_samples),
                          colours = c('black', 'red', 'green'),
                          common_limit = c(-4, 4))

##### Applying Generalised Bayesian Fusion #####

time_choice <- 2
time_mesh <- seq(0, time_choice, 0.005)
particles_to_fuse <- initialise_particle_sets(samples_to_fuse = sub_posteriors_8,
                                              multivariate = TRUE,
                                              number_of_steps = length(time_mesh))

##### Poisson (Hypercube Centre) #####
print('Poisson Fusion (hypercube centre)')
Poisson_GBF_8 <- parallel_generalised_BF_SMC_BLR(particles_to_fuse = particles_to_fuse,
                                                 N = nsamples,
                                                 m = 8,
                                                 time_mesh = time_mesh,
                                                 dim = 5,
                                                 data_split = data_split_8,
                                                 prior_means = rep(0, 5),
                                                 prior_variances = rep(1, 5),
                                                 C = 8,
                                                 precondition_matrices = lapply(sub_posteriors_8, cov),
                                                 resampling_method = 'resid',
                                                 ESS_threshold = 0.5,
                                                 cv_location = 'hypercube_centre',
                                                 diffusion_estimator = 'Poisson',
                                                 seed = seed,
                                                 n_cores = n_cores,
                                                 print_progress_iters = 500)
Poisson_GBF_8$particles <- resample_particle_y_samples(particle_set = Poisson_GBF_8$particles,
                                                       multivariate = TRUE,
                                                       resampling_method = 'resid',
                                                       seed = seed)
print(integrated_abs_distance(full_posterior, Poisson_GBF_8$particles$y_samples))
compare_samples_bivariate(posteriors = list(full_posterior,
                                            Poisson_GBF_8$proposed_samples,
                                            Poisson_GBF_8$particles$y_samples),
                          colours = c('black', 'green', 'red'),
                          common_limit = c(-4, 4))

##### NB (Hypercube Centre) #####
print('NB Fusion (hypercube centre)')
NB_GBF_8 <- parallel_generalised_BF_SMC_BLR(particles_to_fuse = particles_to_fuse,
                                            N = nsamples,
                                            m = 8,
                                            time_mesh = time_mesh,
                                            dim = 5,
                                            data_split = data_split_8,
                                            prior_means = rep(0, 5),
                                            prior_variances = rep(1, 5),
                                            C = 8,
                                            precondition_matrices = lapply(sub_posteriors_8, cov),
                                            resampling_method = 'resid',
                                            ESS_threshold = 0.5,
                                            cv_location = 'hypercube_centre',
                                            diffusion_estimator = 'NB',
                                            seed = seed,
                                            n_cores = n_cores,
                                            print_progress_iters = 500)
NB_GBF_8$particles <- resample_particle_y_samples(particle_set = NB_GBF_8$particles,
                                                  multivariate = TRUE,
                                                  resampling_method = 'resid',
                                                  seed = seed)
print(integrated_abs_distance(full_posterior, NB_GBF_8$particles$y_samples))
compare_samples_bivariate(posteriors = list(full_posterior,
                                            NB_GBF_8$proposed_samples,
                                            NB_GBF_8$particles$y_samples),
                          colours = c('black', 'green', 'red'),
                          common_limit = c(-4, 4))

##### IAD #####

integrated_abs_distance(full_posterior, Poisson_VBF_8$particles$y_samples)
integrated_abs_distance(full_posterior, NB_VBF_8$particles$y_samples)
integrated_abs_distance(full_posterior, Poisson_GBF_8$particles$y_samples)
integrated_abs_distance(full_posterior, NB_GBF_8$particles$y_samples)

##### Save data #####

save.image('SD8_GBF.RData')
