library(DCFusion)
library(HMCBLR)

##### Initialise example #####
seed <- 2021
set.seed(seed)
nsamples <- 10000
ndata <- 1000
C <- 16
time_choice <- 0.5
n_cores <- parallel::detectCores()
true_beta <- c(-3, 1.2, -0.5, 0.8, 3)
frequencies <- c(0.2, 0.3, 0.5, 0.01) # must have length = length(true_beta)-1
# variables for choosing guidance
ESS_threshold <- 0.5
CESS_0_threshold <- 0.2
k1 <- NULL
k2 <- NULL
k3 <- 2
k4 <- 1

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

##### Sampling from sub-posterior C=16 #####

data_split_16 <- split_data(simulated_data, y_col_index = 1, X_col_index = 2:ncol(simulated_data), C = C, as_dataframe = F)
sub_posteriors_16 <- hmc_base_sampler_BLR(nsamples = nsamples,
                                          data_split = data_split_16,
                                          C = C, 
                                          prior_means = rep(0, 5),
                                          prior_variances = rep(1, 5),
                                          warmup = 10000,
                                          seed = seed,
                                          output = T)

##### Applying Vanilla Bayesian Fusion (by passing in identity matrices) #####

# compute estimate for the b for vanilla guidance
data_size <- ceiling(ndata/16)
precondition_matrices <- lapply(sub_posteriors_16, cov)
inv_precondition_matrices <- lapply(precondition_matrices, solve)
Lambda <- inverse_sum_matrices(inv_precondition_matrices)
vanilla_b <- mean(diag(Lambda)*data_size)
# computing guidance for parameters
sub_posterior_means <- t(sapply(sub_posteriors_16, function(sub) apply(sub, 2, mean)))
vanilla_guide <- BF_guidance(condition = 'SH',
                             CESS_0_threshold = CESS_0_threshold,
                             C = C,
                             d = 5,
                             data_size = data_size,
                             b = vanilla_b,
                             sub_posterior_means = sub_posterior_means,
                             k1 = k1,
                             k3 = k3,
                             k4 = k4,
                             vanilla = TRUE)
particles_to_fuse <- initialise_particle_sets(samples_to_fuse = sub_posteriors_16,
                                              multivariate = TRUE,
                                              number_of_steps = length(vanilla_guide$mesh))

##### Poisson (Hypercube Centre) #####
print('Poisson Fusion (hypercube centre)')
Poisson_VBF_16 <- list('reg_mesh' = parallel_GBF_BLR(particles_to_fuse = particles_to_fuse,
                                                     N = nsamples,
                                                     m = 16,
                                                     time_mesh = vanilla_guide$mesh,
                                                     dim = 5,
                                                     data_split = data_split_16,
                                                     prior_means = rep(0, 5),
                                                     prior_variances = rep(1, 5),
                                                     C = C,
                                                     precondition_matrices = rep(list(diag(1, 5)), 16),
                                                     resampling_method = 'resid',
                                                     ESS_threshold = 0.5,
                                                     sub_posterior_means = sub_posterior_means,
                                                     adaptive_mesh = FALSE,
                                                     cv_location = 'hypercube_centre',
                                                     diffusion_estimator = 'Poisson',
                                                     seed = seed,
                                                     n_cores = n_cores,
                                                     print_progress_iters = 500),
                       'adaptive_mesh' = parallel_GBF_BLR(particles_to_fuse = particles_to_fuse,
                                                          N = nsamples,
                                                          m = 16,
                                                          time_mesh = vanilla_guide$mesh,
                                                          dim = 5,
                                                          data_split = data_split_16,
                                                          prior_means = rep(0, 5),
                                                          prior_variances = rep(1, 5),
                                                          C = C,
                                                          precondition_matrices = rep(list(diag(1, 5)), 16),
                                                          resampling_method = 'resid',
                                                          ESS_threshold = 0.5,
                                                          sub_posterior_means = sub_posterior_means,
                                                          adaptive_mesh = TRUE,
                                                          adaptive_mesh_parameters = list('data_size' = data_size,
                                                                                          'b' = vanilla_b,
                                                                                          'k3' = k3,
                                                                                          'k4' = k4,
                                                                                          'vanilla' = TRUE),
                                                          cv_location = 'hypercube_centre',
                                                          diffusion_estimator = 'Poisson',
                                                          seed = seed,
                                                          n_cores = n_cores,
                                                          print_progress_iters = 500))
# regular mesh
Poisson_VBF_16$reg_mesh$particles <- resample_particle_y_samples(particle_set = Poisson_VBF_16$reg_mesh$particles,
                                                                 multivariate = TRUE,
                                                                 resampling_method = 'resid',
                                                                 seed = seed)
print(integrated_abs_distance(full_posterior, Poisson_VBF_16$reg_mesh$particles$y_samples))
compare_samples_bivariate(posteriors = list(full_posterior,
                                            Poisson_VBF_16$reg_mesh$particles$y_samples,
                                            Poisson_VBF_16$reg_mesh$proposed_samples),
                          colours = c('black', 'red', 'green'),
                          common_limit = c(-4, 4))
# adaptive mesh
Poisson_VBF_16$adaptive_mesh$particles <- resample_particle_y_samples(particle_set = Poisson_VBF_16$adaptive_mesh$particles,
                                                                      multivariate = TRUE,
                                                                      resampling_method = 'resid',
                                                                      seed = seed)
print(integrated_abs_distance(full_posterior, Poisson_VBF_16$adaptive_mesh$particles$y_samples))
compare_samples_bivariate(posteriors = list(full_posterior,
                                            Poisson_VBF_16$adaptive_mesh$particles$y_samples,
                                            Poisson_VBF_16$adaptive_mesh$proposed_samples),
                          colours = c('black', 'red', 'green'),
                          common_limit = c(-4, 4))

##### NB (Hypercube Centre) #####
print('NB Fusion (hypercube centre)')
NB_VBF_16 <- list('reg_mesh' = parallel_GBF_BLR(particles_to_fuse = particles_to_fuse,
                                                N = nsamples,
                                                m = 16,
                                                time_mesh = vanilla_guide$mesh,
                                                dim = 5,
                                                data_split = data_split_16,
                                                prior_means = rep(0, 5),
                                                prior_variances = rep(1, 5),
                                                C = C,
                                                precondition_matrices = rep(list(diag(1, 5)), 16),
                                                resampling_method = 'resid',
                                                ESS_threshold = 0.5,
                                                sub_posterior_means = sub_posterior_means,
                                                adaptive_mesh = FALSE,
                                                cv_location = 'hypercube_centre',
                                                diffusion_estimator = 'NB',
                                                seed = seed,
                                                n_cores = n_cores,
                                                print_progress_iters = 500),
                  'adaptive_mesh' = parallel_GBF_BLR(particles_to_fuse = particles_to_fuse,
                                                     N = nsamples,
                                                     m = 16,
                                                     time_mesh = vanilla_guide$mesh,
                                                     dim = 5,
                                                     data_split = data_split_16,
                                                     prior_means = rep(0, 5),
                                                     prior_variances = rep(1, 5),
                                                     C = C,
                                                     precondition_matrices = rep(list(diag(1, 5)), 16),
                                                     resampling_method = 'resid',
                                                     ESS_threshold = 0.5,
                                                     sub_posterior_means = sub_posterior_means,
                                                     adaptive_mesh = TRUE,
                                                     adaptive_mesh_parameters = list('data_size' = data_size,
                                                                                     'b' = vanilla_b,
                                                                                     'k3' = k3,
                                                                                     'k4' = k4,
                                                                                     'vanilla' = TRUE),
                                                     cv_location = 'hypercube_centre',
                                                     diffusion_estimator = 'NB',
                                                     seed = seed,
                                                     n_cores = n_cores,
                                                     print_progress_iters = 500))

# regular mesh
NB_VBF_16$reg_mesh$particles <- resample_particle_y_samples(particle_set = NB_VBF_16$reg_mesh$particles,
                                                            multivariate = TRUE,
                                                            resampling_method = 'resid',
                                                            seed = seed)
print(integrated_abs_distance(full_posterior, NB_VBF_16$reg_mesh$particles$y_samples))
compare_samples_bivariate(posteriors = list(full_posterior,
                                            NB_VBF_16$reg_mesh$particles$y_samples,
                                            NB_VBF_16$reg_mesh$proposed_samples),
                          colours = c('black', 'red', 'green'),
                          common_limit = c(-4, 4))
# adaptive mesh
NB_VBF_16$reg_mesh$particles <- resample_particle_y_samples(particle_set = NB_VBF_16$adaptive_mesh$particles,
                                                            multivariate = TRUE,
                                                            resampling_method = 'resid',
                                                            seed = seed)
print(integrated_abs_distance(full_posterior, NB_VBF_16$adaptive_mesh$particles$y_samples))
compare_samples_bivariate(posteriors = list(full_posterior,
                                            NB_VBF_16$adaptive_mesh$particles$y_samples,
                                            NB_VBF_16$adaptive_mesh$proposed_samples),
                          colours = c('black', 'red', 'green'),
                          common_limit = c(-4, 4))

##### Applying Generalised Bayesian Fusion #####

# computing guidance for parameters
gen_guide <- BF_guidance(condition = 'SH',
                         CESS_0_threshold = CESS_0_threshold,
                         C = C,
                         d = 5,
                         data_size = data_size,
                         sub_posterior_means = sub_posterior_means,
                         precondition_matrices = precondition_matrices,
                         inv_precondition_matrices = inv_precondition_matrices,
                         Lambda = Lambda,
                         k1 = k1,
                         k3 = k3,
                         k4 = k4,
                         vanilla = FALSE)
particles_to_fuse <- initialise_particle_sets(samples_to_fuse = sub_posteriors_16,
                                              multivariate = TRUE,
                                              number_of_steps = length(gen_guide$mesh))

##### Poisson (Hypercube Centre) #####
print('Poisson Fusion (hypercube centre)')
Poisson_GBF_16 <- list('reg_mesh' = parallel_GBF_BLR(particles_to_fuse = particles_to_fuse,
                                                     N = nsamples,
                                                     m = 16,
                                                     time_mesh = gen_guide$mesh,
                                                     dim = 5,
                                                     data_split = data_split_16,
                                                     prior_means = rep(0, 5),
                                                     prior_variances = rep(1, 5),
                                                     C = C,
                                                     precondition_matrices = precondition_matrices,
                                                     resampling_method = 'resid',
                                                     ESS_threshold = 0.5,
                                                     sub_posterior_means = sub_posterior_means,
                                                     adaptive_mesh = FALSE,
                                                     cv_location = 'hypercube_centre',
                                                     diffusion_estimator = 'Poisson',
                                                     seed = seed,
                                                     n_cores = n_cores,
                                                     print_progress_iters = 500),
                       'adaptive_mesh' = parallel_GBF_BLR(particles_to_fuse = particles_to_fuse,
                                                          N = nsamples,
                                                          m = 16,
                                                          time_mesh = gen_guide$mesh,
                                                          dim = 5,
                                                          data_split = data_split_16,
                                                          prior_means = rep(0, 5),
                                                          prior_variances = rep(1, 5),
                                                          C = C,
                                                          precondition_matrices = precondition_matrices,
                                                          resampling_method = 'resid',
                                                          ESS_threshold = 0.5,
                                                          sub_posterior_means = sub_posterior_means,
                                                          adaptive_mesh = TRUE,
                                                          adaptive_mesh_parameters = list('data_size' = data_size,
                                                                                          'k3' = k3,
                                                                                          'k4' = k4,
                                                                                          'vanilla' = FALSE),
                                                          cv_location = 'hypercube_centre',
                                                          diffusion_estimator = 'Poisson',
                                                          seed = seed,
                                                          n_cores = n_cores,
                                                          print_progress_iters = 500))
# regular mesh
Poisson_GBF_16$reg_mesh$particles <- resample_particle_y_samples(particle_set = Poisson_GBF_16$reg_mesh$particles,
                                                                 multivariate = TRUE,
                                                                 resampling_method = 'resid',
                                                                 seed = seed)
print(integrated_abs_distance(full_posterior, Poisson_GBF_16$reg_mesh$particles$y_samples))
compare_samples_bivariate(posteriors = list(full_posterior,
                                            Poisson_GBF_16$reg_mesh$proposed_samples,
                                            Poisson_GBF_16$reg_mesh$particles$y_samples),
                          colours = c('black', 'green', 'red'),
                          common_limit = c(-4, 4))
# adaptive mesh
Poisson_GBF_16$adaptive_mesh$particles <- resample_particle_y_samples(particle_set = Poisson_GBF_16$adaptive_mesh$particles,
                                                                      multivariate = TRUE,
                                                                      resampling_method = 'resid',
                                                                      seed = seed)
print(integrated_abs_distance(full_posterior, Poisson_GBF_16$adaptive_mesh$particles$y_samples))
compare_samples_bivariate(posteriors = list(full_posterior,
                                            Poisson_GBF_16$adaptive_mesh$proposed_samples,
                                            Poisson_GBF_16$adaptive_mesh$particles$y_samples),
                          colours = c('black', 'green', 'red'),
                          common_limit = c(-4, 4))

##### NB (Hypercube Centre) #####
print('NB Fusion (hypercube centre)')
NB_GBF_16 <- list('reg_mesh' = parallel_GBF_BLR(particles_to_fuse = particles_to_fuse,
                                                N = nsamples,
                                                m = 16,
                                                time_mesh = gen_guide$mesh,
                                                dim = 5,
                                                data_split = data_split_16,
                                                prior_means = rep(0, 5),
                                                prior_variances = rep(1, 5),
                                                C = C,
                                                precondition_matrices = precondition_matrices,
                                                resampling_method = 'resid',
                                                ESS_threshold = 0.5,
                                                sub_posterior_means = sub_posterior_means,
                                                adaptive_mesh = FALSE,
                                                diffusion_estimator = 'NB',
                                                seed = seed,
                                                n_cores = n_cores,
                                                print_progress_iters = 500),
                  'adaptive_mesh' = parallel_GBF_BLR(particles_to_fuse = particles_to_fuse,
                                                     N = nsamples,
                                                     m = 16,
                                                     time_mesh = gen_guide$mesh,
                                                     dim = 5,
                                                     data_split = data_split_16,
                                                     prior_means = rep(0, 5),
                                                     prior_variances = rep(1, 5),
                                                     C = C,
                                                     precondition_matrices = precondition_matrices,
                                                     resampling_method = 'resid',
                                                     ESS_threshold = 0.5,
                                                     sub_posterior_means = sub_posterior_means,
                                                     adaptive_mesh = TRUE,
                                                     adaptive_mesh_parameters = list('data_size' = data_size,
                                                                                     'k3' = k3,
                                                                                     'k4' = k4,
                                                                                     'vanilla' = FALSE),
                                                     cv_location = 'hypercube_centre',
                                                     diffusion_estimator = 'NB',
                                                     seed = seed,
                                                     n_cores = n_cores,
                                                     print_progress_iters = 500))
# regular mesh
NB_GBF_16$reg_mesh$particles <- resample_particle_y_samples(particle_set = NB_GBF_16$reg_mesh$particles,
                                                            multivariate = TRUE,
                                                            resampling_method = 'resid',
                                                            seed = seed)
print(integrated_abs_distance(full_posterior, NB_GBF_16$reg_mesh$particles$y_samples))
compare_samples_bivariate(posteriors = list(full_posterior,
                                            NB_GBF_16$reg_mesh$proposed_samples,
                                            NB_GBF_16$reg_mesh$particles$y_samples),
                          colours = c('black', 'green', 'red'),
                          common_limit = c(-4, 4))
# adaptive mesh
NB_GBF_16$adaptive_mesh$particles <- resample_particle_y_samples(particle_set = NB_GBF_16$adaptive_mesh$particles,
                                                                 multivariate = TRUE,
                                                                 resampling_method = 'resid',
                                                                 seed = seed)
print(integrated_abs_distance(full_posterior, NB_GBF_16$adaptive_mesh$particles$y_samples))
compare_samples_bivariate(posteriors = list(full_posterior,
                                            NB_GBF_16$adaptive_mesh$proposed_samples,
                                            NB_GBF_16$adaptive_mesh$particles$y_samples),
                          colours = c('black', 'green', 'red'),
                          common_limit = c(-4, 4))

##### IAD #####

integrated_abs_distance(full_posterior, Poisson_VBF_16$reg_mesh$particles$y_samples)
integrated_abs_distance(full_posterior, Poisson_VBF_16$adaptive_mesh$particles$y_samples)
integrated_abs_distance(full_posterior, NB_VBF_16$reg_mesh$particles$y_samples)
integrated_abs_distance(full_posterior, NB_VBF_16$adaptive_mesh$particles$y_samples)
integrated_abs_distance(full_posterior, Poisson_GBF_16$reg_mesh$particles$y_samples)
integrated_abs_distance(full_posterior, Poisson_GBF_16$adaptive_mesh$particles$y_samples)
integrated_abs_distance(full_posterior, NB_GBF_16$reg_mesh$particles$y_samples)
integrated_abs_distance(full_posterior, NB_GBF_16$adaptive_mesh$particles$y_samples)

##### Save data #####

save.image('SD16_GBF.RData')
