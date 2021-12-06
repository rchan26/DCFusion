library(DCFusion)
library(HMCBLR)

# function to simulate data given frequencies of the coefficients and the coefficients themselves
simulate_data <- function(N,
                          alpha,
                          frequencies,
                          coefficients,
                          G_means = rep(1, length(frequencies)),
                          G_sds = rep(1, length(frequencies)),
                          seed = NULL) {
  if (length(frequencies) != length(coefficients)) {
    stop("length of frequencies and length of coefficients are not the same")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  dim <- length(frequencies)
  X <- matrix(nrow = N, ncol = dim)
  for (d in 1:dim) {
    for (i in 1:N) {
      if (runif(1, 0, 1) < frequencies[d]) {
        X[i,d] <- rnorm(1, G_means[d], G_sds[d])
      } else {
        X[i,d] <- 0
      }
    }
  }
  y <- rep(0, N)
  for (i in 1:N) {
    beta_X <- alpha + sum(X[i,] * coefficients)
    prob <- 1 / (1+exp(-beta_X))
    y[i] <- rbinom(1, 1, prob = prob)
  }
  return(cbind(y, data.frame(X)))
}

# function to check the frequencies of 'active' variables in the dataset
check_activity <- function(data, proportion = T) {
  if (!is.data.frame(data)) {
    stop("check_activity: data is not in a data frame format")
  }
  # create new data frame with same column names
  active_df <- data.frame(matrix(nrow = 1, ncol = ncol(data)))
  colnames(active_df) <- colnames(data)
  # loop through columns and check how many are 'active'
  # i.e. have a 1 in the instance
  for (j in 1:ncol(data)) {
    if (proportion) {
      active_df[,j] <- sum(data[,j]!=0) / nrow(data)
    } else {
      active_df[,j] <- sum(data[,j]!=0)
    }
  }
  print('proportion that each variable is active in the dataset:')
  for (j in 1:ncol(active_df)) {
    print(paste(colnames(active_df)[j], ':', active_df[1,j]))
  }
}

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
simulated_data <- simulate_data(N = ndata,
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

##### Applying Vanilla Bayesian Fusion (by passing in identity matrices) #####

time_choice <- 2
time_mesh <- seq(0, time_choice, 0.005)
particles_to_fuse <- initialise_particle_sets(samples_to_fuse = sub_posteriors_4,
                                              multivariate = TRUE,
                                              number_of_steps = length(time_mesh))

##### Poisson (Hypercube Centre) #####
print('Poisson Fusion (hypercube centre)')
Poisson_GBF_identity_4 <- parallel_generalised_BF_SMC_BLR(particles_to_fuse = particles_to_fuse,
                                                          N = nsamples,
                                                          m = 4,
                                                          time_mesh = time_mesh,
                                                          dim = 5,
                                                          data_split = data_split_4,
                                                          prior_means = rep(0, 5),
                                                          prior_variances = rep(1, 5),
                                                          C = 4,
                                                          precondition_matrices = rep(list(diag(1, 5)), 4),
                                                          resampling_method = 'resid',
                                                          ESS_threshold = 0.5,
                                                          cv_location = 'hypercube_centre',
                                                          diffusion_estimator = 'Poisson',
                                                          seed = seed,
                                                          n_cores = n_cores,
                                                          print_progress_iters = 500)
Poisson_GBF_identity_4$particles <- resample_particle_y_samples(particle_set = Poisson_GBF_identity_4$particles,
                                                                multivariate = TRUE,
                                                                resampling_method = 'resid',
                                                                seed = seed)
print(integrated_abs_distance(full_posterior, Poisson_GBF_identity_4$particles$y_samples))
compare_samples_bivariate(posteriors = list(full_posterior,
                                            Poisson_GBF_identity_4$particles$y_samples,
                                            Poisson_GBF_identity_4$proposed_samples),
                          colours = c('black', 'red', 'green'),
                          common_limit = c(-4, 4))

##### NB (Hypercube Centre) #####
print('NB Fusion (hypercube centre)')
NB_GBF_identity_4 <- parallel_generalised_BF_SMC_BLR(particles_to_fuse = particles_to_fuse,
                                                     N = nsamples,
                                                     m = 4,
                                                     time_mesh = time_mesh,
                                                     dim = 5,
                                                     data_split = data_split_4,
                                                     prior_means = rep(0, 5),
                                                     prior_variances = rep(1, 5),
                                                     C = 4,
                                                     precondition_matrices = rep(list(diag(1, 5)), 4),
                                                     resampling_method = 'resid',
                                                     ESS_threshold = 0.5,
                                                     cv_location = 'hypercube_centre',
                                                     diffusion_estimator = 'NB',
                                                     seed = seed,
                                                     n_cores = n_cores,
                                                     print_progress_iters = 500)
NB_GBF_identity_4$particles <- resample_particle_y_samples(particle_set = NB_GBF_identity_4$particles,
                                                           multivariate = TRUE,
                                                           resampling_method = 'resid',
                                                           seed = seed)
print(integrated_abs_distance(full_posterior, NB_GBF_identity_4$particles$y_samples))
compare_samples_bivariate(posteriors = list(full_posterior,
                                            NB_GBF_identity_4$particles$y_samples,
                                            NB_GBF_identity_4$proposed_samples),
                          colours = c('black', 'red', 'green'),
                          common_limit = c(-4, 4))

##### Applying Generalised Bayesian Fusion #####

time_choice <- 1
time_mesh <- seq(0, time_choice, 0.01)
particles_to_fuse <- initialise_particle_sets(samples_to_fuse = sub_posteriors_4,
                                              multivariate = TRUE,
                                              number_of_steps = length(time_mesh))

##### Poisson (Hypercube Centre) #####
print('Poisson Fusion (hypercube centre)')
Poisson_GBF_4 <- parallel_generalised_BF_SMC_BLR(particles_to_fuse = particles_to_fuse,
                                                 N = nsamples,
                                                 m = 4,
                                                 time_mesh = time_mesh,
                                                 dim = 5,
                                                 data_split = data_split_4,
                                                 prior_means = rep(0, 5),
                                                 prior_variances = rep(1, 5),
                                                 C = 4,
                                                 precondition_matrices = lapply(sub_posteriors_4, cov),
                                                 resampling_method = 'resid',
                                                 ESS_threshold = 0.5,
                                                 cv_location = 'hypercube_centre',
                                                 diffusion_estimator = 'Poisson',
                                                 seed = seed,
                                                 n_cores = n_cores,
                                                 print_progress_iters = 500)
Poisson_GBF_4$particles <- resample_particle_y_samples(particle_set = Poisson_GBF_4$particles,
                                                       multivariate = TRUE,
                                                       resampling_method = 'resid',
                                                       seed = seed)
print(integrated_abs_distance(full_posterior, Poisson_GBF_4$particles$y_samples))
compare_samples_bivariate(posteriors = list(full_posterior,
                                            Poisson_GBF_4$particles$y_samples,
                                            Poisson_GBF_4$proposed_samples),
                          colours = c('black', 'red', 'green'),
                          common_limit = c(-4, 4))

##### NB (Hypercube Centre) #####
print('NB Fusion (hypercube centre)')
NB_GBF_4 <- parallel_generalised_BF_SMC_BLR(particles_to_fuse = particles_to_fuse,
                                            N = nsamples,
                                            m = 4,
                                            time_mesh = time_mesh,
                                            dim = 5,
                                            data_split = data_split_4,
                                            prior_means = rep(0, 5),
                                            prior_variances = rep(1, 5),
                                            C = 4,
                                            precondition_matrices = lapply(sub_posteriors_4, cov),
                                            resampling_method = 'resid',
                                            ESS_threshold = 0.5,
                                            cv_location = 'hypercube_centre',
                                            diffusion_estimator = 'NB',
                                            seed = seed,
                                            n_cores = n_cores,
                                            print_progress_iters = 500)
NB_GBF_4$particles <- resample_particle_y_samples(particle_set = NB_GBF_4$particles,
                                                  multivariate = TRUE,
                                                  resampling_method = 'resid',
                                                  seed = seed)
print(integrated_abs_distance(full_posterior, NB_GBF_4$particles$y_samples))
compare_samples_bivariate(posteriors = list(full_posterior,
                                            NB_GBF_4$particles$y_samples,
                                            NB_GBF_4$proposed_samples),
                          colours = c('black', 'red', 'green'),
                          common_limit = c(-4, 4))
