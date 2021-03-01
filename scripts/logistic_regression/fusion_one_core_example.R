library(hierarchicalFusion)
library(HMCBLR)

# function to simulate data given frequencies of the coefficients and the coefficients themselves
simulate_data <- function(N, alpha, frequencies, coefficients) {
  if (length(frequencies) != length(coefficients)) {
    stop("length of frequencies and length of coefficients are not the same")
  }
  dim <- length(frequencies)
  X <- matrix(nrow = N, ncol = dim)
  for (d in 1:dim) {
    X[,d] <- sample(c(0,1), size = N, replace = T, prob = c((1-frequencies[d]), frequencies[d]))
  }
  y <- rep(0, N)
  for (i in 1:N) {
    beta_X <- alpha + sum(X[i,] * coefficients)
    prob <- 1 / (1+exp(-beta_X))
    y[i] <- rbinom(1, 1, prob = prob)
  }
  data <- data.frame(X)
  return(cbind(y, data))
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
      active_df[,j] <- sum(data[,j]==1) / nrow(data)
    } else {
      active_df[,j] <- sum(data[,j]==1)
    }
  }
  print('proportion that each variable is active in the dataset:')
  for (j in 1:ncol(active_df)) {
    print(paste(colnames(active_df)[j], ':', active_df[1,j]))
  }
}

# set a seed
seed <- 1920
set.seed(seed)

# simulate data set
simulated_data <- simulate_data(N = 1000,
                                alpha = -3,
                                frequencies = c(0.2, 0.5),
                                coefficients = c(1.2, 0.8))

# check activity of the parameters
check_activity(simulated_data)

# organise simulated data
data <- list()
data$y <- simulated_data$y
data$X <- as.matrix(cbind(rep(1, nrow(simulated_data[,2:3])), simulated_data[,2:3]))
colnames(data$X)[1] <- 'intercept'

###### Sampling from full posterier #####

full_posterior <- hmc_sample_BLR(data = data,
                                 C = 1,
                                 prior_means = rep(0, 3),
                                 prior_variances = rep(1, 3),
                                 iterations = 110000,
                                 warmup = 10000,
                                 chains = 1,
                                 seed = seed,
                                 output = T)

##### Fusion on one core #####

particles_to_fuse <- initialise_particle_sets(samples_to_fuse = list(full_posterior),
                                              multivariate = TRUE)

data_split <- split_data(simulated_data, y_col_index = 1, X_col_index = 2:3, C = 1, as_dataframe = F)

# preconditioned fork and join [Poisson estimator]
test_preconditioned_SMC_Poisson <- parallel_fusion_SMC_BLR(particles_to_fuse = particles_to_fuse,
                                                           N = 100000,
                                                           m = 1,
                                                           time = 0.5,
                                                           dim = 3,
                                                           data_split = data_split,
                                                           prior_means = rep(0, 3),
                                                           prior_variances = rep(1, 3),
                                                           C = 1,
                                                           precondition_matrices = list(cov(full_posterior)),
                                                           resampling_method = 'resid',
                                                           ESS_threshold = 0.5,
                                                           diffusion_estimator = 'Poisson',
                                                           bounds_multiplier = 2.5,
                                                           seed = seed)
test_preconditioned_SMC_Poisson$particles <- resample_particle_y_samples(particle_set = test_preconditioned_SMC_Poisson$particles,
                                                                         multivariate = TRUE,
                                                                         resampling_method = 'resid',
                                                                         seed = seed)

# preconditioned fork and join [Negative Binomial estimator]
test_preconditioned_SMC_NB <- parallel_fusion_SMC_BLR(particles_to_fuse = particles_to_fuse,
                                                      N = 100000,
                                                      m = 1,
                                                      time = 0.5,
                                                      dim = 3,
                                                      data_split = data_split,
                                                      prior_means = rep(0, 3),
                                                      prior_variances = rep(1, 3),
                                                      C = 1,
                                                      precondition_matrices = list(cov(full_posterior)),
                                                      resampling_method = 'resid',
                                                      ESS_threshold = 0.5,
                                                      diffusion_estimator = 'NB',
                                                      bounds_multiplier = 2.5,
                                                      seed = seed)
test_preconditioned_SMC_NB$particles <- resample_particle_y_samples(particle_set = test_preconditioned_SMC_NB$particles,
                                                                    multivariate = TRUE,
                                                                    resampling_method = 'resid',
                                                                    seed = seed)

compare_samples_bivariate(posteriors = list(full_posterior,
                                            test_preconditioned_SMC_Poisson$particles$y_samples,
                                            test_preconditioned_SMC_NB$particles$y_samples),
                          colours = c('black', 'red', 'blue'),
                          common_limit = c(-4, 4),
                          title = 'Fusion on one core (Precondition) || N = 100000')
