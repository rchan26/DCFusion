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

##### Applying Fusion #####

##### Poisson (Hypercube Centre) #####
print('Poisson Fusion (hypercube centre)')
Poisson_hc_16 <- bal_binary_fusion_SMC_BLR(N_schedule = rep(nsamples, 4),
                                           m_schedule = rep(2, 4),
                                           time_schedule = rep(time_choice, 4),
                                           base_samples = sub_posteriors_16,
                                           L = 5,
                                           dim = 5,
                                           data_split = data_split_16,
                                           prior_means = rep(0, 5),
                                           prior_variances = rep(1, 5),
                                           C = 16,
                                           precondition = TRUE,
                                           resampling_method = 'resid',
                                           ESS_threshold = 0.5,
                                           cv_location = 'hypercube_centre',
                                           diffusion_estimator = 'Poisson',
                                           seed = seed,
                                           n_cores = n_cores,
                                           print_progress_iters = 500)
Poisson_hc_16$particles <- resample_particle_y_samples(particle_set = Poisson_hc_16$particles[[1]],
                                                       multivariate = TRUE,
                                                       resampling_method = 'resid',
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
                                      dim = 5,
                                      data_split = data_split_16,
                                      prior_means = rep(0, 5),
                                      prior_variances = rep(1, 5),
                                      C = 16,
                                      precondition = TRUE,
                                      resampling_method = 'resid',
                                      ESS_threshold = 0.5,
                                      cv_location = 'hypercube_centre',
                                      diffusion_estimator = 'NB',
                                      seed = seed,
                                      n_cores = n_cores,
                                      print_progress_iters = 500)
NB_hc_16$particles <- resample_particle_y_samples(particle_set = NB_hc_16$particles[[1]],
                                                  multivariate = TRUE,
                                                  resampling_method = 'resid',
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

save.image('SD16.RData')
