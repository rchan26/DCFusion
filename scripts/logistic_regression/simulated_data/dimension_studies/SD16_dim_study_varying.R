library(DCFusion)
library(HMCBLR)

##### Initialise example #####
seed <- 2021
set.seed(seed)
nsamples <- c(10000, 30000, 50000, 100000)
ndata <- 5000
n_cores <- parallel::detectCores()
ESS_threshold <- 0.5
resampling_method <- 'resid'
dim <- c(5,6,7,8,9,10,12,14)
time_choices <- c(0.1, 0.2, 0.25, 0.4, 0.5)
true_beta <- list()
frequencies <- list()
simulated_data <- list()
full_posterior <- list()
data_split_16 <- list()
sub_posteriors_16 <- list()
results <- list()

for (d in 1:length(dim)) {
  print(paste('$$$$$$$$$$ d:', d))
  print(paste('dim:', dim[d]))
  set.seed(seed)
  prior_means <- rep(0, dim[d])
  prior_variances <- rep(1, dim[d])
  true_beta[[d]] <- round(rnorm(dim[d], 0, 1), 2)
  frequencies[[d]] <- round(runif(dim[d]-1, 0.01, 0.99), 2)
  results[[d]] <- list()
  
  # simulate data set
  simulated_data[[d]] <- simulate_LR_data(N = ndata,
                                          alpha = true_beta[[d]][1],
                                          frequencies = frequencies[[d]],
                                          coefficients = true_beta[[d]][2:length(true_beta[[d]])],
                                          seed = seed)
  
  # check activity of the parameters
  check_activity(simulated_data[[d]])
  
  ##### Sampling from full posterior #####
  
  full_data_count <- unique_row_count(y = simulated_data[[d]][,1],
                                      X = cbind('intercept' = rep(1, ndata), simulated_data[[d]][,2:ncol(simulated_data[[d]])]))$full_data_count
  full_posterior[[d]] <- hmc_sample_BLR(full_data_count = full_data_count,
                                        C = 1,
                                        prior_means = prior_means,
                                        prior_variances = prior_variances,
                                        iterations = 110000,
                                        warmup = 10000,
                                        chains = 1,
                                        seed = seed,
                                        output = T)
  
  ##### Sampling from sub-posterior C=16 #####
  
  data_split_16[[d]] <- split_data(simulated_data[[d]],
                                   y_col_index = 1,
                                   X_col_index = 2:ncol(simulated_data[[d]]),
                                   C = 16,
                                   as_dataframe = F)
  sub_posteriors_16[[d]] <- hmc_base_sampler_BLR(nsamples = 100000,
                                                 data_split = data_split_16[[d]],
                                                 C = 16,
                                                 prior_means = prior_means,
                                                 prior_variances = prior_variances,
                                                 warmup = 10000,
                                                 seed = seed,
                                                 output = T)
  
  ##### Applying Fusion #####
  
  ##### NB (Hypercube Centre) #####
  print('NB Fusion (hypercube centre)')
  for (i in 1:length(time_choices)) {
    print(paste('##### i:', i, '|| time:', time_choices[i]))
    results[[d]][[i]] <- list()
    for (j in 1:length(nsamples)) {
      print(paste('!!!!! j:', j, '|| n:', nsamples[j]))
      DC_MCF <- bal_binary_fusion_SMC_BLR(N_schedule = rep(nsamples[j], 4),
                                          m_schedule = rep(2, 4),
                                          time_schedule = rep(time_choices[i], 4),
                                          base_samples = sub_posteriors_16[[d]],
                                          L = 5,
                                          dim = dim[d],
                                          data_split = data_split_16[[d]],
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
                                          print_progress_iters = 500)
      samples <- resample_particle_y_samples(particle_set = DC_MCF$particles[[1]],
                                             multivariate = TRUE,
                                             resampling_method = resampling_method,
                                             seed = seed)$y_samples
      results[[d]][[i]][[j]] <- list('IAD' = integrated_abs_distance(full_post = full_posterior[[d]],
                                                                     fusion_post = samples),
                                     'time' = sum(unlist(DC_MCF$time)),
                                     'samples' = samples,
                                     'CESS' = DC_MCF$CESS,
                                     'ESS' = DC_MCF$ESS,
                                     'resampled' = DC_MCF$resampled,
                                     'phi_bound_intensity' = DC_MCF$phi_bound_intensity,
                                     'phi_kappa' = DC_MCF$phi_kappa)
      print(paste('IAD:', results[[d]][[i]][[j]]$IAD))
      print(paste('time:', results[[d]][[i]][[j]]$time))
      print(paste('log(time):', log(results[[d]][[i]][[j]]$time)))
      print('save_progress')
      save.image('SD16_dim_study_varying.RData')
    }
  }
}

save.image('SD16_dim_study_varying.RData')
