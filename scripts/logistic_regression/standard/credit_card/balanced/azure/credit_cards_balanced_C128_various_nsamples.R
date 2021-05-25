library(hierarchicalFusion)

seed <- 2016

load('credit_cards_balanced_sub_posteriors.RData')

time_choice <- 0.5
test_preconditioned_hierarchical_SMC_Poisson_hc <- list()
test_preconditioned_hierarchical_SMC_NB_hc <- list()
n_samples <- c(2500, 5000, 7500, 10000, 15000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000)
for (i in 1:length(n_samples)) {
  print(paste('Samples:', n_samples[i]))
  ##### Poisson (Hypercube Centre) #####
  print('Poisson Fusion (hypercube centre)')
  test_preconditioned_hierarchical_SMC_Poisson_hc[[i]] <- hierarchical_fusion_SMC_BLR(N_schedule = rep(n_samples[i], 7),
                                                                                      m_schedule = rep(2, 7),
                                                                                      time_schedule = rep(time_choice, 7),
                                                                                      base_samples = sub_posteriors_128,
                                                                                      L = 8,
                                                                                      dim = 5,
                                                                                      data_split = data_split_128,
                                                                                      prior_means = rep(0, 5),
                                                                                      prior_variances = rep(1, 5),
                                                                                      C = 128,
                                                                                      precondition = TRUE,
                                                                                      resampling_method = 'resid',
                                                                                      ESS_threshold = 0.5,
                                                                                      cv_location = 'hypercube_centre',
                                                                                      diffusion_estimator = 'Poisson',
                                                                                      seed = seed,
                                                                                      print_progress_iters = 600)
  test_preconditioned_hierarchical_SMC_Poisson_hc[[i]]$particles <- resample_particle_y_samples(particle_set = test_preconditioned_hierarchical_SMC_Poisson_hc[[i]]$particles[[1]],
                                                                                                multivariate = TRUE,
                                                                                                resampling_method = 'resid',
                                                                                                seed = seed)
  test_preconditioned_hierarchical_SMC_Poisson_hc[[i]]$proposed_samples <- test_preconditioned_hierarchical_SMC_Poisson_hc[[i]]$proposed_samples[[1]]
  print(integrated_abs_distance(full_posterior,
                                test_preconditioned_hierarchical_SMC_Poisson_hc[[i]]$particles$y_samples,
                                bandwidths))
  ##### NB (Hypercube Centre) #####
  print('NB Fusion (hypercube centre)')
  test_preconditioned_hierarchical_SMC_NB_hc[[i]] <- hierarchical_fusion_SMC_BLR(N_schedule = rep(n_samples[i], 7),
                                                                                 m_schedule = rep(2, 7),
                                                                                 time_schedule = rep(time_choice, 7),
                                                                                 base_samples = sub_posteriors_128,
                                                                                 L = 8,
                                                                                 dim = 5,
                                                                                 data_split = data_split_128,
                                                                                 prior_means = rep(0, 5),
                                                                                 prior_variances = rep(1, 5),
                                                                                 C = 128,
                                                                                 precondition = TRUE,
                                                                                 resampling_method = 'resid',
                                                                                 ESS_threshold = 0.5,
                                                                                 cv_location = 'hypercube_centre',
                                                                                 diffusion_estimator = 'NB',
                                                                                 seed = seed,
                                                                                 print_progress_iters = 600)
  test_preconditioned_hierarchical_SMC_NB_hc[[i]]$particles <- resample_particle_y_samples(particle_set = test_preconditioned_hierarchical_SMC_NB_hc[[i]]$particles[[1]],
                                                                                           multivariate = TRUE,
                                                                                           resampling_method = 'resid',
                                                                                           seed = seed)
  print(integrated_abs_distance(full_posterior,
                                test_preconditioned_hierarchical_SMC_NB_hc[[i]]$particles$y_samples,
                                bandwidths))
  test_preconditioned_hierarchical_SMC_NB_hc[[i]]$proposed_samples <- test_preconditioned_hierarchical_SMC_NB_hc[[i]]$proposed_samples[[1]]
  save.image('credit_cards_balanced_C128_various_nsamples.RData')
}

save.image('credit_cards_balanced_C128_various_nsamples.RData')
