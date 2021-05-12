library(hierarchicalFusion)
library(HMCBLR)

seed <- 2016

load("~/OneDrive - The Alan Turing Institute/freefor7/2021-02-12_hierarchical_fusion/results/logistic_regression/credit_cards/credit_cards_unbalanced_sub_posterior_sampling.RData")

time_choices <- c(0.25, 0.5, 1)
test_preconditioned_hierarchical_SMC_Poisson_hc <- list()
test_preconditioned_hierarchical_SMC_NB_hc <- list()
for (i in 1:length(time_choices)) {
  print(paste('Time: ', time_choices[i]))
  ##### Poisson (Hypercube Centre) #####
  print('Poisson Fusion (hypercube centre)')
  test_preconditioned_hierarchical_SMC_Poisson_hc[[i]] <- hierarchical_fusion_SMC_BLR(N_schedule = rep(30000, 4),
                                                                                      m_schedule = rep(2, 4),
                                                                                      time_schedule = rep(time_choices[i], 4),
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
                                                                                      print_progress_iters = 625)
  compare_samples_bivariate(posteriors = list(full_posterior,
                                              test_preconditioned_hierarchical_SMC_Poisson_hc[[i]]$proposed_samples[[1]],
                                              test_preconditioned_hierarchical_SMC_Poisson_hc[[i]]$particles[[1]]$y_samples),
                            colours = c('black', 'darkgreen', 'red'),
                            common_limit = c(-4, 4),
                            title = paste('Credit Cards - C=16 || SMC Hierarchical [Poisson] (h.c.) || Time =', time_choices[i]))
  test_preconditioned_hierarchical_SMC_Poisson_hc[[i]]$particles[[1]] <- resample_particle_y_samples(particle_set = test_preconditioned_hierarchical_SMC_Poisson_hc[[i]]$particles[[1]],
                                                                                                     multivariate = TRUE,
                                                                                                     resampling_method = 'resid',
                                                                                                     seed = seed)
  ##### NB (Hypercube Centre) #####
  print('NB Fusion (hypercube centre)')
  test_preconditioned_hierarchical_SMC_NB_hc[[i]] <- hierarchical_fusion_SMC_BLR(N_schedule = rep(30000, 4),
                                                                                 m_schedule = rep(2, 4),
                                                                                 time_schedule = rep(time_choices[i], 4),
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
                                                                                 print_progress_iters = 625)
  test_preconditioned_hierarchical_SMC_NB_hc[[i]]$particles[[1]] <- resample_particle_y_samples(particle_set = test_preconditioned_hierarchical_SMC_NB_hc[[i]]$particles[[1]],
                                                                                                multivariate = TRUE,
                                                                                                resampling_method = 'resid',
                                                                                                seed = seed)
  compare_samples_bivariate(posteriors = list(full_posterior,
                                              test_preconditioned_hierarchical_SMC_NB_hc[[i]]$proposed_samples[[1]],
                                              test_preconditioned_hierarchical_SMC_NB_hc[[i]]$particles[[1]]$y_samples),
                            colours = c('black', 'darkgreen', 'blue'),
                            common_limit = c(-4, 4),
                            title = paste('Credit Cards - C=16 || SMC Hierarchical [NB] || Time =', time_choices[i]))
}

for (i in 1:length(time_choices)) {
  print(paste('Time: ', time_choices[i]))
  print('Poisson Fusion (h.c.)')
  print(integrated_abs_distance(full_posterior,
                                test_preconditioned_hierarchical_SMC_Poisson_hc[[i]]$particles[[1]]$y_samples,
                                bandwidths))
  print('NB Fusion (h.c.)')
  print(integrated_abs_distance(full_posterior,
                                test_preconditioned_hierarchical_SMC_NB_hc[[i]]$particles[[1]]$y_samples,
                                bandwidths))
}

