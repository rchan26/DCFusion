library(hierarchicalFusion)
library(HMCBLR)

seed <- 2016

load("~/OneDrive - The Alan Turing Institute/freefor7/2021-02-12_hierarchical_fusion/results/logistic_regression/credit_cards/credit_cards_sub_posterior_sampling.RData")

time_choices <- c(0.25, 0.5, 1)
test_preconditioned_hierarchical_SMC_Poisson_mode <- list()
test_preconditioned_hierarchical_SMC_Poisson_hc <- list()
test_preconditioned_hierarchical_SMC_NB_mode <- list()
test_preconditioned_hierarchical_SMC_NB_hc <- list()
for (i in 1:length(time_choices)) {
  print(paste('Time: ', time_choices[i]))
  ##### Poisson (Mode) #####
  print('Poisson Fusion (mode)')
  test_preconditioned_hierarchical_SMC_Poisson_mode[[i]] <- hierarchical_fusion_SMC_BLR(N_schedule = rep(30000, 4),
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
                                                                                        cv_location = 'mode',
                                                                                        diffusion_estimator = 'Poisson',
                                                                                        seed = seed,
                                                                                        print_progress_iters = 1250)
  test_preconditioned_hierarchical_SMC_Poisson_mode[[i]]$particles[[1]] <- resample_particle_y_samples(particle_set = test_preconditioned_hierarchical_SMC_Poisson_mode[[i]]$particles[[1]],
                                                                                                       multivariate = TRUE,
                                                                                                       resampling_method = 'resid',
                                                                                                       seed = seed)
  compare_samples_bivariate(posteriors = list(full_posterior,
                                              test_preconditioned_hierarchical_SMC_Poisson_mode[[i]]$proposed_samples[[1]],
                                              test_preconditioned_hierarchical_SMC_Poisson_mode[[i]]$particles[[1]]$y_samples),
                            colours = c('black', 'darkgreen', 'red'),
                            common_limit = c(-4, 4),
                            title = paste('Credit Cards - C=16 || SMC Hierarchical [Poisson] (mode) || Time =', time_choices[i]))
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
                                                                                      print_progress_iters = 1250)
  test_preconditioned_hierarchical_SMC_Poisson_hc[[i]]$particles[[1]] <- resample_particle_y_samples(particle_set = test_preconditioned_hierarchical_SMC_Poisson_hc[[i]]$particles[[1]],
                                                                                                     multivariate = TRUE,
                                                                                                     resampling_method = 'resid',
                                                                                                     seed = seed)
  compare_samples_bivariate(posteriors = list(full_posterior,
                                              test_preconditioned_hierarchical_SMC_Poisson_hc[[i]]$proposed_samples[[1]],
                                              test_preconditioned_hierarchical_SMC_Poisson_hc[[i]]$particles[[1]]$y_samples),
                            colours = c('black', 'darkgreen', 'red'),
                            common_limit = c(-4, 4),
                            title = paste('Credit Cards - C=16 || SMC Hierarchical [Poisson] (h.c.) || Time =', time_choices[i]))
  ##### NB (Mode) #####
  print('NB Fusion (mode)')
  test_preconditioned_hierarchical_SMC_NB_mode[[i]] <- hierarchical_fusion_SMC_BLR(N_schedule = rep(30000, 4),
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
                                                                                   cv_location = 'mode',
                                                                                   diffusion_estimator = 'NB',
                                                                                   seed = seed,
                                                                                   print_progress_iters = 1250)
  test_preconditioned_hierarchical_SMC_NB_mode[[i]]$particles[[1]] <- resample_particle_y_samples(particle_set = test_preconditioned_hierarchical_SMC_NB_mode[[i]]$particles[[1]],
                                                                                                  multivariate = TRUE,
                                                                                                  resampling_method = 'resid',
                                                                                                  seed = seed)
  compare_samples_bivariate(posteriors = list(full_posterior,
                                              test_preconditioned_hierarchical_SMC_NB_mode[[i]]$proposed_samples[[1]],
                                              test_preconditioned_hierarchical_SMC_NB_mode[[i]]$particles[[1]]$y_samples),
                            colours = c('black', 'darkgreen', 'blue'),
                            common_limit = c(-4, 4),
                            title = paste('Credit Cards - C=16 || SMC Hierarchical [NB] || Time =', time_choices[i]))
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
                                                                                 print_progress_iters = 1250)
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

bandwidths <- rep(0.05, 5)
for (i in 1:length(time_choices)) {
  print(paste('Time: ', time_choices[i]))
  print('Poisson Fusion (mode)')
  print(integrated_abs_distance(full_posterior,
                                test_preconditioned_hierarchical_SMC_Poisson_mode[[i]]$particles[[1]]$y_samples,
                                bandwidths))
  print('NB Fusion (mode)')
  print(integrated_abs_distance(full_posterior,
                                test_preconditioned_hierarchical_SMC_NB_mode[[i]]$particles[[1]]$y_samples,
                                bandwidths))
  print('Poisson Fusion (h.c.)')
  print(integrated_abs_distance(full_posterior,
                                test_preconditioned_hierarchical_SMC_Poisson_hc[[i]]$particles[[1]]$y_samples,
                                bandwidths))
  print('NB Fusion (h.c.)')
  print(integrated_abs_distance(full_posterior,
                                test_preconditioned_hierarchical_SMC_NB_hc[[i]]$particles[[1]]$y_samples,
                                bandwidths))
}

integrated_abs_distance(full_posterior, consensus_mat_16$samples, bandwidths)
integrated_abs_distance(full_posterior, consensus_sca_16$samples, bandwidths)
integrated_abs_distance(full_posterior, neiswanger_16_true$samples, bw = bandwidths)
integrated_abs_distance(full_posterior, neiswanger_16_false$samples, bw = bandwidths)
integrated_abs_distance(full_posterior, weierstrass_16_importance$samples, bandwidths)
integrated_abs_distance(full_posterior, weierstrass_16_rejection$samples, bandwidths)
      
