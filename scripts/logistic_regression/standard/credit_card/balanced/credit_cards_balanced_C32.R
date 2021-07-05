library(hierarchicalFusion)
library(HMCBLR)

seed <- 2016

load("balanced_C32.RData")

time_choices <- c(0.25, 0.5, 1, 1.5)
Poisson_hc_32 <- list()
NB_hc_32 <- list()
for (i in 2:length(time_choices)) {
  print(paste('Time: ', time_choices[i]))
  ##### Poisson (Hypercube Centre) #####
  print('Poisson Fusion (hypercube centre)')
  Poisson_hc_32[[i]] <- hierarchical_fusion_SMC_BLR(N_schedule = rep(30000, 5),
                                                    m_schedule = rep(2, 5),
                                                    time_schedule = rep(time_choices[i], 5),
                                                    base_samples = sub_posteriors_32,
                                                    L = 6,
                                                    dim = 5,
                                                    data_split = data_split_32,
                                                    prior_means = rep(0, 5),
                                                    prior_variances = rep(1, 5),
                                                    C = 32,
                                                    precondition = TRUE,
                                                    resampling_method = 'resid',
                                                    ESS_threshold = 0.5,
                                                    cv_location = 'hypercube_centre',
                                                    diffusion_estimator = 'Poisson',
                                                    seed = seed,
                                                    create_cluster = FALSE,
                                                    print_progress_iters = 1250)
  Poisson_hc_32[[i]]$particles <- resample_particle_y_samples(particle_set = Poisson_hc_32[[i]]$particles[[1]],
                                                              multivariate = TRUE,
                                                              resampling_method = 'resid',
                                                              seed = seed)
  Poisson_hc_32[[i]]$proposed_samples <- Poisson_hc_32[[i]]$proposed_samples[[1]]
  print(integrated_abs_distance(full_posterior,
                                Poisson_hc_32[[i]]$particles$y_samples,
                                bandwidths))
  compare_samples_bivariate(posteriors = list(full_posterior,
                                              Poisson_hc_32[[i]]$proposed_samples,
                                              Poisson_hc_32[[i]]$particles$y_samples),
                            colours = c('black', 'darkgreen', 'red'),
                            common_limit = c(-4, 4),
                            title = paste('Credit Cards - C=32 || SMC Hierarchical [Poisson] (h.c.) || Time =', time_choices[i]))
  ##### NB (Hypercube Centre) #####
  print('NB Fusion (hypercube centre)')
  NB_hc_32[[i]] <- hierarchical_fusion_SMC_BLR(N_schedule = rep(30000, 5),
                                               m_schedule = rep(2, 5),
                                               time_schedule = rep(time_choices[i], 5),
                                               base_samples = sub_posteriors_32,
                                               L = 6,
                                               dim = 5,
                                               data_split = data_split_32,
                                               prior_means = rep(0, 5),
                                               prior_variances = rep(1, 5),
                                               C = 32,
                                               precondition = TRUE,
                                               resampling_method = 'resid',
                                               ESS_threshold = 0.5,
                                               cv_location = 'hypercube_centre',
                                               diffusion_estimator = 'NB',
                                               seed = seed,
                                               create_cluster = FALSE,
                                               print_progress_iters = 1250)
  NB_hc_32[[i]]$particles <- resample_particle_y_samples(particle_set = NB_hc_32[[i]]$particles[[1]],
                                                         multivariate = TRUE,
                                                         resampling_method = 'resid',
                                                         seed = seed)
  NB_hc_32[[i]]$proposed_samples <- NB_hc_32[[i]]$proposed_samples[[1]]
  print(integrated_abs_distance(full_posterior,
                                NB_hc_32[[i]]$particles$y_samples,
                                bandwidths))
  compare_samples_bivariate(posteriors = list(full_posterior,
                                              NB_hc_32[[i]]$proposed_samples,
                                              NB_hc_32[[i]]$particles$y_samples),
                            colours = c('black', 'darkgreen', 'blue'),
                            common_limit = c(-4, 4),
                            title = paste('Credit Cards - C=32 || SMC Hierarchical [NB] || Time =', time_choices[i]))
}

print('Poisson Fusion (h.c.)')
print(sapply(1:length(time_choices), function(i) {
  integrated_abs_distance(full_posterior,
                          Poisson_hc_32[[i]]$particles$y_samples)}))
print('NB Fusion (h.c.)')
print(sapply(1:length(time_choices), function(i) {
  integrated_abs_distance(full_posterior,
                          NB_hc_32[[i]]$particles$y_samples)}))

integrated_abs_distance(full_posterior, consensus_mat_32$samples)
integrated_abs_distance(full_posterior, consensus_sca_32$samples)
integrated_abs_distance(full_posterior, neiswanger_32_true$samples)
integrated_abs_distance(full_posterior, neiswanger_32_false$samples)
integrated_abs_distance(full_posterior, weierstrass_32_importance$samples)
integrated_abs_distance(full_posterior, weierstrass_32_rejection$samples)

