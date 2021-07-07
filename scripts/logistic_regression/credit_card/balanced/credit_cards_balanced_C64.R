library(hierarchicalFusion)
library(HMCBLR)

seed <- 2016

load("balanced_C64.RData")

time_choices <- c(0.5, 1)
Poisson_hc_64 <- list()
NB_hc_64 <- list()
for (i in 1:length(time_choices)) {
  print(paste('Time: ', time_choices[i]))
  ##### Poisson (Hypercube Centre) #####
  print('Poisson Fusion (hypercube centre)')
  Poisson_hc_64[[i]] <- hierarchical_fusion_SMC_BLR(N_schedule = rep(30000, 6),
                                                    m_schedule = rep(2, 6),
                                                    time_schedule = rep(time_choices[i], 6),
                                                    base_samples = sub_posteriors_64,
                                                    L = 7,
                                                    dim = 5,
                                                    data_split = data_split_64,
                                                    prior_means = rep(0, 5),
                                                    prior_variances = rep(1, 5),
                                                    C = 64,
                                                    precondition = TRUE,
                                                    resampling_method = 'resid',
                                                    ESS_threshold = 0.5,
                                                    cv_location = 'hypercube_centre',
                                                    diffusion_estimator = 'Poisson',
                                                    seed = seed,
                                                    print_progress_iters = 1250)
  Poisson_hc_64[[i]]$particles <- resample_particle_y_samples(particle_set = Poisson_hc_64[[i]]$particles[[1]],
                                                                   multivariate = TRUE,
                                                                   resampling_method = 'resid',
                                                                   seed = seed)
  Poisson_hc_64[[i]]$proposed_samples <- Poisson_hc_64[[i]]$proposed_samples[[1]]
  print(integrated_abs_distance(full_posterior,
                                Poisson_hc_64[[i]]$particles$y_samples,
                                bandwidths))
  compare_samples_bivariate(posteriors = list(full_posterior,
                                              Poisson_hc_64[[i]]$proposed_samples,
                                              Poisson_hc_64[[i]]$particles$y_samples),
                            colours = c('black', 'darkgreen', 'red'),
                            common_limit = c(-4, 4),
                            title = paste('Credit Cards - C=64 || SMC Hierarchical [Poisson] (h.c.) || Time =', time_choices[i]))
  ##### NB (Hypercube Centre) #####
  print('NB Fusion (hypercube centre)')
  NB_hc_64[[i]] <- hierarchical_fusion_SMC_BLR(N_schedule = rep(30000, 6),
                                               m_schedule = rep(2, 6),
                                               time_schedule = rep(time_choices[i], 6),
                                               base_samples = sub_posteriors_64,
                                               L = 7,
                                               dim = 5,
                                               data_split = data_split_64,
                                               prior_means = rep(0, 5),
                                               prior_variances = rep(1, 5),
                                               C = 64,
                                               precondition = TRUE,
                                               resampling_method = 'resid',
                                               ESS_threshold = 0.5,
                                               cv_location = 'hypercube_centre',
                                               diffusion_estimator = 'NB',
                                               seed = seed,
                                               create_cluster = FALSE,
                                               print_progress_iters = 1250)
  NB_hc_64[[i]]$particles <- resample_particle_y_samples(particle_set = NB_hc_64[[i]]$particles[[1]],
                                                              multivariate = TRUE,
                                                              resampling_method = 'resid',
                                                              seed = seed)
  NB_hc_64[[i]]$proposed_samples <- NB_hc_64[[i]]$proposed_samples[[1]]
  print(integrated_abs_distance(full_posterior,
                                NB_hc_64[[i]]$particles$y_samples,
                                bandwidths))
  compare_samples_bivariate(posteriors = list(full_posterior,
                                              NB_hc_64[[i]]$proposed_samples,
                                              NB_hc_64[[i]]$particles$y_samples),
                            colours = c('black', 'darkgreen', 'blue'),
                            common_limit = c(-4, 4),
                            title = paste('Credit Cards - C=64 || SMC Hierarchical [NB] || Time =', time_choices[i]))
}

print('Poisson Fusion (h.c.)')
print(sapply(1:length(time_choices), function(i) {
  integrated_abs_distance(full_posterior,
                          Poisson_hc_64[[i]]$particles$y_samples)}))
print('NB Fusion (h.c.)')
print(sapply(1:length(time_choices), function(i) {
  integrated_abs_distance(full_posterior,
                          NB_hc_64[[i]]$particles$y_samples)}))

integrated_abs_distance(full_posterior, consensus_mat_64$samples)
integrated_abs_distance(full_posterior, consensus_sca_64$samples)
integrated_abs_distance(full_posterior, neiswanger_64_true$samples)
integrated_abs_distance(full_posterior, neiswanger_64_false$samples)
integrated_abs_distance(full_posterior, weierstrass_64_importance$samples)
integrated_abs_distance(full_posterior, weierstrass_64_rejection$samples)
