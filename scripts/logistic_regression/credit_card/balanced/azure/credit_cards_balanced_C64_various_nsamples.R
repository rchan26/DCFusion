library(hierarchicalFusion)

seed <- 2016

load('balanced_C64.RData')

time_choice <- 0.5
Poisson_hc_64 <- list()
NB_hc_64 <- list()
n_cores <- parallel::detectCores()
n_samples <- n_cores*c(10, 50, 100, 200, 300, 400, 500, 1000, 2000, 2500, 3000)
for (i in 11:length(n_samples)) {
  print(paste('Samples:', n_samples[i]))
  ##### Poisson (Hypercube Centre) #####
  print('Poisson Fusion (hypercube centre)')
  Poisson_hc_64[[i]] <- hierarchical_fusion_SMC_BLR(N_schedule = rep(n_samples[i], 6),
                                                    m_schedule = rep(2, 6),
                                                    time_schedule = rep(time_choice, 6),
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
                                                    print_progress_iters = 500)
  Poisson_hc_64[[i]]$particles <- resample_particle_y_samples(particle_set = Poisson_hc_64[[i]]$particles[[1]],
                                                              multivariate = TRUE,
                                                              resampling_method = 'resid',
                                                              seed = seed)
  Poisson_hc_64[[i]]$proposed_samples <- Poisson_hc_64[[i]]$proposed_samples[[1]]
  print(integrated_abs_distance(full_posterior,
                                Poisson_hc_64[[i]]$particles$y_samples,
                                bandwidths))
  print(sum(unlist(Poisson_hc_64[[i]]$time))/60/60)
  
  ##### NB (Hypercube Centre) #####
  print('NB Fusion (hypercube centre)')
  NB_hc_64[[i]] <- hierarchical_fusion_SMC_BLR(N_schedule = rep(n_samples[i], 6),
                                               m_schedule = rep(2, 6),
                                               time_schedule = rep(time_choice, 6),
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
                                               print_progress_iters = 500)
  NB_hc_64[[i]]$particles <- resample_particle_y_samples(particle_set = NB_hc_64[[i]]$particles[[1]],
                                                         multivariate = TRUE,
                                                         resampling_method = 'resid',
                                                         seed = seed)
  NB_hc_64[[i]]$proposed_samples <- NB_hc_64[[i]]$proposed_samples[[1]]
  print(integrated_abs_distance(full_posterior,
                                NB_hc_64[[i]]$particles$y_samples,
                                bandwidths))
  print(sum(unlist(NB_hc_64[[i]]$time))/60/60)
  
  save.image('credit_cards_balanced_C64_various_nsamples_v2.RData')
}

##### Time taken (in hours) #####

print('Samples:'); print(n_samples)
print('Poisson time in hours:')
print(sapply(1:length(n_samples), function(i) sum(unlist(Poisson_hc_64[[i]]$time))))
print('NB time in hours:')
print(sapply(1:length(n_samples), function(i) sum(unlist(NB_hc_64[[i]]$time))))

##### IAD #####

print('Poisson Fusion (h.c.)')
print(sapply(1:length(n_samples), function(i) {
  integrated_abs_distance(full_posterior, Poisson_hc_64[[i]]$particles$y_samples)}))
print('NB Fusion (h.c.)')
print(sapply(1:length(n_samples), function(i) {
  integrated_abs_distance(full_posterior, NB_hc_64[[i]]$particles$y_samples)}))

compare_samples_bivariate(posteriors = list(full_posterior,
                                            Poisson_hc_64[[14]]$particles$y_samples,
                                            NB_hc_64[[14]]$particles$y_samples),
                          colours = c('black', 'darkgreen', 'blue'),
                          common_limit = c(-4, 4),
                          title = paste('Credit Cards - C=32 || SMC Hierarchical [NB] || N = 100000'))

save.image('credit_cards_balanced_C64_various_nsamples_v2.RData')

integrated_abs_distance(full_posterior, consensus_mat_64$samples)
integrated_abs_distance(full_posterior, consensus_sca_64$samples)
integrated_abs_distance(full_posterior, neiswanger_64_true$samples)
integrated_abs_distance(full_posterior, neiswanger_64_false$samples)
integrated_abs_distance(full_posterior, weierstrass_64_importance$samples)
integrated_abs_distance(full_posterior, weierstrass_64_rejection$samples)
