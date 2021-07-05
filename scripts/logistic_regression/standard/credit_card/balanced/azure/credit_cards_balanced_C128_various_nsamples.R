library(hierarchicalFusion)

seed <- 2016

load('balanced_C128.RData')

time_choice <- 0.5
Poisson_hc_128 <- list()
NB_hc_128 <- list()
n_samples <- c(2500, 5000, 7500, 10000, 15000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000)
for (i in 1:length(n_samples)) {
  print(paste('Samples:', n_samples[i]))
  ##### Poisson (Hypercube Centre) #####
  print('Poisson Fusion (hypercube centre)')
  Poisson_hc_128[[i]] <- hierarchical_fusion_SMC_BLR(N_schedule = rep(n_samples[i], 7),
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
                                                     create_cluster = TRUE,
                                                     print_progress_iters = 600)
  Poisson_hc_128[[i]]$particles <- resample_particle_y_samples(particle_set = Poisson_hc_128[[i]]$particles[[1]],
                                                               multivariate = TRUE,
                                                               resampling_method = 'resid',
                                                               seed = seed)
  Poisson_hc_128[[i]]$proposed_samples <- Poisson_hc_128[[i]]$proposed_samples[[1]]
  print(integrated_abs_distance(full_posterior,
                                Poisson_hc_128[[i]]$particles$y_samples,
                                bandwidths))
  ##### NB (Hypercube Centre) #####
  print('NB Fusion (hypercube centre)')
  NB_hc_128[[i]] <- hierarchical_fusion_SMC_BLR(N_schedule = rep(n_samples[i], 7),
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
                                                create_cluster = TRUE,
                                                print_progress_iters = 600)
  NB_hc_128[[i]]$particles <- resample_particle_y_samples(particle_set = NB_hc_128[[i]]$particles[[1]],
                                                          multivariate = TRUE,
                                                          resampling_method = 'resid',
                                                          seed = seed)
  NB_hc_128[[i]]$proposed_samples <- NB_hc_128[[i]]$proposed_samples[[1]]
  print(integrated_abs_distance(full_posterior,
                                NB_hc_128[[i]]$particles$y_samples,
                                bandwidths))
  
  save.image('credit_cards_balanced_C128_various_nsamples_v2.RData')
}

##### Time taken (in hours) #####

print('Samples:'); print(n_samples)
print('Poisson time in hours:')
print(sapply(1:length(n_samples), function(i) sum(unlist(Poisson_hc_128[[i]]$time))/60/60))
print('NB time in hours:')
print(sapply(1:length(n_samples), function(i) sum(unlist(NB_hc_128[[i]]$time))/60/60))

##### IAD #####

print('Poisson Fusion (h.c.)')
print(sapply(1:length(n_samples), function(i) {
  integrated_abs_distance(full_posterior, Poisson_hc_128[[i]]$particles$y_samples)}))
print('NB Fusion (h.c.)')
print(sapply(1:length(n_samples), function(i) {
  integrated_abs_distance(full_posterior, NB_hc_128[[i]]$particles$y_samples)}))

compare_samples_bivariate(posteriors = list(full_posterior,
                                            Poisson_hc_128[[14]]$particles$y_samples,
                                            NB_hc_128[[14]]$particles$y_samples),
                          colours = c('black', 'darkgreen', 'blue'),
                          common_limit = c(-4, 4),
                          title = paste('Credit Cards - C=32 || SMC Hierarchical [NB] || N = 100000'))

save.image('credit_cards_balanced_C128_various_nsamples_v2.RData')

integrated_abs_distance(full_posterior, consensus_mat_128$samples)
integrated_abs_distance(full_posterior, consensus_sca_128$samples)
integrated_abs_distance(full_posterior, neiswanger_128_true$samples)
integrated_abs_distance(full_posterior, neiswanger_128_false$samples)
integrated_abs_distance(full_posterior, weierstrass_128_importance$samples)
integrated_abs_distance(full_posterior, weierstrass_128_rejection$samples)
