library(hierarchicalFusion)

seed <- 2016

load('balanced_C128.RData')

time_choice <- 0.5
Poisson_hc_128 <- list()
NB_hc_128 <- list()
n_cores <- parallel::detectCores()
n_samples <- n_cores*c(10, 50, 100, 200, 300, 400, 500, 1000, 2000, 2500, 3000)
for (i in 10:length(n_samples)) {
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

# ##### Time taken (in hours) #####
# 
# print('Samples:'); print(n_samples)
# print('Poisson time in hours:')
# print(sapply(1:length(n_samples), function(i) sum(unlist(Poisson_hc_128[[i]]$time))/60/60))
# print('NB time in hours:')
# print(sapply(1:length(n_samples), function(i) sum(unlist(NB_hc_128[[i]]$time))/60/60))
# 
# ##### IAD #####
# 
# print('Poisson Fusion (h.c.)')
# print(sapply(1:length(n_samples), function(i) {
#   integrated_abs_distance(full_posterior, Poisson_hc_128[[i]]$particles$y_samples)}))
# print('NB Fusion (h.c.)')
# print(sapply(1:length(n_samples), function(i) {
#   integrated_abs_distance(full_posterior, NB_hc_128[[i]]$particles$y_samples)}))
# 
# compare_samples_bivariate(posteriors = list(full_posterior,
#                                             Poisson_hc_128[[14]]$particles$y_samples,
#                                             NB_hc_128[[14]]$particles$y_samples),
#                           colours = c('black', 'darkgreen', 'blue'),
#                           common_limit = c(-4, 4),
#                           title = paste('Credit Cards - C=32 || SMC Hierarchical [NB] || N = 100000'))
# 
# save.image('credit_cards_balanced_C128_various_nsamples.RData')
# 
# integrated_abs_distance(full_posterior, consensus_mat_128$samples)
# integrated_abs_distance(full_posterior, consensus_sca_128$samples)
# integrated_abs_distance(full_posterior, neiswanger_128_true$samples)
# integrated_abs_distance(full_posterior, neiswanger_128_false$samples)
# integrated_abs_distance(full_posterior, weierstrass_128_importance$samples)
# integrated_abs_distance(full_posterior, weierstrass_128_rejection$samples)
