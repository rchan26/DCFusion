library(hierarchicalFusion)

seed <- 2016

load('balanced_C128.RData')

time_choice <- 0.5
NB_hc_128 <- list()
n_samples <- c(1000, 10000, 50000)
for (i in 1:length(n_samples)) {
  print(paste('Samples:', n_samples[i]))
  NB_hc_128[[i]] <- list()
  for (j in 1:10) {
    print(paste('j:', j))
    NB_hc_128[[i]][[j]] <- hierarchical_fusion_SMC_BLR(N_schedule = rep(n_samples[i], 7),
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
                                                       seed = seed*i*j,
                                                       print_progress_iters = 1000)
    NB_hc_128[[i]][[j]]$particles <- resample_particle_y_samples(particle_set = NB_hc_128[[i]][[j]]$particles[[1]],
                                                                 multivariate = TRUE,
                                                                 resampling_method = 'resid',
                                                                 seed = seed)
    NB_hc_128[[i]][[j]]$proposed_samples <- NB_hc_128[[i]][[j]]$proposed_samples[[1]]
    print(integrated_abs_distance(full_posterior,
                                  NB_hc_128[[i]][[j]]$particles$y_samples,
                                  bandwidths))
    print('time_taken:'); print(sum(unlist(NB_hc_128[[i]][[j]]$time))/60/60)
    save.image('C128_replicates.RData')
  }
}

IAD <- list()
for (i in 1:length(n_samples)) {
  IAD[[i]] <- sapply(1:10, function(j) {
    integrated_abs_distance(full_posterior,
                            NB_hc_128[[i]][[j]]$particles$y_samples,
                            bandwidths)})
}

time_taken <- list()
for (i in 1:length(n_samples)) {
  time_taken[[i]] <- sapply(1:10, function(j) sum(unlist(NB_hc_128[[i]][[j]]$time)))
}
