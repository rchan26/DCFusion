library(DCFusion)
library(HMCBLR)

load('NYC64.RData')
# preparing lists to store full posterior and sub-posterior samples
fp <- full_posterior
sp64 <- sub_posteriors_64
full_posterior <- list()
sub_posteriors_64 <- list()
full_posterior[[8]] <- fp
sub_posteriors_64[[8]] <- sp64
rm(fp, sp64)
# preparing lists to store IAD and elapsed run-times
balanced_C64_results <- list('reg' = list(), 'adaptive' = list())
CMC_results <- list('mat' = list(), 'sca' = list())
KDEMC_results <- list('true' = list(), 'false' = list())
WRS_results <- list('importance' = list(), 'rejection' = list())
N_samples <- c(500,1000,2000,5000,10000,15000,20000,30000)
get_approximate_result <- function(result) {
  return(list('IAD' = integrated_abs_distance(full_posterior[[i]], result$samples),
              'IAD2' = integrated_abs_distance(full_posterior[[8]], result$samples),
              'time' = result$time))
}
get_fusion_result <- function(result) {
  return(list('CESS' = result$CESS,
              'ESS' = result$ESS,
              'IAD' = integrated_abs_distance(full_posterior[[i]], result$particles$y_samples),
              'IAD2' = integrated_abs_distance(full_posterior[[8]], result$particles$y_samples),
              'time' = sum(unlist(result$time))))
}
# note that we have already done N = N_samples[8] = 300000, so recording those
i <- 8
CMC_results$mat[[i]] <- get_approximate_result(consensus_mat_64)
CMC_results$sca[[i]] <- get_approximate_result(consensus_sca_64)
KDEMC_results$true[[i]] <- get_approximate_result(neiswanger_true_64)
KDEMC_results$false[[i]] <- get_approximate_result(neiswanger_false_64)
WRS_results$importance[[i]] <- get_approximate_result(weierstrass_importance_64)
WRS_results$rejection[[i]] <- get_approximate_result(weierstrass_rejection_64)
balanced_C64_results$reg[[i]] <- get_fusion_result(balanced_C64$reg)
balanced_C64_results$adaptive[[i]] <- get_fusion_result(balanced_C64$adaptive)
rm(balanced_C64)

for (i in 8) {
  print(paste('##### i:', i))
  ##### sampling from target and sub-posteriors #####

  print('sampling from posterior')
  full_posterior[[i]] <- hmc_sample_BLR(full_data_count = full_data_count,
                                        C = 1,
                                        prior_means = prior_means,
                                        prior_variances = prior_variances,
                                        iterations = N_samples[i] + 10000,
                                        warmup = 10000,
                                        chains = 1,
                                        seed = seed,
                                        output = T)
  print('sampling from sub-posteriors')
  sub_posteriors_64[[i]] <- hmc_base_sampler_BLR(nsamples = N_samples[i],
                                                 data_split = data_split_64,
                                                 C = C,
                                                 prior_means = prior_means,
                                                 prior_variances = prior_variances,
                                                 warmup = 10000,
                                                 seed = seed,
                                                 output = T)

  ##### Applying other methodologies #####

  print('Applying other methodologies')
  consensus_mat_64 <- consensus_scott(S = C, samples_to_combine = sub_posteriors_64[[i]], indep = F)
  consensus_sca_64 <- consensus_scott(S = C, samples_to_combine = sub_posteriors_64[[i]], indep = T)
  neiswanger_true_64 <- neiswanger(S = C,
                                   samples_to_combine = sub_posteriors_64[[i]],
                                   anneal = TRUE)
  neiswanger_false_64 <- neiswanger(S = C,
                                    samples_to_combine = sub_posteriors_64[[i]],
                                    anneal = FALSE)
  weierstrass_importance_64 <- weierstrass(Samples = sub_posteriors_64[[i]],
                                           method = 'importance')
  weierstrass_rejection_64 <- weierstrass(Samples = sub_posteriors_64[[i]],
                                          method = 'reject')

  CMC_results$mat[[i]] <- get_approximate_result(consensus_mat_64)
  CMC_results$sca[[i]] <- get_approximate_result(consensus_sca_64)
  KDEMC_results$true[[i]] <- get_approximate_result(neiswanger_true_64)
  KDEMC_results$false[[i]] <- get_approximate_result(neiswanger_false_64)
  WRS_results$importance[[i]] <- get_approximate_result(weierstrass_importance_64)
  WRS_results$rejection[[i]] <- get_approximate_result(weierstrass_rejection_64)

  rm(consensus_mat_64, consensus_sca_64,
     neiswanger_true_64, neiswanger_false_64,
     weierstrass_importance_64, weierstrass_rejection_64)
  save.image('NYC64_varying_N.RData')

  ##### bal binary combining two sub-posteriors at a time #####

  # regular mesh
  print('regular mesh')
  balanced_C64_reg <- bal_binary_GBF_BLR(N_schedule = rep(N_samples[i], 6),
                                                  m_schedule = rep(2, 6),
                                                  time_mesh = NULL,
                                                  base_samples = sub_posteriors_64[[i]],
                                                  L = 7,
                                                  dim = dim,
                                                  data_split = data_split_64,
                                                  prior_means = prior_means,
                                                  prior_variances = prior_variances,
                                                  C = C,
                                                  precondition = TRUE,
                                                  resampling_method = 'resid',
                                                  ESS_threshold = ESS_threshold,
                                                  adaptive_mesh = FALSE,
                                                  mesh_parameters = list('condition' = 'SH',
                                                                         'CESS_0_threshold' = CESS_0_threshold,
                                                                         'CESS_j_threshold' = CESS_j_threshold,
                                                                         'vanilla' = FALSE),
                                                  diffusion_estimator = diffusion_estimator,
                                                  seed = seed,
                                                  n_cores = n_cores,
                                                  print_progress_iters = 500)
  balanced_C64_reg$particles <- resample_particle_y_samples(particle_set = balanced_C64_reg$particles[[1]],
                                                            multivariate = TRUE,
                                                            resampling_method = 'resid',
                                                            seed = seed)
  balanced_C64_results$reg[[i]] <- get_fusion_result(balanced_C64_reg)
  rm(balanced_C64_reg)
  save.image('NYC64_varying_N.RData')

  # adaptive mesh
  print('adaptive mesh')
  balanced_C64_adpative <- bal_binary_GBF_BLR(N_schedule = rep(N_samples[i], 6),
                                              m_schedule = rep(2, 6),
                                              time_mesh = NULL,
                                              base_samples = sub_posteriors_64[[i]],
                                              L = 7,
                                              dim = dim,
                                              data_split = data_split_64,
                                              prior_means = prior_means,
                                              prior_variances = prior_variances,
                                              C = C,
                                              precondition = TRUE,
                                              resampling_method = 'resid',
                                              ESS_threshold = ESS_threshold,
                                              adaptive_mesh = TRUE,
                                              mesh_parameters = list('condition' = 'SH',
                                                                     'CESS_0_threshold' = CESS_0_threshold,
                                                                     'CESS_j_threshold' = CESS_j_threshold,
                                                                     'vanilla' = FALSE),
                                              diffusion_estimator = diffusion_estimator,
                                              seed = seed,n_cores = n_cores,
                                              print_progress_iters = 500)
  balanced_C64_adpative$particles <- resample_particle_y_samples(particle_set = balanced_C64_adpative$particles[[1]],
                                                                 multivariate = TRUE,
                                                                 resampling_method = 'resid',
                                                                 seed = seed)
  balanced_C64_results$adaptive[[i]] <- get_fusion_result(balanced_C64_adpative)
  rm(balanced_C64_adpative)
  save.image('NYC64_varying_N.RData')
}

plot(x = N_samples, y = sapply(1:8, function(i) balanced_C64_results$adaptive[[i]]$IAD2),
     ylim = c(0, 0.7),
     xlab = '',
     ylab = '',
     xaxt = 'n', lty = 2, lwd = 3, pch = 4, type = 'b')
mtext('N', 1, 2.75, font = 2, cex = 1.5)
mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=seq(0, 30000, 5000), labels = seq(0, 30000, 5000), font = 2, cex = 1.5)
axis(1, at=seq(0, 30000, 2500), labels=rep("", 13), lwd.ticks = 0.5)
axis(2, at=seq(0, 1, 0.1), labels=c("0.0", seq(0.1, 0.9, 0.1), "1.0"), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5)
lines(x = N_samples, y = sapply(1:8, function(i) balanced_C64_results$reg[[i]]$IAD2),
      lty = 3, lwd = 3, type = 'b', pch = 5)
lines(x = N_samples, sapply(1:8, function(i) CMC_results$mat[[i]]$IAD2),
      lty = 4, lwd = 3, type = 'b', pch = 3, col = 'red')
lines(x = N_samples, y = sapply(1:8, function(i) KDEMC_results$false[[i]]$IAD2),
      lty = 5, lwd = 3, type = 'b', pch = 2, col = 'red')
lines(x = N_samples, y = sapply(1:8, function(i) WRS_results$rejection[[i]]$IAD2),
      lty = 6, lwd = 3, type = 'b', pch = 1, col = 'red')
legend(x = 0, y = 0.7,
       legend = c('D&C-GBF (regular mesh)',
                  'D&C-GBF (adaptive mesh)',
                  'CMC',
                  'KDEMC',
                  'WRS'),
       lwd = rep(3, 6),
       lty = c(3,2,4,5,6),
       pch = c(5,4,3,2,1),
       col = c(rep('black', 2), rep('red', 3)),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

plot(x = N_samples, y = log(sapply(1:8, function(i) balanced_C64_results$adaptive[[i]]$time), 2),
     ylim = c(-2, 22),
     xlab = '',
     ylab = '',
     yaxt = 'n',
     xaxt = 'n', lty = 2, lwd = 3, pch = 4, type = 'b')
mtext('N', 1, 2.75, font = 2, cex = 1.5)
mtext('log(Time elapsed in seconds, 2)', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=seq(0, 30000, 5000), labels = seq(0, 30000, 5000), font = 2, cex = 1.5)
axis(1, at=seq(0, 30000, 2500), labels=rep("", 13), lwd.ticks = 0.5)
axis(2, at=seq(-4, 22, 2), labels = seq(-4, 22, 2), font = 2, cex = 1.5)
axis(2, at=seq(-4, 22, 1), labels=rep("", 27), lwd.ticks = 0.5)
lines(x = N_samples, y = log(sapply(1:8, function(i) balanced_C64_results$reg[[i]]$time), 2),
      lty = 3, lwd = 3, type = 'b', pch = 5)
lines(x = N_samples, y = log(sapply(1:8, function(i) CMC_results$mat[[i]]$time), 2),
      lty = 4, lwd = 3, type = 'b', pch = 3, col = 'red')
lines(x = N_samples, y = log(sapply(1:8, function(i) KDEMC_results$false[[i]]$time), 2),
      lty = 5, lwd = 3, type = 'b', pch = 2, col = 'red')
lines(x = N_samples, y = log(sapply(1:8, function(i) WRS_results$rejection[[i]]$time), 2),
      lty = 6, lwd = 3, type = 'b', pch = 1, col = 'red')
legend(x = 2, y = 22,
       legend = c('D&C-GBF (regular mesh)',
                  'D&C-GBF (adaptive mesh)',
                  'CMC',
                  'KDEMC',
                  'WRS'),
       lwd = rep(3, 6),
       lty = c(3,2,4,5,6),
       pch = c(5,4,3,2,1),
       col = c(rep('black', 2), rep('red', 3)),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

plot(x = log(sapply(1:8, function(i) balanced_C64_results$adaptive[[i]]$time), 2),
     y = sapply(1:8, function(i) balanced_C64_results$adaptive[[i]]$IAD2),
     ylim = c(0, 0.7),
     xlim = c(-2, 18),
     xlab = '',
     ylab = '',
     xaxt = 'n', lty = 2, lwd = 3, pch = 4, type = 'b')
mtext('log(Time elapsed in seconds, 2)', 1, 2.75, font = 2, cex = 1.5)
mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=seq(-4, 22, 2), labels = seq(-4, 22, 2), font = 2, cex = 1.5)
axis(1, at=seq(-4, 22, 1), labels=rep("", 27), lwd.ticks = 0.5)
axis(2, at=seq(0, 1, 0.1), labels=c("0.0", seq(0.1, 0.9, 0.1), "1.0"), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5)
lines(x = log(sapply(1:8, function(i) balanced_C64_results$reg[[i]]$time), 2),
      y = sapply(1:8, function(i) balanced_C64_results$reg[[i]]$IAD2),
      lty = 3, lwd = 3, type = 'b', pch = 5)
lines(x = log(sapply(1:8, function(i) CMC_results$mat[[i]]$time), 2),
      y = sapply(1:8, function(i) CMC_results$mat[[i]]$IAD2),
      lty = 4, lwd = 3, type = 'b', pch = 3, col = 'red')
lines(x = log(sapply(1:8, function(i) KDEMC_results$false[[i]]$time), 2),
      y = sapply(1:8, function(i) KDEMC_results$false[[i]]$IAD2),
      lty = 5, lwd = 3, type = 'b', pch = 2, col = 'red')
lines(x = log(sapply(1:8, function(i) WRS_results$rejection[[i]]$time), 2),
      y = sapply(1:8, function(i) WRS_results$rejection[[i]]$IAD2),
      lty = 6, lwd = 3, type = 'b', pch = 1, col = 'red')
legend(x = -2, y = 0.7,
       legend = c('D&C-GBF (regular mesh)',
                  'D&C-GBF (adaptive mesh)',
                  'CMC',
                  'KDEMC',
                  'WRS'),
       lwd = rep(3, 6),
       lty = c(3,2,4,5,6),
       pch = c(5,4,3,2,1),
       col = c(rep('black', 2), rep('red', 3)),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

