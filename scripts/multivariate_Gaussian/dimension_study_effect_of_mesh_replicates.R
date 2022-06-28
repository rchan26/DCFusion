library(DCFusion)

seed <- 1994
set.seed(seed)
nsamples <- 10000
dimension <- 50
C <- 8
corr <- 0.9
opt_bw <- ((4)/(3*nsamples))^(1/5)
diffusion_estimator <- 'NB'
resampling_method <- 'resid'
number_of_replicates <- 10
ESS_threshold <- 0.5
CESS_0_threshold <- 0.05
CESS_j_threshold <- c(0.05, 0.025, 0.01)
vanilla_b <- 1
data_size <- 1
n_cores <- parallel::detectCores()
dc_gbf <- list('regular' = list(), 'adaptive' = list())

collect_results <- function(results, dc, dim, marg_means, marg_sds, seed) {
  if (dc) {
    print(paste('time:', sum(unlist(results$time))))
    print(paste('log(time):', log(sum(unlist(results$time)))))
    return(list('CESS' = results$CESS,
                'time' = results$time,
                'time_mesh' = results$time_mesh,
                'elapsed_time' = results$elapsed_time,
                'resampled' = results$resampled,
                'ESS' = results$ESS,
                'E_nu_j' = results$E_nu_j,
                'chosen' = results$chosen,
                'mesh_terms' = results$mesh_terms,
                'k4_choice' = results$k4_choice,
                'recommended_mesh' = results$recommended_mesh,
                'IAD' = integrated_abs_distance_multiGaussian(fusion_post = resample_particle_y_samples(
                  particle_set = results$particles[[1]],
                  multivariate = TRUE,
                  resampling_method = resampling_method,
                  seed = seed)$y_samples,
                  marg_means = marg_means,
                  marg_sds = marg_sds,
                  bw = rep(opt_bw, dim))))
  } else {
    print(paste('time:', results$time))
    print(paste('log(time):', log(results$time)))
    return(list('CESS' = results$CESS,
                'time_mesh' = results$particles$time_mesh,
                'time' = results$time,
                'elapsed_time' = results$elapsed_time,
                'resampled' = results$resampled,
                'ESS' = results$ESS,
                'E_nu_j' = results$E_nu_j,
                'chosen' = results$chosen,
                'mesh_terms' = results$mesh_terms,
                'k4_choice' = results$k4_choice,
                'IAD' = integrated_abs_distance_multiGaussian(fusion_post = resample_particle_y_samples(
                  particle_set = results$particles,
                  multivariate = TRUE,
                  resampling_method = resampling_method,
                  seed = seed)$y_samples,
                  marg_means = marg_means,
                  marg_sds = marg_sds,
                  bw = rep(opt_bw, dim))))
  }
}

for (i in 1:length(CESS_j_threshold)) {
  print(paste('##### i:', i, '#####'))
  print(paste('%%%%% CESS_j_threshold:', CESS_j_threshold[i], '%%%%%'))
  dc_gbf$regular[[i]] <- list()
  dc_gbf$adaptive[[i]] <- list()
  for (rep in 1:number_of_replicates) {
    print(paste('rep:', rep))
    set.seed(seed*rep*i)
    mean <- rep(0, dimension)
    cov_mat <- matrix(data = corr, nrow = dimension, ncol = dimension)
    diag(cov_mat) <- 1
    cov_mat <- (C/data_size)*cov_mat
    input_samples <- lapply(1:C, function(sub) mvrnormArma(N = nsamples, mu = mean, Sigma = cov_mat))
    if (dimension==1) {
      sub_posterior_means <- as.matrix(sapply(input_samples, function(sub) apply(sub, 2, mean)))
    } else {
      sub_posterior_means <- t(sapply(input_samples, function(sub) apply(sub, 2, mean)))
    }
    ##### Divide-and-Conquer GBF #####
    print('### performing D&C-Generalised Bayesian Fusion (with recommended T, regular mesh)')
    dc_gbf_regular <- bal_binary_GBF_multiGaussian(N_schedule = rep(nsamples, 3),
                                                   m_schedule = rep(2, 3),
                                                   base_samples = input_samples,
                                                   L = 4,
                                                   dim = dimension,
                                                   mean_vecs = rep(list(mean), C),
                                                   Sigmas = rep(list(cov_mat), C),
                                                   C = C,
                                                   precondition = TRUE,
                                                   resampling_method = resampling_method,
                                                   ESS_threshold = ESS_threshold,
                                                   adaptive_mesh = FALSE,
                                                   mesh_parameters = list('condition' = 'SH',
                                                                          'CESS_0_threshold' = CESS_0_threshold,
                                                                          'CESS_j_threshold' = CESS_j_threshold[i],
                                                                          'vanilla' = FALSE),
                                                   diffusion_estimator = diffusion_estimator,
                                                   seed = seed*i*rep)
    dc_gbf$regular[[i]][[rep]] <- collect_results(results = dc_gbf_regular,
                                                  dc = TRUE,
                                                  dim = dimension,
                                                  marg_means = mean,
                                                  marg_sds = rep(sqrt(1/data_size), dimension),
                                                  seed = seed*i*rep)
    print('### performing D&C-Generalised Bayesian Fusion (with recommended T, adaptive mesh)')
    dc_gbf_adaptive <- bal_binary_GBF_multiGaussian(N_schedule = rep(nsamples, 3),
                                                    m_schedule = rep(2, 3),
                                                    base_samples = input_samples,
                                                    L = 4,
                                                    dim = dimension,
                                                    mean_vecs = rep(list(mean), C),
                                                    Sigmas = rep(list(cov_mat), C),
                                                    C = C,
                                                    precondition = TRUE,
                                                    resampling_method = resampling_method,
                                                    ESS_threshold = ESS_threshold,
                                                    adaptive_mesh = TRUE,
                                                    mesh_parameters = list('condition' = 'SH',
                                                                           'CESS_0_threshold' = CESS_0_threshold,
                                                                           'CESS_j_threshold' = CESS_j_threshold[i],
                                                                           'vanilla' = FALSE),
                                                    diffusion_estimator = diffusion_estimator,
                                                    seed = seed*i*rep)
    dc_gbf$adaptive[[i]][[rep]] <- collect_results(results = dc_gbf_adaptive,
                                                   dc = TRUE,
                                                   dim = dimension,
                                                   marg_means = mean,
                                                   marg_sds = rep(sqrt(1/data_size), dimension),
                                                   seed = seed*i*rep)
    
    print('saving progress')
    save.image('dimension_study_effect_of_mesh_replicates.RData')
  }
}

##### Paper: IAD #####
plot(x = CESS_j_threshold,
     y = sapply(1:length(CESS_j_threshold), function(i) {
       mean(sapply(1:number_of_replicates, function(rep) dc_gbf$regular[[i]][[rep]]$IAD))
     }),
     type = 'b', pch = 1, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = CESS_j_threshold,
      y = sapply(1:length(CESS_j_threshold), function(i) {
        mean(sapply(1:number_of_replicates, function(rep) dc_gbf$adaptive[[i]][[rep]]$IAD))
      }),
      pch = 2, lty = 2, lwd = 3, type = 'b', col = 'red')
for (i in 1:length(CESS_j_threshold)) {
  IAD <- sapply(1:number_of_replicates, function(rep) dc_gbf$regular[[i]][[rep]]$IAD)
  points(x = rep(CESS_j_threshold[i], length(IAD)), y = IAD, cex = 0.5, pch = 1, lwd = 1.5)
}
for (i in 1:length(CESS_j_threshold)) {
  IAD <- sapply(1:number_of_replicates, function(rep) dc_gbf$adaptive[[i]][[rep]]$IAD)
  points(x = rep(CESS_j_threshold[i], length(IAD)), y = IAD, cex = 0.5, pch = 2, lwd = 1.5, col = 'red')
}
axis(1, at = log(CESS_j_threshold), labels = log(CESS_j_threshold), font = 2, cex = 1.5)
mtext('CESS_j_threshold', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(-4, 0, 0.1), labels = c("0.0", seq(0.1, 0.9, 0.1), "1.0", seq(1.1, 1.6, 0.1)),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1.6, 0.1), labels=rep("", 17), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
legend(x = 10000, y = 1,
       legend = c('D&C-GBF (regular)',
                  'D&C-GBF (adaptive)'),
       col = c('black', 'red', 'blue', 'green', 'black', 'purple', 'orange', 'pink'),
       lty = c(1, 2, 3, 4, 1, 5, 6, 1),
       pch = c(1, 2, 3, 4, 7, 5, 6, 8),
       lwd = rep(3, 8),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

##### Paper: time #####
plot(x = log(CESS_j_threshold, 2),
     y = sapply(1:length(CESS_j_threshold), function(i) {
       mean(sapply(1:number_of_replicates, function(rep) sum(unlist(dc_gbf$regular[[i]][[rep]]$time))))
     }),
     type = 'b', pch = 1, lty = 1, lwd = 3, ylim = c(2,14), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = log(CESS_j_threshold, 2),
      y = sapply(1:length(CESS_j_threshold), function(i) {
        mean(log(sapply(1:number_of_replicates, function(rep) sum(unlist(dc_gbf$adaptive[[i]][[rep]]$time))), 2))
      }),
      pch = 2, lty = 2, lwd = 3, type = 'b', col = 'red')
axis(1, at = 0:6, labels = 0:6, font = 2, cex = 1.5)
axis(1, at = seq(0, 7, 0.5), labels = rep("", 15), lwd.ticks = 0.5)
mtext('log(CESS_j_threshold, 2)', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 14, 1), labels = seq(0, 14, 1), font = 2, cex = 1.5)
mtext('log(Elapsed time in seconds, 2)', 2, 2.75, font = 2, cex = 1.5)
legend(x = 0, y = 14,
       legend = c('D&C-GBF (regular)',
                  'D&C-GBF (adaptive)',
                  'GBF (regular)',
                  'GBF (adaptive)',
                  'D&C-MCF',
                  'BF (regular)',
                  'BF (adaptive)'),
       col = c('black', 'red', 'blue', 'green', 'black', 'purple', 'orange'),
       lty = c(1, 2, 3, 4, 1, 5, 6),
       pch = c(1, 2, 3, 4, 7, 5, 6),
       lwd = rep(3, 7),
       cex = 1.25,
       text.font = 2,
       bty = 'n')


