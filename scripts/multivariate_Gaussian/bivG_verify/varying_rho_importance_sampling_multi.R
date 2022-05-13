library(DCFusion)

seed <- 1983
set.seed(seed)

# setting parameters
mean <- rep(0, 2)
sd <- rep(1, 2)
correlations <- c(seq(0, 0.9, 0.1), 0.95)
fusion_time <- 1
diffusion_estimator <- 'Poisson'
true_samples <- list()
input_samples <- list()
smc_fusion_standard <- list()
smc_fusion_precondition <- list()

for (i in 1:length(correlations)) {
  print(paste('i:', i))
  cov_mat <- matrix(c(1, correlations[i], correlations[i], 1),
                    nrow = 2, ncol = 2, byrow = T)
  # sample from true target density
  true_samples[[i]] <- mvrnormArma_tempered(N = 10000, mu = mean, Sigma = cov_mat, beta = 2)
  # sampling from the sub-posteriors
  input_samples[[i]] <- lapply(1:2, function(sub) mvrnormArma(N = 10000, mu = mean, Sigma = cov_mat))
  # initialising particle sets 
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples[[i]], multivariate = TRUE)
  
  # standard
  print('### performing standard fusion')
  smc_fusion_standard[[i]] <- parallel_fusion_SMC_multiGaussian(particles_to_fuse = input_particles,
                                                                N = 10000,
                                                                m = 2,
                                                                time = fusion_time,
                                                                dim = 2,
                                                                mu = mean,
                                                                Sigma = cov_mat,
                                                                betas = rep(1, 2),
                                                                precondition_matrices = rep(list(diag(1,2)), 2),
                                                                ESS_threshold = 0,
                                                                diffusion_estimator = diffusion_estimator,
                                                                seed = seed)
  print('ESS:'); print(smc_fusion_standard[[i]]$ESS)
  print('CESS:'); print(smc_fusion_standard[[i]]$CESS)
  
  # precondition
  print('### performing fusion with a preconditioning matrix')
  smc_fusion_precondition[[i]] <- parallel_fusion_SMC_multiGaussian(particles = input_particles,
                                                                    N = 10000,
                                                                    m = 2,
                                                                    time = fusion_time,
                                                                    dim = 2,
                                                                    mu = mean,
                                                                    Sigma = cov_mat,
                                                                    betas = rep(1, 2),
                                                                    precondition_matrices = lapply(input_samples[[i]], cov),
                                                                    ESS_threshold = 0,
                                                                    diffusion_estimator = diffusion_estimator,
                                                                    seed = seed)
  print('ESS:'); print(smc_fusion_precondition[[i]]$ESS)
  print('CESS:'); print(smc_fusion_precondition[[i]]$CESS)
}

save.image('varying_rho_importance_sampling_multi.RData')

plot(correlations, sapply(1:length(correlations), function(i) smc_fusion_standard[[i]]$CESS[1]), ylim = c(0, 10000),
     xlab = 'correlation', ylab = 'rho CESS', col = 'blue', xaxt='n')
axis(1, at=seq(0, 1, 0.25))
lines(correlations, sapply(1:length(correlations), function(i) smc_fusion_standard[[i]]$CESS[1]), col = 'blue')
points(correlations, sapply(1:length(correlations), function(i) smc_fusion_precondition[[i]]$CESS[1]), col = 'red')
lines(correlations, sapply(1:length(correlations), function(i) smc_fusion_precondition[[i]]$CESS[1]), col = 'red')
legend('topleft', legend = c('standard', 'preconditioned'), col = c('blue', 'red'), lty = c(1,1))

plot(correlations, sapply(1:length(correlations), function(i) smc_fusion_standard[[i]]$CESS[2]), ylim = c(0, 10000),
     xlab = 'correlation', ylab = 'Q CESS', col = 'blue')
lines(correlations, sapply(1:length(correlations), function(i) smc_fusion_standard[[i]]$CESS[2]), col = 'blue')
points(correlations, sapply(1:length(correlations), function(i) smc_fusion_precondition[[i]]$CESS[2]), col = 'red')
lines(correlations, sapply(1:length(correlations), function(i) smc_fusion_precondition[[i]]$CESS[2]), col = 'red')
legend('topleft', legend = c('standard', 'preconditioned'), col = c('blue', 'red'), lty = c(1,1))

plot(correlations, sapply(1:length(correlations), function(i) smc_fusion_standard[[i]]$ESS[2]),
     xlab = 'correlation', ylab = 'ESS', ylim = c(0, 10000),
     col = 'black', lwd = 3)
lines(correlations, sapply(1:length(correlations), function(i) smc_fusion_standard[[i]]$ESS[2]), 
      col = 'black', lwd = 3)
points(correlations, sapply(1:length(correlations), function(i) smc_fusion_precondition[[i]]$ESS[2]),
       col = 'black', lwd = 3)
lines(correlations, sapply(1:length(correlations), function(i) smc_fusion_precondition[[i]]$ESS[2]),
      col = 'black', lty = 3, lwd = 3)
legend(x = 0, y = 10000, legend = c('standard', 'preconditioned'), col = c('black', 'black'),
       lty = c(1,3), lwd = c(3,3), bty = 'n')

##### ESS per second

ESS_standard <- sapply(1:length(correlations), function(i) smc_fusion_standard[[i]]$ESS[2])
time_taken_standard <- sapply(1:length(correlations), function(i) smc_fusion_standard[[i]]$time[[1]])
ESS_per_sec_standard <- ESS_standard / time_taken_standard

ESS_precondition <- sapply(1:length(correlations), function(i) smc_fusion_precondition[[i]]$ESS[2])
time_taken_precondition <- sapply(1:length(correlations), function(i) smc_fusion_precondition[[i]]$time[[1]])
ESS_per_sec_precondition <- ESS_precondition / time_taken_precondition

plot(correlations, ESS_per_sec_standard, xlab = 'correlation', ylab = 'ESS/sec', ylim = c(0, 1000),
     col = 'black', lwd = 3)
lines(correlations, ESS_per_sec_standard, col = 'black', lwd = 3)
points(correlations, ESS_per_sec_precondition, col = 'black', lwd = 3)
lines(correlations, ESS_per_sec_precondition, col = 'black', lty = 3, lwd = 3)
legend(x = 0, y = 1000, legend = c('standard', 'preconditioned'), col = c('black', 'black'),
       lty = c(1,3), lwd = c(3,3), bty = 'n')
