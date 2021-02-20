library(hierarchicalFusion)
library(MASS)

seed <- 1983
set.seed(seed)

# setting parameters
mean <- rep(0, 2)
sd <- rep(1, 2)
correlations <- seq(0, 0.8, 0.1)
fusion_time <- 1

true_samples <- list()
true_kde <- list()
input_samples <- list()
fusion_standard <- list()
fusion_precondition <- list()
kde_standard <- list()
kde_precondition <- list()

for (i in 1:length(correlations)) {
  print(paste('i:', i))
  cov_mat <- matrix(c(1, correlations[i], correlations[i], 1),
                    nrow = 2, ncol = 2, byrow = T)
  # sample from true target density
  true_samples[[i]] <- mvrnormArma_tempered(N = 5000, mu = mean, Sigma = cov_mat, beta = 2)
  true_kde[[i]] <- MASS::kde2d(true_samples[[i]][,1], true_samples[[i]][,2], n = 50)
  xlims <- c(-3, 3)
  ylims <- c(-3, 3)
  image(true_kde[[i]], xlim = xlims, ylim = ylims)
  title(paste("correlation:", correlations[i]))
  contour(true_kde[[i]], add = T)
  
  # sampling from the sub-posteriors
  input_samples[[i]] <- lapply(1:2, function(sub) mvrnormArma(N = 5000, mu = mean, Sigma = cov_mat))
  for (sub in 1:2) {
    contour(MASS::kde2d(input_samples[[i]][[sub]][,1],
                        input_samples[[i]][[sub]][,2], n = 50),
            col = 'green', add = T)
  }
  
  # standard
  print('### performing standard fusion')
  fusion_standard[[i]] <- parallel_fusion_biGaussian(N = 10000,
                                                     m = 2,
                                                     time = fusion_time,
                                                     samples_to_fuse = input_samples[[i]],
                                                     mean_vec = mean,
                                                     sd_vec = sd,
                                                     corr = correlations[i],
                                                     betas = rep(1, 2),
                                                     precondition_matrices = rep(list(diag(1,2)), 2),
                                                     seed = seed)
  kde_standard[[i]] <- MASS::kde2d(fusion_standard[[i]]$samples[,1],
                                   fusion_standard[[i]]$samples[,2], 
                                   n = 50)
  print(paste('rho:', fusion_standard[[i]]$rho))
  print(paste('Q:', fusion_standard[[i]]$Q))
  print(paste('time:', fusion_standard[[i]]$time))
  
  # precondition
  print('### performing fusion with a preconditioning matrix')
  fusion_precondition[[i]] <- parallel_fusion_biGaussian(N = 10000,
                                                         m = 2,
                                                         time = fusion_time,
                                                         samples_to_fuse = input_samples[[i]],
                                                         mean_vec = mean,
                                                         sd_vec = sd,
                                                         corr = correlations[i],
                                                         betas = rep(1, 2),
                                                         precondition_matrices = lapply(input_samples[[i]], cov),
                                                         seed = seed)
  kde_precondition[[i]] <- MASS::kde2d(fusion_precondition[[i]]$samples[,1],
                                       fusion_precondition[[i]]$samples[,2], 
                                       n = 50)
  print(paste('rho:', fusion_precondition[[i]]$rho))
  print(paste('Q:', fusion_precondition[[i]]$Q))
  print(paste('time:', fusion_precondition[[i]]$time))
  
  # adding to plot
  contour(kde_standard[[i]], add = T, col = 'blue')
  contour(kde_precondition[[i]], add = T, col = 'red')
  legend('topleft', legend = c('true', 'standard', 'precond'), col = c('black', 'blue', 'red'), lty = c(1,1,1))
  # looking at 1 dimensional kdes
  curve(dnorm(x, 0, 1/sqrt(2)), -4, 4, ylab = 'pdf', main = paste("correlation:", correlations[i]))
  lines(density(fusion_standard[[i]]$samples[,1]), col = 'blue')
  lines(density(fusion_standard[[i]]$samples[,2]), col = 'blue')
  lines(density(fusion_precondition[[i]]$samples[,1]), col = 'red')
  lines(density(fusion_precondition[[i]]$samples[,2]), col = 'red')
}

par(mai = c(1.02, 1, 0.82, 0.42))

plot(correlations, sapply(1:9, function(i) fusion_standard[[i]]$rho), ylim = c(0, 1),
     xlab = 'correlation', ylab = expression(paste('Acceptance Rate for ', rho)), col = 'black', lwd = 3)
lines(correlations, sapply(1:9, function(i) fusion_standard[[i]]$rho), col = 'black', lwd = 3)
points(correlations, sapply(1:9, function(i) fusion_precondition[[i]]$rho), col = 'black', lty = 3, lwd = 3)
lines(correlations, sapply(1:9, function(i) fusion_precondition[[i]]$rho), col = 'black', lwd = 3)
legend(x = 0, y = 1, legend = c('standard', 'preconditioned'), col = c('black', 'black'), lty = c(1,3), lwd = c(3,3), bty = 'n')

plot(correlations, sapply(1:9, function(i) fusion_standard[[i]]$Q), ylim = c(0, 0.4),
     xlab = 'correlation', ylab = expression(paste('Acceptance Rate for ', hat(Q))), col = 'black', lwd = 3)
lines(correlations, sapply(1:9, function(i) fusion_standard[[i]]$Q), col = 'black', lwd = 3)
points(correlations, sapply(1:9, function(i) fusion_precondition[[i]]$Q), col = 'black', lty = 3, lwd = 3)
lines(correlations, sapply(1:9, function(i) fusion_precondition[[i]]$Q), col = 'black', lty = 3, lwd = 3)
legend(x = 0, y = 1, legend = c('standard', 'preconditioned'), col = c('black', 'black'), lty = c(1,3), lwd = c(3,3), bty = 'n')

plot(correlations, log(sapply(1:9, function(i) fusion_standard[[i]]$time)),
     xlab = 'correlation', ylab = 'Time Elapsed in log(seconds)', col = 'black', ylim = c(0, 8), lwd = 3)
lines(correlations, log(sapply(1:9, function(i) fusion_standard[[i]]$time)), col = 'black', lwd = 3)
points(correlations, log(sapply(1:9, function(i) fusion_precondition[[i]]$time)), col = 'black', lty = 3, lwd = 3)
lines(correlations, log(sapply(1:9, function(i) fusion_precondition[[i]]$time)), col = 'black', lty = 3, lwd = 3)
legend(x = 0, y = 8, legend = c('standard', 'preconditioned'), col = c('black', 'black'), lty = c(1,3), lwd = c(3,3), bty = 'n')

