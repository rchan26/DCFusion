library(hierarchicalFusion)

seed <- 1994
set.seed(seed)
different_modes <- seq(0, 8, 0.5)
smc_fusion_samples_precondition <- list()
smc_fusion_samples_standard <- list()

d2product_norm <- function(x, mu1, mu2, sd1, sd2) {
  mu <- (sd1*sd1*mu2 + sd2*sd2*mu1) / (sd1*sd1 + sd2*sd2)
  std_dv <- sqrt(1 / ((1/(sd1*sd1)) + (1/(sd2*sd2))))
  return(dnorm(x, mean = mu, sd = std_dv))
}

for (i in 1:length(different_modes)) {
  print(i)
  samples1 <- rnorm(n = 100000, mean = 0, sd = 1)
  samples2 <- rnorm(n = 100000, mean = different_modes[i], sd = 1)
  input_samples <- list(samples1, samples2)
  input_particles <- initialise_particle_sets(input_samples, multivariate = FALSE)
  
  smc_fusion_samples_precondition[[i]] <- parallel_fusion_SMC_uniGaussian(particles_to_fuse = input_particles, 
                                                                          N = 10000,
                                                                          m = 2, 
                                                                          time = 1,
                                                                          means = c(0, different_modes[i]),
                                                                          sds = c(1, 1), 
                                                                          betas = rep(1, 2),
                                                                          precondition_values = sapply(input_samples, var),
                                                                          ESS_threshold = 0,
                                                                          resampling_method = 'resid',
                                                                          seed = seed)
  
  smc_fusion_samples_standard[[i]] <- parallel_fusion_SMC_uniGaussian(particles_to_fuse = input_particles, 
                                                                      N = 10000,
                                                                      m = 2, 
                                                                      time = 1,
                                                                      means = c(0, different_modes[i]),
                                                                      sds = c(1, 1), 
                                                                      betas = rep(1, 2),
                                                                      precondition_values = rep(1, 2),
                                                                      ESS_threshold = 0,
                                                                      resampling_method = 'resid',
                                                                      seed = seed)
  
  curve(dnorm(x, 0, 1), -5, different_modes[i]+5, ylim = c(0,1), ylab = 'pdf',
        main = paste('mu =', different_modes[i]))
  curve(dnorm(x, different_modes[i], 1), -5, different_modes[i]+5, add = T)
  curve(d2product_norm(x, mu1 = 0, mu2 = different_modes[i], sd1 = 1, sd2 = 1), add = T, lty = 2)
  lines(density(resample_particle_y_samples(particle_set = smc_fusion_samples_standard[[i]]$particles,
                                            multivariate = FALSE,
                                            resampling_method = 'resid',
                                            seed = seed)$y_samples,
                adjust = 2), col = 'green')
  lines(density(resample_particle_y_samples(particle_set = smc_fusion_samples_precondition[[i]]$particles,
                                            multivariate = FALSE,
                                            resampling_method = 'resid',
                                            seed = seed)$y_samples, 
                adjust = 2), col = 'blue')
}

par(mai = c(1.02, 1, 0.82, 0.42))

######################################## running time

plot(x = different_modes, y = sapply(1:length(different_modes), function(i) smc_fusion_samples_precondition[[i]]$time), 
     ylim = c(0, 10), ylab = 'Running time in seconds', xlab = expression(paste('Value of ', mu)),
     col = 'black', lwd = 3, xaxt = "n")
lines(x = different_modes, y = sapply(1:length(different_modes), function(i) smc_fusion_samples_precondition[[i]]$time),
      col = 'black', lwd = 3)
axis(1, at=0:8, labels=0:8)

#################### log

plot(x = different_modes, y = sapply(1:length(different_modes), function(i) log(smc_fusion_samples_precondition[[i]]$time)), 
     ylim = c(-1, 3), ylab = 'Time Elapsed in log(seconds)', xlab = expression(paste('Value of ', mu)), 
     col = 'black', lwd = 3, xaxt = "n")
lines(x = different_modes, y = sapply(1:length(different_modes), function(i) log(smc_fusion_samples_precondition[[i]]$time)),
      col = 'black', lwd = 3)
axis(1, at=0:8, labels=0:8)
axis(1, at=c(0:7)+0.5, labels=rep("", 8), lwd.ticks = 0.5)

######################################## rho CESS

plot(x = different_modes, y = sapply(1:length(different_modes), function(i) smc_fusion_samples_precondition[[i]]$CESS['rho']), ylim = c(0, 10000),
     ylab = expression(paste('CESS after ', rho)), xlab = expression(paste('Value of ', mu)), 
     col = 'black', lwd = 3, xaxt = "n")
lines(x = different_modes, y = sapply(1:length(different_modes), function(i) smc_fusion_samples_precondition[[i]]$CESS['rho']), 
      col = 'black', lwd = 3)
axis(1, at=0:8, labels=0:8)
axis(1, at=c(0:7)+0.5, labels=rep("", 8), lwd.ticks = 0.5)

######################################## Q CESS

plot(x = different_modes, y = sapply(1:length(different_modes), function(i) smc_fusion_samples_precondition[[i]]$CESS['Q']), ylim = c(0, 10000),
     ylab = expression(paste('CESS after ', tilde(Q))), xlab = expression(paste('Value of ', mu)),
     col = 'black', lwd = 3, xaxt = "n")
lines(x = different_modes, y = sapply(1:length(different_modes), function(i) smc_fusion_samples_precondition[[i]]$CESS['Q']),
      col = 'black', lwd = 3)
axis(1, at=0:8, labels=0:8)
axis(1, at=c(0:7)+0.5, labels=rep("", 8), lwd.ticks = 0.5)

######################################## rho*Q ESS (there's no resampling done after rho step)

plot(x = different_modes, y = sapply(1:length(different_modes), function(i) smc_fusion_samples_precondition[[i]]$ESS['Q']), ylim = c(0, 10000),
     ylab = 'ESS', xlab = expression(paste('Value of ', mu)), 
     col = 'black', lwd = 3, xaxt = "n")
lines(x = different_modes, y = sapply(1:length(different_modes), function(i) smc_fusion_samples_precondition[[i]]$ESS['Q']), 
      col = 'black', lwd = 3)
axis(1, at=0:8, labels=0:8)
axis(1, at=c(0:7)+0.5, labels=rep("", 8), lwd.ticks = 0.5)

