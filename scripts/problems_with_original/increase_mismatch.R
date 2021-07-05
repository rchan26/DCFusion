library(hierarchicalFusion)

seed <- 1994
different_modes <- seq(0, 7, 0.5)
fusion_samples_standard <- list()

d2product_norm <- function(x, mu1, mu2, sd1, sd2) {
  mu <- (sd1*sd1*mu2 + sd2*sd2*mu1) / (sd1*sd1 + sd2*sd2)
  std_dv <- sqrt(1 / ((1/(sd1*sd1)) + (1/(sd2*sd2))))
  return(dnorm(x, mean = mu, sd = std_dv))
}

for (i in 1:length(different_modes)) {
  print(i)
  set.seed(seed)
  input_samples <- list(rnorm(n = 100000, mean = 0, sd = 1),
                        rnorm(n = 100000, mean = different_modes[i], sd = 1))
  
  fusion_samples_standard[[i]] <- parallel_fusion_uniGaussian(N = 10000,
                                                              m = 2,
                                                              time = 1,
                                                              samples_to_fuse = input_samples, 
                                                              means = c(0, different_modes[i]),
                                                              sds = c(1, 1), 
                                                              betas = rep(1, 2), 
                                                              precondition_values = rep(1, 2),
                                                              seed = seed)
  
  curve(dnorm(x, 0, 1), -5, different_modes[i]+5, ylim = c(0,1), ylab = 'pdf')
  curve(dnorm(x, different_modes[i], 1), -5, different_modes[i]+5, add = T)
  curve(d2product_norm(x, mu1 = 0, mu2 = different_modes[i], sd1 = 1, sd2 = 1), add = T, lty = 2)
  lines(density(fusion_samples_standard[[i]]$samples), col = 'green')
}

for (i in 1:15) {
  curve(dnorm(x, 0, 1), -5, different_modes[i]+5, ylim = c(0,1), ylab = 'pdf')
  curve(dnorm(x, different_modes[i], 1), -5, different_modes[i]+5, add = T)
  curve(d2product_norm(x, mu1 = 0, mu2 = different_modes[i], sd1 = 1, sd2 = 1), add = T, lty = 2)
  lines(density(fusion_samples_standard[[i]]$samples), col = 'green')
}

par(mai = c(1.02, 1, 0.82, 0.42))

######################################## running time

plot(x = different_modes, y = sapply(1:length(different_modes), function(i) fusion_samples_standard[[i]]$time), 
     ylim = c(0, 2000), ylab = 'Running time in seconds', xlab = expression(paste('Value of ', mu)),
     col = 'black', lwd = 2, xaxt = "n")
lines(x = different_modes, y = sapply(1:length(different_modes), function(i) fusion_samples_standard[[i]]$time),
      col = 'black', lwd = 2)
axis(1, at=0:8, labels=0:8)

#################### log

plot(x = different_modes, y = sapply(1:length(different_modes), function(i) log(fusion_samples_standard[[i]]$time)), 
     ylim = c(-1,12), ylab = 'Time Elapsed in log(seconds)', xlab = expression(paste('Value of ', mu)), 
     col = 'black', lwd = 2, xaxt = "n")
lines(x = different_modes, y = sapply(1:length(different_modes), function(i) log(fusion_samples_standard[[i]]$time)),
      col = 'black', lwd = 2)
axis(1, at=0:8, labels=0:8)
axis(1, at=c(0:7)+0.5, labels=rep("", 8), lwd.ticks = 0.5)

######################################## rho acceptance

plot(x = different_modes, y = sapply(1:length(different_modes), function(i) fusion_samples_standard[[i]]$rho), ylim = c(0, 1),
     ylab = expression(paste('Acceptance Rate for ', rho)), xlab = expression(paste('Value of ', mu)), 
     col = 'black', lwd = 2, xaxt = "n")
lines(x = different_modes, y = sapply(1:length(different_modes), function(i) fusion_samples_standard[[i]]$rho),
      col = 'black', lwd = 2)
axis(1, at=0:8, labels=0:8)
axis(1, at=c(0:7)+0.5, labels=rep("", 8), lwd.ticks = 0.5)

######################################## Q acceptance

plot(x = different_modes, y = sapply(1:length(different_modes), function(i) fusion_samples_standard[[i]]$Q), ylim = c(0, 1),
     ylab = expression(paste('Acceptance Rate for ', hat(Q))), xlab = expression(paste('Value of ', mu)), 
     col = 'black', lwd = 2, xaxt = "n")
lines(x = different_modes, y = sapply(1:length(different_modes), function(i) fusion_samples_standard[[i]]$Q),
      col = 'black', lwd = 2)
axis(1, at=0:8, labels=0:8)
axis(1, at=c(0:7)+0.5, labels=rep("", 8), lwd.ticks = 0.5)

######################################## Overall acceptance

plot(x = different_modes, y = sapply(1:length(different_modes), function(i) fusion_samples_standard[[i]]$rhoQ), ylim = c(0, 1),
     ylab = 'Overall Acceptance Rate', xlab = expression(paste('Value of ', mu)), 
     col = 'black', lwd = 2, xaxt = "n")
lines(x = different_modes, y = sapply(1:length(different_modes), function(i) fusion_samples_standard[[i]]$rhoQ), 
      col = 'black', lwd = 2)
axis(1, at=0:8, labels=0:8)
axis(1, at=c(0:7)+0.5, labels=rep("", 8), lwd.ticks = 0.5)

