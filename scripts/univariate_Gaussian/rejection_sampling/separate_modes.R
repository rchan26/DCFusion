library(hierarchicalFusion)

seed <- 1994
set.seed(seed)
mu_choices <- seq(0, 3.5, 0.25)
sd <- 1
nsamples <- 10000
time_choice <- 1
target_param <- c('mu' = 0, 'sd' = sqrt(1/2))
opt_bw <- ((4/(3*nsamples))^(1/5))*sqrt(1/2)
fusion_samples_precondition <- list()
fusion_samples_standard <- list()

for (i in 1:length(mu_choices)) {
  print(paste('i:', i))
  print(paste('mu:', mu_choices[i]))
  print(paste('difference in means:', 2*mu_choices[i]))
  set.seed(seed)
  input_samples <- list(rnorm(n = nsamples, mean = -mu_choices[i], sd = sd),
                        rnorm(n = nsamples, mean = mu_choices[i], sd = sd))
  # fusion without preconditioning
  print('fusion without preconditioning')
  fusion_samples_standard[[i]] <- parallel_fusion_uniGaussian(N = nsamples,
                                                              m = 2,
                                                              time = time_choice,
                                                              samples_to_fuse = input_samples, 
                                                              means = c(-mu_choices[i], mu_choices[i]),
                                                              sds = c(sd, sd), 
                                                              betas = rep(1, 2),
                                                              precondition_values = rep(1, 2),
                                                              seed = seed)
  print(integrated_abs_distance_uniGaussian(fusion_post = fusion_samples_standard[[i]]$samples,
                                            mean = target_param['mu'],
                                            sd = target_param['sd'],
                                            beta = 1,
                                            bw = opt_bw))
  # fusion with preconditioning
  print('fusion with preconditioning')
  fusion_samples_precondition[[i]] <- parallel_fusion_uniGaussian(N = nsamples,
                                                                  m = 2,
                                                                  time = time_choice,
                                                                  samples_to_fuse = input_samples, 
                                                                  means = c(-mu_choices[i], mu_choices[i]),
                                                                  sds = c(sd, sd), 
                                                                  betas = rep(1, 2),
                                                                  precondition_values = sapply(input_samples, var),
                                                                  seed = seed)
  print(integrated_abs_distance_uniGaussian(fusion_post = fusion_samples_precondition[[i]]$samples,
                                            mean = target_param['mu'],
                                            sd = target_param['sd'],
                                            beta = 1,
                                            bw = opt_bw))
  print('saving progress')
  save.image('separate_modes.RData')
  # curve(dnorm(x, -mu_choices[i], 1), -(mu_choices[i]+5), mu_choices[i]+5,
  #       ylim = c(0,1), ylab = 'pdf', main = paste('Difference in means =', 2*mu_choices[i], '|| rep:', rep))
  # curve(dnorm(x, mu_choices[i], 1), add = T)
  # curve(dnorm(x, 0, sqrt(1/2)), add = T, lty = 2, lwd = 3)
  # lines(density(fusion_samples_standard[[i]]$samples), col = 'green')
  # lines(density(fusion_samples_precondition[[i]]$samples), col = 'blue')
}

######################################## running time

plot(x = 2*mu_choices, y = sapply(1:length(mu_choices), function(i) fusion_samples_precondition[[i]]$time),
     ylim = c(0, 2000), ylab = 'Running time in seconds', xlab = expression(paste('Value of ', mu)),
     col = 'black', lwd = 2, xaxt = "n")
lines(x = 2*mu_choices, y = sapply(1:length(mu_choices), function(i) fusion_samples_precondition[[i]]$time),
      col = 'black', lwd = 2)
axis(1, at=0:8, labels=0:8)

#################### log

plot(x = 2*mu_choices, y = sapply(1:length(mu_choices), function(i) log(fusion_samples_standard[[i]]$time)),
     ylim = c(-1,10), ylab = '', xlab = '',
     col = 'black', lwd = 3, xaxt = "n")
mtext('Difference in sub-posterior means', 1, 2.75, font = 2, cex = 1.5)
mtext('log(Time elapsed in seconds)', 2, 2.75, font = 2, cex = 1.5)
lines(x = 2*mu_choices, y = sapply(1:length(mu_choices), function(i) log(fusion_samples_standard[[i]]$time)),
      col = 'black', lwd = 3)
axis(1, at=0:8, labels=0:8, font = 2, cex = 1.5)
axis(1, at=c(0:7)+0.5, labels=rep("", 8), lwd.ticks = 0.5)
axis(2, at=seq(-2, 12, 2), labels=seq(-2, 12, 2), font = 2, cex = 1.5)
axis(2, at=-1:12, labels=rep("", 14), lwd.ticks = 0.5)

######################################## rho acceptance

plot(x = 2*mu_choices, y = sapply(1:length(mu_choices), function(i) fusion_samples_precondition[[i]]$rho), ylim = c(0, 1),
     ylab = expression(paste('Acceptance Rate for ', rho)), xlab = expression(paste('Value of ', mu)),
     col = 'black', lwd = 2, xaxt = "n")
lines(x = 2*mu_choices, y = sapply(1:length(mu_choices), function(i) fusion_samples_precondition[[i]]$rho),
      col = 'black', lwd = 2)
axis(1, at=0:8, labels=0:8)
axis(1, at=c(0:7)+0.5, labels=rep("", 8), lwd.ticks = 0.5)

######################################## Q acceptance

plot(x = 2*mu_choices, y = sapply(1:length(mu_choices), function(i) fusion_samples_precondition[[i]]$Q), ylim = c(0, 1),
     ylab = expression(paste('Acceptance Rate for ', hat(Q))), xlab = expression(paste('Value of ', mu)),
     col = 'black', lwd = 2, xaxt = "n")
lines(x = 2*mu_choices, y = sapply(1:length(mu_choices), function(i) fusion_samples_precondition[[i]]$Q),
      col = 'black', lwd = 2)
axis(1, at=0:8, labels=0:8)
axis(1, at=c(0:7)+0.5, labels=rep("", 8), lwd.ticks = 0.5)

######################################## Overall acceptance

plot(x = 2*mu_choices, y = sapply(1:length(mu_choices), function(i) fusion_samples_precondition[[i]]$rhoQ), ylim = c(0, 1),
     ylab = 'Overall Acceptance Rate', xlab = expression(paste('Value of ', mu)),
     col = 'black', lwd = 2, xaxt = "n")
lines(x = 2*mu_choices, y = sapply(1:length(mu_choices), function(i) fusion_samples_precondition[[i]]$rhoQ),
      col = 'black', lwd = 2)
axis(1, at=0:8, labels=0:8)
axis(1, at=c(0:7)+0.5, labels=rep("", 8), lwd.ticks = 0.5)

######################################## Overall acceptance (log)

plot(x = 2*mu_choices, y = log(sapply(1:length(mu_choices), function(i) fusion_samples_precondition[[i]]$rhoQ)),
     ylab = 'log(Overall Acceptance Rate)', xlab = expression(paste('Value of ', mu)),
     col = 'black', lwd = 2, xaxt = "n")
lines(x = 2*mu_choices, y = log(sapply(1:length(mu_choices), function(i) fusion_samples_precondition[[i]]$rhoQ)),
      col = 'black', lwd = 2)
axis(1, at=0:8, labels=0:8)
axis(1, at=c(0:7)+0.5, labels=rep("", 8), lwd.ticks = 0.5)
