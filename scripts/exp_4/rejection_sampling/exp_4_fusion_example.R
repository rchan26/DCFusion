library(hierarchicalFusion)
time_choices_exp_4 <- c(0.1, 0.2, 0.5, 1, 1.2)
seed <- 21
set.seed(seed)

# using rejection sampling to obtain input samples
set.seed(seed)
sub_posterior_samples <- base_rejection_sampler_exp_4(beta = 1/4,
                                                      nsamples = 1000000,
                                                      proposal_mean = 0,
                                                      proposal_sd = 1.5,
                                                      dominating_M = 1.4)
curve(exp_4_density(x, beta = 1/4), -4, 4, ylim = c(-0.5,1))
for (i in 1:4) {
  lines(density(sub_posterior_samples[[i]]), col = 'blue')
}

# direct sampling from truth
set.seed(seed)
true_samples <- unlist(base_rejection_sampler_exp_4(beta = 1,
                                                    nsamples = 1000000,
                                                    proposal_mean = 0,
                                                    proposal_sd = 1,
                                                    dominating_M = 1.4))
curve(exp_4_density(x), -4, 4)
curve(1.4*dnorm(x, 0, 1), add = T)
lines(density(true_samples), col = 'blue')

########## USING R CODE [standard, Poisson]
exp_4_fusion_standard_Poisson <- lapply(time_choices_exp_4, function(time) {
  print(paste('time:', time))
  return(parallel_fusion_exp_4(N = 1000000,
                               m = 4,
                               time = time,
                               samples_to_fuse = sub_posterior_samples,
                               mean = 0,
                               betas = rep(1/4, 4),
                               precondition_values = rep(1, 4),
                               diffusion_estimator = 'Poisson',
                               seed = seed))})

plot(x = time_choices_exp_4, 
     y = sapply(1:length(time_choices_exp_4), function(i) exp_4_fusion_standard_Poisson[[i]]$rho),
     xlab = 'Time',
     ylab = 'Prob',
     pch = 4,
     ylim = c(0, 1))
lines(x = time_choices_exp_4, 
      y = sapply(1:length(time_choices_exp_4), function(i) exp_4_fusion_standard_Poisson[[i]]$rho))
points(x = time_choices_exp_4, 
       y = sapply(1:length(time_choices_exp_4), function(i) exp_4_fusion_standard_Poisson[[i]]$Q))
lines(x = time_choices_exp_4, 
      y = sapply(1:length(time_choices_exp_4), function(i) exp_4_fusion_standard_Poisson[[i]]$Q),
      lty = 2)

curve(exp_4_density(x), -4, 4)
for (i in 1:length(time_choices_exp_4)) {
  print(paste('time:', time_choices_exp_4[i]))
  print(ks.test(exp_4_fusion_standard_Poisson[[i]]$samples, true_samples))
  lines(density(exp_4_fusion_standard_Poisson[[i]]$samples), col = 'red')
}

########## USING R CODE [preconditionm Poisson]
exp_4_fusion_precondition_Poisson <- lapply(time_choices_exp_4, function(time) {
  print(paste('time:', time))
  return(parallel_fusion_exp_4(N = 1000000,
                               m = 4,
                               time = time,
                               samples_to_fuse = sub_posterior_samples,
                               mean = 0,
                               betas = rep(1/4, 4),
                               precondition_values = sapply(sub_posterior_samples, var),
                               diffusion_estimator = 'Poisson',
                               seed = seed))})

plot(x = time_choices_exp_4, 
     y = sapply(1:length(time_choices_exp_4), function(i) exp_4_fusion_precondition_Poisson[[i]]$rho),
     xlab = 'Time',
     ylab = 'Prob',
     pch = 4,
     ylim = c(0, 1))
lines(x = time_choices_exp_4, 
      y = sapply(1:length(time_choices_exp_4), function(i) exp_4_fusion_precondition_Poisson[[i]]$rho))
points(x = time_choices_exp_4, 
       y = sapply(1:length(time_choices_exp_4), function(i) exp_4_fusion_precondition_Poisson[[i]]$Q))
lines(x = time_choices_exp_4, 
      y = sapply(1:length(time_choices_exp_4), function(i) exp_4_fusion_precondition_Poisson[[i]]$Q),
      lty = 2)

curve(exp_4_density(x), -4, 4)
for (i in 1:length(time_choices_exp_4)) {
  print(paste('time:', time_choices_exp_4[i]))
  print(ks.test(exp_4_fusion_precondition_Poisson[[i]]$samples, true_samples))
  lines(density(exp_4_fusion_precondition_Poisson[[i]]$samples), col = 'red')
}

########## USING R CODE [standard, NB]
exp_4_fusion_standard_NB <- lapply(time_choices_exp_4, function(time) {
  print(paste('time:', time))
  return(parallel_fusion_exp_4(N = 1000000,
                               m = 4,
                               time = time,
                               samples_to_fuse = sub_posterior_samples,
                               mean = 0,
                               betas = rep(1/4, 4),
                               precondition_values = rep(1, 4),
                               diffusion_estimator = 'NB',
                               seed = seed))})

plot(x = time_choices_exp_4, 
     y = sapply(1:length(time_choices_exp_4), function(i) exp_4_fusion_standard_NB[[i]]$rho),
     xlab = 'Time',
     ylab = 'Prob',
     pch = 4,
     ylim = c(0, 1))
lines(x = time_choices_exp_4, 
      y = sapply(1:length(time_choices_exp_4), function(i) exp_4_fusion_standard_NB[[i]]$rho))
points(x = time_choices_exp_4, 
       y = sapply(1:length(time_choices_exp_4), function(i) exp_4_fusion_standard_NB[[i]]$Q))
lines(x = time_choices_exp_4, 
      y = sapply(1:length(time_choices_exp_4), function(i) exp_4_fusion_standard_NB[[i]]$Q),
      lty = 2)

curve(exp_4_density(x), -4, 4)
for (i in 1:length(time_choices_exp_4)) {
  print(paste('time:', time_choices_exp_4[i]))
  print(ks.test(exp_4_fusion_standard_NB[[i]]$samples, true_samples))
  lines(density(exp_4_fusion_standard_NB[[i]]$samples), col = 'red')
}

sapply(1:length(time_choices_exp_4), function(i) exp_4_fusion_standard_Poisson[[i]]$Q)
sapply(1:length(time_choices_exp_4), function(i) exp_4_fusion_standard_NB[[i]]$Q)

########## USING R CODE [precondition, NB]
exp_4_fusion_precondition_NB <- lapply(time_choices_exp_4, function(time) {
  print(paste('time:', time))
  return(parallel_fusion_exp_4(N = 1000000,
                               m = 4,
                               time = time,
                               samples_to_fuse = sub_posterior_samples,
                               mean = 0,
                               betas = rep(1/4, 4),
                               precondition_values = sapply(sub_posterior_samples, var),
                               diffusion_estimator = 'NB',
                               seed = seed))})

plot(x = time_choices_exp_4, 
     y = sapply(1:length(time_choices_exp_4), function(i) exp_4_fusion_precondition_NB[[i]]$rho),
     xlab = 'Time',
     ylab = 'Prob',
     pch = 4,
     ylim = c(0, 1))
lines(x = time_choices_exp_4, 
      y = sapply(1:length(time_choices_exp_4), function(i) exp_4_fusion_precondition_NB[[i]]$rho))
points(x = time_choices_exp_4, 
       y = sapply(1:length(time_choices_exp_4), function(i) exp_4_fusion_precondition_NB[[i]]$Q))
lines(x = time_choices_exp_4, 
      y = sapply(1:length(time_choices_exp_4), function(i) exp_4_fusion_precondition_NB[[i]]$Q),
      lty = 2)

curve(exp_4_density(x), -4, 4)
for (i in 1:length(time_choices_exp_4)) {
  print(paste('time:', time_choices_exp_4[i]))
  print(ks.test(exp_4_fusion_precondition_NB[[i]]$samples, true_samples))
  lines(density(exp_4_fusion_precondition_NB[[i]]$samples), col = 'red')
}

sapply(1:length(time_choices_exp_4), function(i) exp_4_fusion_precondition_Poisson[[i]]$Q)
sapply(1:length(time_choices_exp_4), function(i) exp_4_fusion_precondition_NB[[i]]$Q)

##### Comparison #####

consensus <- consensus_scott(S = 4,
                             samples_to_combine = lapply(sub_posterior_samples, as.matrix),
                             indep = FALSE)
lines(density(consensus$samples))
neisw <- neiswanger(S = 4,
                    samples_to_combine = lapply(sub_posterior_samples, as.matrix),
                    anneal = FALSE)
lines(density(neisw$samples))
weier <- weierstrass(Samples = sub_posterior_samples,
                     method = 'importance', para.dim = 1)
lines(density(weier$samples))


curve(exp_4_density(x), -4, 4, ylim = c(0, 1), ylab = 'pdf', lwd = 3)
lines(density(exp_4_fusion_precondition[[4]]$samples), col = '#D41159', lty = 5, lwd = 3)
lines(density(consensus$samples), col = '#FFC20A', lty = 4, lwd = 3)
lines(density(neisw$samples), col = '#0C7BDC', lty = 3, lwd = 3) 
lines(density(weier$samples), col = '#994F00', lty = 2, lwd = 3)
legend(x = 2, y = 1,
       legend = c('Target', 'Fusion', 'Consensus', 'Neiswanger', 'Weierstrass'),
       col = c('black', '#D41159', '#FFC20A', '#0C7BDC', '#994F00'),
       lwd = c(3, 3, 3, 3, 3),
       lty = c(1, 5, 4, 3, 2))

