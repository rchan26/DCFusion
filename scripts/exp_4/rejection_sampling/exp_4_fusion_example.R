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

########## USING R CODE [standard]
test_exp_4_fusion_standard <- lapply(time_choices_exp_4, function(time) {
  print(paste('time:', time))
  return(parallel_fusion_exp_4(N = 1000000,
                               m = 4,
                               time = time,
                               samples_to_fuse = sub_posterior_samples,
                               mean = 0,
                               betas = rep(1/4, 4),
                               precondition_values = rep(1, 4),
                               seed = seed))})

curve(exp_4_density(x), -4, 4)
for (i in 1:length(time_choices_exp_4)) {
  print(paste('time:', time_choices_exp_4[i]))
  print(ks.test(test_exp_4_fusion_standard[[i]]$samples, true_samples))
  lines(density(test_exp_4_fusion_standard[[i]]$samples), col = 'red')
}

########## USING R CODE [precondition]
test_exp_4_fusion_precondition <- lapply(time_choices_exp_4, function(time) {
  print(paste('time:', time))
  return(parallel_fusion_exp_4(N = 1000000,
                               m = 4,
                               time = time,
                               samples_to_fuse = sub_posterior_samples,
                               mean = 0,
                               betas = rep(1/4, 4),
                               precondition_values = sapply(sub_posterior_samples, var),
                               seed = seed))})

curve(exp_4_density(x), -4, 4)
for (i in 1:length(time_choices_exp_4)) {
  print(paste('time:', time_choices_exp_4[i]))
  print(ks.test(test_exp_4_fusion_precondition[[i]]$samples, true_samples))
  lines(density(test_exp_4_fusion_precondition[[i]]$samples), col = 'red')
}
