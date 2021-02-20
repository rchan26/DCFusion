library(hierarchicalFusion)
time_choices_exp_4 <- c(0.1, 0.2, 0.5, 1, 1.2)
seed <- 21
set.seed(seed)

# using rejection sampling to obtain input samples
input_samples_exp_4 <- unlist(base_rejection_sampler_exp_4(beta = 1/4,
                                                           nsamples = 250000,
                                                           proposal_mean = 0,
                                                           proposal_sd = 1.5,
                                                           dominating_M = 1.4))
curve(exp_4_density(x, beta = 1/4), -4, 4)
lines(density(input_samples_exp_4), col = 'blue')

########## USING R CODE [standard]
test_exp_4_fusion_1_core_standard <- lapply(time_choices_exp_4, function(time) {
  print(paste('time:', time))
  return(parallel_fusion_exp_4(N = 1000000,
                               m = 1,
                               time = time,
                               samples_to_fuse = list(input_samples_exp_4),
                               mean = 0,
                               betas = 1/4,
                               precondition_values = 1,
                               seed = seed))})

curve(exp_4_density(x, beta = 1/4), -4, 4)
for (i in 1:length(time_choices_exp_4)) {
  print(paste('time:', time_choices_exp_4[i]))
  print(ks.test(test_exp_4_fusion_1_core_standard[[i]]$samples, input_samples_exp_4))
  lines(density(test_exp_4_fusion_1_core_standard[[i]]$samples), col = 'red')
}

########## USING R CODE [precondition]
test_exp_4_fusion_1_core_precondition <- lapply(time_choices_exp_4, function(time) {
  print(paste('time:', time))
  return(parallel_fusion_exp_4(N = 1000000,
                               m = 1,
                               time = time,
                               samples_to_fuse = list(input_samples_exp_4),
                               mean = 0,
                               betas = 1/4,
                               precondition_values = var(input_samples_exp_4),
                               seed = seed))})

curve(exp_4_density(x, beta = 1/4), -4, 4)
for (i in 1:length(time_choices_exp_4)) {
  print(paste('time:', time_choices_exp_4[i]))
  print(ks.test(test_exp_4_fusion_1_core_precondition[[i]]$samples, input_samples_exp_4))
  lines(density(test_exp_4_fusion_1_core_precondition[[i]]$samples), col = 'red')
}
