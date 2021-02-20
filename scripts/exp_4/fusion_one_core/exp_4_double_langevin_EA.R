time_choices_exp_4 <- c(0.1, 0.2, 0.5, 1, 1.2)
set.seed(21)

# using rejection sampling to obtain input samples
input_samples_exp_4 <- unlist(base_rejection_sampler_exp_4(beta = 1/4,
                                                           nsamples = 250000,
                                                           proposal_mean = 0,
                                                           proposal_sd = 1.5,
                                                           dominating_M = 1.4))
curve(exp_4_density(x, beta = 1/4), -4, 4, ylim = c(-0.5,1))
curve(ea_phi_exp_4(x, mean = 0, beta = 1/4), add = T, col = 'red', lty = 2)
lines(density(input_samples_exp_4), col = 'blue')

set.seed(21)
# applying the exact algorithm with double langevin [standard] (USING R CODE)
test_exp_4_DL_R_standard <- lapply(time_choices_exp_4, function(time) {
  print(paste('time:', time))
  return(exp_4_DL(N = 1000000,
                  input_samples = input_samples_exp_4,
                  mean = 0,
                  beta = 1/4,
                  precondition = 1,
                  time = time,
                  code = "R"))})

curve(tempered_target_density_exp_4(x, beta = 1/4), -4, 4)
lines(density(input_samples_exp_4), col = 'blue')
for (i in 1:length(time_choices_exp_4)) {
  print(paste('time:', time_choices_exp_4[i]))
  print(ks.test(test_exp_4_DL_R_standard[[i]], input_samples_exp_4))
  lines(density(test_exp_4_DL_R_standard[[i]]), col = 'green')
}

set.seed(21)
# applying the exact algorithm with double langevin [precondition] (USING R CODE)
test_exp_4_DL_R_precondition <- lapply(time_choices_exp_4, function(time) {
  print(paste('time:', time))
  return(exp_4_DL(N = 1000000,
                  input_samples = input_samples_exp_4,
                  mean = 0,
                  beta = 1/4,
                  precondition = var(input_samples_exp_4),
                  time = time,
                  code = "R"))})

curve(tempered_target_density_exp_4(x, beta = 1/4), -4, 4)
lines(density(input_samples_exp_4), col = 'blue')
for (i in 1:length(time_choices_exp_4)) {
  print(paste('time:', time_choices_exp_4[i]))
  print(ks.test(test_exp_4_DL_R_precondition[[i]], input_samples_exp_4))
  lines(density(test_exp_4_DL_R_precondition[[i]]), col = 'green')
}