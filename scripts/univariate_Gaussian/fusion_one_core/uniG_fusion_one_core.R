mu <- 2.431
sd <- 1.2323
beta <- 0.79834
time_choices_gaussian <- c(0.1, 0.2, 0.5, 1, 2, 3, 4, 5, 7.5)
seed <- 21
set.seed(seed)

# using rejection sampling to obtain input samples
input_samples_uniGaussian <- rnorm_tempered(N = 1000000, mean = mu, sd = sd, beta = beta)
curve(dnorm_tempered(x, mean = mu, sd = sd, beta = beta), -3, 8)
lines(density(input_samples_uniGaussian), col = 'blue')

########## USING R CODE [standard]
test_gaussian_fusion_1_core_standard <- lapply(time_choices_gaussian, function(time) {
  print(paste('time:', time))
  return(parallel_fusion_uniGaussian(N = 1000000,
                                     m = 1,
                                     time = time,
                                     samples_to_fuse = list(input_samples_uniGaussian),
                                     means = mu,
                                     sds = sd,
                                     betas = beta,
                                     precondition_values = 1,
                                     seed = seed))})

curve(dnorm_tempered(x, mean = mu, sd = sd, beta = beta), -3, 8)
lines(density(input_samples_uniGaussian), col = 'blue')
for (i in 1:length(time_choices_gaussian)) {
  print(paste('time:', time_choices_gaussian[i]))
  print(ks.test(test_gaussian_fusion_1_core_standard[[i]]$samples, input_samples_uniGaussian))
  lines(density(test_gaussian_fusion_1_core_standard[[i]]$samples), col = 'red')
}
lines(density(input_samples_uniGaussian), col = 'blue')
curve(dnorm_tempered(x, mean = mu, sd = sd, beta = beta), -3, 8, add = T)

########## USING R CODE [precondition]
test_gaussian_fusion_1_core_precondition <- lapply(time_choices_gaussian, function(time) {
  print(paste('time:', time))
  return(parallel_fusion_uniGaussian(N = 1000000,
                                     m = 1,
                                     time = time,
                                     samples_to_fuse = list(input_samples_uniGaussian),
                                     means = mu,
                                     sds = sd,
                                     betas = beta,
                                     precondition_values = var(input_samples_uniGaussian),
                                     seed = seed))})

curve(dnorm_tempered(x, mean = mu, sd = sd, beta = beta), -3, 8)
for (i in 1:length(time_choices_gaussian)) {
  print(paste('time:', time_choices_gaussian[i]))
  print(ks.test(test_gaussian_fusion_1_core_precondition[[i]]$samples, input_samples_uniGaussian))
  lines(density(test_gaussian_fusion_1_core_precondition[[i]]$samples), col = 'red')
}
lines(density(input_samples_uniGaussian), col = 'blue')
curve(dnorm_tempered(x, mean = mu, sd = sd, beta = beta), -3, 8, add = T)
