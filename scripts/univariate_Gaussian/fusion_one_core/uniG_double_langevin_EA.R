mu <- 2.431
sd <- 1.2323
beta <- 0.79834
time_choices_gaussian <- c(0.1, 0.2, 0.5, 1, 2, 3, 4, 5)
seed <- 21
set.seed(seed)

# using rejection sampling to obtain input samples
input_samples_uniGaussian <- rnorm_tempered(N = 1000000, mean = mu, sd = sd, beta = beta)
curve(dnorm_tempered(x, mean = mu, sd = sd, beta = beta), -3, 8)
lines(density(input_samples_uniGaussian), col = 'blue')

set.seed(seed)
# applying the exact algorithm with double langevin [standard] (USING R CODE)
test_gaussian_DL_standard <- lapply(time_choices_gaussian, function(time) {
  print(paste('time:', time))
  return(ea_uniGaussian_DL(N = 1000000,
                           input_samples = input_samples_uniGaussian,
                           mean = mu,
                           sd = sd,
                           beta = beta,
                           precondition = 1,
                           time = time))})

curve(dnorm_tempered(x, mean = mu, sd = sd, beta = beta), -3, 8)
lines(density(input_samples_uniGaussian), col = 'blue')
for (i in 1:length(time_choices_gaussian)) {
  print(paste('time:', time_choices_gaussian[i]))
  print(ks.test(test_gaussian_DL_standard[[i]], input_samples_uniGaussian))
  lines(density(test_gaussian_DL_standard[[i]]), col = 'red')
}

set.seed(seed)
# applying the exact algorithm with double langevin [precondition] (USING R CODE)
test_gaussian_DL_precondition <- lapply(time_choices_gaussian, function(time) {
  print(paste('time:', time))
  return(ea_uniGaussian_DL(N = 1000000,
                           input_samples = input_samples_uniGaussian,
                           mean = mu,
                           sd = sd,
                           beta = beta,
                           precondition = var(input_samples_uniGaussian),
                           time = time))})

curve(dnorm_tempered(x, mean = mu, sd = sd, beta = beta), -3, 8)
lines(density(input_samples_uniGaussian), col = 'blue')
for (i in 1:length(time_choices_gaussian)) {
  print(paste('time:', time_choices_gaussian[i]))
  print(ks.test(test_gaussian_DL_precondition[[i]], input_samples_uniGaussian))
  lines(density(test_gaussian_DL_precondition[[i]]), col = 'red')
}
