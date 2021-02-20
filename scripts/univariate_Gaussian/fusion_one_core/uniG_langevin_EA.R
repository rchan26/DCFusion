# ----- EA for univariate Gaussian

ea_phi_gaussian <- function(x, mean, sd, beta) {
  y <- x - mean
  return(0.125*(((beta*beta*y*y)/(sd^4))-((2*beta)/(sd*sd))))
}

ea_phi_gaussian_bounds <- function(mean, sd, beta, lower, upper) {
  x <- c(lower, upper)
  if (mean > lower & mean < upper) {
    x <- c(x, mean)
  }
  phi <- ea_phi_gaussian(x = x, mean = mean, sd = sd, beta = beta)
  return(list('LB' = min(phi), 'UB' = max(phi)))
}

ea_phi_gaussian_LB <- function(mean, sd, beta) {
  return(ea_phi_gaussian(x = mean,
                         mean = mean,
                         sd = sd,
                         beta = beta))
}

ea_gaussian_PT <- function(x0,
                           y,
                           s,
                           t,
                           mean,
                           sd,
                           beta) {
  # simulate layer information
  bes_layer <- layeredBB::bessel_layer_simulation(x = x0,
                                                  y = y,
                                                  s = s,
                                                  t = t,
                                                  mult = 0.1)
  lbound_X <- bes_layer$L
  ubound_X <- bes_layer$U
  # calculate lower and upper bounds of phi
  bounds <- ea_phi_gaussian_bounds(mean = mean,
                                   sd = sd,
                                   beta = beta,
                                   lower = lbound_X,
                                   upper = ubound_X)
  LX <- bounds$LB
  UX <- bounds$UB
  # calculate global lower bound of phi
  PHI <- ea_phi_gaussian_LB(mean = mean,
                            sd = sd,
                            beta = beta)
  # simulate the number of points to simulate from Possion distribution
  kap <- rpois(n = 1, lambda = (UX-LX)*(t-s))
  log_acc_prob <- 0
  # if kap > 0, then simulate layered Brownian Bridge
  if (kap > 0) {
    pois_times <- runif(kap, s, t)
    layered_bb <- layeredBB::layered_brownian_bridge(x = x0,
                                                     y = y,
                                                     s = s,
                                                     t = t,
                                                     bessel_layer = bes_layer,
                                                     times = pois_times)
    simulated_path <- layered_bb$simulated_path
    phi <- ea_phi_gaussian(x = simulated_path[1,], mean = mean, sd = sd, beta = beta)
    log_acc_prob <- sum(log(UX-phi))
  }
  return(-(LX-PHI)*(t-s) - kap*log(UX-LX) + log_acc_prob)
}

ea_gaussian_importance <- function(N, input_samples, mean, sd, beta, time) {
  samples <- rep(NA, N)
  weights <- rep(NA, N)
  i <- 0
  while (i < N) {
    # sample x ~ pi
    x <- sample(input_samples, 1)
    # sample y ~ N(x, T)
    y <- rnorm(n = 1, mean = x, sd = sqrt(time))
    # poisson thinning to evaluate acceptance probability
    log_acceptance <- ea_gaussian_PT(x0 = x, 
                                     y = y,
                                     s = 0,
                                     t = time,
                                     mean = mean,
                                     sd = sd,
                                     beta = beta)
    if (log(runif(1, 0, 1)) < log_acceptance) {
      i <- i+1
      samples[i] <- y
      weights[i] <- dnorm_tempered(x = y, mean = mean, sd = sd, beta = beta/2) /
        dnorm_tempered(x = x, mean = mean, sd = sd, beta = beta/2)
    }
  }
  return(list('samples' = samples, 'weights' = weights))
}

mu <- 2.431
sd <- 1.2323
beta <- 0.79834
time_choices_gaussian <- c(0.1, 0.2, 0.5, 1, 2, 3, 4, 5, 7.5)
set.seed(21)

# using rejection sampling to obtain input samples
input_samples_uniGaussian <- rnorm_tempered(N = 1000000, mean = mu, sd = sd, beta = beta)
curve(dnorm_tempered(x, mean = mu, sd = sd, beta = beta), -3, 8)
lines(density(input_samples_uniGaussian), col = 'blue')

set.seed(21)
# applying the exact algorithm (importance)
test_ea_gaussian_importance <- lapply(time_choices_gaussian, function(time) {
  ea_gaussian_importance(N = 1000000,
                         input_samples = input_samples_uniGaussian,
                         mean = mu,
                         sd = sd,
                         beta = beta,
                         time = time)})

curve(dnorm_tempered(x, mean = mu, sd = sd, beta = beta), -3, 8)
for (i in 1:length(time_choices_gaussian)) {
  lines(density(test_ea_gaussian_importance[[i]]$samples), col = 'red')
}

curve(dnorm_tempered(x, mean = mu, sd = sd, beta = beta), -3, 8)
lines(density(input_samples_uniGaussian[1:750000]), col = 'blue')
for (i in 1:length(time_choices_gaussian)) {
  resampled <- sample(x = test_ea_gaussian_importance[[i]]$samples,
                      size = 1000000,
                      replace = TRUE,
                      prob = test_ea_gaussian_importance[[i]]$weights/sum(test_ea_gaussian_importance[[i]]$weights))
  lines(density(resampled), col = 'green')
}

for (i in 1:length(time_choices_gaussian)) {
  norm_weights <- test_ea_gaussian_importance[[i]]$weights / sum(test_ea_gaussian_importance[[i]]$weights)
  ESS <- floor(1/sum(norm_weights^2))
  print(ESS)
  # resampled <- sample(test_ea_gaussian_importance[[i]]$samples, size = ESS, replace = T, prob = norm_weights)
  # samples_to_compare <- sample(input_samples_uniGaussian, size = ESS)
  # print(paste('time:', time_choices_gaussian[i]))
  # print(ks.test(resampled, samples_to_compare))
}