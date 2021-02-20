# ----- EA for exp_4

ea_phi_exp_4 <- function(x, mean, beta) {
  y <- x - mean
  return(0.5*((beta*beta*(y^6)) - (3*beta*y*y)))
}

ea_phi_exp_4_bounds <- function(mean, beta, lower, upper) {
  x <- c(lower, upper)
  if (mean > lower & mean < upper) {
    x <- c(x, mean)
  }
  m1 <- mean - (1/beta)^(0.25)
  m2 <- mean + (1/beta)^(0.25)
  if (m1 > lower & m1 < upper) {
    x <- c(x, m1)
  } 
  if (m2 > lower & m2 < upper) {
    x <- c(x, m2)
  }
  phi <- ea_phi_exp_4(x = x, mean = mean, beta = beta)
  return(list('LB' = min(phi), 'UB' = max(phi)))
}

ea_phi_exp_4_LB <- function(mean, beta) {
  x <- c(mean - (1/beta)^(0.25),
         mean + (1/beta)^(0.25))
  return(min(ea_phi_exp_4(x, mean = mean, beta = beta)))
}

ea_exp_4_PT <- function(x0,
                        y,
                        s,
                        t,
                        mean,
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
  bounds <- ea_phi_exp_4_bounds(mean = mean,
                                beta = beta,
                                lower = lbound_X,
                                upper = ubound_X)
  LX <- bounds$LB
  UX <- bounds$UB
  # calculate global lower bound of phi
  PHI <- ea_phi_exp_4_LB(mean = mean, beta = beta)
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
    phi <- ea_phi_exp_4(x = simulated_path[1,], mean = mean, beta = beta)
    log_acc_prob <- sum(log(UX-phi))
  }
  return(-(LX-PHI)*(t-s) - kap*log(UX-LX) + log_acc_prob)
}

ea_exp_4_importance <- function(N, input_samples, mean, beta, time) {
  samples <- rep(NA, N)
  weights <- rep(NA, N)
  i <- 0
  while (i < N) {
    # sample x ~ pi
    x <- sample(input_samples, 1)
    # sample y ~ N(x, T)
    y <- rnorm(n = 1, mean = x, sd = sqrt(time))
    # poisson thinning to evaluate log acceptance probability
    log_acceptance <- ea_exp_4_PT(x0 = x, 
                                  y = y,
                                  s = 0,
                                  t = time,
                                  mean = mean,
                                  beta = beta)
    if (log(runif(1, 0, 1)) < log_acceptance) {
      i <- i+1
      samples[i] <- y
      weights[i] <- tempered_target_density_exp_4(x = y, mean = mean, beta = beta/2) /
        tempered_target_density_exp_4(x = x, mean = mean, beta = beta/2)
    }
  }
  return(list('samples' = samples, 'weights' = weights))
}

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
# applying the exact algorithm (importance)
test_ea_exp_4_importance <- lapply(time_choices_exp_4, function(time) {
  ea_exp_4_importance(N = 1000000,
                      input_samples = input_samples_exp_4,
                      mean = 0,
                      beta = 1/4,
                      time = time)})

curve(exp_4_density(x, beta = 1/4), -4, 4)
lines(density(input_samples_exp_4), col = 'blue')
for (i in 1:length(time_choices_exp_4)) {
  lines(density(test_ea_exp_4_importance[[i]]$samples), col = 'red')
}

curve(exp_4_density(x, beta = 1/4), -4, 4)
lines(density(input_samples_exp_4[1:900000]), col = 'blue', lwd = 2)
for (i in 1:length(time_choices_exp_4)) {
  resampled <- sample(x = test_ea_exp_4_importance[[i]]$samples,
                      size = 1000000,
                      replace = TRUE,
                      prob = test_ea_exp_4_importance[[i]]$weights/sum(test_ea_exp_4_importance[[i]]$weights))
  lines(density(resampled, bw = 0.08981), col = 'green')
}

for (i in 1:length(time_choices_exp_4)) {
  norm_weights <- test_ea_exp_4_importance[[i]]$weights / sum(test_ea_exp_4_importance[[i]]$weights)
  ESS <- floor(1/sum(norm_weights^2))
  print(ESS)
  # resampled <- sample(test_ea_exp_4_importance[[i]]$samples, size = ESS, replace = T, prob = norm_weights)
  # samples_to_compare <- sample(input_samples_exp_4, size = ESS)
  # print(paste('time:', time_choices_exp_4[i]))
  # print(ks.test(resampled, samples_to_compare))
}
