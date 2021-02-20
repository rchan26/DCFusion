library(hierarchicalFusion)

seed <- 1994
set.seed(seed)

d2product_norm <- function(x, mu1, mu2, sd1, sd2, beta = 1) {
  # density for product of two Gaussians
  sd1 <- sd1/sqrt(beta)
  sd2 <- sd2/sqrt(beta)
  mu <- (sd1*sd1*mu2 + sd2*sd2*mu1) / (sd1*sd1 + sd2*sd2)
  std_dv <- sqrt(1 / ((1/(sd1*sd1)) + (1/(sd2*sd2))))
  return(dnorm(x, mean = mu, sd = std_dv))
}

d3product_norm <- function(x, mu1, mu2, mu3, sd1, sd2, sd3, beta = 1) {
  # density for product of three Gaussians
  sd1 <- sd1/sqrt(beta)
  sd2 <- sd2/sqrt(beta)
  sd3 <- sd3/sqrt(beta)
  mu12 <- (sd1*sd1*mu2 + sd2*sd2*mu1) / (sd1*sd1 + sd2*sd2)
  std_dv12 <- sqrt(1 / ((1/(sd1*sd1)) + (1/(sd2*sd2))))
  mu123 <- (std_dv12*std_dv12*mu3 + sd3*sd3*mu12) / (std_dv12*std_dv12 + sd3*sd3)
  std_dv123 <- sqrt(1 / ((1/(std_dv12*std_dv12)) + (1/(sd3*sd3))))
  return(dnorm(x, mean = mu123, sd = std_dv123))
}

new_parameters <- function(mu1, mu2, sd1, sd2, beta = 1) {
  # parameters when combining two Gaussians
  sd1 <- sd1/sqrt(beta)
  sd2 <- sd2/sqrt(beta)
  mu <- (sd1*sd1*mu2 + sd2*sd2*mu1) / (sd1*sd1 + sd2*sd2)
  std_dv <- sqrt(1 / ((1/(sd1*sd1)) + (1/(sd2*sd2))))
  return(c('mu' = mu, 'sd' = std_dv))
}

time_choice <- 0.5

# setting parameter values for each of the sub-posteriors
mean_values <- c(-1, 1, 8)
sd_values <- sqrt(c(1, 0.5, 4))

# obtain samples from pi1, pi2, pi3 before hand
sub_post_samples <- lapply(X = 1:3, function(i) rnorm(n = 10000, mean = mean_values[i], sd = sd_values[i]))

curve(dnorm(x, mean_values[1], sd_values[1]), -5, 15, ylim = c(0, 0.8), ylab = 'pdf', n = 10000, 
      col = '#003262', lwd = 2, lty = 2)
curve(dnorm(x, mean_values[2], sd_values[2]), add = T, n = 5000, col = '#a41034', lwd = 2, lty = 3)
curve(dnorm(x, mean_values[3], sd_values[3]), add = T, n = 10000, col = '#fcb515', lwd = 2, lty = 1)

lines(density(sub_post_samples[[1]]))
lines(density(sub_post_samples[[2]]))
lines(density(sub_post_samples[[3]]))

all_in_one_precondition <- parallel_fusion_uniGaussian(N = 10000, 
                                                       m = 3,
                                                       time = time_choice,
                                                       samples_to_fuse = sub_post_samples, 
                                                       means = mean_values, 
                                                       sds = sd_values, 
                                                       betas = rep(1, 3),
                                                       precondition_values = sapply(sub_post_samples, var),
                                                       seed = seed)

all_in_one_standard <- parallel_fusion_uniGaussian(N = 10000, 
                                                   m = 3,
                                                   time = time_choice,
                                                   samples_to_fuse = sub_post_samples, 
                                                   means = mean_values, 
                                                   sds = sd_values, 
                                                   betas = rep(1, 3),
                                                   precondition_values = rep(1, 3),
                                                   seed = seed)

curve(d3product_norm(x, mean_values[1], mean_values[2], mean_values[3], sd_values[1], sd_values[2], sd_values[3]), -5, 15, 
      ylim = c(0,1), ylab = 'pdf', n = 10000, col = 'green')
lines(density(all_in_one_precondition$samples))
lines(density(all_in_one_standard$samples))

##################################################
# (pi1*pi2)*pi3
curve(dnorm(x, mean_values[1], sd_values[1]), -5, 15, ylim = c(0,1), ylab = 'pdf', n = 10000, main = 'Method 1')
curve(dnorm(x, mean_values[2], sd_values[2]), add = T, n = 10000)

pi1pi2 <- parallel_fusion_uniGaussian(N = 10000, 
                                      m = 2,
                                      time = time_choice,
                                      samples_to_fuse = list(sub_post_samples[[1]], sub_post_samples[[2]]),
                                      means = c(mean_values[1], mean_values[2]), 
                                      sds = c(sd_values[1], sd_values[2]), 
                                      betas = rep(1, 2), 
                                      precondition_values = c(var(sub_post_samples[[1]]), var(sub_post_samples[[2]])),
                                      seed = seed)
method_1_new_parameters <- new_parameters(mean_values[1], mean_values[2], sd_values[1], sd_values[2])

#####

curve(d2product_norm(x, mean_values[1], mean_values[2], sd_values[1], sd_values[2]), add = T, n = 10000)
curve(dnorm(x, method_1_new_parameters['mu'], method_1_new_parameters['sd']), col = 'red', add = T, n = 10000)
lines(density(pi1pi2$samples), col = '#ffd146')

#####

# trying to perform fusion now on
# (pi1pi2)*pi3
curve(d2product_norm(x, mean_values[1], mean_values[2], sd_values[1], sd_values[2]), -5, 15, ylim = c(0,1), ylab = 'pdf', 
      n = 10000, main = 'Method 1', col = 'red', add = T)
curve(dnorm(x, mean_values[3], sd_values[3]), add = T, n = 10000, col = 'blue')

#####

fused_method_1 <- parallel_fusion_uniGaussian(N = 10000, 
                                              m = 2,
                                              time = time_choice,
                                              samples_to_fuse = list(pi1pi2$samples, sub_post_samples[[3]]),
                                              means = c(method_1_new_parameters['mu'], mean_values[3]),
                                              sds = c(method_1_new_parameters['sd'], sd_values[3]), 
                                              betas = rep(1, 2), 
                                              precondition_values = c(var(pi1pi2$samples), var(sub_post_samples[[3]])),
                                              seed = seed)

##################################################
# (pi1*pi3)*pi2
curve(dnorm(x, mean_values[1], sd_values[1]), -5, 15, ylim = c(0,1), ylab = 'pdf', n = 10000, main = 'Method 2')
curve(dnorm(x, mean_values[3], sd_values[3]), add = T)

pi1pi3 <- parallel_fusion_uniGaussian(N = 10000, 
                                      m = 2,
                                      time = time_choice, 
                                      samples_to_fuse = list(sub_post_samples[[1]], sub_post_samples[[3]]),
                                      means = c(mean_values[1], mean_values[3]), 
                                      sds = c(sd_values[1], sd_values[3]), 
                                      betas = rep(1, 2),
                                      precondition_values = c(var(sub_post_samples[[1]]), var(sub_post_samples[[3]])),
                                      seed = seed)
method_2_new_parameters <- new_parameters(mean_values[1], mean_values[3], sd_values[1], sd_values[3])

#####

curve(d2product_norm(x, mean_values[1], mean_values[3], sd_values[1], sd_values[3]), add = T)
curve(dnorm(x, method_2_new_parameters['mu'], method_2_new_parameters['sd']), col = 'red', add = T)
lines(density(pi1pi3$samples), col = '#ffd146')

#####
# trying to perform fusion now on
# (pi1pi3)*pi2
curve(d2product_norm(x, mean_values[1], mean_values[3], sd_values[1], sd_values[3]), -5, 15, ylim = c(0,1), ylab = 'pdf', 
      n = 10000, main = 'Method 2', col = 'red', add = T)
curve(dnorm(x, mean_values[2], sd_values[2]), add = T, n = 10000, col = 'blue')

#####

fused_method_2 <- parallel_fusion_uniGaussian(N = 10000, 
                                              m = 2,
                                              time = time_choice,
                                              samples_to_fuse = list(pi1pi3$samples, sub_post_samples[[2]]),
                                              means = c(method_2_new_parameters['mu'], mean_values[2]),
                                              sds = c(method_2_new_parameters['sd'], sd_values[2]), 
                                              betas = rep(1, 2),
                                              precondition_values = c(var(pi1pi3$samples), var(sub_post_samples[[2]])),
                                              seed = seed)

##################################################
# (pi2*pi3)*pi1
curve(dnorm(x, mean_values[2], sd_values[2]), -5, 15, ylim = c(0,1), ylab = 'pdf', n = 10000, main = 'Method 3')
curve(dnorm(x, mean_values[3], sd_values[3]), add = T)

pi2pi3 <- parallel_fusion_uniGaussian(N = 10000,
                                      m = 2,
                                      time = time_choice,
                                      samples_to_fuse = list(sub_post_samples[[2]], sub_post_samples[[3]]),
                                      means = c(mean_values[2], mean_values[3]), 
                                      sds = c(sd_values[2], sd_values[3]), 
                                      betas = rep(1, 2),
                                      precondition_values = c(var(sub_post_samples[[2]]), var(sub_post_samples[[3]])),
                                      seed = seed)
method_3_new_parameters <- new_parameters(mean_values[2], mean_values[3], sd_values[2], sd_values[3])

#####

curve(d2product_norm(x, mean_values[2], mean_values[3], sd_values[2], sd_values[3]), add = T)
curve(dnorm(x, method_3_new_parameters['mu'], method_3_new_parameters['sd']), col = 'red', add = T)
lines(density(pi2pi3$samples), col = '#ffd146')

#####

# trying to perform fusion now on
# (pi2pi3)*pi1
curve(d2product_norm(x, mean_values[2], mean_values[3], sd_values[2], sd_values[3]), -5, 15, ylim = c(0,1), ylab = 'pdf',
      n = 10000, main = 'Method 3', col = 'red', add = T)
curve(dnorm(x, mean_values[1], sd_values[1]), add = T, n = 10000, col = 'blue')

#####

fused_method_3 <- parallel_fusion_uniGaussian(N = 10000,
                                              m = 2,
                                              time = time_choice,
                                              samples_to_fuse = list(pi2pi3$samples, sub_post_samples[[1]]),
                                              means = c(method_3_new_parameters['mu'], mean_values[1]),
                                              sds = c(method_3_new_parameters['sd'], sd_values[1]), 
                                              betas = rep(1, 2),
                                              precondition_values = c(var(pi2pi3$samples), var(sub_post_samples[[1]])),
                                              seed = seed)

##################################################

curve(d3product_norm(x, mean_values[1], mean_values[2], mean_values[3], sd_values[1], sd_values[2], sd_values[3]), -5, 15, 
      ylim = c(0,1), ylab = 'pdf', n = 10000, col = 'green')
curve(dnorm(x, mean_values[1], sd_values[1]), add = T, n = 10000, col = '#003262')
curve(dnorm(x, mean_values[2], sd_values[2]), add = T, n = 10000, col = '#a41034')
curve(dnorm(x, mean_values[3], sd_values[3]), add = T, n = 10000, col = '#fcb515')
lines(density(fused_method_1$samples))
lines(density(fused_method_2$samples))
lines(density(fused_method_3$samples))
lines(density(all_in_one_precondition$samples))
lines(density(all_in_one_standard$samples))



#########################

print(paste('All in One Time (preconditioned):', all_in_one_precondition$time))
print(paste('All in One Time (standard):', all_in_one_standard$time))
print(paste('Method 1 Time:', pi1pi2$time, '+', fused_method_1$time, '=', pi1pi2$time + fused_method_1$time))
print(paste('Method 2 Time:', pi1pi3$time, '+', fused_method_2$time, '=', pi1pi3$time + fused_method_2$time))
print(paste('Method 3 Time:', pi2pi3$time, '+', fused_method_3$time, '=', pi2pi3$time + fused_method_3$time))


# combine the black line to get the red line (yellow line is the fusion samples)
# next combine the red line with the blue line



##################################################

# sub-posteriors
curve(dnorm(x, mean_values[1], sd_values[1]), -5, 15, ylim = c(0, 0.8), ylab = 'pdf', n = 10000, 
      col = '#003262', lwd = 2, lty = 2)
curve(dnorm(x, mean_values[2], sd_values[2]), add = T, n = 5000, col = '#a41034', lwd = 2, lty = 3)
curve(dnorm(x, mean_values[3], sd_values[3]), add = T, n = 10000, col = '#fcb515', lwd = 2, lty = 1)
legend(x = -5, y = 0.8, legend = c(expression(f[1]), 
                                   expression(f[2]), 
                                   expression(f[3])),
       col = c('#003262', '#a41034', '#fcb515'), lwd = 2, lty = c(2, 3, 1))

# Method 1
curve(dnorm(x, mean_values[1], sd_values[1]), -5, 15, ylim = c(0, 0.8), ylab = 'pdf', n = 1000, lty = 3, lwd = 2)
curve(dnorm(x, mean_values[2], sd_values[2]), add = T, n = 1000, lty = 3, lwd = 2)
curve(d2product_norm(x, mean_values[1], mean_values[2], sd_values[1], sd_values[2]), add = T, n = 10000, col = '#01653a', lty = 1, lwd = 2)
curve(dnorm(x, mean_values[3], sd_values[3]), add = T, n = 5000, col = '#9e292f', lty = 2, lwd = 2)
legend(x = -5, y = 0.8, legend = c(expression(f[1]), 
                                   expression(f[2]),
                                   expression(f[3]),
                                   expression(f[1]*f[2])),
       col = c('black', 'black', '#9e292f', '#01653a'), lwd = 2, lty = c(3, 3, 2, 1))


# Method 2
curve(dnorm(x, mean_values[1], sd_values[1]), -5, 15, ylim = c(0, 0.8), ylab = 'pdf', n = 1000, lty = 3, lwd = 2)
curve(dnorm(x, mean_values[3], sd_values[3]), add = T, n = 1000, lty = 3, lwd = 2)
curve(d2product_norm(x, mean_values[1], mean_values[3], sd_values[1], sd_values[3]), add = T, n = 10000, col = '#01653a', lty = 1, lwd = 2)
curve(dnorm(x, mean_values[2], sd_values[2]), add = T, n = 5000, col = '#9e292f', lty = 2, lwd = 2)
legend(x = -5, y = 0.8, legend = c(expression(f[1]), 
                                   expression(f[2]),
                                   expression(f[3]),
                                   expression(f[1]*f[3])),
       col = c('black', '#9e292f', 'black', '#01653a'), lwd = 2, lty = c(3, 2, 3, 1))


# Method 3
curve(dnorm(x, mean_values[2], sd_values[2]), -5, 15, ylim = c(0, 0.8), ylab = 'pdf', n = 1000, lty = 3, lwd = 2)
curve(dnorm(x, mean_values[3], sd_values[3]), add = T, n = 1000, lty = 3, lwd = 2)
curve(d2product_norm(x, mean_values[2], mean_values[3], sd_values[2], sd_values[3]), add = T, n = 10000, col = '#01653a', lty = 1, lwd = 2)
curve(dnorm(x, mean_values[1], sd_values[1]), add = T, n = 5000, col = '#9e292f', lty = 2, lwd = 2)
legend(x = -5, y = 0.8, legend = c(expression(f[1]), 
                                   expression(f[2]),
                                   expression(f[3]),
                                   expression(f[2]*f[3])),
       col = c('#9e292f', 'black', 'black', '#01653a'), lwd = 2, lty = c(2, 3, 3, 1))










