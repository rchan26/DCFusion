library(hierarchicalFusion)

d2product_norm <- function(x, mu1, mu2, sd1, sd2) {
  mu <- (sd1*sd1*mu2 + sd2*sd2*mu1) / (sd1*sd1 + sd2*sd2)
  std_dv <- sqrt(1 / ((1/(sd1*sd1)) + (1/(sd2*sd2))))
  return(dnorm(x, mean = mu, sd = std_dv))
}

d3product_norm <- function(x, mu1, mu2, mu3, sd1, sd2, sd3) {
  mu12 <- (sd1*sd1*mu2 + sd2*sd2*mu1) / (sd1*sd1 + sd2*sd2)
  std_dv12 <- sqrt(1 / ((1/(sd1*sd1)) + (1/(sd2*sd2))))
  mu123 <- (std_dv12*std_dv12*mu3 + sd3*sd3*mu12) / (std_dv12*std_dv12 + sd3*sd3)
  std_dv123 <- sqrt(1 / ((1/(std_dv12*std_dv12)) + (1/(sd3*sd3))))
  return(dnorm(x, mean = mu123, sd = std_dv123))
}

d2product_norm_parameters <- function(mu1, mu2, sd1, sd2) {
  mu <- (sd1*sd1*mu2 + sd2*sd2*mu1) / (sd1*sd1 + sd2*sd2)
  std_dv <- sqrt(1 / ((1/(sd1*sd1)) + (1/(sd2*sd2))))
  return(c('mu' = mu, 'sd' = std_dv))
}

d3product_norm_parameters <- function(mu1, mu2, mu3, sd1, sd2, sd3) {
  mu12 <- (sd1*sd1*mu2 + sd2*sd2*mu1) / (sd1*sd1 + sd2*sd2)
  std_dv12 <- sqrt(1 / ((1/(sd1*sd1)) + (1/(sd2*sd2))))
  mu123 <- (std_dv12*std_dv12*mu3 + sd3*sd3*mu12) / (std_dv12*std_dv12 + sd3*sd3)
  std_dv123 <- sqrt(1 / ((1/(std_dv12*std_dv12)) + (1/(sd3*sd3))))
  return(c('mu' = mu123, 'sd' = std_dv123))
}

seed <- 1994
set.seed(seed)
nsamples <- 10000
time_choice <- 0.5
# setting parameter values for each of the sub-posteriors
mean_values <- c(-10, 1, 8)
sd_values <- sqrt(c(1, 0.5, 4))
# target parameters
target_param <- d3product_norm_parameters(mu1 = mean_values[1],
                                          mu2 = mean_values[2],
                                          mu3 = mean_values[3],
                                          sd1 = sd_values[1],
                                          sd2 = sd_values[2],
                                          sd3 = sd_values[3])
# optimal bandwidth for Normal
opt_bw <- ((4/(3*nsamples))^(1/5))*target_param['sd']

# obtain samples from pi1, pi2, pi3 before hand
sub_post_samples <- lapply(X = 1:3, function(i) rnorm(n = nsamples, mean = mean_values[i], sd = sd_values[i]))
curve(dnorm(x, mean_values[1], sd_values[1]), min(mean_values)-4*max(sd_values), max(mean_values)+4*max(sd_values), ylim = c(0, 0.8), ylab = 'pdf', n = 10000, 
      col = '#003262', lwd = 2, lty = 2)
curve(dnorm(x, mean_values[2], sd_values[2]), add = T, n = 10000, col = '#a41034', lwd = 2, lty = 3)
curve(dnorm(x, mean_values[3], sd_values[3]), add = T, n = 10000, col = '#fcb515', lwd = 2, lty = 1)
lines(density(sub_post_samples[[1]]))
lines(density(sub_post_samples[[2]]))
lines(density(sub_post_samples[[3]]))

##### FORK-AND-JOIN #####

# fork-and-join fusion the sub-posterior samples (with preconditioning)
all_in_one <- parallel_fusion_SMC_uniGaussian(particles_to_fuse = initialise_particle_sets(sub_post_samples, FALSE),
                                              N = nsamples,
                                              m = 3,
                                              time = time_choice,
                                              means = mean_values,
                                              sds = sd_values,
                                              betas = rep(1, 3), 
                                              precondition_values = sapply(sub_post_samples, var), 
                                              ESS_threshold = 0.5,
                                              seed = seed)
# fork-and-join fusion the sub-posterior samples (without preconditioning)
all_in_one_standard <- parallel_fusion_SMC_uniGaussian(particles_to_fuse = initialise_particle_sets(sub_post_samples, FALSE),
                                                       N = nsamples,
                                                       m = 3,
                                                       time = time_choice,
                                                       means = mean_values,
                                                       sds = sd_values,
                                                       betas = rep(1, 3), 
                                                       precondition_values = rep(1, 3), 
                                                       ESS_threshold = 0.5,
                                                       seed = seed)

# plot
curve(d3product_norm(x, mean_values[1], mean_values[2], mean_values[3], sd_values[1], sd_values[2], sd_values[3]), min(mean_values)-4*max(sd_values), max(mean_values)+4*max(sd_values), 
      ylim = c(0,1), ylab = 'pdf', n = 10000, col = 'green')
lines(density(resample_particle_y_samples(particle_set = all_in_one$particles,
                                          multivariate = FALSE,
                                          resampling_method = 'resid',
                                          seed = seed)$y_samples,
              bw = opt_bw))
lines(density(resample_particle_y_samples(particle_set = all_in_one_standard$particles,
                                          multivariate = FALSE,
                                          resampling_method = 'resid',
                                          seed = seed)$y_samples,
              bw = opt_bw))

# IAD
print(integrated_abs_distance_uniGaussian(fusion_post = all_in_one$particles$y_samples,
                                          mean = target_param['mu'],
                                          sd = target_param['sd'],
                                          beta = 1,
                                          bw = opt_bw))
print(integrated_abs_distance_uniGaussian(fusion_post = all_in_one_standard$particles$y_samples,
                                          mean = target_param['mu'],
                                          sd = target_param['sd'],
                                          beta = 1,
                                          bw = opt_bw))

##### METHOD 1 #####
# (pi1*pi2)*pi3
curve(dnorm(x, mean_values[1], sd_values[1]), min(mean_values)-4*max(sd_values), max(mean_values)+4*max(sd_values), ylim = c(0,1), ylab = 'pdf', n = 10000, main = 'Method 1')
curve(dnorm(x, mean_values[2], sd_values[2]), add = T, n = 10000)

# part 1
samples_to_fuse <- list(sub_post_samples[[1]], sub_post_samples[[2]])
particles_to_fuse <- initialise_particle_sets(samples_to_fuse = samples_to_fuse,
                                              multivariate = FALSE)
pi1pi2 <- parallel_fusion_SMC_uniGaussian(particles_to_fuse = particles_to_fuse, 
                                          N = nsamples, 
                                          m = 2,
                                          time = time_choice, 
                                          means = c(mean_values[1], mean_values[2]), 
                                          sds = c(sd_values[1], sd_values[2]),
                                          betas = rep(1, 2), 
                                          precondition_values = sapply(samples_to_fuse, var), 
                                          ESS_threshold = 0, 
                                          seed = seed) 
method_1_new_parameters <- d2product_norm_parameters(mean_values[1], mean_values[2], sd_values[1], sd_values[2])

curve(d2product_norm(x, mean_values[1], mean_values[2], sd_values[1], sd_values[2]), add = T, n = 10000)
curve(dnorm(x, method_1_new_parameters['mu'], method_1_new_parameters['sd']), col = 'red', add = T, n = 10000)
lines(density(pi1pi2$particles$y_samples), col = '#ffd146')

# part 2 
# trying to perform fusion now on
# (pi1pi2)*pi3
curve(d2product_norm(x, mean_values[1], mean_values[2], sd_values[1], sd_values[2]), min(mean_values)-4*max(sd_values), max(mean_values)+4*max(sd_values), ylim = c(0,1), ylab = 'pdf', 
      n = 10000, main = 'Method 1', col = 'red', add = T)
curve(dnorm(x, mean_values[3], sd_values[3]), add = T, n = 10000, col = 'blue')

particles_to_fuse <- list(pi1pi2$particles, 
                          create_particle(samples = sub_post_samples[[3]], 
                                          multivariate = FALSE))
fused_method_1 <- parallel_fusion_SMC_uniGaussian(particles_to_fuse = particles_to_fuse,
                                                  N = nsamples, 
                                                  m = 2,
                                                  time = time_choice,
                                                  means = c(method_1_new_parameters['mu'], mean_values[3]),
                                                  sds = c(method_1_new_parameters['sd'], sd_values[3]),
                                                  betas = rep(1, 2), 
                                                  precondition_values = c(pi1pi2$precondition_values[[1]],
                                                                          var(sub_post_samples[[3]])),
                                                  ESS_threshold = 0.5,
                                                  seed = seed)

##### METHOD 2 #####
# (pi1*pi3)*pi2
curve(dnorm(x, mean_values[1], sd_values[1]), min(mean_values)-4*max(sd_values), max(mean_values)+4*max(sd_values), ylim = c(0,1), ylab = 'pdf', n = 10000, main = 'Method 2')
curve(dnorm(x, mean_values[3], sd_values[3]), add = T)

# part 1
samples_to_fuse <- list(sub_post_samples[[1]], sub_post_samples[[3]])
particles_to_fuse <- initialise_particle_sets(samples_to_fuse = samples_to_fuse,
                                              multivariate = FALSE)
pi1pi3 <- parallel_fusion_SMC_uniGaussian(particles_to_fuse = particles_to_fuse,
                                          N = nsamples, 
                                          m = 2,
                                          time = time_choice, 
                                          means = c(mean_values[1], mean_values[3]), 
                                          sds = c(sd_values[1], sd_values[3]), 
                                          betas = rep(1, 2),
                                          precondition_values = sapply(samples_to_fuse, var), 
                                          ESS_threshold = 0.5,
                                          seed = seed)
method_2_new_parameters <- d2product_norm_parameters(mean_values[1], mean_values[3], sd_values[1], sd_values[3])

curve(d2product_norm(x, mean_values[1], mean_values[3], sd_values[1], sd_values[3]), add = T)
curve(dnorm(x, method_2_new_parameters['mu'], method_2_new_parameters['sd']), col = 'red', add = T)
lines(density(pi1pi3$samples), col = '#ffd146')

# part 2
# trying to perform fusion now on
# (pi1pi3)*pi2
curve(d2product_norm(x, mean_values[1], mean_values[3], sd_values[1], sd_values[3]), min(mean_values)-4*max(sd_values), max(mean_values)+4*max(sd_values), ylim = c(0,1), ylab = 'pdf', 
      n = 10000, main = 'Method 2', col = 'red', add = T)
curve(dnorm(x, mean_values[2], sd_values[2]), add = T, n = 10000, col = 'blue')

particles_to_fuse <- list(pi1pi3$particles, 
                          create_particle(samples = sub_post_samples[[2]], 
                                          multivariate = FALSE))
fused_method_2 <- parallel_fusion_SMC_uniGaussian(particles_to_fuse = particles_to_fuse,
                                                  N = nsamples, 
                                                  m = 2,
                                                  time = time_choice,
                                                  means = c(method_2_new_parameters['mu'], mean_values[2]),
                                                  sds = c(method_2_new_parameters['sd'], sd_values[2]), 
                                                  betas = rep(1, 2), 
                                                  precondition_values = c(pi1pi3$precondition_values[[1]], var(sub_post_samples[[2]])),
                                                  ESS_threshold = 0.5,
                                                  seed = seed)

##### METHOD 3 #####
# (pi2*pi3)*pi1
curve(dnorm(x, mean_values[2], sd_values[2]), min(mean_values)-4*max(sd_values), max(mean_values)+4*max(sd_values), ylim = c(0,1), ylab = 'pdf', n = 10000, main = 'Method 3')
curve(dnorm(x, mean_values[3], sd_values[3]), add = T)

# part 1
samples_to_fuse <- list(sub_post_samples[[2]], sub_post_samples[[3]])
particles_to_fuse <- initialise_particle_sets(samples_to_fuse = samples_to_fuse,
                                              multivariate = FALSE)
pi2pi3 <- parallel_fusion_SMC_uniGaussian(particles_to_fuse = particles_to_fuse,
                                          N = nsamples,
                                          m = 2,
                                          time = time_choice,
                                          means = c(mean_values[2], mean_values[3]), 
                                          sds = c(sd_values[2], sd_values[3]), 
                                          betas = rep(1, 2),
                                          precondition_values = sapply(samples_to_fuse, var), 
                                          ESS_threshold = 0.5,
                                          seed = seed)
method_3_new_parameters <- d2product_norm_parameters(mean_values[2], mean_values[3], sd_values[2], sd_values[3])

curve(d2product_norm(x, mean_values[2], mean_values[3], sd_values[2], sd_values[3]), add = T)
curve(dnorm(x, method_3_new_parameters['mu'], method_3_new_parameters['sd']), col = 'red', add = T)
lines(density(pi2pi3$samples), col = '#ffd146')

# part 2
# trying to perform fusion now on
# (pi2pi3)*pi1
curve(d2product_norm(x, mean_values[2], mean_values[3], sd_values[2], sd_values[3]), min(mean_values)-4*max(sd_values), max(mean_values)+4*max(sd_values), ylim = c(0,1), ylab = 'pdf',
      n = 10000, main = 'Method 3', col = 'red', add = T)
curve(dnorm(x, mean_values[1], sd_values[1]), add = T, n = 10000, col = 'blue')

particles_to_fuse <- list(pi2pi3$particles, 
                          create_particle(samples = sub_post_samples[[1]], 
                                          multivariate = FALSE))
fused_method_3 <- parallel_fusion_SMC_uniGaussian(particles_to_fuse = particles_to_fuse,
                                                  N = nsamples,
                                                  m = 2,
                                                  time = time_choice,
                                                  means = c(method_3_new_parameters['mu'], mean_values[1]),
                                                  sds = c(method_3_new_parameters['sd'], sd_values[1]),
                                                  betas = rep(1, 2), 
                                                  precondition_values = c(pi2pi3$precondition_values[[1]], var(sub_post_samples[[1]])),
                                                  ESS_threshold = 0.5,
                                                  seed = seed)

##### PLOTS #####

curve(d3product_norm(x, mean_values[1], mean_values[2], mean_values[3], sd_values[1], sd_values[2], sd_values[3]), min(mean_values)-4*max(sd_values), max(mean_values)+4*max(sd_values), 
      ylim = c(0,1), ylab = 'pdf', n = 10000, col = 'green')
curve(dnorm(x, mean_values[1], sd_values[1]), add = T, n = 10000, col = '#003262')
curve(dnorm(x, mean_values[2], sd_values[2]), add = T, n = 10000, col = '#a41034')
curve(dnorm(x, mean_values[3], sd_values[3]), add = T, n = 10000, col = '#fcb515')
lines(density(resample_particle_y_samples(particle_set = fused_method_1$particles,
                                          multivariate = FALSE,
                                          resampling_method = 'resid',
                                          seed = seed)$y_samples,
              bw = opt_bw))
lines(density(resample_particle_y_samples(particle_set = fused_method_2$particles,
                                          multivariate = FALSE,
                                          resampling_method = 'resid',
                                          seed = seed)$y_samples,
              bw = opt_bw))
lines(density(resample_particle_y_samples(particle_set = fused_method_3$particles,
                                          multivariate = FALSE,
                                          resampling_method = 'resid',
                                          seed = seed)$y_samples,
              bw = opt_bw))
lines(density(resample_particle_y_samples(particle_set = all_in_one$particles,
                                          multivariate = FALSE,
                                          resampling_method = 'resid',
                                          seed = seed)$y_samples,
              bw = opt_bw))

##### TIME #####

print(paste('All in One Time (preconditioned):', all_in_one$time))
print(paste('All in One Time (standard):', all_in_one_standard$time))
print(paste('Method 1 Time:', pi1pi2$time, '+', fused_method_1$time, '=', pi1pi2$time + fused_method_1$time))
print(paste('Method 2 Time:', pi1pi3$time, '+', fused_method_2$time, '=', pi1pi3$time + fused_method_2$time))
print(paste('Method 3 Time:', pi2pi3$time, '+', fused_method_3$time, '=', pi2pi3$time + fused_method_3$time))

##### IAD #####

print(paste('All in One IAD (preconditioned):',
            integrated_abs_distance_uniGaussian(fusion_post = all_in_one$particles$y_samples,
                                                mean = target_param['mu'],
                                                sd = target_param['sd'],
                                                beta = 1,
                                                bw = opt_bw)))
print(paste('All in One IAD (preconditioned):',
            integrated_abs_distance_uniGaussian(fusion_post = all_in_one_standard$particles$y_samples,
                                                mean = target_param['mu'],
                                                sd = target_param['sd'],
                                                beta = 1,
                                                bw = opt_bw)))
print(paste('All in One IAD (preconditioned):',
            integrated_abs_distance_uniGaussian(fusion_post = fused_method_1$particles$y_samples,
                                                mean = target_param['mu'],
                                                sd = target_param['sd'],
                                                beta = 1,
                                                bw = opt_bw)))
print(paste('All in One IAD (preconditioned):',
            integrated_abs_distance_uniGaussian(fusion_post = fused_method_2$particles$y_samples,
                                                mean = target_param['mu'],
                                                sd = target_param['sd'],
                                                beta = 1,
                                                bw = opt_bw)))
print(paste('All in One IAD (preconditioned):',
            integrated_abs_distance_uniGaussian(fusion_post = fused_method_3$particles$y_samples,
                                                mean = target_param['mu'],
                                                sd = target_param['sd'],
                                                beta = 1,
                                                bw = opt_bw)))

##### PAPER PLOTS #####

# sub-posteriors
curve(dnorm(x, mean_values[1], sd_values[1]), min(mean_values)-4*max(sd_values), max(mean_values)+4*max(sd_values), ylim = c(0, 0.8), ylab = 'pdf', n = 10000, 
      col = '#003262', lwd = 2, lty = 2)
curve(dnorm(x, mean_values[2], sd_values[2]), add = T, n = 10000, col = '#a41034', lwd = 2, lty = 3)
curve(dnorm(x, mean_values[3], sd_values[3]), add = T, n = 10000, col = '#fcb515', lwd = 2, lty = 1)
legend(x = -5, y = 0.8, legend = c(expression(f[1]), 
                                   expression(f[2]), 
                                   expression(f[3])),
       col = c('#003262', '#a41034', '#fcb515'), lwd = 2, lty = c(2, 3, 1))

# Method 1
curve(dnorm(x, mean_values[1], sd_values[1]), min(mean_values)-4*max(sd_values), max(mean_values)+4*max(sd_values), ylim = c(0, 0.8), ylab = 'pdf', n = 1000, lty = 3, lwd = 2)
curve(dnorm(x, mean_values[2], sd_values[2]), add = T, n = 1000, lty = 3, lwd = 2)
curve(d2product_norm(x, mean_values[1], mean_values[2], sd_values[1], sd_values[2]), add = T, n = 10000, col = '#01653a', lty = 1, lwd = 2)
curve(dnorm(x, mean_values[3], sd_values[3]), add = T, n = 10000, col = '#9e292f', lty = 2, lwd = 2)
legend(x = -5, y = 0.8, legend = c(expression(f[1]), 
                                   expression(f[2]),
                                   expression(f[3]),
                                   expression(f[1]*f[2])),
       col = c('black', 'black', '#9e292f', '#01653a'), lwd = 2, lty = c(3, 3, 2, 1))


# Method 2
curve(dnorm(x, mean_values[1], sd_values[1]), min(mean_values)-4*max(sd_values), max(mean_values)+4*max(sd_values), ylim = c(0, 0.8), ylab = 'pdf', n = 1000, lty = 3, lwd = 2)
curve(dnorm(x, mean_values[3], sd_values[3]), add = T, n = 1000, lty = 3, lwd = 2)
curve(d2product_norm(x, mean_values[1], mean_values[3], sd_values[1], sd_values[3]), add = T, n = 10000, col = '#01653a', lty = 1, lwd = 2)
curve(dnorm(x, mean_values[2], sd_values[2]), add = T, n = 10000, col = '#9e292f', lty = 2, lwd = 2)
legend(x = -5, y = 0.8, legend = c(expression(f[1]), 
                                   expression(f[2]),
                                   expression(f[3]),
                                   expression(f[1]*f[3])),
       col = c('black', '#9e292f', 'black', '#01653a'), lwd = 2, lty = c(3, 2, 3, 1))


# Method 3
curve(dnorm(x, mean_values[2], sd_values[2]), min(mean_values)-4*max(sd_values), max(mean_values)+4*max(sd_values), ylim = c(0, 0.8), ylab = 'pdf', n = 1000, lty = 3, lwd = 2)
curve(dnorm(x, mean_values[3], sd_values[3]), add = T, n = 1000, lty = 3, lwd = 2)
curve(d2product_norm(x, mean_values[2], mean_values[3], sd_values[2], sd_values[3]), add = T, n = 10000, col = '#01653a', lty = 1, lwd = 2)
curve(dnorm(x, mean_values[1], sd_values[1]), add = T, n = 10000, col = '#9e292f', lty = 2, lwd = 2)
legend(x = -5, y = 0.8, legend = c(expression(f[1]), 
                                   expression(f[2]),
                                   expression(f[3]),
                                   expression(f[2]*f[3])),
       col = c('#9e292f', 'black', 'black', '#01653a'), lwd = 2, lty = c(2, 3, 3, 1))
