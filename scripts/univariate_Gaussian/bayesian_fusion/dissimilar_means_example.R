library(DCFusion)

seed <- 1994
set.seed(seed)
nsamples <- 10000
C <- 2
means <- c(-0.25, 0.25)
beta <- 1/C
a_mesh_vanilla <- seq(0, 0.01, length.out = 6)
a_mesh_gen <- seq(0, 1, length.out = 6)
diffusion_estimator <- 'NB'
ESS_threshold <- 0.5
CESS_0_threshold <- 0.2
vanilla_b <- 1
k1 <- NULL
k2 <- NULL
k3 <- 1
k4 <- 1
data_sizes <- c(250, 500, 1000, 1500, 2000, 2500)
a_results <- list('vanilla' = list(), 'generalised' = list())
b_results <- list('vanilla' = list(), 'generalised' = list())
c_results <- list('vanilla' = list(), 'generalised' = list())
d_results <- list('vanilla' = list(), 'generalised' = list())
SH_adaptive_results <- list('vanilla' = list(), 'generalised' = list())

for (i in 1:length(data_sizes)) {
  set.seed(seed*i)
  sd <- sqrt(1/data_sizes[i])
  curve(dnorm(x, mean = 0, sd = sd), -4*sd, +4*sd,
        main = paste('m:', data_sizes[i]))
  input_samples <- lapply(1:C, function(i) rnorm_tempered(N = nsamples,
                                                          mean = means[i],
                                                          sd = sd,
                                                          beta = beta))
  print(paste('##### m:', data_sizes[i]))
  print(paste('sd:', sd))
  print(paste('var:', sd^2))
  print(paste('sub_posterior var (true):', C*sd^2))
  print(paste('sub_posterior var:', sapply(input_samples, var)))
  ##### Fixed user-specified parameters #####
  print('### performing standard Bayesian Fusion (with standard mesh)')
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = FALSE,
                                              number_of_steps = length(a_mesh_vanilla))
  a_BF_standard <- parallel_GBF_uniGaussian(particles_to_fuse = input_particles,
                                            N = nsamples,
                                            m = C,
                                            time_mesh = a_mesh_vanilla,
                                            means = means,
                                            sds = rep(sd, C),
                                            betas = rep(beta, C),
                                            precondition_values = rep(1, C),
                                            ESS_threshold = ESS_threshold,
                                            diffusion_estimator = diffusion_estimator,
                                            seed = seed*i)
  print('### performing Bayesian Fusion with a preconditioning matrix (with standard mesh)')
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = FALSE,
                                              number_of_steps = length(a_mesh_gen))
  a_BF_generalised <- parallel_GBF_uniGaussian(particles_to_fuse = input_particles,
                                               N = nsamples,
                                               m = C,
                                               time_mesh = a_mesh_gen,
                                               means = means,
                                               sds = rep(sd, C),
                                               betas = rep(beta, C),
                                               precondition_values = sapply(input_samples, var),
                                               ESS_threshold = ESS_threshold,
                                               diffusion_estimator = diffusion_estimator,
                                               seed = seed*i)
  # save results
  a_results$vanilla[[i]] <- list('CESS_0' = a_BF_standard$CESS[1],
                                 'CESS_j' = a_BF_standard$CESS[2:length(a_BF_standard$CESS)],
                                 'CESS_j_avg' = mean(a_BF_standard$CESS[2:length(a_BF_standard$CESS)]),
                                 'n' = length(a_BF_standard$CESS),
                                 'time_mesh' = a_BF_standard$particles$time_mesh,
                                 'time' = a_BF_standard$time)
  a_results$generalised[[i]] <- list('CESS_0' = a_BF_generalised$CESS[1],
                                     'CESS_j' = a_BF_generalised$CESS[2:length(a_BF_generalised$CESS)],
                                     'CESS_j_avg' = mean(a_BF_generalised$CESS[2:length(a_BF_generalised$CESS)]),
                                     'n' = length(a_BF_generalised$CESS),
                                     'time_mesh' = a_BF_generalised$particles$time_mesh,
                                     'time' = a_BF_generalised$time)
  # plot
  lines(density(resample_particle_y_samples(particle_set = a_BF_standard$particles, 
                                            multivariate = FALSE,
                                            resampling_method = 'resid',
                                            seed = seed*i)$y_samples), 
        col = 'red')
  lines(density(resample_particle_y_samples(particle_set = a_BF_generalised$particles, 
                                            multivariate = FALSE,
                                            resampling_method = 'resid',
                                            seed = seed*i)$y_samples), 
        col = 'red', lty = 2)
  
  ##### Recommended scaling of T, fixed n #####
  print('### performing standard Bayesian Fusion (with recommended T, fixed n)')
  vanilla_guide <- BF_guidance(condition = 'SSH',
                               CESS_0_threshold = CESS_0_threshold,
                               C = C,
                               d = 1,
                               data_size = data_sizes[i],
                               b = vanilla_b,
                               sub_posterior_means = sapply(input_samples, mean),
                               k1 = k1,
                               k2 = k2,
                               k3 = k3,
                               k4 = k4,
                               vanilla = TRUE)
  b_mesh_vanilla <- seq(0, vanilla_guide$min_T, length.out = 6)
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = FALSE,
                                              number_of_steps = length(b_mesh_vanilla))
  b_BF_standard <- parallel_GBF_uniGaussian(particles_to_fuse = input_particles,
                                            N = nsamples,
                                            m = C,
                                            time_mesh = b_mesh_vanilla,
                                            means = means,
                                            sds = rep(sd, C),
                                            betas = rep(beta, C),
                                            precondition_values = rep(1, C),
                                            ESS_threshold = ESS_threshold,
                                            diffusion_estimator = diffusion_estimator,
                                            seed = seed*i)
  print('### performing Bayesian Fusion with a preconditioning matrix (with recommended T, fixed n)')
  gen_guide <- BF_guidance(condition = 'SSH',
                           CESS_0_threshold = CESS_0_threshold,
                           C = C,
                           d = 1,
                           data_size = data_sizes[i],
                           sub_posterior_means = sapply(input_samples, mean),
                           precondition_matrices = sapply(input_samples, var),
                           inv_precondition_matrices = 1/sapply(input_samples, var),
                           k1 = k1,
                           k2 = k2,
                           k3 = k3,
                           k4 = k4,
                           vanilla = FALSE)
  b_mesh_gen <- seq(0, gen_guide$min_T, length.out = 6)
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = FALSE,
                                              number_of_steps = length(b_mesh_gen))
  b_BF_generalised <- parallel_GBF_uniGaussian(particles_to_fuse = input_particles,
                                               N = nsamples,
                                               m = C,
                                               time_mesh = b_mesh_gen,
                                               means = means,
                                               sds = rep(sd, C),
                                               betas = rep(beta, C),
                                               precondition_values = sapply(input_samples, var),
                                               ESS_threshold = ESS_threshold,
                                               diffusion_estimator = diffusion_estimator,
                                               seed = seed*i)
  # save results
  b_results$vanilla[[i]] <- list('CESS_0' = b_BF_standard$CESS[1],
                                 'CESS_j' = b_BF_standard$CESS[2:length(b_BF_standard$CESS)],
                                 'CESS_j_avg' = mean(b_BF_standard$CESS[2:length(b_BF_standard$CESS)]),
                                 'n' = length(b_BF_standard$CESS),
                                 'time_mesh' = b_BF_standard$particles$time_mesh,
                                 'time' = b_BF_standard$time)
  b_results$generalised[[i]] <- list('CESS_0' = b_BF_generalised$CESS[1],
                                     'CESS_j' = b_BF_generalised$CESS[2:length(b_BF_generalised$CESS)],
                                     'CESS_j_avg' = mean(b_BF_generalised$CESS[2:length(b_BF_generalised$CESS)]),
                                     'n' = length(b_BF_generalised$CESS),
                                     'time_mesh' = b_BF_generalised$particles$time_mesh,
                                     'time' = b_BF_generalised$time)
  # plot
  lines(density(resample_particle_y_samples(particle_set = b_BF_standard$particles, 
                                            multivariate = FALSE,
                                            resampling_method = 'resid',
                                            seed = seed*i)$y_samples), 
        col = 'blue')
  lines(density(resample_particle_y_samples(particle_set = b_BF_generalised$particles, 
                                            multivariate = FALSE,
                                            resampling_method = 'resid',
                                            seed = seed*i)$y_samples), 
        col = 'blue', lty = 2)
  
  ##### Recommended scaling of T, regular mesh #####
  print('### performing standard Bayesian Fusion (with recommended T, regular mesh)')
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = FALSE,
                                              number_of_steps = length(vanilla_guide$mesh))
  c_BF_standard <- parallel_GBF_uniGaussian(particles_to_fuse = input_particles,
                                            N = nsamples,
                                            m = C,
                                            time_mesh = vanilla_guide$mesh,
                                            means = means,
                                            sds = rep(sd, C),
                                            betas = rep(beta, C),
                                            precondition_values = rep(1, C),
                                            ESS_threshold = ESS_threshold,
                                            diffusion_estimator = diffusion_estimator,
                                            seed = seed*i)
  print('### performing Bayesian Fusion with a preconditioning matrix (with recommended T, regular mesh)')
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = FALSE,
                                              number_of_steps = length(gen_guide$mesh))
  c_BF_generalised <- parallel_GBF_uniGaussian(particles_to_fuse = input_particles,
                                               N = nsamples,
                                               m = C,
                                               time_mesh = gen_guide$mesh,
                                               means = means,
                                               sds = rep(sd, C),
                                               betas = rep(beta, C),
                                               precondition_values = sapply(input_samples, var),
                                               ESS_threshold = ESS_threshold,
                                               diffusion_estimator = diffusion_estimator,
                                               seed = seed*i)
  # save results
  c_results$vanilla[[i]] <- list('CESS_0' = c_BF_standard$CESS[1],
                                 'CESS_j' = c_BF_standard$CESS[2:length(c_BF_standard$CESS)],
                                 'CESS_j_avg' = mean(c_BF_standard$CESS[2:length(c_BF_standard$CESS)]),
                                 'n' = length(c_BF_standard$CESS),
                                 'time_mesh' = c_BF_standard$particles$time_mesh,
                                 'time' = c_BF_standard$time)
  c_results$generalised[[i]] <- list('CESS_0' = c_BF_generalised$CESS[1],
                                     'CESS_j' = c_BF_generalised$CESS[2:length(c_BF_generalised$CESS)],
                                     'CESS_j_avg' = mean(c_BF_generalised$CESS[2:length(c_BF_generalised$CESS)]),
                                     'n' = length(c_BF_generalised$CESS),
                                     'time_mesh' = c_BF_generalised$particles$time_mesh,
                                     'time' = c_BF_generalised$time)
  # plot
  lines(density(resample_particle_y_samples(particle_set = c_BF_standard$particles, 
                                            multivariate = FALSE,
                                            resampling_method = 'resid',
                                            seed = seed*i)$y_samples), 
        col = 'green')
  lines(density(resample_particle_y_samples(particle_set = c_BF_generalised$particles, 
                                            multivariate = FALSE,
                                            resampling_method = 'resid',
                                            seed = seed*i)$y_samples), 
        col = 'green', lty = 2)
  
  ##### Recommended scaling of T, adaptive mesh #####
  print('### performing standard Bayesian Fusion (with recommended T, adaptive mesh)')
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = FALSE,
                                              number_of_steps = length(vanilla_guide$mesh))
  d_BF_standard <- parallel_GBF_uniGaussian(particles_to_fuse = input_particles,
                                            N = nsamples,
                                            m = C,
                                            time_mesh = vanilla_guide$mesh,
                                            means = means,
                                            sds = rep(sd, C),
                                            betas = rep(beta, C),
                                            precondition_values = rep(1, C),
                                            ESS_threshold = ESS_threshold,
                                            sub_posterior_means = sapply(input_samples, mean),
                                            adaptive_mesh = TRUE,
                                            adaptive_mesh_parameters = list('data_size' = data_sizes[i],
                                                                            'b' = vanilla_b,
                                                                            'k3' = k3,
                                                                            'k4' = k4,
                                                                            'vanilla' = TRUE),
                                            diffusion_estimator = diffusion_estimator,
                                            seed = seed*i)
  print('### performing Bayesian Fusion with a preconditioning matrix (with recommended T, adaptive mesh)')
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = FALSE,
                                              number_of_steps = length(gen_guide$mesh))
  d_BF_generalised <- parallel_GBF_uniGaussian(particles_to_fuse = input_particles,
                                               N = nsamples,
                                               m = C,
                                               time_mesh = gen_guide$mesh,
                                               means = means,
                                               sds = rep(sd, C),
                                               betas = rep(beta, C),
                                               precondition_values = sapply(input_samples, var),
                                               ESS_threshold = ESS_threshold,
                                               sub_posterior_means = sapply(input_samples, mean),
                                               adaptive_mesh = TRUE,
                                               adaptive_mesh_parameters = list('data_size' = data_sizes[i],
                                                                               'k3' = k3,
                                                                               'k4' = k4,
                                                                               'vanilla' = FALSE),
                                               diffusion_estimator = diffusion_estimator,
                                               seed = seed*i)
  # save results
  d_results$vanilla[[i]] <- list('CESS_0' = d_BF_standard$CESS[1],
                                 'CESS_j' = d_BF_standard$CESS[2:length(d_BF_standard$CESS)],
                                 'CESS_j_avg' = mean(d_BF_standard$CESS[2:length(d_BF_standard$CESS)]),
                                 'n' = length(d_BF_standard$CESS),
                                 'time_mesh' = d_BF_standard$particles$time_mesh,
                                 'time' = d_BF_standard$time)
  d_results$generalised[[i]] <- list('CESS_0' = d_BF_generalised$CESS[1],
                                     'CESS_j' = d_BF_generalised$CESS[2:length(d_BF_generalised$CESS)],
                                     'CESS_j_avg' = mean(d_BF_generalised$CESS[2:length(d_BF_generalised$CESS)]),
                                     'n' = length(d_BF_generalised$CESS),
                                     'time_mesh' = d_BF_generalised$particles$time_mesh,
                                     'time' = d_BF_generalised$time)
  # plot
  lines(density(resample_particle_y_samples(particle_set = d_BF_standard$particles, 
                                            multivariate = FALSE,
                                            resampling_method = 'resid',
                                            seed = seed*i)$y_samples), 
        col = 'purple')
  lines(density(resample_particle_y_samples(particle_set = d_BF_generalised$particles, 
                                            multivariate = FALSE,
                                            resampling_method = 'resid',
                                            seed = seed*i)$y_samples), 
        col = 'purple', lty = 2)
  
  ##### SH: Recommended scaling of T, adaptive mesh #####
  print('### SH: performing standard Bayesian Fusion (with recommended T, adaptive mesh)')
  vanilla_guide_SH <- BF_guidance(condition = 'SH',
                                  CESS_0_threshold = CESS_0_threshold,
                                  C = C,
                                  d = 1,
                                  data_size = data_sizes[i],
                                  b = vanilla_b,
                                  sub_posterior_means = sapply(input_samples, mean),
                                  k1 = k1,
                                  k3 = k3,
                                  k4 = k4,
                                  vanilla = TRUE)
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = FALSE,
                                              number_of_steps = length(vanilla_guide_SH$mesh))
  SH_adaptive_standard <- parallel_GBF_uniGaussian(particles_to_fuse = input_particles,
                                                   N = nsamples,
                                                   m = C,
                                                   time_mesh = vanilla_guide_SH$mesh,
                                                   means = means,
                                                   sds = rep(sd, C),
                                                   betas = rep(beta, C),
                                                   precondition_values = rep(1, C),
                                                   ESS_threshold = ESS_threshold,
                                                   sub_posterior_means = sapply(input_samples, mean),
                                                   adaptive_mesh = TRUE,
                                                   adaptive_mesh_parameters = list('data_size' = data_sizes[i],
                                                                                   'b' = vanilla_b,
                                                                                   'k3' = k3,
                                                                                   'k4' = k4,
                                                                                   'vanilla' = TRUE),
                                                   diffusion_estimator = diffusion_estimator,
                                                   seed = seed*i)
  print('### SH: performing Bayesian Fusion with a preconditioning matrix (with recommended T, adaptive mesh)')
  gen_guide_SH <- BF_guidance(condition = 'SH',
                              CESS_0_threshold = CESS_0_threshold,
                              C = C,
                              d = 1,
                              data_size = data_sizes[i],
                              sub_posterior_means = sapply(input_samples, mean),
                              precondition_matrices = sapply(input_samples, var),
                              inv_precondition_matrices = 1/sapply(input_samples, var),
                              k1 = k1,
                              k3 = k3,
                              k4 = k4,
                              vanilla = FALSE)
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = FALSE,
                                              number_of_steps = length(gen_guide_SH$mesh))
  SH_adaptive_generalised <- parallel_GBF_uniGaussian(particles_to_fuse = input_particles,
                                                      N = nsamples,
                                                      m = C,
                                                      time_mesh = gen_guide_SH$mesh,
                                                      means = means,
                                                      sds = rep(sd, C),
                                                      betas = rep(beta, C),
                                                      precondition_values = sapply(input_samples, var),
                                                      ESS_threshold = ESS_threshold,
                                                      sub_posterior_means = sapply(input_samples, mean),
                                                      adaptive_mesh = TRUE,
                                                      adaptive_mesh_parameters = list('data_size' = data_sizes[i],
                                                                                      'k3' = k3,
                                                                                      'k4' = k4,
                                                                                      'vanilla' = FALSE),
                                                      diffusion_estimator = diffusion_estimator,
                                                      seed = seed*i)
  # save results
  SH_adaptive_results$vanilla[[i]] <- list('CESS_0' = SH_adaptive_standard$CESS[1],
                                           'CESS_j' = SH_adaptive_standard$CESS[2:length(SH_adaptive_standard$CESS)],
                                           'CESS_j_avg' = mean(SH_adaptive_standard$CESS[2:length(SH_adaptive_standard$CESS)]),
                                           'n' = length(SH_adaptive_standard$CESS),
                                           'time_mesh' = SH_adaptive_standard$particles$time_mesh,
                                           'time' = SH_adaptive_standard$time)
  SH_adaptive_results$generalised[[i]] <- list('CESS_0' = SH_adaptive_generalised$CESS[1],
                                               'CESS_j' = SH_adaptive_generalised$CESS[2:length(SH_adaptive_generalised$CESS)],
                                               'CESS_j_avg' = mean(SH_adaptive_generalised$CESS[2:length(SH_adaptive_generalised$CESS)]),
                                               'n' = length(SH_adaptive_generalised$CESS),
                                               'time_mesh' = SH_adaptive_generalised$particles$time_mesh,
                                               'time' = SH_adaptive_generalised$time)
  # plot
  lines(density(resample_particle_y_samples(particle_set = SH_adaptive_standard$particles, 
                                            multivariate = FALSE,
                                            resampling_method = 'resid',
                                            seed = seed*i)$y_samples), 
        col = 'orange')
  lines(density(resample_particle_y_samples(particle_set = SH_adaptive_generalised$particles, 
                                            multivariate = FALSE,
                                            resampling_method = 'resid',
                                            seed = seed*i)$y_samples), 
        col = 'orange', lty = 2)
}

##### vanilla plots #####

##### Fixed user-specified parameters #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) a_results$vanilla[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) a_results$vanilla[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
axis(1, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('Conditional effective sample size', 2, 2.75, font = 2, cex = 1.5)

##### Recommended scaling of T, fixed n #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) b_results$vanilla[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) b_results$vanilla[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
axis(1, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('Conditional effective sample size', 2, 2.75, font = 2, cex = 1.5)

##### Recommended scaling of T, regular mesh #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) c_results$vanilla[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) c_results$vanilla[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
axis(1, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('Conditional effective sample size', 2, 2.75, font = 2, cex = 1.5)

##### Recommended scaling of T, adaptive mesh #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) d_results$vanilla[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) d_results$vanilla[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
axis(1, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('Conditional effective sample size', 2, 2.75, font = 2, cex = 1.5)

##### SH: Recommended scaling of T, adaptive mesh #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) SH_adaptive_results$vanilla[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) SH_adaptive_results$vanilla[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
axis(1, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('Conditional effective sample size', 2, 2.75, font = 2, cex = 1.5)

##### generalised plots #####

##### Fixed user-specified parameters #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) a_results$generalised[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) a_results$generalised[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
axis(1, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('Conditional effective sample size', 2, 2.75, font = 2, cex = 1.5)

##### Recommended scaling of T, fixed n #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) b_results$generalised[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) b_results$generalised[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
axis(1, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('Conditional effective sample size', 2, 2.75, font = 2, cex = 1.5)

##### Recommended scaling of T, regular mesh #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) c_results$generalised[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) c_results$generalised[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
axis(1, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('Conditional effective sample size', 2, 2.75, font = 2, cex = 1.5)

##### Recommended scaling of T, adaptive mesh #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) d_results$generalised[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) d_results$generalised[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
axis(1, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('Conditional effective sample size', 2, 2.75, font = 2, cex = 1.5)

##### SH: Recommended scaling of T, adaptive mesh #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) SH_adaptive_results$vanilla[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) SH_adaptive_results$vanilla[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
axis(1, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('Conditional effective sample size', 2, 2.75, font = 2, cex = 1.5)

save.image('bf_dissimilar_means_example.RData')
