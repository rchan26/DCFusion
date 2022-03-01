library(DCFusion)

seed <- 1994
set.seed(seed)
nsamples <- 10000
C <- 10
mean <- 0
beta <- 1/C
a_mesh_vanilla <- seq(0, 0.005, length.out = 6)
a_mesh_gen <- seq(0, 1, length.out = 6)
diffusion_estimator <- 'NB'
ESS_threshold <- 0.5
CESS_0_threshold <- 0.2
vanilla_b <- 1
k1 <- NULL
k2 <- NULL
k3 <- 1
k4 <- 1
data_sizes <- c(1000, 5000, 10000, 20000, 30000, 40000, 50000)
a_results <- list('vanilla' = list(), 'generalised' = list())
b_results <- list('vanilla' = list(), 'generalised' = list())
c_results <- list('vanilla' = list(), 'generalised' = list())
d_results <- list('vanilla' = list(), 'generalised' = list())
SSH_adaptive_results <- list('vanilla' = list(), 'generalised' = list())

for (i in 1:length(data_sizes)) {
  set.seed(seed*i)
  sd <- sqrt(1/data_sizes[i])
  opt_bw <- ((4*sd^5)/(3*nsamples))^(1/5)
  curve(dnorm(x, mean = mean, sd = sd), mean-4*sd, mean+4*sd,
        main = paste('m:', data_sizes[i]))
  input_samples <- lapply(1:C, function(i) rnorm_tempered(N = nsamples,
                                                          mean = mean,
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
                                            means = rep(mean, C),
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
                                               means = rep(mean, C),
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
                                 'time' = a_BF_standard$time,
                                 'IAD' = integrated_abs_distance_uniGaussian(
                                   fusion_post = resample_particle_y_samples(particle_set = a_BF_standard$particles, 
                                                                             multivariate = FALSE,
                                                                             resampling_method = 'resid',
                                                                             seed = seed*i)$y_samples,
                                   mean = 0,
                                   sd = sd,
                                   beta = 1,
                                   bw = opt_bw))
  a_results$generalised[[i]] <- list('CESS_0' = a_BF_generalised$CESS[1],
                                     'CESS_j' = a_BF_generalised$CESS[2:length(a_BF_generalised$CESS)],
                                     'CESS_j_avg' = mean(a_BF_generalised$CESS[2:length(a_BF_generalised$CESS)]),
                                     'n' = length(a_BF_generalised$CESS),
                                     'time_mesh' = a_BF_generalised$particles$time_mesh,
                                     'time' = a_BF_generalised$time,
                                     'IAD' = integrated_abs_distance_uniGaussian(
                                       fusion_post = resample_particle_y_samples(particle_set = a_BF_generalised$particles, 
                                                                                 multivariate = FALSE,
                                                                                 resampling_method = 'resid',
                                                                                 seed = seed*i)$y_samples,
                                       mean = 0,
                                       sd = sd,
                                       beta = 1,
                                       bw = opt_bw))
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
  vanilla_guide <- BF_guidance(condition = 'SH',
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
  b_mesh_vanilla <- seq(0, vanilla_guide$min_T, length.out = 6)
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = FALSE,
                                              number_of_steps = length(b_mesh_vanilla))
  b_BF_standard <- parallel_GBF_uniGaussian(particles_to_fuse = input_particles,
                                            N = nsamples,
                                            m = C,
                                            time_mesh = b_mesh_vanilla,
                                            means = rep(mean, C),
                                            sds = rep(sd, C),
                                            betas = rep(beta, C),
                                            precondition_values = rep(1, C),
                                            ESS_threshold = ESS_threshold,
                                            diffusion_estimator = diffusion_estimator,
                                            seed = seed*i)
  print('### performing Bayesian Fusion with a preconditioning matrix (with recommended T, fixed n)')
  gen_guide <- BF_guidance(condition = 'SH',
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
  b_mesh_gen <- seq(0, gen_guide$min_T, length.out = 6)
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = FALSE,
                                              number_of_steps = length(b_mesh_gen))
  b_BF_generalised <- parallel_GBF_uniGaussian(particles_to_fuse = input_particles,
                                               N = nsamples,
                                               m = C,
                                               time_mesh = b_mesh_gen,
                                               means = rep(mean, C),
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
                                 'time' = b_BF_standard$time,
                                 'IAD' = integrated_abs_distance_uniGaussian(
                                   fusion_post = resample_particle_y_samples(particle_set = b_BF_standard$particles, 
                                                                             multivariate = FALSE,
                                                                             resampling_method = 'resid',
                                                                             seed = seed*i)$y_samples,
                                   mean = 0,
                                   sd = sd,
                                   beta = 1,
                                   bw = opt_bw))
  b_results$generalised[[i]] <- list('CESS_0' = b_BF_generalised$CESS[1],
                                     'CESS_j' = b_BF_generalised$CESS[2:length(b_BF_generalised$CESS)],
                                     'CESS_j_avg' = mean(b_BF_generalised$CESS[2:length(b_BF_generalised$CESS)]),
                                     'n' = length(b_BF_generalised$CESS),
                                     'time_mesh' = b_BF_generalised$particles$time_mesh,
                                     'time' = b_BF_generalised$time,
                                     'IAD' = integrated_abs_distance_uniGaussian(
                                       fusion_post = resample_particle_y_samples(particle_set = b_BF_generalised$particles, 
                                                                                 multivariate = FALSE,
                                                                                 resampling_method = 'resid',
                                                                                 seed = seed*i)$y_samples,
                                       mean = 0,
                                       sd = sd,
                                       beta = 1,
                                       bw = opt_bw))
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
                                            means = rep(mean, C),
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
                                               means = rep(mean, C),
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
                                 'time' = c_BF_standard$time,
                                 'IAD' = integrated_abs_distance_uniGaussian(
                                   fusion_post = resample_particle_y_samples(particle_set = c_BF_standard$particles, 
                                                                             multivariate = FALSE,
                                                                             resampling_method = 'resid',
                                                                             seed = seed*i)$y_samples,
                                   mean = 0,
                                   sd = sd,
                                   beta = 1,
                                   bw = opt_bw))
  c_results$generalised[[i]] <- list('CESS_0' = c_BF_generalised$CESS[1],
                                     'CESS_j' = c_BF_generalised$CESS[2:length(c_BF_generalised$CESS)],
                                     'CESS_j_avg' = mean(c_BF_generalised$CESS[2:length(c_BF_generalised$CESS)]),
                                     'n' = length(c_BF_generalised$CESS),
                                     'time_mesh' = c_BF_generalised$particles$time_mesh,
                                     'time' = c_BF_generalised$time,
                                     'IAD' = integrated_abs_distance_uniGaussian(
                                       fusion_post = resample_particle_y_samples(particle_set = c_BF_generalised$particles, 
                                                                                 multivariate = FALSE,
                                                                                 resampling_method = 'resid',
                                                                                 seed = seed*i)$y_samples,
                                       mean = 0,
                                       sd = sd,
                                       beta = 1,
                                       bw = opt_bw))
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
                                            means = rep(mean, C),
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
                                               means = rep(mean, C),
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
                                 'time' = d_BF_standard$time,
                                 'IAD' = integrated_abs_distance_uniGaussian(
                                   fusion_post = resample_particle_y_samples(particle_set = d_BF_standard$particles, 
                                                                             multivariate = FALSE,
                                                                             resampling_method = 'resid',
                                                                             seed = seed*i)$y_samples,
                                   mean = 0,
                                   sd = sd,
                                   beta = 1,
                                   bw = opt_bw))
  d_results$generalised[[i]] <- list('CESS_0' = d_BF_generalised$CESS[1],
                                     'CESS_j' = d_BF_generalised$CESS[2:length(d_BF_generalised$CESS)],
                                     'CESS_j_avg' = mean(d_BF_generalised$CESS[2:length(d_BF_generalised$CESS)]),
                                     'n' = length(d_BF_generalised$CESS),
                                     'time_mesh' = d_BF_generalised$particles$time_mesh,
                                     'time' = d_BF_generalised$time,
                                     'IAD' = integrated_abs_distance_uniGaussian(
                                       fusion_post = resample_particle_y_samples(particle_set = d_BF_generalised$particles, 
                                                                                 multivariate = FALSE,
                                                                                 resampling_method = 'resid',
                                                                                 seed = seed*i)$y_samples,
                                       mean = 0,
                                       sd = sd,
                                       beta = 1,
                                       bw = opt_bw))
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
  
  ##### SSH: Recommended scaling of T, adaptive mesh #####
  print('### SSH: performing standard Bayesian Fusion (with recommended T, adaptive mesh)')
  vanilla_guide_SSH <- BF_guidance(condition = 'SSH',
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
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = FALSE,
                                              number_of_steps = length(vanilla_guide_SSH$mesh))
  SSH_adaptive_standard <- parallel_GBF_uniGaussian(particles_to_fuse = input_particles,
                                                    N = nsamples,
                                                    m = C,
                                                    time_mesh = vanilla_guide_SSH$mesh,
                                                    means = rep(mean, C),
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
  print('### SSH: performing Bayesian Fusion with a preconditioning matrix (with recommended T, adaptive mesh)')
  gen_guide_SSH <- BF_guidance(condition = 'SSH',
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
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = FALSE,
                                              number_of_steps = length(gen_guide_SSH$mesh))
  SSH_adaptive_generalised <- parallel_GBF_uniGaussian(particles_to_fuse = input_particles,
                                                       N = nsamples,
                                                       m = C,
                                                       time_mesh = gen_guide_SSH$mesh,
                                                       means = rep(mean, C),
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
  SSH_adaptive_results$vanilla[[i]] <- list('CESS_0' = SSH_adaptive_standard$CESS[1],
                                            'CESS_j' = SSH_adaptive_standard$CESS[2:length(SSH_adaptive_standard$CESS)],
                                            'CESS_j_avg' = mean(SSH_adaptive_standard$CESS[2:length(SSH_adaptive_standard$CESS)]),
                                            'n' = length(SSH_adaptive_standard$CESS),
                                            'time_mesh' = SSH_adaptive_standard$particles$time_mesh,
                                            'time' = SSH_adaptive_standard$time,
                                            'IAD' = integrated_abs_distance_uniGaussian(
                                              fusion_post = resample_particle_y_samples(particle_set = SSH_adaptive_standard$particles, 
                                                                                        multivariate = FALSE,
                                                                                        resampling_method = 'resid',
                                                                                        seed = seed*i)$y_samples,
                                              mean = 0,
                                              sd = sd,
                                              beta = 1,
                                              bw = opt_bw))
  SSH_adaptive_results$generalised[[i]] <- list('CESS_0' = SSH_adaptive_generalised$CESS[1],
                                                'CESS_j' = SSH_adaptive_generalised$CESS[2:length(SSH_adaptive_generalised$CESS)],
                                                'CESS_j_avg' = mean(SSH_adaptive_generalised$CESS[2:length(SSH_adaptive_generalised$CESS)]),
                                                'n' = length(SSH_adaptive_generalised$CESS),
                                                'time_mesh' = SSH_adaptive_generalised$particles$time_mesh,
                                                'time' = SSH_adaptive_generalised$time,
                                                'IAD' = integrated_abs_distance_uniGaussian(
                                                  fusion_post = resample_particle_y_samples(particle_set = SSH_adaptive_generalised$particles, 
                                                                                            multivariate = FALSE,
                                                                                            resampling_method = 'resid',
                                                                                            seed = seed*i)$y_samples,
                                                  mean = 0,
                                                  sd = sd,
                                                  beta = 1,
                                                  bw = opt_bw))
  # plot
  lines(density(resample_particle_y_samples(particle_set = SSH_adaptive_standard$particles, 
                                            multivariate = FALSE,
                                            resampling_method = 'resid',
                                            seed = seed*i)$y_samples), 
        col = 'orange')
  lines(density(resample_particle_y_samples(particle_set = SSH_adaptive_generalised$particles, 
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
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) a_results$vanilla[[i]]$CESS_j/nsamples)[[ii]]
  points(x = rep(data_sizes[ii], length(cess_j)), y = cess_j, cex = 0.5)
}
axis(1, at = c(1000, seq(10000, 50000, 10000)),
     labels = c(1000, seq(10000, 50000, 10000)), font = 2, cex = 1.5)
axis(1, at = seq(0, 50000, 5000), labels = rep("", 11), lwd.ticks = 0.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('CESS / N', 2, 2.75, font = 2, cex = 1.5)

##### Recommended scaling of T, fixed n #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) b_results$vanilla[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) b_results$vanilla[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) b_results$vanilla[[i]]$CESS_j/nsamples)[[ii]]
  points(x = rep(data_sizes[ii], length(cess_j)), y = cess_j, cex = 0.5)
}
axis(1, at = c(1000, seq(10000, 50000, 10000)),
     labels = c(1000, seq(10000, 50000, 10000)), font = 2, cex = 1.5)
axis(1, at = seq(0, 50000, 5000), labels = rep("", 11), lwd.ticks = 0.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('CESS / N', 2, 2.75, font = 2, cex = 1.5)

##### Recommended scaling of T, regular mesh #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) c_results$vanilla[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) c_results$vanilla[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) c_results$vanilla[[i]]$CESS_j/nsamples)[[ii]]
  points(x = rep(data_sizes[ii], length(cess_j)), y = cess_j, cex = 0.5)
}
axis(1, at = c(1000, seq(10000, 50000, 10000)),
     labels = c(1000, seq(10000, 50000, 10000)), font = 2, cex = 1.5)
axis(1, at = seq(0, 50000, 5000), labels = rep("", 11), lwd.ticks = 0.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('CESS / N', 2, 2.75, font = 2, cex = 1.5)

##### Recommended scaling of T, adaptive mesh #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) d_results$vanilla[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) d_results$vanilla[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) d_results$vanilla[[i]]$CESS_j/nsamples)[[ii]]
  points(x = rep(data_sizes[ii], length(cess_j)), y = cess_j, cex = 0.5)
}
axis(1, at = c(1000, seq(10000, 50000, 10000)),
     labels = c(1000, seq(10000, 50000, 10000)), font = 2, cex = 1.5)
axis(1, at = seq(0, 50000, 5000), labels = rep("", 11), lwd.ticks = 0.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('CESS / N', 2, 2.75, font = 2, cex = 1.5)

##### Compare regular mesh and adaptive mesh times #####
plot(x = c_results$vanilla[[1]]$time_mesh,
     y = rep(data_sizes[1]-500, length(c_results$vanilla[[1]]$time_mesh)),
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,50000),
     xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
for (i in  2:length(data_sizes)) {
  lines(x = c_results$vanilla[[i]]$time_mesh,
        y = rep(data_sizes[i]-500, length(c_results$vanilla[[i]]$time_mesh)),
        type = 'b', pch = 20, lty = 1, lwd = 3)
}
for (i in  1:length(data_sizes)) {
  lines(x = d_results$vanilla[[i]]$time_mesh,
        y = rep(data_sizes[i]+500, length(d_results$vanilla[[i]]$time_mesh)),
        type = 'b', pch = 4, lty = 2, lwd = 3)
}
max(sapply(1:length(data_sizes), function(i) c_results$vanilla[[i]]$time_mesh))
axis(1, at = seq(0, 0.03, 0.01), labels = seq(0, 0.03, 0.01), font = 2, cex = 1.5)
axis(1, at = seq(0, 0.03, 0.005), labels = rep("", 7), lwd.ticks = 0.5)
mtext('Time', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = c(1000, seq(10000, 50000, 10000)),
     labels = c(1000, seq(10000, 50000, 10000)), font = 2, cex = 1.5)
axis(2, at = seq(0, 50000, 5000), labels = rep("", 11), lwd.ticks = 0.5)
mtext('Data Sizes', 2, 2.75, font = 2, cex = 1.5)

# plot(x = d_results$vanilla[[1]]$time_mesh,
#      y = rep(data_sizes[1], length(d_results$vanilla[[1]]$time_mesh)),
#      type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,50000),
#      xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
# for (i in  1:length(data_sizes)) {
#   lines(x = d_results$vanilla[[i]]$time_mesh,
#         y = rep(data_sizes[i], length(d_results$vanilla[[i]]$time_mesh)),
#         type = 'b', pch = 20, lty = 1, lwd = 3)
# }
# # max(sapply(1:length(data_sizes), function(i) max(d_results$vanilla[[i]]$time_mesh)))
# axis(1, at = seq(0, 0.03, 0.01), labels = seq(0, 0.03, 0.01), font = 2, cex = 1.5)
# axis(1, at = seq(0, 0.03, 0.005), labels = rep("", 7), lwd.ticks = 0.5)
# mtext('Time', 1, 2.75, font = 2, cex = 1.5)
# axis(2, at = c(1000, seq(10000, 50000, 10000)),      labels = c(1000, seq(10000, 50000, 10000)), font = 2, cex = 1.5)
# axis(2, at = seq(0, 50000, 5000), labels = rep("", 11), lwd.ticks = 0.5)
# mtext('Data Sizes', 2, 2.75, font = 2, cex = 1.5)

##### SSH: Recommended scaling of T, adaptive mesh #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) SSH_adaptive_results$vanilla[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) SSH_adaptive_results$vanilla[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) SSH_adaptive_results$vanilla[[i]]$CESS_j/nsamples)[[ii]]
  points(x = rep(data_sizes[ii], length(cess_j)), y = cess_j, cex = 0.5)
}
axis(1, at = c(1000, seq(10000, 50000, 10000)),
     labels = c(1000, seq(10000, 50000, 10000)), font = 2, cex = 1.5)
axis(1, at = seq(0, 50000, 5000), labels = rep("", 11), lwd.ticks = 0.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('CESS / N', 2, 2.75, font = 2, cex = 1.5)

##### IAD #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) a_results$vanilla[[i]]$IAD),
     type = 'b', pch = 1, lty = 1, lwd = 3, ylim = c(0,0.6), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) b_results$vanilla[[i]]$IAD),
      pch = 2, lty = 2, lwd = 3, type = 'b')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) c_results$vanilla[[i]]$IAD),
      pch = 3, lty = 3, lwd = 3, type = 'b')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) d_results$vanilla[[i]]$IAD),
      pch = 4, lty = 4, lwd = 3, type = 'b')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) SSH_adaptive_results$vanilla[[i]]$IAD),
      pch = 5, lty = 5, lwd = 3, type = 'b')
axis(1, at = c(1000, seq(10000, 50000, 10000)),
     labels = c(1000, seq(10000, 50000, 10000)), font = 2, cex = 1.5)
axis(1, at = seq(0, 50000, 5000), labels = rep("", 11), lwd.ticks = 0.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = c("0.0", seq(0.1, 0.9, 0.1), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
legend(x = 250, y = 0.6,
       legend = c('Fixed T, fixed n',
                  'SH rec. T, fixed n',
                  'SH rec. T, reg. mesh',
                  'SH rec. T, adapt. mesh',
                  'SSH rec. T, adapt. mesh'),
       lty = 1:5,
       pch = 1:5,
       lwd = rep(3, 5),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

##### time #####
plot(x = data_sizes,
     y = log(sapply(1:length(data_sizes), function(i) a_results$vanilla[[i]]$time)),
     type = 'b', pch = 1, lty = 1, lwd = 3, ylim = c(1,5), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = log(sapply(1:length(data_sizes), function(i) b_results$vanilla[[i]]$time)),
      pch = 2, lty = 2, lwd = 3, type = 'b')
lines(x = data_sizes,
      y = log(sapply(1:length(data_sizes), function(i) c_results$vanilla[[i]]$time)),
      pch = 3, lty = 3, lwd = 3, type = 'b')
lines(x = data_sizes,
      y = log(sapply(1:length(data_sizes), function(i) d_results$vanilla[[i]]$time)),
      pch = 4, lty = 4, lwd = 3, type = 'b')
lines(x = data_sizes,
      y = log(sapply(1:length(data_sizes), function(i) SSH_adaptive_results$vanilla[[i]]$time)),
      pch = 5, lty = 5, lwd = 3, type = 'b')
axis(1, at = c(1000, seq(10000, 50000, 10000)),
     labels = c(1000, seq(10000, 50000, 10000)), font = 2, cex = 1.5)
axis(1, at = seq(0, 50000, 5000), labels = rep("", 11), lwd.ticks = 0.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 5, 1), labels = seq(0, 5, 1), font = 2, cex = 1.5)
mtext('log(Elapsed time in seconds)', 2, 2.75, font = 2, cex = 1.5)
legend(x = 0, y = 5,
       legend = c('Fixed T, fixed n',
                  'SH rec. T, fixed n',
                  'SH rec. T, reg. mesh',
                  'SH rec. T, adapt. mesh',
                  'SSH rec. T, adapt. mesh'),
       lty = 1:5,
       pch = 1:5,
       lwd = rep(3, 5),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

##### generalised plots #####

##### Fixed user-specified parameters #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) a_results$generalised[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) a_results$generalised[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) a_results$generalised[[i]]$CESS_j/nsamples)[[ii]]
  points(x = rep(data_sizes[ii], length(cess_j)), y = cess_j, cex = 0.5)
}
axis(1, at = c(1000, seq(10000, 50000, 10000)),
     labels = c(1000, seq(10000, 50000, 10000)), font = 2, cex = 1.5)
axis(1, at = seq(0, 50000, 5000), labels = rep("", 11), lwd.ticks = 0.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('CESS / N', 2, 2.75, font = 2, cex = 1.5)

##### Recommended scaling of T, fixed n #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) b_results$generalised[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) b_results$generalised[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) b_results$generalised[[i]]$CESS_j/nsamples)[[ii]]
  points(x = rep(data_sizes[ii], length(cess_j)), y = cess_j, cex = 0.5)
}
axis(1, at = c(1000, seq(10000, 50000, 10000)),
     labels = c(1000, seq(10000, 50000, 10000)), font = 2, cex = 1.5)
axis(1, at = seq(0, 50000, 5000), labels = rep("", 11), lwd.ticks = 0.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('CESS / N', 2, 2.75, font = 2, cex = 1.5)

##### Recommended scaling of T, regular mesh #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) c_results$generalised[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) c_results$generalised[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) c_results$generalised[[i]]$CESS_j/nsamples)[[ii]]
  points(x = rep(data_sizes[ii], length(cess_j)), y = cess_j, cex = 0.5)
}
axis(1, at = c(1000, seq(10000, 50000, 10000)),
     labels = c(1000, seq(10000, 50000, 10000)), font = 2, cex = 1.5)
axis(1, at = seq(0, 50000, 5000), labels = rep("", 11), lwd.ticks = 0.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('CESS / N', 2, 2.75, font = 2, cex = 1.5)

##### Recommended scaling of T, adaptive mesh #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) d_results$generalised[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) d_results$generalised[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) d_results$generalised[[i]]$CESS_j/nsamples)[[ii]]
  points(x = rep(data_sizes[ii], length(cess_j)), y = cess_j, cex = 0.5)
}
axis(1, at = c(1000, seq(10000, 50000, 10000)),
     labels = c(1000, seq(10000, 50000, 10000)), font = 2, cex = 1.5)
axis(1, at = seq(0, 50000, 5000), labels = rep("", 11), lwd.ticks = 0.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('CESS / N', 2, 2.75, font = 2, cex = 1.5)

##### Compare regular mesh and adaptive mesh times #####
plot(x = c_results$generalised[[1]]$time_mesh,
     y = rep(data_sizes[1]-500, length(c_results$generalised[[1]]$time_mesh)),
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,50000),
     xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
for (i in  2:length(data_sizes)) {
  lines(x = c_results$generalised[[i]]$time_mesh,
        y = rep(data_sizes[i]-500, length(c_results$generalised[[i]]$time_mesh)),
        type = 'b', pch = 20, lty = 1, lwd = 3)
}
for (i in  1:length(data_sizes)) {
  lines(x = d_results$generalised[[i]]$time_mesh,
        y = rep(data_sizes[i]+500, length(d_results$generalised[[i]]$time_mesh)),
        type = 'b', pch = 4, lty = 2, lwd = 3)
}
max(sapply(1:length(data_sizes), function(i) max(c_results$generalised[[i]]$time_mesh)))
axis(1, at = seq(0, 3, 1), labels = seq(0, 3, 1), font = 2, cex = 1.5)
axis(1, at = seq(0, 3, 0.5), labels = rep("", 7), lwd.ticks = 0.5)
mtext('Time', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = c(1000, seq(10000, 50000, 10000)),
     labels = c(1000, seq(10000, 50000, 10000)), font = 2, cex = 1.5)
axis(2, at = seq(0, 50000, 5000), labels = rep("", 11), lwd.ticks = 0.5)
mtext('Data Sizes', 2, 2.75, font = 2, cex = 1.5)

##### SSH: Recommended scaling of T, adaptive mesh #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) SSH_adaptive_results$vanilla[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) SSH_adaptive_results$vanilla[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) SSH_adaptive_results$generalised[[i]]$CESS_j/nsamples)[[ii]]
  points(x = rep(data_sizes[ii], length(cess_j)), y = cess_j, cex = 0.5)
}
axis(1, at = c(1000, seq(10000, 50000, 10000)),
     labels = c(1000, seq(10000, 50000, 10000)), font = 2, cex = 1.5)
axis(1, at = seq(0, 50000, 5000), labels = rep("", 11), lwd.ticks = 0.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('CESS / N', 2, 2.75, font = 2, cex = 1.5)

##### IAD #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) a_results$generalised[[i]]$IAD),
     type = 'b', pch = 1, lty = 1, lwd = 3, ylim = c(0,0.6), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) b_results$generalised[[i]]$IAD),
      pch = 2, lty = 2, lwd = 3, type = 'b')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) c_results$generalised[[i]]$IAD),
      pch = 3, lty = 3, lwd = 3, type = 'b')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) d_results$generalised[[i]]$IAD),
      pch = 4, lty = 4, lwd = 3, type = 'b')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) SSH_adaptive_results$generalised[[i]]$IAD),
      pch = 5, lty = 5, lwd = 3, type = 'b')
axis(1, at = c(1000, seq(10000, 50000, 10000)),
     labels = c(1000, seq(10000, 50000, 10000)), font = 2, cex = 1.5)
axis(1, at = seq(0, 50000, 5000), labels = rep("", 11), lwd.ticks = 0.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = c("0.0", seq(0.1, 0.9, 0.1), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
legend(x = 0, y = 0.6,
       legend = c('Fixed T, fixed n',
                  'SH rec. T, fixed n',
                  'SH rec. T, reg. mesh',
                  'SH rec. T, adapt. mesh',
                  'SSH rec. T, adapt. mesh'),
       lty = 1:5,
       pch = 1:5,
       lwd = rep(3, 5),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

##### time #####
plot(x = data_sizes,
     y = log(sapply(1:length(data_sizes), function(i) a_results$generalised[[i]]$time)),
     type = 'b', pch = 1, lty = 1, lwd = 3, ylim = c(1,5), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = log(sapply(1:length(data_sizes), function(i) b_results$generalised[[i]]$time)),
      pch = 2, lty = 2, lwd = 3, type = 'b')
lines(x = data_sizes,
      y = log(sapply(1:length(data_sizes), function(i) c_results$generalised[[i]]$time)),
      pch = 3, lty = 3, lwd = 3, type = 'b')
lines(x = data_sizes,
      y = log(sapply(1:length(data_sizes), function(i) d_results$generalised[[i]]$time)),
      pch = 4, lty = 4, lwd = 3, type = 'b')
lines(x = data_sizes,
      y = log(sapply(1:length(data_sizes), function(i) SSH_adaptive_results$generalised[[i]]$time)),
      pch = 5, lty = 5, lwd = 3, type = 'b')
axis(1, at = c(1000, seq(10000, 50000, 10000)),
     labels = c(1000, seq(10000, 50000, 10000)), font = 2, cex = 1.5)
axis(1, at = seq(0, 50000, 5000), labels = rep("", 11), lwd.ticks = 0.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 5, 1), labels = seq(0, 5, 1), font = 2, cex = 1.5)
mtext('log(Elapsed time in seconds)', 2, 2.75, font = 2, cex = 1.5)
legend(x = 0, y = 5,
       legend = c('Fixed T, fixed n',
                  'SH rec. T, fixed n',
                  'SH rec. T, reg. mesh',
                  'SH rec. T, adapt. mesh',
                  'SSH rec. T, adapt. mesh'),
       lty = 1:5,
       pch = 1:5,
       lwd = rep(3, 5),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

save.image('bf_similar_means_example.RData')
