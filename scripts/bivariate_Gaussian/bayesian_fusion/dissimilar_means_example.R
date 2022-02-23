library(DCFusion)

seed <- 1994
set.seed(seed)
nsamples <- 10000
C <- 2
means <- list(rep(-0.25, 2), rep(0.25, 2))
corr <- 0.9
beta <- 1
a_mesh_vanilla <- seq(0, 0.01, length.out = 6)
a_mesh_gen <- seq(0, 2, length.out = 6)
diffusion_estimator <- 'NB'
ESS_threshold <- 0.5
CESS_0_threshold <- 0.2
CESS_j_threshold <- 0.2
vanilla_b <- 1
k1 <- NULL
k2 <- NULL
k3 <- -log(CESS_j_threshold)/2
k4 <- -log(CESS_j_threshold)/2
data_sizes <- c(250, 500, 1000, 1500, 2000, 2500)
a_results <- list('vanilla' = list(), 'generalised' = list())
b_results <- list('vanilla' = list(), 'generalised' = list())
c_results <- list('vanilla' = list(), 'generalised' = list())
d_results <- list('vanilla' = list(), 'generalised' = list())
SH_adaptive_results <- list('vanilla' = list(), 'generalised' = list())

for (i in 1:length(data_sizes)) {
  set.seed(seed*i)
  sd <- sqrt(rep(C, 2)/data_sizes[i])
  # target_cov_mat <- matrix(c(sd[1]^2, sd[1]*sd[2]*corr, sd[1]*sd[2]*corr, sd[2]^2),
  #                          nrow = 2, ncol = 2, byrow = T)/C
  # target_samples <- mvrnormArma(N = nsamples, mu = mean, Sigma = target_cov_mat)
  cov_mat <- matrix(c(sd[1]^2, sd[1]*sd[2]*corr, sd[1]*sd[2]*corr, sd[2]^2),
                    nrow = 2, ncol = 2, byrow = T)
  opt_bw <- ((4*sd^5)/(3*nsamples))^(1/5)
  input_samples <- lapply(1:C, function(sub) mvrnormArma(N = nsamples, mu = means[[sub]], Sigma = cov_mat))
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples, multivariate = TRUE)
  ##### Fixed user-specified parameters #####
  print('### performing standard Bayesian Fusion (with standard mesh)')
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = TRUE,
                                              number_of_steps = length(a_mesh_vanilla))
  a_BF_standard <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                           N = nsamples,
                                           m = C,
                                           time_mesh = a_mesh_vanilla,
                                           mean_vecs = means,
                                           sd_vecs = rep(list(sd), C),
                                           corrs = rep(corr, C),
                                           betas = rep(beta, C),
                                           precondition_matrices = rep(list(diag(1,2)), C),
                                           ESS_threshold = ESS_threshold,
                                           diffusion_estimator = diffusion_estimator,
                                           seed = seed*i)
  print('### performing Bayesian Fusion with a preconditioning matrix (with standard mesh)')
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = TRUE,
                                              number_of_steps = length(a_mesh_gen))
  a_BF_generalised <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                              N = nsamples,
                                              m = C,
                                              time_mesh = a_mesh_gen,
                                              mean_vecs = means,
                                              sd_vecs = rep(list(sd), C),
                                              corrs = rep(corr, C),
                                              betas = rep(beta, C),
                                              precondition_matrices = lapply(input_samples, cov),
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
                                 'elapsed_time' = a_BF_standard$elapsed_time,
                                 'IAD' = integrated_abs_distance_biGaussian(fusion_post = resample_particle_y_samples(
                                   particle_set = a_BF_standard$particles,
                                   multivariate = TRUE,
                                   resampling_method = 'resid',
                                   seed = seed*i)$y_samples,
                                   marg_means = c(0,0),
                                   marg_sds = sqrt(rep(1, 2)/data_sizes[i]),
                                   bw = opt_bw))
  a_results$generalised[[i]] <- list('CESS_0' = a_BF_generalised$CESS[1],
                                     'CESS_j' = a_BF_generalised$CESS[2:length(a_BF_generalised$CESS)],
                                     'CESS_j_avg' = mean(a_BF_generalised$CESS[2:length(a_BF_generalised$CESS)]),
                                     'n' = length(a_BF_generalised$CESS),
                                     'time_mesh' = a_BF_generalised$particles$time_mesh,
                                     'time' = a_BF_generalised$time,
                                     'elapsed_time' = a_BF_generalised$elapsed_time,
                                     'IAD' = integrated_abs_distance_biGaussian(fusion_post = resample_particle_y_samples(
                                       particle_set = a_BF_generalised$particles,
                                       multivariate = TRUE,
                                       resampling_method = 'resid',
                                       seed = seed*i)$y_samples,
                                       marg_means = c(0,0),
                                       marg_sds = sqrt(rep(1, 2)/data_sizes[i]),
                                       bw = opt_bw))
  
  ##### Recommended scaling of T, fixed n #####
  print('### performing standard Bayesian Fusion (with recommended T, fixed n)')
  vanilla_guide <- BF_guidance(condition = 'SSH',
                               CESS_0_threshold = CESS_0_threshold,
                               CESS_j_threshold = CESS_j_threshold,
                               sub_posterior_samples = input_samples,
                               C = C,
                               d = 2,
                               data_size = data_sizes[i],
                               b = vanilla_b,
                               sub_posterior_means = t(sapply(input_samples, function(sub) apply(sub, 2, mean))),
                               k1 = k1,
                               k2 = k2,
                               vanilla = TRUE)
  print(paste('vanilla recommened regular mesh n:', vanilla_guide$n))
  b_mesh_vanilla <- seq(0, vanilla_guide$min_T, length.out = 6)
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = TRUE,
                                              number_of_steps = length(b_mesh_vanilla))
  b_BF_standard <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                           N = nsamples,
                                           m = C,
                                           time_mesh = b_mesh_vanilla,
                                           mean_vecs = means,
                                           sd_vecs = rep(list(sd), C),
                                           corrs = rep(corr, C),
                                           betas = rep(beta, C),
                                           precondition_matrices = rep(list(diag(1,2)), C),
                                           ESS_threshold = ESS_threshold,
                                           diffusion_estimator = diffusion_estimator,
                                           seed = seed*i)
  print('### performing Bayesian Fusion with a preconditioning matrix (with recommended T, fixed n)')
  gen_guide <- BF_guidance(condition = 'SSH',
                           CESS_0_threshold = CESS_0_threshold,
                           CESS_j_threshold = CESS_j_threshold,
                           sub_posterior_samples = input_samples,
                           C = C,
                           d = 2,
                           data_size = data_sizes[i],
                           sub_posterior_means = t(sapply(input_samples, function(sub) apply(sub, 2, mean))),
                           precondition_matrices = lapply(input_samples, cov),
                           inv_precondition_matrices = lapply(input_samples, function(sub) solve(cov(sub))),
                           Lambda = inverse_sum_matrices(lapply(input_samples, function(sub) solve(cov(sub)))),
                           k1 = k1,
                           k2 = k2,
                           vanilla = FALSE)
  print(paste('generalised recommened regular mesh n:', gen_guide$n))
  b_mesh_gen <- seq(0, gen_guide$min_T, length.out = 6)
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = TRUE,
                                              number_of_steps = length(b_mesh_gen))
  b_BF_generalised <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                              N = nsamples,
                                              m = C,
                                              time_mesh = b_mesh_gen,
                                              mean_vecs = means,
                                              sd_vecs = rep(list(sd), C),
                                              corrs = rep(corr, C),
                                              betas = rep(beta, C),
                                              precondition_matrices = lapply(input_samples, cov),
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
                                 'elapsed_time' = b_BF_standard$elapsed_time,
                                 'IAD' = integrated_abs_distance_biGaussian(fusion_post = resample_particle_y_samples(
                                   particle_set = b_BF_standard$particles,
                                   multivariate = TRUE,
                                   resampling_method = 'resid',
                                   seed = seed*i)$y_samples,
                                   marg_means = c(0,0),
                                   marg_sds = sqrt(rep(1, 2)/data_sizes[i]),
                                   bw = opt_bw))
  b_results$generalised[[i]] <- list('CESS_0' = b_BF_generalised$CESS[1],
                                     'CESS_j' = b_BF_generalised$CESS[2:length(b_BF_generalised$CESS)],
                                     'CESS_j_avg' = mean(b_BF_generalised$CESS[2:length(b_BF_generalised$CESS)]),
                                     'n' = length(b_BF_generalised$CESS),
                                     'time_mesh' = b_BF_generalised$particles$time_mesh,
                                     'time' = b_BF_generalised$time,
                                     'elapsed_time' = b_BF_generalised$elapsed_time,
                                     'IAD' = integrated_abs_distance_biGaussian(fusion_post = resample_particle_y_samples(
                                       particle_set = b_BF_generalised$particles,
                                       multivariate = TRUE,
                                       resampling_method = 'resid',
                                       seed = seed*i)$y_samples,
                                       marg_means = c(0,0),
                                       marg_sds = sqrt(rep(1, 2)/data_sizes[i]),
                                       bw = opt_bw))
  
  ##### Recommended scaling of T, regular mesh #####
  print('### performing standard Bayesian Fusion (with recommended T, regular mesh)')
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = TRUE,
                                              number_of_steps = length(vanilla_guide$mesh))
  c_BF_standard <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                           N = nsamples,
                                           m = C,
                                           time_mesh = vanilla_guide$mesh,
                                           mean_vecs = means,
                                           sd_vecs = rep(list(sd), C),
                                           corrs = rep(corr, C),
                                           betas = rep(beta, C),
                                           precondition_matrices = rep(list(diag(1,2)), C),
                                           ESS_threshold = ESS_threshold,
                                           diffusion_estimator = diffusion_estimator,
                                           seed = seed*i)
  print('### performing Bayesian Fusion with a preconditioning matrix (with recommended T, regular mesh)')
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = TRUE,
                                              number_of_steps = length(gen_guide$mesh))
  c_BF_generalised <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                              N = nsamples,
                                              m = C,
                                              time_mesh = gen_guide$mesh,
                                              mean_vecs = means,
                                              sd_vecs = rep(list(sd), C),
                                              corrs = rep(corr, C),
                                              betas = rep(beta, C),
                                              precondition_matrices = lapply(input_samples, cov),
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
                                 'elapsed_time' = c_BF_standard$elapsed_time,
                                 'IAD' = integrated_abs_distance_biGaussian(fusion_post = resample_particle_y_samples(
                                   particle_set = c_BF_standard$particles,
                                   multivariate = TRUE,
                                   resampling_method = 'resid',
                                   seed = seed*i)$y_samples,
                                   marg_means = c(0,0),
                                   marg_sds = sqrt(rep(1, 2)/data_sizes[i]),
                                   bw = opt_bw))
  c_results$generalised[[i]] <- list('CESS_0' = c_BF_generalised$CESS[1],
                                     'CESS_j' = c_BF_generalised$CESS[2:length(c_BF_generalised$CESS)],
                                     'CESS_j_avg' = mean(c_BF_generalised$CESS[2:length(c_BF_generalised$CESS)]),
                                     'n' = length(c_BF_generalised$CESS),
                                     'time_mesh' = c_BF_generalised$particles$time_mesh,
                                     'time' = c_BF_generalised$time,
                                     'elapsed_time' = c_BF_generalised$elapsed_time,
                                     'IAD' = integrated_abs_distance_biGaussian(fusion_post = resample_particle_y_samples(
                                       particle_set = c_BF_generalised$particles,
                                       multivariate = TRUE,
                                       resampling_method = 'resid',
                                       seed = seed*i)$y_samples,
                                       marg_means = c(0,0),
                                       marg_sds = sqrt(rep(1, 2)/data_sizes[i]),
                                       bw = opt_bw))
  
  ##### Recommended scaling of T, adaptive mesh #####
  print('### performing standard Bayesian Fusion (with recommended T, adaptive mesh)')
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = TRUE,
                                              number_of_steps = length(vanilla_guide$mesh))
  d_BF_standard <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                           N = nsamples,
                                           m = C,
                                           time_mesh = vanilla_guide$mesh,
                                           mean_vecs = means,
                                           sd_vecs = rep(list(sd), C),
                                           corrs = rep(corr, C),
                                           betas = rep(beta, C),
                                           precondition_matrices = rep(list(diag(1,2)), C),
                                           ESS_threshold = ESS_threshold,
                                           sub_posterior_means = t(sapply(input_samples, function(sub) apply(sub, 2, mean))),
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
                                              multivariate = TRUE,
                                              number_of_steps = length(gen_guide$mesh))
  d_BF_generalised <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                              N = nsamples,
                                              m = C,
                                              time_mesh = gen_guide$mesh,
                                              mean_vecs = means,
                                              sd_vecs = rep(list(sd), C),
                                              corrs = rep(corr, C),
                                              betas = rep(beta, C),
                                              precondition_matrices = lapply(input_samples, cov),
                                              ESS_threshold = ESS_threshold,
                                              sub_posterior_means = t(sapply(input_samples, function(sub) apply(sub, 2, mean))),
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
                                 'elapsed_time' = d_BF_standard$elapsed_time,
                                 'E_nu_j' = d_BF_standard$E_nu_j,
                                 'E_nu_j_old' = d_BF_standard$E_nu_j_old,
                                 'IAD' = integrated_abs_distance_biGaussian(fusion_post = resample_particle_y_samples(
                                   particle_set = d_BF_standard$particles,
                                   multivariate = TRUE,
                                   resampling_method = 'resid',
                                   seed = seed*i)$y_samples,
                                   marg_means = c(0,0),
                                   marg_sds = sqrt(rep(1, 2)/data_sizes[i]),
                                   bw = opt_bw))
  d_results$generalised[[i]] <- list('CESS_0' = d_BF_generalised$CESS[1],
                                     'CESS_j' = d_BF_generalised$CESS[2:length(d_BF_generalised$CESS)],
                                     'CESS_j_avg' = mean(d_BF_generalised$CESS[2:length(d_BF_generalised$CESS)]),
                                     'n' = length(d_BF_generalised$CESS),
                                     'time_mesh' = d_BF_generalised$particles$time_mesh,
                                     'time' = d_BF_generalised$time,
                                     'elapsed_time' = d_BF_generalised$elapsed_time,
                                     'E_nu_j' = d_BF_generalised$E_nu_j,
                                     'E_nu_j_old' = d_BF_generalised$E_nu_j_old,
                                     'IAD' = integrated_abs_distance_biGaussian(fusion_post = resample_particle_y_samples(
                                       particle_set = d_BF_generalised$particles,
                                       multivariate = TRUE,
                                       resampling_method = 'resid',
                                       seed = seed*i)$y_samples,
                                       marg_means = c(0,0),
                                       marg_sds = sqrt(rep(1, 2)/data_sizes[i]),
                                       bw = opt_bw))
  
  ##### SH: Recommended scaling of T, adaptive mesh #####
  print('### SH: performing standard Bayesian Fusion (with recommended T, adaptive mesh)')
  vanilla_guide_SH <- BF_guidance(condition = 'SH',
                                  CESS_0_threshold = CESS_0_threshold,
                                  CESS_j_threshold = CESS_j_threshold,
                                  sub_posterior_samples = input_samples,
                                  C = C,
                                  d = 2,
                                  data_size = data_sizes[i],
                                  b = vanilla_b,
                                  sub_posterior_means = t(sapply(input_samples, function(sub) apply(sub, 2, mean))),
                                  k1 = k1,
                                  k3 = k3,
                                  k4 = k4,
                                  vanilla = TRUE)
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = TRUE,
                                              number_of_steps = length(vanilla_guide_SH$mesh))
  SH_adaptive_standard <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                                  N = nsamples,
                                                  m = C,
                                                  time_mesh = vanilla_guide_SH$mesh,
                                                  mean_vecs = means,
                                                  sd_vecs = rep(list(sd), C),
                                                  corrs = rep(corr, C),
                                                  betas = rep(beta, C),
                                                  precondition_matrices = rep(list(diag(1,2)), C),
                                                  ESS_threshold = ESS_threshold,
                                                  sub_posterior_means = t(sapply(input_samples, function(sub) apply(sub, 2, mean))),
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
                              CESS_j_threshold = CESS_j_threshold,
                              sub_posterior_samples = input_samples,
                              C = C,
                              d = 2,
                              data_size = data_sizes[i],
                              sub_posterior_means = t(sapply(input_samples, function(sub) apply(sub, 2, mean))),
                              precondition_matrices = lapply(input_samples, cov),
                              inv_precondition_matrices = lapply(input_samples, function(sub) solve(cov(sub))),
                              Lambda = inverse_sum_matrices(lapply(input_samples, function(sub) solve(cov(sub)))),
                              k1 = k1,
                              k3 = k3,
                              k4 = k4,
                              vanilla = FALSE)
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = TRUE,
                                              number_of_steps = length(gen_guide_SH$mesh))
  SH_adaptive_generalised <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                                     N = nsamples,
                                                     m = C,
                                                     time_mesh = gen_guide_SH$mesh,
                                                     mean_vecs = means,
                                                     sd_vecs = rep(list(sd), C),
                                                     corrs = rep(corr, C),
                                                     betas = rep(beta, C),
                                                     precondition_matrices = lapply(input_samples, cov),
                                                     ESS_threshold = ESS_threshold,
                                                     sub_posterior_means = t(sapply(input_samples, function(sub) apply(sub, 2, mean))),
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
                                           'time' = SH_adaptive_standard$time,
                                           'elapsed_time' = SH_adaptive_standard$elapsed_time,
                                           'E_nu_j' = SH_adaptive_standard$E_nu_j,
                                           'E_nu_j_old' = SH_adaptive_standard$E_nu_j_old,
                                           'IAD' = integrated_abs_distance_biGaussian(fusion_post = resample_particle_y_samples(
                                             particle_set = SH_adaptive_standard$particles,
                                             multivariate = TRUE,
                                             resampling_method = 'resid',
                                             seed = seed*i)$y_samples,
                                             marg_means = c(0,0),
                                             marg_sds = sqrt(rep(1, 2)/data_sizes[i]),
                                             bw = opt_bw))
  SH_adaptive_results$generalised[[i]] <- list('CESS_0' = SH_adaptive_generalised$CESS[1],
                                               'CESS_j' = SH_adaptive_generalised$CESS[2:length(SH_adaptive_generalised$CESS)],
                                               'CESS_j_avg' = mean(SH_adaptive_generalised$CESS[2:length(SH_adaptive_generalised$CESS)]),
                                               'n' = length(SH_adaptive_generalised$CESS),
                                               'time_mesh' = SH_adaptive_generalised$particles$time_mesh,
                                               'time' = SH_adaptive_generalised$time,
                                               'elapsed_time' = SH_adaptive_generalised$elapsed_time,
                                               'E_nu_j' = SH_adaptive_generalised$E_nu_j,
                                               'E_nu_j_old' = SH_adaptive_generalised$E_nu_j_old,
                                               'IAD' = integrated_abs_distance_biGaussian(fusion_post = resample_particle_y_samples(
                                                 particle_set = SH_adaptive_generalised$particles,
                                                 multivariate = TRUE,
                                                 resampling_method = 'resid',
                                                 seed = seed*i)$y_samples,
                                                 marg_means = c(0,0),
                                                 marg_sds = sqrt(rep(1, 2)/data_sizes[i]),
                                                 bw = opt_bw))
}

##### vanilla plots #####

##### Fixed user-specified parameters #####
plot(x = data_sizes/100,
     y = sapply(1:length(data_sizes), function(i) a_results$vanilla[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes/100,
      y = sapply(1:length(data_sizes), function(i) a_results$vanilla[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) a_results$vanilla[[i]]$CESS_j/nsamples)[[ii]]
  points(x = rep(data_sizes[ii]/100, length(cess_j)), y = cess_j, cex = 0.5)
}
axis(1, at = seq(0, 2500, 500)/100, labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250)/100, labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('CESS / N', 2, 2.75, font = 2, cex = 1.5)

##### Recommended scaling of T, fixed n #####
plot(x = data_sizes/100,
     y = sapply(1:length(data_sizes), function(i) b_results$vanilla[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes/100,
      y = sapply(1:length(data_sizes), function(i) b_results$vanilla[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) b_results$vanilla[[i]]$CESS_j/nsamples)[[ii]]
  points(x = rep(data_sizes[ii]/100, length(cess_j)), y = cess_j, cex = 0.5)
}
axis(1, at = seq(0, 2500, 500)/100, labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250)/100, labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('CESS / N', 2, 2.75, font = 2, cex = 1.5)

##### Recommended scaling of T, regular mesh #####
plot(x = data_sizes/100,
     y = sapply(1:length(data_sizes), function(i) c_results$vanilla[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes/100,
      y = sapply(1:length(data_sizes), function(i) c_results$vanilla[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) c_results$vanilla[[i]]$CESS_j/nsamples)[[ii]]
  points(x = rep(data_sizes[ii]/100, length(cess_j)), y = cess_j, cex = 0.5)
}
axis(1, at = seq(0, 2500, 500)/100, labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250)/100, labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('CESS / N', 2, 2.75, font = 2, cex = 1.5)

##### Recommended scaling of T, adaptive mesh #####
plot(x = data_sizes/100,
     y = sapply(1:length(data_sizes), function(i) d_results$vanilla[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes/100,
      y = sapply(1:length(data_sizes), function(i) d_results$vanilla[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) d_results$vanilla[[i]]$CESS_j/nsamples)[[ii]]
  points(x = rep(data_sizes[ii]/100, length(cess_j)), y = cess_j, cex = 0.5)
}
axis(1, at = seq(0, 2500, 500)/100, labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250)/100, labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('CESS / N', 2, 2.75, font = 2, cex = 1.5)

##### Compare regular mesh and adaptive mesh times #####
plot(x = c_results$vanilla[[1]]$time_mesh,
     y = rep(data_sizes[1]-25, length(c_results$vanilla[[1]]$time_mesh)),
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,2500),
     xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
for (i in  2:length(data_sizes)) {
  lines(x = c_results$vanilla[[i]]$time_mesh,
        y = rep(data_sizes[i]-25, length(c_results$vanilla[[i]]$time_mesh)),
        type = 'b', pch = 20, lty = 1, lwd = 3)
}
for (i in  1:length(data_sizes)) {
  lines(x = d_results$vanilla[[i]]$time_mesh,
        y = rep(data_sizes[i]+25, length(d_results$vanilla[[i]]$time_mesh)),
        type = 'b', pch = 4, lty = 2, lwd = 3)
}
max(sapply(1:length(data_sizes), function(i) max(c_results$vanilla[[i]]$time_mesh)))
axis(1, at = seq(0, 0.03, 0.01), labels = seq(0, 0.03, 0.01), font = 2, cex = 1.5)
axis(1, at = seq(0, 0.04, 0.005), labels = rep("", 9), lwd.ticks = 0.5)
mtext('Time', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(2, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
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

##### SH: Recommended scaling of T, adaptive mesh #####
plot(x = data_sizes/100,
     y = sapply(1:length(data_sizes), function(i) SH_adaptive_results$vanilla[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes/100,
      y = sapply(1:length(data_sizes), function(i) SH_adaptive_results$vanilla[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) SH_adaptive_results$vanilla[[i]]$CESS_j/nsamples)[[ii]]
  points(x = rep(data_sizes[ii]/100, length(cess_j)), y = cess_j, cex = 0.5)
}
axis(1, at = seq(0, 2500, 500)/100, labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250)/100, labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('CESS / N', 2, 2.75, font = 2, cex = 1.5)

##### IAD #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) a_results$vanilla[[i]]$IAD),
     type = 'b', pch = 1, lty = 1, lwd = 3, ylim = c(0,1.2), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
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
      y = sapply(1:length(data_sizes), function(i) SH_adaptive_results$vanilla[[i]]$IAD),
      pch = 5, lty = 5, lwd = 3, type = 'b')
axis(1, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1.2, 0.1), labels = c("0.0", seq(0.1, 0.9, 0.1), "1.0", seq(1.1, 1.2, 0.1)),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1.2, 0.1), labels=rep("", 13), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
legend(x = 250, y = 1.2,
       legend = c('Fixed T, fixed n',
                  'SSH rec. T, fixed n',
                  'SSH rec. T, reg. mesh',
                  'SSH rec. T, adapt. mesh',
                  'SH rec. T, adapt. mesh'),
       lty = 1:5,
       pch = 1:5,
       lwd = rep(3, 5),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

##### time #####
plot(x = data_sizes,
     y = log(sapply(1:length(data_sizes), function(i) a_results$vanilla[[i]]$time)),
     type = 'b', pch = 1, lty = 1, lwd = 3, ylim = c(0,10), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
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
      y = log(sapply(1:length(data_sizes), function(i) SH_adaptive_results$vanilla[[i]]$time)),
      pch = 5, lty = 5, lwd = 3, type = 'b')
axis(1, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 10, 1), labels = seq(0, 10, 1), font = 2, cex = 1.5)
mtext('log(Elapsed time in seconds)', 2, 2.75, font = 2, cex = 1.5)
legend(x = 250, y = 10,
       legend = c('Fixed T, fixed n',
                  'SSH rec. T, fixed n',
                  'SSH rec. T, reg. mesh',
                  'SSH rec. T, adapt. mesh',
                  'SH rec. T, adapt. mesh'),
       lty = 1:5,
       pch = 1:5,
       lwd = rep(3, 5),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

##### generalised plots #####

##### Fixed user-specified parameters #####
plot(x = data_sizes/100,
     y = sapply(1:length(data_sizes), function(i) a_results$generalised[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes/100,
      y = sapply(1:length(data_sizes), function(i) a_results$generalised[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) a_results$generalised[[i]]$CESS_j/nsamples)[[ii]]
  points(x = rep(data_sizes[ii]/100, length(cess_j)), y = cess_j, cex = 0.5)
}
axis(1, at = seq(0, 2500, 500)/100, labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250)/100, labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('CESS / N', 2, 2.75, font = 2, cex = 1.5)

##### Recommended scaling of T, fixed n #####
plot(x = data_sizes/100,
     y = sapply(1:length(data_sizes), function(i) b_results$generalised[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes/100,
      y = sapply(1:length(data_sizes), function(i) b_results$generalised[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) b_results$generalised[[i]]$CESS_j/nsamples)[[ii]]
  points(x = rep(data_sizes[ii]/100, length(cess_j)), y = cess_j, cex = 0.5)
}
axis(1, at = seq(0, 2500, 500)/100, labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250)/100, labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('CESS / N', 2, 2.75, font = 2, cex = 1.5)

##### Recommended scaling of T, regular mesh #####
plot(x = data_sizes/100,
     y = sapply(1:length(data_sizes), function(i) c_results$generalised[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes/100,
      y = sapply(1:length(data_sizes), function(i) c_results$generalised[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) c_results$generalised[[i]]$CESS_j/nsamples)[[ii]]
  points(x = rep(data_sizes[ii]/100, length(cess_j)), y = cess_j, cex = 0.5)
}
axis(1, at = seq(0, 2500, 500)/100, labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250)/100, labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('CESS / N', 2, 2.75, font = 2, cex = 1.5)

##### Recommended scaling of T, adaptive mesh #####
plot(x = data_sizes/100,
     y = sapply(1:length(data_sizes), function(i) d_results$generalised[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes/100,
      y = sapply(1:length(data_sizes), function(i) d_results$generalised[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) d_results$generalised[[i]]$CESS_j/nsamples)[[ii]]
  points(x = rep(data_sizes[ii]/100, length(cess_j)), y = cess_j, cex = 0.5)
}
axis(1, at = seq(0, 2500, 500)/100, labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250)/100, labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('CESS / N', 2, 2.75, font = 2, cex = 1.5)

##### Compare regular mesh and adaptive mesh times #####
plot(x = c_results$generalised[[1]]$time_mesh,
     y = rep(data_sizes[1]-25, length(c_results$generalised[[1]]$time_mesh)),
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,2500), xlim = c(0, 11),
     xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
for (i in  2:length(data_sizes)) {
  lines(x = c_results$generalised[[i]]$time_mesh,
        y = rep(data_sizes[i]-25, length(c_results$generalised[[i]]$time_mesh)),
        type = 'b', pch = 20, lty = 1, lwd = 3)
}
for (i in  1:length(data_sizes)) {
  lines(x = d_results$generalised[[i]]$time_mesh,
        y = rep(data_sizes[i]+25, length(d_results$generalised[[i]]$time_mesh)),
        type = 'b', pch = 4, lty = 2, lwd = 3)
}
max(sapply(1:length(data_sizes), function(i) max(c_results$generalised[[i]]$time_mesh)))
axis(1, at = seq(0, 11, 1), labels = seq(0, 11, 1), font = 2, cex = 1.5)
axis(1, at = seq(0, 11, 0.5), labels = rep("", 23), lwd.ticks = 0.5)
mtext('Time', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(2, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 2, 2.75, font = 2, cex = 1.5)

##### SH: Recommended scaling of T, adaptive mesh #####
plot(x = data_sizes/100,
     y = sapply(1:length(data_sizes), function(i) SH_adaptive_results$generalised[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes/100,
      y = sapply(1:length(data_sizes), function(i) SH_adaptive_results$generalised[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) SH_adaptive_results$generalised[[i]]$CESS_j/nsamples)[[ii]]
  points(x = rep(data_sizes[ii]/100, length(cess_j)), y = cess_j, cex = 0.5)
}
axis(1, at = seq(0, 2500, 500)/100, labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250)/100, labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('CESS / N', 2, 2.75, font = 2, cex = 1.5)

##### IAD #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) a_results$generalised[[i]]$IAD),
     type = 'b', pch = 1, lty = 1, lwd = 3, ylim = c(0,1.2), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
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
      y = sapply(1:length(data_sizes), function(i) SH_adaptive_results$generalised[[i]]$IAD),
      pch = 5, lty = 5, lwd = 3, type = 'b')
axis(1, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1.2, 0.1), labels = c("0.0", seq(0.1, 0.9, 0.1), "1.0", seq(1.1, 1.2, 0.1)),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1.2, 0.1), labels=rep("", 13), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
legend(x = 250, y = 1.2,
       legend = c('Fixed T, fixed n',
                  'SSH rec. T, fixed n',
                  'SSH rec. T, reg. mesh',
                  'SSH rec. T, adapt. mesh',
                  'SH rec. T, adapt. mesh'),
       lty = 1:5,
       pch = 1:5,
       lwd = rep(3, 5),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

##### time #####
plot(x = data_sizes,
     y = log(sapply(1:length(data_sizes), function(i) a_results$generalised[[i]]$time)),
     type = 'b', pch = 1, lty = 1, lwd = 3, ylim = c(0,10), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
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
      y = log(sapply(1:length(data_sizes), function(i) SH_adaptive_results$generalised[[i]]$time)),
      pch = 5, lty = 5, lwd = 3, type = 'b')
axis(1, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 10, 1), labels = seq(0, 10, 1), font = 2, cex = 1.5)
mtext('log(Elapsed time in seconds)', 2, 2.75, font = 2, cex = 1.5)
legend(x = 250, y = 10,
       legend = c('Fixed T, fixed n',
                  'SSH rec. T, fixed n',
                  'SSH rec. T, reg. mesh',
                  'SSH rec. T, adapt. mesh',
                  'SH rec. T, adapt. mesh'),
       lty = 1:5,
       pch = 1:5,
       lwd = rep(3, 5),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

save.image('bf_bivG_dissimilar_means_example_09.RData')
