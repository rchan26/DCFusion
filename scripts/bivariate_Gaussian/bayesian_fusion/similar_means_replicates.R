library(DCFusion)

seed <- 1994
set.seed(seed)
nsamples <- 10000
C <- 10
mean <- rep(0, 2)
corr <- 0.9
beta <- 1
a_mesh_vanilla <- seq(0, 0.005, length.out = 6)
a_mesh_gen <- seq(0, 1, length.out = 6)
diffusion_estimator <- 'NB'
number_of_replicates <- 10
ESS_threshold <- 0.5
CESS_0_threshold <- 0.2
CESS_j_threshold <- 0.2
vanilla_b <- 1
k1 <- NULL
k2 <- NULL
k3 <- -log(CESS_j_threshold)/2
k4 <- -log(CESS_j_threshold)/2
data_sizes <- c(1000, 5000, 10000, 20000, 30000, 40000)
vanilla_guide <- list()
gen_guide <- list()
vanilla_guide_SSH <- list()
gen_guide_SSH <- list()
a_results <- list('vanilla' = list(), 'generalised' = list())
b_results <- list('vanilla' = list(), 'generalised' = list())
c_results <- list('vanilla' = list(), 'generalised' = list())
d_results <- list('vanilla' = list(), 'generalised' = list())
SSH_adaptive_results <- list('vanilla' = list(), 'generalised' = list())

for (i in 1:length(data_sizes)) {
  print(paste('i:', i))
  print(paste('data_size:', data_sizes[i]))
  vanilla_guide[[i]] <- list()
  gen_guide[[i]] <- list()
  vanilla_guide_SSH[[i]] <- list()
  gen_guide_SSH[[i]] <- list()
  a_results$vanilla[[i]] <- list()
  a_results$generalised[[i]] <- list()
  b_results$vanilla[[i]] <- list()
  b_results$generalised[[i]] <- list()
  c_results$vanilla[[i]] <- list()
  c_results$generalised[[i]] <- list()
  d_results$vanilla[[i]] <- list()
  d_results$generalised[[i]] <- list()
  SSH_adaptive_results$vanilla[[i]] <- list()
  SSH_adaptive_results$generalised[[i]] <- list()
  for (rep in 1:number_of_replicates) {
    print(paste('rep:', rep))
    set.seed(seed*rep*i)
    sd <- sqrt(rep(C, 2)/data_sizes[i])
    cov_mat <- matrix(c(sd[1]^2, sd[1]*sd[2]*corr, sd[1]*sd[2]*corr, sd[2]^2),
                      nrow = 2, ncol = 2, byrow = T)
    opt_bw <- ((4*sd^5)/(3*nsamples))^(1/5)
    input_samples <- lapply(1:C, function(sub) mvrnormArma(N = nsamples, mu = mean, Sigma = cov_mat))
    ##### Fixed user-specified parameters #####
    print('### performing standard Bayesian Fusion (with standard mesh)')
    input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                                multivariate = TRUE,
                                                number_of_steps = length(a_mesh_vanilla))
    a_BF_standard <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                             N = nsamples,
                                             m = C,
                                             time_mesh = a_mesh_vanilla,
                                             mean_vecs = rep(list(mean), C),
                                             sd_vecs = rep(list(sd), C),
                                             corrs = rep(corr, C),
                                             betas = rep(beta, C),
                                             precondition_matrices = rep(list(diag(1,2)), C),
                                             ESS_threshold = ESS_threshold,
                                             diffusion_estimator = diffusion_estimator,
                                             seed = seed*rep*i)
    print('### performing Bayesian Fusion with a preconditioning matrix (with standard mesh)')
    input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                                multivariate = TRUE,
                                                number_of_steps = length(a_mesh_gen))
    a_BF_generalised <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                                N = nsamples,
                                                m = C,
                                                time_mesh = a_mesh_gen,
                                                mean_vecs = rep(list(mean), C),
                                                sd_vecs = rep(list(sd), C),
                                                corrs = rep(corr, C),
                                                betas = rep(beta, C),
                                                precondition_matrices = lapply(input_samples, cov),
                                                ESS_threshold = ESS_threshold,
                                                diffusion_estimator = diffusion_estimator,
                                                seed = seed*rep*i)
    # save results
    a_results$vanilla[[i]][[rep]] <- list('CESS_0' = a_BF_standard$CESS[1],
                                          'CESS_j' = a_BF_standard$CESS[2:length(a_BF_standard$CESS)],
                                          'CESS_j_avg' = mean(a_BF_standard$CESS[2:length(a_BF_standard$CESS)]),
                                          'CESS_j_var' = var(a_BF_standard$CESS[2:length(a_BF_standard$CESS)]),
                                          'n' = length(a_BF_standard$CESS),
                                          'time_mesh' = a_BF_standard$particles$time_mesh,
                                          'time' = a_BF_standard$time,
                                          'elapsed_time' = a_BF_standard$elapsed_time,
                                          'resampled' = a_BF_standard$resampled,
                                          'ESS' = a_BF_standard$ESS,
                                          'IAD' = integrated_abs_distance_biGaussian(fusion_post = resample_particle_y_samples(
                                            particle_set = a_BF_standard$particles,
                                            multivariate = TRUE,
                                            resampling_method = 'resid',
                                            seed = seed*rep*i)$y_samples,
                                            marg_means = mean,
                                            marg_sds = sqrt(rep(1, 2)/data_sizes[i]),
                                            bw = opt_bw))
    a_results$generalised[[i]][[rep]] <- list('CESS_0' = a_BF_generalised$CESS[1],
                                              'CESS_j' = a_BF_generalised$CESS[2:length(a_BF_generalised$CESS)],
                                              'CESS_j_avg' = mean(a_BF_generalised$CESS[2:length(a_BF_generalised$CESS)]),
                                              'CESS_j_var' = var(a_BF_generalised$CESS[2:length(a_BF_generalised$CESS)]),
                                              'n' = length(a_BF_generalised$CESS),
                                              'time_mesh' = a_BF_generalised$particles$time_mesh,
                                              'time' = a_BF_generalised$time,
                                              'elapsed_time' = a_BF_generalised$elapsed_time,
                                              'resampled' = a_BF_generalised$resampled,
                                              'ESS' = a_BF_generalised$ESS,
                                              'IAD' = integrated_abs_distance_biGaussian(fusion_post = resample_particle_y_samples(
                                                particle_set = a_BF_generalised$particles,
                                                multivariate = TRUE,
                                                resampling_method = 'resid',
                                                seed = seed*rep*i)$y_samples,
                                                marg_means = mean,
                                                marg_sds = sqrt(rep(1, 2)/data_sizes[i]),
                                                bw = opt_bw))
    
    ##### Recommended scaling of T, fixed n #####
    print('### performing standard Bayesian Fusion (with recommended T, fixed n)')
    vanilla_guide[[i]][[rep]] <- BF_guidance(condition = 'SH',
                                             CESS_0_threshold = CESS_0_threshold,
                                             CESS_j_threshold = CESS_j_threshold,
                                             sub_posterior_samples = input_samples,
                                             C = C,
                                             d = 2,
                                             data_size = data_sizes[i],
                                             b = vanilla_b,
                                             sub_posterior_means = t(sapply(input_samples, function(sub) apply(sub, 2, mean))),
                                             precondition_matrices = rep(list(diag(1,2)), C),
                                             inv_precondition_matrices = rep(list(diag(1,2)), C),
                                             Lambda = inverse_sum_matrices(rep(list(diag(1,2)), C)),
                                             k1 = k1,
                                             vanilla = TRUE)
    print(paste('vanilla recommened regular mesh n:', vanilla_guide[[i]][[rep]]$n))
    b_mesh_vanilla <- seq(0, vanilla_guide[[i]][[rep]]$min_T, length.out = 6)
    input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                                multivariate = TRUE,
                                                number_of_steps = length(b_mesh_vanilla))
    b_BF_standard <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                             N = nsamples,
                                             m = C,
                                             time_mesh = b_mesh_vanilla,
                                             mean_vecs = rep(list(mean), C),
                                             sd_vecs = rep(list(sd), C),
                                             corrs = rep(corr, C),
                                             betas = rep(beta, C),
                                             precondition_matrices = rep(list(diag(1,2)), C),
                                             ESS_threshold = ESS_threshold,
                                             diffusion_estimator = diffusion_estimator,
                                             seed = seed*rep*i)
    print('### performing Bayesian Fusion with a preconditioning matrix (with recommended T, fixed n)')
    gen_guide[[i]][[rep]] <- BF_guidance(condition = 'SH',
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
                                         vanilla = FALSE)
    print(paste('generalised recommened regular mesh n:', gen_guide[[i]][[rep]]$n))
    b_mesh_gen <- seq(0, gen_guide[[i]][[rep]]$min_T, length.out = 6)
    input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                                multivariate = TRUE,
                                                number_of_steps = length(b_mesh_gen))
    b_BF_generalised <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                                N = nsamples,
                                                m = C,
                                                time_mesh = b_mesh_gen,
                                                mean_vecs = rep(list(mean), C),
                                                sd_vecs = rep(list(sd), C),
                                                corrs = rep(corr, C),
                                                betas = rep(beta, C),
                                                precondition_matrices = lapply(input_samples, cov),
                                                ESS_threshold = ESS_threshold,
                                                diffusion_estimator = diffusion_estimator,
                                                seed = seed*rep*i)
    # save results
    b_results$vanilla[[i]][[rep]] <- list('CESS_0' = b_BF_standard$CESS[1],
                                          'CESS_j' = b_BF_standard$CESS[2:length(b_BF_standard$CESS)],
                                          'CESS_j_avg' = mean(b_BF_standard$CESS[2:length(b_BF_standard$CESS)]),
                                          'CESS_j_var' = var(b_BF_standard$CESS[2:length(b_BF_standard$CESS)]),
                                          'n' = length(b_BF_standard$CESS),
                                          'time_mesh' = b_BF_standard$particles$time_mesh,
                                          'time' = b_BF_standard$time,
                                          'elapsed_time' = b_BF_standard$elapsed_time,
                                          'resampled' = b_BF_standard$resampled,
                                          'ESS' = b_BF_standard$ESS,
                                          'IAD' = integrated_abs_distance_biGaussian(fusion_post = resample_particle_y_samples(
                                            particle_set = b_BF_standard$particles,
                                            multivariate = TRUE,
                                            resampling_method = 'resid',
                                            seed = seed*rep*i)$y_samples,
                                            marg_means = mean,
                                            marg_sds = sqrt(rep(1, 2)/data_sizes[i]),
                                            bw = opt_bw))
    b_results$generalised[[i]][[rep]] <- list('CESS_0' = b_BF_generalised$CESS[1],
                                              'CESS_j' = b_BF_generalised$CESS[2:length(b_BF_generalised$CESS)],
                                              'CESS_j_avg' = mean(b_BF_generalised$CESS[2:length(b_BF_generalised$CESS)]),
                                              'CESS_j_var' = var(b_BF_generalised$CESS[2:length(b_BF_generalised$CESS)]),
                                              'n' = length(b_BF_generalised$CESS),
                                              'time_mesh' = b_BF_generalised$particles$time_mesh,
                                              'time' = b_BF_generalised$time,
                                              'elapsed_time' = b_BF_generalised$elapsed_time,
                                              'resampled' = b_BF_generalised$resampled,
                                              'ESS' = b_BF_generalised$ESS,
                                              'IAD' = integrated_abs_distance_biGaussian(fusion_post = resample_particle_y_samples(
                                                particle_set = b_BF_generalised$particles,
                                                multivariate = TRUE,
                                                resampling_method = 'resid',
                                                seed = seed*rep*i)$y_samples,
                                                marg_means = mean,
                                                marg_sds = sqrt(rep(1, 2)/data_sizes[i]),
                                                bw = opt_bw))
    
    ##### Recommended scaling of T, regular mesh #####
    print('### performing standard Bayesian Fusion (with recommended T, regular mesh)')
    input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                                multivariate = TRUE,
                                                number_of_steps = length(vanilla_guide[[i]][[rep]]$mesh))
    c_BF_standard <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                             N = nsamples,
                                             m = C,
                                             time_mesh = vanilla_guide[[i]][[rep]]$mesh,
                                             mean_vecs = rep(list(mean), C),
                                             sd_vecs = rep(list(sd), C),
                                             corrs = rep(corr, C),
                                             betas = rep(beta, C),
                                             precondition_matrices = rep(list(diag(1,2)), C),
                                             ESS_threshold = ESS_threshold,
                                             diffusion_estimator = diffusion_estimator,
                                             seed = seed*rep*i)
    print('### performing Bayesian Fusion with a preconditioning matrix (with recommended T, regular mesh)')
    input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                                multivariate = TRUE,
                                                number_of_steps = length(gen_guide[[i]][[rep]]$mesh))
    c_BF_generalised <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                                N = nsamples,
                                                m = C,
                                                time_mesh = gen_guide[[i]][[rep]]$mesh,
                                                mean_vecs = rep(list(mean), C),
                                                sd_vecs = rep(list(sd), C),
                                                corrs = rep(corr, C),
                                                betas = rep(beta, C),
                                                precondition_matrices = lapply(input_samples, cov),
                                                ESS_threshold = ESS_threshold,
                                                diffusion_estimator = diffusion_estimator,
                                                seed = seed*rep*i)
    # save results
    c_results$vanilla[[i]][[rep]] <- list('CESS_0' = c_BF_standard$CESS[1],
                                          'CESS_j' = c_BF_standard$CESS[2:length(c_BF_standard$CESS)],
                                          'CESS_j_avg' = mean(c_BF_standard$CESS[2:length(c_BF_standard$CESS)]),
                                          'CESS_j_var' = var(c_BF_standard$CESS[2:length(c_BF_standard$CESS)]),
                                          'n' = length(c_BF_standard$CESS),
                                          'time_mesh' = c_BF_standard$particles$time_mesh,
                                          'time' = c_BF_standard$time,
                                          'elapsed_time' = c_BF_standard$elapsed_time,
                                          'resampled' = c_BF_standard$resampled,
                                          'ESS' = c_BF_standard$ESS,
                                          'IAD' = integrated_abs_distance_biGaussian(fusion_post = resample_particle_y_samples(
                                            particle_set = c_BF_standard$particles,
                                            multivariate = TRUE,
                                            resampling_method = 'resid',
                                            seed = seed*rep*i)$y_samples,
                                            marg_means = mean,
                                            marg_sds = sqrt(rep(1, 2)/data_sizes[i]),
                                            bw = opt_bw))
    c_results$generalised[[i]][[rep]] <- list('CESS_0' = c_BF_generalised$CESS[1],
                                              'CESS_j' = c_BF_generalised$CESS[2:length(c_BF_generalised$CESS)],
                                              'CESS_j_avg' = mean(c_BF_generalised$CESS[2:length(c_BF_generalised$CESS)]),
                                              'CESS_j_var' = var(c_BF_generalised$CESS[2:length(c_BF_generalised$CESS)]),
                                              'n' = length(c_BF_generalised$CESS),
                                              'time_mesh' = c_BF_generalised$particles$time_mesh,
                                              'time' = c_BF_generalised$time,
                                              'elapsed_time' = c_BF_generalised$elapsed_time,
                                              'resampled' = c_BF_generalised$resampled,
                                              'ESS' = c_BF_generalised$ESS,
                                              'IAD' = integrated_abs_distance_biGaussian(fusion_post = resample_particle_y_samples(
                                                particle_set = c_BF_generalised$particles,
                                                multivariate = TRUE,
                                                resampling_method = 'resid',
                                                seed = seed*rep*i)$y_samples,
                                                marg_means = mean,
                                                marg_sds = sqrt(rep(1, 2)/data_sizes[i]),
                                                bw = opt_bw))
    
    ##### Recommended scaling of T, adaptive mesh #####
    print('### performing standard Bayesian Fusion (with recommended T, adaptive mesh)')
    input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                                multivariate = TRUE,
                                                number_of_steps = length(vanilla_guide[[i]][[rep]]$mesh))
    d_BF_standard <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                             N = nsamples,
                                             m = C,
                                             time_mesh = vanilla_guide[[i]][[rep]]$mesh,
                                             mean_vecs = rep(list(mean), C),
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
                                             seed = seed*rep*i)
    print('### performing Bayesian Fusion with a preconditioning matrix (with recommended T, adaptive mesh)')
    input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                                multivariate = TRUE,
                                                number_of_steps = length(gen_guide[[i]][[rep]]$mesh))
    d_BF_generalised <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                                N = nsamples,
                                                m = C,
                                                time_mesh = gen_guide[[i]][[rep]]$mesh,
                                                mean_vecs = rep(list(mean), C),
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
                                                seed = seed*rep*i)
    # save results
    d_results$vanilla[[i]][[rep]] <- list('CESS_0' = d_BF_standard$CESS[1],
                                          'CESS_j' = d_BF_standard$CESS[2:length(d_BF_standard$CESS)],
                                          'CESS_j_avg' = mean(d_BF_standard$CESS[2:length(d_BF_standard$CESS)]),
                                          'CESS_j_var' = var(d_BF_standard$CESS[2:length(d_BF_standard$CESS)]),
                                          'n' = length(d_BF_standard$CESS),
                                          'time_mesh' = d_BF_standard$particles$time_mesh,
                                          'time' = d_BF_standard$time,
                                          'elapsed_time' = d_BF_standard$elapsed_time,
                                          'resampled' = d_BF_standard$resampled,
                                          'ESS' = d_BF_standard$ESS,
                                          'E_nu_j' = d_BF_standard$E_nu_j,
                                          'E_nu_j_old' = d_BF_standard$E_nu_j_old,
                                          'IAD' = integrated_abs_distance_biGaussian(fusion_post = resample_particle_y_samples(
                                            particle_set = d_BF_standard$particles,
                                            multivariate = TRUE,
                                            resampling_method = 'resid',
                                            seed = seed*rep*i)$y_samples,
                                            marg_means = mean,
                                            marg_sds = sqrt(rep(1, 2)/data_sizes[i]),
                                            bw = opt_bw))
    d_results$generalised[[i]][[rep]] <- list('CESS_0' = d_BF_generalised$CESS[1],
                                              'CESS_j' = d_BF_generalised$CESS[2:length(d_BF_generalised$CESS)],
                                              'CESS_j_avg' = mean(d_BF_generalised$CESS[2:length(d_BF_generalised$CESS)]),
                                              'CESS_j_var' = var(d_BF_generalised$CESS[2:length(d_BF_generalised$CESS)]),
                                              'n' = length(d_BF_generalised$CESS),
                                              'time_mesh' = d_BF_generalised$particles$time_mesh,
                                              'time' = d_BF_generalised$time,
                                              'elapsed_time' = d_BF_generalised$elapsed_time,
                                              'resampled' = d_BF_generalised$resampled,
                                              'ESS' = d_BF_generalised$ESS,
                                              'E_nu_j' = d_BF_generalised$E_nu_j,
                                              'E_nu_j_old' = d_BF_generalised$E_nu_j_old,
                                              'IAD' = integrated_abs_distance_biGaussian(fusion_post = resample_particle_y_samples(
                                                particle_set = d_BF_generalised$particles,
                                                multivariate = TRUE,
                                                resampling_method = 'resid',
                                                seed = seed*rep*i)$y_samples,
                                                marg_means = mean,
                                                marg_sds = sqrt(rep(1, 2)/data_sizes[i]),
                                                bw = opt_bw))
    
    ##### SSH: Recommended scaling of T, adaptive mesh #####
    print('### SSH: performing standard Bayesian Fusion (with recommended T, adaptive mesh)')
    vanilla_guide_SSH[[i]][[rep]] <- BF_guidance(condition = 'SSH',
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
    input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                                multivariate = TRUE,
                                                number_of_steps = length(vanilla_guide_SSH[[i]][[rep]]$mesh))
    SSH_adaptive_standard <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                                     N = nsamples,
                                                     m = C,
                                                     time_mesh = vanilla_guide_SSH[[i]][[rep]]$mesh,
                                                     mean_vecs = rep(list(mean), C),
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
                                                     seed = seed*rep*i)
    print('### SSH: performing Bayesian Fusion with a preconditioning matrix (with recommended T, adaptive mesh)')
    gen_guide_SSH[[i]][[rep]] <- BF_guidance(condition = 'SSH',
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
    input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                                multivariate = TRUE,
                                                number_of_steps = length(gen_guide_SSH[[i]][[rep]]$mesh))
    SSH_adaptive_generalised <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                                        N = nsamples,
                                                        m = C,
                                                        time_mesh = gen_guide_SSH[[i]][[rep]]$mesh,
                                                        mean_vecs = rep(list(mean), C),
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
                                                        seed = seed*rep*i)
    # save results
    SSH_adaptive_results$vanilla[[i]][[rep]] <- list('CESS_0' = SSH_adaptive_standard$CESS[1],
                                                     'CESS_j' = SSH_adaptive_standard$CESS[2:length(SSH_adaptive_standard$CESS)],
                                                     'CESS_j_avg' = mean(SSH_adaptive_standard$CESS[2:length(SSH_adaptive_standard$CESS)]),
                                                     'CESS_j_var' = var(SSH_adaptive_standard$CESS[2:length(SSH_adaptive_standard$CESS)]),
                                                     'n' = length(SSH_adaptive_standard$CESS),
                                                     'time_mesh' = SSH_adaptive_standard$particles$time_mesh,
                                                     'time' = SSH_adaptive_standard$time,
                                                     'elapsed_time' = SSH_adaptive_standard$elapsed_time,
                                                     'resampled' = SSH_adaptive_standard$resampled,
                                                     'ESS' = SSH_adaptive_standard$ESS,
                                                     'E_nu_j' = SSH_adaptive_standard$E_nu_j,
                                                     'E_nu_j_old' = SSH_adaptive_standard$E_nu_j_old,
                                                     'IAD' = integrated_abs_distance_biGaussian(fusion_post = resample_particle_y_samples(
                                                       particle_set = SSH_adaptive_standard$particles,
                                                       multivariate = TRUE,
                                                       resampling_method = 'resid',
                                                       seed = seed*rep*i)$y_samples,
                                                       marg_means = mean,
                                                       marg_sds = sqrt(rep(1, 2)/data_sizes[i]),
                                                       bw = opt_bw))
    SSH_adaptive_results$generalised[[i]][[rep]] <- list('CESS_0' = SSH_adaptive_generalised$CESS[1],
                                                         'CESS_j' = SSH_adaptive_generalised$CESS[2:length(SSH_adaptive_generalised$CESS)],
                                                         'CESS_j_avg' = mean(SSH_adaptive_generalised$CESS[2:length(SSH_adaptive_generalised$CESS)]),
                                                         'n' = length(SSH_adaptive_generalised$CESS),
                                                         'time_mesh' = SSH_adaptive_generalised$particles$time_mesh,
                                                         'time' = SSH_adaptive_generalised$time,
                                                         'elapsed_time' = SSH_adaptive_generalised$elapsed_time,
                                                         'resampled' = SSH_adaptive_generalised$resampled,
                                                         'ESS' = SSH_adaptive_generalised$ESS,
                                                         'E_nu_j' = SSH_adaptive_generalised$E_nu_j,
                                                         'E_nu_j_old' = SSH_adaptive_generalised$E_nu_j_old,
                                                         'IAD' = integrated_abs_distance_biGaussian(fusion_post = resample_particle_y_samples(
                                                           particle_set = SSH_adaptive_generalised$particles,
                                                           multivariate = TRUE,
                                                           resampling_method = 'resid',
                                                           seed = seed*rep*i)$y_samples,
                                                           marg_means = mean,
                                                           marg_sds = sqrt(rep(1, 2)/data_sizes[i]),
                                                           bw = opt_bw))
    # save progress
    print('saving progress')
    save.image('bivG_similar_means_replicates.RData')
  }
}

##### IAD #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) {
       mean(sapply(1:number_of_replicates, function(rep) a_results$vanilla[[i]][[rep]]$IAD))
     }),
     type = 'b', pch = 1, lty = 1, lwd = 3, ylim = c(0,1.6), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
for (i in 1:length(data_sizes)) {
  IAD <- sapply(1:number_of_replicates, function(rep) a_results$vanilla[[i]][[rep]]$IAD)
  points(x = rep(data_sizes[i], length(IAD)), y = IAD, cex = 0.5, pch = 1)
}
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) {
        mean(sapply(1:number_of_replicates, function(rep) b_results$vanilla[[i]][[rep]]$IAD))
      }),
      pch = 2, lty = 2, lwd = 3, type = 'b', col = 'blue')
for (i in 1:length(data_sizes)) {
  IAD <- sapply(1:number_of_replicates, function(rep) b_results$vanilla[[i]][[rep]]$IAD)
  points(x = rep(data_sizes[i], length(IAD)), y = IAD, cex = 0.5, pch = 2, col = 'blue')
}
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) {
        mean(sapply(1:number_of_replicates, function(rep) c_results$vanilla[[i]][[rep]]$IAD))
      }),
      pch = 3, lty = 3, lwd = 3, type = 'b', col = 'red')
for (i in 1:length(data_sizes)) {
  IAD <- sapply(1:number_of_replicates, function(rep) c_results$vanilla[[i]][[rep]]$IAD)
  points(x = rep(data_sizes[i], length(IAD)), y = IAD, cex = 0.5, pch = 3, col = 'red')
}
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) {
        mean(sapply(1:number_of_replicates, function(rep) d_results$vanilla[[i]][[rep]]$IAD))
      }),
      pch = 4, lty = 4, lwd = 3, type = 'b', col = 'green')
for (i in 1:length(data_sizes)) {
  IAD <- sapply(1:number_of_replicates, function(rep) d_results$vanilla[[i]][[rep]]$IAD)
  points(x = rep(data_sizes[i], length(IAD)), y = IAD, cex = 0.5, pch = 4, col = 'green')
}
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) {
        mean(sapply(1:number_of_replicates, function(rep) SH_adaptive_results$vanilla[[i]][[rep]]$IAD))
      }),
      pch = 5, lty = 5, lwd = 3, type = 'b', col = 'purple')
for (i in 1:length(data_sizes)) {
  IAD <- sapply(1:number_of_replicates, function(rep) SH_adaptive_results$vanilla[[i]][[rep]]$IAD)
  points(x = rep(data_sizes[i], length(IAD)), y = IAD, cex = 0.5, pch = 5, col = 'purple')
}
axis(1, at = c(1000, seq(10000, 50000, 10000)),
     labels = c(1000, seq(10000, 50000, 10000)), font = 2, cex = 1.5)
axis(1, at = seq(0, 50000, 5000), labels = rep("", 11), lwd.ticks = 0.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1.6, 0.1), labels = c("0.0", seq(0.1, 0.9, 0.1), "1.0", seq(1.1, 1.6, 0.1)),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1.6, 0.1), labels=rep("", 17), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
legend(x = 1000, y = 1.6,
       legend = c('Fixed T, fixed n',
                  'SH rec. T, fixed n',
                  'SH rec. T, reg. mesh',
                  'SH rec. T, adapt. mesh',
                  'SSH rec. T, adapt. mesh'),
       col = c('black', 'blue', 'red', 'green', 'purple'),
       lty = 1:5,
       pch = 1:5,
       lwd = rep(3, 5),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

##### IAD #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) {
       mean(sapply(1:number_of_replicates, function(rep) a_results$generalised[[i]][[rep]]$IAD))
     }),
     type = 'b', pch = 1, lty = 1, lwd = 3, ylim = c(0,1.6), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
for (i in 1:length(data_sizes)) {
  IAD <- sapply(1:number_of_replicates, function(rep) a_results$generalised[[i]][[rep]]$IAD)
  points(x = rep(data_sizes[i], length(IAD)), y = IAD, cex = 0.5, pch = 1)
}
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) {
        mean(sapply(1:number_of_replicates, function(rep) b_results$generalised[[i]][[rep]]$IAD))
      }),
      pch = 2, lty = 2, lwd = 3, type = 'b', col = 'blue')
for (i in 1:length(data_sizes)) {
  IAD <- sapply(1:number_of_replicates, function(rep) b_results$generalised[[i]][[rep]]$IAD)
  points(x = rep(data_sizes[i], length(IAD)), y = IAD, cex = 0.5, pch = 2, col = 'blue')
}
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) {
        mean(sapply(1:number_of_replicates, function(rep) c_results$generalised[[i]][[rep]]$IAD))
      }),
      pch = 3, lty = 3, lwd = 3, type = 'b', col = 'red')
for (i in 1:length(data_sizes)) {
  IAD <- sapply(1:number_of_replicates, function(rep) c_results$generalised[[i]][[rep]]$IAD)
  points(x = rep(data_sizes[i], length(IAD)), y = IAD, cex = 0.5, pch = 3, col = 'red')
}
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) {
        mean(sapply(1:number_of_replicates, function(rep) d_results$generalised[[i]][[rep]]$IAD))
      }),
      pch = 4, lty = 4, lwd = 3, type = 'b', col = 'green')
for (i in 1:length(data_sizes)) {
  IAD <- sapply(1:number_of_replicates, function(rep) d_results$generalised[[i]][[rep]]$IAD)
  points(x = rep(data_sizes[i], length(IAD)), y = IAD, cex = 0.5, pch = 4, col = 'green')
}
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) {
        mean(sapply(1:number_of_replicates, function(rep) SH_adaptive_results$generalised[[i]][[rep]]$IAD))
      }),
      pch = 5, lty = 5, lwd = 3, type = 'b', col = 'purple')
for (i in 1:length(data_sizes)) {
  IAD <- sapply(1:number_of_replicates, function(rep) SH_adaptive_results$generalised[[i]][[rep]]$IAD)
  points(x = rep(data_sizes[i], length(IAD)), y = IAD, cex = 0.5, pch = 5, col = 'purple')
}
axis(1, at = c(1000, seq(10000, 50000, 10000)),
     labels = c(1000, seq(10000, 50000, 10000)), font = 2, cex = 1.5)
axis(1, at = seq(0, 50000, 5000), labels = rep("", 11), lwd.ticks = 0.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1.6, 0.1), labels = c("0.0", seq(0.1, 0.9, 0.1), "1.0", seq(1.1, 1.6, 0.1)),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1.6, 0.1), labels=rep("", 17), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
legend(x = 1000, y = 1.6,
       legend = c('Fixed T, fixed n',
                  'SH rec. T, fixed n',
                  'SH rec. T, reg. mesh',
                  'SH rec. T, adapt. mesh',
                  'SSH rec. T, adapt. mesh'),
       col = c('black', 'blue', 'red', 'green', 'purple'),
       lty = 1:5,
       pch = 1:5,
       lwd = rep(3, 5),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

