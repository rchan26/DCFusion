library(DCFusion)

seed <- 1996
set.seed(seed)
nsamples <- 10000
C <- 2
means <- list(rep(-0.25, 2), rep(0.25, 2))
corr <- 0.9
beta <- 1
a_mesh_vanilla <- seq(0, 0.01, length.out = 6)
a_mesh_gen <- seq(0, 2, length.out = 6)
diffusion_estimator <- 'NB'
resampling_method <- 'resid'
ESS_threshold <- 0.5
CESS_0_threshold <- 0.5
CESS_j_threshold <- 0.5
vanilla_b <- 1
k1 <- NULL
k2 <- NULL
k3 <- -log(CESS_j_threshold)/2
k4 <- -log(CESS_j_threshold)/2
data_sizes <- c(250, 500, 1000, 1500, 2000, 2500)
vanilla_guide <- list()
gen_guide <- list()
vanilla_guide_SH <- list()
gen_guide_SH <- list()
a_results <- list('vanilla' = list(), 'generalised' = list())
b_results <- list('vanilla' = list(), 'generalised' = list())
c_results <- list('vanilla' = list(), 'generalised' = list())
c_check_results <- list('vanilla' = list(), 'generalised' = list())
d1_results <- list('vanilla' = list(), 'generalised' = list())
d2_results <- list('vanilla' = list(), 'generalised' = list())
SH_adaptive_results <- list('vanilla' = list(), 'generalised' = list())

collect_results <- function(results, seed) {
  print(paste('n:', length(results$CESS)-1))
  print(paste('time:', results$time))
  print(paste('log(time):', log(results$time)))
  return(list('CESS_0' = results$CESS[1],
              'CESS_j' = results$CESS[2:length(results$CESS)],
              'CESS_j_avg' = mean(results$CESS[2:length(results$CESS)]),
              'n' = length(results$CESS)-1,
              'time_mesh' = results$particles$time_mesh,
              'time' = results$time,
              'elapsed_time' = results$elapsed_time,
              'resampled' = results$resampled,
              'ESS' = results$ESS,
              'E_nu_j' = results$E_nu_j,
              'chosen' = results$chosen,
              'mesh_terms' = results$mesh_terms,
              'k4_choice' = results$k4_choice,
              'IAD' = integrated_abs_distance_biGaussian(fusion_post = resample_particle_y_samples(
                particle_set = results$particles,
                multivariate = TRUE,
                resampling_method = resampling_method,
                seed = seed)$y_samples,
                marg_means = c(0,0),
                marg_sds = sqrt(rep(1, 2)/data_sizes[i]),
                bw = opt_bw)))
}

for (i in 1:length(data_sizes)) {
  print(paste('i:', i))
  print(paste('data size:', data_sizes[i]))
  set.seed(seed*i)
  sd <- sqrt(rep(C, 2)/data_sizes[i])
  cov_mat <- matrix(c(sd[1]^2, sd[1]*sd[2]*corr, sd[1]*sd[2]*corr, sd[2]^2),
                    nrow = 2, ncol = 2, byrow = T)
  opt_bw <- ((4*sd^5)/(3*nsamples))^(1/5)
  input_samples <- lapply(1:C, function(sub) mvrnormArma(N = nsamples, mu = means[[sub]], Sigma = cov_mat))
  
  #### Fixed user-specified parameters #####
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
  print('### performing Generalised Bayesian Fusion (with standard mesh)')
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
  a_results$vanilla[[i]] <- collect_results(a_BF_standard, seed*i)
  a_results$generalised[[i]] <- collect_results(a_BF_generalised, seed*i)
  
  ##### Recommended scaling of T, fixed n #####
  print('### performing standard Bayesian Fusion (with recommended T, fixed n)')
  vanilla_guide[[i]] <- BF_guidance(condition = 'SSH',
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
  print(paste('vanilla recommened regular mesh n:', vanilla_guide[[i]]$n))
  b_mesh_vanilla <- seq(0, vanilla_guide[[i]]$min_T, length.out = 6)
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
  print('### performing Generalised Bayesian Fusion (with recommended T, fixed n)')
  gen_guide[[i]] <- BF_guidance(condition = 'SSH',
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
  print(paste('generalised recommened regular mesh n:', gen_guide[[i]]$n))
  b_mesh_gen <- seq(0, gen_guide[[i]]$min_T, length.out = 6)
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
  b_results$vanilla[[i]] <- collect_results(b_BF_standard, seed*i)
  b_results$generalised[[i]] <- collect_results(b_BF_generalised, seed*i)
  
  ##### Recommended scaling of T, regular mesh #####
  print('### performing standard Bayesian Fusion (with recommended T, regular mesh)')
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = TRUE,
                                              number_of_steps = length(vanilla_guide[[i]]$mesh))
  c_BF_standard <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                           N = nsamples,
                                           m = C,
                                           time_mesh = vanilla_guide[[i]]$mesh,
                                           mean_vecs = means,
                                           sd_vecs = rep(list(sd), C),
                                           corrs = rep(corr, C),
                                           betas = rep(beta, C),
                                           precondition_matrices = rep(list(diag(1,2)), C),
                                           ESS_threshold = ESS_threshold,
                                           diffusion_estimator = diffusion_estimator,
                                           seed = seed*i)
  print('### performing Generalised Bayesian Fusion (with recommended T, regular mesh)')
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = TRUE,
                                              number_of_steps = length(gen_guide[[i]]$mesh))
  c_BF_generalised <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                              N = nsamples,
                                              m = C,
                                              time_mesh = gen_guide[[i]]$mesh,
                                              mean_vecs = means,
                                              sd_vecs = rep(list(sd), C),
                                              corrs = rep(corr, C),
                                              betas = rep(beta, C),
                                              precondition_matrices = lapply(input_samples, cov),
                                              ESS_threshold = ESS_threshold,
                                              diffusion_estimator = diffusion_estimator,
                                              seed = seed*i)
  # save results
  c_results$vanilla[[i]] <- collect_results(c_BF_standard, seed*i)
  c_results$generalised[[i]] <- collect_results(c_BF_generalised, seed*i)

  ##### Checking regular mesh #####
  print('### performing standard Bayesian Fusion (checking regular mesh)')
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = TRUE,
                                              number_of_steps = length(vanilla_guide[[i]]$mesh))
  c_check_standard <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                              N = nsamples,
                                              m = C,
                                              time_mesh = vanilla_guide[[i]]$mesh,
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
                                                                              'k3' = vanilla_guide[[i]]$k3,
                                                                              'k4' = vanilla_guide[[i]]$k4,
                                                                              'vanilla' = TRUE),
                                              diffusion_estimator = diffusion_estimator,
                                              seed = seed*i)
  print('### performing Generalised Bayesian Fusion (checking regular mesh)')
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = TRUE,
                                              number_of_steps = length(gen_guide[[i]]$mesh))
  c_check_generalised <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                                 N = nsamples,
                                                 m = C,
                                                 time_mesh = gen_guide[[i]]$mesh,
                                                 mean_vecs = means,
                                                 sd_vecs = rep(list(sd), C),
                                                 corrs = rep(corr, C),
                                                 betas = rep(beta, C),
                                                 precondition_matrices = lapply(input_samples, cov),
                                                 ESS_threshold = ESS_threshold,
                                                 sub_posterior_means = t(sapply(input_samples, function(sub) apply(sub, 2, mean))),
                                                 adaptive_mesh = TRUE,
                                                 adaptive_mesh_parameters = list('data_size' = data_sizes[i],
                                                                                 'k3' = gen_guide[[i]]$k3,
                                                                                 'k4' = gen_guide[[i]]$k4,
                                                                                 'vanilla' = FALSE),
                                                 diffusion_estimator = diffusion_estimator,
                                                 seed = seed*i)
  # save results
  c_check_results$vanilla[[i]] <- collect_results(c_check_standard, seed*i)
  c_check_results$generalised[[i]] <- collect_results(c_check_generalised, seed*i)
  
  ##### Recommended scaling of T, adaptive mesh (equal k3,k4) #####
  print('### performing standard Bayesian Fusion (with recommended T, adaptive mesh with equal k3,k4)')
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = TRUE,
                                              number_of_steps = length(vanilla_guide[[i]]$mesh))
  d1_BF_standard <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                            N = nsamples,
                                            m = C,
                                            time_mesh = vanilla_guide[[i]]$mesh,
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
  print('### performing Generalised Bayesian Fusion (with recommended T, adaptive mesh with equal k3,k4)')
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = TRUE,
                                              number_of_steps = length(gen_guide[[i]]$mesh))
  d1_BF_generalised <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                               N = nsamples,
                                               m = C,
                                               time_mesh = gen_guide[[i]]$mesh,
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
  d1_results$vanilla[[i]] <- collect_results(d1_BF_standard, seed*i)
  d1_results$generalised[[i]] <- collect_results(d1_BF_generalised, seed*i)
  
  ##### Recommended scaling of T, adaptive mesh (un-equal k3,k4) #####
  print('### performing standard Bayesian Fusion (with recommended T, adaptive mesh)')
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = TRUE,
                                              number_of_steps = length(vanilla_guide[[i]]$mesh))
  d2_BF_standard <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                            N = nsamples,
                                            m = C,
                                            time_mesh = vanilla_guide[[i]]$mesh,
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
                                                                            'CESS_j_threshold' = CESS_j_threshold,
                                                                            'vanilla' = TRUE),
                                            diffusion_estimator = diffusion_estimator,
                                            seed = seed*i)
  print('### performing Generalised Bayesian Fusion (with recommended T, adaptive mesh)')
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = TRUE,
                                              number_of_steps = length(gen_guide[[i]]$mesh))
  d2_BF_generalised <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                               N = nsamples,
                                               m = C,
                                               time_mesh = gen_guide[[i]]$mesh,
                                               mean_vecs = means,
                                               sd_vecs = rep(list(sd), C),
                                               corrs = rep(corr, C),
                                               betas = rep(beta, C),
                                               precondition_matrices = lapply(input_samples, cov),
                                               ESS_threshold = ESS_threshold,
                                               sub_posterior_means = t(sapply(input_samples, function(sub) apply(sub, 2, mean))),
                                               adaptive_mesh = TRUE,
                                               adaptive_mesh_parameters = list('data_size' = data_sizes[i],
                                                                               'CESS_j_threshold' = CESS_j_threshold,
                                                                               'vanilla' = FALSE),
                                               diffusion_estimator = diffusion_estimator,
                                               seed = seed*i)
  # save results
  d2_results$vanilla[[i]] <- collect_results(d2_BF_standard, seed*i)
  d2_results$generalised[[i]] <- collect_results(d2_BF_generalised, seed*i)
  
  ##### SH: Recommended scaling of T, adaptive mesh #####
  print('### SH: performing standard Bayesian Fusion (with recommended T, adaptive mesh)')
  vanilla_guide_SH[[i]] <- BF_guidance(condition = 'SH',
                                       CESS_0_threshold = CESS_0_threshold,
                                       CESS_j_threshold = CESS_j_threshold,
                                       sub_posterior_samples = input_samples,
                                       C = C,
                                       d = 2,
                                       data_size = data_sizes[i],
                                       b = vanilla_b,
                                       sub_posterior_means = t(sapply(input_samples, function(sub) apply(sub, 2, mean))),
                                       k1 = k1,
                                       vanilla = TRUE)
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = TRUE,
                                              number_of_steps = length(vanilla_guide_SH[[i]]$mesh))
  SH_adaptive_standard <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                                  N = nsamples,
                                                  m = C,
                                                  time_mesh = vanilla_guide_SH[[i]]$mesh,
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
                                                                                  'CESS_j_threshold' = CESS_j_threshold,
                                                                                  'vanilla' = TRUE),
                                                  diffusion_estimator = diffusion_estimator,
                                                  seed = seed*i)
  print('### SH: performing Generalised Bayesian Fusion (with recommended T, adaptive mesh)')
  gen_guide_SH[[i]] <- BF_guidance(condition = 'SH',
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
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = TRUE,
                                              number_of_steps = length(gen_guide_SH[[i]]$mesh))
  SH_adaptive_generalised <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                                     N = nsamples,
                                                     m = C,
                                                     time_mesh = gen_guide_SH[[i]]$mesh,
                                                     mean_vecs = means,
                                                     sd_vecs = rep(list(sd), C),
                                                     corrs = rep(corr, C),
                                                     betas = rep(beta, C),
                                                     precondition_matrices = lapply(input_samples, cov),
                                                     ESS_threshold = ESS_threshold,
                                                     sub_posterior_means = t(sapply(input_samples, function(sub) apply(sub, 2, mean))),
                                                     adaptive_mesh = TRUE,
                                                     adaptive_mesh_parameters = list('data_size' = data_sizes[i],
                                                                                     'CESS_j_threshold' = CESS_j_threshold,
                                                                                     'vanilla' = FALSE),
                                                     diffusion_estimator = diffusion_estimator,
                                                     seed = seed*i)
  # save results
  SH_adaptive_results$vanilla[[i]] <- collect_results(SH_adaptive_standard, seed*i)
  SH_adaptive_results$generalised[[i]] <- collect_results(SH_adaptive_generalised, seed*i)
}

##### vanilla plots #####

##### Paper: Fixed user-specified parameters #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) a_results$vanilla[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) a_results$vanilla[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) a_results$vanilla[[i]]$CESS_j/nsamples)[[ii]]
  points(x = rep(data_sizes[ii], length(cess_j)), y = cess_j, cex = 0.5, pch = 4, lwd = 1.5)
}
axis(1, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('CESS / N', 2, 2.75, font = 2, cex = 1.5)

##### Paper: Recommended scaling of T, fixed n #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) b_results$vanilla[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) b_results$vanilla[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) b_results$vanilla[[i]]$CESS_j/nsamples)[[ii]]
  points(x = rep(data_sizes[ii], length(cess_j)), y = cess_j, cex = 0.5, pch = 4, lwd = 1.5)
}
axis(1, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('CESS / N', 2, 2.75, font = 2, cex = 1.5)

##### Paper: Recommended scaling of T, regular mesh #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) c_results$vanilla[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) c_results$vanilla[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) c_results$vanilla[[i]]$CESS_j/nsamples)[[ii]]
  points(x = rep(data_sizes[ii], length(cess_j)), y = cess_j, cex = 0.5, pch = 4, lwd = 1.5)
}
axis(1, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('CESS / N', 2, 2.75, font = 2, cex = 1.5)

##### Checking regular mesh #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) c_check_results$vanilla[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) c_check_results$vanilla[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) c_check_results$vanilla[[i]]$CESS_j/nsamples)[[ii]]
  points(x = rep(data_sizes[ii], length(cess_j)), y = cess_j, cex = 0.5, pch = 4, lwd = 1.5)
}
axis(1, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('CESS / N', 2, 2.75, font = 2, cex = 1.5)

lapply(1:length(data_sizes), function(i) c_check_results$vanilla[[i]]$chosen)
lapply(1:length(data_sizes), function(i) c_check_results$vanilla[[i]]$E_nu_j < vanilla_guide[[i]]$max_E_nu_j)

##### Recommended scaling of T, adaptive mesh (equal k3,k4) #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) d1_results$vanilla[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) d1_results$vanilla[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) d1_results$vanilla[[i]]$CESS_j/nsamples)[[ii]]
  points(x = rep(data_sizes[ii], length(cess_j)), y = cess_j, cex = 0.5, pch = 4, lwd = 1.5)
}
axis(1, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('CESS / N', 2, 2.75, font = 2, cex = 1.5)

##### Paper: Recommended scaling of T, adaptive mesh (un-equal k3,k4) #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) d2_results$vanilla[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) d2_results$vanilla[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) d2_results$vanilla[[i]]$CESS_j/nsamples)[[ii]]
  points(x = rep(data_sizes[ii], length(cess_j)), y = cess_j, cex = 0.5, pch = 4, lwd = 1.5)
}
axis(1, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('CESS / N', 2, 2.75, font = 2, cex = 1.5)

##### SH: Recommended scaling of T, adaptive mesh #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) SH_adaptive_results$vanilla[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) SH_adaptive_results$vanilla[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) SH_adaptive_results$vanilla[[i]]$CESS_j/nsamples)[[ii]]
  points(x = rep(data_sizes[ii], length(cess_j)), y = cess_j, cex = 0.5, pch = 4, lwd = 1.5)
}
axis(1, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
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
  lines(x = d1_results$vanilla[[i]]$time_mesh,
        y = rep(data_sizes[i]+25, length(d1_results$vanilla[[i]]$time_mesh)),
        type = 'b', pch = 4, lty = 2, lwd = 3)
}
max(sapply(1:length(data_sizes), function(i) max(c_results$vanilla[[i]]$time_mesh)))
axis(1, at = seq(0, 0.03, 0.01), labels = seq(0, 0.03, 0.01), font = 2, cex = 1.5)
axis(1, at = seq(0, 0.04, 0.005), labels = rep("", 9), lwd.ticks = 0.5)
mtext('Time', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(2, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 2, 2.75, font = 2, cex = 1.5)

##### Compare regular mesh and adaptive mesh times (same scale) #####
plot(x = c_results$vanilla[[1]]$time_mesh / tail(c_results$vanilla[[1]]$time_mesh, n = 1),
     y = rep(data_sizes[1]-25, length(c_results$vanilla[[1]]$time_mesh)),
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,2500),
     xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
for (i in  2:length(data_sizes)) {
  lines(x = c_results$vanilla[[i]]$time_mesh / tail(c_results$vanilla[[i]]$time_mesh, n = 1),
        y = rep(data_sizes[i]-25, length(c_results$vanilla[[i]]$time_mesh)),
        type = 'b', pch = 20, lty = 1, lwd = 3)
}
for (i in  1:length(data_sizes)) {
  lines(x = d1_results$vanilla[[i]]$time_mesh / tail(d1_results$vanilla[[i]]$time_mesh, n = 1),
        y = rep(data_sizes[i]+25, length(d1_results$vanilla[[i]]$time_mesh)),
        type = 'b', pch = 4, lty = 2, lwd = 3)
}
# highlight the points where resampling was carried out
for (i in 1:length(data_sizes)) {
  scaled_times <- (c_results$vanilla[[i]]$time_mesh / tail(c_results$vanilla[[i]]$time_mesh, n = 1))
  resampled_times <- scaled_times[c_results$vanilla[[i]]$resampled]
  points(x = resampled_times,
         y = rep(data_sizes[i]-25, length(resampled_times)),
         pch = 20, lty = 1, lwd = 3, col = 'red')
}
# highlight the points where resampling was carried out
for (i in 1:length(data_sizes)) {
  scaled_times <- (d1_results$vanilla[[i]]$time_mesh / tail(d1_results$vanilla[[i]]$time_mesh, n = 1))
  resampled_times <- scaled_times[d1_results$vanilla[[i]]$resampled]
  points(x = resampled_times,
         y = rep(data_sizes[i]+25, length(resampled_times)),
         pch = 4, lty = 1, lwd = 3, col = 'red')
}
axis(1, at = seq(0, 1, 0.1), labels = c("0.0", seq(0.1, 0.9, 0.1), "1.0"), font = 2, cex = 1.5)
axis(1, at = seq(0, 1, 0.05), labels = rep("", 21), lwd.ticks = 0.5)
mtext('Time (scaled)', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(2, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 2, 2.75, font = 2, cex = 1.5)

##### Compare regular mesh and adaptive mesh times (same scale) #####
plot(x = c_results$vanilla[[1]]$time_mesh / tail(c_results$vanilla[[1]]$time_mesh, n = 1),
     y = rep(data_sizes[1]-25, length(c_results$vanilla[[1]]$time_mesh)),
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,2500),
     xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
for (i in  2:length(data_sizes)) {
  lines(x = c_results$vanilla[[i]]$time_mesh / tail(c_results$vanilla[[i]]$time_mesh, n = 1),
        y = rep(data_sizes[i]-25, length(c_results$vanilla[[i]]$time_mesh)),
        type = 'b', pch = 20, lty = 1, lwd = 3)
}
for (i in  1:length(data_sizes)) {
  lines(x = d2_results$vanilla[[i]]$time_mesh / tail(d2_results$vanilla[[i]]$time_mesh, n = 1),
        y = rep(data_sizes[i]+25, length(d2_results$vanilla[[i]]$time_mesh)),
        type = 'b', pch = 4, lty = 2, lwd = 3)
}
# highlight the points where resampling was carried out
for (i in 1:length(data_sizes)) {
  scaled_times <- (c_results$vanilla[[i]]$time_mesh / tail(c_results$vanilla[[i]]$time_mesh, n = 1))
  resampled_times <- scaled_times[c_results$vanilla[[i]]$resampled]
  points(x = resampled_times,
         y = rep(data_sizes[i]-25, length(resampled_times)),
         pch = 20, lty = 1, lwd = 3, col = 'red')
}
# highlight the points where resampling was carried out
for (i in 1:length(data_sizes)) {
  scaled_times <- (d2_results$vanilla[[i]]$time_mesh / tail(d2_results$vanilla[[i]]$time_mesh, n = 1))
  resampled_times <- scaled_times[d2_results$vanilla[[i]]$resampled]
  points(x = resampled_times,
         y = rep(data_sizes[i]+25, length(resampled_times)),
         pch = 4, lty = 1, lwd = 3, col = 'red')
}
axis(1, at = seq(0, 1, 0.1), labels = c("0.0", seq(0.1, 0.9, 0.1), "1.0"), font = 2, cex = 1.5)
axis(1, at = seq(0, 1, 0.05), labels = rep("", 21), lwd.ticks = 0.5)
mtext('Time (scaled)', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(2, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 2, 2.75, font = 2, cex = 1.5)

##### Compare number of mesh points #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) c_results$vanilla[[i]]$n),
     type = 'b', pch = 5, lty = 3, lwd = 3, ylim = c(0,10000), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) d1_results$vanilla[[i]]$n),
      pch = 6, lty = 4, lwd = 3, type = 'b')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) d2_results$vanilla[[i]]$n),
      pch = 4, lty = 2, lwd = 3, type = 'b')
axis(1, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 10000, 2000), labels = seq(0, 10000, 2000),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 10000, 1000), labels=rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('n', 2, 2.75, font = 2, cex = 1.5)
legend(x = 250, y = 10000,
       legend = c('SSH rec. T, reg. mesh',
                  'SSH rec. T, adapt. mesh (k3=k4)',
                  'SSH rec. T, adapt. mesh'),
       lty = c(3,4,2),
       pch = c(5,6,4),
       lwd = rep(3, 6),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

##### Paper: Compare number of mesh points #####
plot(x = data_sizes,
     y = log(sapply(1:length(data_sizes), function(i) c_results$vanilla[[i]]$n), 10),
     type = 'b', pch = 5, lty = 3, lwd = 3, ylim = c(1,4), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = log(sapply(1:length(data_sizes), function(i) d2_results$vanilla[[i]]$n), 10),
      pch = 4, lty = 2, lwd = 3, type = 'b')
axis(1, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(1, 4), labels = seq(1, 4), font = 2, cex = 1.5)
axis(2, at = seq(1, 4, 0.5), labels=rep("", 7), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('log(n, 10)', 2, 2.75, font = 2, cex = 1.5)
legend(x = 250, y = 4,
       legend = c('Regular mesh', 'Adaptive mesh'),
       lty = c(3,2),
       pch = c(5,4),
       lwd = rep(3, 2),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

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
      y = sapply(1:length(data_sizes), function(i) d1_results$vanilla[[i]]$IAD),
      pch = 4, lty = 4, lwd = 3, type = 'b')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) d2_results$vanilla[[i]]$IAD),
      pch = 5, lty = 5, lwd = 3, type = 'b')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) SH_adaptive_results$vanilla[[i]]$IAD),
      pch = 6, lty = 6, lwd = 3, type = 'b')
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
                  'SSH rec. T, adapt. mesh (k3=k4)',
                  'SSH rec. T, adapt. mesh',
                  'SH rec. T, adapt. mesh'),
       lty = 1:6,
       pch = 1:6,
       lwd = rep(3, 6),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

##### time #####
plot(x = data_sizes,
     y = log(sapply(1:length(data_sizes), function(i) a_results$vanilla[[i]]$time)),
     type = 'b', pch = 1, lty = 1, lwd = 3, ylim = c(0,12), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = log(sapply(1:length(data_sizes), function(i) b_results$vanilla[[i]]$time)),
      pch = 2, lty = 2, lwd = 3, type = 'b')
lines(x = data_sizes,
      y = log(sapply(1:length(data_sizes), function(i) c_results$vanilla[[i]]$time)),
      pch = 3, lty = 3, lwd = 3, type = 'b')
lines(x = data_sizes,
      y = log(sapply(1:length(data_sizes), function(i) d1_results$vanilla[[i]]$time)),
      pch = 4, lty = 4, lwd = 3, type = 'b')
lines(x = data_sizes,
      y = log(sapply(1:length(data_sizes), function(i) d2_results$vanilla[[i]]$time)),
      pch = 5, lty = 5, lwd = 3, type = 'b')
lines(x = data_sizes,
      y = log(sapply(1:length(data_sizes), function(i) SH_adaptive_results$vanilla[[i]]$time)),
      pch = 6, lty = 6, lwd = 3, type = 'b')

axis(1, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 12, 1), labels = seq(0, 12, 1), font = 2, cex = 1.5)
mtext('log(Elapsed time in seconds)', 2, 2.75, font = 2, cex = 1.5)
legend(x = 250, y = 12,
       legend = c('Fixed T, fixed n',
                  'SSH rec. T, fixed n',
                  'SSH rec. T, reg. mesh',
                  'SSH rec. T, adapt. mesh (k3=k4)',
                  'SSH rec. T, adapt. mesh',
                  'SH rec. T, adapt. mesh'),
       lty = 1:6,
       pch = 1:6,
       lwd = rep(3, 6),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

##### generalised plots #####

##### Paper: Fixed user-specified parameters #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) a_results$generalised[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) a_results$generalised[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) a_results$generalised[[i]]$CESS_j/nsamples)[[ii]]
  points(x = rep(data_sizes[ii], length(cess_j)), y = cess_j, cex = 0.5, pch = 4, lwd = 1.5)
}
axis(1, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('CESS / N', 2, 2.75, font = 2, cex = 1.5)

##### Paper: Recommended scaling of T, fixed n #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) b_results$generalised[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) b_results$generalised[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) b_results$generalised[[i]]$CESS_j/nsamples)[[ii]]
  points(x = rep(data_sizes[ii], length(cess_j)), y = cess_j, cex = 0.5, pch = 4, lwd = 1.5)
}
axis(1, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('CESS / N', 2, 2.75, font = 2, cex = 1.5)

##### Paper: Recommended scaling of T, regular mesh #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) c_results$generalised[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) c_results$generalised[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) c_results$generalised[[i]]$CESS_j/nsamples)[[ii]]
  points(x = rep(data_sizes[ii], length(cess_j)), y = cess_j, cex = 0.5, pch = 4, lwd = 1.5)
}
axis(1, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('CESS / N', 2, 2.75, font = 2, cex = 1.5)

##### Checking regular mesh #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) c_check_results$generalised[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) c_check_results$generalised[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) c_check_results$generalised[[i]]$CESS_j/nsamples)[[ii]]
  points(x = rep(data_sizes[ii], length(cess_j)), y = cess_j, cex = 0.5, pch = 4, lwd = 1.5)
}
axis(1, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('CESS / N', 2, 2.75, font = 2, cex = 1.5)

lapply(1:length(data_sizes), function(i) c_check_results$generalised[[i]]$chosen)
lapply(1:length(data_sizes), function(i) c_check_results$generalised[[i]]$E_nu_j < gen_guide[[i]]$max_E_nu_j)

##### Recommended scaling of T, adaptive mesh (equal k3,k4) #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) d1_results$generalised[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) d1_results$generalised[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) d1_results$generalised[[i]]$CESS_j/nsamples)[[ii]]
  points(x = rep(data_sizes[ii], length(cess_j)), y = cess_j, cex = 0.5, pch = 4, lwd = 1.5)
}
axis(1, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('CESS / N', 2, 2.75, font = 2, cex = 1.5)

##### Paper: Recommended scaling of T, adaptive mesh (un-equal k3,k4) #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) d2_results$generalised[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) d2_results$generalised[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) d2_results$generalised[[i]]$CESS_j/nsamples)[[ii]]
  points(x = rep(data_sizes[ii], length(cess_j)), y = cess_j, cex = 0.5, pch = 4, lwd = 1.5)
}
axis(1, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('CESS / N', 2, 2.75, font = 2, cex = 1.5)

##### SH: Recommended scaling of T, adaptive mesh #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) SH_adaptive_results$generalised[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) SH_adaptive_results$generalised[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) SH_adaptive_results$generalised[[i]]$CESS_j/nsamples)[[ii]]
  points(x = rep(data_sizes[ii], length(cess_j)), y = cess_j, cex = 0.5, pch = 4, lwd = 1.5)
}
axis(1, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.2), labels = c("0.0", seq(0.2, 0.8, 0.2), "1.0"),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1, 0.1), labels = rep("", 11), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('CESS / N', 2, 2.75, font = 2, cex = 1.5)

##### Compare regular mesh and adaptive mesh times #####
plot(x = c_results$generalised[[1]]$time_mesh,
     y = rep(data_sizes[1]-25, length(c_results$generalised[[1]]$time_mesh)),
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,2500), xlim = c(0, 16),
     xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
for (i in  2:length(data_sizes)) {
  lines(x = c_results$generalised[[i]]$time_mesh,
        y = rep(data_sizes[i]-25, length(c_results$generalised[[i]]$time_mesh)),
        type = 'b', pch = 20, lty = 1, lwd = 3)
}
for (i in  1:length(data_sizes)) {
  lines(x = d1_results$generalised[[i]]$time_mesh,
        y = rep(data_sizes[i]+25, length(d1_results$generalised[[i]]$time_mesh)),
        type = 'b', pch = 4, lty = 2, lwd = 3)
}
max(sapply(1:length(data_sizes), function(i) max(c_results$generalised[[i]]$time_mesh)))
axis(1, at = seq(0, 16, 1), labels = seq(0, 16, 1), font = 2, cex = 1.5)
axis(1, at = seq(0, 16, 0.5), labels = rep("", 33), lwd.ticks = 0.5)
mtext('Time', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(2, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 2, 2.75, font = 2, cex = 1.5)

##### Compare regular mesh and adaptive mesh times (same scale) #####
plot(x = c_results$generalised[[1]]$time_mesh / tail(c_results$generalised[[1]]$time_mesh, n = 1),
     y = rep(data_sizes[1]-25, length(c_results$generalised[[1]]$time_mesh)),
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,2500),
     xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
for (i in  2:length(data_sizes)) {
  lines(x = c_results$generalised[[i]]$time_mesh / tail(c_results$generalised[[i]]$time_mesh, n = 1),
        y = rep(data_sizes[i]-25, length(c_results$generalised[[i]]$time_mesh)),
        type = 'b', pch = 20, lty = 1, lwd = 3)
}
for (i in  1:length(data_sizes)) {
  lines(x = d1_results$generalised[[i]]$time_mesh / tail(d1_results$generalised[[i]]$time_mesh, n = 1),
        y = rep(data_sizes[i]+25, length(d1_results$generalised[[i]]$time_mesh)),
        type = 'b', pch = 4, lty = 2, lwd = 3)
}
# highlight the points where resampling was carried out
for (i in 1:length(data_sizes)) {
  scaled_times <- (c_results$generalised[[i]]$time_mesh / tail(c_results$generalised[[i]]$time_mesh, n = 1))
  resampled_times <- scaled_times[c_results$generalised[[i]]$resampled]
  points(x = resampled_times,
         y = rep(data_sizes[i]-25, length(resampled_times)),
         pch = 20, lty = 1, lwd = 3, col = 'red')
}
# highlight the points where resampling was carried out
for (i in 1:length(data_sizes)) {
  scaled_times <- (d1_results$generalised[[i]]$time_mesh / tail(d1_results$generalised[[i]]$time_mesh, n = 1))
  resampled_times <- scaled_times[d1_results$generalised[[i]]$resampled]
  points(x = resampled_times,
         y = rep(data_sizes[i]+25, length(resampled_times)),
         pch = 4, lty = 1, lwd = 3, col = 'red')
}
axis(1, at = seq(0, 1, 0.1), labels = c("0.0", seq(0.1, 0.9, 0.1), "1.0"), font = 2, cex = 1.5)
axis(1, at = seq(0, 1, 0.05), labels = rep("", 21), lwd.ticks = 0.5)
mtext('Time (scaled)', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(2, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 2, 2.75, font = 2, cex = 1.5)

##### Compare regular mesh and adaptive mesh times (same scale) #####
plot(x = c_results$generalised[[1]]$time_mesh / tail(c_results$generalised[[1]]$time_mesh, n = 1),
     y = rep(data_sizes[1]-25, length(c_results$generalised[[1]]$time_mesh)),
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,2500),
     xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
for (i in  2:length(data_sizes)) {
  lines(x = c_results$generalised[[i]]$time_mesh / tail(c_results$generalised[[i]]$time_mesh, n = 1),
        y = rep(data_sizes[i]-25, length(c_results$generalised[[i]]$time_mesh)),
        type = 'b', pch = 20, lty = 1, lwd = 3)
}
for (i in  1:length(data_sizes)) {
  lines(x = d2_results$generalised[[i]]$time_mesh / tail(d2_results$generalised[[i]]$time_mesh, n = 1),
        y = rep(data_sizes[i]+25, length(d2_results$generalised[[i]]$time_mesh)),
        type = 'b', pch = 4, lty = 2, lwd = 3)
}
# highlight the points where resampling was carried out
for (i in 1:length(data_sizes)) {
  scaled_times <- (c_results$generalised[[i]]$time_mesh / tail(c_results$generalised[[i]]$time_mesh, n = 1))
  resampled_times <- scaled_times[c_results$generalised[[i]]$resampled]
  points(x = resampled_times,
         y = rep(data_sizes[i]-25, length(resampled_times)),
         pch = 20, lty = 1, lwd = 3, col = 'red')
}
# highlight the points where resampling was carried out
for (i in 1:length(data_sizes)) {
  scaled_times <- (d2_results$generalised[[i]]$time_mesh / tail(d2_results$generalised[[i]]$time_mesh, n = 1))
  resampled_times <- scaled_times[d2_results$generalised[[i]]$resampled]
  points(x = resampled_times,
         y = rep(data_sizes[i]+25, length(resampled_times)),
         pch = 4, lty = 1, lwd = 3, col = 'red')
}
axis(1, at = seq(0, 1, 0.1), labels = c("0.0", seq(0.1, 0.9, 0.1), "1.0"), font = 2, cex = 1.5)
axis(1, at = seq(0, 1, 0.05), labels = rep("", 21), lwd.ticks = 0.5)
mtext('Time (scaled)', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(2, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 2, 2.75, font = 2, cex = 1.5)

##### Compare number of mesh points #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) c_results$generalised[[i]]$n),
     type = 'b', pch = 3, lty = 3, lwd = 3, ylim = c(0,4000), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) d1_results$generalised[[i]]$n),
      pch = 4, lty = 4, lwd = 3, type = 'b')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) d2_results$generalised[[i]]$n),
      pch = 5, lty = 5, lwd = 3, type = 'b')
axis(1, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 10000, 1000), labels = seq(0, 10000, 1000),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 10000, 500), labels=rep("", 21), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('n', 2, 2.75, font = 2, cex = 1.5)
legend(x = 250, y = 4000,
       legend = c('SSH rec. T, reg. mesh',
                  'SSH rec. T, adapt. mesh (k3=k4)',
                  'SSH rec. T, adapt. mesh'),
       lty = 3:5,
       pch = 3:5,
       lwd = rep(3, 6),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

##### Paper: Compare number of mesh points #####
plot(x = data_sizes,
     y = log(sapply(1:length(data_sizes), function(i) c_results$generalised[[i]]$n), 10),
     type = 'b', pch = 5, lty = 3, lwd = 3, ylim = c(1,4), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = log(sapply(1:length(data_sizes), function(i) d2_results$generalised[[i]]$n), 10),
      pch = 4, lty = 2, lwd = 3, type = 'b')
axis(1, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(1, 4), labels = seq(1, 4), font = 2, cex = 1.5)
axis(2, at = seq(1, 4, 0.5), labels=rep("", 7), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('log(n, 10)', 2, 2.75, font = 2, cex = 1.5)
legend(x = 250, y = 4,
       legend = c('Regular mesh', 'Adaptive mesh'),
       lty = c(3,2),
       pch = c(5,4),
       lwd = rep(3, 2),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

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
      y = sapply(1:length(data_sizes), function(i) d1_results$generalised[[i]]$IAD),
      pch = 4, lty = 4, lwd = 3, type = 'b')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) d2_results$generalised[[i]]$IAD),
      pch = 5, lty = 5, lwd = 3, type = 'b')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) SH_adaptive_results$generalised[[i]]$IAD),
      pch = 6, lty = 6, lwd = 3, type = 'b')
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
                  'SSH rec. T, adapt. mesh (k3=k4)',
                  'SSH rec. T, adapt. mesh',
                  'SH rec. T, adapt. mesh'),
       lty = 1:6,
       pch = 1:6,
       lwd = rep(3, 6),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

##### time #####
plot(x = data_sizes,
     y = log(sapply(1:length(data_sizes), function(i) a_results$generalised[[i]]$time)),
     type = 'b', pch = 1, lty = 1, lwd = 3, ylim = c(0,12), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = log(sapply(1:length(data_sizes), function(i) b_results$generalised[[i]]$time)),
      pch = 2, lty = 2, lwd = 3, type = 'b')
lines(x = data_sizes,
      y = log(sapply(1:length(data_sizes), function(i) c_results$generalised[[i]]$time)),
      pch = 3, lty = 3, lwd = 3, type = 'b')
lines(x = data_sizes,
      y = log(sapply(1:length(data_sizes), function(i) d1_results$generalised[[i]]$time)),
      pch = 4, lty = 4, lwd = 3, type = 'b')
lines(x = data_sizes,
      y = log(sapply(1:length(data_sizes), function(i) d2_results$generalised[[i]]$time)),
      pch = 5, lty = 5, lwd = 3, type = 'b')
lines(x = data_sizes,
      y = log(sapply(1:length(data_sizes), function(i) SH_adaptive_results$generalised[[i]]$time)),
      pch = 6, lty = 6, lwd = 3, type = 'b')
axis(1, at = seq(0, 2500, 500), labels = seq(0, 2500, 500), font = 2, cex = 1.5)
axis(1, at = seq(0, 2500, 250), labels = rep("", 11), lwd.ticks = 0.5, font = 2, cex = 1.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 12, 1), labels = seq(0, 12, 1), font = 2, cex = 1.5)
mtext('log(Elapsed time in seconds)', 2, 2.75, font = 2, cex = 1.5)
legend(x = 250, y = 12,
       legend = c('Fixed T, fixed n',
                  'SSH rec. T, fixed n',
                  'SSH rec. T, reg. mesh',
                  'SSH rec. T, adapt. mesh (k3=k4)',
                  'SSH rec. T, adapt. mesh',
                  'SH rec. T, adapt. mesh'),
       lty = 1:6,
       pch = 1:6,
       lwd = rep(3, 6),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

save.image('bf_bivG_dissimilar_means_thresh_05.RData')
