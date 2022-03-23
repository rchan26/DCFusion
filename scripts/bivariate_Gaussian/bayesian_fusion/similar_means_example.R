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
ESS_threshold <- 0.5
CESS_0_threshold <- 0.5
CESS_j_threshold <- 0.5
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
c_check_results <- list('vanilla' = list(), 'generalised' = list())
d1_results <- list('vanilla' = list(), 'generalised' = list())
d2_results <- list('vanilla' = list(), 'generalised' = list())
SSH_adaptive_results <- list('vanilla' = list(), 'generalised' = list())

collect_results <- function(results) {
  return(list('CESS_0' = results$CESS[1],
              'CESS_j' = results$CESS[2:length(results$CESS)],
              'CESS_j_avg' = mean(results$CESS[2:length(results$CESS)]),
              'n' = length(results$CESS),
              'time_mesh' = results$particles$time_mesh,
              'time' = results$time,
              'elapsed_time' = results$elapsed_time,
              'resampled' = results$resampled,
              'ESS' = results$ESS,
              'E_nu_j' = results$E_nu_j,
              'chosen' = results$chosen,
              'mesh_terms' = results$mesh_terms,
              'IAD' = integrated_abs_distance_biGaussian(fusion_post = resample_particle_y_samples(
                particle_set = results$particles,
                multivariate = TRUE,
                resampling_method = 'resid',
                seed = seed*i)$y_samples,
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
  input_samples <- lapply(1:C, function(sub) mvrnormArma(N = nsamples, mu = mean, Sigma = cov_mat))
  
  # ##### Fixed user-specified parameters #####
  # print('### performing standard Bayesian Fusion (with standard mesh)')
  # input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
  #                                             multivariate = TRUE,
  #                                             number_of_steps = length(a_mesh_vanilla))
  # a_BF_standard <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
  #                                          N = nsamples,
  #                                          m = C,
  #                                          time_mesh = a_mesh_vanilla,
  #                                          mean_vecs = rep(list(mean), C),
  #                                          sd_vecs = rep(list(sd), C),
  #                                          corrs = rep(corr, C),
  #                                          betas = rep(beta, C),
  #                                          precondition_matrices = rep(list(diag(1,2)), C),
  #                                          ESS_threshold = ESS_threshold,
  #                                          diffusion_estimator = diffusion_estimator,
  #                                          seed = seed*i)
  # print('### performing Generalised Bayesian Fusion (with standard mesh)')
  # input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
  #                                             multivariate = TRUE,
  #                                             number_of_steps = length(a_mesh_gen))
  # a_BF_generalised <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
  #                                             N = nsamples,
  #                                             m = C,
  #                                             time_mesh = a_mesh_gen,
  #                                             mean_vecs = rep(list(mean), C),
  #                                             sd_vecs = rep(list(sd), C),
  #                                             corrs = rep(corr, C),
  #                                             betas = rep(beta, C),
  #                                             precondition_matrices = lapply(input_samples, cov),
  #                                             ESS_threshold = ESS_threshold,
  #                                             diffusion_estimator = diffusion_estimator,
  #                                             seed = seed*i)
  # # save results
  # a_results$vanilla[[i]] <- collect_results(a_BF_standard)
  # a_results$generalised[[i]] <- collect_results(a_BF_generalised)
  
  ##### Recommended scaling of T, fixed n #####
  print('### performing standard Bayesian Fusion (with recommended T, fixed n)')
  vanilla_guide[[i]] <- BF_guidance(condition = 'SH',
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
  print(paste('vanilla recommened regular mesh n:', vanilla_guide[[i]]$n))
  b_mesh_vanilla <- seq(0, vanilla_guide[[i]]$min_T, length.out = 6)
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
                                           seed = seed*i)
  print('### performing Generalised Bayesian Fusion (with recommended T, fixed n)')
  gen_guide[[i]] <- BF_guidance(condition = 'SH',
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
  print(paste('generalised recommened regular mesh n:', gen_guide[[i]]$n))
  b_mesh_gen <- seq(0, gen_guide[[i]]$min_T, length.out = 6)
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
                                              seed = seed*i)
  # save results
  b_results$vanilla[[i]] <- collect_results(b_BF_standard)
  b_results$generalised[[i]] <- collect_results(b_BF_generalised)
  
  ##### Recommended scaling of T, regular mesh #####
  print('### performing standard Bayesian Fusion (with recommended T, regular mesh)')
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = TRUE,
                                              number_of_steps = length(vanilla_guide[[i]]$mesh))
  c_BF_standard <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                           N = nsamples,
                                           m = C,
                                           time_mesh = vanilla_guide[[i]]$mesh,
                                           mean_vecs = rep(list(mean), C),
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
                                              mean_vecs = rep(list(mean), C),
                                              sd_vecs = rep(list(sd), C),
                                              corrs = rep(corr, C),
                                              betas = rep(beta, C),
                                              precondition_matrices = lapply(input_samples, cov),
                                              ESS_threshold = ESS_threshold,
                                              diffusion_estimator = diffusion_estimator,
                                              seed = seed*i)
  # save results
  c_results$vanilla[[i]] <- collect_results(c_BF_standard)
  c_results$generalised[[i]] <- collect_results(c_BF_generalised)

  ##### Checking regular mesh #####
  print('### performing standard Bayesian Fusion (checking regular mesh)')
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = TRUE,
                                              number_of_steps = length(vanilla_guide[[i]]$mesh))
  c_check_standard <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                              N = nsamples,
                                              m = C,
                                              time_mesh = vanilla_guide[[i]]$mesh,
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
                                                 mean_vecs = rep(list(mean), C),
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
  c_check_results$vanilla[[i]] <- collect_results(c_check_standard)
  c_check_results$generalised[[i]] <- collect_results(c_check_generalised)
  
  ##### Recommended scaling of T, adaptive mesh (equal k3,k4) #####
  print('### performing standard Bayesian Fusion (with recommended T, adaptive mesh with equal k3,k4)')
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = TRUE,
                                              number_of_steps = length(vanilla_guide[[i]]$mesh))
  d1_BF_standard <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                            N = nsamples,
                                            m = C,
                                            time_mesh = vanilla_guide[[i]]$mesh,
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
                                            seed = seed*i)
  print('### performing Generalised Bayesian Fusion (with recommended T, adaptive mesh with equal k3,k4)')
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = TRUE,
                                              number_of_steps = length(gen_guide[[i]]$mesh))
  d1_BF_generalised <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                               N = nsamples,
                                               m = C,
                                               time_mesh = gen_guide[[i]]$mesh,
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
                                               seed = seed*i)
  # save results
  d1_results$vanilla[[i]] <- collect_results(d1_BF_standard)
  d1_results$generalised[[i]] <- collect_results(d1_BF_generalised)
  
  ##### Recommended scaling of T, adaptive mesh (un-equal k3,k4) #####
  print('### performing standard Bayesian Fusion (with recommended T, adaptive mesh)')
  input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                              multivariate = TRUE,
                                              number_of_steps = length(vanilla_guide[[i]]$mesh))
  d2_BF_standard <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                            N = nsamples,
                                            m = C,
                                            time_mesh = vanilla_guide[[i]]$mesh,
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
                                                                            'threshold' = CESS_j_threshold,
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
                                               mean_vecs = rep(list(mean), C),
                                               sd_vecs = rep(list(sd), C),
                                               corrs = rep(corr, C),
                                               betas = rep(beta, C),
                                               precondition_matrices = lapply(input_samples, cov),
                                               ESS_threshold = ESS_threshold,
                                               sub_posterior_means = t(sapply(input_samples, function(sub) apply(sub, 2, mean))),
                                               adaptive_mesh = TRUE,
                                               adaptive_mesh_parameters = list('data_size' = data_sizes[i],
                                                                               'threshold' = CESS_j_threshold,
                                                                               'vanilla' = FALSE),
                                               diffusion_estimator = diffusion_estimator,
                                               seed = seed*i)
  # save results
  d2_results$vanilla[[i]] <- collect_results(d2_BF_standard)
  d2_results$generalised[[i]] <- collect_results(d2_BF_generalised)
  
  ##### SSH: Recommended scaling of T, adaptive mesh #####
  print('### SSH: performing standard Bayesian Fusion (with recommended T, adaptive mesh)')
  vanilla_guide_SSH[[i]] <- BF_guidance(condition = 'SSH',
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
                                              number_of_steps = length(vanilla_guide_SSH[[i]]$mesh))
  SSH_adaptive_standard <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                                   N = nsamples,
                                                   m = C,
                                                   time_mesh = vanilla_guide_SSH[[i]]$mesh,
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
                                                                                   'threshold' = CESS_j_threshold,
                                                                                   'vanilla' = TRUE),
                                                   diffusion_estimator = diffusion_estimator,
                                                   seed = seed*i)
  print('### SSH: performing Generalised Bayesian Fusion (with recommended T, adaptive mesh)')
  gen_guide_SSH[[i]] <- BF_guidance(condition = 'SSH',
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
                                              number_of_steps = length(gen_guide_SSH[[i]]$mesh))
  SSH_adaptive_generalised <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                                      N = nsamples,
                                                      m = C,
                                                      time_mesh = gen_guide_SSH[[i]]$mesh,
                                                      mean_vecs = rep(list(mean), C),
                                                      sd_vecs = rep(list(sd), C),
                                                      corrs = rep(corr, C),
                                                      betas = rep(beta, C),
                                                      precondition_matrices = lapply(input_samples, cov),
                                                      ESS_threshold = ESS_threshold,
                                                      sub_posterior_means = t(sapply(input_samples, function(sub) apply(sub, 2, mean))),
                                                      adaptive_mesh = TRUE,
                                                      adaptive_mesh_parameters = list('data_size' = data_sizes[i],
                                                                                      'threshold' = CESS_j_threshold,
                                                                                      'vanilla' = FALSE),
                                                      diffusion_estimator = diffusion_estimator,
                                                      seed = seed*i)
  # save results
  SSH_adaptive_results$vanilla[[i]] <- collect_results(SSH_adaptive_standard)
  SSH_adaptive_results$generalised[[i]] <- collect_results(SSH_adaptive_generalised)
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

##### Checking regular mesh #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) c_check_results$vanilla[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) c_check_results$vanilla[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) c_check_results$vanilla[[i]]$CESS_j/nsamples)[[ii]]
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

##### Recommended scaling of T, adaptive mesh (un-equal k3,k4) #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) d2_results$vanilla[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) d2_results$vanilla[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) d2_results$vanilla[[i]]$CESS_j/nsamples)[[ii]]
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

##### Compare regular mesh and adaptive mesh times #####
plot(x = c_results$vanilla[[1]]$time_mesh,
     y = rep(data_sizes[1]-500, length(c_results$vanilla[[1]]$time_mesh)),
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,40000),
     xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
for (i in  2:length(data_sizes)) {
  lines(x = c_results$vanilla[[i]]$time_mesh,
        y = rep(data_sizes[i]-500, length(c_results$vanilla[[i]]$time_mesh)),
        type = 'b', pch = 20, lty = 1, lwd = 3)
}
for (i in  1:length(data_sizes)) {
  lines(x = d1_results$vanilla[[i]]$time_mesh,
        y = rep(data_sizes[i]+500, length(d1_results$vanilla[[i]]$time_mesh)),
        type = 'b', pch = 4, lty = 2, lwd = 3)
}
# max(sapply(1:length(data_sizes), function(i) c_results$vanilla[[i]]$time_mesh))
axis(1, at = seq(0, 0.03, 0.01), labels = seq(0, 0.03, 0.01), font = 2, cex = 1.5)
axis(1, at = seq(0, 0.03, 0.005), labels = rep("", 7), lwd.ticks = 0.5)
mtext('Time', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = c(1000, seq(10000, 50000, 10000)),
     labels = c(1000, seq(10000, 50000, 10000)), font = 2, cex = 1.5)
axis(2, at = seq(0, 50000, 5000), labels = rep("", 11), lwd.ticks = 0.5)
mtext('Data Sizes', 2, 2.75, font = 2, cex = 1.5)

##### Compare regular mesh and adaptive mesh times (same scale) #####
plot(x = c_results$vanilla[[1]]$time_mesh / tail(c_results$vanilla[[1]]$time_mesh, n = 1),
     y = rep(data_sizes[1]-500, length(c_results$vanilla[[1]]$time_mesh)),
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,40000),
     xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
for (i in  2:length(data_sizes)) {
  lines(x = c_results$vanilla[[i]]$time_mesh / tail(c_results$vanilla[[i]]$time_mesh, n = 1),
        y = rep(data_sizes[i]-500, length(c_results$vanilla[[i]]$time_mesh)),
        type = 'b', pch = 20, lty = 1, lwd = 3)
}
for (i in  1:length(data_sizes)) {
  lines(x = d1_results$vanilla[[i]]$time_mesh / tail(d1_results$vanilla[[i]]$time_mesh, n = 1),
        y = rep(data_sizes[i]+500, length(d1_results$vanilla[[i]]$time_mesh)),
        type = 'b', pch = 4, lty = 2, lwd = 3)
}
# highlight the points where resampling was carried out
for (i in 1:length(data_sizes)) {
  scaled_times <- (c_results$vanilla[[i]]$time_mesh / tail(c_results$vanilla[[i]]$time_mesh, n = 1))
  resampled_times <- scaled_times[c_results$vanilla[[i]]$resampled]
  points(x = resampled_times,
         y = rep(data_sizes[i]-500, length(resampled_times)),
         pch = 20, lty = 1, lwd = 3, col = 'red')
}
# highlight the points where resampling was carried out
for (i in 1:length(data_sizes)) {
  scaled_times <- (d1_results$vanilla[[i]]$time_mesh / tail(d1_results$vanilla[[i]]$time_mesh, n = 1))
  resampled_times <- scaled_times[d1_results$vanilla[[i]]$resampled]
  points(x = resampled_times,
         y = rep(data_sizes[i]+500, length(resampled_times)),
         pch = 4, lty = 1, lwd = 3, col = 'red')
}
axis(1, at = seq(0, 1, 0.1), labels = c("0.0", seq(0.1, 0.9, 0.1), "1.0"), font = 2, cex = 1.5)
axis(1, at = seq(0, 1, 0.05), labels = rep("", 21), lwd.ticks = 0.5)
mtext('Time (scaled)', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = c(1000, seq(10000, 50000, 10000)),
     labels = c(1000, seq(10000, 50000, 10000)), font = 2, cex = 1.5)
axis(2, at = seq(0, 50000, 5000), labels = rep("", 11), lwd.ticks = 0.5)
mtext('Data Sizes', 2, 2.75, font = 2, cex = 1.5)

##### Compare regular mesh and adaptive mesh times (same scale) #####
plot(x = c_results$vanilla[[1]]$time_mesh / tail(c_results$vanilla[[1]]$time_mesh, n = 1),
     y = rep(data_sizes[1]-500, length(c_results$vanilla[[1]]$time_mesh)),
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,40000),
     xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
for (i in  2:length(data_sizes)) {
  lines(x = c_results$vanilla[[i]]$time_mesh / tail(c_results$vanilla[[i]]$time_mesh, n = 1),
        y = rep(data_sizes[i]-500, length(c_results$vanilla[[i]]$time_mesh)),
        type = 'b', pch = 20, lty = 1, lwd = 3)
}
for (i in  1:length(data_sizes)) {
  lines(x = d2_results$vanilla[[i]]$time_mesh / tail(d2_results$vanilla[[i]]$time_mesh, n = 1),
        y = rep(data_sizes[i]+500, length(d2_results$vanilla[[i]]$time_mesh)),
        type = 'b', pch = 4, lty = 2, lwd = 3)
}
# highlight the points where resampling was carried out
for (i in 1:length(data_sizes)) {
  scaled_times <- (c_results$vanilla[[i]]$time_mesh / tail(c_results$vanilla[[i]]$time_mesh, n = 1))
  resampled_times <- scaled_times[c_results$vanilla[[i]]$resampled]
  points(x = resampled_times,
         y = rep(data_sizes[i]-500, length(resampled_times)),
         pch = 20, lty = 1, lwd = 3, col = 'red')
}
# highlight the points where resampling was carried out
for (i in 1:length(data_sizes)) {
  scaled_times <- (d2_results$vanilla[[i]]$time_mesh / tail(d2_results$vanilla[[i]]$time_mesh, n = 1))
  resampled_times <- scaled_times[d2_results$vanilla[[i]]$resampled]
  points(x = resampled_times,
         y = rep(data_sizes[i]+500, length(resampled_times)),
         pch = 4, lty = 1, lwd = 3, col = 'red')
}
axis(1, at = seq(0, 1, 0.1), labels = c("0.0", seq(0.1, 0.9, 0.1), "1.0"), font = 2, cex = 1.5)
axis(1, at = seq(0, 1, 0.05), labels = rep("", 21), lwd.ticks = 0.5)
mtext('Time (scaled)', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = c(1000, seq(10000, 50000, 10000)),
     labels = c(1000, seq(10000, 50000, 10000)), font = 2, cex = 1.5)
axis(2, at = seq(0, 50000, 5000), labels = rep("", 11), lwd.ticks = 0.5)
mtext('Data Sizes', 2, 2.75, font = 2, cex = 1.5)

##### IAD #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) a_results$vanilla[[i]]$IAD),
     type = 'b', pch = 1, lty = 1, lwd = 3, ylim = c(0,1.6), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
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
      y = sapply(1:length(data_sizes), function(i) SSH_adaptive_results$vanilla[[i]]$IAD),
      pch = 6, lty = 6, lwd = 3, type = 'b')
axis(1, at = c(1000, seq(10000, 50000, 10000)),
     labels = c(1000, seq(10000, 50000, 10000)), font = 2, cex = 1.5)
axis(1, at = seq(0, 50000, 5000), labels = rep("", 11), lwd.ticks = 0.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1.6, 0.1), labels = c("0.0", seq(0.1, 0.9, 0.1), "1.0", seq(1.1, 1.6, 0.1)),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1.6, 0.1), labels=rep("", 16), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
legend(x = 1000, y = 1.6,
       legend = c('Fixed T, fixed n',
                  'SH rec. T, fixed n',
                  'SH rec. T, reg. mesh',
                  'SH rec. T, adapt. mesh (k3=k4)',
                  'SH rec. T, adapt. mesh',
                  'SSH rec. T, adapt. mesh'),
       lty = 1:6,
       pch = 1:6,
       lwd = rep(3, 6),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

##### time #####
plot(x = data_sizes,
     y = log(sapply(1:length(data_sizes), function(i) a_results$vanilla[[i]]$time)),
     type = 'b', pch = 1, lty = 1, lwd = 3, ylim = c(1,13), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
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
      y = log(sapply(1:length(data_sizes), function(i) SSH_adaptive_results$vanilla[[i]]$time)),
      pch = 6, lty = 6, lwd = 3, type = 'b')
axis(1, at = c(1000, seq(10000, 50000, 10000)),
     labels = c(1000, seq(10000, 50000, 10000)), font = 2, cex = 1.5)
axis(1, at = seq(0, 50000, 5000), labels = rep("", 11), lwd.ticks = 0.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 13, 1), labels = seq(0, 13, 1), font = 2, cex = 1.5)
mtext('log(Elapsed time in seconds)', 2, 2.75, font = 2, cex = 1.5)
legend(x = 1000, y = 13,
       legend = c('Fixed T, fixed n',
                  'SH rec. T, fixed n',
                  'SH rec. T, reg. mesh',
                  'SH rec. T, adapt. mesh (k3=k4)',
                  'SH rec. T, adapt. mesh',
                  'SSH rec. T, adapt. mesh'),
       lty = 1:6,
       pch = 1:6,
       lwd = rep(3, 6),
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

##### Checking regular mesh #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) c_check_results$generalised[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) c_check_results$generalised[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) c_check_results$generalised[[i]]$CESS_j/nsamples)[[ii]]
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

##### Recommended scaling of T, adaptive mesh (un-equal k3,k4) #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) d2_results$generalised[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) d2_results$generalised[[i]]$CESS_j_avg)/nsamples,
      pch = 20, lty = 2, lwd = 3, type = 'b')
for (ii in 1:length(data_sizes)) {
  cess_j <- lapply(1:length(data_sizes), function(i) d2_results$generalised[[i]]$CESS_j/nsamples)[[ii]]
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

##### SSH: Recommended scaling of T, adaptive mesh #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) SSH_adaptive_results$generalised[[i]]$CESS_0)/nsamples,
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = data_sizes,
      y = sapply(1:length(data_sizes), function(i) SSH_adaptive_results$generalised[[i]]$CESS_j_avg)/nsamples,
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

##### Compare regular mesh and adaptive mesh times #####
plot(x = c_results$generalised[[1]]$time_mesh,
     y = rep(data_sizes[1]-500, length(c_results$generalised[[1]]$time_mesh)),
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,40000),
     xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
for (i in  2:length(data_sizes)) {
  lines(x = c_results$generalised[[i]]$time_mesh,
        y = rep(data_sizes[i]-500, length(c_results$generalised[[i]]$time_mesh)),
        type = 'b', pch = 20, lty = 1, lwd = 3)
}
for (i in  1:length(data_sizes)) {
  lines(x = d1_results$generalised[[i]]$time_mesh,
        y = rep(data_sizes[i]+500, length(d1_results$generalised[[i]]$time_mesh)),
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

##### Compare regular mesh and adaptive mesh times (same scale) #####
plot(x = c_results$generalised[[1]]$time_mesh / tail(c_results$generalised[[1]]$time_mesh, n = 1),
     y = rep(data_sizes[1]-500, length(c_results$generalised[[1]]$time_mesh)),
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,40000),
     xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
for (i in  2:length(data_sizes)) {
  lines(x = c_results$generalised[[i]]$time_mesh / tail(c_results$generalised[[i]]$time_mesh, n = 1),
        y = rep(data_sizes[i]-500, length(c_results$generalised[[i]]$time_mesh)),
        type = 'b', pch = 20, lty = 1, lwd = 3)
}
for (i in  1:length(data_sizes)) {
  lines(x = d1_results$generalised[[i]]$time_mesh / tail(d1_results$generalised[[i]]$time_mesh, n = 1),
        y = rep(data_sizes[i]+500, length(d1_results$generalised[[i]]$time_mesh)),
        type = 'b', pch = 4, lty = 2, lwd = 3)
}
# highlight the points where resampling was carried out
for (i in 1:length(data_sizes)) {
  scaled_times <- (c_results$generalised[[i]]$time_mesh / tail(c_results$generalised[[i]]$time_mesh, n = 1))
  resampled_times <- scaled_times[c_results$generalised[[i]]$resampled]
  points(x = resampled_times,
         y = rep(data_sizes[i]-500, length(resampled_times)),
         pch = 20, lty = 1, lwd = 3, col = 'red')
}
# highlight the points where resampling was carried out
for (i in 1:length(data_sizes)) {
  scaled_times <- (d1_results$generalised[[i]]$time_mesh / tail(d1_results$generalised[[i]]$time_mesh, n = 1))
  resampled_times <- scaled_times[d1_results$generalised[[i]]$resampled]
  points(x = resampled_times,
         y = rep(data_sizes[i]+500, length(resampled_times)),
         pch = 4, lty = 1, lwd = 3, col = 'red')
}
axis(1, at = seq(0, 1, 0.1), labels = c("0.0", seq(0.1, 0.9, 0.1), "1.0"), font = 2, cex = 1.5)
axis(1, at = seq(0, 1, 0.05), labels = rep("", 21), lwd.ticks = 0.5)
mtext('Time (scaled)', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = c(1000, seq(10000, 50000, 10000)),
     labels = c(1000, seq(10000, 50000, 10000)), font = 2, cex = 1.5)
axis(2, at = seq(0, 50000, 5000), labels = rep("", 11), lwd.ticks = 0.5)
mtext('Data Sizes', 2, 2.75, font = 2, cex = 1.5)

##### Compare regular mesh and adaptive mesh times (same scale) #####
plot(x = c_results$generalised[[1]]$time_mesh / tail(c_results$generalised[[1]]$time_mesh, n = 1),
     y = rep(data_sizes[1]-500, length(c_results$generalised[[1]]$time_mesh)),
     type = 'b', pch = 20, lty = 1, lwd = 3, ylim = c(0,40000),
     xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
for (i in  2:length(data_sizes)) {
  lines(x = c_results$generalised[[i]]$time_mesh / tail(c_results$generalised[[i]]$time_mesh, n = 1),
        y = rep(data_sizes[i]-500, length(c_results$generalised[[i]]$time_mesh)),
        type = 'b', pch = 20, lty = 1, lwd = 3)
}
for (i in  1:length(data_sizes)) {
  lines(x = d2_results$generalised[[i]]$time_mesh / tail(d2_results$generalised[[i]]$time_mesh, n = 1),
        y = rep(data_sizes[i]+500, length(d2_results$generalised[[i]]$time_mesh)),
        type = 'b', pch = 4, lty = 2, lwd = 3)
}
# highlight the points where resampling was carried out
for (i in 1:length(data_sizes)) {
  scaled_times <- (c_results$generalised[[i]]$time_mesh / tail(c_results$generalised[[i]]$time_mesh, n = 1))
  resampled_times <- scaled_times[c_results$generalised[[i]]$resampled]
  points(x = resampled_times,
         y = rep(data_sizes[i]-500, length(resampled_times)),
         pch = 20, lty = 1, lwd = 3, col = 'red')
}
# highlight the points where resampling was carried out
for (i in 1:length(data_sizes)) {
  scaled_times <- (d2_results$generalised[[i]]$time_mesh / tail(d2_results$generalised[[i]]$time_mesh, n = 1))
  resampled_times <- scaled_times[d2_results$generalised[[i]]$resampled]
  points(x = resampled_times,
         y = rep(data_sizes[i]+500, length(resampled_times)),
         pch = 4, lty = 1, lwd = 3, col = 'red')
}
axis(1, at = seq(0, 1, 0.1), labels = c("0.0", seq(0.1, 0.9, 0.1), "1.0"), font = 2, cex = 1.5)
axis(1, at = seq(0, 1, 0.05), labels = rep("", 21), lwd.ticks = 0.5)
mtext('Time (scaled)', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = c(1000, seq(10000, 50000, 10000)),
     labels = c(1000, seq(10000, 50000, 10000)), font = 2, cex = 1.5)
axis(2, at = seq(0, 50000, 5000), labels = rep("", 11), lwd.ticks = 0.5)
mtext('Data Sizes', 2, 2.75, font = 2, cex = 1.5)

##### IAD #####
plot(x = data_sizes,
     y = sapply(1:length(data_sizes), function(i) a_results$generalised[[i]]$IAD),
     type = 'b', pch = 1, lty = 1, lwd = 3, ylim = c(0,1.6), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
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
      y = sapply(1:length(data_sizes), function(i) SSH_adaptive_results$generalised[[i]]$IAD),
      pch = 6, lty = 6, lwd = 3, type = 'b')
axis(1, at = c(1000, seq(10000, 50000, 10000)),
     labels = c(1000, seq(10000, 50000, 10000)), font = 2, cex = 1.5)
axis(1, at = seq(0, 50000, 5000), labels = rep("", 11), lwd.ticks = 0.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1.6, 0.1), labels = c("0.0", seq(0.1, 0.9, 0.1), "1.0", seq(1.1, 1.6, 0.1)),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1.6, 0.1), labels=rep("", 16), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
legend(x = 1000, y = 1.6,
       legend = c('Fixed T, fixed n',
                  'SH rec. T, fixed n',
                  'SH rec. T, reg. mesh',
                  'SH rec. T, adapt. mesh (k3=k4)',
                  'SH rec. T, adapt. mesh',
                  'SSH rec. T, adapt. mesh'),
       lty = 1:6,
       pch = 1:6,
       lwd = rep(3, 6),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

##### time #####
plot(x = data_sizes,
     y = log(sapply(1:length(data_sizes), function(i) a_results$generalised[[i]]$time)),
     type = 'b', pch = 1, lty = 1, lwd = 3, ylim = c(1,10), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
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
      y = log(sapply(1:length(data_sizes), function(i) SSH_adaptive_results$generalised[[i]]$time)),
      pch = 6, lty = 6, lwd = 3, type = 'b')
axis(1, at = c(1000, seq(10000, 50000, 10000)),
     labels = c(1000, seq(10000, 50000, 10000)), font = 2, cex = 1.5)
axis(1, at = seq(0, 50000, 5000), labels = rep("", 11), lwd.ticks = 0.5)
mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 10, 1), labels = seq(0, 10, 1), font = 2, cex = 1.5)
mtext('log(Elapsed time in seconds)', 2, 2.75, font = 2, cex = 1.5)
legend(x = 1000, y = 10,
       legend = c('Fixed T, fixed n',
                  'SH rec. T, fixed n',
                  'SH rec. T, reg. mesh',
                  'SH rec. T, adapt. mesh (k3=k4)',
                  'SH rec. T, adapt. mesh',
                  'SSH rec. T, adapt. mesh'),
       lty = 1:6,
       pch = 1:6,
       lwd = rep(3, 6),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

save.image('bf_bivG_similar_means_thresh_05.RData')
