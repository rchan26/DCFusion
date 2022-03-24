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
resampling_method <- 'resid'
number_of_replicates <- 10
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
                resampling_method = resampling_method,
                seed = seed*i*rep)$y_samples,
                marg_means = c(0,0),
                marg_sds = sqrt(rep(1, 2)/data_sizes[i]),
                bw = opt_bw)))
}

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
  d1_results$vanilla[[i]] <- list()
  d1_results$generalised[[i]] <- list()
  d2_results$vanilla[[i]] <- list()
  d2_results$generalised[[i]] <- list()
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
                                             resampling_method = resampling_method,
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
                                                resampling_method = resampling_method,
                                                ESS_threshold = ESS_threshold,
                                                diffusion_estimator = diffusion_estimator,
                                                seed = seed*rep*i)
    # save results
    a_results$vanilla[[i]][[rep]] <- collect_results(a_BF_standard)
    a_results$generalised[[i]][[rep]] <- collect_results(a_BF_generalised)
    
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
                                             resampling_method = resampling_method,
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
                                                resampling_method = resampling_method,
                                                ESS_threshold = ESS_threshold,
                                                diffusion_estimator = diffusion_estimator,
                                                seed = seed*rep*i)
    # save results
    b_results$vanilla[[i]][[rep]] <- collect_results(b_BF_standard)
    b_results$generalised[[i]][[rep]] <- collect_results(b_BF_generalised)
    
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
                                             resampling_method = resampling_method,
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
                                                resampling_method = resampling_method,
                                                ESS_threshold = ESS_threshold,
                                                diffusion_estimator = diffusion_estimator,
                                                seed = seed*rep*i)
    # save results
    c_results$vanilla[[i]][[rep]] <- collect_results(c_BF_standard)
    c_results$generalised[[i]][[rep]] <- collect_results(c_BF_generalised)
    
    ##### Recommended scaling of T, adaptive mesh (equal k3,k4) #####
    print('### performing standard Bayesian Fusion (with recommended T, adaptive mesh with equal k3,k4)')
    input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                                multivariate = TRUE,
                                                number_of_steps = length(vanilla_guide[[i]][[rep]]$mesh))
    d1_BF_standard <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                              N = nsamples,
                                              m = C,
                                              time_mesh = vanilla_guide[[i]][[rep]]$mesh,
                                              mean_vecs = rep(list(mean), C),
                                              sd_vecs = rep(list(sd), C),
                                              corrs = rep(corr, C),
                                              betas = rep(beta, C),
                                              precondition_matrices = rep(list(diag(1,2)), C),
                                              resampling_method = resampling_method,
                                              ESS_threshold = ESS_threshold,
                                              sub_posterior_means = t(sapply(input_samples, function(sub) apply(sub, 2, mean))),
                                              adaptive_mesh = TRUE,
                                              adaptive_mesh_parameters = list('data_size' = data_sizes[i],
                                                                              'b' = vanilla_b,
                                                                              'k3' = k3,
                                                                              'k4' = k4,
                                                                              'vanilla' = TRUE),
                                              diffusion_estimator = diffusion_estimator,
                                              seed = seed*i*rep)
    print('### performing Generalised Bayesian Fusion (with recommended T, adaptive mesh with equal k3,k4)')
    input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                                multivariate = TRUE,
                                                number_of_steps = length(gen_guide[[i]][[rep]]$mesh))
    d1_BF_generalised <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                                 N = nsamples,
                                                 m = C,
                                                 time_mesh = gen_guide[[i]][[rep]]$mesh,
                                                 mean_vecs = rep(list(mean), C),
                                                 sd_vecs = rep(list(sd), C),
                                                 corrs = rep(corr, C),
                                                 betas = rep(beta, C),
                                                 precondition_matrices = lapply(input_samples, cov),
                                                 resampling_method = resampling_method,
                                                 ESS_threshold = ESS_threshold,
                                                 sub_posterior_means = t(sapply(input_samples, function(sub) apply(sub, 2, mean))),
                                                 adaptive_mesh = TRUE,
                                                 adaptive_mesh_parameters = list('data_size' = data_sizes[i],
                                                                                 'k3' = k3,
                                                                                 'k4' = k4,
                                                                                 'vanilla' = FALSE),
                                                 diffusion_estimator = diffusion_estimator,
                                                 seed = seed*i*rep)
    # save results
    d1_results$vanilla[[i]][[rep]] <- collect_results(d1_BF_standard)
    d1_results$generalised[[i]][[rep]] <- collect_results(d1_BF_generalised)
    
    ##### Recommended scaling of T, adaptive mesh (un-equal k3,k4) #####
    print('### performing standard Bayesian Fusion (with recommended T, adaptive mesh)')
    input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                                multivariate = TRUE,
                                                number_of_steps = length(vanilla_guide[[i]][[rep]]$mesh))
    d2_BF_standard <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                              N = nsamples,
                                              m = C,
                                              time_mesh = vanilla_guide[[i]][[rep]]$mesh,
                                              mean_vecs = rep(list(mean), C),
                                              sd_vecs = rep(list(sd), C),
                                              corrs = rep(corr, C),
                                              betas = rep(beta, C),
                                              precondition_matrices = rep(list(diag(1,2)), C),
                                              resampling_method = resampling_method,
                                              ESS_threshold = ESS_threshold,
                                              sub_posterior_means = t(sapply(input_samples, function(sub) apply(sub, 2, mean))),
                                              adaptive_mesh = TRUE,
                                              adaptive_mesh_parameters = list('data_size' = data_sizes[i],
                                                                              'b' = vanilla_b,
                                                                              'threshold' = CESS_j_threshold,
                                                                              'vanilla' = TRUE),
                                              diffusion_estimator = diffusion_estimator,
                                              seed = seed*i*rep)
    print('### performing Generalised Bayesian Fusion (with recommended T, adaptive mesh)')
    input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                                multivariate = TRUE,
                                                number_of_steps = length(gen_guide[[i]][[rep]]$mesh))
    d2_BF_generalised <- parallel_GBF_biGaussian(particles_to_fuse = input_particles,
                                                 N = nsamples,
                                                 m = C,
                                                 time_mesh = gen_guide[[i]][[rep]]$mesh,
                                                 mean_vecs = rep(list(mean), C),
                                                 sd_vecs = rep(list(sd), C),
                                                 corrs = rep(corr, C),
                                                 betas = rep(beta, C),
                                                 precondition_matrices = lapply(input_samples, cov),
                                                 resampling_method = resampling_method,
                                                 ESS_threshold = ESS_threshold,
                                                 sub_posterior_means = t(sapply(input_samples, function(sub) apply(sub, 2, mean))),
                                                 adaptive_mesh = TRUE,
                                                 adaptive_mesh_parameters = list('data_size' = data_sizes[i],
                                                                                 'threshold' = CESS_j_threshold,
                                                                                 'vanilla' = FALSE),
                                                 diffusion_estimator = diffusion_estimator,
                                                 seed = seed*i*rep)
    # save results
    d2_results$vanilla[[i]][[rep]] <- collect_results(d2_BF_standard)
    d2_results$generalised[[i]][[rep]] <- collect_results(d2_BF_generalised)
    
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
                                                     resampling_method = resampling_method,
                                                     ESS_threshold = ESS_threshold,
                                                     sub_posterior_means = t(sapply(input_samples, function(sub) apply(sub, 2, mean))),
                                                     adaptive_mesh = TRUE,
                                                     adaptive_mesh_parameters = list('data_size' = data_sizes[i],
                                                                                     'b' = vanilla_b,
                                                                                     'threshold' = CESS_j_threshold,
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
                                                        resampling_method = resampling_method,
                                                        ESS_threshold = ESS_threshold,
                                                        sub_posterior_means = t(sapply(input_samples, function(sub) apply(sub, 2, mean))),
                                                        adaptive_mesh = TRUE,
                                                        adaptive_mesh_parameters = list('data_size' = data_sizes[i],
                                                                                        'threshold' = CESS_j_threshold,
                                                                                        'vanilla' = FALSE),
                                                        diffusion_estimator = diffusion_estimator,
                                                        seed = seed*rep*i)
    # save results
    SSH_adaptive_results$vanilla[[i]][[rep]] <- collect_results(SSH_adaptive_standard)
    SSH_adaptive_results$generalised[[i]][[rep]] <- collect_results(SSH_adaptive_generalised)
    
    # save progress
    print('saving progress')
    save.image('bivG_similar_means_replicates.RData')
  }
}

# ##### IAD #####
# plot(x = data_sizes,
#      y = sapply(1:length(data_sizes), function(i) {
#        mean(sapply(1:number_of_replicates, function(rep) a_results$vanilla[[i]][[rep]]$IAD))
#      }),
#      type = 'b', pch = 1, lty = 1, lwd = 3, ylim = c(0,1.6), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
# for (i in 1:length(data_sizes)) {
#   IAD <- sapply(1:number_of_replicates, function(rep) a_results$vanilla[[i]][[rep]]$IAD)
#   points(x = rep(data_sizes[i], length(IAD)), y = IAD, cex = 0.5, pch = 1)
# }
# lines(x = data_sizes,
#       y = sapply(1:length(data_sizes), function(i) {
#         mean(sapply(1:number_of_replicates, function(rep) b_results$vanilla[[i]][[rep]]$IAD))
#       }),
#       pch = 2, lty = 2, lwd = 3, type = 'b', col = 'blue')
# for (i in 1:length(data_sizes)) {
#   IAD <- sapply(1:number_of_replicates, function(rep) b_results$vanilla[[i]][[rep]]$IAD)
#   points(x = rep(data_sizes[i], length(IAD)), y = IAD, cex = 0.5, pch = 2, col = 'blue')
# }
# lines(x = data_sizes,
#       y = sapply(1:length(data_sizes), function(i) {
#         mean(sapply(1:number_of_replicates, function(rep) c_results$vanilla[[i]][[rep]]$IAD))
#       }),
#       pch = 3, lty = 3, lwd = 3, type = 'b', col = 'red')
# for (i in 1:length(data_sizes)) {
#   IAD <- sapply(1:number_of_replicates, function(rep) c_results$vanilla[[i]][[rep]]$IAD)
#   points(x = rep(data_sizes[i], length(IAD)), y = IAD, cex = 0.5, pch = 3, col = 'red')
# }
# lines(x = data_sizes,
#       y = sapply(1:length(data_sizes), function(i) {
#         mean(sapply(1:number_of_replicates, function(rep) d_results$vanilla[[i]][[rep]]$IAD))
#       }),
#       pch = 4, lty = 4, lwd = 3, type = 'b', col = 'green')
# for (i in 1:length(data_sizes)) {
#   IAD <- sapply(1:number_of_replicates, function(rep) d_results$vanilla[[i]][[rep]]$IAD)
#   points(x = rep(data_sizes[i], length(IAD)), y = IAD, cex = 0.5, pch = 4, col = 'green')
# }
# lines(x = data_sizes,
#       y = sapply(1:length(data_sizes), function(i) {
#         mean(sapply(1:number_of_replicates, function(rep) e_results$vanilla[[i]][[rep]]$IAD))
#       }),
#       pch = 5, lty = 5, lwd = 3, type = 'b', col = 'orange')
# for (i in 1:length(data_sizes)) {
#   IAD <- sapply(1:number_of_replicates, function(rep) e_results$vanilla[[i]][[rep]]$IAD)
#   points(x = rep(data_sizes[i], length(IAD)), y = IAD, cex = 0.5, pch = 5, col = 'orange')
# }
# lines(x = data_sizes,
#       y = sapply(1:length(data_sizes), function(i) {
#         mean(sapply(1:number_of_replicates, function(rep) SSH_adaptive_results$vanilla[[i]][[rep]]$IAD))
#       }),
#       pch = 6, lty = 6, lwd = 3, type = 'b', col = 'purple')
# for (i in 1:length(data_sizes)) {
#   IAD <- sapply(1:number_of_replicates, function(rep) SSH_adaptive_results$vanilla[[i]][[rep]]$IAD)
#   points(x = rep(data_sizes[i], length(IAD)), y = IAD, cex = 0.5, pch = 6, col = 'purple')
# }
# axis(1, at = c(1000, seq(10000, 50000, 10000)),
#      labels = c(1000, seq(10000, 50000, 10000)), font = 2, cex = 1.5)
# axis(1, at = seq(0, 50000, 5000), labels = rep("", 11), lwd.ticks = 0.5)
# mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
# axis(2, at = seq(0, 1.6, 0.1), labels = c("0.0", seq(0.1, 0.9, 0.1), "1.0", seq(1.1, 1.6, 0.1)),
#      font = 2, cex = 1.5)
# axis(2, at = seq(0, 1.6, 0.1), labels=rep("", 17), lwd.ticks = 0.5,
#      font = 2, cex = 1.5)
# mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
# legend(x = 1000, y = 1.6,
#        legend = c('Fixed T, fixed n',
#                   'SH rec. T, fixed n',
#                   'SH rec. T, reg. mesh',
#                   'SH rec. T, adapt. mesh',
#                   'SH rec. T, reg. mesh (same n as adapt.)',
#                   'SSH rec. T, adapt. mesh'),
#        col = c('black', 'blue', 'red', 'green', 'orange', 'purple'),
#        lty = 1:6,
#        pch = 1:6,
#        lwd = rep(3, 6),
#        cex = 1.25,
#        text.font = 2,
#        bty = 'n')
# 
# ##### IAD #####
# plot(x = data_sizes,
#      y = sapply(1:length(data_sizes), function(i) {
#        mean(sapply(1:number_of_replicates, function(rep) a_results$generalised[[i]][[rep]]$IAD))
#      }),
#      type = 'b', pch = 1, lty = 1, lwd = 3, ylim = c(0,1.6), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
# for (i in 1:length(data_sizes)) {
#   IAD <- sapply(1:number_of_replicates, function(rep) a_results$generalised[[i]][[rep]]$IAD)
#   points(x = rep(data_sizes[i], length(IAD)), y = IAD, cex = 0.5, pch = 1)
# }
# lines(x = data_sizes,
#       y = sapply(1:length(data_sizes), function(i) {
#         mean(sapply(1:number_of_replicates, function(rep) b_results$generalised[[i]][[rep]]$IAD))
#       }),
#       pch = 2, lty = 2, lwd = 3, type = 'b', col = 'blue')
# for (i in 1:length(data_sizes)) {
#   IAD <- sapply(1:number_of_replicates, function(rep) b_results$generalised[[i]][[rep]]$IAD)
#   points(x = rep(data_sizes[i], length(IAD)), y = IAD, cex = 0.5, pch = 2, col = 'blue')
# }
# lines(x = data_sizes,
#       y = sapply(1:length(data_sizes), function(i) {
#         mean(sapply(1:number_of_replicates, function(rep) c_results$generalised[[i]][[rep]]$IAD))
#       }),
#       pch = 3, lty = 3, lwd = 3, type = 'b', col = 'red')
# for (i in 1:length(data_sizes)) {
#   IAD <- sapply(1:number_of_replicates, function(rep) c_results$generalised[[i]][[rep]]$IAD)
#   points(x = rep(data_sizes[i], length(IAD)), y = IAD, cex = 0.5, pch = 3, col = 'red')
# }
# lines(x = data_sizes,
#       y = sapply(1:length(data_sizes), function(i) {
#         mean(sapply(1:number_of_replicates, function(rep) d_results$generalised[[i]][[rep]]$IAD))
#       }),
#       pch = 4, lty = 4, lwd = 3, type = 'b', col = 'green')
# for (i in 1:length(data_sizes)) {
#   IAD <- sapply(1:number_of_replicates, function(rep) d_results$generalised[[i]][[rep]]$IAD)
#   points(x = rep(data_sizes[i], length(IAD)), y = IAD, cex = 0.5, pch = 4, col = 'green')
# }
# lines(x = data_sizes,
#       y = sapply(1:length(data_sizes), function(i) {
#         mean(sapply(1:number_of_replicates, function(rep) e_results$generalised[[i]][[rep]]$IAD))
#       }),
#       pch = 5, lty = 5, lwd = 3, type = 'b', col = 'orange')
# for (i in 1:length(data_sizes)) {
#   IAD <- sapply(1:number_of_replicates, function(rep) e_results$generalised[[i]][[rep]]$IAD)
#   points(x = rep(data_sizes[i], length(IAD)), y = IAD, cex = 0.5, pch = 5, col = 'orange')
# }
# lines(x = data_sizes,
#       y = sapply(1:length(data_sizes), function(i) {
#         mean(sapply(1:number_of_replicates, function(rep) SH_adaptive_results$generalised[[i]][[rep]]$IAD))
#       }),
#       pch = 6, lty = 6, lwd = 3, type = 'b', col = 'purple')
# for (i in 1:length(data_sizes)) {
#   IAD <- sapply(1:number_of_replicates, function(rep) SH_adaptive_results$generalised[[i]][[rep]]$IAD)
#   points(x = rep(data_sizes[i], length(IAD)), y = IAD, cex = 0.5, pch = 6, col = 'purple')
# }
# axis(1, at = c(1000, seq(10000, 50000, 10000)),
#      labels = c(1000, seq(10000, 50000, 10000)), font = 2, cex = 1.5)
# axis(1, at = seq(0, 50000, 5000), labels = rep("", 11), lwd.ticks = 0.5)
# mtext('Data Sizes', 1, 2.75, font = 2, cex = 1.5)
# axis(2, at = seq(0, 1.6, 0.1), labels = c("0.0", seq(0.1, 0.9, 0.1), "1.0", seq(1.1, 1.6, 0.1)),
#      font = 2, cex = 1.5)
# axis(2, at = seq(0, 1.6, 0.1), labels=rep("", 17), lwd.ticks = 0.5,
#      font = 2, cex = 1.5)
# mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
# legend(x = 1000, y = 1.6,
#        legend = c('Fixed T, fixed n',
#                   'SH rec. T, fixed n',
#                   'SH rec. T, reg. mesh',
#                   'SH rec. T, adapt. mesh',
#                   'SH rec. T, reg. mesh (same n as adapt.)',
#                   'SSH rec. T, adapt. mesh'),
#        col = c('black', 'blue', 'red', 'green', 'orange', 'purple'),
#        lty = 1:6,
#        pch = 1:6,
#        lwd = rep(3, 6),
#        cex = 1.25,
#        text.font = 2,
#        bty = 'n')
