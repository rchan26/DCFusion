library(DCFusion)
library(MASS)

seed <- 1983
set.seed(seed)
nsamples <- 10000
fusion_time <- 1
mean <- c(1, 2)
sd <- c(sqrt(0.5), sqrt(2))
cov_mat <- matrix(c(sd[1]^2, 0.9, 0.9, sd[2]^2), nrow = 2, ncol = 2, byrow = T)
corr <- 0.9/(sd[1]*sd[2])
C <- 4
beta <- 1/C
diffusion_estimator <- 'NB'

# sampling from the sub-posteriors (target at inverse temperature 1/8)
input_samples <- lapply(1:C, function(i) mvrnormArma_tempered(1000000,
                                                              mu = mean,
                                                              Sigma = cov_mat,
                                                              beta = beta))

# sample from true target density
true_samples <- mvrnormArma(1000000, mu = mean, Sigma = cov_mat)
true_kde <- MASS::kde2d(true_samples[,1], true_samples[,2], n = 50)
image(true_kde)
contour(true_kde, add = T)

##### Monte Carlo Fusion #####
input_particles_MCF <- initialise_particle_sets(samples_to_fuse = input_samples,
                                                multivariate = TRUE)

print('### performing standard fusion')
MCF_standard <- parallel_fusion_SMC_biGaussian(particles_to_fuse = input_particles_MCF,
                                               N = nsamples,
                                               m = C,
                                               time = fusion_time,
                                               mean_vec = mean,
                                               sd_vec = sd,
                                               corr = corr,
                                               betas = rep(beta, C),
                                               precondition_matrices = rep(list(diag(1,2)), C),
                                               ESS_threshold = 0.5,
                                               diffusion_estimator = diffusion_estimator,
                                               seed = seed)
print('ESS:'); print(MCF_standard$ESS)
print('CESS:'); print(MCF_standard$CESS)

print('### performing fusion with a preconditioning matrix')
MCF_generalised <- parallel_fusion_SMC_biGaussian(particles_to_fuse = input_particles_MCF,
                                                  N = nsamples,
                                                  m = C,
                                                  time = fusion_time,
                                                  mean_vec = mean,
                                                  sd_vec = sd,
                                                  corr = corr,
                                                  betas = rep(beta, C),
                                                  precondition_matrices = lapply(input_samples, cov),
                                                  ESS_threshold = 0.5,
                                                  diffusion_estimator = diffusion_estimator,
                                                  seed = seed)
print('ESS:'); print(MCF_generalised$ESS)
print('CESS:'); print(MCF_generalised$CESS)

# plots
MCF_standard$particles <- resample_particle_y_samples(particle_set = MCF_standard$particles,
                                                      multivariate = TRUE,
                                                      resampling_method = 'resid',
                                                      seed = seed)
MCF_generalised$particles <- resample_particle_y_samples(particle_set = MCF_generalised$particles,
                                                         multivariate = TRUE,
                                                         resampling_method = 'resid',
                                                         seed = seed)
compare_samples_bivariate(list(true_samples,
                               MCF_standard$particles$y_samples,
                               MCF_generalised$particles$y_samples),
                          c('black', 'red', 'blue'),
                          c(-4, 4))

##### Bayesian Fusion (with n=1, so equal to Monte Carlo Fusion) #####
time_mesh_BF_n1 <- seq(0, fusion_time, 1)
input_particles_BF_n1 <- initialise_particle_sets(samples_to_fuse = input_samples,
                                                  multivariate = TRUE,
                                                  number_of_steps = length(time_mesh_BF_n1))

print('### performing standard Bayesian Fusion (with n=1)')
BF_standard_n1 <- parallel_GBF_biGaussian(particles_to_fuse = input_particles_BF_n1,
                                          N = nsamples,
                                          m = C,
                                          time_mesh = time_mesh_BF_n1,
                                          mean_vec = mean,
                                          sd_vec = sd,
                                          corr = corr,
                                          betas = rep(beta, C),
                                          precondition_matrices = rep(list(diag(1,2)), C),
                                          ESS_threshold = 0.5,
                                          diffusion_estimator = diffusion_estimator,
                                          seed = seed)
print('ESS:'); print(BF_standard_n1$ESS)
print('CESS:'); print(BF_standard_n1$CESS)

print('### performing Bayesian Fusion (with n=1) with a preconditioning matrix')
BF_generalised_n1 <- parallel_GBF_biGaussian(particles_to_fuse = input_particles_BF_n1,
                                             N = nsamples,
                                             m = C,
                                             time_mesh = time_mesh_BF_n1,
                                             mean_vec = mean,
                                             sd_vec = sd,
                                             corr = corr,
                                             betas = rep(beta, C),
                                             precondition_matrices = lapply(input_samples, cov),
                                             ESS_threshold = 0.5,
                                             diffusion_estimator = diffusion_estimator,
                                             seed = seed)
print('ESS:'); print(BF_generalised_n1$ESS)
print('CESS:'); print(BF_generalised_n1$CESS)

# plots
BF_standard_n1$particles <- resample_particle_y_samples(particle_set = BF_standard_n1$particles,
                                                        multivariate = TRUE,
                                                        resampling_method = 'resid',
                                                        seed = seed)
BF_generalised_n1$particles <- resample_particle_y_samples(particle_set = BF_generalised_n1$particles,
                                                           multivariate = TRUE,
                                                           resampling_method = 'resid',
                                                           seed = seed)
compare_samples_bivariate(list(true_samples,
                               BF_standard_n1$particles$y_samples,
                               BF_generalised_n1$particles$y_samples),
                          c('black', 'red', 'blue'),
                          c(-4, 4))

##### Bayesian Fusion (with n=20) #####
time_mesh_BF_n20 <- seq(0, fusion_time, fusion_time/20)
input_particles_BF_n20 <- initialise_particle_sets(samples_to_fuse = input_samples,
                                                   multivariate = TRUE,
                                                   number_of_steps = length(time_mesh_BF_n20))

print('### performing standard Bayesian Fusion (with n=20)')
BF_standard_n20 <- parallel_GBF_biGaussian(particles_to_fuse = input_particles_BF_n20,
                                           N = nsamples,
                                           m = C,
                                           time_mesh = time_mesh_BF_n20,
                                           mean_vec = mean,
                                           sd_vec = sd,
                                           corr = corr,
                                           betas = rep(beta, C),
                                           precondition_matrices = rep(list(diag(1,2)), C),
                                           ESS_threshold = 0.5,
                                           diffusion_estimator = diffusion_estimator,
                                           seed = seed)
print('ESS:'); print(BF_standard_n20$ESS)
print('CESS:'); print(BF_standard_n20$CESS)

print('### performing Bayesian Fusion (with n=20) with a preconditioning matrix')
BF_generalised_n20 <- parallel_GBF_biGaussian(particles_to_fuse = input_particles_BF_n20,
                                              N = nsamples,
                                              m = C,
                                              time_mesh = time_mesh_BF_n20,
                                              mean_vec = mean,
                                              sd_vec = sd,
                                              corr = corr,
                                              betas = rep(beta, C),
                                              precondition_matrices = lapply(input_samples, cov),
                                              ESS_threshold = 0.5,
                                              diffusion_estimator = diffusion_estimator,
                                              seed = seed)
print('ESS:'); print(BF_generalised_n20$ESS)
print('CESS:'); print(BF_generalised_n20$CESS)

# plots
# proposals
compare_samples_bivariate(list(true_samples,
                               BF_standard_n20$proposed_samples,
                               BF_generalised_n20$proposed_samples),
                          c('black', 'red', 'blue'),
                          c(-4, 4))
# resampled
BF_standard_n20$particles <- resample_particle_y_samples(particle_set = BF_standard_n20$particles,
                                                         multivariate = TRUE,
                                                         resampling_method = 'resid',
                                                         seed = seed)
BF_generalised_n20$particles <- resample_particle_y_samples(particle_set = BF_generalised_n20$particles,
                                                            multivariate = TRUE,
                                                            resampling_method = 'resid',
                                                            seed = seed)
compare_samples_bivariate(list(true_samples,
                               BF_standard_n20$particles$y_samples,
                               BF_generalised_n20$particles$y_samples),
                          c('black', 'red', 'blue'),
                          c(-4, 4))
