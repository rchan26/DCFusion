library(DCFusion)

seed <- 1994
set.seed(seed)
nsamples <- 10000
fusion_time <- 1
mean <- 1.2
sd <- sqrt(0.1)
C <- 5
beta <- 1/C
diffusion_estimator <- 'NB'

input_samples <- lapply(1:C, function(i) rnorm_tempered(N = nsamples,
                                                        mean = mean,
                                                        sd = sd,
                                                        beta = beta))

##### Monte Carlo Fusion #####
input_particles_MCF <- initialise_particle_sets(samples_to_fuse = input_samples,
                                                multivariate = FALSE)

print('### performing standard fusion')
MCF_standard <- parallel_fusion_SMC_uniGaussian(particles_to_fuse = input_particles_MCF,
                                                N = nsamples,
                                                m = C,
                                                time = fusion_time,
                                                means = rep(mean, C),
                                                sds = rep(sd, C),
                                                betas = rep(beta, C),
                                                precondition_values = rep(1, C),
                                                ESS_threshold = 0.5,
                                                diffusion_estimator = diffusion_estimator,
                                                seed = seed)
print('ESS:'); print(MCF_standard$ESS)
print('CESS:'); print(MCF_standard$CESS)

print('### performing fusion with a preconditioning matrix')
MCF_generalised <- parallel_fusion_SMC_uniGaussian(particles_to_fuse = input_particles_MCF,
                                                   N = nsamples,
                                                   m = C,
                                                   time = fusion_time,
                                                   means = rep(mean, C),
                                                   sds = rep(sd, C),
                                                   betas = rep(beta, C),
                                                   precondition_values = sapply(input_samples, var),
                                                   ESS_threshold = 0.5,
                                                   diffusion_estimator = diffusion_estimator,
                                                   seed = seed)
print('ESS:'); print(MCF_generalised$ESS)
print('CESS:'); print(MCF_generalised$CESS)

# plots
curve(dnorm(x, mean = 1.2, sd = sqrt(0.1)), -3, 5)
lines(density(resample_particle_y_samples(particle_set = MCF_generalised$particles, 
                                          multivariate = FALSE,
                                          resampling_method = 'resid',
                                          seed = seed)$y_samples), 
      col = 'red')
lines(density(resample_particle_y_samples(particle_set = MCF_standard$particles, 
                                          multivariate = FALSE,
                                          resampling_method = 'resid',
                                          seed = seed)$y_samples), 
      col = 'blue')

##### Bayesian Fusion (with n=1, so equal to Monte Carlo Fusion) #####
time_mesh_BF_n1 <- seq(0, fusion_time, 1)
input_particles_BF_n1 <- initialise_particle_sets(samples_to_fuse = input_samples,
                                                  multivariate = FALSE,
                                                  number_of_steps = length(time_mesh_BF_n1))

print('### performing standard Bayesian Fusion (with n=1)')
BF_standard_n1 <- parallel_GBF_uniGaussian(particles_to_fuse = input_particles_BF_n1,
                                           N = nsamples,
                                           m = C,
                                           time_mesh = time_mesh_BF_n1,
                                           means = rep(mean, C),
                                           sds = rep(sd, C),
                                           betas = rep(beta, C),
                                           precondition_values = rep(1, C),
                                           ESS_threshold = 0.5,
                                           diffusion_estimator = diffusion_estimator,
                                           seed = seed)
print('ESS:'); print(BF_standard_n1$ESS)
print('CESS:'); print(BF_standard_n1$CESS)

print('### performing Bayesian Fusion (with n=1) with a preconditioning matrix')
BF_generalised_n1 <- parallel_GBF_uniGaussian(particles_to_fuse = input_particles_BF_n1,
                                              N = nsamples,
                                              m = C,
                                              time_mesh = time_mesh_BF_n1,
                                              means = rep(mean, C),
                                              sds = rep(sd, C),
                                              betas = rep(beta, C),
                                              precondition_values = sapply(input_samples, var),
                                              ESS_threshold = 0.5,
                                              diffusion_estimator = diffusion_estimator,
                                              seed = seed)
print('ESS:'); print(BF_generalised_n1$ESS)
print('CESS:'); print(BF_generalised_n1$CESS)

# plots
curve(dnorm(x, mean = 1.2, sd = sqrt(0.1)), -3, 5)
lines(density(resample_particle_y_samples(particle_set = BF_standard_n1$particles, 
                                          multivariate = FALSE,
                                          resampling_method = 'resid',
                                          seed = seed)$y_samples), 
      col = 'red')
lines(density(resample_particle_y_samples(particle_set = BF_generalised_n1$particles, 
                                          multivariate = FALSE,
                                          resampling_method = 'resid',
                                          seed = seed)$y_samples), 
      col = 'blue')

##### Bayesian Fusion (with n=20) #####
time_mesh_BF_n20 <- seq(0, fusion_time, fusion_time/20)
input_particles_BF_n20 <- initialise_particle_sets(samples_to_fuse = input_samples,
                                                   multivariate = FALSE,
                                                   number_of_steps = length(time_mesh_BF_n20))

print('### performing standard Bayesian Fusion (with n=20)')
BF_standard_n20 <- parallel_GBF_uniGaussian(particles_to_fuse = input_particles_BF_n20,
                                            N = nsamples,
                                            m = C,
                                            time_mesh = time_mesh_BF_n20,
                                            means = rep(mean, C),
                                            sds = rep(sd, C),
                                            betas = rep(beta, C),
                                            precondition_values = rep(1, C),
                                            ESS_threshold = 0.5,
                                            diffusion_estimator = diffusion_estimator,
                                            seed = seed)
print('ESS:'); print(BF_standard_n20$ESS)
print('CESS:'); print(BF_standard_n20$CESS)

print('### performing Bayesian Fusion (with n=20) with a preconditioning matrix')
BF_generalised_n20 <- parallel_GBF_uniGaussian(particles_to_fuse = input_particles_BF_n20,
                                               N = nsamples,
                                               m = C,
                                               time_mesh = time_mesh_BF_n20,
                                               means = rep(mean, C),
                                               sds = rep(sd, C),
                                               betas = rep(beta, C),
                                               precondition_values = sapply(input_samples, var),
                                               ESS_threshold = 0.5,
                                               diffusion_estimator = diffusion_estimator,
                                               seed = seed)
print('ESS:'); print(BF_generalised_n20$ESS)
print('CESS:'); print(BF_generalised_n20$CESS)

# plots
# proposals
curve(dnorm(x, mean = 1.2, sd = sqrt(0.1)), -3, 5)
lines(density(BF_standard_n20$proposed_samples), col = 'red')
lines(density(BF_generalised_n20$proposed_samples), col = 'blue')
# resampled
curve(dnorm(x, mean = 1.2, sd = sqrt(0.1)), -3, 5)
lines(density(resample_particle_y_samples(particle_set = BF_standard_n20$particles, 
                                          multivariate = FALSE,
                                          resampling_method = 'resid',
                                          seed = seed)$y_samples), 
      col = 'red')
lines(density(resample_particle_y_samples(particle_set = BF_generalised_n20$particles, 
                                          multivariate = FALSE,
                                          resampling_method = 'resid',
                                          seed = seed)$y_samples), 
      col = 'blue')
