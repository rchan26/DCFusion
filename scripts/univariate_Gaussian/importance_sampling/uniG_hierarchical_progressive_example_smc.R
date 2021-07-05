library(hierarchicalFusion)

seed <- 1994
set.seed(seed)

input_samples <- lapply(1:8, function(i) rnorm_tempered(N = 50000,
                                                        mean = 1.2, 
                                                        sd = sqrt(0.1),
                                                        beta = 1/8))

######################################## 

hier_smc_standard_Poisson <- hierarchical_fusion_SMC_uniGaussian(N_schedule = rep(50000, 3), 
                                                                 m_schedule = rep(2, 3),
                                                                 time_schedule = rep(1, 3), 
                                                                 base_samples = input_samples, 
                                                                 L = 4, 
                                                                 mean = 1.2,
                                                                 sd = sqrt(0.1), 
                                                                 start_beta = 1/8,
                                                                 precondition = FALSE,
                                                                 diffusion_estimator = 'Poisson',
                                                                 resampling_method = 'resid',
                                                                 ESS_threshold = 0.5,
                                                                 seed = seed)
print(hier_smc_standard_Poisson[c('ESS', 'CESS', 'resampled')])

hier_smc_precondition_Poisson <- hierarchical_fusion_SMC_uniGaussian(N_schedule = rep(50000, 3), 
                                                                     m_schedule = rep(2, 3),
                                                                     time_schedule = rep(1, 3), 
                                                                     base_samples = input_samples, 
                                                                     L = 4, 
                                                                     mean = 1.2,
                                                                     sd = sqrt(0.1), 
                                                                     start_beta = 1/8,
                                                                     precondition = TRUE, 
                                                                     diffusion_estimator = 'Poisson',
                                                                     resampling_method = 'resid',
                                                                     ESS_threshold = 0.5,
                                                                     seed = seed)
print(hier_smc_precondition_Poisson[c('ESS', 'CESS', 'resampled')])

hier_smc_standard_NB <- hierarchical_fusion_SMC_uniGaussian(N_schedule = rep(50000, 3), 
                                                            m_schedule = rep(2, 3),
                                                            time_schedule = rep(1, 3), 
                                                            base_samples = input_samples, 
                                                            L = 4, 
                                                            mean = 1.2,
                                                            sd = sqrt(0.1), 
                                                            start_beta = 1/8,
                                                            precondition = FALSE,
                                                            diffusion_estimator = 'NB',
                                                            resampling_method = 'resid',
                                                            ESS_threshold = 0.5,
                                                            seed = seed)
print(hier_smc_standard_NB[c('ESS', 'CESS', 'resampled')])

hier_smc_precondition_NB <- hierarchical_fusion_SMC_uniGaussian(N_schedule = rep(50000, 3), 
                                                                m_schedule = rep(2, 3),
                                                                time_schedule = rep(1, 3), 
                                                                base_samples = input_samples, 
                                                                L = 4, 
                                                                mean = 1.2,
                                                                sd = sqrt(0.1), 
                                                                start_beta = 1/8,
                                                                precondition = TRUE, 
                                                                diffusion_estimator = 'NB',
                                                                resampling_method = 'resid',
                                                                ESS_threshold = 0.5,
                                                                seed = seed)
print(hier_smc_precondition_NB[c('ESS', 'CESS', 'resampled')])

curve(dnorm(x, mean = 1.2, sd = sqrt(0.1)), -3, 5)
lines(density(resample_particle_y_samples(particle_set = hier_smc_standard_Poisson$particles[[1]], 
                                          multivariate = FALSE,
                                          resampling_method = 'resid',
                                          seed = seed)$y_samples), 
      col = 'green')
lines(density(resample_particle_y_samples(particle_set = hier_smc_precondition_Poisson$particles[[1]], 
                                          multivariate = FALSE,
                                          resampling_method = 'resid',
                                          seed = seed)$y_samples), 
      col = 'blue')
lines(density(resample_particle_y_samples(particle_set = hier_smc_standard_NB$particles[[1]], 
                                          multivariate = FALSE,
                                          resampling_method = 'resid',
                                          seed = seed)$y_samples), 
      col = 'green', lty = 2)
lines(density(resample_particle_y_samples(particle_set = hier_smc_precondition_NB$particles[[1]], 
                                          multivariate = FALSE,
                                          resampling_method = 'resid',
                                          seed = seed)$y_samples), 
      col = 'blue', lty = 2)

######################################## 

prog_smc_standard_Poisson <- progressive_fusion_SMC_uniGaussian(N_schedule = rep(50000, 7),
                                                                time_schedule = rep(1, 7),
                                                                base_samples = input_samples, 
                                                                mean = 1.2, 
                                                                sd = sqrt(0.1), 
                                                                start_beta = 1/8,
                                                                precondition = FALSE, 
                                                                diffusion_estimator = 'Poisson',
                                                                resampling_method = 'resid',
                                                                ESS_threshold = 0.5,
                                                                seed = seed)
print(prog_smc_standard_Poisson[c('ESS', 'CESS', 'resampled')])

prog_smc_precondition_Poisson <- progressive_fusion_SMC_uniGaussian(N_schedule = rep(50000, 7),
                                                                    time_schedule = rep(1, 7),
                                                                    base_samples = input_samples, 
                                                                    mean = 1.2, 
                                                                    sd = sqrt(0.1), 
                                                                    start_beta = 1/8,
                                                                    precondition = TRUE, 
                                                                    diffusion_estimator = 'Poisson',
                                                                    resampling_method = 'resid',
                                                                    ESS_threshold = 0.5,
                                                                    seed = seed)
print(prog_smc_precondition_Poisson[c('ESS', 'CESS', 'resampled')])

prog_smc_standard_NB <- progressive_fusion_SMC_uniGaussian(N_schedule = rep(50000, 7),
                                                           time_schedule = rep(1, 7),
                                                           base_samples = input_samples, 
                                                           mean = 1.2, 
                                                           sd = sqrt(0.1), 
                                                           start_beta = 1/8,
                                                           precondition = FALSE, 
                                                           diffusion_estimator = 'NB',
                                                           resampling_method = 'resid',
                                                           ESS_threshold = 0.5,
                                                           seed = seed)
print(prog_smc_standard_NB[c('ESS', 'CESS', 'resampled')])

prog_smc_precondition_NB <- progressive_fusion_SMC_uniGaussian(N_schedule = rep(50000, 7),
                                                               time_schedule = rep(1, 7),
                                                               base_samples = input_samples, 
                                                               mean = 1.2, 
                                                               sd = sqrt(0.1), 
                                                               start_beta = 1/8,
                                                               precondition = TRUE, 
                                                               diffusion_estimator = 'NB',
                                                               resampling_method = 'resid',
                                                               ESS_threshold = 0.5,
                                                               seed = seed)
print(prog_smc_precondition_NB[c('ESS', 'CESS', 'resampled')])

curve(dnorm(x, mean = 1.2, sd = sqrt(0.1)), -3, 5)
lines(density(resample_particle_y_samples(particle_set = prog_smc_standard_Poisson$particles[[1]], 
                                          multivariate = FALSE,
                                          resampling_method = 'resid',
                                          seed = seed)$y_samples), 
      col = 'green')
lines(density(resample_particle_y_samples(particle_set = prog_smc_precondition_Poisson$particles[[1]], 
                                          multivariate = FALSE,
                                          resampling_method = 'resid',
                                          seed = seed)$y_samples), 
      col = 'blue')
lines(density(resample_particle_y_samples(particle_set = prog_smc_standard_NB$particles[[1]], 
                                          multivariate = FALSE,
                                          resampling_method = 'resid',
                                          seed = seed)$y_samples), 
      col = 'green', lty = 2)
lines(density(resample_particle_y_samples(particle_set = prog_smc_precondition_NB$particles[[1]], 
                                          multivariate = FALSE,
                                          resampling_method = 'resid',
                                          seed = seed)$y_samples), 
      col = 'blue', lty = 2)
