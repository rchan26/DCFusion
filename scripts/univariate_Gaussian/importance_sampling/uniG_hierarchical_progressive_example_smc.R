library(hierarchicalFusion)

seed <- 1994
set.seed(seed)

input_samples <- lapply(1:8, function(i) rnorm_tempered(N = 50000,
                                                        mean = 1.2, 
                                                        sd = sqrt(0.1), 
                                                        beta = 1/8))

######################################## 

hier_smc_standard <- hierarchical_fusion_SMC_uniGaussian(N_schedule = rep(50000, 3), 
                                                         m_schedule = rep(2, 3),
                                                         time_schedule = rep(1, 3), 
                                                         base_samples = input_samples, 
                                                         L = 4, 
                                                         mean = 1.2,
                                                         sd = sqrt(0.1), 
                                                         start_beta = 1/8,
                                                         precondition = FALSE, 
                                                         resampling_method = 'resid',
                                                         ESS_threshold = 0.5,
                                                         seed = seed)

print(hier_smc_standard$ESS)
print(hier_smc_standard$CESS)
print(hier_smc_standard$resampled)

hier_smc_precondition <- hierarchical_fusion_SMC_uniGaussian(N_schedule = rep(50000, 3), 
                                                             m_schedule = rep(2, 3),
                                                             time_schedule = rep(1, 3), 
                                                             base_samples = input_samples, 
                                                             L = 4, 
                                                             mean = 1.2,
                                                             sd = sqrt(0.1), 
                                                             start_beta = 1/8,
                                                             precondition = TRUE, 
                                                             resampling_method = 'resid',
                                                             ESS_threshold = 0.5,
                                                             seed = seed)

print(hier_smc_precondition$ESS)
print(hier_smc_precondition$CESS)
print(hier_smc_precondition$resampled)

curve(dnorm(x, mean = 1.2, sd = sqrt(0.1)), -3, 5)
lines(density(resample_particle_set(hier_smc_standard$particles[[1]], 
                                    multivariate = FALSE,
                                    resampling_method = 'resid',
                                    seed = seed)), 
      col = 'green')
lines(density(resample_particle_set(hier_smc_precondition$particles[[1]], 
                                    multivariate = FALSE,
                                    resampling_method = 'resid',
                                    seed = seed)), 
      col = 'blue')

######################################## 

prog_smc_standard <- progressive_fusion_SMC_uniGaussian(N_schedule = rep(50000, 7),
                                                        time_schedule = rep(1, 7),
                                                        base_samples = input_samples, 
                                                        mean = 1.2, 
                                                        sd = sqrt(0.1), 
                                                        start_beta = 1/8,
                                                        precondition = FALSE, 
                                                        resampling_method = 'resid',
                                                        ESS_threshold = 0.5,
                                                        seed = seed)

print(prog_smc_standard$ESS)
print(prog_smc_standard$CESS)
print(prog_smc_standard$resampled)

prog_smc_precondition <- progressive_fusion_SMC_uniGaussian(N_schedule = rep(50000, 7),
                                                            time_schedule = rep(1, 7),
                                                            base_samples = input_samples, 
                                                            mean = 1.2, 
                                                            sd = sqrt(0.1), 
                                                            start_beta = 1/8,
                                                            precondition = TRUE, 
                                                            resampling_method = 'resid',
                                                            ESS_threshold = 0.5,
                                                            seed = seed)

print(prog_smc_precondition$ESS)
print(prog_smc_precondition$CESS)
print(prog_smc_precondition$resampled)

curve(dnorm(x, mean = 1.2, sd = sqrt(0.1)), -3, 5)
lines(density(resample_particle_set(prog_smc_standard$particles[[1]], 
                                    multivariate = FALSE,
                                    resampling_method = 'resid',
                                    seed = seed)), 
      col = 'green')
lines(density(resample_particle_set(prog_smc_precondition$particles[[1]], 
                                    multivariate = FALSE,
                                    resampling_method = 'resid',
                                    seed = seed)), 
      col = 'blue')
