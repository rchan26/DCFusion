library(hierarchicalFusion)

######################################## examples ########################################

seed <- 21
set.seed(seed)
target_mc <- sample_exp_4(N = 40000,
                          proposal_mean = 0,
                          proposal_sd = 1,
                          dominating_M = 1.35,
                          beta = 1)

######################################## beta = 1/8

# using rejection sampling to obtain input samples
input_samples <- base_rejection_sampler_exp_4(beta = 1/8,
                                              nsamples = 50000,
                                              proposal_mean = 0,
                                              proposal_sd = 1.5,
                                              dominating_M = 1.3)
curve(exp_4_density(x, beta = 1/8), -4, 4)
# check the samples look okay
for (samples in input_samples) {
   lines(density(samples), col = 'blue')
}

########################################

fork_and_join_standard <- hierarchical_fusion_SMC_exp_4(N_schedule = 50000, 
                                                        m_schedule = 8, 
                                                        time_schedule = 1,
                                                        base_samples = input_samples, 
                                                        L = 2, 
                                                        mean = 0,
                                                        start_beta = 1/8, 
                                                        precondition = FALSE, 
                                                        ESS_threshold = 0.5, 
                                                        resampling_method = 'resid',
                                                        seed = seed)

fork_and_join_precondition <- hierarchical_fusion_SMC_exp_4(N_schedule = 50000, 
                                                            m_schedule = 8, 
                                                            time_schedule = 1, 
                                                            base_samples = input_samples, 
                                                            L = 2, 
                                                            mean = 0,
                                                            start_beta = 1/8, 
                                                            precondition = TRUE,
                                                            ESS_threshold = 0.5, 
                                                            resampling_method = 'resid',
                                                            seed = seed)

fork_and_join_standard$ESS
fork_and_join_precondition$ESS
fork_and_join_standard$CESS
fork_and_join_precondition$CESS

curve(exp_4_density(x), -2.5, 2.5)
lines(density(target_mc), col = 'black')
lines(density(resample_particle_y_samples(particle_set = fork_and_join_standard$particles[[1]],
                                          multivariate = FALSE,
                                          resampling_method = 'resid',
                                          seed = seed)$y_samples),
      col = 'green')
lines(density(resample_particle_y_samples(particle_set = fork_and_join_precondition$particles[[1]], 
                                          multivariate = FALSE,
                                          resampling_method = 'resid',
                                          seed = seed)$y_samples), 
      col = 'blue')

########################################

# standard hierarchical fusion
hier_smc_standard <- hierarchical_fusion_SMC_exp_4(N_schedule = rep(50000, 3),
                                                   m_schedule = rep(2, 3),
                                                   time_schedule = rep(1, 3),
                                                   base_samples = input_samples,
                                                   L = 4,
                                                   mean = 0,
                                                   start_beta = 1/8,
                                                   precondition = FALSE,
                                                   ESS_threshold = 0.5, 
                                                   resampling_method = 'resid',
                                                   seed = seed)

print(hier_smc_standard$ESS)
print(hier_smc_standard$CESS)
print(hier_smc_standard$resampled)

# preconditioned hierarchical fusion
hier_smc_precondition <- hierarchical_fusion_SMC_exp_4(N_schedule = rep(50000, 3),
                                                       m_schedule = rep(2, 3),
                                                       time_schedule = rep(1, 3),
                                                       base_samples = input_samples,
                                                       L = 4,
                                                       mean = 0,
                                                       start_beta = 1/8,
                                                       precondition = TRUE,
                                                       ESS_threshold = 0.5, 
                                                       resampling_method = 'resid',
                                                       seed = seed)

print(hier_smc_precondition$ESS)
print(hier_smc_precondition$CESS)
print(hier_smc_precondition$resampled)

curve(exp_4_density(x), -2.5, 2.5)
lines(density(target_mc), col = 'black')
lines(density(resample_particle_y_samples(particle_set = hier_smc_standard$particles[[1]], 
                                          multivariate = FALSE, 
                                          resampling_method = 'resid',
                                          seed = seed)$y_samples), 
      col = 'green')
lines(density(resample_particle_y_samples(particle_set = hier_smc_precondition$particles[[1]], 
                                          multivariate = FALSE, 
                                          resampling_method = 'resid',
                                          seed = seed)$y_samples), 
      col = 'blue')

########################################

# standard progressive fusion
prog_smc_standard <- progressive_fusion_SMC_exp_4(N_schedule = rep(50000, 7),
                                                  time_schedule = rep(1, 7),
                                                  base_samples = input_samples,
                                                  mean = 0,
                                                  start_beta = 1/8,
                                                  precondition = FALSE,
                                                  ESS_threshold = 0.5, 
                                                  resampling_method = 'resid',
                                                  seed = seed)

print(prog_smc_standard$ESS)
print(prog_smc_standard$CESS)
print(prog_smc_standard$resampled)

# preconditioned progressive fusion
prog_smc_precondition <- progressive_fusion_SMC_exp_4(N_schedule = rep(50000, 7),
                                                      time_schedule = rep(1, 7),
                                                      base_samples = input_samples,
                                                      mean = 0,
                                                      start_beta = 1/8,
                                                      precondition = TRUE,
                                                      ESS_threshold = 0.5, 
                                                      resampling_method = 'resid',
                                                      seed = seed)

print(prog_smc_precondition$ESS)
print(prog_smc_precondition$CESS)
print(prog_smc_precondition$resampled)

curve(exp_4_density(x), -2.5, 2.5)
lines(density(target_mc), col = 'black')
lines(density(resample_particle_y_samples(particle_set = prog_smc_standard$particles[[1]], 
                                          multivariate = FALSE, 
                                          resampling_method = 'resid',
                                          seed = seed)$y_samples), 
      col = 'green')
lines(density(resample_particle_y_samples(particle_set = prog_smc_precondition$particles[[1]], 
                                          multivariate = FALSE, 
                                          resampling_method = 'resid',
                                          seed = seed)$y_samples), 
      col = 'blue')
