library(hierarchicalFusion)

######################################## examples ########################################

seed <- 21
set.seed(seed)
target_mc <- sample_exp_4(N = 100000,
                          proposal_mean = 0,
                          proposal_sd = 1,
                          dominating_M = 1.35,
                          beta = 1)

######################################## beta = 1/8

# using rejection sampling to obtain input samples
input_samples <- base_rejection_sampler_exp_4(beta = 1/8,
                                              nsamples = 100000,
                                              proposal_mean = 0,
                                              proposal_sd = 1.5,
                                              dominating_M = 1.3)
curve(exp_4_density(x, beta = 1/8), -4, 4)
# check the samples look okay
for (samples in input_samples) {
   lines(density(samples), col = 'blue')
}

########################################

fork_and_join_standard_Poisson <- hierarchical_fusion_SMC_exp_4(N_schedule = 100000, 
                                                                m_schedule = 8, 
                                                                time_schedule = 1,
                                                                base_samples = input_samples, 
                                                                L = 2, 
                                                                mean = 0,
                                                                start_beta = 1/8, 
                                                                precondition = FALSE, 
                                                                ESS_threshold = 0.5, 
                                                                resampling_method = 'resid',
                                                                diffusion_estimator = 'Poisson',
                                                                seed = seed)

fork_and_join_precondition_Poisson <- hierarchical_fusion_SMC_exp_4(N_schedule = 100000, 
                                                                    m_schedule = 8, 
                                                                    time_schedule = 1, 
                                                                    base_samples = input_samples, 
                                                                    L = 2, 
                                                                    mean = 0,
                                                                    start_beta = 1/8, 
                                                                    precondition = TRUE,
                                                                    ESS_threshold = 0.5, 
                                                                    resampling_method = 'resid',
                                                                    diffusion_estimator = 'Poisson',
                                                                    seed = seed)

fork_and_join_standard_NB <- hierarchical_fusion_SMC_exp_4(N_schedule = 100000, 
                                                           m_schedule = 8, 
                                                           time_schedule = 1,
                                                           base_samples = input_samples, 
                                                           L = 2, 
                                                           mean = 0,
                                                           start_beta = 1/8, 
                                                           precondition = FALSE, 
                                                           ESS_threshold = 0.5, 
                                                           resampling_method = 'resid',
                                                           diffusion_estimator = 'NB',
                                                           gamma_NB_n_points = 2,
                                                           seed = seed)

fork_and_join_precondition_NB <- hierarchical_fusion_SMC_exp_4(N_schedule = 100000, 
                                                               m_schedule = 8, 
                                                               time_schedule = 1, 
                                                               base_samples = input_samples, 
                                                               L = 2, 
                                                               mean = 0,
                                                               start_beta = 1/8, 
                                                               precondition = TRUE,
                                                               ESS_threshold = 0.5, 
                                                               resampling_method = 'resid',
                                                               diffusion_estimator = 'NB',
                                                               gamma_NB_n_points = 2,
                                                               seed = seed)

fork_and_join_standard_Poisson[c('ESS', 'CESS')]
fork_and_join_precondition_Poisson[c('ESS', 'CESS')]
fork_and_join_standard_NB[c('ESS', 'CESS')]
fork_and_join_precondition_NB[c('ESS', 'CESS')]

curve(exp_4_density(x), -2.5, 2.5)
lines(density(target_mc), col = 'black')
lines(density(resample_particle_y_samples(particle_set = fork_and_join_standard_Poisson$particles[[1]],
                                          multivariate = FALSE,
                                          resampling_method = 'resid',
                                          seed = seed)$y_samples),
      col = 'green')
lines(density(resample_particle_y_samples(particle_set = fork_and_join_precondition_Poisson$particles[[1]], 
                                          multivariate = FALSE,
                                          resampling_method = 'resid',
                                          seed = seed)$y_samples), 
      col = 'blue')
lines(density(resample_particle_y_samples(particle_set = fork_and_join_standard_NB$particles[[1]],
                                          multivariate = FALSE,
                                          resampling_method = 'resid',
                                          seed = seed)$y_samples),
      col = 'green', lty = 2)
lines(density(resample_particle_y_samples(particle_set = fork_and_join_precondition_NB$particles[[1]], 
                                          multivariate = FALSE,
                                          resampling_method = 'resid',
                                          seed = seed)$y_samples), 
      col = 'blue', lty = 2)

########################################

# standard hierarchical fusion
hier_smc_standard_Poisson <- hierarchical_fusion_SMC_exp_4(N_schedule = rep(100000, 3),
                                                           m_schedule = rep(2, 3),
                                                           time_schedule = rep(1, 3),
                                                           base_samples = input_samples,
                                                           L = 4,
                                                           mean = 0,
                                                           start_beta = 1/8,
                                                           precondition = FALSE,
                                                           ESS_threshold = 0.5, 
                                                           resampling_method = 'resid',
                                                           diffusion_estimator = 'Poisson',
                                                           seed = seed)
print(hier_smc_standard_Poisson[c('ESS', 'CESS', 'resampled')])

# preconditioned hierarchical fusion
hier_smc_precondition_Poisson <- hierarchical_fusion_SMC_exp_4(N_schedule = rep(100000, 3),
                                                               m_schedule = rep(2, 3),
                                                               time_schedule = rep(1, 3),
                                                               base_samples = input_samples,
                                                               L = 4,
                                                               mean = 0,
                                                               start_beta = 1/8,
                                                               precondition = TRUE,
                                                               ESS_threshold = 0.5, 
                                                               resampling_method = 'resid',
                                                               diffusion_estimator = 'Poisson',
                                                               seed = seed)
print(hier_smc_standard_Poisson[c('ESS', 'CESS', 'resampled')])

# standard hierarchical fusion
hier_smc_standard_NB <- hierarchical_fusion_SMC_exp_4(N_schedule = rep(100000, 3),
                                                      m_schedule = rep(2, 3),
                                                      time_schedule = rep(1, 3),
                                                      base_samples = input_samples,
                                                      L = 4,
                                                      mean = 0,
                                                      start_beta = 1/8,
                                                      precondition = FALSE,
                                                      ESS_threshold = 0.5, 
                                                      resampling_method = 'resid',
                                                      diffusion_estimator = 'NB',
                                                      seed = seed)
print(hier_smc_standard_NB[c('ESS', 'CESS', 'resampled')])

# preconditioned hierarchical fusion
hier_smc_precondition_NB <- hierarchical_fusion_SMC_exp_4(N_schedule = rep(100000, 3),
                                                          m_schedule = rep(2, 3),
                                                          time_schedule = rep(1, 3),
                                                          base_samples = input_samples,
                                                          L = 4,
                                                          mean = 0,
                                                          start_beta = 1/8,
                                                          precondition = TRUE,
                                                          ESS_threshold = 0.5, 
                                                          resampling_method = 'resid',
                                                          diffusion_estimator = 'NB',
                                                          seed = seed)
print(hier_smc_standard_NB[c('ESS', 'CESS', 'resampled')])

curve(exp_4_density(x), -2.5, 2.5)
lines(density(target_mc), col = 'black')
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

# standard progressive fusion
prog_smc_standard_Poisson <- progressive_fusion_SMC_exp_4(N_schedule = rep(100000, 7),
                                                          time_schedule = rep(1, 7),
                                                          base_samples = input_samples,
                                                          mean = 0,
                                                          start_beta = 1/8,
                                                          precondition = FALSE,
                                                          ESS_threshold = 0.5, 
                                                          resampling_method = 'resid',
                                                          diffusion_estimator = 'Poisson',
                                                          seed = seed)
print(prog_smc_standard_Poisson[c('ESS', 'CESS', 'resampled')])

# preconditioned progressive fusion
prog_smc_precondition_Poisson <- progressive_fusion_SMC_exp_4(N_schedule = rep(100000, 7),
                                                              time_schedule = rep(1, 7),
                                                              base_samples = input_samples,
                                                              mean = 0,
                                                              start_beta = 1/8,
                                                              precondition = TRUE,
                                                              ESS_threshold = 0.5, 
                                                              resampling_method = 'resid',
                                                              diffusion_estimator = 'Poisson',
                                                              seed = seed)
print(prog_smc_precondition_Poisson[c('ESS', 'CESS', 'resampled')])

# standard progressive fusion
prog_smc_standard_NB <- progressive_fusion_SMC_exp_4(N_schedule = rep(100000, 7),
                                                     time_schedule = rep(1, 7),
                                                     base_samples = input_samples,
                                                     mean = 0,
                                                     start_beta = 1/8,
                                                     precondition = FALSE,
                                                     ESS_threshold = 0.5, 
                                                     resampling_method = 'resid',
                                                     diffusion_estimator = 'NB',
                                                     seed = seed)
print(prog_smc_standard_NB[c('ESS', 'CESS', 'resampled')])

# preconditioned progressive fusion
prog_smc_precondition_NB <- progressive_fusion_SMC_exp_4(N_schedule = rep(100000, 7),
                                                         time_schedule = rep(1, 7),
                                                         base_samples = input_samples,
                                                         mean = 0,
                                                         start_beta = 1/8,
                                                         precondition = TRUE,
                                                         ESS_threshold = 0.5, 
                                                         resampling_method = 'resid',
                                                         diffusion_estimator = 'NB',
                                                         seed = seed)
print(prog_smc_precondition_NB[c('ESS', 'CESS', 'resampled')])

curve(exp_4_density(x), -2.5, 2.5)
lines(density(target_mc), col = 'black')
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
