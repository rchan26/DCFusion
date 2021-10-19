library(DCFusion)

######################################## examples ########################################

seed <- 21
set.seed(seed)

# setting variables
n_comp <- 2
w_ex <- c(0.35, 0.65)
m_ex <- c(-3, 12)
s_ex <- c(1, 1.5)
b_ex <- 1/2

# sampling from tempered density
nsamples <- 100000
base <- base_rejection_sampler_mixG(nsamples = nsamples,
                                    weights = w_ex,
                                    means = m_ex,
                                    sds = s_ex,
                                    beta = b_ex,
                                    proposal_sds = c(2, 2.5),
                                    dominating_M = 2)
particles_to_fuse <- initialise_particle_sets(samples_to_fuse = base, multivariate = FALSE)

curve(dnorm_mix_tempered(x, n_comp, w_ex, m_ex, s_ex, b_ex), -11, 20, ylab = 'pdf', n = 100000)
for (particles in particles_to_fuse) {
  lines(density(particles$y_samples, adjust = 0.1), col = 'blue')
}
# normal fusion
test_standard_smc <- parallel_fusion_SMC_mixG(particles_to_fuse = particles_to_fuse,
                                              N = 100000,
                                              m = 2,
                                              time = 1,
                                              n_comp = n_comp,
                                              weights = w_ex,
                                              means = m_ex,
                                              sds = s_ex,
                                              betas = rep(b_ex, 2),
                                              precondition_values = rep(1, 2), 
                                              resampling_method = 'resid',
                                              ESS_threshold = 0.5,
                                              bounds_multiplier = 1.2,
                                              seed = seed)

# pre-conditioned fusion
test_precondition_smc <- parallel_fusion_SMC_mixG(particles_to_fuse = particles_to_fuse,
                                                  N = 100000,
                                                  m = 2,
                                                  time = 1/50,
                                                  n_comp = n_comp,
                                                  weights = w_ex,
                                                  means = m_ex,
                                                  sds = s_ex,
                                                  betas = rep(b_ex, 2),
                                                  precondition_values = sapply(base, var),
                                                  resampling_method = 'resid',
                                                  ESS_threshold = 0.5,
                                                  bounds_multiplier = 1.2,
                                                  seed = seed)

curve(dnorm_mix(x, n_comp = n_comp, w_ex, m_ex, s_ex), -11, 20, ylab = 'pdf', n = 100000)
lines(density(test_standard_smc$proposed_samples), col = 'blue', lty = 2)
lines(density(resample_particle_y_samples(particle_set = test_standard_smc$particles,
                                          multivariate = FALSE,
                                          resampling_method = 'resid',
                                          seed = seed)$y_samples,
              adjust = 0.1), col = 'red')
lines(density(test_precondition_smc$proposed_samples), col = 'blue', lty = 2)
lines(density(resample_particle_y_samples(particle_set = test_precondition_smc$particles,
                                          multivariate = FALSE,
                                          resampling_method = 'resid',
                                          seed = seed)$y_samples,
              adjust = 0.1), col = 'red')
