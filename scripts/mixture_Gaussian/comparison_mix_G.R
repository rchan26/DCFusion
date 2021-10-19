library(DCFusion)

######################################## examples ########################################

seed <- 21
set.seed(seed)

# setting variables
n_comp <- 3
w_ex <- c(0.5, 0.2, 0.3)
m_ex <- c(-5, 6, 12)
s_ex <- c(1, 2, 1.5)
b_ex <- 1/4

# sampling from tempered density
nsamples <- 100000
base <- base_rejection_sampler_mixG(nsamples = nsamples,
                                    weights = w_ex,
                                    means = m_ex,
                                    sds = s_ex,
                                    beta = b_ex,
                                    proposal_sds = c(3, 4.5, 4),
                                    dominating_M = 2)
curve(dnorm_mix_tempered(x, n_comp, w_ex, m_ex, s_ex, b_ex), -11, 20, ylab = 'pdf', n = 100000)
# curve(2*dnorm_mix(x, n_comp, w_ex, m_ex, c(3, 4.5, 4)), add = T)
for (samples in base) {
  lines(density(samples, adjust = 0.1), col = 'blue')
}

particles_to_fuse <- initialise_particle_sets(samples_to_fuse = base, multivariate = FALSE)

# normal fusion
test_standard <- parallel_fusion_mixG(N = 100000,
                                      m = 4,
                                      time = 1,
                                      samples_to_fuse = base,
                                      n_comp = n_comp,
                                      weights = w_ex,
                                      means = m_ex,
                                      sds = s_ex,
                                      betas = rep(b_ex, 4),
                                      precondition_values = rep(1, 4),
                                      bounds_multiplier = 1.2,
                                      seed = seed)

# pre-conditioned fusion
test_precondition <- parallel_fusion_mixG(N = 100000,
                                          m = 4,
                                          time = 1/50,
                                          samples_to_fuse = base,
                                          n_comp = n_comp,
                                          weights = w_ex,
                                          means = m_ex,
                                          sds = s_ex,
                                          betas = rep(b_ex, 4),
                                          precondition_values = sapply(base, var),
                                          bounds_multiplier = 1.2,
                                          seed = seed)

# normal fusion
test_standard_smc <- parallel_fusion_SMC_mixG(particles_to_fuse = particles_to_fuse,
                                              N = 100000,
                                              m = 4,
                                              time = 1,
                                              n_comp = n_comp,
                                              weights = w_ex,
                                              means = m_ex,
                                              sds = s_ex,
                                              betas = rep(b_ex, 4),
                                              precondition_values = rep(1, 4),
                                              resampling_method = 'resid',
                                              ESS_threshold = 0.5,
                                              bounds_multiplier = 1.2,
                                              seed = seed)

# pre-conditioned fusion
test_precondition_smc <- parallel_fusion_SMC_mixG(particles_to_fuse = particles_to_fuse,
                                                  N = 100000,
                                                  m = 4,
                                                  time = 1/50,
                                                  n_comp = n_comp,
                                                  weights = w_ex,
                                                  means = m_ex,
                                                  sds = s_ex,
                                                  betas = rep(b_ex, 4),
                                                  precondition_values = sapply(base, var),
                                                  resampling_method = 'resid',
                                                  ESS_threshold = 0.5,
                                                  bounds_multiplier = 1.2,
                                                  seed = seed)

curve(dnorm_mix(x, n_comp = n_comp, w_ex, m_ex, s_ex), -11, 20, ylab = 'pdf', n = 100000)
lines(density(test_standard$samples, adjust = 0.1), col = 'red')
lines(density(test_precondition$samples, adjust = 0.1), col = 'blue')
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

consensus <- consensus_scott(S = 4,
                             samples_to_combine = lapply(base, as.matrix),
                             indep = FALSE)
lines(density(consensus$samples))
neisw <- neiswanger(S = 4,
                    samples_to_combine = lapply(base, as.matrix),
                    anneal = FALSE)
lines(density(neisw$samples))
weier <- weierstrass(Samples = base,
                     method = 'importance', para.dim = 1)
lines(density(weier$samples))

##### plots #####

curve(dnorm_mix(x, n_comp = n_comp, w_ex, m_ex, s_ex), -11, 20, ylab = 'pdf', n = 100000, lwd = 3, ylim = c(0, 0.25))
lines(density(test_precondition$samples, adjust = 0.25), col = '#D41159', lty = 5, lwd = 3)
lines(density(consensus$samples), col = '#FFC20A', lty = 4, lwd = 3)
lines(density(neisw$samples), col = '#0C7BDC', lty = 3, lwd = 3) 
lines(density(weier$samples, adjust = 0.25), col = '#994F00', lty = 2, lwd = 3)
legend(x = 12, y = 0.25,
       legend = c('Target', 'Fusion', 'Consensus', 'Neiswanger', 'Weierstrass'),
       col = c('black', '#D41159', '#FFC20A', '#0C7BDC', '#994F00'),
       lwd = c(3, 3, 3, 3, 3),
       lty = c(1, 5, 4, 3, 2))
