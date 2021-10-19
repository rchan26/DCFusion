library(DCFusion)

######################################## examples ########################################

seed <- 21
set.seed(seed)

# setting variables
w_ex <- c(0.35, 0.65)
m_ex <- c(-3, 2)
s_ex <- c(1, 1.5)
b_ex <- 1/2

# sampling from tempered density
nsamples <- 20000
base <- base_rejection_sampler_mixG(nsamples = nsamples, 
                                    weights = w_ex, 
                                    means = m_ex,
                                    sds = s_ex,
                                    beta = b_ex,
                                    proposal_sds = c(2, 2),
                                    dominating_M = 2)

# check the base samples look okay
curve(dnorm_mix_tempered(x, w_ex, m_ex, s_ex, b_ex), -12.5, 15, ylab = 'pdf')
curve(2*dnorm_mix(x, w_ex, m_ex, c(2, 2)), add = T, lty = 2, col = 2)
for (samples in base) {
  lines(density(samples, adjust = 0.5), col = 'blue')
}

# normal fusion
test_standard <- parallel_fusion_mixG(N = 100000,
                                      m = 2,
                                      time = 1,
                                      samples_to_fuse = base,
                                      weights = w_ex,
                                      means = m_ex,
                                      sds = s_ex,
                                      betas = rep(1/2, 2),
                                      precondition_values = rep(1, 2), 
                                      bounds_multiplier = 1.2,
                                      seed = seed, 
                                      n_cores = 12)

# pre-conditioned fusion
test_precondition <- parallel_fusion_mixG(N = 100000,
                                          m = 2,
                                          time = 0.1,
                                          samples_to_fuse = base,
                                          weights = w_ex,
                                          means = m_ex,
                                          sds = s_ex,
                                          betas = rep(1/2, 2),
                                          precondition_values = sapply(base, var), 
                                          bounds_multiplier = 1.2,
                                          seed = seed, 
                                          n_cores = 12)

curve(dnorm_mix(x, w_ex, m_ex, s_ex), -11, 12, ylab = 'pdf', n = 10000)
lines(density(test_standard$samples, adjust = 0.5), col = 'red')
lines(density(test_precondition$samples, adjust = 0.5), col = 'red')
