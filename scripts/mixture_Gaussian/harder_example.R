library(hierarchicalFusion)

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
curve(dnorm_mix_tempered(x, n_comp, w_ex, m_ex, s_ex, b_ex), -11, 20, ylab = 'pdf', n = 100000)
for (samples in base) {
  lines(density(samples, adjust = 0.1), col = 'blue')
}

# normal fusion
test_standard <- parallel_fusion_mixG(N = 100000,
                                      n_comp = n_comp,
                                      weights = w_ex,
                                      means = m_ex,
                                      sds = s_ex,
                                      time = 1,
                                      m = 2,
                                      samples_to_fuse = base,
                                      betas = rep(1/2, 2),
                                      precondition_values = rep(1, 2), 
                                      bounds_multiplier = 1.2,
                                      seed = seed)

# pre-conditioned fusion
test_precondition <- parallel_fusion_mixG(N = 100000,
                                          n_comp = n_comp,
                                          weights = w_ex,
                                          means = m_ex,
                                          sds = s_ex,
                                          time = 1/50,
                                          m = 2,
                                          samples_to_fuse = base,
                                          betas = rep(1/2, 2),
                                          precondition_values = sapply(base, var),
                                          bounds_multiplier = 1.2,
                                          seed = seed)

curve(dnorm_mix(x, n_comp = n_comp, w_ex, m_ex, s_ex), -11, 20, ylab = 'pdf', n = 100000)
lines(density(test_standard$samples, adjust = 0.1), col = 'red')
lines(density(test_precondition$samples, adjust = 0.1), col = 'blue')

test_standard$rhoQ
test_precondition$rhoQ
