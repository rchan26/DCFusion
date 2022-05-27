library(DCFusion)
time_choices_mixG <- c(0.1, 0.2, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5)
seed <- 2022
set.seed(seed)
C <- 5

# setting variables
C <- 5
w_ex <- c(0.35, 0.2, 0.45)
m_ex <- c(-3, 2, 5)
s_ex <- c(1, 1.5, 0.5)
b_ex <- 1/C

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

# fusion
mixG_fusion <- lapply(time_choices_mixG, function(time) {
  print(paste('time:', time))
  return(parallel_fusion_mixG(N = 100000,
                              m = C,
                              time = time,
                              samples_to_fuse = base,
                              weights = w_ex,
                              means = m_ex,
                              sds = s_ex,
                              betas = rep(b_ex, C),
                              precondition_values = rep(1, C),
                              bounds_multiplier = 1.2,
                              seed = seed))})

plot(x = time_choices_exp_4, 
     y = sapply(1:length(time_choices_exp_4), function(i) exp_4_fusion[[i]]$rho),
     xlab = 'Time',
     ylab = 'Prob',
     pch = 4,
     ylim = c(0, 1))
lines(x = time_choices_exp_4, 
      y = sapply(1:length(time_choices_exp_4), function(i) exp_4_fusion[[i]]$rho))
points(x = time_choices_exp_4, 
       y = sapply(1:length(time_choices_exp_4), function(i) exp_4_fusion[[i]]$Q))
lines(x = time_choices_exp_4, 
      y = sapply(1:length(time_choices_exp_4), function(i) exp_4_fusion[[i]]$Q),
      lty = 2)

curve(dnorm_mix(x, w_ex, m_ex, s_ex), -11, 12, ylab = 'pdf', n = 10000)
for (i in 1:length(time_choices_exp_4)) {
  lines(density(mixG_fusion[[i]]$samples, adjust = 0.1), col = 'red')
}

##### Comparison #####

consensus <- consensus_scott(S = C,
                             samples_to_combine = lapply(sub_posterior_samples, as.matrix),
                             indep = FALSE)
neisw <- neiswanger(S = C,
                    samples_to_combine = lapply(sub_posterior_samples, as.matrix),
                    anneal = FALSE)
weier <- weierstrass(Samples = sub_posterior_samples,
                     method = 'importance', para.dim = 1)

lines(density(neisw$samples))
lines(density(consensus$samples))
lines(density(weier$samples))

curve(dnorm_mix(x, w_ex, m_ex, s_ex), -11, 12, ylab = 'pdf', n = 10000)
lines(density(mixG_fusion[[4]]$samples), col = '#D41159', lty = 5, lwd = 3)
lines(density(consensus$samples), col = '#FFC20A', lty = 4, lwd = 3)
lines(density(neisw$samples), col = '#0C7BDC', lty = 3, lwd = 3) 
lines(density(weier$samples), col = '#994F00', lty = 2, lwd = 3)
legend(x = 2, y = 1,
       legend = c('Target', 'Fusion', 'Consensus', 'Neiswanger', 'Weierstrass'),
       col = c('black', '#D41159', '#FFC20A', '#0C7BDC', '#994F00'),
       lwd = c(3, 3, 3, 3, 3),
       lty = c(1, 5, 4, 3, 2))
