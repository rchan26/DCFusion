library(DCFusion)
time_choices_mixG <- c(0.25, 0.5, 1, 1.5, 2, 2.5, 3)
seed <- 2022
set.seed(seed)
C <- 4

# setting variables
n_components <- 3
w_ex <- c(0.35, 0.2, 0.45)
m_ex <- c(-3, 2, 5)
s_ex <- c(1, 1.5, 0.5)
b_ex <- 1/C

# sampling from tempered density
nsamples <- 10000
base <- base_rejection_sampler_mixG(nsamples = nsamples, 
                                    weights = w_ex, 
                                    means = m_ex,
                                    sds = s_ex,
                                    beta = b_ex,
                                    proposal_sds = c(3, 3, 3),
                                    dominating_M = 2)

# check the base samples look okay
curve(dnorm_mix_tempered(x, n_components, w_ex, m_ex, s_ex, b_ex), -12.5, 15, ylab = 'pdf')
curve(2*dnorm_mix(x, n_components, w_ex, m_ex, c(3, 3, 3)), add = T, lty = 2, col = 2)
for (samples in base) {
  lines(density(samples, adjust = 0.5), col = 'blue')
}

# fusion
mixG_fusion <- lapply(time_choices_mixG, function(time) {
  print(paste('time:', time))
  return(parallel_fusion_mixG(N = nsamples,
                              m = C,
                              time = time,
                              samples_to_fuse = base,
                              n_comp = n_components,
                              weights = w_ex,
                              means = m_ex,
                              sds = s_ex,
                              betas = rep(b_ex, C),
                              precondition_values = rep(1, C),
                              bounds_multiplier = 1.1,
                              seed = seed))})

plot(x = time_choices_mixG,
     y = sapply(1:length(time_choices_mixG), function(i) mixG_fusion[[i]]$rho),
     xlab = '',
     ylab = '',
     xaxt = 'n',
     yaxt = 'n',
     type = 'b',
     pch = 4,
     lty = 2,
     lwd = 3,
     ylim = c(0, 1))
lines(x = time_choices_mixG, 
      y = sapply(1:length(time_choices_mixG), function(i) mixG_fusion[[i]]$Q),
      lty = 3,
      type = 'b',
      lwd = 3)
mtext('T', 1, 2.75, font = 2, cex = 1.5)
mtext('Acceptance Probability', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=c(0.25, seq(0, 3, 0.5)), labels=c(0.25, seq(0, 3, 0.5)), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=c("0.0", c(seq(0.1, 0.9, 0.1), "1.0")), font = 2, cex = 1.5)
legend(x = 0.25, y = 1,
       legend = c(expression(rho^bm), expression(Q^bm)),
       lwd = c(3, 3),
       lty = c(2, 3),
       pch = c(4, 1),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

plot(x = time_choices_mixG, 
     y = log(sapply(1:length(time_choices_mixG), function(i) mixG_fusion[[i]]$time), 2),
     xlab = '',
     ylab = '',
     xaxt = 'n',
     yaxt = 'n',
     type = 'b',
     pch = 4,
     lty = 1,
     lwd = 3,
     ylim = c(0, 10))
mtext('T', 1, 2.75, font = 2, cex = 1.5)
mtext('log(Time elapsed in seconds, 2)', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=c(0.25, seq(0, 3, 0.5)), labels=c(0.25, seq(0, 3, 0.5)), font = 2, cex = 1.5)
axis(2, at=seq(0, 15, 1), labels=seq(0, 15, 1), font = 2, cex = 1.5)

curve(dnorm_mix(x, n_components, w_ex, m_ex, s_ex), -11, 12, ylab = 'pdf', n = 10000)
for (i in 1:length(time_choices_mixG)) {
  lines(density(mixG_fusion[[i]]$samples, adjust = 0.5), col = 'red')
}

##### Comparison #####

consensus <- consensus_scott(S = C,
                             samples_to_combine = lapply(base, as.matrix),
                             indep = FALSE)
neisw <- neiswanger(S = C,
                    samples_to_combine = lapply(base, as.matrix),
                    anneal = FALSE)
weier <- weierstrass(Samples = base,
                     method = 'importance', para.dim = 1)

lines(density(neisw$samples))
lines(density(consensus$samples))
lines(density(weier$samples))

bw <- 0.1
curve(dnorm_mix(x, n_components, w_ex, m_ex, s_ex), -8, 8, n = 10000, 
      ylim = c(0, 0.4), ylab = '', xlab = '', lwd = 5, xaxt = 'n', yaxt = 'n')
mtext('x', 1, 2.75, font = 2, cex = 1.5)
mtext('Density', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=seq(-8, 8, 2), labels=seq(-8, 8, 2), font = 2, cex = 1.5)
axis(1, at=seq(-8, 8, 1), labels=rep("", 17), lwd.ticks = 0.5)
axis(2, at=seq(0, 1, 0.1), labels=c("0.0", c(seq(0.1, 0.9, 0.1), "1.0")), font = 2, cex = 1.5)
lines(density(consensus$samples, bw), col = '#FFC20A', lty = 4, lwd = 3)
lines(density(neisw$samples, bw), col = '#0C7BDC', lty = 3, lwd = 3)
lines(density(weier$samples, bw), col = '#22FF22', lty = 2, lwd = 3)
lines(density(mixG_fusion[[5]]$samples, bw), col = 'red', lty = 5, lwd = 3)
legend(x = -8, y = 0.4,
       legend = c('Target',
                  'MCF',
                  'CMC',
                  'KDEMC',
                  'WRS'),
       lwd = c(3, 3, 3, 3, 3),
       lty = c(1, 5, 4, 3, 2),
       col = c('black', 'red', '#FFC20A', '#0C7BDC', '#22FF22'),
       cex = 1.25,
       text.font = 2,
       bty = 'n')
