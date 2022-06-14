library(DCFusion)
time_choices_exp_4 <- c(0.1, 0.2, 0.5, 1, 1.5, 2, 2.5)
seed <- 2022
set.seed(seed)
C <- 5

# using rejection sampling to obtain input samples
set.seed(seed)
nsamples <- 10000
sub_posterior_samples <- base_rejection_sampler_exp_4(beta = 1/C,
                                                      nsamples = nsamples,
                                                      proposal_mean = 0,
                                                      proposal_sd = 1.5,
                                                      dominating_M = 1.4)
curve(exp_4_density(x, beta = 1/C), -4, 4, ylim = c(-0.5,1))
for (i in 1:5) {
  lines(density(sub_posterior_samples[[i]]), col = 'blue')
}

# direct sampling from truth
set.seed(seed)
true_samples <- unlist(base_rejection_sampler_exp_4(beta = 1,
                                                    nsamples = nsamples,
                                                    proposal_mean = 0,
                                                    proposal_sd = 1,
                                                    dominating_M = 1.4))
curve(exp_4_density(x), -4, 4)
curve(1.4*dnorm(x, 0, 1), add = T)
lines(density(true_samples), col = 'blue')

########## USING R CODE [standard, Poisson]
exp_4_fusion <- lapply(time_choices_exp_4, function(time) {
  print(paste('time:', time))
  return(parallel_fusion_exp_4(N = nsamples,
                               m = C,
                               time = time,
                               samples_to_fuse = sub_posterior_samples,
                               mean = 0,
                               betas = rep(1/C, C),
                               precondition_values = rep(1, C),
                               seed = seed))})

plot(x = time_choices_exp_4, 
     y = sapply(1:length(time_choices_exp_4), function(i) exp_4_fusion[[i]]$rho),
     xlab = '',
     ylab = '',
     xaxt = 'n',
     yaxt = 'n',
     type = 'b',
     pch = 4,
     lty = 2,
     lwd = 3,
     ylim = c(0, 1))
lines(x = time_choices_exp_4, 
      y = sapply(1:length(time_choices_exp_4), function(i) exp_4_fusion[[i]]$Q),
      lty = 3,
      type = 'b',
      lwd = 3)
mtext('T', 1, 2.75, font = 2, cex = 1.5)
mtext('Acceptance Probability', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=c(0.1, seq(0.5, 2.5, 0.5)), labels=c(0.1, 0.5, "1.0", 1.5, "2.0", 2.5), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=c("0.0", c(seq(0.1, 0.9, 0.1), "1.0")), font = 2, cex = 1.5)
legend(x = 0.1, y = 1,
       legend = c(expression(rho^bm), expression(Q^bm)),
       lwd = c(3, 3),
       lty = c(2, 3),
       pch = c(4, 1),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

plot(x = time_choices_exp_4, 
     y = log(sapply(1:length(time_choices_exp_4), function(i) exp_4_fusion[[i]]$time), 2),
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
axis(1, at=c(0.1, seq(0.5, 2.5, 0.5)), labels=c(0.1, 0.5, "1.0", 1.5, "2.0", 2.5), font = 2, cex = 1.5)
axis(2, at=seq(0, 10, 1), labels=seq(0, 10, 1), font = 2, cex = 1.5)

curve(exp_4_density(x), -4, 4)
for (i in 1:length(time_choices_exp_4)) {
  print(paste('time:', time_choices_exp_4[i]))
  print(ks.test(exp_4_fusion[[i]]$samples, true_samples))
  print(integrated_abs_distance_exp_4(exp_4_fusion[[i]]$samples))
  lines(density(exp_4_fusion[[i]]$samples), col = 'red')
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

bw <- 0.1
curve(exp_4_density(x), -3, 3, ylim = c(0, 1), ylab = '', xlab = '', lwd = 5, xaxt = 'n', yaxt = 'n')
mtext('x', 1, 2.75, font = 2, cex = 1.5)
mtext('Density', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=seq(-3, 3, 1), labels=seq(-3, 3, 1), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=c("0.0", c(seq(0.1, 0.9, 0.1), "1.0")), font = 2, cex = 1.5)
lines(density(consensus$samples, bw), col = '#FFC20A', lty = 4, lwd = 3)
lines(density(neisw$samples, bw), col = '#0C7BDC', lty = 3, lwd = 3)
lines(density(weier$samples, bw), col = '#22FF22', lty = 2, lwd = 3)
lines(density(exp_4_fusion[[6]]$samples, bw), col = 'red', lty = 5, lwd = 3)
legend(x = -3, y = 1,
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
