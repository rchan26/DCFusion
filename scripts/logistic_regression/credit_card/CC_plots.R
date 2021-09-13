library(hierarchicalFusion)
load('CC32.RData')
load('CC64.RData')
load('CC128.RData')
save.image('CC.RData')

P_fusion <- c(integrated_abs_distance(full_posterior,
                                      Poisson_hc_32$particles$y_samples),
              integrated_abs_distance(full_posterior,
                                      Poisson_hc_64$particles$y_samples),
              integrated_abs_distance(full_posterior,
                                      Poisson_hc_128$particles$y_samples))
NB_fusion <- c(integrated_abs_distance(full_posterior,
                                       NB_hc_32$particles$y_samples),
               integrated_abs_distance(full_posterior,
                                       NB_hc_64$particles$y_samples),
               integrated_abs_distance(full_posterior,
                                       NB_hc_128$particles$y_samples))
consensus <- c(integrated_abs_distance(full_posterior,
                                       consensus_mat_32$samples),
               integrated_abs_distance(full_posterior,
                                       consensus_mat_64$samples),
               integrated_abs_distance(full_posterior,
                                       consensus_mat_128$samples))
neiswanger <- c(integrated_abs_distance(full_posterior,
                                        neiswanger_false_32$samples),
                integrated_abs_distance(full_posterior,
                                        neiswanger_false_64$samples),
                integrated_abs_distance(full_posterior,
                                        neiswanger_false_128$samples))
weierstrass <- c(integrated_abs_distance(full_posterior,
                                         weierstrass_rejection_32$samples),
                 integrated_abs_distance(full_posterior,
                                         weierstrass_rejection_64$samples),
                 integrated_abs_distance(full_posterior,
                                         weierstrass_rejection_128$samples))

plot(x = c(32, 64, 128), y = NB_fusion,
     ylim = c(0, 1),
     xlab = '',
     ylab = '',
     xaxt = 'n', lwd = 3, pch = 1)
mtext('Number of sub-posteriors (C)', 1, 2.75, font = 2, cex = 1.5)
mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=c(seq(0, 0.9, 0.1), 0.95), labels=c("0.0", c(seq(0.1, 0.9, 0.1), 0.95)), font = 2, cex = 1.5)
axis(1, at=c(32, 64, 128), labels = c(32, 64, 128), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.2), labels=c("0.0", seq(0.2, 0.8, 0.2), "1.0"), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5)
lines(x = c(32, 64, 128), y = NB_fusion,
      lty = 1, lwd = 3)
points(x = c(32, 64, 128), y = consensus,
       pch = 4, lwd = 3)
lines(x = c(32, 64, 128), y = consensus,
      lty = 4, lwd = 3)
points(x = c(32, 64, 128), y = neiswanger,
       pch = 5, lwd = 3)
lines(x = c(32, 64, 128), y = neiswanger,
      lty = 2, lwd = 3)
points(x = c(32, 64, 128), y = weierstrass,
       pch = 6, lwd = 3)
lines(x = c(32, 64, 128), y = weierstrass,
      lty = 3, lwd = 3)
legend(x = 32, y = 1.05,
       legend = c('H. Fusion', 'CMC', 'KDEMC', 'WRS'),
       # col = c('#D41159', '#FFC20A', '#0C7BDC', '#994F00'),
       lwd = c(3, 3, 3),
       pch = c(1, 4, 5, 6),
       lty = c(1, 4, 2, 3),
       bty = 'n')

######################### TIME

P_fusion_time <- c(sum(unlist(Poisson_hc_32$time)),
                   sum(unlist(Poisson_hc_64$time)),
                   sum(unlist(Poisson_hc_128$time)))
NB_fusion_time <- c(sum(unlist(NB_hc_32$time)),
                    sum(unlist(NB_hc_64$time)),
                    sum(unlist(NB_hc_128$time)))
consensus_time <- c(consensus_mat_32$time,
                    consensus_mat_64$time,
                    consensus_mat_128$time)
neiswanger_time <- c(neiswanger_false_32$time,
                     neiswanger_false_64$time,
                     neiswanger_false_128$time)
weierstrass_time <- c(weierstrass_rejection_32$time,
                      weierstrass_rejection_64$time,
                      weierstrass_rejection_128$time)

plot(x = c(32, 64, 128), y = log(NB_fusion_time),
     ylim = c(0, 12),
     xlab = '',
     ylab = '',
     xaxt = 'n', lwd = 3, pch = 1)
mtext('Number of sub-posteriors (C)', 1, 2.75, font = 2, cex = 1.5)
mtext('log(Time elapsed in seconds)', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=c(seq(0, 0.9, 0.1), 0.95), labels=c("0.0", c(seq(0.1, 0.9, 0.1), 0.95)), font = 2, cex = 1.5)
axis(1, at=c(32, 64, 128), labels = c(32, 64, 128), font = 2, cex = 1.5)
axis(2, at=seq(0, 12, 2), labels = seq(0, 12, 2), font = 2, cex = 1.5)
axis(2, at=seq(0, 12, 1), labels=rep("", 13), lwd.ticks = 0.5)
lines(x = c(32, 64, 128), y = log(NB_fusion_time),
      lty = 1, lwd = 3)
points(x = c(32, 64, 128), y = log(consensus_time),
       pch = 4, lwd = 3)
lines(x = c(32, 64, 128), y = log(consensus_time),
      lty = 4, lwd = 3)
points(x = c(32, 64, 128), y = log(neiswanger_time),
       pch = 5, lwd = 3)
lines(x = c(32, 64, 128), y = log(neiswanger_time),
      lty = 2, lwd = 3)
points(x = c(32, 64, 128), y = log(weierstrass_time),
       pch = 6, lwd = 3)
lines(x = c(32, 64, 128), y = log(weierstrass_time),
      lty = 3, lwd = 3)
legend(x = 32, y = 12.5,
       legend = c('H. Fusion', 'CMC', 'KDEMC', 'WRS'),
       # col = c('#D41159', '#FFC20A', '#0C7BDC', '#994F00'),
       lwd = c(3, 3, 3),
       pch = c(1, 4, 5, 6),
       lty = c(1, 4, 2, 3),
       bty = 'n')
