library(DCFusion)
load('CC32.RData')
load('CC64.RData')
load('CC128.RData')

# P_fusion <- c(integrated_abs_distance(full_posterior,
#                                       Poisson_hc_32$particles$y_samples),
#               integrated_abs_distance(full_posterior,
#                                       Poisson_hc_64$particles$y_samples),
#               integrated_abs_distance(full_posterior,
#                                       Poisson_hc_128$particles$y_samples))
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

plot(x = log(c(32, 64, 128), 2), y = NB_fusion,
     ylim = c(0, 0.5),
     xlab = '',
     ylab = '',
     xaxt = 'n', lwd = 3, pch = 20, type = 'b')
mtext('log(C, 2)', 1, 2.75, font = 2, cex = 1.5)
mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=c(seq(0, 0.9, 0.1), 0.95), labels=c("0.0", c(seq(0.1, 0.9, 0.1), 0.95)), font = 2, cex = 1.5)
axis(1, at=log(c(32, 64, 128), 2), labels = log(c(32, 64, 128), 2), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=c("0.0", seq(0.1, 0.9, 0.1), "1.0"), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5)
lines(x = log(c(32, 64, 128), 2), y = consensus,
      lty = 4, lwd = 3, type = 'b', pch = 3, col = 'red')
lines(x = log(c(32, 64, 128), 2), y = neiswanger,
      lty = 5, lwd = 3, type = 'b', pch = 2, col = 'red')
lines(x = log(c(32, 64, 128), 2), y = weierstrass,
      lty = 6, lwd = 3, type = 'b', pch = 1, col = 'red')
legend(x = 5, y = 0.5,
       legend = c('D&C-GMCF', 'CMC', 'KDEMC', 'WRS'),
       lwd = c(3, 3, 3),
       pch = c(20, 3, 2, 1),
       lty = c(1, 4, 5, 6),
       col = c('black', rep('red', 3)),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

######################### TIME

# P_fusion_time <- c(sum(unlist(Poisson_hc_32$time)),
#                    sum(unlist(Poisson_hc_64$time)),
#                    sum(unlist(Poisson_hc_128$time)))
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

plot(x = log(c(32, 64, 128), 2), y = log(NB_fusion_time, 2),
     ylim = c(0, 13),
     xlab = '',
     ylab = '',
     xaxt = 'n', lwd = 3, pch = 20, type = 'b')
mtext('log(C, 2)', 1, 2.75, font = 2, cex = 1.5)
mtext('log(Time elapsed in seconds, 2)', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=c(seq(0, 0.9, 0.1), 0.95), labels=c("0.0", c(seq(0.1, 0.9, 0.1), 0.95)), font = 2, cex = 1.5)
axis(1, at=log(c(32, 64, 128), 2), labels = log(c(32, 64, 128), 2), font = 2, cex = 1.5)
axis(2, at=seq(0, 12, 2), labels = seq(0, 12, 2), font = 2, cex = 1.5)
axis(2, at=seq(0, 13, 1), labels=rep("", 14), lwd.ticks = 0.5)
lines(x = log(c(32, 64, 128), 2), y = log(consensus_time, 2),
      lty = 4, lwd = 3, type = 'b', pch = 3, col = 'red')
lines(x = log(c(32, 64, 128), 2), y = log(neiswanger_time, 2),
      lty = 5, lwd = 3, type = 'b', pch = 2, col = 'red')
lines(x = log(c(32, 64, 128), 2), y = log(weierstrass_time, 2),
      lty = 6, lwd = 3, type = 'b', pch = 1, col = 'red')
legend(x = 5, y = 13,
       legend = c('D&C-GMCF', 'CMC', 'KDEMC', 'WRS'),
       lwd = c(3, 3, 3),
       pch = c(20, 3, 2, 1),
       lty = c(1, 4, 5, 6),
       col = c('black', rep('red', 3)),
       cex = 1.25,
       text.font = 2,
       bty = 'n')
