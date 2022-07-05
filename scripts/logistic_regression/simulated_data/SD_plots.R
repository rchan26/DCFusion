library(DCFusion)
load('SD4.RData')
load('SD8.RData')
load('SD16.RData')
load('SD32.RData')

NB_fusion <- c(integrated_abs_distance(full_posterior,
                                       NB_hc_4$particles$y_samples),
               integrated_abs_distance(full_posterior,
                                       NB_hc_8$particles$y_samples),
               integrated_abs_distance(full_posterior,
                                       NB_hc_16$particles$y_samples),
               integrated_abs_distance(full_posterior,
                                       NB_hc_32$particles$y_samples))
consensus <- c(integrated_abs_distance(full_posterior,
                                       consensus_mat_4$samples),
               integrated_abs_distance(full_posterior,
                                       consensus_mat_8$samples),
               integrated_abs_distance(full_posterior,
                                       consensus_mat_16$samples),
               integrated_abs_distance(full_posterior,
                                       consensus_mat_32$samples))
neiswanger <- c(integrated_abs_distance(full_posterior,
                                        neiswanger_false_4$samples),
                integrated_abs_distance(full_posterior,
                                        neiswanger_false_8$samples),
                integrated_abs_distance(full_posterior,
                                        neiswanger_false_16$samples),
                integrated_abs_distance(full_posterior,
                                        neiswanger_false_32$samples))
weierstrass <- c(integrated_abs_distance(full_posterior,
                                         weierstrass_rejection_4$samples),
                 integrated_abs_distance(full_posterior,
                                         weierstrass_rejection_8$samples),
                 integrated_abs_distance(full_posterior,
                                         weierstrass_rejection_16$samples),
                 integrated_abs_distance(full_posterior,
                                         weierstrass_rejection_32$samples))

plot(x = log(c(4, 8, 16, 32), 2), y = NB_fusion,
     ylim = c(0, 0.5),
     xlab = '',
     ylab = '',
     xaxt = 'n', lwd = 3, pch = 20, type = 'b')
mtext('log(C, 2)', 1, 2.75, font = 2, cex = 1.5)
mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=c(seq(0, 0.9, 0.1), 0.95), labels=c("0.0", c(seq(0.1, 0.9, 0.1), 0.95)), font = 2, cex = 1.5)
axis(1, at=log(c(4, 8, 16, 32), 2), labels = log(c(4, 8, 16, 32), 2), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=c("0.0", seq(0.1, 0.9, 0.1), "1.0"), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5)
lines(x = log(c(4, 8, 16, 32), 2), y = consensus,
      lty = 4, lwd = 3, pch = 3, type = 'b')
lines(x = log(c(4, 8, 16, 32), 2), y = neiswanger,
      lty = 5, lwd = 3, pch = 2, type = 'b')
lines(x = log(c(4, 8, 16, 32), 2), y = weierstrass,
      lty = 6, lwd = 3, pch = 1, type = 'b')
legend(x = 2, y = 0.5,
       legend = c('D&C-MCF', 'CMC', 'KDEMC', 'WRS'),
       lwd = c(3, 3, 3),
       pch = c(20, 3, 2, 1),
       lty = c(1, 4, 5, 6),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

######################### TIME

NB_fusion_time <- c(sum(unlist(NB_hc_4$time)),
                    sum(unlist(NB_hc_8$time)),
                    sum(unlist(NB_hc_16$time)),
                    sum(unlist(NB_hc_32$time)))
consensus_time <- c(consensus_mat_4$time,
                    consensus_mat_8$time,
                    consensus_mat_16$time,
                    consensus_mat_32$time)
neiswanger_time <- c(neiswanger_false_4$time,
                     neiswanger_false_8$time,
                     neiswanger_false_16$time,
                     neiswanger_false_32$time)
weierstrass_time <- c(weierstrass_rejection_4$time,
                      weierstrass_rejection_8$time,
                      weierstrass_rejection_16$time,
                      weierstrass_rejection_32$time)

plot(x = log(c(4, 8, 16, 32), 2), y = log(NB_fusion_time, 2),
     ylim = c(-4, 15),
     xlab = '',
     ylab = '',
     yaxt = 'n',
     xaxt = 'n', lwd = 3, pch = 20, type = 'b')
mtext('log(C, 2)', 1, 2.75, font = 2, cex = 1.5)
mtext('log(Time elapsed in seconds)', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=c(seq(0, 0.9, 0.1), 0.95), labels=c("0.0", c(seq(0.1, 0.9, 0.1), 0.95)), font = 2, cex = 1.5)
axis(1, at=log(c(4, 8, 16, 32), 2), labels = log(c(4, 8, 16, 32), 2), font = 2, cex = 1.5)
axis(2, at=seq(-4, 14, 2), labels = seq(-4, 14, 2), font = 2, cex = 1.5)
axis(2, at=seq(-4, 15, 1), labels=rep("", 20), lwd.ticks = 0.5)
lines(x = log(c(4, 8, 16, 32), 2), y = log(consensus_time, 2),
      lty = 4, lwd = 3, pch = 3, type = 'b')
lines(x = log(c(4, 8, 16, 32), 2), y = log(neiswanger_time, 2),
      lty = 5, lwd = 3, pch = 2, type = 'b')
lines(x = log(c(4, 8, 16, 32), 2), y = log(weierstrass_time, 2),
      lty = 6, lwd = 3, pch = 1, type = 'b')
legend(x = 2, y = 15,
       legend = c('D&C-MCF', 'CMC', 'KDEMC', 'WRS'),
       lwd = c(3, 3, 3),
       pch = c(20, 3, 2, 1),
       lty = c(1, 4, 5, 6),
       cex = 1.25,
       text.font = 2,
       bty = 'n')
