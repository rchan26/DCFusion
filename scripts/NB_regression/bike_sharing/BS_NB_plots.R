library(DCFusion)
load('BS4_NB.RData')
load('BS8_NB.RData')
load('BS16_NB.RData')
load('BS32_NB.RData')
load('BS64_NB.RData')
load('BS128_NB.RData')

balanced <- list('reg' = c(integrated_abs_distance(full_posterior,
                                                   balanced_C4$reg$particles$y_samples),
                           integrated_abs_distance(full_posterior,
                                                   balanced_C8$reg$particles$y_samples),
                           integrated_abs_distance(full_posterior,
                                                   balanced_C16$reg$particles$y_samples),
                           integrated_abs_distance(full_posterior,
                                                   balanced_C32$reg$particles$y_samples),
                           integrated_abs_distance(full_posterior,
                                                   balanced_C64$reg$particles$y_samples),
                           integrated_abs_distance(full_posterior,
                                                   balanced_C128$reg$particles$y_samples)),
                 'adaptive' = c(integrated_abs_distance(full_posterior,
                                                        balanced_C4$adaptive$particles$y_samples),
                                integrated_abs_distance(full_posterior,
                                                        balanced_C8$adaptive$particles$y_samples),
                                integrated_abs_distance(full_posterior,
                                                        balanced_C16$adaptive$particles$y_samples),
                                integrated_abs_distance(full_posterior,
                                                        balanced_C32$adaptive$particles$y_samples),
                                integrated_abs_distance(full_posterior,
                                                        balanced_C64$adaptive$particles$y_samples),
                                integrated_abs_distance(full_posterior,
                                                        balanced_C128$adaptive$particles$y_samples)))
consensus <- c(integrated_abs_distance(full_posterior,
                                       consensus_mat_4$samples),
               integrated_abs_distance(full_posterior,
                                       consensus_mat_8$samples),
               integrated_abs_distance(full_posterior,
                                       consensus_mat_16$samples),
               integrated_abs_distance(full_posterior,
                                       consensus_mat_32$samples),
               integrated_abs_distance(full_posterior,
                                       consensus_mat_64$samples),
               integrated_abs_distance(full_posterior,
                                       consensus_mat_128$samples))
neiswanger <- c(integrated_abs_distance(full_posterior,
                                        neiswanger_false_4$samples),
                integrated_abs_distance(full_posterior,
                                        neiswanger_false_8$samples),
                integrated_abs_distance(full_posterior,
                                        neiswanger_false_16$samples),
                integrated_abs_distance(full_posterior,
                                        neiswanger_false_32$samples),
                integrated_abs_distance(full_posterior,
                                        neiswanger_false_64$samples),
                integrated_abs_distance(full_posterior,
                                        neiswanger_false_128$samples))
weierstrass <- c(integrated_abs_distance(full_posterior,
                                         weierstrass_rejection_4$samples),
                 integrated_abs_distance(full_posterior,
                                         weierstrass_rejection_8$samples),
                 integrated_abs_distance(full_posterior,
                                         weierstrass_rejection_16$samples),
                 integrated_abs_distance(full_posterior,
                                         weierstrass_rejection_32$samples),
                 integrated_abs_distance(full_posterior,
                                         weierstrass_rejection_64$samples),
                 integrated_abs_distance(full_posterior,
                                         weierstrass_rejection_128$samples))

plot(x = log(c(4, 8, 16, 32, 64, 128), 2), y = balanced$adaptive,
     ylim = c(0, 0.3),
     xlab = '',
     ylab = '',
     xaxt = 'n', yaxt = 'n', lty = 2, lwd = 3, pch = 4, type = 'b')
mtext('log(C, 2)', 1, 2.75, font = 2, cex = 1.5)
mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=c(seq(0, 0.9, 0.1), 0.95), labels=c("0.0", c(seq(0.1, 0.9, 0.1), 0.95)), font = 2, cex = 1.5)
axis(1, at=log(c(4, 8, 16, 32, 64, 128), 2), labels = log(c(4, 8, 16, 32, 64, 128), 2), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=c("0.0", seq(0.1, 0.9, 0.1), "1.0"), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.05), labels=rep("", 21), lwd.ticks = 0.5)
lines(x = log(c(4, 8, 16, 32, 64, 128), 2), y = balanced$reg,
      lty = 3, lwd = 3, type = 'b', pch = 5)
lines(x = log(c(4, 8, 16, 32, 64, 128), 2), y = consensus,
      lty = 4, lwd = 3, type = 'b', pch = 3, col = 'red')
lines(x = log(c(4, 8, 16, 32, 64, 128), 2), y = neiswanger,
      lty = 5, lwd = 3, type = 'b', pch = 2, col = 'red')
lines(x = log(c(4, 8, 16, 32, 64, 128), 2), y = weierstrass,
      lty = 6, lwd = 3, type = 'b', pch = 1, col = 'red')
legend(x = 2, y = 0.3,
       legend = c('D&C-GBF (regular mesh)',
                  'D&C-GBF (adaptive mesh)',
                  'CMC',
                  'KDEMC',
                  'WRS'),
       lwd = rep(3, 6),
       lty = c(3,2,4,5,6),
       pch = c(5,4,3,2,1),
       col = c(rep('black', 2), rep('red', 3)),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

######################### TIME

balanced_time <-  list('reg' = c(sum(unlist(balanced_C4$reg$time)),
                                 sum(unlist(balanced_C8$reg$time)),
                                 sum(unlist(balanced_C16$reg$time)),
                                 sum(unlist(balanced_C32$reg$time)),
                                 sum(unlist(balanced_C64$reg$time)),
                                 sum(unlist(balanced_C128$reg$time))),
                       'adaptive' = c(sum(unlist(balanced_C4$adaptive$time)),
                                      sum(unlist(balanced_C8$adaptive$time)),
                                      sum(unlist(balanced_C16$adaptive$time)),
                                      sum(unlist(balanced_C32$adaptive$time)),
                                      sum(unlist(balanced_C64$adaptive$time)),
                                      sum(unlist(balanced_C128$adaptive$time))))
consensus_time <- c(consensus_mat_4$time,
                    consensus_mat_8$time,
                    consensus_mat_16$time,
                    consensus_mat_32$time,
                    consensus_mat_64$time,
                    consensus_mat_128$time)
neiswanger_time <- c(neiswanger_false_4$time,
                     neiswanger_false_8$time,
                     neiswanger_false_16$time,
                     neiswanger_false_32$time,
                     neiswanger_false_64$time,
                     neiswanger_false_128$time)
weierstrass_time <- c(weierstrass_rejection_4$time,
                      weierstrass_rejection_8$time,
                      weierstrass_rejection_16$time,
                      weierstrass_rejection_32$time,
                      weierstrass_rejection_64$time,
                      weierstrass_rejection_128$time)

plot(x = log(c(4, 8, 16, 32, 64, 128), 2), y = log(balanced_time$adaptive, 2),
     ylim = c(-2, 18),
     xlab = '',
     ylab = '',
     yaxt = 'n',
     xaxt = 'n', lty = 2, lwd = 3, pch = 4, type = 'b')
mtext('log(C, 2)', 1, 2.75, font = 2, cex = 1.5)
mtext('log(Time elapsed in seconds, 2)', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=c(seq(0, 0.9, 0.1), 0.95), labels=c("0.0", c(seq(0.1, 0.9, 0.1), 0.95)), font = 2, cex = 1.5)
axis(1, at=log(c(4, 8, 16, 32, 64, 128), 2), labels = log(c(4, 8, 16, 32, 64, 128), 2), font = 2, cex = 1.5)
axis(2, at=seq(-4, 20, 2), labels = seq(-4, 20, 2), font = 2, cex = 1.5)
axis(2, at=seq(-4, 20, 1), labels=rep("", 25), lwd.ticks = 0.5)
lines(x = log(c(4, 8, 16, 32, 64, 128), 2), y = log(balanced_time$reg, 2),
      lty = 3, lwd = 3, type = 'b', pch = 5)
lines(x = log(c(4, 8, 16, 32, 64, 128), 2), y = log(consensus_time, 2),
      lty = 4, lwd = 3, type = 'b', pch = 3, col = 'red')
lines(x = log(c(4, 8, 16, 32, 64, 128), 2), y = log(neiswanger_time, 2),
      lty = 5, lwd = 3, type = 'b', pch = 2, col = 'red')
lines(x = log(c(4, 8, 16, 32, 64, 128), 2), y = log(weierstrass_time, 2),
      lty = 6, lwd = 3, type = 'b', pch = 1, col = 'red')
legend(x = 2, y = 18,
       legend = c('D&C-GBF (regular mesh)',
                  'D&C-GBF (adaptive mesh)',
                  'CMC',
                  'KDEMC',
                  'WRS'),
       lwd = rep(3, 6),
       lty = c(3,2,4,5,6),
       pch = c(5,4,3,2,1),
       col = c(rep('black', 2), rep('red', 3)),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

##### ISBA poster #####

plot(x = log(c(4, 8, 16, 32, 64, 128), 2), y = balanced$adaptive,
     ylim = c(0, 0.4),
     xlab = '',
     ylab = '',
     xaxt = 'n', lty = 1, lwd = 3, pch = 1, type = 'b')
mtext('log(C, 2)', 1, 2.75, font = 2, cex = 1.5)
mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=c(seq(0, 0.9, 0.1), 0.95), labels=c("0.0", c(seq(0.1, 0.9, 0.1), 0.95)), font = 2, cex = 1.5)
axis(1, at=log(c(4, 8, 16, 32, 64, 128), 2), labels = log(c(4, 8, 16, 32, 64, 128), 2), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=c("0.0", seq(0.1, 0.9, 0.1), "1.0"), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5)
lines(x = log(c(4, 8, 16, 32, 64, 128), 2), y = consensus,
      lty = 4, lwd = 3, type = 'b', pch = 4, col = 'red')
lines(x = log(c(4, 8, 16, 32, 64, 128), 2), y = neiswanger,
      lty = 5, lwd = 3, type = 'b', pch = 5, col = 'red')
lines(x = log(c(4, 8, 16, 32, 64, 128), 2), y = weierstrass,
      lty = 6, lwd = 3, type = 'b', pch = 6, col = 'red')
legend(x = 2, y = 0.4,
       legend = c('D&C-GBF',
                  'CMC',
                  'KDEMC',
                  'WRS'),
       lwd = rep(3, 5),
       lty = c(1,4,5,6),
       pch = c(1,4,5,6),
       col = c(rep('black', 1), rep('red', 3)),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

plot(x = log(c(4, 8, 16, 32, 64, 128), 2), y = log(balanced_time$adaptive, 2),
     ylim = c(-4, 22),
     xlab = '',
     ylab = '',
     yaxt = 'n',
     xaxt = 'n', lty = 1, lwd = 3, pch = 1, type = 'b')
mtext('log(C, 2)', 1, 2.75, font = 2, cex = 1.5)
mtext('log(Time elapsed in seconds, 2)', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=c(seq(0, 0.9, 0.1), 0.95), labels=c("0.0", c(seq(0.1, 0.9, 0.1), 0.95)), font = 2, cex = 1.5)
axis(1, at=log(c(4, 8, 16, 32, 64, 128), 2), labels = log(c(4, 8, 16, 32, 64, 128), 2), font = 2, cex = 1.5)
axis(2, at=seq(-4, 22, 4), labels = seq(-4, 22, 4), font = 2, cex = 1.5)
axis(2, at=seq(-4, 22, 1), labels=rep("", 27), lwd.ticks = 0.5)
lines(x = log(c(4, 8, 16, 32, 64, 128), 2), y = log(consensus_time, 2),
      lty = 4, lwd = 3, type = 'b', pch = 4, col = 'red')
lines(x = log(c(4, 8, 16, 32, 64, 128), 2), y = log(neiswanger_time, 2),
      lty = 5, lwd = 3, type = 'b', pch = 5, col = 'red')
lines(x = log(c(4, 8, 16, 32, 64, 128), 2), y = log(weierstrass_time, 2),
      lty = 6, lwd = 3, type = 'b', pch = 6, col = 'red')
legend(x = 2, y = 22,
       legend = c('D&C-GBF',
                  'CMC',
                  'KDEMC',
                  'WRS'),
       lwd = rep(3, 5),
       lty = c(1,4,5,6),
       pch = c(1,4,5,6),
       col = c(rep('black', 1), rep('red', 3)),
       cex = 1.25,
       text.font = 2,
       bty = 'n')
