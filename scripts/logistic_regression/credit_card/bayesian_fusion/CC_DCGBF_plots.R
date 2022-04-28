library(DCFusion)
load('CC32_DCGBF.RData')
load('CC64_DCGBF.RData')
load('CC128_DCGBF.RData')
load('CC256_DCGBF.RData')

balanced <- list('reg' = c(integrated_abs_distance(full_posterior,
                                                   balanced_C32$reg$particles$y_samples),
                           integrated_abs_distance(full_posterior,
                                                   balanced_C64$reg$particles$y_samples),
                           integrated_abs_distance(full_posterior,
                                                   balanced_C128$reg$particles$y_samples)),
                 'adaptive' = c(integrated_abs_distance(full_posterior,
                                                        balanced_C32$adaptive$particles$y_samples),
                                integrated_abs_distance(full_posterior,
                                                        balanced_C64$adaptive$particles$y_samples),
                                integrated_abs_distance(full_posterior,
                                                        balanced_C128$adaptive$particles$y_samples)))
NB_fusion <- c(integrated_abs_distance(full_posterior,
                                       NB_hc_32$particles$y_samples),
               integrated_abs_distance(full_posterior,
                                       NB_hc_64$particles$y_samples),
               integrated_abs_distance(full_posterior,
                                       NB_hc_128$particles$y_samples))
# consensus <- c(integrated_abs_distance(full_posterior,
#                                        consensus_mat_32$samples),
#                integrated_abs_distance(full_posterior,
#                                        consensus_mat_64$samples),
#                integrated_abs_distance(full_posterior,
#                                        consensus_mat_128$samples))
# neiswanger <- c(integrated_abs_distance(full_posterior,
#                                         neiswanger_false_32$samples),
#                 integrated_abs_distance(full_posterior,
#                                         neiswanger_false_64$samples),
#                 integrated_abs_distance(full_posterior,
#                                         neiswanger_false_128$samples))
# weierstrass <- c(integrated_abs_distance(full_posterior,
#                                          weierstrass_rejection_32$samples),
#                  integrated_abs_distance(full_posterior,
#                                          weierstrass_rejection_64$samples),
#                  integrated_abs_distance(full_posterior,
#                                          weierstrass_rejection_128$samples))

plot(x = log(c(32, 64, 128), 2), y = balanced$reg,
     ylim = c(0, 0.2),
     xlab = '',
     ylab = '',
     xaxt = 'n', yaxt = 'n', lwd = 3, pch = 1, type = 'b')
mtext('log(C, 2)', 1, 2.75, font = 2, cex = 1.5)
mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=c(seq(0, 0.9, 0.1), 0.95), labels=c("0.0", c(seq(0.1, 0.9, 0.1), 0.95)), font = 2, cex = 1.5)
axis(1, at=log(c(4, 8, 16, 32), 2), labels = log(c(4, 8, 16, 32), 2), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=c("0.0", seq(0.1, 0.9, 0.1), "1.0"), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.05), labels=rep("", 21), lwd.ticks = 0.5)
lines(x = log(c(32, 64, 128), 2), y = balanced$adaptive,
      lty = 2, lwd = 3, type = 'b', pch = 2)
lines(x = log(c(32, 64, 128), 2), y = NB_fusion,
      lty = 5, lwd = 3, type = 'b', pch = 3)
# lines(x = log(c(32, 64, 128), 2), y = consensus,
#       lty = 6, lwd = 3, type = 'b', pch = 6, col = 'red')
# lines(x = log(c(32, 64, 128), 2), y = neiswanger,
#       lty = 7, lwd = 3, type = 'b', pch = 7, col = 'red')
# lines(x = log(c(32, 64, 128), 2), y = weierstrass,
#       lty = 8, lwd = 3, type = 'b', pch = 8, col = 'red')
legend(x = 5, y = 0.2,
       legend = c('D&C-GBF (reg)', 'D&C-GBF (adaptive)', 'DC-MCF'),
       lwd = c(3, 3, 3),
       lty = c(1, 2, 3),
       pch = c(1, 2, 3),
       col = c(rep('black', 5), rep('red', 3)),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

######################### TIME

balanced_time <-  list('reg' = c(sum(unlist(balanced_C32$reg$time)),
                                 sum(unlist(balanced_C64$reg$time)),
                                 sum(unlist(balanced_C128$reg$time))),
                       'adaptive' = c(sum(unlist(balanced_C32$adaptive$time)),
                                      sum(unlist(balanced_C64$adaptive$time)),
                                      sum(unlist(balanced_C128$adaptive$time))))
NB_fusion_time <- c(sum(unlist(NB_hc_32$time)),
                    sum(unlist(NB_hc_64$time)),
                    sum(unlist(NB_hc_128$time)))
# consensus_time <- c(consensus_mat_32$time,
#                     consensus_mat_64$time,
#                     consensus_mat_128$time)
# neiswanger_time <- c(neiswanger_false_32$time,
#                      neiswanger_false_64$time,
#                      neiswanger_false_128$time)
# weierstrass_time <- c(weierstrass_rejection_32$time,
#                       weierstrass_rejection_64$time,
#                       weierstrass_rejection_128$time)

plot(x = log(c(32, 64, 128), 2), y = log(balanced_time$reg, 2),
     ylim = c(4, 16),
     xlab = '',
     ylab = '',
     xaxt = 'n', lwd = 3, pch = 1, type = 'b', xaxt = 'n', yaxt = 'n')
mtext('log(C, 2)', 1, 2.75, font = 2, cex = 1.5)
mtext('log(Time elapsed in seconds, 2)', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=c(seq(0, 0.9, 0.1), 0.95), labels=c("0.0", c(seq(0.1, 0.9, 0.1), 0.95)), font = 2, cex = 1.5)
axis(1, at=log(c(32, 64, 128), 2), labels = log(c(32, 64, 128), 2), font = 2, cex = 1.5)
axis(2, at=seq(-4, 16, 2), labels = seq(-4, 16, 2), font = 2, cex = 1.5)
axis(2, at=seq(-4, 16, 1), labels=rep("", 21), lwd.ticks = 0.5)
lines(x = log(c(32, 64, 128), 2), y = log(balanced_time$adaptive, 2),
      lty = 2, lwd = 3, type = 'b', pch = 2)
lines(x = log(c(32, 64, 128), 2), y = log(NB_fusion_time, 2),
      lty = 3, lwd = 3, type = 'b', pch = 3)
# lines(x = log(c(32, 64, 128), 2), y = log(consensus_time, 2),
#       lty = 6, lwd = 3, type = 'b', pch = 6, col = 'red')
# lines(x = log(c(32, 64, 128), 2), y = log(neiswanger_time, 2),
#       lty = 7, lwd = 3, type = 'b', pch = 7, col = 'red')
# lines(x = log(c(32, 64, 128), 2), y = log(weierstrass_time, 2),
#       lty = 8, lwd = 3, type = 'b', pch = 8, col = 'red')
legend(x = 5, y = 16,
       legend = c('D&C-GBF (reg)', 'D&C-GBF (adaptive)', 'DC-MCF'),
       lwd = 3,
       lty = 1:3,
       pch = 1:3,
       col = 'black',
       cex = 1.25,
       text.font = 2,
       bty = 'n')
# legend(x = 3.5, y = 16,
#        legend = c('CMC', 'KDEMC', 'WRS'),
#        col = 'red',
#        lwd = c(3, 3, 3, 3),
#        lty = c(4, 2, 3), 
#        pch = c(6, 7, 8),
#        cex = 1.25,
#        text.font = 2,
#        bty = 'n')
