library(DCFusion)
load('PP4.RData')
load('PP8.RData')
load('PP16.RData')
load('PP32.RData')

GBF <- list('reg' = c(integrated_abs_distance(full_posterior,
                                              GBF_4$reg$particles$y_samples),
                      integrated_abs_distance(full_posterior,
                                              GBF_8$reg$particles$y_samples)),
            'adaptive' = c(integrated_abs_distance(full_posterior,
                                                   GBF_4$adaptive$particles$y_samples),
                           integrated_abs_distance(full_posterior,
                                                   GBF_8$adaptive$particles$y_samples)))
balanced <- list('reg' = c(integrated_abs_distance(full_posterior,
                                                   balanced_C4$reg$particles$y_samples),
                           integrated_abs_distance(full_posterior,
                                                   balanced_C8$reg$particles$y_samples),
                           integrated_abs_distance(full_posterior,
                                                   balanced_C16$reg$particles$y_samples),
                           integrated_abs_distance(full_posterior,
                                                   balanced_C32$reg$particles$y_samples)),
                 'adaptive' = c(integrated_abs_distance(full_posterior,
                                                        balanced_C4$adaptive$particles$y_samples),
                                integrated_abs_distance(full_posterior,
                                                        balanced_C8$adaptive$particles$y_samples),
                                integrated_abs_distance(full_posterior,
                                                        balanced_C16$adaptive$particles$y_samples),
                                integrated_abs_distance(full_posterior,
                                                        balanced_C32$adaptive$particles$y_samples)))
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

plot(x = c(4, 8, 16, 32), y = balanced$reg,
     ylim = c(0, 1),
     xlab = '',
     ylab = '',
     xaxt = 'n', lwd = 3, pch = 1, type = 'b')
mtext('Number of sub-posteriors (C)', 1, 2.75, font = 2, cex = 1.5)
mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=c(seq(0, 0.9, 0.1), 0.95), labels=c("0.0", c(seq(0.1, 0.9, 0.1), 0.95)), font = 2, cex = 1.5)
axis(1, at=c(4, 8, 16, 32), labels = c(4, 8, 16, 32), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=c("0.0", seq(0.1, 0.9, 0.1), "1.0"), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5)
lines(x = c(4, 8, 16, 32), y = balanced$adaptive,
      lty = 2, lwd = 3, type = 'b', pch = 4)
lines(x = c(4, 8), y = GBF$reg,
      lty = 3, lwd = 3, type = 'b', pch = 3)
lines(x = c(4, 8), y = GBF$adaptive,
      lty = 4, lwd = 3, type = 'b', pch = 2)
lines(x = c(4, 8, 16, 32), y = NB_fusion,
      lty = 5, lwd = 3, type = 'b', pch = 5)
lines(x = c(4, 8, 16, 32), y = consensus,
      lty = 6, lwd = 3, type = 'b', pch = 6, col = 'red')
lines(x = c(4, 8, 16, 32), y = neiswanger,
      lty = 7, lwd = 3, type = 'b', pch = 7, col = 'red')
lines(x = c(4, 8, 16, 32), y = weierstrass,
      lty = 8, lwd = 3, type = 'b', pch = 8, col = 'red')
legend(x = 4, y = 1,
       legend = c('D&C-GBF (reg)', 'D&C-GBF (adaptive)', 'GBF (reg)', 'GBF (adaptive)', 'DC-MCF', 'CMC', 'KDEMC', 'WRS'),
       lwd = c(3, 3, 3, 3, 3, 3, 3, 3),
       lty = c(1, 2, 3, 4, 5, 6, 7, 8),
       pch = c(1, 4, 3, 2, 5, 6, 7, 8),
       col = c(rep('black', 5), rep('red', 3)),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

######################### TIME

GBF_time <- list('reg' = c(GBF_4$reg$time[[1]],
                           GBF_8$reg$time[[1]]),
                 'adaptive' = c(GBF_4$adaptive$time[[1]],
                                GBF_8$adaptive$time[[1]]))
balanced_time <-  list('reg' = c(sum(unlist(balanced_C4$reg$time)),
                                 sum(unlist(balanced_C8$reg$time)),
                                 sum(unlist(balanced_C16$reg$time)),
                                 sum(unlist(balanced_C32$reg$time))),
                       'adaptive' = c(sum(unlist(balanced_C4$adaptive$time)),
                                      sum(unlist(balanced_C8$adaptive$time)),
                                      sum(unlist(balanced_C16$adaptive$time)),
                                      sum(unlist(balanced_C32$adaptive$time))))
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

plot(x = c(4, 8, 16, 32), y = log(balanced_time$reg),
     ylim = c(-4, 16),
     xlab = '',
     ylab = '',
     xaxt = 'n', lwd = 3, pch = 1, type = 'b', xaxt = 'n', yaxt = 'n')
mtext('Number of sub-posteriors (C)', 1, 2.75, font = 2, cex = 1.5)
mtext('log(Time elapsed in seconds)', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=c(seq(0, 0.9, 0.1), 0.95), labels=c("0.0", c(seq(0.1, 0.9, 0.1), 0.95)), font = 2, cex = 1.5)
axis(1, at=c(4, 8, 16, 32), labels = c(4, 8, 16, 32), font = 2, cex = 1.5)
axis(2, at=seq(-4, 16, 2), labels = seq(-4, 16, 2), font = 2, cex = 1.5)
axis(2, at=seq(-4, 16, 1), labels=rep("", 21), lwd.ticks = 0.5)
lines(x = c(4, 8, 16, 32), y = log(balanced_time$adaptive),
      lty = 2, lwd = 3, type = 'b', pch = 4)
lines(x = c(4, 8), y = log(GBF_time$reg),
      lty = 3, lwd = 3, type = 'b', pch = 3)
lines(x = c(4, 8), y = log(GBF_time$adaptive),
      lty = 4, lwd = 3, type = 'b', pch = 2)
lines(x = c(4, 8, 16, 32), y = log(NB_fusion_time),
      lty = 5, lwd = 3, type = 'b', pch = 5)
lines(x = c(4, 8, 16, 32), y = log(consensus_time),
      lty = 6, lwd = 3, type = 'b', pch = 6, col = 'red')
lines(x = c(4, 8, 16, 32), y = log(neiswanger_time),
      lty = 7, lwd = 3, type = 'b', pch = 7, col = 'red')
lines(x = c(4, 8, 16, 32), y = log(weierstrass_time),
      lty = 8, lwd = 3, type = 'b', pch = 8, col = 'red')
legend(x = 4, y = 16.5,
       legend = c('D&C-GBF (reg)', 'D&C-GBF (adaptive)', 'GBF (reg)', 'GBF (adaptive)', 'DC-MCF'),
       lwd = c(3, 3, 3, 3, 3),
       lty = c(1, 2, 3, 4, 5),
       pch = c(1, 4, 3, 2, 5),
       col = 'black',
       cex = 1.25,
       text.font = 2,
       bty = 'n')
legend(x = 16, y = 16.5,
       legend = c('CMC', 'KDEMC', 'WRS'),
       col = 'red',
       lwd = c(3, 3, 3, 3),
       lty = c(4, 2, 3), 
       pch = c(6, 7, 8),
       cex = 1.25,
       text.font = 2,
       bty = 'n')
