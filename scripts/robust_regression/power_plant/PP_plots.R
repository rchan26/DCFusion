library(DCFusion)
load('PP4.RData')
load('PP8.RData')
load('PP16.RData')
load('PP32.RData')
load('PP64.RData')
load('PP128.RData')

# GBF <- list('reg' = c(integrated_abs_distance(full_posterior,
#                                               GBF_4$reg$particles$y_samples),
#                       integrated_abs_distance(full_posterior,
#                                               GBF_8$reg$particles$y_samples)),
#             'adaptive' = c(integrated_abs_distance(full_posterior,
#                                                    GBF_4$adaptive$particles$y_samples),
#                            integrated_abs_distance(full_posterior,
#                                                    GBF_8$adaptive$particles$y_samples)))
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
NB_fusion <- c(integrated_abs_distance(full_posterior,
                                       NB_hc_4$particles$y_samples),
               integrated_abs_distance(full_posterior,
                                       NB_hc_8$particles$y_samples),
               integrated_abs_distance(full_posterior,
                                       NB_hc_16$particles$y_samples),
               integrated_abs_distance(full_posterior,
                                       NB_hc_32$particles$y_samples),
               integrated_abs_distance(full_posterior,
                                       NB_hc_64$particles$y_samples),
               integrated_abs_distance(full_posterior,
                                       NB_hc_128$particles$y_samples))
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
     ylim = c(0, 0.5),
     xlab = '',
     ylab = '',
     xaxt = 'n', lty = 2, lwd = 3, pch = 4, type = 'b')
mtext('log(C, 2)', 1, 2.75, font = 2, cex = 1.5)
mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=c(seq(0, 0.9, 0.1), 0.95), labels=c("0.0", c(seq(0.1, 0.9, 0.1), 0.95)), font = 2, cex = 1.5)
axis(1, at=log(c(4, 8, 16, 32, 64, 128), 2), labels = log(c(4, 8, 16, 32, 64, 128), 2), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=c("0.0", seq(0.1, 0.9, 0.1), "1.0"), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5)
lines(x = log(c(4, 8, 16, 32, 64, 128), 2), y = balanced$reg,
      lty = 3, lwd = 3, type = 'b', pch = 5)
# lines(x = log(c(4, 8), 2), y = GBF$adaptive,
#       lty = 1, lwd = 3, type = 'b', pch = 1)
lines(x = log(c(4, 8, 16, 32, 64, 128), 2), y = NB_fusion,
      lty = 1, lwd = 3, type = 'b', pch = 20)
lines(x = log(c(4, 8, 16, 32, 64, 128), 2), y = consensus,
      lty = 4, lwd = 3, type = 'b', pch = 3, col = 'red')
lines(x = log(c(4, 8, 16, 32, 64, 128), 2), y = neiswanger,
      lty = 5, lwd = 3, type = 'b', pch = 2, col = 'red')
lines(x = log(c(4, 8, 16, 32, 64, 128), 2), y = weierstrass,
      lty = 6, lwd = 3, type = 'b', pch = 1, col = 'red')
legend(x = 2, y = 0.5,
       legend = c('D&C-GBF (regular mesh)',
                  'D&C-GBF (adaptive mesh)',
                  'D&C-MCF',
                  'CMC',
                  'KDEMC',
                  'WRS'),
       lwd = rep(3, 6),
       lty = c(3,2,1,4,5,6),
       pch = c(5,4,20,3,2,1),
       col = c(rep('black', 3), rep('red', 3)),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

######################### TIME

# GBF_time <- list('reg' = c(GBF_4$reg$time[[1]],
#                            GBF_8$reg$time[[1]]),
#                  'adaptive' = c(GBF_4$adaptive$time[[1]],
#                                 GBF_8$adaptive$time[[1]]))
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
NB_fusion_time <- c(sum(unlist(NB_hc_4$time)),
                    sum(unlist(NB_hc_8$time)),
                    sum(unlist(NB_hc_16$time)),
                    sum(unlist(NB_hc_32$time)),
                    sum(unlist(NB_hc_64$time)),
                    sum(unlist(NB_hc_128$time)))
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
     ylim = c(-2, 20),
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
# lines(x = log(c(4, 8), 2), y = log(GBF_time$adaptive, 2),
#       lty = 1, lwd = 3, type = 'b', pch = 1)
lines(x = log(c(4, 8, 16, 32, 64, 128), 2), y = log(NB_fusion_time, 2),
      lty = 1, lwd = 3, type = 'b', pch = 20)
lines(x = log(c(4, 8, 16, 32, 64, 128), 2), y = log(consensus_time, 2),
      lty = 4, lwd = 3, type = 'b', pch = 3, col = 'red')
lines(x = log(c(4, 8, 16, 32, 64, 128), 2), y = log(neiswanger_time, 2),
      lty = 5, lwd = 3, type = 'b', pch = 2, col = 'red')
lines(x = log(c(4, 8, 16, 32, 64, 128), 2), y = log(weierstrass_time, 2),
      lty = 6, lwd = 3, type = 'b', pch = 1, col = 'red')
legend(x = 2, y = 20,
       legend = c('D&C-GBF (regular mesh)',
                  'D&C-GBF (adaptive mesh)',
                  'D&C-MCF',
                  'CMC',
                  'KDEMC',
                  'WRS'),
       lwd = rep(3, 6),
       lty = c(3,2,1,4,5,6),
       pch = c(5,4,20,3,2,1),
       col = c(rep('black', 3), rep('red', 3)),
       cex = 1.25,
       text.font = 2,
       bty = 'n')
# legend(x = 4.5, y = 20,
#        legend = c('CMC',
#                   'KDEMC',
#                   'WRS'),
#        lwd = rep(3, 3),
#        lty = c(4,5,6),
#        pch = c(3,2,1),
#        col = rep('red', 3),
#        cex = 1.25,
#        text.font = 2,
#        bty = 'n')
