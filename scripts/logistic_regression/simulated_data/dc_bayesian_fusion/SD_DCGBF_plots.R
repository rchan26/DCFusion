library(DCFusion)
load('SD4_DCGBF.RData')
load('SD8_DCGBF.RData')
load('SD16_DCGBF.RData')
load('SD32_DCGBF.RData')
load('SD64_DCGBF.RData')

GBF <- list('adaptive' = c(integrated_abs_distance(full_posterior,
                                                   GBF_4$adaptive$particles$y_samples),
                           integrated_abs_distance(full_posterior,
                                                   GBF_8$adaptive$particles$y_samples),
                           integrated_abs_distance(full_posterior,
                                                   GBF_16$adaptive$particles$y_samples),
                           integrated_abs_distance(full_posterior,
                                                   GBF_32$adaptive$particles$y_samples)))
balanced <- list('reg' = c(integrated_abs_distance(full_posterior,
                                                   balanced_C4$reg$particles$y_samples),
                           integrated_abs_distance(full_posterior,
                                                   balanced_C8$reg$particles$y_samples),
                           integrated_abs_distance(full_posterior,
                                                   balanced_C16$reg$particles$y_samples),
                           integrated_abs_distance(full_posterior,
                                                   balanced_C32$reg$particles$y_samples),
                           integrated_abs_distance(full_posterior,
                                                   balanced_C64$reg$particles$y_samples)),
                 'adaptive' = c(integrated_abs_distance(full_posterior,
                                                        balanced_C4$adaptive$particles$y_samples),
                                integrated_abs_distance(full_posterior,
                                                        balanced_C8$adaptive$particles$y_samples),
                                integrated_abs_distance(full_posterior,
                                                        balanced_C16$adaptive$particles$y_samples),
                                integrated_abs_distance(full_posterior,
                                                        balanced_C32$adaptive$particles$y_samples),
                                integrated_abs_distance(full_posterior,
                                                        balanced_C64$adaptive$particles$y_samples)))
consensus <- c(integrated_abs_distance(full_posterior,
                                       consensus_mat_4$samples),
               integrated_abs_distance(full_posterior,
                                       consensus_mat_8$samples),
               integrated_abs_distance(full_posterior,
                                       consensus_mat_16$samples),
               integrated_abs_distance(full_posterior,
                                       consensus_mat_32$samples),
               integrated_abs_distance(full_posterior,
                                       consensus_mat_64$samples))
neiswanger <- c(integrated_abs_distance(full_posterior,
                                        neiswanger_false_4$samples),
                integrated_abs_distance(full_posterior,
                                        neiswanger_false_8$samples),
                integrated_abs_distance(full_posterior,
                                        neiswanger_false_16$samples),
                integrated_abs_distance(full_posterior,
                                        neiswanger_false_32$samples),
                integrated_abs_distance(full_posterior,
                                        neiswanger_false_64$samples))
weierstrass <- c(integrated_abs_distance(full_posterior,
                                         weierstrass_rejection_4$samples),
                 integrated_abs_distance(full_posterior,
                                         weierstrass_rejection_8$samples),
                 integrated_abs_distance(full_posterior,
                                         weierstrass_rejection_16$samples),
                 integrated_abs_distance(full_posterior,
                                         weierstrass_rejection_32$samples),
                 integrated_abs_distance(full_posterior,
                                         weierstrass_rejection_64$samples))

plot(x = log(c(4, 8, 16, 32, 64), 2), y = balanced$adaptive,
     ylim = c(0, 0.6),
     xlab = '',
     ylab = '',
     xaxt = 'n', lty = 3, lwd = 3, pch = 3, type = 'b')
mtext('log(C, 2)', 1, 2.75, font = 2, cex = 1.5)
mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=c(seq(0, 0.9, 0.1), 0.95), labels=c("0.0", c(seq(0.1, 0.9, 0.1), 0.95)), font = 2, cex = 1.5)
axis(1, at=log(c(4, 8, 16, 32, 64), 2), labels = log(c(4, 8, 16, 32, 64), 2), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=c("0.0", seq(0.1, 0.9, 0.1), "1.0"), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5)
lines(x = log(c(4, 8, 16, 32, 64), 2), y = balanced$reg,
      lty = 2, lwd = 3, type = 'b', pch = 2)
lines(x = log(c(4, 8, 16, 32), 2), y = GBF$adaptive,
      lty = 1, lwd = 3, type = 'b', pch = 1)
lines(x = log(c(4, 8, 16, 32, 64), 2), y = consensus,
      lty = 4, lwd = 3, type = 'b', pch = 4, col = 'red')
lines(x = log(c(4, 8, 16, 32, 64), 2), y = neiswanger,
      lty = 5, lwd = 3, type = 'b', pch = 5, col = 'red')
lines(x = log(c(4, 8, 16, 32, 64), 2), y = weierstrass,
      lty = 6, lwd = 3, type = 'b', pch = 6, col = 'red')
lines(x = log(c(4, 8, 16, 32), 2), y = NB_fusion,
      lty = 7, lwd = 3, type = 'b', pch = 7, col = 'red')
legend(x = 2, y = 0.6,
       legend = c('GBF (adaptive mesh)',
                  'D&C-GBF (regular mesh)',
                  'D&C-GBF (adaptive mesh)',
                  'D&C-MCF',
                  'CMC',
                  'KDEMC',
                  'WRS'),
       lwd = rep(3, 7),
       lty = c(1,2,3,7,4,5,6),
       pch = c(1,2,3,7,4,5,6),
       col = c(rep('black', 3), rep('red', 4)),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

######################### TIME

GBF_time <- list('adaptive' = c(GBF_4$adaptive$time[[1]],
                                GBF_8$adaptive$time[[1]],
                                GBF_16$adaptive$time[[1]],
                                GBF_32$adaptive$time[[1]]))
balanced_time <-  list('reg' = c(sum(unlist(balanced_C4$reg$time)),
                                 sum(unlist(balanced_C8$reg$time)),
                                 sum(unlist(balanced_C16$reg$time)),
                                 sum(unlist(balanced_C32$reg$time)),
                                 sum(unlist(balanced_C64$reg$time))),
                       'adaptive' = c(sum(unlist(balanced_C4$adaptive$time)),
                                      sum(unlist(balanced_C8$adaptive$time)),
                                      sum(unlist(balanced_C16$adaptive$time)),
                                      sum(unlist(balanced_C32$adaptive$time)),
                                      sum(unlist(balanced_C64$adaptive$time))))
consensus_time <- c(consensus_mat_4$time,
                    consensus_mat_8$time,
                    consensus_mat_16$time,
                    consensus_mat_32$time,
                    consensus_mat_64$time)
neiswanger_time <- c(neiswanger_false_4$time,
                     neiswanger_false_8$time,
                     neiswanger_false_16$time,
                     neiswanger_false_32$time,
                     neiswanger_false_64$time)
weierstrass_time <- c(weierstrass_rejection_4$time,
                      weierstrass_rejection_8$time,
                      weierstrass_rejection_16$time,
                      weierstrass_rejection_32$time,
                      weierstrass_rejection_64$time)

plot(x = log(c(4, 8, 16, 32, 64), 2), y = log(balanced_time$adaptive, 2),
     ylim = c(-4, 20),
     xlab = '',
     ylab = '',
     yaxt = 'n',
     xaxt = 'n', lty = 3, lwd = 3, pch = 3, type = 'b')
mtext('log(C, 2)', 1, 2.75, font = 2, cex = 1.5)
mtext('log(Time elapsed in seconds, 2)', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=c(seq(0, 0.9, 0.1), 0.95), labels=c("0.0", c(seq(0.1, 0.9, 0.1), 0.95)), font = 2, cex = 1.5)
axis(1, at=log(c(4, 8, 16, 32, 64), 2), labels = log(c(4, 8, 16, 32, 64), 2), font = 2, cex = 1.5)
axis(2, at=seq(-4, 20, 2), labels = seq(-4, 20, 2), font = 2, cex = 1.5)
axis(2, at=seq(-4, 20, 1), labels=rep("", 25), lwd.ticks = 0.5)
lines(x = log(c(4, 8, 16, 32, 64), 2), y = log(balanced_time$reg, 2),
      lty = 2, lwd = 3, type = 'b', pch = 2)
lines(x = log(c(4, 8, 16, 32), 2), y = log(GBF_time$adaptive, 2),
      lty = 1, lwd = 3, type = 'b', pch = 1)
lines(x = log(c(4, 8, 16, 32, 64), 2), y = log(consensus_time, 2),
      lty = 4, lwd = 3, type = 'b', pch = 4, col = 'red')
lines(x = log(c(4, 8, 16, 32, 64), 2), y = log(neiswanger_time, 2),
      lty = 5, lwd = 3, type = 'b', pch = 5, col = 'red')
lines(x = log(c(4, 8, 16, 32, 64), 2), y = log(weierstrass_time, 2),
      lty = 6, lwd = 3, type = 'b', pch = 6, col = 'red')
lines(x = log(c(4, 8, 16, 32), 2), y = log(NB_fusion_time, 2),
      lty = 7, lwd = 3, type = 'b', pch = 7, col = 'red')
legend(x = 2, y = 20,
       legend = c('GBF (adaptive mesh)',
                  'D&C-GBF (regular mesh)',
                  'D&C-GBF (adaptive mesh)',
                  'D&C-MCF',
                  'CMC',
                  'KDEMC',
                  'WRS'),
       lwd = rep(3, 7),
       lty = c(1,2,3,7,4,5,6),
       pch = c(1,2,3,7,4,5,6),
       col = c(rep('black', 3), rep('red', 4)),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

##### compare to D&C-MCF #####

load('SD4.RData')
load('SD8.RData')
load('SD16.RData')
load('SD32.RData')

P_fusion <- c(integrated_abs_distance(full_posterior,
                                      Poisson_hc_4$particles$y_samples),
              integrated_abs_distance(full_posterior,
                                      Poisson_hc_8$particles$y_samples),
              integrated_abs_distance(full_posterior,
                                      Poisson_hc_16$particles$y_samples))
NB_fusion <- c(integrated_abs_distance(full_posterior,
                                       NB_hc_4$particles$y_samples),
               integrated_abs_distance(full_posterior,
                                       NB_hc_8$particles$y_samples),
               integrated_abs_distance(full_posterior,
                                       NB_hc_16$particles$y_samples),
               integrated_abs_distance(full_posterior,
                                       NB_hc_32$particles$y_samples))
P_fusion_time <- c(sum(unlist(Poisson_hc_4$time)),
                   sum(unlist(Poisson_hc_8$time)),
                   sum(unlist(Poisson_hc_16$time)))
NB_fusion_time <- c(sum(unlist(NB_hc_4$time)),
                    sum(unlist(NB_hc_8$time)),
                    sum(unlist(NB_hc_16$time)),
                    sum(unlist(NB_hc_32$time)))

plot(x = log(c(4, 8, 16, 32, 64), 2), y = balanced$adaptive,
     ylim = c(0, 0.6),
     xlab = '',
     ylab = '',
     xaxt = 'n', lty = 3, lwd = 3, pch = 3, type = 'b')
mtext('log(C, 2)', 1, 2.75, font = 2, cex = 1.5)
mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=c(seq(0, 0.9, 0.1), 0.95), labels=c("0.0", c(seq(0.1, 0.9, 0.1), 0.95)), font = 2, cex = 1.5)
axis(1, at=log(c(4, 8, 16, 32, 64), 2), labels = log(c(4, 8, 16, 32, 64), 2), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=c("0.0", seq(0.1, 0.9, 0.1), "1.0"), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5)
lines(x = log(c(4, 8, 16, 32, 64), 2), y = balanced$reg,
      lty = 2, lwd = 3, type = 'b', pch = 2)
lines(x = log(c(4, 8, 16, 32), 2), y = GBF$adaptive,
      lty = 1, lwd = 3, type = 'b', pch = 1)
lines(x = log(c(4, 8, 16, 32), 2), y = NB_fusion,
      lty = 4, lwd = 3, type = 'b', pch = 4, col = 'red')
legend(x = 2, y = 0.6,
       legend = c('GBF (adaptive)', 'D&C-GBF (reg)', 'D&C-GBF (adaptive)', 'D&C-MCF'),
       lwd = c(3, 3, 3, 3),
       lty = c(1, 2, 3, 4),
       pch = c(1, 2, 3, 4),
       col = c(rep('black', 3), rep('red', 1)),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

plot(x = log(c(4, 8, 16, 32, 64), 2), y = log(balanced_time$adaptive, 2),
     ylim = c(-4, 20),
     xlab = '',
     ylab = '',
     yaxt = 'n',
     xaxt = 'n', lty = 3, lwd = 3, pch = 3, type = 'b')
mtext('log(C, 2)', 1, 2.75, font = 2, cex = 1.5)
mtext('log(Time elapsed in seconds, 2)', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=c(seq(0, 0.9, 0.1), 0.95), labels=c("0.0", c(seq(0.1, 0.9, 0.1), 0.95)), font = 2, cex = 1.5)
axis(1, at=log(c(4, 8, 16, 32, 64), 2), labels = log(c(4, 8, 16, 32, 64), 2), font = 2, cex = 1.5)
axis(2, at=seq(-4, 20, 2), labels = seq(-4, 20, 2), font = 2, cex = 1.5)
axis(2, at=seq(-4, 20, 1), labels=rep("", 25), lwd.ticks = 0.5)
lines(x = log(c(4, 8, 16, 32, 64), 2), y = log(balanced_time$reg, 2),
      lty = 2, lwd = 3, type = 'b', pch = 2)
lines(x = log(c(4, 8, 16, 32), 2), y = log(GBF_time$adaptive, 2),
      lty = 1, lwd = 3, type = 'b', pch = 1)
lines(x = log(c(4, 8, 16, 32), 2), y = log(NB_fusion_time, 2),
      lty = 7, lwd = 3, type = 'b', pch = 4, col = 'red')
legend(x = 2, y = 20,
       legend = c('GBF (adaptive)', 'D&C-GBF (reg)', 'D&C-GBF (adaptive)', 'D&C-MCF'),
       lwd = c(3, 3, 3, 3),
       lty = c(1, 2, 3, 4),
       pch = c(1, 2, 3, 4),
       col = c(rep('black', 3), rep('red', 1)),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

##### ISBA poster #####

plot(x = log(c(4, 8, 16, 32, 64), 2), y = balanced$adaptive,
     ylim = c(0, 0.6),
     xlab = '',
     ylab = '',
     xaxt = 'n', lty = 1, lwd = 3, pch = 1, type = 'b')
mtext('log(C, 2)', 1, 2.75, font = 2, cex = 1.5)
mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=c(seq(0, 0.9, 0.1), 0.95), labels=c("0.0", c(seq(0.1, 0.9, 0.1), 0.95)), font = 2, cex = 1.5)
axis(1, at=log(c(4, 8, 16, 32, 64), 2), labels = log(c(4, 8, 16, 32, 64), 2), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=c("0.0", seq(0.1, 0.9, 0.1), "1.0"), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5)
lines(x = log(c(4, 8, 16, 32), 2), y = GBF$adaptive,
      lty = 3, lwd = 3, type = 'b', pch = 3)
lines(x = log(c(4, 8, 16, 32, 64), 2), y = consensus,
      lty = 4, lwd = 3, type = 'b', pch = 4, col = 'red')
lines(x = log(c(4, 8, 16, 32, 64), 2), y = neiswanger,
      lty = 5, lwd = 3, type = 'b', pch = 5, col = 'red')
lines(x = log(c(4, 8, 16, 32, 64), 2), y = weierstrass,
      lty = 6, lwd = 3, type = 'b', pch = 6, col = 'red')
lines(x = log(c(4, 8, 16, 32), 2), y = NB_fusion,
      lty = 2, lwd = 3, type = 'b', pch = 2, col = 'black')
legend(x = 2, y = 0.6,
       legend = c('D&C-GBF',
                  'D&C-MCF',
                  'GBF',
                  'CMC',
                  'KDEMC',
                  'WRS'),
       lwd = rep(3, 6),
       lty = c(1,2,3,4,5,6),
       pch = c(1,2,3,4,5,6),
       col = c(rep('black', 3), rep('red', 3)),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

plot(x = log(c(4, 8, 16, 32, 64), 2), y = log(balanced_time$adaptive, 2),
     ylim = c(-4, 20),
     xlab = '',
     ylab = '',
     yaxt = 'n',
     xaxt = 'n', lty = 1, lwd = 3, pch = 1, type = 'b')
mtext('log(C, 2)', 1, 2.75, font = 2, cex = 1.5)
mtext('log(Time elapsed in seconds, 2)', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=c(seq(0, 0.9, 0.1), 0.95), labels=c("0.0", c(seq(0.1, 0.9, 0.1), 0.95)), font = 2, cex = 1.5)
axis(1, at=log(c(4, 8, 16, 32, 64), 2), labels = log(c(4, 8, 16, 32, 64), 2), font = 2, cex = 1.5)
axis(2, at=seq(-4, 20, 2), labels = seq(-4, 20, 2), font = 2, cex = 1.5)
axis(2, at=seq(-4, 20, 1), labels=rep("", 25), lwd.ticks = 0.5)
lines(x = log(c(4, 8, 16, 32), 2), y = log(GBF_time$adaptive, 2),
      lty = 3, lwd = 3, type = 'b', pch = 3)
lines(x = log(c(4, 8, 16, 32, 64), 2), y = log(consensus_time, 2),
      lty = 4, lwd = 3, type = 'b', pch = 4, col = 'red')
lines(x = log(c(4, 8, 16, 32, 64), 2), y = log(neiswanger_time, 2),
      lty = 5, lwd = 3, type = 'b', pch = 5, col = 'red')
lines(x = log(c(4, 8, 16, 32, 64), 2), y = log(weierstrass_time, 2),
      lty = 6, lwd = 3, type = 'b', pch = 6, col = 'red')
lines(x = log(c(4, 8, 16, 32), 2), y = log(NB_fusion_time, 2),
      lty = 2, lwd = 3, type = 'b', pch = 2, col = 'black')
legend(x = 2, y = 20,
       legend = c('D&C-GBF',
                  'D&C-MCF',
                  'GBF',
                  'CMC',
                  'KDEMC',
                  'WRS'),
       lwd = rep(3, 6),
       lty = c(1,2,3,4,5,6),
       pch = c(1,2,3,4,5,6),
       col = c(rep('black', 3), rep('red', 3)),
       cex = 1.25,
       text.font = 2,
       bty = 'n')


##### ISBA slides #####

plot(x = log(c(4, 8, 16, 32, 64), 2), y = balanced$adaptive,
     ylim = c(0, 0.6),
     xlab = '',
     ylab = '',
     xaxt = 'n', lty = 1, lwd = 3, pch = 1, type = 'b')
mtext('log(C, 2)', 1, 2.75, font = 2, cex = 1.5)
mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=c(seq(0, 0.9, 0.1), 0.95), labels=c("0.0", c(seq(0.1, 0.9, 0.1), 0.95)), font = 2, cex = 1.5)
axis(1, at=log(c(4, 8, 16, 32, 64), 2), labels = log(c(4, 8, 16, 32, 64), 2), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=c("0.0", seq(0.1, 0.9, 0.1), "1.0"), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5)
lines(x = log(c(4, 8, 16, 32, 64), 2), y = consensus,
      lty = 4, lwd = 3, type = 'b', pch = 4, col = 'red')
lines(x = log(c(4, 8, 16, 32, 64), 2), y = neiswanger,
      lty = 5, lwd = 3, type = 'b', pch = 5, col = 'red')
lines(x = log(c(4, 8, 16, 32, 64), 2), y = weierstrass,
      lty = 6, lwd = 3, type = 'b', pch = 6, col = 'red')
legend(x = 2, y = 0.6,
       legend = c('D&C-GBF',
                  'CMC',
                  'KDEMC',
                  'WRS'),
       lwd = rep(3, 6),
       lty = c(1,4,5,6),
       pch = c(1,4,5,6),
       col = c('black', rep('red', 3)),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

plot(x = log(c(4, 8, 16, 32, 64), 2), y = log(balanced_time$adaptive, 2),
     ylim = c(-4, 16),
     xlab = '',
     ylab = '',
     yaxt = 'n',
     xaxt = 'n', lty = 1, lwd = 3, pch = 1, type = 'b')
mtext('log(C, 2)', 1, 2.75, font = 2, cex = 1.5)
mtext('log(Time elapsed in seconds, 2)', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=c(seq(0, 0.9, 0.1), 0.95), labels=c("0.0", c(seq(0.1, 0.9, 0.1), 0.95)), font = 2, cex = 1.5)
axis(1, at=log(c(4, 8, 16, 32, 64), 2), labels = log(c(4, 8, 16, 32, 64), 2), font = 2, cex = 1.5)
axis(2, at=seq(-4, 20, 2), labels = seq(-4, 20, 2), font = 2, cex = 1.5)
axis(2, at=seq(-4, 20, 1), labels=rep("", 25), lwd.ticks = 0.5)
lines(x = log(c(4, 8, 16, 32, 64), 2), y = log(consensus_time, 2),
      lty = 4, lwd = 3, type = 'b', pch = 4, col = 'red')
lines(x = log(c(4, 8, 16, 32, 64), 2), y = log(neiswanger_time, 2),
      lty = 5, lwd = 3, type = 'b', pch = 5, col = 'red')
lines(x = log(c(4, 8, 16, 32, 64), 2), y = log(weierstrass_time, 2),
      lty = 6, lwd = 3, type = 'b', pch = 6, col = 'red')
legend(x = 2, y = 16,
       legend = c('D&C-GBF',
                  'CMC',
                  'KDEMC',
                  'WRS'),
       lwd = rep(3, 6),
       lty = c(1,4,5,6),
       pch = c(1,4,5,6),
       col = c('black', rep('red', 3)),
       cex = 1.25,
       text.font = 2,
       bty = 'n')
