library(DCFusion)
load('SD4_GBF.RData')
load('SD8_GBF.RData')
load('SD16_GBF.RData')
save.image('SD.RData')

NB_VBF_adapt_fusion <- c(integrated_abs_distance(full_posterior,
                                                 NB_VBF_4$adaptive_mesh$particles$y_samples),
                         integrated_abs_distance(full_posterior,
                                                 NB_VBF_8$adaptive_mesh$particles$y_samples),
                         integrated_abs_distance(full_posterior,
                                                 NB_VBF_16$adaptive_mesh$particles$y_samples))
NB_GBF_reg_fusion <-  c(integrated_abs_distance(full_posterior,
                                                NB_GBF_4$reg_mesh$particles$y_samples),
                        integrated_abs_distance(full_posterior,
                                                NB_GBF_8$reg_mesh$particles$y_samples),
                        integrated_abs_distance(full_posterior,
                                                NB_GBF_16$reg_mesh$particles$y_samples))
NB_GBF_adapt_fusion <-  c(integrated_abs_distance(full_posterior,
                                                  NB_GBF_4$adaptive_mesh$particles$y_samples),
                          integrated_abs_distance(full_posterior,
                                                  NB_GBF_8$adaptive_mesh$particles$y_samples),
                          integrated_abs_distance(full_posterior,
                                                  NB_GBF_16$adaptive_mesh$particles$y_samples))
NB_fusion <- c(integrated_abs_distance(full_posterior,
                                       NB_hc_4$particles$y_samples),
               integrated_abs_distance(full_posterior,
                                       NB_hc_8$particles$y_samples),
               integrated_abs_distance(full_posterior,
                                       NB_hc_16$particles$y_samples))
consensus <- c(integrated_abs_distance(full_posterior,
                                       consensus_mat_4$samples),
               integrated_abs_distance(full_posterior,
                                       consensus_mat_8$samples),
               integrated_abs_distance(full_posterior,
                                       consensus_mat_16$samples))
neiswanger <- c(integrated_abs_distance(full_posterior,
                                        neiswanger_false_4$samples),
                integrated_abs_distance(full_posterior,
                                        neiswanger_false_8$samples),
                integrated_abs_distance(full_posterior,
                                        neiswanger_false_16$samples))
weierstrass <- c(integrated_abs_distance(full_posterior,
                                         weierstrass_rejection_4$samples),
                 integrated_abs_distance(full_posterior,
                                         weierstrass_rejection_8$samples),
                 integrated_abs_distance(full_posterior,
                                         weierstrass_rejection_16$samples))

plot(x = c(4, 8, 16), y = NB_GBF_adapt_fusion,
     ylim = c(0, 1),
     xlab = '',
     ylab = '',
     xaxt = 'n', lwd = 3, pch = 1, type = 'b')
mtext('Number of sub-posteriors (C)', 1, 2.75, font = 2, cex = 1.5)
mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=c(seq(0, 0.9, 0.1), 0.95), labels=c("0.0", c(seq(0.1, 0.9, 0.1), 0.95)), font = 2, cex = 1.5)
axis(1, at=c(4, 8, 16), labels = c(4, 8, 16), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=c("0.0", seq(0.1, 0.9, 0.1), "1.0"), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5)
lines(x = c(4, 8, 16), y = NB_GBF_reg_fusion,
      lty = 5, lwd = 3, type = 'b', pch = 4)
lines(x = c(4, 8, 16), y = NB_VBF_adapt_fusion,
      lty = 6, lwd = 3, type = 'b', pch = 5)
lines(x = c(4, 8, 16), y = NB_fusion,
      lty = 2, lwd = 3, type = 'b', pch = 9)
lines(x = c(4, 8, 16), y = consensus,
      lty = 4, lwd = 3, type = 'b', pch = 6)
lines(x = c(4, 8, 16), y = neiswanger,
      lty = 2, lwd = 3, type = 'b', pch = 7)
lines(x = c(4, 8, 16), y = weierstrass,
      lty = 3, lwd = 3, type = 'b', pch = 8)
legend(x = 4, y = 1,
       legend = c('GBF (adaptive)', 'GBF (regular)', 'VBF (regular)', 'D&C-MCF', 'CMC', 'KDEMC', 'WRS'),
       lwd = c(3, 3, 3, 3, 3, 3, 3),
       lty = c(1, 5, 6, 2, 4, 2, 3),
       pch = c(1, 4, 5, 9, 6, 7, 8),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

######################### TIME

NB_VBF_adapt_time <- c(NB_VBF_4$adaptive_mesh$time,
                       NB_VBF_8$adaptive_mesh$time,
                       NB_VBF_16$adaptive_mesh$time)
NB_GBF_reg_time <-  c(NB_GBF_4$reg_mesh$time,
                      NB_GBF_8$reg_mesh$time,
                      NB_GBF_16$reg_mesh$time)
NB_GBF_adapt_time <-  c(NB_GBF_4$adaptive_mesh$time,
                        NB_GBF_8$adaptive_mesh$time,
                        NB_GBF_16$adaptive_mesh$time)
NB_fusion_time <- c(sum(unlist(NB_hc_4$time)),
                    sum(unlist(NB_hc_8$time)),
                    sum(unlist(NB_hc_16$time)))
consensus_time <- c(consensus_mat_4$time,
                    consensus_mat_8$time,
                    consensus_mat_16$time)
neiswanger_time <- c(neiswanger_false_4$time,
                     neiswanger_false_8$time,
                     neiswanger_false_16$time)
weierstrass_time <- c(weierstrass_rejection_4$time,
                      weierstrass_rejection_8$time,
                      weierstrass_rejection_16$time)

plot(x = c(4, 8, 16), y = log(NB_GBF_adapt_time),
     ylim = c(-2, 15),
     xlab = '',
     ylab = '',
     xaxt = 'n', lwd = 3, pch = 1, type = 'b', xaxt = 'n', yaxt = 'n')
mtext('Number of sub-posteriors (C)', 1, 2.75, font = 2, cex = 1.5)
mtext('log(Time elapsed in seconds)', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=c(seq(0, 0.9, 0.1), 0.95), labels=c("0.0", c(seq(0.1, 0.9, 0.1), 0.95)), font = 2, cex = 1.5)
axis(1, at=c(4, 8, 16), labels = c(4, 8, 16), font = 2, cex = 1.5)
axis(2, at=seq(-2, 16, 2), labels = seq(-2, 16, 2), font = 2, cex = 1.5)
axis(2, at=seq(-2, 16, 1), labels=rep("", 19), lwd.ticks = 0.5)
lines(x = c(4, 8, 16), y = log(NB_GBF_reg_time),
      lty = 5, lwd = 3, type = 'b', pch = 4)
lines(x = c(4, 8, 16), y = log(NB_VBF_adapt_time),
      lty = 6, lwd = 3, type = 'b', pch = 5)
lines(x = c(4, 8, 16), y = log(NB_fusion_time),
      lty = 2, lwd = 3, type = 'b', pch = 9)
lines(x = c(4, 8, 16), y = log(consensus_time),
      lty = 4, lwd = 3, type = 'b', pch = 6)
lines(x = c(4, 8, 16), y = log(neiswanger_time),
      lty = 2, lwd = 3, type = 'b', pch = 7)
lines(x = c(4, 8, 16), y = log(weierstrass_time),
      lty = 3, lwd = 3, type = 'b', pch = 8)
legend(x = 4, y = 16,
       legend = c('GBF (adaptive)', 'GBF (regular)', 'VBF (regular)', 'D&C-MCF', 'CMC', 'KDEMC', 'WRS'),
       lwd = c(3, 3, 3, 3, 3, 3, 3),
       lty = c(1, 5, 6, 2, 4, 2, 3),
       pch = c(1, 4, 5, 9, 6, 7, 8),
       cex = 1.25,
       text.font = 2,
       bty = 'n')
