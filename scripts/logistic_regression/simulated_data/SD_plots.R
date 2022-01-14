library(DCFusion)
load('SD4.RData')
load('SD8.RData')
load('SD16.RData')
save.image('SD.RData')

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

plot(x = c(4, 8, 16), y = NB_fusion,
     ylim = c(0, 0.6),
     xlab = '',
     ylab = '',
     xaxt = 'n', lwd = 3, pch = 1, type = 'l')
mtext('Number of sub-posteriors (C)', 1, 2.75, font = 2, cex = 1.5)
mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=c(seq(0, 0.9, 0.1), 0.95), labels=c("0.0", c(seq(0.1, 0.9, 0.1), 0.95)), font = 2, cex = 1.5)
axis(1, at=c(4, 8, 16), labels = c(4, 8, 16), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=c("0.0", seq(0.1, 0.9, 0.1), "1.0"), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5)
# lines(x = c(4, 8, 16), y = NB_fusion,
#       lty = 1, lwd = 3)
# points(x = c(4, 8, 16), y = consensus,
#        pch = 4, lwd = 3)
lines(x = c(4, 8, 16), y = consensus,
      lty = 4, lwd = 3)
# points(x = c(4, 8, 16), y = neiswanger,
#        pch = 5, lwd = 3)
lines(x = c(4, 8, 16), y = neiswanger,
      lty = 2, lwd = 3)
# points(x = c(4, 8, 16), y = weierstrass,
#        pch = 6, lwd = 3)
lines(x = c(4, 8, 16), y = weierstrass,
      lty = 3, lwd = 3)
legend(x = 4, y = 0.6,
       legend = c('D&C-MCF', 'CMC', 'KDEMC', 'WRS'),
       # col = c('#D41159', '#FFC20A', '#0C7BDC', '#994F00'),
       lwd = c(3, 3, 3),
       # pch = c(1, 4, 5, 6),
       lty = c(1, 4, 2, 3),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

######################### TIME

P_fusion_time <- c(sum(unlist(Poisson_hc_4$time)),
                   sum(unlist(Poisson_hc_8$time)),
                   sum(unlist(Poisson_hc_16$time)))
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

plot(x = c(4, 8, 16), y = log(NB_fusion_time),
     ylim = c(-2, 8),
     xlab = '',
     ylab = '',
     xaxt = 'n', lwd = 3, pch = , type = 'l')
mtext('Number of sub-posteriors (C)', 1, 2.75, font = 2, cex = 1.5)
mtext('log(Time elapsed in seconds)', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=c(seq(0, 0.9, 0.1), 0.95), labels=c("0.0", c(seq(0.1, 0.9, 0.1), 0.95)), font = 2, cex = 1.5)
axis(1, at=c(4, 8, 16), labels = c(4, 8, 16), font = 2, cex = 1.5)
axis(2, at=seq(-2, 12, 2), labels = seq(-2, 12, 2), font = 2, cex = 1.5)
axis(2, at=seq(-2, 12, 1), labels=rep("", 15), lwd.ticks = 0.5)
# lines(x = c(4, 8, 16), y = log(NB_fusion_time),
#       lty = 1, lwd = 3)
# points(x = c(4, 8, 16), y = log(consensus_time),
#        pch = 4, lwd = 3)s
lines(x = c(4, 8, 16), y = log(consensus_time),
      lty = 4, lwd = 3)
# points(x = c(4, 8, 16), y = log(neiswanger_time),
#        pch = 5, lwd = 3)
lines(x = c(4, 8, 16), y = log(neiswanger_time),
      lty = 2, lwd = 3)
# points(x = c(4, 8, 16), y = log(weierstrass_time),
#        pch = 6, lwd = 3)
lines(x = c(4, 8, 16), y = log(weierstrass_time),
      lty = 3, lwd = 3)
legend(x = 4, y = 8,
       legend = c('D&C-MCF', 'CMC', 'KDEMC', 'WRS'),
       # col = c('#D41159', '#FFC20A', '#0C7BDC', '#994F00'),
       lwd = c(3, 3, 3),
       # pch = c(1, 4, 5, 6),
       lty = c(1, 4, 2, 3),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

