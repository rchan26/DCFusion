fusion <- c(integrated_abs_distance(full_posterior,
                                    Poisson_hc_32[[2]]$particles$y_samples),
            integrated_abs_distance(full_posterior,
                                    Poisson_hc_64[[1]]$particles$y_samples),
            integrated_abs_distance(full_posterior,
                                    Poisson_hc_128[[7]]$particles$y_samples))
consensus <- c(integrated_abs_distance(full_posterior,
                                       consensus_mat_32$samples),
               integrated_abs_distance(full_posterior,
                                       consensus_mat_64$samples),
               integrated_abs_distance(full_posterior,
                                       consensus_mat_128$samples))
neiswanger <- c(integrated_abs_distance(full_posterior,
                                        neiswanger_32_false$samples),
                integrated_abs_distance(full_posterior,
                                        neiswanger_64_false$samples),
                integrated_abs_distance(full_posterior,
                                        neiswanger_128_false$samples))
weierstrass <- c(integrated_abs_distance(full_posterior,
                                         weierstrass_32_rejection$samples),
                 integrated_abs_distance(full_posterior,
                                         weierstrass_64_rejection$samples),
                 integrated_abs_distance(full_posterior,
                                         weierstrass_128_rejection$samples))

plot(x = c(32, 64, 128), y = fusion,
     ylim = c(0, 1),
     xlab = 'Number of sub-posteriors (C)',
     ylab = 'Integrated Absolute Distance',
     col = '#D41159', xaxt = 'n', lwd = 3, pch = 1)
lines(x = c(32, 64, 128), y = fusion,
      lty = 5, col = '#D41159', lwd = 3)
axis(1, at = c(32, 64, 128), labels = c(32, 64, 128))
points(x = c(32, 64, 128), y = consensus,
       pch = 4, col = '#FFC20A', lwd = 3)
lines(x = c(32, 64, 128), y = consensus,
       lty = 4, col = '#FFC20A', lwd = 3)
points(x = c(32, 64, 128), y = neiswanger,
       pch = 5, col = '#0C7BDC', lwd = 3)
lines(x = c(32, 64, 128), y = neiswanger,
       lty = 3, col = '#0C7BDC', lwd = 3)
points(x = c(32, 64, 128), y = weierstrass,
       pch = 6, col = '#994F00', lwd = 3)
lines(x = c(32, 64, 128), y = weierstrass,
      lty = 2, col = '#994F00', lwd = 3)
legend(x = 105, y = 1,
       legend = c('Fusion', 'Consensus', 'Neiswanger', 'Weierstrass'),
       col = c('#D41159', '#FFC20A', '#0C7BDC', '#994F00'),
       lwd = c(3, 3, 3),
       pch = c(1, 4, 5, 6),
       lty = c(5, 4, 3, 2))
