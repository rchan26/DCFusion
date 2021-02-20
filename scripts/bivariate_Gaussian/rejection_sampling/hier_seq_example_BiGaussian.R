library(hierarchicalFusion)
library(MASS)

#####

seed <- 1983
set.seed(seed)

# setting parameters
mean <- c(1, 2)
sd <- c(sqrt(0.5), sqrt(2))
cov_mat <- matrix(c(sd[1]^2, 0.9, 0.9, sd[2]^2), nrow = 2, ncol = 2, byrow = T)
corr <- 0.9/(sd[1]*sd[2])
fusion_time <- 1

# sampling from the sub-posteriors (target at inverse temperature 1/8)
input_samples <- lapply(1:8, function(i) mvrnormArma_tempered(10000,
                                                              mu = mean,
                                                              Sigma = cov_mat,
                                                              beta = 1/8))

# sample from true target density
true_samples <- mvrnormArma(100000, mu = mean, Sigma = cov_mat)
true_kde <- MASS::kde2d(true_samples[,1], true_samples[,2], n = 50)
image(true_kde)
contour(true_kde, add = T)

##############################

test_hier_standard <- hierarchical_fusion_biGaussian(N_schedule = rep(10000, 3),
                                                     m_schedule = rep(2, 3),
                                                     time_schedule = rep(fusion_time, 3),
                                                     base_samples = input_samples,
                                                     L = 4,
                                                     mean_vec = mean,
                                                     sd_vec = sd,
                                                     corr = corr,
                                                     start_beta = 1/8,
                                                     precondition = FALSE,
                                                     seed = seed, 
                                                     n_cores = 12)
kde_hier_standard <- MASS::kde2d(test_hier_standard$samples[[1]][,1],
                                 test_hier_standard$samples[[1]][,2], 
                                 n = 50)

test_hier_precondition <- hierarchical_fusion_biGaussian(N_schedule = rep(10000, 3),
                                                         m_schedule = rep(2, 3),
                                                         time_schedule = rep(fusion_time, 3),
                                                         base_samples = input_samples,
                                                         L = 4,
                                                         mean_vec = mean,
                                                         sd_vec = sd,
                                                         corr = corr,
                                                         start_beta = 1/8,
                                                         precondition = TRUE,
                                                         seed = seed, 
                                                         n_cores = 12)
kde_hier_precondition <- MASS::kde2d(test_hier_precondition$samples[[1]][,1],
                                     test_hier_precondition$samples[[1]][,2], 
                                     n = 50)

##############################

test_prog_standard <- progressive_fusion_biGaussian(N_schedule = rep(10000, 7),
                                                    time_schedule = rep(fusion_time, 7),
                                                    mean_vec = mean,
                                                    sd_vec = sd,
                                                    corr = corr,
                                                    start_beta = 1/8,
                                                    precondition = FALSE,
                                                    base_samples = input_samples,
                                                    seed = seed)
kde_prog_standard <- MASS::kde2d(test_prog_standard$samples[[1]][,1],
                                 test_prog_standard$samples[[1]][,2],
                                 n = 50)

test_prog_precondition <- progressive_fusion_biGaussian(N_schedule = rep(10000, 7),
                                                        time_schedule = rep(fusion_time, 7),
                                                        mean_vec = mean,
                                                        sd_vec = sd,
                                                        corr = corr,
                                                        start_beta = 1/8,
                                                        precondition = TRUE,
                                                        base_samples = input_samples,
                                                        seed = seed)
kde_prog_precondition <- MASS::kde2d(test_prog_precondition$samples[[1]][,1],
                                     test_prog_precondition$samples[[1]][,2],
                                     n = 50)

##############################

xlims <- c(-2, 4)
ylims <- c(-4, 8)
image(true_kde, xlim = xlims, ylim = ylims)
contour(true_kde, add = T)
#####
contour(kde_hier_standard, add = T)
contour(kde_hier_precondition, add = T)
contour(kde_prog_standard, add = T)
contour(kde_prog_precondition, add = T)
#####
image(true_kde, xlim = xlims, ylim = ylims)
image(kde_hier_standard, xlim = xlims, ylim = ylims)
image(kde_prog_standard, xlim = xlims, ylim = ylims)
image(kde_hier_precondition, xlim = xlims, ylim = ylims)
image(kde_prog_precondition, xlim = xlims, ylim = ylims)
#####
acceptance_rate_plots(hier1 = test_hier_precondition,
                      hier2 = test_hier_standard,
                      time = fusion_time,
                      hierarchical = TRUE)
acceptance_rate_plots(hier1 = test_prog_precondition,
                      hier2 = test_prog_standard,
                      time = fusion_time,
                      hierarchical = FALSE)

#####

par(oma = c(2, 1, 0, 1))
par(mfrow = c(1,2))

##### hier
# rho
plot(1:length(test_hier_standard$overall_rho), test_hier_standard$overall_rho, ylim = c(0,1), col = 'black',
     ylab = expression(rho), xlab = 'Level', pch = 4, lwd = 1.5, xaxt = 'n')
axis(1, at = c(1,2,3))
lines(1:length(test_hier_standard$overall_rho), test_hier_standard$overall_rho, col = 'black', lwd = 1.5)
points(1:length(test_hier_precondition$overall_rho), test_hier_precondition$overall_rho, col = 'magenta', pch = 2, lwd = 1.5)
lines(1:length(test_hier_precondition$overall_rho), test_hier_precondition$overall_rho, col = 'magenta', lty = 2, lwd = 1.5)
# Q
plot(1:length(test_hier_standard$overall_Q), test_hier_standard$overall_Q, ylim = c(0,1), col = 'black',
     ylab = 'Q', xlab = 'Level', pch = 4, lwd = 1.5, xaxt = 'n')
axis(1, at = c(1,2,3))
lines(1:length(test_hier_standard$overall_Q), test_hier_standard$overall_Q, col = 'black', lwd = 1.5)
points(1:length(test_hier_precondition$overall_Q), test_hier_precondition$overall_Q, col = 'magenta', pch = 2, lwd = 1.5)
lines(1:length(test_hier_precondition$overall_Q), test_hier_precondition$overall_Q, col = 'magenta', lty = 2, lwd = 1.5)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", 
       legend = c('standard', 'time-adapting'), 
       lty = c(1, 2),
       xpd = TRUE, 
       horiz = TRUE, 
       inset = c(0, 0), 
       bty = "n", 
       pch = c(4, 2), 
       col = c('black', 'magenta'))

##### prog
# rho
plot(1:length(test_prog_standard$rho_acc), test_prog_standard$rho_acc, ylim = c(0,1), col = 'black',
     ylab = expression(rho), xlab = 'Level', pch = 4, lwd = 1.5)
lines(1:length(test_prog_standard$rho_acc), test_prog_standard$rho_acc, col = 'black', lwd = 1.5)
points(1:length(test_prog_precondition$rho_acc), test_prog_precondition$rho_acc, col = 'magenta', pch = 2, lwd = 1.5)
lines(1:length(test_prog_precondition$rho_acc), test_prog_precondition$rho_acc, col = 'magenta', lty = 2, lwd = 1.5)
# Q
plot(1:length(test_prog_standard$Q_acc), test_prog_standard$Q_acc, ylim = c(0,1), col = 'black',
     ylab = 'Q', xlab = 'Level', pch = 4, lwd = 1.5)
lines(1:length(test_prog_standard$Q_acc), test_prog_standard$Q_acc, col = 'black')
points(1:length(test_prog_precondition$Q_acc), test_prog_precondition$Q_acc, col = 'magenta', pch = 2, lwd = 1.5)
lines(1:length(test_prog_precondition$Q_acc), test_prog_precondition$Q_acc, col = 'magenta', lty = 2, lwd = 1.5)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", 
       legend = c('standard', 'time-adapting'), 
       lty = c(1, 2),
       xpd = TRUE, 
       horiz = TRUE, 
       inset = c(0, 0), 
       bty = "n", 
       pch = c(4, 2), 
       col = c('black', 'magenta'))

#####

sum(test_hier_standard$overall_time)
sum(test_hier_precondition$overall_time)
sum(test_prog_standard$time)
sum(test_prog_precondition$time)

#####
par(original_par)
layout(matrix(1:4, nrow = 2, byrow = T), widths = c(10,5), heights = c(4,10))

# x1 density plots
par(mar = c(0.25, 5, 1, 0))
plot(density(true_samples[,1]), ylab = 'density', main = NA, axes = F, ylim = c(0, 0.3),
     xlab = NA, col = 'black', lwd = 2)
axis(2, las = 1)
lines(density(test_hier_standard$samples[[1]][,1]), col = 'blue', lty = 2, lwd = 2)
lines(density(test_hier_precondition$samples[[1]][,1]), col = 'violet', lty = 4, lwd = 2)
lines(density(test_prog_standard$samples[[1]][,1]), col = 'green', lty = 3, lwd = 2)
lines(density(test_prog_precondition$samples[[1]][,1]), col = 'deeppink', lty = 5, lwd = 2)

# legend
plot.new()
par(mar = c(0.25,0,0,0))
legend("left", legend = c('true target',
                          'hierarchical (standard)', 
                          'hierarchical (time-adapting)',
                          'progressive (standard)',
                          'progressive (time-adapting)'),
       col = c('black', 'blue', 'violet', 'green', 'deeppink'), 
       lty = c(1, 2, 4, 3, 5), 
       lwd = c(2, 2, 2, 2, 2))

par(mar = c(4, 5, 0, 0))
# bivariate plots
contour(true_kde, lwd = 2, xlab = expression(x[1]), ylab = expression(x[2]))
contour(kde_hier_standard, add = T, col = 'blue', lty = 2, lwd = 2)
contour(kde_hier_precondition, add = T, col = 'violet', lty = 4, lwd = 2)
contour(kde_prog_standard, add = T, col = 'green', lty = 3, lwd = 2)
contour(kde_prog_precondition, add = T, col = 'deeppink', lty = 5, lwd = 2)
axis(3, labels = F, tck = 0.01)
axis(4, labels = F, tck = 0.01)
box()

# x2 density plots
par(mar = c(4, 0.25, 0, 1))
plot(x = density(true_samples[,2])$y, y = density(true_samples[,2])$x, xlab = 'density', 
     main = NA, axes = F, xlim = c(0, 0.25), type = 'l', lwd = 2)
axis(1)
lines(x = density(test_hier_standard$samples[[1]][,2])$y,
      y = density(test_hier_standard$samples[[1]][,2])$x,
      col = 'blue', lty = 2, lwd = 2)
lines(x = density(test_hier_precondition$samples[[1]][,2])$y,
      y = density(test_hier_precondition$samples[[1]][,2])$x,
      col = 'violet', lty = 4, lwd = 2)
lines(x = density(test_prog_standard$samples[[1]][,2])$y, 
      y = density(test_prog_precondition$samples[[1]][,2])$x, 
      col = 'green', lty = 3, lwd = 2)
lines(x = density(test_prog_precondition$samples[[1]][,2])$y, 
      y = density(test_prog_precondition$samples[[1]][,2])$x,
      col = 'deeppink', lty = 5, lwd = 2)
