library(DCFusion)
library(ggplot2)

seed <- 1983
set.seed(seed)

# ---------- script is written to be run on sever since
# ---------- we're replicating experiments many times

# setting parameters
mean <- rep(0, 2)
sd <- rep(1, 2)
correlations <- c(seq(0, 0.9, 0.1), 0.95)
fusion_time <- 1
number_of_replicates <- 100
nsamples <- 10000
opt_bw <- rep((4/(3*nsamples))^(1/5), 2)
bw <- rep(0.1, 2)
smc_fusion_standard <- list()
smc_fusion_precondition <- list()
for (i in 1:length(correlations)) {
  print(paste('i:', i))
  print(paste('correlation:', correlations[i]))
  smc_fusion_standard[[i]] <- list()
  smc_fusion_precondition[[i]] <- list()
  cov_mat <- matrix(c(1, correlations[i], correlations[i], 1), nrow = 2, ncol = 2, byrow = T)
  for (rep in 1:number_of_replicates) {
    print(paste('rep:', rep))
    # sampling from the sub-posteriors
    set.seed(seed*rep*i)
    input_samples <- lapply(1:2, function(sub) mvrnormArma(N = nsamples, mu = mean, Sigma = cov_mat))
    input_particles <- initialise_particle_sets(samples_to_fuse = input_samples, multivariate = TRUE)
    print('### performing standard fusion')
    standard <- parallel_fusion_SMC_biGaussian(particles = input_particles,
                                               N = nsamples,
                                               m = 2,
                                               time = fusion_time,
                                               mean_vec = mean,
                                               sd_vec = sd,
                                               corr = correlations[i],
                                               betas = rep(1, 2),
                                               precondition_matrices = rep(list(diag(1,2)), 2),
                                               ESS_threshold = 0,
                                               seed = seed*rep*i)
    smc_fusion_standard[[i]][[rep]] <- list('ESS' = standard$ESS, 
                                            'CESS' = standard$CESS,
                                            'time' = standard$time, 
                                            'ESS_per_sec' = standard$ESS / standard$time,
                                            'IAD_bw' = integrated_abs_distance_biGaussian(fusion_post = resample_particle_y_samples(
                                              particle_set = standard$particles,
                                              multivariate = TRUE,
                                              resampling_method = 'resid',
                                              seed = seed*rep*i)$y_samples,
                                              marg_means = mean,
                                              marg_sds = rep(1/sqrt(2), 2),
                                              bw = bw),
                                            'IAD_opt_bw' = integrated_abs_distance_biGaussian(fusion_post = resample_particle_y_samples(
                                              particle_set = standard$particles,
                                              multivariate = TRUE,
                                              resampling_method = 'resid',
                                              seed = seed*rep*i)$y_samples,
                                              marg_means = mean,
                                              marg_sds = rep(1/sqrt(2), 2),
                                              bw = opt_bw))
    print('### performing fusion with a preconditioning matrix')
    precondition <- parallel_fusion_SMC_biGaussian(particles = input_particles,
                                                   N = nsamples,
                                                   m = 2,
                                                   time = fusion_time,
                                                   mean_vec = mean,
                                                   sd_vec = sd,
                                                   corr = correlations[i],
                                                   betas = rep(1, 2),
                                                   precondition_matrices = lapply(input_samples, cov),
                                                   ESS_threshold = 0,
                                                   seed = seed*rep*i)
    smc_fusion_precondition[[i]][[rep]] <- list('ESS' = precondition$ESS, 
                                                'CESS' = precondition$CESS,
                                                'time' = precondition$time,
                                                'ESS_per_sec' = precondition$ESS / precondition$time,
                                                'IAD_bw' = integrated_abs_distance_biGaussian(fusion_post = resample_particle_y_samples(
                                                  particle_set = precondition$particles,
                                                  multivariate = TRUE,
                                                  resampling_method = 'resid',
                                                  seed = seed*rep*i)$y_samples,
                                                  marg_means = mean,
                                                  marg_sds = rep(1/sqrt(2), 2),
                                                  bw = bw),
                                                'IAD_opt_bw' = integrated_abs_distance_biGaussian(fusion_post = resample_particle_y_samples(
                                                  particle_set = precondition$particles,
                                                  multivariate = TRUE,
                                                  resampling_method = 'resid',
                                                  seed = seed*rep*i)$y_samples,
                                                  marg_means = mean,
                                                  marg_sds = rep(1/sqrt(2), 2),
                                                  bw = opt_bw))
    print('saving progress')
    save.image('varying_rho_replicates.RData')
  }
}

# ---------- plot ESS for particular simulation

######################################## rho CESS

plot(correlations, sapply(1:length(correlations), function(i) smc_fusion_standard[[i]][[1]]$CESS[1]), ylim = c(0, nsamples),
     xlab = 'correlation', ylab = 'rho CESS', col = 'blue', xaxt='n')
axis(1, at=seq(0, 1, 0.25))
lines(correlations, sapply(1:length(correlations), function(i) smc_fusion_standard[[i]][[1]]$CESS[1]), col = 'blue')
points(correlations, sapply(1:length(correlations), function(i) smc_fusion_precondition[[i]][[1]]$CESS[1]), col = 'red')
lines(correlations, sapply(1:length(correlations), function(i) smc_fusion_precondition[[i]][[1]]$CESS[1]), col = 'red')
legend('topleft', legend = c('MCF', 'GMCF'), col = c('blue', 'red'), lty = c(1,1))

######################################## Q CESS

plot(correlations, sapply(1:length(correlations), function(i) smc_fusion_standard[[i]][[1]]$CESS[2]), ylim = c(0, nsamples),
     xlab = 'correlation', ylab = 'Q CESS', col = 'blue')
lines(correlations, sapply(1:length(correlations), function(i) smc_fusion_standard[[i]][[1]]$CESS[2]), col = 'blue')
points(correlations, sapply(1:length(correlations), function(i) smc_fusion_precondition[[i]][[1]]$CESS[2]), col = 'red')
lines(correlations, sapply(1:length(correlations), function(i) smc_fusion_precondition[[i]][[1]]$CESS[2]), col = 'red')
legend('topleft', legend = c('MCF', 'GMCF'), col = c('blue', 'red'), lty = c(1,1))

######################################## ESS

plot(correlations, sapply(1:length(correlations), function(i) smc_fusion_standard[[i]][[1]]$ESS[2]),
     xlab = 'correlation', ylab = 'ESS', ylim = c(0, 5000),
     col = 'black', lwd = 3)
lines(correlations, sapply(1:length(correlations), function(i) smc_fusion_standard[[i]][[1]]$ESS[2]),
      col = 'black', lwd = 3)
points(correlations, sapply(1:length(correlations), function(i) smc_fusion_precondition[[i]][[1]]$ESS[2]),
       col = 'black', lwd = 3)
lines(correlations, sapply(1:length(correlations), function(i) smc_fusion_precondition[[i]][[1]]$ESS[2]),
      col = 'black', lty = 3, lwd = 3)
legend(x = 0, y = nsamples, legend = c('MCF', 'GMCF'), col = c('black', 'black'),
       lty = c(1,3), lwd = c(3,3), bty = 'n')

######################################## ESS / second

plot(correlations, sapply(1:length(correlations), function(i) smc_fusion_standard[[i]][[1]]$ESS_per_sec[2]),
     xlab = '', ylab = '', xaxt = 'n', ylim = c(0, 1000), col = 'black', lwd = 3)
mtext('Correlation', 1, 2.75, font = 2, cex = 1.5)
mtext('ESS / second', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=correlations, labels=rep("", length(correlations)), lwd.ticks = 0.5)
axis(1, at=c(seq(0, 0.9, 0.1), 0.95), labels=c("0.0", c(seq(0.1, 0.9, 0.1), 0.95)), font = 2, cex = 1.5)
axis(2, at=seq(0, 1000, 200), labels=seq(0, 1000, 200), font = 2, cex = 1.5)
lines(correlations, sapply(1:length(correlations), function(i) smc_fusion_standard[[i]][[1]]$ESS_per_sec[2]),
      col = 'black', lwd = 3)
points(correlations, sapply(1:length(correlations), function(i) smc_fusion_precondition[[i]][[1]]$ESS_per_sec[2]),
       col = 'black', lwd = 3)
lines(correlations, sapply(1:length(correlations), function(i) smc_fusion_precondition[[i]][[1]]$ESS_per_sec[2]),
      col = 'black', lty = 3, lwd = 3)
legend(x = 0, y = 1000, legend = c('MCF', 'GMCF'), col = c('black', 'black'),
       lty = c(1,3), lwd = c(3,3), bty = 'n')

######################################## ESS / second

set.seed(seed)
indices <- sample(100, 50) # take a sample of 10 of the runs
plot(correlations, sapply(1:length(correlations), function(i) {
  mean(sapply(indices, function(j) smc_fusion_standard[[i]][[j]]$ESS_per_sec[2]))
}), xlab = '', ylab = '', xaxt = 'n', ylim = c(0, 1000), col = 'black', lwd = 3, type = 'l')
mtext('Correlation', 1, 2.75, font = 2, cex = 1.5)
mtext('ESS / second', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=correlations, labels=rep("", length(correlations)), lwd.ticks = 0.5)
axis(1, at=c(seq(0, 0.9, 0.1), 0.95), labels=c("0.0", c(seq(0.1, 0.9, 0.1), 0.95)), font = 2, cex = 1.5)
axis(2, at=seq(0, 1000, 200), labels=seq(0, 1000, 200), font = 2, cex = 1.5)
lines(correlations, sapply(1:length(correlations), function(i) {
  mean(sapply(indices, function(j) smc_fusion_standard[[i]][[j]]$ESS_per_sec[2]))
}), col = 'black', lwd = 3)
lines(correlations, sapply(1:length(correlations), function(i) {
  mean(sapply(indices, function(j) smc_fusion_precondition[[i]][[j]]$ESS_per_sec[2]))
}), col = 'black', lty = 3, lwd = 3)
legend(x = 0, y = 1000, legend = c('MCF', 'GMCF'),
       lty = c(1,3), lwd = c(3,3), cex = 1.25, text.font = 2, bty = 'n')

######################################## log(ESS / second)

set.seed(seed)
indices <- sample(100, 10) # take a sample of 10 of the runs
plot(correlations, log(sapply(1:length(correlations), function(i) {
  mean(sapply(indices, function(j) smc_fusion_standard[[i]][[j]]$ESS_per_sec[2]))
})), xlab = '', ylab = '', xaxt = 'n', ylim = c(-5, 10), col = 'black', lwd = 3)
mtext('Correlation', 1, 2.75, font = 2, cex = 1.5)
mtext('log(ESS / second)', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=correlations, labels=rep("", length(correlations)), lwd.ticks = 0.5)
axis(1, at=c(seq(0, 0.9, 0.1), 0.95), labels=c("0.0", c(seq(0.1, 0.9, 0.1), 0.95)), font = 2, cex = 1.5)
axis(2, at=seq(0, 1000, 200), labels=seq(0, 1000, 200), font = 2, cex = 1.5)
lines(correlations, log(sapply(1:length(correlations), function(i) {
  mean(sapply(indices, function(j) smc_fusion_standard[[i]][[j]]$ESS_per_sec[2]))
})), col = 'black', lwd = 3)
points(correlations, log(sapply(1:length(correlations), function(i) {
  mean(sapply(indices, function(j) smc_fusion_precondition[[i]][[j]]$ESS_per_sec[2]))
})), col = 'black', lwd = 3)
lines(correlations, log(sapply(1:length(correlations), function(i) {
  mean(sapply(indices, function(j) smc_fusion_precondition[[i]][[j]]$ESS_per_sec[2]))
})), col = 'black', lty = 3, lwd = 3)
legend(x = 0, y = 1000, legend = c('MCF', 'GMCF'),
       lty = c(1,3), lwd = c(3,3), cex = 1.25, text.font = 2, bty = 'n')

######################################## IAD

set.seed(seed)
indices <- sample(100, 10) # take a sample of 10 of the runs
plot(correlations, sapply(1:length(correlations), function(i) {
  mean(sapply(indices, function(j) smc_fusion_standard[[i]][[j]]$IAD_opt_bw))
}), xlab = '', ylab = '', xaxt = 'n', ylim = c(0, 0.5), col = 'black', lwd = 3)
mtext('Correlation', 1, 2.75, font = 2, cex = 1.5)
mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=correlations, labels=rep("", length(correlations)), lwd.ticks = 0.5)
axis(1, at=c(seq(0, 0.9, 0.1), 0.95), labels=c("0.0", c(seq(0.1, 0.9, 0.1), 0.95)), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=seq(0, 1, 0.1), font = 2, cex = 1.5)
lines(correlations, sapply(1:length(correlations), function(i) {
  mean(sapply(indices, function(j) smc_fusion_standard[[i]][[j]]$IAD_opt_bw))
}), col = 'black', lwd = 3)
points(correlations, sapply(1:length(correlations), function(i) {
  mean(sapply(indices, function(j) smc_fusion_precondition[[i]][[j]]$IAD_opt_bw))
}), col = 'black', lwd = 3)
lines(correlations, sapply(1:length(correlations), function(i) {
  mean(sapply(indices, function(j) smc_fusion_precondition[[i]][[j]]$IAD_opt_bw))
}), col = 'black', lty = 3, lwd = 3)
legend(x = 0, y = 0.5, legend = c('MCF', 'GMCF'),
       lty = c(1,3), lwd = c(3,3), cex = 1.25, text.font = 2, bty = 'n')

######################################## log(IAD)

set.seed(seed)
indices <- sample(100, 10) # take a sample of 10 of the runs
plot(correlations, log(sapply(1:length(correlations), function(i) {
  mean(sapply(indices, function(j) smc_fusion_standard[[i]][[j]]$IAD_opt_bw))
})), xlab = '', ylab = '', xaxt = 'n', ylim = c(-5, 0), col = 'black', lwd = 3)
mtext('Correlation', 1, 2.75, font = 2, cex = 1.5)
mtext('log(Integrated Absolute Distance)', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=correlations, labels=rep("", length(correlations)), lwd.ticks = 0.5)
axis(1, at=c(seq(0, 0.9, 0.1), 0.95), labels=c("0.0", c(seq(0.1, 0.9, 0.1), 0.95)), font = 2, cex = 1.5)
axis(2, at=seq(0, 1000, 200), labels=seq(0, 1000, 200), font = 2, cex = 1.5)
lines(correlations, log(sapply(1:length(correlations), function(i) {
  mean(sapply(indices, function(j) smc_fusion_standard[[i]][[j]]$IAD_opt_bw))
})), col = 'black', lwd = 3)
points(correlations, log(sapply(1:length(correlations), function(i) {
  mean(sapply(indices, function(j) smc_fusion_precondition[[i]][[j]]$IAD_opt_bw))
})), col = 'black', lwd = 3)
lines(correlations, log(sapply(1:length(correlations), function(i) {
  mean(sapply(indices, function(j) smc_fusion_precondition[[i]][[j]]$IAD_opt_bw))
})), col = 'black', lty = 3, lwd = 3)
legend(x = 0, y = 1000, legend = c('MCF', 'GMCF'),
       lty = c(1,3), lwd = c(3,3), cex = 1.25, text.font = 2, bty = 'n')

########################################
# ---------- boxplot ESS for replicates

Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
theme_set(theme_bw())
theme_update(text = element_text(size=12),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank()
)

par(mai = c(1.02, 1, 0.82, 0.42))
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

# create dataframe where each column is the ESS of the replicate simulation for each correlation
ESS <- data.frame()
for (i in 1:length(correlations)) {
  ESS <- rbind(ESS,
               data.frame('ESS' = c(sapply(1:rep, function(rep) smc_fusion_standard[[i]][[rep]]$ESS[2]),
                                    sapply(1:rep, function(rep) smc_fusion_precondition[[i]][[rep]]$ESS[2])),
                          'Correlation' = correlations[i],
                          'Standard' = as.factor(c(rep('MCF', number_of_replicates), rep('GMCF',  number_of_replicates)))))
}
rownames(ESS) <- c()
ESS$Correlation <- as.factor(ESS$Correlation)
ESS$Standard <- as.factor(ESS$Standard)

# plot boxplots for ESS comparison
ggplot(data = ESS, aes(x = Correlation, y = ESS, fill = Standard)) +
  geom_boxplot(outlier.colour = 'white', position = position_dodge(0.9)) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = Okabe_Ito[c(4,5)])

# create dataframe where each column is the ESS/sec of the replication simulation for each correlation
ESS_per_second <- data.frame()
for (i in 1:length(correlations)) {
  ESS_per_second <- rbind(ESS_per_second,
                          data.frame('ESS_per_sec' = c(sapply(1:rep, function(rep) smc_fusion_standard[[i]][[rep]]$ESS_per_sec[2]),
                                                       sapply(1:rep, function(rep) smc_fusion_precondition[[i]][[rep]]$ESS_per_sec[2])),
                                     'Correlation' = correlations[i],
                                     'Standard' = as.factor(c(rep('MCF', number_of_replicates), rep('GMCF',  number_of_replicates)))))
}
rownames(ESS_per_second) <- c()
ESS_per_second$Correlation <- as.factor(ESS_per_second$Correlation)
ESS_per_second$Standard <- as.factor(ESS_per_second$Standard)

# plot boxplots for ESS/sec for ESS comparison
ggplot(data = ESS_per_second, aes(x = Correlation, y = ESS_per_sec, fill = Standard)) +
  geom_boxplot(outlier.colour = 'white', position = position_dodge(0.9)) +
  ylab('ESS / sec') +
  theme(legend.title = element_blank(),
        legend.text = element_text(face = 'bold', size = 12),
        axis.text = element_text(face = 'bold', size = 12),
        axis.title = element_text(face = 'bold', size = 14)) +
  scale_fill_manual(values = Okabe_Ito[c(4,5)])

# create dataframe where each column is the ESS/sec of the replication simulation for each correlation
ESS_over_N <- data.frame()
for (i in 1:length(correlations)) {
  ESS_over_N <- rbind(ESS_over_N,
                      data.frame('ESS_over_N' = c(sapply(1:rep, function(rep) smc_fusion_standard[[i]][[rep]]$ESS[2]/nsamples),
                                                  sapply(1:rep, function(rep) smc_fusion_precondition[[i]][[rep]]$ESS[2]/nsamples)),
                                 'Correlation' = correlations[i],
                                 'Standard' = as.factor(c(rep('MCF', number_of_replicates), rep('GMCF',  number_of_replicates)))))
}
rownames(ESS_over_N) <- c()
ESS_over_N$Correlation <- as.factor(ESS_over_N$Correlation)
ESS_over_N$Standard <- as.factor(ESS_over_N$Standard)

# plot boxplots for ESS/sec for ESS comparison
ggplot(data = ESS_over_N, aes(x = Correlation, y = ESS_over_N, fill = Standard)) +
  geom_boxplot(outlier.colour = 'white', position = position_dodge(0.9)) +
  ylab('ESS / N') +
  theme(legend.title = element_blank(),
        legend.text = element_text(face = 'bold', size = 12),
        axis.text = element_text(face = 'bold', size = 12),
        axis.title = element_text(face = 'bold', size = 14)) +
  scale_fill_manual(values = Okabe_Ito[c(4,5)])
