library(DCFusion)
library(HMCBLR)

load('NYC64_varying_N_approximate_replicates.RData')
# # preparing lists to store full posterior and sub-posterior samples
# fp <- full_posterior
# rm(full_posterior, sub_posteriors_64, balanced_C64,
#    consensus_mat_64, consensus_sca_64,
#    neiswanger_true_64, neiswanger_false_64,
#    weierstrass_importance_64, weierstrass_rejection_64)
# # preparing lists to store IAD and elapsed run-times
# CMC_results <- list()
# KDEMC_results <- list()
# WRS_results <- list()
# N_samples <- c(500,1000,2000,5000,10000,15000,20000,30000,50000,100000,200000)
# get_approximate_result <- function(result) {
#   print(paste('time:', result$time))
#   print(paste('log(time):', log(result$time, 2)))
#   return(list('IAD' = integrated_abs_distance(full_posterior, result$samples),
#               'IAD2' = integrated_abs_distance(fp, result$samples),
#               'time' = result$time))
# }
# number_of_replicates <- 10
# warmup <- 2500

for (i in 10:length(N_samples)) {
  print(paste('########## i:', i))
  ##### sampling from target and sub-posteriors #####
  CMC_results[[i]] <- list()
  # KDEMC_results[[i]] <- list()
  WRS_results[[i]] <- list()
  for (rep in 1:number_of_replicates) {
    print(paste('***** rep:', rep))
    set.seed(seed*rep*i)
    print('sampling from posterior')
    full_posterior <- hmc_sample_BLR(full_data_count = full_data_count,
                                     C = 1,
                                     prior_means = prior_means,
                                     prior_variances = prior_variances,
                                     iterations = N_samples[i] + warmup,
                                     warmup = warmup,
                                     chains = 1,
                                     seed = seed*rep*i,
                                     output = T)
    print('sampling from sub-posteriors')
    sub_posteriors_64 <- hmc_base_sampler_BLR(nsamples = N_samples[i],
                                              data_split = data_split_64,
                                              C = C,
                                              prior_means = prior_means,
                                              prior_variances = prior_variances,
                                              warmup = warmup,
                                              seed = seed*rep*i,
                                              output = T)
    
    ##### Applying approximate methodologies #####
    print('applying approximate methodologies')
    consensus_mat_64 <- consensus_scott(S = C,
                                        samples_to_combine = sub_posteriors_64,
                                        indep = F)
    # neiswanger_false_64 <- neiswanger(S = C,
    #                                   samples_to_combine = sub_posteriors_64,
    #                                   anneal = FALSE)
    weierstrass_rejection_64 <- weierstrass(Samples = sub_posteriors_64,
                                            method = 'reject')
    
    CMC_results[[i]][[rep]] <- get_approximate_result(consensus_mat_64)
    # KDEMC_results[[i]][[rep]] <- get_approximate_result(neiswanger_false_64)
    WRS_results[[i]][[rep]] <- get_approximate_result(weierstrass_rejection_64)
    
    rm(full_posterior,
       sub_posteriors_64,
       consensus_mat_64,
       # neiswanger_false_64,
       weierstrass_rejection_64)
    save.image('NYC64_varying_N_approximate_replicates_2.RData')
  }
}

load('NYC64_varying_N_approximate_replicates_2.RData')
CMC_holder <- CMC_results
KDEMC_holder <- KDEMC_results
WRS_holder <- WRS_results
load('NYC64_varying_N.RData')
CMC_results <- CMC_holder
KDEMC_results <- KDEMC_holder
WRS_results <- WRS_holder
sapply(1:8, function(i) balanced_C64_results$adaptive[[i]]$IAD2)

plot(x = log(sapply(1:8, function(i) balanced_C64_results$adaptive[[i]]$time), 2),
     y = sapply(1:8, function(i) balanced_C64_results$adaptive[[i]]$IAD2),
     ylim = c(0, 0.7),
     xlim = c(-2, 18),
     xlab = '',
     ylab = '',
     xaxt = 'n', lty = 2, lwd = 3, pch = 4, type = 'b')
mtext('log(Time elapsed in seconds, 2)', 1, 2.75, font = 2, cex = 1.5)
mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=seq(-4, 22, 2), labels = seq(-4, 22, 2), font = 2, cex = 1.5)
axis(1, at=seq(-4, 22, 1), labels=rep("", 27), lwd.ticks = 0.5)
axis(2, at=seq(0, 1, 0.1), labels=c("0.0", seq(0.1, 0.9, 0.1), "1.0"), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5)
lines(x = log(sapply(1:8, function(i) balanced_C64_results$reg[[i]]$time), 2),
      y = sapply(1:8, function(i) balanced_C64_results$reg[[i]]$IAD2),
      lty = 3, lwd = 3, type = 'b', pch = 5)
abline(h=min(sapply(1:11, function(i) mean(sapply(1:length(CMC_results[[i]]), function(rep) CMC_results[[i]][[rep]]$IAD2)))), col = 'magenta', lty = 2, lwd = 3)
# abline(h=0.00, col = 'black', lty = 1, lwd = 3)

lines(x = log(sapply(1:9, function(i) mean(sapply(1:length(KDEMC_results[[i]]), function(rep) KDEMC_results[[i]][[rep]]$time))), 2),
      y = sapply(1:9, function(i) mean(sapply(1:length(KDEMC_results[[i]]), function(rep) KDEMC_results[[i]][[rep]]$IAD2))),
      lty = 5, lwd = 3, type = 'b', pch = 2, col = 'blue')
lines(x = log(sapply(1:9, function(i) mean(sapply(1:length(KDEMC_results[[i]]), function(rep) KDEMC_results[[i]][[rep]]$time))), 2),
      y = sapply(1:9, function(i) min(sapply(1:length(KDEMC_results[[i]]), function(rep) KDEMC_results[[i]][[rep]]$IAD2))),
      lty = 1, lwd = 3, type = 'l', pch = 2, col = 'blue')
lines(x = log(sapply(1:9, function(i) mean(sapply(1:length(KDEMC_results[[i]]), function(rep) KDEMC_results[[i]][[rep]]$time))), 2),
      y = sapply(1:9, function(i) max(sapply(1:length(KDEMC_results[[i]]), function(rep) KDEMC_results[[i]][[rep]]$IAD2))),
      lty = 1, lwd = 3, type = 'l', pch = 2, col = 'blue')

lines(x = log(sapply(1:11, function(i) mean(sapply(1:length(WRS_results[[i]]), function(rep) WRS_results[[i]][[rep]]$time))), 2),
      y = sapply(1:11, function(i) mean(sapply(1:length(WRS_results[[i]]), function(rep) WRS_results[[i]][[rep]]$IAD2))),
      lty = 6, lwd = 3, type = 'b', pch = 1, col = 'green')
lines(x = log(sapply(1:11, function(i) mean(sapply(1:length(WRS_results[[i]]), function(rep) WRS_results[[i]][[rep]]$time))), 2),
      y = sapply(1:11, function(i) min(sapply(1:length(WRS_results[[i]]), function(rep) WRS_results[[i]][[rep]]$IAD2))),
      lty = 1, lwd = 3, type = 'l', pch = 1, col = 'green')
lines(x = log(sapply(1:11, function(i) mean(sapply(1:length(WRS_results[[i]]), function(rep) WRS_results[[i]][[rep]]$time))), 2),
      y = sapply(1:11, function(i) max(sapply(1:length(WRS_results[[i]]), function(rep) WRS_results[[i]][[rep]]$IAD2))),
      lty = 1, lwd = 3, type = 'l', pch = 1, col = 'green')

lines(x = log(sapply(1:11, function(i) mean(sapply(1:length(CMC_results[[i]]), function(rep) CMC_results[[i]][[rep]]$time))), 2),
      y = sapply(1:11, function(i) mean(sapply(1:length(CMC_results[[i]]), function(rep) CMC_results[[i]][[rep]]$IAD2))),
      lty = 4, lwd = 3, type = 'b', pch = 3, col = 'red')
lines(x = log(sapply(1:11, function(i) mean(sapply(1:length(CMC_results[[i]]), function(rep) CMC_results[[i]][[rep]]$time))), 2),
      y = sapply(1:11, function(i) min(sapply(1:length(CMC_results[[i]]), function(rep) CMC_results[[i]][[rep]]$IAD2))),
      lty = 1, lwd = 3, type = 'l', pch = 3, col = 'red')
lines(x = log(sapply(1:11, function(i) mean(sapply(1:length(CMC_results[[i]]), function(rep) CMC_results[[i]][[rep]]$time))), 2),
      y = sapply(1:11, function(i) max(sapply(1:length(CMC_results[[i]]), function(rep) CMC_results[[i]][[rep]]$IAD2))),
      lty = 1, lwd = 3, type = 'l', pch = 3, col = 'red')

legend(x = -2, y = 0.7,
       legend = c('D&C-GBF (regular mesh)',
                  'D&C-GBF (adaptive mesh)',
                  'CMC',
                  'KDEMC',
                  'WRS'),
       lwd = rep(3, 6),
       lty = c(3,2,4,5,6),
       pch = c(5,4,3,2,1),
       col = c(rep('black', 2), c('red', 'blue', 'green')),
       cex = 1.25,
       text.font = 2,
       bty = 'n')


##### log(IAD, 2) #####

plot(x = log(sapply(1:8, function(i) balanced_C64_results$adaptive[[i]]$time), 2),
     y = log(sapply(1:8, function(i) balanced_C64_results$adaptive[[i]]$IAD2), 2),
     ylim = c(-6, 1),
     xlim = c(-2, 18),
     xlab = '',
     ylab = '',
     xaxt = 'n', lty = 2, lwd = 3, pch = 4, type = 'b')
mtext('log(Time elapsed in seconds, 2)', 1, 2.75, font = 2, cex = 1.5)
mtext('log(Integrated Absolute Distance, 2)', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=seq(-4, 22, 2), labels = seq(-4, 22, 2), font = 2, cex = 1.5)
axis(1, at=seq(-4, 22, 1), labels=rep("", 27), lwd.ticks = 0.5)
axis(2, at=seq(-6, 1, 1), labels=seq(-6, 1, 1), font = 2, cex = 1.5)
axis(2, at=seq(-6, 1, 0.5), labels=rep("", 15), lwd.ticks = 0.5)
lines(x = log(sapply(1:8, function(i) balanced_C64_results$reg[[i]]$time), 2),
      y = log(sapply(1:8, function(i) balanced_C64_results$reg[[i]]$IAD2), 2),
      lty = 3, lwd = 3, type = 'b', pch = 5)
abline(h=min(log(sapply(1:11, function(i) mean(sapply(1:length(CMC_results[[i]]), function(rep) CMC_results[[i]][[rep]]$IAD2))), 2)),
       col = 'magenta', lty = 2, lwd = 3)
# abline(h=0.00, col = 'black', lty = 1, lwd = 3)

lines(x = log(sapply(1:9, function(i) mean(sapply(1:length(KDEMC_results[[i]]), function(rep) KDEMC_results[[i]][[rep]]$time))), 2),
      y = log(sapply(1:9, function(i) mean(sapply(1:length(KDEMC_results[[i]]), function(rep) KDEMC_results[[i]][[rep]]$IAD2))), 2),
      lty = 5, lwd = 3, type = 'b', pch = 2, col = 'blue')
lines(x = log(sapply(1:9, function(i) mean(sapply(1:length(KDEMC_results[[i]]), function(rep) KDEMC_results[[i]][[rep]]$time))), 2),
      y = log(sapply(1:9, function(i) min(sapply(1:length(KDEMC_results[[i]]), function(rep) KDEMC_results[[i]][[rep]]$IAD2))), 2),
      lty = 1, lwd = 3, type = 'l', pch = 2, col = 'blue')
lines(x = log(sapply(1:9, function(i) mean(sapply(1:length(KDEMC_results[[i]]), function(rep) KDEMC_results[[i]][[rep]]$time))), 2),
      y = log(sapply(1:9, function(i) max(sapply(1:length(KDEMC_results[[i]]), function(rep) KDEMC_results[[i]][[rep]]$IAD2))), 2),
      lty = 1, lwd = 3, type = 'l', pch = 2, col = 'blue')

lines(x = log(sapply(1:11, function(i) mean(sapply(1:length(WRS_results[[i]]), function(rep) WRS_results[[i]][[rep]]$time))), 2),
      y = log(sapply(1:11, function(i) mean(sapply(1:length(WRS_results[[i]]), function(rep) WRS_results[[i]][[rep]]$IAD2))), 2),
      lty = 6, lwd = 3, type = 'b', pch = 1, col = 'green')
lines(x = log(sapply(1:11, function(i) mean(sapply(1:length(WRS_results[[i]]), function(rep) WRS_results[[i]][[rep]]$time))), 2),
      y = log(sapply(1:11, function(i) min(sapply(1:length(WRS_results[[i]]), function(rep) WRS_results[[i]][[rep]]$IAD2))), 2),
      lty = 1, lwd = 3, type = 'l', pch = 1, col = 'green')
lines(x = log(sapply(1:11, function(i) mean(sapply(1:length(WRS_results[[i]]), function(rep) WRS_results[[i]][[rep]]$time))), 2),
      y = log(sapply(1:11, function(i) max(sapply(1:length(WRS_results[[i]]), function(rep) WRS_results[[i]][[rep]]$IAD2))), 2),
      lty = 1, lwd = 3, type = 'l', pch = 1, col = 'green')

lines(x = log(sapply(1:11, function(i) mean(sapply(1:length(CMC_results[[i]]), function(rep) CMC_results[[i]][[rep]]$time))), 2),
      y = log(sapply(1:11, function(i) mean(sapply(1:length(CMC_results[[i]]), function(rep) CMC_results[[i]][[rep]]$IAD2))), 2),
      lty = 4, lwd = 3, type = 'b', pch = 3, col = 'red')
lines(x = log(sapply(1:11, function(i) mean(sapply(1:length(CMC_results[[i]]), function(rep) CMC_results[[i]][[rep]]$time))), 2),
      y = log(sapply(1:11, function(i) min(sapply(1:length(CMC_results[[i]]), function(rep) CMC_results[[i]][[rep]]$IAD2))), 2),
      lty = 1, lwd = 3, type = 'l', pch = 3, col = 'red')
lines(x = log(sapply(1:11, function(i) mean(sapply(1:length(CMC_results[[i]]), function(rep) CMC_results[[i]][[rep]]$time))), 2),
      y = log(sapply(1:11, function(i) max(sapply(1:length(CMC_results[[i]]), function(rep) CMC_results[[i]][[rep]]$IAD2))), 2),
      lty = 1, lwd = 3, type = 'l', pch = 3, col = 'red')

legend(x = -2, y = 1,
       legend = c('D&C-GBF (regular mesh)',
                  'D&C-GBF (adaptive mesh)',
                  'CMC',
                  'KDEMC',
                  'WRS'),
       lwd = rep(3, 6),
       lty = c(3,2,4,5,6),
       pch = c(5,4,3,2,1),
       col = c(rep('black', 2), c('red', 'blue', 'green')),
       cex = 1.25,
       text.font = 2,
       bty = 'n')
