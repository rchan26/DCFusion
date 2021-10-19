library(DCFusion)

seed <- 408
set.seed(seed)
denominator <- 2:32
input_samples <- list()
smc_fnj_results <- list()
smc_bal_results <- list()
smc_prog_results <- list()
target_mc <- sample_exp_4(N = 30000,
                          proposal_mean = 0,
                          proposal_sd = 1,
                          dominating_M = 1.35,
                          beta = 1)
time_choice <- 1

for (i in 1:length(denominator)) {
  print(denominator[i])
  set.seed(seed)
  input_samples[[i]] <- base_rejection_sampler_exp_4(beta = 1/denominator[i],
                                                     nsamples = 30000,
                                                     proposal_mean = 0,
                                                     proposal_sd = 1.5,
                                                     dominating_M = 1.75)
  
  curve(exp_4_density(x, beta = 1/denominator[i]), -4, 4,
        main = denominator[i], ylab = 'tempered pdf')
  for (j in 1:length(input_samples[[i]])) {
    lines(density(input_samples[[i]][[j]]), col = 'blue')
  }
  
  # standard fork and join
  # we have no preconditioning here since we want to compare with the standard MCF approach
  print('performing standard MC fusion')
  smc_fnj_fused <- bal_binary_fusion_SMC_exp_4(N_schedule = 30000,
                                               m_schedule = denominator[i],
                                               time_schedule = time_choice,
                                               base_samples = input_samples[[i]],
                                               L = 2,
                                               mean = 0,
                                               start_beta = 1/denominator[i], 
                                               precondition = FALSE,
                                               resampling_method = 'resid',
                                               ESS_threshold = 0.5,
                                               diffusion_estimator = 'NB',
                                               seed = seed)
  
  smc_fnj_results[[i]] <- list('time' = smc_fnj_fused$time[[1]],
                               'ESS' = smc_fnj_fused$ESS[[1]],
                               'CESS' = smc_fnj_fused$CESS[[1]],
                               'IAD' = integrated_abs_distance_exp_4(fusion_post = resample_particle_y_samples(
                                 particle_set = smc_fnj_fused$particles[[1]],
                                 multivariate = FALSE,
                                 resampling_method = 'resid',
                                 seed = seed)$y_samples))
  
  # balanced binary if denominator[i] is 2, 4, 8, 16 or 32
  if (denominator[i]==2) {
    print('performing balanced binary MC fusion')
    smc_bal_fused <- bal_binary_fusion_SMC_exp_4(N_schedule = 30000,
                                                 m_schedule = 2,
                                                 time_schedule = time_choice,
                                                 base_samples = input_samples[[i]],
                                                 L = 2,
                                                 mean = 0,
                                                 start_beta = 1/2,
                                                 precondition = TRUE, 
                                                 resampling_method = 'resid',
                                                 ESS_threshold = 0.5,
                                                 diffusion_estimator = 'NB',
                                                 seed = seed)
  } else if (denominator[i]==4) {
    print('performing balanced binary MC fusion')
    smc_bal_fused <- bal_binary_fusion_SMC_exp_4(N_schedule = rep(30000, 2),
                                                 m_schedule = rep(2, 2),
                                                 time_schedule = rep(time_choice, 2),
                                                 base_samples = input_samples[[i]],
                                                 L = 3,
                                                 mean = 0,
                                                 start_beta = 1/4,
                                                 precondition = TRUE,
                                                 resampling_method = 'resid',
                                                 ESS_threshold = 0.5,
                                                 diffusion_estimator = 'NB',
                                                 seed = seed)
  } else if (denominator[i]==8) {
    print('performing balanced binary MC fusion')
    smc_bal_fused <- bal_binary_fusion_SMC_exp_4(N_schedule = rep(30000, 3),
                                                 m_schedule = rep(2, 3),
                                                 time_schedule = rep(time_choice, 3),
                                                 base_samples = input_samples[[i]],
                                                 L = 4,
                                                 mean = 0,
                                                 start_beta = 1/8,
                                                 precondition = TRUE,
                                                 resampling_method = 'resid',
                                                 ESS_threshold = 0.5,
                                                 diffusion_estimator = 'NB',
                                                 seed = seed)
  } else if (denominator[i]==16) {
    print('performing balanced binary MC fusion')
    smc_bal_fused <- bal_binary_fusion_SMC_exp_4(N_schedule = rep(30000, 4),
                                                 m_schedule = rep(2, 4),
                                                 time_schedule = rep(time_choice, 4),
                                                 base_samples = input_samples[[i]], 
                                                 L = 5,
                                                 mean = 0,
                                                 start_beta = 1/16,
                                                 precondition = TRUE,
                                                 resampling_method = 'resid',
                                                 ESS_threshold = 0.5,
                                                 diffusion_estimator = 'NB',
                                                 seed = seed)
  } else if (denominator[i]==32) {
    print('performing balanced binary MC fusion')
    smc_bal_fused <- bal_binary_fusion_SMC_exp_4(N_schedule = rep(30000, 5),
                                                 m_schedule = rep(2, 5),
                                                 time_schedule = rep(time_choice, 5),
                                                 base_samples = input_samples[[i]], 
                                                 L = 6,
                                                 mean = 0,
                                                 start_beta = 1/32,
                                                 precondition = TRUE,
                                                 resampling_method = 'resid',
                                                 ESS_threshold = 0.5,
                                                 diffusion_estimator = 'NB',
                                                 seed = seed)
  }
  
  if (denominator[i] %in% c(2, 4, 8, 16, 32)) {
    smc_bal_results[[i]] <- list('time' = smc_bal_fused$time[[1]],
                                 'ESS' = smc_bal_fused$ESS[[1]],
                                 'CESS' = smc_bal_fused$CESS[[1]],
                                 'IAD' = integrated_abs_distance_exp_4(fusion_post = resample_particle_y_samples(
                                   particle_set = smc_bal_fused$particles[[1]],
                                   multivariate = FALSE,
                                   resampling_method = 'resid',
                                   seed = seed)$y_samples))
  } else {
    smc_bal_results[[i]] <- NA
  }
  
  # progressive
  print('performing progressive MC fusion')
  smc_prog_fused <- progressive_fusion_SMC_exp_4(N_schedule = rep(30000, denominator[i]-1),
                                                 time_schedule = rep(time_choice, denominator[i]-1),
                                                 base_samples = input_samples[[i]], 
                                                 mean = 0,
                                                 start_beta = 1/denominator[i],
                                                 precondition = TRUE, 
                                                 resampling_method = 'resid',
                                                 ESS_threshold = 0.5,
                                                 diffusion_estimator = 'NB',
                                                 seed = seed)
  
  smc_prog_results[[i]] <- list('time' = smc_prog_fused$time[[1]],
                                'ESS' = smc_prog_fused$ESS[[1]],
                                'CESS' = smc_prog_fused$CESS[[1]],
                                'IAD' = integrated_abs_distance_exp_4(fusion_post = resample_particle_y_samples(
                                  particle_set = smc_prog_fused$particles[[1]],
                                  multivariate = FALSE,
                                  resampling_method = 'resid',
                                  seed = seed)$y_samples))
  
  ##########
  curve(exp_4_density(x, mean = 0), -4, 4, ylim = c(0, 0.5), main = denominator[i])
  lines(density(target_mc))
  lines(density(smc_fnj_fused$particles[[1]]$y_samples), col = 'orange')
  lines(density(smc_bal_fused$particles[[1]]$y_samples), col = 'green')
  lines(density(smc_prog_fused$particles[[1]]$y_samples), col = 'blue')
}

par(mai = c(1.02, 1, 0.82, 0.42))

######################################## running time

plot(x = 2:16, y = sapply(1:15, function(i) smc_fnj_results[[i]]$time), ylim = c(0, 20),
     ylab = 'Running time in seconds', xlab = 'Number of Subposteriors (C)', col = 'orange', pch = 4)
lines(x = 2:16, y = sapply(1:15, function(i) smc_fnj_results[[i]]$time), col = 'orange')
points(x = c(2, 4, 8, 16), y = c(sum(smc_bal_results[[1]]$time), sum(smc_bal_results[[3]]$time),
                                 sum(smc_bal_results[[7]]$time), sum(smc_bal_results[[15]]$time)), col = 'blue', pch = 4)
lines(x = c(2, 4, 8, 16), y = c(sum(smc_bal_results[[1]]$time), sum(smc_bal_results[[3]]$time),
                                sum(smc_bal_results[[7]]$time), sum(smc_bal_results[[15]]$time)), col = 'blue')
points(x = 2:16, y = sapply(1:15, function(i) sum(smc_prog_results[[i]]$time)), col = 'green', pch = 4)
lines(x = 2:16, y = sapply(1:15, function(i) sum(smc_prog_results[[i]]$time)), col = 'green')

#################### log

plot(x = 2:32, y = sapply(1:31, function(i) log(smc_fnj_results[[i]][[1]])), ylim = c(-1, 10),
     ylab = '(logarithm) Running time in seconds', xlab = 'Number of Subposteriors (C)', col = 'red', pch = 4, lwd = 3)
lines(x = 2:32, y = sapply(1:31, function(i) log(smc_fnj_results[[i]][[1]])), col = 'red')
points(x = c(2, 4, 8, 16, 32), y = log(c(sum(smc_bal_results[[1]]$time), sum(smc_bal_results[[3]]$time),
                                         sum(smc_bal_results[[7]]$time), sum(smc_bal_results[[15]]$time),
                                         sum(smc_bal_results[[31]]$time))), col = 'blue', pch = 0, lwd = 3)
lines(x = c(2, 4, 8, 16, 32), y = log(c(sum(smc_bal_results[[1]]$time), sum(smc_bal_results[[3]]$time),
                                        sum(smc_bal_results[[7]]$time), sum(smc_bal_results[[15]]$time),
                                        sum(smc_bal_results[[31]]$time))), col = 'blue', lty = 2)
points(x = 2:32, y = sapply(1:31, function(i) log(sum(smc_prog_results[[i]]$time))), col = 'darkgreen', pch = 1, lwd = 3)
lines(x = 2:32, y = sapply(1:31, function(i) log(sum(smc_prog_results[[i]]$time))), col = 'darkgreen', lty = 3)
legend(x = 2, y = 10, legend = c('standard', 'balanced binary', 'progressive'),
       lty = c(1, 2, 3), pch = c(4, 0, 1), col = c('red', 'blue', 'darkgreen'))

######################################## ESS (overall)

Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
plot(x = 2:16, y = sapply(1:15, function(i) smc_fnj_results[[i]]$ESS['Q']), ylim = c(0, 30000),
     ylab = 'ESS', xlab = 'Number of Subposteriors (C)',
     col = Okabe_Ito[8], pch = 1, lwd = 3)
lines(x = 2:16, y = sapply(1:15, function(i) smc_fnj_results[[i]]$ESS['Q']),
      col = Okabe_Ito[8], lwd = 3)
points(x = c(2, 4, 8, 16), y = c(smc_bal_results[[1]]$ESS['Q'], smc_bal_results[[3]]$ESS['Q'],
                                 smc_bal_results[[7]]$ESS['Q'], smc_bal_results[[15]]$ESS['Q']), 
       col = Okabe_Ito[5], pch = 0, lwd = 3)
lines(x = c(2, 4, 8, 16), y = c(smc_bal_results[[1]]$ESS['Q'], smc_bal_results[[3]]$ESS['Q'],
                                smc_bal_results[[7]]$ESS['Q'], smc_bal_results[[15]]$ESS['Q']),
      col = Okabe_Ito[5], lty = 2, lwd = 3)
points(x = 2:16, y = sapply(1:15, function(i) smc_prog_results[[i]]$ESS['Q']), 
       col = Okabe_Ito[4], pch = 2, lwd = 3)
lines(x = 2:16, y = sapply(1:15, function(i) smc_prog_results[[i]]$ESS['Q']), 
      col = Okabe_Ito[4], lty = 3, lwd = 3)
legend(x = 2, y = 30000, 
       legend = c('fork-and-join', 'balanced', 'progressive'),
       lty = c(1, 2, 3), 
       lwd = c(3, 3, 3),
       pch = c(1, 0, 2), 
       col = Okabe_Ito[c(8, 5, 4)],
       cex = 1.1,
       bty = 'n')

######################################## IAD (overall)

plot(x = 2:32, y = sapply(1:31, function(i) smc_fnj_results[[i]]$IAD), ylim = c(0, 0.6),
     ylab = '', xlab = '', font.lab = 2, pch = 1, lwd = 3, xaxt = 'n', yaxt = 'n')
axis(1, at = seq(2, 32, 2), labels = seq(2, 32, 2), font = 2, cex = 1.5)
axis(1, at=0:32, labels=rep("", 33), lwd.ticks = 0.5)
mtext('Number of sub-posteriors (C)', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 0.6, 0.05), 
     labels = c("0.0", 0.05, "0.10", 0.15, "0.20", 0.25, "0.30", 0.35, "0.40", 0.45, "0.50", 0.55, 0.6),
     font = 2, cex = 1.5)
mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
lines(x = 2:32, y = sapply(1:31, function(i) smc_fnj_results[[i]]$IAD), lwd = 3)
points(x = c(2, 4, 8, 16, 32), y = c(smc_bal_results[[1]]$IAD, smc_bal_results[[3]]$IAD,
                                     smc_bal_results[[7]]$IAD, smc_bal_results[[15]]$IAD,
                                     smc_bal_results[[31]]$IAD), 
       pch = 0, lwd = 3)
lines(x = c(2, 4, 8, 16, 32), y =c(smc_bal_results[[1]]$IAD, smc_bal_results[[3]]$IAD,
                                   smc_bal_results[[7]]$IAD, smc_bal_results[[15]]$IAD,
                                   smc_bal_results[[31]]$IAD),
      lty = 3, lwd = 3)
points(x = 2:32, y = sapply(1:31, function(i) smc_prog_results[[i]]$IAD), 
       pch = 2, lwd = 3)
lines(x = 2:32, y = sapply(1:31, function(i) smc_prog_results[[i]]$IAD), 
      lty = 2, lwd = 3)
legend(x = 2, y = 0.625, 
       legend = c('fork-and-join', 'balanced', 'progressive'),
       lty = c(1, 3, 2), 
       lwd = c(3, 3, 3),
       pch = c(1, 0, 2), 
       cex = 1.25,
       text.font = 2,
       bty = 'n')

save.image('varying_C_experiments.RData')
