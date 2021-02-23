library(hierarchicalFusion)

seed <- 408
set.seed(seed)
denominator <- 2:16
input_samples <- list()
smc_fnj_results <- list()
smc_hier_results <- list()
smc_prog_results <- list()
time_choice <- 0.5

for (i in 1:length(denominator)) {
  print(denominator[i])
  input_samples[[i]] <- base_rejection_sampler_exp_4(beta = 1/denominator[i],
                                                     nsamples = 10000,
                                                     proposal_mean = 0,
                                                     proposal_sd = 1.5,
                                                     dominating_M = 1.75)
  
  curve(exp_4_density(x, beta = 1/denominator[i]), -4, 4,
        main = denominator[i], ylab = 'tempered pdf')
  for (j in 1:length(input_samples[[i]])) {
    lines(density(input_samples[[i]][[j]]), col = 'blue')
  }
  
  # standard fork and join (with no resampling)
  print('performing standard MC fusion')
  smc_fnj_fused <- hierarchical_fusion_SMC_exp_4(N_schedule = 10000,
                                                 m_schedule = denominator[i],
                                                 time_schedule = time_choice,
                                                 base_samples = input_samples[[i]],
                                                 L = 2,
                                                 mean = 0,
                                                 start_beta = 1/denominator[i], 
                                                 C = denominator[i], 
                                                 precondition = TRUE, 
                                                 ESS_threshold = 0,
                                                 seed = seed)
  
  smc_fnj_results[[i]] <- list('time' = smc_fnj_fused$time[[1]],
                               'ESS' = smc_fnj_fused$ESS[[1]],
                               'CESS' = smc_fnj_fused$CESS[[1]])
  
  # hierarchical if denominator[i] is 2, 4, 8, or 16 (with no resampling)
  if (denominator[i]==2) {
    print('performing hierarchical MC fusion')
    smc_hier_fused <- hierarchical_fusion_SMC_exp_4(N_schedule = 10000,
                                                    m_schedule = 2,
                                                    time_schedule = time_choice,
                                                    base_samples = input_samples[[i]],
                                                    L = 2,
                                                    mean = 0,
                                                    start_beta = 1/2,
                                                    precondition = TRUE, 
                                                    ESS_threshold = 0, 
                                                    seed = seed)
  } else if (denominator[i]==4) {
    print('performing hierarchical MC fusion')
    smc_hier_fused <- hierarchical_fusion_SMC_exp_4(N_schedule = rep(10000, 2),
                                                    m_schedule = rep(2, 2),
                                                    time_schedule = rep(time_choice, 2),
                                                    base_samples = input_samples[[i]],
                                                    L = 3,
                                                    mean = 0,
                                                    start_beta = 1/4,
                                                    precondition = TRUE,
                                                    ESS_threshold = 0, 
                                                    seed = seed)
  } else if (denominator[i]==8) {
    print('performing hierarchical MC fusion')
    smc_hier_fused <- hierarchical_fusion_SMC_exp_4(N_schedule = rep(10000, 3),
                                                    m_schedule = rep(2, 3),
                                                    time_schedule = rep(time_choice, 3),
                                                    base_samples = input_samples[[i]],
                                                    L = 4,
                                                    mean = 0,
                                                    start_beta = 1/8,
                                                    precondition = TRUE,
                                                    ESS_threshold = 0, 
                                                    seed = seed)
  } else if (denominator[i]==16) {
    print('performing hierarchical MC fusion')
    smc_hier_fused <- hierarchical_fusion_SMC_exp_4(N_schedule = rep(10000, 4),
                                                    m_schedule = rep(2, 4),
                                                    time_schedule = rep(time_choice, 4),
                                                    base_samples = input_samples[[i]], 
                                                    L = 5,
                                                    mean = 0,
                                                    start_beta = 1/16,
                                                    precondition = TRUE,
                                                    ESS_threshold = 0, 
                                                    seed = seed)
  }
  
  if (denominator[i] %in% c(2, 4, 8, 16)) {
    smc_hier_results[[i]] <- list('time' = smc_hier_fused$time[[1]],
                                  'ESS' = smc_hier_fused$ESS[[1]],
                                  'CESS' = smc_hier_fused$CESS[[1]])
  } else {
    smc_hier_results[[i]] <- NA
  }
  
  # progressive (with no resampling)
  print('performing progressive MC fusion')
  smc_prog_fused <- progressive_fusion_SMC_exp_4(N_schedule = rep(10000, denominator[i]-1),
                                                 time_schedule = rep(time_choice, denominator[i]-1),
                                                 base_samples = input_samples[[i]], 
                                                 mean = 0,
                                                 start_beta = 1/denominator[i],
                                                 precondition = TRUE, 
                                                 ESS_threshold = 0,
                                                 seed = seed)
  
  smc_prog_results[[i]] <- list('time' = smc_prog_fused$time[[1]],
                                'ESS' = smc_prog_fused$ESS[[1]],
                                'CESS' = smc_prog_fused$CESS[[1]])
  
  ##########
  curve(exp_4_density(x, mean = 0), -4, 4, ylim = c(0, 0.5), main = denominator[i])
  lines(density(resample_particle_y_samples(smc_fnj_fused$particles[[1]], 
                                            multivariate = FALSE, 
                                            seed = seed)$y_samples), 
        col = 'orange')
  lines(density(resample_particle_y_samples(smc_hier_fused$particles[[1]], 
                                            multivariate = FALSE, 
                                            seed = seed)$y_samples), 
        col = 'green')
  lines(density(resample_particle_y_samples(smc_prog_fused$particles[[1]], 
                                            multivariate = FALSE, 
                                            seed = seed)$y_samples), 
        col = 'blue')
}

par(mai = c(1.02, 1, 0.82, 0.42))

######################################## running time

plot(x = 2:16, y = sapply(1:15, function(i) smc_fnj_results[[i]]$time), ylim = c(0, 20),
     ylab = 'Running time in seconds', xlab = 'Number of Subposteriors (C)', col = 'orange', pch = 4)
lines(x = 2:16, y = sapply(1:15, function(i) smc_fnj_results[[i]]$time), col = 'orange')
points(x = c(2, 4, 8, 16), y = c(sum(smc_hier_results[[1]]$time), sum(smc_hier_results[[3]]$time),
                                 sum(smc_hier_results[[7]]$time), sum(smc_hier_results[[15]]$time)), col = 'blue', pch = 4)
lines(x = c(2, 4, 8, 16), y = c(sum(smc_hier_results[[1]]$time), sum(smc_hier_results[[3]]$time),
                                sum(smc_hier_results[[7]]$time), sum(smc_hier_results[[15]]$time)), col = 'blue')
points(x = 2:16, y = sapply(1:15, function(i) sum(smc_prog_results[[i]]$time)), col = 'green', pch = 4)
lines(x = 2:16, y = sapply(1:15, function(i) sum(smc_prog_results[[i]]$time)), col = 'green')

#################### log

plot(x = 2:16, y = sapply(1:15, function(i) log(smc_fnj_results[[i]][[1]])), ylim = c(-1, 10),
     ylab = '(logarithm) Running time in seconds', xlab = 'Number of Subposteriors (C)', col = 'red', pch = 4, lwd = 3)
lines(x = 2:16, y = sapply(1:15, function(i) log(smc_fnj_results[[i]][[1]])), col = 'red')
points(x = c(2, 4, 8, 16), y = log(c(sum(smc_hier_results[[1]]$time), sum(smc_hier_results[[3]]$time),
                                     sum(smc_hier_results[[7]]$time), sum(smc_hier_results[[15]]$time))), col = 'blue', pch = 0, lwd = 3)
lines(x = c(2, 4, 8, 16), y = log(c(sum(smc_hier_results[[1]]$time), sum(smc_hier_results[[3]]$time),
                                    sum(smc_hier_results[[7]]$time), sum(smc_hier_results[[15]]$time))), col = 'blue', lty = 2)
points(x = 2:16, y = sapply(1:15, function(i) log(sum(smc_prog_results[[i]]$time))), col = 'darkgreen', pch = 1, lwd = 3)
lines(x = 2:16, y = sapply(1:15, function(i) log(sum(smc_prog_results[[i]]$time))), col = 'darkgreen', lty = 3)
legend(x = 2, y = 10, legend = c('standard', 'hierarchical', 'progressive'),
       lty = c(1, 2, 3), pch = c(4, 0, 1), col = c('red', 'blue', 'darkgreen'))

######################################## rho.Q acceptance (overall)

plot(x = 2:16, y = sapply(1:15, function(i) smc_fnj_results[[i]]$ESS['Q']), ylim = c(0, 10000),
     ylab = 'ESS', xlab = 'Number of Subposteriors (C)',
     col = Okabe_Ito[8], pch = 1, lwd = 3)
lines(x = 2:16, y = sapply(1:15, function(i) smc_fnj_results[[i]]$ESS['Q']),
      col = Okabe_Ito[8], lwd = 3)
points(x = c(2, 4, 8, 16), y = c(smc_hier_results[[1]]$ESS['Q'], smc_hier_results[[3]]$ESS['Q'],
                                 smc_hier_results[[7]]$ESS['Q'], smc_hier_results[[15]]$ESS['Q']), 
       col = Okabe_Ito[5], pch = 0, lwd = 3)
lines(x = c(2, 4, 8, 16), y = c(smc_hier_results[[1]]$ESS['Q'], smc_hier_results[[3]]$ESS['Q'],
                                smc_hier_results[[7]]$ESS['Q'], smc_hier_results[[15]]$ESS['Q']),
      col = Okabe_Ito[5], lty = 2, lwd = 3)
points(x = 2:16, y = sapply(1:15, function(i) smc_prog_results[[i]]$ESS['Q']), 
       col = Okabe_Ito[4], pch = 2, lwd = 3)
lines(x = 2:16, y = sapply(1:15, function(i) smc_prog_results[[i]]$ESS['Q']), 
      col = Okabe_Ito[4], lty = 3, lwd = 3)
legend(x = 2, y = 10000, 
       legend = c('fork-and-join', 'balanced', 'progressive'),
       lty = c(1, 2, 3), 
       lwd = c(3, 3, 3),
       pch = c(1, 0, 2), 
       col = Okabe_Ito[c(8, 5, 4)],
       cex = 1.1,
       bty = 'n')


save.image('varying_C_experiments.RData')
