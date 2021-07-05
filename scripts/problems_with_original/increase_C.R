library(hierarchicalFusion)

seed <- 408
denominator <- 2:16
input_samples <- list()
fnj_results <- list()
time_choice <- 1

for (i in 1:length(denominator)) {
  print(denominator[i])
  set.seed(seed)
  input_samples[[i]] <- base_rejection_sampler_exp_4(beta = 1/denominator[i],
                                                     nsamples = 10000,
                                                     proposal_mean = 0,
                                                     proposal_sd = 1.5,
                                                     dominating_M = 1.75)
  
  curve(exp_4_density(x, beta = 1/denominator[i]), -4, 4,
        main = denominator[i], ylab = 'tempered pdf')
  for (j in 1:length(input_samples[[i]])) {
    lines(density(input_samples[[i]][[j]]), col = 'black')
  }
  
  # standard fork and join
  print('performing standard fork-and-join MC fusion')
  fnj_fused <- hierarchical_fusion_exp_4(N_schedule = 10000,
                                         m_schedule = denominator[i],
                                         time_schedule = time_choice,
                                         base_samples = input_samples[[i]],
                                         mean = 0,
                                         start_beta = 1/denominator[i],
                                         L = 2,
                                         precondition = FALSE,
                                         seed = seed)
  
  fnj_results[[i]] <- list('time' = fnj_fused$overall_time,
                           'overall_rho' = fnj_fused$overall_rho,
                           'overall_Q' = fnj_fused$overall_Q,
                           'overall_rhoQ' = fnj_fused$overall_rhoQ)
}

par(mai = c(1.02, 1, 0.82, 0.42))

######################################## running time

plot(x = 2:16, y = sapply(1:15, function(i) fnj_results[[i]][[1]]), ylim = c(0, 10),
     ylab = 'Running time in seconds', xlab = 'Number of Subposteriors (C)', col = 'black', pch = 4)
lines(x = 2:16, y = sapply(1:15, function(i) fnj_results[[i]][[1]]), col = 'black')

#################### log

# uses colours defined in colour_scales.R

plot(x = 2:16, y = sapply(1:15, function(i) log(fnj_results[[i]][[1]])), ylim = c(-1, 10),
     ylab = 'Time Elapsed in log(seconds)', xlab = 'Number of Subposteriors (C)',
     col = Okabe_Ito[8], pch = 1, lwd = 3)
lines(x = 2:16, y = sapply(1:15, function(i) log(fnj_results[[i]][[1]])), 
      col = Okabe_Ito[8], lwd = 3)

######################################## rho acceptance

plot(x = 2:16, y = sapply(1:15, function(i) fnj_results[[i]]$overall_rho), ylim = c(0, 1),
     ylab = expression(paste('Acceptance Rate for ', rho)), xlab = 'Number of Subposteriors (C)', 
     col = 'black', pch = 1, lwd = 3)
lines(x = 2:16, y = sapply(1:15, function(i) fnj_results[[i]]$overall_rho), 
      col = 'black', lwd = 3)

######################################## Q acceptance

plot(x = 2:16, y = sapply(1:15, function(i) fnj_results[[i]]$overall_Q), ylim = c(0, 1),
     ylab = expression(paste('Acceptance Rate for ', hat(Q))), xlab = 'Number of Subposteriors (C)', 
     col = 'black', pch = 1, lwd = 3)
lines(x = 2:16, y = sapply(1:15, function(i) fnj_results[[i]]$overall_Q), 
      col = 'black', lwd = 3)

