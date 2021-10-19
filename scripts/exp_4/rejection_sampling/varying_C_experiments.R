library(DCFusion)

seed <- 408
set.seed(seed)
denominator <- 2:16
input_samples <- list()
fnj_results <- list()
bal_results <- list()
prog_results <- list()
time_choice <- 1
nsamples <- 10000

for (i in 1:length(denominator)) {
  print(denominator[i])
  input_samples[[i]] <- base_rejection_sampler_exp_4(beta = 1/denominator[i],
                                                     nsamples = nsamples,
                                                     proposal_mean = 0,
                                                     proposal_sd = 1.5,
                                                     dominating_M = 1.75)
  
  curve(exp_4_density(x, beta = 1/denominator[i]), -4, 4,
        main = denominator[i], ylab = 'tempered pdf')
  for (j in 1:length(input_samples[[i]])) {
    lines(density(input_samples[[i]][[j]]), col = 'black')
  }
  
  # standard fork and join
  print('performing preconditioned fork-and-join MC fusion')
  fnj_fused <- bal_binary_fusion_exp_4(N_schedule = nsamples,
                                       m_schedule = denominator[i],
                                       time_schedule = time_choice,
                                       base_samples = input_samples[[i]],
                                       mean = 0,
                                       start_beta = 1/denominator[i],
                                       L = 2,
                                       precondition = TRUE,
                                       seed = seed)
  
  fnj_results[[i]] <- list('time' = fnj_fused$overall_time,
                           'overall_rho' = fnj_fused$overall_rho,
                           'overall_Q' = fnj_fused$overall_Q,
                           'overall_rhoQ' = fnj_fused$overall_rhoQ)
  
  # balanced binary if denominator[i] is 2, 4, 8, or 16
  if (denominator[i]==2) {
    print('performing preconditioned balanced binary MC fusion')
    bal_fused <- bal_binary_fusion_exp_4(N_schedule = nsamples,
                                         m_schedule = 2,
                                         time_schedule = time_choice,
                                         base_samples = input_samples[[i]],
                                         mean = 0,
                                         start_beta = 1/2,
                                         L = 2,
                                         precondition = TRUE,
                                         seed = seed)
  } else if (denominator[i]==4) {
    print('performing preconditioned balanced binary MC fusion')
    bal_fused <- bal_binary_fusion_exp_4(N_schedule = rep(nsamples, 2),
                                         m_schedule = rep(2, 2),
                                         time_schedule = rep(time_choice, 2),
                                         base_samples = input_samples[[i]],
                                         mean = 0,
                                         start_beta = 1/4,
                                         L = 3,
                                         precondition = TRUE,
                                         seed = seed)
  } else if (denominator[i]==8) {
    print('performing preconditioned balanced binary MC fusion')
    bal_fused <- bal_binary_fusion_exp_4(N_schedule = rep(nsamples, 3),
                                         m_schedule = rep(2, 3),
                                         time_schedule = rep(time_choice, 3),
                                         base_samples = input_samples[[i]],
                                         mean = 0,
                                         start_beta = 1/8,
                                         L = 4,
                                         precondition = TRUE,
                                         seed = seed)
  } else if (denominator[i]==16) {
    print('performing preconditioned balanced binary MC fusion')
    bal_fused <- bal_binary_fusion_exp_4(N_schedule = rep(nsamples, 4),
                                         m_schedule = rep(2, 4),
                                         time_schedule = rep(time_choice, 4),
                                         base_samples = input_samples[[i]], 
                                         mean = 0,
                                         start_beta = 1/16,
                                         L = 5,
                                         precondition = TRUE,
                                         seed = seed)
  }
  
  if (denominator[i] %in% c(2, 4, 8, 16)) {
    bal_results[[i]] <- list('time' = bal_fused$overall_time,
                             'overall_rho' = bal_fused$overall_rho,
                             'overall_Q' = bal_fused$overall_Q,
                             'overall_rhoQ' = bal_fused$overall_rhoQ)
  } else {
    bal_results[[i]] <- NA
  }
  
  # progressive
  print('performing preconditioned progressive MC fusion')
  prog_fused <- progressive_fusion_exp_4(N_schedule = rep(nsamples, denominator[i]-1),
                                         time_schedule = rep(time_choice, denominator[i]-1),
                                         base_samples = input_samples[[i]], 
                                         mean = 0,
                                         start_beta = 1/denominator[i],
                                         precondition = TRUE,
                                         seed = seed)
  
  prog_results[[i]] <- list('time' = prog_fused$time,
                            'overall_rho' = prog_fused$rho_acc,
                            'overall_Q' = prog_fused$Q_acc,
                            'overall_rhoQ' = prog_fused$rhoQ_acc)
  
  curve(exp_4_density(x), -4, 4, ylim = c(0, 0.5), main = denominator[i])
  lines(density(fnj_fused$samples[[1]]), col = 'orange', lty = 2)
  if (!any(is.na(bal_results[[i]]))) {
    lines(density(bal_fused$samples[[1]]), col = 'blue', lty = 2)
  }
  lines(density(prog_fused$samples[[1]]), col = 'darkgreen', lty = 2)
}

par(mai = c(1.02, 1, 0.82, 0.42))

######################################## running time

plot(x = 2:16, y = sapply(1:15, function(i) fnj_results[[i]][[1]]), ylim = c(0, 10),
     ylab = 'Running time in seconds', xlab = 'Number of Subposteriors (C)', col = 'black', pch = 4)
lines(x = 2:16, y = sapply(1:15, function(i) fnj_results[[i]][[1]]), col = 'black')
points(x = c(2, 4, 8, 16), y = c(sum(bal_results[[1]]$time), sum(bal_results[[3]]$time),
                                 sum(bal_results[[7]]$time), sum(bal_results[[15]]$time)), col = 'black', pch = 4)
lines(x = c(2, 4, 8, 16), y = c(sum(bal_results[[1]]$time), sum(bal_results[[3]]$time),
                                sum(bal_results[[7]]$time), sum(bal_results[[15]]$time)), col = 'black')
points(x = 2:16, y = sapply(1:15, function(i) sum(prog_results[[i]]$time)), col = 'black', pch = 4)
lines(x = 2:16, y = sapply(1:15, function(i) sum(prog_results[[i]]$time)), col = 'black')

#################### log

Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
plot(x = 2:16, y = sapply(1:15, function(i) log(fnj_results[[i]][[1]])), ylim = c(-1, 12),
     ylab = '', xlab = '', font.lab = 2,
     col = Okabe_Ito[8], pch = 1, lwd = 3)
axis(1, at = seq(2, 16, 2), labels = seq(2, 16, 2), font = 2, cex = 1.5)
axis(1, at=0:16, labels=rep("", 17), lwd.ticks = 0.5)
mtext('Number of sub-posteriors (C)', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 12, 2), labels = seq(0, 12, 2), font = 2, cex = 1.5)
axis(2, at=0:12, labels=rep("", 13), lwd.ticks = 0.5)
mtext('Time Elapsed in log(seconds)', 2, 2.75, font = 2, cex = 1.5)
lines(x = 2:16, y = sapply(1:15, function(i) log(fnj_results[[i]][[1]])), 
      col = Okabe_Ito[8], lwd = 3)
points(x = 2:16, y = sapply(1:15, function(i) log(sum(prog_results[[i]]$time))), 
       col = Okabe_Ito[4], pch = 2, lwd = 3)
lines(x = c(2, 4, 8, 16), y = log(c(sum(bal_results[[1]]$time), sum(bal_results[[3]]$time),
                                    sum(bal_results[[7]]$time), sum(bal_results[[15]]$time))), 
      col = Okabe_Ito[5], lty = 2, lwd = 3)
lines(x = 2:16, y = sapply(1:15, function(i) log(sum(prog_results[[i]]$time))), 
      col = Okabe_Ito[4], lty = 3, lwd = 3)
points(x = c(2, 4, 8, 16), y = log(c(sum(bal_results[[1]]$time), sum(bal_results[[3]]$time),
                                     sum(bal_results[[7]]$time), sum(bal_results[[15]]$time))), 
       col = Okabe_Ito[5], pch = 0, lwd = 3)
legend(x = 2, y = 12, 
       legend = c('fork-and-join', 'balanced', 'progressive'),
       lty = c(1, 2, 3), 
       lwd = c(3, 3, 3),
       pch = c(1, 0, 2), 
       col = Okabe_Ito[c(8, 5, 4)],
       cex = 1.25,
       text.font = 2,
       bty = 'n')

####################  black

plot(x = 2:16, y = sapply(1:15, function(i) log(fnj_results[[i]][[1]])), ylim = c(-1, 12),
     ylab = '', xlab = '', font.lab = 2, pch = 1, lwd = 3)
axis(1, at = seq(2, 16, 2), labels = seq(2, 16, 2), font = 2, cex = 1.5)
axis(1, at=0:16, labels=rep("", 17), lwd.ticks = 0.5)
mtext('Number of sub-posteriors (C)', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 12, 2), labels = seq(0, 12, 2), font = 2, cex = 1.5)
axis(2, at=0:12, labels=rep("", 13), lwd.ticks = 0.5)
mtext('Time Elapsed in log(seconds)', 2, 2.75, font = 2, cex = 1.5)
lines(x = 2:16, y = sapply(1:15, function(i) log(fnj_results[[i]][[1]])), lwd = 3)
points(x = 2:16, y = sapply(1:15, function(i) log(sum(prog_results[[i]]$time))), 
       pch = 2, lwd = 3)
lines(x = c(2, 4, 8, 16), y = log(c(sum(bal_results[[1]]$time), sum(bal_results[[3]]$time),
                                    sum(bal_results[[7]]$time), sum(bal_results[[15]]$time))), 
      lty = 3, lwd = 3)
lines(x = 2:16, y = sapply(1:15, function(i) log(sum(prog_results[[i]]$time))), 
      lty = 2, lwd = 3)
points(x = c(2, 4, 8, 16), y = log(c(sum(bal_results[[1]]$time), sum(bal_results[[3]]$time),
                                     sum(bal_results[[7]]$time), sum(bal_results[[15]]$time))), 
       pch = 0, lwd = 3)
legend(x = 2, y = 12, 
       legend = c('fork-and-join', 'balanced', 'progressive'),
       lty = c(1, 3, 2), 
       lwd = c(3, 3, 3),
       pch = c(1, 0, 2), 
       cex = 1.25,
       text.font = 2,
       bty = 'n')


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

