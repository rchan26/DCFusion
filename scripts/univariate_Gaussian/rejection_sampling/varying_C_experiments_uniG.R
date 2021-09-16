library(hierarchicalFusion)

seed <- 408
set.seed(seed)
denominator <- 2:16
fnj_standard_results <- list()
fnj_precondition_results <- list()
hier_results <- list()
prog_results <- list()
time_choice <- 1
nsamples <- 10000

for (i in 1:length(denominator)) {
  set.seed(seed)
  print(paste('i:', i))
  print(paste('C:', denominator[i]))
  input_samples <- lapply(1:denominator[i], function(c) {
    rnorm_tempered(N = nsamples, mean = 0, sd = 1, beta = 1/denominator[i])})
  # standard fork and join
  print('performing standard fork-and-join MC fusion')
  fnj_standard_fused <- hierarchical_fusion_uniGaussian(N_schedule = nsamples,
                                                        m_schedule = denominator[i],
                                                        time_schedule = time_choice,
                                                        base_samples = input_samples,
                                                        L = 2,
                                                        mean = 0,
                                                        sd = 1,
                                                        start_beta = 1/denominator[i],
                                                        precondition = FALSE,
                                                        seed = seed)
  fnj_standard_results[[i]] <- list('time' = fnj_standard_fused$overall_time,
                                    'overall_rho' = fnj_standard_fused$overall_rho,
                                    'overall_Q' = fnj_standard_fused$overall_Q,
                                    'overall_rhoQ' = fnj_standard_fused$overall_rhoQ)
  # preconditioned fork and join
  print('performing preconditioned fork-and-join MC fusion')
  fnj_precondition_fused <- hierarchical_fusion_uniGaussian(N_schedule = nsamples,
                                                            m_schedule = denominator[i],
                                                            time_schedule = time_choice,
                                                            base_samples = input_samples,
                                                            L = 2,
                                                            mean = 0,
                                                            sd = 1,
                                                            start_beta = 1/denominator[i],
                                                            precondition = TRUE,
                                                            seed = seed)
  fnj_precondition_results[[i]] <- list('time' = fnj_precondition_fused$overall_time,
                                        'overall_rho' = fnj_precondition_fused$overall_rho,
                                        'overall_Q' = fnj_precondition_fused$overall_Q,
                                        'overall_rhoQ' = fnj_precondition_fused$overall_rhoQ)
  # hierarchical if denominator[i] is 2, 4, 8, or 16
  if (denominator[i]==2) {
    # C=2 is the same for all tree hierarchies
    hier_results[[i]] <- fnj_precondition_results[[i]]
    prog_results[[i]] <- fnj_precondition_results[[i]]
  } else {
    # hierarchical
    if (log(denominator[i], 2)%%1==0) {
      print('performing hierarchical MC fusion')
      hier_fused <- hierarchical_fusion_uniGaussian(N_schedule = rep(nsamples, log(denominator[i], 2)),
                                                    m_schedule = rep(2, log(denominator[i], 2)),
                                                    time_schedule = rep(time_choice, log(denominator[i], 2)),
                                                    base_samples = input_samples,
                                                    L = log(denominator[i], 2)+1,
                                                    mean = 0,
                                                    sd = 1,
                                                    start_beta = 1/denominator[i],
                                                    precondition = TRUE,
                                                    seed = seed)
      hier_results[[i]] <- list('time' = hier_fused$overall_time,
                                'overall_rho' = hier_fused$overall_rho,
                                'overall_Q' = hier_fused$overall_Q,
                                'overall_rhoQ' = hier_fused$overall_rhoQ)
    
    } else {
      hier_results[[i]] <- NA
    }
    # progressive
    print('performing preconditioned progressive MC fusion')
    prog_fused <- progressive_fusion_uniGaussian(N_schedule = rep(nsamples, denominator[i]-1),
                                                 time_schedule = rep(time_choice, denominator[i]-1),
                                                 base_samples = input_samples, 
                                                 mean = 0,
                                                 sd = 1,
                                                 start_beta = 1/denominator[i],
                                                 precondition = TRUE,
                                                 seed = seed)
    prog_results[[i]] <- list('time' = prog_fused$time,
                              'overall_rho' = prog_fused$rho_acc,
                              'overall_Q' = prog_fused$Q_acc,
                              'overall_rhoQ' = prog_fused$rhoQ_acc)
  }
  # curve(dnorm(x), -4, 4, ylim = c(0, 0.5), main = denominator[i])
  # lines(density(fnj_standard_fused$samples[[1]]), col = 'orange', lty = 2)
  # lines(density(fnj_precondition_fused$samples[[1]]), col = 'orange', lty = 2)
  # if (denominator[i]!=2) {
  #   if (!any(is.na(hier_results[[i]]))) {
  #     lines(density(hier_fused$samples[[1]]), col = 'blue', lty = 2)
  #   }
  #   lines(density(prog_fused$samples[[1]]), col = 'darkgreen', lty = 2)
  # }
  print('saving progress')
  save.image('varying_C_experiments_uniG.RData')
}

######################################## running time

plot(x = 2:16, y = sapply(1:15, function(i) fnj_precondition_results[[i]][[1]]), ylim = c(0, 10),
     ylab = 'Running time in seconds', xlab = 'Number of Subposteriors (C)', col = 'black', pch = 4)
lines(x = 2:16, y = sapply(1:15, function(i) fnj_precondition_results[[i]][[1]]), col = 'black')
points(x = c(2, 4, 8, 16), y = c(sum(hier_results[[1]]$time), sum(hier_results[[3]]$time),
                                 sum(hier_results[[7]]$time), sum(hier_results[[15]]$time)), col = 'black', pch = 4)
lines(x = c(2, 4, 8, 16), y = c(sum(hier_results[[1]]$time), sum(hier_results[[3]]$time),
                                sum(hier_results[[7]]$time), sum(hier_results[[15]]$time)), col = 'black')
points(x = 2:16, y = sapply(1:15, function(i) sum(prog_results[[i]]$time)), col = 'black', pch = 4)
lines(x = 2:16, y = sapply(1:15, function(i) sum(prog_results[[i]]$time)), col = 'black')

#################### log

Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
plot(x = 2:16, y = sapply(1:15, function(i) log(fnj_standard_results[[i]][[1]])), ylim = c(-1, 12),
     ylab = '', xlab = '', font.lab = 2,
     col = Okabe_Ito[8], pch = 1, lwd = 3)
axis(1, at = seq(2, 16, 2), labels = seq(2, 16, 2), font = 2, cex = 1.5)
axis(1, at=0:16, labels=rep("", 17), lwd.ticks = 0.5)
mtext('Number of sub-posteriors (C)', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 12, 2), labels = seq(0, 12, 2), font = 2, cex = 1.5)
axis(2, at=0:12, labels=rep("", 13), lwd.ticks = 0.5)
mtext('Time Elapsed in log(seconds)', 2, 2.75, font = 2, cex = 1.5)
points(x = 2:16, y = sapply(1:15, function(i) log(fnj_precondition_results[[i]][[1]])),
      col = Okabe_Ito[8], lwd = 3)
lines(x = 2:16, y = sapply(1:15, function(i) log(fnj_precondition_results[[i]][[1]])),
      col = Okabe_Ito[8], lwd = 3)
lines(x = 2:16, y = sapply(1:15, function(i) log(fnj_standard_results[[i]][[1]])),
      col = Okabe_Ito[8], lwd = 3)
points(x = 2:16, y = sapply(1:15, function(i) log(sum(prog_results[[i]]$time))),
       col = Okabe_Ito[4], pch = 2, lwd = 3)
lines(x = c(2, 4, 8, 16), y = log(c(sum(hier_results[[1]]$time), sum(hier_results[[3]]$time),
                                    sum(hier_results[[7]]$time), sum(hier_results[[15]]$time))),
      col = Okabe_Ito[5], lty = 2, lwd = 3)
lines(x = 2:16, y = sapply(1:15, function(i) log(sum(prog_results[[i]]$time))),
      col = Okabe_Ito[4], lty = 3, lwd = 3)
points(x = c(2, 4, 8, 16), y = log(c(sum(hier_results[[1]]$time), sum(hier_results[[3]]$time),
                                     sum(hier_results[[7]]$time), sum(hier_results[[15]]$time))),
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

plot(x = 2:16, y = sapply(1:15, function(i) log(fnj_precondition_results[[i]][[1]])), ylim = c(-1, 7),
     ylab = '', xlab = '', font.lab = 2, pch = 1, lwd = 3, xaxt='n', yaxt='n')
axis(1, at = seq(2, 16, 2), labels = seq(2, 16, 2), font = 2, cex = 1.5)
axis(1, at=0:16, labels=rep("", 17), lwd.ticks = 0.5)
mtext('Number of sub-posteriors (C)', 1, 2.75, font = 2, cex = 1.5)
axis(2, at=-1:7, labels=-1:7, font = 2, cex = 1.5)
# axis(2, at = seq(0, 8, 2), labels = seq(0, 8, 2), font = 2, cex = 1.5)
# axis(2, at=0:8, labels=rep("", 9), lwd.ticks = 0.5)
mtext('log(Time elapsed in seconds)', 2, 2.75, font = 2, cex = 1.5)
lines(x = 2:16, y = sapply(1:15, function(i) log(fnj_precondition_results[[i]][[1]])), lwd = 3)
points(x = 2:16, y = sapply(1:15, function(i) log(sum(prog_results[[i]]$time))),
       pch = 2, lwd = 3)
lines(x = c(2, 4, 8, 16), y = log(c(sum(hier_results[[1]]$time), sum(hier_results[[3]]$time),
                                    sum(hier_results[[7]]$time), sum(hier_results[[15]]$time))),
      lty = 3, lwd = 3)
lines(x = 2:16, y = sapply(1:15, function(i) log(sum(prog_results[[i]]$time))),
      lty = 2, lwd = 3)
points(x = c(2, 4, 8, 16), y = log(c(sum(hier_results[[1]]$time), sum(hier_results[[3]]$time),
                                     sum(hier_results[[7]]$time), sum(hier_results[[15]]$time))),
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

plot(x = 2:16, y = sapply(1:15, function(i) fnj_precondition_results[[i]]$overall_rho), ylim = c(0, 1),
     ylab = expression(paste('Acceptance Rate for ', rho)), xlab = 'Number of Subposteriors (C)',
     col = 'black', pch = 1, lwd = 3)
lines(x = 2:16, y = sapply(1:15, function(i) fnj_precondition_results[[i]]$overall_rho),
      col = 'black', lwd = 3)

######################################## Q acceptance

plot(x = 2:16, y = sapply(1:15, function(i) fnj_precondition_results[[i]]$overall_Q), ylim = c(0, 1),
     ylab = expression(paste('Acceptance Rate for ', hat(Q))), xlab = 'Number of Subposteriors (C)',
     col = 'black', pch = 1, lwd = 3)
lines(x = 2:16, y = sapply(1:15, function(i) fnj_precondition_results[[i]]$overall_Q),
      col = 'black', lwd = 3)
