library(DCFusion)

seed <- 408
set.seed(seed)
denominator <- c(2, 4, 8, 16, 32, 64, 128, 256)
number_of_replicates <- 50
smc_fnj_results <- list()
smc_bal_results <- list()
smc_prog_results <- list()
nsamples <- 10000
opt_bw <- ((4/(3*nsamples))^(1/5))
bw <- 0.1
time_choice <- 1

for (i in 1:length(denominator)) {
  print(paste('i:', i))
  print(paste('C:', denominator[i]))
  smc_fnj_results[[i]] <- list()
  smc_bal_results[[i]] <- list()
  smc_prog_results[[i]] <- list()
  for (rep in 1:number_of_replicates) {
    print(paste('rep:', rep))
    set.seed(seed*rep*i)
    input_samples <- lapply(1:denominator[i], function(c) {
      rnorm_tempered(N = nsamples, mean = 0, sd = 1, beta = 1/denominator[i])})
    # preconditioned fork and join
    print('performing preconditioned MC fusion')
    smc_fnj_fused <-  hierarchical_fusion_SMC_uniGaussian(N_schedule = nsamples,
                                                          m_schedule = denominator[i],
                                                          time_schedule = time_choice,
                                                          base_samples = input_samples,
                                                          L = 2,
                                                          mean = 0,
                                                          sd = 1,
                                                          start_beta = 1/denominator[i], 
                                                          precondition = TRUE,
                                                          resampling_method = 'resid',
                                                          ESS_threshold = 0.5,
                                                          diffusion_estimator = 'NB',
                                                          seed = seed*rep*i)
    resampled_fnj <- resample_particle_y_samples(particle_set = smc_fnj_fused$particles[[1]],
                                                 multivariate = FALSE,
                                                 resampling_method = 'resid',
                                                 seed = seed*rep*i)$y_samples
    smc_fnj_results[[i]][[rep]] <- list('time' = smc_fnj_fused$time[[1]],
                                        'ESS' = smc_fnj_fused$ESS[[1]],
                                        'CESS' = smc_fnj_fused$CESS[[1]],
                                        'IAD_bw' = integrated_abs_distance_uniGaussian(fusion_post = resampled_fnj,
                                                                                       mean = 0,
                                                                                       sd = 1,
                                                                                       beta = 1,
                                                                                       bw = bw),
                                        'IAD_opt_bw' = integrated_abs_distance_uniGaussian(fusion_post = resampled_fnj,
                                                                                           mean = 0,
                                                                                           sd = 1,
                                                                                           beta = 1,
                                                                                           bw = opt_bw))
    
    if (denominator[i]==2) {
      # C=2 is the same for all tree hierarchies
      smc_bal_results[[i]][[rep]] <- smc_fnj_results[[i]][[rep]]
      smc_prog_results[[i]][[rep]] <- smc_fnj_results[[i]][[rep]]
    } else {
      # hierarchical
      if (log(denominator[i], 2)%%1==0) {
        print('performing hierarchical MC fusion')
        smc_bal_fused <- hierarchical_fusion_SMC_uniGaussian(N_schedule = rep(nsamples, log(denominator[i], 2)),
                                                             m_schedule = rep(2, log(denominator[i], 2)),
                                                             time_schedule = rep(time_choice, log(denominator[i], 2)),
                                                             base_samples = input_samples,
                                                             L = log(denominator[i], 2)+1,
                                                             mean = 0,
                                                             sd = 1,
                                                             start_beta = 1/denominator[i],
                                                             precondition = TRUE,
                                                             resampling_method = 'resid',
                                                             ESS_threshold = 0.5,
                                                             diffusion_estimator = 'NB',
                                                             seed = seed*rep*i)
        resampled_bal <- resample_particle_y_samples(particle_set = smc_bal_fused$particles[[1]],
                                                     multivariate = FALSE,
                                                     resampling_method = 'resid',
                                                     seed = seed*rep*i)$y_samples
        smc_bal_results[[i]][[rep]] <- list('time' = sum(unlist(smc_bal_fused$time)),
                                            'ESS' = smc_bal_fused$ESS,
                                            'CESS' = smc_bal_fused$CESS,
                                            'IAD_bw' = integrated_abs_distance_uniGaussian(fusion_post = resampled_bal,
                                                                                           mean = 0,
                                                                                           sd = 1,
                                                                                           beta = 1,
                                                                                           bw = bw),
                                            'IAD_opt_bw' = integrated_abs_distance_uniGaussian(fusion_post = resampled_bal,
                                                                                               mean = 0,
                                                                                               sd = 1,
                                                                                               beta = 1,
                                                                                               bw = opt_bw))
      } else {
        smc_bal_results[[i]][[rep]] <- NA
      }
      # progressive
      print('performing progressive MC fusion')
      smc_prog_fused <- progressive_fusion_SMC_uniGaussian(N_schedule = rep(nsamples, denominator[i]-1),
                                                           time_schedule = rep(time_choice, denominator[i]-1),
                                                           base_samples = input_samples, 
                                                           mean = 0,
                                                           sd = 1,
                                                           start_beta = 1/denominator[i],
                                                           precondition = TRUE, 
                                                           resampling_method = 'resid',
                                                           ESS_threshold = 0.5,
                                                           diffusion_estimator = 'NB',
                                                           seed = seed*rep*i)
      resampled_prog <- resample_particle_y_samples(particle_set = smc_prog_fused$particles[[1]],
                                                    multivariate = FALSE,
                                                    resampling_method = 'resid',
                                                    seed = seed*rep*i)$y_samples
      smc_prog_results[[i]][[rep]] <- list('time' = sum(unlist(smc_prog_fused$time)),
                                           'ESS' = smc_prog_fused$ESS,
                                           'CESS' = smc_prog_fused$CESS,
                                           'IAD_bw' = integrated_abs_distance_uniGaussian(fusion_post = resampled_prog,
                                                                                          mean = 0,
                                                                                          sd = 1,
                                                                                          beta = 1,
                                                                                          bw = bw),
                                           'IAD_opt_bw' = integrated_abs_distance_uniGaussian(fusion_post = resampled_prog,
                                                                                              mean = 0,
                                                                                              sd = 1,
                                                                                              beta = 1,
                                                                                              bw = opt_bw))
    }
    
    # ##########
    # curve(dnorm(x), -4, 4, ylim = c(0, 0.5), main = paste('Default BW || C =',denominator[i], '|| rep =', rep))
    # lines(density(smc_fnj_fused$particles[[1]]$y_samples), col = 'orange')
    # if (denominator[i]!=2) {
    #   if (!any(is.na(smc_bal_results[[i]][[rep]]))) {
    #     lines(density(smc_bal_fused$particles[[1]]$y_samples), col = 'green')
    #   }
    #   lines(density(smc_prog_fused$particles[[1]]$y_samples), col = 'blue')
    # }
    # ##########
    # curve(dnorm(x), -4, 4, ylim = c(0, 0.5), main = paste('BW =', bw, '|| C =',denominator[i], '|| rep =', rep))
    # lines(density(smc_fnj_fused$particles[[1]]$y_samples, bw = bw), col = 'orange')
    # if (denominator[i]!=2) {
    #   if (!any(is.na(smc_bal_results[[i]][[rep]]))) {
    #     lines(density(smc_bal_fused$particles[[1]]$y_samples, bw = bw), col = 'green')
    #   }
    #   lines(density(smc_prog_fused$particles[[1]]$y_samples, bw = bw), col = 'blue')
    # }

    print('saving progress')
    save.image('varying_C_experiments_uniG_smc_replicates.RData')
  }
}

######################################## running time

plot(x = denominator, y = sapply(1:length(denominator), function(i) {
  mean(sapply(1:number_of_replicates, function(j) smc_fnj_results[[i]][[j]]$time))
}), ylim = c(0, 2000), ylab = '', xlab = '', font.lab = 2, pch = 1, lwd = 3, xaxt = 'n', yaxt = 'n', type = 'l')
axis(1, at = denominator, labels = denominator, font = 2, cex = 1.5)
mtext('Number of sub-posteriors (C)', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 2000, 200), labels = seq(0, 2000, 200), font = 2, cex = 1.5)
mtext('Time elapsed in seconds', 2, 2.75, font = 2, cex = 1.5)
lines(x = denominator, y = sapply(1:length(denominator), function(i) {
  mean(sapply(1:number_of_replicates, function(j) smc_bal_results[[i]][[j]]$time))
}), lty = 3, lwd = 3)
lines(x = denominator, y = sapply(1:length(denominator), function(i) {
  mean(sapply(1:number_of_replicates, function(j) smc_prog_results[[i]][[j]]$time))
}), lty = 2, lwd = 3)
legend(x = 2, y = 2000,
       legend = c('fork-and-join', 'balanced', 'progressive'),
       lty = c(1, 3, 2),
       lwd = c(3, 3, 3),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

#################### log

plot(x = denominator, y = sapply(1:length(denominator), function(i) {
  mean(sapply(1:number_of_replicates, function(j) log(smc_fnj_results[[i]][[j]]$time)))
}), ylim = c(1, 8), ylab = '', xlab = '', font.lab = 2, pch = 1, lwd = 3, xaxt = 'n', yaxt = 'n', type = 'l')
axis(1, at = denominator, labels = denominator, font = 2, cex = 1.5)
mtext('Number of sub-posteriors (C)', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 20, 1), labels = seq(0, 20, 1), font = 2, cex = 1.5)
mtext('log(Time elapsed in seconds)', 2, 2.75, font = 2, cex = 1.5)
lines(x = denominator, y = sapply(1:length(denominator), function(i) {
  mean(sapply(1:number_of_replicates, function(j) log(smc_bal_results[[i]][[j]]$time)))
}), lty = 3, lwd = 3)
lines(x = denominator, y = sapply(1:length(denominator), function(i) {
  mean(sapply(1:number_of_replicates, function(j) log(smc_prog_results[[i]][[j]]$time)))
}), lty = 2, lwd = 3)
legend(x = 2, y = 8,
       legend = c('fork-and-join', 'balanced', 'progressive'),
       lty = c(1, 3, 2),
       lwd = c(3, 3, 3),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

######################################## ESS (overall)

plot(x = denominator, y = sapply(1:length(denominator), function(i) {
  mean(sapply(1:number_of_replicates, function(j) smc_fnj_results[[i]][[j]]$ESS['Q']))
}), ylim = c(0, 10000), ylab = '', xlab = '', font.lab = 2, pch = 1, lwd = 3, xaxt = 'n', yaxt = 'n', type = 'l')
axis(1, at = denominator, labels = denominator, font = 2, cex = 1.5)
mtext('Number of sub-posteriors (C)', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 10000, 1000), labels = seq(0, 10000, 1000), font = 2, cex = 1.5)
mtext('ESS (at top level)', 2, 2.75, font = 2, cex = 1.5)
lines(x = denominator, y = c(mean(sapply(1:number_of_replicates, function(j) smc_bal_results[[1]][[j]]$ESS['Q'])),
                             sapply(2:length(denominator), function(i) {
                               mean(sapply(1:number_of_replicates, function(j) smc_bal_results[[i]][[j]]$ESS[[1]]['Q']))
                             })), lty = 3, lwd = 3)
lines(x = denominator, y = c(mean(sapply(1:number_of_replicates, function(j) smc_prog_results[[1]][[j]]$ESS['Q'])),
                             sapply(2:length(denominator), function(i) {
                               mean(sapply(1:number_of_replicates, function(j) smc_prog_results[[i]][[j]]$ESS[[1]]['Q']))
                             })), lty = 2, lwd = 3)
legend(x = 2, y = 10000,
       legend = c('fork-and-join', 'balanced', 'progressive'),
       lty = c(1, 3, 2),
       lwd = c(3, 3, 3),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

######################################## IAD

plot(x = denominator, y = sapply(1:length(denominator), function(i) {
  mean(sapply(1:number_of_replicates, function(j) smc_fnj_results[[i]][[j]]$IAD_bw))
}), ylim = c(0, 0.6), ylab = '', xlab = '', font.lab = 2, pch = 1, lwd = 3, xaxt = 'n', yaxt = 'n', type = 'l')
axis(1, at = denominator, labels = denominator, font = 2, cex = 1.5)
mtext('Number of sub-posteriors (C)', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 0.8, 0.05),
     labels = c("0.0", 0.05, "0.10", 0.15, "0.20", 0.25, "0.30", 0.35, "0.40", 0.45, "0.50", 0.55, "0.60", 0.65, "0.7", 0.75, "0.80"),
     font = 2, cex = 1.5)
mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
lines(x = denominator, y = sapply(1:length(denominator), function(i) {
  mean(sapply(1:number_of_replicates, function(j) smc_fnj_results[[i]][[j]]$IAD_bw))
}), lwd = 3)
lines(x = denominator, y = sapply(1:length(denominator), function(i) {
  mean(sapply(1:number_of_replicates, function(j) smc_bal_results[[i]][[j]]$IAD_bw))
}), lty = 3, lwd = 3)
lines(x = denominator, y = sapply(1:length(denominator), function(i) {
  mean(sapply(1:number_of_replicates, function(j) smc_prog_results[[i]][[j]]$IAD_bw))
}), lty = 2, lwd = 3)
legend(x = 2, y = 0.6,
       legend = c('fork-and-join', 'balanced', 'progressive'),
       lty = c(1, 3, 2),
       lwd = c(3, 3, 3),
       cex = 1.25,
       text.font = 2,
       bty = 'n')
