library(DCFusion)

seed <- 1994
set.seed(seed)
number_of_replicates <- 50
mu_choices <- c(1, 2, 4, 8, 16, 32, 64, 128)
beta_choices <- c(1/2, 1/4, 1/8, 1/16)
indices <- list(1:7, 1:8, 1:8, 1:8)
sd <- 1
nsamples <- 10000
time_choice <- 1
target_param <- c('mu' = 0, 'sd' = sqrt(1/2))
opt_bw <- ((4/(3*nsamples))^(1/5))*sqrt(1/2)
bw <- 0.1
precondition_bal_results <- list()
# precondition_prog_results <- list()

for (b in 1:length(beta_choices)) {
  print(paste('b:', b))
  print(paste('beta:', beta_choices[b]))
  precondition_bal_results[[b]] <- list()
  # precondition_prog_results[[b]] <- list()
  for (i in indices[[b]]) {
    print(paste('i:', i))
    print(paste('mu:', mu_choices[i]))
    print(paste('difference in means:', 2*mu_choices[i]))
    precondition_bal_results[[b]][[i]] <- list()
    # precondition_prog_results[[b]][[i]] <- list()
    for (rep in 1:number_of_replicates) {
      print(paste('rep:', rep))
      set.seed(seed*rep*b*i)
      # fusion with preconditioning
      print('fusion with preconditioning (with tempering)')
      input_samples <- lapply(1:(1/beta_choices[b]), function(c) {
        list(rnorm_tempered(N = nsamples, mean = -mu_choices[i], sd = sd, beta = beta_choices[b]),
             rnorm_tempered(N = nsamples, mean = mu_choices[i], sd = sd, beta = beta_choices[b]))})
      input_particles <- lapply(input_samples, function(samples) initialise_particle_sets(samples, multivariate = FALSE))
      # part 1: obtaining particle approximation for f^beta using standard fusion
      print('fusion for f^beta')
      precondition_p1 <- lapply(1:(1/beta_choices[b]), function(c) {
        parallel_fusion_SMC_uniGaussian(particles_to_fuse = input_particles[[c]], 
                                        N = nsamples,
                                        m = 2,
                                        time = time_choice,
                                        means = c(-mu_choices[i], mu_choices[i]),
                                        sds = c(sd, sd), 
                                        betas = rep(beta_choices[b], 2),
                                        precondition_values = sapply(input_samples[[c]], var),
                                        ESS_threshold = 0.5,
                                        resampling_method = 'resid',
                                        seed = seed*rep*b*i)})
      # part 2a: obtaining particle approximation for f using D&C fusion (balanced binary)
      print('fusion for f with balanced binary tree')
      precondition_bal_p2 <- bal_binary_fusion_SMC_uniGaussian(N_schedule = rep(nsamples, log((1/beta_choices[b]), 2)),
                                                               m_schedule = rep(2, log((1/beta_choices[b]), 2)),
                                                               time_schedule = rep(time_choice, log((1/beta_choices[b]), 2)),
                                                               base_samples = lapply(1:(1/beta_choices[b]), function(c) precondition_p1[[c]]$particles),
                                                               L = log((1/beta_choices[b]), 2)+1,
                                                               mean = target_param['mu'],
                                                               sd = target_param['sd'],
                                                               start_beta = beta_choices[b],
                                                               precondition = lapply(1:(1/beta_choices[b]), function(c) precondition_p1[[c]]$precondition_values[[1]]),
                                                               resampling_method = 'resid',
                                                               ESS_threshold = 0.5,
                                                               diffusion_estimator = 'NB',
                                                               seed = seed*rep*b*i)
      resampled_precondition_bal <- resample_particle_y_samples(particle_set = precondition_bal_p2$particles[[1]],
                                                                multivariate = FALSE,
                                                                resampling_method = 'resid',
                                                                seed = seed*rep*b*i)$y_samples
      precondition_bal_results[[b]][[i]][[rep]] <- list('time' = sum(sapply(1:(1/beta_choices[b]), function(c) precondition_p1[[c]]$time))+sum(unlist(precondition_bal_p2$time)),
                                                        'ESS' = list(lapply(1:(1/beta_choices[b]), function(c) precondition_p1[[c]]$ESS), precondition_bal_p2$ESS),
                                                        'CESS' = list(lapply(1:(1/beta_choices[b]), function(c) precondition_p1[[c]]$CESS), precondition_bal_p2$CESS),
                                                        'IAD_bw' = integrated_abs_distance_uniGaussian(fusion_post = resampled_precondition_bal,
                                                                                                       mean = target_param['mu'],
                                                                                                       sd = target_param['sd'],
                                                                                                       beta = 1,
                                                                                                       bw = bw),
                                                        'IAD_opt_bw' = integrated_abs_distance_uniGaussian(fusion_post = resampled_precondition_bal,
                                                                                                           mean = target_param['mu'],
                                                                                                           sd = target_param['sd'],
                                                                                                           beta = 1,
                                                                                                           bw = opt_bw))
      # if (beta_choices[b]==1/2) {
      #   # the balanced binary and progressive trees are exactly the same
      #   precondition_prog_results[[b]][[i]][[rep]] <- precondition_bal_results[[b]][[i]][[rep]]
      # } else {
      #   # part 2b: obtaining particle approximation for f using D&C fusion (progressive)
      #   print('fusion for f with progressive tree')
      #   precondition_prog_p2 <- progressive_fusion_SMC_uniGaussian(N_schedule = rep(nsamples, (1/beta_choices[b])-1),
      #                                                              time_schedule = rep(time_choice, (1/beta_choices[b])-1),
      #                                                              base_samples = lapply(1:(1/beta_choices[b]), function(c) precondition_p1[[c]]$particles),
      #                                                              mean = target_param['mu'],
      #                                                              sd = target_param['sd'],
      #                                                              start_beta = beta_choices[b],
      #                                                              precondition = lapply(1:(1/beta_choices[b]), function(c) precondition_p1[[c]]$precondition_values[[1]]),
      #                                                              resampling_method = 'resid',
      #                                                              ESS_threshold = 0.5,
      #                                                              diffusion_estimator = 'NB',
      #                                                              seed = seed*rep*b*i)
      #   resampled_precondition_prog <- resample_particle_y_samples(particle_set = precondition_prog_p2$particles[[1]],
      #                                                              multivariate = FALSE,
      #                                                              resampling_method = 'resid',
      #                                                              seed = seed*rep*b*i)$y_samples
      #   precondition_prog_results[[b]][[i]][[rep]] <- list('time' = sum(sapply(1:(1/beta_choices[b]), function(c) precondition_p1[[c]]$time))+sum(unlist(precondition_prog_p2$time)),
      #                                                      'ESS' = list(lapply(1:(1/beta_choices[b]), function(c) precondition_p1[[c]]$ESS), precondition_prog_p2$ESS),
      #                                                      'CESS' = list(lapply(1:(1/beta_choices[b]), function(c) precondition_p1[[c]]$CESS), precondition_prog_p2$CESS),
      #                                                      'IAD_bw' = integrated_abs_distance_uniGaussian(fusion_post = resampled_precondition_prog,
      #                                                                                                     mean = target_param['mu'],
      #                                                                                                     sd = target_param['sd'],
      #                                                                                                     beta = 1,
      #                                                                                                     bw = bw),
      #                                                      'IAD_opt_bw' = integrated_abs_distance_uniGaussian(fusion_post = resampled_precondition_prog,
      #                                                                                                         mean = target_param['mu'],
      #                                                                                                         sd = target_param['sd'],
      #                                                                                                         beta = 1,
      #                                                                                                         bw = opt_bw))
      # }
      # curve(dnorm(x, -mu_choices[i], sd), -(mu_choices[i]+5), mu_choices[i]+5,
      #       ylim = c(0,1), ylab = 'pdf', main = paste('Difference in means =', 2*mu_choices[i], '|| rep:', rep))
      # curve(dnorm(x, mu_choices[i], sd), add = T)
      # curve(dnorm_tempered(x, -mu_choices[i], sd, beta_choices[b]), add = T)
      # curve(dnorm_tempered(x, mu_choices[i], sd, beta_choices[b]), add = T)
      # curve(dnorm(x, 0, sqrt(1/2)), add = T, lty = 2, lwd = 3)
      # lines(density(resampled_precondition_bal, bw = bw), col = 'green')
      # if (beta_choices[b]!=1/2) {
      #   lines(density(resampled_precondition_prog, bw = bw), col = 'blue')  
      # }
      # for (c in 1:(1/beta_choices[b])) {
      #   lines(density(resample_particle_y_samples(particle_set = precondition_p1[[c]]$particles,
      #                                             multivariate = FALSE,
      #                                             resampling_method = 'resid',
      #                                             seed = seed)$y_samples),
      #         col = 'yellow')
      # }
      print('saving progress')
      save.image('separate_modes_with_tempering_varying_beta.RData')
    }
  }
}

######################################## IAD

# beta = 1/16
plot(x = log(2*mu_choices, 2), y = sapply(1:length(mu_choices), function(i) mean(sapply(1:number_of_replicates, function(j) precondition_bal_results[[4]][[i]][[j]]$IAD_bw))),
     ylim = c(0, 1), xlab = '', ylab = '', col = 'black', lty = 3, lwd = 3, pch = 5, xaxt = "n", yaxt = 'n', type = 'b')
axis(1, at = log(2*mu_choices, 2), labels = log(2*mu_choices, 2),  font = 2, cex = 1.5)
mtext('log(Difference in sub-posterior means, 2)', 1, 2.75, font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.2), labels=c(seq(0, 0.8, 0.2), "1.0"), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5)
mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
# beta = 1/8
lines(x = log(2*mu_choices, 2), y = sapply(1:length(mu_choices), function(i) mean(sapply(1:number_of_replicates, function(j) precondition_bal_results[[3]][[i]][[j]]$IAD_bw))),
      col = 'black', lty = 2, lwd = 3, pch = 4, type = 'b')
# beta = 1/4
lines(x = log(2*mu_choices, 2), y = sapply(1:length(mu_choices), function(i) mean(sapply(1:number_of_replicates, function(j) precondition_bal_results[[2]][[i]][[j]]$IAD_bw))),
      col = 'black', lty = 4, lwd = 3, pch = 3, type = 'b')
# beta = 1/2
lines(x = log(2*mu_choices[1:7], 2), y = c(sapply(1:length(mu_choices[1:6]), function(i) mean(sapply(1:number_of_replicates, function(j) precondition_bal_results[[1]][[i]][[j]]$IAD_bw))), 0.3200),
      col = 'black', lty = 5, lwd = 3, pch = 2, type = 'b')
# beta = 1
lines(x = log(2*mu_choices[1:6], 2), y = sapply(1:length(mu_choices[1:6]), function(i) mean(sapply(1:number_of_replicates, function(j) precondition_results[[i]][[j]]$IAD_bw))),
      col = 'black', lwd = 3, pch = 20, type = 'b')
legend(x = 1, y = 1,
       legend = c(bquote(bold('direct combination')),
                  as.expression(bquote(bold(paste(beta, ' = 1/2')))),
                  as.expression(bquote(bold(paste(beta, ' = 1/4')))),
                  as.expression(bquote(bold(paste(beta, ' = 1/8')))),
                  as.expression(bquote(bold(paste(beta, ' = 1/16'))))),
       lty = c(1, 5, 4, 2, 3),
       lwd = rep(3, 5),
       col = rep('black', 5),
       pch = c(20, 2, 3, 4, 5),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

######################################## time

# beta = 1/16
plot(x = log(2*mu_choices, 2),
     y = log(sapply(1:length(mu_choices), function(i) mean(sapply(1:number_of_replicates, function(j) precondition_bal_results[[4]][[i]][[j]]$time))), 2),
     ylim = c(0, 12), xlab = '', ylab = '', col = 'black', lty = 3, lwd = 3, pch = 5, xaxt = "n", yaxt = 'n', type = 'b')
axis(1, at = log(2*mu_choices, 2), labels = log(2*mu_choices, 2), font = 2, cex = 1.5)
mtext('log(Difference in sub-posterior means, 2)', 1, 2.75, font = 2, cex = 1.5)
axis(2, at=seq(0, 12, 2), labels=seq(0, 12, 2), font = 2, cex = 1.5)
axis(2, at=0:12, labels=rep("", 13), lwd.ticks = 0.5)
mtext('log(Time elapsed in seconds, 2)', 2, 2.75, font = 2, cex = 1.5)
# beta = 1/8
lines(x = log(2*mu_choices, 2),
      y = log(sapply(1:length(mu_choices), function(i) mean(sapply(1:number_of_replicates, function(j) precondition_bal_results[[3]][[i]][[j]]$time))), 2),
      col = 'black', lty = 2, lwd = 3, pch = 4, type = 'b')
# beta = 1/4
lines(x = log(2*mu_choices, 2),
      y = log(sapply(1:length(mu_choices), function(i) mean(sapply(1:number_of_replicates, function(j) precondition_bal_results[[2]][[i]][[j]]$time))), 2),
      col = 'black', lty = 4, lwd = 3, pch = 3, type = 'b')
# beta = 1/2
lines(x = log(2*mu_choices[1:7], 2),
      y = log(sapply(1:length(mu_choices[1:7]), function(i) mean(sapply(1:number_of_replicates, function(j) precondition_bal_results[[1]][[i]][[j]]$time))), 2),
      col = 'black', lty = 5, lwd = 3, pch = 2, type = 'b')
# beta = 1
lines(x = log(2*mu_choices[1:6], 2),
      y = log(sapply(1:length(mu_choices[1:6]), function(i) mean(sapply(1:number_of_replicates, function(j) precondition_results[[i]][[j]]$time))), 2),
      col = 'black', lwd = 3, pch = 20, type = 'b')
legend(x = 1, y = 12,
       legend = c(bquote(bold('direct combination')),
                  as.expression(bquote(bold(paste(beta, ' = 1/2')))),
                  as.expression(bquote(bold(paste(beta, ' = 1/4')))),
                  as.expression(bquote(bold(paste(beta, ' = 1/8')))),
                  as.expression(bquote(bold(paste(beta, ' = 1/16'))))),
       lty = c(1, 5, 4, 2, 3),
       lwd = rep(3, 5),
       col = rep('black', 5),
       pch = c(20, 2, 3, 4, 5),
       cex = 1.25,
       text.font = 2,
       bty = 'n')
