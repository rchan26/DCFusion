library(hierarchicalFusion)

seed <- 1994
set.seed(seed)
number_of_replicates <- 10
mu_choices <- c(1, 2, 4, 8, 16, 32, 64, 128)
beta_choices <- c(1/2, 1/4, 1/8)
indices <- list(1:7, 1:8, 1:8)
sd <- 1
nsamples <- 10000
time_choice <- 1
target_param <- c('mu' = 0, 'sd' = sqrt(1/2))
opt_bw <- ((4/(3*nsamples))^(1/5))*sqrt(1/2)
bw <- 0.1
precondition_hier_results <- list()
precondition_prog_results <- list()

for (b in 1:length(beta_choices)) {
  print(paste('b:', b))
  print(paste('beta:', beta_choices[b]))
  precondition_hier_results[[b]] <- list()
  precondition_prog_results[[b]] <- list()
  for (i in indices[[b]]) {
    print(paste('i:', i))
    print(paste('mu:', mu_choices[i]))
    print(paste('difference in means:', 2*mu_choices[i]))
    precondition_hier_results[[b]][[i]] <- list()
    precondition_prog_results[[b]][[i]] <- list()
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
      # part 2a: obtaining particle approximation for f using hierarchical fusion (balanced binary)
      print('fusion for f with balanced binary tree')
      precondition_hier_p2 <- hierarchical_fusion_SMC_uniGaussian(N_schedule = rep(nsamples, log((1/beta_choices[b]), 2)),
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
      resampled_precondition_hier <- resample_particle_y_samples(particle_set = precondition_hier_p2$particles[[1]],
                                                                 multivariate = FALSE,
                                                                 resampling_method = 'resid',
                                                                 seed = seed*rep*b*i)$y_samples
      precondition_hier_results[[b]][[i]][[rep]] <- list('time' = sum(sapply(1:(1/beta_choices[b]), function(c) precondition_p1[[c]]$time))+sum(unlist(precondition_hier_p2$time)),
                                                         'ESS' = list(lapply(1:(1/beta_choices[b]), function(c) precondition_p1[[c]]$ESS), precondition_hier_p2$ESS),
                                                         'CESS' = list(lapply(1:(1/beta_choices[b]), function(c) precondition_p1[[c]]$CESS), precondition_hier_p2$CESS),
                                                         'IAD_bw' = integrated_abs_distance_uniGaussian(fusion_post = resampled_precondition_hier,
                                                                                                        mean = target_param['mu'],
                                                                                                        sd = target_param['sd'],
                                                                                                        beta = 1,
                                                                                                        bw = bw),
                                                         'IAD_opt_bw' = integrated_abs_distance_uniGaussian(fusion_post = resampled_precondition_hier,
                                                                                                            mean = target_param['mu'],
                                                                                                            sd = target_param['sd'],
                                                                                                            beta = 1,
                                                                                                            bw = opt_bw))
      if (beta_choices[b]==1/2) {
        # the hierarchical and progressive trees are exactly the same
        precondition_prog_results[[b]][[i]][[rep]] <- precondition_hier_results[[b]][[i]][[rep]]
      } else {
        # part 2b: obtaining particle approximation for f using hierarchical fusion (progressive)
        print('fusion for f with progressive tree')
        precondition_prog_p2 <- progressive_fusion_SMC_uniGaussian(N_schedule = rep(nsamples, (1/beta_choices[b])-1),
                                                                   time_schedule = rep(time_choice, (1/beta_choices[b])-1),
                                                                   base_samples = lapply(1:(1/beta_choices[b]), function(c) precondition_p1[[c]]$particles),
                                                                   mean = target_param['mu'],
                                                                   sd = target_param['sd'],
                                                                   start_beta = beta_choices[b],
                                                                   precondition = lapply(1:(1/beta_choices[b]), function(c) precondition_p1[[c]]$precondition_values[[1]]),
                                                                   resampling_method = 'resid',
                                                                   ESS_threshold = 0.5,
                                                                   diffusion_estimator = 'NB',
                                                                   seed = seed*rep*b*i)
        resampled_precondition_prog <- resample_particle_y_samples(particle_set = precondition_prog_p2$particles[[1]],
                                                                   multivariate = FALSE,
                                                                   resampling_method = 'resid',
                                                                   seed = seed*rep*b*i)$y_samples
        precondition_prog_results[[b]][[i]][[rep]] <- list('time' = sum(sapply(1:(1/beta_choices[b]), function(c) precondition_p1[[c]]$time))+sum(unlist(precondition_prog_p2$time)),
                                                           'ESS' = list(lapply(1:(1/beta_choices[b]), function(c) precondition_p1[[c]]$ESS), precondition_prog_p2$ESS),
                                                           'CESS' = list(lapply(1:(1/beta_choices[b]), function(c) precondition_p1[[c]]$CESS), precondition_prog_p2$CESS),
                                                           'IAD_bw' = integrated_abs_distance_uniGaussian(fusion_post = resampled_precondition_prog,
                                                                                                          mean = target_param['mu'],
                                                                                                          sd = target_param['sd'],
                                                                                                          beta = 1,
                                                                                                          bw = bw),
                                                           'IAD_opt_bw' = integrated_abs_distance_uniGaussian(fusion_post = resampled_precondition_prog,
                                                                                                              mean = target_param['mu'],
                                                                                                              sd = target_param['sd'],
                                                                                                              beta = 1,
                                                                                                              bw = opt_bw))
      }
      
      save.image('separate_modes_with_tempering_varying_beta.RData')
      
      curve(dnorm(x, -mu_choices[i], sd), -(mu_choices[i]+5), mu_choices[i]+5,
            ylim = c(0,1), ylab = 'pdf', main = paste('Difference in means =', 2*mu_choices[i], '|| rep:', rep))
      curve(dnorm(x, mu_choices[i], sd), add = T)
      curve(dnorm_tempered(x, -mu_choices[i], sd, beta_choices[b]), add = T)
      curve(dnorm_tempered(x, mu_choices[i], sd, beta_choices[b]), add = T)
      curve(dnorm(x, 0, sqrt(1/2)), add = T, lty = 2, lwd = 3)
      lines(density(resampled_precondition_hier, bw = bw), col = 'green')
      if (beta_choices[b]!=1/2) {
        lines(density(resampled_precondition_prog, bw = bw), col = 'blue')  
      }
      for (c in 1:(1/beta_choices[b])) {
        lines(density(resample_particle_y_samples(particle_set = precondition_p1[[c]]$particles,
                                                  multivariate = FALSE,
                                                  resampling_method = 'resid',
                                                  seed = seed)$y_samples),
              col = 'yellow')
      }
    }
  }
}
curve(dnorm(x, -mu_choices[i], 1), -(mu_choices[i]+5), mu_choices[i]+5,
      ylim = c(0,1), ylab = 'pdf', main = paste('Difference in means =', 2*mu_choices[i]))
curve(dnorm(x, mu_choices[i], 1), add = T)
curve(dnorm_tempered(x, -mu_choices[i], 1, beta = 1/8), add = T)
curve(dnorm_tempered(x, mu_choices[i], 1, beta = 1/8), add = T)

######################################## running time
# also run seperate_modes_smc.R to get standard_results to create the next plots

plot(x = 2*mu_choices[1:7], y = sapply(1:length(mu_choices[1:7]), function(i) mean(sapply(1:number_of_replicates, function(j) precondition_hier_results[[1]][[i]][[j]]$time))),
     ylim = c(0, 600), ylab = 'log(Time elapsed in seconds)', xlab = 'Difference between sub-posterior means',
     col = 'black', lty = 3, lwd = 3, xaxt = "n", type = 'l')
lines(x = 2*mu_choices[1:7], y = sapply(1:length(mu_choices[1:7]), function(i) mean(sapply(1:number_of_replicates, function(j) precondition_prog_results[[1]][[i]][[j]]$time))),
      col = 'black', lty = 2, lwd = 3)
lines(x = 2*mu_choices, y = sapply(1:length(mu_choices), function(i) mean(sapply(1:number_of_replicates, function(j) standard_results[[i]][[j]]$time))),
      col = 'black', lwd = 3)
axis(1, at=0:(2*mu_choices[length(mu_choices)]), labels=0:(2*mu_choices[length(mu_choices)]))
# axis(1, at=2*mu_choices, labels=2*mu_choices)
legend(x = 0, y = 6, 
       legend = c('direct combination', 'balanced', 'progressive'),
       lty = c(1, 3, 2), 
       lwd = c(3, 3, 3),
       # pch = c(1, 0, 2), 
       cex = 1.25,
       text.font = 2,
       bty = 'n')

#################### log running time

plot(x = 2*mu_choices[1:7], y = sapply(1:length(mu_choices[1:7]), function(i) log(mean(sapply(1:number_of_replicates, function(j) precondition_hier_results[[1]][[i]][[j]]$time)))),
     ylim = c(1, 8), ylab = 'log(Time elapsed in seconds)', xlab = 'Difference between sub-posterior means',
     col = 'black', lty = 3, lwd = 3, xaxt = "n", type = 'l')
lines(x = 2*mu_choices[1:7], y = sapply(1:length(mu_choices[1:7]), function(i) log(mean(sapply(1:number_of_replicates, function(j) precondition_prog_results[[1]][[i]][[j]]$time)))),
      col = 'black', lty = 2, lwd = 3)
lines(x = 2*mu_choices, y = sapply(1:length(mu_choices), function(i) log(mean(sapply(1:number_of_replicates, function(j) standard_results[[i]][[j]]$time)))),
      col = 'black', lwd = 3)
axis(1, at=0:(2*mu_choices[length(mu_choices)]), labels=0:(2*mu_choices[length(mu_choices)]))
legend(x = 0, y = 6, 
       legend = c('direct combination', 'balanced', 'progressive'),
       lty = c(1, 3, 2), 
       lwd = c(3, 3, 3),
       # pch = c(1, 0, 2), 
       cex = 1.25,
       text.font = 2,
       bty = 'n')

######################################## IAD

plot(x = 2*mu_choices[1:7], y = sapply(1:length(mu_choices[1:7]), function(i) mean(sapply(1:number_of_replicates, function(j) precondition_hier_results[[1]][[i]][[j]]$IAD))),
     ylim = c(0, 0.6), ylab = 'IAD', xlab = 'Difference between sub-posterior means',
     col = 'black', lty = 3, lwd = 3, xaxt = "n", type = 'l')
# lines(x = 2*mu_choices, y = sapply(1:length(mu_choices), function(i) precondition_hier_results[[i]]$IAD), 
#       col = 'black', lty = 2, lwd = 3)
# points(x = 2*mu_choices, y = sapply(1:length(mu_choices), function(i) precondition_prog_results[[i]]$IAD),
#        col = 'black', lwd = 3)
lines(x = 2*mu_choices[1:7], y = sapply(1:length(mu_choices[1:7]), function(i) mean(sapply(1:number_of_replicates, function(j) precondition_prog_results[[1]][[i]][[j]]$IAD))),
      col = 'black', lty = 2, lwd = 3)
# points(x = 2*mu_choices, y = sapply(1:length(mu_choices), function(i) standard_results[[i]]$IAD),
#        col = 'black', lwd = 3)
lines(x = 2*mu_choices, y = sapply(1:length(mu_choices), function(i) mean(sapply(1:number_of_replicates, function(j) standard_results[[i]][[j]]$IAD))),
      col = 'black', lwd = 3)
# axis(1, at=0:(2*mu_choices[length(mu_choices)]), labels=0:(2*mu_choices[length(mu_choices)]))
axis(1, at=2*mu_choices, labels=2*mu_choices)
legend(x = 0, y = 0.6, 
       legend = c('direct combination', 'balanced', 'progressive'),
       lty = c(1, 3, 2), 
       lwd = c(3, 3, 3),
       # pch = c(1, 0, 2), 
       cex = 1.25,
       text.font = 2,
       bty = 'n')

######################################## log(IAD)

plot(x = 2*mu_choices, y = log(sapply(1:length(mu_choices), function(i) mean(sapply(1:number_of_replicates, function(j) precondition_hier_results[[i]][[j]]$IAD)))),
     ylim = c(-5, 0), ylab = 'log(IAD)', xlab = 'Difference between sub-posterior means',
     col = 'black', lty = 3, lwd = 3, xaxt = "n", type = 'l')
# lines(x = 2*mu_choices, y = sapply(1:length(mu_choices), function(i) precondition_hier_results[[i]]$IAD), 
#       col = 'black', lty = 2, lwd = 3)
# points(x = 2*mu_choices, y = sapply(1:length(mu_choices), function(i) precondition_prog_results[[i]]$IAD),
#        col = 'black', lwd = 3)
lines(x = 2*mu_choices, y = log(sapply(1:length(mu_choices), function(i) mean(sapply(1:number_of_replicates, function(j) precondition_prog_results[[i]][[j]]$IAD)))),
      col = 'black', lty = 2, lwd = 3)
# points(x = 2*mu_choices, y = sapply(1:length(mu_choices), function(i) standard_results[[i]]$IAD),
#        col = 'black', lwd = 3)
lines(x = 2*mu_choices, y = log(sapply(1:length(mu_choices), function(i) mean(sapply(1:number_of_replicates, function(j) standard_results[[i]][[j]]$IAD)))),
      col = 'black', lwd = 3)
axis(1, at=0:(2*mu_choices[length(mu_choices)]), labels=0:(2*mu_choices[length(mu_choices)]))
legend(x = 0, y = 0, 
       legend = c('direct combination', 'balanced', 'progressive'),
       lty = c(1, 3, 2), 
       lwd = c(3, 3, 3),
       # pch = c(1, 0, 2), 
       cex = 1.25,
       text.font = 2,
       bty = 'n')
