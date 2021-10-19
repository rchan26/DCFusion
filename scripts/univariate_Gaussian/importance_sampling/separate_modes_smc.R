library(DCFusion)

seed <- 1994
set.seed(seed)
number_of_replicates <- 50
mu_choices <- c(1, 2, 4, 8, 16, 32)
sd <- 1
nsamples <- 10000
time_choice <- 1
target_param <- c('mu' = 0, 'sd' = sqrt(1/2))
opt_bw <- ((4/(3*nsamples))^(1/5))*sqrt(1/2)
bw <- 0.1
standard_results <- list()
precondition_results <- list()

for (i in 1:length(mu_choices)) {
  print(paste('i:', i))
  print(paste('mu:', mu_choices[i]))
  print(paste('difference in means:', 2*mu_choices[i]))
  standard_results[[i]] <- list()
  precondition_results[[i]] <- list()
  for (rep in 1:number_of_replicates) {
    print(paste('rep:', rep))
    set.seed(seed*rep*i)
    input_samples <- list(rnorm(n = nsamples, mean = -mu_choices[i], sd = sd),
                          rnorm(n = nsamples, mean = mu_choices[i], sd = sd))
    input_particles <- initialise_particle_sets(input_samples, multivariate = FALSE)
    # fusion without preconditioning
    print('fusion without preconditioning')
    standard <- parallel_fusion_SMC_uniGaussian(particles_to_fuse = input_particles, 
                                                N = nsamples,
                                                m = 2, 
                                                time = time_choice,
                                                means = c(-mu_choices[i], mu_choices[i]),
                                                sds = c(sd, sd), 
                                                betas = rep(1, 2),
                                                precondition_values = rep(1, 2),
                                                ESS_threshold = 0.5,
                                                resampling_method = 'resid',
                                                seed = seed*rep*i)
    resampled_standard <- resample_particle_y_samples(particle_set = standard$particles,
                                                      multivariate = FALSE,
                                                      resampling_method = 'resid',
                                                      seed = seed*rep*i)$y_samples
    standard_results[[i]][[rep]] <- list('time' = standard$time,
                                         'ESS' = standard$ESS,
                                         'CESS' = standard$CESS,
                                         'IAD_bw' = integrated_abs_distance_uniGaussian(fusion_post = resampled_standard,
                                                                                        mean = target_param['mu'],
                                                                                        sd = target_param['sd'],
                                                                                        beta = 1,
                                                                                        bw = bw),
                                         'IAD_opt_bw' = integrated_abs_distance_uniGaussian(fusion_post = resampled_standard,
                                                                                            mean = target_param['mu'],
                                                                                            sd = target_param['sd'],
                                                                                            beta = 1,
                                                                                            bw = opt_bw))
    # fusion with preconditioning
    print('fusion with preconditioning')
    precondition <- parallel_fusion_SMC_uniGaussian(particles_to_fuse = input_particles, 
                                                    N = nsamples,
                                                    m = 2, 
                                                    time = time_choice,
                                                    means = c(-mu_choices[i], mu_choices[i]),
                                                    sds = c(sd, sd), 
                                                    betas = rep(1, 2),
                                                    precondition_values = sapply(input_samples, var),
                                                    ESS_threshold = 0.5,
                                                    resampling_method = 'resid',
                                                    seed = seed*rep*i)
    resampled_precondition <- resample_particle_y_samples(particle_set = precondition$particles,
                                                          multivariate = FALSE,
                                                          resampling_method = 'resid',
                                                          seed = seed*rep*i)$y_samples
    precondition_results[[i]][[rep]] <- list('time' = precondition$time,
                                             'ESS' = precondition$ESS,
                                             'CESS' = precondition$CESS,
                                             'IAD_bw' = integrated_abs_distance_uniGaussian(fusion_post = resampled_precondition,
                                                                                            mean = target_param['mu'],
                                                                                            sd = target_param['sd'],
                                                                                            beta = 1,
                                                                                            bw = bw),
                                             'IAD_opt_bw' = integrated_abs_distance_uniGaussian(fusion_post = resampled_precondition,
                                                                                                mean = target_param['mu'],
                                                                                                sd = target_param['sd'],
                                                                                                beta = 1,
                                                                                                bw = opt_bw))
    print('saving progress')
    save.image('separate_modes_smc_uniG.RData')
    # curve(dnorm(x, -mu_choices[i], 1), -(mu_choices[i]+5), mu_choices[i]+5,
    #       ylim = c(0,1), ylab = 'pdf', main = paste('Difference in means =', 2*mu_choices[i], '|| rep:', rep))
    # curve(dnorm(x, mu_choices[i], 1), add = T)
    # curve(dnorm(x, 0, sqrt(1/2)), add = T, lty = 2, lwd = 3)
    # lines(density(resampled_precondition, bw = opt_bw), col = 'green')
    # lines(density(resampled_standard, bw = opt_bw), col = 'blue')
  }
}

######################################## running time

plot(x = 2*mu_choices, y = sapply(1:length(mu_choices), function(i) mean(sapply(1:number_of_replicates, function(j) precondition_results[[i]][[j]]$time))),
     ylim = c(0, 40), ylab = 'Running time in seconds', xlab = 'Difference between sub-posterior means',
     col = 'black', lwd = 3, xaxt = "n")
lines(x = 2*mu_choices, y = sapply(1:length(mu_choices), function(i) mean(sapply(1:number_of_replicates, function(j) precondition_results[[i]][[j]]$time))),
      col = 'black', lwd = 3)
axis(1, at=2*mu_choices, labels=2*mu_choices)

#################### log running time

plot(x = 2*mu_choices, y = log(sapply(1:length(mu_choices), function(i) mean(sapply(1:number_of_replicates, function(j) precondition_results[[i]][[j]]$time)))),
     ylim = c(1, 4), ylab = 'log(Time elapsed in seconds)', xlab = 'Difference between sub-posterior means',
     col = 'black', lwd = 3, xaxt = "n")
lines(x = 2*mu_choices, y = log(sapply(1:length(mu_choices), function(i) mean(sapply(1:number_of_replicates, function(j) precondition_results[[i]][[j]]$time)))),
      col = 'black', lwd = 3)
axis(1, at=2*mu_choices, labels=2*mu_choices)

######################################## IAD (with preconditioning)

plot(x = 2*mu_choices, y = sapply(1:length(mu_choices), function(i) mean(sapply(1:number_of_replicates, function(j) precondition_results[[i]][[j]]$IAD_bw))),
     ylim = c(0, 1), ylab = 'IAD', xlab = 'Difference between sub-posterior means',
     col = 'black', lwd = 3, xaxt = "n")
lines(x = 2*mu_choices, y = sapply(1:length(mu_choices), function(i) mean(sapply(1:number_of_replicates, function(j) precondition_results[[i]][[j]]$IAD_bw))),
      col = 'black', lwd = 3)
points(x = 2*mu_choices, y = sapply(1:length(mu_choices), function(i) mean(sapply(1:number_of_replicates, function(j) precondition_results[[i]][[j]]$IAD_opt_bw))),
       lwd = 3)
lines(x = 2*mu_choices, y = sapply(1:length(mu_choices), function(i) mean(sapply(1:number_of_replicates, function(j) precondition_results[[i]][[j]]$IAD_opt_bw))),
      col = 'black', lwd = 3, lty = 2)
axis(1, at=2*mu_choices, labels=2*mu_choices)
