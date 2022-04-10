library(DCFusion)
library(HMCBLR)

##### Initialise example #####
seed <- 2021
set.seed(seed)
nsamples <- c(10000, 30000, 50000, 100000)
ndata <- 5000
n_cores <- parallel::detectCores()
ESS_threshold <- 0.5
resampling_method <- 'resid'
dim <- c(5,6,7,8,9,10,12,14)
time_choices <- c(0.1, 0.2, 0.25, 0.4, 0.5)
true_beta <- list()
frequencies <- list()
simulated_data <- list()
full_posterior <- list()
data_split_16 <- list()
sub_posteriors_16 <- list()
results <- list()

for (d in 1:length(dim)) {
  print(paste('$$$$$$$$$$ d:', d))
  print(paste('dim:', dim[d]))
  set.seed(seed)
  prior_means <- rep(0, dim[d])
  prior_variances <- rep(1, dim[d])
  true_beta[[d]] <- round(rnorm(dim[d], 0, 1), 2)
  frequencies[[d]] <- round(runif(dim[d]-1, 0.01, 0.99), 2)
  results[[d]] <- list()
  
  # simulate data set
  simulated_data[[d]] <- simulate_LR_data(N = ndata,
                                          alpha = true_beta[[d]][1],
                                          frequencies = frequencies[[d]],
                                          coefficients = true_beta[[d]][2:length(true_beta[[d]])],
                                          seed = seed)
  
  # check activity of the parameters
  check_activity(simulated_data[[d]])
  
  ##### Sampling from full posterior #####
  
  full_data_count <- unique_row_count(y = simulated_data[[d]][,1],
                                      X = cbind('intercept' = rep(1, ndata), simulated_data[[d]][,2:ncol(simulated_data[[d]])]))$full_data_count
  full_posterior[[d]] <- hmc_sample_BLR(full_data_count = full_data_count,
                                        C = 1,
                                        prior_means = prior_means,
                                        prior_variances = prior_variances,
                                        iterations = 110000,
                                        warmup = 10000,
                                        chains = 1,
                                        seed = seed,
                                        output = T)
  
  ##### Sampling from sub-posterior C=16 #####
  
  data_split_16[[d]] <- split_data(simulated_data[[d]],
                                   y_col_index = 1,
                                   X_col_index = 2:ncol(simulated_data[[d]]),
                                   C = 16,
                                   as_dataframe = F)
  sub_posteriors_16[[d]] <- hmc_base_sampler_BLR(nsamples = 100000,
                                                 data_split = data_split_16[[d]],
                                                 C = 16,
                                                 prior_means = prior_means,
                                                 prior_variances = prior_variances,
                                                 warmup = 10000,
                                                 seed = seed,
                                                 output = T)
  
  ##### Applying Fusion #####
  
  ##### NB (Hypercube Centre) #####
  print('NB Fusion (hypercube centre)')
  for (i in 1:length(time_choices)) {
    print(paste('##### i:', i, '|| time:', time_choices[i]))
    results[[d]][[i]] <- list()
    for (j in 1:length(nsamples)) {
      print(paste('!!!!! j:', j, '|| n:', nsamples[j]))
      DC_MCF <- bal_binary_fusion_SMC_BLR(N_schedule = rep(nsamples[j], 4),
                                          m_schedule = rep(2, 4),
                                          time_schedule = rep(time_choices[i], 4),
                                          base_samples = sub_posteriors_16[[d]],
                                          L = 5,
                                          dim = dim[d],
                                          data_split = data_split_16[[d]],
                                          prior_means = prior_means,
                                          prior_variances = prior_variances,
                                          C = 16,
                                          precondition = TRUE,
                                          resampling_method = resampling_method,
                                          ESS_threshold = ESS_threshold,
                                          cv_location = 'hypercube_centre',
                                          diffusion_estimator = 'NB',
                                          seed = seed,
                                          n_cores = n_cores,
                                          print_progress_iters = 500)
      samples <- resample_particle_y_samples(particle_set = DC_MCF$particles[[1]],
                                             multivariate = TRUE,
                                             resampling_method = resampling_method,
                                             seed = seed)$y_samples
      results[[d]][[i]][[j]] <- list('IAD' = integrated_abs_distance(full_post = full_posterior[[d]],
                                                                     fusion_post = samples),
                                     'time' = sum(unlist(DC_MCF$time)),
                                     'samples' = samples,
                                     'CESS' = DC_MCF$CESS,
                                     'ESS' = DC_MCF$ESS,
                                     'resampled' = DC_MCF$resampled)
      print(paste('IAD:', results[[d]][[i]][[j]]$IAD))
      print(paste('time:', results[[d]][[i]][[j]]$time))
      print(paste('log(time):', log(results[[d]][[i]][[j]]$time)))
      print('save_progress')
      save.image('SD16_dim_study_varying.RData')
    }
  }
}

dim <- c(5,6,7,8,9,10,12,14)
plot(x = dim,
     y = sapply(1:length(dim), function(d) results[[d]][[2]][[length(results[[d]][[2]])]]$IAD),
     ylim = c(0, 0.2),
     xlab = '',
     ylab = '',
     xaxt = 'n',
     yaxt = 'n', lwd = 3, pch = 1, type = 'b')
mtext('Dimension', 1, 2.75, font = 2, cex = 1.5)
mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=dim, labels=dim, font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=c("0.0", seq(0.1, 0.9, 0.1), "1.0"), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.05), labels=rep("", 21), lwd.ticks = 0.5)
lines(x = c(5, 10, 12, 14), y = c(0.021, sapply(2:length(c(5, 10, 12, 14)), function(d) {
  integrated_abs_distance(full_posterior[[d]], balanced_C16[[d]]$adaptive$particles$y_samples)})),
  lty = 2, lwd = 3, pch = 4, type = 'b')
lines(x = c(5, 10, 12, 14), y = c(0.02, sapply(2:length(c(5, 10, 12, 14)), function(d) {
  integrated_abs_distance(full_posterior[[d]], balanced_C16[[d]]$reg$particles$y_samples)})),
  lty = 3, lwd = 3, pch = 20, type = 'b')
legend(x = 5, y = 0.2,
       legend = c('D&C-MCF', 'D&C-GBF (reg)', 'D&C-GBF (adaptive)'),
       lwd = c(3, 3, 3),
       lty = c(1, 2, 3),
       col = 'black',
       cex = 1.25,
       text.font = 2,
       bty = 'n')

plot(x = dim, 
     y = log(sapply(1:length(dim), function(d) results[[d]][[2]][[length(results[[d]][[2]])]]$time)),
     ylim = c(-2, 14), xlab = '', ylab = '', yaxt = 'n', xaxt = 'n', lwd = 3, pch = 1, type = 'b')
mtext('Dimension', 1, 2.75, font = 2, cex = 1.5)
mtext('log(Time elapsed in seconds)', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=c(seq(0, 0.9, 0.1), 0.95), labels=c("0.0", c(seq(0.1, 0.9, 0.1), 0.95)), font = 2, cex = 1.5)
axis(1, at=dim, labels = dim, font = 2, cex = 1.5)
axis(2, at=seq(-2, 14, 2), labels = seq(-2, 14, 2), font = 2, cex = 1.5)
axis(2, at=seq(-2, 14, 1), labels=rep("", 17), lwd.ticks = 0.5)
lines(x = c(5, 10, 12, 14), y = c(5.9, log(sapply(2:length(c(5, 10, 12, 14)), function(d) sum(unlist(balanced_C16[[d]]$adaptive$time))))),
      lty = 2, lwd = 3, pch = 4, type = 'b')
lines(x = c(5, 10, 12, 14), y = c(6.1, log(sapply(2:length(c(5, 10, 12, 14)), function(d) sum(unlist(balanced_C16[[d]]$reg$time))))),
      lty = 3, lwd = 3, pch = 20, type = 'b')
legend(x = 5, y = 14,
       legend = c('D&C-MCF', 'D&C-GBF (reg)', 'D&C-GBF (adaptive)'),
       lwd = c(3, 3, 3),
       lty = c(1, 2, 3),
       col = 'black',
       cex = 1.25,
       text.font = 2,
       bty = 'n')

save.image('SD16_dim_study_varying.RData')
