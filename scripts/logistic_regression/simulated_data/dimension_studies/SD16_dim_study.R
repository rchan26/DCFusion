library(DCFusion)
library(HMCBLR)

##### Initialise example #####
seed <- 2021
set.seed(seed)
nsamples <- 20000
ndata <- 5000
time_choice <- 0.5
n_cores <- parallel::detectCores()
ESS_threshold <- 0.5
resampling_method <- 'resid'
dim <- c(5,6,7,8,9,10,15)
true_beta <- list()
frequencies <- list()
simulated_data <- list()
full_posterior <- list()
data_split_16 <- list()
sub_posteriors_16 <- list()
consensus_mat_16 <- list()
consensus_sca_16 <- list()
neiswanger_true_16 <- list()
neiswanger_false_16 <- list()
weierstrass_importance_16 <- list()
weierstrass_rejection_16 <- list()
Poisson_hc_16 <- list()
NB_hc_16 <- list()

for (d in 1:length(dim)) {
  print(paste('$$$$$$$$$$ d:', d))
  print(paste('dim:', dim[d]))
  set.seed(seed)
  prior_means <- rep(0, dim[d])
  prior_variances <- rep(1, dim[d])
  true_beta[[d]] <- round(rnorm(dim[d], 0, 1), 2)
  frequencies[[d]] <- round(runif(dim[d]-1, 0.01, 0.99), 2)
  
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
                                        iterations = nsamples + 10000,
                                        warmup = 10000,
                                        chains = 1,
                                        seed = seed,
                                        output = T)
  
  print(apply(full_posterior[[d]], 2, mean))
  print(true_beta[[d]])
  print(abs(apply(full_posterior[[d]], 2, mean)-true_beta[[d]]))
  
  ##### Sampling from sub-posterior C=16 #####
  
  data_split_16[[d]] <- split_data(simulated_data[[d]], y_col_index = 1, X_col_index = 2:ncol(simulated_data[[d]]), C = 16, as_dataframe = F)
  sub_posteriors_16[[d]] <- hmc_base_sampler_BLR(nsamples = nsamples,
                                                 data_split = data_split_16[[d]],
                                                 C = 16, 
                                                 prior_means = prior_means,
                                                 prior_variances = prior_variances,
                                                 warmup = 10000,
                                                 seed = seed,
                                                 output = T)
  
  ##### Applying other methodologies #####
  
  print('Applying other methodologies')
  consensus_mat_16[[d]] <- consensus_scott(S = 16, samples_to_combine = sub_posteriors_16[[d]], indep = F)
  consensus_sca_16[[d]] <- consensus_scott(S = 16, samples_to_combine = sub_posteriors_16[[d]], indep = T)
  neiswanger_true_16[[d]] <- neiswanger(S = 16,
                                        samples_to_combine = sub_posteriors_16[[d]],
                                        anneal = TRUE)
  neiswanger_false_16[[d]] <- neiswanger(S = 16,
                                         samples_to_combine = sub_posteriors_16[[d]],
                                         anneal = FALSE)
  weierstrass_importance_16[[d]] <- weierstrass(Samples = sub_posteriors_16[[d]],
                                                method = 'importance')
  weierstrass_rejection_16[[d]] <- weierstrass(Samples = sub_posteriors_16[[d]],
                                               method = 'reject')
  
  ##### Applying Fusion #####
  
  ##### Poisson (Hypercube Centre) #####
  print('Poisson Fusion (hypercube centre)')
  Poisson_hc_16[[d]] <- bal_binary_fusion_SMC_BLR(N_schedule = rep(nsamples, 4),
                                                  m_schedule = rep(2, 4),
                                                  time_schedule = rep(time_choice, 4),
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
                                                  diffusion_estimator = 'Poisson',
                                                  seed = seed,
                                                  n_cores = n_cores,
                                                  print_progress_iters = 10)
  Poisson_hc_16[[d]]$particles <- resample_particle_y_samples(particle_set = Poisson_hc_16[[d]]$particles[[1]],
                                                              multivariate = TRUE,
                                                              resampling_method = resampling_method,
                                                              seed = seed)
  Poisson_hc_16[[d]]$proposed_samples <- Poisson_hc_16[[d]]$proposed_samples[[1]]
  print(integrated_abs_distance(full_posterior[[d]], Poisson_hc_16[[d]]$particles$y_samples))
  
  print('save_progress')
  save.image('SD16_dim_study.RData')
  
  ##### NB (Hypercube Centre) #####
  print('NB Fusion (hypercube centre)')
  NB_hc_16[[d]] <- bal_binary_fusion_SMC_BLR(N_schedule = rep(nsamples, 4),
                                             m_schedule = rep(2, 4),
                                             time_schedule = rep(time_choice, 4),
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
  NB_hc_16[[d]]$particles <- resample_particle_y_samples(particle_set = NB_hc_16[[d]]$particles[[1]],
                                                         multivariate = TRUE,
                                                         resampling_method = resampling_method,
                                                         seed = seed)
  NB_hc_16[[d]]$proposed_samples <- NB_hc_16[[d]]$proposed_samples[[1]]
  print(integrated_abs_distance(full_posterior[[d]], NB_hc_16[[d]]$particles$y_samples))
  
  print('save_progress')
  save.image('SD16_dim_study.RData')
  
  ##### IAD #####
  
  print(integrated_abs_distance(full_posterior[[d]], Poisson_hc_16[[d]]$particles$y_samples))
  print(integrated_abs_distance(full_posterior[[d]], NB_hc_16[[d]]$particles$y_samples))
  print(integrated_abs_distance(full_posterior[[d]], consensus_mat_16[[d]]$samples))
  print(integrated_abs_distance(full_posterior[[d]], consensus_sca_16[[d]]$samples))
  print(integrated_abs_distance(full_posterior[[d]], neiswanger_true_16[[d]]$samples))
  print(integrated_abs_distance(full_posterior[[d]], neiswanger_false_16[[d]]$samples))
  print(integrated_abs_distance(full_posterior[[d]], weierstrass_importance_16[[d]]$samples))
  print(integrated_abs_distance(full_posterior[[d]], weierstrass_rejection_16[[d]]$samples))
}

dim <- 5:10
plot(x = dim,
     y = sapply(1:length(dim), function(d) {
       integrated_abs_distance(full_posterior[[d]], Poisson_hc_16[[d]]$particles$y_samples)}),
     ylim = c(0, 0.6),
     xlab = '',
     ylab = '',
     xaxt = 'n', lwd = 3, pch = 1, type = 'l')
mtext('Dimension', 1, 2.75, font = 2, cex = 1.5)
mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=dim, labels=dim, font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=c("0.0", seq(0.1, 0.9, 0.1), "1.0"), font = 2, cex = 1.5)
axis(2, at=seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5)
lines(x = dim, y = sapply(1:6, function(d) {
  integrated_abs_distance(full_posterior[[d]], NB_hc_16[[d]]$particles$y_samples)}),
      lty = 2, lwd = 3)
lines(x = dim, y = sapply(1:6, function(d) {
  integrated_abs_distance(full_posterior[[d]], consensus_mat_16[[d]]$samples)}),
      lty = 3, lwd = 3, col = 'red')
lines(x = dim, y = sapply(1:6, function(d) {
  integrated_abs_distance(full_posterior[[d]], neiswanger_false_16[[d]]$samples)}),
      lty = 4, lwd = 3, col = 'red')
lines(x = dim, y = sapply(1:6, function(d) {
  integrated_abs_distance(full_posterior[[d]], weierstrass_rejection_16[[d]]$samples)}),
      lty = 5, lwd = 3, col = 'red')
legend(x = 5, y = 0.6,
       legend = c('D&C-MCF (Poisson)', 'D&C-MCF (NB)', 'CMC', 'KDEMC', 'WRS'),
       lwd = c(3, 3, 3, 3, 3),
       lty = c(1, 2, 3, 4, 5),
       col = c(rep('black', 2), rep('red', 3)),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

plot(x = dim, y = log(sapply(1:6, function(d) sum(unlist(Poisson_hc_16[[d]]$time)))),
     ylim = c(-2, 10),
     xlab = '',
     ylab = '',
     xaxt = 'n', lwd = 3, pch = , type = 'l')
mtext('Dimension', 1, 2.75, font = 2, cex = 1.5)
mtext('log(Time elapsed in seconds)', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=c(seq(0, 0.9, 0.1), 0.95), labels=c("0.0", c(seq(0.1, 0.9, 0.1), 0.95)), font = 2, cex = 1.5)
axis(1, at=dim, labels = dim, font = 2, cex = 1.5)
axis(2, at=seq(-2, 12, 2), labels = seq(-2, 12, 2), font = 2, cex = 1.5)
axis(2, at=seq(-2, 12, 1), labels=rep("", 15), lwd.ticks = 0.5)
lines(x = dim, y = log(sapply(1:6, function(d) sum(unlist(NB_hc_16[[d]]$time)))),
      lty = 2, lwd = 3)
lines(x = dim, y = log(sapply(1:6, function(d) consensus_mat_16[[d]]$time)),
      lty = 3, lwd = 3, col = 'red')
lines(x = dim, y = log(sapply(1:6, function(d) neiswanger_false_16[[d]]$time)),
      lty = 4, lwd = 3, col = 'red')
lines(x = dim, y = log(sapply(1:6, function(d) weierstrass_rejection_16[[d]]$time)),
      lty = 5, lwd = 3, col = 'red')
legend(x = 5, y = 10,
       legend = c('D&C-MCF (Poisson)', 'D&C-MCF (NB)', 'CMC', 'KDEMC', 'WRS'),
       lwd = c(3, 3, 3, 3, 3),
       lty = c(1, 2, 3, 4, 5),
       col = c(rep('black', 2), rep('red', 3)),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

##### Save data #####

save.image('SD16_dim_study.RData')
