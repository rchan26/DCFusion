library(DCFusion)

seed <- 1994
set.seed(seed)
nsamples <- 10000
C <- 8
corr <- 0.9
opt_bw <- ((4)/(3*nsamples))^(1/5)
diffusion_estimator <- 'NB'
resampling_method <- 'resid'
number_of_replicates <- 10
ESS_threshold <- 0.5
CESS_0_threshold <- 0.05
CESS_j_threshold <- 0.05
vanilla_b <- 1
data_size <- 1
n_cores <- parallel::detectCores()
dimension <- c(1,2,4,6,8,10,12,14,16,18,20,25,30,35,40,45,50,60,70,85,100)
vanilla_guide <- list()
gen_guide <- list()
dc_mcf <- list()
sbf <- list('regular' = list(), 'adaptive' = list())
gbf <- list('regular' = list(), 'adaptive' = list())
dc_gbf <- list('regular' = list(), 'adaptive' = list())
iid_sampling <- list()

collect_results <- function(results, dc, dim, marg_means, marg_sds, seed) {
  if (dc) {
    print(paste('time:', sum(unlist(results$time))))
    print(paste('log(time):', log(sum(unlist(results$time)))))
    return(list('CESS' = results$CESS,
                'time' = results$time,
                'time_mesh' = results$time_mesh,
                'elapsed_time' = results$elapsed_time,
                'resampled' = results$resampled,
                'ESS' = results$ESS,
                'E_nu_j' = results$E_nu_j,
                'chosen' = results$chosen,
                'mesh_terms' = results$mesh_terms,
                'k4_choice' = results$k4_choice,
                'recommended_mesh' = results$recommended_mesh,
                'IAD' = integrated_abs_distance_multiGaussian(fusion_post = resample_particle_y_samples(
                  particle_set = results$particles[[1]],
                  multivariate = TRUE,
                  resampling_method = resampling_method,
                  seed = seed)$y_samples,
                  marg_means = marg_means,
                  marg_sds = marg_sds,
                  bw = rep(opt_bw, dim))))
  } else {
    print(paste('time:', results$time))
    print(paste('log(time):', log(results$time)))
    return(list('CESS' = results$CESS,
                'time_mesh' = results$particles$time_mesh,
                'time' = results$time,
                'elapsed_time' = results$elapsed_time,
                'resampled' = results$resampled,
                'ESS' = results$ESS,
                'E_nu_j' = results$E_nu_j,
                'chosen' = results$chosen,
                'mesh_terms' = results$mesh_terms,
                'k4_choice' = results$k4_choice,
                'IAD' = integrated_abs_distance_multiGaussian(fusion_post = resample_particle_y_samples(
                  particle_set = results$particles,
                  multivariate = TRUE,
                  resampling_method = resampling_method,
                  seed = seed)$y_samples,
                  marg_means = marg_means,
                  marg_sds = marg_sds,
                  bw = rep(opt_bw, dim))))
  }
}

for (d in 1:length(dimension)) {
  print(paste('##### d:', d, '#####'))
  print(paste('%%%%% dimension:', dimension[d], '%%%%%'))
  dc_mcf[[d]] <- list()
  gen_guide[[d]] <- list()
  gbf$regular[[d]] <- list()
  gbf$adaptive[[d]] <- list()
  vanilla_guide[[d]] <- list()
  sbf$regular[[d]] <- list()
  sbf$adaptive[[d]] <- list()
  dc_gbf$regular[[d]] <- list()
  dc_gbf$adaptive[[d]] <- list()
  iid_sampling[[d]] <- list()
  for (rep in 1:number_of_replicates) {
    print(paste('rep:', rep))
    set.seed(seed*rep*d)
    mean <- rep(0, dimension[d])
    cov_mat <- matrix(data = corr, nrow = dimension[d], ncol = dimension[d])
    diag(cov_mat) <- 1
    
    ##### iid sampling #####
    pcm <- proc.time()
    target_samples <- mvrnormArma(N = nsamples, mu = mean, Sigma = cov_mat/data_size)
    time_taken <- (pcm-proc.time())['elapsed']
    iid_sampling[[d]][[rep]] <- list('IAD' = integrated_abs_distance_multiGaussian(fusion_post = target_samples,
                                                                                   marg_means = mean,
                                                                                   marg_sds = rep(sqrt(1/data_size), dimension[d]),
                                                                                   bw = rep(opt_bw, dimension[d])),
                                     'time' = time_taken)
    
    
    cov_mat <- (C/data_size)*cov_mat
    input_samples <- lapply(1:C, function(sub) mvrnormArma(N = nsamples, mu = mean, Sigma = cov_mat))
    if (dimension[d]==1) {
      sub_posterior_means <- as.matrix(sapply(input_samples, function(sub) apply(sub, 2, mean)))
    } else {
      sub_posterior_means <- t(sapply(input_samples, function(sub) apply(sub, 2, mean)))
    }
    
    if (dimension[d] <= 30) {
      ##### Divide-and-Conquer GMCF #####
      print('### performing D&C-Monte Carlo Fusion')
      input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                                  multivariate = TRUE)
      dc_mc_fusion <- bal_binary_GBF_multiGaussian(N_schedule = rep(nsamples, 3),
                                                   m_schedule = rep(2, 3),
                                                   time_mesh = c(0, 1),
                                                   base_samples = input_samples,
                                                   L = 4,
                                                   dim = dimension[d],
                                                   mean_vecs = rep(list(mean), C),
                                                   Sigmas = rep(list(cov_mat), C),
                                                   C = 8,
                                                   precondition = TRUE,
                                                   resampling_method = resampling_method,
                                                   ESS_threshold = ESS_threshold,
                                                   adaptive_mesh = FALSE,
                                                   diffusion_estimator = diffusion_estimator,
                                                   seed = seed*d*rep)
      dc_mcf[[d]][[rep]] <- collect_results(results = dc_mc_fusion,
                                            dc = TRUE,
                                            dim = dimension[d],
                                            marg_means = mean,
                                            marg_sds = rep(sqrt(1/data_size), dimension[d]),
                                            seed = seed*d*rep)

      ##### generalised BF #####
      print('### performing Generalised Bayesian Fusion (with recommended T, regular mesh)')
      gen_guide[[d]][[rep]] <- BF_guidance(condition = 'SH',
                                           CESS_0_threshold = CESS_0_threshold,
                                           CESS_j_threshold = CESS_j_threshold,
                                           sub_posterior_samples = input_samples,
                                           C = C,
                                           d = dimension[d],
                                           data_size = data_size,
                                           sub_posterior_means = sub_posterior_means,
                                           precondition_matrices = lapply(input_samples, cov),
                                           inv_precondition_matrices = lapply(input_samples, function(sub) solve(cov(sub))),
                                           Lambda = inverse_sum_matrices(lapply(input_samples, function(sub) solve(cov(sub)))),
                                           vanilla = FALSE)
      input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                                  multivariate = TRUE,
                                                  number_of_steps = length(gen_guide[[d]][[rep]]$mesh))
      gbf_regular <- parallel_GBF_multiGaussian(particles_to_fuse = input_particles,
                                                N = nsamples,
                                                m = C,
                                                time_mesh = gen_guide[[d]][[rep]]$mesh,
                                                dim = dimension[d],
                                                mean_vecs = rep(list(mean), C),
                                                Sigmas = rep(list(cov_mat), C),
                                                precondition_matrices = lapply(input_samples, cov),
                                                resampling_method = resampling_method,
                                                ESS_threshold = ESS_threshold,
                                                diffusion_estimator = diffusion_estimator,
                                                seed = seed*d*rep)
      gbf$regular[[d]][[rep]] <- collect_results(results = gbf_regular,
                                                 dc = FALSE,
                                                 dim = dimension[d],
                                                 marg_means = mean,
                                                 marg_sds = rep(sqrt(1/data_size), dimension[d]),
                                                 seed = seed*d*rep)
      print('### performing Generalised Bayesian Fusion (with recommended T, adaptive mesh)')
      gbf_adaptive <- parallel_GBF_multiGaussian(particles_to_fuse = input_particles,
                                                 N = nsamples,
                                                 m = C,
                                                 time_mesh = gen_guide[[d]][[rep]]$mesh,
                                                 dim = dimension[d],
                                                 mean_vecs = rep(list(mean), C),
                                                 Sigmas = rep(list(cov_mat), C),
                                                 precondition_matrices = lapply(input_samples, cov),
                                                 resampling_method = resampling_method,
                                                 ESS_threshold = ESS_threshold,
                                                 sub_posterior_means = sub_posterior_means,
                                                 adaptive_mesh = TRUE,
                                                 adaptive_mesh_parameters = list('CESS_j_threshold' = CESS_j_threshold,
                                                                                 'vanilla' = FALSE),
                                                 diffusion_estimator = diffusion_estimator,
                                                 seed = seed*d*rep)
      gbf$adaptive[[d]][[rep]] <- collect_results(results = gbf_adaptive,
                                                  dc = FALSE,
                                                  dim = dimension[d],
                                                  marg_means = mean,
                                                  marg_sds = rep(sqrt(1/data_size), dimension[d]),
                                                  seed = seed*d*rep)
    }

    if (dimension[d] <= 10) {
      ##### standard BF #####
      print('### performing standard Bayesian Fusion (with recommended T, regular mesh)')
      vanilla_guide[[d]][[rep]] <- BF_guidance(condition = 'SH',
                                               CESS_0_threshold = CESS_0_threshold,
                                               CESS_j_threshold = CESS_j_threshold,
                                               sub_posterior_samples = input_samples,
                                               C = C,
                                               d = dimension[d],
                                               data_size = data_size,
                                               b = vanilla_b,
                                               sub_posterior_means = sub_posterior_means,
                                               vanilla = TRUE)
      input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                                  multivariate = TRUE,
                                                  number_of_steps = length(vanilla_guide[[d]][[rep]]$mesh))
      sbf_regular <- parallel_GBF_multiGaussian(particles_to_fuse = input_particles,
                                                N = nsamples,
                                                m = C,
                                                time_mesh = vanilla_guide[[d]][[rep]]$mesh,
                                                dim = dimension[d],
                                                mean_vecs = rep(list(mean), C),
                                                Sigmas = rep(list(cov_mat), C),
                                                precondition_matrices = rep(list(diag(1,dimension[d])), C),
                                                resampling_method = resampling_method,
                                                ESS_threshold = ESS_threshold,
                                                diffusion_estimator = diffusion_estimator,
                                                seed = seed*d*rep)
      sbf$regular[[d]][[rep]] <- collect_results(results = sbf_regular,
                                                 dc = FALSE,
                                                 dim = dimension[d],
                                                 marg_means = mean,
                                                 marg_sds = rep(sqrt(1/data_size), dimension[d]),
                                                 seed = seed*d*rep)
      print('### performing standard Bayesian Fusion (with recommended T, adaptive mesh)')
      sbf_adaptive <- parallel_GBF_multiGaussian(particles_to_fuse = input_particles,
                                                 N = nsamples,
                                                 m = C,
                                                 time_mesh = vanilla_guide[[d]][[rep]]$mesh,
                                                 dim = dimension[d],
                                                 mean_vecs = rep(list(mean), C),
                                                 Sigmas = rep(list(cov_mat), C),
                                                 precondition_matrices = rep(list(diag(1,dimension[d])), C),
                                                 resampling_method = resampling_method,
                                                 ESS_threshold = ESS_threshold,
                                                 sub_posterior_means = sub_posterior_means,
                                                 adaptive_mesh = TRUE,
                                                 adaptive_mesh_parameters = list('data_size' = data_size,
                                                                                 'b' = vanilla_b,
                                                                                 'CESS_j_threshold' = CESS_j_threshold,
                                                                                 'vanilla' = TRUE),
                                                 diffusion_estimator = diffusion_estimator,
                                                 seed = seed*d*rep)
      sbf$adaptive[[d]][[rep]] <- collect_results(results = sbf_adaptive,
                                                  dc = FALSE,
                                                  dim = dimension[d],
                                                  marg_means = mean,
                                                  marg_sds = rep(sqrt(1/data_size), dimension[d]),
                                                  seed = seed*d*rep)
    }

    ##### Divide-and-Conquer GBF #####
    print('### performing D&C-Generalised Bayesian Fusion (with recommended T, regular mesh)')
    dc_gbf_regular <- bal_binary_GBF_multiGaussian(N_schedule = rep(nsamples, 3),
                                                   m_schedule = rep(2, 3),
                                                   base_samples = input_samples,
                                                   L = 4,
                                                   dim = dimension[d],
                                                   mean_vecs = rep(list(mean), C),
                                                   Sigmas = rep(list(cov_mat), C),
                                                   C = C,
                                                   precondition = TRUE,
                                                   resampling_method = resampling_method,
                                                   ESS_threshold = ESS_threshold,
                                                   adaptive_mesh = FALSE,
                                                   mesh_parameters = list('condition' = 'SH',
                                                                          'CESS_0_threshold' = CESS_0_threshold,
                                                                          'CESS_j_threshold' = CESS_j_threshold,
                                                                          'vanilla' = FALSE),
                                                   diffusion_estimator = diffusion_estimator,
                                                   seed = seed*d*rep)
    dc_gbf$regular[[d]][[rep]] <- collect_results(results = dc_gbf_regular,
                                                  dc = TRUE,
                                                  dim = dimension[d],
                                                  marg_means = mean,
                                                  marg_sds = rep(sqrt(1/data_size), dimension[d]),
                                                  seed = seed*d*rep)
    print('### performing D&C-Generalised Bayesian Fusion (with recommended T, adaptive mesh)')
    dc_gbf_adaptive <- bal_binary_GBF_multiGaussian(N_schedule = rep(nsamples, 3),
                                                    m_schedule = rep(2, 3),
                                                    base_samples = input_samples,
                                                    L = 4,
                                                    dim = dimension[d],
                                                    mean_vecs = rep(list(mean), C),
                                                    Sigmas = rep(list(cov_mat), C),
                                                    C = C,
                                                    precondition = TRUE,
                                                    resampling_method = resampling_method,
                                                    ESS_threshold = ESS_threshold,
                                                    adaptive_mesh = TRUE,
                                                    mesh_parameters = list('condition' = 'SH',
                                                                           'CESS_0_threshold' = CESS_0_threshold,
                                                                           'CESS_j_threshold' = CESS_j_threshold,
                                                                           'vanilla' = FALSE),
                                                    diffusion_estimator = diffusion_estimator,
                                                    seed = seed*d*rep)
    dc_gbf$adaptive[[d]][[rep]] <- collect_results(results = dc_gbf_adaptive,
                                                   dc = TRUE,
                                                   dim = dimension[d],
                                                   marg_means = mean,
                                                   marg_sds = rep(sqrt(1/data_size), dimension[d]),
                                                   seed = seed*d*rep)1

    print('saving progress')
    save.image('dimension_study_similar_means_replicates.RData')
  }
}

##### Paper: IAD #####
plot(x = dimension,
     y = sapply(1:length(dimension), function(i) {
       mean(sapply(1:number_of_replicates, function(rep) dc_gbf$regular[[i]][[rep]]$IAD))
     }),
     type = 'b', pch = 1, lty = 1, lwd = 3, ylim = c(0,1), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = dimension,
      y = sapply(1:length(dimension), function(i) {
        mean(sapply(1:number_of_replicates, function(rep) dc_gbf$adaptive[[i]][[rep]]$IAD))
      }),
      pch = 2, lty = 2, lwd = 3, type = 'b', col = 'red')
lines(x = dimension,
      y = sapply(1:length(dimension), function(i) {
        mean(sapply(1:number_of_replicates, function(rep) iid_sampling[[i]][[rep]]$IAD))
      }),
      pch = 8, lty = 1, lwd = 3, type = 'b', col = 'pink')
lines(x = dimension[1:13],
      y = sapply(1:length(dimension[1:13]), function(i) {
        mean(sapply(1:number_of_replicates, function(rep) gbf$regular[[i]][[rep]]$IAD))
      }),
      pch = 3, lty = 3, lwd = 3, type = 'b', col = 'blue')
lines(x = dimension[1:13],
      y = sapply(1:length(dimension[1:13]), function(i) {
        mean(sapply(1:number_of_replicates, function(rep) gbf$adaptive[[i]][[rep]]$IAD))
      }),
      pch = 4, lty = 4, lwd = 3, type = 'b', col = 'green')
lines(x = dimension[1:13],
      y = sapply(1:length(dimension[1:13]), function(i) {
        mean(sapply(1:number_of_replicates, function(rep) dc_mcf[[i]][[rep]]$IAD))
      }),
      pch = 7, lty = 1, lwd = 3, type = 'b', col = 'black')
lines(x = dimension[1:6],
      y = sapply(1:length(dimension[1:6]), function(i) {
        mean(sapply(1:number_of_replicates, function(rep) sbf$regular[[i]][[rep]]$IAD))
      }),
      pch = 5, lty = 5, lwd = 3, type = 'b', col = 'purple')
lines(x = dimension[1:6],
      y = sapply(1:length(dimension[1:6]), function(i) {
        mean(sapply(1:number_of_replicates, function(rep) sbf$adaptive[[i]][[rep]]$IAD))
      }),
      pch = 6, lty = 6, lwd = 3, type = 'b', col = 'orange')
for (i in 1:length(dimension)) {
  IAD <- sapply(1:number_of_replicates, function(rep) dc_gbf$regular[[i]][[rep]]$IAD)
  points(x = rep(dimension[i], length(IAD)), y = IAD, cex = 0.5, pch = 1, lwd = 1.5)
}
for (i in 1:length(dimension)) {
  IAD <- sapply(1:number_of_replicates, function(rep) dc_gbf$adaptive[[i]][[rep]]$IAD)
  points(x = rep(dimension[i], length(IAD)), y = IAD, cex = 0.5, pch = 2, lwd = 1.5, col = 'red')
}
for (i in 1:length(dimension[1:6])) {
  IAD <- sapply(1:number_of_replicates, function(rep) iid_sampling[[i]][[rep]]$IAD)
  points(x = rep(dimension[i], length(IAD)), y = IAD, cex = 0.5, pch = 8, lwd = 1.5, col = 'pink')
}
for (i in 1:length(dimension[1:13])) {
  IAD <- sapply(1:number_of_replicates, function(rep) gbf$regular[[i]][[rep]]$IAD)
  points(x = rep(dimension[i], length(IAD)), y = IAD, cex = 0.5, pch = 3, lwd = 1.5, col = 'blue')
}
for (i in 1:length(dimension[1:13])) {
  IAD <- sapply(1:number_of_replicates, function(rep) gbf$adaptive[[i]][[rep]]$IAD)
  points(x = rep(dimension[i], length(IAD)), y = IAD, cex = 0.5, pch = 4, lwd = 1.5, col = 'green')
}
for (i in 1:length(dimension[1:13])) {
  IAD <- sapply(1:number_of_replicates, function(rep) dc_mcf[[i]][[rep]]$IAD)
  points(x = rep(dimension[i], length(IAD)), y = IAD, cex = 0.5, pch = 7, lwd = 1.5, col = 'black')
}
for (i in 1:length(dimension[1:6])) {
  IAD <- sapply(1:number_of_replicates, function(rep) sbf$regular[[i]][[rep]]$IAD)
  points(x = rep(dimension[i], length(IAD)), y = IAD, cex = 0.5, pch = 5, lwd = 1.5, col = 'purple')
}
for (i in 1:length(dimension[1:6])) {
  IAD <- sapply(1:number_of_replicates, function(rep) sbf$adaptive[[i]][[rep]]$IAD)
  points(x = rep(dimension[i], length(IAD)), y = IAD, cex = 0.5, pch = 6, lwd = 1.5, col = 'orange')
}
axis(1, at = seq(0, dimension[length(dimension)], 10),
     labels = seq(0, dimension[length(dimension)], 10), font = 2, cex = 1.5)
axis(1, at = seq(0, dimension[length(dimension)], 5), labels = rep("", 21), lwd.ticks = 0.5)
mtext('Dimension', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1.6, 0.1), labels = c("0.0", seq(0.1, 0.9, 0.1), "1.0", seq(1.1, 1.6, 0.1)),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1.6, 0.1), labels=rep("", 17), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
legend(x = 30, y = 1,
       legend = c('D&C-GBF (regular)',
                  'D&C-GBF (adaptive)',
                  'GBF (regular)',
                  'GBF (adaptive)',
                  'D&C-MCF',
                  'BF (regular)',
                  'BF (adaptive)',
                  'iid sampling'),
       col = c('black', 'red', 'blue', 'green', 'black', 'purple', 'orange', 'pink'),
       lty = c(1, 2, 3, 4, 1, 5, 6, 1),
       pch = c(1, 2, 3, 4, 7, 5, 6, 8),
       lwd = rep(3, 8),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

##### Paper: IAD #####
plot(x = dimension,
     y = sapply(1:length(dimension), function(i) {
       var(sapply(1:number_of_replicates, function(rep) dc_gbf$regular[[i]][[rep]]$IAD))
     }),
     type = 'b', pch = 1, lty = 1, lwd = 3, ylim = c(0,0.05), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = dimension,
      y = sapply(1:length(dimension), function(i) {
        var(sapply(1:number_of_replicates, function(rep) dc_gbf$adaptive[[i]][[rep]]$IAD))
      }),
      pch = 2, lty = 2, lwd = 3, type = 'b', col = 'red')
lines(x = dimension,
      y = sapply(1:length(dimension), function(i) {
        var(sapply(1:number_of_replicates, function(rep) iid_sampling[[i]][[rep]]$IAD))
      }),
      pch = 8, lty = 1, lwd = 3, type = 'b', col = 'pink')
lines(x = dimension[1:13],
      y = sapply(1:length(dimension[1:13]), function(i) {
        var(sapply(1:number_of_replicates, function(rep) gbf$regular[[i]][[rep]]$IAD))
      }),
      pch = 3, lty = 3, lwd = 3, type = 'b', col = 'blue')
lines(x = dimension[1:13],
      y = sapply(1:length(dimension[1:13]), function(i) {
        var(sapply(1:number_of_replicates, function(rep) gbf$adaptive[[i]][[rep]]$IAD))
      }),
      pch = 4, lty = 4, lwd = 3, type = 'b', col = 'green')
lines(x = dimension[1:13],
      y = sapply(1:length(dimension[1:13]), function(i) {
        var(sapply(1:number_of_replicates, function(rep) dc_mcf[[i]][[rep]]$IAD))
      }),
      pch = 7, lty = 1, lwd = 3, type = 'b', col = 'black')
lines(x = dimension[1:6],
      y = sapply(1:length(dimension[1:6]), function(i) {
        var(sapply(1:number_of_replicates, function(rep) sbf$regular[[i]][[rep]]$IAD))
      }),
      pch = 5, lty = 5, lwd = 3, type = 'b', col = 'purple')
lines(x = dimension[1:6],
      y = sapply(1:length(dimension[1:6]), function(i) {
        var(sapply(1:number_of_replicates, function(rep) sbf$adaptive[[i]][[rep]]$IAD))
      }),
      pch = 6, lty = 6, lwd = 3, type = 'b', col = 'orange')
axis(1, at = seq(0, dimension[length(dimension)], 10),
     labels = seq(0, dimension[length(dimension)], 10), font = 2, cex = 1.5)
axis(1, at = seq(0, dimension[length(dimension)], 5), labels = rep("", 21), lwd.ticks = 0.5)
mtext('Dimension', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 1.6, 0.1), labels = c("0.0", seq(0.1, 0.9, 0.1), "1.0", seq(1.1, 1.6, 0.1)),
     font = 2, cex = 1.5)
axis(2, at = seq(0, 1.6, 0.1), labels=rep("", 17), lwd.ticks = 0.5,
     font = 2, cex = 1.5)
mtext('Integrated Absolute Distance', 2, 2.75, font = 2, cex = 1.5)
legend(x = 30, y = 1,
       legend = c('D&C-GBF (regular)',
                  'D&C-GBF (adaptive)',
                  'GBF (regular)',
                  'GBF (adaptive)',
                  'D&C-MCF',
                  'BF (regular)',
                  'BF (adaptive)',
                  'iid sampling'),
       col = c('black', 'red', 'blue', 'green', 'black', 'purple', 'orange', 'pink'),
       lty = c(1, 2, 3, 4, 1, 5, 6, 1),
       pch = c(1, 2, 3, 4, 7, 5, 6, 8),
       lwd = rep(3, 8),
       cex = 1.25,
       text.font = 2,
       bty = 'n')

##### Paper: time #####
plot(x = log(dimension, 2),
     y = sapply(1:length(dimension), function(i) {
       mean(log(sapply(1:number_of_replicates, function(rep) sum(unlist(dc_gbf$regular[[i]][[rep]]$time))), 2))
     }),
     type = 'b', pch = 1, lty = 1, lwd = 3, ylim = c(2,14), xaxt = 'n', yaxt ='n', xlab = '', ylab = '')
lines(x = log(dimension, 2),
      y = sapply(1:length(dimension), function(i) {
        mean(log(sapply(1:number_of_replicates, function(rep) sum(unlist(dc_gbf$adaptive[[i]][[rep]]$time))), 2))
      }),
      pch = 2, lty = 2, lwd = 3, type = 'b', col = 'red')
lines(x = log(dimension[1:13], 2),
      y = sapply(1:length(dimension[1:13]), function(i) {
        mean(log(sapply(1:number_of_replicates, function(rep) gbf$regular[[i]][[rep]]$time), 2))
      }),
      pch = 3, lty = 3, lwd = 3, type = 'b', col = 'blue')
lines(x = log(dimension[1:13], 2),
      y = sapply(1:length(dimension[1:13]), function(i) {
        mean(log(sapply(1:number_of_replicates, function(rep) gbf$adaptive[[i]][[rep]]$time), 2))
      }),
      pch = 4, lty = 4, lwd = 3, type = 'b', col = 'green')
lines(x = log(dimension[1:13], 2),
      y = sapply(1:length(dimension[1:13]), function(i) {
        mean(log(sapply(1:number_of_replicates, function(rep) sum(unlist(dc_mcf[[i]][[rep]]$time))), 2))
      }),
      pch = 7, lty = 1, lwd = 3, type = 'b', col = 'black')
lines(x = log(dimension[1:6], 2),
      y = sapply(1:length(dimension[1:6]), function(i) {
        mean(log(sapply(1:number_of_replicates, function(rep) sbf$adaptive[[i]][[rep]]$time), 2))
      }),
      pch = 4, lty = 4, lwd = 3, type = 'b', col = 'purple')
lines(x = log(dimension[1:6], 2),
      y = sapply(1:length(dimension[1:6]), function(i) {
        mean(log(sapply(1:number_of_replicates, function(rep) sbf$regular[[i]][[rep]]$time), 2))
      }),
      pch = 3, lty = 3, lwd = 3, type = 'b', col = 'orange')
axis(1, at = 0:6,
     labels = 0:6, font = 2, cex = 1.5)

axis(1, at = seq(0, 7, 0.5), labels = rep("", 15), lwd.ticks = 0.5)
mtext('log(Dimension, 2)', 1, 2.75, font = 2, cex = 1.5)
axis(2, at = seq(0, 14, 1), labels = seq(0, 14, 1), font = 2, cex = 1.5)
mtext('log(Elapsed time in seconds, 2)', 2, 2.75, font = 2, cex = 1.5)
legend(x = 0, y = 14,
       legend = c('D&C-GBF (regular)',
                  'D&C-GBF (adaptive)',
                  'GBF (regular)',
                  'GBF (adaptive)',
                  'D&C-MCF',
                  'BF (regular)',
                  'BF (adaptive)'),
       col = c('black', 'red', 'blue', 'green', 'black', 'purple', 'orange'),
       lty = c(1, 2, 3, 4, 1, 5, 6),
       pch = c(1, 2, 3, 4, 7, 5, 6),
       lwd = rep(3, 7),
       cex = 1.25,
       text.font = 2,
       bty = 'n')
