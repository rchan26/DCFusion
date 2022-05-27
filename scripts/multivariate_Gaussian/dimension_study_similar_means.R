library(DCFusion)

seed <- 1994
set.seed(seed)
nsamples <- 5000
C <- 8
corr <- 0.9
opt_bw <- ((4)/(3*nsamples))^(1/5)
diffusion_estimator <- 'NB'
resampling_method <- 'resid'
ESS_threshold <- 0.5
CESS_0_threshold <- 0.05
CESS_j_threshold <- 0.05
vanilla_b <- 1
data_size <- 1
n_cores <- parallel::detectCores()
dimension <- c(1,2,4,6,8,10,12,14,16,18,20,25,30,35,40,45,50)
# dimension <- c(1,2,4,6,8,10,15,20)
vanilla_guide <- list()
gen_guide <- list()
dc_mcf <- list()
sbf <- list('regular' = list(), 'adaptive' = list())
gbf <- list('regular' = list(), 'adaptive' = list())
dc_gbf <- list('regular' = list(), 'adaptive' = list())

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
  print(paste('##### d:', d))
  print(paste('dimension:', dimension[d]))
  set.seed(seed*d)
  mean <- rep(0, dimension[d])
  cov_mat <- matrix(data = corr, nrow = dimension[d], ncol = dimension[d])
  diag(cov_mat) <- 1
  cov_mat <- (C/data_size)*cov_mat
  input_samples <- lapply(1:C, function(sub) mvrnormArma(N = nsamples, mu = mean, Sigma = cov_mat))
  if (dimension[d]==1) {
    sub_posterior_means <- as.matrix(sapply(input_samples, function(sub) apply(sub, 2, mean)))
  } else {
    sub_posterior_means <- t(sapply(input_samples, function(sub) apply(sub, 2, mean)))
  }
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
                                               seed = seed*d)
  dc_mcf[[d]] <- collect_results(results = dc_mc_fusion,
                                 dc = TRUE,
                                 dim = dimension[d],
                                 marg_means = mean,
                                 marg_sds = rep(sqrt(1/data_size), dimension[d]),
                                 seed = seed*d)
  
  if (dimension[d] <= 10) {
    ##### standard BF #####
    vanilla_guide[[d]] <- BF_guidance(condition = 'SH',
                                      CESS_0_threshold = CESS_0_threshold,
                                      CESS_j_threshold = CESS_j_threshold,
                                      sub_posterior_samples = input_samples,
                                      C = C,
                                      d = dimension[d],
                                      data_size = data_size,
                                      b = vanilla_b,
                                      sub_posterior_means = sub_posterior_means,
                                      vanilla = TRUE)
    print('### performing standard Bayesian Fusion (with recommended T, regular mesh)')
    input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                                multivariate = TRUE,
                                                number_of_steps = length(vanilla_guide[[d]]$mesh))
    sbf_regular <- parallel_GBF_multiGaussian(particles_to_fuse = input_particles,
                                              N = nsamples,
                                              m = C,
                                              time_mesh = vanilla_guide[[d]]$mesh,
                                              dim = dimension[d],
                                              mean_vecs = rep(list(mean), C),
                                              Sigmas = rep(list(cov_mat), C), 
                                              precondition_matrices = rep(list(diag(1,dimension[d])), C),
                                              resampling_method = resampling_method,
                                              ESS_threshold = ESS_threshold,
                                              diffusion_estimator = diffusion_estimator,
                                              seed = seed*d)
    sbf$regular[[d]] <- collect_results(results = sbf_regular,
                                        dc = FALSE,
                                        dim = dimension[d],
                                        marg_means = mean,
                                        marg_sds = rep(sqrt(1/data_size), dimension[d]),
                                        seed = seed*d)
    print('### performing standard Bayesian Fusion (with recommended T, adaptive mesh)')
    sbf_adaptive <- parallel_GBF_multiGaussian(particles_to_fuse = input_particles,
                                               N = nsamples,
                                               m = C,
                                               time_mesh = vanilla_guide[[d]]$mesh,
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
                                               seed = seed*d)
    sbf$adaptive[[d]] <- collect_results(results = sbf_adaptive,
                                         dc = FALSE,
                                         dim = dimension[d],
                                         marg_means = mean,
                                         marg_sds = rep(sqrt(1/data_size), dimension[d]),
                                         seed = seed*d)
    
    ##### generalised BF #####
    gen_guide[[d]] <- BF_guidance(condition = 'SH',
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
    print('### performing Generalised Bayesian Fusion (with recommended T, regular mesh)')
    input_particles <- initialise_particle_sets(samples_to_fuse = input_samples,
                                                multivariate = TRUE,
                                                number_of_steps = length(gen_guide[[d]]$mesh))
    gbf_regular <- parallel_GBF_multiGaussian(particles_to_fuse = input_particles,
                                              N = nsamples,
                                              m = C,
                                              time_mesh = gen_guide[[d]]$mesh,
                                              dim = dimension[d],
                                              mean_vecs = rep(list(mean), C),
                                              Sigmas = rep(list(cov_mat), C),
                                              precondition_matrices = lapply(input_samples, cov),
                                              resampling_method = resampling_method,
                                              ESS_threshold = ESS_threshold,
                                              diffusion_estimator = diffusion_estimator,
                                              seed = seed*d)
    gbf$regular[[d]] <- collect_results(results = gbf_regular,
                                        dc = FALSE,
                                        dim = dimension[d],
                                        marg_means = mean,
                                        marg_sds = rep(sqrt(1/data_size), dimension[d]),
                                        seed = seed*d)
    print('### performing Generalised Bayesian Fusion (with recommended T, adaptive mesh)')
    gbf_adaptive <- parallel_GBF_multiGaussian(particles_to_fuse = input_particles,
                                               N = nsamples,
                                               m = C,
                                               time_mesh = gen_guide[[d]]$mesh,
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
                                               seed = seed*d)
    gbf$adaptive[[d]] <- collect_results(results = gbf_adaptive,
                                         dc = FALSE,
                                         dim = dimension[d],
                                         marg_means = mean,
                                         marg_sds = rep(sqrt(1/data_size), dimension[d]),
                                         seed = seed*d)
  }
  
  ##### Divide-and-Conquer GBF #####
  print('### performing D&C-Generalised Bayesian Fusion (with recommended T, regular mesh)')
  dc_gbf_regular <- bal_binary_GBF_multiGaussian(N_schedule = rep(nsamples, 3),
                                                 m_schedule = rep(2, 3),
                                                 time_mesh = gen_guide[[d]]$mesh,
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
                                                 seed = seed*d)
  dc_gbf$regular[[d]] <- collect_results(results = dc_gbf_regular,
                                         dc = TRUE,
                                         dim = dimension[d],
                                         marg_means = mean,
                                         marg_sds = rep(sqrt(1/data_size), dimension[d]),
                                         seed = seed*d)
  print('### performing D&C-Generalised Bayesian Fusion (with recommended T, adaptive mesh)')
  dc_gbf_adaptive <- bal_binary_GBF_multiGaussian(N_schedule = rep(nsamples, 3),
                                                  m_schedule = rep(2, 3),
                                                  time_mesh = gen_guide[[d]]$mesh,
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
                                                  seed = seed*d)
  dc_gbf$adaptive[[d]] <- collect_results(results = dc_gbf_adaptive,
                                          dc = TRUE,
                                          dim = dimension[d],
                                          marg_means = mean,
                                          marg_sds = rep(sqrt(1/data_size), dimension[d]),
                                          seed = seed*d)
  
  print('saving progress')
  save.image('dimension_study_similar_means_thresh_005_n5000.RData')
}
