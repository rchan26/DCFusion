#' @export
rho_IS_additional <- function(particle_set,
                              time,
                              inv_precondition_matrices,
                              inv_gamma_matrices) {
  log_rho_additional_weights <- sapply(1:particle_set$N, function(i) {
    log_rho_multivariate_additional(y = particle_set$y_samples[i,],
                                    x_mean = particle_set$x_means[i,],
                                    time = time,
                                    inv_precondition_matrices = inv_precondition_matrices,
                                    inv_gamma_matrices = inv_gamma_matrices)})
  # ---------- update particle set
  # update the weights and return updated particle set
  updated_weights <- particle_set$log_weights + log_rho_additional_weights
  # normalise weights
  norm_weights <- particle_ESS(log_weights = updated_weights)
  particle_set$log_weights <- norm_weights$log_weights
  particle_set$lw$rho_additional <- log_rho_additional_weights
  particle_set$normalised_lw$rho_additional <- particle_ESS(log_rho_additional_weights)$log_weights
  particle_set$normalised_weights <- norm_weights$normalised_weights
  particle_set$ESS <- norm_weights$ESS
  # calculate the conditional ESS (i.e. the 1/sum(inc_change^2))
  # where inc_change is the incremental change in weight (= log_rho_additional_weights)
  particle_set$CESS['rho_additional'] <- particle_ESS(log_weights = log_rho_additional_weights)$ESS
  # set the resampled indicator to FALSE
  particle_set$resampled['rho_additional'] <- FALSE
  return(particle_set)
}

#' @export
parallel_fusion_SMC_BLR_alt <- function(particles_to_fuse,
                                        N,
                                        m,
                                        time,
                                        dim,
                                        data_split,
                                        prior_means,
                                        prior_variances,
                                        C,
                                        precondition_matrices,
                                        gamma_matrices,
                                        resampling_method = 'multi',
                                        cv_location = 'hypercube_centre',
                                        diffusion_estimator = 'Poisson',
                                        beta_NB = 10,
                                        gamma_NB_n_points = 2,
                                        seed = NULL,
                                        n_cores = parallel::detectCores(),
                                        cl = NULL,
                                        level = 1,
                                        node = 1,
                                        print_progress_iters = 1000) {
  if (!is.list(particles_to_fuse) | (length(particles_to_fuse)!=m)) {
    stop("parallel_fusion_SMC_BLR_alt: particles_to_fuse must be a list of length m")
  } else if (!all(sapply(particles_to_fuse, function(sub_posterior) ("particle" %in% class(sub_posterior))))) {
    stop("parallel_fusion_SMC_BLR_alt: particles in particles_to_fuse must be \"particle\" objects")
  } else if (!all(sapply(particles_to_fuse, function(sub_posterior) is.matrix(sub_posterior$y_samples)))) {
    stop("parallel_fusion_SMC_BLR_alt: the particles' samples for y should all be matrices")
  } else if (!all(sapply(particles_to_fuse, function(sub_posterior) ncol(sub_posterior$y_samples)==dim))) {
    stop("parallel_fusion_SMC_BLR_alt: the particles' samples for y should all be matrices with dim columns")
  } else if (!is.list(data_split) | length(data_split)!=m) {
    stop("parallel_fusion_SMC_BLR_alt: data_split must be a list of length m")
  } else if (!all(sapply(data_split, function(sub_posterior) (is.list(sub_posterior) & identical(names(sub_posterior), c("y", "X", "full_data_count", "design_count")))))) {
    stop("parallel_fusion_SMC_BLR_alt: each item in data_split must be a list of length 4 with names \'y\', \'X\', \'full_data_count\', \'design_count\'")
  } else if (!all(sapply(1:m, function(i) is.vector(data_split[[i]]$y)))) {
    stop("parallel_fusion_SMC_BLR_alt: for each i in 1:m, data_split[[i]]$y must be a vector")
  } else if (!all(sapply(1:m, function(i) is.matrix(data_split[[i]]$X)))) {
    stop("parallel_fusion_SMC_BLR_alt: for each i in 1:m, data_split[[i]]$X must be a matrix")
  } else if (!all(sapply(1:m, function(i) ncol(data_split[[i]]$X)==dim))) {
    stop("parallel_fusion_SMC_BLR_alt: for each i in 1:m, data_split[[i]]$X must be a matrix with dim columns")
  } else if (!all(sapply(1:m, function(i) length(data_split[[i]]$y)==nrow(data_split[[i]]$X)))) {
    stop("parallel_fusion_SMC_BLR_alt: for each i in 1:m, length(data_split[[i]]$y) and nrow(data_split[[i]]$X) must be equal")
  } else if (!all(sapply(1:m, function(i) is.data.frame(data_split[[i]]$full_data_count)))) {
    stop("parallel_fusion_SMC_BLR_alt: for each i in 1:m, data_split[[i]]$full_data_count must be a data frame")
  } else if (!all(sapply(1:m, function(i) is.data.frame(data_split[[i]]$design_count)))) {
    stop("parallel_fusion_SMC_BLR_alt: for each i in 1:m, data_split[[i]]$design_count must be a data frame")
  } else if (!is.vector(prior_means) | length(prior_means)!=dim) {
    stop("parallel_fusion_SMC_BLR_alt: prior_means must be vectors of length dim")
  } else if (!is.vector(prior_variances) | length(prior_variances)!=dim) {
    stop("parallel_fusion_SMC_BLR_alt: prior_variances must be vectors of length dim")
  } else if (!is.list(precondition_matrices) | (length(precondition_matrices)!=m)) {
    stop("parallel_fusion_SMC_BLR_alt: precondition_matrices must be a list of length m")
  } else if (!(diffusion_estimator %in% c('Poisson', 'NB'))) {
    stop("parallel_fusion_SMC_BLR_alt: diffusion_estimator must be set to either \'Poisson\' or \'NB\'")
  } else if ((cv_location != 'mode') & (cv_location != 'hypercube_centre')) {
    stop("parallel_fusion_SMC_BLR_alt: cv_location must be either \"mode\" or \"hypercube_centre\"")
  }
  # set a seed if one is supplied
  if (!is.null(seed)) {
    set.seed(seed)
  }
  # ---------- resample the particles if they do not have equal weights
  # set a seed if one is supplied
  # check if the resampled indicator if FALSE
  # also check if there are enough samples
  for (c in 1:length(particles_to_fuse)) {
    if ((!particles_to_fuse[[c]]$resampled['Q']) | (particles_to_fuse[[c]]$N!=N)) {
      particles_to_fuse[[c]] <- resample_particle_y_samples(N = N,
                                                            particle_set = particles_to_fuse[[c]],
                                                            multivariate = TRUE,
                                                            resampling_method = resampling_method,
                                                            seed = seed)
    }
  }
  # start time recording
  pcm <- proc.time()
  # ---------- first importance sampling step
  inv_precondition_matrices <- lapply(precondition_matrices, solve)
  inv_gamma_matrices <- lapply(gamma_matrices, solve)
  # importance sampling for rho step
  particles <- rho_IS_multivariate(particles_to_fuse = particles_to_fuse,
                                   dim = dim,
                                   N = N,
                                   m = m,
                                   time = time,
                                   inv_precondition_matrices = inv_gamma_matrices,
                                   inverse_sum_inv_precondition_matrices = inverse_sum_matrices(inv_gamma_matrices),
                                   n_cores = n_cores,
                                   cl = cl)
  # record ESS and CESS after rho step
  ESS <- c('rho' = particles$ESS)
  CESS <- c('rho' = particles$CESS['rho'])
  # ---------- second importance sampling step
  # unbiased estimator for Q
  particles <- Q_IS_BLR(particle_set = particles,
                        m = m,
                        time = time,
                        dim = dim,
                        data_split = data_split,
                        prior_means = prior_means,
                        prior_variances = prior_variances,
                        C = C,
                        proposal_cov = calculate_proposal_cov(time = time, weights = inv_precondition_matrices),
                        precondition_matrices = gamma_matrices,
                        inv_precondition_matrices = inv_gamma_matrices,
                        cv_location = cv_location,
                        diffusion_estimator = diffusion_estimator,
                        beta_NB = beta_NB,
                        gamma_NB_n_points = gamma_NB_n_points,
                        seed = seed,
                        n_cores = n_cores,
                        cl = cl,
                        level = level,
                        node = node,
                        print_progress_iters = print_progress_iters)
  # record ESS and CESS after Q step
  ESS['Q'] <- particles$ESS
  CESS['Q'] <- particles$CESS['Q']
  names(CESS) <- c('rho', 'Q')
  # record proposed samples
  proposed_samples <- particles$y_samples
  # ---------- third importance sampling step
  particles <- rho_IS_additional(particle_set = particles,
                                 time = time,
                                 inv_precondition_matrices = inv_precondition_matrices,
                                 inv_gamma_matrices = inv_gamma_matrices)
  # record ESS and CESS after Q step
  ESS['rho_additional'] <- particles$ESS
  CESS['rho_additional'] <- particles$CESS['rho_additional']
  names(CESS) <- c('rho', 'Q', 'rho_additional')
  if (identical(precondition_matrices, rep(list(diag(1, dim)), m))) {
    return(list('particles' = particles,
                'proposed_samples' = proposed_samples,
                'time' = (proc.time()-pcm)['elapsed'],
                'ESS' = ESS,
                'CESS' = CESS,
                'precondition_matrices' = list(diag(1, dim), precondition_matrices),
                'combined_data' = combine_data(list_of_data = data_split, dim = dim)))
  } else {
    return(list('particles' = particles,
                'proposed_samples' = proposed_samples,
                'time' = (proc.time()-pcm)['elapsed'],
                'ESS' = ESS,
                'CESS' = CESS,
                'precondition_matrices' = list(inverse_sum_matrices(inv_precondition_matrices),
                                               precondition_matrices),
                'combined_data' = combine_data(list_of_data = data_split, dim = dim)))
  }
}
