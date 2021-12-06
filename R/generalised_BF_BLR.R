#' @export
generalised_rho_j_BLR <- function(particle_set,
                                  m,
                                  time_mesh,
                                  dim,
                                  data_split,
                                  prior_means,
                                  prior_variances,
                                  C,
                                  precondition_matrices,
                                  inv_precondition_matrices,
                                  Lambda,
                                  resampling_method = 'multi',
                                  ESS_threshold = 0.5,
                                  cv_location = 'hypercube_centre',
                                  diffusion_estimator,
                                  beta_NB = 10,
                                  gamma_NB_n_points = 2,
                                  local_bounds = TRUE,
                                  seed = NULL,
                                  n_cores = parallel::detectCores(),
                                  cl = NULL,
                                  level = 1,
                                  node = 1,
                                  print_progress_iters = 1000) {
  if (!("particle" %in% class(particle_set))) {
    stop("generalised_rho_j_BLR: particle_set must be a \"particle\" object")
  } else if (!is.list(data_split) | length(data_split)!=m) {
    stop("generalised_rho_j_BLR: data_split must be a list of length m")
  } else if (!all(sapply(data_split, function(sub_posterior) (is.list(sub_posterior) & identical(names(sub_posterior), c("y", "X", "full_data_count", "design_count")))))) {
    stop("generalised_rho_j_BLR: each item in data_split must be a list of length 4 with names \'y\', \'X\', \'full_data_count\', \'design_count\'")
  } else if (!is.vector(time_mesh)) {
    stop("generalised_rho_j_BLR: time_mesh must be an ordered vector of length >= 2")
  } else if (length(time_mesh) < 2) {
    stop("generalised_rho_j_BLR: time_mesh must be an ordered vector of length >= 2")
  } else if (!identical(time_mesh, sort(time_mesh))) {
    stop("generalised_rho_j_BLR: time_mesh must be an ordered vector of length >= 2")
  } else if (!all(sapply(1:m, function(i) is.vector(data_split[[i]]$y)))) {
    stop("generalised_rho_j_BLR: for each i in 1:m, data_split[[i]]$y must be a vector")
  } else if (!all(sapply(1:m, function(i) is.matrix(data_split[[i]]$X)))) {
    stop("generalised_rho_j_BLR: for each i in 1:m, data_split[[i]]$X must be a matrix")
  } else if (!all(sapply(1:m, function(i) ncol(data_split[[i]]$X)==dim))) {
    stop("generalised_rho_j_BLR: for each i in 1:m, ncol(data_split[[i]]$X) must be equal to dim")
  } else if (!all(sapply(1:m, function(i) length(data_split[[i]]$y)==nrow(data_split[[i]]$X)))) {
    stop("generalised_rho_j_BLR: for each i in 1:m, length(data_split[[i]]$y) and nrow(data_split[[i]]$X) must be equal")
  } else if (!all(sapply(1:m, function(i) is.data.frame(data_split[[i]]$full_data_count)))) {
    stop("generalised_rho_j_BLR: for each i in 1:m, data_split[[i]]$full_data_count must be a data frame")
  } else if (!all(sapply(1:m, function(i) is.data.frame(data_split[[i]]$design_count)))) {
    stop("generalised_rho_j_BLR: for each i in 1:m, data_split[[i]]$design_count must be a data frame")
  } else if (!is.vector(prior_means) | length(prior_means)!=dim) {
    stop("generalised_rho_j_BLR: prior_means must be vectors of length dim")
  } else if (!is.vector(prior_variances) | length(prior_variances)!=dim) {
    stop("generalised_rho_j_BLR: prior_variances must be vectors of length dim")
  } else if (!is.list(precondition_matrices) | (length(precondition_matrices)!=m)) {
    stop("Q_IS_BLR: precondition_matrices must be a list of length m")
  } else if (!is.list(inv_precondition_matrices) | (length(inv_precondition_matrices)!=m)) {
    stop("Q_IS_BLR: inv_precondition_matrices must be a list of length m")
  } else if (!(diffusion_estimator %in% c('Poisson', 'NB'))) {
    stop("generalised_rho_j_BLR: diffusion_estimator must be set to either \'Poisson\' or \'NB\'")
  } else if (!any(class(cl)=="cluster") & !is.null(cl)) {
    stop("generalised_rho_j_BLR: cl must be a \"cluster\" object or NULL")
  }
  if (cv_location == 'mode') {
    cv_location <- lapply(1:m, function(c) {
      MLE <- obtain_LR_MLE(dim = dim, data = data_split[[c]])
      X <- as.matrix(subset(data_split[[c]]$full_data_count, select = -c(y, count)))
      list('beta_hat' = MLE,
           'grad_log_hat' = log_BLR_gradient(beta = MLE,
                                             y_labels = data_split[[c]]$full_data_count$y,
                                             X = X,
                                             X_beta = as.vector(X %*% MLE),
                                             count = data_split[[c]]$full_data_count$count,
                                             prior_means = prior_means,
                                             prior_variances = prior_variances,
                                             C = C))})
  } else if (cv_location == 'hypercube_centre') {
    cv_location <- lapply(1:m, function(c) 'hypercube_centre')
  } else {
    stop("generalised_rho_j_BLR: cv_location must be either \"mode\" or \"hypercube_centre\"")
  }
  transform_matrices <- lapply(1:m, function(c) {
    list('to_Z' = expm::sqrtm(inv_precondition_matrices[[c]]),
         'to_X' = expm::sqrtm(precondition_matrices[[c]]))
  })
  N <- particle_set$N
  # ---------- creating parallel cluster
  if (is.null(cl)) {
    cl <- parallel::makeCluster(n_cores, setup_strategy = "sequential", outfile = "VBF_BLR_outfile.txt")
    parallel::clusterExport(cl, varlist = ls("package:layeredBB"))
    close_cluster <- TRUE
  } else {
    close_cluster <- FALSE
  }
  parallel::clusterExport(cl, envir = environment(), varlist = ls("package:DCFusion"))
  if (!is.null(seed)) {
    parallel::clusterSetRNGStream(cl, iseed = seed)
  }
  max_samples_per_core <- ceiling(N/n_cores)
  split_indices <- split(1:N, ceiling(seq_along(1:N)/max_samples_per_core))
  counts <- c('full_data_count', 'design_count')
  ESS <-  c(particle_set$ESS[1], rep(NA, length(time_mesh)-1))
  CESS <- c(particle_set$CESS[1], rep(NA, length(time_mesh)-1))
  resampled <- rep(FALSE, length(time_mesh))
  L_inv <- construct_L_inv(C = m, d = dim, inv_precondition_matrices = inv_precondition_matrices)
  if (is.null(print_progress_iters)) {
    print_progress_iters <- split_N
  }
  # iterative proposals
  for (j in 2:length(time_mesh)) {
    # ----------- resample particle_set (only resample if ESS < N*ESS_threshold)
    if (particle_set$ESS < N*ESS_threshold) {
      resampled[j-1] <- TRUE
      particle_set <- resample_particle_x_samples(N = N,
                                                  particle_set = particle_set,
                                                  multivariate = TRUE,
                                                  step = j-1,
                                                  resampling_method = resampling_method,
                                                  seed = seed)
    } else {
      resampled[j-1] <- FALSE
    }
    # split the x samples from the previous time marginal (and their means) into approximately equal lists
    split_x_samples <- lapply(split_indices, function(indices) particle_set$x_samples[indices])
    split_x_means <- lapply(split_indices, function(indices) particle_set$x_means[indices,,drop = FALSE])
    V <- construct_V(s = time_mesh[j-1],
                     t = time_mesh[j],
                     end_time = time_mesh[length(time_mesh)],
                     C = m,
                     d = dim,
                     inv_precondition_matrices = inv_precondition_matrices,
                     Lambda = Lambda)
    rho_j_weighted_samples <- parallel::parLapply(cl, X = 1:length(split_indices), fun = function(core) {
      split_N <- length(split_indices[[core]])
      x_mean_j <- matrix(data = NA, nrow = split_N, ncol = dim)
      log_rho_j <- rep(0, split_N)
      x_j <- lapply(1:split_N, function(i) {
        M <- construct_M(s = time_mesh[j-1],
                         t = time_mesh[j],
                         end_time = time_mesh[length(time_mesh)],
                         C = m,
                         d = dim,
                         inv_precondition_matrices = inv_precondition_matrices,
                         Lambda = Lambda,
                         V = V,
                         L_inv = L_inv,
                         sub_posterior_samples = split_x_samples[[core]][[i]],
                         sub_posterior_mean = split_x_means[[core]][i,])$M
        if (j!=length(time_mesh)) {
          return(matrix(mvrnormArma(N = 1, mu = M, Sigma = V), nrow = m, ncol = dim, byrow = TRUE))
        } else {
          return(matrix(mvtnorm::rmvnorm(n = 1, mean = M, sigma = V), nrow = m, ncol = dim, byrow = TRUE))
        }
      })
      cat('Level:', level, '|| Step:', j, '|| Node:', node, '|| Core:', core, '|| START \n',
          file = 'generalised_rho_j_BLR_progress.txt', append = T)
      for (i in 1:split_N) {
        x_mean_j[i,] <- weighted_mean_multivariate(matrix = x_j[[i]],
                                                   weights = inv_precondition_matrices,
                                                   inverse_sum_weights = inverse_sum_matrices(inv_precondition_matrices))
        log_rho_j[i] <- sum(sapply(1:m, function(c) {
          tryCatch(expr = ea_BLR_DL_PT(dim = dim,
                                       x0 = as.vector(split_x_samples[[core]][[i]][c,]),
                                       y = as.vector(x_j[[i]][c,]),
                                       s = time_mesh[j-1],
                                       t = time_mesh[j],
                                       data = data_split[[c]][counts],
                                       prior_means = prior_means,
                                       prior_variances = prior_variances,
                                       C = C,
                                       precondition_mat = precondition_matrices[[c]],
                                       transform_mats = transform_matrices[[c]],
                                       cv_location = cv_location[[c]],
                                       diffusion_estimator = diffusion_estimator,
                                       beta_NB = beta_NB,
                                       gamma_NB_n_points = gamma_NB_n_points,
                                       local_bounds = local_bounds,
                                       logarithm = TRUE),
                   error = function(e) {
                     ea_BLR_DL_PT(dim = dim,
                                  x0 = as.vector(split_x_samples[[core]][[i]][c,]),
                                  y = as.vector(x_j[[i]][c,]),
                                  s = time_mesh[j-1],
                                  t = time_mesh[j],
                                  data = data_split[[c]][counts],
                                  prior_means = prior_means,
                                  prior_variances = prior_variances,
                                  C = C,
                                  precondition_mat = precondition_matrices[[c]],
                                  transform_mats = transform_matrices[[c]],
                                  cv_location = cv_location[[c]],
                                  diffusion_estimator = diffusion_estimator,
                                  beta_NB = beta_NB,
                                  gamma_NB_n_points = gamma_NB_n_points,
                                  local_bounds = FALSE,
                                  logarithm = TRUE)})
        }))
        if (i%%print_progress_iters==0) {
          cat('Level:', level, '|| Step:', j, '/', length(time_mesh),
              '|| Node:', node, '|| Core:', core, '||', i, '/', split_N, '\n',
              file = 'generalised_rho_j_BLR_progress.txt', append = T)
        }
      }
      cat('Level:', level, '|| Step:', j, '|| Node:', node, '|| Core:', core, '||', split_N, '/',
          split_N, '\n', file = 'generalised_rho_j_BLR_progress.txt', append = T)
      return(list('x_j' = x_j, 'x_mean_j' = x_mean_j, 'log_rho_j' = log_rho_j))
    })
    # ---------- update particle set
    # update the weights and return updated particle set
    particle_set$x_samples <- unlist(lapply(1:length(split_indices), function(i) {
      rho_j_weighted_samples[[i]]$x_j}), recursive = FALSE)
    particle_set$x_means <- do.call(rbind, lapply(1:length(split_indices), function(i) {
      rho_j_weighted_samples[[i]]$x_mean_j}))
    # update weight and normalise
    log_rho_j <- unlist(lapply(1:length(split_indices), function(i) {
      rho_j_weighted_samples[[i]]$log_rho_j}))
    norm_weights <- particle_ESS(log_weights = particle_set$log_weights + log_rho_j)
    particle_set$log_weights <- norm_weights$log_weights
    particle_set$normalised_weights <- norm_weights$normalised_weights
    particle_set$ESS <- norm_weights$ESS
    ESS[j] <- particle_set$ESS
    # calculate the conditional ESS (i.e. the 1/sum(inc_change^2))
    # where inc_change is the incremental change in weight (= log_rho_j)
    particle_set$CESS[j] <- particle_ESS(log_weights = log_rho_j)$ESS
    CESS[j] <- particle_set$CESS[j]
  }
  if (close_cluster) {
    parallel::stopCluster(cl)
  }
  # set the y samples as the first element of each of the x_samples
  proposed_samples <- t(sapply(1:N, function(i) particle_set$x_samples[[i]][1,]))
  particle_set$y_samples <- proposed_samples
  # ----------- resample particle_set (only resample if ESS < N*ESS_threshold)
  if (particle_set$ESS < N*ESS_threshold) {
    resampled[particle_set$number_of_steps] <- TRUE
    particle_set <- resample_particle_y_samples(N = N,
                                                particle_set = particle_set,
                                                multivariate = TRUE,
                                                resampling_method = resampling_method,
                                                seed = seed)
  } else {
    resampled[particle_set$number_of_steps] <- FALSE
  }
  return(list('particle_set' = particle_set,
              'proposed_samples' = proposed_samples,
              'ESS' = ESS,
              'CESS' = CESS,
              'resampled' = resampled))
}

#' @export
parallel_generalised_BF_SMC_BLR <- function(particles_to_fuse,
                                            N,
                                            m,
                                            time_mesh,
                                            dim,
                                            data_split,
                                            prior_means,
                                            prior_variances,
                                            C,
                                            precondition_matrices,
                                            resampling_method = 'multi',
                                            ESS_threshold = 0.5,
                                            cv_location = 'hypercube_centre',
                                            diffusion_estimator = 'Poisson',
                                            beta_NB = 10,
                                            gamma_NB_n_points = 2,
                                            local_bounds = TRUE,
                                            seed = NULL,
                                            n_cores = parallel::detectCores(),
                                            cl = NULL,
                                            level = 1,
                                            node = 1,
                                            print_progress_iters = 1000) {
  if (!is.list(particles_to_fuse) | (length(particles_to_fuse)!=m)) {
    stop("parallel_generalised_BF_SMC_BLR: particles_to_fuse must be a list of length m")
  } else if (!all(sapply(particles_to_fuse, function(sub_posterior) ("particle" %in% class(sub_posterior))))) {
    stop("parallel_generalised_BF_SMC_BLR: particles in particles_to_fuse must be \"particle\" objects")
  } else if (!all(sapply(particles_to_fuse, function(sub_posterior) is.matrix(sub_posterior$y_samples)))) {
    stop("parallel_generalised_BF_SMC_BLR: the particles' samples for y should all be matrices")
  } else if (!all(sapply(particles_to_fuse, function(sub_posterior) ncol(sub_posterior$y_samples)==dim))) {
    stop("parallel_generalised_BF_SMC_BLR: the particles' samples for y should all be matrices with dim columns")
  } else if (!is.vector(time_mesh)) {
    stop("parallel_generalised_BF_SMC_BLR: time_mesh must be an ordered vector of length >= 2")
  } else if (length(time_mesh) < 2) {
    stop("parallel_generalised_BF_SMC_BLR: time_mesh must be an ordered vector of length >= 2")
  } else if (!identical(time_mesh, sort(time_mesh))) {
    stop("parallel_generalised_BF_SMC_BLR: time_mesh must be an ordered vector of length >= 2")
  } else if (!is.list(data_split) | length(data_split)!=m) {
    stop("parallel_generalised_BF_SMC_BLR: data_split must be a list of length m")
  } else if (!all(sapply(data_split, function(sub_posterior) (is.list(sub_posterior) & identical(names(sub_posterior), c("y", "X", "full_data_count", "design_count")))))) {
    stop("parallel_generalised_BF_SMC_BLR: each item in data_split must be a list of length 4 with names \'y\', \'X\', \'full_data_count\', \'design_count\'")
  } else if (!all(sapply(1:m, function(i) is.vector(data_split[[i]]$y)))) {
    stop("parallel_generalised_BF_SMC_BLR: for each i in 1:m, data_split[[i]]$y must be a vector")
  } else if (!all(sapply(1:m, function(i) is.matrix(data_split[[i]]$X)))) {
    stop("parallel_generalised_BF_SMC_BLR: for each i in 1:m, data_split[[i]]$X must be a matrix")
  } else if (!all(sapply(1:m, function(i) ncol(data_split[[i]]$X)==dim))) {
    stop("parallel_generalised_BF_SMC_BLR: for each i in 1:m, data_split[[i]]$X must be a matrix with dim columns")
  } else if (!all(sapply(1:m, function(i) length(data_split[[i]]$y)==nrow(data_split[[i]]$X)))) {
    stop("parallel_generalised_BF_SMC_BLR: for each i in 1:m, length(data_split[[i]]$y) and nrow(data_split[[i]]$X) must be equal")
  } else if (!all(sapply(1:m, function(i) is.data.frame(data_split[[i]]$full_data_count)))) {
    stop("parallel_generalised_BF_SMC_BLR: for each i in 1:m, data_split[[i]]$full_data_count must be a data frame")
  } else if (!all(sapply(1:m, function(i) is.data.frame(data_split[[i]]$design_count)))) {
    stop("parallel_generalised_BF_SMC_BLR: for each i in 1:m, data_split[[i]]$design_count must be a data frame")
  } else if (!is.vector(prior_means) | length(prior_means)!=dim) {
    stop("parallel_generalised_BF_SMC_BLR: prior_means must be vectors of length dim")
  } else if (!is.vector(prior_variances) | length(prior_variances)!=dim) {
    stop("parallel_generalised_BF_SMC_BLR: prior_variances must be vectors of length dim")
  } else if (!is.list(precondition_matrices) | (length(precondition_matrices)!=m)) {
    stop("parallel_fusion_SMC_BLR: precondition_matrices must be a list of length m")
  } else if (!(diffusion_estimator %in% c('Poisson', 'NB'))) {
    stop("parallel_generalised_BF_SMC_BLR: diffusion_estimator must be set to either \'Poisson\' or \'NB\'")
  } else if ((ESS_threshold < 0) | (ESS_threshold > 1)) {
    stop("parallel_generalised_BF_SMC_BLR: ESS_threshold must be between 0 and 1")
  } else if ((cv_location != 'mode') & (cv_location != 'hypercube_centre')) {
    stop("parallel_generalised_BF_SMC_BLR: cv_location must be either \"mode\" or \"hypercube_centre\"")
  } else if (!any(class(cl)=="cluster") & !is.null(cl)) {
    stop("parallel_generalised_BF_SMC_BLR: cl must be a \"cluster\" object or NULL")
  }
  # set a seed if one is supplied
  if (!is.null(seed)) {
    set.seed(seed)
  }
  # start time recording
  pcm <- proc.time()
  # ---------- first importance sampling step
  # pre-calculating the inverse precondition matrices
  inv_precondition_matrices <- lapply(precondition_matrices, solve)
  Lambda <- inverse_sum_matrices(inv_precondition_matrices)
  particles <- rho_IS_multivariate(particles_to_fuse = particles_to_fuse,
                                   dim = dim,
                                   N = N,
                                   m = m,
                                   time = time_mesh[length(time_mesh)],
                                   inv_precondition_matrices = inv_precondition_matrices,
                                   inverse_sum_inv_precondition_matrices = Lambda,
                                   number_of_steps = length(time_mesh),
                                   time_mesh = time_mesh,
                                   resampling_method = resampling_method,
                                   n_cores = n_cores,
                                   cl = cl)
  # ---------- iterative steps
  rho_j <- generalised_rho_j_BLR(particle_set = particles,
                                 m = m,
                                 time_mesh = time_mesh,
                                 dim = dim,
                                 data_split = data_split,
                                 prior_means = prior_means,
                                 prior_variances = prior_variances,
                                 C = C,
                                 precondition_matrices = precondition_matrices,
                                 inv_precondition_matrices = inv_precondition_matrices,
                                 Lambda = Lambda,
                                 resampling_method = resampling_method,
                                 ESS_threshold = ESS_threshold,
                                 cv_location = cv_location,
                                 diffusion_estimator = diffusion_estimator,
                                 beta_NB = beta_NB,
                                 gamma_NB_n_points = gamma_NB_n_points,
                                 local_bounds = local_bounds,
                                 seed = seed,
                                 n_cores = n_cores,
                                 cl = cl,
                                 level = level,
                                 node = node,
                                 print_progress_iters = print_progress_iters)
  if (identical(precondition_matrices, rep(list(diag(1, dim)), m))) {
    new_precondition_matrices <- list(diag(1, dim), precondition_matrices)
  } else {
    new_precondition_matrices <- list(inverse_sum_matrices(inv_precondition_matrices),
                                      precondition_matrices)
  }
  return(list('particles' = rho_j$particle_set,
              'proposed_samples' = rho_j$proposed_samples,
              'time' = (proc.time()-pcm)['elapsed'],
              'ESS' = rho_j$ESS,
              'CESS' = rho_j$CESS,
              'resampled' = rho_j$resampled,
              'precondition_matrices' = new_precondition_matrices,
              'combined_data' = combine_data(list_of_data = data_split, dim = dim)))  
}
