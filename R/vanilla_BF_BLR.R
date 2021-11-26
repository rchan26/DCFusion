#' @export
vanilla_rho_j_BLR <- function(particle_set,
                     m,
                     time,
                     dim,
                     data_split,
                     prior_means,
                     prior_variances,
                     C,
                     proposal_cov,
                     precondition_matrices,
                     inv_precondition_matrices,
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
    stop("Q_IS_BLR: particle_set must be a \"particle\" object")
  } else if (!is.list(data_split) | length(data_split)!=m) {
    stop("Q_IS_BLR: data_split must be a list of length m")
  } else if (!all(sapply(data_split, function(sub_posterior) (is.list(sub_posterior) & identical(names(sub_posterior), c("y", "X", "full_data_count", "design_count")))))) {
    stop("Q_IS_BLR: each item in data_split must be a list of length 4 with names \'y\', \'X\', \'full_data_count\', \'design_count\'")
  } else if (!all(sapply(1:m, function(i) is.vector(data_split[[i]]$y)))) {
    stop("Q_IS_BLR: for each i in 1:m, data_split[[i]]$y must be a vector")
  } else if (!all(sapply(1:m, function(i) is.matrix(data_split[[i]]$X)))) {
    stop("Q_IS_BLR: for each i in 1:m, data_split[[i]]$X must be a matrix")
  } else if (!all(sapply(1:m, function(i) ncol(data_split[[i]]$X)==dim))) {
    stop("Q_IS_BLR: for each i in 1:m, ncol(data_split[[i]]$X) must be equal to dim")
  } else if (!all(sapply(1:m, function(i) length(data_split[[i]]$y)==nrow(data_split[[i]]$X)))) {
    stop("Q_IS_BLR: for each i in 1:m, length(data_split[[i]]$y) and nrow(data_split[[i]]$X) must be equal")
  } else if (!all(sapply(1:m, function(i) is.data.frame(data_split[[i]]$full_data_count)))) {
    stop("Q_IS_BLR: for each i in 1:m, data_split[[i]]$full_data_count must be a data frame")
  } else if (!all(sapply(1:m, function(i) is.data.frame(data_split[[i]]$design_count)))) {
    stop("Q_IS_BLR: for each i in 1:m, data_split[[i]]$design_count must be a data frame")
  } else if (!is.vector(prior_means) | length(prior_means)!=dim) {
    stop("Q_IS_BLR: prior_means must be vectors of length dim")
  } else if (!is.vector(prior_variances) | length(prior_variances)!=dim) {
    stop("Q_IS_BLR: prior_variances must be vectors of length dim")
  } else if (!is.list(precondition_matrices) | (length(precondition_matrices)!=m)) {
    stop("Q_IS_BLR: precondition_matrices must be a list of length m")
  } else if (!is.list(inv_precondition_matrices) | (length(inv_precondition_matrices)!=m)) {
    stop("Q_IS_BLR: inv_precondition_matrices must be a list of length m")
  } else if (!(diffusion_estimator %in% c('Poisson', 'NB'))) {
    stop("Q_IS_BLR: diffusion_estimator must be set to either \'Poisson\' or \'NB\'")
  } else if (!any(class(cl)=="cluster") & !is.null(cl)) {
    stop("Q_IS_BLR: cl must be a \"cluster\" object or NULL")
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
    stop("Q_IS_BLR: cv_location must be either \"mode\" or \"hypercube_centre\"")
  }
  transform_matrices <- lapply(1:m, function(c) {
    list('to_Z' = expm::sqrtm(inv_precondition_matrices[[c]]),
         'to_X' = expm::sqrtm(precondition_matrices[[c]]))
  })
  N <- particle_set$N
  # ---------- creating parallel cluster
  if (is.null(cl)) {
    cl <- parallel::makeCluster(n_cores, setup_strategy = "sequential", outfile = "SMC_BLR_outfile.txt")
    parallel::clusterExport(cl, varlist = ls("package:layeredBB"))
    close_cluster <- TRUE
  } else {
    close_cluster <- FALSE
  }
  parallel::clusterExport(cl, envir = environment(), 
                          varlist = c(ls(), "ea_phi_BLR_DL_matrix",
                                      "ea_phi_BLR_DL_bounds",
                                      "ea_BLR_DL_PT"))
  if (!is.null(seed)) {
    parallel::clusterSetRNGStream(cl, iseed = seed)
  }
  # split the x samples and their means into approximately equal lists
  max_samples_per_core <- ceiling(N/n_cores)
  split_indices <- split(1:N, ceiling(seq_along(1:N)/max_samples_per_core))
  split_x_samples <- lapply(split_indices, function(indices) particle_set$x_samples[indices])
  split_x_means <- lapply(split_indices, function(indices) particle_set$x_means[indices,,drop = FALSE])
  counts <- c('full_data_count', 'design_count')
  # for each set of x samples, we propose a new value y and assign a weight for it
  # sample for y and importance weight in parallel to split computation
  Q_weighted_samples <- parallel::parLapply(cl, X = 1:length(split_indices), fun = function(core) {
    split_N <- length(split_indices[[core]])
    y_samples <- t(apply(split_x_means[[core]], 1, function(vec) mvrnormArma(N = 1, mu = vec, Sigma = proposal_cov)))
    log_Q_weights <- rep(0, split_N)
    cat('Level:', level, '|| Node:', node, '|| Core:', core, '|| START \n',
        file = 'Q_IS_BLR_progress.txt', append = T)
    if (is.null(print_progress_iters)) {
      print_progress_iters <- split_N
    }
    for (i in 1:split_N) {
      log_Q_weights[i] <- sum(sapply(1:m, function(c) {
        tryCatch(expr = ea_BLR_DL_PT(dim = dim,
                                     x0 = as.vector(split_x_samples[[core]][[i]][c,]),
                                     y = as.vector(y_samples[i,]),
                                     s = 0,
                                     t = time,
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
                                y = as.vector(y_samples[i,]),
                                s = 0,
                                t = time,
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
        cat('Level:', level, '|| Node:', node, '|| Core:', core, '||', i, '/',
            split_N, '\n', file = 'Q_IS_BLR_progress.txt', append = T)
      }
    }
    cat('Completed: Level:', level, '|| Node:', node, '|| Core:', core, '||', split_N, '/',
        split_N, '\n', file = 'Q_IS_BLR_progress.txt', append = T)
    return(list('y_samples' = y_samples, 'log_Q_weights' = log_Q_weights))
  })
  if (close_cluster) {
    parallel::stopCluster(cl)
  }
  # unlist the proposed samples for y and their associated log Q weights
  y_samples <- do.call(rbind, lapply(1:length(split_x_samples), function(i) {
    Q_weighted_samples[[i]]$y_samples}))
  log_Q_weights <- unlist(lapply(1:length(split_x_samples), function(i) {
    Q_weighted_samples[[i]]$log_Q_weights}))
  # ---------- update particle set
  # update the weights and return updated particle set
  particle_set$y_samples <- y_samples
  # normalise weight
  norm_weights <- particle_ESS(log_weights = particle_set$log_weights + log_Q_weights)
  particle_set$log_weights <- norm_weights$log_weights
  particle_set$normalised_weights <- norm_weights$normalised_weights
  particle_set$ESS <- norm_weights$ESS
  # calculate the conditional ESS (i.e. the 1/sum(inc_change^2))
  # where inc_change is the incremental change in weight (= log_Q_weights)
  particle_set$CESS['Q'] <- particle_ESS(log_weights = log_Q_weights)$ESS
  # set the resampled indicator to FALSE
  particle_set$resampled['Q'] <- FALSE
  return(particle_set)
}
