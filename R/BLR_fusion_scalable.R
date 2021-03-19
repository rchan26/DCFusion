#' @export
ea_phi_BLR_DL_scalable <- function(cv_list,
                                   beta,
                                   y_labels,
                                   X,
                                   prior_means,
                                   prior_variances,
                                   C,
                                   precondition_mat) {
  if (is.vector(beta)) {
    return(ea_phi_BLR_DL_vec_scalable(cv_list = cv_list,
                                      beta = beta,
                                      y_labels = y_labels,
                                      X = X,
                                      prior_means = prior_means,
                                      prior_variances = prior_variances,
                                      C = C,
                                      precondition_mat = precondition_mat))
  } else if (is.matrix(beta)) {
    return(ea_phi_BLR_DL_matrix_scalable(cv_list = cv_list,
                                         beta = beta,
                                         y_labels = y_labels,
                                         X = X,
                                         prior_means = prior_means,
                                         prior_variances = prior_variances,
                                         C = C,
                                         precondition_mat = precondition_mat))
  }
  stop("ea_phi_BLR_DL_scalable: beta must be a vector or a matrix")
}

#' @export
control_variates_BLR <- function(dim,
                                 data,
                                 prior_means,
                                 prior_variances,
                                 C,
                                 precondition_mat) {
  beta_hat <- as.vector(unname(glm(formula = data$y ~ data$X[,2:dim], family = 'binomial')$coeff))
  X_beta <- as.vector(data$X %*% beta_hat)
  grad_log_beta_hat <- as.vector(log_BLR_gradient(beta = beta_hat,
                                                  y_labels = data$y,
                                                  X = data$X,
                                                  X_beta = X_beta,
                                                  prior_means = prior_means,
                                                  prior_variances = prior_variances,
                                                  C = C))
  t1 <- as.vector(t(grad_log_beta_hat)%*%precondition_mat%*%grad_log_beta_hat)
  hessian_log_beta_hat <- log_BLR_hessian(X = data$X,
                                          X_beta = X_beta,
                                          prior_variances = prior_variances,
                                          C = C)
  t2 <- sum(precondition_mat * hessian_log_beta_hat)
  return(list('beta_hat' = beta_hat,
              'grad_log_beta_hat' = grad_log_beta_hat,
              'hessian_log_beta_hat' = hessian_log_beta_hat,
              't1' = t1,
              't2' = t2,
              'constant' = 0.5*(t1+t2),
              'data_size' = length(data$y)))
}

#' @export
maximum_distance_hypercube_to_cv <- function(dim,
                                             beta_hat,
                                             transform_mat,
                                             bessel_layers) {
  bounds <- lapply(1:dim, function(d) c(bessel_layers[[d]]$L, bessel_layers[[d]]$U))
  hypercube_Z <- as.matrix(expand.grid(bounds))
  hypercube_X <- t(transform_mat %*% t(hypercube_Z))
  return(max(sapply(1:nrow(hypercube_X), function(i) Euclidean_distance(hypercube_X[i,], beta_hat))))
}

#' @export
ea_phi_BLR_DL_bounds_scalable <- function(cv_list,
                                          dim,
                                          X,
                                          prior_variances,
                                          C,
                                          precondition_mat,
                                          transform_mat,
                                          bessel_layers) {
  dist <- maximum_distance_hypercube_to_cv(dim = dim,
                                           beta_hat = cv_list$beta_hat,
                                           transform_mat = transform_mat,
                                           bessel_layers = bessel_layers)
  P_n <- hessian_bound_BLR(dim = dim,
                           X = X,
                           prior_variances = prior_variances,
                           C = C,
                           precondition_mat = precondition_mat)
  n <- cv_list$data_size
  alpha_bds <- sqrt(sum((2*cv_list$grad_log_beta_hat)^2))
  sum_precond_diag <- sum(diag(precondition_mat))
  bound <- n*P_n*dist*(alpha_bds+n*P_n*dist) + n*P_n*sum_precond_diag
  return(list('LB' = -0.5*bound + cv_list$constant,
              'UB' = 0.5*bound + cv_list$constant,
              'dist' = dist,
              'P_n' = P_n,
              'n' = n,
              'alpha_bds' = alpha_bds,
              'sum_precond_diag' = sum_precond_diag))
}

#' @export
ea_BLR_DL_PT_scalable <- function(dim,
                                  x0,
                                  y,
                                  s,
                                  t,
                                  cv_list,
                                  y_labels,
                                  X,
                                  prior_means,
                                  prior_variances,
                                  C,
                                  precondition_mat,
                                  transform_mats,
                                  diffusion_estimator,
                                  beta_NB = 10,
                                  logarithm) {
  # transform to preconditoned space
  z0 <- transform_mats$to_Z %*% x0
  zt <- transform_mats$to_Z %*% y
  # simulate layer information
  bes_layers <- layeredBB::multi_bessel_layer_simulation(dim = dim,
                                                         x = z0,
                                                         y = zt,
                                                         s = s,
                                                         t = t,
                                                         mult = 0.1)
  bounds <- ea_phi_BLR_DL_bounds_scalable(cv_list = cv_list,
                                          dim = dim,
                                          X = X,
                                          prior_variances = prior_variances,
                                          C = C,
                                          precondition_mat = precondition_mat,
                                          transform_mat = transform_mats$to_X,
                                          bessel_layers = bes_layers)
  LX <- bounds$LB
  UX <- bounds$UB
  if (diffusion_estimator=='Poisson') {
    # simulate the number of points to simulate from Possion distribution
    kap <- rpois(n = 1, lambda = (UX-LX)*(t-s))
    log_acc_prob <- 0
    if (kap > 0) {
      layered_bb <- layeredBB::multi_layered_brownian_bridge(dim = dim,
                                                             x = z0,
                                                             y = zt,
                                                             s = s,
                                                             t = t,
                                                             bessel_layers = bes_layers,
                                                             times = runif(kap, s, t))
      sim_path <- t(transform_mats$to_X %*% layered_bb$simulated_path[1:dim,])
      phi <- ea_phi_BLR_DL_matrix_scalable(cv_list = cv_list,
                                           beta = sim_path,
                                           y_labels = y_labels,
                                           X = X,
                                           prior_means = prior_means,
                                           prior_variances = prior_variances,
                                           C = C,
                                           precondition_mat = precondition_mat)
      terms <- (UX-phi$phi)
      log_acc_prob <- sum(log(terms))
      if (any(terms < 0)) {
        cat('########## \n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        cat('LX:', LX, '\n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        cat('UX:', UX, '\n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        cat('kap:', kap, '\n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        cat('phi:', phi$phi, '\n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        cat('(UX-phi):', terms, '\n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        cat('(phi-LX):', phi$phi-LX, '\n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        cat('##### Some of (UX-phi) are < 0. \n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        cat('phi[which((UX-phi) < 0)]:', phi$phi[which(terms<0)], '\n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        cat('which((UX-phi) < 0):', which(terms < 0), '\n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        cat('phi[which((phi-LX) < 0)]:', phi$phi[which(phi$phi-LX<0)], '\n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        cat('which((phi-LX) < 0):', which((phi$phi-LX) < 0), '\n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        hes_spectral_norm_I <- sapply(1:kap, function(i) spectral_norm_hessian(dim = dim,
                                                                               beta = sim_path[i,],
                                                                               X = X,
                                                                               index = phi$I[i],
                                                                               prior_variances = prior_variances,
                                                                               C = C,
                                                                               precondition_mat = precondition_mat))
        hes_spectral_norm_J <- sapply(1:kap, function(i) spectral_norm_hessian(dim = dim,
                                                                               beta = sim_path[i,],
                                                                               X = X,
                                                                               index = phi$J[i],
                                                                               prior_variances = prior_variances,
                                                                               C = C,
                                                                               precondition_mat = precondition_mat))
        cat('bounds$P_n:', bounds$P_n, '\n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        cat('hes_spectral_norm_I:', hes_spectral_norm_I, '\n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        cat('hes_spectral_norm_J:', hes_spectral_norm_J, '\n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        if (any(hes_spectral_norm_I > bounds$P_n)) {
          cat('##### Some of hes_spectral_norm_I < bounds$P_n \n', file = "SMC_BLR_scalable_bounds.txt", append = T)
          cat('hes_spectral_norm_I[which(hes_spectral_norm_I > bounds$P_n)]:', hes_spectral_norm_I[which(hes_spectral_norm_I > bounds$P_n)], '\n', file = "SMC_BLR_scalable_bounds.txt", append = T)
          cat('which(hes_spectral_norm_I > bounds$P_n):', which(hes_spectral_norm_I > bounds$P_n),  '\n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        } else {
          cat('##### hes_spectral_norm_I < bounds$P_n \n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        }
        if (any(hes_spectral_norm_J > bounds$P_n)) {
          cat('##### Some of hes_spectral_norm_J > bounds$P_n \n', file = "SMC_BLR_scalable_bounds.txt", append = T)
          cat('hes_spectral_norm_J[which(hes_spectral_norm_I > bounds$P_n)]:', hes_spectral_norm_J[which(hes_spectral_norm_J > bounds$P_n)], '\n', file = "SMC_BLR_scalable_bounds.txt", append = T)
          cat('which(hes_spectral_norm_J > bounds$P_n):', which(hes_spectral_norm_J > bounds$P_n),  '\n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        } else {
          cat('##### hes_spectral_norm_J < bounds$P_n \n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        }
        stop('Some of (UX-phi) are < 0.')
      } else if (any((phi$phi - LX) < 0)) {
        cat('########## \n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        cat('LX:', LX, '\n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        cat('UX:', UX, '\n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        cat('kap:', kap, '\n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        cat('phi:', phi$phi, '\n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        cat('(UX-phi):', terms, '\n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        cat('(phi-LX):', phi$phi-LX, '\n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        cat('##### Some of (phi-LX) are < 0. \n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        cat('phi[which((UX-phi) < 0)]:', phi$phi[which(terms<0)], '\n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        cat('which((UX-phi) < 0):', which(terms < 0), '\n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        cat('phi[which((phi-LX) < 0)]:', phi$phi[which(phi$phi-LX<0)], '\n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        cat('which((phi-LX) < 0):', which((phi$phi-LX) < 0), '\n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        hes_spectral_norm_I <- sapply(1:kap, function(i) spectral_norm_hessian(dim = dim,
                                                                               beta = sim_path[i,],
                                                                               X = X,
                                                                               index = phi$I[i],
                                                                               prior_variances = prior_variances,
                                                                               C = C,
                                                                               precondition_mat = precondition_mat))
        hes_spectral_norm_J <- sapply(1:kap, function(i) spectral_norm_hessian(dim = dim,
                                                                               beta = sim_path[i,],
                                                                               X = X,
                                                                               index = phi$J[i],
                                                                               prior_variances = prior_variances,
                                                                               C = C,
                                                                               precondition_mat = precondition_mat))
        cat('bounds$P_n:', bounds$P_n, '\n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        cat('hes_spectral_norm_I:', hes_spectral_norm_I, '\n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        cat('hes_spectral_norm_J:', hes_spectral_norm_J, '\n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        if (any(hes_spectral_norm_I > bounds$P_n)) {
          cat('##### Some of hes_spectral_norm_I < bounds$P_n \n', file = "SMC_BLR_scalable_bounds.txt", append = T)
          cat('hes_spectral_norm_I[which(hes_spectral_norm_I > bounds$P_n)]:', hes_spectral_norm_I[which(hes_spectral_norm_I > bounds$P_n)], '\n', file = "SMC_BLR_scalable_bounds.txt", append = T)
          cat('which(hes_spectral_norm_I > bounds$P_n):', which(hes_spectral_norm_I > bounds$P_n),  '\n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        } else {
          cat('##### hes_spectral_norm_I < bounds$P_n \n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        }
        if (any(hes_spectral_norm_J > bounds$P_n)) {
          cat('##### Some of hes_spectral_norm_J > bounds$P_n \n', file = "SMC_BLR_scalable_bounds.txt", append = T)
          cat('hes_spectral_norm_J[which(hes_spectral_norm_I > bounds$P_n)]:', hes_spectral_norm_J[which(hes_spectral_norm_J > bounds$P_n)], '\n', file = "SMC_BLR_scalable_bounds.txt", append = T)
          cat('which(hes_spectral_norm_J > bounds$P_n):', which(hes_spectral_norm_J > bounds$P_n),  '\n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        } else {
          cat('##### hes_spectral_norm_J < bounds$P_n \n', file = "SMC_BLR_scalable_bounds.txt", append = T)
        }
        stop('Some of (phi-LX) are < 0.')
      }
    }
    if (logarithm) {
      return(-LX*(t-s) - kap*log(UX-LX) + log_acc_prob)
    } else {
      return(exp(-LX*(t-s) - kap*log(UX-LX) + log_acc_prob))
    }
  } else if (diffusion_estimator=="NB") {
    stop("ea_BLR_DL_PT_scalable: NB diffusion_estimator not yet implemented")
  } else {
    stop("ea_BLR_DL_PT_scalable: diffusion_estimator must be set to either \'Poisson\' or \'NB\'")
  }
}

#' @export
Q_IS_BLR_scalable <- function(particle_set,
                              m,
                              time,
                              dim,
                              data_split,
                              prior_means,
                              prior_variances,
                              C,
                              precondition_matrices,
                              inv_precondition_matrices,
                              diffusion_estimator,
                              beta_NB = 10,
                              seed = NULL,
                              n_cores = parallel::detectCores(),
                              level = 1,
                              node = 1) {
  if (!("particle" %in% class(particle_set))) {
    stop("Q_IS_BLR_scalable: particle_set must be a \"particle\" object")
  } else if (!is.list(data_split) | length(data_split)!=m) {
    stop("Q_IS_BLR_scalable: data_split must be a list of length m")
  } else if (!all(sapply(data_split, function(sub_posterior) (is.list(sub_posterior) & identical(names(sub_posterior), c("y", "X")))))) {
    stop("Q_IS_BLR_scalable: each item in data_split must be a list of length 2 with names y and X")
  } else if (!all(sapply(1:m, function(i) is.vector(data_split[[i]]$y)))) {
    stop("Q_IS_BLR_scalable: for each i in 1:m, data_split[[i]]$y must be a vector")
  } else if (!all(sapply(1:m, function(i) is.matrix(data_split[[i]]$X)))) {
    stop("Q_IS_BLR_scalable: for each i in 1:m, data_split[[i]]$X must be a matrix")
  } else if (!all(sapply(1:m, function(i) ncol(data_split[[i]]$X)==dim))) {
    stop("Q_IS_BLR_scalable: for each i in 1:m, ncol(data_split[[i]]$X) must be equal to dim")
  } else if (!all(sapply(1:m, function(i) length(data_split[[i]]$y)==nrow(data_split[[i]]$X)))) {
    stop("Q_IS_BLR_scalable: for each i in 1:m, length(data_split[[i]]$y) and nrow(data_split[[i]]$X) must be equal")
  } else if (!is.vector(prior_means) | length(prior_means)!=dim) {
    stop("Q_IS_BLR_scalable: prior_means must be vectors of length dim")
  } else if (!is.vector(prior_variances) | length(prior_variances)!=dim) {
    stop("Q_IS_BLR_scalable: prior_variances must be vectors of length dim")
  } else if (!is.list(precondition_matrices) | (length(precondition_matrices)!=m)) {
    stop("Q_IS_BLR_scalable: precondition_matrices must be a list of length m")
  } else if (!is.list(inv_precondition_matrices) | (length(inv_precondition_matrices)!=m)) {
    stop("Q_IS_BLR_scalable: inv_precondition_matrices must be a list of length m")
  } else if (!(diffusion_estimator %in% c('Poisson', 'NB'))) {
    stop("Q_IS_BLR_scalable: diffusion_estimator must be set to either \'Poisson\' or \'NB\'")
  }
  transform_matrices <- lapply(1:m, function(c) {
    list('to_Z' = expm::sqrtm(inv_precondition_matrices[[c]]),
         'to_X' = expm::sqrtm(precondition_matrices[[c]]))
  })
  cv_lists <- lapply(1:m, function(c) {
    control_variates_BLR(dim = dim,
                         data = data_split[[c]],
                         prior_means = prior_means,
                         prior_variances = prior_variances,
                         C = C,
                         precondition_mat = precondition_matrices[[c]])
  })
  proposal_cov <- calculate_proposal_cov(time = time, weights = inv_precondition_matrices)
  N <- particle_set$N
  # ---------- creating parallel cluster
  cl <- parallel::makeCluster(n_cores, setup_strategy = "sequential", outfile = "SMC_BLR_outfile.txt")
  varlist <- c(ls(), list("ea_phi_BLR_DL_matrix_scalable",
                          "ea_phi_BLR_DL_bounds_scalable",
                          "ea_BLR_DL_PT_scalable"))
  parallel::clusterExport(cl, envir = environment(), varlist = varlist)
  # exporting functions from layeredBB package to simulate layered Brownian bridges
  parallel::clusterExport(cl, varlist = ls("package:layeredBB"))
  if (!is.null(seed)) {
    parallel::clusterSetRNGStream(cl, iseed = seed)
  }
  # split the x samples and their means into approximately equal lists
  max_samples_per_core <- ceiling(N/n_cores)
  split_indices <- split(1:N, ceiling(seq_along(1:N)/max_samples_per_core))
  split_x_samples <- lapply(split_indices, function(indices) particle_set$x_samples[indices])
  split_x_means <- lapply(split_indices, function(indices) particle_set$x_means[indices,,drop = FALSE])
  # for each set of x samples, we propose a new value y and assign a weight for it
  # sample for y and importance weight in parallel to split computation
  Q_weighted_samples <- parallel::parLapply(cl, X = 1:length(split_indices), fun = function(core) {
    split_N <- length(split_indices[[core]])
    y_samples <- matrix(nrow = split_N, ncol = dim)
    log_Q_weights <- rep(0, split_N)
    cat('Level:', level, '|| Node:', node, '|| Core:', core, '|| START \n',
        file = 'Q_IS_BLR_scalable_progress.txt', append = T)
    for (i in 1:split_N) {
      y_samples[i,] <- mvrnormArma(N = 1, mu = split_x_means[[core]][i,], Sigma = proposal_cov)
      log_Q_weights[i] <- sum(sapply(1:m, function(c) {
        ea_BLR_DL_PT_scalable(dim = dim,
                              x0 = as.vector(split_x_samples[[core]][[i]][c,]),
                              y = as.vector(y_samples[i,]),
                              s = 0,
                              t = time,
                              cv_list = cv_lists[[c]],
                              y_labels = data_split[[c]]$y,
                              X = data_split[[c]]$X,
                              prior_means = prior_means,
                              prior_variances = prior_variances,
                              C = C,
                              precondition_mat = precondition_matrices[[c]],
                              transform_mats = transform_matrices[[c]],
                              diffusion_estimator = diffusion_estimator,
                              beta_NB = beta_NB,
                              logarithm = TRUE)
      }))
      cat('Level:', level, '|| Node:', node, '|| Core:', core, '||', i, '/',
          split_N, '\n', file = 'Q_IS_BLR_scalable_progress.txt', append = T)
    }
    return(list('y_samples' = y_samples, 'log_Q_weights' = log_Q_weights))
  })
  parallel::stopCluster(cl)
  # unlist the proposed samples for y and their associated log Q weights
  y_samples <- do.call(rbind, lapply(1:length(split_x_samples), function(i) {
    Q_weighted_samples[[i]]$y_samples}))
  log_Q_weights <- unlist(lapply(1:length(split_x_samples), function(i) {
    Q_weighted_samples[[i]]$log_Q_weights}))
  # ---------- update particle set
  # update the weights and return updated particle set
  particle_set$y_samples <- y_samples
  particle_set$log_weights <- particle_set$log_weights + log_Q_weights
  # normalise weights
  norm_weights <- particle_ESS(log_weights = particle_set$log_weights)
  particle_set$normalised_weights <- norm_weights$normalised_weights
  particle_set$ESS <- norm_weights$ESS
  # calculate the conditional ESS (i.e. the 1/sum(inc_change^2))
  # where inc_change is the incremental change in weight (= log_Q_weights)
  particle_set$CESS['Q'] <- particle_ESS(log_weights = log_Q_weights)$ESS
  # set the resampled indicator to FALSE
  particle_set$resampled['Q'] <- FALSE
  return(particle_set)
}

#' @export
parallel_fusion_SMC_BLR_scalable <- function(particles_to_fuse,
                                             N,
                                             m,
                                             time,
                                             dim,
                                             data_split,
                                             prior_means,
                                             prior_variances,
                                             C,
                                             precondition_matrices,
                                             resampling_method = 'multi',
                                             ESS_threshold = 0.5,
                                             diffusion_estimator = 'Poisson',
                                             beta_NB = 10,
                                             seed = NULL,
                                             n_cores = parallel::detectCores(),
                                             level = 1,
                                             node = 1) {
  if (!is.list(particles_to_fuse) | (length(particles_to_fuse)!=m)) {
    stop("parallel_fusion_SMC_BLR_scalable: particles_to_fuse must be a list of length m")
  } else if (!all(sapply(particles_to_fuse, function(sub_posterior) ("particle" %in% class(sub_posterior))))) {
    stop("parallel_fusion_SMC_BLR_scalable: particles in particles_to_fuse must be \"particle\" objects")
  } else if (!all(sapply(particles_to_fuse, function(sub_posterior) is.matrix(sub_posterior$y_samples)))) {
    stop("parallel_fusion_SMC_BLR_scalable: the particles' samples for y should all be matrices")
  } else if (!all(sapply(particles_to_fuse, function(sub_posterior) ncol(sub_posterior$y_samples)==dim))) {
    stop("parallel_fusion_SMC_BLR_scalable: the particles' samples for y should all be matrices with dim columns")
  } else if (!is.list(data_split) | length(data_split)!=m) {
    stop("parallel_fusion_SMC_BLR_scalable: data_split must be a list of length m")
  } else if (!all(sapply(data_split, function(sub_posterior) (is.list(sub_posterior) & identical(names(sub_posterior), c("y", "X")))))) {
    stop("parallel_fusion_SMC_BLR_scalable: each item in data_split must be a list of length 2 with names y and X")
  } else if (!all(sapply(1:m, function(i) is.vector(data_split[[i]]$y)))) {
    stop("parallel_fusion_SMC_BLR_scalable: for each i in 1:m, data_split[[i]]$y must be a vector")
  } else if (!all(sapply(1:m, function(i) is.matrix(data_split[[i]]$X)))) {
    stop("parallel_fusion_SMC_BLR_scalable: for each i in 1:m, data_split[[i]]$X must be a matrix")
  } else if (!all(sapply(1:m, function(i) ncol(data_split[[i]]$X)==dim))) {
    stop("parallel_fusion_SMC_BLR_scalable: for each i in 1:m, data_split[[i]]$X must be a matrix with dim columns")
  } else if (!all(sapply(1:m, function(i) length(data_split[[i]]$y)==nrow(data_split[[i]]$X)))) {
    stop("parallel_fusion_SMC_BLR_scalable: for each i in 1:m, length(data_split[[i]]$y) and nrow(data_split[[i]]$X) must be equal")
  } else if (!is.vector(prior_means) | length(prior_means)!=dim) {
    stop("parallel_fusion_SMC_BLR_scalable: prior_means must be vectors of length dim")
  } else if (!is.vector(prior_variances) | length(prior_variances)!=dim) {
    stop("parallel_fusion_SMC_BLR_scalable: prior_variances must be vectors of length dim")
  } else if (!is.list(precondition_matrices) | (length(precondition_matrices)!=m)) {
    stop("parallel_fusion_SMC_BLR_scalable: precondition_matrices must be a list of length m")
  } else if (!(diffusion_estimator %in% c('Poisson', 'NB'))) {
    stop("parallel_fusion_SMC_BLR_scalable: diffusion_estimator must be set to either \'Poisson\' or \'NB\'")
  } else if ((ESS_threshold < 0) | (ESS_threshold > 1)) {
    stop("parallel_fusion_SMC_BLR_scalable: ESS_threshold must be between 0 and 1")
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
  # pre-calculating the inverse precondition matrices
  inv_precondition_matrices <- lapply(precondition_matrices, solve)
  # importance sampling for rho step
  particles <- rho_IS_multivariate(particles_to_fuse = particles_to_fuse,
                                   dim = dim,
                                   N = N,
                                   m = m,
                                   time = time,
                                   inv_precondition_matrices = inv_precondition_matrices,
                                   inverse_sum_inv_precondition_matrices = inverse_sum_matrices(inv_precondition_matrices),
                                   n_cores = n_cores)
  # record ESS and CESS after rho step
  ESS <- c('rho' = particles$ESS)
  CESS <- c('rho' = particles$CESS['rho'])
  # ----------- resample particles
  # only resample if ESS < N*ESS_threshold
  if (particles$ESS < N*ESS_threshold) {
    resampled <- c('rho' = TRUE)
    particles <- resample_particle_x_samples(N = N,
                                             particle_set = particles,
                                             multivariate = TRUE,
                                             resampling_method = resampling_method,
                                             seed = seed)
  } else {
    resampled <- c('rho' = FALSE)
  }
  # ---------- second importance sampling step
  # unbiased estimator for Q
  particles <- Q_IS_BLR_scalable(particle_set = particles,
                                 m = m,
                                 time = time,
                                 dim = dim,
                                 data_split = data_split,
                                 prior_means = prior_means,
                                 prior_variances = prior_variances,
                                 C = C,
                                 precondition_matrices = precondition_matrices,
                                 inv_precondition_matrices = inv_precondition_matrices,
                                 diffusion_estimator = diffusion_estimator,
                                 beta_NB = beta_NB,
                                 seed = seed,
                                 n_cores = n_cores,
                                 level = level,
                                 node = node)
  # record ESS and CESS after Q step
  ESS['Q'] <- particles$ESS
  CESS['Q'] <- particles$CESS['Q']
  names(CESS) <- c('rho', 'Q')
  # record proposed samples
  proposed_samples <- particles$y_samples
  # ----------- resample particles
  # only resample if ESS < N*ESS_threshold
  if (particles$ESS < N*ESS_threshold) {
    resampled['Q'] <- TRUE
    particles <- resample_particle_y_samples(N = N,
                                             particle_set = particles,
                                             multivariate = TRUE,
                                             resampling_method = resampling_method,
                                             seed = seed)
  } else {
    resampled['Q'] <- FALSE
  }
  if (identical(precondition_matrices, rep(list(diag(1, dim)), m))) {
    return(list('particles' = particles,
                'proposed_samples' = proposed_samples,
                'time' = (proc.time()-pcm)['elapsed'],
                'ESS' = ESS,
                'CESS' = CESS,
                'resampled' = resampled,
                'precondition_matrices' = list(diag(1, dim), precondition_matrices),
                'combined_data' = combine_data(list_of_data = data_split, dim = dim)))
  } else {
    return(list('particles' = particles,
                'proposed_samples' = proposed_samples,
                'time' = (proc.time()-pcm)['elapsed'],
                'ESS' = ESS,
                'CESS' = CESS,
                'resampled' = resampled,
                'precondition_matrices' = list(inverse_sum_matrices(inv_precondition_matrices),
                                               precondition_matrices),
                'combined_data' = combine_data(list_of_data = data_split, dim = dim)))
  }
}

#' @export
hierarchical_fusion_SMC_BLR_scalable <- function(N_schedule,
                                                 m_schedule,
                                                 time_schedule,
                                                 base_samples,
                                                 L,
                                                 dim,
                                                 data_split,
                                                 prior_means,
                                                 prior_variances,
                                                 C,
                                                 precondition = TRUE,
                                                 resampling_method = 'multi',
                                                 ESS_threshold = 0.5,
                                                 diffusion_estimator = 'Poisson',
                                                 beta_NB = 10,
                                                 seed = NULL,
                                                 n_cores = parallel::detectCores()) {
  if (!is.vector(N_schedule) | (length(N_schedule)!=(L-1))) {
    stop("hierarchical_fusion_SMC_BLR: N_schedule must be a vector of length (L-1)")
  } else if (!is.vector(m_schedule) | (length(m_schedule)!=(L-1))) {
    stop("hierarchical_fusion_SMC_BLR: m_schedule must be a vector of length (L-1)")
  } else if (!is.vector(time_schedule) | (length(time_schedule)!=(L-1))) {
    stop("hierarchical_fusion_SMC_BLR: time_schedule must be a vector of length (L-1)")
  } else if (!is.list(base_samples) | (length(base_samples)!=C)) {
    stop("hierarchical_fusion_SMC_BLR: base_samples must be a list of length C")
  } else if (!is.list(data_split) | length(data_split)!=C) {
    stop("hierarchical_fusion_SMC_BLR: data_split must be a list of length C")
  } else if (!all(sapply(data_split, function(sub_posterior) (is.list(sub_posterior) & identical(names(sub_posterior), c("y", "X")))))) {
    stop("hierarchical_fusion_SMC_BLR: each item in data_split must be a list of length 2 with names y and X")
  } else if (!all(sapply(1:C, function(i) is.vector(data_split[[i]]$y)))) {
    stop("hierarchical_fusion_SMC_BLR: for each i in 1:C, data_split[[i]]$y must be a vector")
  } else if (!all(sapply(1:C, function(i) is.matrix(data_split[[i]]$X)))) {
    stop("hierarchical_fusion_SMC_BLR: for each i in 1:C, data_split[[i]]$X must be a matrix")
  } else if (!all(sapply(1:C, function(i) ncol(data_split[[i]]$X)==dim))) {
    stop("hierarchical_fusion_SMC_BLR: for each i in 1:C, data_split[[i]]$X must be a matrix with dim columns")
  } else if (!all(sapply(1:m, function(i) length(data_split[[i]]$y)==nrow(data_split[[i]]$X)))) {
    stop("hierarchical_fusion_SMC_BLR: for each i in 1:m, length(data_split[[i]]$y) and nrow(data_split[[i]]$X) must be equal")
  } else if (!all(sapply(base_samples, is.matrix))) {
    stop("hierarchical_fusion_SMC_BLR: the sub-posterior samples in base_samples must be matrices")
  } else if (!all(sapply(base_samples, function(core) ncol(core)==dim))) {
    stop("hierarchical_fusion_SMC_BLR: the sub-posterior samples in base_samples must be matrices with dim columns")
  } else if (!is.vector(prior_means) | length(prior_means)!=dim) {
    stop("hierarchical_fusion_SMC_BLR: prior_means must be vectors of length dim")
  } else if (!is.vector(prior_variances) | length(prior_variances)!=dim) {
    stop("hierarchical_fusion_SMC_BLR: prior_variances must be vectors of length dim")
  } else if (ESS_threshold < 0 | ESS_threshold > 1) {
    stop("hierarchical_fusion_SMC_BLR: ESS_threshold must be between 0 and 1")
  }
  if (is.vector(m_schedule) & (length(m_schedule)==(L-1))) {
    for (l in (L-1):1) {
      if ((C/prod(m_schedule[(L-1):l]))%%1!=0) {
        stop("hierarchical_fusion_SMC_BLR: check that C/prod(m_schedule[(L-1):l])
              is an integer for l=L-1,...,1")
      }
    }
  } else {
    stop("hierarchical_fusion_SMC_BLR: m_schedule must be a vector of length (L-1)")
  }
  # we append 1 to the vector m_schedule to make the indices work later on when we call fusion
  m_schedule <- c(m_schedule, 1)
  # initialising results that we want to keep
  particles <- list()
  particles[[L]] <- initialise_particle_sets(samples_to_fuse = base_samples,
                                             multivariate = TRUE)
  proposed_samples <- list()
  data_inputs <- list()
  data_inputs[[L]] <- data_split
  time <- list()
  ESS <- list()
  CESS <- list()
  resampled <- list()
  precondition_matrices <- list()
  if (precondition) {
    precondition_matrices[[L]] <- lapply(base_samples, cov)
  } else {
    precondition_matrices[[L]] <- lapply(base_samples, function(c) diag(1, dim))
  }
  cat('Starting hierarchical fusion \n', file = 'hierarchical_fusion_SMC_BLR.txt')
  for (k in ((L-1):1)) {
    n_nodes <- max(C/prod(m_schedule[L:k]), 1)
    cat('########################\n', file = 'hierarchical_fusion_SMC_BLR.txt', append = T)
    cat('Starting to fuse', m_schedule[k], 'sub-posteriors for level', k, 'with time',
        time_schedule[k], ', which is using', n_cores, 'cores\n',
        file = 'hierarchical_fusion_SMC_BLR.txt', append = T)
    cat('At this level, the data is split up into', (C/prod(m_schedule[L:(k+1)])), 'subsets\n',
        file = 'hierarchical_fusion_SMC_BLR.txt', append = T)
    cat('There are', n_nodes, 'nodes at the next level each giving', N_schedule[k],
        'samples \n', file = 'hierarchical_fusion_SMC_BLR.txt', append = T)
    cat('########################\n', file = 'hierarchical_fusion_SMC_BLR.txt', append = T)
    fused <- lapply(X = 1:n_nodes, FUN = function(i) {
      previous_nodes <- ((m_schedule[k]*i)-(m_schedule[k]-1)):(m_schedule[k]*i)
      particles_to_fuse <- particles[[k+1]][previous_nodes]
      precondition_mats <- precondition_matrices[[k+1]][previous_nodes]
      parallel_fusion_SMC_BLR_scalable(particles_to_fuse = particles_to_fuse,
                                       N = N_schedule[k],
                                       m = m_schedule[k],
                                       time = time_schedule[k],
                                       dim = dim,
                                       data_split = data_inputs[[k+1]][previous_nodes],
                                       prior_means = prior_means,
                                       prior_variances = prior_variances,
                                       C = (C/prod(m_schedule[L:(k+1)])),
                                       precondition_matrices = precondition_mats,
                                       resampling_method = resampling_method,
                                       ESS_threshold = ESS_threshold,
                                       diffusion_estimator = diffusion_estimator,
                                       beta_NB = beta_NB,
                                       seed = seed,
                                       n_cores = n_cores,
                                       level = k,
                                       node = i)
    })
    # need to combine the correct samples
    particles[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$particles)
    proposed_samples[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$proposed_samples)
    time[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$time)
    ESS[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$ESS)
    CESS[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$CESS)
    resampled[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$resampled)
    precondition_matrices[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$precondition_matrices[[1]])
    data_inputs[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$combined_data)
  }
  cat('Completed hierarchical fusion\n', file = 'hierarchical_fusion_SMC_BLR.txt', append = T)
  if (length(particles[[1]])==1) {
    particles[[1]] <- particles[[1]][[1]]
    proposed_samples[[1]] <- proposed_samples[[1]][[1]]
    time[[1]] <- time[[1]][[1]]
    ESS[[1]] <- ESS[[1]][[1]]
    CESS[[1]] <- CESS[[1]][[1]]
    resampled[[1]] <- resampled[[1]][[1]]
    precondition_matrices[[1]] <- precondition_matrices[[1]][[1]]
    data_inputs[[1]] <- data_inputs[[1]][[1]]
  }
  return(list('particles' = particles,
              'proposed_samples' = proposed_samples,
              'time' = time,
              'ESS' = ESS,
              'CESS' = CESS,
              'resampled' = resampled,
              'precondition_matrices' = precondition_matrices,
              'data_inputs' = data_inputs,
              'diffusion_times' = time_schedule))
}
