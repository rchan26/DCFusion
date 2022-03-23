#' @export
ea_phi_BLR_DL <- function(beta,
                          y_labels,
                          X,
                          count,
                          prior_means,
                          prior_variances,
                          C,
                          precondition_mat) {
  if (is.vector(beta)) {
    return(ea_phi_BLR_DL_vec(beta = beta,
                             y_labels = y_labels,
                             X = X,
                             count = count,
                             prior_means = prior_means,
                             prior_variances = prior_variances,
                             C = C,
                             precondition_mat = precondition_mat))
  } else if (is.matrix(beta)) {
    return(ea_phi_BLR_DL_matrix(beta = beta,
                                y_labels = y_labels,
                                X = X,
                                count = count,
                                prior_means = prior_means,
                                prior_variances = prior_variances,
                                C = C,
                                precondition_mat = precondition_mat))
  }
  stop("ea_phi_BLR_DL: beta must be a vector or a matrix")
}

obtain_LR_MLE <- function(dim, data) {
  return(as.vector(unname(glm(formula = data$y ~ data$X[,2:dim], family = 'binomial')$coeff)))
}

#' Diffusion probability for the Exact Algorithm for Langevin diffusion for
#' Bayesian logistic regression
#' 
#' Simulate Langevin diffusion using the Exact Algorithm where target is the
#' posterior for a logistic regression model with Gaussian priors
#' 
#' @param dim dimension of the predictors (= p+1)
#' @param x0 start value (vector of length dim)
#' @param y end value (vector of length dim)
#' @param s start time
#' @param t end time
#' @param data list of length 4 where data[[c]]$y is the vector for y responses 
#'             and data[[c]]$X is the design matrix for the covariates for
#'             sub-posterior c, data[[c]]$full_data_count is the unique
#'             rows of the full data set with their counts and 
#'             data[[c]]$design_count is the unique rows of the design
#'             matrix and their counts
#' @param prior_means prior for means of predictors
#' @param prior_variances prior for variances of predictors
#' @param C overall number of sub-posteriors
#' @param precondition_mat precondition matrix
#' @param transform_mats list of transformation matrices where 
#'                       transform_mats$to_Z is the transformation matrix to Z space
#'                       and transform_mats$to_X is the transformation matrix to 
#'                       X space
#' @param cv_location string to determine what the location of the control variate
#'                    should be. Must be either 'mode' where the MLE estimator 
#'                    will be used or 'hypercube_centre' (default) to use the centre
#'                    of the simulated hypercube
#' @param diffusion_estimator choice of unbiased estimator for the Exact Algorithm
#'                            between "Poisson" (default) for Poisson estimator
#'                            and "NB" for Negative Binomial estimator
#' @param beta_NB beta parameter for Negative Binomial estimator (default 10)
#' @param gamma_NB_n_points number of points used in the trapezoidal estimation
#'                          of the integral found in the mean of the negative
#'                          binomial estimator (default is 2)
#' @param local_bounds logical value indicating if local bounds for the phi function
#'                     are used (default is TRUE)
#' @param logarithm logical value to determine if log probability is
#'                  returned (TRUE) or not (FALSE)
#' 
#' @return acceptance probability of simulating Langevin diffusion where target
#'         is the posterior for a logistic regression model with Gaussian priors
#' 
#' @export
ea_BLR_DL_PT <- function(dim,
                         x0,
                         y,
                         s,
                         t,
                         data,
                         transformed_design_mat,
                         prior_means,
                         prior_variances,
                         C,
                         precondition_mat,
                         transform_mats,
                         cv_location = 'hypercube_centre',
                         diffusion_estimator,
                         beta_NB = 10,
                         gamma_NB_n_points = 2,
                         local_bounds = TRUE,
                         logarithm) {
  # transform to preconditioned space
  z0 <- transform_mats$to_Z %*% x0
  zt <- transform_mats$to_Z %*% y
  # simulate layer information
  bes_layers <- layeredBB::multi_bessel_layer_simulation(dim = dim,
                                                         x = z0,
                                                         y = zt,
                                                         s = s,
                                                         t = t,
                                                         mult = 0.05)
  lbound_Z <- sapply(1:dim, function(d) bes_layers[[d]]$L)
  ubound_Z <- sapply(1:dim, function(d) bes_layers[[d]]$U)
  # calculate the lower and upper bounds of phi
  if (is.list(cv_location)) {
    if (names(cv_location)==c("beta_hat", "grad_log_hat")) {
      hypercube_vertices <- obtain_hypercube_vertices(bessel_layers = bes_layers,
                                                      dim = dim)
      bounds <- ea_phi_BLR_DL_bounds(beta_hat = as.vector(cv_location$beta_hat),
                                     grad_log_hat = as.vector(cv_location$grad_log_hat),
                                     dim = dim,
                                     transformed_X = transformed_design_mat,
                                     count = data$design_count$count,
                                     prior_variances = prior_variances,
                                     C = C,
                                     transform_mats = transform_mats,
                                     hypercube_vertices = hypercube_vertices,
                                     local_bounds = local_bounds)
    } else {
      stop("ea_BLR_BL_PT: cv_location must be a list or be set to \"hypercube_centre\"")
    }
  } else if (cv_location == 'hypercube_centre') {
    cv_location <- obtain_hypercube_centre_BLR(bessel_layers = bes_layers,
                                               transform_to_X = transform_mats$to_X,
                                               y_labels = data$full_data_count$y,
                                               X = as.matrix(subset(data$full_data_count, select = -c(y, count))),
                                               count = data$full_data_count$count,
                                               prior_means = prior_means,
                                               prior_variances = prior_variances,
                                               C = C)
    hypercube_vertices <- obtain_hypercube_vertices(bessel_layers = bes_layers,
                                                    dim = dim)
    bounds <- ea_phi_BLR_DL_bounds(beta_hat = as.vector(cv_location$beta_hat),
                                   grad_log_hat = as.vector(cv_location$grad_log_hat),
                                   dim = dim,
                                   transformed_X = transformed_design_mat,
                                   count = data$design_count$count,
                                   prior_variances = prior_variances,
                                   C = C,
                                   transform_mats = transform_mats,
                                   hypercube_vertices = hypercube_vertices,
                                   local_bounds = local_bounds)
  } else {
    stop("ea_BLR_BL_PT: cv_location must be a list or be set to \"hypercube_centre\"")
  }
  LX <- bounds$LB
  UX <- bounds$UB
  if (diffusion_estimator=='Poisson') {
    # simulate the number of points to simulate from Poisson distribution
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
      phi <- ea_phi_BLR_DL_matrix(beta = sim_path,
                                  y_labels = data$full_data_count$y,
                                  X = as.matrix(subset(data$full_data_count, select = -c(y, count))),
                                  count = data$full_data_count$count,
                                  prior_means = prior_means,
                                  prior_variances = prior_variances,
                                  C = C,
                                  precondition_mat = precondition_mat)
      terms <- (UX-phi$phi)
      log_acc_prob <- sum(log(terms))
      if (any(terms < 0)) {
        stop('Some of (UX-phi) are < 0. Try local_bounds = FALSE to use global bounds.')
      } else if (any((phi$phi - LX) < 0)) {
        stop('Some of (phi-LX) are < 0. Try local_bounds = FALSE to use global bounds.')
      }
    }
    if (logarithm) {
      return(list('phi' = -LX*(t-s) - kap*log(UX-LX) + log_acc_prob,
                  'LX' = LX, 'UX' = UX, 'kap' = kap, 'bound_intensity' = (UX-LX)*(t-s)))
    } else {
      return(list('phi' = exp(-LX*(t-s) - kap*log(UX-LX) + log_acc_prob),
                  'LX' = LX, 'UX' = UX, 'kap' = kap, 'bound_intensity' = (UX-LX)*(t-s)))
    }
  } else if (diffusion_estimator=="NB") {
    # integral estimate for gamma in NB estimator
    h <- (t-s)/(gamma_NB_n_points-1)
    times_to_eval <- seq(from = s, to = t, by = h)
    integral_estimate <- gamma_NB_BLR(times = times_to_eval,
                                      h = h,
                                      x0 = x0,
                                      y = y,
                                      s = s,
                                      t = t,
                                      y_labels = data$full_data_count$y,
                                      X = as.matrix(subset(data$full_data_count, select = -c(y, count))),
                                      count = data$full_data_count$count,
                                      prior_means = prior_means,
                                      prior_variances = prior_variances,
                                      C = C,
                                      precondition_mat = precondition_mat)
    gamma_NB <- (t-s)*UX - integral_estimate
    kap <- rnbinom(1, size = beta_NB, mu = gamma_NB)
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
      phi <- ea_phi_BLR_DL_matrix(beta = sim_path,
                                  y_labels = data$full_data_count$y,
                                  X = as.matrix(subset(data$full_data_count, select = -c(y, count))),
                                  count = data$full_data_count$count,
                                  prior_means = prior_means,
                                  prior_variances = prior_variances,
                                  C = C,
                                  precondition_mat = precondition_mat)
      terms <- (UX-phi$phi)
      log_acc_prob <- sum(log(terms))
      if (any(terms < 0)) {
        stop('Some of (UX-phi) are < 0. Try local_bounds = FALSE to use global bounds.')
      } else if (any((phi$phi - LX) < 0)) {
        stop('Some of (phi-LX) are < 0. Try local_bounds = FALSE to use global bounds.')
      }
    }
    log_middle_term <- kap*log(t-s) + lgamma(beta_NB) + (beta_NB+kap)*log(beta_NB+gamma_NB) -
      lgamma(beta_NB+kap) - beta_NB*log(beta_NB) - kap*log(gamma_NB)
    if (logarithm) {
      return(list('phi' = -UX*(t-s) + log_middle_term + log_acc_prob,
                  'LX' = LX, 'UX' = UX, 'kap' = kap, 'bound_intensity' = (UX-LX)*(t-s)))
    } else {
      return(list('phi' = exp(-UX*(t-s) + log_middle_term + log_acc_prob),
                  'LX' = LX, 'UX' = UX, 'kap' = kap, 'bound_intensity' = (UX-LX)*(t-s)))
    }
  } else {
    stop("ea_BLR_DL_PT: diffusion_estimator must be set to either \'Poisson\' or \'NB\'")
  }
}

#' Q Importance Sampling Step
#'
#' Q Importance Sampling weighting for Bayesian logistic regression
#'
#' @param particle_set particles set prior to Q importance sampling step
#' @param m number of sub-posteriors to combine
#' @param time time T for fusion algorithm
#' @param dim dimension of the predictors (= p+1)
#' @param data_split list of length m where each item is a list of length 4 where
#'                   for c=1,...,m, data_split[[c]]$y is the vector for y responses and
#'                   data_split[[c]]$X is the design matrix for the covariates for
#'                   sub-posterior c, data_split[[c]]$full_data_count is the unique
#'                   rows of the full data set with their counts and 
#'                   data_split[[c]]$design_count is the unique rows of the design
#'                   matrix and their counts
#' @param prior_means prior for means of predictors
#' @param prior_variances prior for variances of predictors
#' @param C overall number of sub-posteriors
#' @param proposal_cov proposal covariance of Gaussian distribution for Fusion
#' @param precondition_matrices list of length m, where precondition_matrices[[c]]
#'                               is the precondition matrix for sub-posterior c
#' @param inv_precondition_matrices list of length m, where inv_precondition_matrices[[c]]
#'                                  is the inverse precondition matrix for sub-posterior c
#' @param cv_location string to determine what the location of the control variate
#'                    should be. Must be either 'mode' where the MLE estimator 
#'                    will be used or 'hypercube_centre' (default) to use the centre
#'                    of the simulated hypercube
#' @param diffusion_estimator choice of unbiased estimator for the Exact Algorithm
#'                            between "Poisson" (default) for Poisson estimator
#'                            and "NB" for Negative Binomial estimator
#' @param beta_NB beta parameter for Negative Binomial estimator (default 10)
#' @param gamma_NB_n_points number of points used in the trapezoidal estimation
#'                          of the integral found in the mean of the negative
#'                          binomial estimator (default is 2)
#' @param local_bounds logical value indicating if local bounds for the phi function
#'                     are used (default is TRUE)
#' @param seed seed number - default is NULL, meaning there is no seed
#' @param n_cores number of cores to use
#' @param cl an object of class "cluster" for parallel computation in R. If none
#'           is passed, then one is created and used within this function
#' @param level indicates which level this is for the hierarchy (default 1)
#' @param node indicates which node this is for the hierarchy (default 1)
#' @param print_progress_iters number of iterations between each progress update
#'                             (default is 1000). If NULL, progress will only
#'                             be updated when importance sampling is finished
#'
#' @return An updated particle set
#' 
#' @export
Q_IS_BLR <- function(particle_set,
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
  transformed_design_matrices <- lapply(1:m, function(c) {
    as.matrix(subset(data_split[[c]]$design_count, select = -count)) %*% transform_matrices[[c]]$to_X
  })
  N <- particle_set$N
  # ---------- creating parallel cluster
  if (is.null(cl)) {
    cl <- parallel::makeCluster(n_cores, setup_strategy = "sequential", outfile = "SMC_BLR_outfile.txt")
    parallel::clusterExport(cl, varlist = ls("package:layeredBB"))
    parallel::clusterExport(cl, varlist = ls("package:DCFusion"))
    close_cluster <- TRUE
  } else {
    close_cluster <- FALSE
  }
  parallel::clusterExport(cl, envir = environment(), varlist = ls())
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
    bound_intensity <- matrix(nrow = split_N, ncol = m)
    kap <- matrix(nrow = split_N, ncol = m)
    cat('Level:', level, '|| Node:', node, '|| Core:', core, '|| START \n',
        file = 'Q_IS_BLR_progress.txt', append = T)
    if (is.null(print_progress_iters)) {
      print_progress_iters <- split_N
    }
    for (i in 1:split_N) {
      phi <- lapply(1:m, function(c) {
        tryCatch(expr = ea_BLR_DL_PT(dim = dim,
                                     x0 = as.vector(split_x_samples[[core]][[i]][c,]),
                                     y = as.vector(y_samples[i,]),
                                     s = 0,
                                     t = time,
                                     data = data_split[[c]][counts],
                                     transformed_design_mat = transformed_design_matrices[[c]],
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
                                transformed_design_mat = transformed_design_matrices[[c]],
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
                                logarithm = TRUE)})})
      log_Q_weights[i] <- sum(sapply(1:m, function(c) phi[[c]]$phi))
      bound_intensity[i,] <- sapply(1:m, function(c) phi[[c]]$bound_intensity)
      kap[i,] <- sapply(1:m, function(c) phi[[c]]$kap)
    }
    if (i%%print_progress_iters==0) {
      cat('Level:', level, '|| Node:', node, '|| Core:', core, '||', i, '/',
          split_N, '\n', file = 'Q_IS_BLR_progress.txt', append = T)
    }
    cat('Completed: Level:', level, '|| Node:', node, '|| Core:', core, '||', split_N, '/',
        split_N, '\n', file = 'Q_IS_BLR_progress.txt', append = T)
    return(list('y_samples' = y_samples,
                'log_Q_weights' = log_Q_weights,
                'bound_intensity' = bound_intensity,
                'kap' = kap))
  })
  if (close_cluster) {
    parallel::stopCluster(cl)
  }
  # unlist the proposed samples for y and their associated log Q weights
  log_Q_weights <- unlist(lapply(1:length(split_x_samples), function(i) {
    Q_weighted_samples[[i]]$log_Q_weights}))
  # ---------- update particle set
  # update the weights and return updated particle set
  particle_set$y_samples <- do.call(rbind, lapply(1:length(split_x_samples), function(i) {
    Q_weighted_samples[[i]]$y_samples}))
  # normalise weight
  norm_weights <- particle_ESS(log_weights = particle_set$log_weights + log_Q_weights)
  particle_set$log_weights <- norm_weights$log_weights
  particle_set$normalised_weights <- norm_weights$normalised_weights
  particle_set$ESS <- norm_weights$ESS
  # calculate the conditional ESS (i.e. the 1/sum(inc_change^2))
  # where inc_change is the incremental change in weight (= log_Q_weights)
  particle_set$CESS[2] <- particle_ESS(log_weights = log_Q_weights)$ESS
  # set the resampled indicator to FALSE
  particle_set$resampled[2] <- FALSE
  return(list('particle_set' = particle_set,
              'phi_bound_intensity' = do.call(rbind, lapply(1:length(split_x_samples), function(i) {
                Q_weighted_samples[[i]]$bound_intensity})),
              'phi_kappa' = do.call(rbind, lapply(1:length(split_x_samples), function(i) {
                Q_weighted_samples[[i]]$kap}))))
}

#' Generalised Monte Carlo Fusion [parallel]
#' 
#' Generalised Monte Carlo Fusion for Bayesian Logistic Regression
#'
#' @param particles_to_fuse list of length m, where particles_to_fuse[c] contains
#'                          the particles for the c-th sub-posterior. Can
#'                          initialise a this from list of sub-posterior samples
#'                          by using the initialise_particle_sets function
#' @param N number of samples
#' @param m number of sub-posteriors to combine
#' @param time time T for fusion algorithm
#' @param dim dimension of the predictors (= p+1)
#' @param data_split list of length m where each item is a list of length 4 where
#'                   for c=1,...,m, data_split[[c]]$y is the vector for y responses and
#'                   data_split[[c]]$X is the design matrix for the covariates for
#'                   sub-posterior c, data_split[[c]]$full_data_count is the unique
#'                   rows of the full data set with their counts and 
#'                   data_split[[c]]$design_count is the unique rows of the design
#'                   matrix and their counts
#' @param prior_means prior for means of predictors
#' @param prior_variances prior for variances of predictors
#' @param C overall number of sub-posteriors
#' @param precondition_matrices list of length m, where precondition_matrices[[c]]
#'                               is the precondition matrix for sub-posterior c
#' @param resampling_method method to be used in resampling, default is multinomial
#'                          resampling ('multi'). Other choices are stratified
#'                          resampling ('strat'), systematic resampling ('system'),
#'                          residual resampling ('resid')
#' @param ESS_threshold number between 0 and 1 defining the proportion of the
#'                      number of samples that ESS needs to be lower than for
#'                      resampling (i.e. resampling is carried out only when
#'                      ESS < N*ESS_threshold)
#' @param cv_location string to determine what the location of the control variate
#'                    should be. Must be either 'mode' where the MLE estimator 
#'                    will be used or 'hypercube_centre' (default) to use the centre
#'                    of the simulated hypercube
#' @param diffusion_estimator choice of unbiased estimator for the Exact Algorithm
#'                            between "Poisson" (default) for Poisson estimator
#'                            and "NB" for Negative Binomial estimator
#' @param beta_NB beta parameter for Negative Binomial estimator (default 10)
#' @param gamma_NB_n_points number of points used in the trapezoidal estimation
#'                          of the integral found in the mean of the negative
#'                          binomial estimator (default is 2)
#' @param local_bounds logical value indicating if local bounds for the phi function
#'                     are used (default is TRUE)
#' @param seed seed number - default is NULL, meaning there is no seed
#' @param n_cores number of cores to use
#' @param cl an object of class "cluster" for parallel computation in R. If none
#'           is passed, then one is created and used within this function
#' @param level indicates which level this is for the hierarchy (default 1)
#' @param node indicates which node this is for the hierarchy (default 1)
#' @param print_progress_iters number of iterations between each progress update
#'                             (default is 1000). If NULL, progress will only
#'                             be updated when importance sampling is finished
#'
#' @return A list with components:
#' \describe{
#'   \item{particles}{particles returned from fusion sampler}
#'   \item{proposed_samples}{proposal samples from fusion sampler}
#'   \item{time}{run-time of fusion sampler}
#'   \item{ESS}{effective sample size of the particles after each step}
#'   \item{CESS}{conditional effective sample size of the particles after each step}
#'   \item{resampled}{boolean value to indicate if particles were resampled
#'                    after each time step}
#'   \item{precondition_matrices}{list of length 2 where precondition_matrices[[2]]
#'                                are the pre-conditioning matrices that were used
#'                                and precondition_matrices[[1]] are the combined
#'                                precondition matrices}
#'   \item{combined_data}{combined data for the fusion density}
#' }
#'
#' @export
parallel_fusion_SMC_BLR <- function(particles_to_fuse,
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
    stop("parallel_fusion_SMC_BLR: particles_to_fuse must be a list of length m")
  } else if (!all(sapply(particles_to_fuse, function(sub_posterior) ("particle" %in% class(sub_posterior))))) {
    stop("parallel_fusion_SMC_BLR: particles in particles_to_fuse must be \"particle\" objects")
  } else if (!all(sapply(particles_to_fuse, function(sub_posterior) is.matrix(sub_posterior$y_samples)))) {
    stop("parallel_fusion_SMC_BLR: the particles' samples for y should all be matrices")
  } else if (!all(sapply(particles_to_fuse, function(sub_posterior) ncol(sub_posterior$y_samples)==dim))) {
    stop("parallel_fusion_SMC_BLR: the particles' samples for y should all be matrices with dim columns")
  } else if (!is.list(data_split) | length(data_split)!=m) {
    stop("parallel_fusion_SMC_BLR: data_split must be a list of length m")
  } else if (!all(sapply(data_split, function(sub_posterior) (is.list(sub_posterior) & identical(names(sub_posterior), c("y", "X", "full_data_count", "design_count")))))) {
    stop("parallel_fusion_SMC_BLR: each item in data_split must be a list of length 4 with names \'y\', \'X\', \'full_data_count\', \'design_count\'")
  } else if (!all(sapply(1:m, function(i) is.vector(data_split[[i]]$y)))) {
    stop("parallel_fusion_SMC_BLR: for each i in 1:m, data_split[[i]]$y must be a vector")
  } else if (!all(sapply(1:m, function(i) is.matrix(data_split[[i]]$X)))) {
    stop("parallel_fusion_SMC_BLR: for each i in 1:m, data_split[[i]]$X must be a matrix")
  } else if (!all(sapply(1:m, function(i) ncol(data_split[[i]]$X)==dim))) {
    stop("parallel_fusion_SMC_BLR: for each i in 1:m, data_split[[i]]$X must be a matrix with dim columns")
  } else if (!all(sapply(1:m, function(i) length(data_split[[i]]$y)==nrow(data_split[[i]]$X)))) {
    stop("parallel_fusion_SMC_BLR: for each i in 1:m, length(data_split[[i]]$y) and nrow(data_split[[i]]$X) must be equal")
  } else if (!all(sapply(1:m, function(i) is.data.frame(data_split[[i]]$full_data_count)))) {
    stop("parallel_fusion_SMC_BLR: for each i in 1:m, data_split[[i]]$full_data_count must be a data frame")
  } else if (!all(sapply(1:m, function(i) is.data.frame(data_split[[i]]$design_count)))) {
    stop("parallel_fusion_SMC_BLR: for each i in 1:m, data_split[[i]]$design_count must be a data frame")
  } else if (!is.vector(prior_means) | length(prior_means)!=dim) {
    stop("parallel_fusion_SMC_BLR: prior_means must be vectors of length dim")
  } else if (!is.vector(prior_variances) | length(prior_variances)!=dim) {
    stop("parallel_fusion_SMC_BLR: prior_variances must be vectors of length dim")
  } else if (!is.list(precondition_matrices) | (length(precondition_matrices)!=m)) {
    stop("parallel_fusion_SMC_BLR: precondition_matrices must be a list of length m")
  } else if (!(diffusion_estimator %in% c('Poisson', 'NB'))) {
    stop("parallel_fusion_SMC_BLR: diffusion_estimator must be set to either \'Poisson\' or \'NB\'")
  } else if ((ESS_threshold < 0) | (ESS_threshold > 1)) {
    stop("parallel_fusion_SMC_BLR: ESS_threshold must be between 0 and 1")
  } else if ((cv_location != 'mode') & (cv_location != 'hypercube_centre')) {
    stop("parallel_fusion_SMC_BLR: cv_location must be either \"mode\" or \"hypercube_centre\"")
  } else if (!any(class(cl)=="cluster") & !is.null(cl)) {
    stop("parallel_fusion_SMC_BLR: cl must be a \"cluster\" object or NULL")
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
  # importance sampling for rho step
  particles <- rho_IS_multivariate(particles_to_fuse = particles_to_fuse,
                                   dim = dim,
                                   N = N,
                                   m = m,
                                   time = time,
                                   inv_precondition_matrices = inv_precondition_matrices,
                                   inverse_sum_inv_precondition_matrices = inverse_sum_matrices(inv_precondition_matrices),
                                   number_of_steps = 2,
                                   resampling_method = resampling_method,
                                   seed = seed,
                                   n_cores = n_cores,
                                   cl = cl)
  # record ESS and CESS after rho step
  ESS <- c('rho' = particles$ESS)
  CESS <- c('rho' = particles$CESS[1])
  # ----------- resample particles (only resample if ESS < N*ESS_threshold)
  if (particles$ESS < N*ESS_threshold) {
    resampled <- c('rho' = TRUE)
    particles <- resample_particle_x_samples(N = N,
                                             particle_set = particles,
                                             multivariate = TRUE,
                                             step = 1,
                                             resampling_method = resampling_method,
                                             seed = seed)
  } else {
    resampled <- c('rho' = FALSE)
  }
  # ---------- second importance sampling step
  # unbiased estimator for Q
  Q <- Q_IS_BLR(particle_set = particles,
                m = m,
                time = time,
                dim = dim,
                data_split = data_split,
                prior_means = prior_means,
                prior_variances = prior_variances,
                C = C,
                proposal_cov = calculate_proposal_cov(time = time, weights = inv_precondition_matrices),
                precondition_matrices = precondition_matrices,
                inv_precondition_matrices = inv_precondition_matrices,
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
  particles <- Q$particle_set
  # record ESS and CESS after Q step
  ESS['Q'] <- particles$ESS
  CESS['Q'] <- particles$CESS[2]
  names(CESS) <- c('rho', 'Q')
  # record proposed samples
  proposed_samples <- particles$y_samples
  # ----------- resample particles (only resample if ESS < N*ESS_threshold)
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
    new_precondition_matrices <- list(diag(1, dim), precondition_matrices)
  } else {
    new_precondition_matrices <- list(inverse_sum_matrices(inv_precondition_matrices),
                                      precondition_matrices)
  }
  return(list('particles' = particles,
              'proposed_samples' = proposed_samples,
              'time' = (proc.time()-pcm)['elapsed'],
              'ESS' = ESS,
              'CESS' = CESS,
              'resampled' = resampled,
              'phi_bound_intensity' = Q$phi_bound_intensity,
              'phi_kappa' = Q$phi_kappa,
              'precondition_matrices' = new_precondition_matrices,
              'combined_data' = combine_data(list_of_data = data_split, dim = dim)))
}

#' (Balanced Binary) D&C Monte Carlo Fusion using SMC
#' 
#' (Balanced Binary) D&C Monte Carlo Fusion using SMC for Bayesian Logistic Regression
#'
#' @param N_schedule vector of length (L-1), where N_schedule[l] is the
#'                   number of samples per node at level l
#' @param m_schedule vector of length (L-1), where m_schedule[k] is the number
#'                   of samples to fuse for level k
#' @param time_schedule vector of length (L-1), where time_schedule[k] is time
#'                      T for algorithm for level k
#' @param base_samples list of length C, where base_samples[[c]] contains
#'                     the samples for the c-th node in the level
#' @param L total number of levels in the hierarchy
#' @param dim dimension of the predictors (= p+1)
#' @param data_split list of length m where each item is a list of length 4 where
#'                   for c=1,...,m, data_split[[c]]$y is the vector for y responses and
#'                   data_split[[c]]$X is the design matrix for the covariates for
#'                   sub-posterior c, data_split[[c]]$full_data_count is the unique
#'                   rows of the full data set with their counts and 
#'                   data_split[[c]]$design_count is the unique rows of the design
#'                   matrix and their counts
#' @param prior_means prior for means of predictors
#' @param prior_variances prior for variances of predictors
#' @param C number of sub-posteriors at the base level
#' @param precondition either a logical value to determine if preconditioning matrices are
#'                     used (TRUE - and is set to be the variance of the sub-posterior samples)
#'                     or not (FALSE - and is set to be the identity matrix for all sub-posteriors),
#'                     or a list of length (1/start_beta) where precondition[[c]]
#'                     is the preconditioning matrix for sub-posterior c. Default is TRUE
#' @param resampling_method method to be used in resampling, default is
#'                          multinomial resampling ('multi'). Other choices are
#'                          stratified ('strat'), systematic ('system'),
#'                          residual ('resid')
#' @param ESS_threshold number between 0 and 1 defining the proportion
#'                      of the number of samples that ESS needs to be
#'                      lower than for resampling (i.e. resampling is carried
#'                      out only when ESS < N*ESS_threshold)
#' @param cv_location string to determine what the location of the control variate
#'                    should be. Must be either 'mode' where the MLE estimator 
#'                    will be used or 'hypercube_centre' (default) to use the centre
#'                    of the simulated hypercube
#' @param diffusion_estimator choice of unbiased estimator for the Exact Algorithm
#'                            between "Poisson" (default) for Poisson estimator
#'                            and "NB" for Negative Binomial estimator
#' @param beta_NB beta parameter for Negative Binomial estimator (default 10)
#' @param gamma_NB_n_points number of points used in the trapezoidal estimation
#'                          of the integral found in the mean of the negative
#'                          binomial estimator (default is 2)
#' @param local_bounds logical value indicating if local bounds for the phi function
#'                     are used (default is TRUE)
#' @param seed seed number - default is NULL, meaning there is no seed
#' @param n_cores number of cores to use
#' @param print_progress_iters number of iterations between each progress update
#'                             (default is 1000). If NULL, progress will only
#'                             be updated when importance sampling is finished
#'
#' @return A list with components:
#' \describe{
#'   \item{particles}{list of length (L-1), where particles[[l]][[i]] are the
#'                    particles for level l, node i}
#'   \item{proposed_samples}{list of length (L-1), where proposed_samples[[l]][[i]]
#'                           are the proposed samples for level l, node i}
#'   \item{time}{list of length (L-1), where time[[l]][[i]] is the run time
#'               for level l, node i}
#'   \item{ESS}{list of length (L-1), where ESS[[l]][[i]] is the effective
#'              sample size of the particles after each step BEFORE deciding
#'              whether or not to resample for level l, node i}
#'   \item{CESS}{list of length (L-1), where ESS[[l]][[i]] is the conditional
#'               effective sample size of the particles after each step}
#'   \item{resampled}{list of length (L-1), where resampled[[l]][[i]] is a
#'                    boolean value to record if the particles were resampled
#'                    after each step; rho and Q for level l, node i}
#'   \item{precondition_matrices}{pre-conditioning matrices that were used}
#'   \item{data_inputs}{list of length (L-1), where data_inputs[[l]][[i]] is the
#'                      data input for the sub-posterior in level l, node i}
#'   \item{diffusion_times}{vector of length (L-1), where diffusion_times[l]
#'                          are the times for fusion in level l}
#' }
#'
#' @export
bal_binary_fusion_SMC_BLR <- function(N_schedule,
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
                                      cv_location = 'hypercube_centre',
                                      diffusion_estimator = 'Poisson',
                                      beta_NB = 10,
                                      gamma_NB_n_points = 2,
                                      local_bounds = TRUE,
                                      seed = NULL,
                                      n_cores = parallel::detectCores(),
                                      print_progress_iters = 100) {
  if (!is.vector(N_schedule) | (length(N_schedule)!=(L-1))) {
    stop("bal_binary_fusion_SMC_BLR: N_schedule must be a vector of length (L-1)")
  } else if (!is.vector(m_schedule) | (length(m_schedule)!=(L-1))) {
    stop("bal_binary_fusion_SMC_BLR: m_schedule must be a vector of length (L-1)")
  } else if (!is.vector(time_schedule) | (length(time_schedule)!=(L-1))) {
    stop("bal_binary_fusion_SMC_BLR: time_schedule must be a vector of length (L-1)")
  } else if (!is.list(base_samples) | (length(base_samples)!=C)) {
    stop("bal_binary_fusion_SMC_BLR: base_samples must be a list of length C")
  } else if (!is.list(data_split) | length(data_split)!=C) {
    stop("bal_binary_fusion_SMC_BLR: data_split must be a list of length C")
  } else if (!all(sapply(data_split, function(sub_posterior) (is.list(sub_posterior) & identical(names(sub_posterior), c("y", "X", "full_data_count", "design_count")))))) {
    stop("bal_binary_fusion_SMC_BLR: each item in data_split must be a list of length 4 with names \'y\', \'X\', \'full_data_count\', \'design_count\'")
  } else if (!all(sapply(1:C, function(i) is.vector(data_split[[i]]$y)))) {
    stop("bal_binary_fusion_SMC_BLR: for each i in 1:C, data_split[[i]]$y must be a vector")
  } else if (!all(sapply(1:C, function(i) is.matrix(data_split[[i]]$X)))) {
    stop("bal_binary_fusion_SMC_BLR: for each i in 1:C, data_split[[i]]$X must be a matrix")
  } else if (!all(sapply(1:C, function(i) ncol(data_split[[i]]$X)==dim))) {
    stop("bal_binary_fusion_SMC_BLR: for each i in 1:C, data_split[[i]]$X must be a matrix with dim columns")
  } else if (!all(sapply(1:C, function(i) length(data_split[[i]]$y)==nrow(data_split[[i]]$X)))) {
    stop("bal_binary_fusion_SMC_BLR: for each i in 1:C, length(data_split[[i]]$y) and nrow(data_split[[i]]$X) must be equal")
  } else if (!all(sapply(1:C, function(i) is.data.frame(data_split[[i]]$full_data_count)))) {
    stop("bal_binary_fusion_SMC_BLR: for each i in 1:C, data_split[[i]]$full_data_count must be a data frame")
  } else if (!all(sapply(1:C, function(i) is.data.frame(data_split[[i]]$design_count)))) {
    stop("bal_binary_fusion_SMC_BLR: for each i in 1:C, data_split[[i]]$design_count must be a data frame")
  } else if (!is.vector(prior_means) | length(prior_means)!=dim) {
    stop("bal_binary_fusion_SMC_BLR: prior_means must be vectors of length dim")
  } else if (!is.vector(prior_variances) | length(prior_variances)!=dim) {
    stop("bal_binary_fusion_SMC_BLR: prior_variances must be vectors of length dim")
  } else if (ESS_threshold < 0 | ESS_threshold > 1) {
    stop("bal_binary_fusion_SMC_BLR: ESS_threshold must be between 0 and 1")
  }
  if (is.vector(m_schedule) & (length(m_schedule)==(L-1))) {
    for (l in (L-1):1) {
      if ((C/prod(m_schedule[(L-1):l]))%%1!=0) {
        stop("bal_binary_fusion_SMC_BLR: check that C/prod(m_schedule[(L-1):l])
              is an integer for l=L-1,...,1")
      }
    }
  } else {
    stop("bal_binary_fusion_SMC_BLR: m_schedule must be a vector of length (L-1)")
  }
  m_schedule <- c(m_schedule, 1)
  particles <- list()
  if (all(sapply(base_samples, function(sub) class(sub)=='particle'))) {
    particles[[L]] <- base_samples
  } else if (all(sapply(base_samples, is.matrix))) {
    if (!all(sapply(base_samples, function(core) ncol(core)==dim))) {
      stop("bal_binary_fusion_SMC_BLR: the sub-posterior samples in base_samples must be matrices with dim columns")
    }
    particles[[L]] <- initialise_particle_sets(samples_to_fuse = base_samples,
                                               multivariate = TRUE,
                                               number_of_steps = 2)
  } else {
    stop("bal_binary_fusion_SMC_BLR: base_samples must be a list of length C
         containing either items of class \"particle\" (representing particle 
         approximations of the sub-posteriors) or are matrices with dim columns
         (representing un-normalised sample approximations of the sub-posteriors)")
  }
  proposed_samples <- list()
  data_inputs <- list()
  data_inputs[[L]] <- data_split
  time <- list()
  ESS <- list()
  CESS <- list()
  resampled <- list()
  phi_bound_intensity <- list()
  phi_kappa <- list()
  precondition_matrices <- list()
  if (is.logical(precondition)) {
    if (precondition) {
      precondition_matrices[[L]] <- lapply(base_samples, cov)
    } else {
      precondition_matrices[[L]] <- lapply(base_samples, function(c) diag(1, dim))
    }
  } else if (is.list(precondition)) {
    if (length(precondition)==C & all(sapply(precondition, is.matrix))) {
      if (all(sapply(precondition, function(sub) ncol(sub)==dim))) {
        precondition_matrices[[L]] <- precondition  
      }
    }
  } else {
    stop("bal_binary_fusion_SMC_BLR: precondition must be a logical indicating 
          whether or not a preconditioning matrix should be used, or a list of
          length C, where precondition[[c]] is the preconditioning matrix for
          the c-th sub-posterior")
  }
  cl <- parallel::makeCluster(n_cores, setup_strategy = "sequential", outfile = "SMC_BLR_outfile.txt")
  parallel::clusterExport(cl, envir = environment(), varlist = ls())
  parallel::clusterExport(cl, varlist = ls("package:DCFusion"))
  parallel::clusterExport(cl, varlist = ls("package:layeredBB"))
  cat('Starting bal_binary fusion \n', file = 'bal_binary_fusion_SMC_BLR.txt')
  for (k in ((L-1):1)) {
    n_nodes <- max(C/prod(m_schedule[L:k]), 1)
    cat('########################\n', file = 'bal_binary_fusion_SMC_BLR.txt', append = T)
    cat('Starting to fuse', m_schedule[k], 'sub-posteriors for level', k, 'with time',
        time_schedule[k], ', which is using', n_cores, 'cores\n',
        file = 'bal_binary_fusion_SMC_BLR.txt', append = T)
    cat('At this level, the data is split up into', (C/prod(m_schedule[L:(k+1)])), 'subsets\n',
        file = 'bal_binary_fusion_SMC_BLR.txt', append = T)
    cat('There are', n_nodes, 'nodes at the next level each giving', N_schedule[k],
        'samples \n', file = 'bal_binary_fusion_SMC_BLR.txt', append = T)
    cat('########################\n', file = 'bal_binary_fusion_SMC_BLR.txt', append = T)
    fused <- lapply(X = 1:n_nodes, FUN = function(i) {
      previous_nodes <- ((m_schedule[k]*i)-(m_schedule[k]-1)):(m_schedule[k]*i)
      particles_to_fuse <- particles[[k+1]][previous_nodes]
      precondition_mats <- precondition_matrices[[k+1]][previous_nodes]
      parallel_fusion_SMC_BLR(particles_to_fuse = particles_to_fuse,
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
                              cv_location = cv_location,
                              diffusion_estimator = diffusion_estimator,
                              beta_NB = beta_NB,
                              gamma_NB_n_points = gamma_NB_n_points,
                              local_bounds = local_bounds,
                              seed = seed,
                              n_cores = n_cores,
                              cl = cl,
                              level = k,
                              node = i,
                              print_progress_iters = print_progress_iters)
    })
    particles[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$particles)
    proposed_samples[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$proposed_samples)
    time[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$time)
    ESS[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$ESS)
    CESS[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$CESS)
    resampled[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$resampled)
    phi_bound_intensity[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$phi_bound_intensity)
    phi_kappa[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$phi_kappa)
    precondition_matrices[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$precondition_matrices[[1]])
    data_inputs[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$combined_data)
  }
  parallel::stopCluster(cl)
  cat('Completed bal_binary fusion\n', file = 'bal_binary_fusion_SMC_BLR.txt', append = T)
  if (length(particles[[1]])==1) {
    particles[[1]] <- particles[[1]][[1]]
    proposed_samples[[1]] <- proposed_samples[[1]][[1]]
    time[[1]] <- time[[1]][[1]]
    ESS[[1]] <- ESS[[1]][[1]]
    CESS[[1]] <- CESS[[1]][[1]]
    resampled[[1]] <- resampled[[1]][[1]]
    phi_bound_intensity[[1]] <- phi_bound_intensity[[1]][[1]]
    phi_kappa[[1]] <- phi_kappa[[1]][[1]]
    precondition_matrices[[1]] <- precondition_matrices[[1]][[1]]
    data_inputs[[1]] <- data_inputs[[1]][[1]]
  }
  return(list('particles' = particles,
              'proposed_samples' = proposed_samples,
              'time' = time,
              'ESS' = ESS,
              'CESS' = CESS,
              'resampled' = resampled,
              'phi_bound_intensity' = phi_bound_intensity,
              'phi_kappa' = phi_kappa,
              'precondition_matrices' = precondition_matrices,
              'data_inputs' = data_inputs,
              'diffusion_times' = time_schedule))
}
