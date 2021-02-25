#' @export
ea_phi_BLR_DL <- function(beta,
                          y_labels,
                          X,
                          prior_means,
                          prior_variances,
                          C,
                          precondition_mat,
                          transform_mat) {
  if (is.vector(beta)) {
    return(ea_phi_BLR_DL_vec(beta = beta,
                             y_labels = y_labels,
                             X = X,
                             prior_means = prior_means,
                             prior_variances = prior_variances,
                             C = C,
                             precondition_mat = precondition_mat,
                             transform_mat = transform_mat))
  } else if (is.matrix(beta)) {
    return(ea_phi_BLR_DL_matrix(beta = beta,
                                y_labels = y_labels,
                                X = X,
                                prior_means = prior_means,
                                prior_variances = prior_variances,
                                C = C,
                                precondition_mat = precondition_mat,
                                transform_mat = transform_mat))
  }
  stop("ea_phi_BLR_DL: x must be a vector or a matrix")
}

#' @export
ea_phi_BLR_DL_bounds <- function(initial_parameters,
                                 y_labels,
                                 X,
                                 prior_means,
                                 prior_variances,
                                 C,
                                 precondition_mat,
                                 transform_mat,
                                 lower,
                                 upper,
                                 bounds_multiplier) {
  if (bounds_multiplier < 1) {
    stop("ea_phi_BLR_DL_bounds: bounds_multipler should be greater than or equal to 1")
  }
  # checking that initial parameters is between lower and upper bounds
  if (!all(lower <= initial_parameters)) {
    for (d in which(lower > initial_parameters)) {
      initial_parameters[d] <- (lower[d]+upper[d])/2
    }
  }
  if (!all(upper >= initial_parameters)) {
    for (d in which(upper < initial_parameters)) {
      initial_parameters[d] <- (lower[d]+upper[d])/2
    }
  }
  LB <- optim(par = initial_parameters,
              fn = function(beta) {
                ea_phi_BLR_DL_vec(beta = beta,
                                  y_labels = y_labels,
                                  X = X,
                                  prior_means = prior_means,
                                  prior_variances = prior_variances,
                                  C = C,
                                  precondition_mat = precondition_mat,
                                  transform_mat = transform_mat)},
              method = "L-BFGS-B",
              lower = lower,
              upper = upper,
              control = list('fnscale' = 1, 'maxit' = 250))$value
  UB <- optim(par = initial_parameters,
              fn = function(beta) {
                ea_phi_BLR_DL_vec(beta = beta,
                                  y_labels = y_labels,
                                  X = X,
                                  prior_means = prior_means,
                                  prior_variances = prior_variances,
                                  C = C,
                                  precondition_mat = precondition_mat,
                                  transform_mat = transform_mat)},
              method = "L-BFGS-B",
              lower = lower,
              upper = upper,
              control = list('fnscale' = -1, 'maxit' = 250))$value
  # multiply the bounds by bounds_multiplier to compensate the case
  # that the opitimser did not get the bounds exactly correct
  # calculate the Bounds Difference
  BD <- UB-LB
  # by subtracting and adding 0.5*(bounds_multiplier-1)*BD
  # makes the resulting bounds difference be larger by a factor of bounds_multiplier
  return(list('LB' = LB - 0.5*(bounds_multiplier-1)*BD,
              'UB' = UB + 0.5*(bounds_multiplier-1)*BD))
}

#' @export
ea_BLR_DL_PT <- function(dim,
                         x0,
                         y,
                         s,
                         t,
                         y_labels,
                         X,
                         prior_means,
                         prior_variances,
                         C,
                         precondition_mat,
                         transform_mats,
                         diffusion_estimator,
                         beta_NB = 10,
                         bounds_multiplier = 1.1,
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
  lbound_Z <- sapply(1:2, function(dim) bes_layers[[dim]]$L)
  ubound_Z <- sapply(1:2, function(dim) bes_layers[[dim]]$U)
  # calculate the lower and upper bounds of phi
  bounds <- lapply(list(lbound_z, ubound_z), function(init) {
    ea_phi_BLR_DL_bounds(initial_parameters = init,
                         y_labels = y_labels,
                         X = X,
                         prior_means = prior_means,
                         prior_variances = prior_variances,
                         C = C,
                         precondition_mat = precondition_mat,
                         transform_mat = transform_mats$to_X,
                         lower = lbound_Z,
                         upper = ubound_Z,
                         bounds_multiplier = bounds_multiplier)})
  LZ <- min(sapply(1:length(bounds), function(i) bounds[[i]]$LB))
  UZ <- max(sapply(1:length(bounds), function(i) bounds[[i]]$UB))
  if (diffusion_estimator=='Poisson') {
    # simulate the number of points to simulate from Possion distribution
    kap <- rpois(n = 1, lambda = (UZ-LZ)*(t-s))
    log_acc_prob <- 0
    if (kap > 0) {
      layered_bb <- layeredBB::multi_layered_brownian_bridge(dim = dim,
                                                             x = z0,
                                                             y = zt,
                                                             s = s,
                                                             t = t,
                                                             bessel_layers = bes_layers,
                                                             times = runif(kap, s, t))
      phi <- ea_phi_BLR_DL_matrix(beta = t(layered_bb$simulated_path[1:dim,]),
                                  y_labels = y_labels,
                                  X = X,
                                  prior_means = prior_means,
                                  prior_variances = prior_variances,
                                  C = C,
                                  precondition_mat = precondition_mat,
                                  transform_mat = transform_mat)
      terms <- (UZ-phi)
      log_acc_prob <- sum(log(terms))
      if (any(terms < 0)) {
        cat('LZ:', LZ, '\n', file = 'bounds.txt', append = T)
        cat('UZ:', UZ, '\n', file = 'bounds.txt', append = T)
        cat('phi:', phi, '\n', file = 'bounds.txt', append = T)
        cat('(UZ-phi):', terms, '\n', file = 'bounds.txt', append = T)
        cat('(phi-LZ):', phi-LZ, '\n', file = 'bounds.txt', append = T)
        stop('Some of (UZ-phi) are < 0. Try increase bounds_multiplier')
      } else if (any((phi - LZ) < 0)) {
        cat('LZ:', LZ, '\n', file = 'bounds.txt', append = T)
        cat('UZ:', UZ, '\n', file = 'bounds.txt', append = T)
        cat('phi:', phi, '\n', file = 'bounds.txt', append = T)
        cat('(UZ-phi):', terms, '\n', file = 'bounds.txt', append = T)
        cat('(phi-LZ):', phi-LZ, '\n', file = 'bounds.txt', append = T)
        stop('Some of (phi-LZ) are < 0. Try increase bounds_multiplier')
      }
    }
    if (logarithm) {
      return(-LZ*(t-s) - kap*log(UZ-LZ) + log_acc_prob)
    } else {
      return(exp(-LZ*(t-s) - kap*log(UZ-LZ) + log_acc_prob))
    }
  } else if (diffusion_estimator=="NB") {
    # integral estimate for gamma in NB estimator
    integral_estimate <- cubature::adaptIntegrate(f = function(s_) {
      ea_phi_BLR_DL_vec(beta = (z0*(t-s_)+zt*s_)/(t-s),
                        y_labels = y_labels,
                        X = X,
                        prior_means = prior_means,
                        prior_variances = prior_variances,
                        C = C,
                        precondition_mat = precondition_mat,
                        transform_mat = transform_mat)},
      lowerLimit = s,
      upperLimit = t)$integral
    gamma_NB <- (t-s)*UZ - integral_estimate
    # simulate the number of points to simulate from Negative Binomial distribution
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
      phi <- ea_phi_BLR_DL_matrix(beta = t(layered_bb$simulated_path[1:dim,]),
                                  y_labels = y_labels,
                                  X = X,
                                  prior_means = prior_means,
                                  prior_variances = prior_variances,
                                  C = C,
                                  precondition_mat = precondition_mat,
                                  transform_mat = transform_mat)
      terms <- (UZ-phi)
      log_acc_prob <- sum(log(terms))
      if (any(terms < 0)) {
        cat('LZ:', LZ, '\n', file = 'bounds.txt', append = T)
        cat('UZ:', UZ, '\n', file = 'bounds.txt', append = T)
        cat('phi:', phi, '\n', file = 'bounds.txt', append = T)
        cat('(UZ-phi):', terms, '\n', file = 'bounds.txt', append = T)
        cat('(phi-LZ):', phi-LZ, '\n', file = 'bounds.txt', append = T)
        stop('Some of (UZ-phi) are < 0. Try increase bounds_multiplier')
      } else if (any((phi - LZ) < 0)) {
        cat('LZ:', LZ, '\n', file = 'bounds.txt', append = T)
        cat('UZ:', UZ, '\n', file = 'bounds.txt', append = T)
        cat('phi:', phi, '\n', file = 'bounds.txt', append = T)
        cat('(UZ-phi):', terms, '\n', file = 'bounds.txt', append = T)
        cat('(phi-LZ):', phi-LZ, '\n', file = 'bounds.txt', append = T)
        stop('Some of (phi-LZ) are < 0. Try increase bounds_multiplier')
      }
    }
    log_middle_term <- kap*log(t-s) + lgamma(beta_NB) + (beta_NB+kap)*log(beta_NB+gamma_NB) -
      lgamma(beta_NB+kap) - beta_NB*log(beta_NB) - kap*log(gamma_NB)
    if (logarithm) {
      return(-UZ*(t-s) + log_middle_term + log_acc_prob)
    } else {
      return(exp(-UZ*(t-s) + log_middle_term + log_acc_prob))
    }
  } else {
    stop("ea_BLR_DL_PT: diffusion_estimator must be set to either \'Poisson\' or \'NB\'")
  }
}

#' @export
combine_data <- function(list_of_data, dim) {
  if (!is.list(list_of_data)) {
    stop("combine_data: list_of_data must be a list")
  } else if (!all(sapply(list_of_data, function(sub_posterior) (is.list(sub_posterior) & identical(names(data), c("y", "X")))))) {
    stop("combine_data: each item in list_of_data must be a list of length 2 with names y and X")
  } else if (!all(sapply(1:m, function(i) is.vector(list_of_data[[i]]$y)))) {
    stop("combine_data: for each i in 1:length(list_of_data), list_of_data[[i]]$y must be a vector")
  } else if (!all(sapply(1:m, function(i) is.matrix(list_of_data[[i]]$X)))) {
    stop("combine_data: for each i in 1:length(list_of_data), list_of_data[[i]]$X must be a matrix")
  } else if (!all(sapply(1:m, function(i) ncol(list_of_data[[i]]$X)!=dim))) {
    stop("combine_data: for each i in 1:length(list_of_data), ncol(list_of_data[[i]]$X) must be equal to dim")
  }
  y <- unlist(lapply(list_of_data, function(sub_posterior) sub_posterior$y))
  X <- do.call(rbind, lapply(list_of_data, function(sub_posterior) sub_posterior$X))
  return(list('y' = y, 'X' = X))
}

#' @export
Q_IS_BLR <- function(particle_set,
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
                     bounds_multiplier = 1.1,
                     seed = NULL,
                     n_cores = parallel::detectCores(),
                     level = 1,
                     node = 1) {
  if (!("particle" %in% class(particle_set))) {
    stop("Q_IS_BLR: particle_set must be a \"particle\" object")
  } else if (!is.list(data_split) | length(data_split)!=m) {
    stop("Q_IS_BLR: data_split must be a list of length m")
  } else if (!all(sapply(data_split, function(sub_posterior) (is.list(sub_posterior) & identical(names(data), c("y", "X")))))) {
    stop("Q_IS_BLR: each item in data_split must be a list of length 2 with names y and X")
  } else if (!all(sapply(1:m, function(i) is.vector(data_split[[i]]$y)))) {
    stop("Q_IS_BLR: for each i in 1:m, data_split[[i]]$y must be a vector")
  } else if (!all(sapply(1:m, function(i) is.matrix(data_split[[i]]$X)))) {
    stop("Q_IS_BLR: for each i in 1:m, data_split[[i]]$X must be a matrix")
  } else if (!all(sapply(1:m, function(i) ncol(data_split[[i]]$X)!=dim))) {
    stop("Q_IS_BLR: for each i in 1:m, ncol(data_split[[i]]$X) must be equal to dim")
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
  }
  transform_matrices <- lapply(1:m, function(c) {
    list('to_Z' = expm::sqrtm(inv_precondition_matrices[[c]]),
         'to_X' = expm::sqrtm(precondition_matrices[[c]]))
  })
  proposal_cov <- calculate_proposal_cov(time = time, weights = inv_precondition_matrices)
  N <- particle_set$N
  # ---------- creating parallel cluster
  cl <- parallel::makeCluster(n_cores, setup_strategy = "sequential", outfile = "smc_fusion_outfile.txt")
  varlist <- c(ls(), list("ea_phi_BLR_DL_matrix",
                          "ea_phi_BLR_DL_bounds",
                          "ea_BLR_DL_PT"))
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
    y_samples <- matrix(nrow = split_N, ncol = 2)
    log_Q_weights <- rep(0, split_N)
    for (i in 1:split_N) {
      y_samples[i,] <- mvrnormArma(N = 1, mu = split_x_means[[core]][i,], Sigma = proposal_cov)
      log_Q_weights[i] <- sum(sapply(1:m, function(c) {
        ea_BLR_DL_PT(dim = dim,
                     x0 = as.vector(split_x_samples[[core]][[i]][c,]),
                     y = as.vector(y_samples[i,]),
                     s = 0,
                     t = time,
                     y_labels = data_split[[c]]$y,
                     X = data_split[[c]]$X,
                     prior_means = prior_means,
                     prior_variances = prior_variances,
                     C = C,
                     precondition_mat = precondition_matrices[[c]],
                     transform_mats = transform_matrices[[c]],
                     diffusion_estimator = diffusion_estimator,
                     beta_NB = beta_NB,
                     bounds_multiplier = bounds_multiplier,
                     logarithm = TRUE)
      }))
      cat('Level:', level, '|| Node:', node, '|| Core:', core, '||', i, '/',
          split_N, '\n', file = 'Q_IS_progress.txt', append = T)
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

#' Parallel Sequential Monte Carlo Fusion for Bayesian Logistic Regression model
#'
#' @param particles_to_fuse list of length m, where particles_to_fuse[c] contains
#'                          the particles for the c-th sub-posterior. Can
#'                          initialise a this from list of sub-posterior samples
#'                          by using the intialise_particle_sets function
#' @param N number of samples
#' @param m number of sub-posteriors to combine
#' @param time time T for fusion algorithm
#' @param dim dimension of the predictors (= p+1)
#' @param data_split list of length m where each item is a list of length 2 where
#'                   for c=1,...,m, data[[c]]$y is the vector for y responses and
#'                   data[[c]]$x is the design matrix for the covariates for sub-posterior c
#' @param prior_means prior for means of predictors
#' @param prior_variances prior for variances of predictors
#' @param C overall number of sub-posteriors
#' @param precondition_matricies list of length m, where precondition_matrices[[c]]
#'                               is the precondition matrix for sub-posterior c
#' @param resampling_method method to be used in resampling, default is multinomial
#'                          resampling ('multi'). Other choices are stratified
#'                          resampling ('strat'), systematic resampling ('system'),
#'                          residual resampling ('resid')
#' @param ESS_threshold number between 0 and 1 defining the proportion of the
#'                      number of samples that ESS needs to be lower than for
#'                      resampling (i.e. resampling is carried out only when
#'                      ESS < N*ESS_threshold)
#' @param diffusion_estimator choice of unbiased estimator for the Exact Algorithm
#'                            between "Poisson" (default) for Poission estimator
#'                            and "NB" for Negative Binomial estimator
#' @param beta_NB beta parameter for Negative Binomial estimator (default 10). If Poisson
#'                estimator used, then this can be ignored
#' @param bounds_multiplier scalar value to mulitply bounds by
#'                          (should greater than or equal to 1)
#' @param seed seed number - default is NULL, meaning there is no seed
#' @param n_cores number of cores to use
#' @param level indicates which level this is for the hierarchy (default 1)
#' @param node indicates which node this is for the hierarchy (default 1)
#'
#' @return A list with components:
#' \describe{
#'   \item{particles}{particles returned from fusion sampler}
#'   \item{proposed_samples}{proposal samples from fusion sampler}
#'   \item{time}{run-time of fusion sampler}
#'   \item{ESS}{list of length (L-1), where ESS[[l]][[i]] is the effective
#'              sample size of the particles after each step BEFORE deciding
#'              whether or not to resample for level l, node i}
#'   \item{CESS}{list of length (L-1), where CESS[[l]][[i]] is the conditional
#'               effective sample size of the particles after each step}
#'   \item{resampled}{list of length (L-1), where resampled[[l]][[i]] is a
#'                    boolean value to record if the particles were resampled
#'                    after each step; rho and Q for level l, node i}
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
                                    diffusion_estimator = 'Poisson',
                                    beta_NB = 10,
                                    bounds_multiplier = 1.1,
                                    seed = NULL,
                                    n_cores = parallel::detectCores(),
                                    level = 1,
                                    node = 1) {
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
  } else if (!all(sapply(data_split, function(sub_posterior) (is.list(sub_posterior) & identical(names(data), c("y", "X")))))) {
    stop("parallel_fusion_SMC_BLR: each item in data_split must be a list of length 2 with names y and X")
  } else if (!all(sapply(1:m, function(i) is.vector(data_split[[i]]$y)))) {
    stop("parallel_fusion_SMC_BLR: for each i in 1:m, data_split[[i]]$y must be a vector")
  } else if (!all(sapply(1:m, function(i) is.matrix(data_split[[i]]$X)))) {
    stop("parallel_fusion_SMC_BLR: for each i in 1:m, data_split[[i]]$X must be a matrix")
  } else if (!all(sapply(1:m, function(i) ncol(data_split[[i]]$X)==dim))) {
    stop("parallel_fusion_SMC_BLR: for each i in 1:m, data_split[[i]]$X must be a matrix with dim columns")
  } else if (!is.vector(prior_means) | length(prior_means)!=dim) {
    stop("parallel_fusion_SMC_BLR: prior_means must be vectors of length dim")
  } else if (!is.vector(prior_variances) | length(prior_variances)!=dim) {
    stop("parallel_fusion_SMC_BLR: prior_variances must be vectors of length dim")
  } else if (!is.list(precondition_matrices) | (length(precondition_matrices)!=m)) {
    stop("parallel_fusion_SMC_BLR: precondition_matrices must be a list of length m")
  } else if (!is.list(inv_precondition_matrices) | (length(inv_precondition_matrices)!=m)) {
    stop("parallel_fusion_SMC_BLR: inv_precondition_matrices must be a list of length m")
  } else if (!(diffusion_estimator %in% c('Poisson', 'NB'))) {
    stop("parallel_fusion_SMC_BLR: diffusion_estimator must be set to either \'Poisson\' or \'NB\'")
  } else if ((ESS_threshold < 0) | (ESS_threshold > 1)) {
    stop("parallel_fusion_SMC_BLR: ESS_threshold must be between 0 and 1")
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
                                   sum_inv_precondition_matrices = inv_sum_matrices(inv_precondition_matrices))
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
  particles <- Q_IS_BLR(particle_set = particles,
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
                        bounds_multiplier = bounds_multiplier,
                        seed = seed,
                        n_cores = n_cores,
                        level = 1,
                        node = 1)
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
                'ESS' = ESS,
                'CESS' = CESS,
                'resampled' = resampled,
                'precondition_matrices' = list(diag(1, dim), precondition_matrices),
                'combined_data' = combine_data(list_of_data = data_split, dim = dim)))
  } else {
    return(list('particles' = particles,
                'proposed_samples' = proposed_samples,
                'ESS' = ESS,
                'CESS' = CESS,
                'resampled' = resampled,
                'precondition_matrices' = list(inv_sum_matrices(inv_precondition_matrices),
                                               precondition_matrices),
                'combined_data' = combine_data(list_of_data = data_split, dim = dim)))
  }
}

#' Hierarchical Monte Carlo Fusion using SMC for Bayesian Logistic Regression model
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
#' @param data_split list of length C where each item is a list of length 2 where
#'                   for c=1,...,C, data[[c]]$y is the vector for y responses and
#'                   data[[c]]$x is the design matrix for the covariates for sub-posterior c
#' @param prior_means prior for means of predictors
#' @param prior_variances prior for variances of predictors
#' @param C number of sub-posteriors at the base level
#' @param precondition logical value determining whether or not a
#'                     preconditioning matrix is to be used
#' @param resampling_method method to be used in resampling, default is
#'                          multinomial resampling ('multi'). Other choices are
#'                          stratified ('strat'), systematic ('system'),
#'                          residual ('resid')
#' @param ESS_threshold number between 0 and 1 defining the proportion
#'                      of the number of samples that ESS needs to be
#'                      lower than for resampling (i.e. resampling is carried
#'                      out only when ESS < N*ESS_threshold)
#' @param diffusion_estimator choice of unbiased estimator for the Exact Algorithm
#'                            between "Poisson" (default) for Poission estimator
#'                            and "NB" for Negative Binomial estimator
#' @param beta_NB beta parameter for Negative Binomial estimator (default 10). If Poisson
#'                estimator used, then this can be ignored
#' @param bounds_multiplier scalar value to mulitply bounds by
#'                          (should greater than or equal to 1)
#' @param seed seed number - default is NULL, meaning there is no seed
#' @param n_cores number of cores to use
#'
#' @return A list with components:
#' \describe{
#'   \item{particles}{list of length (L-1), where particles[[l]][[i]] are the
#'                    particles for level l, node i}
#'   \item{proposed_samples}{list of length (L-1), where proposed_samples[[l]][[i]]
#'                           are the proposed samples for level l, node i}
#'   \item{time}{list of length (L-1), where time[[l]] is the run time
#'                     for level l, node i}
#'   \item{ESS}{list of length (L-1), where ESS[[l]][[i]] is the effective
#'              sample size of the particles after each step BEFORE deciding
#'              whether or not to resample for level l, node i}
#'   \item{CESS}{list of length (L-1), where ESS[[l]][[i]] is the conditional
#'               effective sample size of the particles after each step}
#'   \item{resampled}{list of length (L-1), where resampled[[l]][[i]] is a
#'                    boolean value to record if the particles were resampled
#'                    after each step; rho and Q for level l, node i}
#'   \item{precondition_matrices}{pre-conditioning matrices that were used}
#'   \item{resampling_method}{method that was used in resampling}
#'   \item{y_inputs}{input y data for each level and node}
#'   \item{X_inputs}{input X data for each level and node}
#'   \item{diffusion_times}{vector of length (L-1), where diffusion_times[l]
#'                          are the times for fusion in level l}
#' }
#'
#' @export
hierarchical_fusion_SMC_BLR <- function(N_schedule,
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
                                        bounds_multiplier = 1.1,
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
  } else if (!all(sapply(data_split, function(sub_posterior) (is.list(sub_posterior) & identical(names(data), c("y", "X")))))) {
    stop("hierarchical_fusion_SMC_BLR: each item in data_split must be a list of length 2 with names y and X")
  } else if (!all(sapply(1:C, function(i) is.vector(data_split[[i]]$y)))) {
    stop("hierarchical_fusion_SMC_BLR: for each i in 1:C, data_split[[i]]$y must be a vector")
  } else if (!all(sapply(1:C, function(i) is.matrix(data_split[[i]]$X)))) {
    stop("hierarchical_fusion_SMC_BLR: for each i in 1:C, data_split[[i]]$X must be a matrix")
  } else if (!all(sapply(1:C, function(i) ncol(data_split[[i]]$X)==dim))) {
    stop("hierarchical_fusion_SMC_BLR: for each i in 1:C, data_split[[i]]$X must be a matrix with dim columns")
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
                              diffusion_estimator = diffusion_estimator,
                              beta_NB = beta_NB,
                              bounds_multiplier = bounds_multiplier,
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
    y_inputs[[1]] <- y_inputs[[1]][[1]]
    X_inputs[[1]] <- X_inputs[[1]][[1]]
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
