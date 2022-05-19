#' phi-function for multivariate tempered Gaussian distribution
#'
#' phi-function for the Exact Algorithm for multivariate tempered Gaussian distribution
#'
#' @param dim dimension
#' @param x vector of length dim
#' @param mu vector of length dim for mean
#' @param inv_Sigma dim x dim inverse covariance matrix
#' @param beta real value
#' @param precondition_mat precondition matrix
#'
#' @return real value
#'
#' @export
ea_phi_multiGaussian_DL <- function(dim,
                                    x,
                                    mu,
                                    inv_Sigma,
                                    beta,
                                    precondition_mat) {
  if (is.vector(x)) {
    if (length(x)==dim) {
      return(ea_phi_multiGaussian_DL_vec(x = x,
                                         mu = mu,
                                         inv_Sigma = inv_Sigma,
                                         beta = beta,
                                         precondition_mat = precondition_mat))
    }
  } else if (is.matrix(x)) {
    if (dim(x)[2]==dim) {
      return(ea_phi_multiGaussian_DL_matrix(x = x,
                                            mu = mu,
                                            inv_Sigma = inv_Sigma,
                                            beta = beta,
                                            precondition_mat = precondition_mat))
    }
  }
  stop("ea_phi_multiGaussian_DL: x must be a vector or length dim or a matrix with dim columns")
}

obtain_hypercube_vertices_multiGaussian <- function(bessel_layers,
                                                    mean = rep(0, dim),
                                                    transform_mat = diag(1, dim),
                                                    dim) {
  if (!is.list(bessel_layers)) {
    stop("obtain_hypercube_vertices_multiGaussian: bessel_layers must be a list of length dim")
  } else if (length(mean)!=dim) {
    stop("obtain_hypercube_vertices_multiGaussian: mean must be a mean of length dim")
  }
  if (dim == 1) {
    if (length(bessel_layers)==1) {
      bessel_layers <- bessel_layers[[1]]
    }
    if ((bessel_layers$L < mean) & (bessel_layers$U > mean)) {
      if (abs(bessel_layers$L-mean) > abs(bessel_layers$U-mean)) {
        V <- matrix(c(bessel_layers$L, mean))  
      } else {
        V <- matrix(c(mean, bessel_layers$U))
      }
    } else {
      V <- matrix(c(bessel_layers$L, bessel_layers$U))
    }
    return(list('vertices' = matrix(c(bessel_layers$L, bessel_layers$U)) %*% transform_mat,
                'V' = V %*% transform_mat))
  } else if (dim > 1) {
    if (length(bessel_layers)!=dim) {
      stop("obtain_hypercube_vertices_multiGaussian: if dim > 1, bessel_layers must be a list of length dim")
    }
    bounds <- lapply(1:dim, function(d) c(bessel_layers[[d]]$L, bessel_layers[[d]]$U))
    B <- lapply(1:dim, function(d) {
      if ((bessel_layers[[d]]$L < mean[d]) & (bessel_layers[[d]]$U > mean[d])) {
        if (abs(bessel_layers[[d]]$L-mean[d]) > abs(bessel_layers[[d]]$U-mean[d])) {
          return(c(bessel_layers[[d]]$L, mean[d]))  
        } else {
          return(c(mean[d], bessel_layers[[d]]$U))  
        }
      } else {
        return(c(bessel_layers[[d]]$L, bessel_layers[[d]]$U))
      }
    })
    vertices <- as.matrix(expand.grid(bounds))
    V <- as.matrix(expand.grid(B))
    colnames(vertices) <- c()
    colnames(V) <- c()
    return(list('vertices' = vertices %*% transform_mat,
                'V' = V %*% transform_mat))
  } else {
    stop("obtain_hypercube_vertices_multiGaussian: dim must be greater than or equal to 1")
  }
}

#' Diffusion probability for the Exact Algorithm for Langevin diffusion for
#' multivariate tempered Gaussian distribution
#' 
#' Simulate Langevin diffusion using the Exact Algorithm where pi =
#' multivariate tempered Gaussian distribution
#' 
#' @param x0 start value (vector of length dim)
#' @param y end value (vector of length dim)
#' @param s start time
#' @param t end time
#' @param dim dimension
#' @param mu vector of length dim for mean
#' @param inv_Sigma dim x dim inverse covariance matrix
#' @param beta real value
#' @param precondition_mat precondition matrix (if non-identity matrix, it should
#'                         be the estimated covariance matrix, i.e. a matrix close 
#'                         to solve(inv_Sigma) - could run into problems if this
#'                         is not the case since a trick is used to compute the
#'                         bounds to avoid evaluating phi at 3^d points)
#' @param transform_mats list of transformation matrices where 
#'                       transform_mats$to_Z is the transformation matrix to Z space
#'                       and transform_mats$to_X is the transformation matrix to 
#'                       X space
#' @param diffusion_estimator choice of unbiased estimator for the Exact Algorithm
#'                            between "Poisson" (default) for Poisson estimator
#'                            and "NB" for Negative Binomial estimator
#' @param beta_NB beta parameter for Negative Binomial estimator (default 10)
#' @param gamma_NB_n_points number of points used in the trapezoidal estimation
#'                          of the integral found in the mean of the negative
#'                          binomial estimator (default is 2)
#' @param logarithm logical value to determine if log probability is
#'                  returned (TRUE) or not (FALSE)
#' 
#' @return acceptance probability of simulating Langevin diffusion with pi =
#'         multivariate tempered Gaussian distribution
#' 
#' @export
ea_multiGaussian_DL_PT <- function(x0,
                                   y,
                                   s,
                                   t,
                                   dim,
                                   mu,
                                   inv_Sigma,
                                   inv_Sigma_Z = transform_mats$to_X %*% inv_Sigma %*% transform_mats$to_X,
                                   beta,
                                   precondition_mat,
                                   transform_mats,
                                   diffusion_estimator = 'Poisson',
                                   beta_NB = 10,
                                   gamma_NB_n_points = 2,
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
                                                         mult = 0.1)
  # calculate the lower and upper bounds of phi
  transformed_mean <- transform_mats$to_Z %*% mu
  if (!identical(precondition_mat, diag(1, dim))) {
    # can find bounds in Z space by considering 2^dim points (corners of hypercube in Z)
    hypercube_vertices_Z <- obtain_hypercube_vertices_multiGaussian(bessel_layers = bes_layers,
                                                                    mean = transformed_mean,
                                                                    dim = dim)
    bounds <- ea_phi_multiGaussian_DL_bounds(mu = transformed_mean,
                                             inv_Sigma = inv_Sigma_Z,
                                             beta = beta,
                                             precondition_mat = diag(1, dim),
                                             hypercube_vertices = hypercube_vertices_Z)  
  } else {
    # alternatively, can find bounds by transforming the vertices and considering
    # potentially up to 3^d points
    hypercube_vertices <- obtain_hypercube_vertices(bessel_layers = bes_layers,
                                                    vector = transformed_mean,
                                                    transform_mat = transform_mats$to_X,
                                                    dim = dim)
    mean_in_bes_layer <- all(sapply(1:dim, function(d) {
      (bes_layers[[d]]$L < transformed_mean[d]) & (bes_layers[[d]]$U > transformed_mean[d])
    }))
    bounds <- ea_phi_multiGaussian_DL_bounds(mu = mu,
                                             inv_Sigma = inv_Sigma,
                                             beta = beta,
                                             precondition_mat = precondition_mat,
                                             hypercube_vertices = hypercube_vertices,
                                             mean_in_bes_layer = mean_in_bes_layer)
  }
  LX <- bounds$LB
  UX <- bounds$UB
  PHI <- ea_phi_multiGaussian_DL_LB(mu = mu,
                                    inv_Sigma = inv_Sigma,
                                    beta = beta,
                                    precondition_mat = precondition_mat)
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
      phi <- ea_phi_multiGaussian_DL_matrix(x = sim_path,
                                            mu = mu,
                                            inv_Sigma = inv_Sigma,
                                            beta = beta,
                                            precondition_mat = precondition_mat)
      log_acc_prob <- sum(log(UX-phi))
    }
    if (logarithm) {
      return(-(LX-PHI)*(t-s) - kap*log(UX-LX) + log_acc_prob)
    } else {
      return(exp(-(LX-PHI)*(t-s) - kap*log(UX-LX) + log_acc_prob))
    }
  } else if (diffusion_estimator=='NB') {
    # integral estimate for gamma in NB estimator
    h <- (t-s)/(gamma_NB_n_points-1)
    times_to_eval <- seq(from = s, to = t, by = h)
    integral_estimate <- gamma_NB_multiGaussian(times = times_to_eval,
                                                h = h,
                                                x0 = x0,
                                                y = y,
                                                s = s,
                                                t = t,
                                                dim = dim,
                                                mu = mu,
                                                inv_Sigma = inv_Sigma,
                                                beta = beta,
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
      phi <- ea_phi_multiGaussian_DL_matrix(x = sim_path,
                                            mu = mu,
                                            inv_Sigma = inv_Sigma,
                                            beta = beta,
                                            precondition_mat = precondition_mat)
      log_acc_prob <- sum(log(UX-phi))
    }
    log_middle_term <- kap*log(t-s) + lgamma(beta_NB) + (beta_NB+kap)*log(beta_NB+gamma_NB) -
      lgamma(beta_NB+kap) - beta_NB*log(beta_NB) - kap*log(gamma_NB)
    if (logarithm) {
      return(-(UX-PHI)*(t-s) + log_middle_term + log_acc_prob)
    } else {
      return(exp(-(UX-PHI)*(t-s) + log_middle_term + log_acc_prob))
    }
  } else {
    stop("ea_multiGaussian_DL_PT: diffusion_estimator must be set to either \'Poisson\' or \'NB\'")
  }
}

#' Generalised Monte Carlo Fusion (rejection sampling) [on a single core]
#'
#' Generalised Monte Carlo Fusion with multivariate Gaussian target
#'
#' @param N number of samples
#' @param m number of sub-posteriors to combine
#' @param time time T for fusion algorithm
#' @param samples_to_fuse list of length m, where samples_to_fuse[c] containing
#'                        the samples for the c-th sub-posterior
#' @param dim dimension
#' @param mu vector of length dim for mean
#' @param inv_Sigma dim x dim inverse covariance matrix
#' @param betas vector of length m, where betas[c] is the inverse temperature (beta)
#'              for c-th sub-posterior (can also pass in one number if they are all
#'              at the same inverse temperature)
#' @param precondition_matrices list of length m, where precondition_matrices[[c]]
#'                               is the precondition matrix for sub-posterior c
#'
#' @return samples: fusion samples
#' @return iters_rho: number of iterations from the first accept/reject step (rho)
#' @return iters_Q: number of iterations from the second (diffusion) accept/reject step (Q)
#'
#' @export
fusion_multiGaussian <- function(N,
                                 m,
                                 time,
                                 samples_to_fuse,
                                 dim,
                                 mu,
                                 inv_Sigma,
                                 betas,
                                 precondition_matrices,
                                 inv_precondition_matrices) {
  if (!is.list(samples_to_fuse) | (length(samples_to_fuse)!=m)) {
    stop("fusion_multiGaussian: samples_to_fuse must be a list of length m")
  } else if (!all(sapply(samples_to_fuse, is.matrix)) | !all(sapply(samples_to_fuse, function(core) ncol(core)==dim))) {
    stop("fusion_multiGaussian: the sub-posteriors in samples_to_fuse must be matrices with two columns")
  } else if (!is.vector(mu) | (length(mu)!=dim)) {
    stop("fusion_multiGaussian: mu must be a vector of length dim")
  } else if (!is.matrix(inv_Sigma)) {
    if (nrow(inv_Sigma)!=dim | ncol(inv_Sigma)!=dim) {
      stop("fusion_multiGaussian: inv_Sigma must be a dim x dim matrix")
    }
  } else if (!is.vector(betas) | (length(betas)!=m)) {
    stop("fusion_multiGaussian: betas must be a vector of length m")
  } else if (!is.list(precondition_matrices) | (length(precondition_matrices)!=m)) {
    stop("fusion_multiGaussian: precondition_matrices must be a list of length m")
  } 
  inv_precondition_matrices <- lapply(precondition_matrices, solve)
  transform_matrices <- lapply(1:m, function(c) {
    list('to_Z' = expm::sqrtm(inv_precondition_matrices[[c]]),
         'to_X' = expm::sqrtm(precondition_matrices[[c]]))
  })
  inv_Sigma_Z <- lapply(1:m, function(c) transform_matrices[[c]]$to_X %*% inv_Sigma %*% transform_matrices[[c]]$to_X)
  fusion_samples <- matrix(data = NA, nrow = N, ncol = dim)
  i <- 0; iters_rho <- 0; iters_Q <- 0
  proposal_cov <- calculate_proposal_cov(time = time, weights = inv_precondition_matrices)
  inverse_sum_inv_precondition_matrices <- inverse_sum_matrices(matrices = inv_precondition_matrices)
  while (i < N) {
    iters_rho <- iters_rho + 1
    x <- t(sapply(samples_to_fuse, function(core) core[sample(nrow(core), 1),]))
    x_mean <- weighted_mean_multivariate(matrix = x,
                                         weights = inv_precondition_matrices,
                                         inverse_sum_weights = inverse_sum_inv_precondition_matrices)
    log_rho_prob <- log_rho_multivariate(x = x,
                                         x_mean = x_mean,
                                         time = time,
                                         inv_precondition_matrices = inv_precondition_matrices)
    if (log(runif(1, 0, 1)) < log_rho_prob) {
      iters_Q <- iters_Q + 1
      y <- as.vector(mvrnormArma(N = 1, mu = x_mean, Sigma = proposal_cov))
      log_Q_prob <- sum(sapply(1:m, function(c) {
        ea_multiGaussian_DL_PT(x0 = x[c,],
                               y = y,
                               s = 0,
                               t = time,
                               dim = dim,
                               mu = mu,
                               inv_Sigma = inv_Sigma,
                               inv_Sigma_Z = inv_Sigma_Z[[c]],
                               beta = betas[c],
                               precondition_mat = precondition_matrices[[c]],
                               transform_mats = transform_matrices[[c]],
                               diffusion_estimator = 'Poisson',
                               logarithm = TRUE)
      }))
      if (log(runif(1, 0, 1)) < log_Q_prob) {
        i <- i+1
        fusion_samples[i,] <- y
      }
    }
  }
  return(list('samples' = fusion_samples,
              'iters_rho' = iters_rho,
              'iters_Q' = iters_Q))
}

#' Generalised Monte Carlo Fusion (rejection sampling) [parallel]
#'
#' Generalised Monte Carlo Fusion with multivariate Gaussian target
#'
#' @param N number of samples
#' @param m number of sub-posteriors to combine
#' @param time time T for fusion algorithm
#' @param samples_to_fuse list of length m, where samples_to_fuse[[c]]
#'                        contains the samples for the c-th sub-posterior
#' @param dim dimension
#' @param mu vector of length dim for mean
#' @param Sigma dim x dim covariance matrix
#' @param betas vector of length c, where betas[c] is the inverse temperature 
#'              value for c-th posterior
#' @param precondition_matrices list of length m, where precondition_matrices[[c]]
#'                               is the precondition matrix for sub-posterior c
#' @param seed seed number - default is NULL, meaning there is no seed
#' @param n_cores number of cores to use
#'
#' @return A list with components:
#' \describe{
#'  \item{samples}{fusion samples}
#'  \item{rho}{acceptance rate for rho step}
#'  \item{Q}{acceptance rate for Q step}
#'  \item{rhoQ}{overall acceptance rate}
#'  \item{time}{run-time of fusion sampler}
#'  \item{rho_iterations}{number of iterations for rho step}
#'  \item{Q_iterations}{number of iterations for Q step}
#'  \item{precondition_matrices}{list of length dim where precondition_matrices[[2]] 
#'                               are the pre-conditioning matrices that were used 
#'                               and precondition_matrices[[1]] are the combined 
#'                               precondition matrices}
#' }
#'
#' @export
parallel_fusion_multiGaussian <- function(N,
                                          m,
                                          time,
                                          samples_to_fuse,
                                          dim,
                                          mu,
                                          Sigma,
                                          betas, 
                                          precondition_matrices,
                                          seed = NULL,
                                          n_cores = parallel::detectCores()) {
  if (!is.list(samples_to_fuse) | (length(samples_to_fuse)!=m)) {
    stop("parallel_fusion_multiGaussian: samples_to_fuse must be a list of length m")
  } else if (!all(sapply(samples_to_fuse, is.matrix))) {
    stop("parallel_fusion_multiGaussian: the sub-posteriors in samples_to_fuse must be matrices")
  } else if (!all(sapply(samples_to_fuse, function(core) ncol(core)==dim))) {
    stop("parallel_fusion_multiGaussian: the sub-posteriors in samples_to_fuse must be matrices with dim columns")
  } else if (!is.vector(mu) | (length(mu)!=dim)) {
    stop("parallel_fusion_multiGaussian: mu must be a vector of length dim")
  } else if (!is.matrix(Sigma)) {
    if (nrow(Sigma)!=dim | ncol(Sigma)!=dim) {
      stop("parallel_fusion_multiGaussian: Sigma must be a dim x dim matrix")
    }
  } else if (!is.vector(betas) | (length(betas)!=m)) {
    stop("parallel_fusion_multiGaussian: betas must be a vector of length m")
  } else if (!is.list(precondition_matrices) | (length(precondition_matrices)!=m)) {
    stop("parallel_fusion_multiGaussian: precondition_matrices must be a list of length m")
  }
  inv_Sigma <- solve(Sigma)
  ######### creating parallel cluster
  cl <- parallel::makeCluster(n_cores)
  parallel::clusterExport(cl, envir = environment(), varlist = ls())
  parallel::clusterExport(cl, varlist = ls("package:DCFusion"))
  parallel::clusterExport(cl, varlist = ls("package:layeredBB"))
  if (!is.null(seed)) {
    parallel::clusterSetRNGStream(cl, iseed = seed)
  }
  # how many samples do we need for each core?
  if (N < n_cores) {
    samples_per_core <- rep(1, N)
  } else {
    samples_per_core <- rep(floor(N/n_cores), n_cores)
    if (sum(samples_per_core)!=N) {
      remainder <- N %% n_cores
      samples_per_core[1:remainder] <- samples_per_core[1:remainder] + 1
    }
  }
  # run fusion in parallel
  pcm <- proc.time()
  fused <- parallel::parLapply(cl, X = 1:length(samples_per_core), fun = function(core) {
    fusion_multiGaussian(N = samples_per_core[core],
                         m = m,
                         time = time,
                         samples_to_fuse = samples_to_fuse,
                         dim = dim,
                         mu = mu,
                         inv_Sigma = inv_Sigma,
                         betas = betas,
                         precondition_matrices = precondition_matrices)
  })
  final <- proc.time() - pcm
  parallel::stopCluster(cl)
  # ---------- return samples and acceptance probabilities
  samples <- do.call(rbind, lapply(1:length(samples_per_core), function(i) fused[[i]]$samples))
  rho_iterations <- sum(sapply(1:length(samples_per_core), function(i) fused[[i]]$iters_rho))
  Q_iterations <- sum(sapply(1:length(samples_per_core), function(i) fused[[i]]$iters_Q))
  rho_acc <- Q_iterations / rho_iterations
  Q_acc <- N / Q_iterations
  rhoQ_acc <- N / rho_iterations
  if (identical(precondition_matrices, rep(list(diag(1, dim)), m))) {
    new_precondition_matrices <- list(diag(1, dim), precondition_matrices)
  } else {
    new_precondition_matrices <- list(inverse_sum_matrices(lapply(precondition_matrices, solve)),
                                      precondition_matrices)
  }
  return(list('samples' = samples,
              'rho' = rho_acc,
              'Q' = Q_acc,
              'rhoQ '= rhoQ_acc,
              'time' = final['elapsed'],
              'rho_iterations' = rho_iterations,
              'Q_iterations' = Q_iterations,
              'precondition_matrices' = new_precondition_matrices))
}

#' (Balanced Binary) D&C Monte Carlo Fusion (rejection sampling)
#'
#' (Balanced Binary) D&C Monte Carlo Fusion with multivariate Gaussian target
#'
#' @param N_schedule vector of length (L-1), where N_schedule[l] is the number 
#'                   of samples per node at level l
#' @param m_schedule vector of length (L-1), where m_schedule[l] is the number 
#'                   of samples to fuse for level l
#' @param time_schedule vector of length(L-1), where time_schedule[l] is the 
#'                      time chosen for Fusion at level l
#' @param base_samples list of length (1/start_beta), where base_samples[[c]] 
#'                     contains the samples for the c-th node in the level
#' @param L total number of levels in the hierarchy
#' @param dim dimension
#' @param mu vector of length dim for mean
#' @param Sigma dim x dim covariance matrix
#' @param start_beta beta for the base level
#' @param precondition either a logical value to determine if preconditioning matrices are
#'                     used (TRUE - and is set to be the variance of the sub-posterior samples)
#'                     or not (FALSE - and is set to be the identity matrix for all sub-posteriors),
#'                     or a list of length (1/start_beta) where precondition[[c]]
#'                     is the preconditioning matrix for sub-posterior c. Default is TRUE
#' @param seed seed number - default is NULL, meaning there is no seed
#' @param n_cores number of cores to use
#' 
#' @return A list with components:
#' \describe{
#'  \item{samples}{list of length (L-1), where samples[[l]][[i]] are the samples 
#'                 for level l, node i}
#'  \item{time}{list of length (L-1), where time[[l]][[i]] is the run time for 
#'              level l, node i}
#'  \item{rho_acc}{list of length (L-1), where rho_acc[[l]][i] is the acceptance 
#'                 rate for first fusion step for level l, node i}
#'  \item{Q_acc}{list of length (L-1), where Q_acc[[l]][i] is the acceptance 
#'               rate for second fusion step for level l, node i}
#'  \item{rhoQ_acc}{list of length (L-1), where rhoQ_acc[[l]][i] is the overall 
#'                  acceptance rate for fusion for level l, node i}
#'  \item{diffusion_times}{vector of length (L-1), where diffusion_times[l] are 
#'                         the times for fusion in level l (= time_schedule)}
#'  \item{precondition_matrices}{preconditioning matrices used in the algorithm 
#'                               for each node}
#'  \item{overall_rho}{vector of length (L-1), where overall_rho[k] is overall 
#'                     acceptance rate for rho of level l}
#'  \item{overall_Q}{vector of length (L-1), where overall_Q[k] is overall 
#'                   acceptance rate for Q of level l}
#'  \item{overall_rhoQ}{vector of length (L-1), where overall_rhoQ[k] is overall 
#'                      acceptance rate for rho*Q of level l}
#'  \item{overall_time}{vector of length (L-1), where overall_time[k] is 
#'                      overall taken for level l}
#' }
#' 
#' @export
bal_binary_fusion_multiGaussian <- function(N_schedule,
                                            m_schedule,
                                            time_schedule,
                                            base_samples,
                                            L,
                                            dim,
                                            mu,
                                            Sigma,
                                            start_beta,
                                            precondition = TRUE,
                                            seed = NULL,
                                            n_cores = parallel::detectCores()) {
  if (!is.vector(N_schedule) | (length(N_schedule)!=(L-1))) {
    stop("bal_binary_fusion_multiGaussian: N_schedule must be a vector of length (L-1)")
  } else if (!is.vector(m_schedule) | (length(m_schedule)!=(L-1))) {
    stop("bal_binary_fusion_multiGaussian: m_schedule must be a vector of length (L-1)")
  } else if (!is.vector(time_schedule) | (length(time_schedule)!=(L-1))) {
    stop("bal_binary_fusion_multiGaussian: time_schedule must be a vector of length (L-1)")
  } else if (!is.list(base_samples) | (length(base_samples)!=(1/start_beta))) {
    stop("bal_binary_fusion_multiGaussian: base_samples must be a list of length (1/start_beta)")
  } else if (!all(sapply(base_samples, is.matrix))) {
    stop("bal_binary_fusion_multiGaussian: the sub-posterior samples in base_samples must be matrices")
  } else if (!all(sapply(base_samples, function(core) ncol(core)==dim))) {
    stop("bal_binary_fusion_multiGaussian: the sub-posterior samples in base_samples must be matrices with dim columns")
  } else if (!is.vector(mu) | (length(mu)!=dim)) {
    stop("bal_binary_fusion_multiGaussian: mu must be a vector of length dim")
  } else if (!is.matrix(Sigma)) {
    if (nrow(Sigma)!=dim | ncol(Sigma)!=dim) {
      stop("bal_binary_fusion_multiGaussian: Sigma must be a dim x dim matrix")
    }
  }
  if (is.vector(m_schedule) & (length(m_schedule)==(L-1))) {
    for (l in (L-1):1) {
      if (((1/start_beta)/prod(m_schedule[(L-1):l]))%%1!=0) {
        stop("bal_binary_fusion_multiGaussian: check that (1/start_beta)/prod(m_schedule[(L-1):l])
              is an integer for l=L-1,...,1")
      }
    }
  } else {
    stop("bal_binary_fusion_multiGaussian: m_schedule must be a vector of length (L-1)")
  }
  m_schedule <- c(m_schedule, 1)
  hier_samples <- list()
  hier_samples[[L]] <- base_samples
  time <- list()
  rho <- list()
  Q <- list()
  rhoQ <- list()
  overall_rho <- rep(0, L-1)
  overall_Q <- rep(0, L-1)
  overall_rhoQ <- rep(0, L-1)
  overall_time <- rep(0, L-1)
  precondition_matrices <- list()
  if (is.logical(precondition)) {
    if (precondition) {
      precondition_matrices[[L]] <- lapply(base_samples, cov)
    } else {
      precondition_matrices[[L]] <- rep(list(diag(1, dim)), (1/start_beta))
    }
  } else if (is.list(precondition)) {
    if (length(precondition)==(1/start_beta) & all(sapply(precondition, is.matrix))) {
      if (all(sapply(precondition, function(sub) ncol(sub)==dim))) {
        precondition_matrices[[L]] <- precondition  
      }
    }
  } else {
    stop("bal_binary_fusion_multiGaussian: precondition must be a logical indicating 
          whether or not a preconditioning matrix should be used, or a list of
          length C, where precondition[[c]] is the preconditioning matrix for
          the c-th sub-posterior")
  }
  cat('Starting bal_binary fusion \n', file = 'bal_binary_fusion_multiGaussian.txt')
  for (k in ((L-1):1)) {
    n_nodes <- max((1/start_beta)/prod(m_schedule[L:k]), 1)
    cat('########################\n', file = 'bal_binary_fusion_multiGaussian.txt', 
        append = T)
    cat('Starting to fuse', m_schedule[k], 'densities of pi^beta, where beta =', 
        prod(m_schedule[L:(k+1)]), '/', (1/start_beta), 'for level', k, 'with time', 
        time_schedule[k], ', which is using', parallel::detectCores(), 'cores\n',
        file = 'bal_binary_fusion_multiGaussian.txt', append = T)
    cat('There are', n_nodes, 'nodes at this level each giving', N_schedule[k],
        'samples for beta =', prod(m_schedule[L:k]), '/', (1/start_beta),
        '\n', file = 'bal_binary_fusion_multiGaussian.txt', append = T)
    cat('########################\n', file = 'bal_binary_fusion_multiGaussian.txt', 
        append = T)
    fused <- lapply(X = 1:n_nodes, FUN = function(i) {
      previous_nodes <- ((m_schedule[k]*i)-(m_schedule[k]-1)):(m_schedule[k]*i)
      precondition_mats <- precondition_matrices[[k+1]][previous_nodes]
      parallel_fusion_multiGaussian(N = N_schedule[k], 
                                    m = m_schedule[k], 
                                    time = time_schedule[k],
                                    samples_to_fuse = hier_samples[[k+1]][previous_nodes],
                                    dim = dim,
                                    mu = mu,
                                    Sigma = Sigma,
                                    betas = rep(prod(m_schedule[L:(k+1)])*(start_beta), m_schedule[k]),
                                    precondition_matrices = precondition_mats,
                                    seed = seed,
                                    n_cores = n_cores)
    })
    cat('Completed fusion for level', k, '\n', file = 'bal_binary_fusion_multiGaussian.txt', 
        append = T)
    cat('########################\n', file = 'bal_binary_fusion_multiGaussian.txt', 
        append = T)
    # need to combine the correct samples
    hier_samples[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$samples)
    rho[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$rho)
    Q[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$Q)
    rhoQ[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$rhoQ)
    time[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$time)
    sum_rho_iterations <- sum(unlist(lapply(1:n_nodes, function(i) fused[[i]]$rho_iterations)))
    sum_Q_iterations <- sum(unlist(lapply(1:n_nodes, function(i) fused[[i]]$Q_iterations)))
    overall_rho[k] <- sum_Q_iterations / sum_rho_iterations
    overall_Q[k] <- N_schedule[k]*n_nodes / sum_Q_iterations
    overall_rhoQ[k] <- N_schedule[k]*n_nodes / sum_rho_iterations
    overall_time[k] <- sum(unlist(time[[k]]))
    precondition_matrices[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$precondition_matrices[[1]])
  }
  cat('Completed bal_binary fusion\n', file = 'bal_binary_fusion_multiGaussian.txt', append = T)
  if (length(hier_samples[[1]])==1) {
    hier_samples[[1]] <- hier_samples[[1]][[1]]
    time[[1]] <- time[[1]][[1]]
    rho[[1]] <- rho[[1]][[1]]
    Q[[1]] <- Q[[1]][[1]]
    rhoQ[[1]] <- rhoQ[[1]][[1]]
    precondition_matrices[[1]] <- precondition_matrices[[1]][[1]]
  }
  return(list('samples' = hier_samples,
              'time' = time,
              'rho_acc' = rho,
              'Q_acc' = Q,
              'rhoQ_acc' = rhoQ,
              'diffusion_times' = time_schedule,
              'precondition_matrices' = precondition_matrices,
              'overall_rho' = overall_rho,
              'overall_Q' = overall_Q,
              'overall_rhoQ' = overall_rhoQ,
              'overall_time' = overall_time))
}

#' (Progressive) D&C Monte Carlo Fusion (rejection sampling)
#'
#' (Progressive) D&C Monte Carlo Fusion with multivariate Gaussian target
#'
#' @param N_schedule vector of length (L-1), where N_schedule[l] is the number 
#'                   of samples per node at level l
#' @param time_schedule vector of length(L-1), where time_schedule[l] is the 
#'                      time chosen for Fusion at level l
#' @param dim dimension
#' @param mu vector of length dim for mean
#' @param Sigma dim x dim covariance matrix
#' @param start_beta beta for the base level
#' @param base_samples list of length (1/start_beta), where base_samples[[c]] 
#'                     contains the samples for the c-th node in the level
#' @param precondition either a logical value to determine if preconditioning matrices are
#'                     used (TRUE - and is set to be the variance of the sub-posterior samples)
#'                     or not (FALSE - and is set to be the identity matrix for all sub-posteriors),
#'                     or a list of length (1/start_beta) where precondition[[c]]
#'                     is the preconditioning matrix for sub-posterior c. Default is TRUE
#' @param seed seed number - default is NULL, meaning there is no seed
#' @param n_cores number of cores to use
#'
#' @return A list with components:
#' \describe{
#'  \item{samples}{list of length (L-1), where samples[[l]][[i]] are the samples 
#'                 for level l, node i}
#'  \item{time}{list of length (L-1), where time[[l]][[i]] is the run time for 
#'              level l, node i}
#'  \item{rho_acc}{list of length (L-1), where rho_acc[[l]][i] is the acceptance 
#'                 rate for first fusion step for level l, node i}
#'  \item{Q_acc}{list of length (L-1), where Q_acc[[l]][i] is the acceptance 
#'               rate for second fusion step for level l, node i}
#'  \item{rhoQ_acc}{list of length (L-1), where rhoQ_acc[[l]][i] is the overall 
#'                  acceptance rate for fusion for level l, node i}
#'  \item{diffusion_times}{vector of length (L-1), where diffusion_times[l] are 
#'                         the times for fusion in level l (= time_schedule)}
#'  \item{precondition_matrices}{preconditioning matrices used in the algorithm 
#'                               for each node}
#' }
#'
#' @export
progressive_fusion_multiGaussian <- function(N_schedule, 
                                             time_schedule,
                                             dim,
                                             mu,
                                             Sigma,
                                             start_beta,
                                             base_samples,
                                             precondition = TRUE,
                                             seed = NULL,
                                             n_cores = parallel::detectCores()) {
  if (!is.vector(N_schedule) | (length(N_schedule)!=(1/start_beta)-1)) {
    stop("progressive_fusion_multiGaussian: N_schedule must be a vector of length ((1/start_beta)-1)")
  } else if (!is.vector(time_schedule) | (length(time_schedule)!=(1/start_beta)-1)) {
    stop("progressive_fusion_multiGaussian: time_schedule must be a vector of length ((1/start_beta)-1)")
  } else if (!is.list(base_samples) | (length(base_samples)!=(1/start_beta))) {
    stop("progressive_fusion_multiGaussian: base_samples must be a list of length (1/start_beta)")
  } else if (!all(sapply(base_samples, is.matrix))) {
    stop("progressive_fusion_multiGaussian: the sub-posterior samples in base_samples must be matrices")
  } else if (!all(sapply(base_samples, function(core) ncol(core)==dim))) {
    stop("progressive_fusion_multiGaussian: the sub-posterior samples in base_samples must be matrices with dim columns")
  } else if (!is.vector(mu) | (length(mu)!=dim)) {
    stop("progressive_fusion_multiGaussian: mu must be a vector of length dim")
  } else if (!is.matrix(Sigma)) {
    if (nrow(Sigma)!=dim | ncol(Sigma)!=dim) {
      stop("progressive_fusion_multiGaussian: Sigma must be a dim x dim matrix")
    }
  }
  prog_samples <- list()
  prog_samples[[(1/start_beta)]] <- base_samples
  time <- rep(0, (1/start_beta)-1)
  rho <- rep(0, (1/start_beta)-1)
  Q <- rep(0, (1/start_beta)-1)
  rhoQ <- rep(0, (1/start_beta)-1)
  precondition_matrices <- list()
  if (is.logical(precondition)) {
    if (precondition) {
      precondition_matrices[[(1/start_beta)]] <- lapply(base_samples, cov)
    } else {
      precondition_matrices[[(1/start_beta)]] <- rep(list(diag(1, dim)), (1/start_beta))
    }
  } else if (is.list(precondition)) {
    if (length(precondition)==(1/start_beta) & all(sapply(precondition, is.matrix))) {
      if (all(sapply(precondition, function(sub) ncol(sub)==dim))) {
        precondition_matrices[[(1/start_beta)]] <- precondition  
      }
    }
  } else {
    stop("progressive_fusion_multiGaussian: precondition must be a logical indicating 
          whether or not a preconditioning matrix should be used, or a list of
          length C, where precondition[[c]] is the preconditioning matrix for
          the c-th sub-posterior")
  }
  index <- 2
  cat('Starting progressive fusion \n', file = 'progressive_fusion_multiGaussian.txt')
  for (k in ((1/start_beta)-1):1) {
    if (k==(1/start_beta)-1) {
      cat('########################\n', file = 'progressive_fusion_multiGaussian.txt', 
          append = T)
      cat('Starting to fuse', 2, 'densities for level', k, 'which is using', 
          parallel::detectCores(), 'cores\n',
          file = 'progressive_fusion_multiGaussian.txt', append = T)
      cat('Fusing samples for beta =', 1, '/', (1/start_beta), 'with time', 
          time_schedule[k], 'to get', N_schedule[k], 'samples for beta =', 2, 
          '/', (1/start_beta), '\n', file = 'progressive_fusion_multiGaussian.txt', 
          append = T)
      cat('########################\n', file = 'progressive_fusion_multiGaussian.txt', 
          append = T)
      samples_to_fuse <- list(base_samples[[1]], base_samples[[2]])
      precondition_mats <- precondition_matrices[[k+1]][1:2]
      fused <- parallel_fusion_multiGaussian(N = N_schedule[k], 
                                             m = 2,
                                             time = time_schedule[k],
                                             samples_to_fuse = samples_to_fuse,
                                             dim = dim,
                                             mu = mu,
                                             Sigma = Sigma,
                                             betas = c(start_beta, start_beta),
                                             precondition_matrices = precondition_mats,
                                             seed = seed,
                                             n_cores = n_cores)
    } else {
      cat('########################\n', file = 'progressive_fusion_multiGaussian.txt',
          append = T)
      cat('Starting to fuse', 2, 'densities for level', k, 'which is using', 
          parallel::detectCores(), 'cores\n', file = 'progressive_fusion_multiGaussian.txt',
          append = T)
      cat('Fusing samples for beta =', index, '/', (1/start_beta), 'and beta =', 
          1, '/', (1/start_beta), 'with time', time_schedule[k], 'to get', 
          N_schedule[k], 'samples for beta =', (index+1), '/', (1/start_beta),
          '\n', file = 'progressive_fusion_multiGaussian.txt', append = T)
      cat('########################\n', file = 'progressive_fusion_multiGaussian.txt', 
          append = T)
      samples_to_fuse <- list(prog_samples[[k+1]], 
                              base_samples[[index+1]])
      precondition_mats <- list(precondition_matrices[[k+1]],
                                precondition_matrices[[(1/start_beta)]][[index+1]]) 
      fused <- parallel_fusion_multiGaussian(N = N_schedule[k], 
                                             m = 2,
                                             time = time_schedule[k],
                                             samples_to_fuse = samples_to_fuse,
                                             dim = dim,
                                             mu = mu,
                                             Sigma = Sigma,
                                             betas = c(index*start_beta, start_beta),
                                             precondition_matrices = precondition_mats,
                                             seed = seed, 
                                             n_cores = n_cores)
      index <- index + 1
    }
    # need to combine the correct samples
    prog_samples[[k]] <- fused$samples
    precondition_matrices[[k]] <- fused$precondition_matrices[[1]]
    rho[k] <- fused$rho
    Q[k] <- fused$Q
    rhoQ[k] <- fused$rhoQ
    time[k] <- fused$time
  }
  cat('Completed progressive fusion\n', file = 'progressive_fusion_multiGaussian.txt', append = T)
  return(list('samples' = prog_samples,
              'time' = time,
              'rho_acc' = rho,
              'Q_acc' = Q,
              'rhoQ_acc' = rhoQ,
              'diffusion_times' = time_schedule,
              'precondition_matrices' = precondition_matrices))
}

#' Q Importance Sampling Step
#'
#' Q Importance Sampling weighting for multivariate Gaussian distributions
#'
#' @param particle_set particles set prior to Q importance sampling step
#' @param m number of sub-posteriors to combine
#' @param time time T for fusion algorithm
#' @param dim dimension
#' @param mu vector of length dim for mean
#' @param inv_Sigma dim x dim inverse covariance matrix
#' @param betas vector of length c, where betas[c] is the inverse temperature 
#'              value for c-th posterior
#' @param precondition_matrices list of length m, where precondition_matrices[[c]]
#'                               is the precondition matrix for sub-posterior c
#' @param inv_precondition_matrices list of length m, where inv_precondition_matrices[[c]]
#'                                  is the inverse precondition matrix for sub-posterior c
#' @param diffusion_estimator choice of unbiased estimator for the Exact Algorithm
#'                            between "Poisson" (default) for Poisson estimator
#'                            and "NB" for Negative Binomial estimator
#' @param beta_NB beta parameter for Negative Binomial estimator (default 10)
#' @param gamma_NB_n_points number of points used in the trapezoidal estimation
#'                          of the integral found in the mean of the negative
#'                          binomial estimator (default is 2)
#' @param seed seed number - default is NULL, meaning there is no seed
#' @param n_cores number of cores to use
#'
#' @return An updated particle set
#' 
#' @export
Q_IS_multiGaussian <- function(particle_set,
                               m,
                               time,
                               dim,
                               mu,
                               inv_Sigma,
                               betas,
                               precondition_matrices,
                               inv_precondition_matrices,
                               diffusion_estimator = 'Poisson',
                               beta_NB = 10,
                               gamma_NB_n_points = 2,
                               seed = NULL,
                               n_cores = parallel::detectCores()) {
  if (!("particle" %in% class(particle_set))) {
    stop("Q_IS_multiGaussian: particle_set must be a \"particle\" object")
  } else if (!is.vector(mu) | (length(mu)!=dim)) {
    stop("Q_IS_multiGaussian: mu must be a vector of length dim")
  } else if (!is.matrix(inv_Sigma)) {
    if (nrow(inv_Sigma)!=dim | ncol(inv_Sigma)!=dim) {
      stop("Q_IS_multiGaussian: inv_Sigma must be a dim x dim matrix")
    }
  } else if (!is.vector(betas) | (length(betas)!=m)) {
    stop("Q_IS_multiGaussian: betas must be a vector of length m")
  } else if (!is.list(precondition_matrices) | (length(precondition_matrices)!=m)) {
    stop("Q_IS_multiGaussian: precondition_matrices must be a list of length m")
  } else if (!is.list(inv_precondition_matrices) | (length(inv_precondition_matrices)!=m)) {
    stop("Q_IS_multiGaussian: inv_precondition_matrices must be a list of length m")
  }
  transform_matrices <- lapply(1:m, function(c) {
    list('to_Z' = expm::sqrtm(inv_precondition_matrices[[c]]),
         'to_X' = expm::sqrtm(precondition_matrices[[c]]))
  })
  inv_Sigma_Z <- lapply(1:m, function(c) transform_matrices[[c]]$to_X %*% inv_Sigma %*% transform_matrices[[c]]$to_X)
  proposal_cov <- calculate_proposal_cov(time = time, weights = inv_precondition_matrices)
  N <- particle_set$N
  # ---------- creating parallel cluster
  cl <- parallel::makeCluster(n_cores)
  parallel::clusterExport(cl, envir = environment(), varlist = ls())
  parallel::clusterExport(cl, varlist = ls("package:DCFusion"))
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
    for (i in 1:split_N) {
      y_samples[i,] <- mvrnormArma(N = 1, mu = split_x_means[[core]][i,], Sigma = proposal_cov)
      log_Q_weights[i] <- sum(sapply(1:m, function(c) {
        ea_multiGaussian_DL_PT(x0 = as.vector(split_x_samples[[core]][[i]][c,]),
                               y = as.vector(y_samples[i,]),
                               s = 0,
                               t = time,
                               dim = dim,
                               mu = mu,
                               inv_Sigma = inv_Sigma,
                               inv_Sigma_Z = inv_Sigma_Z[[c]],
                               beta = betas[c],
                               precondition_mat = precondition_matrices[[c]],
                               transform_mats = transform_matrices[[c]],
                               diffusion_estimator = diffusion_estimator,
                               beta_NB = beta_NB,
                               gamma_NB_n_points = gamma_NB_n_points,
                               logarithm = TRUE)
      }))
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
  return(particle_set)
}

#' Generalised Monte Carlo Fusion [parallel]
#'
#' Generalised Monte Carlo Fusion with multivariate Gaussian target
#'
#' @param particles_to_fuse list of length m, where particles_to_fuse[[c]]
#'                          contains the particles for the c-th sub-posterior
#'                          (a list of particles to fuse can be initialised by 
#'                          initialise_particle_sets() function)
#' @param N number of samples
#' @param m number of sub-posteriors to combine
#' @param time time T for fusion algorithm
#' @param precondition_matrices list of length m, where precondition_matrices[[c]]
#'                               is the precondition matrix for sub-posterior c
#' @param dim dimension
#' @param mu vector of length dim for mean
#' @param Sigma dim x dim covariance matrix
#' @param betas vector of length c, where betas[c] is the inverse temperature 
#'              value for c-th posterior
#' @param resampling_method method to be used in resampling, default is multinomial 
#'                          resampling ('multi'). Other choices are stratified 
#'                          resampling ('strat'), systematic resampling ('system'),
#'                          residual resampling ('resid')
#' @param ESS_threshold number between 0 and 1 defining the proportion 
#'                      of the number of samples that ESS needs to be
#'                      lower than for resampling (i.e. resampling is carried 
#'                      out only when ESS < N*ESS_threshold)
#' @param diffusion_estimator choice of unbiased estimator for the Exact Algorithm
#'                            between "Poisson" (default) for Poisson estimator
#'                            and "NB" for Negative Binomial estimator
#' @param beta_NB beta parameter for Negative Binomial estimator (default 10)
#' @param gamma_NB_n_points number of points used in the trapezoidal estimation
#'                          of the integral found in the mean of the negative
#'                          binomial estimator (default is 2)
#' @param seed seed number - default is NULL, meaning there is no seed
#' @param n_cores number of cores to use
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
#'   \item{precondition_matrices}{list of length dim where precondition_matrices[[2]] 
#'                                are the pre-conditioning matrices that were used 
#'                                and precondition_matrices[[1]] are the combined 
#'                                precondition matrices}
#' }
#' 
#' @export
parallel_fusion_SMC_multiGaussian <- function(particles_to_fuse,
                                              N, 
                                              m,
                                              time,
                                              dim,
                                              mu,
                                              Sigma,
                                              betas,
                                              precondition_matrices, 
                                              resampling_method = 'multi',
                                              ESS_threshold = 0.5,
                                              diffusion_estimator = 'Poisson',
                                              beta_NB = 10,
                                              gamma_NB_n_points = 2,
                                              seed = NULL,
                                              n_cores = parallel::detectCores()) {
  if (!is.list(particles_to_fuse) | (length(particles_to_fuse)!=m)) {
    stop("parallel_fusion_SMC_multiGaussian: particles_to_fuse must be a list of length m")
  } else if (!all(sapply(particles_to_fuse, function(sub_posterior) ("particle" %in% class(sub_posterior))))) {
    stop("parallel_fusion_SMC_multiGaussian: particles in particles_to_fuse must be \"particle\" objects")
  } else if (!all(sapply(particles_to_fuse, function(sub_posterior) is.matrix(sub_posterior$y_samples)))) {
    stop("parallel_fusion_SMC_multiGaussian: the particles' samples for y should all be matrices")
  } else if (!all(sapply(particles_to_fuse, function(sub_posterior) ncol(sub_posterior$y_samples)==dim))) {
    stop("parallel_fusion_SMC_multiGaussian: the particles' samples for y should all be matrices with dim columns")
  } else if (!is.vector(mu) | (length(mu)!=dim)) {
    stop("parallel_fusion_SMC_multiGaussian: mu must be a vector of length dim")
  } else if (!is.matrix(Sigma)) {
    if (nrow(Sigma)!=dim | ncol(Sigma)!=dim) {
      stop("parallel_fusion_SMC_multiGaussian: Sigma must be a dim x dim matrix")
    }
  } else if (!is.vector(betas) | (length(betas)!=m)) {
    stop("parallel_fusion_SMC_multiGaussian: betas must be a vector of length m")
  } else if (!is.list(precondition_matrices) | (length(precondition_matrices)!=m)) {
    stop("parallel_fusion_SMC_multiGaussian: precondition_matrices must be a list of length m")
  } else if ((ESS_threshold < 0) | (ESS_threshold > 1)) {
    stop("parallel_fusion_SMC_multiGaussian: ESS_threshold must be between 0 and 1")
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
                                   n_cores = n_cores)
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
  particles <- Q_IS_multiGaussian(particle_set = particles,
                                  m = m,
                                  time = time,
                                  dim = dim,
                                  mu = mu,
                                  inv_Sigma = solve(Sigma),
                                  betas = betas,
                                  precondition_matrices = precondition_matrices,
                                  inv_precondition_matrices = inv_precondition_matrices,
                                  diffusion_estimator = diffusion_estimator,
                                  beta_NB = beta_NB,
                                  gamma_NB_n_points = gamma_NB_n_points,
                                  seed = seed,
                                  n_cores = n_cores)
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
              'precondition_matrices' = new_precondition_matrices))
}

#' (Balanced Binary) D&C Monte Carlo Fusion using SMC
#'
#' (Balanced Binary) D&C Monte Carlo Fusion with multivariate Gaussian target
#'
#' @param N_schedule vector of length (L-1), where N_schedule[l] is the number 
#'                   of samples per node at level l
#' @param m_schedule vector of length (L-1), where m_schedule[l] is the number 
#'                   of samples to fuse for level l
#' @param time_schedule vector of length(L-1), where time_schedule[l] is the time 
#'                      chosen for Fusion at level l
#' @param base_samples list of length (1/start_beta), where base_samples[[c]] 
#'                     contains the samples for the c-th node in the level
#' @param L total number of levels in the hierarchy
#' @param dim dimension
#' @param mu vector of length dim for mean
#' @param Sigma dim x dim covariance matrix
#' @param start_beta beta for the base level
#' @param precondition either a logical value to determine if preconditioning matrices are
#'                     used (TRUE - and is set to be the variance of the sub-posterior samples)
#'                     or not (FALSE - and is set to be the identity matrix for all sub-posteriors),
#'                     or a list of length (1/start_beta) where precondition[[c]]
#'                     is the preconditioning matrix for sub-posterior c. Default is TRUE
#' @param resampling_method method to be used in resampling, default is multinomial 
#'                          resampling ('multi'). Other choices are stratified 
#'                          resampling ('strat'), systematic resampling ('system'),
#'                          residual resampling ('resid')
#' @param ESS_threshold number between 0 and 1 defining the proportion 
#'                      of the number of samples that ESS needs to be
#'                      lower than for resampling (i.e. resampling is carried 
#'                      out only when ESS < N*ESS_threshold)
#' @param diffusion_estimator choice of unbiased estimator for the Exact Algorithm
#'                            between "Poisson" (default) for Poisson estimator
#'                            and "NB" for Negative Binomial estimator
#' @param beta_NB beta parameter for Negative Binomial estimator (default 10)
#' @param gamma_NB_n_points number of points used in the trapezoidal estimation
#'                          of the integral found in the mean of the negative
#'                          binomial estimator (default is 2)
#' @param seed seed number - default is NULL, meaning there is no seed
#' @param n_cores number of cores to use
#'
#' @return A list with components:
#' \describe{
#'   \item{particles}{list of length (L-1), where particles[[l]][[i]] are the
#'                    particles for level l, node i}
#'   \item{proposed_samples}{list of length (L-1), where proposed_samples[[l]][[i]]
#'                           are the proposed samples for level l, node i}
#'   \item{time}{list of length (L-1), where time[[l]][[i]] is the run time for level l,
#'               node i}
#'   \item{ESS}{list of length (L-1), where ESS[[l]][[i]] is the effective 
#'              sample size of the particles after each step BEFORE deciding 
#'              whether or not to resample for level l, node i}
#'   \item{CESS}{list of length (L-1), where CESS[[l]][[i]] is the conditional
#'               effective sample size of the particles after each step}
#'   \item{resampled}{list of length (L-1), where resampled[[l]][[i]] is a 
#'                    boolean value to record if the particles were resampled
#'                    after each step; rho and Q for level l, node i}
#'   \item{precondition_matrices}{preconditioning matrices used in the algorithm 
#'                               for each node}
#'   \item{diffusion_times}{vector of length (L-1), where diffusion_times[l]
#'                         are the times for fusion in level l}
#' }
#'
#' @export
bal_binary_fusion_SMC_multiGaussian <- function(N_schedule,
                                                m_schedule,
                                                time_schedule,
                                                base_samples,
                                                L,
                                                dim,
                                                mu,
                                                Sigma,
                                                start_beta,
                                                precondition = TRUE,
                                                resampling_method = 'multi',
                                                ESS_threshold = 0.5,
                                                diffusion_estimator = 'Poisson',
                                                beta_NB = 10,
                                                gamma_NB_n_points = 2,
                                                seed = NULL,
                                                n_cores = parallel::detectCores()) {
  if (!is.vector(N_schedule) | (length(N_schedule)!=(L-1))) {
    stop("bal_binary_fusion_SMC_multiGaussian: N_schedule must be a vector of length (L-1)")
  } else if (!is.vector(m_schedule) | (length(m_schedule)!=(L-1))) {
    stop("bal_binary_fusion_SMC_multiGaussian: m_schedule must be a vector of length (L-1)")
  } else if (!is.vector(time_schedule) | (length(time_schedule)!=(L-1))) {
    stop("bal_binary_fusion_SMC_multiGaussian: time_schedule must be a vector of length (L-1)")
  } else if (!is.list(base_samples) | (length(base_samples)!=(1/start_beta))) {
    stop("bal_binary_fusion_SMC_multiGaussian: base_samples must be a list of length (1/start_beta)")
  } else if (!is.vector(mu) | (length(mu)!=dim)) {
    stop("bal_binary_fusion_SMC_multiGaussian: mu must be a vector of length dim")
  } else if (!is.matrix(Sigma)) {
    if (nrow(Sigma)!=dim | ncol(Sigma)!=dim) {
      stop("bal_binary_fusion_SMC_multiGaussian: Sigma must be a dim x dim matrix")
    }
  } else if (ESS_threshold < 0 | ESS_threshold > 1) {
    stop("bal_binary_fusion_SMC_multiGaussian: ESS_threshold must be between 0 and 1")
  }
  if (is.vector(m_schedule) & (length(m_schedule)==(L-1))) {
    for (l in (L-1):1) {
      if (((1/start_beta)/prod(m_schedule[(L-1):l]))%%1!=0) {
        stop("bal_binary_fusion_SMC_multiGaussian: check that (1/start_beta)/prod(m_schedule[(L-1):l])
              is an integer for l=L-1,...,1")
      }
    }
  } else {
    stop("bal_binary_fusion_SMC_multiGaussian: m_schedule must be a vector of length (L-1)")
  }
  m_schedule <- c(m_schedule, 1)
  particles <- list()
  if (all(sapply(base_samples, function(sub) class(sub)=='particle'))) {
    particles[[L]] <- base_samples
  } else if (all(sapply(base_samples, is.matrix))) {
    if (!all(sapply(base_samples, function(core) ncol(core)==dim))) {
      stop("bal_binary_fusion_SMC_multiGaussian: the sub-posterior samples in base_samples must be matrices with dim columns")
    }
    particles[[L]] <- initialise_particle_sets(samples_to_fuse = base_samples,
                                               multivariate = TRUE,
                                               number_of_steps = 2)
  } else {
    stop("bal_binary_fusion_SMC_multiGaussian: base_samples must be a list of length
         (1/start_beta) containing either items of class \"particle\" (representing
         particle approximations of the sub-posteriors) or are matrices with dim columns
         (representing un-normalised sample approximations of the sub-posteriors)")
  }
  proposed_samples <- list()
  time <- list()
  ESS <- list()
  CESS <- list()
  resampled <- list()
  precondition_matrices <- list()
  if (is.logical(precondition)) {
    if (precondition) {
      precondition_matrices[[L]] <- lapply(base_samples, cov)
    } else {
      precondition_matrices[[L]] <- rep(list(diag(1, dim)), (1/start_beta))
    }
  } else if (is.list(precondition)) {
    if (length(precondition)==(1/start_beta) & all(sapply(precondition, is.matrix))) {
      if (all(sapply(precondition, function(sub) ncol(sub)==dim))) {
        precondition_matrices[[L]] <- precondition  
      }
    }
  } else {
    stop("bal_binary_fusion_SMC_multiGaussian: precondition must be a logical indicating 
          whether or not a preconditioning matrix should be used, or a list of
          length C, where precondition[[c]] is the preconditioning matrix for
          the c-th sub-posterior")
  }
  cat('Starting bal_binary fusion \n', file = 'bal_binary_fusion_SMC_multiGaussian.txt')
  for (k in ((L-1):1)) {
    n_nodes <- max((1/start_beta)/prod(m_schedule[L:k]), 1)
    cat('########################\n', file = 'bal_binary_fusion_SMC_multiGaussian.txt',
        append = T)
    cat('Starting to fuse', m_schedule[k], 'densities of pi^beta, where beta =',
        prod(m_schedule[L:(k+1)]), '/', (1/start_beta), 'for level', k, 'with time',
        time_schedule[k], ', which is using', parallel::detectCores(), 'cores\n',
        file = 'bal_binary_fusion_SMC_multiGaussian.txt', append = T)
    cat('There are', n_nodes, 'nodes at this level each giving', N_schedule[k],
        'samples for beta =', prod(m_schedule[L:k]), '/', (1/start_beta),
        '\n', file = 'bal_binary_fusion_SMC_multiGaussian.txt', append = T)
    cat('########################\n', file = 'bal_binary_fusion_SMC_multiGaussian.txt', 
        append = T)
    fused <- lapply(X = 1:n_nodes, FUN = function(i) {
      previous_nodes <- ((m_schedule[k]*i)-(m_schedule[k]-1)):(m_schedule[k]*i)
      particles_to_fuse <- particles[[k+1]][previous_nodes]
      precondition_mats <- precondition_matrices[[k+1]][previous_nodes]
      parallel_fusion_SMC_multiGaussian(particles_to_fuse = particles_to_fuse,
                                        m = m_schedule[k],
                                        time = time_schedule[k],
                                        N = N_schedule[k],
                                        dim = dim,
                                        mu = mu,
                                        Sigma = Sigma,
                                        betas = rep(prod(m_schedule[L:(k+1)])*(start_beta), m_schedule[k]),
                                        precondition_matrices = precondition_mats,
                                        resampling_method = resampling_method,
                                        ESS_threshold = ESS_threshold,
                                        diffusion_estimator = diffusion_estimator,
                                        beta_NB = beta_NB,
                                        gamma_NB_n_points = gamma_NB_n_points,
                                        seed = seed,
                                        n_cores = n_cores)
    })
    # need to combine the correct samples
    particles[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$particles)
    proposed_samples[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$proposed_samples)
    time[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$time)
    ESS[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$ESS)
    CESS[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$CESS)
    resampled[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$resampled)
    precondition_matrices[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$precondition_matrices[[1]])
  }
  cat('Completed bal_binary fusion\n', file = 'bal_binary_fusion_SMC_multiGaussian.txt', append = T)
  if (length(particles[[1]])==1) {
    particles[[1]] <- particles[[1]][[1]]
    proposed_samples[[1]] <- proposed_samples[[1]][[1]]
    time[[1]] <- time[[1]][[1]]
    ESS[[1]] <- ESS[[1]][[1]]
    CESS[[1]] <- CESS[[1]][[1]]
    resampled[[1]] <- resampled[[1]][[1]]
    precondition_matrices[[1]] <- precondition_matrices[[1]][[1]]
  }
  return(list('particles' = particles,
              'proposed_samples' = proposed_samples,
              'time' = time,
              'ESS' = ESS,
              'CESS' = CESS,
              'resampled' = resampled,
              'precondition_matrices' = precondition_matrices,
              'diffusion_times' = time_schedule))
}

#' (Progressive) D&C Monte Carlo Fusion using SMC
#'
#' (Progressive) D&C Monte Carlo Fusion with multivariate Gaussian target
#'
#' @param N_schedule vector of length (L-1), where N_schedule[l] is the number 
#'                   of samples per node at level l
#' @param time_schedule vector of length(L-1), where time_schedule[l] is the time 
#'                      chosen for Fusion at level l
#' @param base_samples list of length (1/start_beta), where base_samples[[c]] 
#'                     contains the samples for the c-th node in the level
#' @param dim dimension
#' @param mu vector of length dim for mean
#' @param Sigma dim x dim covariance matrix
#' @param start_beta beta for the base level
#' @param precondition either a logical value to determine if preconditioning matrices are
#'                     used (TRUE - and is set to be the variance of the sub-posterior samples)
#'                     or not (FALSE - and is set to be the identity matrix for all sub-posteriors),
#'                     or a list of length (1/start_beta) where precondition[[c]]
#'                     is the preconditioning matrix for sub-posterior c. Default is TRUE
#' @param resampling_method method to be used in resampling, default is multinomial 
#'                          resampling ('multi'). Other choices are stratified 
#'                          resampling ('strat'), systematic resampling ('system'),
#'                          residual resampling ('resid')
#' @param ESS_threshold number between 0 and 1 defining the proportion 
#'                      of the number of samples that ESS needs to be
#'                      lower than for resampling (i.e. resampling is carried 
#'                      out only when ESS < N*ESS_threshold)
#' @param diffusion_estimator choice of unbiased estimator for the Exact Algorithm
#'                            between "Poisson" (default) for Poisson estimator
#'                            and "NB" for Negative Binomial estimator
#' @param beta_NB beta parameter for Negative Binomial estimator (default 10)
#' @param gamma_NB_n_points number of points used in the trapezoidal estimation
#'                          of the integral found in the mean of the negative
#'                          binomial estimator (default is 2)
#' @param seed seed number - default is NULL, meaning there is no seed
#' @param n_cores number of cores to use
#'
#' @return A list with components:
#' \describe{
#'   \item{particles}{list of length (L-1), where particles[[l]][[i]] are the
#'                   particles for level l, node i}
#'   \item{proposed_samples}{list of length (L-1), where proposed_samples[[l]][[i]]
#'                          are the proposed samples for level l, node i}
#'   \item{time}{list of length (L-1), where time[[l]][[i]] is the run time for level l,
#'              node i}
#'   \item{ESS}{list of length (L-1), where ESS[[l]][[i]] is the effective 
#'             sample size of the particles after each step BEFORE deciding 
#'             whether or not to resample for level l, node i}
#'   \item{CESS}{list of length (L-1), where CESS[[l]][[i]] is the conditional
#'              effective sample size of the particles after each step}
#'   \item{resampled}{list of length (L-1), where resampled[[l]][[i]] is a 
#'                   boolean value to record if the particles were resampled
#'                   after each step; rho and Q for level l, node i}
#'   \item{precondition_matrices}{preconditioning matrices used in the algorithm 
#'                               for each node}
#'   \item{diffusion_times}{vector of length (L-1), where diffusion_times[l]
#'                         are the times for fusion in level l}
#' }
#'
#' @export
progressive_fusion_SMC_multiGaussian <- function(N_schedule,
                                                 time_schedule,
                                                 base_samples,
                                                 dim,
                                                 mu,
                                                 Sigma,
                                                 start_beta,
                                                 precondition = TRUE,
                                                 resampling_method = 'multi',
                                                 ESS_threshold = 0.5,
                                                 diffusion_estimator = 'Poisson',
                                                 beta_NB = 10,
                                                 gamma_NB_n_points = 2,
                                                 seed = NULL,
                                                 n_cores = parallel::detectCores()) {
  if (!is.vector(N_schedule) | (length(N_schedule)!=(1/start_beta)-1)) {
    stop("progressive_fusion_SMC_multiGaussian: N_schedule must be a vector of length ((1/start_beta)-1)")
  } else if (!is.vector(time_schedule) | (length(time_schedule)!=(1/start_beta)-1)) {
    stop("progressive_fusion_SMC_multiGaussian: time_schedule must be a vector of length ((1/start_beta)-1)")
  } else if (!is.list(base_samples) | (length(base_samples)!=(1/start_beta))) {
    stop("progressive_fusion_SMC_multiGaussian: base_samples must be a list of length (1/start_beta)")
  } else if (!is.vector(mu) | (length(mu)!=dim)) {
    stop("progressive_fusion_SMC_multiGaussian: mu must be a vector of length dim")
  } else if (!is.matrix(Sigma)) {
    if (nrow(Sigma)!=dim | ncol(Sigma)!=dim) {
      stop("progressive_fusion_SMC_multiGaussian: Sigma must be a dim x dim matrix")
    }
  } else if (ESS_threshold < 0 | ESS_threshold > 1) {
    stop("progressive_fusion_SMC_multiGaussian: ESS_threshold must be between 0 and 1")
  }
  particles <- list()
  if (all(sapply(base_samples, function(sub) class(sub)=='particle'))) {
    particles[[(1/start_beta)]] <- base_samples
  } else if (all(sapply(base_samples, is.matrix))) {
    if (!all(sapply(base_samples, function(core) ncol(core)==dim))) {
      stop("progressive_fusion_SMC_multiGaussian: the sub-posterior samples in base_samples must be matrices with dim columns")
    }
    particles[[(1/start_beta)]] <- initialise_particle_sets(samples_to_fuse = base_samples,
                                                            multivariate = TRUE,
                                                            number_of_steps = 2)
  } else {
    stop("progressive_fusion_SMC_multiGaussian: base_samples must be a list of length
         (1/start_beta) containing either items of class \"particle\" (representing
         particle approximations of the sub-posteriors) zor are matrices (representing
         un-normalised sample approximations of the sub-posteriors)")
  }
  proposed_samples <- list()
  time <- list()
  ESS <- list()
  CESS <- list()
  resampled <- list()
  precondition_matrices <- list()
  if (is.logical(precondition)) {
    if (precondition) {
      precondition_matrices[[(1/start_beta)]] <- lapply(base_samples, cov)
    } else {
      precondition_matrices[[(1/start_beta)]] <- rep(list(diag(1, dim)), (1/start_beta))
    }
  } else if (is.list(precondition)) {
    if (length(precondition)==(1/start_beta) & all(sapply(precondition, is.matrix))) {
      if (all(sapply(precondition, function(sub) ncol(sub)==dim))) {
        precondition_matrices[[(1/start_beta)]] <- precondition  
      }
    }
  } else {
    stop("progressive_fusion_SMC_multiGaussian: precondition must be a logical indicating 
          whether or not a preconditioning matrix should be used, or a list of
          length C, where precondition[[c]] is the preconditioning matrix for
          the c-th sub-posterior")
  }
  index <- 2
  cat('Starting progressive fusion \n', file = 'progressive_fusion_SMC_multiGaussian.txt')
  for (k in ((1/start_beta)-1):1) {
    if (k==(1/start_beta)-1) {
      cat('########################\n', file = 'progressive_fusion_SMC_multiGaussian.txt', 
          append = T)
      cat('Starting to fuse', 2, 'densities for level', k, 'which is using', 
          parallel::detectCores(), 'cores\n', 
          file = 'progressive_fusion_SMC_multiGaussian.txt', append = T)
      cat('Fusing samples for beta =', 1, '/', (1/start_beta), 'with time',
          time_schedule[k], 'to get', N_schedule[k], 'samples for beta =', 2, 
          '/', (1/start_beta), '\n', file = 'progressive_fusion_SMC_multiGaussian.txt', 
          append = T)
      cat('########################\n', file = 'progressive_fusion_SMC_multiGaussian.txt',
          append = T)
      particles_to_fuse <- list(particles[[(1/start_beta)]][[1]], 
                                particles[[(1/start_beta)]][[2]])
      precondition_mats <- precondition_matrices[[k+1]][1:2]
      fused <- parallel_fusion_SMC_multiGaussian(particles_to_fuse = particles_to_fuse,
                                                 N = N_schedule[k],
                                                 m = 2,
                                                 time = time_schedule[k],
                                                 dim = dim,
                                                 mu = mu,
                                                 Sigma = Sigma,
                                                 betas = c(start_beta, start_beta),
                                                 precondition_matrices = precondition_mats,
                                                 resampling_method = resampling_method,
                                                 ESS_threshold = ESS_threshold,
                                                 diffusion_estimator = diffusion_estimator,
                                                 beta_NB = beta_NB,
                                                 gamma_NB_n_points = gamma_NB_n_points,
                                                 seed = seed,
                                                 n_cores = n_cores)
      
    } else {
      cat('########################\n', file = 'progressive_fusion_SMC_multiGaussian.txt', 
          append = T)
      cat('Starting to fuse', 2, 'densities for level', k, 'which is using', 
          parallel::detectCores(), 'cores\n',
          file = 'progressive_fusion_SMC_multiGaussian.txt', append = T)
      cat('Fusing samples for beta =', index, '/', (1/start_beta), 'and beta =', 
          1, '/', (1/start_beta), 'with time', time_schedule[k], 'to get',
          N_schedule[k], 'samples for beta =', (index+1), '/', (1/start_beta),
          '\n', file = 'progressive_fusion_SMC_multiGaussian.txt', append = T)
      cat('########################\n', file = 'progressive_fusion_SMC_multiGaussian.txt', 
          append = T)
      particles_to_fuse <- list(particles[[k+1]], 
                                particles[[(1/start_beta)]][[index+1]])
      precondition_mats <- list(precondition_matrices[[k+1]],
                                precondition_matrices[[(1/start_beta)]][[index+1]]) 
      fused <- parallel_fusion_SMC_multiGaussian(particles_to_fuse = particles_to_fuse,
                                                 N = N_schedule[k],
                                                 m = 2,
                                                 time = time_schedule[k],
                                                 dim = dim,
                                                 mu = mu,
                                                 Sigma = Sigma,
                                                 betas = c(index*start_beta, start_beta),
                                                 precondition_matrices = precondition_mats,
                                                 resampling_method = resampling_method,
                                                 ESS_threshold = ESS_threshold,
                                                 diffusion_estimator = diffusion_estimator,
                                                 beta_NB = beta_NB,
                                                 gamma_NB_n_points = gamma_NB_n_points,
                                                 seed = seed,
                                                 n_cores = n_cores)
      index <- index + 1
    }
    # need to combine the correct samples
    particles[[k]] <- fused$particles
    proposed_samples[[k]] <-fused$proposed_samples
    time[[k]] <- fused$time
    ESS[[k]] <- fused$ESS
    CESS[[k]] <- fused$CESS
    resampled[[k]] <- fused$resampled
    precondition_matrices[[k]] <- fused$precondition_matrices[[1]]
  }
  cat('Completed progressive fusion\n', file = 'progressive_fusion_SMC_multiGaussian.txt', append = T)
  return(list('particles' = particles,
              'proposed_samples' = proposed_samples,
              'time' = time,
              'ESS' = ESS,
              'CESS' = CESS,
              'resampled' = resampled,
              'precondition_matrices' = precondition_matrices,
              'diffusion_times' = time_schedule))
}
