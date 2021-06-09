#' Create a particle object
#'
#' Creates an R environment to contain the particles for SMC Fusion
#'
#' @param samples_to_fuse a list of samples (vector if univariate, matrix is multivariate)
#'                        that you wish to perform fusion with
#' @param multivariate logical value indicating if the samples are multivariate
#'                     (TRUE) or not (FALSE)
#'
#' @return A particle environment with components
#' \describe{
#'   \item{y_samples}{samples for y in particle set (initialised as the samples given)}
#'   \item{x_samples}{a list where x_samples[[i]] is the ith x sample for in 
#'                    particle set (all initialised as NA)}
#'   \item{x_mean}{the corresponding means for x_samples (initialised as NA)}
#'   \item{log_weights}{associated logarithm of the weights (initiliased as
#'                      the logarithm of 1/number of samples)}
#'   \item{normalised_weights}{associated normalised weights (initialised as 
#'                             1/number of samples)}
#'   \item{ESS}{effective sample size of particles (initialised as the 
#'              number of samples)}
#'   \item{CESS}{conditional effective sample size of particles after rho and Q step
#'               (initialised as NA for both)}
#'   \item{resampled}{logical value to idicate if particles have been resampled
#'                    after rho and Q step (initialised as FALSE for rho and TRUE for Q)}
#'   \item{N}{Number of particles}
#' }
#'
#' @examples 
#' p <- create_particle(samples = rnorm(10, 0, 1), multivariate = FALSE)
#' 
#' @export
create_particle <- function(samples, multivariate) {
  ps <- new.env(parent = emptyenv())
  if (multivariate) {
    N <- nrow(samples)
    dim <- ncol(samples)
    ps$y_samples <- samples
    ps$x_samples <- rep(list(NA), N)
    ps$x_means <- matrix(data = NA, nrow = N, ncol = dim)
    ps$log_weights <- log(rep(1/N, N))
    ps$normalised_weights <- rep(1/N, N)
    ps$ESS <- N
    ps$CESS <- c('rho' = NA, 'Q' = NA)
    ps$resampled <- c('rho' = FALSE, 'Q' = TRUE)
    ps$N <- N
    class(ps) <- "particle"
    return(ps)
  } else {
    N <- length(samples)
    ps$y_samples <- samples
    ps$x_samples <- rep(list(NA), N)
    ps$x_means <- rep(NA, N)
    ps$log_weights <- log(rep(1/N, N))
    ps$normalised_weights <- rep(1/N, N)
    ps$ESS <- N
    ps$CESS <- c('rho' = NA, 'Q' = NA)
    ps$resampled <- c('rho' = FALSE, 'Q' = TRUE)
    ps$N <- N
    class(ps) <- "particle"
    return(ps)
  }
}

#' Initialise particle sets from a list of samples
#'
#' Function to intialise particle sets from a list of samples
#'
#' @param samples_to_fuse a list of samples that you wish to perform fusion with
#' @param multivariate logical value indicating if the samples are multivariate
#'                     (TRUE) or not (FALSE)
#'
#' @return A list of particles to fuse, where the cth component is the particle 
#'         for sub-posterior c. In particular, each item in the list is an enviroment
#'         with components
#' \describe{
#'   \item{y_samples}{samples for y in particle set (initialised as the samples given)}
#'   \item{x_samples}{a list where x_samples[[i]] is the ith x sample for in 
#'                    particle set (all initialised as NA)}
#'   \item{x_mean}{the corresponding means for x_samples (initialised as NA)}
#'   \item{log_weights}{associated logarithm of the weights (initiliased as
#'                      the logarithm of 1/number of samples)}
#'   \item{normalised_weights}{associated normalised weights (initialised as 
#'                             1/number of samples)}
#'   \item{ESS}{effective sample size of particles (initialised as the 
#'              number of samples)}
#'   \item{CESS}{conditional effective sample size of particles after rho and Q step
#'               (initialised as NA for both)}
#'   \item{resampled}{logical value to idicate if particles have been resampled
#'                    after rho and Q step (initialised as FALSE for rho and TRUE for Q)}
#'   \item{N}{Number of particles}
#' }
#'
#' @examples 
#' # univariate
#' uni_samples <- lapply(1:2, function(i) rnorm(100, 0, 1))
#' particles <- initialise_particle_sets(samples_to_fuse = uni_samples, multivariate = FALSE)
#' 
#' # multivariate
#' multi_samples <- lapply(1:2, function(i) mvrnormArma(100, c(0, 0), diag(2)))
#' particles <- initialise_particle_sets(samples_to_fuse = multi_samples, multivariate = TRUE)
#' 
#' @export
initialise_particle_sets <- function(samples_to_fuse, multivariate) {
  if (!is.list(samples_to_fuse)) {
    stop("initialise_particle_sets: samples_to_fuse must be a list")
  }
  if (multivariate) {
    if (any(!sapply(samples_to_fuse, is.matrix))) {
      stop("initialise_particle_sets: if multivariate is TRUE, each of the samples in samples_to_fuse must be a matrix")
    }
    if (any(!sapply(samples_to_fuse[-1], function(samples) ncol(samples)==ncol(samples_to_fuse[[1]])))) {
      stop("initialise_particle_sets: if multivariate is TRUE, each of the samples in samples_to_fuse must be a matrix with the same number of columns")
    }
  } else {
    if (any(!sapply(samples_to_fuse, is.vector))) {
      stop("initialise_particle_sets: if multivariate is FALSE, each of the samples in samples_to_fuse must be a vector")
    }
  }
  return(lapply(samples_to_fuse, function(samples) create_particle(samples = samples, multivariate = multivariate)))
}

# ---------- functions for resampling schemes from Murray Pollock Github

multi.resamp <- function(normalised_weights, n = length(normalised_weights)) { 
  # Multinomial Resampling
  # Check whether resampling is possible
  if (sum(normalised_weights) > 0 & (n >0)) {
    if (sum(normalised_weights) != 1) {
      normalised_weights <- normalised_weights / sum(normalised_weights)
    }
    return(sample.int(n = length(normalised_weights),
                      size = n,
                      replace = TRUE,
                      prob = normalised_weights))
  } else {
    return(numeric(0))
  }
}

system.resamp <- function(normalised_weights, n = length(normalised_weights)) { 
  # Systematic Resampling
  # Check whether resampling is possible
  if (sum(normalised_weights) > 0 & (n >0)) {
    if (sum(normalised_weights) != 1) {
      normalised_weights <- normalised_weights / sum(normalised_weights)
    }
    cum.wei <- cumsum(normalised_weights)
    samp.vec <- seq(runif(1,0,1/n),1,1/n) 
    return(sapply(1:n, function(i) length(cum.wei[samp.vec[i]>=cum.wei])+1))
  } else {
    return(numeric(0))
  }  
}

strat.resamp <- function(normalised_weights, n = length(normalised_weights)) { 
  # Stratified Resampling
  # Check whether resampling is possible
  if (sum(normalised_weights) > 0 & (n >0)) {
    if (sum(normalised_weights) != 1) {
      normalised_weights <- normalised_weights / sum(normalised_weights)
    }
    cum.wei <- cumsum(normalised_weights)
    samp.vec <- seq(0, (n-1)/n, 1/n) + runif(n, 0, 1/n) 
    return(findInterval(samp.vec,cum.wei)+1) 
  } else {
    return(numeric(0))
  }
}

resid.resamp <- function(normalised_weights, 
                         n = length(normalised_weights), 
                         nest.resamp = strat.resamp) { 
  # Residual Resampling
  # Check whether resampling is possible
  if (sum(normalised_weights) > 0 & (n >0)) {
    if (sum(normalised_weights) != 1) {
      normalised_weights <- normalised_weights / sum(normalised_weights)
    }
    det.resamp <- floor(n*normalised_weights)
    return(c(rep(1:length(normalised_weights), det.resamp),
             nest.resamp((normalised_weights-det.resamp/n),
                         n = n-sum(det.resamp))))
  } else {
    return(numeric(0))
  }
}

#' Resampling methods
#'
#' @param normalised_weights a vector of normalised weights
#' @param method method to be used in resampling, default is multinomial 
#'               resampling ('multi'). Other choices are stratified ('strat'), 
#'               systematic ('system'), residual ('resid')
#' @param n number of samples to resample
#'
#' @return a vector of resampled indicies
#' 
#' @export
resample_indices <- function(normalised_weights,
                             method = 'multi',
                             n = length(normalised_weights)) {
  # double check that weights are normalised
  normalised_weights <- normalised_weights / sum(normalised_weights)
  # check that resampling is possible
  if (sum(normalised_weights) > 0 & (n >0)) {
    if (method == 'multi') {
      return(multi.resamp(normalised_weights = normalised_weights, n = n))
    } else if (method == 'system') {
      return(system.resamp(normalised_weights = normalised_weights, n = n))
    } else if (method == 'strat') {
      return(strat.resamp(normalised_weights = normalised_weights, n = n))
    } else if (method == 'resid') {
      return(resid.resamp(normalised_weights = normalised_weights, n = n))
    } else {
      stop('resample: method must be one of \'mutli\', \'system\', \'strat\' or \'resid\'')
    }
  } else {
    return(numeric(0))
  }
}

#' Resampling a particle set for y
#' 
#' Resamples the particle set for the samples for y
#'
#' @param N number of samples (default to particle_set$N)
#' @param multivariate logical value indicating if the particles are multivariate
#'                     (TRUE) or not (FALSE)
#' @param resampling_method method to be used in resampling, default is multinomial 
#'                          resampling ('multi'). Other choices are stratified 
#'                          resampling ('strat'), systematic resampling ('system'),
#'                          residual resampling ('resid')
#' @param seed seed number - default is NULL, meaning there is no seed
#'
#' @return resampled particle set
#' 
#' @examples
#' particles <- create_particle(samples = rnorm(10, 0, 1), multivariate = FALSE)
#' particles <- resample_particle_y_samples(particle_set = particles,
#'                                  multivariate = FALSE,
#'                                  resampling_method = 'resid',
#'                                  seed = 21)
#' 
#' @export
resample_particle_y_samples <- function(N = particle_set$N,
                                        particle_set, 
                                        multivariate, 
                                        resampling_method = 'multi',
                                        seed = NULL) {
  if (!("particle" %in% class(particle_set))) {
    stop("resample_particle_y_samples: particle_set must be a \"particle\" object")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (particle_set$resampled['Q'] & particle_set$N==N) {
    return(particle_set)
  } else {
    indices <- resample_indices(normalised_weights = particle_set$normalised_weights,
                                method = resampling_method,
                                n = N)
    if (multivariate) {
      particle_set$y_samples <- particle_set$y_samples[indices,]
    } else {
      particle_set$y_samples <- particle_set$y_samples[indices]
    }
    particle_set$log_weights <- rep(log(1/N), N)
    particle_set$normalised_weights <- rep(1/N, N)
    particle_set$ESS <- N
    particle_set$resampled['Q'] <- TRUE
    particle_set$N <- N
    return(particle_set)
  }
}

#' Resampling a particle set for x
#' 
#' Resamples the particle set for the samples for x
#'
#' @param N number of samples (default to particle_set$N)
#' @param particle_set particle object (see initialise_particle_set)
#' @param multivariate logical value indicating if the particles are multivariate
#'                     (TRUE) or not (FALSE)
#' @param resampling_method method to be used in resampling, default is multinomial 
#'                          resampling ('multi'). Other choices are stratified 
#'                          resampling ('strat'), systematic resampling ('system'),
#'                          residual resampling ('resid')
#' @param seed seed number - default is NULL, meaning there is no seed
#'
#' @return resampled particle set
#' 
#' @examples
#' samples_to_fuse <- lapply(1:2, function(i) rnorm(100, 0, 1))
#' particles_to_fuse <- initialise_particle_sets(samples_to_fuse = samples_to_fuse,
#'                                               multivariate = FALSE)
#' precondition_values <- sapply(samples_to_fuse, var)
#' particles <- rho_IS_univariate(particles_to_fuse = particles_to_fuse,
#'                                N = 100,
#'                                m = 2,
#'                                time = 0.5,
#'                                precondition_values = precondition_values)
#' particles <- resample_particle_x_samples(particle_set = particles,
#'                                          multivariate = FALSE,
#'                                          resampling_method = 'resid',
#'                                          seed = 21)
#' 
#' @export
resample_particle_x_samples <- function(N = particle_set$N,
                                        particle_set,
                                        multivariate,
                                        resampling_method = 'multi',
                                        seed = NULL) {
  if (!("particle" %in% class(particle_set))) {
    stop("resample_particle_x_samples: particle_set must be a \"particle\" object")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (particle_set$resampled['rho'] & particle_set$N==N) {
    return(particle_set)
  } else {
    indices <- resample_indices(normalised_weights = particle_set$normalised_weights,
                                method = resampling_method,
                                n = N)
    if (multivariate) {
      particle_set$x_means <- particle_set$x_means[indices,]
    } else {
      particle_set$x_means <- particle_set$x_means[indices]
    }
    particle_set$x_samples <- particle_set$x_samples[indices]
    particle_set$log_weights <- rep(log(1/N), N)
    particle_set$normalised_weights <- rep(1/N, N)
    particle_set$ESS <- N
    particle_set$resampled['Q'] <- TRUE
    particle_set$N <- N
    return(particle_set)
  }
}

#' rho Importance Sampling Step (univariate)
#' 
#' Performs the importance sampling step for rho where target is univariate
#'
#' @param particles_to_fuse list of length m, where particles_to_fuse[[c]]
#'                          contains the particles for the c-th sub-posterior
#'                          (a list of particles to fuse can be initialised by
#'                          initialise_particle_sets() function)
#' @param N number of particles to importance sample
#' @param m number of sub-posteriors to combine
#' @param time time T for fusion algorithm
#' @param precondition_values vector of length m, where precondition_values[[c]]
#'                            is the precondition value for sub-posterior c
#' @param n_cores number of cores to use
#' @param cl an object of class "cluster" for parallel computation in R. If none
#'           is passed, then one is created and used within this function
#' 
#' @return A importance weighted particle set
#' 
#' @examples
#' samples_to_fuse <- lapply(1:2, function(i) rnorm(100, 0, 1))
#' particles_to_fuse <- initialise_particle_sets(samples_to_fuse = samples_to_fuse,
#'                                               multivariate = FALSE)
#' precondition_values <- sapply(samples_to_fuse, var)
#' particles <- rho_IS_univariate(particles_to_fuse = particles_to_fuse,
#'                                N = 100,
#'                                m = 2,
#'                                time = 0.5,
#'                                precondition_values = precondition_values)
#'
#' @export
rho_IS_univariate <- function(particles_to_fuse,
                              N,
                              m,
                              time,
                              precondition_values,
                              n_cores = parallel::detectCores(),
                              cl = NULL) {
  if (!is.list(particles_to_fuse) | (length(particles_to_fuse)!=m)) {
    stop("rho_IS_univariate: particles_to_fuse must be a list of length m")
  } else if (!all(sapply(particles_to_fuse, function(sub_posterior) ("particle" %in% class(sub_posterior))))) {
    stop("rho_IS_univariate: particles in particles_to_fuse must be \"particle\" objects")
  } else if (!all(sapply(particles_to_fuse, function(sub_posterior) sub_posterior$N==N))) {
    stop("rho_IS_univariate: particles in particles_to_fuse must have the same number of samples")
  } else if (!is.vector(precondition_values) | length(precondition_values)!=m) {
    stop("rho_IS_univariate: precondition_values must be a vector of length m")
  } else if (any(class(cl)=="cluster") | is.null(cl)) {
    stop("rho_IS_univariate: cl must be a \"cluster\" object or NULL")
  }
  # ---------- creating parallel cluster
  if (is.null(cl)) {
    print('Creating cluster in rho_IS_univariate')
    cl <- parallel::makeCluster(n_cores, setup_strategy = "sequential")
    parallel::clusterExport(cl, envir = environment(), varlist = c("rho_IS_univariate_"))
    close_cluster <- TRUE
  } else {
    close_cluster <- FALSE
  }
  if (!is.null(seed)) {
    parallel::clusterSetRNGStream(cl, iseed = seed)
  }
  # split the y_samples that we want to combine
  max_samples_per_core <- ceiling(N/n_cores)
  split_indices <- split(1:N, ceiling(seq_along(1:N)/max_samples_per_core))
  split_samples_to_fuse <- lapply(split_indices, function(indices) {
    lapply(1:m, function(i) particles_to_fuse[[i]]$y_samples[indices])
  })
  rho_IS <- parallel::parLapply(cl, X = 1:length(split_indices), fun = function(core) {
    rho_IS_univariate_(samples_to_fuse = split_samples_to_fuse[[core]],
                       N = length(split_indices[[core]]),
                       m = m,
                       time = time,
                       precondition_values = precondition_values)
  })
  if (close_cluster) {
    print('Closing cluster in rho_IS_univariate')
    parallel::stopCluster(cl)
  }
  # ---------- create particle
  ps <- new.env(parent = emptyenv())
  ps$y_samples <- matrix(data = NA, nrow = N, ncol = dim)
  ps$x_samples <- unlist(lapply(1:length(split_indices), function(i) rho_IS[[i]]$x_samples), recursive = FALSE)
  ps$x_means <-  unlist(lapply(1:length(split_indices), function(i) rho_IS[[i]]$x_means))
  lw <- unlist(lapply(1:length(split_indices), function(i) rho_IS[[i]]$log_weights))
  norm_weights <- particle_ESS(log_weights = lw)
  ps$log_weights <- norm_weights$log_weights
  ps$lw$rho <- lw
  ps$normalised_lw$rho <- norm_weights$log_weights
  ps$normalised_weights <- norm_weights$normalised_weights
  ps$ESS <- norm_weights$ESS
  ps$CESS <- c('rho' = norm_weights$ESS, 'Q' = NA)
  ps$resampled <- c('rho' = FALSE, 'Q' = FALSE)
  ps$N <- N
  class(ps) <- "particle"
  return(ps)
}

#' rho Importance Sampling Step (multivariate)
#' 
#' Performs the importance sampling step for rho where target is univariate
#'
#' @param particles_to_fuse list of length m, where particles_to_fuse[[c]]
#'                          contains the particles for the c-th sub-posterior
#'                          (a list of particles to fuse can be initialised by
#'                          initialise_particle_sets() function)
#' @param N number of particles to importance sample
#' @param dim dimension of the particles
#' @param m number of sub-posteriors to combine
#' @param time time T for fusion algorithm
#' @param inv_precondition_matrices list of length m of inverse 
#'                                  preconditioning matrices
#' @param inverse_sum_inv_precondition_matrices the inverse of the sum of the inverse
#'                                              precondition matrices (can be 
#'                                              calculated by passing the inverse 
#'                                              precondition matrices into inverse_sum_matrices())
#' @param n_cores number of cores to use
#' @param cl an object of class "cluster" for parallel computation in R. If none
#'           is passed, then one is created and used within this function
#'
#' @return A importance weighted particle set
#' 
#' @examples
#' samples_to_fuse <- lapply(1:2, function(i) mvrnormArma(100, c(0, 0), diag(2)))
#' particles_to_fuse <- initialise_particle_sets(samples_to_fuse = samples_to_fuse,
#'                                               multivariate = TRUE)
#' precondition_mats <- lapply(samples_to_fuse, cov)
#' inv_precondition_mats <- lapply(precondition_mats, solve)
#' inv_sum_inv_precondition_mats <- inverse_sum_matrices(inv_precondition_mats)
#' particles <- rho_IS_multivariate(particles_to_fuse = particles_to_fuse,
#'                                  N = 100,
#'                                  dim = 2,
#'                                  m = 2,
#'                                  time = 0.5,
#'                                  inv_precondition_matrices = inv_precondition_mats,
#'                                  inverse_sum_inv_precondition_matrices = inv_sum_inv_precondition_mats)
#'
#' @export
rho_IS_multivariate <- function(particles_to_fuse,
                                dim,
                                N,
                                m,
                                time,
                                inv_precondition_matrices,
                                inverse_sum_inv_precondition_matrices,
                                n_cores = parallel::detectCores(),
                                cl = NULL) {
  if (!is.list(particles_to_fuse) | (length(particles_to_fuse)!=m)) {
    stop("rho_IS_multivariate: particles_to_fuse must be a list of length m")
  } else if (!all(sapply(particles_to_fuse, function(sub_posterior) ("particle" %in% class(sub_posterior))))) {
    stop("rho_IS_multivariate: particles in particles_to_fuse must be \"particle\" objects")
  } else if (!all(sapply(particles_to_fuse, function(sub_posterior) sub_posterior$N==N))) {
    stop("rho_IS_multivariate: particles in particles_to_fuse must have the same number of samples")
  } else if (!is.list(inv_precondition_matrices) | length(inv_precondition_matrices)!=m) {
    stop("rho_IS_multivariate: inv_precondition_matrices must be a list of length m")
  } else if (!any(class(cl)=="cluster") & !is.null(cl)) {
    stop("rho_IS_multivariate: cl must be a \"cluster\" object or NULL")
  }
  # ---------- creating parallel cluster
  if (is.null(cl)) {
    print('Creating cluster in rho_IS_multivariate')
    cl <- parallel::makeCluster(n_cores, setup_strategy = "sequential")
    parallel::clusterExport(cl, envir = environment(), varlist = c("rho_IS_multivariate_"))
    close_cluster <- TRUE
  } else {
    close_cluster <- FALSE
  }
  if (!is.null(seed)) {
    parallel::clusterSetRNGStream(cl, iseed = seed)
  }
  # split the y_samples that we want to combine
  max_samples_per_core <- ceiling(N/n_cores)
  split_indices <- split(1:N, ceiling(seq_along(1:N)/max_samples_per_core))
  split_samples_to_fuse <- lapply(split_indices, function(indices) {
    lapply(1:m, function(i) particles_to_fuse[[i]]$y_samples[indices,,drop = FALSE])
  })
  rho_IS <- parallel::parLapply(cl, X = 1:length(split_indices), fun = function(core) {
    rho_IS_multivariate_(samples_to_fuse = split_samples_to_fuse[[core]],
                         dim = dim,
                         N = length(split_indices[[core]]),
                         m = m,
                         time = time,
                         inv_precondition_matrices = inv_precondition_matrices,
                         inverse_sum_inv_precondition_matrices = inverse_sum_inv_precondition_matrices)
  })
  if (close_cluster) {
    print('Closing cluster in rho_IS_multivariate')
    parallel::stopCluster(cl)
  }
  # ---------- create particle
  ps <- new.env(parent = emptyenv())
  ps$y_samples <- matrix(data = NA, nrow = N, ncol = dim)
  ps$x_samples <- unlist(lapply(1:length(split_indices), function(i) rho_IS[[i]]$x_samples), recursive = FALSE)
  ps$x_means <-  do.call(rbind, lapply(1:length(split_indices), function(i) rho_IS[[i]]$x_means))
  lw <- unlist(lapply(1:length(split_indices), function(i) rho_IS[[i]]$log_weights))
  norm_weights <- particle_ESS(log_weights = lw)
  ps$log_weights <- norm_weights$log_weights
  ps$lw$rho <- lw
  ps$normalised_lw$rho <- norm_weights$log_weights
  ps$normalised_weights <- norm_weights$normalised_weights
  ps$ESS <- norm_weights$ESS
  ps$CESS <- c('rho' = norm_weights$ESS, 'Q' = NA)
  ps$resampled <- c('rho' = FALSE, 'Q' = FALSE)
  ps$N <- N
  class(ps) <- "particle"
  return(ps)
}
