#' rho_j Importance Sampling Step
#'
#' rho_j Importance Sampling weighting for bivariate Gaussian distributions
#'
#' @param particle_set particles set prior to Q importance sampling step
#' @param m number of sub-posteriors to combine
#' @param time_mesh time mesh used in Bayesian Fusion
#' @param dim dimension
#' @param mean_vecs list of length m, where mean_vecs[[c]] is a vector of 
#'                  length dim for the mean of sub-posterior c
#' @param inv_Sigmas list of length m, where inv_Sigmas[[c]] is a dim x dim
#'                   inverse covariance matrix for sub-posterior c
#' @param precondition_matrices list of length m, where precondition_matrices[[c]]
#'                               is the precondition matrix for sub-posterior c
#' @param inv_precondition_matrices list of length m, where inv_precondition_matrices[[c]]
#'                                  is the inverse precondition matrix for sub-posterior c
#' @param Lambda inverse of the sum of the inverse precondition matrices (which
#'               can be computed using inverse_sum_matrices(inv_precondition_matrices))
#' @param resampling_method method to be used in resampling, default is multinomial
#'                          resampling ('multi'). Other choices are stratified
#'                          resampling ('strat'), systematic resampling ('system'),
#'                          residual resampling ('resid')
#' @param ESS_threshold number between 0 and 1 defining the proportion of the
#'                      number of samples that ESS needs to be lower than for
#'                      resampling (i.e. resampling is carried out only when
#'                      ESS < N*ESS_threshold)
#' @param sub_posterior_means matrix with m rows and dim columns, where sub_posterior_means[c,]
#'                            is the sub-posterior mean of sub-posterior c
#' @param adaptive_mesh logical value to indicate if an adaptive mesh is used
#'                      (default is FALSE)
#' @param adaptive_mesh_parameters list of parameters used for adaptive mesh
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
#'   \item{particle_set}{updated particle set after the iterative rho_j steps}
#'   \item{proposed_samples}{proposal samples for the last time step}
#'   \item{time}{elapsed time of each step of the algorithm}
#'   \item{ESS}{effective sample size of the particles after each step}
#'   \item{CESS}{conditional effective sample size of the particles after each step}
#'   \item{resampled}{boolean value to indicate if particles were resampled
#'                    after each time step}
#'   \item{E_nu_j}{approximation of the average variation of the trajectories
#'                 at each time step}
#'   \item{chosen}{which term was chosen if using an adaptive mesh at each time step}
#'   \item{mesh_terms}{the evaluated terms in deciding the mesh at each time step}
#' }
#' 
#' @export
rho_j_multiGaussian <- function(particle_set,
                                m,
                                time_mesh,
                                dim,
                                mean_vecs,
                                inv_Sigmas,
                                precondition_matrices,
                                inv_precondition_matrices,
                                Lambda,
                                resampling_method = 'multi',
                                ESS_threshold = 0.5,
                                sub_posterior_means = NULL,
                                adaptive_mesh = FALSE,
                                adaptive_mesh_parameters = NULL,
                                diffusion_estimator = 'Poisson',
                                beta_NB = 10,
                                gamma_NB_n_points = 2,
                                seed = NULL,
                                n_cores = parallel::detectCores()) {
  if (!("particle" %in% class(particle_set))) {
    stop("rho_j_multiGaussian: particle_set must be a \"particle\" object")
  } else if (!is.vector(time_mesh)) {
    stop("rho_j_multiGaussian: time_mesh must be an ordered vector of length >= 2")
  } else if (length(time_mesh) < 2) {
    stop("rho_j_multiGaussian: time_mesh must be an ordered vector of length >= 2")
  } else if (!identical(time_mesh, sort(time_mesh))) {
    stop("rho_j_multiGaussian: time_mesh must be an ordered vector of length >= 2")
  } else if (any(sapply(1:m, function(c) (!is.vector(mean_vecs[[c]]) | (length(mean_vecs[[c]])!=dim))))) {
    stop("rho_j_multiGaussian: mean_vecs[[c]] must be a vector of length dim for each c")
  } else if (any(sapply(1:m, function(c) !is.matrix(inv_Sigmas[[c]])))) {
    stop("rho_j_multiGaussian: inv_Sigmas[[c]] must be a dim x dim matrix for each c")
  } else if (!is.list(precondition_matrices) | (length(precondition_matrices)!=m)) {
    stop("rho_j_multiGaussian: precondition_matrices must be a list of length m")
  } else if (!is.list(inv_precondition_matrices) | (length(inv_precondition_matrices)!=m)) {
    stop("rho_j_multiGaussian: inv_precondition_matrices must be a list of length m")
  } else if ((ESS_threshold < 0) | (ESS_threshold > 1)) {
    stop("rho_j_multiGaussian: ESS_threshold must be between 0 and 1")
  }
  if (adaptive_mesh) {
    if (!is.matrix(sub_posterior_means)) {
      stop("rho_j_multiGaussian: if adaptive_mesh==TRUE, sub_posterior_means must be a (m x dim) matrix")
    } else if (any(dim(sub_posterior_means)!=c(m,dim))) {
      stop("rho_j_multiGaussian: if adaptive_mesh==TRUE, sub_posterior_means must be a (m x dim) matrix")
    }
  }
  transform_matrices <- lapply(1:m, function(c) {
    list('to_Z' = expm::sqrtm(inv_precondition_matrices[[c]]),
         'to_X' = expm::sqrtm(precondition_matrices[[c]]))
  })
  inv_Sigma_Z <- lapply(1:m, function(c) transform_matrices[[c]]$to_X %*% inv_Sigmas[[c]] %*% transform_matrices[[c]]$to_X)
  N <- particle_set$N
  # ---------- creating parallel cluster
  cl <- parallel::makeCluster(n_cores)
  parallel::clusterExport(cl, envir = environment(), varlist = ls())
  parallel::clusterExport(cl, varlist = ls("package:DCFusion"))
  parallel::clusterExport(cl, varlist = ls("package:layeredBB"))
  if (!is.null(seed)) {
    parallel::clusterSetRNGStream(cl, iseed = seed)
  }
  max_samples_per_core <- ceiling(N/n_cores)
  split_indices <- split(1:N, ceiling(seq_along(1:N)/max_samples_per_core))
  elapsed_time <- rep(NA, length(time_mesh)-1)
  ESS <- c(particle_set$ESS[1], rep(NA, length(time_mesh)-1))
  CESS <- c(particle_set$CESS[1], rep(NA, length(time_mesh)-1))
  resampled <- rep(FALSE, length(time_mesh))
  if (adaptive_mesh) {
    E_nu_j <- rep(NA, length(time_mesh))
    chosen <- rep("", length(time_mesh))
    mesh_terms <- rep(list(c(NA,NA)), length(time_mesh))
    k4_choice <- rep(NA, length(time_mesh))
  } else {
    E_nu_j <- NA
    chosen <- NULL
    mesh_terms <- NULL
    k4_choice <- NULL
  }
  # iterative proposals
  end_time <- time_mesh[length(time_mesh)]
  j <- 1
  while (time_mesh[j]!=end_time) {
    pcm <- proc.time()
    j <- j+1
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
    # ----------- if adaptive_mesh==TRUE, find mesh for jth iteration
    if (adaptive_mesh) {
      if (particle_set$number_of_steps < j) {
        particle_set$number_of_steps <- j
        particle_set$CESS[j] <- NA
        particle_set$resampled[j] <- FALSE
      }
      tilde_Delta_j <- mesh_guidance_adaptive(C = m,
                                              d = dim,
                                              data_size = adaptive_mesh_parameters$data_size,
                                              b = adaptive_mesh_parameters$b,
                                              threshold = adaptive_mesh_parameters$CESS_j_threshold,
                                              particle_set = particle_set,
                                              sub_posterior_means = sub_posterior_means,
                                              inv_precondition_matrices = inv_precondition_matrices,
                                              k3 = adaptive_mesh_parameters$k3,
                                              k4 = adaptive_mesh_parameters$k4,
                                              vanilla = adaptive_mesh_parameters$vanilla)
      E_nu_j[j] <- tilde_Delta_j$E_nu_j
      chosen[j] <- tilde_Delta_j$chosen
      mesh_terms[[j]] <- c(tilde_Delta_j$T1, tilde_Delta_j$T2)
      k4_choice[j] <- tilde_Delta_j$k4_choice
      time_mesh[j] <- min(end_time, time_mesh[j-1]+tilde_Delta_j$max_delta_j)
    }
    # split the x samples from the previous time marginal (and their means) into approximately equal lists
    split_x_samples <- lapply(split_indices, function(indices) particle_set$x_samples[indices])
    split_x_means <- lapply(split_indices, function(indices) particle_set$x_means[indices,,drop = FALSE])
    V <- construct_V(s = time_mesh[j-1],
                     t = time_mesh[j],
                     end_time = end_time,
                     C = m,
                     d = dim,
                     precondition_matrices = precondition_matrices,
                     Lambda = Lambda,
                     iteration = j)
    rho_j_weighted_samples <- parallel::parLapply(cl, X = 1:length(split_indices), fun = function(core) {
      split_N <- length(split_indices[[core]])
      x_mean_j <- matrix(data = NA, nrow = split_N, ncol = dim)
      log_rho_j <- rep(0, split_N)
      x_j <- lapply(1:split_N, function(i) {
        M <- construct_M(s = time_mesh[j-1],
                         t = time_mesh[j],
                         end_time = end_time,
                         C = m,
                         d = dim,
                         sub_posterior_samples = split_x_samples[[core]][[i]],
                         sub_posterior_mean = split_x_means[[core]][i,],
                         iteration = j)
        if (time_mesh[j]!=end_time) {
          return(matrix(mvrnormArma(N = 1, mu = M, Sigma = V), nrow = m, ncol = dim, byrow = TRUE))
        } else {
          return(matrix(mvtnorm::rmvnorm(n = 1, mean = M, sigma = V), nrow = m, ncol = dim, byrow = TRUE))
        }
      })
      for (i in 1:split_N) {
        x_mean_j[i,] <- weighted_mean_multivariate(matrix = x_j[[i]],
                                                   weights = inv_precondition_matrices,
                                                   inverse_sum_weights = Lambda)
        log_rho_j[i] <- sum(sapply(1:m, function(c) {
          ea_multiGaussian_DL_PT(x0 = as.vector(split_x_samples[[core]][[i]][c,]),
                                 y = as.vector(x_j[[i]][c,]),
                                 s = time_mesh[j-1],
                                 t = time_mesh[j],
                                 dim = dim,
                                 mu = mean_vecs[[c]],
                                 inv_Sigma = inv_Sigmas[[c]],
                                 inv_Sigma_Z = inv_Sigma_Z[[c]],
                                 beta = 1,
                                 precondition_mat = precondition_matrices[[c]],
                                 transform_mats = transform_matrices[[c]],
                                 diffusion_estimator = diffusion_estimator,
                                 beta_NB = beta_NB,
                                 gamma_NB_n_points = gamma_NB_n_points,
                                 logarithm = TRUE)}))
      }
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
    elapsed_time[j-1] <- (proc.time()-pcm)['elapsed']
  }
  parallel::stopCluster(cl)
  # set the y samples as the first element of each of the x_samples
  if (dim == 1) {
    proposed_samples <- as.matrix(sapply(1:N, function(i) particle_set$x_samples[[i]][1,]))
  } else {
    proposed_samples <- t(sapply(1:N, function(i) particle_set$x_samples[[i]][1,]))
  }
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
  if (adaptive_mesh) {
    CESS <- CESS[1:j]
    ESS <- ESS[1:j]
    resampled <- resampled[1:j]
    particle_set$time_mesh <- time_mesh[1:j]
    elapsed_time <- elapsed_time[1:(j-1)]
    E_nu_j <- E_nu_j[1:j]
    chosen <- chosen[1:j]
    mesh_terms <- mesh_terms[1:j]
    k4_choice <- k4_choice[1:j]
  }
  return(list('particle_set' = particle_set,
              'proposed_samples' = proposed_samples,
              'time' = elapsed_time,
              'ESS' = ESS,
              'CESS' = CESS,
              'resampled' = resampled,
              'E_nu_j' = E_nu_j,
              'chosen' = chosen,
              'mesh_terms' = mesh_terms,
              'k4_choice' = k4_choice))
}

#' Generalised Bayesian Fusion [parallel]
#'
#' Generalised Bayesian Fusion with bivariate Gaussian target
#'
#' @param particles_to_fuse list of length m, where particles_to_fuse[[c]]
#'                          contains the particles for the c-th sub-posterior
#'                          (a list of particles to fuse can be initialised by 
#'                          initialise_particle_sets() function)
#' @param N number of samples
#' @param m number of sub-posteriors to combine
#' @param time_mesh time mesh used in Bayesian Fusion
#' @param precondition_matrices list of length m, where precondition_matrices[[c]]
#'                               is the precondition matrix for sub-posterior c
#' @param dim dimension
#' @param mean_vecs list of length m, where mean_vecs[[c]] is a vector of 
#'                  length dim for the mean of sub-posterior c
#' @param Sigmas list of length m, where Sigmas[[c]] is a dim x dim
#'                   inverse covariance matrix for sub-posterior c
#' @param resampling_method method to be used in resampling, default is multinomial 
#'                          resampling ('multi'). Other choices are stratified 
#'                          resampling ('strat'), systematic resampling ('system'),
#'                          residual resampling ('resid')
#' @param ESS_threshold number between 0 and 1 defining the proportion 
#'                      of the number of samples that ESS needs to be
#'                      lower than for resampling (i.e. resampling is carried 
#'                      out only when ESS < N*ESS_threshold)
#' @param sub_posterior_means matrix with m rows and dim columns, where sub_posterior_means[c,]
#'                            is the sub-posterior mean of sub-posterior c
#' @param adaptive_mesh logical value to indicate if an adaptive mesh is used
#'                      (default is FALSE)
#' @param adaptive_mesh_parameters list of parameters used for adaptive mesh
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
#'   \item{elapsed_time}{elapsed time of each step of the algorithm}
#'   \item{ESS}{effective sample size of the particles after each step}
#'   \item{CESS}{conditional effective sample size of the particles after each step}
#'   \item{resampled}{boolean value to indicate if particles were resampled
#'                    after each time step}
#'   \item{E_nu_j}{approximation of the average variation of the trajectories
#'                 at each time step}
#'   \item{chosen}{which term was chosen if using an adaptive mesh at each time step}
#'   \item{mesh_terms}{the evaluated terms in deciding the mesh at each time step}
#'   \item{precondition_matrices}{list of length 2 where precondition_matrices[[2]] 
#'                                are the pre-conditioning matrices that were used 
#'                                and precondition_matrices[[1]] are the combined 
#'                                precondition matrices}
#'   \item{sub_posterior_means}{list of length 2, where sub_posterior_means[[2]]
#'                              are the sub-posterior means that were used and
#'                              sub_posterior_means[[1]] are the combined
#'                              sub-posterior means}
#' }
#' 
#' @export
parallel_GBF_multiGaussian <- function(particles_to_fuse,
                                       N, 
                                       m,
                                       time_mesh,
                                       dim,
                                       mean_vecs,
                                       Sigmas, 
                                       precondition_matrices, 
                                       resampling_method = 'multi',
                                       ESS_threshold = 0.5,
                                       sub_posterior_means = NULL,
                                       adaptive_mesh = FALSE,
                                       adaptive_mesh_parameters = NULL,
                                       diffusion_estimator = 'Poisson',
                                       beta_NB = 10,
                                       gamma_NB_n_points = 2,
                                       seed = NULL,
                                       n_cores = parallel::detectCores()) {
  if (!is.list(particles_to_fuse) | (length(particles_to_fuse)!=m)) {
    stop("parallel_GBF_multiGaussian: particles_to_fuse must be a list of length m")
  } else if (!all(sapply(particles_to_fuse, function(sub_posterior) ("particle" %in% class(sub_posterior))))) {
    stop("parallel_GBF_multiGaussian: particles in particles_to_fuse must be \"particle\" objects")
  } else if (!all(sapply(particles_to_fuse, function(sub_posterior) is.matrix(sub_posterior$y_samples)))) {
    stop("parallel_GBF_multiGaussian: the particles' samples for y should all be matrices")
  } else if (!all(sapply(particles_to_fuse, function(sub_posterior) ncol(sub_posterior$y_samples)==dim))) {
    stop("parallel_GBF_multiGaussian: the particles' samples for y should all be matrices with dim columns")
  } else if (!is.vector(time_mesh)) {
    stop("parallel_GBF_multiGaussian: time_mesh must be an ordered vector of length >= 2")
  } else if (length(time_mesh) < 2) {
    stop("parallel_GBF_multiGaussian: time_mesh must be an ordered vector of length >= 2")
  } else if (!identical(time_mesh, sort(time_mesh))) {
    stop("parallel_GBF_multiGaussian: time_mesh must be an ordered vector of length >= 2")
  } else if (any(sapply(1:m, function(c) (!is.vector(mean_vecs[[c]]) | (length(mean_vecs[[c]])!=dim))))) {
    stop("parallel_GBF_multiGaussian: mean_vecs[[c]] must be a vector of length dim for each c")
  } else if (any(sapply(1:m, function(c) !is.matrix(Sigmas[[c]])))) {
    stop("parallel_GBF_multiGaussian: Sigmas[[c]] must be a dim x dim matrix for each c")
  } else if (!is.list(precondition_matrices) | (length(precondition_matrices)!=m)) {
    stop("parallel_GBF_multiGaussian: precondition_matrices must be a list of length m")
  } else if ((ESS_threshold < 0) | (ESS_threshold > 1)) {
    stop("parallel_GBF_multiGaussian: ESS_threshold must be between 0 and 1")
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
  inv_Sigmas <- lapply(Sigmas, solve)
  Lambda <- inverse_sum_matrices(inv_precondition_matrices)
  pcm_rho_0 <- proc.time()
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
                                   seed = seed,
                                   n_cores = n_cores)
  elapsed_time_rho_0 <- (proc.time()-pcm_rho_0)['elapsed']
  # ---------- iterative steps
  rho_j <- rho_j_multiGaussian(particle_set = particles,
                               m = m,
                               time_mesh = time_mesh,
                               dim = dim,
                               mean_vecs = mean_vecs,
                               inv_Sigmas = inv_Sigmas,
                               precondition_matrices = precondition_matrices,
                               inv_precondition_matrices = inv_precondition_matrices,
                               Lambda = Lambda,
                               resampling_method = resampling_method,
                               ESS_threshold = ESS_threshold,
                               sub_posterior_means = sub_posterior_means,
                               adaptive_mesh = adaptive_mesh,
                               adaptive_mesh_parameters = adaptive_mesh_parameters,
                               diffusion_estimator = diffusion_estimator,
                               beta_NB = beta_NB,
                               gamma_NB_n_points = gamma_NB_n_points,
                               seed = seed,
                               n_cores = n_cores)
  if (identical(precondition_matrices, rep(list(diag(1, dim)), m))) {
    new_precondition_matrices <- list(diag(1, dim), precondition_matrices)
  } else {
    new_precondition_matrices <- list(Lambda, precondition_matrices)
  }
  if (!is.null(sub_posterior_means)) {
    new_sub_posterior_means <- list(t(weighted_mean_multivariate(matrix = sub_posterior_means,
                                                                 weights = inv_precondition_matrices,
                                                                 inverse_sum_weights = Lambda)),
                                    sub_posterior_means)
  } else {
    new_sub_posterior_means <- list(NULL, sub_posterior_means)
  }
  # combine the means and Sigmas
  new_Sigma <- inverse_sum_matrices(inv_Sigmas)
  new_mean_vec <- as.vector(weighted_mean_multivariate(matrix = do.call(rbind, mean_vecs),
                                                       weights = inv_Sigmas,
                                                       inverse_sum_weights = new_Sigma))
  return(list('particles' = rho_j$particle_set,
              'proposed_samples' = rho_j$proposed_samples,
              'time' = (proc.time()-pcm)['elapsed'],
              'elapsed_time' = c(elapsed_time_rho_0, rho_j$time),
              'ESS' = rho_j$ESS,
              'CESS' = rho_j$CESS,
              'resampled' = rho_j$resampled,
              'E_nu_j' = rho_j$E_nu_j,
              'chosen' = rho_j$chosen,
              'mesh_terms' = rho_j$mesh_terms,
              'k4_choice' = rho_j$k4_choice,
              'precondition_matrices' = new_precondition_matrices,
              'sub_posterior_means' = new_sub_posterior_means,
              'new_mean_vec' = new_mean_vec,
              'new_Sigma' = new_Sigma))
}

#' (Balanced Binary) D&C Generalised Bayesian Fusion using SMC
#'
#' (Balanced Binary) D&C Generalised Bayesian Fusion with multivariate Gaussian target
#'
#' @param N_schedule vector of length (L-1), where N_schedule[l] is the number 
#'                   of samples per node at level l
#' @param m_schedule vector of length (L-1), where m_schedule[l] is the number 
#'                   of samples to fuse for level l
#' @param time_mesh time mesh used in Bayesian Fusion. This can either be a vector
#'                  which will be used for each node in the tree, or it can be
#'                  passed in as NULL, where a recommended mesh will be generated
#'                  using the parameters passed into mesh_parameters
#' @param base_samples list of length C, where base_samples[[c]] 
#'                     contains the samples for the c-th node in the level
#' @param L total number of levels in the hierarchy
#' @param dim dimension
#' @param mu vector of length dim for mean
#' @param Sigma dim x dim covariance matrix
#' @param precondition either a logical value to determine if preconditioning matrices are
#'                     used (TRUE - and is set to be the variance of the sub-posterior samples)
#'                     or not (FALSE - and is set to be the identity matrix for all sub-posteriors),
#'                     or a list of length C where precondition[[c]]
#'                     is the preconditioning matrix for sub-posterior c. Default is TRUE
#' @param resampling_method method to be used in resampling, default is multinomial 
#'                          resampling ('multi'). Other choices are stratified 
#'                          resampling ('strat'), systematic resampling ('system'),
#'                          residual resampling ('resid')
#' @param ESS_threshold number between 0 and 1 defining the proportion 
#'                      of the number of samples that ESS needs to be
#'                      lower than for resampling (i.e. resampling is carried 
#'                      out only when ESS < N*ESS_threshold)
#' @param adaptive_mesh logical value to indicate if an adaptive mesh is used
#'                      (default is FALSE)
#' @param mesh_parameters list of parameters used for mesh
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
#'   \item{time_mesh_used}{list of length (L-1), where time_mesh_used[[l]][[i]]
#'                         is the time_mesh that was used for level l, node i}
#'   \item{ESS}{list of length (L-1), where ESS[[l]][[i]] is the effective 
#'              sample size of the particles after each step BEFORE deciding 
#'              whether or not to resample for level l, node i}
#'   \item{CESS}{list of length (L-1), where CESS[[l]][[i]] is the conditional
#'               effective sample size of the particles after each step}
#'   \item{resampled}{list of length (L-1), where resampled[[l]][[i]] is a 
#'                    boolean value to record if the particles were resampled
#'                    after each step; rho and Q for level l, node i}
#'   \item{E_nu_j}{list of length (L-1), where E_nu_j[[l]][[i]] is the
#'                 approximation of the average variation of the trajectories
#'                 at each time step for level l, node i}
#'   \item{chosen}{list of length (L-1), where chosen[[l]][[i]] indicates
#'                 which term was chosen if using an adaptive mesh at each
#'                 time step for level l, node i}
#'   \item{mesh_terms}{list of length (L-1), where mesh_terms[[l]][[i]] indicates
#'                     the evaluated terms in deciding the mesh at each time step
#'                     for level l, node i}
#'   \item{k4_choice}{list of length (L-1), where k4_choice[[l]][[i]]] indicates 
#'                    which of the roots of k4 were chosen at each time step for
#'                    level l, node i}
#'   \item{precondition_matrices}{preconditioning matrices used in the algorithm 
#'                               for each node}
#'   \item{diffusion_times}{vector of length (L-1), where diffusion_times[l]
#'                         are the times for fusion in level l}
#' }
#'
#' @export
bal_binary_GBF_multiGaussian <- function(N_schedule,
                                         m_schedule,
                                         time_mesh = NULL,
                                         base_samples,
                                         L,
                                         dim,
                                         mean_vecs,
                                         Sigmas, 
                                         C,
                                         precondition = TRUE,
                                         resampling_method = 'multi',
                                         ESS_threshold = 0.5,
                                         adaptive_mesh = FALSE,
                                         mesh_parameters = NULL,
                                         diffusion_estimator = 'Poisson',
                                         beta_NB = 10,
                                         gamma_NB_n_points = 2,
                                         seed = NULL,
                                         n_cores = parallel::detectCores()) {
  if (!is.vector(N_schedule) | (length(N_schedule)!=(L-1))) {
    stop("bal_binary_GBF_multiGaussian: N_schedule must be a vector of length (L-1)")
  } else if (!is.vector(m_schedule) | (length(m_schedule)!=(L-1))) {
    stop("bal_binary_GBF_multiGaussian: m_schedule must be a vector of length (L-1)")
  } else if (!is.list(base_samples) | (length(base_samples)!=C)) {
    stop("bal_binary_GBF_multiGaussian: base_samples must be a list of length C")
  } else if (any(sapply(1:C, function(c) (!is.vector(mean_vecs[[c]]) | (length(mean_vecs[[c]])!=dim))))) {
    stop("parallel_GBF_multiGaussian: mean_vecs[[c]] must be a vector of length dim for each c")
  } else if (any(sapply(1:C, function(c) !is.matrix(Sigmas[[c]])))) {
    stop("parallel_GBF_multiGaussian: Sigmas[[c]] must be a dim x dim matrix for each c")
  } else if (ESS_threshold < 0 | ESS_threshold > 1) {
    stop("bal_binary_GBF_multiGaussian: ESS_threshold must be between 0 and 1")
  }
  if (is.vector(m_schedule) & (length(m_schedule)==(L-1))) {
    for (l in (L-1):1) {
      if ((C/prod(m_schedule[(L-1):l]))%%1!=0) {
        stop("bal_binary_GBF_multiGaussian: check that C/prod(m_schedule[(L-1):l])
              is an integer for l=L-1,...,1")
      }
    }
  } else {
    stop("bal_binary_GBF_multiGaussian: m_schedule must be a vector of length (L-1)")
  }
  if (is.vector(time_mesh)) {
    if (length(time_mesh) < 2) {
      stop("bal_binary_GBF_BLR: time_mesh must be an ordered vector of length >= 2")
    } else if (!identical(time_mesh, sort(time_mesh))) {
      stop("bal_binary_GBF_BLR: time_mesh must be an ordered vector of length >= 2")
    }
  } else if (is.null(time_mesh)) {
    if (!is.list(mesh_parameters)) {
      stop("bal_binary_GBF_BLR: if time_mesh is NULL, mesh_parameters must be a
           list of parameters to obtain guidance for the mesh")
    } 
  } else {
    stop("bal_binary_GBF_BLR: time_mesh must either be an ordered vector of length
           >= 2 or passed as NULL if want to use recommended guidance")
  }
  m_schedule <- c(m_schedule, 1)
  particles <- list()
  if (all(sapply(base_samples, function(sub) class(sub)=='particle'))) {
    particles[[L]] <- base_samples
  } else if (all(sapply(base_samples, is.matrix))) {
    if (!all(sapply(base_samples, function(core) ncol(core)==dim))) {
      stop("bal_binary_GBF_multiGaussian: the sub-posterior samples in base_samples must be matrices with dim columns")
    }
    particles[[L]] <- initialise_particle_sets(samples_to_fuse = base_samples,
                                               multivariate = TRUE,
                                               number_of_steps = 2)
  } else {
    stop("bal_binary_GBF_multiGaussian: base_samples must be a list of length
         C containing either items of class \"particle\" (representing
         particle approximations of the sub-posteriors) or are matrices with dim columns
         (representing un-normalised sample approximations of the sub-posteriors)")
  }
  proposed_samples <- list()
  time <- list()
  time_mesh_used <- list()
  ESS <- list()
  CESS <- list()
  resampled <- list()
  E_nu_j <- list()
  chosen <- list()
  mesh_terms <- list()
  k4_choice <- list()
  new_mean_vecs <- list()
  new_mean_vecs[[L]] <- mean_vecs
  new_Sigmas <- list()
  new_Sigmas[[L]] <- Sigmas
  recommended_mesh <- list()
  precondition_matrices <- list()
  if (is.logical(precondition)) {
    if (precondition) {
      precondition_matrices[[L]] <- lapply(base_samples, cov)
    } else {
      precondition_matrices[[L]] <- rep(list(diag(1, dim)), C)
    }
  } else if (is.list(precondition)) {
    if (length(precondition)==C & all(sapply(precondition, is.matrix))) {
      if (all(sapply(precondition, function(sub) ncol(sub)==dim))) {
        precondition_matrices[[L]] <- precondition  
      }
    }
  } else {
    stop("bal_binary_GBF_multiGaussian: precondition must be a logical indicating 
          whether or not a preconditioning matrix should be used, or a list of
          length C, where precondition[[c]] is the preconditioning matrix for
          the c-th sub-posterior")
  }
  sub_posterior_means <- list()
  if (dimension[d]==1) {
    sub_posterior_means[[L]] <- as.matrix(sapply(base_samples, function(sub) apply(sub, 2, mean)))
  } else {
    sub_posterior_means[[L]] <- t(sapply(base_samples, function(sub) apply(sub, 2, mean)))
  }
  cat('Starting bal_binary fusion \n', file = 'bal_binary_GBF_multiGaussian.txt')
  for (k in ((L-1):1)) {
    n_nodes <- max(C/prod(m_schedule[L:k]), 1)
    cat('########################\n', file = 'bal_binary_GBF_multiGaussian.txt',
        append = T)
    cat('There are', n_nodes, 'nodes at this level each giving', N_schedule[k],
        'samples \n', file = 'bal_binary_GBF_multiGaussian.txt', append = T)
    cat('########################\n', file = 'bal_binary_GBF_multiGaussian.txt', 
        append = T)
    fused <- lapply(X = 1:n_nodes, FUN = function(i) {
      previous_nodes <- ((m_schedule[k]*i)-(m_schedule[k]-1)):(m_schedule[k]*i)
      particles_to_fuse <- particles[[k+1]][previous_nodes]
      precondition_mats <- precondition_matrices[[k+1]][previous_nodes]
      inv_precondition_mats <- lapply(precondition_mats, solve)
      Lambda <- inverse_sum_matrices(inv_precondition_mats)
      sub_post_means <- sub_posterior_means[[k+1]][previous_nodes,,drop=FALSE]
      mean_vecs_to_fuse <- new_mean_vecs[[k+1]][previous_nodes]
      Sigmas_to_fuse <- new_Sigmas[[k+1]][previous_nodes]
      if (is.null(time_mesh)) {
        recommendation <- BF_guidance(condition = mesh_parameters$condition,
                                      CESS_0_threshold = mesh_parameters$CESS_0_threshold,
                                      CESS_j_threshold = mesh_parameters$CESS_j_threshold,
                                      sub_posterior_samples = lapply(1:length(previous_nodes), function(i) {
                                        particles_to_fuse[[i]]$y_samples}),
                                      log_weights = lapply(1:length(previous_nodes), function(i) {
                                        particles_to_fuse[[i]]$log_weights}),
                                      C = m_schedule[k],
                                      d = dim,
                                      data_size = mesh_parameters$data_size,
                                      b = mesh_parameters$b,
                                      sub_posterior_means = sub_post_means,
                                      precondition_matrices = precondition_mats,
                                      inv_precondition_matrices = inv_precondition_mats,
                                      Lambda = Lambda,
                                      lambda = mesh_parameters$lambda,
                                      gamma = mesh_parameters$gamma,
                                      k1 = mesh_parameters$k1,
                                      k2 = mesh_parameters$k2,
                                      k3 = mesh_parameters$k3,
                                      k4 = mesh_parameters$k4,
                                      vanilla = mesh_parameters$vanilla)
      } else {
        recommendation <- list('mesh' = time_mesh)
      }
      return(list('recommendation' = recommendation,
                  'fusion' = parallel_GBF_multiGaussian(particles_to_fuse = particles_to_fuse,
                                                        N = N_schedule[k],
                                                        m = m_schedule[k],
                                                        time_mesh = recommendation$mesh,
                                                        dim = dim,
                                                        mean_vecs = mean_vecs_to_fuse,
                                                        Sigmas = Sigmas_to_fuse,
                                                        precondition_matrices = precondition_mats,
                                                        resampling_method = resampling_method,
                                                        ESS_threshold = ESS_threshold,
                                                        sub_posterior_means = sub_post_means,
                                                        adaptive_mesh = adaptive_mesh,
                                                        adaptive_mesh_parameters = mesh_parameters,
                                                        diffusion_estimator = diffusion_estimator,
                                                        beta_NB = beta_NB,
                                                        gamma_NB_n_points = gamma_NB_n_points,
                                                        seed = seed,
                                                        n_cores = n_cores)))
    })
    # need to combine the correct samples
    particles[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$fusion$particles)
    proposed_samples[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$fusion$proposed_samples)
    time[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$fusion$time)
    time_mesh_used[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$fusion$particles$time_mesh)
    ESS[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$fusion$ESS)
    CESS[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$fusion$CESS)
    resampled[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$fusion$resampled)
    E_nu_j[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$fusion$E_nu_j)
    chosen[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$fusion$chosen)
    mesh_terms[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$fusion$mesh_terms)
    k4_choice[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$fusion$k4_choice)
    precondition_matrices[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$fusion$precondition_matrices[[1]])
    sub_posterior_means[[k]] <- do.call(rbind, lapply(1:n_nodes, function(i) fused[[i]]$fusion$sub_posterior_means[[1]]))
    new_mean_vecs[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$fusion$new_mean_vec)
    new_Sigmas[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$fusion$new_Sigma)
    recommended_mesh[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$recommendation)
  }
  cat('Completed bal_binary fusion\n', file = 'bal_binary_GBF_multiGaussian.txt', append = T)
  if (length(particles[[1]])==1) {
    particles[[1]] <- particles[[1]][[1]]
    proposed_samples[[1]] <- proposed_samples[[1]][[1]]
    time[[1]] <- time[[1]][[1]]
    time_mesh_used[[1]] <- time_mesh_used[[1]][[1]]
    ESS[[1]] <- ESS[[1]][[1]]
    CESS[[1]] <- CESS[[1]][[1]]
    resampled[[1]] <- resampled[[1]][[1]]
    E_nu_j[[1]] <- E_nu_j[[1]][[1]]
    chosen[[1]] <- chosen[[1]][[1]]
    mesh_terms[[1]] <- mesh_terms[[1]][[1]]
    k4_choice[[1]] <- k4_choice[[1]][[1]]
    precondition_matrices[[1]] <- precondition_matrices[[1]][[1]]
    sub_posterior_means[[1]] <- sub_posterior_means[[1]][[1]]
    new_mean_vecs[[1]] <- new_mean_vecs[[1]][[1]]
    new_Sigmas[[1]] <- new_Sigmas[[1]][[1]]
  }
  return(list('particles' = particles,
              'proposed_samples' = proposed_samples,
              'time' = time,
              'time_mesh_used' = time_mesh_used,
              'ESS' = ESS,
              'CESS' = CESS,
              'resampled' = resampled,
              'E_nu_j' = E_nu_j,
              'chosen' = chosen,
              'mesh_terms' = mesh_terms,
              'k4_choice' = k4_choice,
              'precondition_matrices' = precondition_matrices,
              'sub_posterior_means' = sub_posterior_means,
              'recommended_mesh' = recommended_mesh,
              'new_mean_vecs' = new_mean_vecs,
              'new_Sigmas' = new_Sigmas))
}

#' (Progressive) D&C Generalised Bayesian Fusion using SMC
#'
#' (Progressive) D&C Generalised Bayesian Fusion with multivariate Gaussian target
#'
#' @param N_schedule vector of length (L-1), where N_schedule[l] is the number 
#'                   of samples per node at level l
#' @param time_mesh time mesh used in Bayesian Fusion. This can either be a vector
#'                  which will be used for each node in the tree, or it can be
#'                  passed in as NULL, where a recommended mesh will be generated
#'                  using the parameters passed into mesh_parameters
#' @param base_samples list of length C, where base_samples[[c]] 
#'                     contains the samples for the c-th node in the level
#' @param dim dimension
#' @param mu vector of length dim for mean
#' @param Sigma dim x dim covariance matrix
#' @param precondition either a logical value to determine if preconditioning matrices are
#'                     used (TRUE - and is set to be the variance of the sub-posterior samples)
#'                     or not (FALSE - and is set to be the identity matrix for all sub-posteriors),
#'                     or a list of length C where precondition[[c]]
#'                     is the preconditioning matrix for sub-posterior c. Default is TRUE
#' @param resampling_method method to be used in resampling, default is multinomial 
#'                          resampling ('multi'). Other choices are stratified 
#'                          resampling ('strat'), systematic resampling ('system'),
#'                          residual resampling ('resid')
#' @param ESS_threshold number between 0 and 1 defining the proportion 
#'                      of the number of samples that ESS needs to be
#'                      lower than for resampling (i.e. resampling is carried 
#'                      out only when ESS < N*ESS_threshold)
#' @param adaptive_mesh logical value to indicate if an adaptive mesh is used
#'                      (default is FALSE)
#' @param mesh_parameters list of parameters used for mesh
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
#'   \item{time_mesh_used}{list of length (L-1), where time_mesh_used[[l]][[i]]
#'                         is the time_mesh that was used for level l, node i}
#'   \item{ESS}{list of length (L-1), where ESS[[l]][[i]] is the effective 
#'             sample size of the particles after each step BEFORE deciding 
#'             whether or not to resample for level l, node i}
#'   \item{CESS}{list of length (L-1), where CESS[[l]][[i]] is the conditional
#'              effective sample size of the particles after each step}
#'   \item{resampled}{list of length (L-1), where resampled[[l]][[i]] is a 
#'                   boolean value to record if the particles were resampled
#'                   after each step; rho and Q for level l, node i}
#'   \item{E_nu_j}{list of length (L-1), where E_nu_j[[l]][[i]] is the
#'                 approximation of the average variation of the trajectories
#'                 at each time step for level l, node i}
#'   \item{chosen}{list of length (L-1), where chosen[[l]][[i]] indicates
#'                 which term was chosen if using an adaptive mesh at each
#'                 time step for level l, node i}
#'   \item{mesh_terms}{list of length (L-1), where mesh_terms[[l]][[i]] indicates
#'                     the evaluated terms in deciding the mesh at each time step
#'                     for level l, node i}
#'   \item{k4_choice}{list of length (L-1), where k4_choice[[l]][[i]]] indicates 
#'                    which of the roots of k4 were chosen at each time step for
#'                    level l, node i}
#'   \item{precondition_matrices}{preconditioning matrices used in the algorithm 
#'                               for each node}
#'   \item{diffusion_times}{vector of length (L-1), where diffusion_times[l]
#'                         are the times for fusion in level l}
#' }
#'
#' @export
progressive_GBF_multiGaussian <- function(N_schedule,
                                          time_mesh = NULL,
                                          base_samples,
                                          dim,
                                          mean_vecs,
                                          Sigmas, 
                                          C,
                                          precondition = TRUE,
                                          resampling_method = 'multi',
                                          ESS_threshold = 0.5,
                                          adaptive_mesh = FALSE,
                                          mesh_parameters = NULL,
                                          diffusion_estimator = 'Poisson',
                                          beta_NB = 10,
                                          gamma_NB_n_points = 2,
                                          seed = NULL,
                                          n_cores = parallel::detectCores()) {
  if (!is.vector(N_schedule) | (length(N_schedule)!=C-1)) {
    stop("progressive_GBF_multiGaussian: N_schedule must be a vector of length (C-1)")
  } else if (!is.list(base_samples) | (length(base_samples)!=C)) {
    stop("progressive_GBF_multiGaussian: base_samples must be a list of length C")
  } else if (any(sapply(1:C, function(c) (!is.vector(mean_vecs[[c]]) | (length(mean_vecs[[c]])!=dim))))) {
    stop("parallel_GBF_multiGaussian: mean_vecs[[c]] must be a vector of length dim for each c")
  } else if (any(sapply(1:C, function(c) !is.matrix(Sigmas[[c]])))) {
    stop("parallel_GBF_multiGaussian: Sigmas[[c]] must be a dim x dim matrix for each c")
  } else if (ESS_threshold < 0 | ESS_threshold > 1) {
    stop("progressive_GBF_multiGaussian: ESS_threshold must be between 0 and 1")
  }
  if (is.vector(time_mesh)) {
    if (length(time_mesh) < 2) {
      stop("progressive_GBF_multiGaussian: time_mesh must be an ordered vector of length >= 2")
    } else if (!identical(time_mesh, sort(time_mesh))) {
      stop("progressive_GBF_multiGaussian: time_mesh must be an ordered vector of length >= 2")
    }
  } else if (is.null(time_mesh)) {
    if (!is.list(mesh_parameters)) {
      stop("progressive_GBF_multiGaussian: if time_mesh is NULL, mesh_parameters must be a
           list of parameters to obtain guidance for the mesh")
    } 
  } else {
    stop("progressive_GBF_multiGaussian: time_mesh must either be an ordered vector of length
           >= 2 or passed as NULL if want to use recommended guidance")
  }
  particles <- list()
  if (all(sapply(base_samples, function(sub) class(sub)=='particle'))) {
    particles[[C]] <- base_samples
  } else if (all(sapply(base_samples, is.matrix))) {
    if (!all(sapply(base_samples, function(core) ncol(core)==dim))) {
      stop("progressive_GBF_multiGaussian: the sub-posterior samples in base_samples must be matrices with dim columns")
    }
    particles[[C]] <- initialise_particle_sets(samples_to_fuse = base_samples,
                                               multivariate = TRUE,
                                               number_of_steps = 2)
  } else {
    stop("progressive_GBF_multiGaussian: base_samples must be a list of length
         C containing either items of class \"particle\" (representing
         particle approximations of the sub-posteriors) zor are matrices (representing
         un-normalised sample approximations of the sub-posteriors)")
  }
  proposed_samples <- list()
  time <- list()
  time_mesh_used <- list()
  ESS <- list()
  CESS <- list()
  resampled <- list()
  E_nu_j <- list()
  chosen <- list()
  mesh_terms <- list()
  k4_choice <- list()
  new_mean_vecs <- list()
  new_mean_vecs[[C]] <- mean_vecs
  new_Sigmas <- list()
  new_Sigmas[[C]] <- Sigmas
  recommended_mesh <- list()
  precondition_matrices <- list()
  if (is.logical(precondition)) {
    if (precondition) {
      precondition_matrices[[C]] <- lapply(base_samples, cov)
    } else {
      precondition_matrices[[C]] <- rep(list(diag(1, dim)), C)
    }
  } else if (is.list(precondition)) {
    if (length(precondition)==C & all(sapply(precondition, is.matrix))) {
      if (all(sapply(precondition, function(sub) ncol(sub)==dim))) {
        precondition_matrices[[C]] <- precondition  
      }
    }
  } else {
    stop("progressive_GBF_multiGaussian: precondition must be a logical indicating 
          whether or not a preconditioning matrix should be used, or a list of
          length C, where precondition[[c]] is the preconditioning matrix for
          the c-th sub-posterior")
  }
  sub_posterior_means <- list()
  if (dimension[d]==1) {
    sub_posterior_means[[C]] <- as.matrix(sapply(input_samples, function(sub) apply(sub, 2, mean)))
  } else {
    sub_posterior_means[[C]] <- t(sapply(input_samples, function(sub) apply(sub, 2, mean)))
  }
  index <- 2
  cat('Starting progressive fusion \n', file = 'progressive_GBF_multiGaussian.txt')
  for (k in (C-1):1) {
    if (k==C-1) {
      cat('########################\n', file = 'progressive_GBF_multiGaussian.txt', 
          append = T)
      cat('Starting to fuse', 2, 'densities for level', k, 'which is using', 
          parallel::detectCores(), 'cores\n', 
          file = 'progressive_GBF_multiGaussian.txt', append = T)
      cat('########################\n', file = 'progressive_GBF_multiGaussian.txt',
          append = T)
      particles_to_fuse <- list(particles[[C]][[1]], 
                                particles[[C]][[2]])
      precondition_mats <- precondition_matrices[[k+1]][1:2]
      inv_precondition_mats <- lapply(precondition_mats, solve)
      Lambda <- inverse_sum_matrices(inv_precondition_mats)
      sub_post_means <- sub_posterior_means[[k+1]][1:2]
      mean_vecs_to_fuse <- new_mean_vecs[[k+1]][1:2]
      Sigmas_to_fuse <- new_Sigmas[[k+1]][1:2]
      if (is.null(time_mesh)) {
        recommendation <- BF_guidance(condition = mesh_parameters$condition,
                                      CESS_0_threshold = mesh_parameters$CESS_0_threshold,
                                      CESS_j_threshold = mesh_parameters$CESS_j_threshold,
                                      sub_posterior_samples = lapply(1:2, function(i) particles_to_fuse[[i]]$y_samples),
                                      log_weights = lapply(1:2, function(i) particles_to_fuse[[i]]$log_weights),
                                      C = m_schedule[k],
                                      d = dim,
                                      data_size = mesh_parameters$data_size,
                                      b = mesh_parameters$b,
                                      sub_posterior_means = sub_post_means,
                                      precondition_matrices = precondition_mats,
                                      inv_precondition_matrices = inv_precondition_mats,
                                      Lambda = Lambda,
                                      lambda = mesh_parameters$lambda,
                                      gamma = mesh_parameters$gamma,
                                      k1 = mesh_parameters$k1,
                                      k2 = mesh_parameters$k2,
                                      k3 = mesh_parameters$k3,
                                      k4 = mesh_parameters$k4,
                                      vanilla = mesh_parameters$vanilla)
      } else {
        recommendation <- list('mesh' = time_mesh)
      }
      fused <- parallel_GBF_multiGaussian(particles_to_fuse = particles_to_fuse,
                                          N = N_schedule[k],
                                          m = 2,
                                          time_mesh = recommendation$mesh,
                                          dim = dim,
                                          mean_vecs = mean_vecs_to_fuse,
                                          Sigmas = Sigmas_to_fuse,
                                          precondition_matrices = precondition_mats,
                                          resampling_method = resampling_method,
                                          ESS_threshold = ESS_threshold,
                                          sub_posterior_means = sub_post_means,
                                          adaptive_mesh = adaptive_mesh,
                                          adaptive_mesh_parameters = mesh_parameters,
                                          diffusion_estimator = diffusion_estimator,
                                          beta_NB = beta_NB,
                                          gamma_NB_n_points = gamma_NB_n_points,
                                          seed = seed,
                                          n_cores = n_cores)
    } else {
      cat('########################\n', file = 'progressive_GBF_multiGaussian.txt', 
          append = T)
      cat('Starting to fuse', 2, 'densities for level', k, 'which is using', 
          parallel::detectCores(), 'cores\n',
          file = 'progressive_GBF_multiGaussian.txt', append = T)
      cat('########################\n', file = 'progressive_GBF_multiGaussian.txt', 
          append = T)
      particles_to_fuse <- list(particles[[k+1]], 
                                particles[[C]][[index+1]])
      precondition_mats <- list(precondition_matrices[[k+1]],
                                precondition_matrices[[C]][[index+1]])
      inv_precondition_mats <- lapply(precondition_mats, solve)
      Lambda <- inverse_sum_matrices(inv_precondition_mats)
      sub_post_means <- list(sub_posterior_means[[k+1]],
                             sub_posterior_means[[C]][[index+1]])
      mean_vecs_to_fuse <- list(new_mean_vecs[[k+1]], new_mean_vecs[[C]][[index+1]])
      Sigmas_to_fuse <- list(new_Sigmas[[k+1]], new_Sigmas[[C]][[index+1]])
      if (is.null(time_mesh)) {
        recommendation <- BF_guidance(condition = mesh_parameters$condition,
                                      CESS_0_threshold = mesh_parameters$CESS_0_threshold,
                                      CESS_j_threshold = mesh_parameters$CESS_j_threshold,
                                      sub_posterior_samples = lapply(1:2, function(i) particles_to_fuse[[i]]$y_samples),
                                      log_weights = lapply(1:2, function(i) particles_to_fuse[[i]]$log_weights),
                                      C = m_schedule[k],
                                      d = dim,
                                      data_size = mesh_parameters$data_size,
                                      b = mesh_parameters$b,
                                      sub_posterior_means = sub_post_means,
                                      precondition_matrices = precondition_mats,
                                      inv_precondition_matrices = inv_precondition_mats,
                                      Lambda = Lambda,
                                      lambda = mesh_parameters$lambda,
                                      gamma = mesh_parameters$gamma,
                                      k1 = mesh_parameters$k1,
                                      k2 = mesh_parameters$k2,
                                      k3 = mesh_parameters$k3,
                                      k4 = mesh_parameters$k4,
                                      vanilla = mesh_parameters$vanilla)
      } else {
        recommendation <- list('mesh' = time_mesh)
      }
      fused <- parallel_GBF_multiGaussian(particles_to_fuse = particles_to_fuse,
                                          N = N_schedule[k],
                                          m = 2,
                                          time_mesh = recommendation$mesh,
                                          dim = dim,
                                          mean_vecs = mean_vecs_to_fuse,
                                          Sigmas = Sigmas_to_fuse,
                                          precondition_matrices = precondition_mats,
                                          resampling_method = resampling_method,
                                          ESS_threshold = ESS_threshold,
                                          sub_posterior_means = sub_post_means,
                                          adaptive_mesh = adaptive_mesh,
                                          adaptive_mesh_parameters = mesh_parameters,
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
    time_mesh_used[[k]] <- fused$particles$time_mesh
    ESS[[k]] <- fused$ESS
    CESS[[k]] <- fused$CESS
    resampled[[k]] <- fused$resampled
    E_nu_j[[k]] <- fused$E_nu_j
    chosen[[k]] <- fused$chosen
    mesh_terms[[k]] <- fused$chosen
    k4_choice[[k]] <- fused$k4_choice
    precondition_matrices[[k]] <- fused$precondition_matrices[[1]]
    sub_posterior_means[[k]] <- fused$sub_posterior_means[[1]]
    new_mean_vecs[[k]] <- fused$new_mean_vec
    new_Sigmas[[k]] <- fused$new_Sigma
    recommended_mesh[[k]] <- recommendation
  }
  cat('Completed progressive fusion\n', file = 'progressive_GBF_multiGaussian.txt', append = T)
  return(list('particles' = particles,
              'proposed_samples' = proposed_samples,
              'time' = time,
              'time_mesh_used' = time_mesh_used,
              'ESS' = ESS,
              'CESS' = CESS,
              'resampled' = resampled,
              'E_nu_j' = E_nu_j,
              'chosen' = chosen,
              'mesh_terms' = mesh_terms,
              'k4_choice' = k4_choice,
              'precondition_matrices' = precondition_matrices,
              'sub_posterior_means' = sub_posterior_means,
              'recommended_mesh' = recommended_mesh,
              'new_mean_vecs' = new_mean_vecs,
              'new_Sigmas' = new_Sigmas))
}
