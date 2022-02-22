#' rho_j Importance Sampling Step
#'
#' rho_j Importance Sampling weighting for bivariate Gaussian distributions
#'
#' @param particle_set particles set prior to Q importance sampling step
#' @param m number of sub-posteriors to combine
#' @param time_mesh time mesh used in Bayesian Fusion
#' @param mean_vecs list of length m, where mean_vecs[[c]] is a vector of 
#'                  length 2 for the mean of sub-posterior c
#' @param sd_vecs list of length m, where sd_vecs[[c]] is a vector of length 2
#'                for the standard deviation of sub-posterior c
#' @param corrs vector of length m, where corrs[c] give the correlation value
#'              between component 1 and component 2 for sub-posterior c
#' @param betas vector of length m, where betas[c] is the inverse temperature 
#'              value for c-th posterior
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
#' @param sub_posterior_means matrix with m rows and 2 columns, where sub_posterior_means[c,]
#'                            is the sub-posterior mean of sub-posterior c
#' @param adaptive_mesh logical value to indicate if an adaptive mesh is used
#'                      (default is FALSE)
#' @param adaptive_mesh_parameters list of parameters used for adaptive mesh.
#'                                 Items which can be included are data_size,
#'                                 b (population variance if vanilla BF is used),
#'                                 k3, k4 (determines CESS_j threshold), T2
#'                                 (regular mesh recommended), vanilla (logical
#'                                 value to indicate if vanilla BF guidance is used)
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
#'   \item{ESS}{effective sample size of the particles after each step}
#'   \item{CESS}{conditional effective sample size of the particles after each step}
#'   \item{resampled}{boolean value to indicate if particles were resampled
#'                    after each time step}
#' }
#' 
#' @export
rho_j_biGaussian <- function(particle_set,
                             m,
                             time_mesh,
                             mean_vecs,
                             sd_vecs,
                             corrs,
                             betas,
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
    stop("rho_j_biGaussian: particle_set must be a \"particle\" object")
  } else if (!is.vector(time_mesh)) {
    stop("rho_j_biGaussian: time_mesh must be an ordered vector of length >= 2")
  } else if (length(time_mesh) < 2) {
    stop("rho_j_biGaussian: time_mesh must be an ordered vector of length >= 2")
  } else if (!identical(time_mesh, sort(time_mesh))) {
    stop("rho_j_biGaussian: time_mesh must be an ordered vector of length >= 2")
  } else if (any(sapply(1:m, function(c) (!is.vector(mean_vecs[[c]]) | (length(mean_vecs[[c]])!=2))))) {
    stop("rho_j_biGaussian: mean_vecs[[c]] must be a vector of length 2 for each c")
  } else if (any(sapply(1:m, function(c) (!is.vector(sd_vecs[[c]]) | (length(sd_vecs[[c]])!=2))))) {
    stop("rho_j_biGaussian: sd_vecs[[c]] must be a vector of length 2 for each c")
  } else if (!is.vector(corrs) | (length(corrs)!=m)) {
    stop("rho_j_biGaussian: corrs must be a vector of length m")
  } else if (!is.vector(betas) | (length(betas)!=m)) {
    stop("rho_j_biGaussian: betas must be a vector of length m")
  } else if (!is.list(precondition_matrices) | (length(precondition_matrices)!=m)) {
    stop("rho_j_biGaussian: precondition_matrices must be a list of length m")
  } else if (!is.list(inv_precondition_matrices) | (length(inv_precondition_matrices)!=m)) {
    stop("rho_j_biGaussian: inv_precondition_matrices must be a list of length m")
  } else if ((ESS_threshold < 0) | (ESS_threshold > 1)) {
    stop("rho_j_biGaussian: ESS_threshold must be between 0 and 1")
  }
  if (adaptive_mesh) {
    if (!is.matrix(sub_posterior_means)) {
      stop("rho_j_biGaussian: if adaptive_mesh==TRUE, sub_posterior_means must be a (m x 2) matrix")
    } else if (any(dim(sub_posterior_means)!=c(m,2))) {
      stop("rho_j_biGaussian: if adaptive_mesh==TRUE, sub_posterior_means must be a (m x 2) matrix")
    }
  }
  transform_matrices <- lapply(1:m, function(c) {
    list('to_Z' = expm::sqrtm(inv_precondition_matrices[[c]]),
         'to_X' = expm::sqrtm(precondition_matrices[[c]]))
  })
  N <- particle_set$N
  # ---------- creating parallel cluster
  cl <- parallel::makeCluster(n_cores, setup_strategy = "sequential")
  parallel::clusterExport(cl, envir = environment(), varlist = ls())
  parallel::clusterExport(cl, varlist = ls("package:DCFusion"))
  parallel::clusterExport(cl, varlist = ls("package:layeredBB"))
  if (!is.null(seed)) {
    parallel::clusterSetRNGStream(cl, iseed = seed)
  }
  max_samples_per_core <- ceiling(N/n_cores)
  split_indices <- split(1:N, ceiling(seq_along(1:N)/max_samples_per_core))
  counts <- c('full_data_count', 'design_count')
  ESS <-  c(particle_set$ESS[1], rep(NA, length(time_mesh)-1))
  CESS <- c(particle_set$CESS[1], rep(NA, length(time_mesh)-1))
  resampled <- rep(FALSE, length(time_mesh))
  if (adaptive_mesh) {
    E_nu_j <- rep(NA, length(time_mesh))
    E_nu_j_old <- rep(NA, length(time_mesh))
  } else {
    E_nu_j <- NA
    E_nu_j_old <- NA
  }
  # iterative proposals
  end_time <- time_mesh[length(time_mesh)]
  j <- 1
  while (time_mesh[j]!=end_time) {
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
                                              d = 2,
                                              data_size = adaptive_mesh_parameters$data_size,
                                              b = adaptive_mesh_parameters$b,
                                              particle_set = particle_set,
                                              sub_posterior_means = sub_posterior_means,
                                              inv_precondition_matrices = inv_precondition_matrices,
                                              k3 = adaptive_mesh_parameters$k3,
                                              k4 = adaptive_mesh_parameters$k4,
                                              T2 = adaptive_mesh_parameters$T2,
                                              vanilla = adaptive_mesh_parameters$vanilla)
      if (is.null(adaptive_mesh_parameters$T2)) {
        adaptive_mesh_parameters$T2 <- tilde_Delta_j$T2
      }
      E_nu_j[j] <- tilde_Delta_j$E_nu_j
      E_nu_j_old[j] <- tilde_Delta_j$E_nu_j_old
      time_mesh[j] <- min(end_time, time_mesh[j-1]+tilde_Delta_j$max_delta_j)
    }
    # split the x samples from the previous time marginal (and their means) into approximately equal lists
    split_x_samples <- lapply(split_indices, function(indices) particle_set$x_samples[indices])
    split_x_means <- lapply(split_indices, function(indices) particle_set$x_means[indices,,drop = FALSE])
    V <- construct_V(s = time_mesh[j-1],
                     t = time_mesh[j],
                     end_time = end_time,
                     C = m,
                     d = 2,
                     precondition_matrices = precondition_matrices,
                     Lambda = Lambda)
    rho_j_weighted_samples <- parallel::parLapply(cl, X = 1:length(split_indices), fun = function(core) {
      split_N <- length(split_indices[[core]])
      x_mean_j <- matrix(data = NA, nrow = split_N, ncol = 2)
      log_rho_j <- rep(0, split_N)
      x_j <- lapply(1:split_N, function(i) {
        M <- construct_M(s = time_mesh[j-1],
                         t = time_mesh[j],
                         end_time = end_time,
                         C = m,
                         d = 2,
                         precondition_matrices = precondition_matrices,
                         sub_posterior_samples = split_x_samples[[core]][[i]],
                         sub_posterior_mean = split_x_means[[core]][i,])$M
        if (time_mesh[j]!=end_time) {
          return(matrix(mvrnormArma(N = 1, mu = M, Sigma = V), nrow = m, ncol = 2, byrow = TRUE))
        } else {
          return(matrix(mvtnorm::rmvnorm(n = 1, mean = M, sigma = V), nrow = m, ncol = 2, byrow = TRUE))
        }
      })
      for (i in 1:split_N) {
        x_mean_j[i,] <- weighted_mean_multivariate(matrix = x_j[[i]],
                                                   weights = inv_precondition_matrices,
                                                   inverse_sum_weights = Lambda)
        log_rho_j[i] <- sum(sapply(1:m, function(c) {
          ea_biGaussian_DL_PT(x0 = as.vector(split_x_samples[[core]][[i]][c,]),
                              y = as.vector(x_j[[i]][c,]),
                              s = time_mesh[j-1],
                              t = time_mesh[j],
                              mean_vec = mean_vecs[[c]],
                              sd_vec = sd_vecs[[c]],
                              corr = corrs[c],
                              beta = betas[c],
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
  }
  parallel::stopCluster(cl)
  if (adaptive_mesh) {
    CESS <- CESS[1:j]
    ESS <- ESS[1:j]
    resampled <- resampled[1:j]
    particle_set$time_mesh <- time_mesh[1:j]
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
              'resampled' = resampled,
              'E_nu_j' = E_nu_j,
              'E_nu_j_old' = E_nu_j_old))
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
#' @param mean_vecs list of length m, where mean_vecs[[c]] is a vector of 
#'                  length 2 for the mean of sub-posterior c
#' @param sd_vecs list of length m, where sd_vecs[[c]] is a vector of length 2
#'                for the standard deviation of sub-posterior c
#' @param corrs vector of length m, where corrs[c] give the correlation value
#'              between component 1 and component 2 for sub-posterior c
#' @param betas vector of length m, where betas[c] is the inverse temperature 
#'              value for c-th posterior
#' @param resampling_method method to be used in resampling, default is multinomial 
#'                          resampling ('multi'). Other choices are stratified 
#'                          resampling ('strat'), systematic resampling ('system'),
#'                          residual resampling ('resid')
#' @param ESS_threshold number between 0 and 1 defining the proportion 
#'                      of the number of samples that ESS needs to be
#'                      lower than for resampling (i.e. resampling is carried 
#'                      out only when ESS < N*ESS_threshold)
#' @param sub_posterior_means matrix with m rows and 2 columns, where sub_posterior_means[c,]
#'                            is the sub-posterior mean of sub-posterior c
#' @param adaptive_mesh logical value to indicate if an adaptive mesh is used
#'                      (default is FALSE)
#' @param adaptive_mesh_parameters list of parameters used for adaptive mesh.
#'                                 Items which can be included are data_size,
#'                                 b (population variance if vanilla BF is used),
#'                                 k3, k4 (determines CESS_j threshold), T2
#'                                 (regular mesh recommended), vanilla (logical
#'                                 value to indicate if vanilla BF guidance is used)
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
parallel_GBF_biGaussian <- function(particles_to_fuse,
                                    N, 
                                    m,
                                    time_mesh,
                                    mean_vecs,
                                    sd_vecs,
                                    corrs, 
                                    betas,
                                    precondition_matrices, 
                                    resampling_method = 'multi',
                                    ESS_threshold = 0.5,
                                    diffusion_estimator = 'Poisson',
                                    sub_posterior_means = NULL,
                                    adaptive_mesh = FALSE,
                                    adaptive_mesh_parameters = NULL,
                                    beta_NB = 10,
                                    gamma_NB_n_points = 2,
                                    seed = NULL,
                                    n_cores = parallel::detectCores()) {
  if (!is.list(particles_to_fuse) | (length(particles_to_fuse)!=m)) {
    stop("parallel_GBF_biGaussian: particles_to_fuse must be a list of length m")
  } else if (!all(sapply(particles_to_fuse, function(sub_posterior) ("particle" %in% class(sub_posterior))))) {
    stop("parallel_GBF_biGaussian: particles in particles_to_fuse must be \"particle\" objects")
  } else if (!all(sapply(particles_to_fuse, function(sub_posterior) is.matrix(sub_posterior$y_samples)))) {
    stop("parallel_GBF_biGaussian: the particles' samples for y should all be matrices")
  } else if (!all(sapply(particles_to_fuse, function(sub_posterior) ncol(sub_posterior$y_samples)==2))) {
    stop("parallel_GBF_biGaussian: the particles' samples for y should all be matrices with 2 columns")
  } else if (!is.vector(time_mesh)) {
    stop("parallel_GBF_biGaussian: time_mesh must be an ordered vector of length >= 2")
  } else if (length(time_mesh) < 2) {
    stop("parallel_GBF_biGaussian: time_mesh must be an ordered vector of length >= 2")
  } else if (!identical(time_mesh, sort(time_mesh))) {
    stop("parallel_GBF_biGaussian: time_mesh must be an ordered vector of length >= 2")
  } else if (any(sapply(1:m, function(c) (!is.vector(mean_vecs[[c]]) | (length(mean_vecs[[c]])!=2))))) {
    stop("parallel_GBF_biGaussian: mean_vecs[[c]] must be a vector of length 2 for each c")
  } else if (any(sapply(1:m, function(c) (!is.vector(sd_vecs[[c]]) | (length(sd_vecs[[c]])!=2))))) {
    stop("parallel_GBF_biGaussian: sd_vecs[[c]] must be a vector of length 2 for each c")
  } else if (!is.vector(corrs) | (length(corrs)!=m)) {
    stop("parallel_GBF_biGaussian: corrs must be a vector of length m")
  } else if (!is.vector(betas) | (length(betas)!=m)) {
    stop("parallel_GBF_biGaussian: betas must be a vector of length m")
  } else if (!is.list(precondition_matrices) | (length(precondition_matrices)!=m)) {
    stop("parallel_GBF_biGaussian: precondition_matrices must be a list of length m")
  } else if ((ESS_threshold < 0) | (ESS_threshold > 1)) {
    stop("parallel_GBF_biGaussian: ESS_threshold must be between 0 and 1")
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
                                   dim = 2,
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
  # ---------- iterative steps
  rho_j <- rho_j_biGaussian(particle_set = particles,
                            m = m,
                            time_mesh = time_mesh,
                            mean_vecs = mean_vecs,
                            sd_vecs = sd_vecs,
                            corrs = corrs,
                            betas = betas,
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
  if (identical(precondition_matrices, rep(list(diag(1, 2)), m))) {
    new_precondition_matrices <- list(diag(1, 2), precondition_matrices)
  } else {
    new_precondition_matrices <- list(inverse_sum_matrices(inv_precondition_matrices),
                                      precondition_matrices)
  }
  if (!is.null(sub_posterior_means)) {
    new_sub_posterior_means <- list(weighted_mean_multivariate(matrix = sub_posterior_means,
                                                               weights = inv_precondition_matrices,
                                                               inverse_sum_weights = Lambda),
                                    sub_posterior_means)
  } else {
    new_sub_posterior_means <- list(NULL, sub_posterior_means)
  }
  return(list('particles' = rho_j$particle_set,
              'proposed_samples' = rho_j$proposed_samples,
              'time' = (proc.time()-pcm)['elapsed'],
              'ESS' = rho_j$ESS,
              'CESS' = rho_j$CESS,
              'resampled' = rho_j$resampled,
              'E_nu_j' = rho_j$E_nu_j,
              'E_nu_j_old' = rho_j$E_nu_j_old,
              'precondition_matrices' = new_precondition_matrices,
              'sub_posterior_means' = new_sub_posterior_means))
}
