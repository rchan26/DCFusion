#' rho_j Importance Sampling Step
#'
#' rho_j Importance Sampling weighting for univariate Gaussian distributions
#'
#' @param particle_set particles set prior to Q importance sampling step
#' @param m number of sub-posteriors to combine
#' @param time_mesh time mesh used in Bayesian Fusion
#' @param means vector of length m, where means[c] is the mean for c-th 
#'              sub-posterior
#' @param sds vector of length m, where sds[c] is the standard deviation 
#'            for c-th sub-posterior
#' @param betas vector of length c, where betas[c] is the inverse temperature 
#'              value for c-th posterior
#' @param precondition_values vector of length m, where precondition_values[c]
#'                            is the precondition value for sub-posterior c
#' @param resampling_method method to be used in resampling, default is multinomial 
#'                          resampling ('multi'). Other choices are stratified 
#'                          resampling ('strat'), systematic resampling ('system'),
#'                          residual resampling ('resid')
#' @param ESS_threshold number between 0 and 1 defining the proportion 
#'                      of the number of samples that ESS needs to be
#'                      lower than for resampling (i.e. resampling is carried 
#'                      out only when ESS < N*ESS_threshold)
#' @param sub_posterior_means vector of length m, where sub_posterior_means[c]
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
#' }
#' 
#' @export
rho_j_uniGaussian <- function(particle_set,
                              m,
                              time_mesh,
                              means,
                              sds,
                              betas,
                              precondition_values,
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
    stop("rho_j_uniGaussian: particle_set must be a \"particle\" object")
  } else if (!is.vector(time_mesh)) {
    stop("rho_j_uniGaussian: time_mesh must be an ordered vector of length >= 2")
  } else if (length(time_mesh) < 2) {
    stop("rho_j_uniGaussian: time_mesh must be an ordered vector of length >= 2")
  } else if (!identical(time_mesh, sort(time_mesh))) {
    stop("rho_j_uniGaussian: time_mesh must be an ordered vector of length >= 2")
  } else if (!is.vector(means) | (length(means)!=m)) {
    stop("rho_j_uniGaussian: means must be a vector of length m")
  } else if (!is.vector(sds) | (length(sds)!=m)) {
    stop("rho_j_uniGaussian: sds must be a vector of length m")
  } else if (!is.vector(betas) | (length(betas)!=m)) {
    stop("rho_j_uniGaussian: betas must be a vector of length m")
  } else if (!is.vector(precondition_values) | (length(precondition_values)!=m)) {
    stop("rho_j_uniGaussian: precondition_values must be a vector of length m")
  } else if ((ESS_threshold < 0) | (ESS_threshold > 1)) {
    stop("rho_j_uniGaussian: ESS_threshold must be between 0 and 1")
  }
  if (adaptive_mesh) {
    if (!is.vector(sub_posterior_means)) {
      stop("rho_j_uniGaussian: if adaptive_mesh==TRUE, sub_posterior_means must be a vector of length m")
    } else if (length(sub_posterior_means)!=m) {
      stop("rho_j_uniGaussian: if adaptive_mesh==TRUE, sub_posterior_means must be a vector of length m")
    }
  }
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
  } else {
    E_nu_j <- NA
  }
  precondition_matrices <- lapply(precondition_values, as.matrix)
  Lambda <- inverse_sum_matrices(lapply(1/precondition_values, as.matrix))
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
                                                  multivariate = FALSE,
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
                                              d = 1,
                                              data_size = adaptive_mesh_parameters$data_size,
                                              b = adaptive_mesh_parameters$b,
                                              threshold = adaptive_mesh_parameters$CESS_j_threshold,
                                              particle_set = particle_set,
                                              sub_posterior_means = sub_posterior_means,
                                              inv_precondition_matrices = 1/precondition_values,
                                              k3 = adaptive_mesh_parameters$k3,
                                              k4 = adaptive_mesh_parameters$k4,
                                              vanilla = adaptive_mesh_parameters$vanilla)
      E_nu_j[j] <- tilde_Delta_j$E_nu_j
      time_mesh[j] <- min(end_time, time_mesh[j-1]+tilde_Delta_j$max_delta_j)
    }
    # split the x samples from the previous time marginal (and their means) into approximately equal lists
    split_x_samples <- lapply(split_indices, function(indices) particle_set$x_samples[indices])
    split_x_means <- lapply(split_indices, function(indices) particle_set$x_means[indices])
    V <- construct_V(s = time_mesh[j-1],
                     t = time_mesh[j],
                     end_time = end_time,
                     C = m,
                     d = 1,
                     precondition_matrices = precondition_matrices,
                     Lambda = Lambda,
                     iteration = j)
    rho_j_weighted_samples <- parallel::parLapply(cl, X = 1:length(split_indices), fun = function(core) {
      split_N <- length(split_indices[[core]])
      x_mean_j <- rep(0, split_N)
      log_rho_j <- rep(0, split_N)
      x_j <- lapply(1:split_N, function(i) {
        M <- construct_M(s = time_mesh[j-1],
                         t = time_mesh[j],
                         end_time = end_time,
                         C = m,
                         d = 1,
                         sub_posterior_samples = split_x_samples[[core]][[i]],
                         sub_posterior_mean = split_x_means[[core]][i],
                         iteration = j)
        if (time_mesh[j]!=end_time) {
          return(as.vector(mvrnormArma(N = 1, mu = M, Sigma = V)))
        } else {
          return(as.vector(mvtnorm::rmvnorm(n = 1, mean = M, sigma = V)))
        }
      })
      for (i in 1:split_N) {
        x_mean_j[i] <- weighted_mean_univariate(x = x_j[[i]], weights = 1/precondition_values)
        log_rho_j[i] <- sum(sapply(1:m, function(c) {
          ea_uniGaussian_DL_PT(x0 = split_x_samples[[core]][[i]][c],
                               y = x_j[[i]][c],
                               s = time_mesh[j-1],
                               t = time_mesh[j],
                               mean = means[c],
                               sd = sds[c], 
                               beta = betas[c], 
                               precondition = precondition_values[c],
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
    particle_set$x_means <- unlist(lapply(1:length(split_indices), function(i) {
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
  proposed_samples <- sapply(1:N, function(i) particle_set$x_samples[[i]][1])
  particle_set$y_samples <- proposed_samples
  # ----------- resample particle_set (only resample if ESS < N*ESS_threshold)
  if (particle_set$ESS < N*ESS_threshold) {
    resampled[particle_set$number_of_steps] <- TRUE
    particle_set <- resample_particle_y_samples(N = N,
                                                particle_set = particle_set,
                                                multivariate = FALSE,
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
  }
  return(list('particle_set' = particle_set,
              'proposed_samples' = proposed_samples,
              'time' = elapsed_time,
              'ESS' = ESS,
              'CESS' = CESS,
              'resampled' = resampled,
              'E_nu_j' = E_nu_j))
}

#' Generalised Bayesian Fusion [parallel]
#'
#' Generalised Bayesian Fusion with univariate Gaussian target
#'
#' @param particles_to_fuse list of length m, where particles_to_fuse[[c]]
#'                          contains the particles for the c-th sub-posterior
#'                          (a list of particles to fuse can be initialised by 
#'                          initialise_particle_sets() function)
#' @param N number of samples
#' @param m number of sub-posteriors to combine
#' @param time_mesh time mesh used in Bayesian Fusion
#' @param means vector of length m, where means[c] is the mean for c-th 
#'              sub-posterior
#' @param sds vector of length m, where sds[c] is the standard deviation 
#'            for c-th sub-posterior
#' @param betas vector of length c, where betas[c] is the inverse temperature 
#'              value for c-th posterior
#' @param precondition_values vector of length m, where precondition_values[c]
#'                            is the precondition value for sub-posterior c
#' @param resampling_method method to be used in resampling, default is multinomial 
#'                          resampling ('multi'). Other choices are stratified 
#'                          resampling ('strat'), systematic resampling ('system'),
#'                          residual resampling ('resid')
#' @param ESS_threshold number between 0 and 1 defining the proportion 
#'                      of the number of samples that ESS needs to be
#'                      lower than for resampling (i.e. resampling is carried 
#'                      out only when ESS < N*ESS_threshold)
#' @param sub_posterior_means vector of length m, where sub_posterior_means[c]
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
#'   \item{precondition_values}{list of length 2 where precondition_values[[2]] 
#'                              are the pre-conditioning values that were used 
#'                              and precondition_values[[1]] are the combined 
#'                              precondition values}
#'   \item{sub_posterior_means}{list of length 2, where sub_posterior_means[[2]]
#'                              are the sub-posterior means that were used and
#'                              sub_posterior_means[[1]] are the combined
#'                              sub-posterior means}
#' }
#' 
#' @export
parallel_GBF_uniGaussian <- function(particles_to_fuse,
                                     N,
                                     m,
                                     time_mesh,
                                     means,
                                     sds,
                                     betas,
                                     precondition_values,
                                     resampling_method = 'multi',
                                     ESS_threshold = 0.5,
                                     sub_posterior_means = NULL,
                                     adaptive_mesh = FALSE,
                                     adaptive_mesh_parameters = list(),
                                     diffusion_estimator = 'Poisson',
                                     beta_NB = 10,
                                     gamma_NB_n_points = 2,
                                     seed = NULL,
                                     n_cores = parallel::detectCores()) {
  if (!is.list(particles_to_fuse) | (length(particles_to_fuse)!=m)) {
    stop("parallel_GBF_uniGaussian: particles_to_fuse must be a list of length m")
  } else if (!all(sapply(particles_to_fuse, function(sub_posterior) ("particle" %in% class(sub_posterior))))) {
    stop("parallel_GBF_uniGaussian: particles in particles_to_fuse must be \"particle\" objects")
  } else if (!is.vector(time_mesh)) {
    stop("parallel_GBF_uniGaussian: time_mesh must be an ordered vector of length >= 2")
  } else if (length(time_mesh) < 2) {
    stop("parallel_GBF_uniGaussian: time_mesh must be an ordered vector of length >= 2")
  } else if (!identical(time_mesh, sort(time_mesh))) {
    stop("parallel_GBF_uniGaussian: time_mesh must be an ordered vector of length >= 2")
  } else if (!is.vector(means) | (length(means)!=m)) {
    stop("parallel_GBF_uniGaussian: means must be a vector of length m")
  } else if (!is.vector(sds) | (length(sds)!=m)) {
    stop("parallel_GBF_uniGaussian: sds must be a vector of length m")
  } else if (!is.vector(betas) | (length(betas)!=m)) {
    stop("parallel_GBF_uniGaussian: betas must be a vector of length m")
  } else if (!is.vector(precondition_values) | (length(precondition_values)!=m)) {
    stop("parallel_GBF_uniGaussian: precondition_values must be a vector of length m")
  } else if ((ESS_threshold < 0) | (ESS_threshold > 1)) {
    stop("parallel_GBF_uniGaussian: ESS_threshold must be between 0 and 1")
  }
  # set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  # start time recording
  pcm <- proc.time()
  # ---------- first importance sampling step
  pcm_rho_0 <- proc.time()
  particles <- rho_IS_univariate(particles_to_fuse = particles_to_fuse,
                                 N = N,
                                 m = m,
                                 time = time_mesh[length(time_mesh)],
                                 precondition_values = precondition_values,
                                 number_of_steps = length(time_mesh),
                                 time_mesh = time_mesh,
                                 resampling_method = resampling_method,
                                 seed = seed,
                                 n_cores = n_cores)
  elapsed_time_rho_0 <- (proc.time()-pcm_rho_0)['elapsed']
  # ---------- iterative steps
  rho_j <- rho_j_uniGaussian(particle_set = particles,
                             m = m,
                             time_mesh = time_mesh,
                             means = means,
                             sds = sds,
                             betas = betas,
                             precondition_values = precondition_values,
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
  if (identical(precondition_values, rep(1, m))) {
    new_precondition_values <- list(1, precondition_values)
  } else {
    new_precondition_values <- list(1/sum(1/precondition_values), precondition_values)
  }
  if (!is.null(sub_posterior_means)) {
    new_sub_posterior_means <- list(weighted_mean_univariate(x = sub_posterior_means,
                                                             weights = 1/precondition_values),
                                    sub_posterior_means)
  } else {
    new_sub_posterior_means <- list(NULL, sub_posterior_means)
  }
  return(list('particles' = rho_j$particle_set,
              'proposed_samples' = rho_j$proposed_samples,
              'time' = (proc.time()-pcm)['elapsed'],
              'elapsed_time' = c(elapsed_time_rho_0, rho_j$time),
              'ESS' = rho_j$ESS,
              'CESS' = rho_j$CESS,
              'resampled' = rho_j$resampled,
              'E_nu_j' = rho_j$E_nu_j,
              'precondition_matrices' = new_precondition_values,
              'sub_posterior_means' = new_sub_posterior_means))
}
