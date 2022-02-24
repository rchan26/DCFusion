#' rho_j Importance Sampling Step
#'
#' rho_j Importance Sampling weighting for Bayesian logistic regression
#'
#' @param particle_set particles set prior to Q importance sampling step
#' @param m number of sub-posteriors to combine
#' @param time_mesh time mesh used in Bayesian Fusion
#' @param dim dimension of the predictors (= p+1)
#' @param data_split list of length m where each item is a list of length 4 where
#'                   for c=1,...,m, data_split[[c]]$y is the vector for y responses and
#'                   data_split[[c]]$x is the design matrix for the covariates for
#'                   sub-posterior c, data_split[[c]]$full_data_count is the unique
#'                   rows of the full data set with their counts and 
#'                   data_split[[c]]$design_count is the unique rows of the design
#'                   matrix and their counts
#' @param prior_means prior for means of predictors
#' @param prior_variances prior for variances of predictors
#' @param C overall number of sub-posteriors
#' @param precondition_matrices list of length m, where precondition_matrices[[c]]
#'                               is the precondition matrix for sub-posterior c
#' @param inv_precondition_matrices list of length m, where inv_precondition_matrices[[c]]
#'                                  is the inverse precondition matrix for sub-posterior c
#' @param Lambda inverse of the sum of the inverse precondition matrices (which
#'               can be computed using inverse_sum_matrices(inv_precondition_matrices))
#' @param cv_location string to determine what the location of the control variate
#'                    should be. Must be either 'mode' where the MLE estimator 
#'                    will be used or 'hypercube_centre' (default) to use the centre
#'                    of the simulated hypercube
#' @param sub_posterior_means matrix with m rows and d columns, where sub_posterior_means[c,]
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
#'   \item{particle_set}{updated particle set after the iterative rho_j steps}
#'   \item{proposed_samples}{proposal samples for the last time step}
#'   \item{time}{elapsed time of each step of the algorithm}
#'   \item{ESS}{effective sample size of the particles after each step}
#'   \item{CESS}{conditional effective sample size of the particles after each step}
#'   \item{resampled}{boolean value to indicate if particles were resampled
#'                    after each time step}
#'   \item{E_nu_j}{Approximation of the average variation of the trajectories
#'                 at each time step}
#' }
#' 
#' @export
rho_j_BLR <- function(particle_set,
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
                      sub_posterior_means = NULL,
                      adaptive_mesh = FALSE,
                      adaptive_mesh_parameters = NULL,
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
    stop("rho_j_BLR: particle_set must be a \"particle\" object")
  } else if (!is.list(data_split) | length(data_split)!=m) {
    stop("rho_j_BLR: data_split must be a list of length m")
  } else if (!all(sapply(data_split, function(sub_posterior) (is.list(sub_posterior) & identical(names(sub_posterior), c("y", "X", "full_data_count", "design_count")))))) {
    stop("rho_j_BLR: each item in data_split must be a list of length 4 with names \'y\', \'X\', \'full_data_count\', \'design_count\'")
  } else if (!is.vector(time_mesh)) {
    stop("rho_j_BLR: time_mesh must be an ordered vector of length >= 2")
  } else if (length(time_mesh) < 2) {
    stop("rho_j_BLR: time_mesh must be an ordered vector of length >= 2")
  } else if (!identical(time_mesh, sort(time_mesh))) {
    stop("rho_j_BLR: time_mesh must be an ordered vector of length >= 2")
  } else if (!all(sapply(1:m, function(i) is.vector(data_split[[i]]$y)))) {
    stop("rho_j_BLR: for each i in 1:m, data_split[[i]]$y must be a vector")
  } else if (!all(sapply(1:m, function(i) is.matrix(data_split[[i]]$X)))) {
    stop("rho_j_BLR: for each i in 1:m, data_split[[i]]$X must be a matrix")
  } else if (!all(sapply(1:m, function(i) ncol(data_split[[i]]$X)==dim))) {
    stop("rho_j_BLR: for each i in 1:m, ncol(data_split[[i]]$X) must be equal to dim")
  } else if (!all(sapply(1:m, function(i) length(data_split[[i]]$y)==nrow(data_split[[i]]$X)))) {
    stop("rho_j_BLR: for each i in 1:m, length(data_split[[i]]$y) and nrow(data_split[[i]]$X) must be equal")
  } else if (!all(sapply(1:m, function(i) is.data.frame(data_split[[i]]$full_data_count)))) {
    stop("rho_j_BLR: for each i in 1:m, data_split[[i]]$full_data_count must be a data frame")
  } else if (!all(sapply(1:m, function(i) is.data.frame(data_split[[i]]$design_count)))) {
    stop("rho_j_BLR: for each i in 1:m, data_split[[i]]$design_count must be a data frame")
  } else if (!is.vector(prior_means) | length(prior_means)!=dim) {
    stop("rho_j_BLR: prior_means must be vectors of length dim")
  } else if (!is.vector(prior_variances) | length(prior_variances)!=dim) {
    stop("rho_j_BLR: prior_variances must be vectors of length dim")
  } else if (!is.list(precondition_matrices) | (length(precondition_matrices)!=m)) {
    stop("rho_j_BLR: precondition_matrices must be a list of length m")
  } else if (!is.list(inv_precondition_matrices) | (length(inv_precondition_matrices)!=m)) {
    stop("rho_j_BLR: inv_precondition_matrices must be a list of length m")
  } else if (!(diffusion_estimator %in% c('Poisson', 'NB'))) {
    stop("rho_j_BLR: diffusion_estimator must be set to either \'Poisson\' or \'NB\'")
  } else if (!any(class(cl)=="cluster") & !is.null(cl)) {
    stop("rho_j_BLR: cl must be a \"cluster\" object or NULL")
  }
  if (adaptive_mesh) {
    if (!is.matrix(sub_posterior_means)) {
      stop("rho_j_BLR: if adaptive_mesh==TRUE, sub_posterior_means must be a (m x dim) matrix")
    } else if (any(dim(sub_posterior_means)!=c(m,dim))) {
      stop("rho_j_BLR: if adaptive_mesh==TRUE, sub_posterior_means must be a (m x dim) matrix")
    }
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
    stop("rho_j_BLR: cv_location must be either \"mode\" or \"hypercube_centre\"")
  }
  transform_matrices <- lapply(1:m, function(c) {
    list('to_Z' = expm::sqrtm(inv_precondition_matrices[[c]]),
         'to_X' = expm::sqrtm(precondition_matrices[[c]]))
  })
  N <- particle_set$N
  # ---------- creating parallel cluster
  if (is.null(cl)) {
    cl <- parallel::makeCluster(n_cores, setup_strategy = "sequential", outfile = "GBF_BLR_outfile.txt")
    parallel::clusterExport(cl, varlist = ls("package:DCFusion"))
    parallel::clusterExport(cl, varlist = ls("package:layeredBB"))
    close_cluster <- TRUE
  } else {
    close_cluster <- FALSE
  }
  parallel::clusterExport(cl, envir = environment(), varlist = ls())
  if (!is.null(seed)) {
    parallel::clusterSetRNGStream(cl, iseed = seed)
  }
  max_samples_per_core <- ceiling(N/n_cores)
  split_indices <- split(1:N, ceiling(seq_along(1:N)/max_samples_per_core))
  counts <- c('full_data_count', 'design_count')
  elapsed_time <- rep(NA, length(time_mesh)-1)
  ESS <- c(particle_set$ESS[1], rep(NA, length(time_mesh)-1))
  CESS <- c(particle_set$CESS[1], rep(NA, length(time_mesh)-1))
  resampled <- rep(FALSE, length(time_mesh))
  if (adaptive_mesh) {
    E_nu_j <- rep(NA, length(time_mesh))
    E_nu_j_old <- rep(NA, length(time_mesh))
  } else {
    E_nu_j <- NA
    E_nu_j_old <- NA
  }
  if (is.null(print_progress_iters)) {
    print_progress_iters <- split_N
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
                     d = dim,
                     precondition_matrices = precondition_matrices,
                     Lambda = Lambda)
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
                         precondition_matrices = precondition_matrices,
                         sub_posterior_samples = split_x_samples[[core]][[i]],
                         sub_posterior_mean = split_x_means[[core]][i,])$M
        if (time_mesh[j]!=end_time) {
          return(matrix(mvrnormArma(N = 1, mu = M, Sigma = V), nrow = m, ncol = dim, byrow = TRUE))
        } else {
          return(matrix(mvtnorm::rmvnorm(n = 1, mean = M, sigma = V), nrow = m, ncol = dim, byrow = TRUE))
        }
      })
      if (core == 1) {
        cat('##### t_{j-1}:', time_mesh[j-1], '|| t_{j}:', time_mesh[j], '|| T:', 
            end_time, '#####\n', file = 'rho_j_BLR_progress.txt', append = T)  
      }
      cat('Level:', level, '|| Step:', j, '/', length(time_mesh), '|| Node:', node,
          '|| Core:', core, '|| START \n', file = 'rho_j_BLR_progress.txt', append = T)
      for (i in 1:split_N) {
        x_mean_j[i,] <- weighted_mean_multivariate(matrix = x_j[[i]],
                                                   weights = inv_precondition_matrices,
                                                   inverse_sum_weights = Lambda)
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
              file = 'rho_j_BLR_progress.txt', append = T)
        }
      }
      cat('Level:', level, '|| Step:', j, '|| Node:', node, '|| Core:', core, '||', split_N, '/',
          split_N, '\n', file = 'rho_j_BLR_progress.txt', append = T)
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
  if (close_cluster) {
    parallel::stopCluster(cl)
  }
  if (adaptive_mesh) {
    CESS <- CESS[1:j]
    ESS <- ESS[1:j]
    resampled <- resampled[1:j]
    particle_set$time_mesh <- time_mesh[1:j]
    elapsed_time <- elapsed_time[1:(j-1)]
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
              'time' = elapsed_time,
              'ESS' = ESS,
              'CESS' = CESS,
              'resampled' = resampled,
              'E_nu_j' = E_nu_j,
              'E_nu_j_old' = E_nu_j_old))
}

#' Generalised Bayesian Fusion [parallel]
#' 
#' Generalised Bayesian Fusion for Bayesian Logistic Regression
#'
#' @param particles_to_fuse list of length m, where particles_to_fuse[c] contains
#'                          the particles for the c-th sub-posterior. Can
#'                          initialise a this from list of sub-posterior samples
#'                          by using the initialise_particle_sets function
#' @param N number of samples
#' @param m number of sub-posteriors to combine
#' @param time_mesh time mesh used in Bayesian Fusion
#' @param dim dimension of the predictors (= p+1)
#' @param data_split list of length m where each item is a list of length 4 where
#'                   for c=1,...,m, data_split[[c]]$y is the vector for y responses and
#'                   data_split[[c]]$x is the design matrix for the covariates for
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
#' @param sub_posterior_means matrix with m rows and d columns, where sub_posterior_means[c,]
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
#'   \item{elapsed_time}{elapsed time of each step of the algorithm}
#'   \item{ESS}{effective sample size of the particles after each step}
#'   \item{CESS}{conditional effective sample size of the particles after each step}
#'   \item{resampled}{boolean value to indicate if particles were resampled
#'                    after each time step}
#'   \item{E_nu_j}{Approximation of the average variation of the trajectories
#'                 at each time step}
#'   \item{precondition_matrices}{list of length 2 where precondition_matrices[[2]]
#'                                are the pre-conditioning matrices that were used
#'                                and precondition_matrices[[1]] are the combined
#'                                precondition matrices}
#'   \item{sub_posterior_means}{list of length 2, where sub_posterior_means[[2]]
#'                              are the sub-posterior means that were used and
#'                              sub_posterior_means[[1]] are the combined
#'                              sub-posterior means}
#'   \item{combined_data}{combined data for the fusion density}
#' }
#'
#' @export
parallel_GBF_BLR <- function(particles_to_fuse,
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
                             sub_posterior_means = NULL,
                             adaptive_mesh = FALSE,
                             adaptive_mesh_parameters = NULL,
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
    stop("parallel_generalised_BF_SMC_BLR: precondition_matrices must be a list of length m")
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
  # ---------- creating parallel cluster
  if (is.null(cl)) {
    cl <- parallel::makeCluster(n_cores, setup_strategy = "sequential", outfile = "GBF_BLR_outfile.txt")
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
  # ---------- first importance sampling step
  # pre-calculating the inverse precondition matrices
  inv_precondition_matrices <- lapply(precondition_matrices, solve)
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
                                   n_cores = n_cores,
                                   cl = cl)
  elapsed_time_rho_0 <- (proc.time()-pcm_rho_0)['elapsed']
  # ---------- iterative steps
  rho_j <- rho_j_BLR(particle_set = particles,
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
                     sub_posterior_means = sub_posterior_means,
                     adaptive_mesh = adaptive_mesh,
                     adaptive_mesh_parameters = adaptive_mesh_parameters,
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
  if (close_cluster) {
    parallel::stopCluster(cl)
  }
  if (identical(precondition_matrices, rep(list(diag(1, dim)), m))) {
    new_precondition_matrices <- list(diag(1, dim), precondition_matrices)
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
              'elapsed_time' = c(elapsed_time_rho_0, rho_j$time),
              'ESS' = rho_j$ESS,
              'CESS' = rho_j$CESS,
              'resampled' = rho_j$resampled,
              'E_nu_j' = rho_j$E_nu_j,
              'E_nu_j_old' = rho_j$E_nu_j_old,
              'precondition_matrices' = new_precondition_matrices,
              'sub_posterior_means' = new_sub_posterior_means,
              'combined_data' = combine_data(list_of_data = data_split, dim = dim)))
}
