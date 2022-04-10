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
#'                   data_split[[c]]$X is the design matrix for the covariates for
#'                   sub-posterior c
#' @param nu degrees of freedom in t-distribution
#' @param sigma scale parameter in t-distribution
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
#' @param adaptive_mesh_parameters list of parameters used for adaptive mesh
#' @param record logical value indicating if variables such as E[nu_j], chosen,
#'               mesh_terms and k4_choice should be recorded at each iteration
#'               and returned (see return variables for this function) - default
#'               is FALSE
#' @param diffusion_estimator choice of unbiased estimator for the Exact Algorithm
#'                            between "Poisson" (default) for Poisson estimator
#'                            and "NB" for Negative Binomial estimator
#' @param beta_NB beta parameter for Negative Binomial estimator (default 10)
#' @param gamma_NB_n_points number of points used in the trapezoidal estimation
#'                          of the integral found in the mean of the negative
#'                          binomial estimator (default is 2)
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
#' }
#' If record is set to TRUE, additional components will be returned:
#' \describe{
#'   \item{E_nu_j}{approximation of the average variation of the trajectories
#'                 at each time step}
#'   \item{chosen}{which term was chosen if using an adaptive mesh at each time step}
#'   \item{mesh_terms}{the evaluated terms in deciding the mesh at each time step}
#'   \item{k4_choice}{which of the roots of k4 were chosen}
#' }
#' 
#' @export
rho_j_BRR <- function(particle_set,
                      m,
                      time_mesh,
                      dim,
                      data_split,
                      nu,
                      sigma,
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
                      record = FALSE,
                      diffusion_estimator,
                      beta_NB = 10,
                      gamma_NB_n_points = 2,
                      seed = NULL,
                      n_cores = parallel::detectCores(),
                      cl = NULL,
                      level = 1,
                      node = 1,
                      print_progress_iters = 1000) {
  if (!("particle" %in% class(particle_set))) {
    stop("rho_j_BRR: particle_set must be a \"particle\" object")
  } else if (!is.list(data_split) | length(data_split)!=m) {
    stop("rho_j_BRR: data_split must be a list of length m")
  } else if (!is.vector(time_mesh)) {
    stop("rho_j_BRR: time_mesh must be an ordered vector of length >= 2")
  } else if (length(time_mesh) < 2) {
    stop("rho_j_BRR: time_mesh must be an ordered vector of length >= 2")
  } else if (!identical(time_mesh, sort(time_mesh))) {
    stop("rho_j_BRR: time_mesh must be an ordered vector of length >= 2")
  } else if (!all(sapply(1:m, function(i) is.vector(data_split[[i]]$y)))) {
    stop("rho_j_BRR: for each i in 1:m, data_split[[i]]$y must be a vector")
  } else if (!all(sapply(1:m, function(i) is.matrix(data_split[[i]]$X)))) {
    stop("rho_j_BRR: for each i in 1:m, data_split[[i]]$X must be a matrix")
  } else if (!all(sapply(1:m, function(i) ncol(data_split[[i]]$X)==dim))) {
    stop("rho_j_BRR: for each i in 1:m, ncol(data_split[[i]]$X) must be equal to dim")
  } else if (!all(sapply(1:m, function(i) length(data_split[[i]]$y)==nrow(data_split[[i]]$X)))) {
    stop("rho_j_BRR: for each i in 1:m, length(data_split[[i]]$y) and nrow(data_split[[i]]$X) must be equal")
  } else if (!is.vector(prior_means) | length(prior_means)!=dim) {
    stop("rho_j_BRR: prior_means must be vectors of length dim")
  } else if (!is.vector(prior_variances) | length(prior_variances)!=dim) {
    stop("rho_j_BRR: prior_variances must be vectors of length dim")
  } else if (!is.list(precondition_matrices) | (length(precondition_matrices)!=m)) {
    stop("rho_j_BRR: precondition_matrices must be a list of length m")
  } else if (!is.list(inv_precondition_matrices) | (length(inv_precondition_matrices)!=m)) {
    stop("rho_j_BRR: inv_precondition_matrices must be a list of length m")
  } else if (!(diffusion_estimator %in% c('Poisson', 'NB'))) {
    stop("rho_j_BRR: diffusion_estimator must be set to either \'Poisson\' or \'NB\'")
  } else if (!any(class(cl)=="cluster") & !is.null(cl)) {
    stop("rho_j_BRR: cl must be a \"cluster\" object or NULL")
  }
  if (adaptive_mesh) {
    if (!is.matrix(sub_posterior_means)) {
      stop("rho_j_BRR: if adaptive_mesh==TRUE, sub_posterior_means must be a (m x dim) matrix")
    } else if (any(dim(sub_posterior_means)!=c(m,dim))) {
      stop("rho_j_BRR: if adaptive_mesh==TRUE, sub_posterior_means must be a (m x dim) matrix")
    }
  }
  if (cv_location == 'mode') {
    cv_location <- lapply(1:m, function(c) {
      MLE <- obtain_LR_MLE(dim = dim, data = data_split[[c]])
      list('beta_hat' = MLE,
           'grad_log_hat' = log_BRR_gradient(beta = MLE,
                                             y_resp = data_split[[c]]$y,
                                             X = data_split[[c]]$X,
                                             X_beta = as.vector(X %*% MLE),
                                             nu = nu,
                                             sigma = sigma,
                                             prior_means = prior_means,
                                             prior_variances = prior_variances,
                                             C = C))})
  } else if (cv_location == 'hypercube_centre') {
    cv_location <- lapply(1:m, function(c) 'hypercube_centre')
  } else {
    stop("rho_j_BRR: cv_location must be either \"mode\" or \"hypercube_centre\"")
  }
  transform_matrices <- lapply(1:m, function(c) {
    list('to_Z' = expm::sqrtm(inv_precondition_matrices[[c]]),
         'to_X' = expm::sqrtm(precondition_matrices[[c]]))
  })
  transformed_design_matrices <- lapply(1:m, function(c) data_split[[c]]$X %*% transform_matrices[[c]]$to_X)
  N <- particle_set$N
  # ---------- creating parallel cluster
  if (is.null(cl)) {
    cl <- parallel::makeCluster(n_cores, setup_strategy = "sequential", outfile = "GBF_BRR_outfile.txt")
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
  elapsed_time <- rep(NA, length(time_mesh)-1)
  ESS <- c(particle_set$ESS[1], rep(NA, length(time_mesh)-1))
  CESS <- c(particle_set$CESS[1], rep(NA, length(time_mesh)-1))
  resampled <- rep(FALSE, length(time_mesh))
  if (record) {
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
  }
  if (is.null(print_progress_iters)) {
    print_progress_iters <- split_N
  }
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
      if (record) {
        E_nu_j[j] <- tilde_Delta_j$E_nu_j
        chosen[j] <- tilde_Delta_j$chosen
        mesh_terms[[j]] <- c(tilde_Delta_j$T1, tilde_Delta_j$T2)
        k4_choice[j] <- tilde_Delta_j$k4_choice  
      }
      time_mesh[j] <- min(end_time, time_mesh[j-1]+tilde_Delta_j$max_delta_j)
    }
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
      if (core == 1) {
        cat('##### t_{j-1}:', time_mesh[j-1], '|| t_{j}:', time_mesh[j], '|| T:', 
            end_time, '#####\n', file = 'rho_j_BRR_progress.txt', append = T)  
      }
      cat('Level:', level, '|| Step:', j, '/', length(time_mesh), '|| Node:', node,
          '|| Core:', core, '|| START \n', file = 'rho_j_BRR_progress.txt', append = T)
      for (i in 1:split_N) {
        x_mean_j[i,] <- weighted_mean_multivariate(matrix = x_j[[i]],
                                                   weights = inv_precondition_matrices,
                                                   inverse_sum_weights = Lambda)
        phi <- lapply(1:m, function(c) {
          ea_BRR_DL_PT(dim = dim,
                       x0 = as.vector(split_x_samples[[core]][[i]][c,]),
                       y = as.vector(x_j[[i]][c,]),
                       s = time_mesh[j-1],
                       t = time_mesh[j],
                       data = data_split[[c]],
                       transformed_design_mat = transformed_design_matrices[[c]],
                       nu = nu,
                       sigma = sigma,
                       prior_means = prior_means,
                       prior_variances = prior_variances,
                       C = C,
                       precondition_mat = precondition_matrices[[c]],
                       transform_mats = transform_matrices[[c]],
                       cv_location = cv_location[[c]],
                       diffusion_estimator = diffusion_estimator,
                       beta_NB = beta_NB,
                       gamma_NB_n_points = gamma_NB_n_points,
                       logarithm = TRUE)})
        log_rho_j[i] <- sum(sapply(1:m, function(c) phi[[c]]$phi))
        if (i%%print_progress_iters==0) {
          cat('Level:', level, '|| Step:', j, '/', length(time_mesh),
              '|| Node:', node, '|| Core:', core, '||', i, '/', split_N, '\n',
              file = 'rho_j_BRR_progress.txt', append = T)
        }
      }
      cat('Level:', level, '|| Step:', j, '/', length(time_mesh),
          '|| Node:', node, '|| Core:', core, '|| DONE ||', split_N, '/',
          split_N, '\n', file = 'rho_j_BRR_progress.txt', append = T)
      return(list('x_j' = x_j, 'x_mean_j' = x_mean_j, 'log_rho_j' = log_rho_j))
    })
    # ---------- update particle set
    particle_set$x_samples <- unlist(lapply(1:length(split_indices), function(i) {
      rho_j_weighted_samples[[i]]$x_j}), recursive = FALSE)
    particle_set$x_means <- do.call(rbind, lapply(1:length(split_indices), function(i) {
      rho_j_weighted_samples[[i]]$x_mean_j}))
    log_rho_j <- unlist(lapply(1:length(split_indices), function(i) {
      rho_j_weighted_samples[[i]]$log_rho_j}))
    norm_weights <- particle_ESS(log_weights = particle_set$log_weights + log_rho_j)
    particle_set$log_weights <- norm_weights$log_weights
    particle_set$normalised_weights <- norm_weights$normalised_weights
    particle_set$ESS <- norm_weights$ESS
    ESS[j] <- particle_set$ESS
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
  if (record) {
    E_nu_j <- E_nu_j[1:j]
    chosen <- chosen[1:j]
    mesh_terms <- mesh_terms[1:j]
    k4_choice <- k4_choice[1:j]
  }
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
  if (record) {
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
  } else {
    return(list('particle_set' = particle_set,
                'proposed_samples' = proposed_samples,
                'time' = elapsed_time,
                'ESS' = ESS,
                'CESS' = CESS,
                'resampled' = resampled))
  }
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
#'                   data_split[[c]]$X is the design matrix for the covariates for
#'                   sub-posterior c
#' @param nu degrees of freedom in t-distribution
#' @param sigma scale parameter in t-distribution
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
#' @param adaptive_mesh_parameters list of parameters used for adaptive mesh
#' @param record logical value indicating if variables such as E[nu_j], chosen,
#'               mesh_terms and k4_choice should be recorded at each iteration
#'               and returned (see return variables for this function) - default
#'               is FALSE
#' @param diffusion_estimator choice of unbiased estimator for the Exact Algorithm
#'                            between "Poisson" (default) for Poisson estimator
#'                            and "NB" for Negative Binomial estimator
#' @param beta_NB beta parameter for Negative Binomial estimator (default 10)
#' @param gamma_NB_n_points number of points used in the trapezoidal estimation
#'                          of the integral found in the mean of the negative
#'                          binomial estimator (default is 2)
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
#'   \item{time_mesh}{time_mesh used}
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
#'   \item{combined_data}{combined data for the fusion density}
#' }
#' If record is set to TRUE, additional components will be returned:
#' \describe{
#'   \item{E_nu_j}{approximation of the average variation of the trajectories
#'                 at each time step}
#'   \item{chosen}{which term was chosen if using an adaptive mesh at each time step}
#'   \item{mesh_terms}{the evaluated terms in deciding the mesh at each time step}
#'   \item{k4_choice}{which of the roots of k4 were chosen}
#' }
#' 
#' @export
parallel_GBF_BRR <- function(particles_to_fuse,
                             N,
                             m,
                             time_mesh,
                             dim,
                             data_split,
                             nu,
                             sigma,
                             prior_means,
                             prior_variances,
                             C,
                             precondition_matrices,
                             inv_precondition_matrices = NULL,
                             Lambda = NULL,
                             resampling_method = 'multi',
                             ESS_threshold = 0.5,
                             cv_location = 'hypercube_centre',
                             sub_posterior_means = NULL,
                             adaptive_mesh = FALSE,
                             adaptive_mesh_parameters = NULL,
                             record = FALSE,
                             diffusion_estimator = 'Poisson',
                             beta_NB = 10,
                             gamma_NB_n_points = 2,
                             seed = NULL,
                             n_cores = parallel::detectCores(),
                             cl = NULL,
                             level = 1,
                             node = 1,
                             print_progress_iters = 1000) {
  if (!is.list(particles_to_fuse) | (length(particles_to_fuse)!=m)) {
    stop("parallel_generalised_BF_SMC_BRR: particles_to_fuse must be a list of length m")
  } else if (!all(sapply(particles_to_fuse, function(sub_posterior) ("particle" %in% class(sub_posterior))))) {
    stop("parallel_generalised_BF_SMC_BRR: particles in particles_to_fuse must be \"particle\" objects")
  } else if (!all(sapply(particles_to_fuse, function(sub_posterior) is.matrix(sub_posterior$y_samples)))) {
    stop("parallel_generalised_BF_SMC_BRR: the particles' samples for y should all be matrices")
  } else if (!all(sapply(particles_to_fuse, function(sub_posterior) ncol(sub_posterior$y_samples)==dim))) {
    stop("parallel_generalised_BF_SMC_BRR: the particles' samples for y should all be matrices with dim columns")
  } else if (!is.vector(time_mesh)) {
    stop("parallel_generalised_BF_SMC_BRR: time_mesh must be an ordered vector of length >= 2")
  } else if (length(time_mesh) < 2) {
    stop("parallel_generalised_BF_SMC_BRR: time_mesh must be an ordered vector of length >= 2")
  } else if (!identical(time_mesh, sort(time_mesh))) {
    stop("parallel_generalised_BF_SMC_BRR: time_mesh must be an ordered vector of length >= 2")
  } else if (!is.list(data_split) | length(data_split)!=m) {
    stop("parallel_generalised_BF_SMC_BRR: data_split must be a list of length m")
  } else if (!all(sapply(1:m, function(i) is.vector(data_split[[i]]$y)))) {
    stop("parallel_generalised_BF_SMC_BRR: for each i in 1:m, data_split[[i]]$y must be a vector")
  } else if (!all(sapply(1:m, function(i) is.matrix(data_split[[i]]$X)))) {
    stop("parallel_generalised_BF_SMC_BRR: for each i in 1:m, data_split[[i]]$X must be a matrix")
  } else if (!all(sapply(1:m, function(i) ncol(data_split[[i]]$X)==dim))) {
    stop("parallel_generalised_BF_SMC_BRR: for each i in 1:m, data_split[[i]]$X must be a matrix with dim columns")
  } else if (!all(sapply(1:m, function(i) length(data_split[[i]]$y)==nrow(data_split[[i]]$X)))) {
    stop("parallel_generalised_BF_SMC_BRR: for each i in 1:m, length(data_split[[i]]$y) and nrow(data_split[[i]]$X) must be equal")
  } else if (!is.vector(prior_means) | length(prior_means)!=dim) {
    stop("parallel_generalised_BF_SMC_BRR: prior_means must be vectors of length dim")
  } else if (!is.vector(prior_variances) | length(prior_variances)!=dim) {
    stop("parallel_generalised_BF_SMC_BRR: prior_variances must be vectors of length dim")
  } else if (!is.list(precondition_matrices) | (length(precondition_matrices)!=m)) {
    stop("parallel_generalised_BF_SMC_BRR: precondition_matrices must be a list of length m")
  } else if (!(diffusion_estimator %in% c('Poisson', 'NB'))) {
    stop("parallel_generalised_BF_SMC_BRR: diffusion_estimator must be set to either \'Poisson\' or \'NB\'")
  } else if ((ESS_threshold < 0) | (ESS_threshold > 1)) {
    stop("parallel_generalised_BF_SMC_BRR: ESS_threshold must be between 0 and 1")
  } else if ((cv_location != 'mode') & (cv_location != 'hypercube_centre')) {
    stop("parallel_generalised_BF_SMC_BRR: cv_location must be either \"mode\" or \"hypercube_centre\"")
  } else if (!any(class(cl)=="cluster") & !is.null(cl)) {
    stop("parallel_generalised_BF_SMC_BRR: cl must be a \"cluster\" object or NULL")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  pcm <- proc.time()
  # ---------- creating parallel cluster
  if (is.null(cl)) {
    cl <- parallel::makeCluster(n_cores, setup_strategy = "sequential", outfile = "GBF_BRR_outfile.txt")
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
  if (is.null(inv_precondition_matrices)) {
    inv_precondition_matrices <- lapply(precondition_matrices, solve)  
  }
  if (is.null(Lambda)) {
    Lambda <- inverse_sum_matrices(inv_precondition_matrices)  
  }
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
  rho_j <- rho_j_BRR(particle_set = particles,
                     m = m,
                     time_mesh = time_mesh,
                     dim = dim,
                     data_split = data_split,
                     nu = nu,
                     sigma = sigma,
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
                     record = record,
                     diffusion_estimator = diffusion_estimator,
                     beta_NB = beta_NB,
                     gamma_NB_n_points = gamma_NB_n_points,
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
    new_sub_posterior_means <- list(t(weighted_mean_multivariate(matrix = sub_posterior_means,
                                                                 weights = inv_precondition_matrices,
                                                                 inverse_sum_weights = Lambda)),
                                    sub_posterior_means)
  } else {
    new_sub_posterior_means <- list(NULL, sub_posterior_means)
  }
  if (record) {
    return(list('particles' = rho_j$particle_set,
                'proposed_samples' = rho_j$proposed_samples,
                'time' = (proc.time()-pcm)['elapsed'],
                'elapsed_time' = c(elapsed_time_rho_0, rho_j$time),
                'time_mesh' = rho_j$particle_set$time_mesh,
                'ESS' = rho_j$ESS,
                'CESS' = rho_j$CESS,
                'resampled' = rho_j$resampled,
                'E_nu_j' = rho_j$E_nu_j,
                'chosen' = rho_j$chosen,
                'mesh_terms' = rho_j$mesh_terms,
                'k4_choice' = rho_j$k4_choice,
                'precondition_matrices' = new_precondition_matrices,
                'sub_posterior_means' = new_sub_posterior_means,
                'combined_data' = combine_data(list_of_data = data_split, dim = dim)))
  } else {
    return(list('particles' = rho_j$particle_set,
                'proposed_samples' = rho_j$proposed_samples,
                'time' = (proc.time()-pcm)['elapsed'],
                'elapsed_time' = c(elapsed_time_rho_0, rho_j$time),
                'time_mesh' = rho_j$particle_set$time_mesh,
                'ESS' = rho_j$ESS,
                'CESS' = rho_j$CESS,
                'resampled' = rho_j$resampled,
                'precondition_matrices' = new_precondition_matrices,
                'sub_posterior_means' = new_sub_posterior_means,
                'combined_data' = combine_data(list_of_data = data_split, dim = dim)))
  }
}

#' (Balanced Binary) D&C Monte Carlo Fusion using SMC
#' 
#' (Balanced Binary) D&C Monte Carlo Fusion using SMC for Bayesian Logistic Regression
#'
#' @param N_schedule vector of length (L-1), where N_schedule[l] is the
#'                   number of samples per node at level l
#' @param m_schedule vector of length (L-1), where m_schedule[k] is the number
#'                   of samples to fuse for level k
#' @param time_mesh time mesh used in Bayesian Fusion. This can either be a vector
#'                  which will be used for each node in the tree, or it can be
#'                  passed in as NULL, where a recommended mesh will be generated
#'                  using the parameters passed into mesh_parameters
#' @param base_samples list of length C, where base_samples[[c]] contains
#'                     the samples for the c-th node in the level
#' @param L total number of levels in the hierarchy
#' @param dim dimension of the predictors (= p+1)
#' @param data_split list of length m where each item is a list of length 4 where
#'                   for c=1,...,m, data_split[[c]]$y is the vector for y responses and
#'                   data_split[[c]]$X is the design matrix for the covariates for
#'                   sub-posterior c
#' @param nu degrees of freedom in t-distribution
#' @param sigma scale parameter in t-distribution
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
#' @param adaptive_mesh logical value to indicate if an adaptive mesh is used
#'                      (default is FALSE)
#' @param adaptive_mesh_parameters list of parameters used for adaptive mesh
#' @param record logical value indicating if variables such as E[nu_j], chosen,
#'               mesh_terms and k4_choice should be recorded at each iteration
#'               and returned (see return variables for this function) - default
#'               is FALSE
#' @param diffusion_estimator choice of unbiased estimator for the Exact Algorithm
#'                            between "Poisson" (default) for Poisson estimator
#'                            and "NB" for Negative Binomial estimator
#' @param beta_NB beta parameter for Negative Binomial estimator (default 10)
#' @param gamma_NB_n_points number of points used in the trapezoidal estimation
#'                          of the integral found in the mean of the negative
#'                          binomial estimator (default is 2)
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
#'   \item{elapsed_time}{list of length (L-1), where elapsed_time[[l]][[i]]
#'                       is the elapsed time of each step of the algorithm for
#'                       level l, node i}
#'   \item{time_mesh}{list of length (L-1), where time_mesh[[l]][[i]]
#'                    is the time_mesh used for level l, node i}
#'   \item{ESS}{list of length (L-1), where ESS[[l]][[i]] is the effective
#'              sample size of the particles after each step BEFORE deciding
#'              whether or not to resample for level l, node i}
#'   \item{CESS}{list of length (L-1), where ESS[[l]][[i]] is the conditional
#'               effective sample size of the particles after each step}
#'   \item{resampled}{list of length (L-1), where resampled[[l]][[i]] is a
#'                    boolean value to record if the particles were resampled
#'                    after each step; rho and Q for level l, node i}
#'   \item{precondition_matrices}{pre-conditioning matrices that were used}
#'   \item{sub_posterior_means}{sub-posterior means that were used}
#'   \item{recommended_mesh}{list of length (L-1), where recommended_mesh[[l]][[i]]
#'                           is the recommended mesh for level l, node i}
#'   \item{data_inputs}{list of length (L-1), where data_inputs[[l]][[i]] is the
#'                      data input for the sub-posterior in level l, node i}
#' }
#' If record is set to TRUE, additional components will be returned:
#' \describe{
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
#' }
#'
#' @export
bal_binary_GBF_BRR <- function(N_schedule,
                               m_schedule,
                               time_mesh = NULL,
                               base_samples,
                               L,
                               dim,
                               data_split,
                               nu,
                               sigma,
                               prior_means,
                               prior_variances,
                               C,
                               precondition = TRUE,
                               resampling_method = 'multi',
                               ESS_threshold = 0.5,
                               cv_location = 'hypercube_centre',
                               adaptive_mesh = FALSE,
                               mesh_parameters = NULL,
                               record = FALSE,
                               diffusion_estimator = 'Poisson',
                               beta_NB = 10,
                               gamma_NB_n_points = 2,
                               seed = NULL,
                               n_cores = parallel::detectCores(),
                               print_progress_iters = 1000) {
  if (!is.vector(N_schedule) | (length(N_schedule)!=(L-1))) {
    stop("bal_binary_GBF_BRR: N_schedule must be a vector of length (L-1)")
  } else if (!is.vector(m_schedule) | (length(m_schedule)!=(L-1))) {
    stop("bal_binary_GBF_BRR: m_schedule must be a vector of length (L-1)")
  } else if (!is.list(base_samples) | (length(base_samples)!=C)) {
    stop("bal_binary_GBF_BRR: base_samples must be a list of length C")
  } else if (!is.list(data_split) | length(data_split)!=C) {
    stop("bal_binary_GBF_BRR: data_split must be a list of length C")
  } else if (!all(sapply(1:C, function(i) is.vector(data_split[[i]]$y)))) {
    stop("bal_binary_GBF_BRR: for each i in 1:C, data_split[[i]]$y must be a vector")
  } else if (!all(sapply(1:C, function(i) is.matrix(data_split[[i]]$X)))) {
    stop("bal_binary_GBF_BRR: for each i in 1:C, data_split[[i]]$X must be a matrix")
  } else if (!all(sapply(1:C, function(i) ncol(data_split[[i]]$X)==dim))) {
    stop("bal_binary_GBF_BRR: for each i in 1:C, data_split[[i]]$X must be a matrix with dim columns")
  } else if (!all(sapply(1:C, function(i) length(data_split[[i]]$y)==nrow(data_split[[i]]$X)))) {
    stop("bal_binary_GBF_BRR: for each i in 1:C, length(data_split[[i]]$y) and nrow(data_split[[i]]$X) must be equal")
  } else if (!is.vector(prior_means) | length(prior_means)!=dim) {
    stop("bal_binary_GBF_BRR: prior_means must be vectors of length dim")
  } else if (!is.vector(prior_variances) | length(prior_variances)!=dim) {
    stop("bal_binary_GBF_BRR: prior_variances must be vectors of length dim")
  } else if (ESS_threshold < 0 | ESS_threshold > 1) {
    stop("bal_binary_GBF_BRR: ESS_threshold must be between 0 and 1")
  }
  if (is.vector(m_schedule) & (length(m_schedule)==(L-1))) {
    for (l in (L-1):1) {
      if ((C/prod(m_schedule[(L-1):l]))%%1!=0) {
        stop("bal_binary_GBF_BRR: check that C/prod(m_schedule[(L-1):l])
              is an integer for l=L-1,...,1")
      }
    }
  } else {
    stop("bal_binary_GBF_BRR: m_schedule must be a vector of length (L-1)")
  }
  if (is.vector(time_mesh)) {
    if (length(time_mesh) < 2) {
      stop("bal_binary_GBF_BRR: time_mesh must be an ordered vector of length >= 2")
    } else if (!identical(time_mesh, sort(time_mesh))) {
      stop("bal_binary_GBF_BRR: time_mesh must be an ordered vector of length >= 2")
    }
  } else if (is.null(time_mesh)) {
    if (!is.list(mesh_parameters)) {
      stop("bal_binary_GBF_BRR: if time_mesh is NULL, mesh_parameters must be a
           list of parameters to obtain guidance for the mesh")
    } 
  } else {
    stop("bal_binary_GBF_BRR: time_mesh must either be an ordered vector of length
           >= 2 or passed as NULL if want to use recommended guidance")
  }
  m_schedule <- c(m_schedule, 1)
  particles <- list()
  if (all(sapply(base_samples, function(sub) class(sub)=='particle'))) {
    particles[[L]] <- base_samples
  } else if (all(sapply(base_samples, is.matrix))) {
    if (!all(sapply(base_samples, function(core) ncol(core)==dim))) {
      stop("bal_binary_GBF_BRR: the sub-posterior samples in base_samples must be matrices with dim columns")
    }
    particles[[L]] <- initialise_particle_sets(samples_to_fuse = base_samples,
                                               multivariate = TRUE,
                                               number_of_steps = 2)
  } else {
    stop("bal_binary_GBF_BRR: base_samples must be a list of length C
         containing either items of class \"particle\" (representing particle 
         approximations of the sub-posteriors) or are matrices with dim columns
         (representing un-normalised sample approximations of the sub-posteriors)")
  }
  proposed_samples <- list()
  data_inputs <- list()
  data_inputs[[L]] <- data_split
  time <- list()
  elapsed_time <- list()
  used_time_mesh <- list()
  ESS <- list()
  CESS <- list()
  resampled <- list()
  if (record) {
    E_nu_j <- list()
    chosen <- list()
    mesh_terms <- list()
    k4_choice <- list()  
  }
  recommended_mesh <- list()
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
    stop("bal_binary_GBF_BRR: precondition must be a logical indicating 
          whether or not a preconditioning matrix should be used, or a list of
          length C, where precondition[[c]] is the preconditioning matrix for
          the c-th sub-posterior")
  }
  sub_posterior_means <- list()
  sub_posterior_means[[L]] <- t(sapply(base_samples, function(sub) apply(sub, 2, mean)))
  cl <- parallel::makeCluster(n_cores, setup_strategy = "sequential", outfile = "SMC_BRR_outfile.txt")
  parallel::clusterExport(cl, envir = environment(), varlist = ls())
  parallel::clusterExport(cl, varlist = ls("package:DCFusion"))
  parallel::clusterExport(cl, varlist = ls("package:layeredBB"))
  cat('Starting bal_binary fusion \n', file = 'bal_binary_GBF_BRR.txt')
  for (k in ((L-1):1)) {
    n_nodes <- max(C/prod(m_schedule[L:k]), 1)
    cat('########################\n', file = 'bal_binary_GBF_BRR.txt', append = T)
    cat('Starting to fuse', m_schedule[k], 'sub-posteriors for level', k, 
        'using', n_cores, 'cores\n', file = 'bal_binary_GBF_BRR.txt', append = T)
    cat('At this level, the data is split up into', (C/prod(m_schedule[L:(k+1)])), 'subsets\n',
        file = 'bal_binary_GBF_BRR.txt', append = T)
    cat('There are', n_nodes, 'nodes at the next level each giving', N_schedule[k],
        'samples \n', file = 'bal_binary_GBF_BRR.txt', append = T)
    cat('########################\n', file = 'bal_binary_GBF_BRR.txt', append = T)
    fused <- lapply(X = 1:n_nodes, FUN = function(i) {
      previous_nodes <- ((m_schedule[k]*i)-(m_schedule[k]-1)):(m_schedule[k]*i)
      particles_to_fuse <- particles[[k+1]][previous_nodes]
      precondition_mats <- precondition_matrices[[k+1]][previous_nodes]
      inv_precondition_mats <- lapply(precondition_mats, solve)
      Lambda <- inverse_sum_matrices(inv_precondition_mats)
      sub_post_means <- sub_posterior_means[[k+1]][previous_nodes,]
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
                                      data_size = length(data_inputs[[k+1]][[1]]$y),
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
                                      trial_k3_by = mesh_parameters$trial_k3_by,
                                      vanilla = mesh_parameters$vanilla)
      } else {
        recommendation <- list('mesh' = time_mesh)
      }
      return(list('recommendation' = recommendation,
                  'fusion' = parallel_GBF_BRR(particles_to_fuse = particles_to_fuse,
                                              N = N_schedule[k],
                                              m = m_schedule[k],
                                              time_mesh = recommendation$mesh,
                                              dim = dim,
                                              data_split = data_inputs[[k+1]][previous_nodes],
                                              nu = nu,
                                              sigma = sigma,
                                              prior_means = prior_means,
                                              prior_variances = prior_variances,
                                              C = (C/prod(m_schedule[L:(k+1)])),
                                              precondition_matrices = precondition_mats,
                                              inv_precondition_matrices = inv_precondition_mats,
                                              Lambda = Lambda,
                                              resampling_method = resampling_method,
                                              ESS_threshold = ESS_threshold,
                                              cv_location = cv_location,
                                              sub_posterior_means = sub_post_means,
                                              adaptive_mesh = adaptive_mesh,
                                              adaptive_mesh_parameters = mesh_parameters,
                                              record = record,
                                              diffusion_estimator = diffusion_estimator,
                                              beta_NB = beta_NB,
                                              gamma_NB_n_points = gamma_NB_n_points,
                                              seed = seed,
                                              n_cores = n_cores,
                                              cl = cl,
                                              level = k,
                                              node = i,
                                              print_progress_iters = print_progress_iters)))
    })
    particles[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$fusion$particles)
    proposed_samples[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$fusion$proposed_samples)
    time[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$fusion$time)
    elapsed_time[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$fusion$elapsed_time)
    used_time_mesh[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$fusion$time_mesh)
    ESS[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$fusion$ESS)
    CESS[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$fusion$CESS)
    resampled[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$fusion$resampled)
    if (record) {
      E_nu_j[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$fusion$E_nu_j)
      chosen[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$fusion$chosen)
      mesh_terms[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$fusion$mesh_terms)
      k4_choice[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$fusion$k4_choice)  
    }
    precondition_matrices[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$fusion$precondition_matrices[[1]])
    sub_posterior_means[[k]] <- do.call(rbind, lapply(1:n_nodes, function(i) fused[[i]]$fusion$sub_posterior_means[[1]]))
    data_inputs[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$fusion$combined_data)
    recommended_mesh[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$recommendation)
  }
  parallel::stopCluster(cl)
  cat('Completed bal_binary fusion\n', file = 'bal_binary_GBF_BRR.txt', append = T)
  if (length(particles[[1]])==1) {
    particles[[1]] <- particles[[1]][[1]]
    proposed_samples[[1]] <- proposed_samples[[1]][[1]]
    time[[1]] <- time[[1]][[1]]
    elapsed_time[[1]] <- elapsed_time[[1]][[1]]
    used_time_mesh[[1]] <- used_time_mesh[[1]][[1]]
    ESS[[1]] <- ESS[[1]][[1]]
    CESS[[1]] <- CESS[[1]][[1]]
    resampled[[1]] <- resampled[[1]][[1]]
    if (record) {
      E_nu_j[[1]] <- E_nu_j[[1]][[1]]
      chosen[[1]] <- chosen[[1]][[1]]
      mesh_terms[[1]] <- mesh_terms[[1]][[1]]
      k4_choice[[1]] <- k4_choice[[1]][[1]]  
    }
    precondition_matrices[[1]] <- precondition_matrices[[1]][[1]]
    sub_posterior_means[[1]] <- sub_posterior_means[[1]][[1]]
    recommended_mesh[[1]] <- recommended_mesh[[1]][[1]]
    data_inputs[[1]] <- data_inputs[[1]][[1]]
  }
  if (record) {
    return(list('particles' = particles,
                'proposed_samples' = proposed_samples,
                'time' = time,
                'elapsed_time' = elapsed_time,
                'time_mesh' = used_time_mesh,
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
                'data_inputs' = data_inputs))    
  } else {
    return(list('particles' = particles,
                'proposed_samples' = proposed_samples,
                'time' = time,
                'elapsed_time' = elapsed_time,
                'time_mesh' = used_time_mesh,
                'ESS' = ESS,
                'CESS' = CESS,
                'resampled' = resampled,
                'precondition_matrices' = precondition_matrices,
                'sub_posterior_means' = sub_posterior_means,
                'recommended_mesh' = recommended_mesh,
                'data_inputs' = data_inputs))
  }
}
