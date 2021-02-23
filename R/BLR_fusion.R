#' @export
ea_phi_BLR_DL_bounds <- function(initial_parameters,
                                 y_labels,
                                 X,
                                 prior_means,
                                 prior_variances,
                                 C,
                                 lower,
                                 upper,
                                 precondition_mat, 
                                 transform_mat,
                                 bounds_multiplier = 1.2) {
  # checking that initial parameters is between lower and upper bounds
  # might not be if we use a transformation to z space
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
  # obtain lower bound using 'optim' function in R's stats package
  LB <- optim(par = initial_parameters,
              fn = function(beta) ISFusion:::phi_BLR(beta = beta,
                                                     y_labels = y_labels,
                                                     X = X,
                                                     prior_means = prior_means,
                                                     prior_variances = prior_variances,
                                                     C = C,
                                                     precondition_mat = precondition_mat,
                                                     transform_mat = transform_mat),
              method = "L-BFGS-B",
              lower = lower,
              upper = upper,
              control = list('fnscale' = 1, 'maxit' = 250))
  # obtain upper bound using 'optim' function in R's stats package
  UB <- optim(par = initial_parameters,
              fn = function(beta) ISFusion:::phi_BLR(beta = beta,
                                                     y_labels = y_labels,
                                                     X = X,
                                                     prior_means = prior_means,
                                                     prior_variances = prior_variances,
                                                     C = C,
                                                     precondition_mat = precondition_mat,
                                                     transform_mat = transform_mat),
              method = "L-BFGS-B",
              lower = lower,
              upper = upper,
              control = list('fnscale' = -1, 'maxit' = 250))
  # multiply the bounds by bounds_multiplier to compensate the case 
  # if the opitimser did not get the bounds correct
  # calculate the Bounds Difference
  BD <- UB$value - LB$value
  # subtract or add 0.5*(bounds_multiplier-1)*BD
  # makes the resulting bounds difference be bounds_multipler larger
  return(list('LB' = LB$value - 0.5*(bounds_multiplier-1)*BD, 
              'UB' = UB$value + 0.5*(bounds_multiplier-1)*BD))
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
                         inv_precondition_mat, 
                         diffusion_estimator,
                         beta_NB = 10,
                         bounds_multiplier = 1.2, 
                         logarithm) {
  # create variables to store matricies for transforming between spaces
  transform_to_z <- expm::sqrtm(inv_precondition_mat)
  transform_to_x <- expm::sqrtm(precondition_mat)
  # transform start and end point to preconditioned space
  z0 <- as.vector(transform_to_z %*% x0)
  zt <- as.vector(transform_to_z %*% y)
  
  # ----- simulate layer information to bound sample path
  # simulate layer information to bound sample path
  # create sequence vector
  a_seq <- seq(from = 0, to = ceiling(t-s), by = (t-s)/20)
  # simulate layer information (Bessel layer)
  bes_layers_z <- layeredBB::multi_bessel_layer_simulation(dim = dim,
                                                           x = z0,
                                                           y = zt,
                                                           s = s,
                                                           t = t,
                                                           a = a_seq)
  # ----- compute bounds for phi_BLR
  # set lower and upper bounds of the Brownian bridge
  lbound_z <- sapply(X = 1:dim, function(d) min(z0[d], zt[d]) - bes_layers_z[[d]]$a[bes_layers_z[[d]]$l])
  ubound_z <- sapply(X = 1:dim, function(d) max(z0[d], zt[d]) + bes_layers_z[[d]]$a[bes_layers_z[[d]]$l])
  
  # calculate upper and lower bounds of phi given the simulated sample layer information
  # we try different initial parameters and take the smallest and largest bounds obtained
  bounds <- lapply(list(lbound_z, ubound_z), function(init) {
    bound_phi_BLR(initial_parameters = init,
                  y_labels = y_labels,
                  X = X,
                  prior_means = prior_means,
                  prior_variances = prior_variances,
                  C = C,
                  lower = lbound_z,
                  upper = ubound_z, 
                  precondition_mat = precondition_mat,
                  transform_mat = transform_to_x,
                  bounds_multiplier = bounds_multiplier)})
  LZ <- min(sapply(1:length(bounds), function(i) bounds[[i]]$LB))
  UZ <- max(sapply(1:length(bounds), function(i) bounds[[i]]$UB))
  
  if (diffusion_estimator == 'Poisson') {
    # ----- construct Poisson estimator
    # simulate skeletal points for Poisson thinning
    # simulate number of points to simulate from Poisson distribution
    kap <- rpois(1, (UZ-LZ)*(t-s))
    log_acc_prob <- 0
    if (kap > 0) {
      # simulating time points from a Uniform(0, (t-s)) distribution
      # simulating sample path at simulated times
      layered_BB <- layeredBB::multi_layered_brownian_bridge(dim = dim,
                                                             x = z0,
                                                             y = zt,
                                                             s = s,
                                                             t = t,
                                                             layers = bes_layers_z,
                                                             times = runif(kap, s, t))
      # layered_bb includes points for times s and t - but we only want those simulated at simulated times
      pois_points <- as.matrix(layered_BB[1:(nrow(layered_BB)-1), 2:(ncol(layered_BB)-1)])
      phi_eval <- rep(0, ncol(pois_points))
      for (j in 1:ncol(pois_points)) {
        phi_eval[j] <- ISFusion:::phi_BLR(beta = as.vector(pois_points[,j]),
                                          y_labels = y_labels,
                                          X = X,
                                          prior_means = prior_means, 
                                          prior_variances = prior_variances, 
                                          C = C, 
                                          precondition_mat = precondition_mat, 
                                          transform_mat = transform_to_x)
      }
      terms <- (UZ - phi_eval)
      if (any(terms < 0)) {
        cat('LZ:', LZ, '\n', file = 'bounds.txt', append = T)
        cat('UZ:', UZ, '\n', file = 'bounds.txt', append = T)
        cat('phi_eval:', phi_eval, '\n', file = 'bounds.txt', append = T)
        cat('(UZ-phi_eval):', terms, '\n', file = 'bounds.txt', append = T)
        cat('(phi_eval-LZ):', phi_eval-LZ, '\n', file = 'bounds.txt', append = T)
        stop('Some of (UZ-phi_eval) are < 0')
      }
      log_acc_prob <- sum(log(terms))
      # checking lower bounds
      if (any((phi_eval - LZ) < 0)) {
        cat('LZ:', LZ, '\n', file = 'bounds.txt', append = T)
        cat('UZ:', UZ, '\n', file = 'bounds.txt', append = T)
        cat('phi_eval:', phi_eval, '\n', file = 'bounds.txt', append = T)
        cat('(UZ-phi_eval):', terms, '\n', file = 'bounds.txt', append = T)
        cat('(phi_eval-LZ):', phi_eval-LZ, '\n', file = 'bounds.txt', append = T)
        stop('Some of (phi_eval-LZ) are < 0')
      }
    } 
    
    if (logarithm) {
      return(-LZ*(t-s) - kap*log(UZ-LZ) + log_acc_prob)
    } else {
      return(exp(-LZ*(t-s) - kap*log(UZ-LZ) + log_acc_prob))
    }
  } else if (diffusion_estimator == "NB") {
    # ----- construct Negative Binomial estimator
    integral_estimate <- cubature::adaptIntegrate(f = function(s_) {
      phi_BLR(beta = (z0*(t-s_) + zt*s_)/(t-s),
              y_labels = y_labels,
              X = X,
              prior_means = prior_means,
              prior_variances = prior_variances,
              C = C,
              precondition_mat = precondition_mat,
              transform_mat = transform_to_x)},
      lowerLimit = s,
      upperLimit = t)$integral
    gamma_NB <- (t-s)*UZ - integral_estimate
    
    # simulate skeletal points for NB estimator
    # simulate number of points to simulate from NB distribution
    kap <- rnbinom(1, size = beta_NB, mu = gamma_NB)
    log_acc_prob <- 0
    if (kap > 0) {
      # simulating time points from a Uniform(0, (t-s)) distribution
      # simulating sample path at simulated times
      layered_BB <- layeredBB::multi_layered_brownian_bridge(dim = dim,
                                                             x = z0,
                                                             y = zt,
                                                             s = s,
                                                             t = t,
                                                             layers = bes_layers_z,
                                                             times = runif(kap, s, t))
      # layered_bb includes points for times s and t - but we only want those simulated at simulated times
      kappa_points <- as.matrix(layered_BB[1:(nrow(layered_BB)-1), 2:(ncol(layered_BB)-1)])
      phi_eval <- rep(0, ncol(kappa_points))
      for (j in 1:ncol(kappa_points)) {
        phi_eval[j] <- ISFusion:::phi_BLR(beta = as.vector(kappa_points[,j]),
                                          y_labels = y_labels,
                                          X = X,
                                          prior_means = prior_means, 
                                          prior_variances = prior_variances, 
                                          C = C, 
                                          precondition_mat = precondition_mat, 
                                          transform_mat = transform_to_x)
      }
      terms <- (UZ - phi_eval)
      if (any(terms < 0)) {
        cat('LZ:', LZ, '\n', file = 'bounds.txt', append = T)
        cat('UZ:', UZ, '\n', file = 'bounds.txt', append = T)
        cat('phi_eval:', phi_eval, '\n', file = 'bounds.txt', append = T)
        cat('(UZ-phi_eval):', terms, '\n', file = 'bounds.txt', append = T)
        cat('(phi_eval-LZ):', phi_eval-LZ, '\n', file = 'bounds.txt', append = T)
        stop('Some of (UZ-phi) are < 0')
      }
      log_acc_prob <- sum(log(terms))
      # checking lower bounds
      if (any((phi_eval - LZ) < 0)) {
        cat('LZ:', LZ, '\n', file = 'bounds.txt', append = T)
        cat('UZ:', UZ, '\n', file = 'bounds.txt', append = T)
        cat('phi_eval:', phi_eval, '\n', file = 'bounds.txt', append = T)
        cat('(UZ-phi_eval):', terms, '\n', file = 'bounds.txt', append = T)
        cat('(phi_eval-LZ):', phi_eval-LZ, '\n', file = 'bounds.txt', append = T)
        stop('Some of (phi-LZ) are < 0')
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
    stop("log_diffusion_probability_BLR: diffusion_estimator must be set to 
    either \'Poisson\' or \'NB\'")
  }
}

# rho_IS_R <- function(particles_to_fuse,
#                      N,
#                      dim,
#                      time,
#                      m,
#                      inv_precondition_matrices,
#                      sum_inv_precondition_matrices) {
#   # right now in this function we do not update the weights
#   # we set them at the result from assigning the weights with rho
#   x_samples <- rep(list(matrix(nrow = m, ncol = dim)), N)
#   x_means <- matrix(nrow = N, ncol = dim)
#   log_rho_weights <- rep(0, N)
#   for (i in 1:N) {
#     for (c in 1:m) {
#       x_samples[[i]][c,] <- particles_to_fuse[[c]]$y_samples[i,]
#     }
#     x_means[i,] <- weighted_mean_matrix(matrix = x_samples[[i]],
#                                         weights = inv_precondition_matrices,
#                                         inv_weights_sum = sum_inv_precondition_matrices)
#     # could the particles_to_fuse weights to this in order to update weights from previous level
#     log_rho_weights[i] <- log_rho(x = x_samples[[i]],
#                                   x_mean = x_means[i,],
#                                   time = time,
#                                   inv_precondition_matrices = inv_precondition_matrices)
#   }
#   # normalise log weights and calculate the ESS
#   norm_weights <- particle_ESS(log_weights = log_rho_weights)
#   return(list('y_samples' = matrix(nrow = N, ncol = dim),
#               'x_samples' = x_samples, 
#               'x_means' = x_means, 
#               'log_weights' = log_rho_weights,
#               'normalised_weights' = norm_weights$normalised_weights,
#               'ESS' = norm_weights$ESS,
#               'CESS' = norm_weights$ESS,
#               'resampled' = FALSE,
#               'N' = N))
# }


Q_IS_BLR <- function(particle_set,
                     dim,
                     y_split,
                     X_split,
                     prior_means,
                     prior_variances,
                     time,
                     m,
                     C,
                     precondition_matrices,
                     inv_precondition_matrices,
                     diffusion_estimator,
                     beta_NB = 10,
                     bounds_multiplier = 1.2,
                     seed = NULL,
                     level = 1,
                     node = 1,
                     n_cores = parallel::detectCores()) {
  # calculate proposal covariance matrix
  proposal_cov <- calculate_proposal_cov(time = time, weights = inv_precondition_matrices)
  N <- particle_set$N
  # ---------- creating parallel cluster
  cl <- parallel::makeCluster(n_cores, 
                              setup_strategy = "sequential",
                              outfile = "smc_fusion_outfile.txt")
  # creating variable and functions list to pass into cluster using clusterExport()
  varlist <- c(ls(), list("phi_BLR", 
                          "bound_phi_BLR",
                          "diffusion_probability_BLR"))
  parallel::clusterExport(cl, envir = environment(), varlist = varlist)
  # exporting functions from layeredBB package to simulate layered Brownian bridges
  parallel::clusterExport(cl, varlist = ls("package:layeredBB"))
  if (!is.null(seed)) {
    # setting seed for the cores in the cluster
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
    # for each particle (each x), propose a value for y and update weights
    for (i in 1:split_N) {
      # simulate proposed value y from a Gaussian distribution
      y_samples[i,] <- ISFusion:::mvrnormArma(N = 1,
                                              mu = split_x_means[[core]][i,],
                                              Sigma = proposal_cov)
      # simulate m diffusion and obtain weight for particle set
      for (c in 1:m) {
        log_Q <- diffusion_probability_BLR(dim = dim,
                                           x0 = as.vector(split_x_samples[[core]][[i]][c,]),
                                           y = as.vector(y_samples[i,]),
                                           s = 0,
                                           t = time,
                                           y_labels = y_split[[c]],
                                           X = X_split[[c]],
                                           prior_means = prior_means,
                                           prior_variances = prior_variances,
                                           C = C,
                                           precondition_mat = precondition_matrices[[c]],
                                           inv_precondition_mat = inv_precondition_matrices[[c]],
                                           diffusion_estimator = diffusion_estimator, 
                                           beta_NB = beta_NB, 
                                           bounds_multiplier = bounds_multiplier,
                                           logarithm = TRUE)
        # update particle weights
        log_Q_weights[i] <- log_Q_weights[i] + log_Q
      }
      cat('Level:', level, '|| Node:', node, '|| Core:', core, '||', i, '/', 
          split_N, '\n', file = 'Q_IS_progress.txt', append = T)
    }
    return(list('y_samples' = y_samples, 'log_Q_weights' = log_Q_weights))
  })
  # stopping cluster
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
  particle_set$CESS <- particle_ESS(log_weights = log_Q_weights)$ESS
  # set the resampled indicator to FALSE
  particle_set$resampled <- FALSE
  return(particle_set)
}

# -------------------- parallelised SMC fusion -------------------- #

#' Parallel Sequential Monte Carlo Fusion for Bayesian Logistic Regression model

#' @param particles_to_fuse list of length m, where particles_to_fuse[c] contains 
#'                          the particles for the c-th sub-posterior. Can 
#'                          initialise a this from list of sub-posterior samples 
#'                          by using the intialise_particle_sets function
#' @param N number of samples
#' @param dim dimension of the predictors (= p+1)
#' @param y_split list of length m, where y_split[[c]] is the y responses for
#'                sub-posterior c
#' @param X_split list of length m, where X_split[[c]] is the design matrix for
#'                sub-posterior c
#' @param prior_means prior for means of predictors
#' @param prior_variances prior for variances of predictors
#' @param time time T for fusion algorithm
#' @param m number of sub-posteriors to combine
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
#' @param seed seed number - default is NULL, meaning there is no seed
#' @param level indicates which level this is for the hierarchy (default 1)
#' @param node indicates which node this is for the hierarchy (default 1)
#' @param n_cores number of cores to use
#'
#' @return A list with components:
#' \describe{
#'   \item{particles}{weighted particle set from fusion}
#'   \item{proposed_samples}{proposal samples for y}
#'   \item{ESS}{effective sample size after rho and Q step}
#'   \item{resampled}{boolean value to record if the particles were resampled
#'                    after each step; rho and Q}
#'   \item{time}{user time elapsed for algorithm to run}
#'   \item{precondition_matrices}{list of length 2 where precondition_matrices[[2]] 
#'                                are the pre-conditioning matrices that were used 
#'                                and precondition_matrices[[1]] are the combined 
#'                                precondition matrices}
#'   \item{combined_y}{combined y responses after fusion}
#'   \item{combined_X}{combined design matrix after fusion}
#' }
#'
#' @export
parallel_fusion_SMC_BLR <- function(particles_to_fuse,
                                    N, 
                                    dim,
                                    y_split,
                                    X_split,
                                    prior_means,
                                    prior_variances,
                                    time,
                                    m,
                                    C,
                                    precondition_matrices, 
                                    resampling_method = 'multi',
                                    ESS_threshold = 0.5,
                                    diffusion_estimator = 'Poisson',
                                    beta_NB = 10,
                                    bounds_multiplier = 1.2,
                                    seed = NULL,
                                    level = 1,
                                    node = 1,
                                    n_cores = parallel::detectCores()) {
  if (length(particles_to_fuse)!=m) {
    stop("parallel_fusion_SMC_BLR: particles_to_fuse must be a list of length m")
  } else if (length(y_split)!=m) {
    stop("parallel_fusion_SMC_BLR: y_split must be a list of length m")
  } else if (length(X_split)!=m) {
    stop("parallel_fusion_SMC_BLR: X_split must be a list of length m")
  } else if (length(precondition_matrices)!=m) {
    stop("parallel_fusion_SMC_BLR: precondtion_matrices must be a list of length m")
  } else if (diffusion_estimator != 'Poisson' & diffusion_estimator != 'NB') {
    stop("parallel_fusion_SMC_BLR: diffusion estimator must be either \'Poisson\' or \'NB\'")
  } 
  # check that all the samples in particles_to_fuse are matrices with dim columns
  for (c in 1:length(particles_to_fuse)) {
    if (ncol(particles_to_fuse[[c]]$y_samples)!=dim) {
      stop("parallel_fusion_SMC_BLR: check that particles_to_fuse contains 
           matrices with dim columns")
    }
  }
  # set a seed if one is supplied
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # ---------- resample the particles if they do not have equal weights
  # check if the resampled indicator if FALSE
  # ---------- also check if there are enough samples
  for (c in 1:length(particles_to_fuse)) {
    if (!particles_to_fuse[[c]]$resampled) {
      indices <- resample_indices(normalised_weights = particles_to_fuse[[c]]$normalised_weights,
                                  method = resampling_method,
                                  n = N)
      particles_to_fuse[[c]]$y_samples <- particles_to_fuse[[c]]$y_samples[indices,]
      # reset weights and ESS
      particles_to_fuse[[c]]$normalised_weights[] <- 1/N
      particles_to_fuse[[c]]$log_weights[] <- log(1/N)
      particles_to_fuse[[c]]$ESS <- N
    } else if (particles_to_fuse[[c]]$N < N) {
      indices <- resample_indices(normalised_weights = particles_to_fuse[[c]]$normalised_weights,
                                  method = resampling_method,
                                  n = N)
      particles_to_fuse[[c]]$y_samples <- particles_to_fuse[[c]]$y_samples[indices,]
      # reset weights and ESS
      particles_to_fuse[[c]]$normalised_weights[] <- 1/N
      particles_to_fuse[[c]]$log_weights[] <- log(1/N)
      particles_to_fuse[[c]]$ESS <- N
    }
  }
  # ---------- first importance sampling step 
  # pre-calculating the inverse precondition matrices
  inv_precondition_matrices <- lapply(precondition_matrices, solve)
  # importance samping for rho step
  particles <- rho_IS_R(particles_to_fuse = particles_to_fuse,
                        N = N,
                        dim = dim,
                        time = time,
                        m = m,
                        inv_precondition_matrices = inv_precondition_matrices, 
                        sum_inv_precondition_matrices = inv_sum_matrices(inv_precondition_matrices))
  # record ESS and CESS after rho step 
  ESS <- c('rho' = particles$ESS)
  CESS <- c('rho' = particles$CESS)
  # ----------- resample particles
  # only resample if ESS < N*ESS_threshold
  if (particles$ESS < N*ESS_threshold) {
    particles$resampled <- TRUE
    resampled <- c('rho' = TRUE)
    indices <- resample_indices(normalised_weights = particles$normalised_weights,
                                method = resampling_method,
                                n = N)
    particles$x_samples <- particles$x_samples[indices]
    particles$x_means <- particles$x_means[indices,]
    # reset weights and ESS
    particles$normalised_weights[] <- 1/N
    particles$log_weights[] <- log(1/N)
    particles$ESS <- N
  } else {
    particles$resampled <- FALSE
    resampled <- c('rho' = FALSE)
  }
  # ---------- second importance sampling step
  # # need to find the lower bound K for the Exact Algorithm (second importance 
  # # sampling step for Q) - cancels out but calculate anyway to pass 
  # # into diffusion_probability_biGaussian()
  # K <- sapply(1:m, function(c) {
  #   BayesLogitFusion:::phi_LB_BLR(X = X_split[[c]],
  #                                 prior_variances = prior_variances,
  #                                 C = C,
  #                                 precondition_mat = precondition_matrices[[c]])
  # })
  # unbiased estimator for Q
  particles <- Q_IS_BLR(particle_set = particles, 
                        dim = dim,
                        y_split = y_split,
                        X_split = X_split, 
                        prior_means = prior_means,
                        prior_variances = prior_variances, 
                        time = time,
                        m = m, 
                        C = C,
                        precondition_matrices = precondition_matrices,
                        inv_precondition_matrices = inv_precondition_matrices,
                        diffusion_estimator = diffusion_estimator,
                        beta_NB = beta_NB,
                        bounds_multiplier = bounds_multiplier,
                        seed = seed,
                        level = level,
                        node = node,
                        n_cores = n_cores)
  # record ESS and CESS after Q step
  ESS['Q'] <- particles$ESS
  CESS['Q'] <- particles$CESS
  # record proposed samples
  proposed_samples <- particles$y_samples
  # ---------- resample particles to return an equally weighted particle set 
  # only resample if ESS < N*ESS_threshold
  if (particles$ESS < N*ESS_threshold) {
    particles$resampled <- TRUE
    resampled['Q'] <- TRUE
    indices <- resample_indices(normalised_weights = particles$normalised_weights,
                                method = resampling_method,
                                n = N)
    particles$y_samples <- particles$y_samples[indices,]
    # reset weights and ESS
    particles$normalised_weights[] <- 1/N
    particles$log_weights[] <- log(1/N)
    particles$ESS <- N
  } else {
    particles$resampled <- FALSE
    resampled['Q'] <- FALSE
  }
  if (identical(precondition_matrices, rep(list(diag(1, dim)), m))) {
    return(list('particles' = particles,
                'proposed_samples' = proposed_samples,
                'ESS' = ESS,
                'CESS' = CESS,
                'resampled' = resampled,
                'precondition_matrices' = list(diag(1, dim), 
                                               precondition_matrices),
                'combined_y' = do.call(append, y_split),
                'combined_X' = do.call(rbind, X_split)))
  } else {
    return(list('particles' = particles,
                'proposed_samples' = proposed_samples,
                'ESS' = ESS,
                'CESS' = CESS,
                'resampled' = resampled,
                'precondition_matrices' = list(inv_sum_matrices(inv_precondition_matrices),
                                               precondition_matrices),
                'combined_y' = do.call(append, y_split),
                'combined_X' = do.call(rbind, X_split)))
  }
}

# -------------------- hierarchical SMC fusion -------------------- #

#' Hierarchical Monte Carlo Fusion using SMC for Bayesian Logistic Regression model
#'
#' @param N_schedule vector of length (L-1), where N_schedule[l] is the 
#'                   number of samples per node at level l
#' @param dim dimension of the predictors (= p+1)
#' @param y_split list of length C, where y_split[[c]] is the y responses 
#'                for sub-posterior c
#' @param X_split list of length C, where X_split[[c]] is the design matrix 
#'                for sub-posterior c
#' @param prior_means prior for means of predictors
#' @param prior_variances prior for variances of predictors
#' @param time_schedule vector of length (L-1), where time_schedule[k] is time 
#'                      T for algorithm for level k
#' @param m_schedule vector of length (L-1), where m_schedule[k] is the number 
#'                   of samples to fuse for level k
#' @param C number of sub-posteriors at the base level
#' @param precondition logical value determining whether or not a 
#'                     preconditioning matrix is to be used
#' @param L total number of levels in the hierarchy
#' @param base_samples list of length C, where base_samples[[c]] contains
#'                     the samples for the c-th node in the level
#' @param resampling_method method to be used in resampling, default is 
#'                          multinomial resampling ('multi'). Other choices are 
#'                          stratified ('strat'), systematic ('system'), 
#'                          residual ('resid')
#' @param ESS_threshold number between 0 and 1 defining the proportion 
#'                      of the number of samples that ESS needs to be
#'                      lower than for resampling (i.e. resampling is carried 
#'                      out only when ESS < N*ESS_threshold)
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
                                        dim,
                                        y_split,
                                        X_split,
                                        prior_means,
                                        prior_variances,
                                        time_schedule,
                                        m_schedule,
                                        C,
                                        precondition = FALSE,
                                        L,
                                        base_samples,
                                        resampling_method = 'multi',
                                        ESS_threshold = 0.5,
                                        diffusion_estimator = 'Poisson',
                                        beta_NB = 10,
                                        bounds_multiplier = 1.2,
                                        seed = NULL,
                                        n_cores = parallel::detectCores()) {
  # check variables are of the right length
  if (length(N_schedule) != (L-1)) {
    stop('hierarchical_fusion_SMC_BLR: length of N_schedule must be equal to (L-1)')
  } else if (!is.list(y_split) || length(y_split)!= C) {
    stop('hierarchical_fusion_SMC_BLR: check that y_split is a list of length C')
  } else if (!is.list(X_split) || length(X_split)!= C) {
    stop('hierarchical_fusion_SMC_BLR: check that X_split is a list of length C')
  } else if (length(time_schedule) != (L-1)) {
    stop('hierarchical_fusion_SMC_BLR: length of time_schedule must be equal to (L-1)')
  } else if (length(base_samples) != C) {
    stop('hierarchical_fusion_SMC_BLR: length of base_samples must be equal to C')
  }
  # output warning to say that top level does not have C=1
  if (C != prod(m_schedule)) {
    warning('hierarchical_fusion_SMC_BLR: C != prod(m_schedule) - the top level 
            does not have C=1')
  }
  # check that at each level, we are fusing a suitable number
  if (length(m_schedule) == (L-1)) {
    for (l in (L-1):1) {
      if ((C/prod(m_schedule[(L-1):l]))%%1 != 0) {
        stop('hierarchical_fusion_SMC_BLR: check that C/prod(m_schedule[(L-1):l]) 
             is an integer for l=L-1,...,1')
      }
    }
  } else {
    stop('hierarchical_fusion_SMC_BLR: m_schedule must be a vector of length (L-1)')
  }
  # check ESS_threshold is strictly between 0 and 1
  if (ESS_threshold < 0 || ESS_threshold > 1) {
    stop('hierarchical_fusion_SMC_BLR: ESS_threshold must be between 0 and 1')
  }
  
  # we append 1 to the vector m_schedule to make the indices work later on when we call fusion
  m_schedule <- c(m_schedule, 1)
  # initialising results that we want to keep
  particles <- list()
  particles[[L]] <- initialise_particle_sets(samples_to_fuse = base_samples, 
                                             multivariate = TRUE)
  proposed_samples <- list()
  y_inputs <- list()
  y_inputs[[L]] <- y_split # base level input samples for y
  X_inputs <- list()
  X_inputs[[L]] <- X_split # base level input samples for X
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
  
  # add to output file that starting hierarchical fusion
  cat('Starting hierarchical fusion \n', file = 'hierarchical_fusion_SMC_BLR.txt')
  # parallelising tasks for each level going up the hiearchy
  for (k in ((L-1):1)) {
    # since previous level has C/prod(m_schedule[L:(k-1)]) nodes and we fuse m_schedule[k] of these
    n_nodes <- C/prod(m_schedule[L:k])
    # performing Fusion for this level
    # printing out some stuff to log file to track the progress
    cat('########################\n', file = 'hierarchical_fusion_SMC_BLR.txt', append = T)
    cat('Starting to fuse', m_schedule[k], 'sub-posteriors for level', k, 'with time',
        time_schedule[k], ', which is using', n_cores, 'cores\n',
        file = 'hierarchical_fusion_SMC_BLR.txt', append = T)
    cat('At this level, the data is split up into', (C / prod(m_schedule[L:(k+1)])), 'subsets\n',
        file = 'hierarchical_fusion_SMC_BLR.txt', append = T)
    cat('There are', n_nodes, 'nodes at the next level each giving', N_schedule[k],
        'samples \n', file = 'hierarchical_fusion_SMC_BLR.txt', append = T)
    cat('########################\n', file = 'hierarchical_fusion_SMC_BLR.txt', append = T)
    # starting fusion
    fused <- lapply(X = 1:n_nodes, FUN = function(i) {
      previous_nodes <- ((m_schedule[k]*i)-(m_schedule[k]-1)):(m_schedule[k]*i)
      particles_to_fuse <- particles[[k+1]][previous_nodes]
      precondition_mats <- precondition_matrices[[k+1]][previous_nodes]
      pcm <- proc.time()
      return(c(parallel_fusion_SMC_BLR(particles_to_fuse = particles_to_fuse,
                                       N = N_schedule[k],
                                       dim = dim,
                                       y_split = y_inputs[[k+1]][previous_nodes],
                                       X_split = X_inputs[[k+1]][previous_nodes],
                                       prior_means = prior_means,
                                       prior_variances = prior_variances,
                                       time = time_schedule[k],
                                       m = m_schedule[k],
                                       C = (C / prod(m_schedule[L:(k+1)])),
                                       precondition_matrices = precondition_mats,
                                       resampling_method = resampling_method,
                                       ESS_threshold = ESS_threshold,
                                       diffusion_estimator = diffusion_estimator,
                                       beta_NB = beta_NB,
                                       bounds_multiplier = bounds_multiplier,
                                       seed = seed,
                                       level = k,
                                       node = i,
                                       n_cores = n_cores),
               'time' = (proc.time()-pcm)['elapsed']))
    })
    # need to combine the correct samples
    particles[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$particles)
    proposed_samples[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$proposed_samples)
    time[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$time)
    ESS[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$ESS)
    CESS[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$CESS)
    # in a hierarchical setting, samples are always resampled after Q for intermediate levels
    if (k != 1) {
      resampled[[k]] <- lapply(1:n_nodes, function(i) {
        c(fused[[i]]$resampled['rho'], 'Q' = TRUE)})
    } else {
      resampled[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$resampled)
    }
    precondition_matrices[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$precondition_matrices[[1]])
    y_inputs[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$combined_y)
    X_inputs[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$combined_X)
  }
  
  # print completion
  cat('Completed hierarchical fusion\n', file = 'hierarchical_fusion_SMC_BLR.txt',
      append = T)
  if (C == prod(m_schedule)) {
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
              'resampling_method' = resampling_method,
              'y_inputs' = y_inputs,
              'X_inputs' = X_inputs,
              'diffusion_times' = time_schedule))
}

