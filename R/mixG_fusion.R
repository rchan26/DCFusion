#' Obtain bounds for phi function
#'
#' Finds the lower and upper bounds of the phi function between an interval
#'
#' @param n_comp integer number of components of mixture Gaussian
#' @param weights vector: weights of mixture Gaussian
#' @param means vector: means of mixture Gassuan
#' @param sds vector: st.devs of mixture Gaussian
#' @param beta real value
#' @param precondition precondition value
#' @param PHI global lower bound of phi
#' @param lower lower end of interval
#' @param upper upper end of interval
#' @param bounds_multiplier scalar value to mulitply bounds by 
#'                          (should greater than or equal to 1)
#' @param seq_length desired length of the sequence that we compute phi in the
#'                   interval [lower, upper]. Default is 1000
#'
#' @return vector of lower and upper bound of phi-function between [lower, upper]
#'
#' @examples
#' weights <- c(0.4, 0.6)
#' means <- c(-3, 6)
#' sds <- c(1, 2)
#' beta <- 0.5
#' precondition <- 1.5
#' lower <- -10
#' upper <- 4
#' curve(ea_phi_mixG_DL(x,
#'                      n_comp = 2,
#'                      weights = weights,
#'                      means = means,
#'                      sds = sds,
#'                      beta = beta,
#'                      precondition = precondition),
#'       from = -20, to = 25, n = 10000,
#'       ylab = 'phi')
#' PHI <- ea_phi_mixG_DL_LB(n_comp = 2,
#'                          weights = weights,
#'                          means = means,
#'                          sds = sds,
#'                          beta = beta,
#'                          precondition = precondition,
#'                          bounds_multiplier = 1)
#' bounds <- ea_phi_mixG_DL_bounds(n_comp = 2,
#'                                 weights = weights,
#'                                 means = means,
#'                                 sds = sds,
#'                                 beta = beta,
#'                                 precondition = precondition,
#'                                 PHI = PHI,
#'                                 lower = lower,
#'                                 upper = upper,
#'                                 bounds_multiplier = 1)
#' abline(v=c(lower, upper))
#' abline(h=PHI, col = 'red')
#' abline(h=c(bounds$UB, bounds$LB), col = 'blue', lty = 2)
#'
#' @export
ea_phi_mixG_DL_bounds <- function(n_comp,
                                  weights,
                                  means,
                                  sds,
                                  beta,
                                  precondition,
                                  PHI,
                                  lower,
                                  upper,
                                  bounds_multiplier,
                                  seq_length = 1000) {
  if (bounds_multiplier < 1) {
    stop("ea_phi_mixG_DL_bounds: bounds_multipler should be greater than or equal to 1")
  }
  # ---------- finding maxima and minimma in the interval
  # first just evaluate the phi function along a sequence
  sequence_phi_eval <- ea_phi_mixG_DL(x = seq(lower, upper, length.out = seq_length),
                                      n_comp = n_comp,
                                      weights = weights,
                                      means = means,
                                      sds = sds,
                                      beta = beta,
                                      precondition = precondition)
  # also calculate bounds using an optimiser implemented in base R
  LB <- min(c(sequence_phi_eval, 
              optimise(f = function(x) ea_phi_mixG_DL(x, n_comp = n_comp, weights, means, sds, beta, precondition),
                       interval = c(lower, upper),
                       maximum = FALSE, 
                       tol = .Machine$double.eps^0.8)$objective))
  UB <- max(c(sequence_phi_eval,
              optimise(f = function(x) ea_phi_mixG_DL(x, n_comp = n_comp, weights, means, sds, beta, precondition),
                       interval = c(lower, upper),
                       maximum = TRUE, 
                       tol = .Machine$double.eps^0.8)$objective))
  # multiply the bounds by bounds_multiplier to compensate the case 
  # that the opitimser did not get the bounds exactly correct
  # calculate the Bounds Difference
  BD <- UB-LB
  # by subtracting and adding 0.5*(bounds_multiplier-1)*BD
  # makes the resulting bounds difference be larger by a factor of bounds_multiplier
  # for lower bound, make sure its not less than the global lower bound PHI
  return(list('LB' = max(LB - 0.5*(bounds_multiplier-1)*BD, PHI),
              'UB' = UB + 0.5*(bounds_multiplier-1)*BD))
}

#' Exact Algorithm for langevin diffusion with pi as a tempered mixture Gaussian
#'
#' Simulate langevin diffusion using the Exact Algorithm with pi = exp(-(beta*x^4)/2)
#'
#' @param x0 start value
#' @param y end value
#' @param s start time
#' @param t end time
#' @param n_comp integer number of components of mixture Gaussian
#' @param weights vector: weights of mixture Gaussian
#' @param means vector: means of mixture Gassuan
#' @param sds vector: st.devs of mixture Gaussian
#' @param beta real value
#' @param precondition precondition value (i.e the covariance for 
#'                     the Langevin diffusion)
#' @param bounds_multiplier scalar value to mulitply bounds by 
#'                          (should greater than or equal to 1)
#' @param logarithm logical value to determine if log probability is 
#'                  returned (TRUE) or not (FALSE)
#'                  
#' @examples 
#' weights <- c(0.4, 0.6)
#' means <- c(-8, 15)
#' sds <- c(1, 2)
#' beta <- 1/4
#' precondition <- 1.5
#' x0 <- -9
#' y <- -8
#' s <- 0
#' t <- 1
#' # simulate event probability of diffusion from -9 to -8 between [0,1]
#' ea_mixG_DL_PT(x0 = x0,
#'               y = y,
#'               s = s,
#'               t = t,
#'               n_comp = 2,
#'               weights = weights,
#'               means = means,
#'               sds = sds,
#'               beta = beta,
#'               precondition = precondition,
#'               bounds_multiplier = 1,
#'               logarithm = FALSE)
#' # simulate event probability of diffusion from -9 to 15 between [0,1]
#' x0 <- -9
#' y <- 15
#' ea_mixG_DL_PT(x0 = x0,
#'               y = y,
#'               s = s,
#'               t = t,
#'               n_comp = 2,
#'               weights = weights,
#'               means = means,
#'               sds = sds,
#'               beta = beta,
#'               precondition = precondition,
#'               bounds_multiplier = 1,
#'               logarithm = FALSE)
#'
#' @export
ea_mixG_DL_PT <- function(x0,
                          y,
                          s,
                          t,
                          n_comp,
                          weights,
                          means,
                          sds,
                          beta,
                          precondition,
                          bounds_multiplier = 1.1,
                          logarithm) {
  # transform to preconditoned space
  z0 <- x0 / sqrt(precondition)
  zt <- y / sqrt(precondition)
  # simulate layer information
  bes_layer <- layeredBB::bessel_layer_simulation(x = z0,
                                                  y = zt,
                                                  s = s,
                                                  t = t,
                                                  mult = 0.1)
  lbound_X <- sqrt(precondition) * bes_layer$L
  ubound_X <- sqrt(precondition) * bes_layer$U
  # calculate upper and lower bounds of phi given the simulated sample layer information
  PHI <- ea_phi_mixG_DL_LB(n_comp = n_comp,
                           weights = weights,
                           means = means,
                           sds = sds,
                           beta = beta,
                           precondition = precondition,
                           bounds_multiplier = bounds_multiplier)
  bounds <- ea_phi_mixG_DL_bounds(n_comp = n_comp,
                                  weights = weights, 
                                  means = means,
                                  sds = sds,
                                  beta = beta,
                                  precondition = precondition,
                                  PHI = PHI,
                                  lower = lbound_X,
                                  upper = ubound_X,
                                  bounds_multiplier = bounds_multiplier,
                                  seq_length = 1000)
  LX <- bounds$LB
  UX <- bounds$UB
  # simulate number of points to simulate from Poisson distribution
  kap <- rpois(n = 1, lambda = (UX-LX)*(t-s))
  log_acc_prob <- 0
  if (kap > 0) {
    layered_bb <- layeredBB::layered_brownian_bridge(x = z0,
                                                     y = zt,
                                                     s = s,
                                                     t = t,
                                                     bessel_layer = bes_layer,
                                                     times = runif(kap, s, t))
    phi <- ea_phi_mixG_DL(x = sqrt(precondition) * layered_bb$simulated_path[1,],
                          n_comp = n_comp,
                          weights = weights,
                          means = means,
                          sds = sds,
                          beta = beta,
                          precondition = precondition)
    log_acc_prob <- sum(log(UX-phi))
  }
  if (logarithm) {
    return(-(LX-PHI)*(t-s) - kap*log(UX-LX) + log_acc_prob)
  } else {
    return(exp(-(LX-PHI)*(t-s) - kap*log(UX-LX) + log_acc_prob))
  }
}

#' Preconditioned Monte Carlo Fusion [on a single core]
#'
#' Monte Carlo Fusion for sub-posteriors which are tempered mixture Gaussians
#'
#' @param N number of samples
#' @param m number of sub-posteriors to combine
#' @param time time T for fusion algorithm
#' @param samples_to_fuse list of length m, where samples_to_fuse[c] containing
#'                        the samples for the c-th sub-posterior
#' @param n_comp integer number of components of mixture Gaussian
#' @param weights vector: weights of mixture Gaussian
#' @param means vector: means of mixture Gassuan
#' @param sds vector: st.devs of mixture Gaussian
#' @param betas vector of length m, where betas[c] is the inverse temperature
#'              (beta) for c-th sub-posterior (can also pass in one number if
#'              they are all at the same inverse temperature)
#' @param precondition_values vector of length m, where precondition_values[c]
#'                            is the precondition value for sub-posterior c
#' @param bounds_multiplier scalar value to mulitply bounds by 
#'                          (should greater than or equal to 1)
#' @param seed seed number - default is NULL, meaning there is no seed
#' @param level integer for level in hierarchy - default 1
#' @param node integer for node in level in hierarchy - default 1
#' @param n_cores number of cores to use
#'
#' @return A list with components:
#' \describe{
#'   \item{samples}{fusion samples}
#'   \item{iters_rho}{number of iterations for rho step}
#'   \item{iters_Q}{number of iterations for Q step}
#' }
#'
#' @export
fusion_mixG <- function(N,
                        m,
                        time,
                        samples_to_fuse,
                        n_comp,
                        weights,
                        means,
                        sds,
                        betas,
                        precondition_values,
                        bounds_multiplier = 1.1,
                        level = 1,
                        node = 1,
                        core = 1) {
  if (length(weights)!=n_comp) {
    stop("fusion_mixG: weights must be a vector of length n_comp")
  } else if (length(means)!=n_comp) {
    stop("fusion_mixG: means must be a vector of length n_comp")
  } else if (length(sds)!=n_comp) {
    stop("fusion_mixG: sds must be a vector of length n_comp")
  } else if (!is.list(samples_to_fuse) | length(samples_to_fuse)!=m) {
    stop("fusion_mixG: samples_to_fuse must be a list of length m")
  } else if (!is.vector(betas) | (length(betas)!=m)) {
    stop("fusion_mixG: betas must be a vector of length m")
  } else if (!is.vector(precondition_values) | length(precondition_values)!=m) {
    stop("fusion_mixG: precondition_values must be a vector of length m")
  }
  fusion_samples <- rep(NA, N); i <- 0; iters_rho <- 0; iters_Q <- 0
  proposal_sd <- sqrt(time / sum(1/precondition_values))
  while (i < N) {
    iters_rho <- iters_rho+1
    x <- sapply(samples_to_fuse, function(core) sample(x = core, size = 1))
    weighted_avg <- weighted_mean_univariate(x = x, weights = 1/precondition_values)
    log_rho_prob <- log_rho_univariate(x = x,
                                       weighted_mean = weighted_avg,
                                       time = time,
                                       precondition_values = precondition_values)
    if (log(runif(1, 0, 1)) < log_rho_prob) {
      iters_Q <- iters_Q + 1
      y <- rnorm(n = 1, mean = weighted_avg, sd = proposal_sd)
      log_Q_prob <- sum(sapply(1:m, function(c) {
        ea_mixG_DL_PT(x0 = x[c],
                      y = y,
                      s = 0,
                      t = time,
                      n_comp = n_comp,
                      weights = weights,
                      means = means,
                      sds = sds,
                      beta = betas[c],
                      precondition = precondition_values[c],
                      bounds_multiplier = bounds_multiplier,
                      logarithm = TRUE)}))
      if (log(runif(1,0,1)) < log_Q_prob) {
        i <- i+1
        fusion_samples[i] <- y
        cat('Level:', level, '||', 'Node:', node, '|| Core:', core, '|| Iteration:', i, '/', N, '\n',
            file = 'fusion_mixG_progress.txt', append = T)
      } 
    }
  }
  cat('Node', node, '|| Core:', core, '|| Completed fusion for level', level, '\n',
      file = 'fusion_mixG_progress.txt', append = T)
  return(list('samples' = fusion_samples,
              'iters_rho' = iters_rho,
              'iters_Q' = iters_Q))
}

#' Preconditioned Monte Carlo Fusion [parallel]
#'
#' Monte Carlo Fusion for sub-posteriors which are tempered mixture Gaussians
#'
#' @param N number of samples
#' @param m number of sub-posteriors to combine
#' @param time time T for fusion algorithm
#' @param samples_to_fuse list of length m, where samples_to_fuse[c] containing
#'                        the samples for the c-th sub-posterior
#' @param n_comp integer number of components of mixture Gaussian
#' @param weights vector: weights of mixture Gaussian
#' @param means vector: means of mixture Gassuan
#' @param sds vector: st.devs of mixture Gaussian
#' @param betas vector of length m, where betas[c] is the inverse temperature
#'              (beta) for c-th sub-posterior (can also pass in one number if
#'              they are all at the same inverse temperature)
#' @param precondition_values vector of length m, where precondition_values[c]
#'                            is the precondition value for sub-posterior c
#' @param bounds_multiplier scalar value to mulitply bounds by 
#'                          (should greater than or equal to 1)
#' @param seed seed number - default is NULL, meaning there is no seed
#' @param n_cores number of cores to use
#' @param level integer for level in hierarchy - default 1
#' @param node integer for node in level in hierarchy - default 1
#'
#' @return A list with components:
#' \describe{
#'   \item{samples}{fusion samples}
#'   \item{rho}{acceptance rate for rho step}
#'   \item{Q}{acceptance rate for Q step}
#'   \item{rhoQ}{overall acceptance rate}
#'   \item{time}{run-time of fusion sampler}
#'   \item{rho_iterations}{number of iterations for rho step}
#'   \item{Q_iterations}{number of iterations for Q step}
#'   \item{precondition_values}{list of length 2 where precondition_values[[2]]
#'                              are the pre-conditioning values that were used
#'                              and precondition_values[[1]] are the combined
#'                              precondition values}
#' }
#'
#' @export
parallel_fusion_mixG <- function(N,
                                 m,
                                 time,
                                 samples_to_fuse,
                                 n_comp,
                                 weights,
                                 means,
                                 sds,
                                 betas,
                                 precondition_values,
                                 bounds_multiplier = 1.1,
                                 seed = NULL,
                                 n_cores = parallel::detectCores(),
                                 level = 1,
                                 node = 1) {
  if (length(weights)!=n_comp) {
    stop("parallel_fusion_mixG: weights must be a vector of length n_comp")
  } else if (length(means)!=n_comp) {
    stop("parallel_fusion_mixG: means must be a vector of length n_comp")
  } else if (length(sds)!=n_comp) {
    stop("parallel_fusion_mixG: sds must be a vector of length n_comp")
  } else if (!is.list(samples_to_fuse) | length(samples_to_fuse)!=m) {
    stop("parallel_fusion_mixG: samples_to_fuse must be a list of length m")
  } else if (!is.vector(betas) | (length(betas)!=m)) {
    stop("parallel_fusion_mixG: betas must be a vector of length m")
  } else if (!is.vector(precondition_values) | length(precondition_values)!=m) {
    stop("parallel_fusion_mixG: precondition_values must be a vector of length m")
  }
  # ---------- creating parallel cluster
  cl <- parallel::makeCluster(n_cores, setup_strategy = "sequential")
  parallel::clusterExport(cl, envir = environment(),
                          varlist = c(ls(), list("ea_phi_mixG_DL",
                                                 "ea_phi_mixG_DL_LB",
                                                 "ea_phi_mixG_DL_bounds",
                                                 "ea_mixG_DL_PT",
                                                 "log_rho_univariate",
                                                 "fusion_mixG")))
  # exporting functions from layeredBB package to simulate layered Brownian bridges
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
  fused <- parallel::parLapply(cl, X = 1:length(samples_per_core), fun = function(i) {
    fusion_mixG(N = samples_per_core[i],
                m = m,
                time = time,
                samples_to_fuse = samples_to_fuse,
                n_comp = n_comp,
                weights = weights,
                means = means,
                sds = sds,
                betas = betas,
                precondition_values = precondition_values,
                bounds_multiplier = bounds_multiplier,
                level = level,
                node = node,
                core = i)
  })
  final <- proc.time() - pcm
  parallel::stopCluster(cl)
  # ---------- return samples and acceptance probabilities
  samples <- unlist(lapply(1:length(samples_per_core), function(i) fused[[i]]$samples))
  rho_iterations <- sum(sapply(1:length(samples_per_core), function(i) fused[[i]]$iters_rho))
  Q_iterations <- sum(sapply(1:length(samples_per_core), function(i) fused[[i]]$iters_Q))
  rho_acc <- Q_iterations / rho_iterations
  Q_acc <- N / Q_iterations
  rhoQ_acc <- N / rho_iterations
  if (identical(precondition_values, rep(1, m))) {
    return(list('samples' = samples,
                'rho' = rho_acc,
                'Q' = Q_acc,
                'rhoQ '= rhoQ_acc,
                'time' = final['elapsed'],
                'rho_iterations' = rho_iterations,
                'Q_iterations' = Q_iterations,
                'precondition_values' = list(1, precondition_values)))
  } else {
    return(list('samples' = samples,
                'rho' = rho_acc,
                'Q' = Q_acc,
                'rhoQ '= rhoQ_acc,
                'time' = final['elapsed'],
                'rho_iterations' = rho_iterations,
                'Q_iterations' = Q_iterations,
                'precondition_values' = list(1/sum(1/precondition_values),
                                             precondition_values)))
  }
}

#' Hierarchical Monte Carlo Fusion
#'
#' Hierarchical Monte Carlo Fusion with base level with nodes that are
#' tempered mixture Gaussians
#'
#' @param N_schedule vector of length (L-1), where N_schedule[l] is the number
#'                   of samples per node at level l
#' @param m_schedule vector of length (L-1), where m_schedule[l] is the number
#'                   of samples to fuse for level l
#' @param time_schedule vector of legnth(L-1), where time_schedule[l] is the
#'                      time chosen for Fusion at level l
#' @param base_samples list of length (1/start_beta), where samples_to_fuse[c]
#'                     contains the samples for the c-th node in the level
#' @param L total number of levels in the hierarchy
#' @param n_comp integer number of components of mixture Gaussian
#' @param weights vector: weights of mixture Gaussian
#' @param means vector: means of mixture Gaussian
#' @param sds vector: st.devs of mixture Gaussian
#' @param start_beta beta for the base level
#' @param precondition either a logical value to determine if preconditioning values are
#'                     used (TRUE - and is set to be the variance of the sub-posterior samples)
#'                     or not (FALSE - and is set to be 1 for all sub-posteriors),
#'                     or a list of length (1/start_beta) where precondition[[c]]
#'                     is the preconditioning value for sub-posterior c. Default is TRUE
#' @param bounds_multiplier scalar value to mulitply bounds by 
#'                          (should greater than or equal to 1)
#' @param seed seed number - default is NULL, meaning there is no seed
#' @param n_cores number of cores to use
#'
#' @return A list with components:
#' \describe{
#'   \item{samples}{list of length (L-1), where samples[[l]][[i]] are the samples
#'                  for level l, node i}
#'   \item{time}{list of length (L-1), where time[[l]] is the run time for level
#'               l, node i}
#'   \item{rho_acc}{list of length (L-1), where rho_acc[[l]][i] is the acceptance
#'                  rate for first fusion step for level l, node i}
#'   \item{Q_acc}{list of length (L-1), where Q_acc[[l]][i] is the acceptance
#'                rate for second fusion step for level l, node i}
#'   \item{rhoQ_acc}{list of length (L-1), where rhoQ_acc[[l]][i] is the overall
#'                   acceptance rate for fusion for level l, node i}
#'   \item{diffusion_times}{vector of length (L-1), where diffusion_times[l] are
#'                          the times for fusion in level l (= time_schedule)}
#'   \item{overall_rho}{vector of length (L-1), where overall_rho[k] is overall
#'                      acceptance rate for rho of level l}
#'   \item{overall_Q}{vector of length (L-1), where overall_Q[k] is overall
#'                    acceptance rate for Q of level l}
#'   \item{overall_rhoQ}{vector of length (L-1), where overall_rhoQ[k] is overall
#'                       acceptance rate for rho*Q of level l}
#'   \item{overall_time}{vector of length (L-1), where overall_time[k] is
#'                       overall taken for level l}
#'   \item{precondition_values}{preconditioning values used in the algorithm
#'                              for each node}
#' }
#'
#' @export
hierarchical_fusion_mixG <- function(N_schedule,
                                     m_schedule,
                                     time_schedule,
                                     base_samples,
                                     L,
                                     n_comp,
                                     weights,
                                     means,
                                     sds,
                                     start_beta,
                                     precondition = TRUE,
                                     bounds_multiplier = 1.1,
                                     seed = NULL,
                                     n_cores = parallel::detectCores()) {
  if (length(weights)!=n_comp) {
    stop("hierarchical_fusion_mixG: weights must be a vector of length n_comp")
  } else if (length(means)!=n_comp) {
    stop("hierarchical_fusion_mixG: means must be a vector of length n_comp")
  } else if (length(sds)!=n_comp) {
    stop("hierarchical_fusion_mixG: sds must be a vector of length n_comp")
  } else if (!is.vector(N_schedule) | (length(N_schedule)!=(L-1))) {
    stop("hierarchical_fusion_mixG: N_schedule must be a vector of length (L-1)")
  } else if (!is.vector(m_schedule) | (length(m_schedule)!=(L-1))) {
    stop("hierarchical_fusion_mixG: m_schedule must be a vector of length (L-1)")
  } else if (!is.vector(time_schedule) | (length(time_schedule)!=(L-1))) {
    stop("hierarchical_fusion_mixG: time_schedule must be a vector of length (L-1)")
  } else if (!is.list(base_samples) | (length(base_samples)!=(1/start_beta))) {
    stop("hierarchical_fusion_mixG: base_samples must be a list of length (1/start_beta)")
  }
  if (is.vector(m_schedule) & (length(m_schedule)==(L-1))) {
    for (l in (L-1):1) {
      if (((1/start_beta)/prod(m_schedule[(L-1):l]))%%1!=0) {
        stop("hierarchical_fusion_mixG: check that (1/start_beta)/prod(m_schedule[(L-1):l])
             is an integer for l=L-1,...,1")
      }
    }
  } else {
    stop("hierarchical_fusion_mixG: m_schedule must be a vector of length (L-1)")
  }
  # we append 1 to the vector m_schedule to make the indices work later on when we call fusion
  # we need this so that we can set the right value for beta when fusing up the levels
  m_schedule <- c(m_schedule, 1)
  # initialising study results
  hier_samples <- list()
  hier_samples[[L]] <- base_samples # base level
  time <- list()
  rho <- list()
  Q <- list()
  rhoQ <- list()
  overall_rho <- rep(0, L-1)
  overall_Q <- rep(0, L-1)
  overall_rhoQ <- rep(0, L-1)
  overall_time <- rep(0, L-1)
  precondition_values <- list()
  if (is.logical(precondition)) {
    if (precondition) {
      precondition_values[[L]] <- lapply(base_samples, var)
    } else {
      precondition_values[[L]] <- lapply(1:length(base_samples), function(i) 1)
    }
  } else if (is.list(precondition)) {
    if (length(precondition)==(1/start_beta)) {
      precondition_values[[L]] <- precondition
    }
  } else {
    stop("hierarchical_fusion_mixG: precondition must be a logical indicating 
          whether or not a preconditioning value should be used, or a list of
          length C, where precondition[[c]] is the preconditioning value for
          the c-th sub-posterior")
  }
  cat('Starting hierarchical fusion \n', file = 'hierarchical_fusion_mixG.txt')
  for (k in ((L-1):1)) {
    # since previous level has (1/beta)/prod(m_schedule[L:(k-1)]) nodes and we 
    # fuse m_schedule[k] of these
    n_nodes <- max((1/start_beta)/prod(m_schedule[L:k]), 1)
    cat('########################\n', file = 'hierarchical_fusion_mixG.txt',
        append = T)
    cat('Starting to fuse', m_schedule[k], 'densities for level', k, 'with time',
        time_schedule[k], '- using', parallel::detectCores(), 'cores\n',
        file = 'hierarchical_fusion_mixG.txt', append = T)
    cat('There are', n_nodes, 'nodes at this level each giving', N_schedule[k],
        'samples for beta =', prod(m_schedule[L:k]), '/', (1/start_beta),
        '\n', file = 'hierarchical_fusion_mixG.txt', append = T)
    cat('########################\n', file = 'hierarchical_fusion_mixG.txt',
        append = T)
    fused <- lapply(X = 1:n_nodes, FUN = function(i) {
      previous_nodes <- ((m_schedule[k]*i)-(m_schedule[k]-1)):(m_schedule[k]*i)
      precondition_vals <- unlist(precondition_values[[k+1]][previous_nodes])
      parallel_fusion_mixG(N = N_schedule[k],
                           m = m_schedule[k],
                           time = time_schedule[k],
                           samples_to_fuse = hier_samples[[k+1]][previous_nodes],
                           n_comp = n_comp,
                           weights = weights,
                           means = means,
                           sds = sds,
                           betas = prod(m_schedule[L:(k+1)])*(start_beta),
                           precondition_values = precondition_vals,
                           bounds_multiplier = bounds_multiplier,
                           seed = seed,
                           n_cores = n_cores,
                           level = k,
                           node = i)
    })
    # need to combine the correct samples
    hier_samples[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$samples)
    # obtaining the acceptance rates for all nodes in the current level
    rho[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$rho)
    Q[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$Q)
    rhoQ[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$rhoQ)
    time[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$time)
    # get number of iterations in this level
    sum_rho_iterations <- sum(unlist(lapply(1:n_nodes, function(i) fused[[i]]$rho_iterations)))
    sum_Q_iterations <- sum(unlist(lapply(1:n_nodes, function(i) fused[[i]]$Q_iterations)))
    # acceptance rate for whole level
    overall_rho[k] <- sum_Q_iterations / sum_rho_iterations
    overall_Q[k] <- N_schedule[k]*n_nodes / sum_Q_iterations
    overall_rhoQ[k] <- N_schedule[k]*n_nodes / sum_rho_iterations
    overall_time[k] <- sum(unlist(time[[k]]))
    precondition_values[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$precondition_values[[1]])
  }
  cat('Completed hierarchical fusion\n', file = 'hierarchical_fusion_mixG.txt',
      append = T)
  if (length(hier_samples[[1]])==1) {
    hier_samples[[1]] <- hier_samples[[1]][[1]]
    time[[1]] <- time[[1]][[1]]
    rho[[1]] <- rho[[1]][[1]]
    Q[[1]] <- Q[[1]][[1]]
    rhoQ[[1]] <- rhoQ[[1]][[1]]
  }
  return(list('samples' = hier_samples,
              'time' = time,
              'rho_acc' = rho,
              'Q_acc' = Q,
              'rhoQ_acc' = rhoQ,
              'diffusion_times' = time_schedule,
              'precondition_values' = precondition_values,
              'overall_rho' = overall_rho,
              'overall_Q' = overall_Q,
              'overall_rhoQ' = overall_rhoQ,
              'overall_time' = overall_time))
}

#' Progressive Monte Carlo Fusion
#'
#' Progressive Monte Carlo Fusion with base level with nodes that are
#' tempered mixture Gaussians
#'
#' @param N_schedule vector of length (L-1), where N_schedule[l] is the number
#'                   of samples per node at level l
#' @param time_schedule vector of legnth(L-1), where time_schedule[l] is the
#'                      time chosen for Fusion at level l
#' @param base_samples list of length (1/start_beta), where samples_to_fuse[c]
#'                     contains the samples for the c-th node in the level
#' @param n_comp integer number of components of mixture Gaussian
#' @param weights vector: weights of mixture Gaussian
#' @param means vector: means of mixture Gassuan
#' @param sds vector: st.devs of mixture Gaussian
#' @param start_beta beta for the base level
#' @param precondition either a logical value to determine if preconditioning values are
#'                     used (TRUE - and is set to be the variance of the sub-posterior samples)
#'                     or not (FALSE - and is set to be 1 for all sub-posteriors),
#'                     or a list of length (1/start_beta) where precondition[[c]]
#'                     is the preconditioning value for sub-posterior c. Default is TRUE
#' @param bounds_multiplier scalar value to mulitply bounds by 
#'                          (should greater than or equal to 1)
#' @param seed seed number - default is NULL, meaning there is no seed
#' @param n_cores number of cores to use
#'
#' @return A list with components:
#' \describe{
#'   \item{samples}{list of length (L-1), where samples[[l]][[i]] are the samples
#'                  for level l, node i}
#'   \item{time}{list of length (L-1), where time[[l]] is the run time for level
#'               l, node i}
#'   \item{rho_acc}{list of length (L-1), where rho_acc[[l]][i] is the acceptance
#'                  rate for first fusion step for level l, node i}
#'   \item{Q_acc}{list of length (L-1), where Q_acc[[l]][i] is the acceptance
#'                rate for second fusion step for level l, node i}
#'   \item{rhoQ_acc}{list of length (L-1), where rhoQ_acc[[l]][i] is the overall
#'                   acceptance rate for fusion for level l, node i}
#'   \item{diffusion_times}{vector of length (L-1), where diffusion_times[l] are
#'                          the times for fusion in level l (= time_schedule)}
#'   \item{precondition_values}{preconditioning values used in the algorithm
#'                              for each node}
#' }
#'
#' @export
progressive_fusion_mixG <- function(N_schedule,
                                    time_schedule,
                                    base_samples,
                                    n_comp,
                                    weights,
                                    means,
                                    sds,
                                    start_beta,
                                    precondition = TRUE,
                                    bounds_multiplier = 1.1,
                                    seed = NULL,
                                    n_cores = parallel::detectCores()) {
  if (length(weights)!=n_comp) {
    stop("progressive_fusion_mixG: weights must be a vector of length n_comp")
  } else if (length(means)!=n_comp) {
    stop("progressive_fusion_mixG: means must be a vector of length n_comp")
  } else if (length(sds)!=n_comp) {
    stop("progressive_fusion_mixG: sds must be a vector of length n_comp")
  } else if (!is.vector(N_schedule) | (length(N_schedule)!=(1/start_beta)-1)) {
    stop("progressive_fusion_mixG: N_schedule must be a vector of length ((1/start_beta)-1)")
  } else if (!is.vector(time_schedule) | (length(time_schedule)!=(1/start_beta)-1)) {
    stop("progressive_fusion_mixG: time_schedule must be a vector of length ((1/start_beta)-1)")
  } else if (!is.list(base_samples) | (length(base_samples)!=(1/start_beta))) {
    stop("progressive_fusion_mixG: base_samples must be a list of length (1/start_beta)")
  }
  # initialising results
  prog_samples <- list()
  prog_samples[[(1/start_beta)]] <- base_samples # base level
  time <- rep(0, (1/start_beta)-1)
  rho <- rep(0, (1/start_beta)-1)
  Q <- rep(0, (1/start_beta)-1)
  rhoQ <- rep(0, (1/start_beta)-1)
  precondition_values <- list()
  if (is.logical(precondition)) {
    if (precondition) {
      precondition_values[[(1/start_beta)]] <- lapply(base_samples, var)
    } else {
      precondition_values[[(1/start_beta)]] <- lapply(1:length(base_samples), function(i) 1)
    }
  } else if (is.list(precondition)) {
    if (length(precondition)==(1/start_beta)) {
      precondition_values[[(1/start_beta)]] <- precondition
    }
  } else {
    stop("progressive_fusion_mixG: precondition must be a logical indicating 
          whether or not a preconditioning value should be used, or a list of
          length C, where precondition[[c]] is the preconditioning value for
          the c-th sub-posterior")
  }
  index <- 2
  cat('Starting progressive fusion \n', file = 'progressive_fusion_mixG.txt')
  for (k in ((1/start_beta)-1):1) {
    if (k==(1/start_beta)-1) {
      cat('########################\n', file = 'progressive_fusion_mixG.txt',
          append = T)
      cat('Starting to fuse', 2, 'densities for level', k, 'which is using',
          parallel::detectCores(), 'cores\n',
          file = 'progressive_fusion_mixG.txt', append = T)
      cat('Fusing samples for beta =', 1, '/', (1/start_beta), 'with time',
          time_schedule[k], 'to get', N_schedule[k], 'samples for beta =', 2,
          '/', (1/start_beta), '\n', file = 'progressive_fusion_mixG.txt',
          append = T)
      cat('########################\n', file = 'progressive_fusion_mixG.txt',
          append = T)
      samples_to_fuse <- list(base_samples[[1]], base_samples[[2]])
      precondition_vals <- unlist(precondition_values[[k+1]][1:2])
      fused <- parallel_fusion_mixG(N = N_schedule[k],
                                    m = 2,
                                    time = time_schedule[k],
                                    samples_to_fuse = samples_to_fuse,
                                    n_comp = n_comp,
                                    weights = weights,
                                    means = means,
                                    sds = sds,
                                    betas = c(start_beta, start_beta),
                                    precondition_values = precondition_vals,
                                    bounds_multiplier = bounds_multiplier,
                                    seed = seed,
                                    n_cores = n_cores,
                                    level = k)
    } else {
      # printing out some stuff to log file to track the progress
      cat('########################\n', file = 'progressive_fusion_mixG.txt',
          append = T)
      cat('Starting to fuse', 2, 'densities for level', k, 'which is using',
          parallel::detectCores(), 'cores\n',
          file = 'progressive_fusion_mixG.txt', append = T)
      cat('Fusing samples for beta =', index, '/', (1/start_beta), 'and beta =',
          1, '/', (1/start_beta), 'with time', time_schedule[k], 'to get',
          N_schedule[k], 'samples for beta =', (index+1), '/', (1/start_beta),
          '\n', file = 'progressive_fusion_mixG.txt', append = T)
      cat('########################\n', file = 'progressive_fusion_mixG.txt',
          append = T)
      # starting fusion
      samples_to_fuse <- list(prog_samples[[k+1]],
                              base_samples[[index+1]])
      precondition_vals <- c(precondition_values[[k+1]],
                             precondition_values[[(1/start_beta)]][[index+1]])
      fused <- parallel_fusion_mixG(N = N_schedule[k],
                                    m = 2,
                                    time = time_schedule[k],
                                    samples_to_fuse = samples_to_fuse,
                                    n_comp = n_comp,
                                    weights = weights,
                                    means = means,
                                    sds = sds,
                                    betas = c(index*start_beta, start_beta),
                                    precondition_values = precondition_vals,
                                    bounds_multiplier = bounds_multiplier,
                                    seed = seed,
                                    n_cores = n_cores,
                                    level = k)
      index <- index + 1
    }
    # need to combine the correct samples
    prog_samples[[k]] <- fused$samples
    precondition_values[[k]] <- fused$precondition_values[[1]]
    # obtaining the acceptance rates for all nodes in the current level
    rho[k] <- fused$rho
    Q[k] <- fused$Q
    rhoQ[k] <- fused$rhoQ
    time[k] <- fused$time
  }
  cat('Completed progressive fusion\n',
      file = 'progressive_fusion_mixG.txt', append = T)
  return(list('samples' = prog_samples,
              'time' = time,
              'rho_acc' = rho,
              'Q_acc' = Q,
              'rhoQ_acc' = rhoQ,
              'diffusion_times' = time_schedule,
              'precondition_values' = precondition_values))
}

#' Q Importance Sampling Step for sub-posteriors of the form exp(-((x^4)*beta)/2)
#'
#' Q Importance Sampling weighting for sub-posteriors of the form exp(-((x^4)*beta)/2)
#'
#' @param particle_set particles set prior to Q importance sampling step
#' @param m number of sub-posteriors to combine
#' @param time time T for fusion algorithm
#' @param n_comp integer number of components of mixture Gaussian
#' @param weights vector: weights of mixture Gaussian
#' @param means vector: means of mixture Gaussian
#' @param sds vector: st.devs of mixture Gaussian
#' @param betas vector of length c, where betas[c] is the inverse temperature 
#'              value for c-th posterior
#' @param precondition_values vector of length m, where precondition_values[c]
#'                            is the precondition value for sub-posterior c
#' @param bounds_multiplier scalar value to mulitply bounds by 
#'                          (should greater than or equal to 1)
#' @param seed seed number - default is NULL, meaning there is no seed
#' @param n_cores number of cores to use
#'
#' @return An updated particle set
#' 
#' @export
Q_IS_mixG <- function(particle_set,
                      m,
                      time,
                      n_comp,
                      weights,
                      means,
                      sds,
                      betas,
                      precondition_values,
                      bounds_multiplier = 1.1,
                      seed = NULL,
                      n_cores = parallel::detectCores()) {
  if (!("particle" %in% class(particle_set))) {
    stop("Q_IS_mixG: particle_set must be a \"particle\" object")
  } else if (length(weights)!=n_comp) {
    stop("Q_IS_mixG: weights must be a vector of length n_comp")
  } else if (length(means)!=n_comp) {
    stop("Q_IS_mixG: means must be a vector of length n_comp")
  } else if (length(sds)!=n_comp) {
    stop("Q_IS_mixG: sds must be a vector of length n_comp")
  } else if (!is.vector(betas) | (length(betas)!=m)) {
    stop("Q_IS_mixG: betas must be a vector of length m")
  } else if (!is.vector(precondition_values) | (length(precondition_values)!=m)) {
    stop("Q_IS_mixG: precondition_values must be a vector of length m")
  }
  proposal_sd <- sqrt(time / sum(1/precondition_values))
  N <- particle_set$N
  # ---------- creating parallel cluster
  cl <- parallel::makeCluster(n_cores, setup_strategy = "sequential")
  parallel::clusterExport(cl, envir = environment(),
                          varlist = c(ls(), list("ea_phi_mixG_DL",
                                                 "ea_phi_mixG_DL_LB",
                                                 "ea_phi_mixG_DL_bounds",
                                                 "ea_mixG_DL_PT")))
  # exporting functions from layeredBB package to simulate layered Brownian bridges
  parallel::clusterExport(cl, varlist = ls("package:layeredBB"))
  if (!is.null(seed)) {
    parallel::clusterSetRNGStream(cl, iseed = seed)
  }
  # split the x samples and their means into approximately equal lists
  max_samples_per_core <- ceiling(N/n_cores)
  split_indices <- split(1:N, ceiling(seq_along(1:N)/max_samples_per_core))
  split_x_samples <- lapply(split_indices, function(indices) particle_set$x_samples[indices])
  split_x_means <- lapply(split_indices, function(indices) particle_set$x_means[indices])
  # for each set of x samples, we propose a new value y and assign a weight for it
  # sample for y and importance weight in parallel to split computation
  Q_weighted_samples <- parallel::parLapply(cl, X = 1:length(split_indices), fun = function(core) {
    split_N <- length(split_indices[[core]])
    y_samples <- rep(0, split_N)
    log_Q_weights <- rep(0, split_N)
    for (i in 1:split_N) {
      y_samples[i] <- rnorm(1, mean = split_x_means[[core]][i], sd = proposal_sd)
      log_Q_weights[i] <- sum(sapply(1:m, function(c) {
        ea_mixG_DL_PT(x0 = split_x_samples[[core]][[i]][c],
                      y = y_samples[i],
                      s = 0,
                      t = time,
                      n_comp = n_comp,
                      weights = weights,
                      means = means,
                      sds = sds,
                      beta = betas[c],
                      precondition = precondition_values[c],
                      bounds_multiplier = bounds_multiplier,
                      logarithm = TRUE)
      }))
    }
    return(list('y_samples' = y_samples, 'log_Q_weights' = log_Q_weights))
  })
  parallel::stopCluster(cl)
  # unlist the proposed samples for y and their associated log Q weights
  y_samples <- unlist(lapply(1:length(split_x_samples), function(i) {
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

#' SMC Fusion for sub-posteriors of the form exp(-((x^4)*beta)/2)
#'
#' SMC Fusion for sub-posteriors of the form exp(-((x^4)*beta)/2)
#'
#' @param particles_to_fuse list of length m, where particles_to_fuse[[c]]
#'                          contains the particles for the c-th sub-posterior
#'                          (a list of particles to fuse can be initialised by 
#'                          initialise_particle_sets() function)
#' @param N number of samples
#' @param m number of sub-posteriors to combine
#' @param time time T for fusion algorithm
#' @param n_comp integer number of components of mixture Gaussian
#' @param weights vector: weights of mixture Gaussian
#' @param means vector: means of mixture Gaussian
#' @param sds vector: st.devs of mixture Gaussian
#' @param betas vector of length c, where betas[c] is the inverse temperature 
#'              value for c-th posterior
#' @param precondition_values vector of length m, where precondition_values[c]
#'                            is the precondition value for sub-posterior c
#' @param bounds_multiplier scalar value to mulitply bounds by 
#'                          (should greater than or equal to 1)
#' @param resampling_method method to be used in resampling, default is multinomial 
#'                          resampling ('multi'). Other choices are stratified 
#'                          resampling ('strat'), systematic resampling ('system'),
#'                          residual resampling ('resid')
#' @param ESS_threshold number between 0 and 1 defining the proportion 
#'                      of the number of samples that ESS needs to be
#'                      lower than for resampling (i.e. resampling is carried 
#'                      out only when ESS < N*ESS_threshold)
#' @param seed seed number - default is NULL, meaning there is no seed
#' @param n_cores number of cores to use
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
#'   \item{precondition_values}{list of length 2 where precondition_values[[2]] 
#'                              are the pre-conditioning values that were used 
#'                              and precondition_values[[1]] are the combined 
#'                              precondition values}
#' }
#' 
#' @export
parallel_fusion_SMC_mixG <- function(particles_to_fuse,
                                     N, 
                                     m,
                                     time,
                                     n_comp,
                                     weights,
                                     means,
                                     sds,
                                     betas,
                                     precondition_values,
                                     bounds_multiplier = 1.1,
                                     resampling_method = 'multi',
                                     ESS_threshold = 0.5,
                                     seed = NULL,
                                     n_cores = parallel::detectCores()) {
  if (!is.list(particles_to_fuse) | (length(particles_to_fuse)!=m)) {
    stop("parallel_fusion_SMC_mixG: particles_to_fuse must be a list of length m")
  } else if (!all(sapply(particles_to_fuse, function(sub_posterior) ("particle" %in% class(sub_posterior))))) {
    stop("parallel_fusion_SMC_mixG: particles in particles_to_fuse must be \"particle\" objects")
  } else if (length(weights)!=n_comp) {
    stop("parallel_fusion_SMC_mixG: weights must be a vector of length n_comp")
  } else if (length(means)!=n_comp) {
    stop("parallel_fusion_SMC_mixG: means must be a vector of length n_comp")
  } else if (length(sds)!=n_comp) {
    stop("parallel_fusion_SMC_mixG: sds must be a vector of length n_comp")
  } else if (!is.vector(betas) | (length(betas)!=m)) {
    stop("parallel_fusion_SMC_mixG: betas must be a vector of length m")
  } else if (!is.vector(precondition_values) | (length(precondition_values)!=m)) {
    stop("parallel_fusion_SMC_mixG: precondition_values must be a vector of length m")
  } else if ((ESS_threshold < 0) | (ESS_threshold > 1)) {
    stop("parallel_fusion_SMC_mixG: ESS_threshold must be between 0 and 1")
  }
  # set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  # start time recording
  pcm <- proc.time()
  # ---------- first importance sampling step 
  particles <- rho_IS_univariate(particles_to_fuse = particles_to_fuse,
                                 N = N,
                                 m = m,
                                 time = time,
                                 precondition_values = precondition_values,
                                 resampling_method = resampling_method,
                                 n_cores = n_cores)
  # record ESS and CESS after rho step 
  ESS <- c('rho' = particles$ESS)
  CESS <- c('rho' = particles$CESS['rho'])
  # ----------- resample particles
  # only resample if ESS < N*ESS_threshold
  if (particles$ESS < N*ESS_threshold) {
    resampled <- c('rho' = TRUE)
    particles <- resample_particle_x_samples(N = N,
                                             particle_set = particles,
                                             multivariate = FALSE,
                                             resampling_method = resampling_method,
                                             seed = seed)
  } else {
    resampled <- c('rho' = FALSE)
  }
  # ---------- second importance sampling step
  # unbiased estimator for Q
  particles <- Q_IS_mixG(particle_set = particles,
                         m = m,
                         time = time,
                         n_comp = n_comp,
                         weights = weights,
                         means = means,
                         sds = sds,
                         betas = betas,
                         precondition_values = precondition_values,
                         bounds_multiplier = bounds_multiplier,
                         seed = seed,
                         n_cores = n_cores)
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
                                             multivariate = FALSE,
                                             resampling_method = resampling_method,
                                             seed = seed)
  } else {
    resampled['Q'] <- FALSE
  }
  if (identical(precondition_values, rep(1, m))) {
    return(list('particles' = particles,
                'proposed_samples' = proposed_samples,
                'time' = (proc.time()-pcm)['elapsed'],
                'ESS' = ESS,
                'CESS' = CESS,
                'resampled' = resampled,
                'precondition_values' = list(1, precondition_values)))
  } else {
    return(list('particles' = particles,
                'proposed_samples' = proposed_samples,
                'time' = (proc.time()-pcm)['elapsed'],
                'ESS' = ESS,
                'CESS' = CESS,
                'resampled' = resampled,
                'precondition_values' = list(1/sum(1/precondition_values), 
                                             precondition_values)))
  }
}

#' Hierarchical SMC Fusion 
#'
#' Hierarchical SMC Fusion with base level of the form exp(-((x^4)*beta)/2)
#'
#' @param N_schedule vector of length (L-1), where N_schedule[l] is the number 
#'                   of samples per node at level l
#' @param m_schedule vector of length (L-1), where m_schedule[l] is the number 
#'                   of samples to fuse for level l
#' @param time_schedule vector of legnth(L-1), where time_schedule[l] is the time 
#'                      chosen for Fusion at level l
#' @param base_samples list of length (1/start_beta), where samples_to_fuse[c] 
#'                     contains the samples for the c-th node in the level
#' @param L total number of levels in the hierarchy
#' @param n_comp integer number of components of mixture Gaussian
#' @param weights vector: weights of mixture Gaussian
#' @param means vector: means of mixture Gaussian
#' @param sds vector: st.devs of mixture Gaussian
#' @param start_beta beta for the base level
#' @param precondition either a logical value to determine if preconditioning values are
#'                     used (TRUE - and is set to be the variance of the sub-posterior samples)
#'                     or not (FALSE - and is set to be 1 for all sub-posteriors),
#'                     or a list of length (1/start_beta) where precondition[[c]]
#'                     is the preconditioning value for sub-posterior c. Default is TRUE
#' @param bounds_multiplier scalar value to mulitply bounds by 
#'                          (should greater than or equal to 1)
#' @param resampling_method method to be used in resampling, default is multinomial 
#'                          resampling ('multi'). Other choices are stratified 
#'                          resampling ('strat'), systematic resampling ('system'),
#'                          residual resampling ('resid')
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
#'                   particles for level l, node i}
#'   \item{proposed_samples}{list of length (L-1), where proposed_samples[[l]][[i]]
#'                          are the proposed samples for level l, node i}
#'   \item{time}{list of length (L-1), where time[[l]] is the run time for level l,
#'              node i}
#'   \item{ESS}{list of length (L-1), where ESS[[l]][[i]] is the effective 
#'             sample size of the particles after each step BEFORE deciding 
#'             whether or not to resample for level l, node i}
#'   \item{CESS}{list of length (L-1), where CESS[[l]][[i]] is the conditional
#'              effective sample size of the particles after each step}
#'   \item{resampled}{list of length (L-1), where resampled[[l]][[i]] is a 
#'                   boolean value to record if the particles were resampled
#'                   after each step; rho and Q for level l, node i}
#'   \item{precondition_values}{preconditioning values used in the algorithm 
#'                             for each node}
#'   \item{diffusion_times}{vector of length (L-1), where diffusion_times[l]
#'                         are the times for fusion in level l}
#' }
#'
#' @export
hierarchical_fusion_SMC_mixG <- function(N_schedule,
                                         m_schedule,
                                         time_schedule,
                                         base_samples,
                                         L,
                                         n_comp,
                                         weights,
                                         means,
                                         sds,
                                         start_beta,
                                         precondition = TRUE,
                                         bounds_multiplier = 1.1,
                                         resampling_method = 'multi',
                                         ESS_threshold = 0.5,
                                         seed = NULL,
                                         n_cores = parallel::detectCores()) {
  if (length(weights)!=n_comp) {
    stop("hierarchical_fusion_SMC_mixG: weights must be a vector of length n_comp")
  } else if (length(means)!=n_comp) {
    stop("hierarchical_fusion_SMC_mixG: means must be a vector of length n_comp")
  } else if (length(sds)!=n_comp) {
    stop("hierarchical_fusion_SMC_mixG: sds must be a vector of length n_comp")
  } else if (!is.vector(N_schedule) | (length(N_schedule)!=(L-1))) {
    stop("hierarchical_fusion_SMC_mixG: N_schedule must be a vector of length (L-1)")
  } else if (!is.vector(m_schedule) | (length(m_schedule)!=(L-1))) {
    stop("hierarchical_fusion_SMC_mixG: m_schedule must be a vector of length (L-1)")
  } else if (!is.vector(time_schedule) | (length(time_schedule)!=(L-1))) {
    stop("hierarchical_fusion_SMC_mixG: time_schedule must be a vector of length (L-1)")
  } else if (!is.list(base_samples) | (length(base_samples)!=(1/start_beta))) {
    stop("hierarchical_fusion_SMC_mixG: base_samples must be a list of length (1/start_beta)")
  } else if ((ESS_threshold < 0) | (ESS_threshold > 1)) {
    stop("hierarchical_fusion_SMC_mixG: ESS_threshold must be between 0 and 1")
  }
  if (is.vector(m_schedule) & (length(m_schedule)==(L-1))) {
    for (l in (L-1):1) {
      if (((1/start_beta)/prod(m_schedule[(L-1):l]))%%1!=0) {
        stop("hierarchical_fusion_SMC_mixG: check that (1/start_beta)/prod(m_schedule[(L-1):l])
             is an integer for l=L-1,...,1")
      }
    }
  } else {
    stop("hierarchical_fusion_SMC_mixG: m_schedule must be a vector of length (L-1)")
  }
  # we append 1 to the vector m_schedule to make the indices work later on when we call fusion
  m_schedule <- c(m_schedule, 1)
  # initialising results
  particles <- list()
  if (all(sapply(base_samples, function(sub) class(sub)=='particle'))) {
    particles[[L]] <- base_samples
  } else if (all(sapply(base_samples, is.vector))) {
    particles[[L]] <- initialise_particle_sets(samples_to_fuse = base_samples, multivariate = FALSE)
  } else {
    stop("hierarchical_fusion_SMC_mixG: base_samples must be a list of length 
         (1/start_beta) containing either items of class \"particle\" (representing
         particle approximations of the sub-posteriors) or are vectors
         (representing un-normalised sample approximations of the sub-posteriors)")
  }
  proposed_samples <- list()
  time <- list()
  ESS <- list()
  CESS <- list()
  resampled <- list()
  precondition_values <- list()
  if (is.logical(precondition)) {
    if (precondition) {
      precondition_values[[L]] <- lapply(base_samples, var)
    } else {
      precondition_values[[L]] <- lapply(1:length(base_samples), function(i) 1)
    }
  } else if (is.list(precondition)) {
    if (length(precondition)==(1/start_beta)) {
      precondition_values[[L]] <- precondition
    }
  } else {
    stop("hierarchical_fusion_SMC_mixG: precondition must be a logical indicating 
          whether or not a preconditioning value should be used, or a list of
          length C, where precondition[[c]] is the preconditioning value for
          the c-th sub-posterior")
  }
  cat('Starting hierarchical fusion \n', file = 'hierarchical_fusion_SMC_mixG.txt')
  for (k in ((L-1):1)) {
    # since previous level has (1/beta)/prod(m_schedule[L:(k-1)]) nodes and we 
    # fuse m_schedule[k] of these
    n_nodes <- max((1/start_beta)/prod(m_schedule[L:k]), 1)
    cat('########################\n', file = 'hierarchical_fusion_SMC_mixG.txt',
        append = T)
    cat('Starting to fuse', m_schedule[k], 'densities of pi^beta, where beta =',
        prod(m_schedule[L:(k+1)]), '/', (1/start_beta), 'for level', k, 'with time',
        time_schedule[k], ', which is using', parallel::detectCores(), 'cores\n',
        file = 'hierarchical_fusion_SMC_mixG.txt', append = T)
    cat('There are', n_nodes, 'nodes at this level each giving', N_schedule[k],
        'samples for beta =', prod(m_schedule[L:k]), '/', (1/start_beta),
        '\n', file = 'hierarchical_fusion_SMC_mixG.txt', append = T)
    cat('########################\n', file = 'hierarchical_fusion_SMC_mixG.txt', 
        append = T)
    fused <- lapply(X = 1:n_nodes, FUN = function(i) {
      previous_nodes <- ((m_schedule[k]*i)-(m_schedule[k]-1)):(m_schedule[k]*i)
      precondition_vals <- unlist(precondition_values[[k+1]][previous_nodes])
      parallel_fusion_SMC_mixG(particles_to_fuse = particles[[k+1]][previous_nodes],
                               N = N_schedule[k],
                               m = m_schedule[k],
                               time = time_schedule[k],
                               n_comp = n_comp,
                               weights = weights,
                               means = means,
                               sds = sds,
                               betas = prod(m_schedule[L:(k+1)])*(start_beta),
                               precondition_values = precondition_vals, 
                               bounds_multiplier = bounds_multiplier,
                               resampling_method = resampling_method,
                               ESS_threshold = ESS_threshold,
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
    precondition_values[[k]] <- lapply(1:n_nodes, function(i) fused[[i]]$precondition_values[[1]])
  }
  cat('Completed hierarchical fusion\n', file = 'hierarchical_fusion_SMC_mixG.txt', 
      append = T)
  if (length(particles[[1]])==1) {
    particles[[1]] <- particles[[1]][[1]]
    proposed_samples[[1]] <- proposed_samples[[1]][[1]]
    time[[1]] <- time[[1]][[1]]
    ESS[[1]] <- ESS[[1]][[1]]
    CESS[[1]] <- CESS[[1]][[1]]
    resampled[[1]] <- resampled[[1]][[1]]
    precondition_values[[1]] <- precondition_values[[1]][[1]]
  }
  return(list('particles' = particles,
              'proposed_samples' = proposed_samples,
              'time' = time,
              'ESS' = ESS,
              'CESS' = CESS,
              'resampled' = resampled,
              'precondition_values' = precondition_values,
              'diffusion_times' = time_schedule))
}

#' Progressive SMC Fusion
#'
#' Progressive SMC Fusion with base level of the form exp(-((x^4)*beta)/2)
#'
#' @param N_schedule vector of length (L-1), where N_schedule[l] is the number 
#'                   of samples per node at level l
#' @param time_schedule vector of legnth(L-1), where time_schedule[l] is the time 
#'                      chosen for Fusion at level l
#' @param base_samples list of length (1/start_beta), where samples_to_fuse[c] 
#'                     contains the samples for the c-th node in the level
#' @param n_comp integer number of components of mixture Gaussian
#' @param weights vector: weights of mixture Gaussian
#' @param means vector: means of mixture Gaussian
#' @param sds vector: st.devs of mixture Gaussian
#' @param start_beta beta for the base level
#' @param precondition either a logical value to determine if preconditioning values are
#'                     used (TRUE - and is set to be the variance of the sub-posterior samples)
#'                     or not (FALSE - and is set to be 1 for all sub-posteriors),
#'                     or a list of length (1/start_beta) where precondition[[c]]
#'                     is the preconditioning value for sub-posterior c. Default is TRUE
#' @param bounds_multiplier scalar value to mulitply bounds by 
#'                          (should greater than or equal to 1)
#' @param resampling_method method to be used in resampling, default is multinomial 
#'                          resampling ('multi'). Other choices are stratified 
#'                          resampling ('strat'), systematic resampling ('system'),
#'                          residual resampling ('resid')
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
#'                   particles for level l, node i}
#'   \item{proposed_samples}{list of length (L-1), where proposed_samples[[l]][[i]]
#'                          are the proposed samples for level l, node i}
#'   \item{time}{list of length (L-1), where time[[l]] is the run time for level l,
#'              node i}
#'   \item{ESS}{list of length (L-1), where ESS[[l]][[i]] is the effective 
#'             sample size of the particles after each step BEFORE deciding 
#'             whether or not to resample for level l, node i}
#'   \item{CESS}{list of length (L-1), where CESS[[l]][[i]] is the conditional
#'              effective sample size of the particles after each step}
#'   \item{resampled}{list of length (L-1), where resampled[[l]][[i]] is a 
#'                   boolean value to record if the particles were resampled
#'                   after each step; rho and Q for level l, node i}
#'   \item{precondition_values}{preconditioning values used in the algorithm 
#'                             for each node}
#'   \item{diffusion_times}{vector of length (L-1), where diffusion_times[l]
#'                         are the times for fusion in level l}
#' }
#'
#' @export
progressive_fusion_SMC_mixG <- function(N_schedule,
                                        time_schedule,
                                        base_samples,
                                        n_comp,
                                        weights,
                                        means,
                                        sds,
                                        start_beta,
                                        precondition = TRUE,
                                        bounds_multiplier = 1.1,
                                        resampling_method = 'multi',
                                        ESS_threshold = 0.5,
                                        seed = NULL,
                                        n_cores = parallel::detectCores()) {
  if (length(weights)!=n_comp) {
    stop("progressive_fusion_SMC_mixG: weights must be a vector of length n_comp")
  } else if (length(means)!=n_comp) {
    stop("progressive_fusion_SMC_mixG: means must be a vector of length n_comp")
  } else if (length(sds)!=n_comp) {
    stop("progressive_fusion_SMC_mixG: sds must be a vector of length n_comp")
  } else if (!is.vector(N_schedule) | (length(N_schedule)!=(1/start_beta)-1)) {
    stop("progressive_fusion_SMC_mixG: N_schedule must be a vector of length ((1/start_beta)-1)")
  } else if (!is.vector(time_schedule) | (length(time_schedule)!=(1/start_beta)-1)) {
    stop("progressive_fusion_SMC_mixG: time_schedule must be a vector of length ((1/start_beta)-1")
  } else if (!is.list(base_samples) | (length(base_samples)!=(1/start_beta))) {
    stop("progressive_fusion_SMC_mixG: base_samples must be a list of length (1/start_beta)")
  } else if (ESS_threshold < 0 | ESS_threshold > 1) {
    stop("progressive_fusion_SMC_mixG: ESS_threshold must be between 0 and 1")
  }
  # initialising results
  particles <- list()
  if (all(sapply(base_samples, function(sub) class(sub)=='particle'))) {
    particles[[(1/start_beta)]] <- base_samples
  } else if (all(sapply(base_samples, is.vector))) {
    particles[[(1/start_beta)]] <- initialise_particle_sets(samples_to_fuse = base_samples, multivariate = FALSE)
  } else {
    stop("progressive_fusion_SMC_mixG: base_samples must be a list of length
         (1/start_beta) containing either items of class \"particle\" (representing
         particle approximations of the sub-posteriors) or are vectors (representing
         un-normalised sample approximations of the sub-posteriors)")
  }
  proposed_samples <- list()
  time <- list()
  ESS <- list()
  CESS <- list()
  resampled <- list()
  precondition_values <- list()
  if (is.logical(precondition)) {
    if (precondition) {
      precondition_values[[(1/start_beta)]] <- lapply(base_samples, var)
    } else {
      precondition_values[[(1/start_beta)]] <- lapply(1:length(base_samples), function(i) 1)
    }
  } else if (is.list(precondition)) {
    if (length(precondition)==(1/start_beta)) {
      precondition_values[[(1/start_beta)]] <- precondition
    }
  } else {
    stop("progressive_fusion_SMC_mixG: precondition must be a logical indicating 
          whether or not a preconditioning value should be used, or a list of
          length C, where precondition[[c]] is the preconditioning value for
          the c-th sub-posterior")
  }
  index <- 2
  cat('Starting progressive fusion \n', file = 'progressive_fusion_SMC_mixG.txt')
  for (k in ((1/start_beta)-1):1) {
    if (k==(1/start_beta)-1) {
      cat('########################\n', file = 'progressive_fusion_SMC_mixG.txt', 
          append = T)
      cat('Starting to fuse', 2, 'densities for level', k, 'which is using', 
          parallel::detectCores(), 'cores\n', 
          file = 'progressive_fusion_SMC_mixG.txt', append = T)
      cat('Fusing samples for beta =', 1, '/', (1/start_beta), 'with time',
          time_schedule[k], 'to get', N_schedule[k], 'samples for beta =', 2, 
          '/', (1/start_beta), '\n', file = 'progressive_fusion_SMC_mixG.txt', 
          append = T)
      cat('########################\n', file = 'progressive_fusion_SMC_mixG.txt',
          append = T)
      particles_to_fuse <- list(particles[[(1/start_beta)]][[1]], 
                                particles[[(1/start_beta)]][[2]])
      precondition_vals <- unlist(precondition_values[[k+1]][1:2])
      fused <- parallel_fusion_SMC_mixG(particles_to_fuse = particles_to_fuse,
                                        N = N_schedule[k],
                                        m = 2,
                                        time = time_schedule[k],
                                        n_comp = n_comp,
                                        weights = weights,
                                        means = means,
                                        sds = sds,
                                        betas = c(start_beta, start_beta),
                                        precondition_values = precondition_vals,
                                        bounds_multiplier = bounds_multiplier,
                                        resampling_method = resampling_method,
                                        ESS_threshold = ESS_threshold,
                                        seed = seed,
                                        n_cores = n_cores)
    } else {
      cat('########################\n', file = 'progressive_fusion_SMC_mixG.txt', 
          append = T)
      cat('Starting to fuse', 2, 'densities for level', k, 'which is using', 
          parallel::detectCores(), 'cores\n',
          file = 'progressive_fusion_SMC_mixG.txt', append = T)
      cat('Fusing samples for beta =', index, '/', (1/start_beta), 'and beta =', 
          1, '/', (1/start_beta), 'with time', time_schedule[k], 'to get',
          N_schedule[k], 'samples for beta =', (index+1), '/', (1/start_beta),
          '\n', file = 'progressive_fusion_SMC_mixG.txt', append = T)
      cat('########################\n', file = 'progressive_fusion_SMC_mixG.txt', 
          append = T)
      particles_to_fuse <- list(particles[[k+1]], 
                                particles[[(1/start_beta)]][[index+1]])
      precondition_vals <- c(precondition_values[[k+1]],
                             precondition_values[[(1/start_beta)]][[index+1]]) 
      fused <- parallel_fusion_SMC_mixG(particles_to_fuse = particles_to_fuse, 
                                        N = N_schedule[k],
                                        m = 2,
                                        time = time_schedule[k],
                                        n_comp = n_comp,
                                        weights = weights,
                                        means = means,
                                        sds = sds,
                                        betas = c(index*start_beta, start_beta),
                                        precondition_values = precondition_vals,
                                        bounds_multiplier = bounds_multiplier,
                                        resampling_method = resampling_method,
                                        ESS_threshold = ESS_threshold,
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
    precondition_values[[k]] <- fused$precondition_values[[1]]
  }
  cat('Completed progressive fusion\n', 
      file = 'progressive_fusion_SMC_mixG.txt', append = T)
  return(list('particles' = particles,
              'proposed_samples' = proposed_samples,
              'time' = time,
              'ESS' = ESS,
              'CESS' = CESS,
              'resampled' = resampled,
              'precondition_values' = precondition_values,
              'diffusion_times' = time_schedule))
}
