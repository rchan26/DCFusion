#' Initialise particle sets from a list of samples
#'
#' Function to intialise particle sets from a list of samples
#'
#' @param samples_to_fuse a list of samples that you wish to perform fusion with
#' @param multivariate logical value indicating if the samples are multivariate
#'                     (TRUE) or not (FALSE)
#'
#' @return A list with components
#' \describe{
#'   \item{y_samples}{samples for y in particle set (initialised as the samples given)}
#'   \item{x_samples}{a list where x_samples[[i]] is the ith x sample for in 
#'                    particle set (initialised as NA)}
#'   \item{x_mean}{the corresponding means for x_samples (initialised as NA)}
#'   \item{log_weights}{associated logarithm of the weights (initiliased as
#'                      the logarithm of 1 / number of samples)}
#'   \item{normalised_weights}{associated normalised weights (initialised as 
#'                             1 / number of samples)}
#'   \item{ESS}{effective sample size of particles (initialised as the 
#'              number of samples)}
#'   \item{CESS}{conditional effective sample size of particles (initialised 
#'               as the number of samples)}
#'   \item{resampled}{logical value to idicate if particles have been resampled
#'                    (initialised as TRUE)}
#'   \item{N}{Number of particles}
#' }
#'
#' @examples 
#' samples <- rnorm(100, 0, 1)
#' initialise_particle_sets(samples_to_fuse = samples, multivariate = FALSE)
#' 
#' @export
initialise_particle_sets <- function(samples_to_fuse, multivariate) {
  # create a list for the particles where y_samples are the sub-posterior samples
  # initialise weights
  if (multivariate) {
    return(lapply(samples_to_fuse, function(sub_posterior) {
      N <- nrow(sub_posterior)
      dim <- ncol(sub_posterior)
      return(list('y_samples' = sub_posterior, 
                  'x_samples' = rep(list(NA), N),
                  'x_means' = matrix(nrow = N, ncol = dim),
                  'log_weights' = log(rep(1/N, N)),
                  'normalised_weights' = rep(1/N, N), 
                  'ESS' = N,
                  'CESS' = NA,
                  'resampled' = TRUE,
                  'N' = N))
    }))
  } else {
    return(lapply(samples_to_fuse, function(sub_posterior) {
      N <- length(sub_posterior)
      return(list('y_samples' = sub_posterior, 
                  'x_samples' = rep(list(NA), N),
                  'x_means' = rep(NA, N),
                  'log_weights' = log(rep(1/N, N)),
                  'normalised_weights' = rep(1/N, N), 
                  'ESS' = N,
                  'CESS' = NA,
                  'resampled' = TRUE,
                  'N' = N))
    }))
  }
}

#' Resampling a particle set
#'
#' @param particle_set set of samples for x and y with corresponding weights
#' @param multivariate logical value indicating if the particles are multivariate
#'                     (TRUE) or not (FALSE)
#' @param resampling_method method to be used in resampling, default is multinomial 
#'                          resampling ('multi'). Other choices are stratified 
#'                          resampling ('strat'), systematic resampling ('system'),
#'                          residual resampling ('resid')
#' @param seed seed number - default is NULL, meaning there is no seed
#'
#' @return y parameter samples from resampled particle set
#' 
#' @export
resample_particle_set <- function(particle_set, 
                                  multivariate = TRUE, 
                                  resampling_method = 'multi',
                                  seed = NULL) {
  # set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  # check if the particles have not already been re-sampled
  if (particle_set$resampled) {
    return(particle_set$y_samples)
  } else {
    indices <- resample_indices(normalised_weights = particle_set$normalised_weights,
                                method = resampling_method,
                                n = particle_set$N)
    if (multivariate) {
      return(particle_set$y_samples[indices,])
    } else {
      return(particle_set$y_samples[indices])
    }
  }
}

# ---------- functions for resampling schemes from Murray Pollock Github

multi.resamp <- function(normalised_weights, 
                         n = length(normalised_weights)) { 
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

system.resamp <- function(normalised_weights, 
                          n = length(normalised_weights)) { 
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

strat.resamp <- function(normalised_weights, 
                         n = length(normalised_weights)) { 
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
