#' Tempered target density exp(-(((x-mean)^4)*beta)/2)
#'
#' Normalised tempered target density
#'
#' @param x real value
#' @param mean mean value, defaults to 0
#' @param beta real value between 0 and 1, defaults to 1
#'
#' @return real value
#'
#' @examples
#' curve(exp_4_density(x, 1/2), from = -4, to = 4)
#'
#' @export
exp_4_density <- function(x, mean = 0, beta = 1) {
  norm_constant <- integrate(function(y) exp(-(((y-mean)^4)*beta)/2), 
                             lower = -Inf, upper = Inf)$value
  return(exp(-(((x-mean)^4)*beta)/2) / norm_constant)
}

#' Rejection sampler for tempered target exp(-(beta*(x-mean)^4)/2)
#'
#' Sample from tempered target using rejection sampling with normal proposals
#'
#' @param N number of samples
#' @param mean mean value, defaults to 0
#' @param proposal_mean mean of proposal density for rejection sampler
#' @param proposal_sd st.dev of proposal density for rejection sampler
#' @param dominiating_M constant M to bound the target density
#' @param beta real value between 0 and 1
#'
#' @return samples from tempered target density
#'
#' @examples
#' dominating_dnorm <- function(x) 1.3*dnorm(x, mean = 0, sd = 1)
#' curve(dominating_dnorm(x), -5, 5, col = 'red')
#' curve(dnorm(x, mean = 0, sd = 1), -5, 5, add = T, col = 'green')
#' curve(exp_4_density(x, 1/2), -5, 5, add = T, col = 'blue')
#' sample_from_fc(N = 100000, proposal_mean = 0, proposal_sd = 1, dominating_M = 1.3)
#'
#' @export
sample_exp_4 <- function(N, 
                         mean = 0, 
                         proposal_mean, 
                         proposal_sd,
                         dominating_M, 
                         beta) {
  dominating_function <- function(x) {
    return(dominating_M*dnorm(x, mean = proposal_mean, sd = proposal_sd))
    
  }
  target_fc <- function(x) {
    norm_constant <- integrate(function(y) exp(-(((y-mean)^4)*beta)/2),
                               lower = -Inf, upper = Inf)$value
    return(exp(-(((x-mean)^4)*beta)/2) / norm_constant)
  }
  i <- 0
  samples <- rep(NA, N)
  iters <- 0
  while (i < N) {
    iters <- iters+1
    prop <- rnorm(1, mean = proposal_mean, sd = proposal_sd)
    if (runif(1,0,1) < (target_fc(prop) / dominating_function(prop))) {
      i <- i+1
      samples[i] <- prop
    }
  }
  print(paste('acceptance rate:', N, '/', iters, '=', N/iters))
  return(samples)
}

#' Rejection sampler for base level
#'
#' Sample for base level (samples of exp(-(beta*(x-mean)^4)/2))
#'
#' @param beta real value between 0 and 1
#' @param mean mean value, defaults to 0
#' @param nsamples number of samples per node
#' @param proposal_mean mean of proposal density for rejection sampler
#' @param proposal_sd st.dev of proposal density for rejection sampler
#' @param dominiating_M constant M to bound the target density
#'
#' @return return_name return description
#'
#' @examples
#' base_rejection_sampler_exp_4(beta = 1/2, 
#'                              nsamples = 100000, 
#'                              proposal_mean = 0, 
#'                              proposal_sd = 1, 
#'                              dominating_M = 1.3)
#'                              
#' @export
base_rejection_sampler_exp_4 <- function(beta, 
                                         mean = 0, 
                                         nsamples, 
                                         proposal_mean, 
                                         proposal_sd, 
                                         dominating_M) {
  print(paste('sampling from tempered target density with beta =', 1, '/', (1/beta)))
  return(lapply(1:(1/beta), function(i) {
    sample_exp_4(N = nsamples,
                 mean = mean,
                 proposal_mean = proposal_mean,
                 proposal_sd = proposal_sd,
                 dominating_M = dominating_M,
                 beta = beta)}))
}
