check_equal_dim <- function(list_of_vectors) {
  return(all(sapply(list_of_vectors[-1], function(x) length(x)==length(list_of_vectors[[1]]))))
}

#' Density of a mixture Gaussian
#'
#' Returns the value of a mixture Gaussian at given point x
#'
#' @param x real value
#' @param n_comp integer number of components of mixture Gaussian
#' @param weights vector: weights of mixture Gaussian
#' @param means vector: means of mixture Gassuan
#' @param sds vector: st.devs of mixture Gaussian
#'
#' @return value: value of a mixture Gaussian at given point x
#'
#' @examples
#' weights <- c(0.4, 0.6)
#' means <- c(-3, 6)
#' sds <- c(1 ,2)
#' curve(dnorm_mix(x,
#'                 n_comp = n_comp,
#'                 weights = weights,
#'                 means = means,
#'                 sds = sds),
#'       -10, 15,
#'       ylim = c(0, 0.2), ylab = 'pdf')
#'
#' @export
dnorm_mix <- function(x, n_comp, weights, means, sds) {
  if (length(weights)!=n_comp) {
    stop("dnorm_mix: weights must be a vector of length n_comp")
  } else if (length(means)!=n_comp) {
    stop("dnorm_mix: means must be a vector of length n_comp")
  } else if (length(sds)!=n_comp) {
    stop("dnorm_mix: sds must be a vector of length n_comp")
  }
  value <- rep(0, length(x))
  for (i in 1:length(weights)) {
    value <- value + weights[i]*dnorm(x=x, mean=means[i], sd=sds[i])
  }
  return(value)
}

#' Density of tempered mixture Gaussian
#'
#' Returns the value of a normalised tempered mixture Gaussian at given point x
#'
#' @param x real value
#' @param n_comp integer number of components of mixture Gaussian
#' @param weights vector: weights of mixture Gaussian
#' @param means vector: means of mixture Gassuan
#' @param sds vector: st.devs of mixture Gaussian
#' @param beta real value
#' @param normalised boolean value to determine if normalisation constant is calculated
#'
#' @return value: value of a normalised tempered mixture Gaussian at given point x
#'
#' @examples
#' weights <- c(0.4, 0.6)
#' means <- c(-3, 6)
#' sds <- c(1, 2)
#' beta <- 1/4
#' curve(dnorm_mix_tempered(x,
#'                          n_comp = 2,
#'                          weights = weights,
#'                          means = means,
#'                          sds = sds,
#'                          beta = beta,
#'                          normalised = FALSE),
#'       -15, 20, ylim = c(0, 0.8), ylab = 'pdf')
#' curve(dnorm_mix_tempered(x,
#'                          n_comp = 2,
#'                          weights = weights,
#'                          means = means,
#'                          sds = sds,
#'                          beta = beta,
#'                          normalised = TRUE),
#'       add = T, lty = 2)
#'
#' @export
dnorm_mix_tempered <- function(x, n_comp, weights, means, sds, beta, normalised = TRUE) {
  if (length(weights)!=n_comp) {
    stop("dnorm_mix_tempered: weights must be a vector of length n_comp")
  } else if (length(means)!=n_comp) {
    stop("dnorm_mix_tempered: means must be a vector of length n_comp")
  } else if (length(sds)!=n_comp) {
    stop("dnorm_mix_tempered: sds must be a vector of length n_comp")
  }
  if (normalised) {
    Z <- integrate(function(x) {
      dnorm_mix_tempered_unnormalised(x=x, w=weights, m=means, s=sds, b=beta)
    },
    lower = -Inf, upper = Inf, subdivisions = 500L, rel.tol = .Machine$double.eps^0.8)$value
    return(dnorm_mix_tempered_unnormalised(x=x, w=weights, m=means, s=sds, b=beta) / Z)
  } else {
    return(dnorm_mix_tempered_unnormalised(x=x, w=weights, m=means, s=sds, b=beta))
  }
}

#' Simulate exactly a mixture Gaussian distribution
#'
#' Simulate exactly a mixture Gaussian distribution using component sampling
#'
#' @param n number of samples required
#' @param weights vector: weights of mixture Gaussian
#' @param means vector: means of mixture Gassuan
#' @param sds vector: st.devs of mixture Gaussian
#'
#' @return samples from a mixture Gaussian distribution
#'
#' @examples
#' rmix_norm(N = 10000, weights = c(0.4, 0.6), means = c(-3, 9), sds = c(2, 0.4))
#'
#' @export
rmix_norm <- function(N, weights, means, sds) {
  if (sum(weights)!=1) {
    stop("rmix_norm: weights do not add up to 1")
  } else if (!check_equal_dim(list(weights, means, sds))) {
    stop("rmix_norm: dimensions of input vectors do not match")
  }
  w <- sample(1:length(weights), prob=weights, size=N, replace=T)
  samples <- rnorm(n = N, mean=means[w], sd=sds[w])
  return(samples)
}

#' Rejection sampler for tempered mixture Gaussian
#'
#' Sample from tempered target using rejection sampling with mixture Gaussian proposals
#'
#' @param N number of samples
#' @param weights vector: weights of mixture Gaussian
#' @param means vector: means of mixture Gassuan
#' @param sds vector: st.devs of mixture Gaussian
#' @param beta real value between 0 and 1
#' @param proposal_sds vector: st.devs of proposal mixture Gaussian
#' @param dominating_M constant M to bound the target density
#'
#' @return samples from tempered target density
#'
#' @export
rmix_norm_tempered_rejection_sampler <- function(N,
                                                 weights,
                                                 means,
                                                 sds,
                                                 beta,
                                                 proposal_sds,
                                                 dominating_M) {
  if (sum(weights)!=1) {
    stop("rmix_norm_tempered_rejection_sampler: weights do not add up to 1")
  } else if (!check_equal_dim(list(weights, means, sds))) {
    stop("rmix_norm_tempered_rejection_sampler: dimensions of input vectors do not match")
  }
  dominating_function <- function(x) {
    return(dominating_M*dnorm_mix(x = x, 
                                  n_comp = length(weights),
                                  weights = weights,
                                  means = means,
                                  sds = proposal_sds))
  }
  target_fc <- function(x) {
    return(dnorm_mix_tempered(x = x,
                              n_comp = length(weights),
                              weights = weights,
                              means = means,
                              sds = sds,
                              beta = beta))
  }
  i <- 0
  samples <- rep(NA, N)
  iters <- 0
  while (i < N) {
    iters <- iters+1
    prop <- rmix_norm(N = 1,
                      weights = weights,
                      means = means,
                      sds = proposal_sds)
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
#' Sample for base level tempered mixture Gaussian
#'
#' @param nsamples number of samples per node
#' @param weights vector: weights of mixture Gaussian
#' @param means vector: means of mixture Gassuan
#' @param sds vector: st.devs of mixture Gaussian
#' @param beta real value between 0 and 1
#' @param proposal_sds vector: st.devs of proposal mixture Gaussian
#' @param dominating_M constant M to bound the target density
#'
#' @return list of length (1/beta) with samples from tempered mixture Gaussian
#' 
#' @export
base_rejection_sampler_mixG <- function(nsamples,
                                        weights,
                                        means,
                                        sds,
                                        beta,
                                        proposal_sds,
                                        dominating_M) {
  if (sum(weights)!=1) {
    stop("base_rejection_sampler_mixG: weights do not add up to 1")
  } else if (!check_equal_dim(list(weights, means, sds))) {
    stop("base_rejection_sampler_mixG: dimensions of input vectors do not match")
  }
  print(paste('sampling from tempered target density with beta =', 1, '/', (1/beta)))
  return(lapply(1:(1/beta), function(i) {
    rmix_norm_tempered_rejection_sampler(N = nsamples,
                                         weights = weights,
                                         means = means,
                                         sds = sds,
                                         beta = beta,
                                         proposal_sds = proposal_sds,
                                         dominating_M = dominating_M)}))
}