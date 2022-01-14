#' @export
construct_V_vanilla <- function(s,
                                t,
                                end_time,
                                C,
                                d) {
  if (t <= s) {
    stop("construct_V_vanilla: must have s < t <= end_time")
  } else if (end_time < t) {
    stop("construct_V_vanilla: must have s < t <= end_time")
  }
  C1 <- ((t-s)*(end_time-t)) / (end_time-s)
  C2 <- ((t-s)^2)/(C*(end_time-s))
  Sigma <- matrix(data = NA,  nrow = C, ncol = C)
  for (i in 1:C) {
    for (j in 1:C) {
      Sigma[i,j] <- C2
      if (i==j) {
        Sigma[i,j] <- Sigma[i,j] + C1
      }
    }
  }
  return(kronecker(Sigma, diag(1, d)))
}

#' @export
construct_M_vanilla <- function(s,
                                t,
                                end_time,
                                C,
                                d,
                                sub_posterior_samples,
                                sub_posterior_mean) {
  if (t <= s) {
    stop("construct_M_vanilla: must have s < t <= end_time")
  } else if (end_time < t) {
    stop("construct_M_vanilla: must have s < t <= end_time")
  }
  if (d==1) {
    if (!is.vector(sub_posterior_samples)) {
      stop("construct_M_vanilla: if d==1, sub_posterior_samples must be a vector of length C")
    } else if (length(sub_posterior_samples)!=C) {
      stop("construct_M_vanilla: if d==1, sub_posterior_samples must be a vector of length C")
    } else if (!is.double(sub_posterior_mean) | length(sub_posterior_mean)!=1) {
      stop("construct_M_vanilla: if d==1, sub_posterior_mean must be a double")
    }
  } else if (d > 1)  {
    if (!is.matrix(sub_posterior_samples)){
      stop("construct_M_vanilla: if d>1, sub_posterior_samples must be a (C x d) matrix")
    } else if (any(dim(sub_posterior_samples)!=c(C,d))) {
      stop("construct_M_vanilla: if d>1, sub_posterior_samples must be a (C x d) matrix")
    } else if (!is.vector(sub_posterior_mean)) {
      stop("construct_M_vanilla: if d>1, sub_posterior_mean must be a vector of length d")
    } else if (length(sub_posterior_mean)!=d) {
      stop("construct_M_vanilla: if d>1, sub_posterior_mean must be a vector of length d")
    }
  } else {
    stop("construct_M_vanilla: d must be greater than or equal to 1")
  }
  C1 <- (end_time-t)/(end_time-s)
  C2 <- (t-s)/(end_time-s)
  M <- rep(NA, C*d)
  M_list <- rep(list(rep(NA, d)), C)
  for (i in 1:C) {
    i_fill <- d*(i-1)+(1:d)
    M[i_fill] <- C1*sub_posterior_samples[i,] + C2*sub_posterior_mean
    M_list[[i]] <- M[i_fill]
  }
  return(list('M' = M, 'M_list' = M_list))
}

#' @export
simulate_xj_vanilla <- function(s,
                                t,
                                end_time,
                                C,
                                d,
                                M,
                                xi = NULL,
                                eta = NULL) {
  if (t <= s) {
    stop("simulate_xj_vanilla: must have s < t <= end_time")
  } else if (end_time < t) {
    stop("simulate_xj_vanilla: must have s < t <= end_time")
  } else if (!is.double(M) | length(M)!=C*d) {
    stop("simulate_xj_vanilla: M must be a vector of length C*d")
  }
  if (is.null(xi)) {
    xi <- mvrnormArma(N = 1, mu = rep(0, d), Sigma = diag(1, d))
  } else {
    if (!is.vector(xi)) {
      stop("simulate_xj_vanilla: if xi is not NULL, xi must be a vector of length d")
    } else if (length(xi)!=d) {
      stop("simulate_xj_vanilla: if xi is not NULL, xi must be a vector of length d")
    }
  }
  if (is.null(eta)) {
    eta <- lapply(1:C, function(c) {
      mvrnormArma(N = 1, mu = rep(0, d), Sigma = diag(1, d))})
  } else {
    if (!is.list(eta)) {
      stop("simulate_xj_vanilla: if eta is not NULL, eta must be a list of length C")
    } else if (length(eta)!=C) {
      stop("simulate_xj_vanilla: if eta is not NULL, eta must be a list of length C")
    } else if (!all(sapply(1:C, function(c) is.vector(eta[[c]])))) {
      stop("simulate_xj_vanilla: if eta is not NULL, eta[[c]] must be a vector of length d for all c=1,...,C")
    } else if (!all(sapply(1:C, function(c) length(eta[[c]])==d))) {
      stop("simulate_xj_vanilla: if eta is not NULL, eta[[c]] must be a vector of length d for all c=1,...,C")
    }
  }
  Delta_j <- t-s
  C1 <- Delta_j / sqrt(C*(end_time-s))
  C2 <- sqrt(Delta_j*(end_time-t)/(end_time-s))
  xj <- rep(NA, C*d)
  xj_list <- rep(list(rep(NA, d)), C)
  for (i in 1:C) {
    i_fill <- d*(i-1)+(1:d)
    xj[i_fill] <- C1*xi + C2*eta[[i]] + M[i_fill]
    xj_list[[i]] <- xj[i_fill]
  }
  return(list('xj' = xj, 'xj_list' = xj_list))
}
