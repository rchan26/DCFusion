#' @export
construct_V <- function(s,
                        t,
                        end_time,
                        C,
                        d,
                        precondition_matrices,
                        Lambda) {
  if (t <= s) {
    stop("construct_V: must have s < t <= end_time")
  } else if (end_time < t) {
    stop("construct_V: must have s < t <= end_time")
  } else if (!is.list(precondition_matrices)) {
    stop("construct_V: precondition_matrices must be a list of length C")
  } else if (length(precondition_matrices)!=C) {
    stop("construct_V: precondition_matrices must be a list of length C")
  }
  if (d==1) {
    if (!is.double(Lambda) | length(Lambda)!=1) {
      stop("construct_V: if d==1, Lambda must be a double")
    } else if (!all(sapply(1:C, function(c) is.double(precondition_matrices[[c]]) & length(precondition_matrices[[c]])==1))) {
      stop("construct_V: if d==1, precondition_matrices[[c]] must a double")
    }
  } else if (d > 1)  {
    if (!is.matrix(Lambda)) {
      stop("construct_V: if d>1, Lambda must be a (d x d) matrix")
    } else if (!all(dim(Lambda)==d)) {
      stop("construct_V: if d>1, Lambda must be a (d x d) matrix")
    } else if (!all(sapply(1:C, function(c) is.matrix(precondition_matrices[[c]])))) {
      stop("construct_V: if d>1, precondition_matrices[[c]] must be a (d x d) matrix for all c=1,...,C")
    } else if (!all(sapply(1:C, function(c) all(dim(precondition_matrices)==d)))) {
      stop("construct_V: if d>1, precondition_matrices[[c]] must be a (d x d) matrix for all c=1,...,C")
    }
  } else {
    stop("construct_V: d must be greater than or equal to 1")
  }
  C1 <- ((t-s)*(end_time-t)) / (end_time-s)
  C2 <- ((t-s)^2) / (end_time-s)
  V <- matrix(data = NA, nrow = C*d, ncol = C*d)
  for (i in 1:C) {
    for (j in 1:C) {
      i_fill <- d*(i-1)+(1:d)
      j_fill <- d*(j-1)+(1:d)
      V[i_fill,j_fill] <- C2*Lambda
      if (i==j) {
        V[i_fill,j_fill] <- V[i_fill,j_fill] + C1*precondition_matrices[[i]]
      }
    }
  }
  return(V)
}

#' @export
construct_M <- function(s,
                        t,
                        end_time,
                        C,
                        d,
                        precondition_matrices,
                        sub_posterior_samples,
                        sub_posterior_mean) {
  if (t <= s) {
    stop("construct_M: must have s < t <= end_time")
  } else if (end_time < t) {
    stop("construct_M: must have s < t <= end_time")
  } else if (!is.list(precondition_matrices)) {
    stop("construct_M: precondition_matrices must be a list of length C")
  } else if (length(precondition_matrices)!=C) {
    stop("construct_M: precondition_matrices must be a list of length C")
  }
  if (d==1) {
    if (!is.vector(sub_posterior_samples)) {
      stop("construct_M: if d==1, sub_posterior_samples must be a vector of length C")
    } else if (length(sub_posterior_samples)!=C) {
      stop("construct_M: if d==1, sub_posterior_samples must be a vector of length C")
    } else if (!all(sapply(1:C, function(c) is.double(precondition_matrices[[c]]) & length(precondition_matrices[[c]]==1)))) {
      stop("construct_M: if d==1, precondition_matrices[[c]] must a double")
    } else if (length(sub_posterior_mean)!=1) {
      stop("construct_M: if d==1, sub_posterior_mean must be a numeric of length 1")
    }
  } else if (d > 1)  {
    if (!is.matrix(sub_posterior_samples)){
      stop("construct_M: if d>1, sub_posterior_samples must be a (C x d) matrix")
    } else if (any(dim(sub_posterior_samples)!=c(C,d))) {
      stop("construct_M: if d>1, sub_posterior_samples must be a (C x d) matrix")
    } else if (!all(sapply(1:C, function(c) is.matrix(precondition_matrices[[c]])))) {
      stop("construct_M: if d>1, precondition_matrices[[c]] must be a (d x d) matrix for all c=1,...,C")
    } else if (!all(sapply(1:C, function(c) all(dim(precondition_matrices)==d)))) {
      stop("construct_M: if d>1, precondition_matrices[[c]] must be a (d x d) matrix for all c=1,...,C")
    } else if (length(sub_posterior_mean)!=d) {
      stop("construct_M: if d>1, sub_posterior_mean must be a vector of length d")
    }
  } else {
    stop("construct_M: d must be greater than or equal to 1")
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
simulate_xj <- function(s,
                        t,
                        end_time,
                        C,
                        d,
                        precondition_matrices,
                        Lambda,
                        M,
                        xi = NULL,
                        eta = NULL) {
  if (t <= s) {
    stop("simulate_xj: must have s < t <= end_time")
  } else if (end_time < t) {
    stop("simulate_xj: must have s < t <= end_time")
  } else if (!is.list(precondition_matrices)) {
    stop("simulate_xj: precondition_matrices must be a list of length C")
  } else if (length(precondition_matrices)!=C) {
    stop("simulate_xj: precondition_matrices must be a list of length C")
  } else if (!is.double(M) | length(M)!=C*d) {
    stop("simulate_xj: M must be a vector of length C*d")
  }
  if (d==1) {
    if (!is.double(Lambda) | length(Lambda)!=1) {
      stop("simulate_xj: if d==1, Lambda must be a double")
    } else if (!all(sapply(1:C, function(c) is.double(precondition_matrices[[c]]) & length(precondition_matrices[[c]])==1))) {
      stop("simulate_xj: if d==1, precondition_matrices[[c]] must a double")
    }
  } else if (d > 1)  {
    if (!is.matrix(Lambda)) {
      stop("simulate_xj: if d>1, Lambda must be a (d x d) matrix")
    } else if (!all(dim(Lambda)==d)) {
      stop("simulate_xj: if d>1, Lambda must be a (d x d) matrix")
    } else if (!all(sapply(1:C, function(c) is.matrix(precondition_matrices[[c]])))) {
      stop("simulate_xj: if d>1, precondition_matrices[[c]] must be a (d x d) matrix for all c=1,...,C")
    } else if (!all(sapply(1:C, function(c) all(dim(precondition_matrices)==d)))) {
      stop("simulate_xj: if d>1, precondition_matrices[[c]] must be a (d x d) matrix for all c=1,...,C")
    }
  }
  if (is.null(xi)) {
    xi <- mvrnormArma(N = 1, mu = rep(0, d), Sigma = Lambda)  
  } else {
    if (!is.vector(xi)) {
      stop("simulate_xj: if xi is not NULL, xi must be a vector of length d")
    } else if (length(xi)!=d) {
      stop("simulate_xj: if xi is not NULL, xi must be a vector of length d")
    }
  }
  if (is.null(eta)) {
    eta <- lapply(1:C, function(c) {
      mvrnormArma(N = 1, mu = rep(0, d), Sigma = precondition_matrices[[c]])})
  } else {
    if (!is.list(eta)) {
      stop("simulate_xj: if eta is not NULL, eta must be a list of length C")
    } else if (length(eta)!=C) {
      stop("simulate_xj: if eta is not NULL, eta must be a list of length C")
    } else if (!all(sapply(1:C, function(c) is.vector(eta[[c]])))) {
      stop("simulate_xj: if eta is not NULL, eta[[c]] must be a vector of length d for all c=1,...,C")
    } else if (!all(sapply(1:C, function(c) length(eta[[c]])==d))) {
      stop("simulate_xj: if eta is not NULL, eta[[c]] must be a vector of length d for all c=1,...,C")
    }
  }
  Delta_j <- t-s
  C1 <- Delta_j / sqrt(end_time-s)
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
