#' @export
construct_V_vanilla <- function(s,
                                t,
                                end_time,
                                C,
                                d,
                                inv = FALSE) {
  if (t <= s) {
    stop("construct_V_vanilla: must have s < t < end_time")
  } else if (end_time <= t) {
    stop("construct_V_vanilla: must have s < t < end_time")
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
  V <- kronecker(Sigma, diag(1, d))
  if (inv) {
    return(list('V' = V, 'V_inv' = solve(V)))  
  } else {
    return(list('V' = V))
  }
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
    stop("construct_M_vanilla: must have s < t < end_time")
  } else if (end_time <= t) {
    stop("construct_M_vanilla: must have s < t < end_time")
  }
  if (d==1) {
    if (!is.vector(sub_posterior_samples)) {
      stop("construct_M_vanilla: if d==1, sub_posterior_samples must be a vector of length C")
    } else if (length(sub_posterior_samples)!=C) {
      stop("construct_M_vanilla: if d==1, sub_posterior_samples must be a vector of length C")
    } else if (!is.double(sub_posterior_mean)) {
      stop("construct_M_vanilla: if d==1, sub_posterior_mean must be a double")
    }
  } else if (d > 1)  {
    if (!is.matrix(sub_posterior_samples)){
      stop("construct_M_vanilla: if d>1, sub_posterior_samples must be a (C x d) matrix")
    } else if (dim(sub_posterior_samples)!=c(C, d)) {
      stop("construct_M_vanilla: if d>1, sub_posterior_samples must be a (C x d) matrix")
    } else if (!is.vector(sub_posterior_mean)) {
      stop("construct_M_vanilla: if d==1, sub_posterior_mean must be a vector of length d")
    } else if (length(sub_posterior_mean)==d) {
      stop("construct_M_vanilla: if d==1, sub_posterior_mean must be a vector of length d")
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
  return(list('M' = M,
              'M_list' = M_list))
}

#' #' @export
#' vanilla_BF_proposal <- function(s,
#'                                 t,
#'                                 end_time,
#'                                 C,
#'                                 d,
#'                                 sub_posterior_samples,
#'                                 sub_posterior_mean) {
#'   V <- construct_V_vanilla(s, t, end_time, C, d)
#'   M <- construct_M_vanilla(s, t, end_time, C, d, sub_posterior_samples, sub_posterior_mean)
#'   return(list('M' = M$M,
#'               'M_list' = M$M_list,
#'               'V' = V$V,
#'               'V_inv' = V$V_inv))
#' }
