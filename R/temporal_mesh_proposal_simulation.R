#' @export
construct_V_vanilla <- function(s, t, end_time, C, d) {
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
  return(list('V' = V, 'V_inv' = solve(V)))
}

#' @export
construct_M_vanilla <- function(s, t, end_time, C, d, sub_posterior_samples) {
  sum_vec <- sub_posterior_samples[[C]]
  for (i in 1:(C-1)) {
    sum_vec <- sum_vec + sub_posterior_samples[[i]]
  }
  x_bar <- sum_vec / C
  C1 <- (end_time-t)/(end_time-s)
  C2 <- (t-s)/(end_time-s)
  M <- rep(NA, C*d)
  M_list <- rep(list(rep(NA, d)), C)
  for (i in 1:C) {
    i_fill <- d*(i-1)+(1:d)
    M[i_fill] <- C1*sub_posterior_samples[[i]] + C2*x_bar
    M_list[[i]] <- M[i_fill]
  }
  return(list('M' = M, 'M_list' = M_list))
}

#' @export
long_vec_to_list_of_vec <- function(C, d, vec) {
  list_of_vec <- rep(list(rep(NA, d)), C)
  for (c in 1:C) {
    indices <- d*(c-1)+(1:d)
    list_of_vec[[c]] <- vec[indices]
  }
  return(list_of_vec)
}

#' @export
construct_V <- function(s, t, end_time, C, d, Lambda, inv_precondition_matrices) {
  C1 <- (end_time-s) / ((t-s)*(end_time-t))
  C2 <- -1/(end_time-t)
  V_inv <- matrix(data = NA, nrow = C*d, ncol = C*d)
  for (i in 1:C) {
    for (j in 1:C) {
      i_fill <- d*(i-1)+(1:d)
      j_fill <- d*(j-1)+(1:d)
      V_inv[i_fill,j_fill] <- C2*(inv_precondition_matrices[[i]]%*%Lambda%*%inv_precondition_matrices[[j]])
      if (i==j) {
        V_inv[i_fill,j_fill] <- V_inv[i_fill,j_fill] + C1*inv_precondition_matrices[[i]]
      }
    }
  }
  return(list('V' = solve(V_inv), 'V_inv' = V_inv))
}

#' @export
construct_L_inv <- function(inv_precondition_matrices) {
  L_inv <- matrix(data = 0, nrow = C*d, ncol = C*d)
  for (i in 1:C) {
    i_fill <- d*(i-1)+(1:d)
    L_inv[i_fill,i_fill] <- inv_precondition_matrices[[i]]
  }
  return(L_inv)
}

#' @export
construct_M <- function(s, t, end_time, C, d, Lambda, inv_precondition_matrices, sub_posterior_samples) {
  covariance_mat <- construct_V(s, t, end_time, C, d, Lambda, inv_precondition_matrices)
  V <- covariance_mat$V
  L_inv <- construct_L_inv(inv_precondition_matrices)
  x <- unlist(sub_posterior_samples)
  M <- as.vector(V%*%L_inv%*%x / (t-s))
  return(list('M' = M, 'M_list' = long_vec_to_list_of_vec(C, d, M)))
}
