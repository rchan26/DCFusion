#' @export
long_vec_to_list_of_vec <- function(C,
                                    d,
                                    vec) {
  list_of_vec <- rep(list(rep(NA, d)), C)
  for (c in 1:C) {
    indices <- d*(c-1)+(1:d)
    list_of_vec[[c]] <- vec[indices]
  }
  return(list_of_vec)
}

construct_V <- function(s,
                        t,
                        end_time,
                        C,
                        d,
                        Lambda,
                        inv_precondition_matrices) {
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
  return(list('V' = solve(V_inv),
              'V_inv' = V_inv))
}

construct_L_inv <- function(C, d, inv_precondition_matrices) {
  L_inv <- matrix(data = 0, nrow = C*d, ncol = C*d)
  for (i in 1:C) {
    i_fill <- d*(i-1)+(1:d)
    L_inv[i_fill,i_fill] <- inv_precondition_matrices[[i]]
  }
  return(L_inv)
}

construct_M <- function(s,
                        t,
                        end_time,
                        C,
                        d,
                        Lambda,
                        inv_precondition_matrices,
                        sub_posterior_samples) {
  covariance_mat <- construct_V(s, t, end_time, C, d, Lambda, inv_precondition_matrices)
  L_inv <- construct_L_inv(inv_precondition_matrices)
  M <- as.vector(covariance_mat$V%*%L_inv%*%unlist(sub_posterior_samples)/(t-s))
  return(list('M' = M,
              'M_list' = long_vec_to_list_of_vec(C, d, M),
              'V' = covariance_mat$V,
              'V_inv' = covariance_mat$V_inv))
}

#' @export
generalised_BF_proposal <- function(s,
                                    t,
                                    end_time,
                                    C,
                                    d,
                                    Lambda,
                                    inv_precondition_matrices,
                                    sub_posterior_samples) {
  if (t <= s) {
    stop("generalised_BF_proposal: must have s < t < end_time")
  } else if (end_time <= t) {
    stop("generalised_BF_proposal: must have s < t < end_time")
  } else if (!is.list(sub_posterior_samples)) {
    stop("generalised_BF_proposal: sub_posterior_samples must be a list of length C")
  } else if (length(sub_posterior_samples)!=C) {
    stop("generalised_BF_proposal: sub_posterior_samples must be a list of length C")
  } else if (!is.list(inv_precondition_matrices)) {
    stop("generalised_BF_proposal: inv_precondition_matrices must be a list of length C")
  } else if (length(inv_precondition_matrices)!=C) {
    stop("generalised_BF_proposal: inv_precondition_matrices must be a list of length C")
  } 
  if (d==1) {
    if (!all(sapply(1:length(sub_posterior_samples), function(i) is.double(sub_posterior_samples[[i]])))) {
      stop("generalised_BF_proposal: if d==1, sub_posterior_samples[[i]] must be a double for all i=1,...,C")
    } else if (!is.double(Lambda)) {
      stop("generalised_BF_proposal: if d==1, Lambda must be a double")
    } else if (!all(sapply(1:length(sub_posterior_samples), function(i) is.double(inv_precondition_matrices[[i]])))) {
      stop("generalised_BF_proposal: if d==1, inv_precondition_matrices[[i]] must be a double for all i=1,...,C")
    }
  } else if (d > 1) {
    if (!all(sapply(1:length(sub_posterior_samples), function(i) length(sub_posterior_samples[[i]])==d))) {
      stop("vanilla_BF_proposal: if d>1, sub_posterior_samples[[i]] must be a vector of length d for all i=1,...,C")
    } else if (!all(dim(Lambda)==d)) {
      stop("generalised_BF_proposal: if d>1, Lambda must matrix with dimension dxd")
    } else if (!all(sapply(1:length(sub_posterior_samples), function(i) all(dim(inv_precondition_matrices[[i]])==d)))) {
      stop("generalised_BF_proposal: if d==1, inv_precondition_matrices[[i]] must matrix with dimension dxd for all i=1,...,C")
    }
  } else {
    stop("generalised_BF_proposal: d must be greater than or equal to 1")
  }
  M <- construct_M(s, t, end_time, C, d, Lambda, inv_precondition_matrices, sub_posterior_samples)
  return(list('M' = M$M,
              'M_list' = M$M_list,
              'V' = M$V,
              'V_inv' = M$V_inv))
}
