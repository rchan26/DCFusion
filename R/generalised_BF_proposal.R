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

#' @export
construct_V <- function(s,
                        t,
                        end_time,
                        C,
                        d,
                        inv_precondition_matrices,
                        Lambda) {
  if (t <= s) {
    stop("construct_V: must have s < t <= end_time")
  } else if (end_time < t) {
    stop("construct_V: must have s < t <= end_time")
  } else if (!is.list(inv_precondition_matrices)) {
    stop("construct_V: inv_precondition_matrices must be a list of length C")
  } else if (length(inv_precondition_matrices)!=C) {
    stop("construct_V: inv_precondition_matrices must be a list of length C")
  }
  if (d==1) {
    if (!is.double(Lambda)) {
      stop("construct_V: if d==1, Lambda must be a double")
    } else if (!all(sapply(1:C, function(c) is.double(inv_precondition_matrices)))) {
      stop("construct_V: if d==1, inv_precondition_matrices[[c]] must a double")
    }
  } else if (d > 1)  {
    if (!is.matrix(Lambda)) {
      stop("construct_V: if d>1, Lambda must be a (d x d) matrix")
    } else if (!all(dim(Lambda)==d)) {
      stop("construct_V: if d>1, Lambda must be a (d x d) matrix")
    } else if (!all(sapply(1:C, function(c) is.matrix(inv_precondition_matrices[[c]])))) {
      stop("construct_V: if d>1, inv_precondition_matrices[[c]] must be a (d x d) matrix for all c=1,...,C")
    } else if (!all(sapply(1:C, function(c) all(dim(inv_precondition_matrices)==d)))) {
      stop("construct_V: if d>1, inv_precondition_matrices[[c]] must be a (d x d) matrix for all c=1,...,C")
    }
  } else {
    stop("construct_V: d must be greater than or equal to 1")
  }
  if (t != end_time) {
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
    return(solve(V_inv))
  } else {
    proposal_cov <- calculate_proposal_cov(time = t-s,
                                           weights = inv_precondition_matrices)
    # create (C x C) block matrix with elements proposal_cov
    block_cols <- do.call(cbind, replicate(C, proposal_cov, simplify = FALSE))
    return(do.call(rbind, replicate(C, block_cols, simplify = FALSE)))
  }
}

#' @export
construct_L_inv <- function(C, d, inv_precondition_matrices) {
  L_inv <- matrix(data = 0, nrow = C*d, ncol = C*d)
  for (i in 1:C) {
    i_fill <- d*(i-1)+(1:d)
    L_inv[i_fill,i_fill] <- inv_precondition_matrices[[i]]
  }
  return(L_inv)
}

#' @export
construct_M <- function(s,
                        t,
                        end_time,
                        C,
                        d,
                        inv_precondition_matrices,
                        Lambda,
                        V,
                        L_inv,
                        sub_posterior_samples,
                        sub_posterior_mean) {
  if (t <= s) {
    stop("construct_M: must have s < t <= end_time")
  } else if (end_time < t) {
    stop("construct_M: must have s < t <= end_time")
  } else if (!is.list(inv_precondition_matrices)) {
    stop("construct_M: inv_precondition_matrices must be a list of length C")
  } else if (length(inv_precondition_matrices)!=C) {
    stop("construct_M: inv_precondition_matrices must be a list of length C")
  }
  if (d==1) {
    if (!is.vector(sub_posterior_samples)) {
      stop("construct_M: if d==1, sub_posterior_samples must be a vector of length C")
    } else if (length(sub_posterior_samples)!=C) {
      stop("construct_M: if d==1, sub_posterior_samples must be a vector of length C")
    } else if (!is.double(Lambda)) {
      stop("construct_M: if d==1, Lambda must be a double")
    } else if (!all(sapply(1:C, function(c) is.double(inv_precondition_matrices)))) {
      stop("construct_M: if d==1, inv_precondition_matrices[[c]] must a double")
    }
  } else if (d > 1)  {
    if (!is.matrix(sub_posterior_samples)){
      stop("construct_M: if d>1, sub_posterior_samples must be a (C x d) matrix")
    } else if (any(dim(sub_posterior_samples)!=c(C,d))) {
      stop("construct_M: if d>1, sub_posterior_samples must be a (C x d) matrix")
    } else if (!is.matrix(Lambda)) {
      stop("construct_M: if d>1, Lambda must be a (d x d) matrix")
    } else if (!all(dim(Lambda)==d)) {
      stop("construct_M: if d>1, Lambda must be a (d x d) matrix")
    } else if (!all(sapply(1:C, function(c) is.matrix(inv_precondition_matrices[[c]])))) {
      stop("construct_M: if d>1, inv_precondition_matrices[[c]] must be a (d x d) matrix for all c=1,...,C")
    } else if (!all(sapply(1:C, function(c) all(dim(inv_precondition_matrices)==d)))) {
      stop("construct_M: if d>1, inv_precondition_matrices[[c]] must be a (d x d) matrix for all c=1,...,C")
    }
  } else {
    stop("construct_M: d must be greater than or equal to 1")
  }
  if (t != end_time) {
    if (!is.matrix(V)) {
      stop("construct_M: if t != end_time, V must be a (C*d x C*d) marix")
    } else if (!all(dim(V)==C*d)) {
      stop("construct_M: if t != end_time, V must be a (C*d x C*d) marix")
    } else if (!is.matrix(L_inv)) {
      stop("construct_M: if t != end_time, L-inv must be a (C*d x C*d) marix")
    } else if (!all(dim(L_inv)==C*d)) {
      stop("construct_M: if t != end_time, L_inv must be a (C*d x C*d) marix")
    }
    M <- as.vector(V%*%L_inv%*%as.vector(t(sub_posterior_samples))/(t-s))
    return(list('M' = M,
                'M_list' = long_vec_to_list_of_vec(C, d, M)))    
  } else {
    return(list('M' = rep(as.vector(sub_posterior_mean), C),
                'M_list' = rep(list(as.vector(sub_posterior_mean)), C)))
  }
}
