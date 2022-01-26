#' @export
T_guidance <- function(condition,
                       threshold = NULL,
                       C,
                       d,
                       data_size = NULL,
                       b = NULL,
                       sub_posterior_means,
                       precondition_matrices = NULL,
                       inv_precondition_matrices = NULL,
                       Lambda = NULL,
                       lambda = NULL,
                       gamma = NULL,
                       k1 = NULL,
                       k2 = NULL,
                       vanilla = NULL) {
  if (!(condition %in% c('SH', 'SSH'))) {
    stop("T_guidance: condition must be either \"SH\" or \"SSH\"")
  }
  if (is.null(threshold)) {
    if (is.null(k1)) {
      stop("T_guidance: if threshold is not passed, k1 must be passed in")
    } else if (is.null(k2) & condition=="SSH") {
      stop("T_guidance: if threshold is not passed and condition==\"SSH\", k1 and k2 must be passed in")
    }
  } else if ((threshold < 0) | (threshold > 1)) {
    stop("T_guidance: threshold must be between 0 and 1")
  }
  if (d==1) {
    if (!is.vector(sub_posterior_means)) {
      stop("T_guidance: if d==1, sub_posterior_means must be a vector of length C")
    } else if (length(sub_posterior_means)!=C) {
      stop("T_guidance: if d==1, sub_posterior_means must be a vector of length C")
    }
  } else if (d > 1)  {
    if (!is.matrix(sub_posterior_samples)){
      stop("T_guidance: if d>1, sub_posterior_samples must be a (C x d) matrix")
    } else if (any(dim(sub_posterior_samples)!=c(C,d))) {
      stop("T_guidance: if d>1, sub_posterior_samples must be a (C x d) matrix")
    } 
  } else {
    stop("T_guidance: d must be greater than or equal to 1")
  }
  if (is.null(b)) {
    b <- 1
  } else {
    if (b < 0) {
      stop("T_guidance: b must be greater than 0")
    }
  }
  if (is.null(data_size)) {
    data_size <- C
  } else {
    if (data_size < 0) {
      stop("T_guidance: data_size must be greater than 0")
    }
  }
  if (is.null(vanilla)) {
    vanilla <- FALSE
  }
  if (vanilla) {
    if (d==1) {
      precondition_matrices <- rep(1, C)
      inv_precondition_matrices <- rep(1, C)
    } else {
      precondition_matrices <- rep(list(diag(1,d)), C)
      inv_precondition_matrices <- rep(list(diag(1,d)), C)
      Lambda <- diag(1,d)/C
    }
  } else {
    # set b to 1 if vanilla == FALSE
    b <- 1
    if (d==1) {
      if (!all(sapply(1:C, function(c) is.double(precondition_matrices[[c]]) & length(precondition_matrices[[c]]==1)))) {
        stop("T_guidance: if d==1, precondition_matrices[[c]] must a double for all c=1,...,C")
      } else if (!all(sapply(1:C, function(c) is.double(inv_precondition_matrices[[c]]) & length(inv_precondition_matrices[[c]]==1)))) {
        stop("T_guidance: if d==1, inv_precondition_matrices[[c]] must a double for all c=1,...,C")
      }
      precondition_matrices <- unlist(precondition_matrices)
      inv_precondition_matrices <- unlist(inv_precondition_matrices)
    } else if (d>1) {
      if (!is.matrix(Lambda)) {
        stop("T_guidance: if d>1, Lambda must be a (d x d) matrix")
      } else if (!all(dim(Lambda)==d)) {
        stop("T_guidance: if d>1, Lambda must be a (d x d) matrix")
      } else if (!all(sapply(1:C, function(c) is.matrix(precondition_matrices[[c]])))) {
        stop("T_guidance: if d>1, precondition_matrices[[c]] must be a (d x d) matrix for all c=1,...,C")
      } else if (!all(sapply(1:C, function(c) all(dim(precondition_matrices)==d)))) {
        stop("T_guidance: if d>1, precondition_matrices[[c]] must be a (d x d) matrix for all c=1,...,C")
      } else if (!all(sapply(1:C, function(c) is.matrix(inv_precondition_matrices[[c]])))) {
        stop("T_guidance: if d>1, inv_precondition_matrices[[c]] must be a (d x d) matrix for all c=1,...,C")
      } else if (!all(sapply(1:C, function(c) all(dim(inv_precondition_matrices)==d)))) {
        stop("T_guidance: if d>1, inv_precondition_matrices[[c]] must be a (d x d) matrix for all c=1,...,C")
      }
    }
  }
  ##### compute mean and variance of sub-posterior means #####
  if (d==1) {
    a_tilde <- weighted_mean_univariate(x = sub_posterior_means,
                                        weights = inv_precondition_matrices)
    sigma_a_2 <- weighted_variance_univariate(x = sub_posterior_means,
                                              x_mean = a_tilde,
                                              precondition_values = precondition_matrices)
  } else {
    a_tilde <- weighted_mean_multivariate(matrix = sub_posterior_means,
                                          weights = inv_precondition_matrices,
                                          inverse_sum_weights = Lambda)
    sigma_a_2 <- weighted_variance_multivariate(x = sub_posterior_means,
                                                x_mean = a_tilde,
                                                inv_precondition_matrices = inv_precondition_matrices)
  }
  # if vanilla == FALSE, b is set to 1 above
  if (condition == 'SH') {
    ##### SH(lambda) #####
    if (is.null(lambda)) {
      lambda <- data_size*sigma_a_2/(b*(C-1))
    } else {
      if (lambda < 0) {
        stop("T_guidance: lambda must be greater than 0")
      }
    }
    if (is.null(k1)) {
      k1 <- sqrt(-(lambda+d/2)/log(threshold))  
    } else {
      if (k1 < 0) {
        stop("T_guidance: k1 must be greater than 0")
      }
    }
    CESS_0_threshold <- exp(-(lambda+d/2)/(k1^2))
    min_T <- b*(C^(1.5))*k1/data_size
    return(list('min_T' = min_T,
                'CESS_0_threshold' = CESS_0_threshold,
                'sigma_a_2' = sigma_a_2,
                'lambda' = lambda,
                'k1' = k1))
  } else if (condition == 'SSH') {
    ##### SSH(gamma) #####
    if (is.null(gamma)) {
      gamma <- b*sigma_a_2  
    } else {
      if (gamma < 0) {
        stop("T_guidance: gamma must be greater than 0")
      }
    }
    if (is.null(k1)) {
      k1 <- sqrt(-((data_size*gamma/C)+(d/2))/log(threshold))  
    } else {
      if (k1 < 0) {
        stop("T_guidance: k1 must be greater than 0")
      }
    }
    if (is.null(k2)) {
      k2 <- b*C*k1/data_size  
    } else {
      if (k2 < 0) {
        stop("T_guidance: k2 must be greater than 0")
      }
    }
    k1_order <- sqrt(data_size*gamma/C)
    k2_order <- b*sqrt(C*gamma/data_size)
    T1 <- b*(C^(1.5))*k1/data_size
    T2 <- (C^(0.5))*k2
    CESS_0_threshold <- exp(-(gamma/(k1*k2))-(d/(2*k1^2)))
    return(list('min_T' = max(T1, T2),
                'CESS_0_threshold' = CESS_0_threshold,
                'T1' = T1,
                'T2' = T2,
                'sigma_a_2' = sigma_a_2,
                'gamma' = gamma,
                'k1' = k1,
                'k2' = k2,
                'k1_order' = k1_order,
                'k2_order' = k2_order))
  }
}

#' @export
mesh_guidance_adaptive <- function(C,
                                   d,
                                   data_size = NULL,
                                   b = NULL,
                                   trajectories,
                                   sub_posterior_means,
                                   inv_precondition_matrices = NULL,
                                   k3 = NULL,
                                   k4 = NULL,
                                   T2 = NULL,
                                   vanilla = NULL) {
  if (d==1) {
    if (!is.vector(sub_posterior_means)) {
      stop("T_guidance: if d==1, sub_posterior_means must be a vector of length C")
    } else if (length(sub_posterior_means)!=C) {
      stop("T_guidance: if d==1, sub_posterior_means must be a vector of length C")
    }
  } else if (d > 1)  {
    if (!is.matrix(sub_posterior_samples)){
      stop("T_guidance: if d>1, sub_posterior_samples must be a (C x d) matrix")
    } else if (any(dim(sub_posterior_samples)!=c(C,d))) {
      stop("T_guidance: if d>1, sub_posterior_samples must be a (C x d) matrix")
    } 
  } else {
    stop("T_guidance: d must be greater than or equal to 1")
  }
  if (is.null(b)) {
    b <- 1
  } else {
    if (b < 0) {
      stop("T_guidance: b must be greater than 0")
    }
  }
  if (is.null(data_size)) {
    data_size <- C
  } else {
    if (data_size < 0) {
      stop("T_guidance: data_size must be greater than 0")
    }
  }
  if (is.null(vanilla)) {
    vanilla <- FALSE
  }
  if (vanilla) {
    if (d==1) {
      inv_precondition_matrices <- rep(1, C)
    } else {
      inv_precondition_matrices <- rep(list(diag(1,d)), C)
    }
  } else {
    # set b to 1 if vanilla == FALSE
    b <- 1
    if (d==1) {
      if (!all(sapply(1:C, function(c) is.double(inv_precondition_matrices[[c]]) & length(inv_precondition_matrices[[c]]==1)))) {
        stop("T_guidance: if d==1, inv_precondition_matrices[[c]] must a double for all c=1,...,C")
      }
      inv_precondition_matrices <- unlist(inv_precondition_matrices)
    } else if (d>1) {
      if (!all(sapply(1:C, function(c) is.matrix(inv_precondition_matrices[[c]])))) {
        stop("T_guidance: if d>1, inv_precondition_matrices[[c]] must be a (d x d) matrix for all c=1,...,C")
      } else if (!all(sapply(1:C, function(c) all(dim(inv_precondition_matrices)==d)))) {
        stop("T_guidance: if d>1, inv_precondition_matrices[[c]] must be a (d x d) matrix for all c=1,...,C")
      }
    }
  }
  ##### compute average variation of the trajectories #####
  if (d==1) {
    E_nu_j <- weighted_trajectory_variation_univariate(x_samples = trajectories,
                                                       sub_posterior_means = sub_posterior_means,
                                                       precondition_values = 1/inv_precondition_matrices)
  } else {
    E_nu_j <- weighted_trajectory_variation_multivariate(x_samples = trajectories,
                                                         sub_posterior_means = sub_posterior_means,
                                                         inv_precondition_matrices = inv_precondition_matrices)
  }
  if (is.null(k3)) {
    k3 <- 0.5
  } else {
    if (k3 < 0) {
      stop("T_guidance: k3 must be greater than 0")
    }
  }
  if (is.null(k4)) {
    k4 <- 0.5
  } else {
    if (k4< 0) {
      stop("T_guidance: k4 must be greater than 0")
    }
  }
  T1 <- (b^2)*C*k3/(E_nu_j*(data_size^2))
  if (is.null(T2)) {
    T2 <- sqrt((b^2)*C*k4/(2*(data_size^2)*d))
  }
  return(list('max_delta_j' = min(T1, T2),
              'CESS_j_treshold' = exp(-k3-k4),
              'T1' = T1,
              'T2' = T2,
              'E_nu_j' = E_nu_j))
}

#' @export
mesh_guidance_regular <- function(C,
                                  d,
                                  data_size = NULL,
                                  b = NULL,
                                  k3 = NULL,
                                  k4 = NULL,
                                  vanilla = NULL) {
  if (is.null(b)) {
    b <- 1
  } else {
    if (b < 0) {
      stop("mesh_guidance_regular: b must be greater than 0")
    }
  }
  if (is.null(data_size)) {
    data_size <- C
  } else {
    if (data_size < 0) {
      stop("mesh_guidance_regular: data_size must be greater than 0")
    }
  }
  if (is.null(vanilla)) {
    vanilla <- FALSE
  }
  if (!vanilla) {
    b <- 1
  }
  if (is.null(k3)) {
    k3 <- 0.5
  } else {
    if (k3 < 0) {
      stop("mesh_guidance_regular: k3 must be greater than 0")
    }
  }
  if (is.null(k4)) {
    k4 <- 0.5
  } else {
    if (k4< 0) {
      stop("mesh_guidance_regular: k4 must be greater than 0")
    }
  }
  return(list('max_delta_j' = sqrt((b^2)*C*k4/(2*(data_size^2)*d)),
              'CESS_j_treshold' = exp(-k3-k4)))
}

#' @export
BF_guidance <- function(condition,
                        CESS_0_threshold = NULL,
                        C,
                        d,
                        data_size = NULL,
                        b = NULL,
                        sub_posterior_means,
                        precondition_matrices = NULL,
                        inv_precondition_matrices = NULL,
                        Lambda = NULL,
                        lambda = NULL,
                        gamma = NULL,
                        k1 = NULL,
                        k2 = NULL,
                        k3 = NULL,
                        k4 = NULL,
                        vanilla = NULL) {
  T_guide <- T_guidance(condition = condition,
                        threshold = CESS_0_threshold,
                        C = C,
                        d = d,
                        data_size = data_size,
                        b = b,
                        sub_posterior_means = sub_posterior_means,
                        precondition_matrices = precondition_matrices,
                        inv_precondition_matrices = inv_precondition_matrices,
                        Lambda = Lambda,
                        lambda = lambda,
                        gamma = gamma,
                        k1 = k1,
                        k2 = k2,
                        vanilla = vanilla)
  mesh_guide <- mesh_guidance_regular(C = C,
                                      d = d,
                                      data_size = data_size,
                                      b = b,
                                      k3 = k3,
                                      k4 = k4,
                                      vanilla = vanilla)
  rec_mesh <- seq(from = 0, 
                  to = T_guide$min_T,
                  length.out = ceiling(T_guide$min_T/mesh_guide$max_delta_j)+1)
  return(c(T_guide, mesh_guide, list('n' = length(rec_mesh), 'mesh' = rec_mesh)))
}
