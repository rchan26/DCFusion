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
  } else if (d > 1) {
    if (!is.matrix(sub_posterior_means)){
      stop("T_guidance: if d>1, sub_posterior_means must be a (C x d) matrix")
    } else if (any(dim(sub_posterior_means)!=c(C,d))) {
      stop("T_guidance: if d>1, sub_posterior_means must be a (C x d) matrix")
    } 
  } else {
    stop("T_guidance: d must be greater than or equal to 1")
  }
  if (is.null(vanilla)) {
    warning('T_guidance: vanilla is set to FALSE by default')
    vanilla <- FALSE
  }
  if (is.null(data_size)) {
    warning('T_guidance: data_size is set to 1 by default')
    data_size <- 1
  } else {
    if (data_size < 0) {
      stop("T_guidance: data_size must be greater than 0")
    }
  }
  if (is.null(b)) {
    if (vanilla) {
      warning('T_guidance: b is set to 1 by default')
      b <- 1
    } else {
      warning('T_guidance: b is set to data_size/C by default')
      b <- data_size/C
    }
  } else {
    if (b < 0) {
      stop("T_guidance: b must be greater than 0")
    }
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
  if (condition == 'SH') {
    if (is.null(lambda)) {
      warning("T_guidance: lambda is set to 1 by default")
      lambda <- 1
    } else {
      if (lambda < 0) {
        stop("T_guidance: lambda must be greater than 0")
      }
    }
    if (is.null(k1)) {
      warning("T_guidance: k1 set to sqrt(-(lambda+d/2)/log(threshold)) = ",
              sqrt(-(lambda+d/2)/log(threshold)), " by default")
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
                'b' = b,
                'k1' = k1,
                'precondition_matrices' = precondition_matrices,
                'inv_precondition_matrices' = inv_precondition_matrices,
                'Lambda' = Lambda))
  } else if (condition == 'SSH') {
    if (is.null(gamma)) {
      warning("T_guidance: gamma is set to sigma_a_2/b = ", sigma_a_2/b, " by default")
      gamma <- sigma_a_2/b
    } else {
      if (gamma < 0) {
        stop("T_guidance: gamma must be greater than 0")
      }
    }
    if (is.null(k1)) {
      warning("T_guidance: k1 set to sqrt(-((data_size*gamma/C)+(d/2))/log(threshold)) = ",
              sqrt(-((data_size*gamma/C)+(d/2))/log(threshold)), " by default")
      k1 <- sqrt(-((data_size*gamma/C)+(d/2))/log(threshold))
    } else {
      if (k1 < 0) {
        stop("T_guidance: k1 must be greater than 0")
      }
    }
    if (is.null(k2)) {
      warning("T_guidance: k2 set to b*C*k1/data_size = ",
              b*C*k1/data_size, " by default")
      k2 <- b*C*k1/data_size
    } else {
      if (k2 < 0) {
        stop("T_guidance: k2 must be greater than 0")
      }
    }
    T1 <- b*(C^(1.5))*k1/data_size
    T2 <- (C^(0.5))*k2
    CESS_0_threshold <- exp(-(gamma/(k1*k2))-(d/(2*k1^2)))
    return(list('min_T' = max(T1, T2),
                'CESS_0_threshold' = CESS_0_threshold,
                'T1' = T1,
                'T2' = T2,
                'sigma_a_2' = sigma_a_2,
                'gamma' = gamma,
                'b' = b,
                'k1' = k1,
                'k2' = k2,
                'precondition_matrices' = precondition_matrices,
                'inv_precondition_matrices' = inv_precondition_matrices,
                'Lambda' = Lambda))
  }
}

mesh_T1 <- function(k3, b, C, E_nu_j, data_size, logarithm = FALSE) {
  if (logarithm) {
    return(2*log(b)+log(C)+log(k3)-log(E_nu_j)-2*log(data_size))
  } else {
    return(b*b*C*k3/(E_nu_j*data_size*data_size))  
  }
}

mesh_T2 <- function(k4, b, C, data_size, d, logarithm = FALSE) {
  if (logarithm) {
    return(0.5*(2*log(b)+log(C)+log(k4)-log(2)-2*log(data_size)-log(d)))
  } else {
    return(sqrt(b*b*C*k4/(2*data_size*data_size*d)))  
  }
}

#' @export
mesh_guidance_adaptive <- function(C,
                                   d,
                                   data_size = NULL,
                                   b = NULL,
                                   particle_set,
                                   sub_posterior_means,
                                   inv_precondition_matrices = NULL,
                                   k3 = NULL,
                                   k4 = NULL,
                                   T2 = NULL,
                                   vanilla = NULL) {
  if (!("particle" %in% class(particle_set))) {
    stop("mesh_guidance_adaptive: particle_set must be a \"particle\" object")
  }
  if (d==1) {
    if (!is.vector(sub_posterior_means)) {
      stop("mesh_guidance_adaptive: if d==1, sub_posterior_means must be a vector of length C")
    } else if (length(sub_posterior_means)!=C) {
      stop("mesh_guidance_adaptive: if d==1, sub_posterior_means must be a vector of length C")
    }
  } else if (d > 1) {
    if (!is.matrix(sub_posterior_means)){
      stop("mesh_guidance_adaptive: if d>1, sub_posterior_means must be a (C x d) matrix")
    } else if (any(dim(sub_posterior_means)!=c(C,d))) {
      stop("mesh_guidance_adaptive: if d>1, sub_posterior_means must be a (C x d) matrix")
    }
  } else {
    stop("mesh_guidance_adaptive: d must be greater than or equal to 1")
  }
  if (is.null(vanilla)) {
    warning('mesh_guidance_adaptive: vanilla is set to FALSE by default')
    vanilla <- FALSE
  }
  if (is.null(data_size)) {
    warning('mesh_guidance_adaptive: data_size is set to 1 by default')
    data_size <- 1
  } else {
    if (data_size < 0) {
      stop("mesh_guidance_adaptive: data_size must be greater than 0")
    }
  }
  if (is.null(b)) {
    if (vanilla) {
      warning('mesh_guidance_adaptive: b is set to 1 by default')
      b <- 1
    } else {
      warning('mesh_guidance_adaptive: b is set to data_size/C by default')
      b <- data_size/C
    }
  } else {
    if (b < 0) {
      stop("mesh_guidance_adaptive: b must be greater than 0")
    }
  }
  if (vanilla) {
    if (d==1) {
      inv_precondition_matrices <- rep(1, C)
    } else {
      inv_precondition_matrices <- rep(list(diag(1,d)), C)
    }
  } else {
    if (d==1) {
      if (!all(sapply(1:C, function(c) is.double(inv_precondition_matrices[[c]]) & length(inv_precondition_matrices[[c]]==1)))) {
        stop("mesh_guidance_adaptive: if d==1, inv_precondition_matrices[[c]] must a double for all c=1,...,C")
      }
      inv_precondition_matrices <- unlist(inv_precondition_matrices)
    } else if (d>1) {
      if (!all(sapply(1:C, function(c) is.matrix(inv_precondition_matrices[[c]])))) {
        stop("mesh_guidance_adaptive: if d>1, inv_precondition_matrices[[c]] must be a (d x d) matrix for all c=1,...,C")
      } else if (!all(sapply(1:C, function(c) all(dim(inv_precondition_matrices)==d)))) {
        stop("mesh_guidance_adaptive: if d>1, inv_precondition_matrices[[c]] must be a (d x d) matrix for all c=1,...,C")
      }
    }
  }
  if (d==1) {
    E_nu_j <- weighted_trajectory_variation_univariate(x_samples = particle_set$x_samples,
                                                       normalised_weights = particle_set$normalised_weights,
                                                       sub_posterior_means = sub_posterior_means,
                                                       precondition_values = 1/inv_precondition_matrices)
    E_nu_j_old <- weighted_trajectory_variation_univariate(x_samples = particle_set$x_samples,
                                                           normalised_weights = rep(1/length(particle_set$x_samples), length(particle_set$x_samples)),
                                                           sub_posterior_means = sub_posterior_means,
                                                           precondition_values = 1/inv_precondition_matrices)
  } else {
    E_nu_j <- weighted_trajectory_variation_multivariate(x_samples = particle_set$x_samples,
                                                         normalised_weights = particle_set$normalised_weights,
                                                         sub_posterior_means = sub_posterior_means,
                                                         inv_precondition_matrices = inv_precondition_matrices)
    E_nu_j_old <- weighted_trajectory_variation_multivariate(x_samples = particle_set$x_samples,
                                                             normalised_weights = rep(1/length(particle_set$x_samples), length(particle_set$x_samples)),
                                                             sub_posterior_means = sub_posterior_means,
                                                             inv_precondition_matrices = inv_precondition_matrices)
  }
  if (is.null(k3)) {
    warning('mesh_guidance_adaptive: k3 is set to 0.5 by default')
    k3 <- 0.5
  } else {
    if (k3 < 0) {
      stop("mesh_guidance_adaptive: k3 must be greater than 0")
    }
  }
  if (is.null(k4)) {
    warning('mesh_guidance_adaptive: k4 is set to 0.5 by default')
    k4 <- 0.5
  } else {
    if (k4< 0) {
      stop("mesh_guidance_adaptive: k4 must be greater than 0")
    }
  }
  T1 <- mesh_T1(k3, b, C, E_nu_j, data_size)
  if (is.null(T2)) {
    T2 <- mesh_T2(k4, b, C, data_size, d)
  }
  return(list('max_delta_j' = min(T1, T2),
              'CESS_j_treshold' = exp(-k3-k4),
              'T1' = T1,
              'T2' = T2,
              'E_nu_j' = E_nu_j,
              'E_nu_j_old' = E_nu_j_old))
}

#' @export
mesh_guidance_regular <- function(C,
                                  d,
                                  data_size = NULL,
                                  b = NULL,
                                  threshold = NULL,
                                  k3 = NULL,
                                  k4 = NULL,
                                  max_E_nu_j = NULL,
                                  trial_k3_by = 0.000001,
                                  vanilla = NULL) {
  if (is.null(vanilla)) {
    warning('mesh_guidance_regular: vanilla is set to FALSE by default')
    vanilla <- FALSE
  }
  if (is.null(data_size)) {
    warning('mesh_guidance_regular: data_size is set to 1 by default')
    data_size <- 1
  } else {
    if (data_size < 0) {
      stop("mesh_guidance_regular: data_size must be greater than 0")
    }
  }
  if (is.null(b)) {
    if (vanilla) {
      warning('mesh_guidance_regular: b is set to 1 by default')
      b <- 1
    } else {
      warning('mesh_guidance_regular: b is set to data_size/C by default')
      b <- data_size/C
    }
  } else {
    if (b < 0) {
      stop("mesh_guidance_regular: b must be greater than 0")
    }
  }
  if (is.null(threshold)) {
    if (is.null(k3)) {
      stop("mesh_guidance_regular: if threshold is not passed, k3 must be passed in")
    } else if (is.null(k4)) {
      stop("mesh_guidance_regular: if threshold is not passed, k4 must be passed in")
    }
  } else if (is.numeric(threshold)) {
    if ((threshold < 0) | (threshold > 1)) {
      stop("mesh_guidance_regular: threshold must be between 0 and 1")
    }
    if (!is.numeric(max_E_nu_j)) {
      stop("mesh_guidance_regular: if threshold is not NULL, max_E_nu_j must be a numeric")
    }
    i <- 1
    trial_k3 <- c(-log(.Machine$double.xmin)+1, 100, 10, seq(5, trial_k3_by, -trial_k3_by))
    k3 <- trial_k3[i]
    k4 <- -log(threshold)-k3
    # find the first k3 in trial_k3 such that k4 > 0
    while (k4 < 0) {
      i <- i+1
      k3 <- trial_k3[i]
      k4 <- -log(threshold)-k3
    }
    # compute T1 and T2 for current k3 and k4
    T1 <- mesh_T1(k3 = k3, b = b, C = C, E_nu_j = max_E_nu_j, data_size = data_size)
    T2 <- mesh_T2(k4 = k4, b = b, C = C, data_size = data_size, d = d)
    # print(paste('k3:', k3))
    # print(paste('k4:', k4))
    # print(paste('exp(-k3-k4):', exp(-k3-k4)))
    # print(paste('T1:', T1))
    # print(paste('T2:', T2))
    while (T1 < T2) {
      # loop through possible k3, and compute the corresponding k4 such that the
      # user-specified lower bound is satisfied
      if (i==length(trial_k3)) {
        trial_k3[i+1] <- trial_k3[i]/10
      }
      i <- i+1
      k3 <- trial_k3[i]
      k4 <- -log(threshold)-k3
      T1 <- mesh_T1(k3 = k3, b = b, C = C, E_nu_j = max_E_nu_j, data_size = data_size)
      T2 <- mesh_T2(k4 = k4, b = b, C = C, data_size = data_size, d = d)
      # print(paste('k3:', k3))
      # print(paste('k4:', k4))
      # print(paste('T1:', T1))
      # print(paste('T2:', T2))
      if (T1 == 0) {
        # occurs when k3 is becoming too small and so T1 is 0
        # we cannot make T2 less than 0, so we haven't been able to find k3, k4
        # this can be solved sometimes by making the sequence of numbers that are
        # trialed for k3 finer (i.e. making trial_k3_by smaller - it is 0.00001 by default)
        stop("mesh_guidance_regular: couldn't find a k3 and k4. Try again with a smaller value for trial_k3_by (default is 0.00001)")
      }
    }
    # If we are able to find a suitable k3 and k4 that satisfies the bound,
    # (the previous loop tries finds the first small enough k3 which where T1 > T2)
    # we can continue to push k3 smaller and k4 larger while T1 > T2 to get a
    # computationally more efficient algorithm, since we'd want k4 as large as
    # possible until we cannot guarantee it will be smaller than T1 on average
    while (T1 > T2) {
      # print('pushing k3 lower')
      if (i==length(trial_k3)) {
        trial_k3[i+1] <- trial_k3[i]/10
      }
      i <- i+1
      k3 <- trial_k3[i]
      k4 <- -log(threshold)-k3
      T1 <- mesh_T1(k3 = k3, b = b, C = C, E_nu_j = max_E_nu_j, data_size = data_size)
      T2 <- mesh_T2(k4 = k4, b = b, C = C, data_size = data_size, d = d)
      # print(paste('k3:', k3))
      # print(paste('k4:', k4))
      # print(paste('T1:', T1))
      # print(paste('T2:', T2))
      if (T1 < T2) {
        k3 <- trial_k3[i-1]
        k4 <- -log(threshold)-k3
      }
    }
    warning("mesh_guidance_regular: k3 set to ", k3, " and k4 set to ", k4, " || exp(-k3-k4) = ", exp(-k3-k4))
  } else {
    stop("mesh_guidance_regular: threshold either is NULL or is a numeric value between 0 and 1")
  }
  return(list('max_delta_j' = mesh_T2(k4 = k4, b = b, C = C, data_size = data_size, d = d),
              'CESS_j_treshold' = exp(-k3-k4),
              'delta_j_T1' = mesh_T1(k3 = k3, b = b, C = C, E_nu_j = max_E_nu_j, data_size = data_size),
              'delta_j_T2' = mesh_T2(k4 = k4, b = b, C = C, data_size = data_size, d = d),
              'delta_j_T2<delta_j_T1' = mesh_T2(k4 = k4, b = b, C = C, data_size = data_size, d = d)<
                mesh_T1(k3 = k3, b = b, C = C, E_nu_j = max_E_nu_j, data_size = data_size),
              'k3' = k3,
              'k4' = k4,
              'max_E_nu_j' = max_E_nu_j))
}

#' @export
BF_guidance <- function(condition,
                        CESS_0_threshold = NULL,
                        CESS_j_threshold = NULL,
                        sub_posterior_samples,
                        log_weights = NULL,
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
                        trial_k3_by = 0.000001,
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
  if (d == 1) {
    N <- min(sapply(sub_posterior_samples, length))
    if (is.null(log_weights)) {
      log_weights <- rep(list(rep(log(1/N), N)), C)
    }
    max_E_nu_j <- compute_max_E_nu_j_univariate(N = N,
                                                sub_posterior_samples = sub_posterior_samples,
                                                log_weights = log_weights,
                                                time = T_guide$min_T,
                                                sub_posterior_means = sub_posterior_means,
                                                precondition_values = T_guide$precondition_matrices)
  } else {
    N <- min(sapply(sub_posterior_samples, nrow))
    if (is.null(log_weights)) {
      log_weights <- rep(list(rep(log(1/N), N)), C)
    }
    max_E_nu_j <- compute_max_E_nu_j_multivariate(N = N,
                                                  dim = d,
                                                  sub_posterior_samples = sub_posterior_samples,
                                                  log_weights = log_weights,
                                                  time = T_guide$min_T,
                                                  sub_posterior_means = sub_posterior_means,
                                                  inv_precondition_matrices = T_guide$inv_precondition_matrices,
                                                  Lambda = T_guide$Lambda)
  }
  mesh_guide <- mesh_guidance_regular(C = C,
                                      d = d,
                                      data_size = data_size,
                                      b = b,
                                      threshold = CESS_j_threshold,
                                      k3 = k3,
                                      k4 = k4,
                                      max_E_nu_j = max_E_nu_j$sumed,
                                      trial_k3_by = trial_k3_by,
                                      vanilla = vanilla)
  rec_mesh <- seq(from = 0, to = T_guide$min_T, by = mesh_guide$max_delta_j)
  if (rec_mesh[length(rec_mesh)]!=T_guide$min_T) {
    rec_mesh[length(rec_mesh)+1] <- T_guide$min_T
  }
  return(c(T_guide, mesh_guide, list('n' = length(rec_mesh), 'mesh' = rec_mesh)))
}
