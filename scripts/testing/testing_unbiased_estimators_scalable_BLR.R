cv_list = control_variates_BLR(dim = 3,
                               data = data,
                               prior_means = rep(0, 3),
                               prior_variances = rep(1, 3),
                               C = 1,
                               precondition_mat = cov(full_posterior))

test_scalable_phi <- function(beta, 
                              y_labels,
                              X,
                              prior_means,
                              prior_variances,
                              C,
                              precondition_mat,
                              transform_mats,
                              cv_list,
                              dim,
                              n_evaluate) {
  phi <- ea_phi_BLR_DL(beta = beta,
                       y_labels = y_labels,
                       X = X,
                       prior_means = prior_means,
                       prior_variances = prior_variances,
                       C = C,
                       precondition_mat = precondition_mat,
                       transform_mat = diag(1, dim))
  phi_transformed <- ea_phi_BLR_DL(beta = as.vector(transform_mats$to_Z %*% beta),
                                   y_labels = y_labels,
                                   X = X,
                                   prior_means = prior_means,
                                   prior_variances = prior_variances,
                                   C = C,
                                   precondition_mat = precondition_mat,
                                   transform_mat = transform_mats$to_X)
  scalable_phi <- sapply(1:n_evaluate, function(i) ea_phi_BLR_DL_vec_scalable(cv_list = cv_list,
                                                                              beta = beta,
                                                                              y_labels = y_labels,
                                                                              X = X,
                                                                              prior_means = prior_means,
                                                                              prior_variances = prior_variances,
                                                                              C = C,
                                                                              precondition_mat = precondition_mat))
  print(paste('phi:', phi))
  print(paste('phi_transformed:', phi_transformed))
  print(paste('mean(scalable_phi)', mean(scalable_phi)))
  return(list('phi' = phi,
              'phi_transformed' = phi_transformed,
              'scalable_phi' = scalable_phi))
}

test_phi_1 <- test_scalable_phi(beta = full_posterior[1,], 
                                y_labels = data$y,
                                X = data$X,
                                prior_means = rep(0, 3), 
                                prior_variances = rep(1, 3),
                                C = 1,
                                precondition_mat = cov(full_posterior),
                                transform_mats = list('to_Z' = expm::sqrtm(solve(cov(full_posterior))),
                                                      'to_X' = expm::sqrtm(cov(full_posterior))),
                                cv_list = cv_list,
                                dim = 3,
                                n_evaluate = 100000)

test_phi_1 <- test_scalable_phi(beta = full_posterior[50,], 
                                y_labels = data$y,
                                X = data$X,
                                prior_means = rep(0, 3), 
                                prior_variances = rep(1, 3),
                                C = 1,
                                precondition_mat = cov(full_posterior),
                                transform_mats = list('to_Z' = expm::sqrtm(solve(cov(full_posterior))),
                                                      'to_X' = expm::sqrtm(cov(full_posterior))),
                                cv_list = cv_list,
                                dim = 3,
                                n_evaluate = 100000)

test_scalable_alpha <- function(beta,
                                y_labels,
                                X,
                                prior_means,
                                prior_variances,
                                C,
                                cv_list,
                                n_evaluate) {
  X_beta <- X %*% beta
  X_beta_hat <- X %*% cv_list$beta_hat
  log_beta <- log_BLR_gradient(beta = beta,
                               y_labels = y_labels,
                               X = X,
                               X_beta = X_beta,
                               prior_means = prior_means,
                               prior_variances = prior_variances,
                               C = C) 
  log_beta_hat <- log_BLR_gradient(beta = cv_list$beta_hat,
                                   y_labels = y_labels,
                                   X = X,
                                   X_beta = X_beta_hat,
                                   prior_means = prior_means,
                                   prior_variances = prior_variances,
                                   C = C) 
  alpha <- log_beta - log_beta_hat
  scalable_alpha <- sapply(1:n_evaluate, function(i) alpha_tilde(index = sample(0:(cv_list$data_size-1), 1),
                                                                 beta = beta,
                                                                 beta_hat = cv_list$beta_hat,
                                                                 y_labels = y_labels,
                                                                 X = X,
                                                                 data_size = cv_list$data_size,
                                                                 prior_means = prior_means,
                                                                 prior_variances = prior_variances,
                                                                 C = C))
  print('alpha:'); print(as.vector(alpha))
  print('mean(scalable_alpha):'); print(apply(scalable_alpha, 1, mean))
  return(list('alpha' = alpha,
              'scalable_alpha' = scalable_alpha))
}

test_alpha_1 <- test_scalable_alpha(beta = full_posterior[1,],
                                    y_labels = data$y,
                                    X = data$X,
                                    prior_means = rep(0, 3),
                                    prior_variances = rep(1, 3),
                                    C = 1,
                                    cv_list = cv_list,
                                    n_evaluate = 100000)

test_alpha_2 <- test_scalable_alpha(beta = full_posterior[50,],
                                    y_labels = data$y,
                                    X = data$X,
                                    prior_means = rep(0, 3),
                                    prior_variances = rep(1, 3),
                                    C = 1,
                                    cv_list = cv_list,
                                    n_evaluate = 100000)

test_scalable_div_alpha <- function(beta,
                                    X,
                                    prior_variances,
                                    C,
                                    precondition_mat,
                                    cv_list,
                                    n_evaluate) {
  X_beta <- X %*% beta
  X_beta_hat <- X %*% cv_list$beta_hat
  div_log_beta <- div_log_BLR_gradient(X = X,
                                       X_beta = X_beta,
                                       prior_variances = prior_variances,
                                       C = C,
                                       precondition_mat = precondition_mat)
  div_log_beta_hat <- div_log_BLR_gradient(X = X,
                                           X_beta = X_beta_hat,
                                           prior_variances = prior_variances,
                                           C = C,
                                           precondition_mat = precondition_mat) 
  div_alpha <- div_log_beta - div_log_beta_hat
  scalable_div_alpha <- sapply(1:n_evaluate, function(i) div_alpha_tilde(index = sample(0:(cv_list$data_size-1), 1),
                                                                         beta = beta,
                                                                         beta_hat = cv_list$beta_hat,
                                                                         X = X,
                                                                         data_size = cv_list$data_size,
                                                                         prior_variances = prior_variances,
                                                                         C = C,
                                                                         precondition_mat = precondition_mat))
  print(paste('div_alpha:', div_alpha))
  print(paste('mean(scalable_div_alpha):', mean(scalable_div_alpha)))
  return(list('div_alpha' = div_alpha,
              'scalable_div_alpha' = scalable_div_alpha))
}

test_div_alpha_1 <- test_scalable_div_alpha(beta = full_posterior[1,],
                                            X = data$X,
                                            prior_variances = rep(1, 3),
                                            C = 1,
                                            precondition_mat = precondition_mat,
                                            cv_list = cv_list,
                                            n_evaluate = 100000)

test_div_alpha_2 <- test_scalable_div_alpha(beta = full_posterior[50,],
                                            X = data$X,
                                            prior_variances = rep(1, 3),
                                            C = 1,
                                            precondition_mat = precondition_mat,
                                            cv_list = cv_list,
                                            n_evaluate = 100000)
