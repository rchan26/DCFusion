library(hierarchicalFusion)

beta <- full_posterior[15,]

ea_phi_BLR_DL_vec(beta = beta,
                  y_labels = data$y,
                  X = data$X,
                  prior_means = rep(0, 3),
                  prior_variances = rep(1, 3),
                  C = 1,
                  precondition_mat = cov(full_posterior),
                  transform_mat = diag(1, 3))
ea_phi_BLR_DL_vec_Z(beta = beta,
                    y_labels = data$y,
                    X = data$X, 
                    prior_means = rep(0, 3),
                    prior_variances = rep(1, 3),
                    C = 1,
                    precondition_mat = cov(full_posterior),
                    transform_mat = expm::sqrtm(cov(full_posterior)))

# equal with identity matrix
ea_phi_BLR_DL_vec(beta = beta,
                  y_labels = data$y,
                  X = data$X,
                  prior_means = rep(0, 3),
                  prior_variances = rep(1, 3),
                  C = 1,
                  precondition_mat = diag(1, 3),
                  transform_mat = diag(1, 3))
ea_phi_BLR_DL_vec_Z(beta = beta,
                    y_labels = data$y,
                    X = data$X,
                    prior_means = rep(0, 3),
                    prior_variances = rep(1, 3),
                    C = 1,
                    precondition_mat = diag(1, 3),
                    transform_mat = diag(1, 3))

# equal with covariance preconditioning matrix
grad <- log_BLR_gradient(beta = beta,
                         y_labels = data$y,
                         X = data$X,
                         X_beta = data$X %*% beta,
                         prior_means = rep(0, 3), 
                         prior_variances = rep(1, 3),
                         C = 1)
t(grad) %*% cov(full_posterior) %*% grad
grad_Z <- log_BLR_gradient_Z(beta = beta,
                             y_labels = data$y,
                             transformed_X = data$X %*% expm::sqrtm(cov(full_posterior)),
                             X_beta = data$X %*% beta,
                             prior_means = rep(0, 3), 
                             prior_variances = rep(1, 3),
                             C = 1,
                             precondition_mat = cov(full_posterior),
                             transform_mat = expm::sqrtm(cov(full_posterior)))
sum(grad_Z^2)

# equal with identity matrix
grad <- log_BLR_gradient(beta = beta,
                         y_labels = data$y,
                         X = data$X,
                         X_beta = data$X %*% beta,
                         prior_means = rep(0, 3), 
                         prior_variances = rep(1, 3),
                         C = 1)
t(grad) %*% diag(1,3) %*% grad
grad_Z <- log_BLR_gradient_Z(beta = beta,
                             y_labels = data$y,
                             transformed_X = data$X %*% diag(1, 3),
                             X_beta = data$X %*% beta,
                             prior_means = rep(0, 3), 
                             prior_variances = rep(1, 3),
                             C = 1,
                             precondition_mat = diag(1, 3),
                             transform_mat = diag(1, 3))
sum(grad_Z^2)

# not equal with covariance preconditioning matrix
div_log_BLR_gradient(X = data$X,
                     X_beta = data$X %*% beta,
                     prior_variances = rep(1, 3),
                     C = 1,
                     precondition_mat = cov(full_posterior))
div_log_BLR_gradient_Z(X_beta = data$X %*% beta,
                       transformed_X = data$X %*% expm::sqrtm(cov(full_posterior)),
                       prior_variances = rep(1, 3),
                       C = 1,
                       precondition_mat = cov(full_posterior),
                       transform_mat = expm::sqrtm(cov(full_posterior)))
term2(X = data$X,
      X_beta = data$X %*% beta,
      prior_variances = rep(1, 3),
      C = 1,
      precondition_mat = cov(full_posterior))

# equal with identity matrix
div_log_BLR_gradient(X = data$X,
                     X_beta = data$X %*% beta,
                     prior_variances = rep(1, 3),
                     C = 1,
                     precondition_mat = diag(1, 3))
div_log_BLR_gradient_Z(X_beta = data$X %*% beta,
                       transformed_X = data$X %*% diag(1, 3),
                       prior_variances = rep(1, 3),
                       C = 1,
                       precondition_mat = diag(1, 3),
                       transform_mat = diag(1, 3))
term2(X = data$X,
      X_beta = data$X %*% beta,
      prior_variances = rep(1, 3),
      C = 1,
      precondition_mat = diag(1, 3))

