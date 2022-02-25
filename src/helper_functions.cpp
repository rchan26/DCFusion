#include "../inc/helper_functions.hpp"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

double find_max(const Rcpp::NumericVector &vect) {
  double current_max = vect.at(0);
  for (const auto &element: vect) {
    if (element > current_max) {
      current_max = element;
    }
  }
  return current_max;
}

double find_min(const Rcpp::NumericVector &vect) {
  double current_min = vect.at(0);
  for (const auto &element: vect) {
    if (element < current_min) {
      current_min = element;
    }
  }
  return current_min;
}

//' Calculate the weighted mean (univariate)
//' 
//' Calculation of weighted mean when the target is univariate
//'
//' @param x vector of values
//' @param weights vector of weights
//'
//' @return the weighted mean of x
//' 
//' @examples
//' x <- c(0.4, 0.2, 0.5, 0.9, 1.4)
//' w <- c(3, 2, 5, 1, 2)
//' weighted_mean_univariate(x = x,
//'                          weights = w)
//' # returns the same using weighted.mean(x, w) function in R
//' weighted.mean(x, w)
//'
//' weighted_mean_univariate(x = x,
//'                          weights = rep(1, 5))
//' #returns the same using standard mean(x) function in R
//' mean(x)
// [[Rcpp::export]]
double weighted_mean_univariate(const Rcpp::NumericVector &x,
                                const Rcpp::NumericVector &weights) {
  if (x.size() != weights.size()) {
    stop("weighted_mean: x and weights must be of the same size");
  }
  double w_mean = 0.0;
  for (int i=0; i < x.size(); ++i) {
    w_mean += weights[i]*x[i];
  }
  return w_mean / Rcpp::sum(weights);
}

//' Calculate the logarithm of rho (univariate)
//' 
//' Calculation of the log of rho acceptance probability or weight when target is
//' univariate
//'
//' @param x vector of sampled sub-posterior values
//' @param x_mean weighted mean of sampled sub-posterior values
//' @param time time T for fusion algorithm
//' @param precondition_values precondition values associated to each sub-posterior
//'
//' @return the logarithm of rho
//' 
//' @examples
//' x <- rnorm(4, 0, 1)
//' precondition_vals <- c(1, 2, 3, 4)
//' x_mean <- weighted_mean_univariate(x = x,
//'                                    weights = 1/precondition_vals)
//' log_rho_univariate(x = x,
//'                    x_mean = x_mean,
//'                    time = 0.5,
//'                    precondition_values = precondition_vals)
// [[Rcpp::export]]
double log_rho_univariate(const Rcpp::NumericVector &x,
                          const double &x_mean,
                          const double &time,
                          const Rcpp::NumericVector &precondition_values) {
  if (x.size() != precondition_values.size()) {
    stop("log_rho_univariate: x and precondition_values must be of the same size");
  }
  double numerator = 0.0;
  for (int i=0; i < x.size(); ++i) {
    double D = (x_mean - x[i]);
    double num = D*D/precondition_values[i];
    numerator += num;
  }
  return (-numerator/(2*time));
}

//' Calculate the variance of numbers (univariate)
//' 
//' Calculation of the weighted variance of numbers
//'
//' @param x vector of sampled sub-posterior values
//' @param x_mean weighted mean of sampled sub-posterior values
//' @param precondition_values precondition values associated to each sub-posterior
//'
//' @return the weighted variance of numbers
//'
//' @examples
//' x <- rnorm(4, 0, 1)
//' precondition_vals <- c(1, 2, 3, 4)
//' x_mean <- weighted_mean_univariate(x = x,
//'                                    weights = 1/precondition_vals)
//' weighted_variance_univariate(x = x,
//'                              x_mean = x_mean,
//'                              precondition_values = precondition_vals)
// [[Rcpp::export]]
double weighted_variance_univariate(const arma::vec &x,
                                    const double &x_mean,
                                    const arma::vec &precondition_values) {
  if (x.size() != precondition_values.size()) {
    stop("weighted_variance_univariate: x and precondition_values must be of the same size");
  }
  double w_variance = 0.0;
  for (int i=0; i < x.size(); ++i) {
    double D = x_mean - x[i];
    double num = D*D/precondition_values[i];
    w_variance += num;
  }
  return (w_variance/x.size());
}

//' Calculate approximation to expectation of nu_j (univariate)
//' 
//' Calculation of the scaled/weighted average variation of the C trajectories
//' with respect to their individual sub-posterior means
//'
//' @param list where x_samples[[i]] ith collection of the C trajectories
//' @param sub_posterior_means vector of length C of sub-posterior means
//' @param precondition_values precondition values associated to each sub-posterior
//'
//' @return the approximated expectation of nu_j
//'
//' @examples
//' # x_samples has 5 samples and C=4
//' N <- 10
//' C <- 4
//' x_samples <- lapply(1:N, function(i) rnorm(C))
//' normalised_weights <- rep(1/N, N)
//' sub_posterior_means <- rnorm(C)
//' precond <- 1:C
//' weighted_trajectory_variation_univariate(x_samples = x_samples,
//'                                          normalised_weights = normalised_weights,
//'                                          sub_posterior_means = sub_posterior_means,
//'                                          precondition_values = precond)
//' # should be equal to the result of this:
//' sum(sapply(1:N, function(i) {
//'   sum((((x_samples[[i]]-sub_posterior_means)^2)/precond))/C
//' }))/N
// [[Rcpp::export]]
double weighted_trajectory_variation_univariate(const Rcpp::List &x_samples,
                                                const arma::vec &normalised_weights,
                                                const arma::vec &sub_posterior_means,
                                                const arma::vec &precondition_values) {
  const double C = precondition_values.size();
  if (sub_posterior_means.size()!=C) {
    stop("weighted_trajectory_variation_univariate: sub_posterior_means must be a vector of length C=length(precondition_values)");
  } else if (x_samples.size()!=normalised_weights.size()) {
    stop("weighted_trajectory_variation_univariate: x_samples and normalised_weights must have the same length");
  }
  double variation = 0.0;
  for (int i=0; i < x_samples.size(); ++i) {
    const arma::vec &x_i = x_samples[i];
    if (x_i.size()!=C) {
      stop("weighted_trajectory_variation_univariate: x_samples[i] must be a list of vectors of length C=length(precondition_values)");
    }
    double ith_variation = 0.0;
    for (int c=0; c < C; ++c) {
      const double D = x_i[c] - sub_posterior_means[c];
      ith_variation += D*D/precondition_values[c];
    }
    variation += normalised_weights.at(i)*ith_variation/C;
  }
  return variation;
}

// [[Rcpp::export]]
Rcpp::List compute_max_E_nu_j_univariate(const double &N,
                                         const Rcpp::List &sub_posterior_samples,
                                         const Rcpp::List &log_weights,
                                         const double &time,
                                         const arma::vec &sub_posterior_means,
                                         const Rcpp::NumericVector &precondition_values) {
  const double C = precondition_values.size();
  if (sub_posterior_samples.size()!=C) {
    stop("compute_max_E_nu_j_univariate: sub_posterior_samples must be a list of length C=length(precondition_values)");
  } else if (sub_posterior_means.size()!=C) {
    stop("compute_max_E_nu_j_univariate: sub_posterior_means must be a vector of length C=length(precondition_values)");
  } else if (log_weights.size()!=C) {
    stop("compute_max_E_nu_j_univariate: log_weights must be a list of length C=length(precondition_values)");
  }
  // computing the weights for the initialisation step
  Rcpp::NumericVector x_means(N);
  Rcpp::NumericVector log_rho_weights(N);
  for (int i=0; i < N; ++i) {
    Rcpp::checkUserInterrupt();
    Rcpp::NumericVector particle(C);
    for (int c=0; c < C; ++c) {
      const Rcpp::NumericVector &sub_post_samples = sub_posterior_samples[c];
      particle[c] = sub_post_samples.at(i);
      const Rcpp::NumericVector &lw = log_weights[c];
      log_rho_weights[i] += lw.at(i);
    }
    x_means[i] = weighted_mean_univariate(particle, 1/precondition_values);
    log_rho_weights[i] += log_rho_univariate(particle, x_means[i], time, precondition_values);
  }
  const Rcpp::List ESS = particle_ESS(log_rho_weights);
  const Rcpp::NumericVector &normalised_weights = ESS["normalised_weights"];
  // using those normalised weights to compute the maximum E[nu_j]
  Rcpp::NumericVector variation(N);
  for (int i=0; i < N; ++i) {
    Rcpp::checkUserInterrupt();
    double ith_variation = 0.0;
    for (int c=0; c < C; ++c) {
      const double D = x_means[i] - sub_posterior_means[c];
      ith_variation += D*D/precondition_values[c];
    }
    variation.at(i) += normalised_weights.at(i)*ith_variation/C;
  }
  const double sumed = Rcpp::sum(variation);
  const double maxed = N*Rcpp::max(variation);
  return Rcpp::List::create(Named("sumed", sumed), Named("maxed", maxed));
}

//' Calculate the inverse of a sum of matrices
//'
//' Calculation of the inverse of a sum of a list of matrices
//'
//' @param matrices list of matrices (of same dimension)
//'
//' @return the inverse of the sum of the matrices
//'
//' @examples
//' m1 <- matrix(c(1,2,3,4), nrow = 2, ncol = 2)
//' m2 <- matrix(c(5,6,7,8), nrow = 2, ncol = 2)
//' m3 <- matrix(c(9,10,11,12), nrow = 2, ncol = 2)
//' inverse_sum_matrices(list(m1, m2, m3))
//' # returns the same as using solve() in R
//' solve(m1+m2+m3)
// [[Rcpp::export]]
arma::mat inverse_sum_matrices(const Rcpp::List &matrices) {
  const arma::mat &matrix_1 = matrices[0];
  arma::mat matrix_sum(matrix_1.n_cols, matrix_1.n_cols, arma::fill::zeros);
  for (int c=0; c < matrices.size(); ++c) {
    const arma::mat &matrix = matrices[c];
    matrix_sum += matrix;
  }
  return(arma::inv(matrix_sum));
}

//' Calculate the weighted mean (multivariate)
//'
//' Calculation of weighted mean when the target is multivariate
//' 
//' @param matrix an (m x n) matrix where the ith row is the ith sample
//' @param weights list of m matrices of the same dimension (n x n)
//' @param inverse_sum_weights the inverse of the sum of the weights (can be 
//'                        calculated by passing in weights to inverse_sum_matrices)
//'
//' @return proposal mean vector
//' 
//' @examples
//' m1 <- matrix(c(1,2,3,4), nrow = 2, ncol = 2)
//' m2 <- matrix(c(5,6,7,8), nrow = 2, ncol = 2)
//' m3 <- matrix(c(9,10,11,12), nrow = 2, ncol = 2)
//' X <- matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2)
//' weighted_mean_multivariate(matrix = X,
//'                            weights = list(m1, m2, m3),
//'                            inverse_sum_weights = inverse_sum_matrices(list(m1, m2, m3)))
// [[Rcpp::export]]
arma::vec weighted_mean_multivariate(const arma::mat &matrix,
                                     const Rcpp::List &weights,
                                     const arma::mat &inverse_sum_weights) {
  arma::vec weighted_mean(matrix.n_cols, arma::fill::zeros);
  for (int c=0; c < matrix.n_rows; ++c) {
    const arma::mat &weight = weights[c];
    weighted_mean += weight*arma::trans(matrix.row(c));
  }
  return(inverse_sum_weights*weighted_mean);
}

//' Calculate the proposal covariance matrix
//' 
//' Calculation of the proposal covariance matrix for Monte Carlo Fusion
//'
//' @param time T for fusion algorithm
//' @param weights list of m matrices of the same dimension
//'
//' @return proposal covariance matrix
//' 
//' @examples
//' m1 <- matrix(c(1,2,3,4), nrow = 2, ncol = 2)
//' m2 <- matrix(c(5,6,7,8), nrow = 2, ncol = 2)
//' m3 <- matrix(c(9,10,11,12), nrow = 2, ncol = 2)
//' calculate_proposal_cov(time = 0.5, weights = list(m1, m2, m3))
// [[Rcpp::export]]
arma::mat calculate_proposal_cov(const double &time,
                                 const Rcpp::List &weights) {
  const arma::mat &matrix_1 = weights[0];
  arma::mat lambda_inv(matrix_1.n_cols, matrix_1.n_cols, arma::fill::zeros);
  for (int c=0; c < weights.size(); ++c) {
    const arma::mat &weight = weights[c];
    lambda_inv += weight;
  }
  return(arma::inv(lambda_inv/time));
}

//' Row-wise subtraction of a vector to rows of a matrix
//' 
//' Calculates the subtraction of a vector to each row of a matrix
//'
//' @param X matrix (n x m)
//' @param vect vector of length m
//'
//' @return an (n x m) matrix, Y, where Y[i,] = X[i,]-vect
//' 
//' @examples
//' X <- matrix(c(1,2,3,4,5,6,7,8), nrow = 4, ncol = 2, byrow = T)
//' row_wise_subtraction(X = X, vect = c(1,2))
// [[Rcpp::export]]
arma::mat row_wise_subtraction(const arma::mat &X,
                               const arma::vec &vect) {
  arma::mat new_mat(size(X), arma::fill::zeros);
  const arma::rowvec row_vect = trans(vect);
  for (int row=0; row < X.n_rows; ++row) {
    new_mat.row(row) = X.row(row) - row_vect;
  }
  return(new_mat);
}

//' Calculate the logarithm of rho (multivariate)
//' 
//' Calculation of the log of rho acceptance probability or weight when target is
//' multivariate
//'
//' @param x an (m x n) matrix where the ith row is the ith sample
//' @param x_mean a vector of length n (the weighted mean of x samples)
//' @param time time T for fusion algorithm
//' @param inv_precondition_matrices list of length m of inverse 
//'                                  preconditioning matrices
//'
//' @return the logarithm of rho
//'
//' @examples
//' # set covariance matrices
//' Sig1 <- diag(2)
//' Sig2 <- matrix(c(2, 0.5, 0.5, 2), nrow = 2, ncol = 2)
//' Sig3 <- matrix(c(4, -3.2, -3.2, 4), nrow = 2, ncol = 2)
//' # sample some x values and store in the rows
//' x <- matrix(nrow = 3, ncol = 2)
//' x[1,] <- mvrnormArma(N = 1, mu = c(0, 0), Sigma = Sig1)
//' x[2,] <- mvrnormArma(N = 1, mu = c(0, 0), Sigma = Sig2)
//' x[3,] <- mvrnormArma(N = 1, mu = c(0, 0), Sigma = Sig3)
//' # calcualte precondition matrices and their inverses
//' precondition_matrices <- list(Sig1, Sig2, Sig3)
//' inv_precondition_matrices <- lapply(precondition_matrices, solve)
//' inverse_sum_weights <- inverse_sum_matrices(precondition_matrices)
//' # calculate the weighted mean where weights are the inverse precondition matrices
//' x_mean <- weighted_mean_multivariate(matrix = x,
//'                                      weights = precondition_matrices,
//'                                      inverse_sum_weights = inverse_sum_weights)
//' # calculate logarithm of rho with time T = 0.5
//' log_rho_multivariate(x = x,
//'                      x_mean = x_mean,
//'                      time = 0.5,
//'                      inv_precondition_matrices = inv_precondition_matrices)
// [[Rcpp::export]]
double log_rho_multivariate(const arma::mat &x,
                            const arma::vec &x_mean,
                            const double &time,
                            const Rcpp::List &inv_precondition_matrices) {
  if (x.n_cols != x_mean.size()) {
    stop("log_rho_multivariate: ncol(x) and x_mean.size() must be equal");
  } else if (x.n_rows != inv_precondition_matrices.size()) {
    stop("log_rho_multivariate: nrow(x) and length(inv_precondition_matrices) must be equal");
  }
  const arma::mat D = row_wise_subtraction(x, x_mean);
  double numerator = 0.0;
  for (int c=0; c < D.n_rows; ++c) {
    const arma::mat &inv_precond = inv_precondition_matrices[c];
    numerator += as_scalar((D.row(c)*inv_precond)*trans(D.row(c)));
  }
  return (-numerator/(2*time));
}

//' Calculate the variance of vectors (multivariate)
//' 
//' Calculation of the weighted variance of vectors
//'
//' @param x an (m x n) matrix where the ith row is the ith sample
//' @param x_mean a vector of length n (the weighted mean of x samples)
//' @param inv_precondition_matrices list of length m of inverse 
//'                                  preconditioning matrices
//'
//' @return the weighted variance of vectors
//'
//' @examples
//' # set covariance matrices
//' Sig1 <- diag(2)
//' Sig2 <- matrix(c(2, 0.5, 0.5, 2), nrow = 2, ncol = 2)
//' Sig3 <- matrix(c(4, -3.2, -3.2, 4), nrow = 2, ncol = 2)
//' # sample some x values and store in the rows
//' x <- matrix(nrow = 3, ncol = 2)
//' x[1,] <- mvrnormArma(N = 1, mu = c(0, 0), Sigma = Sig1)
//' x[2,] <- mvrnormArma(N = 1, mu = c(0, 0), Sigma = Sig2)
//' x[3,] <- mvrnormArma(N = 1, mu = c(0, 0), Sigma = Sig3)
//' # calculate precondition matrices and their inverses
//' precondition_matrices <- list(Sig1, Sig2, Sig3)
//' inv_precondition_matrices <- lapply(precondition_matrices, solve)
//' inverse_sum_weights <- inverse_sum_matrices(precondition_matrices)
//' # calculate the weighted mean where weights are the inverse precondition matrices
//' x_mean <- weighted_mean_multivariate(matrix = x,
//'                                      weights = precondition_matrices,
//'                                      inverse_sum_weights = inverse_sum_weights)
//' weighted_variance_multivariate(x = x,
//'                                x_mean = x_mean,
//'                                inv_precondition_matrices = inv_precondition_matrices)
// [[Rcpp::export]]
double weighted_variance_multivariate(const arma::mat &x,
                                      const arma::vec &x_mean,
                                      const Rcpp::List &inv_precondition_matrices) {
  if (x.n_cols != x_mean.size()) {
    stop("weighted_variance_multivariate: ncol(x) and x_mean.size() must be equal");
  } else if (x.n_rows != inv_precondition_matrices.size()) {
    stop("weighted_variance_multivariate: nrow(x) and length(inv_precondition_matrices) must be equal");
  }
  const arma::mat D = row_wise_subtraction(x, x_mean);
  double w_variance = 0.0;
  for (int c=0; c < D.n_rows; ++c) {
    const arma::mat &inv_precond = inv_precondition_matrices[c];
    w_variance += as_scalar((D.row(c)*inv_precond)*trans(D.row(c)));
  }
  return (w_variance/x.n_rows);
}

//' Calculate approximation to expectation of nu_j (multivariate)
//' 
//' Calculation of the scaled/weighted average variation of the C trajectories
//' with respect to their individual sub-posterior means
//'
//' @param list where x_samples[[i]] ith collection of the C trajectories
//' @param sub_posterior_means matrix with C rows of sub-posterior means
//' @param inv_precondition_matrices list of length m of inverse 
//'                                  preconditioning matrices
//'
//' @return the approximated expectation of nu_j
//'
//' @examples
//' N <- 10
//' C <- 4
//' d <- 3
//' x_samples <- lapply(1:N, function(i) mvrnormArma(C, rep(0,d), diag(1,d)))
//' normalised_weights <- rep(1/N, N)
//' sub_posterior_means <- mvrnormArma(C, rep(0,d), diag(1,d))
//' precond <- lapply(1:C, function(c) diag(c, d))
//' inv_precond <- lapply(precond, solve)
//' weighted_trajectory_variation_multivariate(x_samples = x_samples,
//'                                            normalised_weights = normalised_weights,
//'                                            sub_posterior_means = sub_posterior_means,
//'                                            inv_precondition_matrices = inv_precond)
//' # should be equal to the result of this:
//' sum(sapply(1:N, function(i) {
//'   sum(sapply(1:C, function(c) {
//'     diff <- x_samples[[i]][c,]-sub_posterior_means[c,]
//'     return(t(diff) %*% inv_precond[[c]] %*% diff)
//'   }))/C
//' }))/N
// [[Rcpp::export]]
double weighted_trajectory_variation_multivariate(const Rcpp::List &x_samples,
                                                  const arma::vec &normalised_weights,
                                                  const arma::mat &sub_posterior_means,
                                                  const Rcpp::List &inv_precondition_matrices) {
  const double C = inv_precondition_matrices.size();
  if (sub_posterior_means.n_rows!=C) {
    stop("weighted_trajectory_variation_multivariate: sub_posterior_means must be a matrix with C=length(inv_precondition_matrices) rows");
  } else if (x_samples.size()!=normalised_weights.size()) {
    stop("weighted_trajectory_variation_multivariate: x_samples and normalised_weights must have the same length");
  }
  double variation = 0.0;
  for (int i=0; i < x_samples.size(); ++i) {
    const arma::mat &x_i = x_samples[i];
    if (x_i.n_rows!=C) {
      stop("weighted_trajectory_variation_multivariate: x_samples must be a list of matrices with C=length(inv_precondition_matrices) rows");
    }
    double ith_variation = 0.0;
    for (int c=0; c < C; ++c) {
      const arma::rowvec &x_i_c = x_i.row(c);
      const arma::rowvec &a_c = sub_posterior_means.row(c);
      const arma::rowvec D = x_i_c - a_c;
      const arma::mat &inv_precond = inv_precondition_matrices[c];
      ith_variation += as_scalar((D*inv_precond)*trans(D));
    }
    variation += normalised_weights.at(i)*ith_variation/C;
  }
  return variation;
}

// [[Rcpp::export]]
Rcpp::List compute_max_E_nu_j_multivariate(const double &N,
                                           const int &dim,
                                           const Rcpp::List &sub_posterior_samples,
                                           const Rcpp::List &log_weights,
                                           const double &time,
                                           const arma::mat &sub_posterior_means,
                                           const Rcpp::List &inv_precondition_matrices,
                                           const arma::mat &Lambda) {
  const double C = inv_precondition_matrices.size();
  if (sub_posterior_samples.size()!=C) {
    stop("compute_max_E_nu_j_multivariate: sub_posterior_samples must be a list of length C=length(inv_precondition_matrices)");
  } if (sub_posterior_means.n_rows!=C) {
    stop("compute_max_E_nu_j_multivariate: sub_posterior_means must be a matrix with C=length(inv_precondition_matrices) rows");
  } else if (log_weights.size()!=C) {
    stop("compute_max_E_nu_j_multivariate: log_weights must be a list of length C=length(precondition_values)");
  }
  // computing the weights for the initialisation step
  arma::mat x_means(N, dim, arma::fill::zeros);
  Rcpp::NumericVector log_rho_weights(N);
  for (int i=0; i < N; ++i) {
    Rcpp::checkUserInterrupt();
    arma::mat particle(C, dim);
    for (int c=0; c < C; ++c) {
      const arma::mat &sub_post_samples = sub_posterior_samples[c];
      particle.row(c) = sub_post_samples.row(i);
      const Rcpp::NumericVector &lw = log_weights[c];
      log_rho_weights[i] += lw.at(i);
    }
    arma::vec particle_mean = weighted_mean_multivariate(particle,
                                                         inv_precondition_matrices,
                                                         Lambda);
    log_rho_weights[i] += log_rho_multivariate(particle,
                                               particle_mean,
                                               time,
                                               inv_precondition_matrices);
    x_means.row(i) = arma::trans(particle_mean);
  }
  const Rcpp::List ESS = particle_ESS(log_rho_weights);
  const Rcpp::NumericVector &normalised_weights = ESS["normalised_weights"];
  // using those normalised weights to compute the maximum E[nu_j]
  Rcpp::NumericVector variation(N);
  for (int i=0; i < N; ++i) {
    Rcpp::checkUserInterrupt();
    double ith_variation = 0.0;
    for (int c=0; c < C; ++c) {
      const arma::rowvec &a_c = sub_posterior_means.row(c);
      const arma::rowvec D = x_means.row(i) - a_c;
      const arma::mat &inv_precond = inv_precondition_matrices[c];
      ith_variation += as_scalar((D*inv_precond)*trans(D));
    }
    variation.at(i) += normalised_weights.at(i)*ith_variation/C;
  }
  const double sumed = Rcpp::sum(variation);
  const double maxed = N*Rcpp::max(variation);
  return Rcpp::List::create(Named("sumed", sumed), Named("maxed", maxed));
}

//' Calculate the logarithm of the sum of the exponentials of the arguments
//' 
//' Calculation of the log of the sum of exponential of x, but avoids computational 
//' underflow / overflow
//' 
//' @param x vector
//'
//' @return the logarithm of the sum of the exponentials of elements of x
//'         i.e. returns log(sum(exp(x)))
//'
//' @examples
//' x <- c(1000.01, 1000.02)
//' y1 <- log(sum(exp(x)))
//' print(y1)
//' y2 <- logsumexp(x)
//' print(y2) 
// [[Rcpp::export]]
double logsumexp(const Rcpp::NumericVector &x) {
  const double x_star = Rcpp::max(x);
  return (x_star+log(sum(exp(x-x_star))));
}

//' Calculate the ESS
//'
//' Calculates the ESS and the normalised weights from the logarithm of weihgts
//'
//' @param log_weights the logarithm of the particle weights
//'
//' @return A list with components
//' \describe{
//'  \item{log_weights}{the logarithm of the particle weights}
//'  \item{normalised_weights}{the normalised particle weights}
//'  \item{ESS}{the effective sample size of particles (the inverse of the
//'             sum of the weights squared)}
//' }
//'
//' @examples
//' particle_ESS(log_weights = c(1000.01, 1000.02, 1000.03))
// [[Rcpp::export]]
Rcpp::List particle_ESS(const Rcpp::NumericVector &log_weights) {
  const Rcpp::NumericVector normalised_weights = exp(log_weights-logsumexp(log_weights));
  return Rcpp::List::create(Named("log_weights", Rcpp::log(normalised_weights)),
                            Named("normalised_weights", normalised_weights),
                            Named("ESS", 1/Rcpp::sum(Rcpp::pow(normalised_weights, 2))),
                            Named("input", log_weights));
}

// [[Rcpp::export]]
Rcpp::List rho_IS_univariate_(const Rcpp::List &samples_to_fuse,
                              const int &N,
                              const int &m,
                              const double &time,
                              const Rcpp::NumericVector &precondition_values) {
  Rcpp::NumericVector x(m);
  Rcpp::NumericVector x_means(N);
  Rcpp::List x_samples = rep(Rcpp::List::create(x), N);
  Rcpp::NumericVector log_rho_weights(N);
  for (int i=0; i < N; ++i) {
    Rcpp::checkUserInterrupt();
    Rcpp::NumericVector particle(m);
    for (int c=0; c < m; ++c) {
      const Rcpp::NumericVector &sub_post_samples = samples_to_fuse[c];
      particle[c] = sub_post_samples.at(i);
    }
    x_means[i] = weighted_mean_univariate(particle, 1/precondition_values);
    log_rho_weights[i] = log_rho_univariate(particle, x_means[i], time, precondition_values);
    x_samples[i] = particle;
  }
  return Rcpp::List::create(Named("x_samples", x_samples),
                            Named("x_means", x_means),
                            Named("log_weights", log_rho_weights));
}

// [[Rcpp::export]]
Rcpp::List rho_IS_multivariate_(const Rcpp::List &samples_to_fuse,
                                const int &dim,
                                const int &N,
                                const int &m,
                                const double &time,
                                const Rcpp::List &inv_precondition_matrices,
                                const arma::mat &inverse_sum_inv_precondition_matrices) {
  arma::mat x(m, dim, arma::fill::zeros);
  arma::mat x_means(N, dim, arma::fill::zeros);
  Rcpp::List x_samples = rep(Rcpp::List::create(x), N);
  Rcpp::NumericVector log_rho_weights(N);
  for (int i=0; i < N; ++i) {
    Rcpp::checkUserInterrupt();
    arma::mat particle(m, dim);
    for (int c=0; c < m; ++c) {
      const arma::mat &sub_post_samples = samples_to_fuse[c];
      particle.row(c) = sub_post_samples.row(i);
    }
    arma::vec particle_mean = weighted_mean_multivariate(particle,
                                                         inv_precondition_matrices,
                                                         inverse_sum_inv_precondition_matrices);
    log_rho_weights[i] = log_rho_multivariate(particle,
                                              particle_mean,
                                              time,
                                              inv_precondition_matrices);
    x_samples[i] = particle;
    x_means.row(i) = arma::trans(particle_mean);
  }
  return Rcpp::List::create(Named("x_samples", x_samples),
                            Named("x_means", x_means),
                            Named("log_weights", log_rho_weights));
}

//' Simulate from a Multivariate Gaussian Distribution
//'
//' Produces samples from the specified multivariate Gaussian distribution
//'
//' @param N the number of samples required
//' @param mu a vector of mean values
//' @param Sigma positive-definite symmetric matrix specifying the covariance of variables
//'
//' @return samples from a multivariate Gaussian distribution
//' 
//' @examples
//' samples <- mvrnormArma(N = 10000, mu = c(0, 0), Sigma = diag(2))
// [[Rcpp::export]]
arma::mat mvrnormArma(const int &N,
                      const arma::vec &mu,
                      const arma::mat &Sigma) {
  const arma::mat Y = arma::randn(N, Sigma.n_cols);
  return arma::repmat(mu, 1, N).t() + Y * arma::chol(Sigma);
}

//' Simulate from a tempered Multivariate Gaussian Distribution
//'
//' Produces samples from the specified tempered multivariate Gaussian distribution
//'
//' @param N the number of samples required
//' @param mu a vector of mean values
//' @param Sigma positive-definite symmetric matrix specifying the covariance of variables
//' @param beta inverse temperature  
//'
//' @return samples from a tempered multivariate Gaussian distribution
//'
//' @examples
//' samples <- mvrnormArma(N = 10000,
//'                        mu = c(0, 0),
//'                        Sigma = diag(2))
//' tempered_samples <- mvrnormArma_tempered(N = 10000,
//'                                          mu = c(0, 0),
//'                                          Sigma = diag(2),
//'                                          beta = 2)
// [[Rcpp::export]]
arma::mat mvrnormArma_tempered(const int &N,
                               const arma::vec &mu,
                               const arma::mat &Sigma,
                               const double &beta) {
  const arma::mat new_Sigma = Sigma/beta;
  const arma::mat Y = arma::randn(N, new_Sigma.n_cols);
  return arma::repmat(mu, 1, N).t() + Y * arma::chol(new_Sigma);
}

//' Scaled distance between two vectors
//' 
//' Calculates the scaled distance between two vectors, i.e. calculates the norm of matrix*(x-y)
//' If matrix == identity matrix, this is just the Euclidean distance
//'
//' @param x vector
//' @param y vector
//' @param matrix matrix
//'
//' @return the scaled distance between vectors x and y with matrix 
//' 
//' @examples
//' x <- c(0.3, 0.2, 0.5, 1.2)
//' y <- c(-0.5, 0.8, 1.4, 0.9)
//' scaled_distance(x, y, diag(1, 4))
//' # should equal to the Euclidean distance:
//' sqrt(0.8^2 + 0.6^2 + 0.9^2 + 0.3^2)
// [[Rcpp::export]]
double scaled_distance(const arma::vec &x,
                       const arma::vec &y,
                       const arma::mat &matrix) {
  if (x.size()!=y.size()) {
    stop("scaled_distance: x and y must be the same size");
  } else if (matrix.n_rows!=matrix.n_cols) {
    stop("scaled_distance: matrix must be a square matrix");
  } else if (matrix.n_rows!=x.size()) {
    stop("scaled_distance: dimensions of x, y and matrix are not correct");
  }
  const double sum_dist = arma::sum(arma::square(matrix*(x-y)));
  return(std::sqrt(sum_dist));
}

//' Spectral radius of a symmetric matrix
//' 
//' Calculates the spectral radius of a symmetric matrix A (the largest absolute eigenvalue)
//' 
//' @param A matrix
//' 
//' @return The spectral radius (largest absolute eigenvalue) of A
//' 
//' @examples
//' # symmetric matrix
//' # should equal 2.5
//' spectral_radius(matrix(c(2, 0.5, 0.5, 2), nrow = 2, ncol = 2))
//' # non-symmetrix matrix
//' # should equal 10
//' spectral_radius(matrix(c(9, -1, 2, -2, 8, 4, 1, 1, 8), nrow = 3, ncol = 3, byrow = T))
// [[Rcpp::export]]
double spectral_radius(const arma::mat &A) {
  if (A.is_symmetric()) {
    arma::vec abs_eigenvals = arma::abs(arma::eig_sym(A));
    return(abs_eigenvals.max());
  } else {
    arma::vec abs_eigenvals = arma::abs(arma::eig_gen(A));
    return(abs_eigenvals.max());
  }
}

//' Absolute eigenvalues of a matrix
//' 
//' Calculates the absolute eigenvalues of a matrix
//' 
//' @param A matrix matrix
//' 
//' @return The absolute eigenvalues of A
//' 
//' @examples
//' # symmetric matrix
//' # should equal 2.5, 1.5
//' abs_eigenvals(matrix(c(2, 0.5, 0.5, 2), nrow = 2, ncol = 2))
//' # non-symmetrix matrix
//' # should equal 10, 10, 5
//' abs_eigenvals(matrix(c(9, -1, 2, -2, 8, 4, 1, 1, 8), nrow = 3, ncol = 3, byrow = T))
// [[Rcpp::export]]
arma::vec abs_eigenvals(const arma::mat &A) {
  if (A.is_symmetric()) {
    return(arma::abs(arma::eig_sym(A)));
  } else {
    return(arma::abs(arma::eig_gen(A)));
  }
}

//' Maximal distance to hypercube
//' 
//' Calculates the maximal distance from a vector to a hypercube
//' 
//' @param beta_hat vector
//' @param hypercube vertices matrix of hypercube vertices
//' @param transform_to_X transformation matrix to X-space (original space)
//' @param transform_to_X transformation matrix to Z-space (transformed space)
//' 
//' @return The maximal distance from beta_hat to point in hypercube
// [[Rcpp::export]]
double maximal_distance_hypercube_to_cv(const arma::vec &beta_hat,
                                        const arma::mat &hypercube_vertices,
                                        const arma::mat &transform_to_X,
                                        const arma::mat &transform_to_Z) {
  const arma::mat hypercube_X = arma::trans(transform_to_X * arma::trans(hypercube_vertices));
  arma::vec distances(hypercube_vertices.n_rows, arma::fill::zeros);
  for (int i=0; i < hypercube_X.n_rows; ++i) {
    distances.at(i) = scaled_distance(arma::trans(hypercube_X.row(i)), beta_hat, transform_to_Z);
  }
  return(distances.max());
}
