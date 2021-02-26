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
  double w_mean = 0;
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
//' @param weighted_mean weighted mean of sampled sub-posterior values
//' @param time time T for fusion algorithm
//' @param precondition_values precondition values associated to each sub-posterior
//'
//' @return the logarithm of rho
//' 
//' @examples
//' x <- rnorm(4, 0, 1)
//' precondition_vals <- c(1, 2, 3, 4)
//' weighted_mean <- weighted_mean_univariate(x = x,
//'                                           weights = 1/precondition_vals)
//' log_rho_univariate(x = x,
//'                    weighted_mean = weighted_mean,
//'                    time = 0.5,
//'                    precondition_values = precondition_vals)
// [[Rcpp::export]] 
double log_rho_univariate(const Rcpp::NumericVector &x,
                          const double &weighted_mean,
                          const double &time,
                          const Rcpp::NumericVector &precondition_values) {
  if (x.size() != precondition_values.size()) {
    stop("log_rho_univariate: x and precondition_values must be of the same size");
  }
  double log_rho = 0;
  for (int i=0; i < x.size(); ++i) {
    double D = (weighted_mean - x[i]);
    double num = D*D/precondition_values[i];
    log_rho -= num / (2*time);
  }
  return log_rho;
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
//' @param weights list of m matricies of the same dimension (n x n)
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
arma::mat calculate_proposal_cov(const double &time, const Rcpp::List &weights) {
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
arma::mat row_wise_subtraction(const arma::mat &X, const arma::vec &vect) {
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
  const arma::mat x_minus_x_mean = row_wise_subtraction(x, x_mean);
  double numerator = 0.0;
  for (int c=0; c < x_minus_x_mean.n_rows; ++c) {
    const arma::mat &inv_precond = inv_precondition_matrices[c];
    numerator += as_scalar((x_minus_x_mean.row(c)*inv_precond)*trans(x_minus_x_mean.row(c)));
  }
  return (-numerator/(2*time));
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
  return Rcpp::List::create(Named("log_weights", log_weights),
                            Named("normalised_weights", normalised_weights),
                            Named("ESS", 1/Rcpp::sum(Rcpp::pow(normalised_weights, 2))));
}


// [[Rcpp::export]]
Rcpp::List rho_IS_univariate_(const Rcpp::List &particles_to_fuse,
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
      const Rcpp::Environment &particles = particles_to_fuse[c];
      const Rcpp::NumericVector &sub_post_samples = particles["y_samples"];
      particle[c] = sub_post_samples.at(i);
    }
    x_means[i] = weighted_mean_univariate(particle, 1/precondition_values);
    log_rho_weights[i] = log_rho_univariate(particle, x_means[i], time, precondition_values);
    x_samples[i] = particle;
  }
  return Rcpp::List::create(Named("x_samples", x_samples),
                            Named("x_means", x_means),
                            Named("log_weights", log_rho_weights),
                            Named("norm_weights", particle_ESS(log_rho_weights)));
}

// [[Rcpp::export]]
Rcpp::List rho_IS_multivariate_(const Rcpp::List &particles_to_fuse,
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
      const Rcpp::Environment &particles = particles_to_fuse[c];
      const arma::mat &sub_post_samples = particles["y_samples"];
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
                            Named("log_weights", log_rho_weights),
                            Named("norm_weights", particle_ESS(log_rho_weights)));
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
