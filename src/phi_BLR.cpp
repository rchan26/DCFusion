#include "../inc/phi_BLR.hpp"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

arma::vec log_BLR_gradient(const arma::vec &beta,
                           const arma::vec &y_labels,
                           const arma::mat &X,
                           const arma::vec &X_beta,
                           const arma::vec &prior_means,
                           const arma::vec &prior_variances,
                           const double &C) {
  arma::vec gradient(beta.size(), arma::fill::zeros);
  for (int k=0; k < X.n_cols; ++k) {
    for (int i=0; i < X.n_rows; ++i) {
      gradient.at(k) += X.at(i, k)*(y_labels.at(i)-(1/(1+exp(-X_beta.at(i)))));
    }
    gradient.at(k) += -((beta.at(k)-prior_means.at(k)) / (C*prior_variances.at(k)));
  }
  return(gradient);
}

double div_log_BLR_gradient(const arma::mat &X,
                            const arma::vec &X_beta,
                            const arma::vec &prior_variances,
                            const double &C,
                            const arma::mat &precondition_mat) {
  double divergence = 0;
  for (int k=0; k < X.n_cols; ++k) {
    double diver = 0;
    for (int i=0; i < X.n_rows; ++i) {
      const double &exp_X_beta = exp(X_beta.at(i));
      diver -= ((X.at(i, k)*X.at(i, k)*exp_X_beta) / ((1+exp_X_beta)*(1+exp_X_beta)));
    }
    diver -= 1/(C*prior_variances.at(k));
    diver *= precondition_mat.at(k,k);
    divergence += diver;
  }
  return(divergence);
}

// [[Rcpp::export]]
double ea_phi_BLR_DL(const arma::vec &beta,
                     const arma::vec &y_labels,
                     const arma::mat &X,
                     const arma::vec &prior_means,
                     const arma::vec &prior_variances,
                     const double &C,
                     const arma::mat &precondition_mat,
                     const arma::mat &transform_mat) {
  // transform beta
  const arma::vec transformed_beta = transform_mat * beta;
  // computing the matrix multiplication of X*transformed_beta, which is a (n x 1) vector
  arma::vec X_beta = X * transformed_beta;
  // obtain the gradient of log(fc), where fc is the posterior distribution
  // for Bayesian Logistic regression with Gaussian priors
  arma::vec gradient = log_BLR_gradient(transformed_beta,
                                        y_labels,
                                        X,
                                        X_beta,
                                        prior_means,
                                        prior_variances,
                                        C);
  // define a variable for the norm of the gradient squared
  // (the first term in phi)
  double first_term = as_scalar((arma::trans(gradient)*precondition_mat)*gradient);
  // obtain the divergence of the gradient of log(fc)
  // (the second term in phi)
  double divergence = div_log_BLR_gradient(X,
                                           X_beta,
                                           prior_variances,
                                           C,
                                           precondition_mat);
  return(0.5*(first_term+divergence));
}

// [[Rcpp::export]]
double ea_phi_BLR_DL_LB(const arma::mat &X,
                        const arma::vec &prior_variances,
                        const double &C,
                        const arma::mat &precondition_mat) {
  // calculating the sums in the lower bound
  arma::vec design_sum(X.n_cols, arma::fill::zeros);
  double prior_variances_sum = 0;
  // loop through the dimensions to add relevant parts for the lower bound
  for (int k=0; k < X.n_cols; ++k) {
    for (int i=0; i < X.n_rows; ++i) {
      design_sum.at(k) += X.at(i,k)*X.at(i,k);
    }
    design_sum.at(k) *= precondition_mat.at(k,k);
    prior_variances_sum += precondition_mat.at(k,k)/(C*prior_variances.at(k));
  }
  // calculating lower bound for phi
  double LB = -(arma::sum(design_sum)/4) - prior_variances_sum;
  return(0.5*LB);
}