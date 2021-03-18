#include "../inc/phi_BLR.hpp"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
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
      gradient.at(k) += X.at(i,k)*(y_labels.at(i)-(1/(1+exp(-X_beta.at(i)))));
    }
    gradient.at(k) -= (beta.at(k)-prior_means.at(k))/(C*prior_variances.at(k));
  }
  return(gradient);
}

// [[Rcpp::export]]
double div_log_BLR_gradient(const arma::mat &X,
                            const arma::vec &X_beta,
                            const arma::vec &prior_variances,
                            const double &C,
                            const arma::mat &precondition_mat) {
  double divergence = 0;
  for (int k=0; k < X.n_cols; ++k) {
    double diver = 0;
    for (int i=0; i < X.n_rows; ++i) {
      const double exp_X_beta = exp(X_beta.at(i));
      const double ratio = exp_X_beta/((1+exp_X_beta)*(1+exp_X_beta));
      diver -= X.at(i,k)*X.at(i,k)*ratio;
    }
    diver -= 1/(C*prior_variances.at(k));
    diver *= precondition_mat.at(k,k);
    divergence += diver;
  }
  return(divergence);
}

// [[Rcpp::export]]
double ea_phi_BLR_DL_vec(const arma::vec &beta,
                         const arma::vec &y_labels,
                         const arma::mat &X,
                         const arma::vec &prior_means,
                         const arma::vec &prior_variances,
                         const double &C,
                         const arma::mat &precondition_mat,
                         const arma::mat &transform_mat) {
  const arma::vec transformed_beta = transform_mat * beta;
  const arma::vec X_beta = X * transformed_beta;
  const arma::vec gradient = log_BLR_gradient(transformed_beta,
                                              y_labels,
                                              X,
                                              X_beta,
                                              prior_means,
                                              prior_variances,
                                              C);
  const double t1 = as_scalar((arma::trans(gradient)*precondition_mat)*gradient);
  const double t2 = div_log_BLR_gradient(X,
                                         X_beta,
                                         prior_variances,
                                         C,
                                         precondition_mat);
  return(0.5*(t1+t2));
}

// [[Rcpp::export]]
Rcpp::NumericVector ea_phi_BLR_DL_matrix(const arma::mat &beta,
                                         const arma::vec &y_labels,
                                         const arma::mat &X,
                                         const arma::vec &prior_means,
                                         const arma::vec &prior_variances,
                                         const double &C,
                                         const arma::mat &precondition_mat,
                                         const arma::mat &transform_mat) {
  Rcpp::NumericVector phi(beta.n_rows);
  for (int i=0; i < beta.n_rows; ++i) {
    phi[i] = ea_phi_BLR_DL_vec(arma::trans(beta.row(i)),
                               y_labels,
                               X,
                               prior_means,
                               prior_variances,
                               C,
                               precondition_mat,
                               transform_mat);
  }
  return(phi);
}

// // [[Rcpp::export]]
// double ea_phi_BLR_DL_LB(const arma::mat &X,
//                         const arma::vec &prior_variances,
//                         const double &C,
//                         const arma::mat &precondition_mat) {
//   arma::vec design_sum(X.n_cols, arma::fill::zeros);
//   double prior_variances_sum = 0;
//   for (int k=0; k < X.n_cols; ++k) {
//     for (int i=0; i < X.n_rows; ++i) {
//       design_sum.at(k) += X.at(i,k)*X.at(i,k);
//     }
//     design_sum.at(k) *= precondition_mat.at(k,k);
//     prior_variances_sum += precondition_mat.at(k,k)/(C*prior_variances.at(k));
//   }
//   double LB = -(arma::sum(design_sum)/4) - prior_variances_sum;
//   return(0.5*LB);
// }

// [[Rcpp::export]]
arma::vec log_BLR_gradient_Z(const arma::vec &beta,
                             const arma::vec &y_labels,
                             const arma::vec &X_beta,
                             const arma::mat &transformed_X,
                             const arma::vec &prior_means,
                             const arma::vec &prior_variances,
                             const double &C,
                             const arma::mat &precondition_mat,
                             const arma::mat &transform_mat) {
  arma::vec gradient(beta.size(), arma::fill::zeros);
  for (int k=0; k < transformed_X.n_cols; ++k) {
    for (int i=0; i < transformed_X.n_rows; ++i) {
      gradient.at(k) += transformed_X.at(i,k)*(y_labels.at(i)-(1/(1+exp(-X_beta.at(i)))));
    }
    for (int j=0; j < transformed_X.n_cols; ++j) {
      gradient.at(k) -= (transform_mat.at(j,k)*(beta.at(j)-prior_means.at(j)))/(C*prior_variances.at(j));
    }
  }
  return(gradient);
}

// [[Rcpp::export]]
double div_log_BLR_gradient_Z(const arma::vec &X_beta,
                              const arma::mat &transformed_X,
                              const arma::vec &prior_variances,
                              const double &C,
                              const arma::mat &precondition_mat,
                              const arma::mat &transform_mat) {
  double divergence = 0;
  for (int k=0; k < transformed_X.n_cols; ++k) {
    for (int i=0; i < transformed_X.n_rows; ++i) {
      const double exp_X_beta = exp(X_beta.at(i));
      const double ratio = exp_X_beta/((1+exp_X_beta)*(1+exp_X_beta));
      divergence -= transformed_X.at(i,k)*transformed_X.at(i,k)*ratio;
    }
    for (int j=0; j < transformed_X.n_cols; ++j) {
      divergence -= (transform_mat.at(j,k)*transform_mat.at(j,k))/(C*prior_variances.at(k));
    }
  }
  return(divergence);
  // double divergence = 0;
  // for (int k=0; k < transformed_X.n_cols; ++k) {
  //   for (int i=0; i < transformed_X.n_rows; ++i) {
  //     const double exp_X_beta = exp(X_beta.at(i));
  //     const double ratio = exp_X_beta/((1+exp_X_beta)*(1+exp_X_beta));
  //     divergence  -= transformed_X.at(i,k)*transformed_X.at(i,k)*ratio;
  //   }
  //   for (int j=0; j < transformed_X.n_cols; ++j) {
  //     divergence -= (transform_mat.at(j,k)*transform_mat.at(j,k))/(C*prior_variances.at(j));
  //   }
  // }
  // return(divergence);
}

// [[Rcpp::export]]
double ea_phi_BLR_DL_vec_Z(const arma::vec &beta,
                           const arma::vec &y_labels,
                           const arma::mat &X,
                           const arma::vec &prior_means,
                           const arma::vec &prior_variances,
                           const double &C,
                           const arma::mat &precondition_mat,
                           const arma::mat &transform_mat) {
  const arma::vec X_beta = X * beta;
  const arma::mat transformed_X = X * transform_mat;
  const arma::vec gradient = log_BLR_gradient_Z(beta,
                                                y_labels,
                                                X_beta,
                                                transformed_X,
                                                prior_means,
                                                prior_variances,
                                                C,
                                                precondition_mat,
                                                transform_mat);
  const double t1 = arma::dot(gradient, gradient);
  const double t2 = div_log_BLR_gradient_Z(X_beta,
                                           transformed_X,
                                           prior_variances,
                                           C,
                                           precondition_mat,
                                           transform_mat);
  return(0.5*(t1+t2));
}