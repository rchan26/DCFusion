#include "../inc/phi_BLR_scalable.hpp"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

arma::vec datum_log_BLR_gradient(const arma::vec &beta,
                                 const double &y,
                                 const arma::rowvec &X,
                                 const int &data_size,
                                 const arma::vec &prior_means,
                                 const arma::vec &prior_variances,
                                 const double &C) {
  arma::vec gradient(beta.size(), arma::fill::zeros);
  const double exp_minus_X_beta = exp(-arma::dot(X, beta));
  for (int k=0; k < beta.size(); ++k) {
    gradient.at(k) += X.at(k)*(y-(1/(1+exp_minus_X_beta)));
    gradient.at(k) -= (beta.at(k)-prior_means.at(k))/(data_size*C*prior_variances.at(k));
  }
  return(gradient);
}

arma::mat datum_log_BLR_hessian(const arma::vec &beta,
                                const arma::rowvec &X,
                                const int &data_size,
                                const arma::vec &prior_variances,
                                const double &C) {
  arma::mat hessian(X.size(), X.size(), arma::fill::zeros);
  const double exp_X_beta = exp(arma::dot(X, beta));
  const double ratio = exp_X_beta/((1+exp_X_beta)*(1+exp_X_beta));
  for (int j=0; j < X.size(); ++j) {
    for (int k=0; k < X.size(); ++k) {
      hessian.at(j,k) -= X.at(j)*X.at(k)*ratio;
    }
  }
  for (int k=0; k < X.n_cols; ++k) {
    hessian.at(k,k) -= 1/(data_size*C*prior_variances.at(k));
  }
  return(hessian);
}

// [[Rcpp::export]]
arma::vec alpha_tilde(const int &index,
                      const arma::vec &beta,
                      const arma::vec &beta_hat,
                      const arma::vec &y_labels,
                      const arma::mat &X,
                      const int &data_size,
                      const arma::vec &prior_means,
                      const arma::vec &prior_variances,
                      const double &C) {
  const double &y = y_labels.at(index);
  const arma::rowvec &X_index = X.row(index);
  const arma::vec grad_log = datum_log_BLR_gradient(beta,
                                                    y,
                                                    X_index,
                                                    data_size,
                                                    prior_means,
                                                    prior_variances,
                                                    C);
  const arma::vec grad_log_hat = datum_log_BLR_gradient(beta_hat,
                                                        y,
                                                        X_index,
                                                        data_size,
                                                        prior_means,
                                                        prior_variances,
                                                        C);
  return(data_size*(grad_log-grad_log_hat));
}

// [[Rcpp::export]]
arma::mat log_BLR_hessian_tilde(const int &index,
                                const arma::vec &beta,
                                const arma::vec &beta_hat,
                                const arma::mat &X,
                                const int &data_size,
                                const arma::vec &prior_variances,
                                const double &C) {
  const arma::rowvec &X_index = X.row(index);
  const arma::mat hessian = datum_log_BLR_hessian(beta,
                                                  X_index,
                                                  data_size,
                                                  prior_variances,
                                                  C);
  const arma::mat hessian_hat = datum_log_BLR_hessian(beta_hat,
                                                      X_index,
                                                      data_size,
                                                      prior_variances,
                                                      C);
  return(data_size*(hessian-hessian_hat));
}

// [[Rcpp::export]]
Rcpp::List ea_phi_BLR_DL_vec_scalable(const Rcpp::List &cv_list,
                                      const arma::vec &beta,
                                      const arma::vec &y_labels,
                                      const arma::mat &X,
                                      const arma::vec &prior_means,
                                      const arma::vec &prior_variances,
                                      const double &C,
                                      const arma::mat &precondition_mat) {
  const arma::vec &beta_hat = cv_list["beta_hat"];
  const int &data_size = cv_list["data_size"];
  const arma::vec &grad_log_beta_hat = cv_list["grad_log_beta_hat"];
  const int I = as<int>(Rcpp::sample(X.n_rows, 1))-1;
  const int J = as<int>(Rcpp::sample(X.n_rows, 1))-1;
  const arma::vec alpha_I = alpha_tilde(I,
                                        beta,
                                        beta_hat,
                                        y_labels,
                                        X,
                                        data_size,
                                        prior_means,
                                        prior_variances,
                                        C);
  const arma::vec alpha_J = alpha_tilde(J,
                                        beta,
                                        beta_hat,
                                        y_labels,
                                        X,
                                        data_size,
                                        prior_means,
                                        prior_variances,
                                        C);
  const double t1 = as_scalar(arma::trans(alpha_I)*precondition_mat*(2*grad_log_beta_hat+alpha_J));
  const arma::mat hessian = log_BLR_hessian_tilde(I,
                                                  beta,
                                                  beta_hat,
                                                  X,
                                                  data_size,
                                                  prior_variances,
                                                  C);
  const double t2 = arma::trace(precondition_mat*hessian);
  const double &constant = cv_list["constant"];
  return(Rcpp::List::create(Rcpp::Named("phi", 0.5*(t1+t2)+constant),
                            Rcpp::Named("I", I+1),
                            Rcpp::Named("J", J+1)));
}

// [[Rcpp::export]]
Rcpp::List ea_phi_BLR_DL_matrix_scalable(const Rcpp::List &cv_list,
                                         const arma::mat &beta,
                                         const arma::vec &y_labels,
                                         const arma::mat &X,
                                         const arma::vec &prior_means,
                                         const arma::vec &prior_variances,
                                         const double &C,
                                         const arma::mat &precondition_mat) {
  Rcpp::NumericVector phi(beta.n_rows);
  Rcpp::NumericVector I(beta.n_rows);
  Rcpp::NumericVector J(beta.n_rows);
  for (int i=0; i < beta.n_rows; ++i) {
    Rcpp::List phi_eval = ea_phi_BLR_DL_vec_scalable(cv_list,
                                                     arma::trans(beta.row(i)),
                                                     y_labels,
                                                     X,
                                                     prior_means,
                                                     prior_variances,
                                                     C,
                                                     precondition_mat);
    phi[i] = phi_eval["phi"];
    I[i] = phi_eval["I"];
    J[i] = phi_eval["J"];
  }
  return(Rcpp::List::create(Rcpp::Named("phi", phi),
                            Rcpp::Named("I", I),
                            Rcpp::Named("J", J)));
}

// [[Rcpp::export]]
double spectral_norm_bound_BLR(const int &dim,
                               const arma::mat &X,
                               const arma::vec &prior_variances,
                               const double &C,
                               const arma::mat &sqrt_precondition_mat,
                               const arma::mat &precondition_mat) {
  // ----- compute the Hessian of the transformed log sub-posterior
  const arma::mat transformed_X = X * sqrt_precondition_mat;
  arma::mat hessian(dim, dim, arma::fill::zeros);
  const int data_size = X.n_rows;
  for (int i=0; i < dim; ++i) {
    for (int j=0; j < dim; ++j) {
      const double col_i_max = abs(transformed_X.col(i)).max();
      const double col_j_max = abs(transformed_X.col(j)).max();
      hessian.at(i,j) -= col_i_max*col_j_max/4;
      for (int k=0; k < dim; ++k) {
        hessian.at(i,j) -= sqrt_precondition_mat.at(k,i)*sqrt_precondition_mat.at(k,j)/(data_size*C*prior_variances.at(i));
      }
    }
  }
  // ----- calculate spectral norm of A = precondition_mat * hessian
  arma::mat A = precondition_mat * hessian;
  // find eigenvalues of (A^{*} * A), where A^{*} is the (conjugate) transpose of A
  arma::vec eigenvals = arma::eig_sym(arma::trans(A)*A);
  // obtain the largest eigenvalue
  double max_eigenval = (arma::abs(eigenvals)).max();
  // return square root of the largest eigenvalue
  return(sqrt(max_eigenval));
}

// [[Rcpp::export]]
double spectral_norm_BLR(const int &dim,
                         const arma::vec &beta,
                         const arma::mat &X,
                         const int &index,
                         const arma::vec &prior_variances,
                         const double &C,
                         const arma::mat &precondition_mat) {
  const arma::rowvec &X_index = X.row(index-1);
  arma::mat hessian(dim, dim, arma::fill::zeros);
  const int data_size = X.n_rows;
  const double exp_X_beta = exp(arma::dot(X_index, beta));
  const double ratio = exp_X_beta/((1+exp_X_beta)*(1+exp_X_beta));
  for (int i=0; i < dim; ++i) {
    for (int j=0; j < dim; ++j) {
      if (i==j) {
        hessian.at(i,i) = -(X_index.at(i)*X_index.at(i)*ratio)-(1/(data_size*C*prior_variances.at(i)));
      } else {
        hessian.at(i,j) = -(X_index.at(i)*X_index.at(j)*ratio);
      }
    }
  }
  arma::mat A = precondition_mat * hessian;
  arma::vec eigenvals = arma::eig_sym(arma::trans(A)*A);
  double max_eigenval = (arma::abs(eigenvals)).max();
  return(sqrt(max_eigenval));
}
