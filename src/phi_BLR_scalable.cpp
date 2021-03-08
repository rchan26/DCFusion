#include "../inc/phi_BLR_scalable.hpp"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

arma::vec datum_log_BLR_gradient(const arma::vec &beta,
                                 const double &y,
                                 const arma::rowvec &X,
                                 const double &X_beta,
                                 const int &data_size,
                                 const arma::vec &prior_means,
                                 const arma::vec &prior_variances,
                                 const double &C) {
  arma::vec gradient(beta.size(), arma::fill::zeros);
  for (int k=0; k < beta.size(); ++k) {
    gradient.at(k) += X.at(k)*(y-(1/(1+exp(-X_beta))));
    gradient.at(k) -= (beta.at(k)-prior_means.at(k))/(data_size*C*prior_variances.at(k));
  }
  return(gradient);
}

double datum_div_log_BLR_gradient(const arma::rowvec &X,
                                  const double &X_beta,
                                  const int &data_size,
                                  const arma::vec &prior_variances,
                                  const double &C,
                                  const arma::mat &precondition_mat) {
  double divergence = 0;
  const double exp_X_beta = exp(X_beta);
  for (int k=0; k < X.size(); ++k) {
    double diver = 0;
    diver -= (X.at(k)*X.at(k)*exp_X_beta)/((1+exp_X_beta)*(1+exp_X_beta));
    diver -= 1/(data_size*C*prior_variances.at(k));
    diver *= precondition_mat.at(k,k);
    divergence += diver;
  }
  return(divergence);
}

arma::vec alpha_tilde(const int &index,
                      const arma::vec &beta,
                      const arma::vec &beta_hat,
                      const arma::vec &y_labels,
                      const arma::mat &X,
                      const int &data_size,
                      const arma::vec &prior_means,
                      const arma::vec &prior_variances,
                      const double &C) {
  const double y = y_labels.at(index);
  const rowvec X_index = X.row(index);
  const double X_beta = arma::dot(X_index, beta);
  arma::vec grad_log = datum_log_BLR_gradient(beta,
                                              y,
                                              X_index,
                                              X_beta,
                                              data_size,
                                              prior_means,
                                              prior_variances,
                                              C);
  arma::vec grad_log_hat = datum_log_BLR_gradient(beta_hat,
                                                  y,
                                                  X_index,
                                                  X_beta,
                                                  data_size,
                                                  prior_means,
                                                  prior_variances,
                                                  C);
  return(data_size*(grad_log-grad_log_hat));
}

double div_alpha_tilde(const int &index,
                       const arma::vec &beta,
                       const arma::vec &beta_hat,
                       const arma::mat &X,
                       const int &data_size,
                       const arma::vec &prior_variances,
                       const double &C,
                       const arma::mat &precondition_mat) {
  const rowvec X_index = X.row(index);
  const double X_beta = arma::dot(X_index, beta);
  double div_grad_log = datum_div_log_BLR_gradient(X_index,
                                                   X_beta,
                                                   data_size,
                                                   prior_variances,
                                                   C,
                                                   precondition_mat);
  double div_grad_log_hat = datum_div_log_BLR_gradient(X_index,
                                                       X_beta,
                                                       data_size,
                                                       prior_variances,
                                                       C,
                                                       precondition_mat);
  
  return(data_size*(div_grad_log-div_grad_log_hat));
}

// [[Rcpp::export]]
double ea_phi_BLR_DL_vec_scalable(const Rcpp::List &cv_list,
                                  const arma::vec &beta,
                                  const arma::vec &y_labels,
                                  const arma::mat &X,
                                  const arma::vec &prior_means,
                                  const arma::vec &prior_variances,
                                  const double &C,
                                  const arma::mat &precondition_mat,
                                  const arma::mat &transform_mat) {
  const arma::vec transformed_beta = transform_mat * beta;
  const double beta_hat = cv_list["beta_hat"];
  const int data_size = cv_list["data_size"]
  const double grad_log_beta_hat = cv_list["grad_log_beta_hat"];
  const int I = as<int>(Rcpp::sample(X.n_rows, 1))-1;
  const int J = as<int>(Rcpp::sample(X.n_rows, 1))-1;
  const arma::vec alpha_I = alpha_tilde(I,
                                        transformed_beta,
                                        beta_hat,
                                        y_labels,
                                        X,
                                        data_size,
                                        prior_means,
                                        prior_variances,
                                        C);
  const arma::vec alpha_J = alpha_tilde(J,
                                        transformed_beta,
                                        beta_hat,
                                        y_labels,
                                        X,
                                        data_size,
                                        prior_means,
                                        prior_variances,
                                        C);
  const double t1 = as_scalar(arma::trans(alpha_I)*precondition_mat*(2*grad_log_beta_hat+alpha_J));
  const double t2 = div_alpha_tilde(I,
                                    transformed_beta,
                                    beta_hat,
                                    X,
                                    data_size,
                                    prior_variances,
                                    C,
                                    precondition_mat);
  return(0.5*(t1+t2));
}

// [[Rcpp::export]]
Rcpp::NumericVector ea_phi_BLR_DL_matrix_scalable(const Rcpp::List &cv_list,
                                                  const arma::mat &beta,
                                                  const arma::vec &y_labels,
                                                  const arma::mat &X,
                                                  const arma::vec &prior_means,
                                                  const arma::vec &prior_variances,
                                                  const double &C,
                                                  const arma::mat &precondition_mat,
                                                  const arma::mat &transform_mat) {
  Rcpp::NumericVector phi(beta.n_rows);
  for (int i=0; i < beta.n_rows; ++i) {
    phi[i] = ea_phi_BLR_DL_vec_scalable(cv_list,
                                        arma::trans(beta.row(i)),
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

double hessian_bound(const int &dim,
                     const arma::mat &X,
                     const arma::vec &prior_variances,
                     const double &C,
                     const arma::mat &precondition_mat) {
  // ----- compute hessian matrix
  arma::mat hessian(dim, dim, arma::fill::zeros);
  int data_size = X.n_rows;
  for (int i=0; i < dim; ++i) {
    for(int j=0; j < dim; ++j) {
      if (i==j) {
        double design_mat_max = abs(X.col(i)).max();
        hessian.at(i,i) = -(design_mat_max*design_mat_max/4)-(1/(data_size*C*prior_variances.at(i)));
      } else {
        double design_mat_max_i = abs(X.col(i)).max();
        double design_mat_max_j = abs(X.col(j)).max();
        hessian.at(i,j) = -(design_mat_max_i*design_mat_max_j/4);
      }
    }
  }
  // ----- calculate spectral norm of A = precondition_mat * hessian
  arma::mat A = precondition_mat * hessian;
  // find eigenvalues of (A^{*} * A), where A^{*} is the conjugate transpose of A
  arma::vec eigenvals = arma::eig_sym(arma::trans(A)*A);
  // obtain the largest eigenvalue
  double max_eigenval = (arma::abs(eigenvals)).max();
  // return square root of the largest eigenvalue
  return(sqrt(max_eigenval));
}

// double sqrt_norm(const arma::vec &vect) {
//   return(sqrt(arma::sum(arma::pow(vect, 2))));
// }
// 
// // [[Rcpp::export]]
// Rcpp::List ea_phi_BLR_DL_scalable_bounds(const Rcpp::List &cv_list,
//                                          const double &dim,
//                                          const arma::mat &X,
//                                          const arma::vec &prior_variances,
//                                          const double &C,
//                                          const arma::mat &precondition_mat,
//                                          const Rcpp::NumericVector &lower,
//                                          const Rcpp::NumericVector &upper) {
//   // ----- calcuating values needed to compute the bounds given by Pollock et al. 2020 eq. 20 and 21
//   // find the maximal distance from a point in the hypercube to the centre point beta_hat
//   const double beta_hat = cv_list["beta_hat"];
//   double dist = maximum_distance_to_beta_hat(dim,
//                                              beta_hat,
//                                              lower,
//                                              upper);
//   // calculate spectral norm of precondition_mat * Hessian matrix for the log per datum posterior
//   double hes_bds = hessian_bound(dim,
//                                  X,
//                                  prior_variances,
//                                  C,
//                                  precondition_mat);
//   // ----- calcuate bounds
//   const int dsz = cv_list["data_size"]
//   const double grad_log_beta_hat = cv_list["grad_log_beta_hat"];
//   double grad_log_beta_hat_bds = sqrt_norm(2*grad_log_beta_hat);
//   double precond_diag_sum = arma::sum(precondition_mat.diag());
//   Rcout << "grad_log_beta_hat_bds: " << grad_log_beta_hat_bds << "\n";
//   Rcout << "dsz: " << dsz << "\n";
//   Rcout << "dist: " << dist << "\n";
//   Rcout << "hes_bds: " << hes_bds << "\n";
//   Rcout << "precond_diag_sum: " << precond_diag_sum << "\n";
//   // calculate bounds given by Pollock et al. 2020 eq. 20 and 21
//   double LB = -dsz*hes_bds*(dist*grad_log_beta_hat_bds+precond_diag_sum);
//   double UB = dsz*hes_bds*dist*(grad_log_beta_hat_bds+dsz*hes_bds*dist) + dsz*precond_diag_sum*hes_bds;
//   return Rcpp::List::create(Rcpp::Named("LB") = 0.5*LB, Rcpp::Named("UB") = 0.5*UB);
// }
// 
// // NOTE: a problem is that my bounds [lower, upper] which are the Bessel layers are in Z space
// // beta_hat is in X-space so the distance could be very large
// // I need to get a function that 
