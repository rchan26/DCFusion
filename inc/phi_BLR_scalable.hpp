#ifndef PHI_SCALABLE_BLR
#define PHI_SCALABLE_BLR

#include <RcppArmadillo.h>

arma::vec datum_log_BLR_gradient(const arma::vec &beta,
                                 const double &y,
                                 const arma::rowvec &X,
                                 const int &data_size,
                                 const arma::vec &prior_means,
                                 const arma::vec &prior_variances,
                                 const double &C);

double datum_div_log_BLR_gradient(const arma::vec &beta,
                                  const arma::rowvec &X,
                                  const int &data_size,
                                  const arma::vec &prior_variances,
                                  const double &C,
                                  const arma::mat &precondition_mat);

arma::vec alpha_tilde(const int &index,
                      const arma::vec &beta,
                      const arma::vec &beta_hat,
                      const arma::vec &y_labels,
                      const arma::mat &X,
                      const int &data_size,
                      const arma::vec &prior_means,
                      const arma::vec &prior_variances,
                      const double &C);

double div_alpha_tilde(const int &index,
                       const arma::vec &beta,
                       const arma::vec &beta_hat,
                       const arma::mat &X,
                       const int &data_size,
                       const arma::vec &prior_variances,
                       const double &C,
                       const arma::mat &precondition_mat);

double ea_phi_BLR_DL_vec_scalable(const Rcpp::List &cv_list,
                                  const arma::vec &beta,
                                  const arma::vec &y_labels,
                                  const arma::mat &X,
                                  const arma::vec &prior_means,
                                  const arma::vec &prior_variances,
                                  const double &C,
                                  const arma::mat &precondition_mat,
                                  const arma::mat &transform_mat);

Rcpp::NumericVector ea_phi_BLR_DL_matrix_scalable(const Rcpp::List &cv_list,
                                                  const arma::mat &beta,
                                                  const arma::vec &y_labels,
                                                  const arma::mat &X,
                                                  const arma::vec &prior_means,
                                                  const arma::vec &prior_variances,
                                                  const double &C,
                                                  const arma::mat &precondition_mat,
                                                  const arma::mat &transform_mat);

double hessian_bound_BLR(const int &dim,
                         const arma::mat &X,
                         const arma::vec &prior_variances,
                         const double &C,
                         const arma::mat &precondition_mat);

// double sqrt_norm(const arma::vec &vect);
// 
// Rcpp::List ea_phi_BLR_DL_scalable_bounds(const Rcpp::List cv_list,
//                                          const double &dim,
//                                          const arma::mat &X,
//                                          const arma::vec &prior_variances,
//                                          const double &C,
//                                          const arma::mat &precondition_mat,
//                                          const Rcpp::NumericVector &lower,
//                                          const Rcpp::NumericVector &upper);

#endif
