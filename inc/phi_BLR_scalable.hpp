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

Rcpp::List ea_phi_BLR_DL_vec_scalable(const Rcpp::List &cv_list,
                                      const arma::vec &beta,
                                      const arma::vec &y_labels,
                                      const arma::mat &X,
                                      const arma::vec &prior_means,
                                      const arma::vec &prior_variances,
                                      const double &C,
                                      const arma::mat &precondition_mat);

Rcpp::List ea_phi_BLR_DL_matrix_scalable(const Rcpp::List &cv_list,
                                         const arma::mat &beta,
                                         const arma::vec &y_labels,
                                         const arma::mat &X,
                                         const arma::vec &prior_means,
                                         const arma::vec &prior_variances,
                                         const double &C,
                                         const arma::mat &precondition_mat);

double spectral_norm_bound_BLR(const int &dim,
                               const arma::mat &X,
                               const arma::vec &prior_variances,
                               const double &C,
                               const arma::mat &sqrt_precondition_mat,
                               const arma::mat &precondition_mat);

#endif
