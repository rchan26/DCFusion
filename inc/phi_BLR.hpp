#ifndef PHI_BLR
#define PHI_BLR

#include <RcppArmadillo.h>

arma::vec log_BLR_gradient(const arma::vec &beta,
                           const arma::vec &y_labels,
                           const arma::mat &X, 
                           const arma::vec &X_beta,
                           const arma::vec &prior_means,
                           const arma::vec &prior_variances,
                           const double &C);

double div_log_BLR_gradient(const arma::mat &X, 
                            const arma::vec &X_beta,
                            const arma::vec &prior_variances,
                            const double &C,
                            const arma::mat &precondition_mat);

Rcpp::List ea_phi_BLR_DL_vec(const arma::vec &beta,
                             const arma::vec &y_labels,
                             const arma::mat &X,
                             const arma::vec &prior_means,
                             const arma::vec &prior_variances,
                             const double &C,
                             const arma::mat &precondition_mat,
                             const arma::mat &transform_mat);

Rcpp::List ea_phi_BLR_DL_matrix(const arma::mat &beta,
                                const arma::vec &y_labels,
                                const arma::mat &X,
                                const arma::vec &prior_means,
                                const arma::vec &prior_variances,
                                const double &C,
                                const arma::mat &precondition_mat,
                                const arma::mat &transform_mat);

double spectral_radius_BLR(const arma::vec &beta,
                           const int &dim,
                           const arma::mat &X,
                           const arma::vec &prior_variances,
                           const double &C,
                           const arma::mat &Lambda);

Rcpp::List max_multiplication(const arma::mat &matrix,
                              const Rcpp::List &bessel_layers);

Rcpp::List spectral_radius_bound_BLR_Z(const int &dim,
                                       const Rcpp::List &bessel_layers,
                                       const arma::mat &X,
                                       const arma::vec &prior_variances,
                                       const double &C,
                                       const arma::mat &sqrt_Lambda);

Rcpp::List spectral_radius_global_bound_BLR_Z(const int &dim,
                                              const arma::mat &X,
                                              const arma::vec &prior_variances,
                                              const double &C,
                                              const arma::mat &sqrt_Lambda);

Rcpp::List obtain_hypercube_centre(const Rcpp::List &bessel_layers,
                                   const arma::mat &transform_to_X,
                                   const arma::vec &y_labels,
                                   const arma::mat &X,
                                   const arma::vec &prior_means,
                                   const arma::vec &prior_variances,
                                   const double &C);

double maximal_distance_hypercube_to_cv(const arma::vec &beta_hat,
                                        const arma::mat &hypercube_vertices,
                                        const arma::mat &transform_to_X,
                                        const arma::mat &transform_to_Z);

Rcpp::List ea_phi_BLR_DL_bounds(const arma::vec &beta_hat,
                                const arma::vec &grad_log_hat,
                                const int &dim,
                                const arma::mat &X,
                                const arma::vec &prior_variances,
                                const double &C,
                                const Rcpp::List &transform_mats,
                                const Rcpp::List &bessel_layers,
                                const arma::mat &hypercube_vertices);

#endif