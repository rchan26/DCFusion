#ifndef PHI_BRR
#define PHI_BRR

#include <RcppArmadillo.h>

arma::vec log_BPR_gradient(const arma::vec &beta,
                           const arma::vec &y_count,
                           const arma::mat &X,
                           const arma::vec &X_beta,
                           const arma::vec &count,
                           const arma::vec &prior_means,
                           const arma::vec &prior_variances,
                           const double &C);

arma::mat log_BPR_hessian(const arma::mat &X,
                          const arma::vec &X_beta,
                          const arma::vec &count,
                          const arma::vec &prior_variances,
                          const double &C);

Rcpp::List ea_phi_BPR_DL_vec(const arma::vec &beta,
                             const arma::vec &y_count,
                             const arma::mat &X,
                             const arma::vec &count,
                             const arma::vec &prior_means,
                             const arma::vec &prior_variances,
                             const double &C,
                             const arma::mat &precondition_mat);

Rcpp::List ea_phi_BPR_DL_matrix(const arma::mat &beta,
                                const arma::vec &y_count,
                                const arma::mat &X,
                                const arma::vec &count,
                                const arma::vec &prior_means,
                                const arma::vec &prior_variances,
                                const double &C,
                                const arma::mat &precondition_mat);

double spectral_radius_BPR(const arma::vec &beta,
                           const arma::mat &X,
                           const arma::vec &count,
                           const arma::vec &prior_variances,
                           const double &C,
                           const arma::mat &Lambda);

Rcpp::List obtain_hypercube_centre_BPR(const Rcpp::List &bessel_layers,
                                       const arma::mat &transform_to_X,
                                       const arma::vec &y_count,
                                       const arma::mat &X,
                                       const arma::vec &count,
                                       const arma::vec &prior_means,
                                       const arma::vec &prior_variances,
                                       const double &C);

Rcpp::List spectral_radius_bound_BPR_Z(const int &dim,
                                       const arma::mat &V,
                                       const arma::mat &transformed_X,
                                       const arma::vec &count,
                                       const arma::vec &prior_variances,
                                       const double &C,
                                       const arma::mat &sqrt_Lambda);

Rcpp::List ea_phi_BPR_DL_bounds(const arma::vec &beta_hat,
                                const arma::vec &grad_log_hat,
                                const int &dim,
                                const arma::mat &transformed_X,
                                const arma::vec &count,
                                const arma::vec &prior_variances,
                                const double &C,
                                const Rcpp::List &transform_mats,
                                const Rcpp::List &hypercube_vertices);

double gamma_NB_BPR(const arma::vec &times,
                    const double &h,
                    const arma::vec &x0,
                    const arma::vec &y,
                    const double &s,
                    const double &t,
                    const arma::vec &y_count,
                    const arma::mat &X,
                    const arma::vec &count,
                    const arma::vec &prior_means,
                    const arma::vec &prior_variances,
                    const double &C,
                    const arma::mat &precondition_mat);

#endif