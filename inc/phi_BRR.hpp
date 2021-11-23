#ifndef PHI_BRR
#define PHI_BRR

#include <RcppArmadillo.h>

arma::vec log_BRR_gradient(const arma::vec &beta,
                           const arma::vec &y_resp,
                           const arma::mat &X,
                           const arma::vec &X_beta,
                           const double &nu,
                           const double &sigma,
                           const arma::vec &prior_means,
                           const arma::vec &prior_variances,
                           const double &C);

arma::mat log_BRR_hessian(const arma::vec &y_resp,
                          const arma::mat &X,
                          const arma::vec &X_beta,
                          const double &nu,
                          const double &sigma,
                          const arma::vec &prior_variances,
                          const double &C);

Rcpp::List ea_phi_BRR_DL_vec(const arma::vec &beta,
                             const arma::vec &y_resp,
                             const arma::mat &X,
                             const double &nu,
                             const double &sigma,
                             const arma::vec &prior_means,
                             const arma::vec &prior_variances,
                             const double &C,
                             const arma::mat &precondition_mat);

Rcpp::List ea_phi_BRR_DL_matrix(const arma::mat &beta,
                                const arma::vec &y_resp,
                                const arma::mat &X,
                                const double &nu,
                                const double &sigma,
                                const arma::vec &prior_means,
                                const arma::vec &prior_variances,
                                const double &C,
                                const arma::mat &precondition_mat);

double spectral_radius_BRR(const arma::vec &beta,
                           const int &dim,
                           const arma::vec &y_resp,
                           const arma::mat &X,
                           const double &nu,
                           const double &sigma,
                           const arma::vec &prior_variances,
                           const double &C,
                           const arma::mat &Lambda);

Rcpp::List obtain_hypercube_centre_BRR(const Rcpp::List &bessel_layers,
                                       const arma::mat &transform_to_X,
                                       const arma::vec &y_resp,
                                       const arma::mat &X,
                                       const double &nu,
                                       const double &sigma,
                                       const arma::vec &prior_means,
                                       const arma::vec &prior_variances,
                                       const double &C);

double gamma_NB_BRR(const arma::vec &times,
                    const double &h,
                    const arma::vec &x0,
                    const arma::vec &y,
                    const double &s,
                    const double &t,
                    const arma::vec &y_resp,
                    const arma::mat &X,
                    const double &nu,
                    const double &sigma,
                    const arma::vec &prior_means,
                    const arma::vec &prior_variances,
                    const double &C,
                    const arma::mat &precondition_mat);

#endif