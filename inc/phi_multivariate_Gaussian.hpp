#ifndef PHI_MULTIVARIATE_GAUSSIAN
#define PHI_MULTIVARIATE_GAUSSIAN

#include <RcppArmadillo.h>

double ea_phi_multiGaussian_DL_vec(const arma::vec &x,
                                   const arma::vec &mu,
                                   const arma::mat &inv_Sigma,
                                   const double &beta,
                                   const arma::mat &precondition_mat);

Rcpp::NumericVector ea_phi_multiGaussian_DL_matrix(const arma::mat &x,
                                                   const arma::vec &mu,
                                                   const arma::mat &inv_Sigma,
                                                   const double &beta,
                                                   const arma::mat &precondition_mat);

Rcpp::List ea_phi_multiGaussian_DL_bounds(const arma::vec &mu,
                                          const arma::mat &inv_Sigma,
                                          const double &beta,
                                          const arma::mat &precondition_mat,
                                          const Rcpp::List &hypercube_vertices,
                                          const bool &mean_in_bes_layer);

double ea_phi_multiGaussian_DL_LB(const arma::vec &mu,
                                  const arma::mat &inv_Sigma,
                                  const double &beta,
                                  const arma::mat &precondition_mat);

double gamma_NB_multiGaussian(const arma::vec &times,
                              const double &h,
                              const arma::vec &x0,
                              const arma::vec &y,
                              const double &s,
                              const double &t,
                              const int &dim,
                              const arma::vec &mu,
                              const arma::mat &inv_Sigma,
                              const double &beta,
                              const arma::mat &precondition_mat);

#endif