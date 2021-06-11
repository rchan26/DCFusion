#ifndef PHI_BIVARIATE_GAUSSIAN
#define PHI_BIVARIATE_GAUSSIAN

#include <RcppArmadillo.h>

double ea_phi_biGaussian_DL_vec(const arma::vec &x,
                                const arma::vec &mean_vec,
                                const arma::vec &sd_vec,
                                const double &corr,
                                const double &beta,
                                const arma::mat &precondition_mat,
                                const arma::mat &transform_mat);

Rcpp::NumericVector ea_phi_biGaussian_DL_matrix(const arma::mat &x,
                                                const arma::vec &mean_vec,
                                                const arma::vec &sd_vec,
                                                const double &corr,
                                                const double &beta,
                                                const arma::mat &precondition_mat,
                                                const arma::mat &transform_mat);

Rcpp::List ea_phi_biGaussian_DL_bounds(const arma::vec &mean_vec,
                                       const arma::vec &sd_vec,
                                       const double &corr,
                                       const double &beta,
                                       const arma::mat &precondition_mat,
                                       const arma::mat &transform_to_Z,
                                       const arma::mat &transform_to_X,
                                       const Rcpp::NumericVector &lower,
                                       const Rcpp::NumericVector &upper);

double ea_phi_biGaussian_DL_LB(const arma::vec &mean_vec,
                               const arma::vec &sd_vec,
                               const double &corr,
                               const double &beta,
                               const arma::mat &precondition_mat);

double gamma_NB_biGaussian(const arma::vec &times,
                           const double &h,
                           const arma::vec &x0,
                           const arma::vec &y,
                           const double &s,
                           const double &t,
                           const arma::vec &mean_vec,
                           const arma::vec &sd_vec,
                           const double &corr,
                           const double &beta,
                           const arma::mat &precondition_mat,
                           const arma::mat &transform_mat);

#endif