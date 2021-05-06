#ifndef HELPER_FUNCTIONS
#define HELPER_FUNCTIONS

#include <RcppArmadillo.h>

double find_max(const Rcpp::NumericVector &vect);

double find_min(const Rcpp::NumericVector &vect);

double weighted_mean_univariate(const Rcpp::NumericVector &x,
                                const Rcpp::NumericVector &weights);

double log_rho_univariate(const Rcpp::NumericVector &x,
                          const double &weighted_mean,
                          const double &time,
                          const Rcpp::NumericVector &precondition_values);

arma::mat inverse_sum_matrices(const Rcpp::List &matrices);

arma::vec weighted_mean_multivariate(const arma::mat &matrix,
                                     const Rcpp::List &weights,
                                     const arma::mat &inverse_sum_weights);

arma::mat calculate_proposal_cov(const double &time, const Rcpp::List &weights);

double scaled_distance(const arma::vec &x, const arma::vec &y, const arma::mat &matrix);

double spectral_radius(const arma::mat &A);

arma::vec abs_eigenvals(const arma::mat &A);

arma::mat row_wise_subtraction(const arma::mat &X, const arma::vec &vect);

double log_rho_multivariate(const arma::mat &x,
                            const arma::vec &x_mean,
                            const double &time,
                            const Rcpp::List &inv_precondition_matrices);

double logsumexp(const Rcpp::NumericVector &x);

Rcpp::List particle_ESS(const Rcpp::NumericVector &log_weights);

Rcpp::List rho_IS_univariate_(const Rcpp::List &samples_to_fuse,
                              const int &N,
                              const int &m,
                              const double &time,
                              const Rcpp::NumericVector &precondition_values);

Rcpp::List rho_IS_multivariate_(const Rcpp::List &samples_to_fuse,
                                const int &dim,
                                const int &N,
                                const int &m,
                                const double &time,
                                const Rcpp::List &inv_precondition_matrices,
                                const arma::mat &inverse_sum_inv_precondition_matrices);

arma::mat mvrnormArma(const int &N,
                      const arma::vec &mu,
                      const arma::mat &Sigma);

arma::mat mvrnormArma_tempered(const int &N,
                               const arma::vec &mu,
                               const arma::mat &Sigma,
                               const double &beta);

#endif