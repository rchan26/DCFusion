#ifndef GENERALISED_BF_PROPOSAL
#define GENERALISED_BF_PROPOSAL

#include <RcppArmadillo.h>

const arma::mat construct_V_cpp(const double &s,
                                const double &t,
                                const double &end_time,
                                const int &C,
                                const int &d,
                                const Rcpp::List &precondition_matrices,
                                const arma::mat &Lambda);

const arma::vec construct_M_cpp(const double &s,
                                const double &t,
                                const double &end_time,
                                const int &C,
                                const int &d,
                                const arma::mat &sub_posterior_samples,
                                const arma::rowvec &sub_posterior_mean);

#endif