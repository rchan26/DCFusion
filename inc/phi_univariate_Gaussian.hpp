#ifndef PHI_UNIVARIATE_GAUSSIAN
#define PHI_UNIVARIATE_GAUSSIAN

#include <RcppArmadillo.h>
#include <cmath>

Rcpp::NumericVector ea_phi_uniGaussian_DL(const Rcpp::NumericVector &x,
                                          const double &mean,
                                          const double &sd,
                                          const double &beta,
                                          const double &precondition);

Rcpp::List ea_phi_uniGaussian_DL_bounds(const double &mean,
                                        const double &sd,
                                        const double &beta,
                                        const double &precondition,
                                        const double &lower,
                                        const double &upper);

double ea_phi_uniGaussian_DL_LB(const double &mean,
                                const double &sd,
                                const double &beta,
                                const double &precondition);

#endif