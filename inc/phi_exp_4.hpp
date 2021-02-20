#ifndef PHI_EXP_4
#define PHI_EXP_4

#include <RcppArmadillo.h>
#include <cmath>

Rcpp::NumericVector ea_phi_exp_4_DL(const Rcpp::NumericVector &x,
                                    const double &mean,
                                    const double &beta,
                                    const double &precondition);

Rcpp::List ea_phi_exp_4_DL_bounds(const double &mean,
                                  const double &beta,
                                  const double &precondition,
                                  const double &lower,
                                  const double &upper);

double ea_phi_exp_4_DL_LB(const double &mean,
                          const double &beta,
                          const double &precondition);

#endif