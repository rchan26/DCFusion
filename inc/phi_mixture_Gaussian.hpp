#ifndef PHI_MIXTURE_GAUSSIAN
#define PHI_MIXTURE_GAUSSIAN

#include <RcppArmadillo.h>

Rcpp::NumericVector dnorm_mix_tempered_unnormalised(const Rcpp::NumericVector &x,
                                                    const Rcpp::NumericVector &w,
                                                    const Rcpp::NumericVector &m,
                                                    const Rcpp::NumericVector &s,
                                                    const double &b);

Rcpp::NumericVector dnorm_deriv1(const Rcpp::NumericVector &x,
                                 const double &mean,
                                 const double &sd);

Rcpp::NumericVector dnorm_deriv2(const Rcpp::NumericVector &x,
                                 const double &mean,
                                 const double &sd);

Rcpp::NumericVector dnorm_mix_deriv1(const Rcpp::NumericVector &x,
                                     const Rcpp::NumericVector &weights,
                                     const Rcpp::NumericVector &means,
                                     const Rcpp::NumericVector &sds);

Rcpp::NumericVector dnorm_mix_deriv2(const Rcpp::NumericVector &x,
                                     const Rcpp::NumericVector &weights,
                                     const Rcpp::NumericVector &means,
                                     const Rcpp::NumericVector &sds);

Rcpp::List dnorm_mix_tempered_derivatives(const Rcpp::NumericVector &x,
                                          const Rcpp::NumericVector &weights,
                                          const Rcpp::NumericVector &means,
                                          const Rcpp::NumericVector &sds,
                                          const double &beta);

Rcpp::NumericVector ea_phi_mixG_DL(const Rcpp::NumericVector &x,
                                   const int &n_comp,
                                   const Rcpp::NumericVector &weights,
                                   const Rcpp::NumericVector &means,
                                   const Rcpp::NumericVector &sds,
                                   const double &beta,
                                   const double &precondition);

double ea_phi_mixG_DL_LB(const int &n_comp,
                         const Rcpp::NumericVector &weights,
                         const Rcpp::NumericVector &means,
                         const Rcpp::NumericVector &sds,
                         const double &beta,
                         const double &precondition,
                         const double &bounds_multiplier);

#endif