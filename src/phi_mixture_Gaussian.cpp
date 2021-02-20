#include "../inc/phi_mixture_Gaussian.hpp"
#include "../inc/helper_functions.hpp"

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericVector dnorm_mix_tempered_unnormalised(const Rcpp::NumericVector &x,
                                                    const Rcpp::NumericVector &w,
                                                    const Rcpp::NumericVector &m,
                                                    const Rcpp::NumericVector &s,
                                                    const double &b) {
  Rcpp::NumericVector value(x.size());
  for (int i=0; i < w.size(); ++i) {
    value += w[i]*Rcpp::dnorm(x, m[i], s[i]);
  }
  return Rcpp::pow(value, b);
}

Rcpp::NumericVector dnorm_deriv1(const Rcpp::NumericVector &x,
                                 const double &mean,
                                 const double &sd) {
  return Rcpp::dnorm(x, mean, sd)*(mean-x)/(sd*sd);
}

Rcpp::NumericVector dnorm_deriv2(const Rcpp::NumericVector &x,
                                 const double &mean,
                                 const double &sd) {;
  return Rcpp::dnorm(x, mean, sd)*((x-mean)*(x-mean)-(sd*sd))/(sd*sd*sd*sd);
}

Rcpp::NumericVector dnorm_mix_deriv1(const Rcpp::NumericVector &x,
                                     const Rcpp::NumericVector &weights,
                                     const Rcpp::NumericVector &means,
                                     const Rcpp::NumericVector &sds) {
  Rcpp::NumericVector value(x.size());
  for (int i=0; i < weights.size(); ++i) {
    value += weights[i]*dnorm_deriv1(x, means[i], sds[i]);
  }
  return value;
}

Rcpp::NumericVector dnorm_mix_deriv2(const Rcpp::NumericVector &x,
                                     const Rcpp::NumericVector &weights,
                                     const Rcpp::NumericVector &means,
                                     const Rcpp::NumericVector &sds) {
  Rcpp::NumericVector value(x.size());
  for (int i=0; i < weights.size(); ++i) {
    value += weights[i]*dnorm_deriv2(x, means[i], sds[i]);
  }
  return value;
}

Rcpp::List dnorm_mix_tempered_derivatives(const Rcpp::NumericVector &x,
                                          const Rcpp::NumericVector &weights,
                                          const Rcpp::NumericVector &means,
                                          const Rcpp::NumericVector &sds,
                                          const double &beta) {
  const Rcpp::NumericVector q = dnorm_mix_tempered_unnormalised(x, weights, means, sds, 1);
  const Rcpp::NumericVector q_prime = dnorm_mix_deriv1(x, weights, means, sds);
  const Rcpp::NumericVector q_dprime = dnorm_mix_deriv2(x, weights, means, sds);
  const Rcpp::NumericVector first_derivative = (beta*q_prime*Rcpp::pow(q, beta-1));
  const Rcpp::NumericVector t1 = beta*(beta-1)*(q_prime*q_prime)*Rcpp::pow(q, beta-2);
  const Rcpp::NumericVector t2 = beta*q_dprime*Rcpp::pow(q, beta-1);
  const Rcpp::NumericVector second_derivative = t1+t2;
  return Rcpp::List::create(Named("first_derivative", first_derivative),
                            Named("second_derivative", second_derivative));
}

//' phi-function for a tempered mixture Gaussian
//'
//' phi-function for the Exact Algorithm for a tempered mixture Gaussian
//'
//' @param x real value
//' @param n_comp integer number of components of mixture Gaussian
//' @param weights vector: weights of mixture Gaussian
//' @param means vector: means of mixture Gassuan
//' @param sds vector: st.devs of mixture Gaussian
//' @param beta real value
//' @param precondition precondition value
//'
//' @return real value
//'
//' @examples
//' weights <- c(0.4, 0.6)
//' means <- c(-3, 6)
//' sds <- c(1, 2)
//' beta <- 1/4
//' precondition <- 1
//' curve(dnorm_mix_tempered(x,
//'                          n_comp = 2,
//'                          weights = weights,
//'                          means = means,
//'                          sds = sds,
//'                          beta = beta),
//'       -20, 40,  n = 100000, ylim = c(-0.5, 2), ylab = 'y')
//' curve(ea_phi_mixG_DL(x,
//'                      n_comp = 2,
//'                      weights = weights,
//'                      means = means,
//'                      sds = sds,
//'                      beta = beta,
//'                      precondition = precondition),
//'       add = T, n = 100000, col = 'red')
// [[Rcpp::export]]
Rcpp::NumericVector ea_phi_mixG_DL(const Rcpp::NumericVector &x,
                                   const int &n_comp,
                                   const Rcpp::NumericVector &weights,
                                   const Rcpp::NumericVector &means,
                                   const Rcpp::NumericVector &sds,
                                   const double &beta,
                                   const double &precondition) {
  if (weights.size()!=n_comp) {
    stop("ea_phi_mixG_DL: weights must be a vector of length n_comp");
  } else if (means.size()!=n_comp) {
    stop("ea_phi_mixG_DL: means must be a vector of length n_comp");
  } else if (sds.size()!=n_comp) {
    stop("ea_phi_mixG_DL: sds must be a vector of length n_comp");
  }
  Rcpp::NumericVector fc = dnorm_mix_tempered_unnormalised(x, weights, means, sds, beta);
  Rcpp::List fc_deriv = dnorm_mix_tempered_derivatives(x, weights, means, sds, beta);
  Rcpp::NumericVector fc_prime = fc_deriv["first_derivative"];
  Rcpp::NumericVector fc_dprime = fc_deriv["second_derivative"];
  Rcpp::NumericVector dlogfc_squared = Rcpp::pow(fc_prime/fc, 2);
  Rcpp::NumericVector d2logfc = ((fc*fc_dprime)-(fc_prime*fc_prime))/(fc*fc);
  return 0.5*precondition*(dlogfc_squared+d2logfc);
}

//' Obtain the global lower bound for phi function
//'
//' Finds the global bound of the phi function between a given interval
//'
//' @param n_comp integer number of components of mixture Gaussian
//' @param weights vector: weights of mixture Gaussian
//' @param means vector: means of mixture Gassuan
//' @param sds vector: st.devs of mixture Gaussian
//' @param beta real value
//' @param normalised boolean value to determine if normalisation constant is calculated
//'
//' @return The global lower bound of phi
//'
//' @examples
//' weights <- c(0.4, 0.6)
//' means <- c(-3, 6)
//' sds <- c(1, 2)
//' beta <- 1
//' precondition <- 1
//' curve(ea_phi_mixG_DL(x,
//'                      n_comp = 2,
//'                      weights = weights,
//'                      means = means,
//'                      sds = sds,
//'                      beta = beta,
//'                      precondition = precondition),
//'       from = -20, to = 25, n = 10000,
//'       ylab = 'phi')
//' PHI <- ea_phi_mixG_DL_LB(n_comp = 2,
//'                          weights = weights,
//'                          means = means,
//'                          sds = sds,
//'                          beta = beta,
//'                          precondition = precondition,
//'                          bounds_multiplier = 1)
//' abline(h = PHI, col = 'red')
// [[Rcpp::export]]
double ea_phi_mixG_DL_LB(const int &n_comp,
                         const Rcpp::NumericVector &weights,
                         const Rcpp::NumericVector &means,
                         const Rcpp::NumericVector &sds,
                         const double &beta,
                         const double &precondition,
                         const double &bounds_multiplier) {
  if (bounds_multiplier < 1) {
    stop("ea_phi_mixG_DL_LB: bounds_multipler should be greater than or equal to 1");
  }
  return bounds_multiplier*find_min(ea_phi_mixG_DL(means, n_comp, weights, means, sds, beta, precondition));
}
