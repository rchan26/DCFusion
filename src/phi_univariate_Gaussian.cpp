#include "../inc/phi_univariate_Gaussian.hpp"
#include "../inc/helper_functions.hpp"

using namespace Rcpp;

//' phi-function for tempered Gaussian distribution
//'
//' phi-function for the Exact Algorithm for tempered Gaussian distribution
//'
//' @param x real value
//' @param mean mean value
//' @param sd standard deviation value
//' @param beta real value
//' @param precondition precondition value
//'
//' @return real value
//'
//' @examples
//' curve(ea_phi_uniGaussian_DL(x, 0, 1, 1, 1), from = -4, to = 4)
// [[Rcpp::export]]
Rcpp::NumericVector ea_phi_uniGaussian_DL(const Rcpp::NumericVector &x,
                                          const double &mean,
                                          const double &sd,
                                          const double &beta,
                                          const double &precondition) {
  Rcpp::NumericVector value(x.size());
  const double term1 = (beta*beta)/(sd*sd*sd*sd);
  const double t2 = beta/(sd*sd);
  for (int i=0; i < x.size(); ++i) {
    double y = x.at(i) - mean;
    double t1 = term1*y*y;
    value.at(i) = 0.5*precondition*(t1+t2);
  }
  return(value);
}

//' Obtain bounds for phi function
//'
//' Finds the lower and upper bounds of the phi function between an interval
//'
//' @param mean mean value
//' @param sd standard deviation value
//' @param beta real value
//' @param lower lower end of interval
//' @param upper upper end of interval
//' @param precondition precondition value
//'
//' @return A list of components
//' \describe{
//'   \item{LB}{lower bound of phi}
//'   \item{UB}{upper bound of phi}
//' }
//'
//' @examples
//' mu <- 0.423
//' sd <- 3.231
//' beta <- 0.8693
//' precondition <- 1.243
//' lower <- -2.823
//' upper <- 4.322
//' curve(ea_phi_uniGaussian_DL(x, mu, sd, beta, precondition), lower, upper)
//' abline(h=ea_phi_uniGaussian_DL_bounds(mean = mu,
//'                                       sd = sd,
//'                                       beta = beta,
//'                                       lower = lower,
//'                                       upper = upper,
//'                                       precondition = precondition),
//'        col = 'red', lty = 2)
//' 
//' # another example where the mean is not in the interval
//' mu <- 0.423
//' sd <- 3.231
//' beta <- 0.8693
//' precondition <- 1.243
//' lower <- 2.823
//' upper <- 5.322
//' curve(ea_phi_uniGaussian_DL(x, mu, sd, beta, precondition), lower, upper)
//' abline(h=ea_phi_uniGaussian_DL_bounds(mean = mu,
//'                                       sd = sd,
//'                                       beta = beta,
//'                                       precondition = precondition,
//'                                       lower = lower,
//'                                       upper = upper),
//'        col = 'red', lty = 2)
// [[Rcpp::export]]
Rcpp::List ea_phi_uniGaussian_DL_bounds(const double &mean,
                                        const double &sd,
                                        const double &beta,
                                        const double &precondition,
                                        const double &lower,
                                        const double &upper) {
  Rcpp::NumericVector x = Rcpp::NumericVector::create(lower, upper);
  if (mean > lower & mean < upper) {
    x.push_back(mean);
  } 
  Rcpp::NumericVector phi = ea_phi_uniGaussian_DL(x, mean, sd, beta, precondition);
  return Rcpp::List::create(Named("LB", find_min(phi)), Named("UB", find_max(phi)));
}

//' Obtain the global lower bound for phi function
//'
//' Finds the global bound of the phi function between a given interval
//'
//' @param mean mean value
//' @param sd standard deviation value
//' @param beta real value
//' @param precondition precondition value
//'
//' @return The global lower bound of phi
//'
//' @examples
//' mu <- 0.423
//' sd <- 3.231
//' beta <- 0.8693
//' precondition <- 1.243
//' lower <- -2.823
//' upper <- 4.322
//' curve(ea_phi_uniGaussian_DL(x, mu, sd, beta, precondition), lower, upper)
//' abline(h=ea_phi_uniGaussian_DL_bounds(mean = mu,
//'                                       sd = sd,
//'                                       beta = beta,
//'                                       precondition = precondition,
//'                                       lower = lower,
//'                                       upper = upper),
//'        col = 'red', lty = 2)
//' abline(h=ea_phi_uniGaussian_DL_LB(mean = mu,
//'                                   sd = sd,
//'                                   beta = beta, 
//'                                   precondition = precondition),
//'        col = 'blue', lty = 3)
//' 
//' # another example where the mean is not in the interval
//' mu <- 0.423
//' sd <- 3.231
//' beta <- 0.8693
//' precondition <- 1.243
//' lower <- 2.823
//' upper <- 5.322
//' curve(ea_phi_uniGaussian_DL(x, mu, sd, beta, precondition), lower, upper)
//' abline(h=ea_phi_uniGaussian_DL_bounds(mean = mu,
//'                                       sd = sd,
//'                                       beta = beta,
//'                                       lower = lower,
//'                                       upper = upper,
//'                                       precondition = precondition),
//'        col = 'red', lty = 2)
//' abline(h=ea_phi_uniGaussian_DL_LB(mean = mu,
//'                                   sd = sd,
//'                                   beta = beta, 
//'                                   precondition = precondition),
//'        col = 'blue', lty = 3)
// [[Rcpp::export]]
double ea_phi_uniGaussian_DL_LB(const double &mean,
                                const double &sd,
                                const double &beta,
                                const double &precondition) {
  Rcpp::NumericVector x = Rcpp::NumericVector::create(mean);
  return ea_phi_uniGaussian_DL(x, mean, sd, beta, precondition)[0];
}
