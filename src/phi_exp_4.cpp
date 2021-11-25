#include "../inc/phi_exp_4.hpp"
#include "../inc/helper_functions.hpp"

using namespace Rcpp;

//' phi-function for exp(-(((x-mean)^4)*beta)/2)
//'
//' phi-function for the Exact Algorithm for exp(-(((x-mean)^4)*beta)/2)
//'
//' @param x real value
//' @param beta beta value
//' @param mean mean value
//' @param precondition precondition value
//'
//' @return real value
//'
//' @examples
//' curve(ea_phi_exp_4_DL(x, 0, 1, 1), from = -4, to = 4)
// [[Rcpp::export]]
Rcpp::NumericVector ea_phi_exp_4_DL(const Rcpp::NumericVector &x,
                                    const double &mean,
                                    const double &beta,
                                    const double &precondition) {
  Rcpp::NumericVector phi(x.size());
  for (int i=0; i < x.size(); ++i) {
    double y = x[i]-mean;
    phi[i] = precondition*((2*beta*beta*std::pow(y, 6))-(3*beta*y*y));
  }
  return phi;
}

//' Obtain bounds for phi function
//'
//' Finds the lower and upper bounds of the phi function between a given interval
//' 
//' @param mean mean value
//' @param beta beta value
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
//' mu <- 0.435
//' beta <- 0.583
//' precondition <- 1.243
//' lower <- 0
//' upper <- 1.2
//' 
//' curve(ea_phi_exp_4_DL(x, mu, beta, precondition), lower, upper, ylab = 'phi')
//' abline(v=c(lower, upper))
//' abline(h=ea_phi_exp_4_DL_bounds(mean = mu,
//'                                 beta = beta,
//'                                 precondition = precondition,
//'                                 lower = lower,
//'                                 upper = upper)
//'        col = 'red', lty = 2)
// [[Rcpp::export]]
Rcpp::List ea_phi_exp_4_DL_bounds(const double &mean,
                                  const double &beta,
                                  const double &precondition,
                                  const double &lower,
                                  const double &upper) {
  Rcpp::NumericVector x = Rcpp::NumericVector::create(lower, upper);
  if (mean > lower & mean < upper) {
    x.push_back(mean);
  }
  double m1 = mean - std::pow(1/(2*beta), 0.25);
  double m2 = mean + std::pow(1/(2*beta), 0.25);
  if (m1 > lower & m1 < upper) {
    x.push_back(m1);
  }
  if (m2 > lower & m2 < upper) {
    x.push_back(m2);
  }
  Rcpp::NumericVector phi = ea_phi_exp_4_DL(x, mean, beta, precondition);
  return Rcpp::List::create(Named("LB", find_min(phi)), Named("UB", find_max(phi)));
}

//' Obtain the global lower bound for phi function
//'
//' Finds the global bound of the phi function between a given interval
//' 
//' @param mean mean value
//' @param beta beta value
//' @param precondition precondition value
//'
//' @return The global lower bound of phi
//'
//' @examples
//' mu <- 0.435
//' beta <- 0.583
//' precondition <- 1.243
//' lower <- 0
//' upper <- 1.6
//' 
//' curve(ea_phi_exp_4_DL(x, mu, beta, precondition), lower, upper, ylab = 'phi')
//' abline(v=c(lower, upper))
//' abline(h=ea_phi_exp_4_DL_LB(mean = mu,
//'                             beta = beta,
//'                             precondition = precondition))
//' abline(h=ea_phi_exp_4_DL_bounds(mean = mu,
//'                                 beta = beta,
//'                                 lower = lower,
//'                                 upper = upper,
//'                                 precondition = precondition), 
//'       col = 'red', lty = 2)
// [[Rcpp::export]]
double ea_phi_exp_4_DL_LB(const double &mean,
                          const double &beta,
                          const double &precondition) {
  Rcpp::NumericVector x = Rcpp::NumericVector::create(mean - std::pow(1/(2*beta), 0.25), 
                                                      mean + std::pow(1/(2*beta), 0.25));
  return find_min(ea_phi_exp_4_DL(x, mean, beta, precondition));
}

// [[Rcpp::export]]
double gamma_NB_exp_4(const Rcpp::NumericVector &times,
                      const double &h,
                      const double &x0,
                      const double &y,
                      const double &s,
                      const double &t,
                      const double &mean,
                      const double &beta,
                      const double &precondition) {
  if (times.size() < 2) {
    stop("gamma_NB_exp_4: length of times must be at least 2"); 
  }
  Rcpp::NumericVector eval = (x0*(t-times)+y*(times-s))/(t-s);
  Rcpp::NumericVector phi = ea_phi_exp_4_DL(eval, mean, beta, precondition);
  double sum_phi_eval = 0;
  for (int i=0; i < phi.size(); ++i) {
    if (i==0 || i==phi.size()-1) {
      sum_phi_eval += phi[i];
    } else {
      sum_phi_eval += 2*phi[i];
    }
  }
  return(h*sum_phi_eval/2);
}
