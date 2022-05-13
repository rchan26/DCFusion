#include "../inc/phi_multivariate_Gaussian.hpp"
#include "../inc/helper_functions.hpp"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double ea_phi_multiGaussian_DL_vec(const arma::vec &x,
                                   const arma::vec &mu,
                                   const arma::mat &inv_Sigma,
                                   const double &beta,
                                   const arma::mat &precondition_mat) {
  const arma::mat hessian = -beta*inv_Sigma;
  const arma::vec gradient = hessian*(x-mu);
  const double t1 = as_scalar((arma::trans(gradient)*precondition_mat)*gradient);
  const double t2 = arma::trace(precondition_mat*hessian);
  return(0.5*(t1+t2));
}

// [[Rcpp::export]]
Rcpp::NumericVector ea_phi_multiGaussian_DL_matrix(const arma::mat &x,
                                                   const arma::vec &mu,
                                                   const arma::mat &inv_Sigma,
                                                   const double &beta,
                                                   const arma::mat &precondition_mat) {
  Rcpp::NumericVector phi(x.n_rows);
  for (int i=0; i < x.n_rows; ++i) {
    phi[i] = ea_phi_multiGaussian_DL_vec(arma::trans(x.row(i)),
                                         mu,
                                         inv_Sigma,
                                         beta,
                                         precondition_mat);
  }
  return(phi);
}

//' Obtain bounds for phi function
//'
//' Finds the lower and upper bounds of the phi function between an interval
//'
//' @param mu vector of length dim for mean
//' @param inv_Sigma dim x dim inverse covariance matrix
//' @param beta real value
//' @param precondition_mat dim x dim precondition matrix
//' @param hypercube_vertices list with item named "V" to determine the points
//'                           to evaluate phi which give the bounds of phi
//'
//' @return A list of components
//' \describe{
//'   \item{LB}{lower bound of phi}
//'   \item{UB}{upper bound of phi}
//' }
// [[Rcpp::export]]
Rcpp::List ea_phi_multiGaussian_DL_bounds(const arma::vec &mu,
                                          const arma::mat &inv_Sigma,
                                          const double &beta,
                                          const arma::mat &precondition_mat,
                                          const Rcpp::List &hypercube_vertices) {
  const arma::mat &V = hypercube_vertices["V"];
  Rcpp::NumericVector values = ea_phi_multiGaussian_DL_matrix(V,
                                                              mu,
                                                              inv_Sigma,
                                                              beta,
                                                              precondition_mat);
  return Rcpp::List::create(Rcpp::Named("LB", find_min(values)),
                            Rcpp::Named("UB", find_max(values)));
}

//' Obtain the global lower bound for phi function
//'
//' Finds the global bound of the phi function between a given interval
//'
//' @param mu vector of length dim for mean
//' @param inv_Sigma dim x dim inverse covariance matrix
//' @param beta real value
//' @param precondition_mat dim x dim precondition matrix
//'
//' @return The global lower bound of phi
// [[Rcpp::export]]
double ea_phi_multiGaussian_DL_LB(const arma::vec &mu,
                                  const arma::mat &inv_Sigma,
                                  const double &beta,
                                  const arma::mat &precondition_mat) {
  return ea_phi_multiGaussian_DL_vec(mu,
                                     mu,
                                     inv_Sigma,
                                     beta,
                                     precondition_mat);
}

// [[Rcpp::export]]
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
                              const arma::mat &precondition_mat) {
  if (times.size() < 2) {
    stop("gamma_NB_multiGaussian: length of times must be at least 2"); 
  }
  double sum_phi_eval = 0;
  for (int i=0; i < times.size(); ++i) {
    const arma::vec eval = (x0*(t-times.at(i))+y*(times.at(i)-s))/(t-s);
    const double phi_eval = ea_phi_multiGaussian_DL_vec(eval,
                                                        mu,
                                                        inv_Sigma,
                                                        beta,
                                                        precondition_mat);
    if (i==0 || i==times.size()-1) {
      sum_phi_eval += phi_eval;
    } else {
      sum_phi_eval += 2*phi_eval;
    }
  }
  return(h*sum_phi_eval/2);
}