#include "../inc/phi_bivariate_Gaussian.hpp"
#include "../inc/helper_functions.hpp"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double ea_phi_biGaussian_DL_vec(const arma::vec &x,
                                const arma::vec &mean_vec,
                                const arma::vec &sd_vec,
                                const double &corr,
                                const double &beta,
                                const arma::mat &precondition_mat,
                                const arma::mat &transform_mat) {
  // transform point
  const arma::vec transformed_x = transform_mat * x;
  const arma::vec y = transformed_x - mean_vec;
  // calculate the gradient for both dimensions
  arma::vec grad_log_fc(2, arma::fill::zeros);
  const double mult_term = -(beta/(1-corr*corr));
  grad_log_fc.at(0) = mult_term*((y.at(0)/(sd_vec.at(0)*sd_vec.at(0)))-(corr*y.at(1)/(sd_vec.at(0)*sd_vec.at(1))));
  grad_log_fc.at(1) = mult_term*((y.at(1)/(sd_vec.at(1)*sd_vec.at(1)))-(corr*y.at(0)/(sd_vec.at(0)*sd_vec.at(1))));
  // calculate the first term in phi
  double t1 = as_scalar((arma::trans(grad_log_fc)*precondition_mat)*grad_log_fc);
  // calculate second term in phi
  double t2 = - (precondition_mat.at(0,0)*beta/(sd_vec.at(0)*sd_vec.at(0)*(1-corr*corr)))
    - (precondition_mat.at(1,1)*beta/(sd_vec.at(1)*sd_vec.at(1)*(1-corr*corr)));
  return(0.5*(t1+t2));
}

// [[Rcpp::export]]
Rcpp::NumericVector ea_phi_biGaussian_DL_matrix(const arma::mat &x,
                                                const arma::vec &mean_vec,
                                                const arma::vec &sd_vec,
                                                const double &corr,
                                                const double &beta,
                                                const arma::mat &precondition_mat,
                                                const arma::mat &transform_mat) {
  Rcpp::NumericVector phi(x.n_rows);
  const double mult_term = -(beta/(1-corr*corr));
  for (int i=0; i < x.n_rows; ++i) {
    const arma::vec transformed_x = transform_mat * arma::trans(x.row(i));
    const arma::vec y = transformed_x - mean_vec;
    arma::vec grad_log_fc(2, arma::fill::zeros);
    grad_log_fc.at(0) = mult_term*((y.at(0)/(sd_vec.at(0)*sd_vec.at(0)))-(corr*y.at(1)/(sd_vec.at(0)*sd_vec.at(1))));
    grad_log_fc.at(1) = mult_term*((y.at(1)/(sd_vec.at(1)*sd_vec.at(1)))-(corr*y.at(0)/(sd_vec.at(0)*sd_vec.at(1))));
    double t1 = as_scalar((arma::trans(grad_log_fc)*precondition_mat)*grad_log_fc);
    double t2 = - (precondition_mat.at(0,0)*beta/(sd_vec.at(0)*sd_vec.at(0)*(1-corr*corr)))
      - (precondition_mat.at(1,1)*beta/(sd_vec.at(1)*sd_vec.at(1)*(1-corr*corr)));
    phi[i] = 0.5*(t1+t2);
  }
  return(phi);
}

//' Obtain bounds for phi function
//'
//' Finds the lower and upper bounds of the phi function between an interval
//'
//' @param mean_vec vector of length 2 for mean
//' @param sd_vec vector of length 2 for standard deviation
//' @param corr correlation value between component 1 and component 2
//' @param beta real value
//' @param lower vector of length 2 for the lower end of interval
//' @param upper vector of length 2 for the upper end of interval
//' @param precondition_mat precondition matrix
//' @param transform_to_Z the transformation matrix to Z-space
//' @param transform_to_X the transformation matrix to X-space
//'
//' @return A list of components
//' \describe{
//'   \item{LB}{lower bound of phi}
//'   \item{UB}{upper bound of phi}
//' }
// [[Rcpp::export]]
Rcpp::List ea_phi_biGaussian_DL_bounds(const arma::vec &mean_vec,
                                       const arma::vec &sd_vec,
                                       const double &corr,
                                       const double &beta,
                                       const arma::mat &precondition_mat,
                                       const arma::mat &transform_to_Z,
                                       const arma::mat &transform_to_X,
                                       const Rcpp::NumericVector &lower,
                                       const Rcpp::NumericVector &upper) {
  if (lower.size() != 2) {
    stop("ea_phi_biGaussian_DL_bounds: lower is not a vector of length 2");
  } else if (upper.size() != 2) {
    stop("ea_phi_biGaussian_DL_bounds: upper is not a vector of length 2");
  }
  // check corner values
  arma::vec evaluate(2);
  Rcpp::NumericVector values(4);
  // first corner
  evaluate = {lower[0], lower[1]};
  values[0] = ea_phi_biGaussian_DL_vec(evaluate, mean_vec, sd_vec, corr, beta, precondition_mat, transform_to_X);
  // second corner
  evaluate = {lower[0], upper[1]};
  values[1] = ea_phi_biGaussian_DL_vec(evaluate, mean_vec, sd_vec, corr, beta,precondition_mat, transform_to_X);
  // third corner
  evaluate = {upper[0], lower[1]};
  values[2] = ea_phi_biGaussian_DL_vec(evaluate, mean_vec, sd_vec, corr, beta, precondition_mat, transform_to_X);
  // fourth corner
  evaluate = {upper[0], upper[1]};
  values[3] = ea_phi_biGaussian_DL_vec(evaluate, mean_vec, sd_vec, corr, beta, precondition_mat, transform_to_X);
  // check if mean value in the bounds
  const arma::vec mean_in_z_space = transform_to_Z * mean_vec;
  if (mean_in_z_space.at(0) >= lower.at(0) && mean_in_z_space.at(0) <= upper.at(0)) {
    if (mean_in_z_space.at(1) >= lower.at(1) && mean_in_z_space.at(1) <= upper.at(1)) {
      evaluate = {mean_vec[0], mean_vec[1]};
      values.push_back(ea_phi_biGaussian_DL_vec(evaluate, mean_vec, sd_vec, corr, beta, precondition_mat, arma::eye(2, 2)));
    }
  }
  return Rcpp::List::create(Rcpp::Named("LB", find_min(values)),
                            Rcpp::Named("UB", find_max(values)));
}

//' Obtain the global lower bound for phi function
//'
//' Finds the global bound of the phi function between a given interval
//'
//' @param x vector of length 2
//' @param mean_vec vector of length 2 for mean
//' @param sd_vec vector of length 2 for standard deviation
//' @param corr correlation value between component 1 and component 2
//' @param beta real value
//' @param precondition_mat precondition matrix
//' @param transform_mat the transformation matrix
//'
//' @return The global lower bound of phi
// [[Rcpp::export]]
double ea_phi_biGaussian_DL_LB(const arma::vec &mean_vec,
                               const arma::vec &sd_vec,
                               const double &corr,
                               const double &beta,
                               const arma::mat &precondition_mat) {
  return ea_phi_biGaussian_DL_vec(mean_vec,
                                  mean_vec,
                                  sd_vec,
                                  corr,
                                  beta,
                                  precondition_mat,
                                  arma::eye(2,2));
}
