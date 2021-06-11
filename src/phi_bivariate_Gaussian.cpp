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
  const double var1 = sd_vec.at(0)*sd_vec.at(0);
  const double var2 = sd_vec.at(1)*sd_vec.at(1);
  const double sd1sd2 = sd_vec.at(0)*sd_vec.at(1);
  const double mult_term = -beta/(1-corr*corr);
  const arma::vec transformed_x = transform_mat*x;
  const arma::vec y = transformed_x-mean_vec;
  arma::vec grad_log_fc = {mult_term*((y.at(0)/var1)-(corr*y.at(1)/sd1sd2)),
                           mult_term*((y.at(1)/var2)-(corr*y.at(0)/sd1sd2))};
  const double t1 = as_scalar((arma::trans(grad_log_fc)*precondition_mat)*grad_log_fc);
  arma::mat hessian = {{mult_term/var1, -mult_term*corr/sd1sd2},
                       {-mult_term*corr/sd1sd2, mult_term/var2}};
  const double t2 = arma::trace(precondition_mat*hessian);
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
  for (int i=0; i < x.n_rows; ++i) {
    phi[i] = ea_phi_biGaussian_DL_vec(arma::trans(x.row(i)),
                                      mean_vec,
                                      sd_vec,
                                      corr,
                                      beta,
                                      precondition_mat,
                                      transform_mat);
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
  arma::mat evaluate = {{lower[0], lower[1]},
                        {lower[0], upper[1]},
                        {upper[0], lower[1]},
                        {upper[0], upper[1]}};
  Rcpp::NumericVector values = ea_phi_biGaussian_DL_matrix(evaluate,
                                                           mean_vec,
                                                           sd_vec,
                                                           corr,
                                                           beta,
                                                           precondition_mat,
                                                           transform_to_X);
  const arma::vec mean_in_z_space = transform_to_Z * mean_vec;
  if (mean_in_z_space.at(0) >= lower.at(0) && mean_in_z_space.at(0) <= upper.at(0)) {
    if (mean_in_z_space.at(1) >= lower.at(1) && mean_in_z_space.at(1) <= upper.at(1)) {
      values.push_back(ea_phi_biGaussian_DL_vec(mean_vec, mean_vec, sd_vec, corr, beta, precondition_mat, arma::eye(2, 2)));
    }
  }
  return Rcpp::List::create(Rcpp::Named("LB", find_min(values)), Rcpp::Named("UB", find_max(values)));
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

// [[Rcpp::export]]
double gamma_NB_biGaussian(const arma::vec &times,
                           const double &h,
                           const arma::vec &x0,
                           const arma::vec &y,
                           const double &s,
                           const double &t,
                           const arma::vec &mean_vec,
                           const arma::vec &sd_vec,
                           const double &corr,
                           const double &beta,
                           const arma::mat &precondition_mat,
                           const arma::mat &transform_mat) {
  if (times.size() < 2) {
    stop("gamma_NB_biGaussian: length of times must be at least 2"); 
  }
  double sum_phi_eval = 0;
  for (int i=0; i < times.size(); ++i) {
    const arma::vec eval = x0+times.at(i)*(-x0+y)/(t-s);
    Rcpp::List phi = ea_phi_biGaussian_DL_vec(eval,
                                              mean_vec,
                                              sd_vec,
                                              corr,
                                              beta,
                                              precondition_mat,
                                              transform_mat);
    const double &phi_eval = phi["phi"];
    if (i==0 || i==times.size()-1) {
      sum_phi_eval += phi_eval;
    } else {
      sum_phi_eval += 2*phi_eval;
    }
  }
  return(h*sum_phi_eval/2);
}

// // [[Rcpp::export]]
// Rcpp::List terms_biGaussian_X(const arma::vec &x,
//                               const arma::vec &mean_vec,
//                               const arma::vec &sd_vec,
//                               const double &corr,
//                               const double &beta,
//                               const arma::mat &precondition_mat,
//                               const arma::mat &sqrt_precondition_mat) {
//   const double var1 = sd_vec.at(0)*sd_vec.at(0);
//   const double var2 = sd_vec.at(1)*sd_vec.at(1);
//   const double sd1sd2 = sd_vec.at(0)*sd_vec.at(1);
//   const double mult_term = -beta/(1-corr*corr);
//   const arma::vec y = x - mean_vec;
//   arma::vec grad_log_fc(2, arma::fill::zeros);
//   grad_log_fc.at(0) = mult_term*((y.at(0)/var1)-(corr*y.at(1)/sd1sd2));
//   grad_log_fc.at(1) = mult_term*((y.at(1)/var2)-(corr*y.at(0)/sd1sd2));
//   const double t1 = as_scalar((arma::trans(grad_log_fc)*precondition_mat)*grad_log_fc);
//   // const double t2 = (precondition_mat.at(0,0)*mult_term/var1) + (precondition_mat.at(1,1)*mult_term/var2);
//   const double t2_long = mult_term*((sqrt_precondition_mat.at(0,0)*sqrt_precondition_mat.at(0,0)/var1) -
//                                (2*corr*(sqrt_precondition_mat.at(0,0)*sqrt_precondition_mat.at(1,0))/sd1sd2) +
//                                (sqrt_precondition_mat.at(1,0)*sqrt_precondition_mat.at(1,0)/var2) +
//                                (sqrt_precondition_mat.at(0,1)*sqrt_precondition_mat.at(0,1)/var1) -
//                                (2*corr*(sqrt_precondition_mat.at(0,1)*sqrt_precondition_mat.at(1,1))/sd1sd2) +
//                                (sqrt_precondition_mat.at(1,1)*sqrt_precondition_mat.at(1,1)/var2));
//   arma::mat hessian = {{mult_term/var1, -mult_term*corr/sd1sd2}, {-mult_term*corr/sd1sd2, mult_term/var2}};
//   // hessian.at(0,0) = mult_term/var1;
//   // hessian.at(0,1) = -mult_term*corr/sd1sd2;
//   // hessian.at(1,0) = -mult_term*corr/sd1sd2;
//   // hessian.at(1,1) = mult_term/var2;
//   const double t2 = arma::accu(precondition_mat % hessian);
//   Rcout << "t1: " << t1 << "\n";
//   Rcout << "t2 (long): " << t2_long << "\n";
//   Rcout << "t2: " << t2 << "\n";
//   Rcout << "phi: " << 0.5*(t1+t2) << "\n";
//   return(Rcpp::List::create(Named("t1", t1), Named("t2", t2), Named("phi", 0.5*(t1+t2))));
// }

// // [[Rcpp::export]]
// Rcpp::List terms_biGaussian_Z(const arma::vec &x,
//                               const arma::vec &mean_vec,
//                               const arma::vec &sd_vec,
//                               const double &corr,
//                               const double &beta,
//                               const arma::mat &sqrt_precondition_mat) {
//   const double &sqrt_Lam_11 = sqrt_precondition_mat.at(0,0);
//   const double &sqrt_Lam_12 = sqrt_precondition_mat.at(0,1);
//   const double &sqrt_Lam_21 = sqrt_precondition_mat.at(1,0);
//   const double &sqrt_Lam_22 = sqrt_precondition_mat.at(1,1);
//   const double var1 = sd_vec.at(0)*sd_vec.at(0);
//   const double var2 = sd_vec.at(1)*sd_vec.at(1);
//   const double sd1sd2 = sd_vec.at(0)*sd_vec.at(1);
//   const double mult_term = -beta/(1-corr*corr);
//   const arma::vec y = x-mean_vec;
//   arma::vec grad_log_fc(2, arma::fill::zeros);
//   grad_log_fc.at(0) = mult_term*((sqrt_Lam_11*y.at(0)/var1) - (corr*(sqrt_Lam_11*y.at(1) + sqrt_Lam_21*y.at(0))/sd1sd2) + (sqrt_Lam_21*y.at(1)/var2));
//   grad_log_fc.at(1) = mult_term*((sqrt_Lam_12*y.at(0)/var1) - (corr*(sqrt_Lam_12*y.at(1) + sqrt_Lam_22*y.at(0))/sd1sd2) + (sqrt_Lam_22*y.at(1)/var2));
//   const double t1 = arma::dot(grad_log_fc, grad_log_fc);
//   const double t2 = mult_term*((sqrt_Lam_11*sqrt_Lam_11/var1) - (2*corr*(sqrt_Lam_11*sqrt_Lam_21)/sd1sd2) + (sqrt_Lam_21*sqrt_Lam_21/var2)) 
//     + mult_term*((sqrt_Lam_12*sqrt_Lam_12/var1) - (2*corr*(sqrt_Lam_12*sqrt_Lam_22)/sd1sd2) + (sqrt_Lam_22*sqrt_Lam_22/var2));
//   Rcout << "t1: " << t1 << "\n";
//   Rcout << "t2: " << t2 << "\n";
//   Rcout << "phi: " << 0.5*(t1+t2) << "\n";
//   return(Rcpp::List::create(Named("t1", t1), Named("t2", t2), Named("phi", 0.5*(t1+t2))));
// }
