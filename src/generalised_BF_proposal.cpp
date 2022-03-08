#include "../inc/generalised_BF_proposal.hpp"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
const arma::mat construct_V_cpp(const double &s,
                                const double &t,
                                const double &end_time,
                                const int &C,
                                const int &d,
                                const Rcpp::List &precondition_matrices,
                                const arma::mat &Lambda) {
  const double C1 = (t-s)*(end_time-t)/(end_time-s);
  const double C2 = (t-s)*(t-s)/(end_time-s);
  const arma::mat T2 = C2*Lambda;
  arma::mat V(C*d, C*d, arma::fill::zeros);
  for (int i=0; i < C; ++i) {
    const arma::mat &precond_c = precondition_matrices[i];
    const arma::mat T1 = C1*precond_c;
    const int first_row = d*i;
    const int last_row = d*(i+1)-1;
    for (int j=0; j < C; ++j) {
      const int first_col = d*j;
      const int last_col = d*(j+1)-1;
      if (i==j) {
        V.submat(first_row, first_col, last_row, last_col) = T1+T2;
      } else {
        V.submat(first_row, first_col, last_row, last_col) = T2;
      }
    }
  }
  return(V);
}

// [[Rcpp::export]]
const arma::vec construct_M_cpp(const double &s,
                                const double &t,
                                const double &end_time,
                                const int &C,
                                const int &d,
                                const arma::mat &sub_posterior_samples,
                                const arma::rowvec &sub_posterior_mean) {
  const double C1 = (end_time-t)/(end_time-s);
  const double C2 = (t-s)/(end_time-s);
  const arma::rowvec T2 = C2*sub_posterior_mean;
  arma::rowvec M(C*d, arma::fill::zeros);
  for (int i=0; i < C; ++i) {
    const arma::rowvec &sub_posterior_c = sub_posterior_samples.row(i);
    const arma::rowvec T1 = C1*sub_posterior_c;
    const int first_index = d*i;
    const int last_index = d*(i+1)-1;
    M.subvec(first_index, last_index) = T1+T2;
  }
  return(arma::trans(M));
}
