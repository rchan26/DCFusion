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
//' @param hypercube_vertices matrix which determines the points to evaluate
//'                           phi which give the bounds of phi
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
                                          const arma::mat &V) {
  Rcpp::NumericVector values = ea_phi_multiGaussian_DL_matrix(V,
                                                              mu,
                                                              inv_Sigma,
                                                              beta,
                                                              precondition_mat);;
  return Rcpp::List::create(Rcpp::Named("LB", find_min(values)),
                            Rcpp::Named("UB", find_max(values)));
}

// [[Rcpp::export]]
arma::mat POE_multiGaussian_DL(const Rcpp::List &bessel_layers,
                               const arma::vec &mean,
                               const int &dim) {
  if (mean.size()!=dim) {
    stop("POE_multiGaussian_DL: mean must be a mean of length dim");
  }
  if (dim == 1) {
    double L;
    double U;
    double lower_x;
    double upper_x;
    if (bessel_layers.size()==1) {
      const Rcpp::List &bes_layer = bessel_layers[0];
      L = bes_layer["L"];
      U = bes_layer["U"];
    } else {
      L = bessel_layers["L"];
      U = bessel_layers["U"];
    }
    if ((L<mean[0]) && (U>mean[0])) {
      // mean is in Bessel layer, so lower bound is attained at mean
      // upper bound attained at point furthest away from mean
      lower_x = mean[0];
      if (std::abs(L-mean[0]) > std::abs(U-mean[0])) {
        upper_x = L;
      } else {
        upper_x = U;
      }
    } else {
      // mean is not in Bessel layer, so lower bound is attained at point closest to mean
      // upper bound is attained at point furthest away from the mean
      if (std::abs(L-mean[0]) > std::abs(U-mean[0])) {
        lower_x = U;
        upper_x = L;
      } else {
        lower_x = L;
        upper_x = U;
      }
    }
    arma::mat V(2, dim, arma::fill::zeros);
    V.row(0) = lower_x;
    V.row(1) = upper_x;
    return V;
  } else if (dim > 1) {
    if (bessel_layers.size()!=dim) {
      stop("POE_multiGaussian_DL: if dim > 1, bessel_layers must be a list of length dim");
    }
    arma::vec L(dim, arma::fill::zeros);
    arma::vec U(dim, arma::fill::zeros);
    arma::rowvec lower_x(dim, arma::fill::zeros);
    arma::rowvec upper_x(dim, arma::fill::zeros);
    for (int d=0; d < dim; ++d) {
      const Rcpp::List &bes_layer = bessel_layers[d];
      L[d] = bes_layer["L"];
      U[d] = bes_layer["U"];
      if ((L[d]<mean[d]) && (U[d]>mean[d])) {
        // mean is in Bessel layer, so lower bound is attained at mean
        // upper bound attained at point furthest away from mean
        lower_x[d] = mean[d];
        if (std::abs(L[d]-mean[d]) > std::abs(U[d]-mean[d])) {
          upper_x[d] = L[d];
        } else {
          upper_x[d] = U[d];
        }
      } else {
        // mean is not in Bessel layer, so lower bound is attained at point closest to mean
        // upper bound is attained at point furthest away from the mean
        if (std::abs(L[d]-mean[d]) > std::abs(U[d]-mean[d])) {
          lower_x[d] = U[d];
          upper_x[d] = L[d];
        } else {
          lower_x[d] = L[d];
          upper_x[d] = U[d];
        }
      }
    }
    arma::mat V(2, dim, arma::fill::zeros);
    V.row(0) = lower_x;
    V.row(1) = upper_x;
    return V;
  } else {
    stop("POE_multiGaussian_DL: dim must be greater than or equal to 1");
  }
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