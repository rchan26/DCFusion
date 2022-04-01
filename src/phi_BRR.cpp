#include "../inc/helper_functions.hpp"
#include "../inc/phi_BRR.hpp"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec log_BRR_gradient(const arma::vec &beta,
                           const arma::vec &y_resp,
                           const arma::mat &X,
                           const arma::vec &X_beta,
                           const double &nu,
                           const double &sigma,
                           const arma::vec &prior_means,
                           const arma::vec &prior_variances,
                           const double &C) {
  arma::vec gradient(beta.size(), arma::fill::zeros);
  const double nu_sigma_sq = nu*sigma*sigma;
  for (int k=0; k < X.n_cols; ++k) {
    double sum = 0;
    for (int i=0; i < X.n_rows; ++i) {
      const double diff = y_resp.at(i)-X_beta.at(i);
      sum += X.at(i,k)*diff / (nu_sigma_sq+(diff*diff));
    }
    gradient.at(k) = ((nu+1)*sum) - ((beta.at(k)-prior_means.at(k))/(C*prior_variances.at(k)));
  }
  return(gradient);
}

// [[Rcpp::export]]
arma::mat log_BRR_hessian(const arma::vec &y_resp,
                          const arma::mat &X,
                          const arma::vec &X_beta,
                          const double &nu,
                          const double &sigma,
                          const arma::vec &prior_variances,
                          const double &C) {
  arma::mat hessian(X.n_cols, X.n_cols, arma::fill::zeros);
  const double nu_sigma_sq = nu*sigma*sigma;
  for (int i=0; i < X.n_rows; ++i) {
    const double diff_sq = (y_resp.at(i)-X_beta.at(i))*(y_resp.at(i)-X_beta.at(i));
    const double ratio = (diff_sq-nu_sigma_sq) /
      ((diff_sq+nu_sigma_sq)*(diff_sq+nu_sigma_sq));
    for (int j=0; j < X.n_cols; ++j) {
      for (int k=0; k <= j; ++k) {
        hessian.at(j,k) += X.at(i,j)*X.at(i,k)*ratio;
      }
    }
  }
  for (int j=0; j < X.n_cols; ++j) {
    hessian.at(j,j) = (nu+1)*hessian.at(j,j) - (1/(C*prior_variances.at(j)));
    for (int k=0; k < j; ++k) {
      hessian.at(j,k) = (nu+1)*hessian.at(j,k);
      hessian.at(k,j) = hessian.at(j,k);
    }
  }
  return(hessian);
}

// [[Rcpp::export]]
Rcpp::List ea_phi_BRR_DL_vec(const arma::vec &beta,
                             const arma::vec &y_resp,
                             const arma::mat &X,
                             const double &nu,
                             const double &sigma,
                             const arma::vec &prior_means,
                             const arma::vec &prior_variances,
                             const double &C,
                             const arma::mat &precondition_mat) {
  const arma::vec X_beta = X * beta;
  const arma::vec gradient = log_BRR_gradient(beta,
                                              y_resp,
                                              X,
                                              X_beta,
                                              nu,
                                              sigma,
                                              prior_means,
                                              prior_variances,
                                              C);
  const double t1 = as_scalar((arma::trans(gradient)*precondition_mat)*gradient);
  const arma::mat hessian = log_BRR_hessian(y_resp,
                                            X,
                                            X_beta,
                                            nu,
                                            sigma,
                                            prior_variances,
                                            C);
  const double t2 = arma::trace(precondition_mat*hessian);
  return(Rcpp::List::create(Named("phi", 0.5*(t1+t2)),
                            Named("t1", t1),
                            Named("t2", t2)));
}

// [[Rcpp::export]]
Rcpp::List ea_phi_BRR_DL_matrix(const arma::mat &beta,
                                const arma::vec &y_resp,
                                const arma::mat &X,
                                const double &nu,
                                const double &sigma,
                                const arma::vec &prior_means,
                                const arma::vec &prior_variances,
                                const double &C,
                                const arma::mat &precondition_mat) {
  Rcpp::NumericVector phi(beta.n_rows);
  Rcpp::NumericVector t1(beta.n_rows);
  Rcpp::NumericVector t2(beta.n_rows);
  for (int i=0; i < beta.n_rows; ++i) {
    Rcpp::List phi_eval = ea_phi_BRR_DL_vec(arma::trans(beta.row(i)),
                                            y_resp,
                                            X,
                                            nu,
                                            sigma,
                                            prior_means,
                                            prior_variances,
                                            C,
                                            precondition_mat);
    phi[i] = phi_eval["phi"];
    t1[i] = phi_eval["t1"];
    t2[i] = phi_eval["t2"];
  }
  return(Rcpp::List::create(Rcpp::Named("phi", phi),
                            Rcpp::Named("t1", t1),
                            Rcpp::Named("t2", t2)));
}

// [[Rcpp::export]]
double spectral_radius_BRR(const arma::vec &beta,
                           const arma::vec &y_resp,
                           const arma::mat &X,
                           const double &nu,
                           const double &sigma,
                           const arma::vec &prior_variances,
                           const double &C,
                           const arma::mat &Lambda) {
  const arma::vec X_beta = X * beta;
  return(spectral_radius(Lambda * log_BRR_hessian(y_resp,
                                                  X,
                                                  X_beta,
                                                  nu,
                                                  sigma,
                                                  prior_variances,
                                                  C)));
}

// [[Rcpp::export]]
Rcpp::List obtain_hypercube_centre_BRR(const Rcpp::List &bessel_layers,
                                       const arma::mat &transform_to_X,
                                       const arma::vec &y_resp,
                                       const arma::mat &X,
                                       const double &nu,
                                       const double &sigma,
                                       const arma::vec &prior_means,
                                       const arma::vec &prior_variances,
                                       const double &C) {
  arma::vec centre(bessel_layers.size(), arma::fill::zeros);
  for (int i=0; i < bessel_layers.size(); ++i) {
    const Rcpp::List &b_layer = bessel_layers[i];
    const double &L = b_layer["L"];
    const double &U = b_layer["U"];
    centre.at(i) = 0.5*(L+U);
  }
  const arma::vec beta_hat = transform_to_X * centre;
  const arma::vec X_beta = X * beta_hat;
  return(Rcpp::List::create(Named("beta_hat", beta_hat),
                            Named("grad_log_hat", log_BRR_gradient(beta_hat,
                                                       y_resp,
                                                       X,
                                                       X_beta,
                                                       nu,
                                                       sigma,
                                                       prior_means,
                                                       prior_variances,
                                                       C))));
}

// // [[Rcpp::export]]
// Rcpp::List spectral_radius_bound_BRR_Z(const int &dim,
//                                        const arma::mat &V,
//                                        const arma::vec &y_resp,
//                                        const arma::mat &transformed_X,
//                                        const double &nu,
//                                        const double &sigma,
//                                        const arma::vec &prior_variances,
//                                        const double &C,
//                                        const arma::mat &sqrt_Lambda) {
//   arma::mat hessian(dim, dim, arma::fill::zeros);
//   arma::vec terms(V.n_rows, arma::fill::zeros);
//   const double nu_sigma_sq = nu*sigma*sigma;
//   for (int i=0; i < transformed_X.n_rows; ++i) {
//     for (int v=0; v < V.n_rows; ++v) {
//       const double product = arma::dot(transformed_X.row(i), V.row(v));
//       const double denominator = (y_resp.at(i)-product)*(y_resp.at(i)-product) + nu_sigma_sq;
//       terms.at(v) = (1/denominator)-((2*nu_sigma_sq)/(denominator*denominator));
//     }
//     const double max_term = terms.max();
//     for (int k=0; k < dim; ++k) {
//       for (int l=0; l <= k; ++l) {
//         hessian.at(k,l) -= transformed_X.at(i,k)*transformed_X.at(i,l)*max_term;
//       }
//     }
//   }
//   for (int k=0; k < dim; ++k) {
//     for (int l=0; l <= k; ++l) {
//       for (int j=0; j < dim; ++j) {
//         hessian.at(k,l) -= sqrt_Lambda.at(j,k)*sqrt_Lambda(j,l)/(C*prior_variances.at(j));
//       }
//       if (l!=k) {
//         hessian.at(l,k) = hessian.at(k,l);
//       }
//     }
//   }
//   const arma::vec abs_eigen = abs_eigenvals(hessian);
//   return(Rcpp::List::create(Named("spectral_radius", abs_eigen.max()),
//                             Named("abs_eigenvals", abs_eigen)));
// }

// [[Rcpp::export]]
Rcpp::List spectral_radius_global_bound_BRR_Z(const int &dim,
                                              const arma::mat &transformed_X,
                                              const double &nu,
                                              const double &sigma,
                                              const arma::vec &prior_variances,
                                              const double &C,
                                              const arma::mat &sqrt_Lambda) {
  arma::mat hessian(dim, dim, arma::fill::zeros);
  const double ratio = nu*sigma*sigma/8;
  for (int i=0; i < transformed_X.n_rows; ++i) {
    for (int k=0; k < dim; ++k) {
      for (int l=0; l <= k; ++l) {
        hessian.at(k,l) += transformed_X.at(i,k)*transformed_X.at(i,l)*ratio;
      }
    }
  }
  for (int k=0; k < dim; ++k) {
    for (int l=0; l <= k; ++l) {
      for (int j=0; j < dim; ++j) {
        hessian.at(k,l) -= sqrt_Lambda.at(j,k)*sqrt_Lambda(j,l)/(C*prior_variances.at(j));
      }
      if (l!=k) {
        hessian.at(l,k) = hessian.at(k,l);
      }
    }
  }
  const arma::vec abs_eigen = abs_eigenvals(hessian);
  return(Rcpp::List::create(Named("spectral_radius", abs_eigen.max()),
                            Named("abs_eigenvals", abs_eigen)));
}

// [[Rcpp::export]]
Rcpp::List ea_phi_BRR_DL_bounds(const arma::vec &beta_hat,
                                const arma::vec &grad_log_hat,
                                const int &dim,
                                const arma::vec &y_resp,
                                const arma::mat &transformed_X,
                                const double &nu,
                                const double &sigma,
                                const arma::vec &prior_variances,
                                const double &C,
                                const Rcpp::List &transform_mats,
                                const Rcpp::List &hypercube_vertices) {
  const arma::mat &transform_to_X = transform_mats["to_X"];
  const arma::mat &transform_to_Z = transform_mats["to_Z"];
  const double vec_norm = std::sqrt(arma::sum(arma::square(transform_to_X*grad_log_hat)));
  const arma::mat &vertices = hypercube_vertices["vertices"];
  const double dist = maximal_distance_hypercube_to_cv(beta_hat,
                                                       vertices,
                                                       transform_to_X,
                                                       transform_to_Z);
  Rcpp::List spectral_radius_bds = spectral_radius_global_bound_BRR_Z(dim,
                                                                      transformed_X,
                                                                      nu,
                                                                      sigma,
                                                                      prior_variances,
                                                                      C,
                                                                      transform_to_X);;
  double P_n_Lambda = spectral_radius_bds["spectral_radius"];
  return(Rcpp::List::create(Named("LB", -0.5*dim*P_n_Lambda),
                            Named("UB", 0.5*((vec_norm+dist*P_n_Lambda)*(vec_norm+dist*P_n_Lambda)+dim*P_n_Lambda)),
                            Named("dist", dist),
                            Named("P_n_Lambda", P_n_Lambda),
                            Named("t1_bds", (vec_norm+dist*P_n_Lambda)*(vec_norm+dist*P_n_Lambda)),
                            Named("t2_bds", dim*P_n_Lambda)));
}

// [[Rcpp::export]]
double gamma_NB_BRR(const arma::vec &times,
                    const double &h,
                    const arma::vec &x0,
                    const arma::vec &y,
                    const double &s,
                    const double &t,
                    const arma::vec &y_resp,
                    const arma::mat &X,
                    const double &nu,
                    const double &sigma,
                    const arma::vec &prior_means,
                    const arma::vec &prior_variances,
                    const double &C,
                    const arma::mat &precondition_mat) {
  if (times.size() < 2) {
    stop("gamma_NB_BRR: length of times must be at least 2");
  }
  double sum_phi_eval = 0;
  for (int i=0; i < times.size(); ++i) {
    const arma::vec eval = (x0*(t-times.at(i))+y*(times.at(i)-s))/(t-s);
    Rcpp::List phi = ea_phi_BRR_DL_vec(eval,
                                       y_resp,
                                       X,
                                       nu,
                                       sigma,
                                       prior_means,
                                       prior_variances,
                                       C,
                                       precondition_mat);
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
// arma::vec log_BRR_gradient_Z(const arma::vec &beta,
//                              const arma::vec &y_resp,
//                              const arma::vec &X_beta,
//                              const arma::mat &transformed_X,
//                              const double &nu,
//                              const double &sigma,
//                              const arma::vec &prior_means,
//                              const arma::vec &prior_variances,
//                              const double &C,
//                              const arma::mat &sqrt_precondition_mat) {
//   arma::vec gradient(beta.size(), arma::fill::zeros);
//   const double nu_sigma_sq = nu*sigma*sigma;
//   for (int k=0; k < transformed_X.n_cols; ++k) {
//     double sum = 0;
//     for (int i=0; i < transformed_X.n_rows; ++i) {
//       const double diff = y_resp.at(i)-X_beta.at(i);
//       sum += (transformed_X.at(i,k)*diff) / (nu_sigma_sq + (diff*diff));
//     }
//     gradient.at(k) = (nu+1)*sum;
//     for (int j=0; j < transformed_X.n_cols; ++j) {
//       gradient.at(k) -= (sqrt_precondition_mat.at(j,k)*(beta.at(j)-prior_means.at(j)))/(C*prior_variances.at(j));
//     }
//   }
//   return(gradient);
// }
// 
// // [[Rcpp::export]]
// double term2_Z(const arma::vec &y_resp,
//                const arma::vec &X_beta,
//                const arma::mat &transformed_X,
//                const double &nu,
//                const double &sigma,
//                const arma::vec &prior_variances,
//                const double &C,
//                const arma::mat &precondition_mat,
//                const arma::mat &sqrt_precondition_mat) {
//   double divergence = 0;
//   const double nu_sigma_sq = nu*sigma*sigma;
//   for (int k=0; k < transformed_X.n_cols; ++k) {
//     for (int i=0; i < transformed_X.n_rows; ++i) {
//       double diff_sq = (y_resp.at(i)-X_beta.at(i))*(y_resp.at(i)-X_beta.at(i));
//       const double ratio = (diff_sq-nu_sigma_sq) /
//         ((diff_sq+nu_sigma_sq)*(diff_sq+nu_sigma_sq));
//       divergence += (nu+1)*(transformed_X.at(i,k)*transformed_X.at(i,k)*ratio);
//     }
//     for (int j=0; j < transformed_X.n_cols; ++j) {
//       divergence -= (sqrt_precondition_mat.at(j,k)*sqrt_precondition_mat.at(j,k))/(C*prior_variances.at(k));
//     }
//   }
//   return(divergence);
// }
// 
// // [[Rcpp::export]]
// Rcpp::List ea_phi_BRR_DL_vec_Z(const arma::vec &beta,
//                                const arma::vec &y_resp,
//                                const arma::mat &X,
//                                const double &nu,
//                                const double &sigma,
//                                const arma::vec &prior_means,
//                                const arma::vec &prior_variances,
//                                const double &C,
//                                const arma::mat &precondition_mat,
//                                const arma::mat &sqrt_precondition_mat) {
//   const arma::vec X_beta = X * beta;
//   const arma::mat transformed_X = X * sqrt_precondition_mat;
//   const arma::vec gradient = log_BRR_gradient_Z(beta,
//                                                 y_resp,
//                                                 X_beta,
//                                                 transformed_X,
//                                                 nu,
//                                                 sigma,
//                                                 prior_means,
//                                                 prior_variances,
//                                                 C,
//                                                 sqrt_precondition_mat);
//   const double t1 = arma::dot(gradient, gradient);
//   const double t2 = term2_Z(y_resp,
//                             X_beta,
//                             transformed_X,
//                             nu,
//                             sigma,
//                             prior_variances,
//                             C,
//                             precondition_mat,
//                             sqrt_precondition_mat);
//   return(Rcpp::List::create(Named("phi", 0.5*(t1+t2)),
//                             Named("t1", t1),
//                             Named("t2", t2)));
// }
// 
// // [[Rcpp::export]]
// double term2_X(const arma::vec &y_resp,
//                const arma::mat &X,
//                const arma::vec &X_beta,
//                const double &nu,
//                const double &sigma,
//                const arma::vec &prior_variances,
//                const double &C,
//                const arma::mat &precondition_mat) {
//   return(arma::trace(precondition_mat * log_BRR_hessian(y_resp,
//                                                         X,
//                                                         X_beta,
//                                                         nu,
//                                                         sigma,
//                                                         prior_variances,
//                                                         C)));
// }
