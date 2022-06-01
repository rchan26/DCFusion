#include "../inc/helper_functions.hpp"
#include "../inc/phi_BNBR.hpp"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec log_BNBR_gradient(const arma::vec &beta,
                            const arma::vec &y_count,
                            const arma::mat &X,
                            const arma::vec &X_beta,
                            const arma::vec &count,
                            const double &phi_rate,
                            const arma::vec &prior_means,
                            const arma::vec &prior_variances,
                            const double &C) {
  arma::vec gradient(beta.size(), arma::fill::zeros);
  for (int k=0; k < X.n_cols; ++k) {
    for (int i=0; i < X.n_rows; ++i) {
      const double exp_X_beta = exp(X_beta.at(i));
      gradient.at(k) += count.at(i)*X.at(i,k)*(y_count.at(i)-((y_count.at(i)+phi_rate)*exp_X_beta/(exp_X_beta+phi_rate)));
    }
    gradient.at(k) -= (beta.at(k)-prior_means.at(k))/(C*prior_variances.at(k));
  }
  return(gradient);
}

// [[Rcpp::export]]
arma::mat log_BNBR_hessian(const arma::vec &y_count,
                           const arma::mat &X,
                           const arma::vec &X_beta,
                           const arma::vec &count,
                           const double &phi_rate,
                           const arma::vec &prior_variances,
                           const double &C) {
  arma::mat hessian(X.n_cols, X.n_cols, arma::fill::zeros);
  for (int i=0; i < X.n_rows; ++i) {
    const double exp_X_beta = exp(X_beta.at(i));
    const double ratio = count.at(i)*(y_count.at(i)+phi_rate)*phi_rate*exp_X_beta/((exp_X_beta+phi_rate)*(exp_X_beta+phi_rate));
    for (int j=0; j < X.n_cols; ++j) {
      for (int k=0; k <= j; ++k) {
        hessian.at(j,k) -= X.at(i,j)*X.at(i,k)*ratio;
      }
    }
  }
  for (int j=0; j < X.n_cols; ++j) {
    hessian.at(j,j) -= 1/(C*prior_variances.at(j));
    for (int k=0; k < j; ++k) {
      hessian.at(k,j) = hessian.at(j,k);
    }
  }
  return(hessian);
}

// [[Rcpp::export]]
Rcpp::List ea_phi_BNBR_DL_vec(const arma::vec &beta,
                              const arma::vec &y_count,
                              const arma::mat &X,
                              const arma::vec &count,
                              const double &phi_rate,
                              const arma::vec &prior_means,
                              const arma::vec &prior_variances,
                              const double &C,
                              const arma::mat &precondition_mat) {
  const arma::vec X_beta = X * beta;
  const arma::vec gradient = log_BNBR_gradient(beta,
                                               y_count,
                                               X,
                                               X_beta,
                                               count,
                                               phi_rate,
                                               prior_means,
                                               prior_variances,
                                               C);
  const double t1 = as_scalar((arma::trans(gradient)*precondition_mat)*gradient);
  const arma::mat hessian = log_BNBR_hessian(y_count,
                                             X,
                                             X_beta,
                                             count,
                                             phi_rate,
                                             prior_variances,
                                             C);
  const double t2 = arma::trace(precondition_mat*hessian);
  return(Rcpp::List::create(Named("phi", 0.5*(t1+t2)),
                            Named("t1", t1),
                            Named("t2", t2)));
}

// [[Rcpp::export]]
Rcpp::List ea_phi_BNBR_DL_matrix(const arma::mat &beta,
                                 const arma::vec &y_count,
                                 const arma::mat &X,
                                 const arma::vec &count,
                                 const double &phi_rate,
                                 const arma::vec &prior_means,
                                 const arma::vec &prior_variances,
                                 const double &C,
                                 const arma::mat &precondition_mat) {
  Rcpp::NumericVector phi(beta.n_rows);
  Rcpp::NumericVector t1(beta.n_rows);
  Rcpp::NumericVector t2(beta.n_rows);
  for (int i=0; i < beta.n_rows; ++i) {
    Rcpp::List phi_eval = ea_phi_BNBR_DL_vec(arma::trans(beta.row(i)),
                                             y_count,
                                             X,
                                             count,
                                             phi_rate,
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
double spectral_radius_BNBR(const arma::vec &y_count,
                            const arma::vec &beta,
                            const arma::mat &X,
                            const arma::vec &count,
                            const double &phi_rate,
                            const arma::vec &prior_variances,
                            const double &C,
                            const arma::mat &Lambda) {
  const arma::vec X_beta = X * beta;
  return(spectral_radius(Lambda * log_BNBR_hessian(y_count,
                                                   X,
                                                   X_beta,
                                                   count,
                                                   phi_rate,
                                                   prior_variances,
                                                   C)));
}

// [[Rcpp::export]]
Rcpp::List obtain_hypercube_centre_BNBR(const Rcpp::List &bessel_layers,
                                        const arma::mat &transform_to_X,
                                        const arma::vec &y_count,
                                        const arma::mat &X,
                                        const arma::vec &count,
                                        const arma::vec &prior_means,
                                        const arma::vec &prior_variances,
                                        const double &C,
                                        const double &phi_rate) {
  arma::vec centre(bessel_layers.size(), arma::fill::zeros);
  for (int i=0; i < bessel_layers.size(); ++i) {
    const Rcpp::List &b_layer = bessel_layers[i];
    const double &L = b_layer["L"];
    const double &U = b_layer["U"];
    centre.at(i) = 0.5*(L+U);
  }
  const arma::vec beta_hat = transform_to_X * centre;
  const arma::vec X_beta = X * beta_hat;
  return(Rcpp::List::create(Named("centre", centre),
                            Named("beta_hat", beta_hat),
                            Named("grad_log_hat", log_BNBR_gradient(beta_hat,
                                                        y_count,
                                                        X,
                                                        X_beta,
                                                        count,
                                                        phi_rate,
                                                        prior_means,
                                                        prior_variances,
                                                        C))));
}

// [[Rcpp::export]]
double obtain_G_max(const int &dim,
                    const arma::vec &transformed_X_vec,
                    const Rcpp::List &bessel_layers,
                    const arma::vec &z_hat,
                    const double &phi_rate) {
  const double F_hat = arma::dot(transformed_X_vec, z_hat);
  const double log_phi_rate = log(phi_rate);
  if (F_hat > log_phi_rate) {
    // perform minimisation of F = transform_X_rowvec * z for z in bessel_layers
    const double F_min = optimise_vector_product(dim,
                                                 transformed_X_vec,
                                                 bessel_layers,
                                                 true);
    if (F_min < log_phi_rate) {
      return (0.25/phi_rate);
    } else {
      const double exp_u = exp(F_min);
      return exp_u/((1+exp_u)*(1+exp_u));
    }
  } else {
    // perform maximisation of F = transform_X_rowvec * z for z in bessel_layers
    const double F_max = optimise_vector_product(dim,
                                                 transformed_X_vec,
                                                 bessel_layers,
                                                 false);
    if (F_max > log_phi_rate) {
      return (0.25/phi_rate);
    } else {
      const double exp_u = exp(F_max);
      return exp_u/((1+exp_u)*(1+exp_u));
    }
  }
}

// [[Rcpp::export]]
Rcpp::List spectral_radius_bound_BNBR_Z(const int &dim,
                                        const Rcpp::List &bessel_layers,
                                        const arma::vec &z_hat,
                                        const arma::vec &y_count,
                                        const arma::mat &transformed_X,
                                        const arma::vec &count,
                                        const double &phi_rate,
                                        const arma::vec &prior_variances,
                                        const double &C,
                                        const arma::mat &sqrt_Lambda) {
  arma::mat hessian(dim, dim, arma::fill::zeros);
  for (int i=0; i < transformed_X.n_rows; ++i) {
    const double constant = count.at(i)*(y_count.at(i)+phi_rate)*phi_rate;
    const double G_max = obtain_G_max(dim, arma::trans(transformed_X.row(i)), bessel_layers, z_hat, phi_rate);
    for (int k=0; k < dim; ++k) {
      for (int l=0; l <= k; ++l) {
        hessian.at(k,l) += constant*std::abs(transformed_X.at(i,k))*std::abs(transformed_X.at(i,l))*G_max;
      }
    }
  }
  for (int k=0; k < dim; ++k) {
    for (int l=0; l <= k; ++l) {
      for (int j=0; j < dim; ++j) {
        hessian.at(k,l) += sqrt_Lambda.at(j,k)*sqrt_Lambda(j,l)/(C*prior_variances.at(j));
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
Rcpp::List spectral_radius_global_bound_BNBR_Z(const int &dim,
                                               const arma::vec &y_count,
                                               const arma::mat &transformed_X,
                                               const arma::vec &count,
                                               const double &phi_rate,
                                               const arma::vec &prior_variances,
                                               const double &C,
                                               const arma::mat &sqrt_Lambda) {
  arma::mat hessian(dim, dim, arma::fill::zeros);
  const double ratio = 0.25/phi_rate;
  for (int i=0; i < transformed_X.n_rows; ++i) {
    const double constant = count.at(i)*(y_count.at(i)+phi_rate)*phi_rate;
    for (int k=0; k < dim; ++k) {
      for (int l=0; l <= k; ++l) {
        hessian.at(k,l) += constant*std::abs(transformed_X.at(i,k))*std::abs(transformed_X.at(i,l))*ratio;
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
Rcpp::List ea_phi_BNBR_DL_bounds(const arma::vec &beta_hat,
                                 const arma::vec &grad_log_hat,
                                 const arma::vec &hypercube_centre_Z,
                                 const int &dim,
                                 const arma::vec &y_count,
                                 const arma::mat &transformed_X,
                                 const arma::vec &count,
                                 const double &phi_rate,
                                 const arma::vec &prior_variances,
                                 const double &C,
                                 const Rcpp::List &transform_mats,
                                 const Rcpp::List &hypercube_vertices,
                                 const Rcpp::List &bessel_layers,
                                 const bool &local_bounds) {
  const arma::mat &transform_to_X = transform_mats["to_X"];
  const arma::mat &transform_to_Z = transform_mats["to_Z"];
  const double vec_norm = std::sqrt(arma::sum(arma::square(transform_to_X*grad_log_hat)));
  const arma::mat &vertices = hypercube_vertices["vertices"];
  const double dist = maximal_distance_hypercube_to_cv(beta_hat,
                                                       vertices,
                                                       transform_to_X,
                                                       transform_to_Z,
                                                       true);
  Rcpp::List spectral_radius_bds;
  double P_n_Lambda;
  if (local_bounds) {
    spectral_radius_bds = spectral_radius_bound_BNBR_Z(dim,
                                                       bessel_layers,
                                                       hypercube_centre_Z,
                                                       y_count,
                                                       transformed_X,
                                                       count,
                                                       phi_rate,
                                                       prior_variances,
                                                       C,
                                                       transform_to_X);
    P_n_Lambda = spectral_radius_bds["spectral_radius"];
  } else {
    spectral_radius_bds = spectral_radius_global_bound_BNBR_Z(dim,
                                                              y_count,
                                                              transformed_X,
                                                              count,
                                                              phi_rate,
                                                              prior_variances,
                                                              C,
                                                              transform_to_X);
    P_n_Lambda = spectral_radius_bds["spectral_radius"];
  }
  return(Rcpp::List::create(Named("LB", -0.5*dim*P_n_Lambda),
                            Named("UB", 0.5*((vec_norm+dist*P_n_Lambda)*(vec_norm+dist*P_n_Lambda)+dim*P_n_Lambda)),
                            Named("dist", dist),
                            Named("P_n_Lambda", P_n_Lambda),
                            Named("t1_bds", (vec_norm+dist*P_n_Lambda)*(vec_norm+dist*P_n_Lambda)),
                            Named("t2_bds", dim*P_n_Lambda)));
}

// [[Rcpp::export]]
double gamma_NB_BNBR(const arma::vec &times,
                     const double &h,
                     const arma::vec &x0,
                     const arma::vec &y,
                     const double &s,
                     const double &t,
                     const arma::vec &y_count,
                     const arma::mat &X,
                     const arma::vec &count,
                     const double &phi_rate,
                     const arma::vec &prior_means,
                     const arma::vec &prior_variances,
                     const double &C,
                     const arma::mat &precondition_mat) {
  if (times.size() < 2) {
    stop("gamma_NB_BNBR: length of times must be at least 2"); 
  } 
  double sum_phi_eval = 0;
  for (int i=0; i < times.size(); ++i) {
    const arma::vec eval = (x0*(t-times.at(i))+y*(times.at(i)-s))/(t-s);
    Rcpp::List phi = ea_phi_BNBR_DL_vec(eval,
                                        y_count,
                                        X,
                                        count,
                                        phi_rate,
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
// arma::vec log_BNBR_gradient_Z(const arma::vec &beta,
//                               const arma::vec &y_count,
//                               const arma::vec &X_beta,
//                               const arma::mat &transformed_X,
//                               const arma::vec &count,
//                               const double &phi_rate,
//                               const arma::vec &prior_means,
//                               const arma::vec &prior_variances,
//                               const double &C,
//                               const arma::mat &sqrt_precondition_mat) {
//   arma::vec gradient(beta.size(), arma::fill::zeros);
//   for (int k=0; k < transformed_X.n_cols; ++k) {
//     for (int i=0; i < transformed_X.n_rows; ++i) {
//       const double exp_X_beta = exp(X_beta.at(i));
//       gradient.at(k) += count.at(i)*transformed_X.at(i,k)*(y_count.at(i)-((y_count.at(i)+phi_rate)*exp_X_beta/(exp_X_beta+phi_rate)));
//     }
//     for (int j=0; j < transformed_X.n_cols; ++j) {
//       gradient.at(k) -= (sqrt_precondition_mat.at(j,k)*(beta.at(j)-prior_means.at(j)))/(C*prior_variances.at(j));
//     }
//   }
//   return(gradient);
// }

// // [[Rcpp::export]]
// double BNBR_term2_Z(const arma::vec &y_count,
//                     const arma::vec &X_beta,
//                     const arma::mat &transformed_X,
//                     const arma::vec &count,
//                     const double &phi_rate,
//                     const arma::vec &prior_variances,
//                     const double &C,
//                     const arma::mat &transform_mat) {
//   double divergence = 0;
//   for (int k=0; k < transformed_X.n_cols; ++k) {
//     for (int i=0; i < transformed_X.n_rows; ++i) {
//       const double exp_X_beta = exp(X_beta.at(i));
//       const double ratio = count.at(i)*(y_count.at(i)+phi_rate)*phi_rate*exp_X_beta/((exp_X_beta+phi_rate)*(exp_X_beta+phi_rate));
//       divergence -= transformed_X.at(i,k)*transformed_X.at(i,k)*ratio;
//     }
//     for (int j=0; j < transformed_X.n_cols; ++j) {
//       divergence -= (transform_mat.at(j,k)*transform_mat.at(j,k))/(C*prior_variances.at(j));
//     }
//   }
//   return(divergence);
// }

// // [[Rcpp::export]]
// Rcpp::List ea_phi_BNBR_DL_vec_Z(const arma::vec &beta,
//                                 const arma::vec &y_count,
//                                 const arma::mat &X,
//                                 const arma::vec &count,
//                                 const double &phi_rate,
//                                 const arma::vec &prior_means,
//                                 const arma::vec &prior_variances,
//                                 const double &C,
//                                 const arma::mat &transform_mat) {
//   const arma::vec X_beta = X * beta;
//   const arma::mat transformed_X = X * transform_mat;
//   const arma::vec gradient = log_BNBR_gradient_Z(beta,
//                                                  y_count,
//                                                  X_beta,
//                                                  transformed_X,
//                                                  count,
//                                                  phi_rate,
//                                                  prior_means,
//                                                  prior_variances,
//                                                  C,
//                                                  transform_mat);
//   const double t1 = arma::dot(gradient, gradient);
//   const double t2 = BNBR_term2_Z(y_count,
//                                  X_beta,
//                                  transformed_X,
//                                  count,
//                                  phi_rate,
//                                  prior_variances,
//                                  C,
//                                  transform_mat);
//   return(Rcpp::List::create(Named("phi", 0.5*(t1+t2)),
//                             Named("t1", t1),
//                             Named("t2", t2)));
// }

// // [[Rcpp::export]]
// double BNBR_term2_X(const arma::vec &y_count,
//                     const arma::mat &X,
//                     const arma::vec &X_beta,
//                     const arma::vec &count,
//                     const double &phi_rate,
//                     const arma::vec &prior_variances,
//                     const double &C,
//                     const arma::mat &precondition_mat) {
//   const arma::mat hessian = log_BNBR_hessian(y_count,
//                                              X,
//                                              X_beta,
//                                              count,
//                                              phi_rate,
//                                              prior_variances,
//                                              C);
//   return(arma::trace(precondition_mat * hessian));
// }
