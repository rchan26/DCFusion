#include "../inc/helper_functions.hpp"
#include "../inc/phi_BLR.hpp"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec log_BLR_gradient(const arma::vec &beta,
                           const arma::vec &y_labels,
                           const arma::mat &X,
                           const arma::vec &X_beta,
                           const arma::vec &prior_means,
                           const arma::vec &prior_variances,
                           const double &C) {
  arma::vec gradient(beta.size(), arma::fill::zeros);
  for (int k=0; k < X.n_cols; ++k) {
    for (int i=0; i < X.n_rows; ++i) {
      gradient.at(k) += X.at(i,k)*(y_labels.at(i)-(1/(1+exp(-X_beta.at(i)))));
    }
    gradient.at(k) -= (beta.at(k)-prior_means.at(k))/(C*prior_variances.at(k));
  }
  return(gradient);
}

// [[Rcpp::export]]
arma::mat log_BLR_hessian(const arma::mat &X,
                          const arma::vec &X_beta,
                          const arma::vec &prior_variances,
                          const double &C) {
  arma::mat hessian(X.n_cols, X.n_cols, arma::fill::zeros);
  for (int i=0; i < X.n_rows; ++i) {
    const double exp_X_beta = exp(X_beta.at(i));
    const double ratio = exp_X_beta/((1+exp_X_beta)*(1+exp_X_beta));
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
Rcpp::List ea_phi_BLR_DL_vec(const arma::vec &beta,
                             const arma::vec &y_labels,
                             const arma::mat &X,
                             const arma::vec &prior_means,
                             const arma::vec &prior_variances,
                             const double &C,
                             const arma::mat &precondition_mat,
                             const arma::mat &transform_mat) {
  const arma::vec transformed_beta = transform_mat * beta;
  const arma::vec X_beta = X * transformed_beta;
  const arma::vec gradient = log_BLR_gradient(transformed_beta,
                                              y_labels,
                                              X,
                                              X_beta,
                                              prior_means,
                                              prior_variances,
                                              C);
  const double t1 = as_scalar((arma::trans(gradient)*precondition_mat)*gradient);
  const arma::mat hessian = log_BLR_hessian(X, X_beta, prior_variances, C);
  const double t2 = arma::trace(precondition_mat*hessian);
  return(Rcpp::List::create(Named("phi", 0.5*(t1+t2)),
                            Named("t1", t1),
                            Named("t2", t2)));
}

// [[Rcpp::export]]
Rcpp::List ea_phi_BLR_DL_matrix(const arma::mat &beta,
                                const arma::vec &y_labels,
                                const arma::mat &X,
                                const arma::vec &prior_means,
                                const arma::vec &prior_variances,
                                const double &C,
                                const arma::mat &precondition_mat,
                                const arma::mat &transform_mat) {
  Rcpp::NumericVector phi(beta.n_rows);
  Rcpp::NumericVector t1(beta.n_rows);
  Rcpp::NumericVector t2(beta.n_rows);
  for (int i=0; i < beta.n_rows; ++i) {
    Rcpp::List phi_eval = ea_phi_BLR_DL_vec(arma::trans(beta.row(i)),
                                            y_labels,
                                            X,
                                            prior_means,
                                            prior_variances,
                                            C,
                                            precondition_mat,
                                            transform_mat);
    phi[i] = phi_eval["phi"];
    t1[i] = phi_eval["t1"];
    t2[i] = phi_eval["t2"];
  }
  return(Rcpp::List::create(Rcpp::Named("phi", phi),
                            Rcpp::Named("t1", t1),
                            Rcpp::Named("t2", t2)));
}

// [[Rcpp::export]]
double spectral_radius_BLR(const arma::vec &beta,
                           const int &dim,
                           const arma::mat &X,
                           const arma::vec &prior_variances,
                           const double &C,
                           const arma::mat &Lambda) {
  arma::mat hessian(dim, dim, arma::fill::zeros);
  const arma::vec X_beta = X * beta;
  for (int i=0; i < X.n_rows; ++i) {
    const double exp_X_beta = exp(X_beta.at(i));
    const double ratio = exp_X_beta/((1+exp_X_beta)*(1+exp_X_beta));
    for (int j=0; j < dim; ++j) {
      for (int k=0; k <= j; ++k) {
        hessian.at(j,k) -= X.at(i,j)*X.at(i,k)*ratio;
      }
    }
  }
  for (int j=0; j < dim; ++j) {
    hessian.at(j,j) -= 1/(C*prior_variances.at(j));
    for (int k=0; k < j; ++k) {
      hessian.at(k,j) = hessian.at(j,k);
    }
  }
  return(spectral_radius(Lambda * hessian));
}

// [[Rcpp::export]]
Rcpp::List max_multiplication(const arma::mat &matrix,
                              const Rcpp::List &bessel_layers) {
  if (matrix.n_cols!=bessel_layers.size()) {
    stop("max_multiplication: number of columns of matrix must equal the length of bessel_layers");
  }
  arma::vec max_mult(matrix.n_rows, arma::fill::zeros);
  arma::vec min_mult(matrix.n_rows, arma::fill::zeros);
  for (int j=0; j < matrix.n_cols; ++j) {
    const Rcpp::List &bes_layer = bessel_layers[j];
    const double &lower = bes_layer["L"];
    const double &upper = bes_layer["U"];
    for (int i=0; i < matrix.n_rows; ++i) {
      const double mat_lb = matrix.at(i,j) * lower;
      const double mat_ub = matrix.at(i,j) * upper;
      if (mat_ub > mat_lb) {
        max_mult.at(i) += mat_ub;
        min_mult.at(i) += mat_lb;
      } else {
        max_mult.at(i) += mat_lb;
        min_mult.at(i) += mat_ub;
      }
    }
  }
  return(Rcpp::List::create(Named("max", max_mult),
                            Named("min", min_mult),
                            Named("max_abs", arma::abs(max_mult)),
                            Named("min_abs", arma::abs(min_mult))));
}

// [[Rcpp::export]]
Rcpp::List spectral_radius_bound_BLR_Z(const int &dim,
                                       const Rcpp::List &bessel_layers,
                                       const arma::mat &X,
                                       const arma::vec &prior_variances,
                                       const double &C,
                                       const arma::mat &sqrt_Lambda) {
  arma::mat hessian(dim, dim, arma::fill::zeros);
  // obtain the lower and upper bound on X_beta
  const arma::mat transformed_X = X * sqrt_Lambda;
  const Rcpp::List u_bound = max_multiplication(transformed_X, bessel_layers);
  const arma::vec LB = u_bound["min_abs"];
  const arma::vec UB = u_bound["max_abs"];
  for (int i=0; i < X.n_rows; ++i) {
    // e^u/((1+e^u)^2) is largest when x is closest to 0, hence take the smaller of the bounds
    const double exp_u = exp(std::min(LB.at(i), UB.at(i)));
    const double ratio = exp_u/((1+exp_u)*(1+exp_u));
    for (int k=0; k < dim; ++k) {
      for (int l=0; l <= k; ++l) {
        hessian.at(k,l) -= transformed_X.at(i,k)*transformed_X.at(i,l)*ratio;
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
Rcpp::List spectral_radius_global_bound_BLR_Z(const int &dim,
                                              const arma::mat &X,
                                              const arma::vec &prior_variances,
                                              const double &C,
                                              const arma::mat &sqrt_Lambda) {
  // ----- compute the matrix A = Hessian of the transformed log sub-posterior
  const arma::mat transformed_X = X * sqrt_Lambda;
  arma::mat hessian(dim, dim, arma::fill::zeros);
  for (int i=0; i < X.n_rows; ++i) {
    for (int k=0; k < dim; ++k) {
      for (int l=0; l <= k; ++l) {
        hessian.at(k,l) -= transformed_X.at(i,k)*transformed_X.at(i,l)/4;
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
Rcpp::List obtain_hypercube_centre(const Rcpp::List &bessel_layers,
                                   const arma::mat &transform_to_X,
                                   const arma::vec &y_labels,
                                   const arma::mat &X,
                                   const arma::vec &prior_means,
                                   const arma::vec &prior_variances,
                                   const double &C) {
  // calculate the hypercube centre
  arma::vec centre(bessel_layers.size(), arma::fill::zeros);
  for (int i=0; i < bessel_layers.size(); ++i) {
    const Rcpp::List &b_layer = bessel_layers[i];
    const double &L = b_layer["L"];
    const double &U = b_layer["U"];
    centre.at(i) = 0.5*(L+U);
  }
  const arma::vec beta_hat = transform_to_X*centre;
  // evaluate the gradient of the log posterior at the hypercube centre
  const arma::vec X_beta = X * beta_hat;
  return(Rcpp::List::create(Named("beta_hat", beta_hat),
                            Named("grad_log_hat", log_BLR_gradient(beta_hat,
                                                       y_labels,
                                                       X,
                                                       X_beta,
                                                       prior_means,
                                                       prior_variances,
                                                       C))));
}

// [[Rcpp::export]]
double maximal_distance_hypercube_to_cv(const arma::vec &beta_hat,
                                        const arma::mat &hypercube_vertices,
                                        const arma::mat &transform_to_X,
                                        const arma::mat &transform_to_Z) {
  const arma::mat hypercube_X = arma::trans(transform_to_X * arma::trans(hypercube_vertices));
  arma::vec distances(hypercube_vertices.n_rows, arma::fill::zeros);
  for (int i=0; i < hypercube_X.n_rows; ++i) {
    distances.at(i) = scaled_distance(arma::trans(hypercube_X.row(i)), beta_hat, transform_to_Z);
  }
  return(distances.max());
}

// [[Rcpp::export]]
Rcpp::List ea_phi_BLR_DL_bounds(const arma::vec &beta_hat,
                                const arma::vec &grad_log_hat,
                                const int &dim,
                                const arma::mat &X,
                                const arma::vec &prior_variances,
                                const double &C,
                                const Rcpp::List &transform_mats,
                                const Rcpp::List &bessel_layers,
                                const arma::mat &hypercube_vertices) {
  const arma::mat &transform_to_X = transform_mats["to_X"];
  const arma::mat &transform_to_Z = transform_mats["to_Z"];
  const double vec_norm = std::sqrt(arma::sum(arma::square(transform_to_X*grad_log_hat)));
  const double dist = maximal_distance_hypercube_to_cv(beta_hat,
                                                       hypercube_vertices,
                                                       transform_to_X,
                                                       transform_to_Z);
  const Rcpp::List spectral_radius_bds = spectral_radius_bound_BLR_Z(dim,
                                                                     bessel_layers,
                                                                     X,
                                                                     prior_variances,
                                                                     C,
                                                                     transform_to_X);
  const double P_n_Lambda = spectral_radius_bds["spectral_radius"];
  // const arma::vec abs_eigenvalues = spectral_radius_bds["abs_eigenvals"];
  return(Rcpp::List::create(Named("LB", -0.5*dim*P_n_Lambda),
                            Named("UB", 0.5*((vec_norm+dist*P_n_Lambda)*(vec_norm+dist*P_n_Lambda)+dim*P_n_Lambda)),
                            Named("dist", dist),
                            Named("P_n_Lambda", P_n_Lambda),
                            Named("t1_bds", (vec_norm+dist*P_n_Lambda)*(vec_norm+dist*P_n_Lambda)),
                            Named("t2_bds", dim*P_n_Lambda)));
}

// [[Rcpp::export]]
double gamma_NB_estimate_BLR(const arma::vec &times,
                             const double &h,
                             const double &s,
                             const double &t,
                             const arma::vec &x0,
                             const arma::vec &y,
                             const arma::vec &y_labels,
                             const arma::mat &X,
                             const arma::vec &prior_means,
                             const arma::vec &prior_variances,
                             const double &C,
                             const arma::mat &precondition_mat,
                             const arma::mat &transform_mat) {
  if (times.size() < 2) {
    stop("gamma_NB_estimate_BLR: length of times must be at least 2"); 
  } 
  Rcpp::NumericVector phi(times.size());
  double sum_phi_eval = 0;
  for (int i=0; i < times.size(); ++i) {
    const arma::vec eval = (x0*(t-times.at(i)) + y*times.at(i))/(t-s);
    Rcpp::List phi = ea_phi_BLR_DL_vec(eval,
                                       y_labels,
                                       X,
                                       prior_means,
                                       prior_variances,
                                       C,
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
// double ea_phi_BLR_DL_LB(const arma::mat &X,
//                         const arma::vec &prior_variances,
//                         const double &C,
//                         const arma::mat &precondition_mat) {
//   arma::vec design_sum(X.n_cols, arma::fill::zeros);
//   double prior_variances_sum = 0;
//   for (int k=0; k < X.n_cols; ++k) {
//     for (int i=0; i < X.n_rows; ++i) {
//       design_sum.at(k) += X.at(i,k)*X.at(i,k);
//     }
//     design_sum.at(k) *= precondition_mat.at(k,k);
//     prior_variances_sum += precondition_mat.at(k,k)/(C*prior_variances.at(k));
//   }
//   double LB = -(arma::sum(design_sum)/4) - prior_variances_sum;
//   return(0.5*LB);
// }

// // [[Rcpp::export]]
// arma::vec log_BLR_gradient_Z(const arma::vec &beta,
//                              const arma::vec &y_labels,
//                              const arma::vec &X_beta,
//                              const arma::mat &transformed_X,
//                              const arma::vec &prior_means,
//                              const arma::vec &prior_variances,
//                              const double &C,
//                              const arma::mat &sqrt_precondition_mat) {
//   arma::vec gradient(beta.size(), arma::fill::zeros);
//   for (int k=0; k < transformed_X.n_cols; ++k) {
//     for (int i=0; i < transformed_X.n_rows; ++i) {
//       gradient.at(k) += transformed_X.at(i,k)*(y_labels.at(i)-(1/(1+exp(-X_beta.at(i)))));
//     }
//     for (int j=0; j < transformed_X.n_cols; ++j) {
//       gradient.at(k) -= (sqrt_precondition_mat.at(j,k)*(beta.at(j)-prior_means.at(j)))/(C*prior_variances.at(j));
//     }
//   }
//   return(gradient);
// }

// // [[Rcpp::export]]
// double term2_Z(const arma::vec &X_beta,
//                const arma::mat &transformed_X,
//                const arma::vec &prior_variances,
//                const double &C,
//                const arma::mat &precondition_mat,
//                const arma::mat &transform_mat) {
//   double divergence = 0;
//   for (int k=0; k < transformed_X.n_cols; ++k) {
//     for (int i=0; i < transformed_X.n_rows; ++i) {
//       const double exp_X_beta = exp(X_beta.at(i));
//       const double ratio = exp_X_beta/((1+exp_X_beta)*(1+exp_X_beta));
//       divergence -= transformed_X.at(i,k)*transformed_X.at(i,k)*ratio;
//     }
//     for (int j=0; j < transformed_X.n_cols; ++j) {
//       divergence -= (transform_mat.at(j,k)*transform_mat.at(j,k))/(C*prior_variances.at(k));
//     }
//   }
//   return(divergence);
// }

// // [[Rcpp::export]]
// double ea_phi_BLR_DL_vec_Z(const arma::vec &beta,
//                            const arma::vec &y_labels,
//                            const arma::mat &X,
//                            const arma::vec &prior_means,
//                            const arma::vec &prior_variances,
//                            const double &C,
//                            const arma::mat &precondition_mat,
//                            const arma::mat &transform_mat) {
//   const arma::vec X_beta = X * beta;
//   const arma::mat transformed_X = X * transform_mat;
//   const arma::vec gradient = log_BLR_gradient_Z(beta,
//                                                 y_labels,
//                                                 X_beta,
//                                                 transformed_X,
//                                                 prior_means,
//                                                 prior_variances,
//                                                 C,
//                                                 precondition_mat,
//                                                 transform_mat);
//   const double t1 = arma::dot(gradient, gradient);
//   const double t2 = term2_Z(X_beta,
//                             transformed_X,
//                             prior_variances,
//                             C,
//                             precondition_mat,
//                             transform_mat);
//   return(0.5*(t1+t2));
// }

// // [[Rcpp::export]]
// double term2_X(const arma::mat &X,
//                const arma::vec &X_beta,
//                const arma::vec &prior_variances,
//                const double &C,
//                const arma::mat &precondition_mat) {
//   const arma::mat hessian = log_BLR_hessian(X, X_beta, prior_variances, C);
//   return(arma::trace(precondition_mat * hessian));
// }
