% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BNBR_fusion.R
\name{ea_BNBR_DL_PT}
\alias{ea_BNBR_DL_PT}
\title{Diffusion probability for the Exact Algorithm for Langevin diffusion for
Bayesian Negative Binomial regression}
\usage{
ea_BNBR_DL_PT(
  dim,
  x0,
  y,
  s,
  t,
  data,
  transformed_design_mat,
  phi_rate,
  prior_means,
  prior_variances,
  C,
  precondition_mat,
  transform_mats,
  diffusion_estimator,
  beta_NB = 10,
  gamma_NB_n_points = 2,
  local_bounds = TRUE,
  logarithm
)
}
\arguments{
\item{dim}{dimension of the predictors (= p+1)}

\item{x0}{start value (vector of length dim)}

\item{y}{end value (vector of length dim)}

\item{s}{start time}

\item{t}{end time}

\item{data}{list of length 4 where data[[c]]$y is the vector for y responses 
and data[[c]]$X is the design matrix for the covariates for
sub-posterior c, data[[c]]$full_data_count is the unique
rows of the full data set with their counts and 
data[[c]]$design_count is the unique rows of the design
matrix and their counts}

\item{phi_rate}{rate parameter for NB distribution}

\item{prior_means}{prior for means of predictors}

\item{prior_variances}{prior for variances of predictors}

\item{C}{overall number of sub-posteriors}

\item{precondition_mat}{precondition matrix}

\item{transform_mats}{list of transformation matrices where 
transform_mats$to_Z is the transformation matrix to Z space
and transform_mats$to_X is the transformation matrix to 
X space}

\item{diffusion_estimator}{choice of unbiased estimator for the Exact Algorithm
between "Poisson" (default) for Poisson estimator
and "NB" for Negative Binomial estimator}

\item{beta_NB}{beta parameter for Negative Binomial estimator (default 10)}

\item{gamma_NB_n_points}{number of points used in the trapezoidal estimation
of the integral found in the mean of the negative
binomial estimator (default is 2)}

\item{local_bounds}{logical value indicating if local bounds for the phi function
are used (default is TRUE)}

\item{logarithm}{logical value to determine if log probability is
returned (TRUE) or not (FALSE)}
}
\value{
acceptance probability of simulating Langevin diffusion where target
        is the posterior for a Negative Binomial regression model with Gaussian priors
}
\description{
Simulate Langevin diffusion using the Exact Algorithm where target is the
posterior for a Negative Binomial regression model with Gaussian priors
}
