% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{ea_phi_biGaussian_DL_bounds}
\alias{ea_phi_biGaussian_DL_bounds}
\title{Obtain bounds for phi function}
\usage{
ea_phi_biGaussian_DL_bounds(
  mean_vec,
  sd_vec,
  corr,
  beta,
  precondition_mat,
  transform_to_Z,
  transform_to_X,
  lower,
  upper
)
}
\arguments{
\item{mean_vec}{vector of length 2 for mean}

\item{sd_vec}{vector of length 2 for standard deviation}

\item{corr}{correlation value between component 1 and component 2}

\item{beta}{real value}

\item{precondition_mat}{precondition matrix}

\item{transform_to_Z}{the transformation matrix to Z-space}

\item{transform_to_X}{the transformation matrix to X-space}

\item{lower}{vector of length 2 for the lower end of interval}

\item{upper}{vector of length 2 for the upper end of interval}
}
\value{
A list of components
\describe{
  \item{LB}{lower bound of phi}
  \item{UB}{upper bound of phi}
}
}
\description{
Finds the lower and upper bounds of the phi function between an interval
}
