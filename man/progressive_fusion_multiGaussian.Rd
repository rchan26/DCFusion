% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multivariate_Gaussian_fusion.R
\name{progressive_fusion_multiGaussian}
\alias{progressive_fusion_multiGaussian}
\title{(Progressive) D&C Monte Carlo Fusion (rejection sampling)}
\usage{
progressive_fusion_multiGaussian(
  N_schedule,
  time_schedule,
  dim,
  mu,
  Sigma,
  start_beta,
  base_samples,
  precondition = TRUE,
  seed = NULL,
  n_cores = parallel::detectCores()
)
}
\arguments{
\item{N_schedule}{vector of length (L-1), where N_schedule[l] is the number 
of samples per node at level l}

\item{time_schedule}{vector of length(L-1), where time_schedule[l] is the 
time chosen for Fusion at level l}

\item{dim}{dimension}

\item{mu}{vector of length dim for mean}

\item{Sigma}{dim x dim covariance matrix}

\item{start_beta}{beta for the base level}

\item{base_samples}{list of length (1/start_beta), where base_samples[[c]] 
contains the samples for the c-th node in the level}

\item{precondition}{either a logical value to determine if preconditioning matrices are
used (TRUE - and is set to be the variance of the sub-posterior samples)
or not (FALSE - and is set to be the identity matrix for all sub-posteriors),
or a list of length (1/start_beta) where precondition[[c]]
is the preconditioning matrix for sub-posterior c. Default is TRUE}

\item{seed}{seed number - default is NULL, meaning there is no seed}

\item{n_cores}{number of cores to use}
}
\value{
A list with components:
\describe{
 \item{samples}{list of length (L-1), where samples[[l]][[i]] are the samples 
                for level l, node i}
 \item{time}{list of length (L-1), where time[[l]][[i]] is the run time for 
             level l, node i}
 \item{rho_acc}{list of length (L-1), where rho_acc[[l]][i] is the acceptance 
                rate for first fusion step for level l, node i}
 \item{Q_acc}{list of length (L-1), where Q_acc[[l]][i] is the acceptance 
              rate for second fusion step for level l, node i}
 \item{rhoQ_acc}{list of length (L-1), where rhoQ_acc[[l]][i] is the overall 
                 acceptance rate for fusion for level l, node i}
 \item{diffusion_times}{vector of length (L-1), where diffusion_times[l] are 
                        the times for fusion in level l (= time_schedule)}
 \item{precondition_matrices}{preconditioning matrices used in the algorithm 
                              for each node}
}
}
\description{
(Progressive) D&C Monte Carlo Fusion with multivariate Gaussian target
}
