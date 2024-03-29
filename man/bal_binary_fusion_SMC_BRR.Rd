% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BRR_fusion.R
\name{bal_binary_fusion_SMC_BRR}
\alias{bal_binary_fusion_SMC_BRR}
\title{(Balanced Binary) D&C Monte Carlo Fusion using SMC}
\usage{
bal_binary_fusion_SMC_BRR(
  N_schedule,
  m_schedule,
  time_schedule,
  base_samples,
  L,
  dim,
  data_split,
  nu,
  sigma,
  prior_means,
  prior_variances,
  C,
  precondition = TRUE,
  resampling_method = "multi",
  ESS_threshold = 0.5,
  diffusion_estimator = "Poisson",
  beta_NB = 10,
  gamma_NB_n_points = 2,
  seed = NULL,
  n_cores = parallel::detectCores(),
  print_progress_iters = 1000
)
}
\arguments{
\item{N_schedule}{vector of length (L-1), where N_schedule[l] is the
number of samples per node at level l}

\item{m_schedule}{vector of length (L-1), where m_schedule[k] is the number
of samples to fuse for level k}

\item{time_schedule}{vector of length (L-1), where time_schedule[k] is time
T for algorithm for level k}

\item{base_samples}{list of length C, where base_samples[[c]] contains
the samples for the c-th node in the level}

\item{L}{total number of levels in the hierarchy}

\item{dim}{dimension of the predictors (= p+1)}

\item{data_split}{list of length m where each item is a list of length 4 where
for c=1,...,m, data_split[[c]]$y is the vector for y responses and
data_split[[c]]$X is the design matrix for the covariates for
sub-posterior c}

\item{nu}{degrees of freedom in t-distribution}

\item{sigma}{scale parameter in t-distribution}

\item{prior_means}{prior for means of predictors}

\item{prior_variances}{prior for variances of predictors}

\item{C}{number of sub-posteriors at the base level}

\item{precondition}{either a logical value to determine if preconditioning matrices are
used (TRUE - and is set to be the variance of the sub-posterior samples)
or not (FALSE - and is set to be the identity matrix for all sub-posteriors),
or a list of length (1/start_beta) where precondition[[c]]
is the preconditioning matrix for sub-posterior c. Default is TRUE}

\item{resampling_method}{method to be used in resampling, default is
multinomial resampling ('multi'). Other choices are
stratified ('strat'), systematic ('system'),
residual ('resid')}

\item{ESS_threshold}{number between 0 and 1 defining the proportion
of the number of samples that ESS needs to be
lower than for resampling (i.e. resampling is carried
out only when ESS < N*ESS_threshold)}

\item{diffusion_estimator}{choice of unbiased estimator for the Exact Algorithm
between "Poisson" (default) for Poisson estimator
and "NB" for Negative Binomial estimator}

\item{beta_NB}{beta parameter for Negative Binomial estimator (default 10)}

\item{gamma_NB_n_points}{number of points used in the trapezoidal estimation
of the integral found in the mean of the negative
binomial estimator (default is 2)}

\item{seed}{seed number - default is NULL, meaning there is no seed}

\item{n_cores}{number of cores to use}

\item{print_progress_iters}{number of iterations between each progress update
(default is 1000). If NULL, progress will only
be updated when importance sampling is finished}
}
\value{
A list with components:
\describe{
  \item{particles}{list of length (L-1), where particles[[l]][[i]] are the
                   particles for level l, node i}
  \item{proposed_samples}{list of length (L-1), where proposed_samples[[l]][[i]]
                          are the proposed samples for level l, node i}
  \item{time}{list of length (L-1), where time[[l]][[i]] is the run time
              for level l, node i}
  \item{ESS}{list of length (L-1), where ESS[[l]][[i]] is the effective
             sample size of the particles after each step BEFORE deciding
             whether or not to resample for level l, node i}
  \item{CESS}{list of length (L-1), where ESS[[l]][[i]] is the conditional
              effective sample size of the particles after each step}
  \item{resampled}{list of length (L-1), where resampled[[l]][[i]] is a
                   boolean value to record if the particles were resampled
                   after each step; rho and Q for level l, node i}
  \item{precondition_matrices}{pre-conditioning matrices that were used}
  \item{data_inputs}{list of length (L-1), where data_inputs[[l]][[i]] is the
                     data input for the sub-posterior in level l, node i}
  \item{diffusion_times}{vector of length (L-1), where diffusion_times[l]
                         are the times for fusion in level l}
}
}
\description{
(Balanced Binary) D&C Monte Carlo Fusion using SMC for Bayesian Logistic Regression
}
