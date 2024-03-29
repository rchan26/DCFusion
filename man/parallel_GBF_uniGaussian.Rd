% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/univariate_Gaussian_generalised_BF.R
\name{parallel_GBF_uniGaussian}
\alias{parallel_GBF_uniGaussian}
\title{Generalised Bayesian Fusion [parallel]}
\usage{
parallel_GBF_uniGaussian(
  particles_to_fuse,
  N,
  m,
  time_mesh,
  means,
  sds,
  betas,
  precondition_values,
  resampling_method = "multi",
  ESS_threshold = 0.5,
  sub_posterior_means = NULL,
  adaptive_mesh = FALSE,
  adaptive_mesh_parameters = list(),
  diffusion_estimator = "Poisson",
  beta_NB = 10,
  gamma_NB_n_points = 2,
  seed = NULL,
  n_cores = parallel::detectCores()
)
}
\arguments{
\item{particles_to_fuse}{list of length m, where particles_to_fuse[[c]]
contains the particles for the c-th sub-posterior
(a list of particles to fuse can be initialised by 
initialise_particle_sets() function)}

\item{N}{number of samples}

\item{m}{number of sub-posteriors to combine}

\item{time_mesh}{time mesh used in Bayesian Fusion}

\item{means}{vector of length m, where means[c] is the mean for c-th 
sub-posterior}

\item{sds}{vector of length m, where sds[c] is the standard deviation 
for c-th sub-posterior}

\item{betas}{vector of length c, where betas[c] is the inverse temperature 
value for c-th posterior}

\item{precondition_values}{vector of length m, where precondition_values[c]
is the precondition value for sub-posterior c}

\item{resampling_method}{method to be used in resampling, default is multinomial 
resampling ('multi'). Other choices are stratified 
resampling ('strat'), systematic resampling ('system'),
residual resampling ('resid')}

\item{ESS_threshold}{number between 0 and 1 defining the proportion 
of the number of samples that ESS needs to be
lower than for resampling (i.e. resampling is carried 
out only when ESS < N*ESS_threshold)}

\item{sub_posterior_means}{vector of length m, where sub_posterior_means[c]
is the sub-posterior mean of sub-posterior c}

\item{adaptive_mesh}{logical value to indicate if an adaptive mesh is used
(default is FALSE)}

\item{adaptive_mesh_parameters}{list of parameters used for adaptive mesh}

\item{diffusion_estimator}{choice of unbiased estimator for the Exact Algorithm
between "Poisson" (default) for Poisson estimator
and "NB" for Negative Binomial estimator}

\item{beta_NB}{beta parameter for Negative Binomial estimator (default 10)}

\item{gamma_NB_n_points}{number of points used in the trapezoidal estimation
of the integral found in the mean of the negative
binomial estimator (default is 2)}

\item{seed}{seed number - default is NULL, meaning there is no seed}

\item{n_cores}{number of cores to use}
}
\value{
A list with components:
\describe{
  \item{particles}{particles returned from fusion sampler}
  \item{proposed_samples}{proposal samples from fusion sampler}
  \item{time}{run-time of fusion sampler}
  \item{elapsed_time}{elapsed time of each step of the algorithm}
  \item{ESS}{effective sample size of the particles after each step}
  \item{CESS}{conditional effective sample size of the particles after each step}
  \item{resampled}{boolean value to indicate if particles were resampled
                   after each time step}
  \item{E_nu_j}{approximation of the average variation of the trajectories
                at each time step}
  \item{precondition_values}{list of length 2 where precondition_values[[2]] 
                             are the pre-conditioning values that were used 
                             and precondition_values[[1]] are the combined 
                             precondition values}
  \item{sub_posterior_means}{list of length 2, where sub_posterior_means[[2]]
                             are the sub-posterior means that were used and
                             sub_posterior_means[[1]] are the combined
                             sub-posterior means}
}
}
\description{
Generalised Bayesian Fusion with univariate Gaussian target
}
