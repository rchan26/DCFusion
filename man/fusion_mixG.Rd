% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mixG_fusion.R
\name{fusion_mixG}
\alias{fusion_mixG}
\title{Generalised Monte Carlo Fusion (rejection sampling) [on a single core]}
\usage{
fusion_mixG(
  N,
  m,
  time,
  samples_to_fuse,
  n_comp,
  weights,
  means,
  sds,
  betas,
  precondition_values,
  bounds_multiplier = 1.1,
  level = 1,
  node = 1,
  core = 1
)
}
\arguments{
\item{N}{number of samples}

\item{m}{number of sub-posteriors to combine}

\item{time}{time T for fusion algorithm}

\item{samples_to_fuse}{list of length m, where samples_to_fuse[c] containing
the samples for the c-th sub-posterior}

\item{n_comp}{integer number of components of mixture Gaussian}

\item{weights}{vector: weights of mixture Gaussian}

\item{means}{vector: means of mixture Gaussian}

\item{sds}{vector: st.devs of mixture Gaussian}

\item{betas}{vector of length m, where betas[c] is the inverse temperature
(beta) for c-th sub-posterior (can also pass in one number if
they are all at the same inverse temperature)}

\item{precondition_values}{vector of length m, where precondition_values[c]
is the precondition value for sub-posterior c}

\item{bounds_multiplier}{scalar value to multiply bounds by 
(should greater than or equal to 1)}

\item{level}{integer for level in hierarchy - default 1}

\item{node}{integer for node in level in hierarchy - default 1}

\item{seed}{seed number - default is NULL, meaning there is no seed}

\item{n_cores}{number of cores to use}
}
\value{
A list with components:
\describe{
  \item{samples}{fusion samples}
  \item{iters_rho}{number of iterations for rho step}
  \item{iters_Q}{number of iterations for Q step}
}
}
\description{
Generalised Monte Carlo Fusion with mixture Gaussian target
}
