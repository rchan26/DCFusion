#' Tempered Gaussian Distribution
#'
#' Density of tempered Gaussian distribution
#'
#' @param x real value
#' @param mean mean value
#' @param sd standard deviation value
#' @param beta real value
#'
#' @return density of tempered Gaussian distribution
#'
#' @examples
#' # standard normal
#' curve(dnorm(x), -5, 5)
#' # standard normal to the power 1/2
#' curve(dnorm_tempered(x, mean = 0, sd = 1, beta = 1/2), -5, 5, add = T, col = 'blue')
#'
#' @export
dnorm_tempered <- function(x, mean = 0, sd = 1, beta) {
  return(dnorm(x, mean = mean, sd = sd/sqrt(beta)))
}

#' Tempered Gaussian Distribution Sampler
#'
#' Samples from a tempered Gaussian distribution
#'
#' @param x real value
#' @param mean mean value
#' @param sd standard deviation value
#' @param beta real value
#'
#' @return N samples from a tempered Gaussian distribution
#'
#' @examples
#' # standard normal
#' curve(dnorm(x), -5, 5)
#' # standard normal to the power 1/2
#' curve(dnorm_tempered(x, beta = 1/2), -5, 5, add = T, col = 'blue')
#' samples <- rnorm_tempered(N = 10000, mean = 0, sd = 1, beta = 1/2)
#' lines(density(samples), col = 'red')
#'
#' @export
rnorm_tempered <- function(N, mean, sd, beta) {
  return(rnorm(n = N, mean = mean, sd = sd/sqrt(beta)))
}
