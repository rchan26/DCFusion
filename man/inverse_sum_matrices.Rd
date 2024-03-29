% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{inverse_sum_matrices}
\alias{inverse_sum_matrices}
\title{Calculate the inverse of a sum of matrices}
\usage{
inverse_sum_matrices(matrices)
}
\arguments{
\item{matrices}{list of matrices (of same dimension)}
}
\value{
the inverse of the sum of the matrices
}
\description{
Calculation of the inverse of a sum of a list of matrices
}
\examples{
m1 <- matrix(c(1,2,3,4), nrow = 2, ncol = 2)
m2 <- matrix(c(5,6,7,8), nrow = 2, ncol = 2)
m3 <- matrix(c(9,10,11,12), nrow = 2, ncol = 2)
inverse_sum_matrices(list(m1, m2, m3))
# returns the same as using solve() in R
solve(m1+m2+m3)
}
