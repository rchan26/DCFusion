% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{abs_eigenvals}
\alias{abs_eigenvals}
\title{Absolute eigenvalues of a matrix}
\usage{
abs_eigenvals(A)
}
\arguments{
\item{A}{matrix matrix}
}
\value{
The absolute eigenvalues of A
}
\description{
Calculates the absolute eigenvalues of a matrix
}
\examples{
# symmetric matrix
# should equal 2.5, 1.5
abs_eigenvals(matrix(c(2, 0.5, 0.5, 2), nrow = 2, ncol = 2))
# non-symmetric matrix
# should equal 10, 10, 5
abs_eigenvals(matrix(c(9, -1, 2, -2, 8, 4, 1, 1, 8), nrow = 3, ncol = 3, byrow = T))
}
