% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transformations.R
\name{center_matrix}
\alias{center_matrix}
\title{Simple function to center data}
\usage{
center_matrix(mat, centers = NULL)
}
\arguments{
\item{mat}{a transformed data matrix.}

\item{centers}{a vector of centering values. Length must be equal to the
number of columns in \code{mat}. If \code{NULL} (default) then
centering values are determined by the function (see Details).}
}
\value{
Transformed data matrix. Has an attribute \code{scaled:center} that stores
values used for centering.
}
\description{
Simple function to center data
}
\details{
Centering is done by means for continuous variables (i.e. more than 2
distinct values), and the minimum for binary variables.

It is assumed all categorical variables in the data are represented by
binary dummy variables.
}
\examples{
mat = matrix(1:100, nrow = 10)
center_matrix(mat)

}
