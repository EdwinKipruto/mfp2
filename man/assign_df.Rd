% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mfp2.R
\name{assign_df}
\alias{assign_df}
\title{Helper to assign degrees of freedom}
\usage{
assign_df(x, df_default = 4)
}
\arguments{
\item{x}{input matrix.}

\item{df_default}{default df to be used. Default is 4.}
}
\value{
Vector of length \code{ncol(x)} with degrees of freedom for each variable in \code{x}.
}
\description{
Determine the number of unique values in a variable. To be used in \code{\link[=mfp2]{mfp2()}}.
}
\details{
Variables with fewer than or equal to three unique values, for example,
will be assigned df = 1. df = 2 will be assigned to variables with 4-5
unique values, and df = 4 will be assigned to variables with unique values
greater than or equal to 6.
}
\examples{
x <- matrix(1:100, nrow = 10)
assign_df(x)

}
