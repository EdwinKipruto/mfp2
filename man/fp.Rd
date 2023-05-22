% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fp.R
\name{fp}
\alias{fp}
\title{Helper to assign attributes to a variable undergoing FP-transformation}
\usage{
fp(
  x,
  df = 4,
  alpha = 0.05,
  select = 0.05,
  shift = NULL,
  scale = NULL,
  center = TRUE,
  acd = FALSE,
  powers = NULL
)
}
\arguments{
\item{x}{a vector representing a continuous variable undergoing
fp-transformation.}

\item{df, alpha, select, shift, scale, center, acd}{See \code{\link[=mfp2]{mfp2()}}) for details.}

\item{powers}{a vector of powers to be evaluated for \code{x}. Default is \code{NULL}
and \code{powers = c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)} will be used.}
}
\value{
The vector \code{x} with new attributes relevant for fp-transformation. All
arguments passed to this function will be stored as attributes.
}
\description{
Used in formula interface to \code{mfp2()}.
}