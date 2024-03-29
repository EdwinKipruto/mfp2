% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_metrics.R
\name{deviance_gaussian}
\alias{deviance_gaussian}
\title{Deviance computations as used in mfp in stata}
\usage{
deviance_gaussian(rss, weights, n)
}
\arguments{
\item{rss}{residual sum of squares.}

\item{weights}{numeric vector of weights used in computation of \code{rss}.}

\item{n}{number of observations used to compute \code{rss}.}
}
\value{
A numeric value representing the deviance of a Gaussian model.
}
\description{
Deviance computations as used in mfp in stata
}
\details{
Note that this is not the usual formula of deviance used in R, but
uses the formula found here https://www.stata.com/manuals/rfp.pdf.

It can be applied for normal error models, but should not be used for other
kinds of glms.
}
