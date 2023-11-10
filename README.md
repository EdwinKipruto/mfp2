# `mfp2`

## Overview

`mfp2` implements multivariable fractional polynomial (MFP) models and various
extensions. It allows the selection of variables and functional forms when
modelling the relationship of a data matrix `x` and some outcome `y`. Currently, it
supports generalized linear models and Cox proportional hazards models.
Additionally, it has the ability to model a sigmoid relationship between covariate `x` and an outcome variable `y`
using approximate cumulative distribution (ACD) transformation- a feature that a standard fractional polynomial function cannot achieve. 

## Compatibility with existing software packages

`mfp2` closely emulates the functionality of the `mfp` and `mfpa` package in Stata.

It augments the functionality of the existing `mfp` R
package by:

-   a matrix and a formula interface for input
-   sigmoid transformations via the ACD transformation
-   estimation and plotting of contrasts and partial linear predictors to
    investigate and visualize non-linear effects
-   various optimizations to increase speed and user friendliness

## Installation

``` r
# Install the development version from GitHub
# install.packages("pak")
pak::pak("EdwinKipruto/mfp2")

# or 
# install.packages("remotes")
remotes::install_github("EdwinKipruto/mfp2")
```
<!-- badges: start -->
[![R-CMD-check](https://github.com/EdwinKipruto/mfp2/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EdwinKipruto/mfp2/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## References

To learn more about the MFP algorithm, a good place to start is the book by
Royston, P. and Sauerbrei, W., 2008. *Multivariable Model - Building: A
Pragmatic Approach to Regression Analysis based on Fractional Polynomials for
Modelling Continuous Variables.* John Wiley & Sons.

For insights into the ACD transformation, please refer to Royston (2014). *A smooth covariate rank transformation for use in regression
models with a sigmoid dose–response function.* The Stata Journal
