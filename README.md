# `mfp2`

## Overview

`mfp2` implements multivariable fractional polynomial (MFP) models and various
extensions. It allows the selection of variables and functional forms when
modelling the relationship of a data matrix `x` and some outcome `y`. Currently
supports generalized linear models and Cox proportional hazards models.

## Compatibility with existing software packages

`mfp2` closely emulates the functionality of the [`mfp` package in
stata](https://www.stata.com/manuals13/rmfp.pdf).

It augments the functionality of the existing [`mfp` R
package](https://cran.r-project.org/web/packages/mfp/index.html) by:

-   a matrix and a formula interface for input
-   sigmoid transformations via the approximate cumulative distribution (ACD)
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

## References

To learn more about the MFP algorithm, a good place to start is the book by
Royston, P. and Sauerbrei, W., 2008. *Multivariable Model - Building: A
Pragmatic Approach to Regression Analysis based on Fractional Polynomials for
Modelling Continuous Variables.* John Wiley & Sons.
