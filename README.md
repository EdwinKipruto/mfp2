# `mfp2`

## Overview

`mfp2` implements multivariable fractional polynomial (MFP) models and their 
extensions. It performs variable selection and functional form selection when 
modeling the relationship between a covariate matrix `x` and an outcome `y`. 
The function currently supports generalized linear models with families 
`"gaussian"`, `"binomial"`, and `"poisson"`, as well as Cox proportional hazards 
models (`"cox"`). 

In addition to the standard MFP, it allows modeling of a sigmoid-shaped covariate 
effects using the approximate cumulative distribution (ACD) transformation, 
which cannot be achieved with standard fractional polynomial functions. 
It also supports semi-continuous covariates with a “spike at zero” through 
a dedicated two-stage selection procedure. 

## Compatibility with existing software packages

 `mfp2` closely emulates the functionality of the `mfp` and `mfpa` packages in Stata. 
 It extends the functionality of the existing `mfp` package in R by providing:

 - both matrix and formula interfaces for input data,  
 - sigmoid transformations via the approximate cumulative distribution (ACD) transformation,  
 - support for covariates with a spike at zero,  
 - estimation and plotting of contrasts and partial linear predictors to 
   investigate and visualize non-linear effects, and  
 - various computational optimizations to improve speed and usability.  

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

To learn more about the MFP algorithm, see Royston and Sauerbrei (2008), 
*Multivariable Model-Building: A Pragmatic Approach to Regression Analysis 
based on Fractional Polynomials for Modelling Continuous Variables*. 
John Wiley & Sons.

For details on the ACD transformation, see Royston (2014), 
*A smooth covariate rank transformation for use in regression models with a 
sigmoid dose–response function*. The Stata Journal.

For the spike-at-zero algorithm, see Becher et al. (2012), 
*Analysing covariates with spike at zero: a modified FP procedure and conceptual issues*. 
Biometrical Journal, 54(5), 686–700.



