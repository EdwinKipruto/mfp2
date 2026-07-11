# `mfp2`

## Overview

`mfp2` implements multivariable fractional polynomial (MFP) models and related
extensions. It performs variable selection and functional form selection when
modeling the relationship between a covariate matrix `x` and an outcome `y`.
The package supports generalized linear models with families `"gaussian"`,
`"binomial"`, and `"poisson"`, as well as Cox proportional hazards models
(`"cox"`).

In addition to standard MFP modeling, `mfp2` supports sigmoid-shaped covariate
effects using the approximate cumulative distribution (ACD) transformation,
which cannot be represented by standard fractional polynomial functions. It also
supports semi-continuous covariates with a spike at zero through a dedicated
two-stage selection procedure.

The package also provides `mfpi()`, which investigates interactions between a
categorical grouping variable, such as treatment group, and continuous
covariates using fractional polynomial submodels. This allows treatment-by-
covariate or group-by-covariate interactions to be assessed while allowing
nonlinear covariate effects.

## Compatibility with existing software packages

`mfp2` closely emulates the functionality of the `mfp` and `mfpa` packages in
Stata. It extends the functionality of the existing `mfp` package in R by
providing:

- both matrix and formula interfaces for input data,
- sigmoid transformations via the approximate cumulative distribution (ACD)
  transformation,
- support for covariates with a spike at zero,
- fractional-polynomial interaction analyses between categorical grouping
  variables and continuous covariates,
- estimation and plotting of contrasts and partial linear predictors to
  investigate and visualize nonlinear effects, and
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
sigmoid dose-response function*. The Stata Journal.

For the spike-at-zero algorithm, see Becher et al. (2012),
*Analysing covariates with spike at zero: a modified FP procedure and conceptual
issues*. Biometrical Journal, 54(5), 686-700.

For fractional-polynomial interaction analyses between treatment and continuous
covariates, see Royston and Sauerbrei (2004), *A new approach to modelling
interactions between treatment and continuous covariates in clinical trials by
using fractional polynomial submodels*. Statistics in Medicine, 23(16),
2509-2525. doi: [10.1002/sim.1815](https://doi.org/10.1002/sim.1815)

See also Royston and Sauerbrei (2013), *Interaction of treatment with a
continuous variable: simulation study of significance level for several methods
of analysis*. Statistics in Medicine, 32(22), 3788-3803. doi:
[10.1002/sim.5813](https://doi.org/10.1002/sim.5813)

For simulation results on power, see Royston and Sauerbrei (2014), *Interaction
of treatment with a continuous variable: simulation study of power for several
methods of analysis*. Statistics in Medicine, 33(27), 4695-4708. doi:
[10.1002/sim.6308](https://doi.org/10.1002/sim.6308)



