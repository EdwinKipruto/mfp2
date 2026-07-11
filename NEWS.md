# mfp2 1.1.0

## New functionality

* Added `mfpi()` for fractional polynomial interaction analysis between a categorical grouping variable and continuous covariates. The function supports four flexibility levels, from `flex = "flex1"` to `flex = "flex4"`, which determine how fractional polynomial powers are estimated and constrained across groups. Both matrix and formula interfaces are available through `mfpi.default()` and `mfpi.formula()`, respectively.

* Added `predict.mfpi()` for generating predictions from `mfpi` objects. The method supports group-specific fitted functions, between-group difference curves, and pointwise standard errors and confidence intervals calculated using the delta method.

* Added `plot.mfpi()` for visualizing `mfpi` model results. The method supports fitted-function plots, group-difference plots, and combined side-by-side displays using the `patchwork` package.

* Added support in `mfp2()` for modelling semicontinuous covariates using a spike-at-zero approach.

* Added the `zero` and `catzero` arguments to `mfp2()` for specifying semicontinuous covariates and modelling their positive components.

* Extended the `family` argument in `mfp2()` to accept GLM family objects and alternative link functions, including specifications such as `stats::binomial(link = "probit")`.

* Added `plot.mfp2()`, an S3 `plot()` method for `mfp2` objects. The method provides partial predictor and contrast plots and is now the recommended interface for visualizing fitted `mfp2` models. The existing `fracplot()` function is retained as an alias for backward compatibility.

## Bug fixes

* Fixed an issue in `predict.mfp2()` when `type = "response"` was used with binomial-family models.

## Performance and compatibility

* Improved the computational efficiency of `mfp2()`.

* Improved support for namespace-qualified fractional polynomial terms, including `mfp2::fp()` and `mfp2::fp2()`.

# mfp2 1.0.1

## Changes and improvements

* Updated the final regression coefficients so that they are reported on the original scale of the data.

* Added the `nseq` argument to `predict.mfp2()`.

* Improved the package documentation and added references concerning the influence of influential observations on fractional polynomial functions.

## Bug fixes

* Fixed an issue in `predict.mfp2()` when `type = "terms"` was used.

* Fixed an issue in `mfp2()` where the `keep` argument was not applied when `criterion = "aic"` or `criterion = "bic"`.

* Corrected the BIC calculation for Cox models so that the effective sample size is based on the number of observed events rather than the total number of observations.

# mfp2 1.0.0

* Initial CRAN release.