# mfp2 1.1.0

## New functionality

* Added Gamma, inverse-Gaussian, and negative-binomial regression. GLM families can now be supplied as names, functions, or family objects, including alternative links.

* Added ordinal regression through `ordinal_family()` and multinomial regression through `multinomial_family()`, with prediction and plotting support.

* Added parametric survival regression through `survreg_family()` and Fine–Gray competing-risks regression through `finegray_family()`.

* Added generalized estimating equations through `gee_family()`, with cluster identifiers supplied through `id` and optional measurement order through `waves`. GEE models default to QICu-based selection.

* Added `mfpi()` for interaction analysis with a categorical grouping variable, including prediction and plotting methods for group-specific functions and differences.

* Added spike-at-zero modelling and options for modelling positive predictor values with or without a separate zero indicator.

* Added support for assessing categorical predictors as complete terms, with `term_groups` available in the default interface.

* Added `force_max_fp_vars` to retain specified predictors at their maximum permitted FP degree.

* Added `plot.mfp2()`. `fracplot()` is deprecated and will be removed in a future
  release; use `plot()` with a fitted `mfp2` object instead.

* Added the `fitter` argument to choose between `stats::glm.fit()` (`"base"`)
  and the optional `fastglm` backend for GLM families.

## Interface changes

* Variable-specific `df`, `select`, `alpha`, `shift`, and `scale` settings in the default methods now use named vectors matched to predictor columns. Partial named vectors are supported; unnamed vectors containing more than one value are rejected.

* Cox models now default to `ties = "efron"` rather than `"breslow"`, and to `nocenter = c(-1, 0, 1)` rather than `NULL`. The `"exact"` tie option is no longer accepted.

* Observation weights must be finite and strictly positive; zero weights are no longer accepted.

* ACD settings now use `acd_vars` in the default interface and `acd = TRUE` within `fp()` or `fp2()`. The former `acdx` spelling remains available as a deprecated alias.

* The direct `strata` argument to `mfp2.formula()` is deprecated in favour of `strata(...)` in the formula. For formula-stratified models, prediction data must contain the original stratification variables.

## Fixes and improvements

* Fixed binomial response predictions for new observations.

* Corrected formula evaluation during subsetting and improved handling of unused factor levels. Removing levels from a factor with custom contrasts now produces an informative error.

* Improved support for namespace-qualified terms such as `mfp2::fp()` and `mfp2::fp2()`.


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
