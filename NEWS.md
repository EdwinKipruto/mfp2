# mfp2 1.1.0

## New functionality

* Expanded the `family` interface to all likelihood-based base `glm()`
  families and their supported links. Quasi families remain unsupported because
  MFP selection requires a likelihood. The existing Cox path is unchanged, and
  the package now adds `survreg_family()` for parametric survival models and
  `finegray_family()` for Fine--Gray subdistribution-hazards models. Repeated
  survival candidate fits use matrix-level hot paths; final fits retain native
  `survreg` or `coxph` objects.

* Added `mfpi()` for fractional polynomial interaction analysis between a categorical grouping variable and continuous covariates. The function supports four flexibility levels, from `flex = "flex1"` to `flex = "flex4"`, which determine how fractional polynomial powers are estimated and constrained across groups. Both matrix and formula interfaces are available through `mfpi.default()` and `mfpi.formula()`, respectively.

* Added `predict.mfpi()` for generating predictions from `mfpi` objects. The method supports group-specific fitted functions, between-group difference curves, and pointwise standard errors and confidence intervals calculated using the delta method.

* Added `plot.mfpi()` for visualizing `mfpi` model results. The method supports fitted-function plots, group-difference plots, and combined side-by-side displays using the `patchwork` package.

* Added support in `mfp2()` for modelling semicontinuous covariates using a spike-at-zero approach.

* Added `prop_zero` to retained spike-at-zero term metadata and to the SAZ sections printed by `mfp2()` and `mfpi()`. The value is the proportion of finite observations in the fitted sample that belong to the structural-zero component.


* Added the `zero` and `catzero` arguments to `mfp2()` for specifying semicontinuous covariates and modelling their positive components.

* Extended the `family` argument in `mfp2()` to accept GLM family objects and alternative link functions, including specifications such as `stats::binomial(link = "probit")`.

* Added `plot.mfp2()`, an S3 `plot()` method for `mfp2` objects. The method provides partial predictor and contrast plots and is now the recommended interface for visualizing fitted `mfp2` models. The existing `fracplot()` function is retained as an alias for backward compatibility.

* Changed per-variable `shift` and `scale` settings in `mfp2.default()` from positional, column-index-based vectors to named vectors matched to `colnames(x)`. Named vectors may specify only a subset of predictors, with unspecified values estimated automatically; their order is irrelevant. Unnamed vectors with more than one value are now rejected to prevent settings from being assigned to the wrong columns.

* Changed per-variable `df`, `select`, and `alpha` settings in `mfp2.default()` and `mfpi.default()` from positional vectors to named overrides matched to `colnames(x)`. Named vectors may specify only a subset of predictors; omitted predictors use the package defaults (`df = 4`, `select = 0.05`, and `alpha = 0.05`). Unnamed vectors with more than one value are now rejected.


## Subsetting and factor handling

* Updated `mfp2.formula()` so that `subset` expressions are evaluated exactly once using standard formula semantics: names are resolved from `data` first and then from the formula environment.

* Formula calls now use the same lookup rule for observation-level `weights`,
  `offset`, and Fine--Gray `id` expressions. The deprecated
  `mfp2.formula()` `strata` compatibility argument follows the same convention
  during its transition period. Complete vectors are still validated before
  subsetting.

* Formula methods now construct the complete model frame once and retain requested rows directly. This avoids non-standard-evaluation failures involving internal row-index variables and prevents formula expressions, offsets, and strata terms from being evaluated a second time.

* For formula fits with a non-`NULL` subset, unused factor levels are dropped and default unordered or ordered contrasts are regenerated from the retained levels. Factor-specific custom contrasts are preserved when all levels remain; if subsetting removes a level from a factor with an attached custom contrast, fitting stops with a targeted error rather than silently changing the coding.

## Bug fixes

* Formula-stratified `mfp2` and `mfpi` models now require prediction strata to
  be supplied through the original `strata(...)` variables in `newdata`.
  Supplying a separate `strata` argument is rejected with targeted guidance,
  preventing raw values such as `0` and `1` from being confused with fitted
  formula labels such as `pf=0` and `pf=1`. Matrix-interface prediction retains
  its existing `strata` argument.

* Corrected centering for zero-handled FP and ACD bases. Their uncentered basis
  remains zero at exact-zero observations, after which the complete
  fitted-sample column mean is subtracted from every row. Centering therefore
  no longer changes fitted values or creates an implicit zero effect when no
  zero indicator was requested. Plot documentation and titles now distinguish
  this case from `catzero` and spike-at-zero models.

* Fixed an issue in `predict.mfp2()` when `type = "response"` was used with binomial-family models.

* Corrected the negative-binomial null deviance for models with offsets. The
  reporting-only null fit now follows `MASS::glm.nb()` by retaining the fitted
  theta and fitting an intercept plus offset; repeated candidate fits remain on
  the existing `fastglm` hot path.

* Added the native `type = "linear"` alias for complete-model `survreg`
  predictions. It is equivalent to `type = "link"` and `type = "lp"`.

## Performance and compatibility

* Aligned survival stratification with formula conventions. Formula calls now
  use `strata(...)`: the direct `strata` argument to `mfp2.formula()` remains as
  a deprecated compatibility path, while `mfpi.formula()` accepts strata only
  through the formula. The `strata` argument remains unchanged in both matrix
  interfaces.

* Renamed the public ACD options for clarity: matrix interfaces now use
  `acd_vars`, while formula terms use `fp(..., acd = TRUE)` or
  `fp2(..., acd = TRUE)`. The former `acdx` spelling remains available as a
  deprecated compatibility alias in both interfaces.

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
