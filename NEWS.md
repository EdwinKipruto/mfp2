# mfp2 1.1.0

* `mfpi()`: Added fractional polynomial interaction analysis between a categorical grouping variable and continuous covariates. Supports four flexibility levels (`flex = "flex1"` through `flex = "flex4"`) controlling how FP powers are estimated and constrained across groups. Both matrix (`mfpi.default()`) and formula (`mfpi.formula()`) interfaces are provided.
* `predict.mfpi()`: Added prediction method for `mfpi` objects, including fitted functions per group, group-difference curves, and pointwise delta-method standard errors and confidence intervals.
* `plot.mfpi()`: Added plotting method for `mfpi` objects. Supports fitted-function plots, difference plots, and combined side-by-side displays via the `patchwork` package.
* `mfp2()`: Added spike-at-zero algorithm for modelling semicontinuous covariates.
* `mfp2()`: Introduced `zero` and `catzero` options to model only the positive part of a semicontinuous covariate.
* `mfp2()`: The `family` argument now supports GLM family functions, such as `stats::binomial(link = "probit")`, and different link functions.
* `predict.mfp2()`: Fixed a bug when using `type = "response"` with the binomial family.
* `mfp2()`: Improved computational performance.
* `plot.mfp2()`: Added an S3 `plot()` method for `mfp2` objects, superseding `fracplot()`. `plot(fit)` is now the recommended way to produce partial predictor and contrast plots. `fracplot()` is retained as an alias for backward compatibility.

# mfp2 1.0.1

* Final regression coefficients are now expressed on the original scale of the data.
* `predict.mfp2()`: Fixed a bug when using `type = "terms"`.
* `predict.mfp2()`: Added the `nseq` argument.
* `mfp2()`: Fixed an issue where the `keep` argument was inactive when `criterion = "aic"` or `criterion = "bic"`.
* Fixed a bug in BIC calculation for Cox models: `nobs` now correctly uses the number of events rather than the total number of observations.
* Improved documentation for clarity and added references regarding the effects of influential points in FP functions.

# mfp2 1.0.0

* Initial CRAN submission.