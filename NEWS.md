# mfp2 1.0.0

* Initial CRAN submission.

# mfp2 1.0.1

* Final regression coefficients are now expressed on the original scale of the data.
* `predict.mfp2()`: Fixed a bug when using `type = "terms"`.
* `predict.mfp2()`: Added the `nseq` argument.
* `mfp2()`: Fixed an issue where the `keep` argument was inactive when `criterion = "aic"` or `criterion = "bic"`.
* Fixed a bug in BIC calculation for Cox models: `nobs` now correctly uses the number of events rather than the total number of observations.
* Improved documentation for clarity and added references regarding the effects of influential points in FP functions.
* `DESCRIPTION`: Removed `NeedsCompilation` field as the package contains no native code.

# mfp2 1.0.2
* Added spike-at-zero algorithm for modeling semicontinuous covariates.  
* `mfp2()`: Introduced `zero` and `catzero` options to model only the positive part of a semicontinuous covariate.  
* `mfp2()`: `family` argument now supports GLM family functions (e.g., `stats::binomial(link = "probit")`) and different link functions.  
* `predict.mfp2()`: Fixed bug when using `type = "response"` with the binomial family.  
* `transform_data_step()`: Improved efficiency by transforming a variable only when its fractional polynomial power or spike decision changes.
* `fracplot()`: Extended to allow plotting variables with a spike at zero. 
