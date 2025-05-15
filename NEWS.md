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
* 