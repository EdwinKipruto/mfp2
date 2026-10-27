## R CMD check results

0 errors | 0 warnings | 0 notes

## Release summary

This updates mfp2 from version 1.0.1 to version 1.1.0.

The main changes are:

* Added `mfpi()` for interaction analysis with a categorical grouping variable, with prediction and plotting methods for group-specific functions and differences.
* Added spike-at-zero modelling and options for modelling positive predictor values with or without a separate zero indicator.
* Expanded model support to include multinomial regression, generalized estimating equations, parametric survival regression, and Fine–Gray models.
* Expanded support for GLM families and alternative link functions.
* Added `plot.mfp2()`; `fracplot()` remains available for backward compatibility.
* Corrected prediction, confidence-interval, centring, subsetting, and factor-handling behaviour, and improved computational performance.

Full details are provided in `NEWS.md`.

## Interface changes

Variable-specific `df`, `select`, `alpha`, `shift`, and `scale` settings in the default methods now use named vectors matched to predictor columns. Partial named vectors are supported; unnamed vectors containing more than one value are rejected.

ACD settings now use `acd_vars` in the default interface and `acd = TRUE` within formula terms. The former `acdx` spelling remains available as a deprecated alias.

Formula-based survival models use `strata(...)` in the formula. The direct `strata` argument to `mfp2.formula()` remains available as a deprecated compatibility option.

## Reverse dependencies

All three identified reverse dependencies were checked using `revdepcheck`, comparing the released CRAN version of mfp2 with the local development version:

* ctlr 0.1.0
* midoc 1.0.0
* test2norm 0.3.0.1

No new errors, warnings, or notes were detected.

3 packages checked | 3 passed | 0 broken
