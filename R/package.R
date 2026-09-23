#' mfp2: Multivariable Fractional Polynomial Models and Extensions
#'
#' The `mfp2` package fits multivariable fractional polynomial models with
#' simultaneous selection of variables and functional forms for continuous
#' predictors.
#'
#' @section Main functions:
#' - [mfp2()] fits multivariable fractional polynomial models.
#' - [fp()] specifies fractional-polynomial and related modelling options
#'   within a formula.
#' - [mfpi()] investigates interactions between a categorical grouping
#'   variable and continuous, binary, or categorical variables.
#'
#' Fitted models can be examined using `print()`, `summary()`, and `coef()`,
#' used for prediction with `predict()`, and visualised with `plot()`.
#'
#' @section Supported models:
#' [mfp2()] and [mfpi()] support every likelihood-based family supplied by base
#' [stats::glm()]: Gaussian, binomial, Poisson, Gamma, and inverse-Gaussian
#' regression, including the links accepted by their family objects. Quasi
#' families are not supported because MFP selection requires a likelihood.
#' Optional negative-binomial fitting requires the `fastglm` package.
#' Multinomial logistic models use [multinomial_family()] and the `nnet`
#' package. Their selected FP transformations are common across logits and
#' their regression coefficients are logit specific.
#' Proportional-odds ordinal models use [ordinal_family()] and the optional
#' `rms` package. They estimate a common predictor effect across response
#' cut-points.
#'
#' Survival models are specified as separate families: `family = "cox"` for
#' Cox proportional hazards, [survreg_family()] for the parametric distributions
#' implemented by [survival::survreg()], and [finegray_family()] for Fine--Gray
#' proportional subdistribution hazards.
#'
#' @section Modelling extensions:
#' In addition to standard fractional-polynomial modelling, `mfp2` provides:
#'
#' - approximate cumulative distribution (ACD) transformations for selected
#'   continuous predictors;
#' - spike-at-zero (SAZ) modelling for semi-continuous predictors that contain
#'   a distinct zero component and a positive continuous component;
#' - interaction analysis between continuous, binary, or categorical predictors
#'   and categorical groups through [mfpi()].
#'
#' Spike-at-zero modelling can be requested in the formula interface with
#' `fp(x, spike = TRUE)` or in the matrix interface with
#' `spike_vars = "x"`. The procedure assesses whether the final model should
#' retain both the zero-component indicator and the positive continuous
#' component, or only one of them.
#' Covariates requested through `zero`, `catzero`, or `spike` must be
#' nonnegative: `x = 0` defines the zero component and `x > 0` defines the
#' positive component. Negative values are rejected and must be explicitly
#' recoded if they should represent the zero group.
#'
#' @section Getting started:
#' See `vignette("mfp2_introduction", package = "mfp2")` for an introduction
#' to model fitting, interpretation, prediction, and plotting.
#'
#' Further guidance is available in:
#'
#' - `vignette("MFP_Introduction", package = "mfp2")` for MFP methodology;
#' - `vignette("mfp2_ACD", package = "mfp2")` for ACD transformations;
#' - the [mfp2()] help page for the complete spike-at-zero algorithm;
#' - `vignette("mfpi", package = "mfp2")` for interaction analysis.
#'
#' @seealso
#' [mfp2()], [fp()], [mfpi()], [multinomial_family()], [ordinal_family()],
#' [survreg_family()], [finegray_family()], [predict.mfp2()], [plot.mfp2()]
#'
#' @useDynLib mfp2, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @keywords package
"_PACKAGE"
