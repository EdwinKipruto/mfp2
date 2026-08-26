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
#'   variable and continuous covariates.
#'
#' Fitted models can be examined using `print()`, `summary()`, and `coef()`,
#' used for prediction with `predict()`, and visualised with `plot()`.
#'
#' @section Supported models:
#' [mfp2()] supports Gaussian, binomial, Poisson, and optional
#' negative-binomial regression, together with Cox proportional hazards models
#' for right-censored survival outcomes. Negative-binomial fitting requires the
#' optional `fastglm` package.
#'
#' @section Modelling extensions:
#' In addition to standard fractional-polynomial modelling, `mfp2` provides:
#'
#' - approximate cumulative distribution (ACD) transformations for selected
#'   continuous predictors;
#' - spike-at-zero (SAZ) modelling for semi-continuous predictors that contain
#'   a distinct zero component and a positive continuous component;
#' - interaction analysis between continuous predictors and categorical groups
#'   through [mfpi()].
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
#' - `vignette("mfp2_spike", package = "mfp2")` for spike-at-zero modelling;
#' - `vignette("mfpi", package = "mfp2")` for interaction analysis.
#'
#' @seealso
#' [mfp2()], [fp()], [mfpi()], [predict.mfp2()], [plot.mfp2()]
#'
#' @useDynLib mfp2, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @keywords package
"_PACKAGE"
