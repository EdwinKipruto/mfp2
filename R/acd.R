#' Helper to reset acd transformation for variables with few values
#'
#' To be used in \code{fit_mfp()}.
#' This function resets the `acdx` parameter (logical vector) of variables with
#' less than 5 distinct values to `FALSE`.
#'
#' @param x a design matrix of dimension nobs x nvars where nvars is the number
#' of predictors excluding an intercept.
#' @param acdx a named logical vector of length nvars indicating which continuous
#' variables should undergo the approximate cumulative distribution (ACD)
#' transformation. May be ordered differently than the columns of `x`.
#'
#' @return
#' Logical vector of same length as `acdx`.
#' @keywords internal
#' @noRd
reset_acd <- function(x, acdx) {

  if (!is.logical(acdx) || is.null(names(acdx))) {
    stop("`acdx` must be a named logical vector.", call. = FALSE)
  }

  if (anyNA(acdx)) {
    stop("`acdx` must not contain missing values.", call. = FALSE)
  }

  # exit early if all acdx values are FALSE
  if (!any(acdx)) {
    return(acdx)
  }

  names_acd <- names(acdx)[acdx]

  # number of unique values of each column in acdx
  n_unique <- apply(
    x[, names_acd, drop = FALSE],
    2,
    function(col) length(unique(col))
  )

  ind_reset <- which(n_unique < 5)

  if (length(ind_reset) > 0L) {
    vars_reset <- names_acd[ind_reset]
    acdx[vars_reset] <- FALSE

    warning(
      "i For any variable with fewer than 5 unique values no acd transformation can be performed.\n",
      sprintf(
        "i The requested acd transform has been reset to FALSE for the following variables: %s.",
        paste0(vars_reset, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  acdx
}


#' Fit an Approximate Cumulative Distribution Transformation
#'
#' Estimates the approximate cumulative distribution (ACD) transformation
#' described by Royston (2014). The transformation maps the observed
#' distribution of a continuous covariate approximately onto the interval
#' `(0, 1)`.
#'
#' @details
#' The observed values are first shifted and scaled according to `shift` and
#' `scale`. The function then calculates the normal-score rank transformation
#'
#' \deqn{
#' z_i = \Phi^{-1}\left\{\frac{\operatorname{rank}(x_i)-0.5}{n}\right\},
#' }
#'
#' where tied observations receive average ranks.
#'
#' An FP1 model is fitted for each candidate power:
#'
#' \deqn{
#' E(z_i) = \beta_0 + \beta_1 x_i^p.
#' }
#'
#' The power giving the best fit is selected, and the resulting ACD
#' transformation is:
#'
#' \deqn{
#' \operatorname{ACD}(x_i) =
#' \Phi\left(\hat{\beta}_0 + \hat{\beta}_1 x_i^p\right).
#' }
#'
#' The transformed values lie between zero and one and can be used to represent
#' a smooth sigmoid-shaped covariate effect in a regression model.
#'
#' Most users do not need to call `fit_acd()` directly because ACD modelling
#' can be requested through [mfp2()]. This function is useful for inspecting,
#' reproducing, or applying an ACD transformation in a custom workflow.
#'
#' If `shift = NULL`, an appropriate shift is estimated automatically using
#' \code{find_shift_factor()}. If `scale = NULL`, an appropriate scale is estimated
#' automatically using \code{find_scale_factor()}.
#'
#' When `zero = TRUE`, only positive values are used in the continuous
#' transformation. Nonpositive values are assigned a transformed value of zero,
#' and no shift is applied.
#'
#' @param x Numeric vector containing the covariate values. Missing values are
#'   not allowed.
#'
#' @param powers Optional numeric vector containing the candidate FP1 powers.
#'   If `NULL`, the standard set
#'   `c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)` is used. Power zero represents
#'   the logarithmic transformation.
#'
#' @param shift Numeric scalar used to shift `x` before transformation.
#'   The default, `0`, applies no shift. Use `NULL` to estimate a suitable
#'   shift automatically.
#'
#' @param scale Positive numeric scalar used to scale the shifted values.
#'   The default, `1`, applies no scaling. Use `NULL` to estimate a suitable
#'   scale automatically.
#'
#' @param zero Logical scalar. If `TRUE`, the ACD transformation is fitted to
#'   the positive values and nonpositive values are assigned zero. If `FALSE`,
#'   all observations are included after any requested shifting and scaling.
#'
#' @return
#' A list containing:
#'
#' \describe{
#'   \item{\code{acd}}{
#'     Numeric vector containing the fitted ACD-transformed values.
#'   }
#'
#'   \item{\code{beta0}}{
#'     Estimated intercept of the selected FP1 approximation.
#'   }
#'
#'   \item{\code{beta1}}{
#'     Estimated coefficient of the selected FP1 approximation.
#'   }
#'
#'   \item{\code{power}}{
#'     Selected FP1 power.
#'   }
#'
#'   \item{\code{shift}}{
#'     Shift applied before fitting the transformation.
#'   }
#'
#'   \item{\code{scale}}{
#'     Scale applied before fitting the transformation.
#'   }
#' }
#'
#' @references
#' Royston, P. (2014). A smooth covariate rank transformation for use in
#' regression models with a sigmoid dose-response function.
#' \emph{The Stata Journal}, 14(2), 329--341.
#'
#' Royston, P. and Sauerbrei, W. (2016). mfpa: Extension of mfp using the ACD
#' covariate transformation for enhanced parametric multivariable modeling.
#' \emph{The Stata Journal}, 16(1), 72--87.
#'
#' @examples
#' set.seed(42)
#'
#' # Positive-valued covariate requiring no shift.
#' x <- rgamma(100, shape = 2, rate = 0.5)
#' acd_fit <- fit_acd(x)
#'
#' acd_fit$power
#' head(acd_fit$acd)
#'
#' # Estimate the shift and scale automatically.
#' x_unscaled <- rnorm(100, mean = 10, sd = 20)
#' acd_auto <- fit_acd(
#'   x_unscaled,
#'   shift = NULL,
#'   scale = NULL
#' )
#'
#' \dontrun{
#' # ACD can ordinarily be requested directly during MFP fitting.
#' fit <- mfp2(
#'   outcome ~ fp(exposure, acd = TRUE) + fp(age),
#'   data = analysis_data,
#'   verbose = FALSE
#' )
#' }
#'
#' @keywords internal
#' @noRd
fit_acd <- function(x, powers = NULL, shift = 0, scale = 1, zero = FALSE,
                    fitter = c("base", "fastglm")) {
  fitter <- match.arg(fitter)

  # --- Input checks ---
  if (!is.numeric(x)) {
    stop("`x` must be a numeric vector.")
  }
  if (!is.vector(x)) {
    stop("`x` must be a vector (not a matrix, data.frame, or list).")
  }
  if (length(x) < 2) {
    stop("`x` must contain at least two values.")
  }
  if (anyNA(x)) {
    stop("`x` contains missing values. Please remove or impute them before calling fit_acd().")
  }

  if (!is.null(powers)) {
    if (!is.numeric(powers) || any(!is.finite(powers))) {
      stop("`powers` must be a numeric vector of finite values or NULL.")
    }
  } else {
    powers <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
  }

  if (!is.null(shift) && (!is.numeric(shift) || length(shift) != 1 || is.na(shift))) {
    stop("`shift` must be a single numeric value or NULL.")
  }

  if (!is.null(scale) && (!is.numeric(scale) || length(scale) != 1 || scale <= 0 || is.na(scale))) {
    stop("`scale` must be a single positive numeric value or NULL.")
  }

  if (!is.logical(zero) || length(zero) != 1 || is.na(zero)) {
    stop("`zero` must be a single logical value (TRUE or FALSE).")
  }

  # --- Preprocessing ---
  if (is.null(shift)) {
    shift <- find_shift_factor(x)
  }

  if (zero) {
    x[x <= 0] <- 0
    # zero transformation overrides shift
    shift <- 0
  }

  if (is.null(scale)) {
    scale <- find_scale_factor(x)
  }

  x <- (x + shift) / scale

  # check whether acd is estimable
  #if (!all(x > 0) && !zero) {
  #  warning("All values of `x` must be positive after shifting. ",
  #       "Try specifying a larger `shift` or use `shift = NULL` to estimate it automatically.")
  #}
  # --- ACD transformation ---
  # see here for details: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5663339/pdf/emss-59479.pdf)
  n <- length(x)
  z <- stats::qnorm((rank(x, ties.method = "average") - 0.5) / n)

  # estimate the best p in model E(z) = beta0 + beta1*x^p using the data
  fit <- find_best_fp1_for_acd(
    y = z, x = x, powers = powers, zero = zero, fitter = fitter
  )

  coefx <- fit$fit$coefficients
  zhat <- fit$fit$fitted.values

  list(acd = stats::pnorm(zhat),
       beta0 = coefx[1],
       beta1 = coefx[2],
       power = fit$power,
       shift = shift,
       scale = scale)
}

#' Function to apply Approximate Cumulative Distribution (ACD)
#'
#' Applies the acd transformation as outlined in Royston (2014) and Royston and
#' Sauerbrei (2016).
#' Designed to work with the output of \code{fit_acd()}, Please refer to the corresponding
#' documentation for more details.
#'
#' @param x a numeric vector.
#' @param beta0,beta1 each a numeric value, representing the coefficients of
#' the FP1 model for the ACD transformation.
#' @param power a numeric value, estimated power to be used in the FP1 model for
#' the ACD transformation.
#' @param shift a numeric value that is used to shift the values of `x` to
#' positive values.
#' @param scale a numeric value used to scale `x`.
#' @param zero Logical indicating whether only positive values of the variable
#' should be transformed, with nonpositive values (zero or negative) set to zero.
#' If \code{TRUE}, transformation is applied only to positive values; nonpositive values
#' are replaced with zero before transformation.
#' @param ... not used.
#'
#' @return
#' The transformed input vector `x`.
#' @keywords internal
#' @noRd
apply_acd <- function(x, beta0, beta1, power, shift, scale, zero, ...) {

  if (length(power) != 1) {
    stop("! `power` must be a single numeric value.")
  }

  x_power <- transform_vector_fp(
    x            = x,
    power        = power,
    shift        = shift,
    scale        = scale,
    zero         = zero,
    check_binary = FALSE
  )

  zhat <- beta0 + beta1 * x_power[, 1L]

  stats::pnorm(zhat)
}

#' Function to fit univariable FP1 models for acd transformation
#'
#' To be used in \code{fit_acd()}.
#'
#' @inheritParams fit_acd
#' @param y normal cdf of rank transform of `x`.
#'
#' @return
#' The best FP power with smallest deviance and the fitted model.
#' @keywords internal
#' @noRd
find_best_fp1_for_acd <- function(x,
                                  y,
                                  powers,
                                  zero,
                                  fitter = c("base", "fastglm")) {
  fitter <- match.arg(fitter)

  if (!is.null(dim(x))) {
    stop("! `x` must be a vector.")
  }

  # Generate all possible FP1 transformations.
  trafo <- generate_transformations_fp(
    x = x,
    degree = 1L,
    powers = powers,
    zero = zero
  )$data

  # Fit linear Gaussian models for each FP1 function.
  n_powers <- length(powers)
  family_gaussian <- stats::gaussian()
  family_string_gaussian <- family_gaussian$family

  # Reuse one intercept-augmented design matrix across candidate fits.
  #
  # Only the intercept name is required by fit_glm(..., x_has_intercept = TRUE).
  # The FP column does not need candidate-specific names because the selected
  # power is returned separately. A stable "fp1" name is enough for readable
  # coefficient names and avoids relying on colnames(trafo[[i]]), which may be
  # NULL.
  first_trafo <- trafo[[1L]]

  design_mat <- cbind(
    "(Intercept)" = rep.int(1, NROW(first_trafo)),
    "fp1" = first_trafo[, 1L]
  )

  fp_col <- 2L

  # Store deviance and model object for each candidate power.
  devs <- vector("numeric", n_powers)
  fits <- vector("list", n_powers)

  for (i in seq_len(n_powers)) {
    data_xi <- trafo[[i]]

    # Replace only the numeric contents of the FP column.
    # Do not update the column name: "fp1" is intentionally stable.
    design_mat[, fp_col] <- data_xi[, 1L]

    fit <- fit_model(
      x = design_mat,
      y = y,
      family = family_gaussian,
      family_string = family_string_gaussian,
      fitter = fitter,
      x_has_intercept = TRUE
    )

    devs[i] <- -2 * fit$logl
    fits[[i]] <- fit$fit
  }

  index <- which.min(devs)

  list(
    power = powers[[index]],
    fit = fits[[index]]
  )
}
