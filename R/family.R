#' Normalize and validate a model family
#'
#' @param family Character family name, GLM family function, or GLM family object.
#' @param family_arg Character label used in error messages for function inputs.
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{family}}{A GLM family object, or character \code{"negbin"} or \code{"cox"}.}
#'   \item{\code{family_string}}{The normalized character family name.}
#' }
#'
#' @keywords internal
#' @noRd
normalize_family_argument <- function(family, family_arg = deparse(substitute(family))) {
  allowed_families <- c("gaussian", "binomial", "poisson", "negbin", "cox")
  family_arg <- paste(family_arg, collapse = " ")

  if (is.character(family)) {
    if (length(family) != 1L) {
      stop(
        sprintf(
          "! `family` must be a single character string; got %d values: %s.",
          length(family),
          paste(family, collapse = ", ")
        ),
        "\ni Supported character families are: gaussian, binomial, poisson, negbin, cox.",
        call. = FALSE
      )
    }

    if (!family %in% allowed_families) {
      stop(
        sprintf("! Invalid family: '%s'.", family),
        sprintf(
          "\ni Supported character families are: %s.",
          paste(allowed_families, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    family_obj <- if (family %in% c("cox", "negbin")) {
      family
    } else {
      switch(
        family,
        gaussian = stats::gaussian(),
        binomial = stats::binomial(),
        poisson  = stats::poisson()
      )
    }

    return(list(
      family = family_obj,
      family_string = family
    ))
  }

  if (is.function(family)) {
    family_obj <- tryCatch(
      family(),
      error = function(e) {
        stop(
          sprintf(
            "! Could not create a family object from the provided function `%s`. Error: %s",
            family_arg,
            conditionMessage(e)
          ),
          call. = FALSE
        )
      }
    )

    if (!inherits(family_obj, "family")) {
      stop(
        sprintf(
          "! The provided function `%s` did not return a valid GLM family object.",
          family_arg
        ),
        call. = FALSE
      )
    }

    family_string <- family_obj$family

    if (identical(family_string, "cox")) {
      stop(
        "! `cox` must be specified as a character string, not as a function.",
        call. = FALSE
      )
    }

    if (!family_string %in% setdiff(allowed_families, "negbin")) {
      stop(
        sprintf(
          "! Invalid family returned by `%s`: '%s'. Supported families are: %s.",
          family_arg,
          family_string,
          paste(allowed_families, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    return(list(
      family = family_obj,
      family_string = family_string
    ))
  }

  if (inherits(family, "family")) {
    family_string <- family$family

    if (identical(family_string, "cox")) {
      stop(
        "! `cox` must be specified as a character string, not as a family object.",
        call. = FALSE
      )
    }

    if (!family_string %in% setdiff(allowed_families, "negbin")) {
      stop(
        sprintf(
          "! Invalid family: '%s'. Supported families are: %s.",
          family_string,
          paste(allowed_families, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    return(list(
      family = family,
      family_string = family_string
    ))
  }

  stop(
    "! `family` must be a character string, a GLM family function, or a GLM family object.",
    call. = FALSE
  )
}

#' Extract family name for the formula interface
#'
#' Internal helper used by `mfp2.formula()` to determine the family name before
#' calling `mfp2.default()`. This is needed because the formula interface must
#' know whether Cox-specific formula terms such as `strata()` are allowed, but
#' family validation and normalization are still handled by `mfp2.default()`.
#'
#' @param family A character family name, a GLM family function, or a GLM family
#'   object.
#'
#' @return A single character string giving the family name.
#'
#' @keywords internal
#' @noRd
get_family_string_formula <- function(family, family_arg = deparse(substitute(family))) {
  normalize_family_argument(
    family,
    family_arg = family_arg
  )$family_string
}


#' Resolve a model family once for repeated internal fits
#'
#' @param family Character family name, family function, family object, or "cox".
#' @return A resolved GLM family object, or the character string "cox".
#' @keywords internal
#' @noRd
resolve_fit_model_family <- function(family) {
  normalize_family_argument(family)$family
}

#' Validate response against model family
#'
#' Validates the response object `y` after the family has been normalized.
#' This function checks only response shape and family-specific admissibility;
#' it does not modify `y`.
#'
#' @param y Response vector, matrix, factor, or `survival::Surv()` object.
#' @param family_string Normalized family name: `"gaussian"`, `"binomial"`,
#'   `"poisson"`, `"negbin"`, or `"cox"`.
#' @param nobs Expected number of observations.
#'
#' @return Invisibly returns `TRUE`.
#'
#' @keywords internal
#' @noRd
validate_family_response <- function(y, family_string, nobs) {

  if (!is.character(family_string) || length(family_string) != 1L) {
    stop("! `family_string` must be a single character string.", call. = FALSE)
  }

  if (!is.numeric(nobs) || length(nobs) != 1L ||
      anyNA(nobs) || !is.finite(nobs) || nobs < 1L) {
    stop("! `nobs` must be a positive finite scalar.", call. = FALSE)
  }

  if (family_string == "cox") {
    if (!survival::is.Surv(y)) {
      stop(
        "! For `family = 'cox'`, `y` must be a `survival::Surv()` object.",
        call. = FALSE
      )
    }

    if (NROW(y) != nobs) {
      stop(
        paste0(
          "! `y` has ", NROW(y), " rows but `x` has ", nobs,
          " rows; they must match."
        ),
        call. = FALSE
      )
    }

    type <- attr(y, "type", exact = TRUE)
    if (!identical(type, "right")) {
      stop(
        paste0(
          "! Only right-censored survival data are currently supported; ",
          "`y` has censoring type '", type, "'."
        ),
        call. = FALSE
      )
    }

    if (anyNA(y) || any(!is.finite(as.matrix(y)))) {
      stop(
        "! For `family = 'cox'`, `y` must contain only finite, non-missing values.",
        call. = FALSE
      )
    }

    return(invisible(TRUE))
  }

  if (survival::is.Surv(y)) {
    stop(
      paste0(
        "! Response is a `survival::Surv()` object but family = '",
        family_string, "'. Set `family = 'cox'`."
      ),
      call. = FALSE
    )
  }

  if (NROW(y) != nobs) {
    stop(
      paste0(
        "! `y` has ", NROW(y), " observations but `x` has ", nobs,
        " rows; they must match."
      ),
      call. = FALSE
    )
  }

  if (is.data.frame(y)) {
    stop(
      "! `y` must not be a data frame.",
      call. = FALSE
    )
  }

  if (is.matrix(y)) {
    if (family_string != "binomial") {
      stop(
        "! Matrix responses are only supported for `family = 'binomial'`.",
        "i For grouped binomial counts, use `y = cbind(successes, failures)`.",
        call. = FALSE
      )
    }

    if (!is.numeric(y) || ncol(y) != 2L) {
      stop(
        "! For `family = 'binomial'`, matrix `y` must be a numeric two-column matrix.",
        "i Use `y = cbind(successes, failures)`.",
        call. = FALSE
      )
    }

    if (anyNA(y) || any(!is.finite(y)) || any(y < 0)) {
      stop(
        "! For `family = 'binomial'`, matrix `y` counts must be finite, non-missing, and non-negative.",
        call. = FALSE
      )
    }

    if (any(rowSums(y) <= 0)) {
      stop(
        "! For `family = 'binomial'`, each row of matrix `y` must contain at least one trial.",
        call. = FALSE
      )
    }

    return(invisible(TRUE))
  }

  if (family_string == "gaussian") {
    if (!is.numeric(y)) {
      stop(
        "! For `family = 'gaussian'`, `y` must be a numeric vector.",
        call. = FALSE
      )
    }

    if (anyNA(y) || any(!is.finite(y))) {
      stop(
        "! For `family = 'gaussian'`, `y` must contain only finite, non-missing values.",
        call. = FALSE
      )
    }

    return(invisible(TRUE))
  }

  if (family_string == "negbin") {
    if (!is.numeric(y)) {
      stop(
        "! For `family = 'negbin'`, `y` must be a numeric vector of counts.",
        call. = FALSE
      )
    }

    if (anyNA(y) || any(!is.finite(y)) || any(y < 0) || any(y != floor(y))) {
      stop(
        "! For `family = 'negbin'`, `y` must contain only finite, non-missing, non-negative integer counts.",
        call. = FALSE
      )
    }

    return(invisible(TRUE))
  }

  if (family_string == "poisson") {
    if (!is.numeric(y)) {
      stop(
        "! For `family = 'poisson'`, `y` must be a numeric vector.",
        call. = FALSE
      )
    }

    if (anyNA(y) || any(!is.finite(y)) || any(y < 0)) {
      stop(
        "! For `family = 'poisson'`, `y` must contain only finite, non-missing, non-negative values.",
        call. = FALSE
      )
    }

    return(invisible(TRUE))
  }

  if (family_string == "binomial") {
    if (is.factor(y)) {
      if (nlevels(y) != 2L) {
        stop(
          "! For `family = 'binomial'`, factor `y` must have exactly two levels.",
          call. = FALSE
        )
      }

      if (anyNA(y)) {
        stop(
          "! For `family = 'binomial'`, factor `y` must not contain missing values.",
          call. = FALSE
        )
      }

      return(invisible(TRUE))
    }

    if (is.numeric(y)) {
      if (anyNA(y) || any(!is.finite(y))) {
        stop(
          "! For `family = 'binomial'`, numeric `y` must contain only finite, non-missing values.",
          call. = FALSE
        )
      }

      if (any(y < 0 | y > 1)) {
        stop(
          "! For `family = 'binomial'`, numeric `y` must contain values in [0, 1].",
          call. = FALSE
        )
      }

      return(invisible(TRUE))
    }

    stop(
      "! For `family = 'binomial'`, `y` must be a numeric vector, two-level factor, or numeric two-column matrix.",
      call. = FALSE
    )
  }

  stop(
    paste0("! Unsupported family: '", family_string, "'."),
    call. = FALSE
  )
}