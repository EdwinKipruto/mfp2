#' Validate a logical scalar or vector argument
#'
#' Checks that an argument is logical, contains no missing values, and has one
#' of the allowed lengths.
#'
#' @param arg Object to validate.
#' @param name Character scalar giving the argument name used in error messages.
#' @param allowed_lengths Integer vector of permitted lengths.
#' @param allow_null Logical; whether \code{NULL} is permitted. Default
#'   \code{FALSE}.
#' @param hint Optional character scalar appended (on its own line) to any
#'   error message, e.g. to point the user to an alternative, per-variable way
#'   of supplying the same option. Default \code{NULL} (no hint added).
#'
#' @return Invisibly returns \code{TRUE}.
#'
#' @keywords internal
#' @noRd
validate_logical_vector <- function(arg, name, allowed_lengths,
                                    allow_null = FALSE, hint = NULL) {
  # Appended to every error message below; empty string when no hint is given,
  # so it can be passed unconditionally without altering messages that have
  # none.
  hint_line <- if (!is.null(hint)) paste0("\n", hint) else ""

  if (is.null(arg)) {
    if (allow_null) return(invisible(TRUE))

    stop(
      sprintf("! `%s` must not be NULL.", name),
      hint_line,
      call. = FALSE
    )
  }

  if (!is.logical(arg)) {
    stop(
      sprintf("! `%s` must be logical.", name),
      sprintf("i Current type is: %s.", typeof(arg)),
      hint_line,
      call. = FALSE
    )
  }

  if (anyNA(arg)) {
    stop(
      sprintf("! `%s` must not contain NA values.", name),
      hint_line,
      call. = FALSE
    )
  }

  if (!length(arg) %in% allowed_lengths) {
    stop(
      sprintf(
        "! `%s` must have length %s; got length %d.",
        name,
        paste(allowed_lengths, collapse = " or "),
        length(arg)
      ),
      hint_line,
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Resolve the supported Cox tie-handling method
#'
#' MFP and MFPI selection fit many Cox candidate models through the low-level
#' [survival::coxph.fit()] path. That fitter does not implement the exact
#' partial likelihood, so accepting `ties = "exact"` would make selection use
#' a different tie calculation from an exact final Cox fit. Reject it here
#' rather than silently changing the likelihood used during selection.
#'
#' @param ties Character vector supplied to the public `ties` argument.
#'
#' @return One of `"breslow"` or `"efron"`.
#'
#' @keywords internal
#' @noRd
resolve_mfp_ties <- function(ties = c("breslow", "efron")) {
  # Exact Cox ties are intentionally unsupported for MFP-based selection.
  if (identical(ties, "exact")) {
    stop(
      "'ties = \"exact\"' is not supported for MFP selection; ",
      "use 'efron' or 'breslow'.",
      call. = FALSE
    )
  }

  match.arg(ties)
}

#' Validate Observation Weights
#'
#' Applies the common observation-weight contract used by the public MFP and
#' MFPI interfaces. All supplied weights must be numeric, finite, non-missing,
#' strictly positive, and aligned with the original observations. Zero weights
#' are deliberately not supported because likelihood-based MFP selection can
#' become undefined for some model families even when the underlying fitter
#' accepts zero prior weights.
#'
#' @param weights `NULL` or the user-supplied observation weights.
#' @param nobs Number of observations before `subset` is applied.
#'
#' @return Invisibly returns `TRUE`.
#'
#' @keywords internal
#' @noRd
validate_model_weights <- function(weights, nobs) {
  if (is.null(weights)) {
    return(invisible(TRUE))
  }

  if (!is.numeric(nobs) || length(nobs) != 1L || is.na(nobs) ||
      !is.finite(nobs) || nobs < 0 || nobs != floor(nobs)) {
    stop(
      "Internal error: `nobs` must be a non-negative finite integer.",
      call. = FALSE
    )
  }
  nobs <- as.integer(nobs)

  if (!is.numeric(weights)) {
    stop(
      "! `weights` must be numeric.",
      sprintf("i Current type is: %s.", typeof(weights)),
      call. = FALSE
    )
  }

  if (length(weights) != nobs) {
    stop(
      sprintf(
        "! `weights` must have one value per observation.\ni `weights` has length %d, but %d observations were expected.",
        length(weights), nobs
      ),
      call. = FALSE
    )
  }

  # Diagnose missing values before the more general finite-value check so the
  # user sees the actual input defect rather than a generic non-finite error.
  if (anyNA(weights)) {
    stop("! `weights` must not contain missing values.", call. = FALSE)
  }

  if (any(!is.finite(weights))) {
    stop("! `weights` must contain only finite values.", call. = FALSE)
  }

  # Use one strict rule for every family. In particular, Gaussian glm() can
  # estimate coefficients with zero prior weights, but its likelihood/AIC can
  # become infinite, which makes likelihood-ratio MFP selection undefined.
  # Rejecting zero weights at the public boundary avoids family-specific
  # special cases and keeps all candidate-model comparisons well-defined.
  if (any(weights <= 0)) {
    stop(
      "! `weights` must be strictly positive; zero and negative weights are not supported.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Validate Predictor Column Names
#'
#' Applies the common predictor-name contract before variable-specific settings
#' are matched or internally generated formulas are constructed. Predictor
#' names must be present, unique, non-missing, non-empty, and must not contain
#' the backtick character. Other non-syntactic names, including spaces and
#' hyphens, remain valid and are handled by the existing formula-quoting code.
#'
#' @param nms Character vector of predictor column names, or `NULL`.
#' @param object Character scalar used to identify the object in error messages.
#'
#' @return Invisibly returns `TRUE`.
#'
#' @keywords internal
#' @noRd
validate_predictor_names <- function(nms, object = "`x`") {
  if (is.null(nms)) {
    stop(
      sprintf("! %s must have column names.", object),
      call. = FALSE
    )
  }

  if (anyNA(nms)) {
    stop(
      sprintf("! Predictor names in %s must not be missing.", object),
      call. = FALSE
    )
  }

  if (any(!nzchar(nms))) {
    stop(
      sprintf("! Predictor names in %s must not be empty.", object),
      call. = FALSE
    )
  }

  if (anyDuplicated(nms)) {
    duplicated_names <- unique(nms[duplicated(nms) | duplicated(nms, fromLast = TRUE)])
    stop(
      sprintf(
        "! Predictor names in %s must be unique. Duplicated name(s): %s.",
        object,
        paste(duplicated_names, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  has_backtick <- grepl("`", nms, fixed = TRUE)
  if (any(has_backtick)) {
    bad_names <- unique(nms[has_backtick])
    stop(
      sprintf(
        "! Predictor names in %s must not contain backticks (`). Rename: %s.",
        object,
        paste(bad_names, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Normalize Cox Stratification Input
#'
#' Converts user-facing Cox stratification input to one factor with one value
#' per observation before it reaches the low-level Cox fitting code. A matrix
#' or data frame is interpreted as multiple stratification variables and its
#' columns are combined with [survival::strata()], matching the high-level
#' treatment used by [survival::coxph()].
#'
#' @param strata Stratification input supplied by a public fitting interface or
#'   reconstructed from a formula `strata()` term.
#' @param nobs Expected number of observations.
#'
#' @return `NULL` or a factor of length `nobs`.
#'
#' @keywords internal
#' @noRd
normalize_cox_strata <- function(strata, nobs) {
  if (is.null(strata)) {
    return(NULL)
  }

  if (!is.numeric(nobs) || length(nobs) != 1L || is.na(nobs) ||
      !is.finite(nobs) || nobs < 0 || nobs != floor(nobs)) {
    stop("Internal error: `nobs` must be a non-negative finite integer.", call. = FALSE)
  }
  nobs <- as.integer(nobs)

  # Higher-dimensional arrays have no unambiguous row-wise stratification
  # semantics. Two-dimensional matrices/data frames are handled explicitly
  # below as one stratification variable per column.
  if (is.array(strata) && !is.matrix(strata)) {
    stop(
      "! `strata` must be a vector, factor, matrix, or data frame.",
      call. = FALSE
    )
  }

  if (is.matrix(strata) || is.data.frame(strata)) {
    if (NROW(strata) != nobs) {
      stop(
        sprintf(
          "! `strata` must have one row per observation.\ni `strata` has %d rows, but %d observations were expected.",
          NROW(strata), nobs
        ),
        call. = FALSE
      )
    }

    if (NCOL(strata) < 1L) {
      stop("! `strata` must contain at least one stratification variable.", call. = FALSE)
    }

    strata_df <- as.data.frame(strata, stringsAsFactors = FALSE)
    bad_col <- vapply(
      strata_df,
      function(z) !(is.atomic(z) || is.factor(z)),
      logical(1L)
    )
    if (any(bad_col)) {
      stop(
        "! Each column of `strata` must be an atomic vector or factor.",
        call. = FALSE
      )
    }

    # Missing strata are a distinct input error and should be reported before
    # the broader numeric-finiteness check below (because is.finite(NA) is
    # FALSE). This keeps NA diagnostics specific while still rejecting Inf.
    if (anyNA(strata_df)) {
      stop("! `strata` must not contain missing values.", call. = FALSE)
    }

    # Numeric strata are labels, but infinite numeric labels are invalid. Keep
    # this check in the shared normalizer so fitting and prediction enforce the
    # same contract. NA values have already been handled explicitly above.
    bad_numeric <- vapply(
      strata_df,
      function(z) is.numeric(z) && any(!is.finite(z)),
      logical(1L)
    )
    if (any(bad_numeric)) {
      stop("! Numeric `strata` values must be finite.", call. = FALSE)
    }

    # survival::coxph() combines multiple formula strata into a single factor
    # before converting it to the integer codes required by coxph.fit().
    out <- do.call(
      survival::strata,
      c(as.list(strata_df), list(shortlabel = TRUE))
    )
  } else {
    if (!(is.atomic(strata) || is.factor(strata))) {
      stop(
        "! `strata` must be a vector, factor, matrix, or data frame.",
        call. = FALSE
      )
    }

    if (length(strata) != nobs) {
      stop(
        sprintf(
          "! `strata` must have one value per observation.\ni `strata` has length %d, but %d observations were expected.",
          length(strata), nobs
        ),
        call. = FALSE
      )
    }

    # Check missing values first so NA is reported as missing rather than as a
    # generic non-finite numeric value. Infinite numeric labels remain invalid.
    if (anyNA(strata)) {
      stop("! `strata` must not contain missing values.", call. = FALSE)
    }

    if (is.numeric(strata) && any(!is.finite(strata))) {
      stop("! Numeric `strata` values must be finite.", call. = FALSE)
    }

    # Treat the supplied values as categorical labels irrespective of their
    # storage mode. This is essential for character strata, which cannot be
    # safely coerced with as.integer() inside the low-level fitter.
    out <- factor(strata)
  }

  if (anyNA(out)) {
    stop("! `strata` must not contain missing values.", call. = FALSE)
  }

  # Subsetting can leave unused factor levels. They have no role in the Cox
  # likelihood and are dropped so stored strata metadata reflects fitted data.
  droplevels(out)
}


#' Validate a probability scalar or vector argument
#'
#' Checks that an argument is numeric, finite, non-missing, has length one or
#' the number of variables, and lies in the closed interval \code{[0, 1]}.
#'
#' @param arg Object to validate.
#' @param name Character scalar giving the argument name used in error messages.
#' @param nvars Number of variables expected when a vector is supplied.
#' @param hint Optional character scalar appended (on its own line) to any
#'   error message. Default \code{NULL} (no hint added).
#'
#' @return Invisibly returns \code{TRUE}.
#'
#' @keywords internal
#' @noRd
validate_probability_vector <- function(arg, name, nvars, hint = NULL) {
  hint_line <- if (!is.null(hint)) paste0("\n", hint) else ""

  if (!is.numeric(arg)) {
    stop(
      sprintf("! `%s` must be numeric.", name),
      sprintf("i Current type is: %s.", typeof(arg)),
      hint_line,
      call. = FALSE
    )
  }

  if (anyNA(arg)) {
    stop(
      sprintf("! `%s` must not contain NA values.", name),
      hint_line,
      call. = FALSE
    )
  }

  if (any(!is.finite(arg))) {
    stop(
      sprintf("! `%s` must contain only finite values.", name),
      hint_line,
      call. = FALSE
    )
  }

  if (!length(arg) %in% c(1L, nvars)) {
    stop(
      sprintf(
        "! `%s` must be a single number or a numeric vector of length %d; got length %d.",
        name, nvars, length(arg)
      ),
      hint_line,
      call. = FALSE
    )
  }

  if (any(arg < 0 | arg > 1)) {
    stop(
      sprintf("! `%s` must contain values between 0 and 1.", name),
      hint_line,
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Validate an open-interval probability scalar
#'
#' Checks that an argument is a single finite numeric value strictly between
#' zero and one. This is used for prediction confidence-level controls where
#' endpoint probabilities are not meaningful.
#'
#' @param arg Object to validate.
#' @param name Character scalar giving the argument name used in error messages.
#'
#' @return Invisibly returns \code{TRUE}.
#'
#' @keywords internal
#' @noRd
validate_open_probability_scalar <- function(arg, name) {
  # Reuse the package-wide numeric/finite/missing/length checks, then tighten
  # the allowed range from [0, 1] to the open interval (0, 1).
  validate_probability_vector(arg = arg, name = name, nvars = 1L)

  if (arg <= 0 || arg >= 1) {
    stop(
      sprintf("! `%s` must be a single numeric value strictly between 0 and 1.", name),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Validate a numeric scalar or vector argument
#'
#' Checks that an argument is numeric, has length one or the number of variables,
#' optionally allows \code{NULL}, optionally allows \code{NA}, and optionally
#' requires strictly positive finite values.
#'
#' @param arg Object to validate.
#' @param name Character scalar giving the argument name used in error messages.
#' @param nvars Number of variables expected when a vector is supplied.
#' @param allow_null Logical; whether \code{NULL} is permitted.
#' @param allow_na Logical; whether \code{NA} values are permitted.
#' @param strictly_positive Logical; whether non-missing values must be greater
#'   than zero.
#' @param hint Optional character scalar appended (on its own line) to any
#'   error message, e.g. to point the user to an alternative, per-variable way
#'   of supplying the same option. Default \code{NULL} (no hint added).
#'
#' @return Invisibly returns \code{TRUE}.
#'
#' @keywords internal
#' @noRd
validate_numeric_vector <- function(arg,
                                    name,
                                    nvars,
                                    allow_null = TRUE,
                                    allow_na = FALSE,
                                    strictly_positive = FALSE,
                                    hint = NULL) {
  hint_line <- if (!is.null(hint)) paste0("\n", hint) else ""

  if (is.null(arg)) {
    if (allow_null) return(invisible(TRUE))

    stop(
      sprintf("! `%s` must not be NULL.", name),
      hint_line,
      call. = FALSE
    )
  }

  if (!is.numeric(arg)) {
    stop(
      sprintf("! `%s` must be numeric.", name),
      sprintf("i Current type is: %s.", typeof(arg)),
      hint_line,
      call. = FALSE
    )
  }

  if (!length(arg) %in% c(1L, nvars)) {
    # nvars == 1L is the "scalar option" case (e.g. formula-interface globals
    # or fp() term attributes): phrase the allowed shape as a plain single
    # number rather than the slightly odd "vector of length 1".
    allowed_desc <- if (nvars == 1L) {
      if (allow_null) "NULL or a single number" else "a single number"
    } else if (allow_null) {
      sprintf("NULL, a single number, or a numeric vector of length %d", nvars)
    } else {
      sprintf("a single number or a numeric vector of length %d", nvars)
    }

    stop(
      sprintf("! `%s` must be %s; got length %d.", name, allowed_desc, length(arg)),
      hint_line,
      call. = FALSE
    )
  }

  if (!allow_na && anyNA(arg)) {
    stop(
      sprintf("! `%s` must not contain NA values.", name),
      hint_line,
      call. = FALSE
    )
  }

  finite_values <- if (allow_na) arg[!is.na(arg)] else arg

  if (length(finite_values) > 0L && any(!is.finite(finite_values))) {
    stop(
      sprintf("! `%s` must contain only finite values.", name),
      hint_line,
      call. = FALSE
    )
  }

  if (strictly_positive &&
      length(finite_values) > 0L &&
      any(finite_values <= 0)) {
    stop(
      sprintf("! `%s` must contain only positive non-missing values.", name),
      hint_line,
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Normalize a Per-Column Numeric Setting
#'
#' Validate and align a numeric preprocessing setting such as `shift` or
#' `scale` to the columns of a matrix. Unnamed scalars are recycled. Named
#' inputs are matched by name rather than by position and may optionally specify
#' only a subset of columns.
#'
#' When `allow_partial_named = TRUE`, a named scalar is deliberately treated as
#' a one-column specification rather than as a global scalar. The result is a
#' complete vector in `column_names` order: supplied entries retain their values
#' and unspecified entries are `NA_real_` sentinels for downstream automatic
#' estimation. Missing values present in the caller's input are validated before
#' those internal sentinels are introduced.
#'
#' @param value Object to normalize.
#' @param column_names Character vector of matrix column names.
#' @param argument_name Character scalar used in error messages.
#' @param allow_null Logical; whether `NULL` requests automatic values.
#' @param scalar_recycle Logical; whether a scalar is recycled to every column.
#' @param strictly_positive Logical; whether supplied non-missing values must be
#'   strictly positive.
#' @param allow_na Logical; whether missing values may be retained as internal
#'   automatic-value sentinels. Public matrix interfaces use `FALSE`; formula
#'   interfaces use `TRUE` only for vectors constructed internally.
#' @param allow_partial_named Logical; whether a named vector may specify only a
#'   subset of columns. Unspecified columns are returned as `NA_real_` so their
#'   values can be estimated by the existing preprocessing pipeline.
#'
#' @return A numeric vector named and ordered exactly like `column_names`.
#'   For partial named inputs, automatic entries are represented by `NA_real_`.
#'   For an unnamed scalar, every entry contains the recycled scalar value.
#' @keywords internal
#' @noRd
normalize_named_numeric_setting <- function(value,
                                            column_names,
                                            argument_name,
                                            allow_null = TRUE,
                                            scalar_recycle = TRUE,
                                            strictly_positive = FALSE,
                                            allow_na = FALSE,
                                            allow_partial_named = FALSE) {
  n_columns <- length(column_names)
  shape_message <- if (allow_partial_named) {
    sprintf(
      "`%s` must be a single unnamed numeric value or a named numeric vector for one or more columns of `x`.",
      argument_name
    )
  } else {
    sprintf(
      "`%s` must be a single numeric value or a named numeric vector with one value for each column of `x`.",
      argument_name
    )
  }

  if (is.null(value)) {
    if (!allow_null) {
      stop(sprintf("`%s` must not be NULL.", argument_name), call. = FALSE)
    }

    return(stats::setNames(rep(NA_real_, n_columns), column_names))
  }

  # Intercept-only formula models legitimately delegate a zero-column matrix
  # together with zero-length internally constructed setting vectors.
  if (n_columns == 0L && length(value) == 0L) {
    return(stats::setNames(numeric(0L), column_names))
  }

  # typeof() deliberately excludes logical and complex values. In particular,
  # TRUE/FALSE must not be accepted as numeric shift/scale settings.
  if (!is.numeric(value) || !typeof(value) %in% c("integer", "double")) {
    stop(
      sprintf("`%s` must contain numeric values, not %s values.",
              argument_name, typeof(value)),
      call. = FALSE
    )
  }

  # Missing values supplied at the public interface remain invalid. Missing
  # values introduced below for unspecified partial settings are internal
  # sentinels and are intentionally retained for automatic estimation.
  if (!allow_na && anyNA(value)) {
    stop(sprintf("`%s` must not contain missing values.", argument_name),
         call. = FALSE)
  }

  value_names <- names(value)
  named_partial <- allow_partial_named &&
    !is.null(value_names) &&
    length(value) > 0L

  if (length(value) == 1L && !named_partial) {
    if (!scalar_recycle && n_columns != 1L) {
      stop(shape_message, call. = FALSE)
    }

    value <- rep(as.numeric(value), n_columns)
    names(value) <- column_names
  } else {
    if (is.null(value_names) ||
        length(value_names) != length(value) ||
        anyNA(value_names) ||
        any(!nzchar(value_names))) {
      stop(shape_message, call. = FALSE)
    }

    if (anyDuplicated(value_names)) {
      duplicated_names <- unique(value_names[duplicated(value_names)])
      stop(
        sprintf(
          "`%s` names must be unique; duplicated name(s): %s.",
          argument_name,
          paste(duplicated_names, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    unknown_names <- setdiff(value_names, column_names)

    if (allow_partial_named) {
      if (length(unknown_names) > 0L) {
        stop(
          sprintf(
            "`%s` contains unknown column name(s): %s.",
            argument_name,
            paste(unknown_names, collapse = ", ")
          ),
          call. = FALSE
        )
      }

      normalized <- stats::setNames(rep(NA_real_, n_columns), column_names)
      normalized[value_names] <- as.numeric(value)
      value <- normalized
    } else {
      missing_names <- setdiff(column_names, value_names)

      if (length(value) != n_columns ||
          length(missing_names) > 0L ||
          length(unknown_names) > 0L) {
        details <- character(0L)
        if (length(missing_names) > 0L) {
          details <- c(
            details,
            sprintf("missing: %s", paste(missing_names, collapse = ", "))
          )
        }
        if (length(unknown_names) > 0L) {
          details <- c(
            details,
            sprintf("unknown: %s", paste(unknown_names, collapse = ", "))
          )
        }

        detail_suffix <- if (length(details) > 0L) {
          paste0(" (", paste(details, collapse = "; "), ")")
        } else {
          ""
        }

        stop(
          sprintf(
            "`%s` names must match `colnames(x)` exactly%s.",
            argument_name,
            detail_suffix
          ),
          call. = FALSE
        )
      }

      value <- as.numeric(value[column_names])
      names(value) <- column_names
    }
  }

  finite_values <- value[!is.na(value)]
  if (length(finite_values) > 0L && any(!is.finite(finite_values))) {
    stop(sprintf("`%s` must contain only finite values.", argument_name),
         call. = FALSE)
  }

  if (strictly_positive &&
      length(finite_values) > 0L &&
      any(finite_values <= 0)) {
    stop(sprintf("`%s` must contain only strictly positive values.",
                 argument_name),
         call. = FALSE)
  }

  value
}


#' Normalize a scalar or named partial override setting
#'
#' Matrix interfaces use this helper for settings such as `df`, `select`, and
#' `alpha`. An unnamed scalar is a global value and is recycled to every column.
#' A named vector is matched to `column_names`, may specify only a subset of
#' columns, and fills omitted columns from `default`. Unnamed vectors with more
#' than one value are rejected so that settings cannot be assigned positionally.
#'
#' @param value Numeric scalar or named numeric vector to normalize.
#' @param column_names Character vector of matrix column names.
#' @param default Numeric scalar or complete per-column numeric vector used for
#'   columns omitted from a named input.
#' @param argument_name Character scalar used in error messages.
#'
#' @return A list with three elements: `value`, a complete numeric vector named
#'   and ordered exactly like `column_names`; `supplied`, a named logical vector
#'   identifying entries explicitly supplied by name; and `global_scalar`, a
#'   logical scalar indicating that the caller supplied an unnamed scalar.
#' @keywords internal
#' @noRd
normalize_named_override_setting <- function(value,
                                             column_names,
                                             default,
                                             argument_name) {
  n_columns <- length(column_names)
  shape_message <- sprintf(
    "`%s` must be a single unnamed numeric value or a named numeric vector for one or more columns of `x`.",
    argument_name
  )

  # Intercept-only formula models can delegate a zero-column matrix together
  # with zero-length internally constructed setting vectors.
  if (n_columns == 0L && length(value) == 0L) {
    return(list(
      value = stats::setNames(numeric(0L), column_names),
      supplied = stats::setNames(logical(0L), column_names),
      global_scalar = FALSE
    ))
  }

  if (!is.numeric(value) || !typeof(value) %in% c("integer", "double")) {
    stop(
      sprintf("`%s` must contain numeric values, not %s values.",
              argument_name, typeof(value)),
      call. = FALSE
    )
  }

  if (length(value) == 0L) {
    stop(shape_message, call. = FALSE)
  }

  if (anyNA(value) || any(!is.finite(value))) {
    stop(
      sprintf("`%s` must contain only finite, non-missing values.",
              argument_name),
      call. = FALSE
    )
  }

  normalize_default <- function(x) {
    if (!is.numeric(x) || !typeof(x) %in% c("integer", "double") ||
        anyNA(x) || any(!is.finite(x))) {
      stop(
        sprintf("Internal default for `%s` is malformed.", argument_name),
        call. = FALSE
      )
    }

    if (length(x) == 1L) {
      return(stats::setNames(rep(as.numeric(x), n_columns), column_names))
    }

    if (length(x) != n_columns) {
      stop(
        sprintf("Internal default for `%s` has the wrong length.",
                argument_name),
        call. = FALSE
      )
    }

    x_names <- names(x)
    if (!is.null(x_names)) {
      if (anyNA(x_names) || any(!nzchar(x_names)) || anyDuplicated(x_names) ||
          !setequal(x_names, column_names)) {
        stop(
          sprintf("Internal default names for `%s` are malformed.",
                  argument_name),
          call. = FALSE
        )
      }
      x <- x[column_names]
    }

    stats::setNames(as.numeric(x), column_names)
  }

  value_names <- names(value)
  named_input <- !is.null(value_names) && length(value_names) > 0L

  if (length(value) == 1L && !named_input) {
    normalized <- stats::setNames(
      rep(as.numeric(value), n_columns),
      column_names
    )
    return(list(
      value = normalized,
      supplied = stats::setNames(rep(TRUE, n_columns), column_names),
      global_scalar = TRUE
    ))
  }

  if (is.null(value_names) ||
      length(value_names) != length(value) ||
      anyNA(value_names) ||
      any(!nzchar(value_names))) {
    stop(shape_message, call. = FALSE)
  }

  if (anyDuplicated(column_names)) {
    stop(
      sprintf(
        "`colnames(x)` must be unique when `%s` is supplied as a named vector.",
        argument_name
      ),
      call. = FALSE
    )
  }

  if (anyDuplicated(value_names)) {
    duplicated_names <- unique(value_names[duplicated(value_names)])
    stop(
      sprintf(
        "`%s` names must be unique; duplicated name(s): %s.",
        argument_name,
        paste(duplicated_names, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  unknown_names <- setdiff(value_names, column_names)
  if (length(unknown_names) > 0L) {
    stop(
      sprintf(
        "`%s` contains unknown column name(s): %s.",
        argument_name,
        paste(unknown_names, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  normalized <- normalize_default(default)
  normalized[value_names] <- as.numeric(value)
  supplied <- stats::setNames(column_names %in% value_names, column_names)

  list(
    value = normalized,
    supplied = supplied,
    global_scalar = FALSE
  )
}


#' Normalize a scalar or named partial logical setting
#'
#' Matrix interfaces use this helper for logical settings such as `center`.
#' An unnamed scalar is a global value and is recycled to every column. A named
#' vector is matched to `column_names`, may specify only a subset of columns,
#' and fills omitted columns from `default`. Unnamed vectors with more than one
#' value are rejected so that settings cannot be assigned positionally.
#'
#' @param value Logical scalar or named logical vector to normalize.
#' @param column_names Character vector of matrix column names.
#' @param default Logical scalar or complete per-column logical vector used for
#'   columns omitted from a named input.
#' @param argument_name Character scalar used in error messages.
#'
#' @return A logical vector named and ordered exactly like `column_names`.
#' @keywords internal
#' @noRd
normalize_named_logical_setting <- function(value,
                                            column_names,
                                            default,
                                            argument_name) {
  n_columns <- length(column_names)
  shape_message <- paste0(
    "`", argument_name, "` must be a single unnamed logical value or a named ",
    "logical vector for one or more columns of `x`. Unnamed logical vectors ",
    "with more than one value are not allowed; supply a scalar or name each value."
  )

  # Intercept-only formula models can delegate a zero-column matrix together
  # with zero-length internally constructed setting vectors.
  if (n_columns == 0L && length(value) == 0L) {
    return(stats::setNames(logical(0L), column_names))
  }

  if (!is.logical(value)) {
    stop(
      sprintf("`%s` must contain logical values, not %s values.",
              argument_name, typeof(value)),
      call. = FALSE
    )
  }

  if (length(value) == 0L) {
    stop(shape_message, call. = FALSE)
  }

  if (anyNA(value)) {
    stop(
      sprintf("`%s` must contain only non-missing logical values.",
              argument_name),
      call. = FALSE
    )
  }

  normalize_default <- function(x) {
    if (!is.logical(x) || anyNA(x)) {
      stop(
        sprintf("Internal default for `%s` is malformed.", argument_name),
        call. = FALSE
      )
    }

    if (length(x) == 1L) {
      return(stats::setNames(rep(x, n_columns), column_names))
    }

    if (length(x) != n_columns) {
      stop(
        sprintf("Internal default for `%s` has the wrong length.",
                argument_name),
        call. = FALSE
      )
    }

    x_names <- names(x)
    if (!is.null(x_names)) {
      if (anyNA(x_names) || any(!nzchar(x_names)) || anyDuplicated(x_names) ||
          !setequal(x_names, column_names)) {
        stop(
          sprintf("Internal default names for `%s` are malformed.",
                  argument_name),
          call. = FALSE
        )
      }
      x <- x[column_names]
    }

    stats::setNames(as.logical(x), column_names)
  }

  value_names <- names(value)
  named_input <- !is.null(value_names) && length(value_names) > 0L

  if (length(value) == 1L && !named_input) {
    return(stats::setNames(rep(value, n_columns), column_names))
  }

  if (is.null(value_names) ||
      length(value_names) != length(value) ||
      anyNA(value_names) ||
      any(!nzchar(value_names))) {
    stop(shape_message, call. = FALSE)
  }

  if (anyDuplicated(column_names)) {
    stop(
      sprintf(
        "`colnames(x)` must be unique when `%s` is supplied as a named vector.",
        argument_name
      ),
      call. = FALSE
    )
  }

  if (anyDuplicated(value_names)) {
    duplicated_names <- unique(value_names[duplicated(value_names)])
    stop(
      sprintf(
        "`%s` names must be unique; duplicated name(s): %s.",
        argument_name,
        paste(duplicated_names, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  unknown_names <- setdiff(value_names, column_names)
  if (length(unknown_names) > 0L) {
    stop(
      sprintf(
        "`%s` contains unknown column name(s): %s.",
        argument_name,
        paste(unknown_names, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  normalized <- normalize_default(default)
  normalized[value_names] <- as.logical(value)
  normalized
}


#' Validate a positive integer scalar argument
#'
#' Checks that an argument is a single finite positive integer.
#'
#' @param arg Object to validate.
#' @param name Character scalar giving the argument name used in error messages.
#'
#' @return Invisibly returns \code{TRUE}.
#'
#' @keywords internal
#' @noRd
validate_positive_integer_scalar <- function(arg, name) {
  if (!is.numeric(arg) || length(arg) != 1L || anyNA(arg) || !is.finite(arg)) {
    stop(
      sprintf("! `%s` must be a single positive integer.", name),
      call. = FALSE
    )
  }

  # Check integer-valuedness without coercing first: as.integer() warns and
  # returns NA for finite doubles outside the representable integer range.
  if (arg < 1 || arg > .Machine$integer.max || arg != floor(arg)) {
    stop(
      sprintf("! `%s` must be a single positive integer.", name),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Validate variable names against available columns
#'
#' Checks that an argument is a character vector of known column names.
#'
#' @param arg Object to validate.
#' @param name Character scalar giving the argument name used in error messages.
#' @param vnames Character vector of valid variable names.
#' @param allow_null Logical; whether \code{NULL} is permitted.
#'
#' @return Invisibly returns \code{TRUE}.
#'
#' @keywords internal
#' @noRd
validate_variable_names <- function(arg, name, vnames, allow_null = TRUE) {
  if (is.null(arg)) {
    if (allow_null) return(invisible(TRUE))

    stop(
      sprintf("! `%s` must not be NULL.", name),
      call. = FALSE
    )
  }

  if (!is.character(arg)) {
    stop(
      sprintf("! `%s` must be a character vector of column names in `x`.", name),
      sprintf("i Current type is: %s.", typeof(arg)),
      call. = FALSE
    )
  }

  if (anyNA(arg)) {
    stop(
      sprintf("! `%s` must not contain NA values.", name),
      call. = FALSE
    )
  }

  if (any(arg == "")) {
    stop(
      sprintf("! `%s` must not contain empty strings.", name),
      call. = FALSE
    )
  }

  unknown <- setdiff(unique(arg), vnames)

  if (length(unknown) > 0L) {
    stop(
      sprintf(
        "! Unknown variable name(s) in `%s`: %s.",
        name,
        paste(unknown, collapse = ", ")
      ),
      sprintf(
        "i Available column names in `x` are: %s.",
        paste(vnames, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Validate a single formula-interface scalar default
#'
#' Thin wrapper around \code{validate_numeric_vector()} /
#' \code{validate_logical_vector()} for the scalar (\code{nvars = 1L}) options
#' used as global defaults in \code{mfp2.formula()} (\code{df}, \code{alpha},
#' \code{select}, \code{shift}, \code{scale}, \code{center}). Adds a `hint`
#' pointing the user to the per-variable \code{fp()} override syntax.
#'
#' @param arg Object to validate.
#' @param arg_name Character scalar giving the argument name used in error
#'   messages.
#' @param type Either \code{"numeric"} or \code{"logical"}.
#' @param allow_null Logical; whether \code{NULL} is permitted. Default
#'   \code{FALSE}.
#' @param allow_na Logical; whether \code{NA} is permitted. Default
#'   \code{FALSE}.
#' @param strictly_positive Logical; whether a non-missing numeric value must
#'   be greater than zero (ignored for \code{type = "logical"}). Default
#'   \code{FALSE}.
#' @param hint Optional character scalar appended (on its own line) to any
#'   error message. Default \code{NULL}.
#'
#' @return Invisibly returns \code{TRUE}.
#'
#' @keywords internal
#' @noRd
validate_formula_scalar <- function(arg,
                                    arg_name,
                                    type,
                                    allow_null = FALSE,
                                    allow_na = FALSE,
                                    strictly_positive = FALSE,
                                    hint = NULL) {
  switch(
    type,
    numeric = validate_numeric_vector(
      arg = arg, name = arg_name, nvars = 1L,
      allow_null = allow_null, allow_na = allow_na,
      strictly_positive = strictly_positive, hint = hint
    ),
    logical = validate_logical_vector(
      arg = arg, name = arg_name, allowed_lengths = 1L,
      allow_null = allow_null, hint = hint
    ),
    stop("Internal error: unsupported validator type.", call. = FALSE)
  )
}


#' Validate a single formula-interface probability default
#'
#' Thin wrapper around \code{validate_probability_vector()} for the scalar
#' (\code{nvars = 1L}) significance-level options used as global defaults in
#' \code{mfp2.formula()} (\code{alpha}, \code{select}). Adds a `hint` pointing
#' the user to the per-variable \code{fp()} override syntax.
#'
#' @param arg Object to validate.
#' @param arg_name Character scalar giving the argument name used in error
#'   messages.
#' @param hint Optional character scalar appended (on its own line) to any
#'   error message. Default \code{NULL}.
#'
#' @return Invisibly returns \code{TRUE}.
#'
#' @keywords internal
#' @noRd
validate_formula_probability <- function(arg, arg_name, hint = NULL) {
  validate_probability_vector(arg = arg, name = arg_name, nvars = 1L, hint = hint)
}


#' Validate a single scalar attribute supplied inside an `fp()` term
#'
#' Thin wrapper around \code{validate_numeric_vector()} /
#' \code{validate_logical_vector()} for the scalar (\code{nvars = 1L})
#' attributes attached to an individual \code{fp()} term (\code{df},
#' \code{alpha}, \code{select}, \code{shift}, \code{scale}, \code{center},
#' \code{acdx}, \code{zero}, \code{catzero}, \code{spike},
#' \code{force_max_fp}). \code{fp()} only validates shape/type issues that
#' could corrupt attribute extraction; full range and semantic validation is
#' performed later by \code{mfp2.default()} once formula options have been
#' expanded to per-variable vectors.
#'
#' @param arg Object to validate.
#' @param arg_name Character scalar giving the `fp()` argument name used in
#'   error messages (e.g. \code{"df"}).
#' @param var_name Character scalar giving the name of the variable passed to
#'   \code{fp()} (e.g. \code{"age"} for \code{fp(age)}), used to identify
#'   which `fp()` call raised the error.
#' @param type Either \code{"numeric"} or \code{"logical"}.
#' @param allow_null Logical; whether \code{NULL} is permitted. Default
#'   \code{FALSE}.
#' @param allow_na Logical; whether \code{NA} is permitted. Default
#'   \code{FALSE}.
#'
#' @return Invisibly returns \code{TRUE}.
#'
#' @keywords internal
#' @noRd
validate_scalar_fp <- function(arg, arg_name, var_name, type,
                               allow_null = FALSE, allow_na = FALSE) {
  hint <- sprintf("i This is the `%s` argument inside fp(%s).", arg_name, var_name)

  switch(
    type,
    numeric = validate_numeric_vector(
      arg = arg, name = arg_name, nvars = 1L,
      allow_null = allow_null, allow_na = allow_na, hint = hint
    ),
    logical = validate_logical_vector(
      arg = arg, name = arg_name, allowed_lengths = 1L,
      allow_null = allow_null, hint = hint
    ),
    stop("Internal error: unsupported validator type.", call. = FALSE)
  )
}


#' Normalize and Bind Formula Specials
#'
#' Rewrites formula-special calls that base R does not recognize when
#' namespace-qualified, e.g. survival::strata(x) -> strata(x) and
#' stats::offset(x) -> offset(x). It also binds the bare `fp()` and `fp2()`
#' names to this package's implementations so another attached package cannot
#' intercept their evaluation. The original call object should still be stored
#' on the fitted object; this normalized formula is for internal parsing.
#'
#' @param formula A model formula.
#'
#' @return A formula with selected namespace-qualified specials rewritten to
#'   their bare names.
#'
#' @keywords internal
#' @noRd
normalize_formula_special_namespaces <- function(formula) {
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula.", call. = FALSE)
  }

  is_namespaced_symbol <- function(x, namespace, name) {
    is.call(x) &&
      length(x) == 3L &&
      identical(x[[1L]], as.name("::")) &&
      identical(as.character(x[[2L]]), namespace) &&
      identical(as.character(x[[3L]]), name)
  }

  rewrite_call <- function(expr) {
    if (!is.call(expr)) {
      return(expr)
    }

    fun <- expr[[1L]]

    if (is_namespaced_symbol(fun, "survival", "strata")) {
      expr[[1L]] <- as.name("strata")
    } else if (is_namespaced_symbol(fun, "stats", "offset")) {
      expr[[1L]] <- as.name("offset")
    }

    if (length(expr) > 1L) {
      for (i in seq.int(2L, length(expr))) {
        expr[[i]] <- rewrite_call(expr[[i]])
      }
    }

    expr
  }

  out <- formula

  if (length(out) >= 2L) {
    out[[2L]] <- rewrite_call(out[[2L]])
  }

  if (length(out) >= 3L) {
    out[[3L]] <- rewrite_call(out[[3L]])
  }

  # Give model.frame()/terms() a reliable lookup path for package-owned formula
  # helpers while preserving all other user-scope lookups through the parent.
  # Ensure fp() in the formula resolves to mfp2's version even if mfp is loaded.
  # Both packages export fp(); without this, model.frame() evaluates fp() via
  # the search path and may find mfp::fp(), which lacks the attributes mfp2 expects.
  env <- new.env(parent = environment(formula))
  env$fp <- fp2
  env$fp2 <- fp2
  env$strata <- survival::strata
  env$offset <- stats::offset

  environment(out) <- env
  out
}
