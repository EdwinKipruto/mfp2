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

  if (arg != as.integer(arg) || arg < 1L) {
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


#' Normalize Namespace-Qualified Formula Specials
#'
#' Rewrites formula-special calls that base R does not recognize when
#' namespace-qualified, e.g. survival::strata(x) -> strata(x) and
#' stats::offset(x) -> offset(x). The original call object should still be
#' stored on the fitted object; this normalized formula is for internal parsing.
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

  # Give model.frame()/terms() a reliable lookup path for the bare specials
  # after rewriting, while preserving all user-scope lookups through the parent.
  env <- new.env(parent = environment(formula))
  env$strata <- survival::strata
  env$offset <- stats::offset

  environment(out) <- env
  out
}