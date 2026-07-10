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