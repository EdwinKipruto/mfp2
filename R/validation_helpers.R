#' Validate a logical scalar or vector argument
#'
#' Checks that an argument is logical, contains no missing values, and has one
#' of the allowed lengths.
#'
#' @param arg Object to validate.
#' @param name Character scalar giving the argument name used in error messages.
#' @param allowed_lengths Integer vector of permitted lengths.
#'
#' @return Invisibly returns \code{TRUE}.
#'
#' @keywords internal
#' @noRd
validate_logical_vector <- function(arg, name, allowed_lengths) {
  if (!is.logical(arg)) {
    stop(
      sprintf("! `%s` must be logical.", name),
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
  
  if (!length(arg) %in% allowed_lengths) {
    stop(
      sprintf(
        "! `%s` must have length %s; got length %d.",
        name,
        paste(allowed_lengths, collapse = " or "),
        length(arg)
      ),
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
#'
#' @return Invisibly returns \code{TRUE}.
#'
#' @keywords internal
#' @noRd
validate_probability_vector <- function(arg, name, nvars) {
  if (!is.numeric(arg)) {
    stop(
      sprintf("! `%s` must be numeric.", name),
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
  
  if (any(!is.finite(arg))) {
    stop(
      sprintf("! `%s` must contain only finite values.", name),
      call. = FALSE
    )
  }
  
  if (!length(arg) %in% c(1L, nvars)) {
    stop(
      sprintf(
        "! `%s` must be a single number or a numeric vector of length %d; got length %d.",
        name, nvars, length(arg)
      ),
      call. = FALSE
    )
  }
  
  if (any(arg < 0 | arg > 1)) {
    stop(
      sprintf("! `%s` must contain values between 0 and 1.", name),
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
                                    strictly_positive = FALSE) {
  if (is.null(arg)) {
    if (allow_null) return(invisible(TRUE))
    
    stop(
      sprintf("! `%s` must not be NULL.", name),
      call. = FALSE
    )
  }
  
  if (!is.numeric(arg)) {
    stop(
      sprintf("! `%s` must be numeric.", name),
      sprintf("i Current type is: %s.", typeof(arg)),
      call. = FALSE
    )
  }
  
  if (!length(arg) %in% c(1L, nvars)) {
    stop(
      sprintf(
        "! `%s` must be NULL, a single number, or a numeric vector of length %d; got length %d.",
        name, nvars, length(arg)
      ),
      call. = FALSE
    )
  }
  
  if (!allow_na && anyNA(arg)) {
    stop(
      sprintf("! `%s` must not contain NA values.", name),
      call. = FALSE
    )
  }
  
  finite_values <- if (allow_na) arg[!is.na(arg)] else arg
  
  if (length(finite_values) > 0L && any(!is.finite(finite_values))) {
    stop(
      sprintf("! `%s` must contain only finite values.", name),
      call. = FALSE
    )
  }
  
  if (strictly_positive &&
      length(finite_values) > 0L &&
      any(finite_values <= 0)) {
    stop(
      sprintf("! `%s` must contain only positive non-missing values.", name),
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