#' Normalize One Candidate FP Power Vector
#'
#' Validates and normalizes one vector of candidate fractional-polynomial base
#' powers. The input represents candidate base powers, not a selected FP power
#' vector. Duplicate values are removed because repeated selected powers are
#' generated later by the FP candidate generator for degree-2 and higher FP
#' models.
#'
#' @param powers Numeric vector of candidate FP base powers.
#' @param context Character scalar used in error messages.
#'
#' @return Sorted unique numeric vector of finite candidate powers.
#'
#' @keywords internal
#' @noRd
normalize_fp_power_vector <- function(powers, context = "`powers`") {
  if (!is.numeric(powers)) {
    stop(
      paste0(context, " must be a numeric vector."),
      call. = FALSE
    )
  }
  
  if (length(powers) == 0L || anyNA(powers) || any(!is.finite(powers))) {
    stop(
      paste0(
        context,
        " must contain at least one finite, non-missing numeric value."
      ),
      call. = FALSE
    )
  }
  
  sort(unique(as.numeric(powers)))
}

#' Validate and Normalize a Named Candidate FP Power List
#'
#' Validates a named list of candidate fractional-polynomial base-power sets and
#' returns a complete named list, one element per predictor. The \code{powers}
#' argument defines candidate base powers, not selected FP power vectors.
#' Duplicate values are removed.
#'
#' The optional \code{df} argument is accepted for compatibility with older
#' internal callers, but this helper does not perform df-dependent closed-test
#' validation. That validation must happen after early df modifications such as
#' ACD forcing and SAZ positive-component df capping, via
#' \code{validate_mfp_candidate_powers()} in \code{fit_mfp()}.
#'
#' @param powers User-supplied powers argument. Must be \code{NULL} or a named
#'   list of numeric vectors.
#' @param vnames Character vector of valid predictor names.
#' @param default_powers Numeric vector of default candidate powers.
#' @param arg_name Character scalar used in error messages.
#' @param df Optional numeric scalar or vector aligned with \code{vnames}.
#'   Accepted for compatibility; df-dependent validation is deferred to
#'   \code{validate_mfp_candidate_powers()}.
#'
#' @return Named list of normalized candidate-power vectors, one element per
#'   predictor in \code{vnames}.
#'
#' @keywords internal
#' @noRd
validate_fp_power_list <- function(powers,
                                   vnames,
                                   default_powers = c(-2, -1, -0.5, 0, 0.5, 1, 2, 3),
                                   arg_name = "powers",
                                   df = NULL) {
  if (!is.character(vnames) || length(vnames) == 0L ||
      anyNA(vnames) || any(!nzchar(vnames))) {
    stop(
      "`vnames` must be a non-empty character vector of predictor names.",
      call. = FALSE
    )
  }
  
  # Do not perform df-dependent candidate-power validation here. The effective
  # df can still change later (for example ACD forces df = 4 and retained SAZ
  # variables may have df capped using only their positive component). The closed-
  # test check is therefore deferred to validate_mfp_candidate_powers() inside
  # fit_mfp(), after those modifications are complete.
  if (!is.null(df)) {
    invisible(df)
  }
  
  default_powers <- normalize_fp_power_vector(
    default_powers,
    context = "`default_powers`"
  )
  
  power_list <- stats::setNames(
    replicate(length(vnames), default_powers, simplify = FALSE),
    vnames
  )
  
  if (is.null(powers)) {
    return(power_list)
  }
  
  if (!is.list(powers)) {
    stop(
      paste0("! `", arg_name, "` must be a named list."),
      call. = FALSE
    )
  }
  
  pow_names <- names(powers)
  
  if (is.null(pow_names) || length(pow_names) != length(powers) ||
      anyNA(pow_names) || any(!nzchar(pow_names))) {
    stop(
      paste0("! Every element of `", arg_name, "` must have a non-empty name."),
      call. = FALSE
    )
  }
  
  if (anyDuplicated(pow_names)) {
    dup <- unique(pow_names[duplicated(pow_names)])
    stop(
      paste0(
        "! `", arg_name, "` must not contain duplicated names. ",
        "Duplicated name(s): ",
        paste(dup, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }
  
  unknown_names <- setdiff(pow_names, vnames)
  
  if (length(unknown_names) > 0L) {
    stop(
      paste0(
        "! The following names in `", arg_name,
        "` do not match any predictor column: ",
        paste(unknown_names, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }
  
  powers_norm <- stats::setNames(
    vector("list", length(powers)),
    pow_names
  )
  
  for (nm in pow_names) {
    powers_norm[[nm]] <- normalize_fp_power_vector(
      powers[[nm]],
      context = paste0("`", arg_name, "[[\"", nm, "\"]]`")
    )
  }
  
  power_list[names(powers_norm)] <- powers_norm
  
  
  power_list
}