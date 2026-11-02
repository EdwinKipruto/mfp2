#' Calculates the total number of estimated fractional polynomial powers.
#'
#' This function calculates the number of estimated fractional polynomial power
#' parameters in a list of selected powers. It is intended for use when adding
#' extra FP-power penalties to model-selection criteria, not for calculating the
#' total number of regression coefficients.
#'
#' A single power \code{1} fitted as the fixed-linear row contributes zero
#' FP-power parameters. If the same numerical power was selected from the FP1
#' candidate class, set the corresponding element of \code{searched_fp} to
#' \code{TRUE}; it then contributes one estimated-power parameter. Repeated
#' powers, including \code{c(1, 1)}, are treated as FP terms and contribute
#' their length. ACD encodings such as \code{c(2, NA)} or
#' \code{c(NA, 2)} contribute one FP power after removing the structural
#' \code{NA}.
#'
#' If \code{spike_decision} is supplied, variables with decision
#' \code{saz_decision_codes[["binary_only"]]} are treated as binary-only spike
#' terms and contribute zero FP powers, regardless of their stored powers.
#'
#' @param x A named list of fractional polynomial powers.
#' @param spike_decision Optional named integer vector with spike decisions.
#'   Codes are defined by \code{saz_decision_codes}: \code{cont_binary} for
#'   FP/linear plus binary spike indicator, \code{continuous_only} for FP/linear
#'   only, and \code{binary_only} for binary spike indicator only.
#' @param searched_fp Optional named logical vector indicating whether each
#'   selected form came from an FP candidate class rather than the fixed-linear
#'   row. Defaults to \code{FALSE} for every term.
#'
#' @return Integer value denoting the total number of estimated FP powers.
#'
#' @keywords internal
#' @noRd
calculate_number_fp_powers <- function(x, spike_decision = NULL,
                                       searched_fp = NULL) {
  if (is.null(x) || length(x) == 0L) {
    return(0L)
  }

  if (!is.list(x)) {
    stop("`x` must be a list of fractional polynomial powers.", call. = FALSE)
  }

  if (!is.null(spike_decision)) {
    if (is.null(names(x)) || is.null(names(spike_decision))) {
      stop(
        "`x` and `spike_decision` must both be named when `spike_decision` is supplied.",
        call. = FALSE
      )
    }

    missing_spike <- setdiff(names(x), names(spike_decision))
    if (length(missing_spike) > 0L) {
      stop(
        sprintf(
          "`spike_decision` is missing entries for: %s.",
          paste0(missing_spike, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    spike_decision <- spike_decision[names(x)]

    if (anyNA(spike_decision) || any(!spike_decision %in% unname(saz_decision_codes))) {
      stop(
        "`spike_decision` must contain only values 1L, 2L, or 3L.",
        call. = FALSE
      )
    }
  } else {
    spike_decision <- setNames(
      rep(saz_decision_codes[["continuous_only"]], length(x)),
      names(x)
    )
  }

  if (is.null(searched_fp)) {
    searched_fp <- setNames(rep(FALSE, length(x)), names(x))
  } else {
    if (is.null(names(x)) || is.null(names(searched_fp))) {
      stop(
        "`x` and `searched_fp` must both be named when `searched_fp` is supplied.",
        call. = FALSE
      )
    }

    missing_searched <- setdiff(names(x), names(searched_fp))
    if (length(missing_searched) > 0L) {
      stop(
        sprintf(
          "`searched_fp` is missing entries for: %s.",
          paste0(missing_searched, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    searched_fp <- searched_fp[names(x)]
    if (!is.logical(searched_fp) || anyNA(searched_fp)) {
      stop("`searched_fp` must contain only TRUE or FALSE.", call. = FALSE)
    }
  }

  count_one <- function(p, spike_decision, searched_fp) {
    # Binary-only spike: no FP power parameter.
    if (identical(as.integer(spike_decision), saz_decision_codes[["binary_only"]])) {
      return(0L)
    }

    if (is.null(p) || length(p) == 0L) {
      return(0L)
    }

    # Eliminated variable, e.g. NA or c(NA, NA).
    if (all(is.na(p))) {
      return(0L)
    }

    # Remove structural NA values from ACD encodings, e.g. c(2, NA), c(NA, 2).
    p <- p[!is.na(p)]

    if (length(p) == 0L) {
      return(0L)
    }

    # Linear term: no estimated FP power parameter.
    if (length(p) == 1L && p == 1 && !searched_fp) {
      return(0L)
    }

    # Non-linear FP term. Repeated powers such as c(1, 1) count as FP2.
    length(p)
  }

  sum(
    mapply(
      count_one,
      p = x,
      spike_decision = spike_decision,
      searched_fp = searched_fp,
      SIMPLIFY = TRUE,
      USE.NAMES = FALSE
    )
  )
}
