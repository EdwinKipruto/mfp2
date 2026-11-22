#' Restricted cubic spline term
#'
#' `rcs()` specifies a nested restricted cubic spline in an `mfp2()` formula.
#' `df` is the maximum number of coefficients, including the linear term.
#' Knots include both boundary knots; hence `length(knots) == df + 1`.
#' When knots are omitted, knots follow the placement rules used by
#' `rms::rcs()`, without requiring rms or Hmisc. The fixed basis is used in
#' column order so simpler candidate models are nested in the maximum.
#' For `catzero` or `spike`, the predictor must be nonnegative and automatic
#' knots are calculated from positive fitting observations. Spike eligibility
#' follows `min_saz_prop` in `mfp2()`; an ineligible spike request retains
#' `catzero` only if it was explicitly set.
#' If tied values yield fewer knots, the effective maximum `df` is reduced.
#' The returned matrix also stores the predictor expression, requested `df`,
#' original knot argument, and zero-handling options as attributes so `mfp2()`
#' can resolve the fitted representation.
#'
#' @param x A numeric predictor.
#' @param df Maximum spline degrees of freedom, at least 2.
#' @param knots Optional strictly increasing knot locations, including boundaries.
#' @param center Center retained spline basis columns using the fitted-sample
#'   mean for continuous columns and the minimum for binary columns. Defaults
#'   to `TRUE`, as for `fp()`. Set to `FALSE` to leave the basis uncentered.
#'   Centering changes the intercept, not fitted values.
#' @param catzero Add an exact-zero indicator whenever the term is retained.
#' @param spike Select the spline and zero-indicator components separately.
#' @return A matrix with one row per observation, a linear column followed by
#'   nonlinear columns (up to `df - 1`), and an optional `catzero` column. The
#'   formula interface may rebuild it using the knots selected for the fit.
#' @export
rcs <- function(x, df = 4L, knots = NULL, catzero = FALSE, spike = FALSE,
                center = TRUE) {
  predictor_name <- deparse(substitute(x))
  if (!is.numeric(x) || !is.null(dim(x)) || is.factor(x)) {
    stop("`rcs()` requires a numeric vector.", call. = FALSE)
  }
  if (length(df) != 1L || !is.numeric(df) || is.na(df) ||
      !is.finite(df) || df < 2L || df > .Machine$integer.max - 1L ||
      df != as.integer(df)) {
    stop("`rcs()` requires an integer `df` of at least 2.", call. = FALSE)
  }
  if (!is.logical(catzero) || length(catzero) != 1L || is.na(catzero) ||
      !is.logical(spike) || length(spike) != 1L || is.na(spike)) {
    stop("`catzero` and `spike` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(center) || length(center) != 1L || is.na(center)) {
    stop("`center` must be TRUE or FALSE.", call. = FALSE)
  }
  if ((catzero || spike) && any(x < 0, na.rm = TRUE)) {
    stop("`catzero` and `spike` require a nonnegative predictor.", call. = FALSE)
  }
  supplied_knots <- knots
  # Zero is handled by its own indicator, so it does not influence the
  # automatically chosen spline knots when either zero option is requested.
  knot_observations <- if (catzero || spike) x[x > 0 & !is.na(x)] else x
  if (is.null(knots)) {
    # model.frame() evaluates this wrapper on prediction rows too. If those
    # rows cannot support automatic knots, make a provisional valid basis;
    # predict.mfp2() replaces its columns using the knots saved at fitting.
    knots <- tryCatch(rcs_default_knots(knot_observations, df), error = function(e) {
      finite_observations <- x[is.finite(x)]
      fallback_centre <- if (length(finite_observations)) {
        finite_observations[[1L]]
      } else 0
      seq(fallback_centre - 1, fallback_centre + 1, length.out = df + 1L)
    })
  } else {
    rcs_validate_knots(knots, df)
  }
  basis <- rcs_basis(x, knots, indicator = catzero || spike)
  attr(basis, "rcs_name") <- predictor_name
  attr(basis, "rcs_df") <- as.integer(df)
  attr(basis, "rcs_supplied_knots") <- supplied_knots
  attr(basis, "rcs_catzero") <- catzero
  attr(basis, "rcs_spike") <- spike
  attr(basis, "rcs_center") <- center
  basis
}

#' Check knot count, ordering, and finiteness
#'
#' The vector contains both boundary knots and must have one more entry than
#' the number of spline basis columns. An invalid vector raises an error.
#'
#' @param knots Numeric knot locations, including both boundary knots.
#' @param df Number of spline basis columns, including the linear column.
#' @return The validated knot vector, invisibly.
#' @keywords internal
#' @noRd
rcs_validate_knots <- function(knots, df) {
  if (!is.numeric(knots) || length(knots) != df + 1L ||
      anyNA(knots) || any(!is.finite(knots)) ||
      any(diff(knots) <= 0)) {
    stop("`rcs()` needs `df + 1` finite, strictly increasing knots, including boundaries.",
         call. = FALSE)
  }
  invisible(knots)
}

#' Choose automatic restricted cubic spline knots
#'
#' Uses finite observations and follows the knot-placement rules used by
#' `rms::rcs()`: choose inner distinct values when there are few, otherwise
#' use quantiles. The `fractied` option controls special handling of frequent
#' extreme values. After duplicate knots are removed, there may be fewer than
#' the requested `df + 1` knots; fewer than three is an error.
#'
#' @param x Observations used to place knots; callers exclude zero when a
#'   separate zero indicator is requested.
#' @param df Requested maximum number of spline basis columns.
#' @return A strictly increasing numeric knot vector, including boundaries.
#' @keywords internal
#' @noRd
rcs_default_knots <- function(x, df) {
  placement_observations <- x[is.finite(x)]
  observation_count <- length(placement_observations)
  requested_knot_count <- df + 1L
  if (observation_count < 6L) {
    stop("At least six observations are required for automatic RCS knots.",
         call. = FALSE)
  }
  distinct_observations <- sort(unique(placement_observations))
  # When few distinct values are available, use every interior value and
  # leave the observed minimum and maximum outside the knot vector.
  if (length(distinct_observations) - 2L <= requested_knot_count) {
    knots <- distinct_observations[-c(1L, length(distinct_observations))]
  } else {
    outer_probability <- if (requested_knot_count > 6L) {
      0.025
    } else if (requested_knot_count > 3L) {
      0.05
    } else 0.10
    fractied <- getOption("fractied", 0.05)
    if (!is.numeric(fractied) || length(fractied) != 1L ||
        is.na(fractied) || !is.finite(fractied)) {
      stop("`fractied` must be one finite numeric option.", call. = FALSE)
    }
    lower_reserved_knot <- upper_reserved_knot <- numeric(0L)
    remaining_knot_count <- requested_knot_count
    # When an extreme value is common but interior values are not, reserve
    # a knot just inside that extreme and exclude it from quantile placement.
    if (fractied > 0 && fractied < 1) {
      value_proportions <- as.numeric(table(placement_observations)) / observation_count
      if (max(value_proportions[-c(1L, length(value_proportions))]) < fractied) {
        if (value_proportions[[1L]] >= fractied) {
          lower_reserved_knot <- min(placement_observations[
            placement_observations > min(placement_observations)
          ])
          placement_observations <- placement_observations[
            placement_observations > lower_reserved_knot
          ]
          remaining_knot_count <- remaining_knot_count - 1L
        }
        if (value_proportions[[length(value_proportions)]] >= fractied) {
          upper_reserved_knot <- max(placement_observations[
            placement_observations < max(placement_observations)
          ])
          placement_observations <- placement_observations[
            placement_observations < upper_reserved_knot
          ]
          remaining_knot_count <- remaining_knot_count - 1L
        }
      }
    }
    if (remaining_knot_count == 1L) {
      knots <- stats::median(placement_observations)
    } else {
      # Quantiles place the remaining knots across the observations left
      # after any extreme-tie adjustment.
      quantile_probs <- if (remaining_knot_count == 2L) {
        seq(0.5, 1 - outer_probability, length.out = remaining_knot_count)
      } else {
        seq(outer_probability, 1 - outer_probability,
            length.out = remaining_knot_count)
      }
      knots <- as.numeric(stats::quantile(
        placement_observations, probs = quantile_probs, names = FALSE
      ))
      # Tied quantiles can leave too few locations. Try a finer grid; if
      # that still fails with a reserved lower knot, construct three knots.
      if (length(unique(knots)) < min(remaining_knot_count, 3L)) {
        knots <- as.numeric(stats::quantile(
          placement_observations,
          probs = seq(outer_probability, 1 - outer_probability,
                      length.out = 2L * remaining_knot_count),
          names = FALSE
        ))
        if (length(lower_reserved_knot) && length(unique(knots)) < 3L) {
          midpoint <- if (length(upper_reserved_knot)) {
            (lower_reserved_knot + upper_reserved_knot) / 2
          } else stats::median(placement_observations)
          knots <- c(lower_reserved_knot, midpoint,
                     if (length(upper_reserved_knot)) upper_reserved_knot else
                       as.numeric(stats::quantile(placement_observations,
                                                  1 - outer_probability)))
        }
      }
      # For fewer than 100 remaining observations, use the
      # fifth values from the ends unless an extreme-tie knot was reserved.
      if (length(placement_observations) < 100L) {
        sorted_observations <- sort(placement_observations)
        if (!length(lower_reserved_knot)) knots[[1L]] <- sorted_observations[[5L]]
        if (!length(upper_reserved_knot)) {
          knots[[remaining_knot_count]] <- sorted_observations[[
            length(sorted_observations) - 4L
          ]]
        }
      }
    }
    knots <- c(lower_reserved_knot, knots, upper_reserved_knot)
  }
  # Ties can lower the effective df; validation uses the number that remains.
  knots <- sort(unique(as.numeric(knots)))
  if (length(knots) < 3L) {
    stop("Automatic knots yielded fewer than three distinct locations.",
         call. = FALSE)
  }
  if (length(knots) > df + 1L) {
    stop("Automatic knot placement produced more knots than requested; supply explicit knots.",
         call. = FALSE)
  }
  rcs_validate_knots(knots, length(knots) - 1L)
  knots
}

#' Construct the restricted cubic spline design matrix
#'
#' The first column is the original predictor. Each subsequent column uses a
#' truncated cubic at one knot, adjusted with the final two knots so the
#' resulting spline is linear beyond the boundary knots. Columns use the same
#' normalization as the `rms`/`Hmisc` basis. Taking the first columns gives
#' smaller candidate models from the same fixed knot vector.
#'
#' @param x Numeric predictor values, with one matrix row per value.
#' @param knots Strictly increasing knots, including both boundary knots.
#' @param indicator Whether to append a column for exact zeros.
#' @return A numeric matrix with `length(knots) - 1` spline columns named
#'   `linear`, `nonlinear1`, and so on; when `indicator` is true, a final
#'   `catzero` column is appended. Missing predictor values yield missing rows.
#' @keywords internal
#' @noRd
rcs_basis <- function(x, knots, indicator = FALSE) {
  df <- length(knots) - 1L
  rcs_validate_knots(knots, df)
  # The last two knots anchor the correction that makes each nonlinear
  # component linear beyond the boundary knots. The scale normalizes its
  # cubed distances to match the rms/Hmisc basis.
  last_knot <- knots[[length(knots)]]
  penultimate_knot <- knots[[length(knots) - 1L]]
  knot_scale <- (last_knot - knots[[1L]])^(2 / 3)
  basis <- matrix(NA_real_, length(x), df)
  basis[, 1L] <- x
  # Add one nonlinear column for each knot before the last two. The same
  # columns remain in place when a candidate model uses fewer of them.
  for (j in seq_len(df - 1L)) {
    basis[, j + 1L] <- pmax((x - knots[[j]]) / knot_scale, 0)^3 +
      ((penultimate_knot - knots[[j]]) * pmax((x - last_knot) / knot_scale, 0)^3 -
       (last_knot - knots[[j]]) * pmax((x - penultimate_knot) / knot_scale, 0)^3) /
      (last_knot - penultimate_knot)
  }
  colnames(basis) <- c("linear", paste0("nonlinear", seq_len(df - 1L)))
  if (indicator) basis <- cbind(basis, catzero = as.integer(x == 0))
  basis
}

#' Identify an RCS call in a formula term label
#'
#' Matches a whole `rcs(...)` label with an optional `mfp2::` or `mfp2:::`
#' prefix. This checks the label's shape; it does not parse its arguments.
#'
#' @param z Character vector of formula term labels.
#' @return A logical vector indicating which labels match.
#' @keywords internal
#' @noRd
is_rcs_term <- function(z) {
  grepl("^(mfp2:::?)?rcs\\(.*\\)$", z)
}

#' Describe the selected spline and zero-indicator components
#'
#' A missing `df` means the term was excluded; zero means only its zero
#' indicator was retained; one means a linear predictor. Larger values
#' describe a restricted cubic spline with that many columns.
#'
#' @param df Selected spline column count; `NA` for exclusion and zero for
#'   the zero-indicator-only form.
#' @param indicator Whether to add the zero-indicator suffix when `df > 0`.
#' @return A human-readable label for the selected form.
#' @keywords internal
#' @noRd
rcs_form_label <- function(df, indicator = FALSE) {
  if (is.na(df)) return("excluded")
  if (df == 0L) return("zero indicator only")
  spline_label <- if (df == 1L) "linear" else sprintf("RCS(df=%d)", df)
  if (indicator) paste0(spline_label, " + zero indicator") else spline_label
}

#' Label a fitted spline coefficient by its basis column
#'
#' The linear column keeps the predictor name. Nonlinear columns include
#' their position and the fitted knots; the optional final column denotes
#' the exact-zero indicator.
#'
#' @param term Original predictor name.
#' @param position One-based column position in the fitted design matrix.
#' @param knots Full knot vector used in the fit; position `length(knots)`
#'   identifies the optional zero-indicator column.
#' @return A human-readable column label, including knots for nonlinear terms.
#' @keywords internal
#' @noRd
rcs_column_label <- function(term, position, knots) {
  if (position == length(knots)) return(sprintf("I(%s = 0)", term))
  if (position == 1L) return(term)
  sprintf("RCS_%d(%s; knots = c(%s))", position - 1L, term,
          paste(format(knots, digits = 15, trim = TRUE), collapse = ", "))
}
