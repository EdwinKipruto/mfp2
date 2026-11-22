#' Cubic B-spline term with nested, fixed-knot candidates
#'
#' `bspline()` adds a cubic B-spline term to an `mfp2()` formula. `df` is the
#' maximum number of spline coefficients (excluding the model intercept), so
#' there are `df - 3` interior knots. Knots are selected once from fitting
#' observations and reused for every candidate and prediction. Candidate
#' knot sets grow symmetrically from the middle: for `df = 6`, the candidates
#' are linear, a one-knot spline with four coefficients, and the three-knot
#' spline with six coefficients. There is no knotless cubic candidate.
#' For the familiar `bs()` spelling, use [bs()]; `bspline()` remains available
#' for existing formulas. Both wrappers select the same nested models.
#'
#' @param x Numeric predictor.
#' @param df Maximum number of spline coefficients, at least 4.
#' @param knots Optional vector of `df - 3` increasing interior knots.
#' @param boundary_knots Optional two increasing boundary knots. By default
#'   the range of the fitting observations is used.
#' @param center Center retained spline basis columns using the fitted-sample
#'   mean for continuous columns and the minimum for binary columns. Defaults
#'   to `TRUE`, as for `fp()`. Set to `FALSE` to leave the basis uncentered.
#'   Centering changes the intercept, not fitted values.
#' @param catzero Retain a separate exact-zero indicator with the spline.
#' @param spike Select spline and zero-indicator components separately.
#' @return A design matrix with one row per observation and an optional
#'   `catzero` column. `mfp2()` replaces its provisional columns with the
#'   basis built from fitted knots.
#' @export
bspline <- function(x, df = 6L, knots = NULL, boundary_knots = NULL,
                    catzero = FALSE, spike = FALSE, center = TRUE) {
  predictor_name <- deparse(substitute(x))
  if (!is.numeric(x) || !is.null(dim(x)) || is.factor(x)) {
    stop("`bspline()` requires a numeric vector.", call. = FALSE)
  }
  if (length(df) != 1L || !is.numeric(df) || is.na(df) ||
      !is.finite(df) || df < 4L || df > .Machine$integer.max - 3L ||
      df != as.integer(df)) {
    stop("`bspline()` requires an integer `df` of at least 4.", call. = FALSE)
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
  if (!is.null(boundary_knots)) bspline_validate_boundaries(boundary_knots)
  if (!is.null(knots)) bspline_validate_knots(knots, df, boundary_knots)
  finite_x <- x[is.finite(x)]
  fitted_boundaries <- if (is.null(boundary_knots)) {
    if (!is.null(knots)) {
      # Provisional evaluation on one prediction row must also accommodate
      # supplied knots; fitted boundaries are resolved from the training rows.
      limits <- range(c(finite_x, knots))
      margin <- max(diff(limits), 1)
      limits + c(-margin, margin)
    } else if (length(finite_x) && min(finite_x) < max(finite_x)) {
      range(finite_x)
    } else {
      centre <- if (length(finite_x)) finite_x[[1L]] else 0
      c(centre - 1, centre + 1)
    }
  } else boundary_knots
  knot_observations <- if (catzero || spike) finite_x[finite_x > 0] else finite_x
  provisional_knots <- if (is.null(knots)) {
    tryCatch(bspline_default_knots(knot_observations, df, fitted_boundaries),
             error = function(e) {
               seq(fitted_boundaries[[1L]], fitted_boundaries[[2L]],
                   length.out = df - 1L)[-c(1L, df - 1L)]
             })
  } else knots
  basis <- bspline_basis(x, provisional_knots, fitted_boundaries,
                         indicator = catzero || spike)
  attr(basis, "bspline_name") <- predictor_name
  attr(basis, "bspline_df") <- as.integer(df)
  attr(basis, "bspline_supplied_knots") <- knots
  attr(basis, "bspline_boundary_knots") <- boundary_knots
  attr(basis, "bspline_catzero") <- catzero
  attr(basis, "bspline_spike") <- spike
  attr(basis, "bspline_center") <- center
  basis
}

#' Select a cubic B-spline term with `bs()`
#'
#' Use `bs()` inside an `mfp2()` formula to select among exclusion, a linear
#' effect, and nested cubic B-splines. It evaluates the same basis and uses
#' the same default maximum of six coefficients as [bspline()]. Unlike
#' [splines::bs()], this wrapper participates in `mfp2()` model selection.
#' The spline degree is fixed at three and the model supplies its own
#' intercept. If both packages are attached, use `mfp2::bs()` for selection
#' and `splines::bs()` to construct a basis directly.
#'
#' @param x Numeric predictor.
#' @param df Maximum number of spline coefficients; defaults to 6.
#' @param knots Optional vector of `df - 3` increasing interior knots.
#' @param Boundary.knots Optional two increasing boundary knots. Defaults to
#'   the range of fitting observations. This spelling matches `splines::bs()`.
#' @param catzero Retain an exact-zero indicator with the spline.
#' @param spike Select spline and zero-indicator components separately.
#' @param center Center spline basis columns; defaults to `TRUE`.
#' @return A provisional cubic spline basis for use as a formula term;
#'   `mfp2()` resolves and saves the fitted knots for prediction.
#' @export
bs <- function(x, df = 6L, knots = NULL, Boundary.knots = NULL,
               catzero = FALSE, spike = FALSE, center = TRUE) {
  predictor_name <- deparse(substitute(x))
  basis <- bspline(x, df = df, knots = knots,
                   boundary_knots = Boundary.knots,
                   catzero = catzero, spike = spike, center = center)
  attr(basis, "bspline_name") <- predictor_name
  basis
}

#' Validate B-spline boundary knots
#' @param boundaries Pair of boundary values.
#' @return The pair, invisibly; invalid pairs raise an error.
#' @keywords internal
#' @noRd
bspline_validate_boundaries <- function(boundaries) {
  if (!is.numeric(boundaries) || length(boundaries) != 2L ||
      anyNA(boundaries) || any(!is.finite(boundaries)) ||
      boundaries[[1L]] >= boundaries[[2L]]) {
    stop("`bspline()` requires two finite, increasing boundary knots.", call. = FALSE)
  }
  invisible(boundaries)
}

#' Validate B-spline interior knots
#' @param knots Interior knot values.
#' @param df Requested number of spline coefficients.
#' @param boundaries Optional boundary knot pair.
#' @return The knots, invisibly; invalid locations raise an error.
#' @keywords internal
#' @noRd
bspline_validate_knots <- function(knots, df, boundaries = NULL) {
  if (!is.numeric(knots) || length(knots) != df - 3L || anyNA(knots) ||
      any(!is.finite(knots)) || any(diff(knots) <= 0) ||
      (!is.null(boundaries) &&
       (knots[[1L]] <= boundaries[[1L]] ||
        tail(knots, 1L) >= boundaries[[2L]]))) {
    stop("`bspline()` needs `df - 3` increasing interior knots inside the boundaries.",
         call. = FALSE)
  }
  invisible(knots)
}

#' Choose fixed interior knots from fitting observations
#' @param x Finite observations eligible to place interior knots.
#' @param df Requested maximum number of spline coefficients.
#' @param boundaries Boundary knot pair, based on all fitting values.
#' @return Sorted, distinct interior knots (possibly fewer than requested).
#' @keywords internal
#' @noRd
bspline_default_knots <- function(x, df, boundaries) {
  bspline_validate_boundaries(boundaries)
  if (length(x) < 6L) {
    stop("At least six observations are required for automatic B-spline knots.",
         call. = FALSE)
  }
  n <- df - 3L
  probs <- seq_len(n) / (n + 1L)
  knots <- sort(unique(as.numeric(stats::quantile(x, probs, names = FALSE))))
  knots <- knots[knots > boundaries[[1L]] & knots < boundaries[[2L]]]
  if (!length(knots)) {
    stop("Automatic B-spline knots yielded no distinct interior location.",
         call. = FALSE)
  }
  bspline_validate_knots(knots, length(knots) + 3L, boundaries)
  knots
}

#' Specify the centered, nested interior-knot sets
#' @param knots Sorted vector of fixed interior knots.
#' @return A list of knot subsets, smallest first, ending in all knots.
#' @keywords internal
#' @noRd
bspline_knot_sets <- function(knots) {
  count <- length(knots)
  sizes <- seq.int(if (count %% 2L) 1L else 2L, count, by = 2L)
  lapply(sizes, function(size) {
    start <- (count - size) %/% 2L + 1L
    knots[seq.int(start, length.out = size)]
  })
}

#' Find B-spline columns that extend each smaller spline space
#'
#' Knot insertion makes the smaller spline space a subspace of the larger
#' one. On a dense fixed grid, choose independent B-spline columns after
#' the intercept and linear term. This yields one prefix per candidate without
#' treating raw `bs()` column prefixes as smaller-knot spline models.
#' @param knots Full set of fixed interior knots.
#' @param boundaries Fixed boundary knot pair.
#' @return A list of selected knot-set and column pairs.
#' @keywords internal
#' @noRd
bspline_basis_terms <- function(knots, boundaries) {
  knot_positions <- c(boundaries[[1L]], knots, boundaries[[2L]])
  # Sample every knot interval, including intervals too narrow to be hit by
  # a regular grid over the full range.
  grid <- sort(unique(c(
    seq(boundaries[[1L]], boundaries[[2L]], length.out = 101L),
    unlist(lapply(seq_len(length(knot_positions) - 1L), function(i) {
      seq(knot_positions[[i]], knot_positions[[i + 1L]], length.out = 7L)
    }))
  )))
  # Rescaling the linear column keeps rank checks stable when x is far from 0.
  independent <- cbind(1, (grid - boundaries[[1L]]) /
                           (boundaries[[2L]] - boundaries[[1L]]))
  rank <- qr(independent, tol = 1e-10)$rank
  terms <- list()
  for (subset in bspline_knot_sets(knots)) {
    candidate <- splines::bs(grid, knots = subset, degree = 3L,
                             Boundary.knots = boundaries, intercept = FALSE)
    for (column in seq_len(ncol(candidate))) {
      extended <- cbind(independent, candidate[, column])
      next_rank <- qr(extended, tol = 1e-10)$rank
      if (next_rank > rank) {
        independent <- extended
        rank <- next_rank
        terms[[length(terms) + 1L]] <- list(knots = subset, column = column)
      }
    }
    if (length(terms) != length(subset) + 2L) {
      stop("Unable to construct a nested B-spline basis from these knots.",
           call. = FALSE)
    }
  }
  terms
}

#' Evaluate the nested cubic B-spline basis
#' @param x Numeric predictor values.
#' @param knots Fixed full interior-knot vector.
#' @param boundaries Fixed boundary knots.
#' @param indicator Whether to append an exact-zero indicator.
#' @param terms Optional preselected basis-column descriptors.
#' @return A matrix with a linear column, hierarchical spline columns, and
#'   optionally a zero-indicator column.
#' @keywords internal
#' @noRd
bspline_basis <- function(x, knots, boundaries, indicator = FALSE, terms = NULL) {
  bspline_validate_boundaries(boundaries)
  bspline_validate_knots(knots, length(knots) + 3L, boundaries)
  if (is.null(terms)) terms <- bspline_basis_terms(knots, boundaries)
  basis <- matrix(NA_real_, length(x), length(terms) + 1L)
  basis[, 1L] <- x
  # Descriptors sharing a knot set select columns of the same full basis.
  # Evaluate that basis once, then release it before the next knot set.
  pending <- seq_along(terms)
  while (length(pending)) {
    knots_j <- terms[[pending[[1L]]]]$knots
    same <- vapply(terms[pending], function(term) identical(term$knots, knots_j),
                   logical(1L))
    positions <- pending[same]
    columns <- vapply(terms[positions], function(term) term$column, integer(1L))
    evaluated <- splines::bs(x, knots = knots_j, degree = 3L,
                             Boundary.knots = boundaries, intercept = FALSE,
                             warn.outside = FALSE)
    basis[, positions + 1L] <- evaluated[, columns, drop = FALSE]
    pending <- pending[!same]
  }
  colnames(basis) <- c("linear", paste0("nonlinear", seq_along(terms)))
  if (indicator) basis <- cbind(basis, catzero = as.integer(x == 0))
  basis
}

#' Recognize a standalone B-spline formula term
#' @param z Formula term labels.
#' @return A logical vector of matches.
#' @keywords internal
#' @noRd
is_bspline_term <- function(z) grepl("^(mfp2:::?)?(bspline|bs)\\(.*\\)$", z)

#' Label the selected B-spline form
#' @param df Selected coefficient count, or `NA` for exclusion.
#' @param indicator Whether the form includes a zero indicator.
#' @return A human-readable form label.
#' @keywords internal
#' @noRd
bspline_form_label <- function(df, indicator = FALSE) {
  if (is.na(df)) return("excluded")
  if (df == 0L) return("zero indicator only")
  label <- if (df == 1L) "linear" else sprintf("B-spline(df=%d)", df)
  if (indicator) paste0(label, " + zero indicator") else label
}

#' Label a fitted B-spline design column
#' @param term Original predictor name.
#' @param position Position in the full spline basis.
#' @param spec Fitted B-spline specification.
#' @return A human-readable column label.
#' @keywords internal
#' @noRd
bspline_column_label <- function(term, position, spec) {
  if (position == spec$df + 1L) return(sprintf("I(%s = 0)", term))
  if (position == 1L) return(term)
  basis_term <- spec$basis_terms[[position - 1L]]
  sprintf("BS_%d(%s; knots = c(%s); column = %d)", position - 1L, term,
          paste(format(basis_term$knots, digits = 15, trim = TRUE), collapse = ", "),
          basis_term$column)
}
