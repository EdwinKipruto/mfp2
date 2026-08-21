# Variable-ordering helpers for the MFP backfitting algorithm
#
# fit_mfp() supplies the resolved linear reference representation:
#
#   ordinary term      -> x
#   zero term          -> x+ (nonpositive values already recoded to zero)
#   catzero/SAZ term   -> x+ plus a prebuilt structural-zero indicator
#
# This file does not recode zero values or construct catzero indicators. It only
# assembles the supplied reusable blocks, fits the full reference once, and uses
# that same fit for leave-one-conceptual-term-out significance ordering.


# -----------------------------------------------------------------------------
# order_variables() -----------------------------------------------------------
# -----------------------------------------------------------------------------

#' Fit the Full Linear Reference Model and Determine Predictor Visiting Order
#'
#' Fits the resolved full linear reference model and determines the order in
#' which conceptual predictors are visited by the MFP backfitting algorithm.
#'
#' @details
#' The caller supplies the continuous reference matrix in its final starting
#' form. Thus ordinary terms remain unchanged, while terms with \code{zero}
#' handling have already had nonpositive values recoded to zero. Optional
#' \code{catzero_blocks} contain the already-built binary structural-zero
#' indicators. These indicators are appended to the corresponding conceptual
#' terms only for fitting the reference model; they are not recalculated here.
#'
#' The full reference model is fitted unconditionally because its null and
#' fitted-model statistics are retained by \code{fit_mfp()}. For GLMs, the null
#' deviance reported by the full fit is retained directly; no separate
#' intercept-only likelihood fit is required. For Cox models the reference
#' remains intercept-free.
#'
#' Significance-based ranking is performed only when more than one conceptual
#' predictor is present and \code{xorder} is \code{"ascending"} or
#' \code{"descending"}. One reduced model is fitted per conceptual predictor by
#' removing that predictor's complete reference block. Consequently a
#' \code{catzero} or retained spike-at-zero term is tested jointly through its
#' positive-part continuous column and binary indicator. Test degrees of freedom
#' are the fitted rank difference between the full and reduced models.
#'
#' \describe{
#'   \item{\code{"ascending"}}{Smallest p-value first.}
#'   \item{\code{"descending"}}{Largest p-value first.}
#'   \item{\code{"original"}}{Preserve the supplied conceptual-term order and
#'     fit no reduced ordering models.}
#' }
#'
#' @section Deviance convention:
#' Deviance is family-specific and is taken directly from the full reference
#' fit returned by \code{fit_model()}. For GLMs, \code{null_deviance} and
#' \code{linear_deviance} are the \code{null.deviance} and \code{deviance}
#' values returned by the GLM fit. For Cox models, they are minus twice the null
#' and fitted partial log-likelihoods.
#'
#' @param xorder Character scalar controlling the conceptual-predictor visiting
#'   order. Supported values are \code{"ascending"}, \code{"descending"}, and
#'   \code{"original"}.
#' @param x Numeric matrix containing the continuous reference columns, excluding
#'   the intercept. Structural-zero recoding must already have been applied by
#'   the caller.
#' @param term_to_columns Named list mapping each conceptual term to its columns
#'   in \code{x}. Grouped fixed terms may map to multiple columns.
#' @param catzero_blocks Optional named list, indexed by conceptual term, of
#'   prebuilt one-column binary structural-zero indicator matrices. NULL entries
#'   are ignored. The matrices are reused as supplied and are not recomputed.
#' @param y Response used to fit the models.
#' @param family GLM family object used by \code{fit_model()}, or character
#'   \code{"cox"} for a Cox proportional-hazards model.
#' @param family_string Normalized character family name.
#' @param weights Optional observation weights passed to \code{fit_model()}.
#' @param offset Optional linear-predictor offset passed to \code{fit_model()}.
#' @param strata Optional Cox stratification object. Ignored for GLMs.
#' @param method Cox tie-handling method. Ignored for GLMs.
#' @param control Model-fitting control object.
#' @param nocenter Cox centering-suppression argument. Ignored for GLMs.
#' @param fitter GLM fitting backend; ignored for Cox models.
#'
#' @return A list containing \code{variables_ordered}, \code{null_deviance},
#'   \code{linear_deviance}, \code{linear_logl}, \code{linear_df}, and
#'   \code{null_logl}.
#'
#' @keywords internal
#' @noRd
order_variables <- function(xorder = "ascending",
                            x,
                            y,
                            family,
                            family_string,
                            weights = NULL,
                            offset = NULL,
                            strata = NULL,
                            method = NULL,
                            control = NULL,
                            nocenter = NULL,
                            term_to_columns = NULL,
                            catzero_blocks = NULL,
                            fitter = "base") {
  if (is.null(term_to_columns)) {
    term_to_columns <- stats::setNames(as.list(colnames(x)), colnames(x))
  }

  predictor_names <- names(term_to_columns)
  n_predictors <- length(term_to_columns)

  # Cox fast fits require integer stratum identifiers. Convert once for the
  # complete reference/ordering operation rather than once per reduced model.
  strata_ordering <- if (identical(family_string, "cox") &&
                         !is.null(strata) && !is.integer(strata)) {
    as.integer(strata)
  } else {
    strata
  }

  use_glm_intercept_template <- !identical(family_string, "cox")

  # Fast path: when no catzero/SAZ indicator blocks are present, the supplied x
  # is already the complete reference design. No structural-zero assembly or
  # extra predictor matrix is created. Plain zero terms need no special branch
  # here because fit_mfp() has already recoded their continuous column in x.
  if (is.null(catzero_blocks)) {
    x_fit <- if (use_glm_intercept_template) {
      assemble_design_matrix(
        blocks = list(x),
        nobs = NROW(x),
        intercept = TRUE
      )
    } else {
      x
    }
    reference_term_to_columns <- term_to_columns
  } else {
    reference_design <- assemble_linear_reference_design(
      x = x,
      term_to_columns = term_to_columns,
      catzero_blocks = catzero_blocks,
      intercept = use_glm_intercept_template
    )
    x_fit <- reference_design$x
    reference_term_to_columns <- reference_design$term_to_columns
  }

  # Fit the resolved full reference once. The same likelihood/df are reused by
  # every leave-one-term-out comparison, so the full model is never refitted for
  # significance ordering.
  full_reference <- fit_full_linear_reference(
    x = x_fit,
    y = y,
    family = family,
    family_string = family_string,
    fitter = fitter,
    weights = weights,
    offset = offset,
    strata = strata_ordering,
    method = method,
    control = control,
    nocenter = nocenter,
    x_has_intercept = use_glm_intercept_template
  )

  rank_predictors <- n_predictors > 1L && !identical(xorder, "original")

  variables_ordered <- if (rank_predictors) {
    order_variables_by_significance(
      xorder = xorder,
      x = x_fit,
      term_to_columns = reference_term_to_columns,
      y = y,
      family = family,
      family_string = family_string,
      fitter = fitter,
      weights = weights,
      offset = offset,
      strata = strata_ordering,
      method = method,
      control = control,
      nocenter = nocenter,
      full_reference = full_reference,
      x_has_intercept = use_glm_intercept_template
    )
  } else {
    predictor_names
  }

  list(
    variables_ordered = variables_ordered,
    null_deviance = full_reference$null_deviance,
    linear_deviance = full_reference$model_deviance,
    linear_logl = full_reference$logl,
    linear_df = full_reference$df,
    null_logl = if (identical(family_string, "cox")) {
      full_reference$null_logl
    } else {
      NA_real_
    }
  )
}


# -----------------------------------------------------------------------------
# assemble_linear_reference_design() -----------------------------------------
# -----------------------------------------------------------------------------

#' Assemble the Linear Reference Design from Reusable Blocks
#'
#' Appends precomputed catzero/SAZ indicator blocks to the continuous reference
#' matrix and extends the conceptual term-to-column mapping. This helper does not
#' recode x and does not calculate any indicator values.
#'
#' @inheritParams order_variables
#' @param intercept Logical; prepend a GLM intercept when TRUE.
#'
#' @return A list containing the assembled fit matrix \code{x} and extended
#'   \code{term_to_columns} mapping.
#' @keywords internal
#' @noRd
assemble_linear_reference_design <- function(x,
                                             term_to_columns,
                                             catzero_blocks,
                                             intercept = FALSE) {
  if (is.null(names(catzero_blocks))) {
    stop("Internal error: catzero reference blocks must be named.", call. = FALSE)
  }

  active <- names(catzero_blocks)[vapply(
    catzero_blocks,
    function(block) !is.null(block) && NCOL(block) > 0L,
    logical(1L)
  )]

  # A non-NULL all-empty list is accepted for defensive internal calls, but the
  # normal fit_mfp() path passes NULL and therefore bypasses this helper entirely.
  if (length(active) == 0L) {
    x_fit <- if (isTRUE(intercept)) {
      assemble_design_matrix(list(x), NROW(x), intercept = TRUE)
    } else {
      x
    }
    return(list(x = x_fit, term_to_columns = term_to_columns))
  }

  missing_terms <- setdiff(active, names(term_to_columns))
  if (length(missing_terms) > 0L) {
    stop(
      "Internal error: catzero reference blocks do not match conceptual terms.",
      call. = FALSE
    )
  }

  for (term in active) {
    block <- catzero_blocks[[term]]
    if (NCOL(block) != 1L || NROW(block) != NROW(x)) {
      stop(
        sprintf(
          "Internal error: catzero reference block '%s' must have one column and match x rows.",
          term
        ),
        call. = FALSE
      )
    }
  }

  # Indicator blocks are stored for backfitting with the generic column name
  # "catzero". In the assembled linear-reference design, name each appended
  # indicator directly from the conceptual term that generated it, following
  # mfp2's existing "<term>_bin" convention. The stored backfitting matrices
  # themselves are left unchanged.
  indicator_names <- paste0(active, "_bin")

  # Do not silently repair a collision with make.unique(). A generated *_bin
  # name identifies a specific structural-zero indicator, so an existing design
  # column with the same name would make the conceptual term mapping ambiguous.
  # Check before allocating the augmented reference matrix.
  if (anyDuplicated(indicator_names)) {
    stop(
      "Internal error: duplicate catzero indicator names were generated.",
      call. = FALSE
    )
  }

  name_conflicts <- indicator_names %in% colnames(x)
  if (any(name_conflicts)) {
    stop(
      sprintf(
        "Generated catzero indicator name(s) already exist in the model matrix: %s.",
        paste(indicator_names[name_conflicts], collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # assemble_design_matrix() allocates the destination once and copies x and the
  # existing indicator blocks directly into it. No cbind(x, indicators)
  # intermediate is created.
  x_fit <- assemble_design_matrix(
    blocks = c(list(x), unname(catzero_blocks[active])),
    nobs = NROW(x),
    intercept = intercept
  )

  # The active indicator blocks are appended in the same order as `active`, and
  # each block has exactly one column (validated above). Rename only those
  # assembled columns; catzero_blocks continues to use the generic "catzero"
  # column name required by backfitting.
  indicator_positions <- NCOL(x_fit) - length(active) + seq_along(active)
  colnames(x_fit)[indicator_positions] <- indicator_names

  reference_term_to_columns <- term_to_columns
  for (index in seq_along(active)) {
    term <- active[[index]]
    reference_term_to_columns[[term]] <- c(
      reference_term_to_columns[[term]],
      indicator_names[[index]]
    )
  }

  list(
    x = x_fit,
    term_to_columns = reference_term_to_columns
  )
}


# -----------------------------------------------------------------------------
# fit_full_linear_reference() -------------------------------------------------
# -----------------------------------------------------------------------------

#' Fit the Full Linear Reference Model
#'
#' Fits the resolved full linear reference design supplied by \code{order_variables()}.
#' The design may contain positive-part continuous columns and structural-zero
#' indicators in addition to ordinary linear columns. The fitted model supplies
#' the null/full-reference deviances and the likelihood/df used by the
#' leave-one-conceptual-term-out ordering tests.
#'
#' @inheritParams order_variables
#'
#' @return A lightweight model-fit wrapper returned by \code{fit_model()},
#'   including \code{logl}, \code{df}, \code{rank}, \code{null_deviance},
#'   \code{model_deviance}, coefficients, and degrees of freedom. Cox
#'   fits additionally include \code{null_logl}. The underlying fast-fit object is
#'   not retained.
#'
#' @keywords internal
#' @noRd
fit_full_linear_reference <- function(x,
                                      y,
                                      family,
                                      family_string,
                                      weights,
                                      offset,
                                      strata,
                                      method,
                                      control,
                                      nocenter,
                                      fitter = "base",
                                      x_has_intercept = FALSE) {
  fit_model(
    x = x,
    y = y,
    family = family,
    family_string = family_string,
    fitter = fitter,
    weights = weights,
    offset = offset,
    method = method,
    strata = strata,
    control = control,
    rownames = rownames(x),
    nocenter = nocenter,
    fast = TRUE,
    calculate_fit_statistics = TRUE,
    keep_fit = FALSE,
    x_has_intercept = x_has_intercept
  )
}


# -----------------------------------------------------------------------------
# order_variables_by_significance() -------------------------------------------
# -----------------------------------------------------------------------------

#' Order Predictors by Leave-One-Out Significance
#'
#' Orders predictors using likelihood-ratio tests comparing the full linear
#' reference model with models obtained by removing one predictor at a time.
#'
#' @details
#' This helper does not fit the full linear model. The caller supplies
#' \code{full_reference}, ensuring that the invariant full-model fit is computed
#' only once and reused for every reduced-model comparison.
#'
#' A predictor may correspond to one or more design-matrix columns. All columns
#' belonging to one conceptual term are removed together, and the test degrees
#' of freedom equal the fitted rank difference between the full and reduced
#' models. If that difference is not positive, or if a valid likelihood-ratio
#' statistic cannot be formed, the predictor receives
#' \code{NA} as its ordering p-value and is placed after predictors with valid
#' p-values. No predictor is dropped from the returned order.
#'
#' Ties retain the original conceptual-term order.
#'
#' @inheritParams order_variables
#' @param full_reference Model-fit wrapper returned by
#'   \code{fit_full_linear_reference()} for the resolved full reference model.
#'
#' @return Character vector containing all predictor names in significance-based
#'   visiting order.
#'
#' @keywords internal
#' @noRd
order_variables_by_significance <- function(xorder,
                                            x,
                                            y,
                                            family,
                                            family_string,
                                            weights,
                                            offset,
                                            strata,
                                            method,
                                            control,
                                            nocenter,
                                            full_reference,
                                            term_to_columns = NULL,
                                            fitter = "base",
                                            x_has_intercept = FALSE) {
  if (is.null(term_to_columns)) {
    predictor_columns <- if (isTRUE(x_has_intercept)) {
      colnames(x)[-1L]
    } else {
      colnames(x)
    }
    term_to_columns <- stats::setNames(
      as.list(predictor_columns),
      predictor_columns
    )
  }

  predictor_names <- names(term_to_columns)
  n_predictors <- length(term_to_columns)
  has_mapped_terms <- any(vapply(
    predictor_names,
    function(term) {
      cols <- term_to_columns[[term]]
      length(cols) != 1L || !identical(cols[[1L]], term)
    },
    logical(1L)
  ))

  # Cache matrix metadata and full-reference scalars once. These values are
  # invariant across all leave-one-term-out fits and may otherwise be looked up
  # repeatedly inside the hot loop.
  x_colnames <- colnames(x)
  x_rownames <- rownames(x)
  full_df <- full_reference$df
  full_logl <- full_reference$logl
  descending <- identical(xorder, "descending")

  # Precompute integer column positions to drop. Integer subsetting avoids
  # rebuilding character keep-vectors with setdiff() on every iteration and
  # avoids repeated name-to-position matching during matrix subsetting.
  if (has_mapped_terms) {
    drop_indices <- lapply(term_to_columns, function(columns) {
      indices <- match(columns, x_colnames)

      if (anyNA(indices)) {
        stop(
          "Internal error: term-to-column mapping does not match design matrix.",
          call. = FALSE
        )
      }

      indices
    })
  } else {
    intercept_offset <- as.integer(isTRUE(x_has_intercept))
    drop_indices <- lapply(
      seq_len(n_predictors),
      function(index) index + intercept_offset
    )
  }

  # Initialize with NA rather than zero. A failed or non-identifiable comparison
  # must not be interpreted as overwhelming evidence against the predictor.
  p_values <- stats::setNames(
    rep(NA_real_, n_predictors),
    predictor_names
  )

  for (predictor_index in seq_len(n_predictors)) {
    # Materialize only the reduced matrix needed for the current test. Matrix
    # subsetting necessarily allocates here, so keep this large temporary alive
    # for the shortest possible interval to reduce peak memory and GC pressure.
    reduced_x <- x[
      ,
      -drop_indices[[predictor_index]],
      drop = FALSE
    ]

    reduced_fit <- fit_model(
      x = reduced_x,
      y = y,
      family = family,
      family_string = family_string,
      fitter = fitter,
      weights = weights,
      offset = offset,
      method = method,
      strata = strata,
      control = control,
      rownames = x_rownames,
      nocenter = nocenter,
      fast = TRUE,
      x_has_intercept = x_has_intercept
    )

    # fast = TRUE does not retain the supplied design matrix in the returned
    # wrapper. Drop our reference immediately so the previous reduced matrix is
    # collectible before the next large subset is allocated. Do not force gc()
    # here; allowing R to collect naturally avoids a full GC on every predictor.
    reduced_x <- NULL

    # Copy the two scalar results that are needed below, then release the fit
    # wrapper as well. This shortens the lifetime of any backend temporaries
    # reachable from it without changing the likelihood-ratio calculation.
    reduced_df <- reduced_fit$df
    reduced_logl <- reduced_fit$logl
    reduced_fit <- NULL

    # Use the fitted rank contribution of the omitted conceptual term. Raw
    # design-block width can exceed this difference when one or more grouped
    # columns are aliased or otherwise non-estimable in the fitted model.
    lrt_df <- full_df - reduced_df
    lrt_statistic <- 2 * (full_logl - reduced_logl)

    # A valid nested-model likelihood-ratio test requires a positive df
    # difference and finite likelihoods. Small negative statistics can occur
    # from numerical rounding, so truncate such values to zero.
    if (is.finite(lrt_df) && lrt_df > 0L && is.finite(lrt_statistic)) {
      p_values[predictor_index] <- stats::pchisq(
        q = max(0, lrt_statistic),
        df = lrt_df,
        lower.tail = FALSE
      )
    }
  }

  # Radix ordering is stable, so tied p-values retain the original conceptual-
  # term order without allocating an additional seq_along() tie-break vector.
  # Direct decreasing order also avoids allocating a negated copy of p_values.
  ordering_index <- order(
    p_values,
    decreasing = descending,
    na.last = TRUE,
    method = "radix"
  )

  predictor_names[ordering_index]
}