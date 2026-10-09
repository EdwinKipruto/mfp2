#' Estimate a Common GEE Dispersion for Model Selection
#'
#' Constructs a full, maximum-degree reference mean model and returns its
#' estimated GEE scale. Called before MFP backfitting when a GEE AIC or BIC fit
#' has no user-supplied `gee_family(reference_dispersion = ...)`.
#'
#' @details
#' Each conceptual predictor is included. Terms mapped to multiple raw columns
#' (grouped terms), and terms with `df = 1`, are assigned power 1, even when 1
#' is not in their FP candidate set. For other terms, FP degree is `df / 2`. The
#' function generates the permitted power combinations with
#' `generate_powers_fp()` or, for ACD terms, `generate_powers_acd()`. It chooses
#' the initial combination without fitting candidates. To avoid transformations
#' that diverge near zero, it prefers powers at least 1, followed by positive
#' powers below 1, zero (a logarithm), and negative powers. Combinations are
#' ranked by their most sensitive power, then total sensitivity. For FP1, a
#' stable nonlinear power is preferred to power 1; the target is power 2.
#' For higher degrees, distinct powers are preferred to repeated powers, then
#' combinations are ordered by squared distance from (1, 2) for FP2 or
#' (1, 2, 3, ...) for higher degrees. Ties retain the
#' generated candidate order. ACD degree-1 rows have an inactive first power
#' (`NA`), so only their second power is ranked. If a term has only one allowed
#' FP power, its FP2 uses the existing repeated-power transformation.
#'
#' The function first tries the highest-ranked combination for every term. If
#' that attempt fails during transformation, the design check, fitting, or
#' scale validation, it makes simple fallback attempts in the original term
#' order. Each fallback changes just one term to its second-ranked allowed
#' power combination, leaving every other term at its original power. It stops
#' at the first valid scale. Thus at most one initial fit and one fallback per
#' eligible term are attempted; power combinations are not fitted to optimize
#' the reference model.
#'
#' For each attempt, `expand_term_metadata_to_columns()` maps the powers and
#' resolved term settings to design columns. `transform_matrix()` constructs
#' the full reference design using the existing centering, ACD, and
#' structural-zero transformations. Spike-at-zero terms include both their
#' continuous and binary components. The prepared GEE family supplies the
#' fitted response, effective weights (including grouped-binomial totals),
#' offset, cluster IDs, and working correlation. The design must have full
#' column rank after adding an intercept and at least one residual degree of
#' freedom. A design passing these checks receives one lightweight `fit_model()`
#' GEE fit; its fitted `geese.fit()` scale (`gamma[1]`) must be positive and
#' finite. Only the first successful fit supplies the selection dispersion.
#'
#' The returned scale is computed once and reused for every QICu and QBIC
#' comparison at every selection step and backfitting cycle. It does not fix
#' geepack's scale during candidate or final fits, and does not alter robust
#' Wald p-value selection. If no attempt succeeds, the function reports the
#' number of attempts and the last failure, then stops with an error. It does
#' not drop terms, lower FP degree, choose powers outside the permitted set,
#' or substitute a default scale. The user can instead supply a positive
#' `reference_dispersion` in `gee_family()` to skip the reference fit.
#'
#' @param x Shifted and scaled training design matrix.
#' @param y Training response.
#' @param family Prepared GEE family, containing the response, effective weights,
#'   clusters, offset, and working correlation.
#' @param df Named maximum MFP degrees of freedom for conceptual terms.
#' @param powers Named list of permitted FP powers for each conceptual term.
#' @param term_to_columns Named mapping from conceptual terms to raw columns.
#' @param center,acdx,zero,catzero,spike Resolved transformation settings.
#' @param acd_parameter Previously estimated ACD parameters.
#' @param center_method Optional centering methods for raw design columns.
#' @param control GEE fit controls.
#' @return A list containing `dispersion` (positive estimated GEE scale) and
#'   `powers` (the successful reference powers for every conceptual term).
#' @keywords internal
estimate_reference_dispersion <- function(x, y, family, df, powers,
                                          term_to_columns, center, acdx, zero,
                                          catzero, spike, acd_parameter,
                                          center_method, control) {
  terms <- names(term_to_columns)
  candidate_powers <- setNames(vector("list", length(terms)), terms)

  for (term in terms) {
    if (length(term_to_columns[[term]]) > 1L || df[[term]] == 1L) {
      candidate_powers[[term]] <- matrix(1, nrow = 1L)
      next
    }
    degree <- as.integer(df[[term]] / 2L)
    combinations <- if (isTRUE(acdx[[term]])) {
      generate_powers_acd(degree, powers[[term]])
    } else {
      generate_powers_fp(degree, powers[[term]])
    }
    active <- if (isTRUE(acdx[[term]]) && degree == 1L) 2L else
      seq_len(ncol(combinations))
    active_powers <- combinations[, active, drop = FALSE]
    targets <- if (length(active) == 1L) 2 else seq_along(active)
    distance <- rowSums(sweep(active_powers, 2L, targets, "-")^2)
    sensitivity <- (active_powers < 1) + (active_powers <= 0) +
      (active_powers < 0)
    worst_sensitivity <- apply(sensitivity, 1L, max)
    total_sensitivity <- rowSums(sensitivity)
    repeated <- if (length(active) > 1L) {
      apply(active_powers, 1L,
            function(p) anyDuplicated(p) > 0L)
    } else {
      rep.int(FALSE, nrow(combinations))
    }
    linear_fp1 <- if (length(active) == 1L) {
      active_powers[, 1L] == 1
    } else {
      rep.int(FALSE, nrow(combinations))
    }
    candidate_powers[[term]] <- combinations[
      order(worst_sensitivity, total_sensitivity, linear_fp1, repeated, distance,
            seq_len(nrow(combinations))), , drop = FALSE
    ]
  }

  try_reference <- function(indices) {
    reference_powers <- setNames(lapply(seq_along(terms), function(i) {
      unname(candidate_powers[[i]][indices[i], ])
    }), terms)
    settings <- expand_term_metadata_to_columns(
      term_to_columns = term_to_columns, powers = reference_powers,
      center = center, acdx = acdx, zero = zero, catzero = catzero,
      spike = spike, spike_decision = setNames(
        rep.int(saz_decision_codes[["cont_binary"]], length(terms)), terms
      ), acd_parameter = acd_parameter
    )
    transformed <- transform_matrix(
      x = x, power_list = settings$powers, center = settings$center,
      acdx = settings$acdx, acd_parameter_list = settings$acd_parameter,
      zero = settings$zero, catzero = settings$catzero,
      spike = settings$spike, spike_decision = settings$spike_decision,
      reset_zero = FALSE, center_method = center_method
    )
    design <- transformed$x_transformed
    if (is.null(design) || nrow(design) <= ncol(design) + 1L ||
        qr(cbind(1, design))$rank != ncol(design) + 1L) {
      stop("The full GEE reference design is too large or rank deficient.")
    }
    fitted <- fit_model(
      x = design, y = y, family = family, family_string = "gee",
      control = control, fast = TRUE, keep_fit = FALSE,
      gee_selection_criterion = "aic"
    )
    dispersion <- fitted$gee_scale
    if (!is.numeric(dispersion) || length(dispersion) != 1L ||
        !is.finite(dispersion) || dispersion <= 0) {
      stop("The full GEE reference model did not yield a positive finite dispersion.")
    }
    list(dispersion = dispersion, powers = reference_powers)
  }

  original <- rep.int(1L, length(terms))
  result <- tryCatch(try_reference(original), error = function(e) e)
  if (!inherits(result, "error")) return(result)
  last_error <- conditionMessage(result)
  attempts <- 1L
  for (i in seq_along(terms)) {
    if (nrow(candidate_powers[[i]]) < 2L) next
    alternative <- original
    alternative[i] <- 2L
    attempts <- attempts + 1L
    result <- tryCatch(try_reference(alternative), error = function(e) e)
    if (!inherits(result, "error")) return(result)
    last_error <- conditionMessage(result)
  }
  stop("The full GEE reference model failed after ", attempts,
       " power configuration(s). Last failure: ", last_error,
       " Supply `reference_dispersion` in `gee_family()` or simplify the model.",
       call. = FALSE)
}
