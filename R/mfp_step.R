#' Decide Whether an MFP Comparison Supports Simplification
#'
#' MFP removes a variable or simplifies a functional form only when the
#' comparison p-value is a valid finite probability that is strictly greater
#' than the relevant selection threshold. A missing/non-finite p-value is not
#' evidence in favour of the simpler model, so the current (more complex)
#' representation is retained. This also prevents NA p-values from propagating
#' into model indices during selection.
#'
#' @param pvalue Numeric scalar p-value from an MFP comparison.
#' @param threshold Numeric scalar selection threshold.
#'
#' @return Logical scalar.
#' @keywords internal
#' @noRd
mfp_pvalue_exceeds <- function(pvalue, threshold) {
  is.numeric(pvalue) &&
    length(pvalue) == 1L &&
    !is.na(pvalue) &&
    is.finite(pvalue) &&
    pvalue >= 0 && pvalue <= 1 &&
    is.numeric(threshold) &&
    length(threshold) == 1L &&
    !is.na(threshold) &&
    is.finite(threshold) &&
    pvalue > threshold
}

#' Select the best FP function for one variable
#'
#' See \code{mfp2()} for the notation used here and \code{fit_mfp()} for the
#' overall fitting procedure.
#'
#' @param x Numeric design matrix with one row per observation and expanded
#'   columns, excluding the intercept. Ordinary FP variables have been shifted
#'   and scaled; active ACD variables are shifted but use scale 1.
#' @param y Response vector or \code{Surv} object.
#' @param xi Name of the term being assessed. A continuous term may take an FP
#'   form; a term represented by several columns is assessed as one fixed
#'   linear group.
#' @param term_to_columns Named list identifying the columns belonging to each
#'   term. Grouped terms are tested jointly.
#' @param weights Observation weights.
#' @param offset Observation offsets.
#' @param df Numeric scalar giving the maximum selection degrees of freedom
#'   for \code{xi}. A value of 1 restricts the term to a linear form.
#' @param powers_current Named list containing the current form of each term.
#'   The forms of other terms are used when assessing \code{xi}.
#' @param family Resolved model family supplied by the fitting procedure.
#' @param family_string Name of the selected family, such as
#'   \code{"gaussian"}.
#' @param criterion Selection criterion: \code{"pvalue"}, \code{"aic"}, or
#'   \code{"bic"}.
#' @param select Significance level for retaining \code{xi} when using
#'   \code{"pvalue"}. A grouped term is tested as a whole. A value of 1
#'   forces the term to remain.
#' @param alpha Significance level for choosing among forms when using
#'   \code{"pvalue"}. For an eligible spike-at-zero term, it also applies
#'   to the Stage-2 tests of its separate components.
#' @param keep Names of terms that must remain in the model.
#' @param powers Named list of permitted FP powers for each variable.
#' @param method Method for handling tied event times in Cox regression.
#' @param strata Factor identifying combinations of stratification variables.
#' @param nocenter Values exempted from centering in Cox models; see
#'   \code{survival::coxph()}.
#' @param acdx Named logical vector indicating which terms use ACD.
#'   Terms represented by several columns must have value \code{FALSE}.
#' @param ftest Whether to use F-tests for Gaussian models.
#' @param control Model-fitting controls.
#' @param rownames Row names used when fitting Cox models.
#' @param catzero Named list of exact-zero indicators. Each entry is
#'   \code{NULL} or a one-column matrix indicating which values equal zero.
#' @param zero Named logical vector indicating exact-zero handling for
#'   eligible nonnegative terms.
#' @param spike Named logical vector indicating spike-at-zero handling for
#'   eligible nonnegative terms. Zero and positive values form the two parts.
#' @param acd_parameter Named list of ACD parameters, with one entry per term.
#'   An entry is \code{NULL} when ACD is not active for that term.
#' @param spike_decision Named vector describing the selected form for each
#'   spike-at-zero term: \code{1} retains both the positive-value effect and
#'   the zero indicator; \code{2} retains the positive-value effect only; and
#'   \code{3} retains the zero indicator only.
#' @param prev_adj_params Information from earlier steps that can be reused
#'   when assessing other terms.
#' @param transform_cache Previously prepared forms of individual terms that
#'   can be reused during the current fit.
#' @param force_max_fp Named logical vector indicating terms that must use
#'   their fullest permitted form. A term restricted to \code{df = 1} remains
#'   linear. For other terms, the permitted powers are searched within the
#'   maximum form, without comparing that form with simpler alternatives.
#'   An eligible spike-at-zero term retains both its positive-value effect
#'   and zero indicator.
#' @param has_offset Whether an offset was supplied before missing offsets
#'   were replaced with zeros internally.
#' @param n_obs Size used for model comparisons. This is usually the number
#'   of observations. It is the number of observed events for Cox models,
#'   selected events for Fine--Gray models, effective observations for
#'   multinomial models, and independent groups for GEE models.
#' @param fitter Model-fitting method supplied by the calling procedure.
#' @param precomputed_adj Previously prepared information about the other
#'   terms, used by lower-level fitting functions. This is not an argument
#'   to \code{find_best_fp_step()} itself.
#' @param focal_basis_cache Previously prepared values for the term being
#'   assessed, used by lower-level searches. They must correspond to the
#'   current observations, zero handling, and permitted powers.
#' @param verbose Whether to print details of the selection step.
#'
#' @details
#' Each term is assessed while accounting for the current forms of the other
#' terms. Terms are reconsidered in later cycles, so a term removed in one
#' cycle may be selected again. Each term starts with a linear form, but its
#' configured maximum \code{df} is available from the first cycle.
#'
#' With \code{criterion = "pvalue"}, the function selection procedure tests
#' whether to retain a term and whether its form can be simplified. A term
#' restricted to \code{df = 1} is assessed by comparing its linear form with
#' its exclusion. The columns of a grouped term are assessed together.
#'
#' With \code{"aic"} or \code{"bic"}, eligible forms are compared using the
#' chosen criterion. Eligible spike-at-zero terms are assessed by comparing
#' their available combinations of positive-value and zero-indicator effects.
#' Under \code{"pvalue"}, a selected, non-forced spike-at-zero term instead
#' has a second stage of tests for those effects.
#'
#' A term marked by \code{force_max_fp} keeps its fullest permitted form under
#' any criterion. The permitted powers within that form are still searched.
#'
#' @return A list containing \code{power_best}, \code{searched_fp},
#'   \code{spike_decision}, \code{current_adj_params}, and
#'   \code{transform_cache}. A missing value in \code{power_best} may mean
#'   that the term was removed or that only its zero indicator was retained.
#'   Consult \code{spike_decision} to distinguish these cases.
#'
#' @references
#' Royston, P. and Sauerbrei, W. (2008). \emph{Multivariable Model-Building:
#' A Pragmatic Approach to Regression Analysis Based on Fractional
#' Polynomials for Modelling Continuous Variables}. John Wiley & Sons.
#'
#' Royston, P. and Sauerbrei, W. (2016). mfpa: Extension of mfp using the ACD
#' covariate transformation for enhanced parametric multivariable modeling.
#' \emph{The Stata Journal}, 16(1), 72--87.
#' @keywords internal
#' @noRd
find_best_fp_step <- function(x,
                              y,
                              xi,
                              weights,
                              offset,
                              df,
                              powers_current,
                              family,
                              family_string,
                              criterion,
                              select,
                              alpha,
                              keep,
                              powers,
                              method,
                              strata,
                              nocenter,
                              acdx,
                              ftest,
                              control,
                              rownames,
                              zero,
                              catzero,# a named list of binary variables
                              spike,
                              spike_decision,
                              acd_parameter,
                              prev_adj_params,
                              transform_cache = NULL,
                              force_max_fp,
                              has_offset,
                              n_obs,
                              verbose,
                              term_to_columns,
                              fitter = "base",
                              step_memo = NULL) {

  # `fit_mfp()` constructs the complete conceptual-term lookup once and the
  # backfitting cycle passes it explicitly through every selection routine.
  # find_best_fp_step() dispatches xi's model selection to the path implied by
  # force_max_fp, spike-at-zero eligibility, and the selection criterion.
  # P-value SAZ selection keeps its two closed-testing stages. AIC/BIC SAZ
  # selection instead compares the complete candidate family jointly, because
  # there is no closed-test error-control argument requiring two stages.

  degree <- as.numeric(df / 2)

  criterion_lower <- tolower(criterion)
  force_max_active <- isTRUE(force_max_fp[[xi]])
  saz_active <- isTRUE(spike[[xi]])

  # Shared arguments for every path below. The conceptual-term lookup is
  # passed to every selector because an adjustment term can be a grouped term
  # spanning multiple raw design-matrix columns.
  # Pass the term lookup to every selector, not only to select_linear(). Even
  # when xi is a continuous FP term, its adjustment set can contain grouped
  # categorical terms whose raw columns must be assembled as one block.
  selector_args <- list(
    x = x, xi = xi, keep = keep, degree = degree, acdx = acdx,
    y = y, family = family, family_string = family_string,
    fitter = fitter,
    weights = weights, offset = offset, force_max_fp = force_max_fp,
    powers_current = powers_current, powers = powers,
    criterion = criterion, ftest = ftest, select = select, alpha = alpha,
    method = method, strata = strata, nocenter = nocenter, n_obs = n_obs,
    control = control, rownames = rownames, zero = zero, catzero = catzero,
    spike = spike, spike_decision = spike_decision, has_offset = has_offset,
    acd_parameter = acd_parameter, prev_adj_params = prev_adj_params,
    transform_cache = transform_cache,
    term_to_columns = term_to_columns,
    calculate_gaussian_deviance = isTRUE(ftest)
  )

  # Path 1: force_max_fp fixes the final form before any selection criterion.
  # For df = 1 the maximum form is linear; adding xi to the selector's local
  # keep set forces that row. For df > 1, select_force_max_fp() fits only the
  # requested highest FP form. An eligible SAZ variable retains its binary
  # component as well, so neither the joint-IC path nor SAZ Stage 2 may undo
  # the explicit force request.
  if (force_max_active) {
    if (df == 1) {
      forced_args <- selector_args
      forced_args$keep <- unique(c(keep, xi))
      fit1 <- do.call(select_linear, forced_args)
      # `keep` reports the user's keep setting, not the local implementation
      # device used above to force the only available continuous form.
      fit1$keep <- xi %in% keep
    } else {
      fit1 <- do.call(select_force_max_fp, selector_args)
    }

    if (verbose) {
      print_mfp_step(xi = xi, criterion = criterion, fit = fit1)
    }

    power_best <- normalize_selected_step_powers(fit1, xi, acdx)
    searched_fp <- selected_step_is_searched_fp(fit1)
    if (saz_active) {
      spike_decision[[xi]] <- saz_decision_codes[["cont_binary"]]
    }

    return(list(
      power_best = power_best,
      searched_fp = searched_fp,
      spike_decision = spike_decision,
      current_adj_params = fit1$current_adj_params,
      transform_cache = fit1$transform_cache
    ))
  }

  # Path 2: an eligible SAZ term under AIC/BIC is one joint model-selection
  # problem. The positive-only and positive-plus-binary branches each perform
  # their own power search, then null, binary-only, and both sets of positive
  # candidates compete by the requested information criterion. This avoids
  # retaining powers optimized conditionally on the binary component.
  if (saz_active && criterion_lower %in% c("aic", "bic")) {
    fit1 <- do.call(select_saz_ic, selector_args)

    if (verbose) {
      print_mfp_step(xi = xi, criterion = criterion, fit = fit1)
    }

    return(list(
      power_best = normalize_selected_step_powers(fit1, xi, acdx),
      searched_fp = selected_step_is_searched_fp(fit1),
      spike_decision = fit1$spike_decision,
      current_adj_params = fit1$current_adj_params,
      transform_cache = fit1$transform_cache
    ))
  }

  # Path 3: ordinary FP/ACD selection, plus Stage 1 of p-value SAZ selection.
  # P-value SAZ deliberately retains the historical two-stage split because
  # it is part of the closed-test construction. Non-SAZ AIC/BIC selection also
  # uses these ordinary selectors and is unaffected by the joint SAZ path.
  if (df == 1) {
    select_fct <- select_linear
  } else if (acdx[xi]) {
    select_fct <- if (criterion_lower == "pvalue") select_ra2_acd else select_ic_acd
  } else {
    select_fct <- if (criterion_lower == "pvalue") select_ra2 else select_ic
  }

  # Reuse the previous search for xi when nothing that defines its adjustment
  # model has changed since then (see fit_mfp()). The search depends only on
  # run-constant inputs plus the current forms of the other terms and the SAZ
  # decisions, so the selector result is identical. SAZ terms are excluded
  # because their stage 2 needs the uncompacted adjustment matrices.
  fit1 <- NULL
  memo_key <- NULL
  if (is.environment(step_memo) && !isTRUE(spike[[xi]])) {
    other_powers <- powers_current
    other_powers[[xi]] <- NULL
    memo_key <- list(
      df = df,
      powers = lapply(other_powers, unname),
      spike_decision = spike_decision
    )
    cached <- step_memo[[xi]]
    if (!is.null(cached) && identical(cached$key, memo_key)) {
      fit1 <- cached$fit1
      fit1$transform_cache <- transform_cache
    }
  }
  if (is.null(fit1)) {
    fit1 <- do.call(select_fct, selector_args)
    if (!is.null(memo_key)) {
      stored <- fit1
      stored$transform_cache <- NULL
      if (!is.null(stored$current_adj_params[[xi]])) {
        stored$current_adj_params[[xi]] <-
          compact_prev_adj_cache_entry(stored$current_adj_params[[xi]])
      }
      assign(xi, list(key = memo_key, fit1 = stored), envir = step_memo)
    }
  }

  if (verbose) {
    print_mfp_step(xi = xi, criterion = criterion, fit = fit1)
  }

  power_best <- normalize_selected_step_powers(fit1, xi, acdx)
  searched_fp <- selected_step_is_searched_fp(fit1)

  # Step 4: If xi was eliminated (stage 1 selected null), stop here ----------
  # Stage 1 selected null.
  # Nothing more to do.
  # For spike variables, clear stale binary-only decision.
  if (all(is.na(power_best))) {
    # Stage 2 is conditional on stage 1 selecting the full
    # continuous + binary spike model. If stage 1 selects null,
    # any previous binary-only decision is stale.
    if (isTRUE(spike[[xi]])) {
      # continuous_only is the neutral non-binary-only state here because the
      # variable will be eliminated.
      spike_decision[[xi]] <- saz_decision_codes[["continuous_only"]]
    }

    return(list(
      power_best = power_best,
      searched_fp = searched_fp,
      spike_decision = spike_decision,
      current_adj_params = fit1$current_adj_params,
      transform_cache = fit1$transform_cache
    ))
  }

  # Step 5: If xi is not a spike variable, stop here (no stage 2 needed) -----
  if (!isTRUE(spike[[xi]])) {
    # Stage 1 selected a non-null model,
    # but xi is not a spike variable.
    # Nothing more to do.
    return(list(
      power_best = power_best,
      searched_fp = searched_fp,
      spike_decision = spike_decision,
      current_adj_params = fit1$current_adj_params,
      transform_cache = fit1$transform_cache
    ))
  }

  # Path 4: Stage 2 is reached only by a selected, non-forced p-value SAZ term.
  # If we get here:
  # power_best is not NA
  # and xi is a spike variable
  # Therefore run SAZ stage 2.
  # ----------------------------------------------------------------------------
  # Evaluate spike at zero (SAZ) variables to update spike_decision.
  # This is stage 2 of the SAZ algorithm, computed only when the variable was
  # selected in stage 1 and is not force_max_fp. See evaluate_saz_stage2()
  # (spike_at_zero.R) for reduced-model fitting, the decision rule, and
  # stage-2 printing.
  # ----------------------------------------------------------------------------
  stage2 <- evaluate_saz_stage2(
    fit1 = fit1,
    xi = xi,
    power_best = power_best,
    y = y,
    weights = weights,
    offset = offset,
    family = family,
    family_string = family_string,
    fitter = fitter,
    method = method,
    strata = strata,
    nocenter = nocenter,
    control = control,
    rownames = rownames,
    has_offset = has_offset,
    n_obs = n_obs,
    criterion = criterion,
    alpha = alpha,
    ftest = ftest,
    spike_decision = spike_decision,
    verbose = verbose
  )

  return(list(
    power_best = power_best,
    searched_fp = if (identical(
      as.integer(stage2$spike_decision[[xi]]),
      saz_decision_codes[["binary_only"]]
    )) FALSE else searched_fp,
    spike_decision = stage2$spike_decision,
    current_adj_params = fit1$current_adj_params,
    transform_cache = fit1$transform_cache
  ))

}

#' Normalise and Name the Selected Power Vector From a Step Selector
#'
#' Collects the selected FP power vector returned by a single MFP or SAZ
#' step and normalises it to the naming convention used by the surrounding
#' engine. Centralising this in one helper is important for joint information-
#' criterion selection because the binary-only and null models both carry
#' `NA` powers but must be distinguished by their `spike_decision` value.
#'
#' @param fit A step-selection result carrying `power_best` and, for SAZ
#'   candidates, additional decision metadata.
#' @param xi Character scalar naming the predictor.
#' @param acdx Logical flag indicating whether the predictor is being
#'   assessed with the ACD extension.
#'
#' @return A named numeric vector of selected FP powers. Names identify the
#'   term (and its ACD companion, when applicable); missing values remain
#'   `NA_real_` so downstream code can distinguish binary-only from null
#'   representations. Searched-FP identity is returned separately by
#'   `selected_step_is_searched_fp()` so public power vectors remain plain.
#'
#' @keywords internal
#' @noRd
normalize_selected_step_powers <- function(fit, xi, acdx) {
  power_best <- as.numeric(fit$power_best)

  if (!fit$acd) {
    power_best <- power_best[!is.na(power_best)]
    if (length(power_best) == 0L) {
      power_best <- NA_real_
    }
  }

  names(power_best) <- name_transformed_variables(
    xi, length(power_best), acd = acdx[xi]
  )
  power_best
}

#' Identify Whether a Step Selected a Searched FP Form
#'
#' Keeps model-class identity separate from the numeric power vector. This is
#' necessary for FP1(power 1), which has the same numerical basis as the fixed
#' linear row but carries the searched-FP selection degrees of freedom.
#'
#' @param fit A step-selection result with `metrics` and `model_best`.
#' @return A single logical value.
#' @keywords internal
#' @noRd
selected_step_is_searched_fp <- function(fit) {
  metric_names <- if (is.null(fit$metrics)) NULL else rownames(fit$metrics)
  if (is.null(metric_names) || length(fit$model_best) != 1L ||
      is.na(fit$model_best) || fit$model_best < 1L ||
      fit$model_best > length(metric_names)) {
    return(FALSE)
  }

  # Candidate row labels are generated internally as FP1, FP2, ... (possibly
  # followed by an ACD or binary suffix). A fixed prefix check is sufficient
  # and avoids regex-based discovery of model identity.
  startsWith(trimws(metric_names[fit$model_best]), "FP")
}


#' Find the best FP function of a given degree for one variable
#'
#' Searches FP1 or higher-degree FP candidates. For parameter definitions, see
#' \code{find_best_fp_step()}.
#'
#' @details
#' Within a given FP degree, the candidate power combinations have the same
#' selection degrees of freedom. Likelihood-based models therefore choose the
#' candidate with the largest log-likelihood. For GEE models selected by AIC or
#' BIC, the function chooses the largest quasi-likelihood \eqn{Q}. Within a
#' fixed degree, this also minimizes QICu or QBIC because the penalty is
#' constant. Stata-compatible GEE p-value selection is the exception: it
#' chooses the smallest negative overall robust Wald chi-square.
#'
#' Each selected FP power adds one selection degree of freedom beyond its
#' fitted coefficient. The FP part of FP1 therefore has 2 selection degrees
#' of freedom, FP2 has 4, and FPm has \eqn{2m}. A separate zero indicator,
#' when included, adds its own fitted coefficient.
#'
#' When \code{degree = 1}, the function searches every power in the supplied
#' set. Power 1 is included only if that set contains it. A fixed linear
#' model and a searched FP1 model with power 1 have the same fitted mean
#' function, but different selection degrees of freedom: the fixed linear
#' form has 1, while the searched FP1 form has 2. If a matching linear fit is
#' supplied, its fitted-model results can be reused with the additional FP
#' selection degree of freedom. This distinction also applies when the
#' function is used in GEE, ACD, spike-at-zero, forced-FP, or MFPI selection.
#'
#' Candidates that cannot be fitted are skipped. Selection is made among the
#' candidates with valid results.
#'
#' @section ACD transformation:
#' The function also searches ACD forms when \code{acdx} is \code{TRUE} for
#' \code{xi}. With the default set of eight powers, degree 1 generates 8
#' candidates and degree 2 generates 64 candidates. Without ACD, the same
#' default set generates 36 degree-2 candidates. A restricted power set
#' changes these counts. ACD degrees above 2 are not supported.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{acd}: whether ACD was active for the variable.
#'   \item \code{powers}: a matrix of the candidate power combinations.
#'   \item \code{power_best}: the selected power or powers.
#'   \item \code{metrics}: a matrix of results for the candidates; entries
#'     for candidates that could not be fitted remain missing.
#'   \item \code{model_best}: the row number of the selected candidate in
#'     \code{metrics}.
#'   \item \code{zero}: whether exact-zero handling was active.
#'   \item \code{catzero}: whether a zero indicator was included.
#'   \item \code{current_adj_params}: information about the selected term
#'     and adjustment terms that can be reused in later steps.
#' }
#'
#' @inheritParams find_best_fp_step
#' @param degree FP degree \code{m}, equal to the number of transformed power
#'   terms.
#' @param linear_fit An existing linear fit for the same variable and
#'   adjustment terms, if available. Its results are reused for the power-1
#'   candidate when that candidate is present.
#' @param ... Arguments passed to \code{fit_model()}.
#' @keywords internal
#' @noRd
find_best_fpm_step <- function(x,
                               xi,
                               degree,
                               y,
                               powers_current,
                               powers,
                               acdx,
                               family,
                               family_string,
                               zero,
                               catzero, # a list of binary variables or null
                               spike,
                               spike_decision, # a numeric vector
                               acd_parameter,
                               prev_adj_params,
                               has_offset,
                               precomputed_adj = NULL,
                               focal_basis_cache = NULL,
                               n_obs,
                               term_to_columns,
                               criterion = "pvalue",
                               linear_fit = NULL,
                               ...
) {
  # fit_mfp() prepares the mapping from terms to design-matrix columns and
  # passes it through the selection steps. This keeps grouped terms together
  # when fitting the focal term and the adjustment terms.
  #
  # fit_mfp() also calculates n_obs once. It is usually the number of
  # observations, but may be an event count, an effective observation count,
  # or a cluster count, depending on the model.

  # Search the powers supplied for xi. A degree-1 search includes power 1
  # only when it is present in that set; it is not added automatically.

  # Prepare candidate transformations for xi and the current adjustment
  # terms. Reuse precomputed adjustment data when available. Ordinary FP
  # and ACD searches normally share a compact basis across candidates;
  # the selected candidate is assembled in full after the search.
  x_transformed <- transform_data_step(
    x = x, xi = xi, df = 2 * degree, powers_current = powers_current,
    powers = powers, acdx = acdx, zero = zero, catzero = catzero,
    spike = spike, spike_decision = spike_decision,
    acd_parameter = acd_parameter,
    prev_adj_params = prev_adj_params,
    precomputed_adj = precomputed_adj,
    focal_basis_cache = focal_basis_cache,
    term_to_columns = term_to_columns,
    compact_fp = !isTRUE(acdx[[xi]]),
    compact_acd = isTRUE(acdx[[xi]])
  )

  # Use exact list extraction: persistent cache entries intentionally retain
  # `data_adj = NULL` beside `data_adj_list`, and partial `$` matching must
  # never substitute the per-variable list for the assembled numeric matrix.
  data_adj <- x_transformed[["data_adj", exact = TRUE]]
  # data_fp is NULL by design during an FP or ACD search. find_best_fpm_step()
  # requests a compact representation, so transform_data_step() puts the
  # candidates in fp_basis or acd_basis instead of building a separate full
  # matrix for every candidate.
  data_fp <- x_transformed$data_fp
  fp_basis <- x_transformed$fp_basis
  acd_basis <- x_transformed$acd_basis

  # The fitting loop is agnostic to whether the shared columns came from an
  # ordinary FP basis or an ACD basis. Both expose the same basis/candidate_map
  # contract and use the same in-place C++ column-copy helper.
  has_fp_basis <- !is.null(fp_basis)
  has_acd_basis <- !is.null(acd_basis)

  # A focal variable can use either the ordinary FP compact representation
  # or the ACD compact representation, never both. Both may be NULL when
  # the historical materialized-candidate path is used. Enforce this
  # invariant explicitly so a future control-flow change cannot silently
  # select the wrong compact basis.
  if (has_fp_basis && has_acd_basis) {
    stop(
      "Internal error: both `fp_basis` and `acd_basis` are non-NULL.",
      call. = FALSE
    )
  }

  compact_basis <- if (has_fp_basis) {
    fp_basis
  } else if (has_acd_basis) {
    acd_basis
  } else {
    NULL
  }
  use_compact_basis <- !is.null(compact_basis)
  has_adj <- !is.null(data_adj) && NCOL(data_adj) > 0L

  # The repeated candidate matrix carries an explicit intercept for GLMs,
  # multinomial models, and survreg. Ordinal thresholds are added internally by
  # orm.fit(), so adding then removing an all-ones column would copy every
  # ordinal candidate design unnecessarily. Other ordinal paths retain their
  # existing matrix convention for now.
  x_has_intercept <- mfp2_candidate_matrix_has_intercept(family_string)


  # Build a reusable design-matrix template ---------------------------
  # Compact FP and ACD candidates no longer exist as lists of complete n x
  # degree matrices. All unique transformed terms live in one shared basis, and
  # each candidate is represented by one row of a small integer candidate_map.
  # The design matrix is allocated once and its focal columns are overwritten
  # directly from that basis in C++ for each candidate.
  if (use_compact_basis) {
    basis <- compact_basis$basis
    candidate_map <- compact_basis$candidate_map
    n_candidates <- nrow(candidate_map)
    n_candidate_cols <- ncol(candidate_map)
    use_catzero <- !is.null(compact_basis$catzero)
    n_xi_cols <- n_candidate_cols + as.integer(use_catzero)

    # Placeholder focal block used only to allocate the reusable design matrix.
    # catzero is invariant across candidates, so copy it once here rather than
    # once per candidate.
    xi_template <- matrix(
      0,
      nrow = nrow(basis),
      ncol = n_xi_cols
    )

    if (use_catzero) {
      xi_template[, 1L] <- compact_basis$catzero[, 1L]
      colnames(xi_template) <- c(
        "catzero",
        paste0("V", seq_len(n_candidate_cols))
      )
    }

    xi_cols <- seq_len(n_xi_cols) + as.integer(x_has_intercept)
    basis_target_cols <- if (use_catzero) {
      xi_cols[-1L]
    } else {
      xi_cols
    }

    if (x_has_intercept || has_adj) {
      design_mat <- assemble_design_matrix(
        blocks = list(xi_template, data_adj),
        nobs = nrow(xi_template),
        intercept = x_has_intercept
      )
    } else {
      # Cox without adjustment variables still needs one reusable focal working
      # matrix because compact candidates are not separately materialized.
      design_mat <- xi_template
    }

  } else {
    # Historical materialized representation retained only for defensive
    # fallback paths that explicitly opt out of compact candidate generation.
    first_xi <- data_fp[[1L]]
    n_xi_cols <- NCOL(first_xi)
    n_candidates <- length(data_fp)

    # The xi block starts after the intercept for GLMs and in column 1 for Cox.
    xi_cols <- seq_len(n_xi_cols) + as.integer(x_has_intercept)

    design_mat <- NULL
    if (x_has_intercept || has_adj) {
      design_mat <- assemble_design_matrix(
        blocks = list(first_xi, data_adj),
        nobs = nrow(first_xi),
        intercept = x_has_intercept
      )
    }
  }

  # GEE fits need an additional numerical column scale after the FP transform.
  # The adjustment columns and every shared basis column remain unchanged
  # throughout this search, so prepare their scaled copies only once. Linear,
  # null, and non-compact GEE fits continue using fit_gee()'s ordinary path.
  gee_scaling <- NULL
  if (mfp2_family_is_gee(family_string) && use_compact_basis) {
    gee_scaling <- prepare_gee_candidate_scaling(design_mat, basis)
    design_mat <- gee_scaling$design
    basis <- gee_scaling$basis
  }

  # Step 4: Fit one model per candidate FP power set and score it ------------
  metric_names <- c(
    "logl",
    "df",
    "deviance_rs",
    "deviance_gaussian",
    "aic",
    "bic",
    "df_resid"
  )
  metrics <- matrix(
    NA_real_,
    nrow = n_candidates,
    ncol = length(metric_names),
    dimnames = list(NULL, metric_names)
  )
  # Power 1 has the same mean model as the fixed linear rung. Keep its
  # candidate row (and its FP1 selection penalty), but reuse the fitted model.
  linear_candidate <- integer(0)
  if (degree == 1L && !is.null(linear_fit)) {
    power_rows <- x_transformed$powers_fp
    if (isTRUE(acdx[[xi]])) {
      linear_candidate <- which(is.na(power_rows[, 1L]) &
                                  power_rows[, 2L] == 1)
    } else {
      linear_candidate <- which(power_rows[, 1L] == 1)
    }
  }
  # All candidates in this batch have the same design width and class count.
  # Build nnet's mask and fixed offset-weight layout once rather than once per
  # transformed FP candidate.
  multinomial_optimizer <- if (identical(family_string, "multinomial")) {
    mfp2_multinomial_optimizer_structure(
      p = ncol(design_mat),
      c_classes = family$prepared$n_classes,
      has_offset = family$prepared$has_offset
    )
  } else {
    NULL
  }
  # Correlation is re-estimated for every distinct GEE candidate. Carry only
  # its last successful one-parameter estimate as the next candidate's start;
  # the cache is local to this power search and cannot cross focal variables.
  gee_alpha_start <- NULL

  # A single candidate that cannot be fitted (no convergence, numerically
  # rank-deficient design) is skipped instead of stopping the whole analysis.
  # The null, linear and final fits are not wrapped and remain strict.
  failed_candidates <- character(0)
  fit_candidate <- function(...) {
    tryCatch(
      fit_model(...),
      mfp2_candidate_failure = function(e) {
        failed_candidates[[length(failed_candidates) + 1L]] <<- conditionMessage(e)
        NULL
      }
    )
  }

  for (i in seq_len(n_candidates)) {
    if (i %in% linear_candidate) {
      metrics[i, ] <- mfp2_shift_selection_df(
        linear_fit$metrics[1L, ], 1, n_obs,
        qbic_n = if (mfp2_family_is_gee(family_string))
          family$prepared$qbic_n else n_obs
      )
      next
    }
    if (use_compact_basis) {
      # Copy only the degree-sized column map. The n-length transformed values
      # remain in the shared FP/ACD basis and are copied directly into the
      # reusable design matrix in C++, avoiding an n x degree temporary matrix.
      design_mat <- copy_fp_basis_candidate_cpp(
        target = design_mat,
        basis = basis,
        source_cols = as.integer(candidate_map[i, ]),
        target_cols = as.integer(basis_target_cols)
      )

      # Only the focal columns change between candidates. The remaining
      # factors were computed once along with the scaled adjustment design.
      gee_column_scales <- NULL
      if (!is.null(gee_scaling)) {
        gee_column_scales <- gee_scaling$column_scales
        gee_column_scales[basis_target_cols] <-
          gee_scaling$basis_scales[as.integer(candidate_map[i, ])]
      }

      fit <- fit_candidate(
        x               = design_mat,
        y               = y,
        family          = family,
        family_string   = family_string,
        has_offset      = has_offset,
        x_has_intercept = x_has_intercept,
        keep_coefficients = FALSE,
        multinomial_optimizer = multinomial_optimizer,
        gee_alpha_start = gee_alpha_start,
        gee_column_scales = gee_column_scales,
        gee_selection_criterion = criterion,
        gee_allow_failed_candidate = mfp2_family_is_gee(family_string),
        ...
      )

    } else {
      # Defensive materialized-candidate fallback.
      data_xi <- data_fp[[i]]

      if (!is.null(design_mat) && NCOL(data_xi) == n_xi_cols) {
        design_mat[, xi_cols] <- data_xi
        fit <- fit_candidate(
          x               = design_mat,
          y               = y,
          family          = family,
          family_string   = family_string,
          has_offset      = has_offset,
          x_has_intercept = x_has_intercept,
          keep_coefficients = FALSE,
          multinomial_optimizer = multinomial_optimizer,
          gee_alpha_start = gee_alpha_start,
          gee_selection_criterion = criterion,
          gee_allow_failed_candidate = mfp2_family_is_gee(family_string),
          ...
        )
      } else {
        # Defensive fallback: candidate generation should produce a fixed
        # number of xi columns within one degree, but preserve behaviour if it
        # ever does not. Build the destination in one allocation rather than
        # chaining cbind() calls.
        if (!x_has_intercept && !has_adj) {
          x_fit <- data_xi
        } else {
          x_fit <- assemble_design_matrix(
            blocks = list(data_xi, data_adj),
            nobs = nrow(data_xi),
            intercept = x_has_intercept
          )
        }
        fit <- fit_candidate(
          x               = x_fit,
          y               = y,
          family          = family,
          family_string   = family_string,
          has_offset      = has_offset,
          x_has_intercept = x_has_intercept,
          keep_coefficients = FALSE,
          multinomial_optimizer = if (
            !is.null(multinomial_optimizer) &&
            ncol(x_fit) == multinomial_optimizer$p
          ) multinomial_optimizer else NULL,
          gee_alpha_start = gee_alpha_start,
          gee_selection_criterion = criterion,
          gee_allow_failed_candidate = mfp2_family_is_gee(family_string),
          ...
        )
      }
    }

    if (is.null(fit)) {
      next
    }

    if (mfp2_family_is_gee(family_string) && !is.null(fit$gee_alpha)) {
      gee_alpha_start <- fit$gee_alpha
    }

    # A failed GEE power candidate has no valid quasi-likelihood or Wald
    # score. Leave its row missing and choose only among successful candidates.
    # Required null/linear and final fits do not opt into this exclusion.
    if (mfp2_family_is_gee(family_string) && identical(fit$converged, FALSE)) {
      next
    }

    metrics[i, ] <- calculate_model_metrics(
      fit, # fit includes catzero df when used, so this part does not change
      n_obs,
      degree
    )

  }

  # Step 5: Pick the best candidate within this fixed FP degree. Ordinary
  # likelihood models maximize log-likelihood, while GEE information-criterion
  # paths maximize quasi-likelihood Q. GEE p-value selection is Stata-compatible
  # and instead minimizes its negative-global-Wald selection deviance.
  gee_pvalue <- mfp2_family_is_gee(family_string) &&
    identical(tolower(criterion), "pvalue")
  selection_score <- if (gee_pvalue) "deviance_rs" else "logl"
  valid_candidates <- which(is.finite(metrics[, selection_score]))
  if (length(failed_candidates) > 0L) {
    warning(
      sprintf(
        "i %d of %d FP%s candidate model(s) for '%s' could not be fitted and were skipped.\ni First problem: %s\ni Setting an explicit shift for '%s' may avoid this.",
        length(failed_candidates), n_candidates, degree, xi,
        failed_candidates[[1L]], xi
      ),
      call. = FALSE
    )
  }
  if (length(valid_candidates) == 0L) {
    stop(
      sprintf("No converged %s FP candidates are available for '%s' at degree %s.",
              family_string, xi, degree),
      call. = FALSE
    )
  }
  model_best <- if (gee_pvalue) {
    valid_candidates[[which.min(metrics[valid_candidates, selection_score])]]
  } else {
    valid_candidates[[which.max(metrics[valid_candidates, selection_score])]]
  }

  # SAZ stage 2 needs the complete winning xi matrix. Under either compact
  # representation, materialize exactly that one candidate after model selection.
  x_transformed$current_params[[xi]]$data_xi <- if (!is.null(fp_basis)) {
    materialize_fp_basis_candidate(fp_basis, model_best)
  } else if (!is.null(acd_basis)) {
    materialize_acd_basis_candidate(acd_basis, model_best)
  } else {
    data_fp[[model_best]]
  }

  list(
    acd = acdx[xi],
    powers = x_transformed$powers_fp,
    power_best = x_transformed$powers_fp[model_best, , drop = TRUE],
    metrics = metrics,
    model_best = model_best,
    zero = zero[xi],
    catzero = ifelse(!is.null(catzero[[xi]]), TRUE, FALSE),
    current_adj_params = x_transformed$current_params
  )
}

#' Fit a model without the current term
#'
#' Fits a model using the current adjustment terms while excluding \code{xi}
#' and any columns belonging to it. The adjustment data are built from the
#' current forms of the other terms, or reused when already available.
#'
#' @return A list with three entries:
#' \itemize{
#'   \item \code{powers}: \code{NA}, indicating that \code{xi} is absent.
#'   \item \code{metrics}: a one-row matrix, labelled \code{null}, containing
#'     the fitted model's performance measures.
#'   \item \code{current_adj_params}: adjustment-term powers, spike decisions,
#'     and transformed data, stored under \code{xi} for reuse in later steps.
#' }
#'
#' @inheritParams find_best_fp_step
#' @param precomputed_adj Previously prepared adjustment data. If \code{NULL},
#'   the function builds the adjustment data for this step.
#' @param ... Additional arguments passed to \code{fit_model()}.
#' @keywords internal
#' @noRd
fit_null_step <- function(x,
                          xi,
                          y,
                          powers_current,
                          powers,
                          acdx,
                          family,
                          family_string,
                          zero,
                          catzero,
                          spike,
                          spike_decision,
                          acd_parameter,
                          prev_adj_params,
                          has_offset,
                          precomputed_adj = NULL,
                          n_obs,
                          term_to_columns,
                          ...
) {

  # The lookup is passed explicitly even when a precomputed adjustment block
  # is available, so this helper remains consistent when it builds the block.
  # Step 1: Build (or reuse) the adjustment-variable design matrix ----------
  # Null model uses adjustment variables only. Do not call
  # transform_data_step(df = 1) here: that would generate the current-variable
  # linear candidate and then discard it.
  adj <- if (is.null(precomputed_adj)) {
    build_adjustment_step(
      x = x,
      xi = xi,
      powers_current = powers_current,
      powers = powers,
      acdx = acdx,
      zero = zero,
      catzero = catzero,
      spike = spike,
      spike_decision = spike_decision,
      acd_parameter = acd_parameter,
      prev_adj_params = prev_adj_params,
      term_to_columns = term_to_columns
    )
  } else {
    precomputed_adj
  }

  current_params <- list()
  current_params[[xi]] <- list(
    powers_adj = adj$powers_adj,
    spike_decision_adj = adj$spike_decision_adj,
    data_adj_list = adj$data_adj_list,
    data_adj = adj[["data_adj", exact = TRUE]]
  )

  # Step 2: Fit the null model (adjustment variables only, xi excluded) -----
  # An empty matrix can be returned which creates problem for cox model
  # convert it to NULL
  x_tran <- adj[["data_adj", exact = TRUE]]
  if (!is.null(x_tran) && !is.matrix(x_tran)) {
    stop("Internal error: adjustment data must be a matrix or NULL.", call. = FALSE)
  }
  if (is.null(x_tran) || NCOL(x_tran) == 0L) {
    x_tran <- NULL
  }

  # fit null model
  # i.e. a model that does not contain xi but only adjustment variables
  # In addition, adjustment model can be NULL, so we have an intercept only
  model_null <- fit_model(x = x_tran,
                          y = y,
                          family = family,
                          family_string   = family_string,
                          has_offset = has_offset,
                          keep_coefficients = FALSE,
                          ...
  )

  # Step 3: Score the null model and return, with xi's power fixed at NA ----
  list(
    powers = NA,
    metrics = rbind(null = calculate_model_metrics(model_null, n_obs)),
    current_adj_params = current_params
  )
}

#' Change the selection degrees of freedom without refitting a model
#'
#' A fixed linear term and a searched FP1 term with power 1 produce the same
#' fitted model, but the FP1 term receives one additional selection degree of
#' freedom because its power was chosen from a candidate set. This function
#' adjusts the scores for that difference while leaving the fitted likelihood
#' and deviance unchanged. It can also reverse the adjustment.
#'
#' The same calculation applies to GEE scores. Its robust Wald selection
#' deviance is left unchanged; only the degree-of-freedom count and associated
#' penalties are adjusted.
#'
#' @param metrics Named numeric vector of model scores, including \code{df},
#'   \code{aic}, \code{bic}, and \code{df_resid}.
#' @param additional Number of selection degrees of freedom to add. Use
#'   \code{1} to change a fixed linear score to a searched FP1 score, or
#'   \code{-1} for the reverse change.
#' @param n_obs Sample size used as the default BIC penalty count when
#'   \code{qbic_n} is not supplied.
#' @param qbic_n Count used in the BIC penalty. This may differ from
#'   \code{n_obs} for GEE.
#'
#' @return The updated named score vector. \code{df} increases by
#'   \code{additional}, \code{df_resid} decreases by the same amount,
#'   \code{aic} changes by twice that amount, and \code{bic} changes by
#'   \code{log(qbic_n)} times that amount. Other entries are unchanged.
#' @keywords internal
#' @noRd
mfp2_shift_selection_df <- function(metrics, additional, n_obs,
                                    qbic_n = n_obs) {
  metrics["df"] <- metrics["df"] + additional
  metrics["aic"] <- metrics["aic"] + 2 * additional
  metrics["bic"] <- metrics["bic"] + log(qbic_n) * additional
  metrics["df_resid"] <- metrics["df_resid"] - additional
  metrics
}

#' Fit the current term in its fixed linear form
#'
#' Fits \code{xi} alongside the current adjustment terms. An ordinary numeric
#' term enters with power 1. When ACD is active, the linear form uses
#' \code{A(x)}. A grouped term enters as its complete fixed linear block.
#' If an exact-zero indicator is active for \code{xi}, that indicator is
#' included alongside the term.
#'
#' If a supplied FP1 search already fitted the matching power-1 candidate,
#' its results are reused. The FP1 selection degree of freedom is removed
#' from the reported scores because this function represents a fixed form.
#' Otherwise, the function fits and scores the linear model directly.
#'
#' @return A list with three entries:
#' \itemize{
#'   \item \code{powers}: a matrix describing the fixed form of \code{xi}.
#'   \item \code{metrics}: a one-row matrix of model scores.
#'   \item \code{current_adj_params}: the transformed focal and adjustment
#'     data and related information for reuse in later steps.
#' }
#'
#' @inheritParams find_best_fp_step
#' @param precomputed_adj Previously prepared adjustment data. If \code{NULL},
#'   the adjustment data are prepared during this step.
#' @param fp1_fit Results from an FP1 search for the same term and adjustment
#'   set. Its power-1 candidate is reused when present.
#' @param ... Additional arguments passed to \code{fit_model()} when a new
#'   model fit is needed.
#' @keywords internal
#' @noRd
fit_linear_step <- function(x,
                            xi,
                            y,
                            powers_current,
                            powers,
                            acdx,
                            family,
                            family_string,
                            zero,
                            catzero,
                            spike,
                            spike_decision,
                            acd_parameter,
                            prev_adj_params,
                            has_offset,
                            n_obs,
                            precomputed_adj = NULL,
                            fp1_fit = NULL,
                            term_to_columns,
                            ...
                            ) {
  # n_obs (number of observations, or target events for PH models) is computed once
  # in fit_mfp() and passed down as a parameter, rather than recomputed here.

  # Step 1: Transform xi as a linear term (power 1) plus the adjustment set -
  # transform all data as given by current working model
  # set variable of interest to linear term only
  x_transformed <- transform_data_step(
    x = x, xi = xi, df = 1, powers_current = powers_current, acdx = acdx,
    powers = powers, zero = zero, catzero = catzero, spike = spike,
    spike_decision = spike_decision, acd_parameter = acd_parameter,
    prev_adj_params = prev_adj_params,
    precomputed_adj = precomputed_adj,
    term_to_columns = term_to_columns
  )

  x_transformed$current_params[[xi]]$data_xi <- x_transformed$data_fp[[1]]

  # An FP1 search may have preceded this rung (the df = 2 closed test).
  # Reuse its power-1 candidate even when another power won the search.
  fp1_row <- integer(0)
  if (!is.null(fp1_fit)) {
    candidate_powers <- fp1_fit$powers
    if (isTRUE(acdx[[xi]])) {
      fp1_row <- which(is.na(candidate_powers[, 1L]) &
                         candidate_powers[, 2L] == 1)
    } else {
      fp1_row <- which(candidate_powers[, 1L] == 1)
    }
  }
  if (length(fp1_row) == 1L) {
    metrics <- rbind(linear = mfp2_shift_selection_df(
      fp1_fit$metrics[fp1_row, ], -1, n_obs,
      qbic_n = if (mfp2_family_is_gee(family_string))
        family$prepared$qbic_n else n_obs
    ))
    if (isTRUE(acdx[[xi]])) rownames(metrics) <- "linear(., A(x))"
    return(list(
      powers = x_transformed$powers_fp,
      metrics = metrics,
      current_adj_params = x_transformed$current_params
    ))
  }

  # Step 2: Assemble the design matrix and fit the linear model -------------
  # Fit a model based on the assumption that xi is linear.
  # If catzero or spike-at-zero is active for xi, the corresponding binary
  # component is already included in data_xi.
  data_xi <- x_transformed$data_fp[[1L]]
  data_adj <- x_transformed[["data_adj", exact = TRUE]]

  has_adj <- !is.null(data_adj) && NCOL(data_adj) > 0L
  use_glm_intercept_template <- mfp2_family_has_intercept(family_string)

  # Build the final matrix in one allocation. GLM fits receive a leading
  # intercept and set x_has_intercept below; Cox fits remain intercept-free.
  # When a Cox fit has no adjustment block, reuse data_xi without copying it.
  if (!use_glm_intercept_template && !has_adj) {
    x_fit <- data_xi
  } else {
    x_fit <- assemble_design_matrix(
      blocks = list(data_xi, data_adj),
      nobs = nrow(data_xi),
      intercept = use_glm_intercept_template
    )
  }

  model_linear <- fit_model(
    x = x_fit,
    y = y,
    family = family,
    family_string = family_string,
    has_offset = has_offset,
    x_has_intercept = use_glm_intercept_template,
    keep_coefficients = FALSE,
    ...
  )

  # Step 3: Score the linear model (0 additional df beyond the 1 already
  # counted for the linear term itself) and label the metrics row for ACD.
  metrics <- rbind(
    linear = calculate_model_metrics(
      obj = model_linear,
      n_obs = n_obs,
      df_additional = 0
    )
  )

  if (acdx[xi]) {
    rownames(metrics) <- "linear(., A(x))"
  }

  list(
    powers = x_transformed$powers_fp,
    metrics = metrics,
    current_adj_params = x_transformed$current_params
  )
}

#' Choose between excluding a term and fitting its fixed linear form
#'
#' Used when \code{xi} is restricted to a linear form. The function fits two
#' models with the same adjustment terms: one excludes \code{xi}, and the
#' other includes it. A term represented by several columns is included or
#' excluded as a whole. If an exact-zero indicator is active for \code{xi},
#' the included model contains both the term and that indicator.
#'
#' The adjustment terms are prepared once and shared by both fits. The
#' function also supports a fixed linear ACD form.
#'
#' @details
#' If \code{xi} is in \code{keep}, the included model is selected. Otherwise,
#' \code{"pvalue"} selection uses the null-versus-linear test and
#' \code{select}; \code{"aic"} and \code{"bic"} select the model with the
#' smaller corresponding score. The test statistic and p-value are calculated
#' for all three criteria, but do not determine AIC or BIC selection.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{keep}: whether \code{xi} is in \code{keep}.
#'   \item \code{acd}: whether ACD is active for \code{xi}.
#'   \item \code{powers}: a two-row matrix describing the excluded and
#'     included models.
#'   \item \code{power_best}: the powers in the selected row. An entirely
#'     missing result indicates exclusion; an ACD form may contain a missing
#'     entry for an unused component.
#'   \item \code{metrics}: a two-row matrix of model scores, in the same
#'     order as \code{powers}.
#'   \item \code{model_best}: the row number of the selected model.
#'   \item \code{pvalue} and \code{statistic}: results of the
#'     null-versus-linear comparison.
#'   \item \code{is_gee}: whether the model uses GEE.
#'   \item \code{zero}, \code{catzero}, and \code{spike}: the settings for
#'     exact-zero handling and spike-at-zero treatment.
#'   \item \code{current_adj_params}: information from the selected fit
#'     that can be reused in later steps.
#'   \item \code{transform_cache}: reusable transformed adjustment terms.
#' }
#'
#' @inheritParams find_best_fp_step
#' @param degree Not used; the included form is fixed as linear.
#' @param force_max_fp Not used directly. The caller handles this setting
#'   before invoking this function.
#' @param ... Additional arguments passed through the two fitting helpers
#'   to \code{fit_model()}.
#' @keywords internal
#' @noRd
select_linear <- function(x,
                          xi,
                          keep,
                          degree,
                          acdx,
                          y,
                          powers_current,
                          powers,
                          criterion,
                          ftest,
                          select,
                          alpha,
                          family,
                          family_string,
                          zero,
                          catzero,
                          spike,
                          spike_decision,
                          acd_parameter,
                          prev_adj_params,
                          transform_cache = NULL,
                          force_max_fp,
                          has_offset,
                          n_obs,
                          term_to_columns,
                          ...) {

  is_gee <- mfp2_family_is_gee(family_string)
  # `term_to_columns` is required here because the focal term or any of its
  # adjustment terms may correspond to more than one raw design-matrix column.
  # select_linear() is only used when df = 1 for xi (see find_best_fp_step()'s
  # dispatch logic): this covers variables that are restricted to a linear
  # effect, e.g. dummy/categorical variables or any continuous variable the
  # user explicitly restricted via df = 1. Because no FP alternative exists
  # for such variables, the closed-test procedure collapses to a single
  # null-vs-linear comparison

  # Step 1: Build the adjustment matrix (all variables except xi, transformed
  # at their *current* powers_current) once, and reuse it for both candidate
  # fits below. This is valid because neither the null nor the linear model
  # for xi changes what the adjustment variables look like - only xi's own
  # representation differs between the two candidates.
  precomputed_adj <- build_adjustment_step(
    x = x,
    xi = xi,
    powers_current = powers_current,
    powers = powers,
    acdx = acdx,
    zero = zero,
    catzero = catzero,
    spike = spike,
    spike_decision = spike_decision,
    acd_parameter = acd_parameter,
    prev_adj_params = prev_adj_params,
    transform_cache = transform_cache,
    term_to_columns = term_to_columns
  )

  # If xi itself contributes a structural-zero binary indicator (because it is
  # an active spike variable, or was passed through catzero_vars), the linear
  # model for xi always includes that indicator alongside the continuous term.
  # The " + Binary" suffix on row/test names makes this visible in printed
  # output and in the metrics/powers matrices below, without changing which
  # models are actually compared (there is still only null vs. linear).
  has_binary_xi <- isTRUE(spike[[xi]]) || !is.null(catzero[[xi]])
  binary_suffix <- if (has_binary_xi) " + Binary" else ""
  linear_name <- paste0("linear", binary_suffix)

  # Step 2: Fit the two candidate models -------------------------------------
  # Model 1: Null model (xi excluded entirely, adjustment variables only).
  fit_null <- fit_null_step(
    x = x, xi = xi, y = y,
    powers_current = powers_current, powers = powers, acdx = acdx,
    family = family, family_string = family_string, zero = zero,
    catzero = catzero, spike_decision = spike_decision, has_offset = has_offset,
    spike = spike, acd_parameter = acd_parameter,
    prev_adj_params = prev_adj_params, precomputed_adj = precomputed_adj,
    n_obs = n_obs, term_to_columns = term_to_columns, ...
  )

  # Model 2: Linear model (xi included with FP power fixed at 1).
  fit_linear <- fit_linear_step(
    x = x, xi = xi, y = y,
    powers_current = powers_current, powers = powers, acdx = acdx,
    family = family, family_string = family_string, zero = zero,
    catzero = catzero, spike_decision = spike_decision, spike = spike,
    has_offset = has_offset, n_obs = n_obs,
    acd_parameter = acd_parameter, prev_adj_params = prev_adj_params,
    precomputed_adj = precomputed_adj,
    term_to_columns = term_to_columns, ...
  )

  # Stack both candidates' powers/metrics into 2-row matrices so the rest of
  # this function (and the caller) can treat them uniformly, the same way
  # select_ra2()/select_ic() stack an arbitrary number of candidate rows.
  # Extract powers and metrics of interest
  powers <- rbind(fit_null$powers, fit_linear$powers)
  metrics <- rbind(fit_null$metrics, fit_linear$metrics)

  rownames(powers) <- c("null", linear_name)
  rownames(metrics) <- c("null", linear_name)

  # Step 3: Test null vs. linear jointly -----------------------------------
  # A singleton contributes one regression coefficient. A grouped categorical
  # term contributes one coefficient per raw design column, so the metric df
  # difference automatically gives the joint block-test degrees of freedom.
  if (ftest) {
    # note that ftest is only TRUE if model is gaussian
    stats <- calculate_f_test(
      deviances = metrics[, "deviance_gaussian"],
      dfs_resid = metrics[, "df_resid"],
      n_obs = n_obs
    )
  } else if (is_gee) {
    stats <- calculate_deviance_test(
      deviances = metrics[, "deviance_rs"],
      dfs = metrics[, "df"]
    )
  } else {
    # Ordinary likelihood-ratio chi-square test (2 * difference in log-lik),
    # compared with a chi-square distribution using the fitted-model df
    # difference (one for a singleton, group size for a grouped term).
    stats <- calculate_lr_test(metrics[, "logl"], metrics[, "df"])
  }

  # Compute the corresponding p-value
  pvalue <- stats$pvalue
  test_name <- paste0("null vs ", linear_name)

  names(pvalue) <- test_name
  statistic <- stats$statistic
  names(statistic) <- test_name

  # Step 4: Select null or linear, or force linear if xi is in `keep` --------
  # `keep` overrides ordinary selection entirely: a kept variable is never
  # dropped, so it is fixed at the linear model regardless of significance or
  # information criterion.
  # Check whether the variable should be forced into the model; index 1 denotes
  # a null, while 2 denotes a linear model
  if (xi %in% keep) {
    model_best <- 2
  } else {
    # For "pvalue": drop xi only when the null-vs-linear p-value strictly
    # exceeds `select`. The strict `>` is intentional and matches the original
    # Stata/CRAN MFP endpoint convention: select = 1 is a forcing value, so even
    # an exact p-value of 1 does not remove the term. This mirrors Test 1 of the
    # FP closed-test procedure (see select_ra2()), just with only one candidate
    # non-null model instead of several FP degrees.
    # For "aic"/"bic": simply pick whichever of the two rows has the smaller
    # criterion value; there's no separate functional-form step because a
    # linear-only variable has no functional form to choose between.
    model_best <- switch(
      tolower(criterion),
      "pvalue" = if (mfp_pvalue_exceeds(pvalue, select)) 1 else 2,
      "aic" = which.min(metrics[, "aic", drop = TRUE]),
      "bic" = which.min(metrics[, "bic", drop = TRUE])
    )
  }

  list(
    keep = xi %in% keep,
    acd = acdx[xi],
    powers = powers,
    power_best = powers[model_best, ],
    metrics = metrics,
    model_best = model_best,
    pvalue = pvalue,
    statistic = statistic,
    is_gee = is_gee,
    zero = zero[xi],
    catzero = ifelse(!is.null(catzero[[xi]]), TRUE, FALSE),
    spike = spike[xi],
    # Return the adjustment-variable cache belonging to whichever model was
    # actually selected (not always fit_linear's), so that the cache stored in
    # prev_adj_params for the *next* cycle correctly reflects xi's current
    # state (excluded vs. included) rather than always assuming inclusion.
    current_adj_params = if (model_best == 1) fit_null$current_adj_params else fit_linear$current_adj_params,
    transform_cache = precomputed_adj$transform_cache
  )
}

#' Select a fractional-polynomial form using closed tests
#'
#' Used by \code{find_best_fp_step()} for a non-ACD term when
#' \code{criterion = "pvalue"}. The function compares the best FP model at
#' the requested maximum degree with successively simpler forms.
#'
#' @details
#' The adjustment terms keep their current forms throughout these comparisons.
#' Within each FP degree, \code{find_best_fpm_step()} first selects the best
#' permitted power combination. For GEE p-value selection, that search uses
#' the robust Wald score; other supported models use likelihood.
#'
#' The procedure tests models in this order:
#' \enumerate{
#'   \item Compare the best FP model at the maximum degree with a model that
#'     excludes \code{xi}. If the comparison p-value is greater than
#'     \code{select}, exclude \code{xi}, unless it is in \code{keep}.
#'   \item Compare the maximum-degree FP model with a fixed linear form.
#'     If the p-value is greater than \code{alpha}, select the linear form.
#'   \item If the maximum degree exceeds 1, compare it with the best FP1,
#'     FP2, and subsequent lower-degree models in that order. Select the
#'     first lower-degree model whose comparison p-value is greater than
#'     \code{alpha}. If none qualifies, retain the maximum-degree model.
#' }
#'
#' A simpler form is selected only when its comparison p-value is strictly
#' greater than the relevant threshold. Thus \code{select = 1} prevents
#' exclusion, and \code{alpha = 1} prevents simplification, even when a
#' comparison p-value equals 1. Putting \code{xi} in \code{keep} prevents
#' exclusion but still allows selection of a simpler form.
#'
#' For an eligible spike-at-zero term, these comparisons form the first
#' selection stage. The calling function handles any subsequent choice
#' between its continuous and zero-indicator components.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{keep}: whether \code{xi} is in \code{keep}.
#'   \item \code{acd}: \code{FALSE}.
#'   \item \code{powers}: a matrix of powers for the candidate models reached
#'     before selection stops. Its rows begin with the maximum-degree FP
#'     model and the null model, followed when needed by the linear form
#'     and lower-degree FP models.
#'   \item \code{power_best}: the selected powers, or \code{NA} if
#'     \code{xi} was excluded.
#'   \item \code{metrics}: model scores in the same row order as
#'     \code{powers}.
#'   \item \code{model_best}: the selected row of \code{metrics}.
#'   \item \code{pvalue} and \code{statistic}: named results for each
#'     comparison performed before selection stopped.
#'   \item \code{is_gee}: whether the model uses GEE.
#'   \item \code{spike}: whether \code{xi} is being treated as a
#'     spike-at-zero term.
#'   \item \code{current_adj_params}: information from the selected fit
#'     for reuse in later steps.
#'   \item \code{transform_cache}: reusable transformed adjustment terms.
#' }
#'
#' @param degree Maximum FP degree considered for \code{xi}.
#' @param ... Additional arguments passed through the fitting helpers to
#'   \code{fit_model()}.
#' @inheritParams find_best_fp_step
#'
#' @references
#' Royston, P. and Sauerbrei, W. (2008). \emph{Multivariable Model-Building:
#' A Pragmatic Approach to Regression Analysis Based on Fractional
#' Polynomials for Modelling Continuous Variables}. John Wiley & Sons.
#'
#' @seealso \code{select_ra2_acd()}
#' @keywords internal
#' @noRd
select_ra2 <- function(x,
                       xi,
                       keep,
                       degree,
                       acdx,
                       y,
                       powers_current,
                       powers,
                       criterion,
                       ftest,
                       select,
                       alpha,
                       family,
                       family_string,
                       zero,
                       catzero,
                       spike,
                       spike_decision,
                       acd_parameter,
                       prev_adj_params,
                       transform_cache = NULL,
                       force_max_fp,
                       has_offset,
                       n_obs,
                       term_to_columns,
                       ...) {

  is_gee <- mfp2_family_is_gee(family_string)

  # select_ra2() implements the RA2 closed-test function-selection procedure
  # (Royston & Sauerbrei 2008, Ch. 4 & 6) for a single continuous variable xi
  # at a given candidate FP degree. Called once per FP degree from
  # find_best_fp_step() (criterion = "pvalue", non-ACD variables only).
  # `degree` is the *maximum* degree under consideration (df / 2 for xi); the
  # tests below successively fall back to simpler models (null, linear, then
  # FP1, FP2, ... up to degree - 1) until one is not significantly worse than
  # FPm(degree), or none are, in which case FPm(degree) itself is retained.
  #
  # degree = 1 (linear-only variables) is handled by select_linear() instead;
  # this guard is defensive and should not normally trigger from
  # find_best_fp_step()'s dispatch logic.
  if (degree < 1) {
    return(NULL)
  }

  # Step 0: Set up the test statistic helper and shared naming -------------
  # simplify testing by defining test helper function
  if (ftest) {
    calculate_test <- function(metrics, n_obs, reduced_fit = NULL, full_fit = NULL) {
      calculate_f_test(
        deviances = metrics[, "deviance_gaussian", drop = TRUE],
        dfs_resid = metrics[, "df_resid", drop = TRUE],
        n_obs = n_obs
      )
    }
  } else if (is_gee) {
    calculate_test <- function(metrics, n_obs, reduced_fit = NULL, full_fit = NULL) {
      calculate_deviance_test(
        deviances = metrics[, "deviance_rs", drop = TRUE],
        dfs = metrics[, "df", drop = TRUE]
      )
    }
  } else {
    calculate_test <- function(metrics, n_obs, reduced_fit = NULL, full_fit = NULL) {
      calculate_lr_test(
        logl = metrics[, "logl", drop = TRUE],
        dfs = metrics[, "df", drop = TRUE]
      )
    }
  }

  # step 1: setup output list
  # Model labels get a binary suffix when the focal variable contributes a
  # structural-zero indicator, either through spike handling or catzero.
  # Compute this once to keep row names, test labels, and candidate names
  # consistent throughout the selection function.
  has_binary_xi <- isTRUE(spike[[xi]]) || !is.null(catzero[[xi]])
  binary_suffix <- if (has_binary_xi) " + Binary" else ""

  fpmax <- paste0("FP", degree, binary_suffix)

  # output list
  # `res` accumulates one row per candidate model tested (in `metrics` and
  # `powers`) and one entry per pairwise test performed (in `statistic` and
  # `pvalue`) as the procedure works through Tests 1-3 below. Each early
  # return fixes `model_best` as a row index into `res$metrics`/`res$powers`
  # at that point.
  res <- list(
    keep = xi %in% keep,
    acd = FALSE,
    powers = NULL,
    power_best = NULL,
    metrics = NULL,
    model_best = NULL,
    statistic = NULL,
    pvalue = NULL,
    is_gee = is_gee,
    spike = spike[xi],
    current_adj_params = NULL,
    transform_cache = NULL
  )

  # Step 2: Build the adjustment matrix (all variables except xi, transformed
  # at their *current* powers_current) once, and reuse it for all candidate
  # fits below. This is valid because neither the null nor the FPm models
  # for xi changes what the adjustment variables look like - only xi's own
  # representation differs between the candidate models.
  # Thread the term lookup into adjustment construction so grouped terms are
  # included as complete blocks while xi undergoes ordinary FP selection.
  precomputed_adj <- build_adjustment_step(
    x = x,
    xi = xi,
    powers_current = powers_current,
    powers = powers,
    acdx = acdx,
    zero = zero,
    catzero = catzero,
    spike = spike,
    spike_decision = spike_decision,
    acd_parameter = acd_parameter,
    prev_adj_params = prev_adj_params,
    transform_cache = transform_cache,
    term_to_columns = term_to_columns
  )
  res$transform_cache <- precomputed_adj$transform_cache

  # Step 3: Test null vs. FPm. Asks "is xi associated with the outcome at all,
  # in its most flexible form at this degree?" If not significant (at the
  # `select` level), xi is dropped, unless forced via `keep`.

  # fit highest fp (FPm) and null model
  fit_fpmax <- find_best_fpm_step(
    x = x, xi = xi, degree = degree, y = y,
    powers_current = powers_current, powers = powers, acdx = acdx,
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter, spike = spike,
    prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
    precomputed_adj = precomputed_adj, term_to_columns = term_to_columns,
    criterion = criterion, ...
  )

  fit_null <- fit_null_step(
    x = x, xi = xi, y = y,
    powers_current = powers_current, powers = powers, acdx = acdx,
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter, spike = spike,
    prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
    precomputed_adj = precomputed_adj, term_to_columns = term_to_columns, ...
  )
  res$metrics <- rbind(
    fit_fpmax$metrics[fit_fpmax$model_best, ],
    fit_null$metrics
  )
  rownames(res$metrics) <- c(fpmax, "null")
  # Row 1 = FPm, the reference model every later test compares
  # against; row 2 = null model, whose power is NA since xi is excluded.
  res$powers <- rbind(fit_fpmax$power_best, NA)

  # return also the adjustment parameters
  res$current_adj_params <- fit_fpmax$current_adj_params

  # Test 1: test for overall significance (null vs best FPm)
  # df for tests are 2 * degree
  stats <- calculate_test(
    res$metrics[c("null", fpmax), ], n_obs,
    reduced_fit = fit_null, full_fit = fit_fpmax
  )
  res$statistic <- stats$statistic
  names(res$statistic) <- sprintf("%s vs null", fpmax)
  res$pvalue <- stats$pvalue
  names(res$pvalue) <- names(res$statistic)

  # Ensure keep is not NULL
  current_keep <- if (is.null(keep)) character(0) else keep

  # MFP uses a strict upper-tail boundary for simplification/removal. Therefore
  # select = 1 truly forces inclusion: p = 1 is retained, and only p > select
  # can trigger removal. The same convention is used for alpha below.
  if (mfp_pvalue_exceeds(stats$pvalue, select) && !(xi %in% current_keep)) {
    # Test 1 not significant and xi not forced: eliminate xi (model_best = 2,
    # the null row).
    # not selected and not forced into model
    res$power_best = NA
    res$model_best = 2
    res$current_adj_params <- fit_null$current_adj_params
    return(res)
  }

  # Step 4: Test linear vs. FPm. xi is significant overall (or forced in); now
  # ask "is the extra flexibility of FPm actually needed, or would a simple
  # linear term do just as well?" This test has one fewer df than Test 1 because
  # the linear model already spends 1 df on xi (vs. 0 for the null model).

  fit_lin <- fit_linear_step(
    x = x, xi = xi, y = y,
    powers_current = powers_current, powers = powers, acdx = acdx,
    family = family,family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter, spike = spike,
    prev_adj_params = prev_adj_params, has_offset = has_offset,n_obs = n_obs,
    fp1_fit = if (degree == 1L) fit_fpmax else NULL,
    precomputed_adj = precomputed_adj, term_to_columns = term_to_columns, ...
  )

  old_names <- rownames(res$metrics)
  res$metrics <- rbind(
    res$metrics,
    fit_lin$metrics
  )

  lin_names <- paste0("linear", binary_suffix)
  rownames(res$metrics) <- c(old_names, lin_names)
  # The linear candidate has 1 power (fixed at 1); pad with NA so its row
  # matches res$powers' width (degree columns, one per FPm(degree) power).
  res$powers = rbind(res$powers,
                     ensure_length(fit_lin$powers, ncol(res$powers)))

  stats <- calculate_test(
    res$metrics[c(lin_names, fpmax), ], n_obs,
    reduced_fit = fit_lin, full_fit = fit_fpmax
  )
  old_names <- names(res$statistic)
  res$statistic <- c(res$statistic, stats$statistic)
  names(res$statistic) <- c(old_names, sprintf("%s vs %s", fpmax, lin_names))
  res$pvalue <- c(res$pvalue, stats$pvalue)
  names(res$pvalue) <- names(res$statistic)

  if (mfp_pvalue_exceeds(stats$pvalue, alpha)) {
    # Test 2 not significant: the extra flexibility beyond linear is not
    # needed. Accept linear and stop (model_best = 3: rows were inserted in
    # order FPm(degree), null, linear).
    # no non-linearity detected
    res$power_best = 1
    res$model_best = 3
    res$current_adj_params <- fit_lin$current_adj_params
    return(res)
  }

  # Step 5: Test FPx vs. FPm, where x = 1,2,..m-1. Non-linearity was detected
  # (Test 2 significant), but FPm may still be more complex than necessary.
  # Walk up through increasingly complex lower degrees (FP1, FP2, ...) and
  # stop as soon as one is not significantly worse than FPm; that
  # lower-degree FP becomes the final choice. If every lower degree is
  # significantly worse, FPm itself is retained (Step 6, after the loop).
  # do this for all fps with lower degrees. dfs for tests are decreasing

  if (degree > 1) {
    # Vector of lower degrees starting with FP1.
    lower_degrees <- seq_len(degree - 1L)

    # Construct FP names once
    fp_names <- paste0("FP", lower_degrees, binary_suffix)

    for (i in seq_along(lower_degrees)) {
      current_degree <- lower_degrees[i]
      fpm <- fp_names[i]

      fit_fpm <- find_best_fpm_step(
        x = x, xi = xi, degree = current_degree, y = y,
        powers_current = powers_current, powers = powers, acdx = acdx,
        family = family, family_string = family_string, zero = zero, catzero = catzero,
        spike_decision = spike_decision, acd_parameter = acd_parameter, spike = spike,
        prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
        precomputed_adj = precomputed_adj, term_to_columns = term_to_columns,
        criterion = criterion,
        linear_fit = if (current_degree == 1L) fit_lin else NULL, ...
      )
      # Append metrics
      old_names = rownames(res$metrics)
      res$metrics <- rbind(
        res$metrics,
        fit_fpm$metrics[fit_fpm$model_best, ]
      )
      rownames(res$metrics) <- c(old_names, fpm)

      # Append powers, padded to match the fpmax column width (lower-degree
      # candidates have fewer powers than FPm(degree)).
      # Append powers
      res$powers <- rbind(res$powers,
                          ensure_length(fit_fpm$power_best, ncol(res$powers)))

      # Calculate test statistics
      stats <- calculate_test(
        res$metrics[c(fpm, fpmax), ], n_obs,
        reduced_fit = fit_fpm, full_fit = fit_fpmax
      )
      old_names <- names(res$statistic)
      res$statistic <- c(res$statistic, stats$statistic)
      names(res$statistic) <- c(old_names, sprintf("%s vs %s", fpmax, fpm))

      # Append p-values
      res$pvalue <- c(res$pvalue, stats$pvalue)
      names(res$pvalue) <- names(res$statistic)

      # Stop early if non-linearity detected
      if (mfp_pvalue_exceeds(stats$pvalue, alpha)) {
        # FPm(current_degree) is not significantly worse than FPm(degree):
        # accept it and stop climbing further; model_best is the row just
        # appended (the last row of res$metrics at this point).
        # non-linearity detected, but lower than maximum degree
        res$power_best = fit_fpm$powers[fit_fpm$model_best, , drop = TRUE]
        res$model_best = nrow(res$metrics)
        res$current_adj_params <- fit_fpm$current_adj_params
        return(res)
      }
    }

  }

  # Step 6: Every lower-degree FP tested was significantly worse than
  # FPm (or degree == 1, so there were no lower degrees to test at
  # all): retain the highest-degree FPm itself (row 1, inserted first).
  res$power_best <- fit_fpmax$powers[fit_fpmax$model_best, , drop = TRUE]
  res$model_best <- 1
  res$current_adj_params <- fit_fpmax$current_adj_params

  res
}

#' Select an ACD form using closed tests
#'
#' Used by \code{find_best_fp_step()} when \code{criterion = "pvalue"}
#' and ACD is active for \code{xi}. The procedure compares the full
#' \code{FP1(x, A(x))} form with simpler forms of \code{x} or \code{A(x)}.
#'
#' @details
#' The function prepares the adjustment terms once and keeps their current
#' forms throughout the comparisons. It performs up to five tests, stopping
#' when a form is selected:
#' \enumerate{
#'   \item Compare the full form with a model that excludes \code{xi}.
#'     Exclude \code{xi} if the p-value exceeds \code{select}, unless
#'     \code{xi} is in \code{keep}.
#'   \item Compare the full form with a fixed linear effect of \code{x}.
#'     Select the linear effect if the p-value exceeds \code{alpha}.
#'   \item Compare the full form with \code{FP1(x, .)}. Select that form
#'     if the p-value exceeds \code{alpha}.
#'   \item Compare the full form with \code{FP1(., A(x))}. Retain the
#'     full form unless the p-value exceeds \code{alpha}.
#'   \item If \code{FP1(., A(x))} remains eligible, compare it with a
#'     fixed linear effect of \code{A(x)}. Select the linear effect if
#'     the p-value exceeds \code{alpha}; otherwise select
#'     \code{FP1(., A(x))}.
#' }
#'
#' Exclusion or simplification requires a valid p-value strictly greater
#' than the relevant threshold. Thus \code{select = 1} prevents exclusion
#' and \code{alpha = 1} prevents simplification, including when a p-value
#' equals 1. Putting \code{xi} in \code{keep} prevents exclusion but still
#' allows a simpler form to be selected.
#'
#' For an eligible spike-at-zero term, this is the first selection stage.
#' The calling function handles any later choice between its continuous
#' and zero-indicator components.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{keep}: whether \code{xi} is in \code{keep}.
#'   \item \code{acd}: \code{TRUE}, indicating that ACD forms were assessed.
#'     The selected form need not contain an ACD component.
#'   \item \code{powers}: a two-column matrix describing the models reached
#'     before selection stops. Rows are added in this order:
#'     \code{FP1(x, A(x))}, null, linear \code{x}, \code{FP1(x, .)},
#'     \code{FP1(., A(x))}, and linear \code{A(x)}.
#'   \item \code{power_best}: the selected pair of powers. A missing entry
#'     can mark an unused component; two missing entries indicate exclusion.
#'   \item \code{metrics}: model scores in the same row order as
#'     \code{powers}.
#'   \item \code{model_best}: the selected row of \code{metrics}.
#'   \item \code{pvalue} and \code{statistic}: named results for each test
#'     performed before selection stopped.
#'   \item \code{is_gee}: whether the model uses GEE.
#'   \item \code{spike}: whether \code{xi} is treated as a spike-at-zero term.
#'   \item \code{current_adj_params}: information from the selected fit
#'     for reuse in later steps.
#'   \item \code{transform_cache}: reusable transformed adjustment terms.
#' }
#'
#' @param degree Not used directly. This procedure searches the fixed ACD
#'   forms described above.
#' @param ... Additional arguments passed through the fitting helpers to
#'   \code{fit_model()}.
#' @inheritParams find_best_fp_step
#'
#' @references
#' Royston, P. and Sauerbrei, W. (2016). mfpa: Extension of mfp using the
#' ACD covariate transformation for enhanced parametric multivariable
#' modeling. \emph{The Stata Journal}, 16(1), 72--87.
#'
#' @seealso \code{select_ra2()}
#' @keywords internal
#' @noRd
select_ra2_acd <- function(x,
                           xi,
                           keep,
                           degree,
                           acdx,
                           y,
                           powers_current,
                           powers,
                           criterion,
                           ftest,
                           select,
                           alpha,
                           family,
                           family_string,
                           zero,
                           catzero,
                           spike,
                           spike_decision,
                           acd_parameter,
                           prev_adj_params,
                           transform_cache = NULL,
                           force_max_fp,
                           has_offset,
                           n_obs,
                           term_to_columns,
                           ...) {

  is_gee <- mfp2_family_is_gee(family_string)

  # This implements the FSPA (function-selection procedure for ACD) 5-test
  # closed procedure for the 6 ACD sub-models (M1-M6), documented in detail in
  # mfp2()'s "Closed test procedure to choose final model" section:
  #   M1 = FP1(x, A(x))   [fpmax, the most flexible/reference model]
  #   M2 = FP1(x, .)      (regular FP1 in x, ACD term dropped)
  #   M3 = FP1(., A(x))   (FP1 in A(x), x term dropped)
  #   M4 = linear(x)      (linear in x, ACD term dropped)
  #   M5 = linear(., A(x)) (linear in A(x), x term dropped)
  #   M6 = null           (xi omitted entirely)
  # Tests run in the order M6 vs M1, M4 vs M1, M2 vs M1, M3 vs M1, M5 vs M3.

  # simplify testing by defining test helper function
  if (ftest) {
    calculate_test <- function(metrics, n_obs, reduced_fit = NULL, full_fit = NULL) {
      calculate_f_test(
        deviances = metrics[, "deviance_gaussian", drop = TRUE],
        dfs_resid = metrics[, "df_resid", drop = TRUE],
        n_obs = n_obs
      )
    }
  } else if (is_gee) {
    calculate_test <- function(metrics, n_obs, reduced_fit = NULL, full_fit = NULL) {
      calculate_deviance_test(
        deviances = metrics[, "deviance_rs", drop = TRUE],
        dfs = metrics[, "df", drop = TRUE]
      )
    }
  } else {
    calculate_test <- function(metrics, n_obs, reduced_fit = NULL, full_fit = NULL) {
      calculate_lr_test(
        logl = metrics[, "logl", drop = TRUE],
        dfs = metrics[, "df", drop = TRUE]
      )
    }
  }

  # Model labels get a binary suffix when the focal variable contributes a
  # structural-zero indicator, either through spike handling or catzero.
  # Compute this once to keep row names, test labels, and candidate names
  # consistent throughout the selection function.
  has_binary_xi <- isTRUE(spike[[xi]]) || !is.null(catzero[[xi]])
  binary_suffix <- if (has_binary_xi) " + Binary" else ""
  fpmax <- paste0("FP1(x, A(x))", binary_suffix)

  # acdx_reset_xi disables the ACD flag for xi only, used when fitting
  # candidates (e.g. linear(x), FP1(x, .)) that do not involve A(x).
  acdx_reset_xi <- acdx
  acdx_reset_xi[xi] = FALSE

  # output list
  res <- list(
    keep = xi %in% keep,
    acd = TRUE,
    powers = NULL,
    power_best = NULL,
    metrics = NULL,
    model_best = NULL,
    statistic = NULL,
    pvalue = NULL,
    is_gee = is_gee,
    spike = spike[xi],
    current_adj_params = NULL,
    transform_cache = NULL
  )

  # Build the adjustment matrix (all variables except xi, transformed
  # at their *current* powers_current) once, and reuse it for all candidate
  # fits below.
  # Thread the term lookup into adjustment construction so grouped terms are
  # included as complete blocks while xi undergoes ACD selection.
  precomputed_adj <- build_adjustment_step(
    x = x,
    xi = xi,
    powers_current = powers_current,
    powers = powers,
    acdx = acdx,
    zero = zero,
    catzero = catzero,
    spike = spike,
    spike_decision = spike_decision,
    acd_parameter = acd_parameter,
    prev_adj_params = prev_adj_params,
    transform_cache = transform_cache,
    term_to_columns = term_to_columns
  )
  res$transform_cache <- precomputed_adj$transform_cache

  # Test 1 (M6 vs M1): null vs. FP1(x, A(x)), df = 4. Asks whether xi (in its
  # most flexible ACD-augmented form) is associated with the outcome at all.
  # Not significant (at `select`) => drop xi entirely (M6), unless kept.

  # fit highest fp and null model for initial step
  fit_fpmax <- find_best_fpm_step(
    x = x, xi = xi, degree = 2, y = y,
    powers_current = powers_current, powers = powers, acdx = acdx,
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter, spike = spike,
    prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
    precomputed_adj = precomputed_adj, term_to_columns = term_to_columns,
    criterion = criterion, ...
  )

  fit_null <- fit_null_step(
    x = x, xi = xi, y = y,
    powers_current = powers_current, powers = powers, acdx = acdx,
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter, spike = spike,
    prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
    precomputed_adj = precomputed_adj, term_to_columns = term_to_columns, ...
  )
  # Row 1 = M1 (FP1(x, A(x))), the reference model for Tests 1-4; row 2 = M6
  # (null), whose power is a 2-column NA since neither x nor A(x) is included.
  res$metrics <- rbind(
    fit_fpmax$metrics[fit_fpmax$model_best, ],
    fit_null$metrics
  )
  rownames(res$metrics) <- c(fpmax, "null")
  res$powers <- rbind(fit_fpmax$power_best, fit_null$powers)

  res$current_adj_params <- fit_fpmax$current_adj_params

  # test for overall significance
  # df for tests are degree * 2 = 4
  stats <- calculate_test(
    res$metrics[c("null", fpmax), ], n_obs,
    reduced_fit = fit_null, full_fit = fit_fpmax
  )
  res$statistic <- stats$statistic
  names(res$statistic) <- sprintf("%s vs null", fpmax)
  res$pvalue <- stats$pvalue
  names(res$pvalue) <- names(res$statistic)

  # Ensure keep is not NULL
  current_keep <- if (is.null(keep)) character(0) else keep

  # Preserve the same strict MFP endpoint convention as select_ra2(): select = 1
  # forces inclusion even at p = 1, and alpha = 1 below prevents simplification
  # at the exact upper boundary.
  if (mfp_pvalue_exceeds(stats$pvalue, select) && !(xi %in% current_keep)) {
    # Test 1 not significant and xi not forced: eliminate xi (model_best = 2,
    # the null/M6 row).
    # not selected and not forced into model
    res$power_best = matrix(c(NA, NA), ncol = 2)
    res$model_best = 2
    res$current_adj_params <- fit_null$current_adj_params
    return(res)
  }

  # Test 2 (M4 vs M1): linear(x) vs. FP1(x, A(x)), df = 3. Asks whether the
  # simplest possible model (plain linear in x, no ACD term at all) is
  # already as good as the full M1 model. acdx_reset_xi is used here (and in
  # Test 3) because M4/M2 do not involve A(x) for xi.
  fit_lin <- fit_linear_step(
    x = x, xi = xi, y = y,
    powers_current = powers_current, powers = powers, acdx = acdx_reset_xi,
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter, spike = spike,
    prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
    precomputed_adj = precomputed_adj, term_to_columns = term_to_columns, ...
  )

  old_names = rownames(res$metrics)
  res$metrics <- rbind(
    res$metrics,
    fit_lin$metrics
  )

  lin_names <- paste0("linear", binary_suffix)
  rownames(res$metrics) <- c(old_names, lin_names)
  res$powers <- rbind(res$powers, c(1, NA))

  stats <- calculate_test(
    res$metrics[c(lin_names, fpmax), ], n_obs,
    reduced_fit = fit_lin, full_fit = fit_fpmax
  )
  old_names <- names(res$statistic)
  res$statistic <- c(res$statistic, stats$statistic)
  names(res$statistic) <- c(old_names, sprintf("%s vs %s", fpmax, lin_names))
  res$pvalue <- c(res$pvalue, stats$pvalue)
  names(res$pvalue) <- names(res$statistic)

  if (mfp_pvalue_exceeds(stats$pvalue, alpha)) {
    # Test 2 not significant: plain linear(x) (M4) is not significantly worse
    # than M1. Accept it and stop (model_best = 3: rows inserted in order M1,
    # null/M6, linear/M4).
    # no non-linearity detected
    res$power_best = matrix(c(1, NA), ncol = 2)
    res$model_best = 3
    res$current_adj_params <- fit_lin$current_adj_params
    return(res)
  }

  # Test 3 (M2 vs M1): FP1(x, .) vs. FP1(x, A(x)), df = 2. M4 was rejected
  # (Test 2 significant); now ask whether an ordinary (non-ACD) FP1 in x
  # alone captures the relationship as well as the full M1 model.
  # test for functional form, comparison with FP1(x, .)
  fit <- find_best_fpm_step(
    x = x, xi = xi, degree = 1, y = y,
    powers_current = powers_current, powers = powers, acdx = acdx_reset_xi,
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter, spike = spike,
    prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
    precomputed_adj = precomputed_adj, term_to_columns = term_to_columns,
    criterion = criterion, linear_fit = fit_lin, ...
  )

  old_names = rownames(res$metrics)
  res$metrics <- rbind(
    res$metrics,
    fit$metrics[fit$model_best, ]
  )

  FP_names <- paste0("FP1(x, .)", binary_suffix)

  rownames(res$metrics) <- c(old_names, FP_names)
  res$powers <- rbind(res$powers, c(fit$power_best, NA))

  stats <- calculate_test(
    res$metrics[c(FP_names, fpmax), ], n_obs,
    reduced_fit = fit, full_fit = fit_fpmax
  )
  old_names <- names(res$statistic)
  res$statistic <- c(res$statistic, stats$statistic)
  names(res$statistic) <- c(old_names, sprintf("%s vs %s", fpmax, FP_names))
  res$pvalue <- c(res$pvalue, stats$pvalue)
  names(res$pvalue) <- names(res$statistic)

  if (mfp_pvalue_exceeds(stats$pvalue, alpha)) {
    # Test 3 not significant: ordinary FP1(x, .) (M2) is not significantly
    # worse than M1. Accept it and stop (model_best = 4).
    # FP1(x, .) is good enough
    res$power_best = matrix(c(fit$power_best, NA), ncol = 2)
    res$model_best = 4
    res$current_adj_params <- fit$current_adj_params
    return(res)
  }

  # Test 4 (M3 vs M1): FP1(., A(x)) vs. FP1(x, A(x)), df = 2. M2 was rejected
  # too; now ask the complementary question - does dropping x and keeping
  # only the ACD term, FP1(., A(x)), capture the relationship as well as M1?
  #
  # a *significant* result now means M3 is rejected
  # and M1 must be kept as final (no further simplification possible), so we
  # stop with model_best = 1. A *non-significant* result means M3 survives
  # this test and we must go on to Test 5 to compare it against M5.
  # test for functional form, comparison with FP1(., A(x))
  fit_fp1a <- find_best_fpm_step(
    x = x, xi = xi, degree = 1, y = y,
    powers_current = powers_current, powers = powers, acdx = acdx,
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter, spike = spike,
    prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
    precomputed_adj = precomputed_adj, term_to_columns = term_to_columns,
    criterion = criterion, ...
  )

  old_names = rownames(res$metrics)
  res$metrics <- rbind(
    res$metrics,
    fit_fp1a$metrics[fit_fp1a$model_best, ]
  )

  FP1_names <- paste0("FP1(., A(x))", binary_suffix)
  rownames(res$metrics) <- c(old_names, FP1_names)
  res$powers <- rbind(res$powers, fit_fp1a$power_best)

  stats <- calculate_test(
    res$metrics[c(FP1_names, fpmax), ], n_obs,
    reduced_fit = fit_fp1a, full_fit = fit_fpmax
  )
  old_names <- names(res$statistic)
  res$statistic <- c(res$statistic, stats$statistic)
  names(res$statistic) <- c(old_names, sprintf("%s vs %s", fpmax, FP1_names))
  res$pvalue <- c(res$pvalue, stats$pvalue)
  names(res$pvalue) <- names(res$statistic)

  if (!mfp_pvalue_exceeds(stats$pvalue, alpha)) {
    # Test 4 significant at the package-wide inclusive significance boundary
    # (p <= alpha), or indeterminate because the p-value is invalid: M3 has not
    # supplied valid evidence for simplification. Retain the full M1 model.
    # FP1(x, A(x)) is the best
    res$power_best = fit_fpmax$power_best
    res$model_best = 1
    res$current_adj_params <- fit_fpmax$current_adj_params
    return(res)
  }

  # Test 5 (M5 vs M3): linear(., A(x)) vs. FP1(., A(x)), df = 1. M3 survived
  # Test 4 (was not significantly worse than M1), so now decide whether M3
  # itself can be simplified further down to a plain linear term in A(x).
  # Note this test compares against M3 (not M1) as the reference model -
  # the only test in this function that does so.
  # return best model between FP1(., A(x)) and linear(., A(x))
  fit_lineara <- fit_linear_step(
    x = x, xi = xi, y = y,
    powers_current = powers_current, powers = powers, acdx = acdx,
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter,spike = spike,
    prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
    fp1_fit = fit_fp1a,
    precomputed_adj = precomputed_adj, term_to_columns = term_to_columns, ...
  )

  old_names = rownames(res$metrics)
  res$metrics <- rbind(
    res$metrics,
    fit_lineara$metrics
  )

  linx_names <- paste0("linear(., A(x))", binary_suffix)
  rownames(res$metrics) <- c(old_names, linx_names)
  res$powers <- rbind(res$powers, fit_lineara$powers)

  stats <- calculate_test(
    res$metrics[c(linx_names, FP1_names), ], n_obs,
    reduced_fit = fit_lineara, full_fit = fit_fp1a
  )
  old_names <- names(res$statistic)
  res$statistic <- c(res$statistic, stats$statistic)
  names(res$statistic) <- c(old_names, sprintf("%s vs %s", FP1_names, linx_names))
  res$pvalue <- c(res$pvalue, stats$pvalue)
  names(res$pvalue) <- names(res$statistic)

  if (!mfp_pvalue_exceeds(stats$pvalue, alpha)) {
    # Test 5 significant at the inclusive boundary (p <= alpha), or
    # indeterminate because the p-value is invalid: M5 has not supplied valid
    # evidence for simplification. Retain M3, FP1(., A(x)).
    # use FP1(., A(x))
    res$power_best = fit_fp1a$power_best
    res$model_best = 5
    res$current_adj_params <- fit_fp1a$current_adj_params
    return(res)
  }

  # Test 5 not significant: M5 is not significantly worse than M3, so the
  # simplest surviving model is linear(., A(x)) - the end of the procedure,
  # since M5 has no simpler ACD sub-model to fall back to.
  # use linear(., A(x))
  res$power_best = matrix(c(NA, 1), ncol = 2)
  res$model_best = 6
  res$current_adj_params <- fit_lineara$current_adj_params

  res
}

#' Select a forced maximum FP form
#'
#' Used when \code{force_max_fp} is active for a term whose permitted form
#' is more complex than linear. The term must remain in the model at its
#' maximum permitted form under p-value, AIC, or BIC selection.
#'
#' @details
#' For an ordinary FP term, the function searches power combinations at the
#' requested \code{degree}. For an ACD term, it searches the full
#' \code{FP1(x, A(x))} form, represented internally by an ACD degree-2 search.
#' It fits the permitted power combinations within that form and selects the
#' best one. It does not fit null, linear, lower-degree FP, or reduced ACD
#' alternatives.
#'
#' The adjustment terms are prepared once and reused throughout the power
#' search. The returned \code{powers} and \code{metrics} each contain only
#' the winning candidate, even though multiple power combinations may have
#' been fitted.
#'
#' When an eligible spike-at-zero term is forced, the selected form retains
#' both its continuous effect and zero indicator. The caller therefore skips
#' the later comparisons that could remove either component.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{keep}: whether \code{xi} is also in \code{keep}.
#'   \item \code{acd}: whether ACD is active for \code{xi}.
#'   \item \code{powers} and \code{metrics}: one-row matrices for the
#'     winning power combination.
#'   \item \code{power_best}: the winning power or powers.
#'   \item \code{model_best}: \code{1}, the only returned row.
#'   \item \code{statistic} and \code{pvalue}: empty numeric vectors under
#'     p-value selection, or \code{NA} under AIC or BIC. No comparison
#'     with a simpler form is performed.
#'   \item \code{is_gee}: whether the model uses GEE.
#'   \item \code{spike}: whether \code{xi} is treated as a spike-at-zero term.
#'   \item \code{current_adj_params}: information from the winning fit
#'     for reuse in later steps.
#'   \item \code{transform_cache}: reusable transformed adjustment terms.
#' }
#'
#' @inheritParams find_best_fp_step
#' @param degree Requested ordinary FP degree. An ACD term instead uses
#'   its fixed maximum form.
#' @param ... Additional arguments passed through the fitting helper to
#'   \code{fit_model()}.
#' @keywords internal
#' @noRd
select_force_max_fp <- function(x,
                                xi,
                                keep,
                                degree,
                                acdx,
                                y,
                                powers_current,
                                powers,
                                criterion,
                                ftest,
                                select,
                                alpha,
                                family,
                                family_string,
                                zero,
                                catzero,
                                spike,
                                spike_decision,
                                acd_parameter,
                                prev_adj_params,
                                transform_cache = NULL,
                                force_max_fp,
                                has_offset,
                                n_obs,
                                term_to_columns,
                                ...) {

  # Keep the guard explicit because direct internal calls should never skip
  # model comparisons unless the caller has actually requested force_max_fp.
  # The criterion itself is deliberately unrestricted: forcing has the same
  # documented meaning for p-value, AIC, and BIC selection.
  if (!isTRUE(force_max_fp[xi])) {
    stop(
      "select_force_max_fp() requires force_max_fp = TRUE for `xi`.",
      call. = FALSE
    )
  }

  is_acd <- isTRUE(acdx[xi])
  has_binary_xi <- isTRUE(spike[[xi]]) || !is.null(catzero[[xi]])
  binary_suffix <- if (has_binary_xi) " + Binary" else ""

  # Adjustment-variable transformations are invariant across candidate forms.
  # Build/reuse them once, exactly as the ordinary IC selectors do.
  precomputed_adj <- build_adjustment_step(
    x = x,
    xi = xi,
    powers_current = powers_current,
    powers = powers,
    acdx = acdx,
    zero = zero,
    catzero = catzero,
    spike = spike,
    spike_decision = spike_decision,
    acd_parameter = acd_parameter,
    prev_adj_params = prev_adj_params,
    transform_cache = transform_cache,
    term_to_columns = term_to_columns
  )

  # ACD's maximum model is FP1(x, A(x)); internally it occupies two power
  # slots, so degree = 2 is the fixed search used to identify its best powers.
  # Ordinary FP uses the requested maximum degree directly. The package-wide
  # FP1 rule retains p = 1 whenever it is part of the supplied candidate set.
  forced_degree <- if (is_acd) 2 else degree

  fit_max <- find_best_fpm_step(
    x = x, xi = xi, degree = forced_degree, y = y,
    powers_current = powers_current, powers = powers, acdx = acdx,
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter, spike = spike,
    prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
    precomputed_adj = precomputed_adj, term_to_columns = term_to_columns,
    criterion = criterion, ...
  )

  # Keep only the winning fixed-form row.  No null, linear, lower-degree FP,
  # or reduced ACD models are fitted or materialized by this selector.
  model_label <- if (is_acd) {
    paste0("FP1(x, A(x))", binary_suffix)
  } else {
    paste0("FP", degree, binary_suffix)
  }

  powers_best <- fit_max$powers[fit_max$model_best, , drop = FALSE]
  metrics_best <- fit_max$metrics[fit_max$model_best, , drop = FALSE]
  rownames(powers_best) <- model_label
  rownames(metrics_best) <- model_label

  list(
    keep = xi %in% keep,
    acd = is_acd,
    powers = powers_best,
    power_best = fit_max$power_best,
    metrics = metrics_best,
    model_best = 1L,
    # No pairwise test is performed on the forced path. For p-value printing,
    # use a zero-length p-value vector so the one-row metrics table remains
    # one row rather than recycling an artificial NA comparison row.
    statistic = if (tolower(criterion) == "pvalue") numeric(0) else NA_real_,
    pvalue = if (tolower(criterion) == "pvalue") numeric(0) else NA_real_,
    is_gee = mfp2_family_is_gee(family_string),
    spike = spike[xi],
    current_adj_params = fit_max$current_adj_params,
    transform_cache = precomputed_adj$transform_cache
  )
}

#' Select an ordinary FP model using AIC or BIC
#'
#' Used by \code{find_best_fp_step()} when \code{criterion = "aic"} or
#' \code{"bic"} for a non-ACD term. Forced terms and eligible spike-at-zero
#' terms are handled by separate selectors.
#'
#' @details
#' The function compares a model excluding \code{xi}, a fixed linear model,
#' and the best FP model at each degree from 1 through \code{degree}. The
#' adjustment terms keep their current forms throughout these comparisons.
#'
#' \code{find_best_fpm_step()} searches the permitted power combinations
#' within each FP degree. Because combinations at the same degree have the
#' same selection penalty, it chooses the one with the highest likelihood,
#' or highest quasi-likelihood for GEE. This function then selects the form
#' with the lowest AIC or BIC. It does not use \code{select} or \code{alpha}.
#'
#' If \code{xi} is in \code{keep}, the null model is excluded from the final
#' comparison; the linear and FP forms remain eligible.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{keep}: whether \code{xi} is in \code{keep}.
#'   \item \code{acd}: \code{FALSE}.
#'   \item \code{powers}: a matrix describing the compared models, ordered
#'     null, linear, FP1, FP2, and so on through the requested degree.
#'   \item \code{power_best}: the selected row of \code{powers}. Its power
#'     is missing when the null model is selected.
#'   \item \code{metrics}: model scores in the same row order as
#'     \code{powers}.
#'   \item \code{model_best}: the selected row of \code{metrics}.
#'   \item \code{pvalue} and \code{statistic}: \code{NA}, because selection
#'     is based on AIC or BIC rather than a significance test.
#'   \item \code{is_gee}: whether the model uses GEE.
#'   \item \code{spike}: whether \code{xi} is marked as a spike-at-zero term.
#'   \item \code{current_adj_params}: information from the selected fit
#'     for reuse in later steps.
#'   \item \code{transform_cache}: reusable transformed adjustment terms.
#' }
#'
#' @param degree Maximum FP degree considered for \code{xi}.
#' @param ... Additional arguments passed through the fitting helpers to
#'   \code{fit_model()}.
#' @inheritParams find_best_fp_step
#' @seealso \code{select_ra2()}, \code{select_ic_acd()}
#' @keywords internal
#' @noRd
select_ic <- function(x,
                      xi,
                      keep,
                      degree,
                      acdx,
                      y,
                      powers_current,
                      powers,
                      criterion,
                      ftest,
                      select,
                      alpha,
                      family,
                      family_string,
                      zero,
                      catzero,
                      spike,
                      spike_decision,
                      acd_parameter,
                      prev_adj_params,
                      transform_cache = NULL,
                      force_max_fp,
                      has_offset,
                      n_obs,
                      term_to_columns,
                      ...) {

  # select_ic() implements ordinary (non-forced) AIC/BIC-based selection
  # (criterion = "aic"/"bic"), the simpler alternative to select_ra2()'s
  # sequential closed-test
  # procedure: instead of a chain of pairwise significance tests, every
  # relevant model (null, linear, FP1, FP2, ..., FPm) is fit once and
  # the one with the smallest AIC/BIC is chosen directly. There is no
  # `select`/`alpha` threshold involved here - the criterion value alone
  # decides, which is why `statistic`/`pvalue` below are just placeholders
  # (NA) rather than real test results.
  #
  # degree = 1 (linear-only variables) is handled by select_linear() instead;
  # this guard is defensive and should not normally trigger from
  # find_best_fp_step()'s dispatch logic.
  if (degree < 1) {
    return(NULL)
  }

  # Model labels get a binary suffix when the focal variable contributes a
  # structural-zero indicator, either through spike handling or catzero.
  # Compute this once to keep row names, test labels, and candidate names
  # consistent throughout the selection function.
  has_binary_xi <- isTRUE(spike[[xi]]) || !is.null(catzero[[xi]])
  binary_suffix <- if (has_binary_xi) " + Binary" else ""
  fpmax <- paste0("FP", degree, binary_suffix)

  # output list
  res <- list(
    keep = xi %in% keep,
    acd = FALSE,
    powers = NULL,
    power_best = NULL,
    metrics = NULL,
    model_best = NULL,
    statistic = NA,
    pvalue = NA,
    is_gee = mfp2_family_is_gee(family_string),
    spike = spike[xi],
    current_adj_params = NULL,
    transform_cache = NULL
  )


  # Thread the term lookup into adjustment construction so grouped terms are
  # included as complete blocks during AIC/BIC FP selection.
  precomputed_adj <- build_adjustment_step(
    x = x,
    xi = xi,
    powers_current = powers_current,
    powers = powers,
    acdx = acdx,
    zero = zero,
    catzero = catzero,
    spike = spike,
    spike_decision = spike_decision,
    acd_parameter = acd_parameter,
    prev_adj_params = prev_adj_params,
    transform_cache = transform_cache,
    term_to_columns = term_to_columns
  )
  res$transform_cache <- precomputed_adj$transform_cache
  # Step 1: Fit all relevant models (null, linear, FP1..FPm) ------------------
  # Unlike select_ra2(), there's no early-stopping here: every candidate is
  # fit up front and compared at the end, since the decision rule (smallest
  # AIC/BIC) needs all of them regardless of any individual model's fit.
  # Null Model
  fit_null <- fit_null_step(
    x = x, xi = xi, y = y,
    powers_current = powers_current, powers = powers, acdx = acdx,
    family = family,family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter, spike = spike,
    prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
    precomputed_adj = precomputed_adj, term_to_columns = term_to_columns,
    gee_selection_criterion = criterion, ...
  )
  # Linear model
  fit_lin <- fit_linear_step(
    x = x, xi = xi, y = y,
    powers_current = powers_current, powers = powers, acdx = acdx,
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter, spike = spike,
    prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
    precomputed_adj = precomputed_adj, term_to_columns = term_to_columns,
    gee_selection_criterion = criterion, ...
  )

  res$current_adj_params <- fit_null$current_adj_params


  # Build the focal numerical basis once at the maximum requested degree.
  # select_ic() necessarily evaluates every degree, so this basis will be reused
  # and there is no RA2-style early-stopping tradeoff. Build from the complete
  # user-supplied power set, including power 1; every degree-specific view uses
  # the same candidate semantics.
  xi_cols <- term_to_columns[[xi]]
  if (is.null(xi_cols) || length(xi_cols) != 1L) {
    stop(
      "Internal error: IC FP basis reuse requires one numeric focal column.",
      call. = FALSE
    )
  }

  focal_basis_cache <- build_shared_focal_fp_basis(
    x = x[, xi_cols, drop = TRUE],
    max_degree = degree,
    powers = powers[[xi]],
    zero = zero[xi],
    catzero = catzero[[xi]]
  )

  # All FPm models: fit the best FP1, FP2...FPm candidate at each
  # degree in turn (find_best_fpm_step() already picks the best power
  # combination within each fixed degree by deviance, which is equivalent to
  # AIC/BIC minimisation at fixed df since the penalty term is the same for
  # all power combinations sharing that degree).
  fits_fpm <- vector("list", degree)

  for (m in seq_len(degree)) {
    fits_fpm[[m]] <- find_best_fpm_step(
      x = x, xi = xi, degree = m, y = y,
      powers_current = powers_current, powers = powers, acdx = acdx,
      family = family, family_string = family_string, zero = zero, catzero = catzero,
      spike_decision = spike_decision,acd_parameter = acd_parameter,spike = spike,
      prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
      precomputed_adj = precomputed_adj, focal_basis_cache = focal_basis_cache,
      term_to_columns = term_to_columns, criterion = criterion,
      linear_fit = if (m == 1L) fit_lin else NULL, ...
    )
  }

  names(fits_fpm) <- paste0("FP", seq_len(degree), binary_suffix)

  # Step 2: Assemble powers/metrics for null, linear, and each best FPm into a
  # single indexed summary (candidate_names gives the row order/labels) -----
  # output summary - only output best fpm models
  len_max <- ncol(fits_fpm[[fpmax]]$powers)

  candidate_names <- c(
    "null",
    paste0("linear", binary_suffix),
    names(fits_fpm)
  )

  res$powers <- lapply(fits_fpm, function(x) {
    ensure_length(x$powers[x$model_best, , drop = FALSE], len_max)
  })

  res$powers <- do.call(
    rbind,
    c(
      list(
        ensure_length(fit_null$powers, len_max),
        ensure_length(fit_lin$powers, len_max)
      ),
      res$powers
    )
  )
  rownames(res$powers) <- candidate_names

  res$metrics <- lapply(fits_fpm, function(x) {
    x$metrics[x$model_best, , drop = FALSE]
  })

  res$metrics <- do.call(
    rbind,
    c(
      list(
        fit_null$metrics,
        fit_lin$metrics
      ),
      res$metrics
    )
  )
  rownames(res$metrics) <- candidate_names

  # Step 3: Select the best model by AIC/BIC -----------------------------
  # force_max_fp is handled before this function is dispatched, by
  # select_force_max_fp(). Consequently every model fitted above is a genuine
  # competitor here; the only special restriction is whether `keep` excludes
  # the null model from consideration.
  if (xi %in% keep) {
    # keep: xi must remain in the model, so the null-model row (row 1) is
    # excluded from the comparison before taking which.min(); the resulting
    # index is then shifted by +1 to account for the excluded row and land
    # back on the correct row of the full res$metrics/res$powers.
    # Prevent selection of null model; choose best among linear through FPm.
    # Build the full safe row sequence and drop row 1 rather than using 2:n,
    # which can create a descending sequence when bounds change unexpectedly.
    ind_select <- seq_len(nrow(res$metrics))[-1L]
    res$model_best <- which.min(
      res$metrics[ind_select, tolower(criterion), drop = TRUE]
    )
    # shift by 1 since null row was excluded
    res$model_best <- res$model_best + 1
  } else {
    # Unrestricted: every candidate (including null) is eligible; simply
    # take whichever row has the smallest AIC/BIC.
    # Unrestricted: choose best among null through FPm.
    ind_select <- seq_len(nrow(res$metrics))
    res$model_best <- which.min(
      res$metrics[ind_select, tolower(criterion), drop = TRUE]
    )
  }

  # Map the winning row index back to the adjustment-parameter cache from
  # whichever underlying fit produced it: row 1 = null, row 2 = linear, rows
  # 3+ = fits_fpm[[row - 2]] (FP1, FP2, ... in the same order they were
  # appended to res$metrics/res$powers above).
  res$power_best <- res$powers[res$model_best, , drop = FALSE]
  res$current_adj_params <- if (res$model_best == 1L) {
    fit_null$current_adj_params
  } else if (res$model_best == 2L) {
    fit_lin$current_adj_params
  } else {
    fits_fpm[[res$model_best - 2L]]$current_adj_params
  }

  res
}

#' Select an ACD model using AIC or BIC
#'
#' Used by \code{find_best_fp_step()} when ACD is active for \code{xi}
#' and \code{criterion = "aic"} or \code{"bic"}.
#'
#' @details
#' The function compares six forms: null, linear \code{x}, linear
#' \code{A(x)}, \code{FP1(x, .)}, \code{FP1(., A(x))}, and
#' \code{FP1(x, A(x))}. For each FP form, it first searches the permitted
#' power combinations. It then selects the form with the lowest AIC or BIC.
#' The adjustment terms keep their current forms throughout the search.
#'
#' If \code{xi} is in \code{keep}, the null form is excluded from the final
#' comparison. The function does not use \code{select} or \code{alpha}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{keep}: whether \code{xi} is in \code{keep}.
#'   \item \code{acd}: \code{TRUE}, indicating that ACD forms were assessed;
#'     the selected form need not contain an ACD component.
#'   \item \code{powers}: a two-column matrix describing the six forms in
#'     the order listed above.
#'   \item \code{power_best}: the selected row of \code{powers}. Missing
#'     entries mark components absent from that form.
#'   \item \code{metrics}: model scores in the same row order as
#'     \code{powers}.
#'   \item \code{model_best}: the selected row of \code{metrics}.
#'   \item \code{pvalue} and \code{statistic}: \code{NA}, because selection
#'     is based on AIC or BIC.
#'   \item \code{is_gee}: whether the model uses GEE.
#'   \item \code{spike}: whether \code{xi} is marked as a spike-at-zero term.
#'   \item \code{current_adj_params}: information from the selected fit
#'     for reuse in later steps.
#'   \item \code{transform_cache}: reusable transformed adjustment terms.
#' }
#'
#' @param degree Not used directly; the function compares its six fixed
#'   candidate forms.
#' @param ... Additional arguments passed through the fitting helpers to
#'   \code{fit_model()}.
#' @inheritParams find_best_fp_step
#' @seealso \code{select_ic()}
#' @keywords internal
#' @noRd
select_ic_acd <- function(x,
                          xi,
                          keep,
                          degree,
                          acdx,
                          y,
                          powers_current,
                          powers,
                          criterion,
                          ftest,
                          select,
                          alpha,
                          family,
                          family_string,
                          zero,
                          catzero,
                          spike,
                          spike_decision,
                          acd_parameter,
                          prev_adj_params,
                          transform_cache = NULL,
                          force_max_fp,
                          has_offset,
                          n_obs,
                          term_to_columns,
                          ...) {

  # select_ic_acd() is the ordinary (non-forced) AIC/BIC counterpart of
  # select_ra2_acd(): instead of
  # the 5-test closed sequence over M1-M6, it fits the 6 relevant ACD
  # sub-models directly (null, linear(x), linear(., A(x)), FP1(x, .),
  # FP1(., A(x)), FP1(x, A(x))) and picks whichever has the smallest AIC/BIC.
  # As in select_ic(), `statistic`/`pvalue` are just NA placeholders since no
  # significance testing is involved.
  acdx_reset_xi <- acdx
  acdx_reset_xi[xi] <- FALSE
  # Model labels get a binary suffix when the focal variable contributes a
  # structural-zero indicator, either through spike handling or catzero.
  # Compute this once to keep row names, test labels, and candidate names
  # consistent throughout the selection function.
  has_binary_xi <- isTRUE(spike[[xi]]) || !is.null(catzero[[xi]])
  binary_suffix <- if (has_binary_xi) " + Binary" else ""

  # output list
  res <- list(
    keep = xi %in% keep,
    acd = TRUE,
    powers = NULL,
    power_best = NULL,
    metrics = NULL,
    model_best = NULL,
    statistic = NA,
    pvalue = NA,
    is_gee = mfp2_family_is_gee(family_string),
    spike = spike[xi],
    current_adj_params = NULL,
    transform_cache = NULL
  )


  # Thread the term lookup into adjustment construction so grouped terms are
  # included as complete blocks during AIC/BIC ACD selection.
  precomputed_adj <- build_adjustment_step(
    x = x,
    xi = xi,
    powers_current = powers_current,
    powers = powers,
    acdx = acdx,
    zero = zero,
    catzero = catzero,
    spike = spike,
    spike_decision = spike_decision,
    acd_parameter = acd_parameter,
    prev_adj_params = prev_adj_params,
    transform_cache = transform_cache,
    term_to_columns = term_to_columns
  )
  res$transform_cache <- precomputed_adj$transform_cache
  # Step 1: Fit all relevant models: null, linear(x), linear(., A(x)), and the
  # three best FP1 candidates (FP1(x, .), FP1(., A(x)), FP1(x, A(x))) --------
  # As in select_ic(), every candidate is fit up front (no early stopping),
  # since AIC/BIC selection needs all of them for comparison. acdx_reset_xi
  # (ACD disabled for xi only) is used for the null/linear(x)/FP1(x, .)
  # candidates, which do not involve A(x); the ordinary acdx is used for the
  # linear(., A(x))/FP1(., A(x))/FP1(x, A(x)) candidates, which do.
  # fit all relevant models
  fit_null <- fit_null_step(
    x = x, xi = xi, y = y,
    powers_current = powers_current, powers = powers, acdx = acdx,
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter, spike = spike,
    prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
    precomputed_adj = precomputed_adj, term_to_columns = term_to_columns,
    gee_selection_criterion = criterion, ...
  )

  # linear(x, .)
  fit_lin <- fit_linear_step(
    x = x, xi = xi, y = y,
    powers_current = powers_current, powers = powers, acdx = acdx_reset_xi,
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter, spike = spike,
    prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
    precomputed_adj = precomputed_adj, term_to_columns = term_to_columns,
    gee_selection_criterion = criterion, ...
  )

  # linear(., A(x))
  fit_lina <- fit_linear_step(
    x = x, xi = xi, y = y,
    powers_current = powers_current, powers = powers, acdx = acdx,
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter,spike = spike,
    prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
    precomputed_adj = precomputed_adj, term_to_columns = term_to_columns,
    gee_selection_criterion = criterion, ...
  )

  res$current_adj_params <- fit_null$current_adj_params

  # Build one joint degree-2 x/A(x) numerical basis for all three nonlinear ACD
  # competitors below. Its x columns are also valid for FP1(x, .), while its
  # A(x) columns serve FP1(., A(x)); the joint candidate uses one of each. The
  # cache is constructed from this variable's actual allowed powers, so custom
  # and variable-specific power sets do not rely on default column positions.
  xi_cols <- term_to_columns[[xi]]
  if (is.null(xi_cols) || length(xi_cols) != 1L) {
    stop(
      "Internal error: IC ACD basis reuse requires one numeric focal column.",
      call. = FALSE
    )
  }

  acd_parameter_xi <- acd_parameter[[xi]]
  acd_training_values_xi <- if (!is.null(acd_parameter_xi)) {
    acd_parameter_xi[["acd", exact = TRUE]]
  } else {
    NULL
  }

  focal_basis_cache <- build_shared_focal_acd_basis(
    x = x[, xi_cols, drop = TRUE],
    powers = powers[[xi]],
    zero = zero[xi],
    catzero = catzero[[xi]],
    acd_parameter = acd_parameter_xi,
    acd_training_values = acd_training_values_xi
  )

  # The three FP1 variants (in x only, in A(x) only, and jointly FP1(x, A(x)))
  # are computed via find_best_fpm_step() at the appropriate degree (degree 1
  # for the single-transform variants, degree 2 for the joint FP1(x, A(x))
  # since it involves two power slots) and collected in a named list so Step 2
  # can iterate over them generically.


  fits <- setNames(
    list(
      # FP1(x, .)
      find_best_fpm_step(
        x = x, xi = xi, degree = 1, y = y,
        powers_current = powers_current, powers = powers, acdx = acdx_reset_xi,
        family = family, family_string = family_string, zero = zero, catzero = catzero,
        spike_decision = spike_decision,acd_parameter = acd_parameter,spike = spike,
        prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
        precomputed_adj = precomputed_adj, focal_basis_cache = focal_basis_cache,
        term_to_columns = term_to_columns, criterion = criterion,
        linear_fit = fit_lin, ...
      ),
      # FP1(., A(x))
      find_best_fpm_step(
        x = x, xi = xi, degree = 1, y = y,
        powers_current = powers_current, powers = powers, acdx = acdx,
        family = family, family_string = family_string, zero = zero, catzero = catzero,
        spike_decision = spike_decision, acd_parameter = acd_parameter,spike = spike,
        prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
        precomputed_adj = precomputed_adj, focal_basis_cache = focal_basis_cache,
        term_to_columns = term_to_columns, criterion = criterion,
        linear_fit = fit_lina, ...
      ),
      # FP1(x, A(x))
      find_best_fpm_step(
        x = x, xi = xi, degree = 2, y = y,
        powers_current = powers_current, powers = powers, acdx = acdx,
        family = family, family_string = family_string, zero = zero, catzero = catzero,
        spike_decision = spike_decision,acd_parameter = acd_parameter,spike = spike,
        prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
        precomputed_adj = precomputed_adj, focal_basis_cache = focal_basis_cache,
        term_to_columns = term_to_columns, criterion = criterion, ...
      )
    ),
    c(
      paste0("FP1(x, .)", binary_suffix),
      paste0("FP1(., A(x))", binary_suffix),
      paste0("FP1(x, A(x))", binary_suffix)
    )
  )

  # Step 2: Assemble powers/metrics for null, the two linear variants, and the
  # three best FP1 candidates into a single indexed summary -----------------
  # Assemble output summary ----------------------------------------------------
  # output summary - only output best fpm models
  len_max <- 2L

  candidate_names <- c(
    "null",
    paste0("linear", binary_suffix),
    paste0("linear(., A(x))", binary_suffix),
    names(fits)
  )

  res$powers <- lapply(fits, function(x) {
    ensure_length(x$powers[x$model_best, , drop = FALSE], len_max)
  })

  res$powers <- do.call(
    rbind,
    c(
      list(
        ensure_length(fit_null$powers, len_max),
        ensure_length(fit_lin$powers, len_max),
        ensure_length(fit_lina$powers, len_max)
      ),
      res$powers
    )
  )
  rownames(res$powers) <- candidate_names

  res$metrics <- lapply(fits, function(x) {
    x$metrics[x$model_best, , drop = FALSE]
  })

  res$metrics <- do.call(
    rbind,
    c(
      list(
        fit_null$metrics,
        fit_lin$metrics,
        fit_lina$metrics
      ),
      res$metrics
    )
  )
  rownames(res$metrics) <- candidate_names
  # Step 3: Select the best model by AIC/BIC -----------------------------
  # force_max_fp is handled by select_force_max_fp() before this function is
  # dispatched. All six models above therefore remain genuine competitors;
  # `keep` only determines whether the null row is eligible.
  if (xi %in% keep) {
    # keep: exclude the null row, take which.min() over the rest, then
    # shift the index by +1 to land back on the correct row of the full
    # res$metrics/res$powers.
    # Prevent selection of null model; choose best among linear through FP1(x, A(x)).
    # As above, avoid a programmatic 2:n sequence.
    ind_select <- seq_len(nrow(res$metrics))[-1L]
    res$model_best <- which.min(
      res$metrics[ind_select, tolower(criterion), drop = TRUE]
    )
    res$model_best <- res$model_best + 1L
  } else {
    # Unrestricted: every candidate (including null) is eligible.
    # Unrestricted: choose best among null through FP1(x, A(x))
    ind_select = seq_len(nrow(res$metrics))
    res$model_best <- which.min(
      res$metrics[ind_select, tolower(criterion), drop = TRUE]
    )
  }

  # Map the winning row index back to the adjustment-parameter cache from
  # whichever underlying fit produced it: row 1 = null, row 2 = linear(x),
  # row 3 = linear(., A(x)), rows 4+ = fits[[row - 3]] (FP1(x,.), FP1(.,A(x)),
  # FP1(x,A(x)), in the same order they were inserted into `fits` above).
  res$power_best = res$powers[res$model_best, , drop = FALSE]
  res$current_adj_params <- if (res$model_best == 1L) {
    fit_null$current_adj_params
  } else if (res$model_best == 2L) {
    fit_lin$current_adj_params
  } else if (res$model_best == 3L) {
    fit_lina$current_adj_params
  } else {
    fits[[res$model_best - 3L]]$current_adj_params
  }

  res
}

#' Build adjustment data for single-column terms in C++
#'
#' Prepares aligned inputs for \code{build_adjustment_step_loop_cpp()}.
#' Within \code{build_adjustment_step()}, this path handles adjustment terms
#' that each map to one column with the same name. Terms mapped to different
#' column names or to several columns are assembled separately.
#'
#' The C++ function reuses valid transformed blocks where possible. Otherwise,
#' it transforms the current term, applies any stored ACD parameters, and
#' includes a zero indicator when required. It does not estimate or refit
#' ACD parameters.
#'
#' @param x Numeric working design matrix, or a numeric data frame that can
#'   be converted to a matrix. Its columns must include \code{vars_adj}.
#' @param vars_adj Names of the adjustment terms, in the order used for the
#'   returned data.
#' @param powers_adj Named list of current selected powers for those terms.
#' @param acdx_adj Named logical vector or list indicating active ACD terms.
#' @param zero_adj Named logical vector or list indicating exact-zero handling.
#' @param catzero Named list of zero-indicator matrices. Each relevant entry
#'   is \code{NULL} or a one-column matrix with \code{nrow(x)} rows.
#' @param spike_adj Named logical vector or list indicating spike-at-zero terms.
#' @param spike_decision_int_adj Named integer vector of current spike decisions.
#' @param acd_parameter_adj Named list of stored ACD parameters. An active
#'   ACD adjustment term must have a non-\code{NULL} entry.
#' @param eliminated Named logical vector indicating terms whose current
#'   powers are all missing.
#' @param spike_binary_only_flags Named logical vector identifying spike terms
#'   represented only by their zero indicator. Such terms are handled before
#'   the ordinary eliminated-term check.
#' @param current_power_keys_adj Named list of normalized current power keys
#'   used to check whether stored transformations can be reused.
#' @param prev_power_keys_adj Named list of stored power keys, or \code{NULL}
#'   when no previous keys are available.
#' @param prev_xi Previous cache data containing transformed blocks and spike
#'   decisions, or \code{NULL}. This may come from the focal term's earlier
#'   step or from the shared per-term cache.
#' @param has_prev Whether \code{prev_xi} contains cache data to check.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{data_adj_list}: a named list with one transformed matrix per
#'     adjustment term. A term contributing no columns has a matrix with
#'     zero columns.
#'   \item \code{data_adj}: those matrices joined by columns. It is
#'     \code{NULL} when there are no adjustment terms, or a matrix with
#'     zero columns when terms exist but none contributes columns.
#' }
#' @keywords internal
#' @noRd
mfp2_build_adjustment_step_loop <- function(x,
                                            vars_adj,
                                            powers_adj,
                                            acdx_adj,
                                            zero_adj,
                                            catzero,
                                            spike_adj,
                                            spike_decision_int_adj,
                                            acd_parameter_adj,
                                            eliminated,
                                            spike_binary_only_flags,
                                            current_power_keys_adj,
                                            prev_power_keys_adj = NULL,
                                            prev_xi = NULL,
                                            has_prev = FALSE) {
  # Step 1: Handle the trivial case of no adjustment variables --------------
  # The zero-adjustment-variable case is handled here and returns data_adj = NULL.
  # If adjustment variables exist but all contribute zero columns, the C++ helper
  # returns an n x 0 matrix.
  if (length(vars_adj) == 0L) {
    return(list(
      data_adj_list = list(),
      data_adj      = NULL
    ))
  }

  # Step 2: Normalize/validate inputs for the C++ hot loop -------------------
  # The C++ interface expects a numeric matrix. In normal MFP internals `x`
  # should already be a numeric matrix, but this keeps the bridge robust.
  if (!is.matrix(x)) {
    x <- data.matrix(x)
  }

  if (!is.numeric(x)) {
    stop("Internal error: `x` must be numeric in build_adjustment_step().",
         call. = FALSE)
  }

  # Convert list-like logical inputs to plain aligned logical vectors. This
  # avoids relying on Rcpp's coercion of named one-element list entries.
  acdx_flag <- vapply(acdx_adj, isTRUE, logical(1L))
  zero_flag <- vapply(zero_adj, isTRUE, logical(1L))
  spike_flag <- vapply(spike_adj, isTRUE, logical(1L))
  eliminated_flag <- vapply(eliminated, isTRUE, logical(1L))
  spike_binary_only_flag <- vapply(spike_binary_only_flags, isTRUE, logical(1L))

  # build_adjustment_step() should only apply stored ACD parameters. It should
  # not fit/refit ACD inside the MFP step loop.
  if (any(acdx_flag)) {
    missing_acd_parameter <- vars_adj[
      acdx_flag &
        vapply(acd_parameter_adj, is.null, logical(1L))
    ]

    if (length(missing_acd_parameter) > 0L) {
      stop(
        paste0(
          "Internal error: ACD adjustment variable(s) are missing stored ",
          "acd_parameter values in build_adjustment_step(): ",
          paste(missing_acd_parameter, collapse = ", "),
          "."
        ),
        call. = FALSE
      )
    }
  }

  # Step 3: Prepare previous-cycle cache inputs (for cache-hit reuse) -------
  # Previous per-variable matrices are used for cache hits. If no previous cache
  # exists, pass aligned NULL placeholders.
  prev_data_adj_list <- vector("list", length(vars_adj))
  names(prev_data_adj_list) <- vars_adj

  prev_spike_decision_int_adj <- rep(NA_integer_, length(vars_adj))
  names(prev_spike_decision_int_adj) <- vars_adj

  if (has_prev) {
    if (is.null(prev_xi)) {
      stop("Internal error: `has_prev = TRUE` but `prev_xi` is NULL.",
           call. = FALSE)
    }

    prev_data_adj_list <- prev_xi$data_adj_list[vars_adj]

    prev_spike_decision_int_adj <- as.integer(
      prev_xi$spike_decision_adj[vars_adj]
    )
    names(prev_spike_decision_int_adj) <- vars_adj
  }

  # If there is no previous cache, pass an aligned list of NULL power keys.
  if (!has_prev || is.null(prev_power_keys_adj)) {
    prev_power_keys_adj <- vector("list", length(vars_adj))
    names(prev_power_keys_adj) <- vars_adj
  }

  # Pass column positions instead of slicing x. C++ uses zero-based indices.
  x_col_index <- match(vars_adj, colnames(x))

  if (anyNA(x_col_index)) {
    stop(
      paste0(
        "Internal error: adjustment variable(s) are missing from `x`: ",
        paste(vars_adj[is.na(x_col_index)], collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  # Slice catzero here so the C++ inputs are all aligned to vars_adj.
  catzero_adj <- catzero[vars_adj]

  # Step 4: Call the C++ hot loop with all inputs aligned to vars_adj -------
  cpp_adj <- build_adjustment_step_loop_cpp(
    x                            = x,
    x_col_index                  = as.integer(x_col_index - 1L),
    vars_adj                     = vars_adj,
    powers_adj                   = unname(powers_adj),
    acdx_adj                     = unname(acdx_flag),
    zero_adj                     = unname(zero_flag),
    catzero_adj                  = unname(catzero_adj),
    spike_adj                    = unname(spike_flag),
    spike_decision_int_adj       = unname(as.integer(spike_decision_int_adj)),
    acd_parameter_adj            = unname(acd_parameter_adj),
    eliminated                   = unname(eliminated_flag),
    spike_binary_only_flags      = unname(spike_binary_only_flag),
    current_power_keys_adj       = unname(current_power_keys_adj),
    prev_power_keys_adj          = unname(prev_power_keys_adj),
    prev_data_adj_list           = unname(prev_data_adj_list),
    prev_spike_decision_int_adj  = unname(prev_spike_decision_int_adj),
    has_prev                     = has_prev
  )

  # Step 5: Restore names on the returned per-variable matrices --------------
  # Restore names defensively. The C++ helper also assigns names, but keeping
  # this in R makes the cache contract explicit.
  data_adj_list <- cpp_adj$data_adj_list
  names(data_adj_list) <- vars_adj

  list(
    data_adj_list = data_adj_list,
    data_adj      = cpp_adj[["data_adj", exact = TRUE]]
  )
}

#' Build transformed adjustment data for one MFP step
#'
#' Builds the adjustment data used when assessing the focal term \code{xi}.
#' The adjustment terms are all current terms other than \code{xi}; each uses
#' its current powers and spike-at-zero decision.
#'
#' @details
#' The function can reuse an unchanged transformed block from the shared
#' \code{transform_cache} or from an earlier step for the same focal term.
#' It rebuilds blocks whose powers or spike decisions have changed, then
#' joins the blocks into the adjustment matrix for this step. The complete
#' adjustment matrix is assembled for the current step rather than retained
#' in the cache between cycles.
#'
#' Terms mapped to one column with the same name use the C++ adjustment
#' helper. Terms represented by a mapped column or a group of columns are
#' handled as fixed linear blocks.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{powers_adj}: named list of current powers for adjustment
#'     terms, or \code{NULL} if there are no adjustment terms.
#'   \item \code{spike_decision_adj}: named vector of their current spike
#'     decisions, or \code{NULL} if there are no adjustment terms.
#'   \item \code{data_adj_list}: named list with one transformed matrix per
#'     adjustment term, or an empty list if there are none. A term that
#'     contributes no columns has a matrix with zero columns.
#'   \item \code{data_adj}: the matrices in \code{data_adj_list} joined by
#'     columns. It is \code{NULL} if there are no adjustment terms, or an
#'     \code{n}-row matrix with zero columns if terms exist but none
#'     contributes columns.
#'   \item \code{transform_cache}: updated shared cache of transformed
#'     blocks for reuse in later steps.
#' }
#'
#' @inheritParams transform_data_step
#' @keywords internal
#' @noRd
build_adjustment_step <- function(x,
                                  xi,
                                  powers_current,
                                  powers,
                                  acdx,
                                  zero,
                                  catzero,
                                  spike,
                                  spike_decision,
                                  acd_parameter,
                                  prev_adj_params,
                                  transform_cache = NULL,
                                  term_to_columns) {
  # `term_to_columns` has already been normalized by fit_mfp(); this function
  # only consumes the mapping when assembling singleton and grouped blocks.
  # ---------------------------------------------------------------------------
  # Purpose
  # ---------------------------------------------------------------------------
  # Build the adjustment matrix for the current focal variable `xi`.
  #
  # In the MFP cycle, each variable is updated conditional on all other variables.
  # For focal variable `xi`, the adjustment variables are all variables except
  # `xi`. Their transformed columns are collected into `data_adj`.
  #
  # Two cache layers are deliberately kept separate:
  #   * prev_adj_params[[xi]] is focal-variable-specific and retains metadata
  #     plus per-variable transformed blocks as a fallback on a later cycle.
  #     It intentionally does not retain the complete assembled data_adj matrix.
  #   * transform_cache[[v]] is keyed by adjustment variable and can reuse v's
  #     transformed block across different focal variables in the same fit.
  #
  # The shared cache is an ordinary list threaded explicitly through the MFP
  # call chain. It is created once per fit_mfp() call, so there is no mutable
  # global/package state and no environment-based side effect.
  #
  # Important internal invariant:
  #   catzero[[varname]] is either NULL or an n x 1 numeric/integer matrix.
  #
  # It is NOT a vector inside this inner fitting code.

  # Step 1: Identify the adjustment variable set (everything except xi) -----
  # Use powers_current names for the canonical variable order.
  # Accessing x columns by name avoids copying/reordering the full x matrix.
  names_powers_current <- names(powers_current)

  # Normalize the run-scoped per-variable cache. Keeping all conceptual terms
  # in a named list makes ownership explicit and lets R share unchanged matrix
  # objects under copy-on-modify semantics.
  if (is.null(transform_cache)) {
    transform_cache <- setNames(
      vector("list", length(names_powers_current)),
      names_powers_current
    )
  } else {
    if (!is.list(transform_cache) || is.null(names(transform_cache))) {
      stop(
        "Internal error: `transform_cache` must be a named list.",
        call. = FALSE
      )
    }

    missing_cache_names <- setdiff(names_powers_current, names(transform_cache))
    if (length(missing_cache_names) > 0L) {
      transform_cache[missing_cache_names] <- vector(
        "list",
        length(missing_cache_names)
      )
    }

    transform_cache <- transform_cache[names_powers_current]
  }

  # Adjustment variables are all variables except the current focal variable.
  vars_adj <- setdiff(names_powers_current, xi)

  # Return contract:
  #   * If there are no adjustment variables, data_adj is NULL.
  #   * If adjustment variables exist but all contribute zero columns,
  #     data_adj is an n x 0 matrix.
  #   * If at least one adjustment variable contributes columns,
  #     data_adj is the column-bound adjustment matrix.
  #
  # Downstream code should test NCOL(data_adj) > 0L rather than relying on NULL
  # versus n x 0 matrix semantics.

  # These are returned even when there are no adjustment variables.
  powers_adj <- NULL
  spike_decision_adj <- NULL

  # Defaults for the no-adjustment-variable case.
  data_adj <- NULL
  data_adj_list <- list()

  if (length(vars_adj) > 0L) {
    # Step 2: Slice adjustment-variable metadata and compute cache keys -----
    # Previous cached adjustment information for this focal variable.
    # prev_adj_params is indexed by xi.
    prev_xi <- prev_adj_params[[xi]]
    has_prev <- !is.null(prev_xi)

    # Slice all adjustment-variable metadata once.
    # This avoids repeated indexing into the full objects inside the loop.
    powers_adj <- powers_current[vars_adj]

    # A variable is considered eliminated when all stored powers are NA.
    # Eliminated variables contribute zero adjustment columns unless they are
    # binary-only spike variables, which are handled before this branch.
    eliminated <- vapply(
      powers_adj,
      function(p) all(is.na(as.numeric(p))),
      logical(1L)
    )

    acdx_adj <- acdx[vars_adj]
    zero_adj <- zero[vars_adj]
    acd_parameter_adj <- acd_parameter[vars_adj]
    spike_adj <- spike[vars_adj]
    spike_decision_adj <- spike_decision[vars_adj]

    # Convert spike decisions once to integer.
    # This is important because identical(3, 3L) is FALSE in R.
    # Cache comparisons below use identical(), so current and previous spike
    # decisions must have the same type.
    spike_decision_int_adj <- as.integer(spike_decision_adj)
    names(spike_decision_int_adj) <- vars_adj

    # -------------------------------------------------------------------------
    # Precompute binary-only spike flags
    # -------------------------------------------------------------------------
    # spike_decision == 3 means the continuous FP part is ignored and the
    # variable is represented only by the structural-zero indicator.
    #
    # This branch must be checked before the eliminated/all-NA-powers branch,
    # because for binary-only spike variables the FP powers are intentionally
    # irrelevant.
    spike_binary_only_flags <- vapply(
      vars_adj,
      function(v) {
        isTRUE(spike_adj[[v]]) &&
          identical(
            spike_decision_int_adj[[v]],
            saz_decision_codes[["binary_only"]]
          )
      },
      logical(1L)
    )
    names(spike_binary_only_flags) <- vars_adj

    # -------------------------------------------------------------------------
    # Precompute normalized power keys for cache comparison
    # -------------------------------------------------------------------------
    # For ordinary variables, the cache key is the current FP powers.
    #
    # For binary-only spike variables, the FP powers are irrelevant. The helper
    # normalize_powers_for_convergence() converts their power key to NA_real_
    # whenever spike_decision == 3. That prevents unnecessary recomputation if
    # only the ignored continuous powers changed.
    current_power_keys_adj <- normalize_powers_for_convergence(
      powers_adj,
      spike_decision_int_adj
    )

    # Previous normalized power keys are needed only if a previous cache exists.
    prev_power_keys_adj <- NULL

    if (has_prev) {
      # Keep previous spike decisions as integer for consistent comparison.
      prev_spike_decision_int_adj <- as.integer(
        prev_xi$spike_decision_adj[vars_adj]
      )
      names(prev_spike_decision_int_adj) <- vars_adj

      prev_power_keys_adj <- normalize_powers_for_convergence(
        prev_xi$powers_adj[vars_adj],
        prev_spike_decision_int_adj
      )
    }

    # -------------------------------------------------------------------------
    # No focal-level whole-matrix cache
    # -------------------------------------------------------------------------
    # Older versions retained prev_xi$data_adj and returned that assembled
    # n x adjustment matrix directly when every adjustment transformation was
    # unchanged. Although this avoided one matrix assembly on a full cache hit,
    # keeping one such matrix for every focal variable makes retained memory
    # grow approximately as O(n * p^2). The persistent cache now keeps only
    # per-variable blocks/metadata, and the adjustment matrix is assembled once
    # for the current focal evaluation from reusable blocks below.

    # -------------------------------------------------------------------------
    # Shared per-variable cache: reuse blocks across different focal variables
    # -------------------------------------------------------------------------
    # Prefer transform_cache[[v]], which follows v across focal-variable
    # evaluations. If an entry is not available, fall back to the historical
    # prev_adj_params[[xi]] per-variable block so direct/internal callers retain
    # the old cache behavior as well. The C++ loop performs the actual validity
    # check using normalized power keys plus spike decisions.
    reuse_power_keys_adj <- setNames(
      vector("list", length(vars_adj)),
      vars_adj
    )
    reuse_data_adj_list <- setNames(
      vector("list", length(vars_adj)),
      vars_adj
    )
    reuse_spike_decision_adj <- setNames(
      rep(NA_integer_, length(vars_adj)),
      vars_adj
    )

    for (v in vars_adj) {
      cached <- transform_cache[[v]]

      if (
        !is.null(cached) &&
        !is.null(cached$data) &&
        !is.null(cached$power_key) &&
        length(cached$spike_decision) == 1L
      ) {
        # Use single-bracket list assignment for the key so the aligned list
        # cannot shrink if a malformed cache entry ever contains NULL.
        reuse_power_keys_adj[v] <- list(cached$power_key)
        reuse_data_adj_list[[v]] <- cached$data
        reuse_spike_decision_adj[[v]] <- as.integer(cached$spike_decision)
      } else if (has_prev && !is.null(prev_xi$data_adj_list[[v]])) {
        reuse_power_keys_adj[v] <- list(prev_power_keys_adj[[v]])
        reuse_data_adj_list[[v]] <- prev_xi$data_adj_list[[v]]
        reuse_spike_decision_adj[[v]] <- prev_spike_decision_int_adj[[v]]
      }
    }

    has_reuse_cache <- any(!vapply(reuse_data_adj_list, is.null, logical(1L)))
    reuse_cache <- if (has_reuse_cache) {
      list(
        data_adj_list = reuse_data_adj_list,
        spike_decision_adj = reuse_spike_decision_adj
      )
    } else {
      NULL
    }

    # Step 3: Preserve the historical C++ path only for identity-mapped
    # singleton terms. A categorical term may have one dummy column but a
    # different conceptual name, for example `svi` -> `svi1`; such terms must
    # be assembled through their raw-column lookup rather than matched to x by
    # the conceptual name.
    mapped_as_block <- vapply(
      vars_adj,
      function(term) {
        cols <- term_to_columns[[term]]
        length(cols) != 1L || !identical(cols[[1L]], term)
      },
      logical(1L)
    )
    names(mapped_as_block) <- vars_adj

    if (!any(mapped_as_block)) {
      cpp_adj <- mfp2_build_adjustment_step_loop(
        x                         = x,
        vars_adj                  = vars_adj,
        powers_adj                = powers_adj,
        acdx_adj                  = acdx_adj,
        zero_adj                  = zero_adj,
        catzero                   = catzero,
        spike_adj                 = spike_adj,
        spike_decision_int_adj    = spike_decision_int_adj,
        acd_parameter_adj         = acd_parameter_adj,
        eliminated                = eliminated,
        spike_binary_only_flags   = spike_binary_only_flags,
        current_power_keys_adj    = current_power_keys_adj,
        prev_power_keys_adj       = reuse_power_keys_adj,
        prev_xi                   = reuse_cache,
        has_prev                  = has_reuse_cache
      )

      data_adj_list <- cpp_adj$data_adj_list
      data_adj      <- cpp_adj[["data_adj", exact = TRUE]]
    } else {
      direct_vars <- vars_adj[!mapped_as_block]
      block_vars  <- vars_adj[mapped_as_block]

      data_adj_list <- vector("list", length(vars_adj))
      names(data_adj_list) <- vars_adj

      if (length(direct_vars) > 0L) {
        cpp_adj <- mfp2_build_adjustment_step_loop(
          x                         = x,
          vars_adj                  = direct_vars,
          powers_adj                = powers_adj[direct_vars],
          acdx_adj                  = acdx_adj[direct_vars],
          zero_adj                  = zero_adj[direct_vars],
          catzero                   = catzero,
          spike_adj                 = spike_adj[direct_vars],
          spike_decision_int_adj    = spike_decision_int_adj[direct_vars],
          acd_parameter_adj         = acd_parameter_adj[direct_vars],
          eliminated                = eliminated[direct_vars],
          spike_binary_only_flags   = spike_binary_only_flags[direct_vars],
          current_power_keys_adj    = current_power_keys_adj[direct_vars],
          prev_power_keys_adj       = reuse_power_keys_adj[direct_vars],
          prev_xi                   = reuse_cache,
          has_prev                  = has_reuse_cache
        )
        data_adj_list[direct_vars] <- cpp_adj$data_adj_list[direct_vars]
      }

      # Fixed categorical blocks are already on their required linear design
      # scale. Reuse a shared block when its key and spike state still match;
      # otherwise rebuild it from the mapped raw columns. This covers both
      # multi-column factors and binary factors whose single dummy column has
      # a different name from the conceptual term.
      for (term in block_vars) {
        block_cache_hit <-
          !is.null(reuse_data_adj_list[[term]]) &&
          identical(
            reuse_power_keys_adj[[term]],
            current_power_keys_adj[[term]]
          ) &&
          !is.na(reuse_spike_decision_adj[[term]]) &&
          reuse_spike_decision_adj[[term]] == spike_decision_int_adj[[term]]

        if (block_cache_hit) {
          data_adj_list[[term]] <- reuse_data_adj_list[[term]]
        } else if (isTRUE(eliminated[[term]])) {
          data_adj_list[[term]] <- matrix(
            numeric(0L),
            nrow = nrow(x),
            ncol = 0L,
            dimnames = list(rownames(x), character(0L))
          )
        } else {
          data_adj_list[[term]] <- x[, term_to_columns[[term]], drop = FALSE]
        }
      }

      contributes <- vapply(data_adj_list, NCOL, integer(1L)) > 0L
      if (any(contributes)) {
        data_adj <- do.call(cbind, data_adj_list[contributes])
      } else {
        data_adj <- matrix(
          numeric(0L),
          nrow = nrow(x),
          ncol = 0L,
          dimnames = list(rownames(x), character(0L))
        )
      }
    }
  }

  # Store the current per-variable blocks after assembly. These entries are
  # independent of the focal variable, so the next focal evaluation can reuse
  # them whenever the normalized power key and spike decision still match.
  if (length(vars_adj) > 0L) {
    for (v in vars_adj) {
      transform_cache[[v]] <- list(
        power_key = current_power_keys_adj[[v]],
        spike_decision = spike_decision_int_adj[[v]],
        data = data_adj_list[[v]]
      )
    }
  }

  list(
    powers_adj = powers_adj,
    spike_decision_adj = spike_decision_adj,
    data_adj_list = data_adj_list,
    data_adj = data_adj,
    transform_cache = transform_cache
  )
}

#' Prepare candidate forms for a term and its adjustment terms
#'
#' Prepares the candidate forms of the focal term \code{xi} and the current
#' forms of all other terms. Adjustment terms may use FP or ACD
#' transformations, fixed linear blocks, or zero indicators. Unchanged
#' adjustment-term blocks can be reused from earlier steps.
#'
#' @param x Working design matrix containing the columns needed for the
#'   focal and adjustment terms. Ordinary FP variables have already been
#'   shifted and scaled; active ACD variables use scale 1.
#' @param xi Name of the focal conceptual term.
#' @param powers_current Named list of current powers for all conceptual terms.
#' @param df Numeric scalar controlling the candidate form of \code{xi}.
#'   A value of 1 requests its fixed linear form; FP searches use twice
#'   the requested FP degree.
#' @param powers Named list of permitted FP powers for each term.
#' @param acdx Named logical vector indicating active ACD terms.
#' @param zero Named logical vector indicating exact-zero handling.
#' @param catzero Named list of exact-zero indicator matrices. Each entry
#'   is \code{NULL} or a one-column matrix with \code{nrow(x)} rows.
#' @param spike Named logical vector indicating spike-at-zero terms.
#' @param spike_decision Named vector of spike-at-zero decisions for the
#'   conceptual terms. For an adjustment term, \code{1} includes its
#'   transformed and zero-indicator components, \code{2} includes only
#'   its transformed component, and \code{3} includes only its indicator.
#' @param acd_parameter Named list of stored ACD parameters, with one
#'   entry per conceptual term.
#' @param prev_adj_params Earlier adjustment information that may allow
#'   unchanged transformed blocks to be reused.
#' @param precomputed_adj Optional adjustment data from
#'   \code{build_adjustment_step()}. If supplied, the function reuses it
#'   instead of building the adjustment data again.
#' @param focal_basis_cache Optional shared numerical basis for \code{xi},
#'   used to prepare candidates at the requested degree without repeating
#'   the same transformations.
#' @param term_to_columns Named list mapping each conceptual term to its
#'   column or columns in \code{x}.
#' @param compact_fp Whether to represent ordinary FP candidates with
#'   a shared basis and candidate map instead of separate full matrices.
#' @param compact_acd Whether to use the corresponding compact
#'   representation for ACD candidates.
#'
#' @details
#' The function first builds or reuses the adjustment data, then prepares
#' the candidates for \code{xi}. Terms represented by several columns and
#' variables with at most three distinct values use a fixed, materialized
#' form. For an eligible numeric FP or ACD search, compact mode stores
#' shared transformed columns and a map describing each candidate.
#'
#' In compact mode, \code{data_fp} is \code{NULL}; the candidates are
#' represented by \code{fp_basis} or \code{acd_basis}. Otherwise,
#' \code{data_fp} contains a list of candidate matrices.
#'
#' With \code{df = 1}, an ordinary numeric term uses its fixed linear form.
#' For an active ACD term, the fixed form uses \code{A(x)}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{powers_fp}: a matrix describing the candidate powers
#'     for \code{xi}.
#'   \item \code{data_fp}: a list of candidate matrices for \code{xi},
#'     or \code{NULL} when a compact basis is used.
#'   \item \code{fp_basis}: the compact ordinary-FP basis and candidate
#'     map, or \code{NULL}.
#'   \item \code{acd_basis}: the compact ACD basis and candidate map,
#'     or \code{NULL}.
#'   \item \code{powers_adj}: named list of current powers for adjustment
#'     terms, or \code{NULL} if there are none.
#'   \item \code{data_adj}: assembled adjustment matrix; \code{NULL} if
#'     there are no adjustment terms. It may have zero columns if
#'     adjustment terms exist but none contributes a column.
#'   \item \code{current_params}: information stored under \code{xi}
#'     for reuse, including adjustment powers, spike decisions,
#'     per-term transformed blocks, and the assembled adjustment matrix.
#' }
#' @keywords internal
#' @noRd
transform_data_step <- function(x,
                                xi,
                                powers_current,
                                df,
                                powers,
                                acdx,
                                zero,
                                catzero,
                                spike,
                                spike_decision,
                                acd_parameter,
                                prev_adj_params,
                                precomputed_adj = NULL,
                                focal_basis_cache = NULL,
                                term_to_columns,
                                compact_fp = FALSE,
                                compact_acd = FALSE
) {
  # `term_to_columns` is supplied by the fitting call chain and maps the focal
  # conceptual term to the raw column block used below.
  # Step 1: Build (or reuse) the adjustment matrix for every variable except xi
  # Access x columns by name; avoid copying/reordering x for every candidate
  # transformation call.
  if (is.null(precomputed_adj)) {
    adj <- build_adjustment_step(
      x = x,
      xi = xi,
      powers_current = powers_current,
      powers = powers,
      acdx = acdx,
      zero = zero,
      catzero = catzero,
      spike = spike,
      spike_decision = spike_decision,
      acd_parameter = acd_parameter,
      prev_adj_params = prev_adj_params,
      term_to_columns = term_to_columns
    )
  } else {
    adj <- precomputed_adj
  }

  # Step 2: Generate the FP/ACD candidate transformations for xi itself -----
  # `fp_basis` and `acd_basis` are compact hot-loop representations. Linear,
  # categorical, low-cardinality, and callers that do not request compact mode
  # retain the historical materialized `data_fp` representation.
  fp_basis <- NULL
  acd_basis <- NULL

  xi_cols <- term_to_columns[[xi]]
  if (is.null(xi_cols)) {
    stop(
      sprintf("Internal error: term '%s' is missing from `term_to_columns`.", xi),
      call. = FALSE
    )
  }

  if (length(xi_cols) > 1L) {
    # A multi-column term is categorical by construction. Treat its complete
    # design block as one fixed linear candidate and skip all FP/cardinality
    # logic that applies only to a single numeric predictor.
    data_fp <- list(x[, xi_cols, drop = FALSE])
    powers_fp <- matrix(1, nrow = 1L, ncol = length(xi_cols))
  } else {
    # Preserve the historical single-column path exactly.
    data_xi <- x[, xi_cols, drop = TRUE]

    # A selector-scoped focal basis is valid only for this exact training
    # vector and structural-zero representation. Keep these checks O(1): the
    # optimization is intended to remove repeated n-length transformations, so
    # we deliberately do not rescan/compare every cached value here.
    if (!is.null(focal_basis_cache)) {
      cache_has_catzero <- !is.null(focal_basis_cache$catzero)
      step_has_catzero <- !is.null(catzero[[xi]])

      if (is.null(focal_basis_cache$basis) ||
          !is.matrix(focal_basis_cache$basis) ||
          nrow(focal_basis_cache$basis) != length(data_xi) ||
          !identical(isTRUE(focal_basis_cache$zero), isTRUE(zero[xi])) ||
          cache_has_catzero != step_has_catzero) {
        stop(
          paste0(
            "Internal error: `focal_basis_cache` does not match the current ",
            "focal training data/zero handling."
          ),
          call. = FALSE
        )
      }
    }

    if (length(unique(data_xi)) <= 3) {
      # Variables with <= 3 distinct values are treated as linear/binary: no FP
      # candidates are generated. The original nonnegative values are used as
      # supplied, and a catzero indicator column is added if present (df = 1 is
      # enforced upstream for such variables; see assign_df()).

      # add catzero variable if not null; if null cbind will remove it
      data_xi_mat <- matrix(data_xi, ncol = 1L)

      if (!is.null(catzero[[xi]])) {
        data_xi_mat <- cbind(catzero[[xi]], data_xi_mat)
        colnames(data_xi_mat)[1L] <- "catzero"
      }

      data_fp <- list(data_xi_mat)
      powers_fp <- matrix(1, nrow = 1L, ncol = 1L)

    } else {
      if (acdx[xi]) {
        # fit_mfp() estimates ACD once per variable and retains the fitted A(x)
        # vector in acd_parameter[[xi]]$acd for the training sample. Pass that
        # vector explicitly to the candidate generator so it can be reused
        # without calling apply_acd() again. Legacy/direct callers whose stored
        # parameter list has no $acd value simply fall back to the old apply path.
        acd_parameter_xi <- acd_parameter[[xi]]
        acd_training_values_xi <- if (!is.null(acd_parameter_xi)) {
          acd_parameter_xi[["acd", exact = TRUE]]
        } else {
          NULL
        }

        # ACD searches can use the same compact-basis strategy as ordinary FP:
        # store each unique transform of x and A(x) once and keep only a small
        # integer candidate map. Degree 0/1 naturally has one active A(x) term.
        if (isTRUE(compact_acd)) {
          if (!is.null(focal_basis_cache)) {
            # IC selectors may already hold the joint degree-2 x/A(x) basis.
            # Generate only this degree's ACD candidate specification and map
            # it onto those existing numerical columns. Custom power sets and
            # candidate ordering remain controlled by generate_powers_acd().
            fpd <- view_shared_focal_acd_basis(
              shared_basis = focal_basis_cache,
              degree = floor(df / 2),
              powers = powers[[xi]]
            )
          } else {
            fpd <- generate_transformations_acd_basis(
              data_xi,
              degree = floor(df / 2),
              powers = powers[[xi]],
              zero = zero[xi],
              catzero = catzero[[xi]],
              acd_parameter = acd_parameter_xi,
              acd_training_values = acd_training_values_xi
            )
          }
          acd_basis <- fpd
          data_fp <- NULL
        } else {
          fpd <- generate_transformations_acd(
            data_xi,
            degree = floor(df / 2),
            powers = powers[[xi]],
            zero = zero[xi],
            catzero = catzero[[xi]],
            acd_parameter = acd_parameter_xi,
            acd_training_values = acd_training_values_xi
          )
          data_fp <- fpd$data
        }
      } else {
        # note that degree is df / 2
        if (isTRUE(compact_fp)) {
          if (!is.null(focal_basis_cache)) {
            # Reuse the selector's maximum-degree numerical basis. The current
            # call still generates its own degree-specific candidates, including
            # power 1 for FP1, and maps them by explicit power/repetition
            # metadata rather than fixed column positions.
            fpd <- view_shared_focal_fp_basis(
              shared_basis = focal_basis_cache,
              degree = floor(df / 2),
              powers = powers[[xi]]
            )
          } else {
            fpd <- generate_transformations_fp_basis(
              data_xi,
              degree = floor(df / 2),
              powers = powers[[xi]],
              zero = zero[xi],
              catzero = catzero[[xi]]
            )
          }
          fp_basis <- fpd
          data_fp <- NULL
        } else {
          fpd <- generate_transformations_fp(
            data_xi,
            degree = floor(df / 2),
            powers = powers[[xi]],
            zero = zero[xi],
            catzero = catzero[[xi]]
          )
          data_fp <- fpd$data
        }
      }

      powers_fp <- fpd$powers
    }
  }

  # Step 3: Cache this step's adjustment-variable state for the next step ---
  # Store everything under xi in current_params
  current_params <- list()
  current_params[[xi]] <- list(
    powers_adj = adj$powers_adj,
    spike_decision_adj = adj$spike_decision_adj,
    data_adj_list = adj$data_adj_list,
    data_adj = adj[["data_adj", exact = TRUE]]
  )

  # Return results and current parameters for next step
  list(
    powers_fp = powers_fp,
    data_fp = data_fp,
    fp_basis = fp_basis,
    acd_basis = acd_basis,
    powers_adj = adj$powers_adj,
    data_adj = adj[["data_adj", exact = TRUE]],
    current_params = current_params
  )
}

#' Pad a vector to a fixed length
#'
#' Pads `x` with `fill` until it has length `size`.
#'
#' Used internally to make sure matrix rows have matching dimensions. This
#' helper assumes that `length(x) <= size`; it is intended for padding, not
#' truncation.
#'
#' @param x Input vector.
#' @param size Desired length of `x`.
#' @param fill Value used to pad `x` when it is shorter than `size`.
#'   Defaults to `NA`.
#'
#' @return `x` unchanged if `length(x) == size`; otherwise `x` followed by
#'   enough `fill` values to reach length `size`.
#'
#' @keywords internal
#' @noRd
ensure_length <- function(x, size, fill = NA) {
  if (length(x) == size) {
    return(x)
  }

  x_new <- rep(fill, size)
  x_new[seq_along(x)] <- x
  x_new
}
