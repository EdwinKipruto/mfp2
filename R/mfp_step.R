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


#' Function to estimate the best FP functions for a single variable
#'
#' See \code{mfp2()} for a brief summary on the notation used here and
#' \code{fit_mfp()} for an overview of the fitting procedure.
#'
#' @param x Numeric design matrix with one row per observation and raw columns
#'   already expanded as needed. It excludes the intercept and is assumed to be
#'   shifted and scaled.
#' @param y a vector for the response variable or a `Surv` object.
#' @param xi Character scalar naming the current conceptual term. For a
#'   singleton continuous term its FP form is assessed; for a multi-column term
#'   the complete fixed linear block is assessed jointly.
#' @param term_to_columns Named list mapping each conceptual term to the raw
#'   design-matrix columns that represent it. Multi-column entries are treated
#'   as fixed linear blocks and tested jointly. The lookup is threaded
#'   explicitly through every selection routine so grouped adjustment
#'   terms are available during both linear and FP model searches.
#' @param weights a vector of observation weights of length nobs.
#' @param offset a vector of length nobs of offsets.
#' @param df a numeric vector indicating the maximum degrees of freedom for the
#' variable of interest `xi`.
#' @param powers_current Named list with one selected-power state per
#'   conceptual term, used for all adjustment terms in the current step.
#' @param family Either a character string naming the family (e.g., "gaussian", "binomial", "cox")
#'   or a function that returns a GLM family object (e.g., stats::gaussian).
#'   For Cox models, only a character string "cox" is allowed.
#' @param family_string A character string representing the selected family,
#'   e.g., "gaussian".
#' @param criterion a character string defining the criterion used to select
#' variables and FP models of different degrees.
#' @param select Numeric scalar giving the nominal significance level used to
#'   decide whether the current term `xi` is retained during MFP backfitting.
#'   For a grouped categorical term, this value applies to the joint test of
#'   its complete design-matrix block. A value of 1 forces the term to remain.
#' @param alpha Numeric significance level for tests between FP models of
#'   different degrees for `xi` and, for an eligible spike-at-zero term, for
#'   the Stage-2 component-removal tests.
#' @param keep a character vector with names of variables to be kept
#' in the model.
#' @param powers a named list of numeric values that sets the permitted FP
#' powers for each covariate.
#' @param method a character string specifying the method for tie handling in
#' Cox regression.
#' @param strata a factor of all possible combinations of stratification
#' variables. Returned from [survival::strata()].
#' @param nocenter a numeric vector with a list of values for fitting Cox
#' models. See [survival::coxph()] for details.
#' @param acdx Named logical vector with one value per conceptual term;
#'   multi-column terms must be \code{FALSE}.
#' @param ftest a logical indicating the use of the F-test for Gaussian models.
#' @param control a list with parameters for model fit.
#' @param rownames a parameter for Cox models.
#' @param catzero A named list of exact-zero indicators. Each element is either
#' `NULL` or an n x 1 integer/numeric matrix containing `I(x == 0)`. Non-NULL
#' elements indicate variables for which that binary column is available.
#' @param zero Named logical vector with one value per conceptual term,
#'   indicating exact-zero handling for nonnegative singleton terms. Negative
#'   values are rejected before this function is called.
#' @param spike Named logical vector with one value per conceptual term,
#'   indicating spike-at-zero handling for nonnegative singleton terms, with
#'   `x == 0` defining the zero component and `x > 0` the positive component.
#' @param acd_parameter Named list of ACD parameters produced by \code{fit_acd()},
#' with length equal to \code{ncol(x)}. Each list element corresponds to a variable;
#' if an element is \code{NULL}, the variable was not specified in the
#' \code{acdx} argument of \code{fit_mfp}.
#' @param spike_decision Named vector indicating how spike-at-zero (SAZ)
#' variables are handled. Each element corresponds to a variable and encodes
#' the selected strategy: `1` = include FP for positive values plus binary SAZ,
#' `2` = treat as continuous FP only, `3` = include binary SAZ only.
#' @param prev_adj_params Named list storing adjustment metadata and reusable
#' per-variable transformed blocks from previous steps for each focal variable.
#' Complete assembled adjustment matrices are not retained between cycles.
#' @param transform_cache Named list keyed by variable name. Each populated
#'   entry stores that variable's normalized power key, spike decision, and
#'   transformed adjustment block so it can be reused across focal variables
#'   within one fit_mfp() call.
#' @param force_max_fp A logical vector of length \code{nvars}, named by
#'   variable name. If \code{TRUE} for a non-linear variable \code{xi},
#'   \code{select_force_max_fp()} fits only the most complex functional form
#'   allowed by \code{df}: the highest-likelihood FP model at the requested
#'   degree for ordinary FP variables, or \code{FP1(x, A(x))} for ACD
#'   variables. Null, linear, and lower-degree alternatives are not fitted
#'   because variable selection and functional-form simplification are
#'   explicitly bypassed for forced terms under p-value, AIC, and BIC
#'   selection. For an eligible spike-at-zero term, forcing applies to the
#'   complete maximum SAZ representation: the maximum continuous component
#'   plus its binary zero indicator. SAZ stage 2 is therefore skipped because
#'   its only purpose is to compare reduced component representations. The best
#'   power combination within the forced form is still determined by
#'   \code{find_best_fpm_step()}.
#' @param retain_linear_fp1 Internal logical used specifically for MFPI's
#'   forced fixed-degree power-selection fits. If \code{TRUE}, a degree-1 FP
#'   search retains power \code{1} when it is present in the supplied candidate
#'   set. It does not add power \code{1} when absent. Ordinary MFP closed
#'   testing keeps the default \code{FALSE}, fitting the linear model separately
#'   and excluding power 1 from its non-linear FP1 candidate search.
#' @param has_offset logical indicating whether an offset was specified before
#' missing offsets were replaced by zeros internally.
#' @param n_obs Numeric; number of observations (or observed events, for Cox
#' models), used for AIC/BIC-based selection and small-sample F-tests.
#' @param precomputed_adj Optional internal adjustment object returned by
#' \code{build_adjustment_step()}. When supplied, adjustment-variable
#' transformations are reused and only current-variable candidate
#' transformations are generated. \code{find_best_fp_step()} itself does not
#' take this argument; it is documented here because several sibling
#' functions below (\code{find_best_fpm_step()}, \code{fit_null_step()},
#' \code{fit_linear_step()}, \code{select_linear()}, \code{select_ra2()},
#' \code{select_ra2_acd()}, \code{select_force_max_fp()}, \code{select_ic()},
#' \code{select_ic_acd()}) use
#' \code{@@inheritParams find_best_fp_step} and rely on this entry.
#' @param focal_basis_cache Optional selector-scoped numerical basis for the
#' focal variable. Information-criterion selectors build this once at their
#' maximum required degree and pass it to lower-degree searches. Candidate
#' powers/order are still generated separately for each degree; only the
#' n-length transformed basis columns are reused. This cache is internal and
#' must refer to the same training observations, structural-zero handling, and
#' variable-specific allowed powers as the current focal term.
#' @param verbose a logical; run in verbose mode.
#'
#' @details
#' The function selection procedure (FSP) is used if the p-value criterion is
#' chosen, whereas the criteria AIC and BIC select the model with the smallest
#' AIC and BIC, respectively.
#'
#' It uses transformations for all other variables to assess the FP form of
#' the current variable of interest. This function covers three main use cases:
#'
#' * the linear case (`df = 1`) to test between null and linear models (see
#' \code{select_linear()}). This step differs from the mfp case because
#' linear models only use 1 df, while estimation of (every) fp power adds
#' another df. This is also the case applied for categorical variables for
#' which `df` are set to 1.
#' * the case that an acd transformation is requested (`acdx` is `TRUE`
#' for `xi`) for the variable of interest (see \code{find_best_fpm_step()}).
#' * the (usual) case of the normal mfp algorithm to assess non-linear
#' functional forms (see \code{find_best_fpm_step()}).
#'
#' Note that these cases do not encompass the setting that a variable is not
#' selected, because the evaluation is done for each variable in each cycle.
#' A variable which was de-selected in earlier cycles may be added to the
#' working model again. Also see \code{find_best_fp_cycle()}.
#'
#' The adjustment in each step uses the current fp powers given in
#' `powers_current` for all other variables to determine the adjustment set
#' and transformations in the  working model.
#'
#' Note that the algorithm starts by setting all `df = 1`, and higher fps
#' are evaluated in turn starting from the first step in the first cycle.
#'
#' @section Functional form selection:
#' There are 3 criteria to decide for the current best functional form of a
#' continuous variable.
#'
#' The first option for `criterion = "pvalue"` is the function selection
#' procedure as outlined in e.g. Chapters 4 and 6 of Royston and
#' Sauerbrei (2008), also abbreviated as "RA2".
#' It is a closed testing procedure and is implemented in \code{select_ra2()} and
#' extended for ACD transformation in \code{select_ra2_acd()} according to
#' Royston and Sauerbrei (2016).
#'
#' For the other criteria `aic` and `bic`, ordinary selection fits FP models
#' up to the desired degree and chooses the model with the lowest information
#' criterion via \code{select_ic()} (or \code{select_ic_acd()} for ACD).
#' Across all three criteria, `force_max_fp` is handled before the ordinary
#' selectors are dispatched: \code{select_force_max_fp()} fits only the
#' predetermined maximum functional form and therefore avoids model fits that
#' cannot affect the result. For an eligible spike-at-zero variable, the forced
#' result is the complete maximum SAZ representation (maximum continuous form
#' plus binary zero indicator), so the reduced-component comparisons of SAZ
#' stage 2 are also bypassed.
#'
#' @return
#' A numeric vector indicating the best powers for `xi`. Entries can be
#' `NA` if variable is to be removed from the working model. Note that this
#' vector may include up to two `NA` entries when ACD transformation is
#' requested, but otherwise is either a vector with all numeric entries, or a
#' single `NA`.
#'
#' @references
#' Royston, P. and Sauerbrei, W., 2008. \emph{Multivariable Model - Building:
#' A Pragmatic Approach to Regression Anaylsis based on Fractional Polynomials
#' for Modelling Continuous Variables. John Wiley & Sons.}\cr
#'
#' Royston, P. and Sauerbrei, W., 2016. \emph{mfpa: Extension of mfp using the
#' ACD covariate transformation for enhanced parametric multivariable modeling.
#' The Stata Journal, 16(1), pp.72-87.}
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
                              retain_linear_fp1 = FALSE,
                              has_offset,
                              n_obs,
                              verbose,
                              term_to_columns,
                              fitter = "base") {

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
      selector_args$retain_linear_fp1 <- retain_linear_fp1
      fit1 <- do.call(select_force_max_fp, selector_args)
    }

    if (verbose) {
      print_mfp_step(xi = xi, criterion = criterion, fit = fit1)
    }

    power_best <- normalize_selected_step_powers(fit1, xi, acdx)
    if (saz_active) {
      spike_decision[[xi]] <- saz_decision_codes[["cont_binary"]]
    }

    return(list(
      power_best = power_best,
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

  fit1 <- do.call(select_fct, selector_args)

  if (verbose) {
    print_mfp_step(xi = xi, criterion = criterion, fit = fit1)
  }

  power_best <- normalize_selected_step_powers(fit1, xi, acdx)

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
    spike_decision = stage2$spike_decision,
    current_adj_params = fit1$current_adj_params,
    transform_cache = fit1$transform_cache
  ))

}

# Normalize and name the selected power vector returned by a step selector.
# Keeping this in one helper is important for joint IC because binary-only and
# null models both carry NA powers but have different spike_decision values.
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


#' Function to find the best FP functions of given degree for a single variable
#'
#' Handles the FP1 and the higher order FP cases. For parameter definitions, see
#' \code{find_best_fp_step()}.
#'
#' @details
#' The "best" model is determined by the highest likelihood (or smallest
#' deviance by our definition as minus twice the log-likelihood). This is also
#' the case for the use of information criteria, as all models investigated in
#' this function have the same df, so the penalization term is equal for all
#' models and only their likelihoods differ.
#'
#' Note that the estimation of each fp power adds a degree of freedom. Thus,
#' all fp1s have 2 df, all fp2s have 4 df and so on.
#'
#' In the case that `degree = 1`, the linear model (fp power of 1) is NOT
#' returned, as it is not considered to be a fractional polynomial in this
#' algorithm.
#' A linear model has only one df, whereas the same function regarded as fp
#' would have 2 fp.
#'
#' @section ACD transformation:
#' This function also handles the case of ACD transformations if `acdx` is set
#' to `TRUE` for `xi`. In this case, if `degree = 1`, then 7 models are
#' assessed (like for the non-acd case it excludes the linear case),
#' and if `degree = 2`, then 64 models are assessed (unlike the 36 models
#' for non-acd transformation). Other settings for `degree` are currently not
#' supported when used with ACD transformations.
#'
#' @return
#' A list with several components:
#'
#' * `acd`: logical indicating if an ACD transformation was applied for `xi`.
#' * `powers`: fp powers investigated in step.
#' * `power_best`: the best power found. `power_best` will always be a
#' two-column matrix when an ACD transformation is used, otherwise the number
#' of columns will depend on `degree`.
#' * `metrics`: a matrix with performance indices for all models investigated.
#' Same number of rows as, and indexed by, `powers`.
#' * `model_best`: row index of best model in `metrics`.
#' * `zero`: Logical indicating whether a zero transformation was applied to \code{xi}.
#'   In this case, exact-zero values of \code{xi} remain zero before transformation,
#'   and only positive values were transformed.
#' * `catzero`: Logical indicating whether a combination of a zero transformation
#'   and a binary indicator variable was applied to \code{xi}. This means that
#'   exact-zero values of \code{xi} remain zero, only positive values are
#'   transformed, and an additional binary variable was created to indicate
#'   whether \code{xi} was exactly zero or positive.
#' @inheritParams find_best_fp_step
#' @param degree degrees of freedom for fp transformation of `xi`.
#' @param ... parameters passed to `fit_model()`.
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
                               retain_linear_fp1 = FALSE,
                               ...
) {
  # The conceptual-term lookup is normalized once in fit_mfp() and passed
  # explicitly so every FP candidate uses the same focal and adjustment
  # raw-column blocks.
  # n_obs (number of observations, or events for Cox models) is computed once
  # in fit_mfp() and passed down as a parameter, rather than recomputed here.

  # Step 1: Define the degree-1 candidate set for xi -------------------------
  if (degree == 1L &&
      (!isTRUE(retain_linear_fp1) || isTRUE(acdx[xi]))) {
    # Ordinary MFP closed testing fits the linear model separately, so power 1
    # is excluded from its FP1 candidate search. MFPI differs deliberately:
    # when the user prespecifies FP1, the conventional FP1 class is searched
    # and p = 1 is retained if supplied, so it may legitimately win the FP1
    # power search. This switch never inserts p = 1 into a restricted user
    # candidate set. ACD keeps its historical behaviour because its degree-1
    # p = 1 candidate duplicates an ACD-linear model handled separately by the
    # ACD selection procedure.
    #
    # Do not remove power 1 from the full powers list: transform_data_step()
    # also receives powers and may pass powers[[v]] to adjustment-variable
    # transformations. Adjustment variables must keep their original allowed
    # power sets.
    powers[[xi]] <- setdiff(powers[[xi]], 1)
  }

  # Step 2: Generate candidate FP/ACD transformations for xi and the current
  # adjustment set (or reuse precomputed_adj if the caller already built it).
  # Both ordinary FP and ACD searches now use compact shared bases. Only the
  # winning candidate is materialized after the search for downstream SAZ/cache
  # consumers that genuinely need a complete focal matrix.
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

  # GLMs use an explicit intercept in the design matrix; Cox models do not.
  # Keep this in one flag so matrix construction and model fitting stay aligned.
  x_has_intercept <- !identical(family_string, "cox")


  # Step 3: Build a reusable design-matrix template ---------------------------
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

  for (i in seq_len(n_candidates)) {
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

      fit <- fit_model(
        x               = design_mat,
        y               = y,
        family          = family,
        family_string   = family_string,
        has_offset      = has_offset,
        x_has_intercept = x_has_intercept,
        ...
      )

    } else {
      # Defensive materialized-candidate fallback.
      data_xi <- data_fp[[i]]

      if (!is.null(design_mat) && NCOL(data_xi) == n_xi_cols) {
        design_mat[, xi_cols] <- data_xi
        fit <- fit_model(
          x               = design_mat,
          y               = y,
          family          = family,
          family_string   = family_string,
          has_offset      = has_offset,
          x_has_intercept = x_has_intercept,
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
        fit <- fit_model(
          x               = x_fit,
          y               = y,
          family          = family,
          family_string   = family_string,
          has_offset      = has_offset,
          x_has_intercept = x_has_intercept,
          ...
        )
      }
    }

    metrics[i, ] <- calculate_model_metrics(
      fit, # fit includes catzero df when used, so this part does not change
      n_obs,
      degree
    )
  }

  # Step 5: Pick the candidate with the highest log-likelihood ---------------
  # (equivalent to lowest deviance/AIC/BIC here, since all candidates in this
  # call share the same degrees of freedom; see @details above).
  model_best <- as.numeric(which.max(metrics[, "logl"]))

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

#' Function to fit a null model excluding variable of interest
#'
#' "Null" model here refers to a model which does not include the variable
#' of interest `xi`.
#' For parameter definitions, see \code{find_best_fp_step()}. All parameters
#' captured by `...` are passed on to \code{fit_model()}.
#'
#' @return
#' A list with three entries:
#'
#' * `powers`: FP power(s) of `xi` in fitted model - in this case `NA`.
#' * `metrics`: A matrix with performance indices for fitted model.
#' * `current_adj_params`: Adjustment-variable transformations for `xi`,
#'   cached for reuse in later steps (see \code{prev_adj_params} in
#'   \code{find_best_fp_step()}).
#'
#' @inheritParams find_best_fp_step
#' @param ... Parameters passed to `fit_model()`.
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
  # In addition, adjustment model can be NULL, so we have intercept only
  model_null <- fit_model(x = x_tran,
                          y = y,
                          family = family,
                          family_string   = family_string,
                          has_offset = has_offset,
                          ...
  )

  # Step 3: Score the null model and return, with xi's power fixed at NA ----
  list(
    powers = NA,
    metrics = rbind(null = calculate_model_metrics(model_null, n_obs)),
    current_adj_params = current_params
  )
}

#' Function to fit linear model for variable of interest
#'
#' "Linear" model here refers to a model that includes the variable
#' of interest \code{xi} with an FP (fractional polynomial) power of 1.
#' Note that \code{xi} may be ACD-transformed if indicated by \code{acdx[xi]}.
#' If the variable was passed through the \code{catzero} argument in \code{mfp2()},
#' both the continuous variable and its corresponding binary indicator
#' will be included in the model as linear terms.
#' For parameter definitions, see \code{find_best_fp_step}.
#' All parameters captured by \code{...} are passed to \code{fit_model}.
#'
#' @return
#' A list with three entries:
#'
#' * `powers`: FP power(s) of `xi` (or its ACD transformation) in fitted model.
#' * `metrics`: A matrix with performance indices for fitted model.
#' * `current_adj_params`: Adjustment-variable transformations for `xi`,
#'   cached for reuse in later steps (see \code{prev_adj_params} in
#'   \code{find_best_fp_step()}).
#' @inheritParams find_best_fp_step
#' @param ... Parameters passed to `fit_model()`.
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
                            ...,
                            precomputed_adj = NULL,
                            term_to_columns) {
  # n_obs (number of observations, or events for Cox models) is computed once
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

  # Step 2: Assemble the design matrix and fit the linear model -------------
  # Fit a model based on the assumption that xi is linear.
  # If catzero or spike-at-zero is active for xi, the corresponding binary
  # component is already included in data_xi.
  data_xi <- x_transformed$data_fp[[1L]]
  data_adj <- x_transformed[["data_adj", exact = TRUE]]

  has_adj <- !is.null(data_adj) && NCOL(data_adj) > 0L
  use_glm_intercept_template <- !identical(family_string, "cox")

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

  if (acdx[xi])
    rownames(metrics) <- "linear(., A(x))"

  list(
    powers = x_transformed$powers_fp,
    metrics = metrics,
    current_adj_params = x_transformed$current_params
  )
}

#' Helper function to select between null and linear term for a single variable
#'
#' To be used in \code{find_best_fp_step()}. Only used if `df = 1` for a variable.
#' Handles all criteria for selection.
#' For parameter explanations, see \code{find_best_fp_step()}. All parameters
#' captured by `...` are passed to \code{fit_model()}.
#'
#' @details
#' This function assesses a single variable of interest `xi` regarding its
#' functional form in the current working model as indicated by
#' `powers_current`, with the choice between excluding `xi` ("null model") and
#' including a linear term ("linear fp") for `xi`.
#'
#' Note that this function handles an ACD transformation for `xi` as well.
#'
#' When a variable is forced into the model by including it in `keep`, then
#' this function will not exclude it from the model (by setting its power to
#' `NA`), but will only choose the linear model.
#'
#' @return
#' A list with several components:
#'
#' * `keep`: logical indicating if `xi` is forced into model.
#' * `acd`: logical indicating if an ACD transformation was applied for `xi`.
#' * `powers`: fp powers investigated in step, indexing `metrics`.
#' * `power_best`: a numeric vector with the best power found. The returned
#' best power may be `NA`, indicating the variable has been removed from the
#' model.
#' * `metrics`: a matrix with performance indices for all models investigated.
#' Same number of rows as, and indexed by, `powers`.
#' * `model_best`: row index of best model in `metrics`.
#' * `pvalue`: p-value for comparison of linear and null model.
#' * `statistic`: test statistic used, depends on `ftest`.
#' * `zero`: Logical indicating whether a zero transformation was applied to \code{xi}.
#'   In this case, exact-zero values of \code{xi} remain zero before transformation,
#'   and only positive values were transformed.
#' * `catzero`: Logical indicating whether a combination of a zero transformation
#'   and a binary indicator variable was applied to \code{xi}. This means that
#'   exact-zero values of \code{xi} remain zero, only positive values are
#'   transformed, and an additional binary variable was created to indicate
#'   whether \code{xi} was exactly zero or positive.
#' * `spike`: Logical; whether `xi` is (still) treated as a spike-at-zero
#'   variable, carried through from the `spike` argument.
#' * `current_adj_params`: Adjustment-variable transformations for `xi` from
#'   whichever of the null/linear models was selected, cached for reuse in
#'   later steps (see \code{prev_adj_params} in \code{find_best_fp_step()}).
#' @param degree not used.
#' @param force_max_fp not used
#' @param ... passed to fitting functions.
#' @inheritParams find_best_fp_step
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
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, has_offset = has_offset, spike = spike,
    acd_parameter = acd_parameter, prev_adj_params = prev_adj_params,
    precomputed_adj = precomputed_adj, n_obs = n_obs, term_to_columns = term_to_columns, ...
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

#' Function selection procedure based on closed testing procedure
#'
#' Used in \code{find_best_fp_step()} when `criterion = "pvalue"`.
#' For parameter explanations, see \code{find_best_fp_step()}. All parameters
#' captured by `...` are passed to \code{fit_model()}.
#'
#' @details
#' In case `criterion = "pvalue"` the function selection procedure as outlined
#' in Chapters 4 and 6 of Royston and Sauerbrei (2008) is used.
#'
#' * \emph{Step 1}: Test the best FP\emph{m} function against a null model at the
#' significance level specified by \code{select}, using 2\emph{m} degrees of freedom.
#' If the test is not significant, the variable is excluded. Otherwise, proceed
#' to Step 2.
#' * \emph{Step 2}: Test the best FP\emph{m} function against a linear model at
#' the significance level specified by \code{alpha}, using 2\emph{m}-1 degrees
#' of freedom. If the test is not significant, select the linear model.
#' Otherwise, proceed to Step 3.
#' * \emph{Step 3}: Test the best FP\emph{m} function against the best FP1 model
#' at the significance level specified by \code{alpha}, using 2\emph{m}-2 degrees
#' of freedom. If the test is not significant, retain the best FP1 model. Otherwise,
#' repeat this step by comparing FP\emph{m} to all remaining lower-order FP models,
#' down to \eqn{\mathrm{FP}_{m-1}}, which is tested with 2 degrees of freedom.
#' If the final test is not significant, retain the best \eqn{\mathrm{FP}_{m-1}}
#' model; otherwise, retain the best FP\emph{m} model.
#'
#' MFP uses a strict boundary when deciding that a simpler model is adequate:
#' simplification/removal occurs only when the comparison p-value is strictly
#' greater than the relevant threshold. Consequently, `select = 1` forces
#' inclusion and `alpha = 1` prevents functional-form simplification, including
#' the exact boundary case where a comparison yields `p = 1`.
#'
#' Note that the "best" FP\emph{x} model used in each step refers to the model
#' that applies an FP\emph{x} transformation to the variable of interest and
#' achieves the highest likelihood among all such models, given the current
#' power transformations for all other variables. This procedure is described
#' in Section 4.8 of Royston and Sauerbrei (2008). The best FP\emph{x} models
#' are computed by \code{find_best_fpm_step}.
#'
#' When a variable is forced into the model by including it in the \code{keep}
#' argument of \code{mfp2()}, this function will not exclude it (i.e., will not
#' set its power to \code{NA}), but will instead select its functional form.
#'
#' @return
#' A list with several components:
#'
#' * `keep`: logical indicating if `xi` is forced into model.
#' * `acd`: logical indicating if an ACD transformation was applied for `xi`,
#' i.e. `FALSE` in this case.
#' * `powers`: (best) fp powers investigated in step, indexing `metrics`.
#' Always starts with highest power, then null, then linear, then FP in
#' increasing degree (e.g. FP2, null, linear, FP1).
#' * `power_best`: a numeric vector with the best power found. The returned
#' best power may be `NA`, indicating the variable has been removed from the
#' model.
#' * `metrics`: a matrix with performance indices for all models investigated.
#' Same number of rows as, and indexed by, `powers`.
#' * `model_best`: row index of best model in `metrics`.
#' * `pvalue`: p-value for comparison of linear and null model.
#' * `statistic`: test statistic used, depends on `ftest`.
#' * `spike`: Logical; whether `xi` is (still) treated as a spike-at-zero
#'   variable, carried through from the `spike` argument.
#' * `current_adj_params`: Adjustment-variable transformations for `xi` from
#'   the selected model, cached for reuse in later steps (see
#'   \code{prev_adj_params} in \code{find_best_fp_step()}).
#' @references
#' Royston, P. and Sauerbrei, W., 2008. \emph{Multivariable Model - Building:
#' A Pragmatic Approach to Regression Anaylsis based on Fractional Polynomials
#' for Modelling Continuous Variables. John Wiley & Sons.}
#'
#' @seealso
#' \code{select_ra2_acd()}
#'
#' @param degree integer > 0 giving the degree for the FP transformation.
#' @param ... passed to fitting functions \code{fit_model()}.
#' @inheritParams find_best_fp_step
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
    calculate_test <- function(metrics, n_obs) {
      calculate_f_test(
        deviances = metrics[, "deviance_gaussian", drop = TRUE],
        dfs_resid = metrics[, "df_resid", drop = TRUE],
        n_obs = n_obs
      )
    }
  } else {
    calculate_test <- function(metrics, n_obs) {
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
    precomputed_adj = precomputed_adj, term_to_columns = term_to_columns, ...
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
  stats <- calculate_test(res$metrics[c("null", fpmax), ], n_obs)
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

  stats <- calculate_test(res$metrics[c(lin_names, fpmax), ], n_obs)
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
        precomputed_adj = precomputed_adj,  term_to_columns = term_to_columns, ...
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
      stats <- calculate_test(res$metrics[c(fpm, fpmax), ], n_obs)
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
        res$power_best = fit_fpm$powers[fit_fpm$model_best, , drop = FALSE]
        res$model_best = nrow(res$metrics)
        res$current_adj_params <- fit_fpm$current_adj_params
        return(res)
      }
    }

  }

  # Step 6: Every lower-degree FP tested was significantly worse than
  # FPm (or degree == 1, so there were no lower degrees to test at
  # all): retain the highest-degree FPm itself (row 1, inserted first).
  res$power_best <- fit_fpmax$powers[fit_fpmax$model_best, , drop = FALSE]
  res$model_best <- 1
  res$current_adj_params <- fit_fpmax$current_adj_params

  res
}

#' Function selection procedure for ACD based on closed testing procedure
#'
#' Used in \code{find_best_fp_step()} when `criterion = "pvalue"` and an
#' ACD transformation is requested for `xi`.
#' For parameter explanations, see \code{find_best_fp_step()}. All parameters
#' captured by `...` are passed on to \code{fit_model()}.
#'
#' @details
#' This function extends the algorithm used in \code{select_ra2()} to allow the
#' usage of ACD transformations. The implementation follows the description
#' in Royston and Sauerbrei (2016). The procedure is outlined in detail in
#' the corresponding section in the documentation of \code{mfp2()}.
#'
#' When a variable is forced into the model by including it in `keep`, then
#' this function will not exclude it from the model (by setting its power to
#' `NA`), but will only choose its functional form.
#'
#' @return
#' A list with several components:
#'
#' * `keep`: logical indicating if `xi` is forced into model.
#' * `acd`: logical indicating if an ACD transformation was applied for `xi`,
#' i.e. `FALSE` in this case.
#' * `powers`: (best) fp powers investigated in step, indexing `metrics`.
#' Ordering: FP1(x, A(x)), null, linear, FP1(x, .), linear(., A(x)),
#' FP1(., A(x)).
#' * `power_best`: a numeric vector with the best power found. The returned
#' best power may be `NA`, indicating the variable has been removed from the
#' model.
#' * `metrics`: a matrix with performance indices for all models investigated.
#' Same number of rows as, and indexed by, `powers`.
#' * `model_best`: row index of best model in `metrics`.
#' * `pvalue`: p-value for comparison of linear and null model.
#' * `statistic`: test statistic used, depends on `ftest`.
#' * `spike`: Logical; whether `xi` is (still) treated as a spike-at-zero
#'   variable, carried through from the `spike` argument.
#' * `current_adj_params`: Adjustment-variable transformations for `xi` from
#'   the selected model, cached for reuse in later steps (see
#'   \code{prev_adj_params} in \code{find_best_fp_step()}).
#'
#' @references
#' Royston, P. and Sauerbrei, W., 2016. \emph{mfpa: Extension of mfp using the
#' ACD covariate transformation for enhanced parametric multivariable modeling.
#' The Stata Journal, 16(1), pp.72-87.}
#'
#' @seealso
#' \code{select_ra2()}
#'
#' @param degree integer > 0 giving the degree for the FP transformation.
#' @param ... passed to fitting functions.
#' @inheritParams find_best_fp_step
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
    calculate_test <- function(metrics, n_obs) {
      calculate_f_test(
        deviances = metrics[, "deviance_gaussian", drop = TRUE],
        dfs_resid = metrics[, "df_resid", drop = TRUE],
        n_obs = n_obs
      )
    }
  } else {
    calculate_test <- function(metrics, n_obs) {
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
    precomputed_adj = precomputed_adj, term_to_columns = term_to_columns, ...
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
  stats <- calculate_test(res$metrics[c("null", fpmax), ], n_obs)
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

  stats <- calculate_test(res$metrics[c(lin_names, fpmax), ], n_obs)
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
    precomputed_adj = precomputed_adj,  term_to_columns = term_to_columns, ...
  )

  old_names = rownames(res$metrics)
  res$metrics <- rbind(
    res$metrics,
    fit$metrics[fit$model_best, ]
  )

  FP_names <- paste0("FP1(x, .)", binary_suffix)

  rownames(res$metrics) <- c(old_names, FP_names)
  res$powers <- rbind(res$powers, c(fit$power_best, NA))

  stats <- calculate_test(res$metrics[c(FP_names, fpmax), ], n_obs)
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
    precomputed_adj = precomputed_adj, term_to_columns = term_to_columns, ...
  )

  old_names = rownames(res$metrics)
  res$metrics <- rbind(
    res$metrics,
    fit_fp1a$metrics[fit_fp1a$model_best, ]
  )

  FP1_names <- paste0("FP1(., A(x))", binary_suffix)
  rownames(res$metrics) <- c(old_names, FP1_names)
  res$powers <- rbind(res$powers, fit_fp1a$power_best)

  stats <- calculate_test(res$metrics[c(FP1_names, fpmax), ], n_obs)
  old_names <- names(res$statistic)
  res$statistic <- c(res$statistic, stats$statistic)
  names(res$statistic) <- c(old_names, sprintf("%s vs %s", fpmax, FP1_names))
  res$pvalue <- c(res$pvalue, stats$pvalue)
  names(res$pvalue) <- names(res$statistic)

  if (stats$pvalue < alpha) {
    # Test 4 significant: M3 is significantly worse than M1, so M1 cannot be
    # simplified any further. Stop with the full FP1(x, A(x)) model.
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

  stats <- calculate_test(res$metrics[c(linx_names, FP1_names), ], n_obs)
  old_names <- names(res$statistic)
  res$statistic <- c(res$statistic, stats$statistic)
  names(res$statistic) <- c(old_names, sprintf("%s vs %s", FP1_names, linx_names))
  res$pvalue <- c(res$pvalue, stats$pvalue)
  names(res$pvalue) <- names(res$statistic)

  if (stats$pvalue < alpha) {
    # Test 5 significant: M5 (linear in A(x)) is significantly worse than M3,
    # so M3 cannot be simplified further. Stop with FP1(., A(x)).
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
#' Internal selector used when `force_max_fp` is active for a non-linear term
#' under p-value, AIC, or BIC selection. Because the public option explicitly
#' bypasses variable selection and functional-form simplification, the maximum
#' functional form is predetermined. Fitting the null, linear, lower-degree FP,
#' or reduced ACD alternatives therefore cannot change the result. This helper
#' builds the adjustment set once and fits only the required maximum form.
#'
#' For an ordinary FP term, the required model is the best FP of the requested
#' `degree`.  For an ACD term, the required model is the full
#' `FP1(x, A(x))` form, represented internally by the degree-2 ACD search.
#' The best power combination within that fixed form is still chosen by
#' `find_best_fpm_step()`.
#'
#' The returned object deliberately has the same shape as the other selectors,
#' but `metrics` and `powers` contain only the single model that was actually
#' fitted. This is sufficient for verbose printing. When the forced term is an
#' eligible spike-at-zero variable, \code{find_best_fp_step()} retains this full
#' continuous-plus-binary representation directly and skips SAZ stage 2 because
#' no reduced component representation is permitted to replace it.
#'
#' @inheritParams find_best_fp_step
#' @param degree integer > 0 giving the requested ordinary FP degree.
#' @param ... passed to fitting functions.
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
                                retain_linear_fp1 = FALSE,
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
  # Ordinary FP uses the requested maximum degree directly. For MFPI FP1,
  # retain_linear_fp1 is passed through unchanged so p = 1 remains available
  # when it is part of the supplied candidate set.
  forced_degree <- if (is_acd) 2 else degree

  fit_max <- find_best_fpm_step(
    x = x, xi = xi, degree = forced_degree, y = y,
    powers_current = powers_current, powers = powers, acdx = acdx,
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter, spike = spike,
    prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
    precomputed_adj = precomputed_adj, term_to_columns = term_to_columns,
    retain_linear_fp1 = retain_linear_fp1, ...
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
    spike = spike[xi],
    current_adj_params = fit_max$current_adj_params,
    transform_cache = precomputed_adj$transform_cache
  )
}

#' Function selection procedure based on information criteria
#'
#' Used in \code{find_best_fp_step()} when `criterion = "aic"` or `"bic"`.
#' For parameter explanations, see \code{find_best_fp_step()}. All parameters
#' captured by `...` are passed on to \code{fit_model()}.
#'
#' @details
#' In case an information criterion is used to select the best model the
#' selection procedure simply fits all relevant models and selects the best
#' one according to the given criterion.
#'
#' "Relevant" models for a given degree are the null model excluding the
#' variable of interest, the linear model and all best FP models up to the
#' specified degree.
#'
#' In case an ACD transformation is requested, then the models assessed
#' are the null model, the linear model in x and A(x), the best FP1 models in
#' x and A(x), and the best FP1(x, A(x)) model.
#'
#' Note that the "best" FPx model used in this function are given by the models
#' using a FPx transformation for the variable of interest and having the
#' highest likelihood of all such models given the current powers for all other
#' variables, as outlined in Section 4.8 of Royston and Sauerbrei (2008).
#' These best FPx models are computed in \code{find_best_fpm_step()}.
#' Keep in mind that for a fixed number of degrees of freedom (i.e. fixed m),
#' the model with the highest likelihood is the same as the model with the best
#' information criterion of any kind since all the models share the same
#' penalty term.
#'
#' When a variable is forced into the model by including it in `keep`, then
#' this function will not exclude it from the model (by setting its power to
#' `NA`), but will only choose its functional form.
#'
#' @return
#' A list with several components:
#'
#' * `keep`: logical indicating if `xi` is forced into model.
#' * `acd`: logical indicating if an ACD transformation was applied for `xi`,
#' i.e. `FALSE` in this case.
#' * `powers`: (best) fp powers investigated in step, indexing `metrics`.
#' Ordered by increasing complexity, i.e. null, linear, FP1, FP2 and so on.
#' For ACD transformation, it is null, linear, linear(., A(x)), FP1(x, .),
#' FP1(., A(x)) and FP1(x, A(x)).
#' * `power_best`: a numeric vector with the best power found. The returned
#' best power may be `NA`, indicating the variable has been removed from the
#' model.
#' * `metrics`: a matrix with performance indices for all best models
#' investigated. Same number of rows as, and indexed by, `powers`.
#' * `model_best`: row index of best model in `metrics`.
#' * `pvalue`: p-value for comparison of linear and null model, `NA` in this
#' case..
#' * `statistic`: test statistic used, depends on `ftest`, `NA` in this
#' case.
#' * `spike`: Logical; whether `xi` is (still) treated as a spike-at-zero
#'   variable, carried through from the `spike` argument.
#' * `current_adj_params`: Adjustment-variable transformations for `xi` from
#'   the selected model, cached for reuse in later steps (see
#'   \code{prev_adj_params} in \code{find_best_fp_step()}).
#'
#' @seealso
#' \code{select_ra2()}
#'
#' @param degree integer > 0 giving the degree for the FP transformation.
#' @param ... passed to fitting functions.
#' @inheritParams find_best_fp_step
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
    precomputed_adj = precomputed_adj, term_to_columns = term_to_columns, ...
  )
  # Linear model
  fit_lin <- fit_linear_step(
    x = x, xi = xi, y = y,
    powers_current = powers_current, powers = powers, acdx = acdx,
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter, spike = spike,
    prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
    precomputed_adj = precomputed_adj, term_to_columns = term_to_columns, ...
  )

  res$current_adj_params <- fit_null$current_adj_params


  # Build the focal numerical basis once at the maximum requested degree.
  # select_ic() necessarily evaluates every degree, so this basis will be reused
  # and there is no RA2-style early-stopping tradeoff. Build from the complete
  # user-supplied power set (including power 1); find_best_fpm_step() continues
  # to remove power 1 only from its degree-1 candidate list, exactly as before.
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
      term_to_columns = term_to_columns, ...
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

#' @describeIn select_ic Function to select ACD based transformation.
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
    precomputed_adj = precomputed_adj, term_to_columns = term_to_columns, ...
  )

  # linear(x, .)
  fit_lin <- fit_linear_step(
    x = x, xi = xi, y = y,
    powers_current = powers_current, powers = powers, acdx = acdx_reset_xi,
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter, spike = spike,
    prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
    precomputed_adj = precomputed_adj, term_to_columns = term_to_columns, ...
  )

  # linear(., A(x))
  fit_lina <- fit_linear_step(
    x = x, xi = xi, y = y,
    powers_current = powers_current, powers = powers, acdx = acdx,
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter,spike = spike,
    prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
    precomputed_adj = precomputed_adj, term_to_columns = term_to_columns, ...
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
        term_to_columns = term_to_columns, ...
      ),
      # FP1(., A(x))
      find_best_fpm_step(
        x = x, xi = xi, degree = 1, y = y,
        powers_current = powers_current, powers = powers, acdx = acdx,
        family = family, family_string = family_string, zero = zero, catzero = catzero,
        spike_decision = spike_decision, acd_parameter = acd_parameter,spike = spike,
        prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
        precomputed_adj = precomputed_adj, focal_basis_cache = focal_basis_cache,
        term_to_columns = term_to_columns, ...
      ),
      # FP1(x, A(x))
      find_best_fpm_step(
        x = x, xi = xi, degree = 2, y = y,
        powers_current = powers_current, powers = powers, acdx = acdx,
        family = family, family_string = family_string, zero = zero, catzero = catzero,
        spike_decision = spike_decision,acd_parameter = acd_parameter,spike = spike,
        prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
        precomputed_adj = precomputed_adj, focal_basis_cache = focal_basis_cache,
        term_to_columns = term_to_columns, ...
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

#' Build Adjustment-Step Matrices Using the C++ Hot Loop
#'
#' Prepares aligned inputs for `build_adjustment_step_loop_cpp()` and calls the
#' C++ implementation of the per-variable adjustment loop used inside
#' `build_adjustment_step()`.
#'
#' This helper is intentionally scoped to the MFP adjustment-step path. It does
#' not estimate ACD parameters. ACD adjustment variables are expected to have
#' stored `acd_parameter` values already; the C++ code applies those stored
#' parameters when recomputing ACD adjustment columns.
#'
#' @param x Numeric matrix or numeric data frame containing the original
#'   predictors.
#' @param vars_adj Character vector of adjustment-variable names.
#' @param powers_adj Named list of current selected powers for `vars_adj`.
#' @param acdx_adj Named logical vector/list indicating ACD variables.
#' @param zero_adj Named logical vector/list indicating zero handling.
#' @param catzero Named list of structural-zero indicator matrices. Elements
#'   are either `NULL` or `nrow(x) x 1` matrices.
#' @param spike_adj Named logical vector/list indicating spike-at-zero variables.
#' @param spike_decision_int_adj Named integer vector of spike decisions.
#' @param acd_parameter_adj Named list of stored ACD parameters for `vars_adj`.
#' @param eliminated Named logical vector indicating variables with all-NA powers.
#' @param spike_binary_only_flags Named logical vector indicating
#'   `spike_decision == 3`.
#' @param current_power_keys_adj Named list of normalized current power keys.
#' @param prev_power_keys_adj Named list of normalized previous power keys, or
#'   `NULL` if there is no previous cache.
#' @param prev_xi Previous per-variable cache object containing aligned
#'   transformed blocks and spike decisions, or `NULL`. Historically this was
#'   the cache for the current focal variable; it may now be a shared cache.
#' @param has_prev Logical scalar; whether `prev_xi` is available.
#'
#' @return A list with:
#' \describe{
#'   \item{data_adj_list}{Named list of per-variable adjustment matrices.}
#'   \item{data_adj}{The final column-bound adjustment matrix.}
#' }
#'
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
#' Internal helper used by transform_data_step() and selection functions to
#' avoid recomputing the adjustment matrix repeatedly for the same focal
#' variable, powers_current, and spike_decision state.
#' @inheritParams transform_data_step
#' @return
#' A list with five entries:
#' * `powers_adj`: named list of FP powers used for each adjustment variable
#'   (i.e. \code{powers_current} restricted to variables other than `xi`), or
#'   `NULL` if there are no adjustment variables.
#' * `spike_decision_adj`: named vector of spike-at-zero decisions for each
#'   adjustment variable, or `NULL` if there are no adjustment variables.
#' * `data_adj_list`: named list of per-adjustment-variable transformed
#'   matrices (used as a cache for the next call), or an empty list if there
#'   are no adjustment variables.
#' * `data_adj`: the column-bound adjustment design matrix, `NULL` if there are
#'   no adjustment variables, or an `n x 0` matrix if adjustment variables
#'   exist but all currently contribute zero columns (e.g. all eliminated).
#' * `transform_cache`: updated run-scoped per-variable transformation cache.
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
#' Function to extract and transform adjustment variables
#'
#' This function prepares transformed data for a focal predictor `xi` and its
#' adjustment variables. Adjustment variables are transformed using either
#' fractional polynomials or acd transformations, depending on their assigned
#' powers and parameters. Spike-at-zero effects can be incorporated using
#' binary indicators. Previously computed adjustment variables can be reused
#' if their parameters have not changed.
#' @param x a matrix of predictors that includes the variable of interest `xi`.
#' It is assumed that continuous variables have already been shifted and scaled.
#' @param xi Name of the conceptual term being assessed. It may be a
#'   singleton FP-eligible predictor or a fixed multi-column linear block. All
#'   other conceptual terms are adjustment terms.
#' @param powers_current a named list of FP powers of all variables of interest,
#' including `xi`. Note that these powers are updated during backfitting or MFP
#' cycles.
#' @param df a numeric vector of degrees of freedom for `xi`.
#' @param powers a set of allowed FP powers.
#' @param acdx a logical vector indicating the use of acd transformation.
#' @param zero named logical vector of length ncol(x)
#' @param spike A logical vector indicating which columns of \code{x} contain
#' a spike at zero. The length and order of \code{spike} must match those of
#' the columns in \code{x}.
#' @param catzero A named list of binary indicator variables of length \code{ncol(x)}
#' for exact-zero values, created when specific variables are passed to the
#' \code{catzero} argument of \code{fit_mfp}. If an element of the list is
#' \code{NULL}, it indicates that the corresponding variable was not specified by
#' the user in the \code{catzero} argument of \code{fit_mfp}. Here, \code{catzero}
#' is a list of binary variables, not a named logical vector as in \code{fit_mfp}.
#' @param acd_parameter Named list of ACD parameters produced by \code{fit_acd()},
#' with length equal to \code{ncol(x)}.
#' @param spike_decision a named numeric vector with the same names as the
#' `x` matrix. Each element controls how the corresponding adjustment
#' variable contributes to the adjustment matrix:
#'
#' * `1`: combine both the transformed adjustment variable and its binary
#'   indicator (`catzero`).
#' * `2`: include only the transformed adjustment variable.
#' * `3`: include only the binary indicator (`catzero`).
#' This applies only to adjustment variables, not the focal predictor `xi`.
#' @param prev_adj_params Named list containing results from a previous call to
#' this function. Used to avoid recomputing adjustment variables if both powers
#' and `spike_decision` remain unchanged for a certain variable.
#' @param precomputed_adj Optional internal adjustment object returned by
#' \code{build_adjustment_step()}. When supplied, adjustment-variable
#' transformations are reused and only current-variable candidate
#' transformations are generated.
#' @param compact_fp Logical; use the shared-basis representation for ordinary
#'   FP candidates instead of materializing every candidate matrix.
#' @param compact_acd Logical; use the shared-basis representation for ACD
#'   candidates instead of materializing every candidate matrix.
#' @details
#' This function extracts the adjustment variables and applies the corresponding
#' FP or ACD transformations based on `powers_current`. When evaluating the variable
#' of interest `xi`, it is necessary to account for other variables in the model,
#' which may be transformed or untransformed depending on their individual powers.
#' Some powers may be `NA`, indicating that the corresponding variable has been
#' excluded from the adjustment set.
#'
#' To improve efficiency, the function avoids recomputing adjustment variables if
#' their parameters (`powers` and `spike_decision`) have not changed from a
#' previous step, as provided in `prev_adj_params`.
#'
#' The role of `spike_decision` is to determine how each adjustment variable is
#' represented in the presence of potential spike-at-zero effects. For every
#' adjustment variable, the function can include the transformed variable, the
#' binary indicator for exact-zero values (`catzero`), or both. This makes it
#' possible to model spike-at-zero behavior directly in the adjustment matrix,
#' while preserving flexibility in how variables are included.
#'
#' The function also returns candidate information for predictor `xi`. In the
#' historical materialized mode, `data_fp` contains one matrix per candidate.
#' In compact mode, candidate values are represented by a shared basis plus a
#' small integer map. For default ACD degree 2 this means 16 unique basis
#' columns plus a 64 x 2 map instead of 64 separate n x 2 matrices.
#'
#' When `df = 1`, this function returns data unchanged, i.e. a "linear"
#' transformation with power equal to 1. In case `acdx[xi] = TRUE`, the
#' acd transformation is applied.
#'
#' @return
#' A list containing `power_best` (numeric vector of best FP powers for `xi`,
#' possibly including `NA`), `spike_decision` (updated spike-at-zero strategy),
#' and `current_adj_params` (adjustment variable transformations used in this step).
#' \item{powers_fp}{Numeric vector of FP powers used for `xi`.}
#' \item{data_fp}{List of transformed data for `xi`, or `NULL` in compact mode.}
#' \item{fp_basis}{Compact ordinary-FP basis, or `NULL`.}
#' \item{acd_basis}{Compact ACD basis, or `NULL`.}
#' \item{powers_adj}{Named list of FP powers used for adjustment variables.}
#' \item{data_adj}{Matrix of transformed adjustment variables, or `NULL` if none.}
#' \item{current_params}{Named list containing adjustment variable parameters
#'  for reuse in later steps, including \code{powers_adj},
#'   \code{spike_decision_adj}, \code{data_adj_list}, and \code{data_adj}.}
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
            # call still generates its own degree-specific candidates (including
            # FP1's existing removal of power 1) and maps them by explicit
            # power/repetition metadata rather than fixed column positions.
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
