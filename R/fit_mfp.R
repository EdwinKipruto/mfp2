#' Validate MFP Candidate Powers Against Closed-Test Requirements
#'
#' Checks that each variable with \code{df > 1} has at least one candidate
#' fractional-polynomial power other than \code{1}. In the MFP selection engine,
#' power \code{1} is fitted separately as the ordinary linear model and is
#' therefore excluded from degree-1 FP candidate fitting. The closed-test
#' procedure requires an available FP1 candidate whenever \code{df > 1}, so a
#' candidate-power set containing only \code{1} is invalid in that case.
#'
#' Repeated selected powers remain valid for FP2 and higher-degree models. For
#' example, \code{c(1, 1)} is a valid selected FP2 power vector, corresponding
#' to the basis \eqn{x} and \eqn{x \log(x)}. This helper does not reject repeated
#' selected powers; it only rejects candidate-power sets that collapse to the
#' single value \code{1} when the algorithm needs a non-linear FP candidate.
#'
#' @param powers Named list of numeric candidate-power vectors, one element per
#'   variable. Each element is treated as a candidate set: values are coerced to
#'   numeric, deduplicated, sorted, and checked after removing missing and
#'   non-finite values.
#' @param df Named numeric or integer vector of degrees of freedom, one value
#'   per variable. Values greater than \code{1} indicate that FP model selection
#'   may require degree-1 candidate fitting.
#'
#' @return Invisibly returns \code{TRUE} if all candidate-power sets are
#'   compatible with the closed-test procedure.
#'
#' @details
#' This is an internal selection-engine validation. It should be called after
#' \code{powers} and \code{df} have been normalized to named per-variable
#' objects and before calls to \code{select_ra2()}, \code{select_ic()}, or
#' \code{find_best_fpm_step()}.
#'
#' @keywords internal
#' @noRd
validate_mfp_candidate_powers <- function(powers, df) {
  # Step 1: Flag variables whose candidate-power set collapses to just `1`.
  # A candidate set of exactly {1} after cleaning (dedup, drop NA/non-finite)
  # gives the closed-test procedure no non-linear FP1 candidate to test against
  # the linear model, which is required whenever df > 1.
  only_linear_candidate <- vapply(
    powers,
    function(v) {
      v <- sort(unique(as.numeric(v)))
      v <- v[!is.na(v) & is.finite(v)]
      length(v) == 1L && identical(v, 1)
    },
    logical(1L)
  )
  
  # Step 2: Report an error listing every variable where this is a problem
  # (only-linear candidate set combined with df > 1), rather than failing on
  # just the first offending variable.
  vars_invalid <- names(which(only_linear_candidate & df > 1L))
  
  if (length(vars_invalid) > 0L) {
    stop(
      paste0(
        "The following variable(s) have `df > 1` but their candidate-power ",
        "set contains only power 1: ",
        paste(vars_invalid, collapse = ", "),
        ". Power 1 is fitted separately as the ordinary linear model and is ",
        "excluded from FP1 candidate fitting. The closed-test procedure requires ",
        "at least one available FP1 candidate when `df > 1`. Use `df = 1` for a ",
        "purely linear effect, or include at least one non-1 candidate power."
      ),
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}


#' Function for fitting a model using the MFP, MFPA or spike-at-zero algorithm
#'
#' This internal function implements the Multivariable Fractional Polynomial (MFP),
#' MFP with Approximate Cumulative Distribution (MFPA), and spike-at-zero (SAZ)
#' algorithms. It is not exported and is intended to be called from \code{mfp2()}.
#' While most parameters are documented in \code{mfp2()}, their form may differ here.
#' The function does not perform argument checks and expects that all inputs
#' have been properly prepared by \code{mfp2()}.
#'
#' @param x A numeric design matrix with one row per observation and one or
#'   more raw columns per conceptual term. It excludes the intercept; categorical
#'   terms may already be expanded into dummy columns. Data are assumed to be
#'   shifted and scaled.
#' @param term_to_columns Named list mapping every conceptual term to the raw
#'   columns of \code{x} that represent it. Explicitly mapped entries are fixed
#'   linear design blocks, including a one-column mapping whose raw column name
#'   differs from the conceptual term. Identity singleton entries retain
#'   ordinary MFP behavior.
#' @param y A vector for the response variable or a \code{Surv} object.
#' @param weights A vector of observation weights of length nobs.
#' @param offset A vector of length nobs of offsets.
#' @param cycles An integer representing the maximum number of iteration cycles
#'   during which FP powers for all predictors are updated.
#' @param scale Named numeric vector with one scaling factor per conceptual
#'   term. It is not applied during power selection cycles, but is used to
#'   backscale the corresponding raw column or column block before the final fit.
#' @param shift Named numeric vector with one shift per conceptual term. Shifting
#'   has already been applied upstream; the values are reordered to conform to
#'   \code{xorder} and stored in the returned object.
#' @param df Named numeric vector with one degrees-of-freedom setting per
#'   conceptual term. Explicitly mapped terms must have \code{df = 1}.
#' @param center Named logical vector with one centering setting per conceptual
#'   term.
#' @param family Either a character string specifying the model family
#'   (e.g., \code{"gaussian"}, \code{"binomial"}, \code{"poisson"},
#'   \code{"cox"}) or a function that returns a GLM family object such as
#'   \code{stats::gaussian(link = "identity")}. For Cox models, only the
#'   character string \code{"cox"} is allowed.
#' @param family_string A character string representing the selected family,
#'   e.g., \code{"gaussian"}.
#' @param criterion A character string defining the criterion used to select
#'   variables and FP models of different degrees.
#' @param select Named numeric vector with one significance threshold per
#'   conceptual term, used during backfitting to decide whether the term is
#'   retained.
#' @param alpha Named numeric vector with one functional-form significance
#'   threshold per conceptual term.
#' @param keep A character vector with names of variables to be kept in the
#'   model regardless of selection criteria.
#' @param xorder A string determining the order of entry of the covariates into
#'   the model-selection algorithm.
#' @param powers Named list of permitted FP powers, one element per conceptual
#'   term. Explicitly mapped terms use the fixed linear power \code{1}.
#' @param method A character string specifying the method for tie handling in
#'   Cox regression.
#' @param strata A factor of all possible combinations of stratification
#'   variables. Returned from \code{survival::strata()}.
#' @param nocenter A numeric vector with a list of values for fitting Cox
#'   models. See \code{survival::coxph()} for details.
#' @param acdx Named logical vector with one value per conceptual term,
#'   indicating which singleton continuous terms undergo the approximate
#'   cumulative distribution (ACD) transformation.
#' @param ftest Logical. If \code{TRUE} and \code{family = "gaussian"}, use an
#'   F-test rather than a chi-square likelihood-ratio test.
#' @param control A list with parameters for model fit. See
#'   \code{survival::coxph()} or \code{stats::glm()} for details.
#' @param zero Named logical vector with one value per conceptual term,
#'   indicating which singleton terms treat non-positive values as zero before
#'   FP transformation.
#' @param catzero Named logical vector with one value per conceptual term,
#'   indicating which singleton terms also receive a binary zero indicator.
#'   Internally, values \code{x <= 0} are recoded to zero and the indicator is
#'   equivalent to \code{I(original x <= 0)}.
#' @param spike Named logical vector with one value per conceptual term,
#'   indicating which singleton terms are assessed using the SAZ algorithm.
#' @param min_saz_component_prop Numeric in \eqn{(0, 0.5)}. Minimum required
#'   proportion in each component of a spike-at-zero covariate: the
#'   zero component and the positive component. A requested
#'   spike-at-zero variable is retained for SAZ modelling only if both component
#'   proportions are at least this value.
#' @param saz_pre_resolved Logical. If \code{TRUE}, spike-at-zero eligibility
#'   has already been resolved before shift/scale preprocessing by the caller.
#'   In that case, \code{fit_mfp()} does not call \code{reset_spike()} again.
#'   If \code{FALSE}, \code{fit_mfp()} performs a defensive late eligibility
#'   reset for internal callers that have not yet been updated.
#' @param force_max_fp Named logical vector with one value per conceptual term.
#'   If \code{TRUE} for a term, forces selection of the most complex functional form at the
#'   degree specified by \code{df}, bypassing AIC/BIC comparison against
#'   simpler forms. Has no effect when \code{criterion = "pvalue"}.
#' @param has_offset logical indicating whether an offset was specified before
#' missing offsets were replaced by zeros internally.
#' @param verbose Logical. If \code{TRUE}, progress information is printed
#'   during fitting. Default \code{FALSE}.
#'
#' @section Algorithm:
#' \enumerate{
#'   \item \strong{Full linear reference and variable ordering.} A model
#'     containing all candidate predictors as ordinary linear terms is fitted
#'     once. Its null and fitted-model deviances are retained as
#'     \code{null_deviance} and \code{linear_deviance}. When
#'     \code{xorder} is \code{"ascending"} or \code{"descending"} and
#'     more than one predictor is present, leave-one-predictor-out likelihood-
#'     ratio tests determine the visiting order. For \code{xorder = "original"}
#'     or a single predictor, no reduced ordering models are fitted.
#'   \item \strong{Pre-processing.} Initial FP powers are set to 1. ACD
#'     transformation setup, zero/catzero/spike handling, and spike eligibility
#'     checks are performed. See the \emph{Spike-at-zero handling} section.
#'   \item \strong{MFP backfitting cycles.} FP powers and spike decisions are
#'     updated iteratively via \code{find_best_fp_cycle()} until convergence or
#'     the maximum number of cycles is reached. Previously computed adjustment
#'     transformations are cached in \code{prev_adj_params} to avoid
#'     recomputation across cycles.
#'   \item \strong{Final transformation.} After convergence, \code{x} is
#'     backscaled (if \code{scale != 1}) to restore the shifted-but-not-scaled
#'     variable, then FP-transformed using the selected powers. Centering and
#'     binary indicator creation for \code{catzero} variables are applied.
#'   \item \strong{Model fitting.} The final model is fitted on the transformed
#'     design matrix and returned as an \code{mfp2} object.
#' }
#'
#' @section Spike-at-zero handling:
#' In normal public calls, spike-at-zero eligibility should already have been
#' resolved before shift/scale preprocessing. This is necessary because variables
#' reset from SAZ to ordinary FP may need ordinary FP shifting rather than the
#' shift-to-zero behaviour used for active zero/catzero/spike variables.
#'
#' Inside \code{fit_mfp()}, the following operations are applied:
#' \enumerate{
#'   \item The cascade is enforced: \code{catzero[spike] <- TRUE} and
#'     \code{zero[catzero] <- TRUE}.
#'   \item A temporary recoded copy of \code{x} is created for spike-specific
#'     checks and positive-part df capping.
#'   \item If \code{saz_pre_resolved = FALSE}, \code{reset_spike()} is called as
#'     a defensive fallback for internal callers that have not resolved SAZ
#'     eligibility before preprocessing.
#'   \item For retained spike variables, the maximum FP degrees of freedom may
#'     be reduced according to the number of distinct positive values. This
#'     mirrors the ordinary MFP df-capping rule, but applies it to the positive
#'     continuous component rather than to the full variable including zeros.
#' }
#'
#' @return An object of class \code{"mfp2"} built from the final fitted
#'   GLM or Cox model and augmented with MFP-specific metadata. In addition to
#'   the selected powers, transformation settings, and convergence information,
#'   the returned object contains three family-specific deviance values:
#'   \code{null_deviance}, \code{linear_deviance}, and
#'   \code{mfp_deviance}. For GLMs these are the null and residual deviances
#'   returned by the fitted GLM objects. For Cox models they are minus twice the
#'   corresponding null and fitted partial log-likelihoods. See \code{mfp2()}
#'   for the remaining components.
#'
#' @references
#' Sauerbrei, W. and Royston, P. (1999). Building multivariable prognostic
#' and diagnostic models: transformation of the predictors by using fractional
#' polynomials. \emph{Journal of the Royal Statistical Society: Series A},
#' 162, 71--94.
#'
#' Royston, P. and Sauerbrei, W. (2016). mfpa: Extension of mfp using the ACD
#' covariate transformation for enhanced parametric multivariable modeling.
#' \emph{The Stata Journal}, 16(1), 72--87.
#'
#' Lorenz, E., Jenkner, C., Sauerbrei, W., and Becher, H. (2017). Modeling
#' variables with a spike at zero: examples and practical recommendations.
#' \emph{American Journal of Epidemiology}, 185(8), 650--660.
#'
#' Lorenz, E., Jenkner, C., Sauerbrei, W., and Becher, H. (2019). Modeling
#' exposures with a spike at zero: simulation study and practical application
#' to survival data. \emph{Biostatistics & Epidemiology}, 3(1), 23--37.
#'
#' @seealso \code{mfp2()}, \code{find_best_fp_cycle()}, \code{reset_spike()}
#' @importFrom utils modifyList
#' @keywords internal
#' @noRd
fit_mfp <- function(x,
                    y,
                    weights,
                    offset,
                    cycles,
                    scale,
                    shift,
                    df,
                    center,
                    family,
                    family_string,
                    criterion,
                    select,
                    alpha,
                    keep,
                    xorder,
                    powers,
                    method,
                    strata,
                    nocenter,
                    acdx,
                    ftest,
                    control,
                    zero,
                    catzero,
                    spike,
                    min_saz_component_prop,
                    saz_pre_resolved = FALSE,
                    force_max_fp,
                    has_offset,
                    verbose,
                    term_to_columns = NULL) {
  
  # fit_mfp() implements the full MFP/MFPA/SAZ algorithm described in the
  # "Algorithm" section above: order variables, pre-process ACD/zero/catzero/
  # spike settings, run backfitting cycles to select FP powers, then transform
  # and fit the final model. mfp2.default()/mfp2.formula() are expected to have
  # already validated and shifted/scaled all inputs before calling this function.
  
  if (is.null(term_to_columns)) {
    term_to_columns <- stats::setNames(as.list(colnames(x)), colnames(x))
  }
  
  variables_x <- names(term_to_columns)
  
  # A non-trivial mapping exists when a conceptual term either spans multiple
  # raw columns or maps to a single raw column with a different name. The latter
  # occurs for two-level factors and binary group dummies (for example,
  # `treatment` -> `treatment1`). Both cases require raw-column expansion.
  has_mapped_terms <- any(vapply(
    variables_x,
    function(term) {
      cols <- term_to_columns[[term]]
      length(cols) != 1L || !identical(cols[[1L]], term)
    },
    logical(1L)
  ))
  
  # Step 1: Resolve the family object and optionally report initial df --------
  # Resolve GLM family objects once for all repeated internal model fits.
  # Public mfp2.default() already does this, but keeping it here makes direct
  # internal calls to fit_mfp() avoid repeated stats::gaussian()/binomial()/
  # poisson() construction as well. Cox remains the character string "cox".
  family_fit <- resolve_fit_model_family(family)
  
  # Print the starting df for each variable before the backfitting cycles
  # begin, so users can see what df each predictor was assigned prior to any
  # cardinality-based or SAZ-based capping performed later in this function.
  if (isTRUE(verbose)) {
    df_text <- utils::capture.output(
      print(
        matrix(df, nrow = 1, dimnames = list("df", variables_x)),
        quote = FALSE
      )
    )
    
    message(
      "Initial degrees of freedom:\n",
      paste(df_text, collapse = "\n")
    )
  }
  
  # Step 2: Fit the full linear reference and determine visiting order -------
  # order_variables() always fits the model containing every candidate
  # predictor as an ordinary linear term. That fit supplies both the null-model
  # deviance and the full-linear-model deviance. Its log-likelihood is separately
  # reused for leave-one-predictor-out likelihood-ratio tests when
  # significance-based ordering is requested.
  #
  # The call is deliberately unconditional: with one predictor, or with
  # xorder = "original", no reduced ordering models are needed, but the full
  # linear reference deviance is still required for the fitted mfp2 object.
  ordering_result <- order_variables(
    xorder        = xorder, 
    x             = x,
    term_to_columns = term_to_columns,
    y             = y,
    family        = family_fit,
    family_string = family_string,
    weights       = weights,
    offset        = offset,
    strata        = strata, 
    method        = method,
    control       = control,
    nocenter      = nocenter
  )
  
  variables_ordered <- ordering_result$variables_ordered
  null_deviance <- ordering_result$null_deviance
  linear_deviance <- ordering_result$linear_deviance
  
  # Log-likelihood-scale quantities for the Model Fit block of summary.mfp2().
  # The linear-reference numbers come from the full-linear fit that
  # order_variables() computed once at the top of the algorithm. `null_logl`
  # is populated for Cox (returned free by coxph.fit); for GLMs it is filled
  # in below, once, after backfitting has converged.
  linear_logl <- ordering_result$linear_logl
  linear_df   <- ordering_result$linear_df
  null_logl   <- ordering_result$null_logl
  
  if (verbose) {
    message(sprintf(
      "Visiting order: %s",
      paste0(variables_ordered, collapse = ", ")
    ))
  }
  
  # Step 3: Initialize FP powers and align all per-variable vectors/x to order
  # Every variable starts the first cycle as linear (power = 1); the
  # backfitting cycles below will update powers_current as better-fitting FP
  # forms are found. All per-variable option vectors (and x itself) are then
  # reordered to variables_ordered so the rest of the function can index them
  # consistently by name/position.
  powers_current <- stats::setNames(
    as.list(rep(1, length(variables_ordered))),
    variables_ordered
  )
  
  # Re-order all per-variable vectors to match variables_ordered
  alpha        <- setNames(alpha,   variables_x)[variables_ordered]
  select       <- setNames(select,  variables_x)[variables_ordered]
  df           <- setNames(df,      variables_x)[variables_ordered]
  center       <- setNames(center,  variables_x)[variables_ordered]
  shift        <- if (!is.null(names(shift))) {
    shift[variables_ordered]
  } else {
    setNames(shift, variables_x)[variables_ordered]
  }
  scale        <- if (!is.null(names(scale))) {
    scale[variables_ordered]
  } else {
    setNames(scale, variables_x)[variables_ordered]
  }
  acdx         <- setNames(acdx,    variables_x)[variables_ordered]
  zero         <- setNames(zero,    variables_x)[variables_ordered]
  catzero      <- setNames(catzero, variables_x)[variables_ordered]
  spike        <- setNames(spike,   variables_x)[variables_ordered]
  force_max_fp <- setNames(force_max_fp, variables_x)[variables_ordered]
  powers       <- powers[variables_ordered]
  
  # Reorder the raw design columns by conceptual-term visiting order. For the
  # historical singleton case this is identical to the previous column reorder.
  term_to_columns <- term_to_columns[variables_ordered]
  raw_columns_ordered <- unlist(term_to_columns, use.names = FALSE)
  if (!identical(colnames(x), raw_columns_ordered)) {
    x <- x[, raw_columns_ordered, drop = FALSE]
  }
  if (!identical(colnames(x), raw_columns_ordered)) {
    stop("Internal error: x column order does not match term_to_columns.",
         call. = FALSE)
  }
  
  
  # `keep` forces variables into the final model regardless of their
  # backfitting significance, by setting their selection threshold to 1
  # (i.e. always significant) for the p-value criterion.
  if (!is.null(keep)) {
    select[which(names(select) %in% keep)] <- 1
  }
  
  # Step 4: Configure the ACD (approximate cumulative distribution) transform -
  # Requesting ACD for a variable forces df = 4, since the FSPA function-
  # selection procedure for ACD needs the full FP1(p1, p2) model space (see
  # "Details on approximate cumulative distribution transformation" in mfp2()).
  # Initial powers_current for ACD variables are set to c(1, NA): a linear
  # term for x and no ACD term yet, to be updated in the first backfitting cycle.
  if (any(acdx)) {
    acdx           <- reset_acd(x, acdx)
    variables_acd  <- names(acdx)[acdx]
    powers_current <- utils::modifyList(
      powers_current,
      sapply(variables_acd, function(v) c(1, NA), simplify = FALSE)
    )
    df[variables_ordered %in% variables_acd] <- 4
  }
  
  # Step 5: Resolve the spike/catzero/zero cascade and SAZ eligibility --------
  # Public callers should resolve SAZ eligibility before shift/scale
  # preprocessing. This prevents variables reset from SAZ to ordinary FP from
  # retaining a shift = 0 that was chosen only because spike was temporarily TRUE.
  #
  # fit_mfp() still enforces the spike/catzero/zero cascade. It also keeps a
  # defensive reset path for internal callers that have not yet pre-resolved SAZ
  # eligibility.
  
  user_catzero <- catzero   # pre-cascade user intent
  user_zero    <- zero      # pre-cascade user intent
  
  catzero[spike] <- TRUE    # spike implies catzero
  zero[catzero]  <- TRUE    # catzero implies zero
  
  # Temporary recoded data used only for spike eligibility checks.
  # Do not mutate the real x here.
  x_for_spike <- x
  
  zero_aligned_for_spike <- if (has_mapped_terms) {
    stats::setNames(
      rep(unname(zero[names(term_to_columns)]), lengths(term_to_columns)),
      unlist(term_to_columns, use.names = FALSE)
    )[colnames(x_for_spike)]
  } else {
    zero[colnames(x_for_spike)]
  }
  cols_to_zero_for_spike <- which(zero_aligned_for_spike)
  
  if (length(cols_to_zero_for_spike) > 0L) {
    for (j in cols_to_zero_for_spike) {
      x_for_spike[x_for_spike[, j] <= 0, j] <- 0
    }
  }
  
  # Defensive reset for ineligible spike variables.
  #
  # Public callers such as mfp2.default() should resolve SAZ eligibility before
  # shift/scale preprocessing and call fit_mfp(..., saz_pre_resolved = TRUE).
  # That avoids the unsafe situation where a spike variable is reset to ordinary
  # FP after shift = 0 has already been chosen.
  #
  # This late reset is kept only for internal callers that have not yet
  # pre-resolved SAZ eligibility.
  if (!isTRUE(saz_pre_resolved) && any(spike)) {
    result <- reset_spike(
      x                      = x_for_spike,
      spike                  = spike,
      user_catzero           = user_catzero,
      user_zero              = user_zero,
      min_saz_component_prop = min_saz_component_prop
    )
    
    spike   <- result$spike
    catzero <- result$catzero
    zero    <- result$zero
  }
  
  # For retained spike-at-zero variables, cap the maximum FP df using only the
  # positive component. Ordinary assign_df() uses the full variable,
  # but for SAZ the relevant information for the FP part is x > 0 after
  # temporary zero recoding.
  df <- cap_spike_df(
    x     = x_for_spike,
    df    = df,
    spike = spike
  )
  
  # Validate candidate powers after all early df modifications. This includes
  # ACD forcing df = 4 and SAZ-specific df capping for retained spike variables.
  validate_mfp_candidate_powers(
    powers = powers,
    df = df
  )
  
  # Spike decision initialisation.
  # continuous_only means standard FP algorithm by default.
  spike_decision        <- rep(
    saz_decision_codes[["continuous_only"]],
    length(variables_ordered)
  )
  names(spike_decision) <- variables_ordered
  
  # Step 6: Recode real x for zero handling and build catzero indicator matrices
  # Now that reset_spike() has produced the final zero/catzero/spike vectors, we
  # can safely mutate the actual x used by the MFP cycles.
  #
  # Only variables with final zero == TRUE are recoded. Therefore, variables whose
  # spike request was rejected and whose user-specified zero/catzero status was
  # FALSE remain untouched.
  
  zero_x       <- zero
  zero_aligned <- if (has_mapped_terms) {
    stats::setNames(
      rep(unname(zero[names(term_to_columns)]), lengths(term_to_columns)),
      unlist(term_to_columns, use.names = FALSE)
    )[colnames(x)]
  } else {
    zero[colnames(x)]
  }
  cols_to_zero <- which(zero_aligned)
  
  if (length(cols_to_zero) > 0L) {
    for (j in cols_to_zero) {
      x[x[, j] <= 0, j] <- 0
    }
  }
  
  # zero_x is intentionally kept identical to zero (not reset to FALSE), even
  # though x has already been physically recoded to 0 above. The backfitting
  # cycles below (find_best_fp_cycle() / transform_data_step()) still need
  # this flag at cycle time to know which variables have a zero component, so
  # that FP transformations during variable/degree selection are computed only
  # over the positive part while the recoded zeros are left untouched.
  
  # Step 7: Build catzero binary-indicator matrices and cache ACD parameters -
  # Binary zero indicators for catzero variables.
  #
  # Because non-positive values have already been recoded to zero above,
  # x == 0 here means:
  #
  #   I(original x <= 0)
  #
  # not merely:
  #
  #   I(original x == 0)
  #
  # This is the intended interpretation of catzero/spike variables.
  catzero_mat_list <- lapply(names(catzero), function(v) {
    if (!isTRUE(catzero[[v]])) {
      return(NULL)
    }
    
    matrix(
      as.integer(x[, v] <= 0),
      ncol = 1L,
      dimnames = list(rownames(x), "catzero")
    )
  })
  
  names(catzero_mat_list) <- names(catzero)
  
  # ACD parameters (beta0, beta1, power, shift, scale for the rank-based
  # power-linear approximation) are estimated once per ACD variable here and
  # cached, rather than being refit on every backfitting cycle/step that needs
  # to apply the ACD transformation.
  acd_parameter <- lapply(names(acdx), function(v) {
    if (isTRUE(acdx[[v]])) {
      fit_acd(
        x      = x[, v],
        powers = powers[[v]],
        zero   = zero[[v]]
      )
    } else {
      NULL
    }
  })
  
  names(acd_parameter) <- names(acdx)
  
  # Number of events (Cox) or observations (all other families), used for
  # AIC/BIC-based selection and to guard against fitting a Cox model with no
  # observed events.
  if (family_string == "cox") {
    status <- y[, ncol(y)]
    n_obs <- sum(!is.na(status) & status > 0)
    
    if (!is.finite(n_obs) || n_obs <= 0L) {
      stop(
        "Cox selection requires at least one observed event.",
        call. = FALSE
      )
    }
  } else {
    n_obs <- nrow(x)
  }
  
  # Step 8: Run MFP backfitting cycles until convergence -----------------------
  # A cycle is one complete pass through all variables (see find_best_fp_cycle()
  # documentation); convergence means neither the selected powers nor the SAZ
  # stage-2 decisions changed compared to the previous cycle.
  j         <- 1L
  converged <- FALSE
  
  prev_adj_params       <- vector("list", length = length(variables_ordered))
  names(prev_adj_params) <- variables_ordered
  
  while (j <= cycles) {
    if (verbose) {
      message(sprintf(
        "%s\nRunning MFP Cycle %d\n%s",
        strrep("-", 21), j, strrep("-", 21)
      ))
    }
    
    fit_best_cycle <- find_best_fp_cycle(
      x               = x,
      term_to_columns = term_to_columns,
      y               = y,
      powers_current  = powers_current,
      df              = df,
      weights         = weights,
      offset          = offset,
      family          = family_fit,
      family_string   = family_string,
      criterion       = criterion,
      select          = select,
      alpha           = alpha,
      keep            = keep,
      powers          = powers,
      ftest           = ftest,
      control         = control,
      rownames        = rownames(x),
      strata          = strata,
      nocenter        = nocenter,
      method          = method,
      acdx            = acdx,
      zero            = zero_x,        
      catzero         = catzero_mat_list,  # named list of binary indicators
      spike           = spike,
      spike_decision  = spike_decision,
      acd_parameter   = acd_parameter,
      prev_adj_params = prev_adj_params,
      force_max_fp    = force_max_fp,
      has_offset      = has_offset,
      n_obs           = n_obs,
      verbose         = verbose
    )
    
    powers_updated        <- fit_best_cycle$powers_current
    spike_decision_updated <- fit_best_cycle$spike_decision
    prev_adj_params       <- fit_best_cycle$prev_adj_params
    
    # Compare powers after normalizing away the stored-but-inactive continuous
    # power of binary-only spike variables (spike_decision = 3), since that
    # power does not affect the fitted adjustment matrix and would otherwise
    # cause spurious non-convergence.
    powers_same <- identical(
      normalize_powers_for_convergence(powers_current, spike_decision),
      normalize_powers_for_convergence(powers_updated, spike_decision_updated)
    )
    
    spike_same <- identical(spike_decision, spike_decision_updated)
    
    if (powers_same && spike_same) {
      converged <- TRUE
      if (verbose)
        message(sprintf(
          "Fractional polynomial fitting algorithm converged after %d cycle(s).",
          j
        ))
      break
    } else {
      powers_current <- powers_updated
      spike_decision <- spike_decision_updated
      j <- j + 1L
    }
  }
  
  if (!converged) {
    warning(
      sprintf("i No convergence after %d cycles.", cycles),
      " Results of the last iteration are reported.",
      call. = FALSE
    )
  }
  
  # Step 9: Apply the final FP/ACD transformation and centering ---------------
  # Backscale x before final FP transformation so that coefficients are on
  # the phi(x + shift) scale, matching what a user would expect from mfp2().
  # Scaling was applied upstream for numerical stability during power selection
  # and has no effect on power selection itself.
  if (any(scale != 1)) {
    scale_for_columns <- if (has_mapped_terms) {
      stats::setNames(
        rep(unname(scale[names(term_to_columns)]), lengths(term_to_columns)),
        unlist(term_to_columns, use.names = FALSE)
      )
    } else {
      scale
    }
    x <- backscale_matrix(x, scale_for_columns)
  }
  
  # ACD parameters are estimated on the scaled working x.
  # Final fitting and prediction later pass shifted-but-not-scaled x into
  # transform_matrix(). Store the training scale so apply_acd() reconstructs
  # the same scaled ACD input before using beta0, beta1, and power.
  acd_parameter_final <- acd_parameter
  
  for (v in names(acd_parameter_final)) {
    if (!is.null(acd_parameter_final[[v]])) {
      # Remove training-data ACD values. fit_acd() returns $acd as the
      # transformed training vector, but only beta0/beta1/power/shift/scale
      # are needed for apply_acd() during prediction. Keeping $acd wastes
      # memory and can cause confusion.
      acd_parameter_final[[v]]$acd <- NULL
      acd_parameter_final[[v]]$shift <- 0
      acd_parameter_final[[v]]$scale <- scale[[v]]
    }
  }
  
  if (has_mapped_terms) {
    expanded_final <- expand_term_metadata_to_columns(
      term_to_columns = term_to_columns,
      powers = powers_current,
      center = center,
      acdx = acdx,
      zero = zero,
      catzero = catzero,
      spike = spike,
      spike_decision = spike_decision,
      acd_parameter = acd_parameter_final
    )
    power_list_final <- expanded_final$powers
    center_final <- expanded_final$center
    acdx_final <- expanded_final$acdx
    zero_final <- expanded_final$zero
    catzero_final <- expanded_final$catzero
    spike_final <- expanded_final$spike
    spike_decision_final <- expanded_final$spike_decision
    acd_parameter_final_columns <- expanded_final$acd_parameter
  } else {
    power_list_final <- powers_current
    center_final <- center
    acdx_final <- acdx
    zero_final <- zero
    catzero_final <- catzero
    spike_final <- spike
    spike_decision_final <- spike_decision
    acd_parameter_final_columns <- acd_parameter_final
  }
  
  data_transformed <- transform_matrix(
    x                  = x,
    power_list         = power_list_final,
    center             = center_final,
    acdx               = acdx_final,
    acd_parameter_list = acd_parameter_final_columns,
    zero               = zero_final,
    catzero            = catzero_final,
    spike              = spike_final,
    reset_zero         = FALSE,
    spike_decision     = spike_decision_final
  )
  
  # Update catzero for final metadata.
  # Final-time `catzero` is logical metadata. It should reflect whether a
  # *_bin column is present in the final design matrix.
  catzero_effective <- catzero
  
  for (v in names(catzero_effective)) {
    if (isTRUE(spike[[v]])) {
      if (as.integer(spike_decision[[v]]) ==
          saz_decision_codes[["continuous_only"]]) {
        # Spike continuous-only or null: no spike *_bin column.
        catzero_effective[[v]] <- FALSE
      } else if (as.integer(spike_decision[[v]]) ==
                 saz_decision_codes[["binary_only"]]) {
        # Spike binary-only: *_bin is the selected term.
        catzero_effective[[v]] <- TRUE
      }
    } else if (all(is.na(powers_current[[v]]))) {
      # Ordinary catzero is tied to the parent variable.
      # If the parent is eliminated, its ordinary *_bin column is absent too.
      catzero_effective[[v]] <- FALSE
    }
  }
  
  # Update catzero_list to match final effective catzero status
  for (v in names(catzero_effective)) {
    if (!isTRUE(catzero_effective[[v]])) {
      catzero_mat_list[[v]] <- NULL
    }
  }
  
  # Binary-only spike variables (spike_decision = 3) contribute no continuous
  # FP/ACD term to the final model, only the *_bin indicator handled above.
  # Their stored continuous power is no longer meaningful past this point, so
  # it is set to NA here to keep fp_powers/fp_terms metadata consistent with
  # the actual final design matrix.
  for (v in names(spike_decision)) {
    if (spike[v] && spike_decision[v] == saz_decision_codes[["binary_only"]]) {
      powers_current[[v]] <- NA
    }
  }
  
  # Step 10: Fit the final model and assemble the returned mfp2 object --------
  modelfit <- fit_model(
    x             = data_transformed$x_transformed,
    y             = y,
    family        = family_fit,
    family_string = family_string,
    weights       = weights,
    offset        = offset,
    method        = method,
    strata        = strata,
    control       = control,
    rownames      = rownames(data_transformed$x_transformed),
    nocenter      = nocenter,
    fast          = FALSE,
    has_offset    = has_offset
  )
  
  # predict.coxph() expresses prediction offsets relative to the mean training
  # offset, independently of the covariate reference. This is final-model
  # metadata, so calculate it once here after the MFP backfitting cycles have
  # terminated and the final Cox model has been fitted. Internal candidate fits
  # do not need to calculate or store it.
  cox_offset_reference <- if (identical(family_string, "cox")) {
    if (isTRUE(has_offset)) unname(mean(offset)) else 0
  } else {
    NULL
  }
  
  # fit_model() stores the exact relationship between transformed source
  # columns and fitted coefficient names for every full formula-based fit.
  transformed_to_model_columns <- modelfit$transformed_to_model_columns
  if (is.null(transformed_to_model_columns)) {
    stop(
      "Internal error: final model lacks transformed-to-model column metadata.",
      call. = FALSE
    )
  }
  
  # The final fit supplies the family-specific deviance for the selected MFP
  # model. For GLMs this is the residual deviance; for Cox models it is minus
  # twice the fitted partial log-likelihood.
  mfp_deviance <- modelfit$model_deviance
  
  # -- Model Fit display quantities -------------------------------------------
  # Populate the log-likelihood-scale numbers used by summary.mfp2()'s Model
  # Fit block. All computation is done here, once, after backfitting has
  # converged -- never inside the MFP candidate-fitting loop.
  #
  # For Cox, null_logl was returned free by coxph.fit() during the
  # linear-reference fit. For GLMs, fit_model() does not compute it (that
  # would cost an intercept-only refit inside every interior candidate fit),
  # so we do it here on demand: a single stats::glm.fit() call on an
  # intercept-only design, with the same weights and offset used by the MFP
  # fit. A rare failure yields NA rather than aborting the summary.
  if (identical(family_string, "cox")) {
    # already populated from ordering_result$null_logl above
  } else if (is.null(null_logl) || !is.finite(null_logl)) {
    null_logl <- tryCatch({
      null_x <- matrix(
        rep.int(1, NROW(y)),
        ncol = 1L,
        dimnames = list(NULL, "(Intercept)")
      )
      null_weights <- if (is.null(weights)) rep.int(1, NROW(y)) else weights
      null_offset  <- if (is.null(offset))  rep.int(0, NROW(y)) else offset
      null_fit <- stats::glm.fit(
        x = null_x,
        y = y,
        family = family_fit,
        weights = null_weights,
        offset = null_offset
      )
      null_df_glm <- if (null_fit$family$family == "gaussian") 2L else 1L
      unname(null_df_glm - null_fit$aic / 2)
    }, error = function(e) NA_real_)
  }
  
  mfp_logl <- modelfit$logl
  mfp_df   <- modelfit$df
  
  # Assemble the mfp2 object: start from the glm/coxph fit object returned by
  # fit_model() and layer on MFP-specific metadata (selected powers, shift/
  # scale/center, zero/catzero/spike status, convergence, etc.) expected by
  # print.mfp2()/summary.mfp2()/predict.mfp2().
  fit <- utils::modifyList(
    modelfit$fit,
    list(
      centers         = data_transformed$centers,
      acd_parameter   = data_transformed$acd_parameter,
      convergence_mfp = converged,
      null_deviance   = null_deviance,
      linear_deviance = linear_deviance,
      mfp_deviance    = mfp_deviance,
      # Log-likelihood-scale quantities used by summary.mfp2()'s Model Fit
      # block. All three -2 log L values are on the same convention, so the
      # summary can render them directly without any per-family adjustment.
      null_logl       = null_logl,
      linear_logl     = linear_logl,
      linear_df       = linear_df,
      mfp_logl        = mfp_logl,
      mfp_df          = mfp_df,
      x_original = if (has_mapped_terms) {
        selected_terms <- names(powers_current)[
          !vapply(powers_current, function(p) all(is.na(p)), logical(1L)) |
            (spike & spike_decision == saz_decision_codes[["binary_only"]])
        ]
        x[, unlist(term_to_columns[selected_terms], use.names = FALSE), drop = FALSE]
      } else {
        x[, names(powers_current[
          !sapply(powers_current, function(p) all(is.na(p))) |
            (spike & spike_decision == saz_decision_codes[["binary_only"]])
        ]), drop = FALSE]
      },
      y_original = y,
      fp_terms        = create_fp_terms(
        powers_current, acdx, df, select, alpha, criterion, zero,
        catzero_effective, spike, spike_decision,
        term_to_columns = term_to_columns,
        transformed_to_model_columns = transformed_to_model_columns,
        coefficients = modelfit$fit$coefficients
      ),
      transformations = data.frame(shift = shift, scale = scale, center = center),
      fp_powers       = powers_current,
      acd             = acdx,
      zero            = zero,
      catzero         = catzero_effective,
      catzero_list = catzero_mat_list,
      spike              = spike,
      spike_dec       = spike_decision,
      transformed_to_model_columns = transformed_to_model_columns,
      cox_offset_reference = cox_offset_reference
    )
  )
  
  # Store the complete conceptual-term lookup for every fit. Identity mappings
  # are required by prediction just as explicit categorical mappings are, and
  # keeping one invariant avoids reconstructing term structure downstream.
  fit$term_to_columns <- term_to_columns
  
  class(fit) <- c("mfp2", class(fit))
  fit
}

#' Helper to run cycles of the mfp algorithm 
#' 
#' This function estimates the best FP functions for all predictors in the 
#' current cycle. To be used in \code{fit_mfp()}.
#' 
#' @details 
#' A cycle is defined as a complete pass through all the predictors in the input
#' matrix `x`, while a step is defined as the assessment of a single predictor. 
#' This algorithm is described in Sauerbrei et al. (2006) and given in detail
#' in Royston and Sauerbrei (2008), in particular chapter 6.
#' 
#' Briefly, a cycle works as follows: it takes as input the data matrix along with
#' a set of current best fp powers for each variable. In each step, the fp
#' powers of a single covariate are assessed, while adjusting for other
#' covariates. Adjustment variables are transformed using their current
#' fp powers (this is done in \code{transform_data_step()} and the fp powers 
#' of the variable of interest are tested using the closed test procedure
#' (conducted in \code{find_best_fp_step()}).
#' Some of the adjustment variables may have their fp power set to `NA`, 
#' which means they were not selected from the working model and are not used
#' in that step. The results from all steps are returned, completing a cycle.
#' 
#' Note that in each cycle every variable is evaluated.This includes variables
#' that may have been eliminated in previous cycles. They will re-enter each
#' new cycle for potential inclusion in the working model or to be re-evaluated
#' for elimination.
#' 
#' The current adjustment set is always given through the current fp powers, 
#' which are updated in each step (denoted as `powers_current`). 
#'
#' If \code{catzero} variables are supplied, the algorithm will automatically create 
#' the corresponding binary variables and include them in the model. Additionally, 
#' each binary variable and its associated continuous variable will be treated as 
#' one predictor, and they will be tested jointly for inclusion in the model.
#'  
#' 
#' @references 
#' Royston, P. and Sauerbrei, W., 2008. \emph{Multivariable Model - Building: 
#' A Pragmatic Approach to Regression Anaylsis based on Fractional Polynomials 
#' for Modelling Continuous Variables. John Wiley & Sons.}\cr
#' Sauerbrei, W., Meier-Hirmer, C., Benner, A. and Royston, P., 2006. 
#' \emph{Multivariable regression model building by using fractional 
#' polynomials: Description of SAS, STATA and R programs. 
#' Comput Stat Data Anal, 50(12): 3464-85.}
#' Sauerbrei, W. and Royston, P., 1999. \emph{Building multivariable prognostic 
#' and diagnostic models: transformation of the predictors by using fractional 
#' polynomials. J Roy Stat Soc a Sta, 162:71-94.}
#' 
#' @inheritParams fit_mfp
#' @param powers_current a list of length equal to the number of variables, 
#' indicating the fp powers to be used in the current step for all variables 
#' (except `xi`). 
#' @param catzero A named list of binary indicator variables of length \code{ncol(x)} 
#' for nonpositive values, created when specific variables are passed to the 
#' \code{catzero} argument of \code{fit_mfp}. If an element of the list is 
#' \code{NULL}, it indicates that the corresponding variable was not specified by
#' the user in the \code{catzero} argument of \code{fit_mfp}. Here, \code{catzero}
#' is a list of binary variables, not a named logical vector as in \code{fit_mfp}.
#' @param acd_parameter Named list of ACD parameters produced by `fit_acd()`, 
#' with length equal to \code{ncol(x)}. Each list element corresponds to a variable; 
#' if an element is \code{NULL}, the variable was not specified in the 
#' \code{acdx} argument of \code{fit_mfp}.
#' @param spike_decision Named vector indicating how spike-at-zero (SAZ) 
#' variables are handled. Each element corresponds to a variable and encodes 
#' the selected strategy: `1` = include FP for positive values plus binary SAZ, 
#' `2` = treat as continuous FP only, `3` = include binary SAZ only.
#' @param rownames passed to \code{survival::coxph.fit()}.
#' @param prev_adj_params Named list used to store previously computed adjustment 
#' variable transformations. This is updated at each step and reused in the next 
#' cycle to avoid recomputation.
#' @param has_offset logical indicating whether an offset was specified before
#' missing offsets were replaced by zeros internally.
#' @param n_obs Numeric; number of observations (or observed events, for Cox
#' models), used for AIC/BIC-based selection and small-sample F-tests. Computed
#' once in \code{fit_mfp()} and passed through unchanged across cycles.
#' 
#' @return 
#' A list with updated components `powers_current` (current FP powers for all 
#' variables), `spike_decision` (updated spike-at-zero decisions), and
#' `prev_adj_params` (adjustment variable transformations to be used in the next
#' cycle).
#' @keywords internal
#' @noRd
find_best_fp_cycle <- function(x,
                               y, 
                               powers_current, 
                               df, 
                               weights, 
                               offset, 
                               family, 
                               family_string, 
                               criterion,
                               select, 
                               alpha, 
                               keep, 
                               powers, 
                               method, 
                               strata, 
                               verbose, 
                               ftest, 
                               control,
                               rownames,
                               nocenter,
                               zero,
                               catzero, 
                               spike,
                               spike_decision,
                               acd_parameter,
                               acdx,
                               prev_adj_params,
                               force_max_fp,
                               has_offset,
                               n_obs,
                               term_to_columns = NULL
) {
  
  if (is.null(term_to_columns)) {
    term_to_columns <- stats::setNames(as.list(colnames(x)), colnames(x))
  }
  
  # Variable visiting order within the cycle is fixed by the names of
  # powers_current (set once in fit_mfp() according to `xorder`); it does not
  # change from cycle to cycle.
  names_x <- names(powers_current)
  
  # Step through every variable once, updating its FP power/spike decision in
  # turn while adjusting for the *current* powers of all other variables.
  # Because powers_current is updated in place after each step, later
  # variables in this same pass are adjusted using already-updated powers of
  # earlier variables (Gauss-Seidel-style backfitting, not Jacobi-style).
  for (xi in names_x) {
    # Assess xi's best FP form (NA / linear / FP1 / FP2 / ...) via the closed
    # test procedure, holding the other variables' current powers fixed as
    # the adjustment set.
    fit_best_fp_step <- find_best_fp_step(
      x = x, # raw columns are resolved through term_to_columns
      term_to_columns = term_to_columns,
      y = y,
      xi = xi,
      powers_current = powers_current,
      weights = weights,
      offset = offset,
      df = df[xi], 
      select = select[xi], 
      alpha = alpha[xi],
      keep = keep,
      family = family,
      family_string = family_string,
      criterion = criterion,
      powers = powers,
      method = method,
      strata = strata,
      ftest = ftest,
      control = control,
      rownames = rownames,
      nocenter = nocenter,
      acdx = acdx,
      zero = zero,
      catzero = catzero,
      spike = spike,
      spike_decision = spike_decision,
      acd_parameter = acd_parameter,
      prev_adj_params = prev_adj_params,
      force_max_fp = force_max_fp,
      has_offset   = has_offset,
      n_obs        = n_obs,
      verbose = verbose
    )
    # Explicitly mapped terms have one conceptual linear/null state even when
    # a binary factor contributes only one non-identity dummy column.
    power_best <- fit_best_fp_step$power_best
    if (term_uses_column_mapping(xi, term_to_columns[[xi]])) {
      power_best <- if (all(is.na(power_best))) NA_real_ else 1
    }
    
    # Update powers_current and spike_decision in place so that the next
    # variable in this loop (and the next cycle) sees xi's latest result.
    powers_current[[xi]] <- power_best
    spike_decision       <- fit_best_fp_step$spike_decision
    # Cache xi's adjustment-variable transformations, keyed by xi, so the next
    # cycle can reuse them instead of recomputing from scratch.
    prev_adj_params[[xi]] <- fit_best_fp_step$current_adj_params[[xi]]
  }
  
  list(powers_current  = powers_current,
       spike_decision  = spike_decision, 
       prev_adj_params = prev_adj_params)
}

#' Normalize selected powers before checking convergence
#'
#' Helper used during the MFP backfitting cycle to compare selected powers
#' between successive cycles. For spike-at-zero variables selected as
#' binary-only (`spike_decision = 3`), the stored continuous power is ignored
#' because it does not affect the fitted adjustment matrix. Such powers are
#' normalized to `NA` before comparison.
#'
#' This function is only for convergence checking. It should not be used to
#' mutate `powers_current` during the fitting cycle, because cycle-time
#' transformation code may still require the stored stage-1 power to route
#' binary-only spike variables correctly.
#'
#' @param powers Named list of selected power vectors, keyed by parent variable
#'   name.
#' @param spike_decision Named integer vector/list of spike-at-zero decisions,
#'   using the same variable namespace as `powers`. Decision `3` means
#'   binary-only spike.
#'
#' @return Named list with the same names as `powers`, where powers for
#'   `spike_decision = 3` variables are replaced by `NA_real_`, and all other
#'   power vectors are unnamed.
#'
#' @keywords internal
#' @noRd
normalize_powers_for_convergence <- function(powers, spike_decision) {
  # For each variable, replace the stored power with NA if it is binary-only
  # spike (decision 3, where the continuous power is inert metadata), and
  # strip names from the power vector otherwise so that comparisons via
  # identical() aren't tripped up by incidental name differences.
  out <- Map(
    function(power, decision) {
      if (identical(
        as.integer(decision),
        saz_decision_codes[["binary_only"]]
      )) {
        NA_real_
      } else {
        unname(power)
      }
    },
    powers,
    spike_decision[names(powers)]
  )
  
  names(out) <- names(powers)
  out
}

#' Calculate degrees of freedom for a transformed variable
#'
#' Helper function used in `fit_mfp()` to determine the final number of
#' degrees of freedom (df) contributed by a variable, depending on its selected
#' powers, final/effective catzero status, and spike-at-zero decision.
#'
#' @param powers Numeric vector of selected powers for the variable. Can
#' contain `NA` values, for example for inactive ACD components. If all values
#' are `NA`, the continuous transformed component is inactive.
#' @param spike_decision Integer scalar (1, 2, or 3) specifying spike-at-zero
#'   handling:
#'   * `1` - include both the continuous transformed term(s) and the binary
#'   * `2` - include only the continuous transformed term(s).
#'   * `3` - include only the binary spike-at-zero indicator.
#' @param catzero Logical scalar indicating whether the final/effective
#'   catzero binary indicator is included for this variable. This should be the
#'   final metadata after applying spike-at-zero decisions, not the cycle-time
#'   `catzero` list of binary vectors.
#'
#' @details
#' The package uses the following df convention:
#' * If `spike_decision = 3`, df = 1 because the variable contributes only the
#'   binary spike-at-zero indicator.
#' * Otherwise, if all entries in `powers` are `NA`, df = 0 because the variable
#'   contributes no continuous component and no binary-only spike component.
#' * If the variable is modeled linearly, exactly `powers = 1`, df = 1.
#' * Otherwise, for fractional polynomials of degree *m*, where *m* is the
#'   number of non-`NA` powers, df = 2 * m.
#' * If `catzero = TRUE`, one additional df is added for the binary catzero
#'   indicator. This covers both ordinary non-spike catzero variables and
#'   spike-at-zero variables with `spike_decision = 1`. Binary-only spike
#'   variables are handled by the first rule and are not double-counted.
#'
#' Examples: if `powers = c(1, 2)` and `spike_decision = 2`, then df = 4.
#' If `powers = NA` and `spike_decision = 2`, then df = 0. If
#' `spike_decision = 3`, then df = 1.
#'
#' @return Integer scalar giving the degrees of freedom for the variable.
#'
#' @examples
#' \dontrun{
#' calculate_df(c(1, 2), 2, FALSE)  # df = 4
#' calculate_df(1, 1, TRUE)         # df = 2: linear + binary indicator
#' calculate_df(c(NA, NA), 2, FALSE)# df = 0: unselected variable
#' calculate_df(2, 3, TRUE)         # df = 1: binary spike only
#' calculate_df(0.5, 2, TRUE)       # df = 3: nonlinear FP1 + catzero
#' }
#' @keywords internal
#' @noRd
calculate_df <- function(powers, spike_decision, catzero = FALSE) {
  
  # Step 1: Validate scalar inputs -----------------------------------------
  if (length(spike_decision) != 1L ||
      is.na(spike_decision) ||
      !(as.integer(spike_decision) %in% unname(saz_decision_codes))) {
    stop("`spike_decision` must be a single integer 1, 2, or 3.")
  }
  
  if (length(catzero) != 1L || is.na(catzero) || !is.logical(catzero)) {
    stop("`catzero` must be a single TRUE/FALSE value.")
  }
  
  spike_decision <- as.integer(spike_decision)
  
  # Step 2: Apply the df rules in the order documented above --------------
  # Binary-only spike: exactly one binary indicator column.
  if (spike_decision == saz_decision_codes[["binary_only"]]) {
    return(1L)
  }
  
  powers <- as.numeric(powers)
  
  # Unselected variable.
  if (all(is.na(powers))) {
    return(0L)
  }
  
  # Linear model contributes 1 df; an FP of degree m (m non-NA powers)
  # contributes 2*m df, regardless of whether any powers are repeated.
  p <- as.numeric(powers[!is.na(powers)])
  
  df <- if (length(p) == 1L && p == 1) {
    1L
  } else {
    2L * length(p)
  }
  
  # `catzero` is the final/effective binary-indicator flag.
  # It already covers ordinary catzero variables and spike_decision == 1.
  if (isTRUE(catzero)) {
    df <- df + 1L
  }
  
  df
}

#' Helper to convert a nested list with same or different length into a matrix
#' 
#' Converts the per-variable list of selected FP powers (which may have a
#' different length for each variable, e.g. 1 power for FP1, 2 for FP2) into a
#' single rectangular matrix suitable for use as columns in
#' \code{create_fp_terms()}'s output data frame.
#' 
#' To be used in \code{fit_mfp()}.
#' 
#' @param power_list Named list of numeric power vectors, one element per
#'   variable, as stored in \code{fit_mfp()}'s \code{powers_current}. Elements
#'   may have different lengths (e.g. length 1 for a linear/FP1 variable,
#'   length 2 for FP2).
#' 
#' @return 
#' A numeric matrix with one row per variable (in the order of
#' \code{power_list}, row names taken from its names) and one column per power
#' slot up to the largest FP degree present (columns named \code{"power1"},
#' \code{"power2"}, ...). For variables whose selected FP degree is lower than
#' the matrix's number of columns, the extra entries are \code{NA}.
#' @keywords internal
#' @noRd
convert_powers_list_to_matrix <- function(power_list) {
  # Step 1: Determine how many power columns are needed, i.e. the largest
  # number of selected powers across all variables (FP2 has 2, FP1 has 1).
  psize <- sapply(power_list, length)
  maxp <- max(psize)
  
  # Step 2: For each power slot i = 1..maxp, extract the i-th power of every
  # variable. Indexing past a shorter vector's length yields NA automatically,
  # which is exactly the desired padding for lower-degree variables (e.g. the
  # second power of an FP1 variable becomes NA).
  new_list_powers <- vector(mode = "list", length = length(power_list))
  for (i in 1:maxp) {
    new_list_powers[[i]] <- sapply(power_list, function(x) x[i])
  }
  
  # Step 3: Column-bind the per-slot vectors into a matrix and label the
  # columns power1, power2, ... .
  matp           <- do.call(cbind, new_list_powers)
  colnames(matp) <- paste0("power", 1:maxp)
  
  matp
}

#' Helper to create overview table of fp terms
#' 
#' To be used in \code{fit_mfp()}.
#' 
#' @param spike_decision Integer vector indicating the modeling decision for
#' spike-at-zero variables.
#' 
#' @return 
#' Dataframe with overview of all fp terms. Each row represents a variable, 
#' with rownames giving the name of the variable. Variables with acd 
#' transformation are prefixed by `A_` by the `print` and `summary` methods. 
#' The dataframe comprises the following columns: 
#' 
#' * `df_initial`: initial degrees of freedom.
#' * `select`: significance level used for backward elimination (or criterion name if not "pvalue").
#' * `alpha`: significance level for FP terms (or criterion name if not "pvalue").
#' * `acd`: logical, whether an ACD transformation was applied.
#' * `zero`: logical, indicates whether only the positive values of the variable
#'  are transformed (i.e., whether the FP function is applied exclusively to 
#'  values greater than zero).
#' * `catzero`: logical, whether a binary variable for zero values was created.
#' * `spike`: logical, indicates presence of a spike-at-zero variable.
#' * `spike_decision`: integer code describing how the spike-at-zero variable is modeled.
#' * `selected`: logical, whether the FP term is included in the final model.
#' * `df_final`: final estimated degrees of freedom for the variable.
#' * `power1, power2, ...`: final estimated FP powers (as many columns as needed).
#' 
#' @inheritParams fit_mfp
#' @param fp_powers Named list of final selected FP powers, one numeric vector
#' per variable (as stored in \code{fit_mfp()}'s \code{powers_current}). An
#' all-\code{NA} vector means the variable was not selected into the final
#' model. Used both to derive the \code{selected}/\code{df_final} columns and,
#' via \code{convert_powers_list_to_matrix()}, the \code{power1, power2, ...}
#' columns.
#' @param term_to_columns Optional complete conceptual-term-to-raw-column
#'   mapping. When supplied, explicitly mapped retained terms report their
#'   fitted rank contribution rather than one df per stored linear power.
#' @param transformed_to_model_columns Optional named mapping from final
#'   transformed design columns to exact fitted coefficient names.
#' @param coefficients Optional named coefficient vector from the final fitted
#'   model. Missing coefficients are treated as non-estimable and do not
#'   contribute to a grouped term's final df.
#' @keywords internal
#' @noRd
create_fp_terms <- function(fp_powers, 
                            acdx, 
                            df,
                            select, 
                            alpha, 
                            criterion,
                            zero,
                            catzero,
                            spike, 
                            spike_decision,
                            term_to_columns = NULL,
                            transformed_to_model_columns = NULL,
                            coefficients = NULL) {
  
  # Step 1: Align every input vector/list to the same variable order/subset
  # as fp_powers, since callers may pass vectors covering a different (e.g.
  # unordered) set of variable names.
  vars <- names(fp_powers)
  
  acdx <- acdx[vars]
  df <- df[vars]
  select <- select[vars]
  alpha <- alpha[vars]
  zero <- zero[vars]
  catzero <- catzero[vars]
  spike <- spike[vars]
  spike_decision <- spike_decision[vars]
  
  # Step 2: Calculate final degrees of freedom. Continuous terms use the MFP
  # convention implemented by calculate_df(). Explicitly mapped fixed-linear
  # terms are different: one conceptual term can represent several fitted
  # design columns, so its final df is the number of estimable coefficients in
  # that block rather than the single power value stored for selection.
  df_final <- mapply(
    calculate_df,
    fp_powers,
    spike_decision,
    catzero,
    SIMPLIFY = TRUE
  )
  names(df_final) <- vars
  
  if (!is.null(term_to_columns)) {
    missing_terms <- setdiff(vars, names(term_to_columns))
    if (length(missing_terms) > 0L) {
      stop(
        sprintf(
          "Internal error: term-to-column metadata is missing term(s): %s.",
          paste(missing_terms, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    
    mapped <- mapped_term_flags(term_to_columns[vars])
    mapped_selected <- mapped & !vapply(
      fp_powers,
      function(p) all(is.na(p)),
      logical(1L)
    )
    
    if (any(mapped_selected)) {
      if (is.null(transformed_to_model_columns) || is.null(coefficients)) {
        stop(
          paste0(
            "Internal error: fitted column and coefficient metadata are ",
            "required to calculate grouped-term degrees of freedom."
          ),
          call. = FALSE
        )
      }
      
      coefficient_names <- names(coefficients)
      if (is.null(coefficient_names)) {
        stop(
          "Internal error: final fitted coefficients must have names.",
          call. = FALSE
        )
      }
      
      for (term in vars[mapped_selected]) {
        raw_columns <- term_to_columns[[term]]
        transformed_columns <- paste0(raw_columns, ".1")
        
        missing_transformed <- setdiff(
          transformed_columns,
          names(transformed_to_model_columns)
        )
        if (length(missing_transformed) > 0L) {
          stop(
            sprintf(
              paste0(
                "Internal error: transformed-column metadata for grouped ",
                "term '%s' is missing: %s."
              ),
              term,
              paste(missing_transformed, collapse = ", ")
            ),
            call. = FALSE
          )
        }
        
        model_columns <- unname(
          transformed_to_model_columns[transformed_columns]
        )
        missing_coefficients <- setdiff(model_columns, coefficient_names)
        if (length(missing_coefficients) > 0L) {
          stop(
            sprintf(
              paste0(
                "Internal error: coefficient metadata for grouped term '%s' ",
                "is missing: %s."
              ),
              term,
              paste(missing_coefficients, collapse = ", ")
            ),
            call. = FALSE
          )
        }
        
        df_final[[term]] <- sum(!is.na(coefficients[model_columns]))
      }
    }
  }
  
  # Step 3: Assemble the one-row-per-variable overview table.
  fp_terms <- data.frame(
    # initial degrees of freedom
    df_initial = df, 
    select = select, 
    alpha = alpha, 
    acd = acdx, 
    zero = zero,
    catzero = catzero,
    spike = spike,
    # Spike decision
    spike_dec = spike_decision,
    
    # A variable is "selected" if it has a non-NA continuous power, or if it
    # is a binary-only spike variable (spike_decision = 3): such variables
    # have fp_powers all NA by design (see calculate_df()), yet still
    # contribute a term (the binary indicator) to the final model.
    selected = mapply(function(p, sd) {
      if (sd == saz_decision_codes[["binary_only"]]) TRUE  # binary-only spike is still selected
      else !all(is.na(p))
    }, fp_powers, spike_decision),
    df_final = unname(df_final),
    # Adds power1, power2, ... columns (NA-padded to the largest selected FP
    # degree across all variables).
    convert_powers_list_to_matrix(fp_powers)
  )
  
  rownames(fp_terms) <- names(fp_powers)
  
  # Step 4: For non-p-value criteria, `select`/`alpha` no longer represent
  # significance levels, so replace them with the criterion name (e.g. "AIC")
  # for a clearer summary/print display.
  if (criterion != "pvalue") {
    fp_terms$select <- toupper(criterion)
    fp_terms$alpha <- toupper(criterion)
  }
  
  fp_terms
}


#' Backscale Columns of a Matrix (Internal)
#'
#' Multiplies each column of a numeric matrix by a corresponding scalar value 
#' from a named vector. Typically used to reverse prior scaling (i.e., backscaling).
#' This is an internal helper function and not intended for direct use by package 
#' users.
#'
#' @param x A numeric matrix with column names, or `NULL`.
#' @param scalex A named numeric vector. Each name must match a column name of `x`.
#'
#' @return A matrix with backscaled columns, or `NULL` if `x` is `NULL`.
#' @keywords internal
#' @noRd
backscale_matrix <- function(x, scalex) {
  # Step 1: Validate inputs (NULL passthrough, matrix/numeric/names checks).
  # If x is NULL, return NULL
  if (is.null(x)) {
    return(NULL)
  }
  
  # Check: x must be a matrix
  if (!is.matrix(x)) {
    stop("`x` must be a matrix.")
  }
  
  # Check: x must be numeric
  if (!is.numeric(x)) {
    stop("`x` must be a numeric matrix.")
  }
  
  # Check: column names must be present
  vnames <- colnames(x)
  if (is.null(vnames)) {
    stop("`x` must have column names.")
  }
  
  # Check: scalex must be a named numeric vector
  if (!is.numeric(scalex) || is.null(names(scalex))) {
    stop("`scalex` must be a named numeric vector.")
  }
  
  # Check: all columns in x must have matching names in scalex
  missing_cols <- setdiff(vnames, names(scalex))
  if (length(missing_cols) > 0) {
    stop("Missing scaling values for column(s): ", paste(missing_cols, collapse = ", "))
  }
  
  # Step 2: Backscale, i.e. undo `x / scale` (applied upstream before FP power
  # selection) by multiplying each column back by its scale factor.
  unscale_x <- sweep(x, 2, scalex[vnames], FUN = "*")
  return(unscale_x)
}