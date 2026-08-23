#' Validate MFP Candidate Powers Against Closed-Test Requirements
#'
#' Checks that each variable with \code{df > 1} has at least one candidate
#' fractional-polynomial power other than \code{1}, unless that variable has
#' been explicitly exempted from this closed-test requirement. In the ordinary
#' MFP selection engine, power \code{1} is fitted separately as the linear model
#' and is therefore excluded from the degree-1 FP candidate search. The
#' closed-test procedure consequently requires a non-linear FP1 candidate when
#' \code{df > 1}.
#'
#' The exemption is used by MFPI's internally forced FP1 power search. MFPI
#' prespecifies the FP degree and searches the conventional FP1 class, in which
#' the linear power \eqn{p = 1} remains a legitimate candidate. This exemption
#' affects only validation; it does not alter \code{df}, add power \code{1} to
#' a candidate set, or change ordinary \code{mfp2()} closed-test behaviour.
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
#' @param allow_linear_only_fp1 Optional named logical vector, aligned with
#'   \code{df}. \code{TRUE} exempts the corresponding variable from the
#'   closed-test error when its candidate-power set contains only power
#'   \code{1}. This is an internal exception for MFPI's forced FP1 search; the
#'   default is \code{FALSE} for every variable.
#'
#' @return Invisibly returns \code{TRUE} if all candidate-power sets are
#'   compatible with the requested selection procedure.
#'
#' @details
#' This is an internal selection-engine validation. It should be called after
#' \code{powers} and \code{df} have been normalized to named per-variable
#' objects and before calls to \code{select_ra2()}, \code{select_ic()}, or
#' \code{find_best_fpm_step()}.
#'
#' @keywords internal
#' @noRd
validate_mfp_candidate_powers <- function(powers,
                                          df,
                                          allow_linear_only_fp1 = NULL) {
  if (is.null(allow_linear_only_fp1)) {
    allow_linear_only_fp1 <- stats::setNames(
      rep(FALSE, length(df)),
      names(df)
    )
  } else {
    allow_linear_only_fp1 <- stats::setNames(
      as.logical(allow_linear_only_fp1),
      names(allow_linear_only_fp1)
    )[names(df)]
    allow_linear_only_fp1[is.na(allow_linear_only_fp1)] <- FALSE
  }

  # Flag variables whose candidate-power set collapses to just p = 1 after
  # cleaning. In ordinary MFP that leaves no non-linear FP1 candidate for the
  # closed test against the separately fitted linear model.
  only_linear_candidate <- vapply(
    powers,
    function(v) {
      v <- sort(unique(as.numeric(v)))
      v <- v[!is.na(v) & is.finite(v)]
      length(v) == 1L && identical(v, 1)
    },
    logical(1L)
  )

  # MFPI may explicitly exempt its forced FP1 focal terms because there the
  # FP degree is prespecified and p = 1 is a valid member of the FP1 search.
  # All non-exempt variables retain the ordinary MFP validation rule.
  vars_invalid <- names(which(
    only_linear_candidate &
      df > 1L &
      !allow_linear_only_fp1
  ))

  if (length(vars_invalid) > 0L) {
    stop(
      paste0(
        "The following variable(s) have `df > 1` but their candidate-power ",
        "set contains only power 1: ",
        paste(vars_invalid, collapse = ", "),
        ". Power 1 is fitted separately as the ordinary linear model and is ",
        "excluded from FP1 candidate fitting in the ordinary MFP closed-test ",
        "procedure. Use `df = 1` for a purely linear effect, or include at ",
        "least one non-1 candidate power."
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Warn About Low Information Relative to Initial MFP Complexity
#'
#' Replaces the former absolute five-observation failure with a diagnostic that
#' accounts for the model actually requested. The denominator is the sum of the
#' resolved initial term degrees of freedom: ordinary FP terms use their final
#' search df, mapped fixed-linear blocks use their number of design columns, and
#' each active catzero/SAZ representation adds its binary indicator. The
#' intercept is intentionally excluded because this is an MFP term-complexity
#' diagnostic rather than an algebraic rank calculation.
#'
#' @param x Fitting design matrix after row subsetting.
#' @param y Validated fitting response.
#' @param weights Positive fitting weights.
#' @param family_string Normalized family name.
#' @param df Named, resolved initial df vector by conceptual term.
#' @param term_to_columns Complete conceptual-term-to-design-column mapping.
#' @param catzero Named logical vector after the spike-to-catzero cascade.
#' @param threshold Maximum information units per initial df that triggers the
#'   warning.
#'
#' @return Invisibly, a list containing the information count, initial df, and
#'   their ratio. The return value is useful for focused internal tests.
#'
#' @keywords internal
#' @noRd
warn_mfp_information_ratio <- function(x,
                                       y,
                                       weights,
                                       family_string,
                                       df,
                                       term_to_columns,
                                       catzero,
                                       threshold = 5) {
  initial_df <- df

  # A grouped/factor term is fixed linear in the MFP search (`df = 1`) but may
  # occupy several estimable design columns. Count those columns so the warning
  # reflects the complete starting block rather than its search setting.
  mapped <- mapped_term_flags(term_to_columns)
  if (any(mapped)) {
    initial_df[mapped] <- lengths(term_to_columns[mapped])
  }

  # catzero already includes retained spike terms because fit_mfp() applies the
  # spike -> catzero cascade before calling this helper. Its binary indicator is
  # an additional fitted component and therefore contributes one initial df.
  indicator_terms <- names(catzero)[catzero]
  if (length(indicator_terms) > 0L) {
    initial_df[indicator_terms] <- initial_df[indicator_terms] + 1L
  }

  model_df <- sum(initial_df)

  # Use the information unit most closely tied to estimation for each family.
  # Cox models use observed events. For binomial models, the smaller outcome
  # total is used; this also supports grouped success/failure responses and
  # prior trial weights. Other supported families use fitted observations.
  if (identical(family_string, "cox")) {
    status <- y[, ncol(y)]
    information <- sum(status > 0)
    unit_label <- "events"
  } else if (identical(family_string, "binomial")) {
    if (is.matrix(y)) {
      successes <- sum(weights * y[, 1L])
      failures <- sum(weights * y[, 2L])
    } else {
      response <- if (is.factor(y)) {
        as.integer(y) - 1L
      } else {
        as.numeric(y)
      }
      successes <- sum(weights * response)
      failures <- sum(weights * (1 - response))
    }
    information <- min(successes, failures)
    unit_label <- "minority outcome units"
  } else {
    information <- nrow(x)
    unit_label <- "observations"
  }

  ratio <- information / model_df

  if (is.finite(ratio) && ratio <= threshold) {
    warning(
      sprintf(
        paste0(
          "The model has a low information-to-complexity ratio ",
          "(%.2f %s per initial model degree of freedom; ",
          "%.2f information units / %d initial df).\n",
          "Estimates and fractional-polynomial selection may be unstable."
        ),
        ratio,
        unit_label,
        information,
        as.integer(model_df)
      ),
      call. = FALSE
    )
  }

  invisible(list(
    information = information,
    initial_df = model_df,
    ratio = ratio
  ))
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
#'   terms may already be expanded into dummy columns. Ordinary FP variables are
#'   assumed to be shifted and scaled; active ACD variables are shifted but use
#'   scale 1.
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
#'   term. For ordinary FP terms it is used to backscale the corresponding raw
#'   column or column block before the final fit. Active ACD terms use scale 1.
#' @param shift Named numeric vector with one shift per conceptual term. Shifting
#'   has already been applied upstream; the values are reordered to conform to
#'   \code{xorder} and stored in the returned object.
#' @param df Named numeric vector with one degrees-of-freedom setting per
#'   conceptual term. Explicitly mapped terms must have \code{df = 1}.
#' @param center Named logical vector with one centering setting per conceptual
#'   term.
#' @param family Either a character string specifying the model family
#'   (e.g., \code{"gaussian"}, \code{"binomial"}, \code{"poisson"},
#'   \code{"negbin"}, \code{"cox"}) or a function that returns a GLM family object such as
#'   \code{stats::gaussian(link = "identity")}. For Cox models, only the
#'   character string \code{"cox"} is allowed.
#' @param family_string A character string representing the selected family,
#'   e.g., \code{"gaussian"}.
#' @param criterion One character value defining the selection criterion:
#'   \code{"pvalue"}, \code{"aic"}, or \code{"bic"}. It is validated and
#'   normalized at the fitting-engine boundary because \code{fit_mfp()} is also
#'   called directly by MFPI internals rather than exclusively through the
#'   public \code{mfp2()} methods.
#' @param select Named numeric vector with one significance threshold per
#'   conceptual term, used during backfitting to decide whether the term is
#'   retained.
#' @param alpha Named numeric vector with one functional-form significance
#'   threshold per conceptual term.
#' @param keep A character vector with names of variables to be kept in the
#'   model regardless of selection criteria.
#' @param xorder A string determining the visiting order of conceptual terms.
#'   For significance-based ordering, terms are ranked by leave-one-term-out
#'   likelihood-ratio tests from the resolved full linear reference model. A
#'   zero term contributes its positive-part column, while a catzero or retained
#'   spike term contributes the positive-part column plus its binary indicator.
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
#'   If \code{TRUE} for a non-linear term, forces selection of the most complex
#'   functional form allowed by \code{df} under p-value, AIC, or BIC selection,
#'   bypassing variable elimination and simpler-form comparisons. For an
#'   eligible spike-at-zero term, the forced result is the complete maximum SAZ
#'   representation (maximum positive-component form plus binary zero indicator),
#'   so the reduced-component comparisons in SAZ Stage 2 are not run.
#' @param retain_linear_fp1 Internal logical specifically for MFPI. MFPI uses
#'   \code{fit_mfp()} to perform forced fixed-degree FP power searches, but its
#'   FP1 search follows the conventional FP1 class and therefore retains the
#'   linear power \eqn{p = 1} as a candidate. Set this to \code{TRUE} only for
#'   those MFPI power-selection fits. The default \code{FALSE} preserves the
#'   ordinary \code{mfp2()} closed-test behaviour, where the linear model is
#'   fitted separately from the non-linear FP1 candidates. This option never
#'   adds power \code{1} when it is absent from the supplied candidate set.
#' @param has_offset logical indicating whether an offset was specified before
#' missing offsets were replaced by zeros internally.
#' @param verbose Logical. If \code{TRUE}, progress information is printed
#'   during fitting. Default \code{FALSE}.
#' @param warn_low_information Internal logical. If \code{TRUE}, emit one
#'   information-to-complexity warning after the effective initial df have been
#'   resolved. Public \code{mfp2()} calls enable it; auxiliary internal fits use
#'   the default \code{FALSE}.
#'
#' @section Algorithm:
#' \enumerate{
#'   \item \strong{Pre-processing and structural-zero resolution.} ACD
#'     eligibility is resolved on the unre-coded design. The
#'     spike-to-catzero-to-zero hierarchy and any defensive SAZ reset are then
#'     applied, followed by SAZ df capping.
#'   \item \strong{Reference representation.} Final zero terms are recoded to
#'     their positive part and catzero/retained spike indicators are built once.
#'     These same objects are reused by the reference fit and backfitting.
#'   \item \strong{Full linear reference and variable ordering.} The full
#'     reference contains each conceptual term in its resolved linear starting
#'     form: ordinary x, positive-part x for zero terms, and positive-part x plus
#'     the structural-zero indicator for catzero/retained spike terms. Its fit
#'     statistics are retained as \code{null_deviance},
#'     \code{linear_deviance}, \code{linear_logl}, and \code{linear_df}. For
#'     ascending/descending ordering, each conceptual term is removed as one
#'     block and compared with the full reference by a likelihood-ratio test.
#'     With \code{xorder = "original"} or one term, no reduced ordering fits
#'     are required.
#'   \item \strong{MFP backfitting cycles.} FP powers and spike decisions are
#'     updated iteratively via \code{find_best_fp_cycle()} until convergence or
#'     the maximum number of cycles is reached. Previously computed adjustment
#'     transformations are cached in \code{prev_adj_params} to avoid
#'     recomputation across cycles.
#'   \item \strong{Final transformation.} After convergence, \code{x} is
#'     backscaled (if \code{scale != 1}) to restore the shifted-but-not-scaled
#'     variable, then FP-transformed using the selected powers. Centering is
#'     applied and the already-built catzero indicators are retained or removed
#'     according to the final SAZ decision.
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
#'   \item If \code{saz_pre_resolved = FALSE}, \code{reset_spike()} is called as
#'     a defensive fallback for internal callers that have not resolved SAZ
#'     eligibility before preprocessing.
#'   \item Spike-specific proportions and positive-part df capping are evaluated
#'     directly with \code{x <= 0} / \code{x > 0} predicates before physical
#'     zero recoding.
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
                    retain_linear_fp1 = FALSE,
                    has_offset,
                    verbose,
                    term_to_columns = NULL,
                    fitter = "base",
                    warn_low_information = FALSE) {

  # fit_mfp() is the boundary of the reusable fitting engine. Public mfp2()
  # methods already use match.arg(), but MFPI and future internal callers can
  # invoke this function directly. Validate before touching any other required
  # argument so an invalid criterion cannot fall through a downstream switch()
  # to NULL and cause an unrelated indexing failure.
  if (!is.character(criterion) ||
      length(criterion) != 1L ||
      is.na(criterion) ||
      !nzchar(criterion)) {
    stop(
      "`criterion` must be one character value: 'pvalue', 'aic', or 'bic'.",
      call. = FALSE
    )
  }

  criterion <- tolower(criterion)

  if (!criterion %in% c("pvalue", "aic", "bic")) {
    stop(
      "`criterion` must be one of 'pvalue', 'aic', or 'bic'.",
      call. = FALSE
    )
  }

  # fit_mfp() implements the full MFP/MFPA/SAZ algorithm described in the
  # "Algorithm" section above: resolve ACD/zero/catzero/spike settings, build
  # the linear reference and visiting order, run backfitting cycles, then transform
  # and fit the final model. mfp2.default()/mfp2.formula() are expected to have
  # already validated and preprocessed all inputs: ordinary FP variables are
  # shifted/scaled, whereas active ACD variables are shifted with scale 1.

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

  # Step 1: Resolve the family object -----------------------------------------
  # Resolve GLM family objects once for all repeated internal model fits.
  # Public mfp2.default() already does this, but keeping it here makes direct
  # internal calls to fit_mfp() avoid repeated stats::gaussian()/binomial()/
  # poisson() construction as well. Cox remains the character string "cox".
  family_fit <- resolve_fit_model_family(family)

  # Initial df are reported after the reference model and visiting order have
  # been resolved. Before that point, keep all per-variable settings aligned to
  # `variables_x`, the conceptual-term order supplied to fit_mfp().

  # Step 2: Align settings and configure ACD on the unre-coded design ---------
  # Ordering is intentionally deferred until the structural-zero representation
  # is final. ACD eligibility must be resolved first because reset_acd() inspects
  # the unre-coded covariate values. Naming these vectors here is alignment only;
  # it does not change their conceptual-term order.
  df           <- stats::setNames(df,      variables_x)
  scale        <- if (!is.null(names(scale))) {
    scale[variables_x]
  } else {
    stats::setNames(scale, variables_x)
  }
  acdx         <- stats::setNames(acdx,    variables_x)
  zero         <- stats::setNames(zero,    variables_x)
  catzero      <- stats::setNames(catzero, variables_x)
  spike        <- stats::setNames(spike,   variables_x)
  force_max_fp <- stats::setNames(force_max_fp, variables_x)
  powers       <- powers[variables_x]

  variables_acd <- character(0L)

  # Requesting ACD forces df = 4 because FSPA evaluates the full FP1(p1, p2)
  # model space. Keep x unre-coded here so reset_acd() sees the same values used
  # to determine ACD eligibility.
  if (any(acdx)) {
    acdx          <- reset_acd(x, acdx)
    variables_acd <- names(acdx)[acdx]

    # Public mfp2() callers set active ACD variables to scale 1 before entering
    # fit_mfp(). Keep a defensive path for direct internal callers that supplied
    # a scaled working matrix: restore shifted, unscaled values before any
    # FP/ACD candidate or linear reference term is generated.
    # ACD should not be scaled but can be shifted, see mfpa paper
    acd_scaled <- variables_acd[scale[variables_acd] != 1]
    if (length(acd_scaled) > 0L) {
      x[, acd_scaled] <- sweep(
        x[, acd_scaled, drop = FALSE],
        2,
        scale[acd_scaled],
        "*"
      )
      scale[acd_scaled] <- 1
    }

    df[variables_acd] <- 4
  }

  # Step 3: Resolve spike/catzero/zero hierarchy and SAZ eligibility ----------
  # Preserve explicit user intent before applying the hierarchy. reset_spike()
  # uses these copies to restore zero/catzero settings when a requested spike
  # term is ineligible.
  user_catzero <- catzero
  user_zero    <- zero

  catzero[spike] <- TRUE    # spike implies catzero
  zero[catzero]  <- TRUE    # catzero implies zero

  # Public callers normally resolve SAZ eligibility before shift/scale
  # preprocessing. Keep the defensive reset for direct internal callers.
  if (!isTRUE(saz_pre_resolved) && any(spike)) {
    result <- reset_spike(
      x                      = x,
      spike                  = spike,
      user_catzero           = user_catzero,
      user_zero              = user_zero,
      min_saz_component_prop = min_saz_component_prop
    )

    spike   <- result$spike
    catzero <- result$catzero
    zero    <- result$zero
  }

  # Retain the structural-zero share before x is physically recoded. This is
  # descriptive SAZ metadata, not a model-selection result.
  prop_zero <- calculate_saz_prop_zero(
    x = x,
    spike = spike,
    term_to_columns = term_to_columns
  )

  # For retained spike-at-zero variables, cap the maximum FP df using only the
  # positive component. The positive/nonpositive split must also be evaluated
  # before physical recoding.
  df <- cap_spike_df(
    x     = x,
    df    = df,
    spike = spike
  )

  # The public mfp2() interfaces request this diagnostic only after all rules
  # that determine the starting model complexity have run. In particular, ACD
  # has forced df = 4, SAZ eligibility has been resolved, and the positive-part
  # cardinality cap has been applied. This is therefore more informative than
  # either a fixed minimum row count or n / ncol(x).
  if (isTRUE(warn_low_information)) {
    warn_mfp_information_ratio(
      x               = x,
      y               = y,
      weights         = weights,
      family_string   = family_string,
      df               = df,
      term_to_columns  = term_to_columns,
      catzero          = catzero,
      threshold        = 5
    )
  }

  # Validate candidate powers after ACD forcing and SAZ-specific df capping.
  # MFPI may deliberately retain p = 1 as an FP1 candidate when that degree is
  # prespecified; ordinary mfp2() still requires a non-linear FP1 candidate when
  # df > 1 because its linear model is fitted separately.
  allow_linear_only_fp1 <- stats::setNames(
    rep(FALSE, length(df)),
    names(df)
  )
  if (isTRUE(retain_linear_fp1)) {
    forced_fp1 <- force_max_fp & df == 2L
    allow_linear_only_fp1[forced_fp1] <- TRUE
  }

  validate_mfp_candidate_powers(
    powers = powers,
    df = df,
    allow_linear_only_fp1 = allow_linear_only_fp1
  )

  # Step 4: Recode zero components in x once ----------------------------------
  # Apply the final zero flags after cascade/reset. The same recoded x is used
  # by the full linear reference, variable ordering, and all backfitting cycles;
  # no ordering-specific copy is created.
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

  # zero_x remains TRUE for zero-component terms even though x has already been
  # recoded. The cycle-level transformation code still needs this flag to leave
  # structural zeros untouched when FP transformations are evaluated.

  # Step 5: Build catzero indicators once -------------------------------------
  # Build the binary structural-zero columns after zero recoding, then reuse the
  # same matrices in the full linear reference and in backfitting. This avoids
  # calculating the indicator twice. Each active catzero term is a conceptual
  # two-component block: positive-part x plus I(original x <= 0).
  catzero_mat_list <- stats::setNames(
    vector("list", length(catzero)),
    names(catzero)
  )

  catzero_terms <- names(catzero)[catzero]
  if (length(catzero_terms) > 0L) {
    for (v in catzero_terms) {
      columns <- term_to_columns[[v]]
      if (length(columns) != 1L) {
        stop(
          sprintf(
            "Internal error: catzero term '%s' must map to one design column.",
            v
          ),
          call. = FALSE
        )
      }

      column <- columns[[1L]]
      catzero_mat_list[[v]] <- matrix(
        as.integer(x[, column] <= 0),
        ncol = 1L,
        dimnames = list(rownames(x), "catzero")
      )
    }
  }

  # Step 6: Fit the full linear reference and determine visiting order --------
  # The reference model now matches the resolved starting representation:
  #   ordinary term        -> x
  #   zero term            -> x+ (already recoded in Step 4)
  #   catzero/spike term   -> x+ + I(x <= 0)
  # catzero indicators are passed as already-built blocks; order_variables()
  # only assembles the fit matrix and never recreates zero/catzero data.
  #
  # The full reference is fitted for every xorder because its fit statistics are
  # retained in the mfp2 object. Reduced leave-one-term-out fits are needed only
  # for ascending/descending significance ordering with more than one term.
  ordering_result <- order_variables(
    xorder          = xorder,
    x               = x,
    term_to_columns = term_to_columns,
    catzero_blocks  = if (length(catzero_terms) > 0L) {
      catzero_mat_list
    } else {
      NULL
    },
    y               = y,
    family          = family_fit,
    family_string   = family_string,
    fitter          = fitter,
    weights         = weights,
    offset          = offset,
    strata          = strata,
    method          = method,
    control         = control,
    nocenter        = nocenter
  )

  variables_ordered <- ordering_result$variables_ordered
  null_deviance     <- ordering_result$null_deviance
  linear_deviance   <- ordering_result$linear_deviance
  linear_logl       <- ordering_result$linear_logl
  linear_df         <- ordering_result$linear_df
  null_logl         <- ordering_result$null_logl

  # Step 7: Initialize powers and align all working objects to visiting order --
  # Only now is the visiting order known. Reorder every per-variable object and
  # the raw design columns together so subsequent code can index consistently.
  powers_current <- stats::setNames(
    as.list(rep(1, length(variables_ordered))),
    variables_ordered
  )

  alpha        <- stats::setNames(alpha,  variables_x)[variables_ordered]
  select       <- stats::setNames(select, variables_x)[variables_ordered]
  df           <- df[variables_ordered]
  center       <- stats::setNames(center, variables_x)[variables_ordered]
  shift        <- if (!is.null(names(shift))) {
    shift[variables_ordered]
  } else {
    stats::setNames(shift, variables_x)[variables_ordered]
  }
  scale        <- scale[variables_ordered]
  acdx         <- acdx[variables_ordered]
  zero         <- zero[variables_ordered]
  zero_x       <- zero_x[variables_ordered]
  catzero      <- catzero[variables_ordered]
  spike        <- spike[variables_ordered]
  force_max_fp <- force_max_fp[variables_ordered]
  powers       <- powers[variables_ordered]
  prop_zero    <- prop_zero[variables_ordered]
  catzero_mat_list <- catzero_mat_list[variables_ordered]

  term_to_columns <- term_to_columns[variables_ordered]
  raw_columns_ordered <- unlist(term_to_columns, use.names = FALSE)
  if (!identical(colnames(x), raw_columns_ordered)) {
    x <- x[, raw_columns_ordered, drop = FALSE]
  }
  if (!identical(colnames(x), raw_columns_ordered)) {
    stop(
      "Internal error: x column order does not match term_to_columns.",
      call. = FALSE
    )
  }

  # Active ACD terms start as a linear x component with no ACD component yet.
  # Their df was already forced to 4 in Step 2.
  if (any(acdx)) {
    variables_acd <- names(acdx)[acdx]
    powers_current <- utils::modifyList(
      powers_current,
      sapply(variables_acd, function(v) c(1, NA), simplify = FALSE)
    )
  }

  # `keep` forces variables into the final model regardless of their
  # backfitting significance under the p-value criterion.
  if (!is.null(keep)) {
    select[which(names(select) %in% keep)] <- 1
  }

  if (isTRUE(verbose)) {
    # `df` is the effective continuous FP/ACD search setting after all early
    # capping. For display, convert it to the df represented when selection
    # starts. Mapped fixed terms contribute one df per raw design column, and a
    # retained SAZ term contributes one additional structural-zero indicator df.
    df_initial_display <- df

    mapped_display <- mapped_term_flags(term_to_columns)
    if (any(mapped_display)) {
      df_initial_display[mapped_display] <- lengths(
        term_to_columns[mapped_display]
      )
    }

    if (any(spike)) {
      df_initial_display[spike] <- df_initial_display[spike] + 1L
    }

    message(sprintf(
      "Visiting order: %s",
      paste0(variables_ordered, collapse = ", ")
    ))

    df_text <- utils::capture.output(
      print(
        matrix(
          df_initial_display,
          nrow = 1,
          dimnames = list("df", variables_ordered)
        ),
        quote = FALSE
      )
    )

    message(
      "Initial degrees of freedom:\n",
      paste(df_text, collapse = "\n")
    )
  }

  # Spike decision initialisation. continuous_only means standard FP handling
  # unless the SAZ procedure later retains a structural-zero component.
  spike_decision <- stats::setNames(
    rep(saz_decision_codes[["continuous_only"]], length(variables_ordered)),
    variables_ordered
  )

  # Step 8: Cache ACD parameters and selection sample size --------------------
  # ACD parameters are estimated once per active variable and reused throughout
  # the backfitting cycles.
  acd_parameter <- lapply(names(acdx), function(v) {
    if (isTRUE(acdx[[v]])) {
      fit_acd(
        x      = x[, v],
        powers = powers[[v]],
        zero   = zero[[v]],
        fitter = fitter
      )
    } else {
      NULL
    }
  })
  names(acd_parameter) <- names(acdx)

  # Number of events (Cox) or observations (other families), used for AIC/BIC
  # selection and to guard against a Cox fit with no observed events.
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

  # Step 9: Run MFP backfitting cycles until convergence ----------------------
  # A cycle is one complete pass through all variables (see find_best_fp_cycle()
  # documentation); convergence means neither the selected powers nor the SAZ
  # stage-2 decisions changed compared to the previous cycle.
  j         <- 1L
  converged <- FALSE

  prev_adj_params        <- vector("list", length = length(variables_ordered))
  names(prev_adj_params) <- variables_ordered

  # Cache individual transformed adjustment-variable blocks across focal
  # variables. This is separate from prev_adj_params: the latter is keyed by
  # focal variable and retains only the metadata/per-variable blocks needed as
  # a fallback on the next cycle. It deliberately does NOT retain the complete
  # assembled data_adj matrix, because doing so for every focal variable grows
  # retained memory approximately as O(n * p^2). transform_cache is keyed by
  # the variable being transformed and remains the primary reuse mechanism.
  # The ordinary list is scoped to this fit_mfp() call and threaded explicitly
  # through the cycle, avoiding global state and environment side effects.
  transform_cache <- setNames(
    vector("list", length(variables_ordered)),
    variables_ordered
  )

  while (!converged && j <= cycles) {
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
      fitter          = fitter,
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
      transform_cache = transform_cache,
      force_max_fp    = force_max_fp,
      retain_linear_fp1 = retain_linear_fp1,
      has_offset      = has_offset,
      n_obs           = n_obs,
      verbose         = verbose
    )

    powers_updated         <- fit_best_cycle$powers_current
    spike_decision_updated <- fit_best_cycle$spike_decision
    prev_adj_params        <- fit_best_cycle$prev_adj_params
    transform_cache        <- fit_best_cycle$transform_cache

    # Compare powers after normalizing away the stored-but-inactive continuous
    # power of binary-only spike variables (spike_decision = 3), since that
    # power does not affect the fitted adjustment matrix and would otherwise
    # cause spurious non-convergence.
    powers_same <- identical(
      normalize_powers_for_convergence(
        powers_current,
        spike_decision
      ),
      normalize_powers_for_convergence(
        powers_updated,
        spike_decision_updated
      )
    )

    spike_same <- identical(
      spike_decision,
      spike_decision_updated
    )

    # Store the results from the completed cycle before deciding whether another
    # cycle is required. This also ensures that the final cycle's results are
    # retained when the maximum number of cycles is reached.
    powers_current <- powers_updated
    spike_decision <- spike_decision_updated

    converged <- powers_same && spike_same

    if (converged) {
      if (verbose) {
        message(sprintf(
          "Fractional polynomial fitting algorithm converged after %d cycle(s).",
          j
        ))
      }
    } else {
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

  # Step 10: Apply the final FP/ACD transformation and centering --------------
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

  # ACD parameters are estimated on shifted, unscaled predictors. New fits
  # therefore store scale = 1. The field is retained for compatibility with the
  # existing transformation and prediction helpers and with legacy model objects.
  acd_parameter_final <- acd_parameter

  for (v in names(acd_parameter_final)) {
    if (!is.null(acd_parameter_final[[v]])) {
      # Remove training-data ACD values. fit_acd() returns $acd as the
      # transformed training vector, but only beta0/beta1/power/shift/scale
      # are needed for apply_acd() during prediction. Keeping $acd wastes
      # memory and can cause confusion.
      acd_parameter_final[[v]]$acd <- NULL
      acd_parameter_final[[v]]$shift <- 0
      acd_parameter_final[[v]]$scale <- 1
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

  # Step 11: Fit the final model and assemble the returned mfp2 object --------
  modelfit <- fit_model(
    x             = data_transformed$x_transformed,
    y             = y,
    family        = family_fit,
    family_string = family_string,
    fitter        = fitter,
    weights       = weights,
    offset        = offset,
    method        = method,
    strata        = strata,
    control       = control,
    rownames      = rownames(data_transformed$x_transformed),
    nocenter      = nocenter,
    fast          = FALSE,
    calculate_fit_statistics = TRUE,
    keep_fit      = TRUE,
    has_offset    = has_offset,
    # Reserve original predictor names in addition to transformed final-model
    # columns so package-created response/offset/strata helpers cannot reuse a
    # user-facing name such as y, offset_, or strata_.
    reserved_names = colnames(x)
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

  # The full linear reference likelihood was computed once in Step 6 by order_variables().
  # Cox also supplies its null partial likelihood from that same fit.
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
      # Likelihood-scale quantities retained for inference and compatibility.
      # GLM Model Fit reporting uses the deviance fields above; Cox reporting
      # uses the equivalent -2 partial-log-likelihood deviances.
      null_logl       = null_logl,
      linear_logl     = linear_logl,
      linear_df       = linear_df,
      mfp_logl        = mfp_logl,
      mfp_df          = mfp_df,
      family_string   = family_string,
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
        prop_zero = prop_zero,
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
      transformed_column_to_source = data_transformed$transformed_column_to_source,
      transformed_column_component = data_transformed$transformed_column_component,
      transformed_column_zero_handled =
        data_transformed$transformed_column_zero_handled,
      transformed_column_centered = data_transformed$transformed_column_centered,
      cox_offset_reference = cox_offset_reference,
      fitter          = fitter
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
#' @param prev_adj_params Named focal-variable cache of previously computed
#' adjustment metadata and per-variable transformed blocks. Complete assembled
#' adjustment matrices are deliberately excluded to limit retained memory. The
#' cache is updated at each step and reused in the next cycle.
#' @param transform_cache Named run-scoped list of individual transformed
#'   variable blocks, shared across focal-variable evaluations and updated
#'   explicitly as the MFP cycle proceeds.
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
                               transform_cache = NULL,
                               force_max_fp,
                               retain_linear_fp1 = FALSE,
                               has_offset,
                               n_obs,
                               term_to_columns = NULL,
                               fitter = "base"
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
      fitter = fitter,
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
      transform_cache = transform_cache,
      force_max_fp = force_max_fp,
      retain_linear_fp1 = retain_linear_fp1,
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
    # Carry the run-scoped per-variable cache forward immediately so later
    # focal variables in this same cycle can reuse unchanged transformations.
    transform_cache      <- fit_best_fp_step$transform_cache
    # Cache xi's adjustment-variable transformations, keyed by xi, so the next
    # cycle can reuse unchanged per-variable blocks instead of recomputing them.
    #
    # Important memory rule: current_adj_params must contain data_adj/data_xi
    # while find_best_fp_step() is running because SAZ stage 2 reuses those
    # matrices. Once the focal-variable step is complete they are disposable.
    # Persisting the complete n x adjustment matrix for every xi makes
    # prev_adj_params scale approximately as O(n * p^2), so strip those two
    # focal-step matrices before carrying the cache into the next cycle.
    prev_adj_params[[xi]] <- compact_prev_adj_cache_entry(
      fit_best_fp_step$current_adj_params[[xi]]
    )
  }

  list(powers_current  = powers_current,
       spike_decision  = spike_decision,
       prev_adj_params = prev_adj_params,
       transform_cache = transform_cache)
}

#' Compact a focal-variable adjustment cache entry
#'
#' The active focal-variable step temporarily needs the assembled adjustment
#' matrix (`data_adj`) and, for SAZ stage 2, the selected focal design (`data_xi`).
#' Neither matrix is needed after that step has finished. Keeping them inside
#' `prev_adj_params` for every focal variable retains an n-by-p-scale matrix p
#' times and can therefore make cache memory grow approximately as O(n * p^2).
#'
#' This helper releases only those disposable whole-step matrices. The smaller
#' metadata and `data_adj_list` blocks are retained so the historical per-focal
#' fallback cache remains available when the shared `transform_cache` does not
#' contain a reusable entry.
#'
#' @param params_xi Adjustment-cache entry for one focal variable.
#' @return `params_xi` with `data_adj` and `data_xi` retained as explicit
#'   named `NULL` entries, so their large matrix payloads are released without
#'   exposing the cache to R's partial `$` matching.
#' @keywords internal
#' @noRd
compact_prev_adj_cache_entry <- function(params_xi) {
  if (is.null(params_xi)) {
    return(NULL)
  }

  # Keep explicit named NULL slots rather than removing the elements.
  #
  # This is important because `$` performs partial matching on lists. If the
  # `data_adj` element were removed while `data_adj_list` remained, an access
  # such as `params_xi$data_adj` could silently resolve to `data_adj_list` and
  # leak a list into code that requires a numeric matrix. Single-bracket list
  # assignment preserves the exact names while releasing the large matrices.
  params_xi["data_adj"] <- list(NULL)
  params_xi["data_xi"] <- list(NULL)
  params_xi
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
  # This helper is called only after per-variable powers have been assembled.
  # An empty list or zero-length power vector indicates an inconsistent
  # internal state; fail explicitly rather than letting 1:maxp create a
  # descending or otherwise unintended sequence.
  if (!is.list(power_list) || length(power_list) == 0L) {
    stop(
      "Internal error: `power_list` must contain at least one variable.",
      call. = FALSE
    )
  }

  # Step 1: Determine how many power columns are needed, i.e. the largest
  # number of selected powers across all variables (FP2 has 2, FP1 has 1).
  psize <- vapply(power_list, length, integer(1L))
  if (any(psize == 0L)) {
    stop(
      "Internal error: every element of `power_list` must contain at least one power.",
      call. = FALSE
    )
  }
  maxp <- max(psize)

  # Step 2: For each power slot i = 1..maxp, extract the i-th power of every
  # variable. `seq_len()` is intentional here: unlike `1:maxp`, it remains
  # empty when its upper bound is zero and therefore cannot create 1, 0.
  # Indexing past a shorter vector's length yields NA automatically, which is
  # exactly the desired padding for lower-degree variables.
  new_list_powers <- vector(mode = "list", length = maxp)
  for (i in seq_len(maxp)) {
    new_list_powers[[i]] <- sapply(power_list, function(x) x[i])
  }

  # Step 3: Column-bind the per-slot vectors into a matrix and label the
  # columns power1, power2, ... .
  matp           <- do.call(cbind, new_list_powers)
  colnames(matp) <- paste0("power", seq_len(maxp))

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
#' * `df_setting`: MFP complexity setting supplied to the selection algorithm.
#'   For explicitly mapped grouped terms, `1` denotes a fixed linear block.
#' * `df_initial`: degrees of freedom represented by the term in the initial
#'   MFP model. For grouped terms this is the number of member design columns,
#'   rather than the fixed-linear `df_setting` value. For an eligible SAZ term,
#'   this also includes the one-df structural-zero binary indicator.
#' * `select`: significance level used for backward elimination (or criterion name if not "pvalue").
#' * `alpha`: significance level for FP terms (or criterion name if not "pvalue").
#' * `acd`: logical, whether an ACD transformation was applied.
#' * `zero`: logical, indicates whether only the positive values of the variable
#'  are transformed (i.e., whether the FP function is applied exclusively to
#'  values greater than zero).
#' * `catzero`: logical, whether a binary variable for zero values was created.
#' * `spike`: logical, indicates presence of a spike-at-zero variable.
#' * `prop_zero`: proportion of finite fitting-sample observations in the
#'   structural-zero component for retained SAZ terms; `NA` otherwise.
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
#' @param prop_zero Optional named numeric vector containing the structural-zero
#'   proportion for each retained SAZ term. Missing or non-SAZ terms are stored
#'   as `NA_real_`.
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
                            prop_zero = NULL,
                            term_to_columns = NULL,
                            transformed_to_model_columns = NULL,
                            coefficients = NULL) {

  # Step 1: Align every input vector/list to the same variable order/subset
  # as fp_powers, since callers may pass vectors covering a different (e.g.
  # unordered) set of variable names.
  vars <- names(fp_powers)

  acdx <- acdx[vars]
  df <- df[vars]
  df_setting <- df
  select <- select[vars]
  alpha <- alpha[vars]
  zero <- zero[vars]
  catzero <- catzero[vars]
  spike <- spike[vars]
  spike_decision <- spike_decision[vars]
  if (is.null(prop_zero)) {
    prop_zero <- stats::setNames(rep(NA_real_, length(vars)), vars)
  } else {
    prop_zero <- prop_zero[vars]
  }
  prop_zero[!spike] <- NA_real_

  # Step 2: Separate the MFP search setting from the initial model degrees of
  # freedom. For ordinary continuous terms these values coincide. Explicitly
  # mapped terms are fixed linear blocks (`df_setting = 1`) but can represent
  # several raw design columns, all of which contribute to `df_initial`. The
  # package validates grouped design rank before fitting, so the member-column
  # count is the initial rank contribution for fitted package objects.
  df_initial <- df_setting
  mapped <- stats::setNames(rep(FALSE, length(vars)), vars)

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
    if (any(mapped)) {
      df_initial[mapped] <- lengths(term_to_columns[vars][mapped])
    }
  }

  # Every eligible SAZ term enters the MFP procedure with the structural-zero
  # binary indicator in addition to its continuous FP/ACD component. `spike`
  # has already been eligibility-resolved upstream, whereas `catzero` here is
  # final/effective metadata and may be FALSE after Stage 2 removes the binary
  # component. Count the indicator from `spike` so `df_initial` describes the
  # model that was actually entered into selection (for example, df = 4 starts
  # as 5 df: FP2 contributes 4 df and the zero indicator contributes 1 df).
  if (any(spike)) {
    df_initial[spike] <- df_initial[spike] + 1L
  }

  # Step 3: Calculate final degrees of freedom. Continuous terms use the MFP
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

  # Step 4: Assemble the one-row-per-variable overview table.
  fp_terms <- data.frame(
    # MFP search setting and initial model degrees of freedom
    df_setting = unname(df_setting),
    df_initial = unname(df_initial),
    select = select,
    alpha = alpha,
    acd = acdx,
    zero = zero,
    catzero = catzero,
    spike = spike,
    prop_zero = unname(prop_zero),
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

  # Step 5: For non-p-value criteria, `select`/`alpha` no longer represent
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
