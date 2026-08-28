# -----------------------------------------------------------------------------
# Spike-at-zero decision helpers ----------------------------------------------
# -----------------------------------------------------------------------------

#' Spike-at-Zero Decision Codes
#'
#' Internal constants for final spike-at-zero (SAZ) decisions. These values are
#' stored in fitted objects and must remain stable.
#'
#' \describe{
#'   \item{\code{cont_binary = 1L}}{
#'     Continuous FP/linear/ACD component plus binary zero-indicator retained.
#'   }
#'   \item{\code{continuous_only = 2L}}{
#'     Continuous FP/linear/ACD component retained; binary zero-indicator dropped.
#'   }
#'   \item{\code{binary_only = 3L}}{
#'     Binary zero-indicator retained; continuous component dropped.
#'   }
#' }
#'
#' @keywords internal
#' @noRd
saz_decision_codes <- c(
  cont_binary     = 1L,
  continuous_only = 2L,
  binary_only     = 3L
)


#' Render Spike-at-Zero Decision Labels
#'
#' Maps SAZ decision codes to print, plot, and verbose-output labels.
#' Plot titles use the same component terminology as printed SAZ decisions
#' while explicitly identifying spike-at-zero terms.
#'
#' @param decision Numeric or integer SAZ decision code.
#' @param style Character scalar. One of \code{"print"}, \code{"plot_title"},
#'   or \code{"stage2_model"}.
#' @param continuous_label Character scalar used by plot and stage-2 labels.
#' @param full_label Character scalar used as the stage-2 full-model label.
#' @param binary_label Character scalar used as the stage-2 binary-only label.
#' @param unknown Character scalar used for invalid codes.
#'
#' @return Character vector for \code{style = "print"}; character scalar for
#'   other styles.
#'
#' @keywords internal
#' @noRd
saz_decision_label <- function(decision,
                               style = c("print", "plot_title", "stage2_model"),
                               continuous_label = "",
                               full_label = NULL,
                               binary_label = "Binary",
                               unknown = "unknown") {
  style <- match.arg(style)
  decision_int <- as.integer(decision)

  # "print" is the only vectorized style: it is used to label a whole column
  # of spike_dec values in fp_terms (one row per variable), so it must accept
  # and return a vector rather than a single decision. The other two styles
  # only ever label a single variable's decision at a time (a plot title, or
  # one row in the stage-2 summary table), so they short-circuit below to a
  # scalar `unknown` if handed anything other than exactly one, non-NA value.
  if (style == "print") {
    out <- rep(unknown, length(decision_int))
    out[decision_int == saz_decision_codes[["cont_binary"]]] <- "continuous + binary"
    out[decision_int == saz_decision_codes[["continuous_only"]]] <- "continuous only"
    out[decision_int == saz_decision_codes[["binary_only"]]] <- "binary only"
    return(out)
  }

  if (length(decision_int) != 1L || is.na(decision_int)) {
    return(unknown)
  }

  # binary_only (decision 3): only the structural-zero indicator survives
  # stage 2; there is no continuous component to name, so continuous_label is
  # never referenced in this branch.
  if (decision_int == saz_decision_codes[["binary_only"]]) {
    if (style == "plot_title") {
      return("spike at zero \u2014 binary only")
    }

    return(binary_label)
  }

  # cont_binary (decision 1): both components survive. Note that `full_label`
  # has no default value of its own kind (default NULL) - callers using
  # style = "stage2_model" are expected to always supply it explicitly (see
  # print_mfp_summary.R), since there is no sensible generic fallback text for
  # "the label of the full two-component model" the way there is for the
  # continuous-only/binary-only cases.
  if (decision_int == saz_decision_codes[["cont_binary"]]) {
    if (style == "plot_title") {
      return(paste0(
        "spike at zero \u2014 continuous + binary: ",
        continuous_label
      ))
    }

    return(full_label)
  }

  # continuous_only (decision 2): only the continuous FP/linear/ACD component
  # survives; the binary zero indicator is dropped, so binary_label is never
  # referenced in this branch.
  if (decision_int == saz_decision_codes[["continuous_only"]]) {
    if (style == "plot_title") {
      return(paste0("spike at zero \u2014 continuous only: ", continuous_label))
    }

    return(continuous_label)
  }

  # decision_int matched none of the three known codes (e.g. corrupted data);
  # fall back to the generic placeholder rather than erroring, since this
  # function is used in display/labeling contexts, not validation contexts.
  unknown
}

#' Calculate Spike-at-Zero Structural-Zero Proportions
#'
#' Internal helper that reports, for each retained spike-at-zero term, the
#' proportion of finite observations belonging to the structural-zero
#' component. Exact-zero values are counted directly as structural zeros;
#' callers reject negative values before this helper is reached.
#' Non-SAZ terms receive
#' `NA_real_`.
#'
#' @param x Numeric design matrix on the actual fitting sample.
#' @param spike Named logical vector indexed by conceptual term.
#' @param term_to_columns Optional complete conceptual-term-to-raw-column
#'   mapping. Active SAZ terms must map to exactly one raw column.
#'
#' @return Named numeric vector aligned with `spike`.
#' @keywords internal
#' @noRd
calculate_saz_prop_zero <- function(x, spike, term_to_columns = NULL) {
  if (is.null(names(spike))) {
    stop("Internal error: `spike` must be a named logical vector.", call. = FALSE)
  }

  if (is.null(term_to_columns)) {
    term_to_columns <- stats::setNames(as.list(names(spike)), names(spike))
  }

  out <- stats::setNames(rep(NA_real_, length(spike)), names(spike))
  active_terms <- names(spike)[!is.na(spike) & spike]

  for (term in active_terms) {
    columns <- term_to_columns[[term]]

    if (is.null(columns) || length(columns) != 1L) {
      stop(
        sprintf(
          "Internal error: active SAZ term '%s' must map to exactly one design column.",
          term
        ),
        call. = FALSE
      )
    }

    column <- columns[[1L]]
    if (!column %in% colnames(x)) {
      stop(
        sprintf(
          "Internal error: design column '%s' for SAZ term '%s' is unavailable.",
          column, term
        ),
        call. = FALSE
      )
    }

    values <- x[, column]
    finite <- is.finite(values)
    if (any(finite)) {
      out[[term]] <- mean(values[finite] == 0)
    }
  }

  out
}

#' Reset Spike-at-Zero Indicators and Undo Cascade for Ineligible Variables
#'
#' Evaluates whether each variable flagged as a spike-at-zero (\code{spike})
#' is eligible for the SAZ algorithm based on representation of both components:
#' the zero component and the positive continuous component. For
#' ineligible variables, \code{spike} is reset to \code{FALSE} and the
#' \code{catzero} and \code{zero} flags are restored to the user's original
#' pre-cascade values. This correctly undoes the implicit cascade
#' \code{catzero[spike] <- TRUE}, \code{zero[catzero] <- TRUE} that
#' \code{fit_mfp()} applies before calling this function, ensuring that only
#' what the user explicitly requested is preserved.
#'
#' @section Background -- the cascade problem:
#' In \code{fit_mfp()}, the following cascade is applied before calling
#' \code{reset_spike()}:
#' \preformatted{
#'   catzero[spike]  <- TRUE   # spike implies catzero
#'   zero[catzero]   <- TRUE   # catzero implies zero
#' }
#' If \code{reset_spike()} only resets \code{spike}, \code{catzero} and
#' \code{zero} remain \code{TRUE} for ineligible variables even though the user
#' never asked for them. This function corrects that by accepting the
#' pre-cascade user values and restoring them for any variable where
#' \code{spike} is reset.
#'
#' @param x A numeric matrix or data frame with column names. Only columns named
#'   in \code{spike} are examined. Values are interpreted on the raw structural-
#'   zero scale: finite values \code{== 0} belong to the zero component and
#'   finite values \code{> 0} belong to the positive component. The caller does
#'   not need to materialize a temporarily zero-recoded copy.
#' @param spike A named logical vector. \code{TRUE} indicates the column was
#'   flagged as a spike-at-zero variable. Should be the post-cascade version
#'   after \code{catzero[spike] <- TRUE} has been applied in \code{fit_mfp()}.
#' @param user_catzero A named logical vector of the user's original
#'   \code{catzero} specification before the cascade
#'   \code{catzero[spike] <- TRUE} was applied. Must have the same names as
#'   \code{spike}.
#' @param user_zero A named logical vector of the user's original \code{zero}
#'   specification before the cascade \code{zero[catzero] <- TRUE} was applied.
#'   Must have the same names as \code{spike}.
#' @param min_saz_prop Numeric in \eqn{(0, 0.5)}. Minimum required
#'   proportion in each component of a spike-at-zero covariate: the
#'   structural-zero component and the positive continuous component. A
#'   requested spike-at-zero variable is retained only if both component
#'   proportions are at least this value.
#'
#' @return A named list with three elements, each a named logical vector of the
#'   same length as \code{spike}:
#' \describe{
#'   \item{\code{spike}}{Updated spike vector. Entries for ineligible variables
#'     are set to \code{FALSE}; all other entries are unchanged.}
#'   \item{\code{catzero}}{Updated catzero vector. For ineligible variables
#'     where \code{spike} was just reset, restored to \code{user_catzero}. All
#'     other entries retain their post-cascade values, i.e. \code{TRUE} for
#'     variables that remain eligible spike variables.}
#'   \item{\code{zero}}{Updated zero vector. For ineligible variables, restored
#'     to \code{user_zero}. All other entries retain their post-cascade values.}
#' }
#'
#' @section Cascade restoration examples:
#' For a variable whose \code{spike} is reset, the outcome depends on what the
#' user originally specified:
#' \describe{
#'   \item{User specified \code{spike} only, not \code{catzero} or
#'     \code{zero}}{After reset: \code{spike = FALSE},
#'     \code{catzero = FALSE}, \code{zero = FALSE}. The variable is treated as
#'     a standard continuous predictor -- no binary indicator or exact-zero handling.}
#'   \item{User specified \code{spike} and \code{zero = TRUE}}{After reset:
#'     \code{spike = FALSE}, \code{catzero = FALSE}, \code{zero = TRUE}.
#'     Exact-zero values are still structural zeros before FP transformation,
#'     but no binary indicator is added and the SAZ algorithm is not run.}
#'   \item{User specified \code{spike} and \code{catzero = TRUE}}{After reset:
#'     \code{spike = FALSE}, \code{catzero = TRUE}, \code{zero = TRUE}. The
#'     binary structural-zero indicator \code{I(x == 0)} is still added as explicitly
#'     requested via \code{catzero}, but the SAZ selection algorithm is not run.}
#' }
#'
#' @details
#' Two conditions cause a variable's spike indicator to be reset:
#' \enumerate{
#'   \item \strong{Insufficient component representation} -- either the
#'     structural-zero proportion or the positive-observation proportion is
#'     below \code{min_saz_prop}.
#'   \item \strong{Binary variable} -- exactly two effective finite values
#'     when exact zero and one distinct positive value are present.
#'     The positive part would then contain a single unique value, making FP
#'     transformation degenerate.
#' }
#' A warning is issued for each reset reason, identifying the affected variables.
#'
#' @keywords internal
#' @noRd
reset_spike <- function(x, spike, user_catzero, user_zero,
                        min_saz_prop = 0.10) {

  # Early exit: no spike variables to evaluate. In this case there was no
  # spike-implied cascade to preserve, so return the user's original zero and
  # catzero settings unchanged.
  if (!any(spike)) {
    return(list(spike = spike, catzero = user_catzero, zero = user_zero))
  }

  names_spike <- names(spike)[spike]

  # Component proportions for each requested spike-at-zero variable.
  # Work directly on the unre-coded values: x == 0 defines the structural-zero
  # component, while x > 0 defines the positive component. Avoiding a physical
  # rewrite prevents copy-on-modify from
  # duplicating the full design matrix solely for this eligibility check.
  #
  # is.finite() excludes NA/NaN/Inf from both counts, so n_observed below can
  # be smaller than nrow(x) if any values are non-finite; prop_zero/
  # prop_positive are proportions of the *observed* (finite) values, not of
  # all rows. This matters for the eligibility check further below: a variable
  # with many non-finite values could still be judged eligible or ineligible
  # based only on the finite subset.
  n_zero <- vapply(
    names_spike,
    function(v) {
      sum(is.finite(x[, v]) & x[, v] == 0)
    },
    integer(1L)
  )

  n_positive <- vapply(
    names_spike,
    function(v) {
      sum(is.finite(x[, v]) & x[, v] > 0)
    },
    integer(1L)
  )

  n_observed <- n_zero + n_positive

  prop_zero <- n_zero / n_observed
  prop_positive <- n_positive / n_observed

  # Binary variables are not eligible for spike-at-zero modelling. Match the
  # exact-zero definition without constructing a recoded vector: the finite
  # zero value contributes one level, while each distinct positive value
  # remains a separate level.
  is_binary <- vapply(
    names_spike,
    function(v) {
      values <- x[, v]
      values <- values[is.finite(values)]

      has_zero_component <- any(values == 0)
      n_positive_levels <- length(unique(values[values > 0]))
      n_effective_levels <- as.integer(has_zero_component) + n_positive_levels

      n_effective_levels == 2L
    },
    logical(1L)
  )

  # Identify variables whose spike option must be reset, by reason. A
  # variable can fail for either or both reasons; to_reset is the union so
  # each ineligible variable is only reset once regardless of how many
  # reasons apply. The two reason sets are kept separate only so the warning
  # messages below can report the actual cause(s) to the user.
  to_reset_component <- names_spike[
    is.na(prop_zero) |
      is.na(prop_positive) |
      prop_zero < min_saz_prop |
      prop_positive < min_saz_prop
  ]

  to_reset_binary <- names_spike[is_binary]
  to_reset <- union(to_reset_component, to_reset_binary)

  # Warn about component-representation resets. This walks through several
  # specific, mutually-exclusive causes (in order of how "extreme" the
  # shortfall is: no data at all, all-one-component, then falling just short
  # of the threshold on one side or the other) rather than always emitting the
  # same generic message, so the user gets a precise, actionable explanation
  # of *why* eligibility failed for each variable.
  if (length(to_reset_component) > 0L) {
    component_reason <- vapply(
      to_reset_component,
      function(v) {
        p0 <- prop_zero[[v]]
        pp <- prop_positive[[v]]
        nz <- n_zero[[v]]
        np <- n_positive[[v]]
        no <- n_observed[[v]]

        if (is.na(no) || no == 0L) {
          return("no finite zero or positive observations are available")
        }

        if (nz == 0L && np > 0L) {
          return(sprintf(
            "all observed values are positive; zero proportion is 0 < `min_saz_prop = %s`",
            min_saz_prop
          ))
        }

        if (np == 0L && nz > 0L) {
          return(sprintf(
            "all observed values are structural zeros; positive observation proportion is 0 < `min_saz_prop = %s`",
            min_saz_prop
          ))
        }

        if (is.na(p0) || p0 < min_saz_prop) {
          return(sprintf(
            "zero proportion is %.4g < `min_saz_prop = %s`",
            p0,
            min_saz_prop
          ))
        }

        if (is.na(pp) || pp < min_saz_prop) {
          return(sprintf(
            "positive observation proportion is %.4g < `min_saz_prop = %s`",
            pp,
            min_saz_prop
          ))
        }

        "component representation is insufficient"
      },
      character(1L)
    )

    warning(
      paste(
        c(
          "The spike-at-zero option has been reset for the following variable(s):",
          sprintf("- %s: %s", to_reset_component, component_reason)
        ),
        collapse = "\n"
      ),
      call. = FALSE
    )
  }

  # IMPORTANT: initialise from the post-cascade state, not from the raw user
  # inputs. For every spike variable that remains eligible, the internal state
  # must satisfy:
  #   spike   => catzero
  #   catzero => zero
  # This preserves the implicit cascade applied by fit_mfp():
  #   catzero[spike] <- TRUE
  #   zero[catzero]  <- TRUE
  #
  catzero <- user_catzero
  zero    <- user_zero
  catzero[spike] <- TRUE
  zero[catzero]  <- TRUE

  if (length(to_reset) > 0L) {
    # Reset spike only for variables that failed eligibility checks.
    spike[to_reset] <- FALSE

    # For reset variables only, restore catzero and zero to the user's explicit
    # pre-cascade choices. This undoes the spike-implied cascade only where the
    # spike option was rejected. Eligible spike variables keep catzero = TRUE
    # and zero = TRUE.
    catzero[to_reset] <- user_catzero[to_reset]
    zero[to_reset]    <- user_zero[to_reset]
  }

  # Re-apply the ordinary catzero -> zero cascade after restoring user choices.
  # This second application is necessary (not redundant with the one above):
  # for a variable in to_reset, catzero was just overwritten with
  # user_catzero[v], which could itself be TRUE if the user explicitly passed
  # that variable through catzero_vars independently of spike_vars. In that
  # case zero must be forced back to TRUE for it too, because catzero_vars
  # always implies zero_vars regardless of what happens to the spike request.
  zero[catzero] <- TRUE

  # Defensive internal consistency checks. These should never fail if the
  # cascade above has been applied correctly.
  if (any(spike & !catzero)) {
    stop("Internal error: spike variables must also have catzero = TRUE.", call. = FALSE)
  }
  if (any(catzero & !zero)) {
    stop("Internal error: catzero variables must also have zero = TRUE.", call. = FALSE)
  }

  list(spike = spike, catzero = catzero, zero = zero)
}

#' Resolve Spike-at-Zero Eligibility Before Preprocessing
#'
#' Internal helper used before shift and scale parameters are chosen.
#'
#' The spike-at-zero cascade implies that `spike` variables are temporarily
#' treated as `catzero`, and `catzero` variables are temporarily treated as
#' `zero`. This affects preprocessing because zero-handled variables should keep
#' the structural-zero component at zero.
#'
#' If a requested spike variable is later found ineligible for SAZ, it should
#' revert to the user's explicit `catzero` and `zero` choices before ordinary
#' shift/scale preprocessing is performed. This helper performs that early
#' eligibility resolution.
#'
#' @param x Raw numeric design matrix or data frame, before shift/scale
#'   preprocessing.
#' @param spike Named logical vector indicating requested spike-at-zero
#'   variables.
#' @param catzero Named logical vector indicating requested categorical-zero
#'   variables.
#' @param zero Named logical vector indicating requested zero-handled variables.
#' @param min_saz_prop Numeric in `(0, 0.5)`. Minimum required
#'   proportion in both the structural-zero and positive components.
#'
#' @return A named list with updated `spike`, `catzero`, and `zero` logical
#'   vectors.
#'
#' @keywords internal
#' @noRd
resolve_saz_eligibility <- function(x,
                                    spike,
                                    catzero,
                                    zero,
                                    min_saz_prop = 0.10) {
  # This is the *early* eligibility check, run by mfp2.default() before
  # shift/scale preprocessing (see the "Resolve spike-at-zero eligibility"
  # section of fit_mfp()'s documentation). Its whole purpose is to decide
  # spike/catzero/zero status up front, on the raw data, so that downstream
  # shift estimation (find_shift_factor()) sees the *final* zero/catzero
  # status for each variable and never picks a shift assuming spike-at-zero
  # handling that later turns out to be ineligible.
  if (!any(spike)) {
    return(list(
      spike = spike,
      catzero = catzero,
      zero = zero
    ))
  }

  # Save the user's pre-cascade choices so reset_spike() can restore them
  # for any variable that turns out ineligible.
  user_catzero <- catzero
  user_zero <- zero

  # Delegate directly on the raw matrix. reset_spike() interprets finite exact
  # zeros as the structural-zero component without allocating and rewriting a
  # full copy of x. Public callers reject negative values first.
  reset_spike(
    x = x,
    spike = spike,
    user_catzero = user_catzero,
    user_zero = user_zero,
    min_saz_prop = min_saz_prop
  )
}

#' Cap FP Degrees of Freedom for Spike-at-Zero Positive Components
#'
#' Internal helper used by \code{fit_mfp()} after \code{reset_spike()}.
#'
#' For ordinary variables, \code{assign_df()} caps the maximum FP degrees of
#' freedom according to the number of distinct values in the full variable. For
#' spike-at-zero variables, the continuous FP part is fitted to the positive
#' component, while exact-zero observations are represented structurally
#' through the binary zero indicator. Therefore the SAZ-specific df cap must be
#' based on the number of distinct positive values, not on the number of distinct
#' values in the full variable including zero.
#'
#' The rules mirror \code{assign_df()}:
#' \itemize{
#'   \item at most 3 distinct positive values: force linear, \code{df = 1};
#'   \item 4 or 5 distinct positive values: cap at FP1, \code{df = min(df, 2)};
#'   \item 6 or more distinct positive values: keep the requested df.
#' }
#'
#' @param x Numeric matrix or data frame on the fitting scale. For spike
#'   variables, only finite positive values are used to determine the positive-
#'   component cardinality; exact-zero values need not be physically recoded.
#' @param df Named integer vector of current maximum df values.
#' @param spike Named logical vector indicating final retained spike-at-zero
#'   variables after \code{reset_spike()}.
#'
#' @return Updated named integer vector of df values.
#'
#' @keywords internal
#' @noRd
cap_spike_df <- function(x, df, spike) {
  if (!any(spike)) {
    return(df)
  }

  # Only cap df for variables that are both flagged as spike (post-eligibility
  # resolution) and present in df's names; the intersect() guards against
  # spike containing names that df doesn't track (shouldn't normally happen,
  # but keeps this defensive rather than erroring on a name mismatch).
  spike_vars <- intersect(names(spike)[spike], names(df))

  if (length(spike_vars) == 0L) {
    return(df)
  }

  # Accumulate a record of every variable whose df was actually reduced, so a
  # single combined warning can be issued at the end (rather than one warning
  # per variable), and so the message can report the before/after df and the
  # positive-value count that triggered the change.
  changed <- character(0L)
  old_values <- integer(0L)
  new_values <- integer(0L)
  unique_values <- integer(0L)

  for (v in spike_vars) {
    # Cardinality is assessed only over the *positive* values of x[, v], not
    # the full column, because for a spike variable the FP/linear component is
    # only ever fitted to the positive part - the structural-zero part is
    # represented separately by the binary indicator. Using the full column's
    # cardinality here (as ordinary assign_df() does for non-spike variables)
    # would be wrong: it could, for example, see a variable with only 3
    # distinct positive values but many zeros, and wrongly conclude there's
    # enough variation for a high-degree FP, when the FP is actually only ever
    # evaluated on those 3 distinct positive values.
    positive_values <- x[, v]
    positive_values <- positive_values[
      is.finite(positive_values) & positive_values > 0
    ]

    u_positive <- length(unique(positive_values))

    old_df <- as.integer(df[[v]])
    new_df <- old_df

    # Same three-tier rule as assign_df(), just applied to the positive-only
    # count: <=3 distinct positive values forces linear (df=1); 4-5 caps at
    # FP1 (df=2, but never *raising* df, hence min(2L, old_df)); >=6 leaves
    # the requested df untouched.
    if (u_positive <= 3L) {
      new_df <- 1L
    } else if (u_positive <= 5L) {
      new_df <- min(2L, old_df)
    }

    if (!identical(new_df, old_df)) {
      changed <- c(changed, v)
      old_values <- c(old_values, old_df)
      new_values <- c(new_values, new_df)
      unique_values <- c(unique_values, u_positive)
      df[[v]] <- new_df
    }
  }

  if (length(changed) > 0L) {
    warning(
      "For the following spike-at-zero variable(s), the maximum FP df was ",
      "reduced because the positive component has few distinct values: ",
      paste0(
        changed,
        " (",
        unique_values,
        " distinct positive values: df ",
        old_values,
        " -> ",
        new_values,
        ")",
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }

  df
}


#' Rank Joint SAZ Information-Criterion Candidates by Simplicity
#'
#' Builds explicit tie-breaking metadata for the one-stage SAZ AIC/BIC
#' comparison. For ordinary FP candidates, positive-component complexity is
#' ordered as none, linear, FP1, FP2, and so on. ACD candidates use the same
#' functional-form hierarchy as their ordinary AIC/BIC selector. Component
#' count is retained as a final substantive tie-break after functional-form
#' complexity.
#'
#' Keeping this metadata separate from candidate row positions ensures that
#' reordering the joint metrics table cannot change the selected model.
#'
#' @param candidate_names Character vector of joint SAZ candidate names.
#' @param acd Logical scalar indicating whether the positive component uses
#'   the ACD candidate family.
#'
#' @return A data frame with integer columns `positive_complexity` and
#'   `component_count`, indexed by `candidate_names`.
#'
#' @keywords internal
#' @noRd
saz_ic_candidate_simplicity <- function(candidate_names, acd = FALSE) {
  if (!is.character(candidate_names) || anyNA(candidate_names) ||
      anyDuplicated(candidate_names)) {
    stop(
      "Internal error: joint SAZ candidate names must be unique strings.",
      call. = FALSE
    )
  }

  has_binary <- grepl(" \\+ Binary$", candidate_names)
  base_names <- sub(" \\+ Binary$", "", candidate_names)
  positive_complexity <- rep(NA_integer_, length(candidate_names))

  no_positive <- base_names %in% c("null", "Binary")
  positive_complexity[no_positive] <- 0L

  if (isTRUE(acd)) {
    # The ACD selector's established order runs from ordinary linear through
    # the two single-transformation FP1 forms to the joint FP1(x, A(x)) form.
    acd_forms <- c(
      "linear",
      "linear(., A(x))",
      "FP1(x, .)",
      "FP1(., A(x))",
      "FP1(x, A(x))"
    )
    is_positive <- !no_positive
    positive_complexity[is_positive] <- match(
      base_names[is_positive],
      acd_forms
    )
  } else {
    positive_complexity[base_names == "linear"] <- 1L
    is_fp <- grepl("^FP[0-9]+$", base_names)
    positive_complexity[is_fp] <-
      as.integer(sub("^FP", "", base_names[is_fp])) + 1L
  }

  if (anyNA(positive_complexity)) {
    stop(
      "Internal error: unrecognised joint SAZ candidate form: ",
      paste(candidate_names[is.na(positive_complexity)], collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  component_count <- ifelse(
    base_names == "null",
    0L,
    ifelse(base_names == "Binary", 1L, ifelse(has_binary, 2L, 1L))
  )

  out <- data.frame(
    positive_complexity = as.integer(positive_complexity),
    component_count = as.integer(component_count),
    row.names = candidate_names
  )
  out
}


#' Select the Winning Joint SAZ Information-Criterion Candidate
#'
#' Applies an explicit, order-independent hierarchy to the eligible one-stage
#' SAZ candidates: minimum AIC/BIC, smaller adjusted model df, simpler
#' positive-component form, and fewer retained components. Candidate order is
#' used only as a defensive fallback if every substantive quantity is equal.
#'
#' @param metrics Numeric candidate metrics matrix with row names and columns
#'   `df` plus the requested information criterion.
#' @param criterion Character scalar, `"aic"` or `"bic"`.
#' @param eligible Integer indices of eligible rows in `metrics`.
#' @param acd Logical scalar indicating whether the positive component uses
#'   the ACD candidate family.
#'
#' @return Integer row index of the selected candidate.
#'
#' @keywords internal
#' @noRd
select_saz_ic_winner <- function(metrics, criterion, eligible, acd = FALSE) {
  criterion <- tolower(criterion)
  if (!criterion %in% c("aic", "bic")) {
    stop("Internal error: joint SAZ selection requires AIC or BIC.",
         call. = FALSE)
  }
  if (is.null(rownames(metrics)) || !all(c("df", criterion) %in% colnames(metrics))) {
    stop("Internal error: incomplete joint SAZ candidate metrics.",
         call. = FALSE)
  }

  eligible <- as.integer(eligible)
  ic_values <- metrics[eligible, criterion, drop = TRUE]
  valid_ic <- is.finite(ic_values)
  if (!any(valid_ic)) {
    stop("All eligible SAZ information-criterion values are non-finite.",
         call. = FALSE)
  }

  best_ic <- min(ic_values[valid_ic])
  tied <- eligible[valid_ic & ic_values == best_ic]

  # First prefer lower adjusted model df. These are the same df used in the
  # AIC/BIC penalty, including the extra allowance for FP-power selection.
  if (length(tied) > 1L) {
    tied_df <- metrics[tied, "df", drop = TRUE]
    finite_df <- is.finite(tied_df)
    if (any(finite_df)) {
      best_df <- min(tied_df[finite_df])
      tied <- tied[finite_df & tied_df == best_df]
    }
  }

  simplicity <- saz_ic_candidate_simplicity(
    rownames(metrics),
    acd = acd
  )

  # Equal adjusted df can still represent different kinds of complexity. For
  # ordinary FP SAZ selection, the only such substantive pairs are binary-only
  # versus positive-only linear (1 df), and both-components linear versus
  # positive-only FP1 (2 df). Prefer the simpler positive-component form.
  if (length(tied) > 1L) {
    tied_positive <- simplicity$positive_complexity[tied]
    tied <- tied[tied_positive == min(tied_positive)]
  }

  # This normally resolves no additional ordinary-FP ties, but makes the
  # intended hierarchy complete for any candidate family with equal criterion,
  # adjusted df, and positive-form complexity.
  if (length(tied) > 1L) {
    tied_components <- simplicity$component_count[tied]
    tied <- tied[tied_components == min(tied_components)]
  }

  tied[[1L]]
}


#' Joint AIC/BIC Selection for a Spike-at-Zero Term
#'
#' Fits the complete information-criterion candidate family for one eligible
#' spike-at-zero term. Positive-only and positive-plus-binary models perform
#' independent functional-form searches; the null and binary-only models are
#' then added before one global AIC/BIC minimisation.
#'
#' For an ordinary FP term with maximum degree `d`, the compressed comparison
#' contains `2*d + 4` rows: null, binary, linear with and without binary, and
#' the best FP1 through FPd models with and without binary. "Best" within a
#' fixed degree means minimum deviance, which is also minimum AIC/BIC because
#' all power combinations at that degree have the same complexity.
#'
#' Exact AIC/BIC ties are resolved by smaller adjusted model df, simpler
#' positive-component form, and then fewer retained components. The explicit
#' hierarchy prevents candidate row order from determining a substantive tie.
#'
#' @inheritParams find_best_fp_step
#' @param ... Additional arguments passed through to the ordinary selectors
#'   and ultimately to `fit_model()`.
#'
#' @return A selector result compatible with `find_best_fp_step()`, including
#'   the joint metrics table and the final `spike_decision`.
#'
#' @keywords internal
#' @noRd
select_saz_ic <- function(x,
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
  criterion_lower <- tolower(criterion)
  if (!criterion_lower %in% c("aic", "bic")) {
    stop("select_saz_ic() requires criterion = 'aic' or 'bic'.", call. = FALSE)
  }
  if (!isTRUE(spike[[xi]]) || is.null(catzero[[xi]])) {
    stop("select_saz_ic() requires an eligible spike-at-zero term.", call. = FALSE)
  }

  # SAZ decision codes are an integer-valued package contract (1L, 2L, 3L).
  # Normalize the complete named vector here because assigning an integer code
  # into a caller-supplied double vector would otherwise silently retain
  # storage mode "double" in this selector's result.
  spike_decision <- stats::setNames(
    as.integer(spike_decision),
    names(spike_decision)
  )

  # Select the appropriate functional-form engine. degree = 0.5 corresponds
  # to df = 1, where linear is the maximum and only continuous form.
  branch_selector <- if (degree <= 0.5) {
    select_linear
  } else if (isTRUE(acdx[[xi]])) {
    select_ic_acd
  } else {
    select_ic
  }

  branch_args <- list(
    x = x, xi = xi, keep = keep, degree = degree, acdx = acdx,
    y = y, powers_current = powers_current, powers = powers,
    criterion = criterion, ftest = ftest, select = select, alpha = alpha,
    family = family, family_string = family_string, zero = zero,
    catzero = catzero, spike = spike, spike_decision = spike_decision,
    acd_parameter = acd_parameter, prev_adj_params = prev_adj_params,
    transform_cache = transform_cache, force_max_fp = force_max_fp,
    has_offset = has_offset, n_obs = n_obs,
    term_to_columns = term_to_columns, ...
  )

  # Joint path A: optimize every positive-plus-binary candidate. This is the
  # historical Stage-1 family, but here it is only one branch of a joint IC
  # comparison and does not constrain the powers used by the next branch.
  fit_bp <- do.call(branch_selector, branch_args)

  # Joint path B: remove the binary component first, then repeat the complete
  # continuous-form search. Assigning through `[` preserves xi's named NULL
  # list element; `[[<- NULL` would delete it and misalign later variables.
  positive_args <- branch_args
  positive_args$catzero[xi] <- list(NULL)
  positive_args$spike[[xi]] <- FALSE
  positive_args$spike_decision[[xi]] <-
    saz_decision_codes[["continuous_only"]]
  positive_args$transform_cache <- fit_bp$transform_cache
  fit_p <- do.call(branch_selector, positive_args)

  # Joint path C: fit binary-only using the same adjustment-variable design.
  # Its single indicator coefficient needs no extra FP-power df penalty.
  params_binary <- fit_bp$current_adj_params
  adjustment_matrix <- params_binary[[xi]]$data_adj
  if (is.null(adjustment_matrix) || NCOL(adjustment_matrix) == 0L) {
    adjustment_matrix <- NULL
  }

  binary_matrix <- catzero[[xi]]
  if (!is.matrix(binary_matrix)) {
    binary_matrix <- as.matrix(binary_matrix)
  }
  colnames(binary_matrix) <- "catzero"

  is_cox <- identical(family_string, "cox")
  x_binary <- if (is_cox) {
    if (is.null(adjustment_matrix)) binary_matrix else
      cbind(binary_matrix, adjustment_matrix)
  } else {
    assemble_design_matrix(
      blocks = list(binary_matrix, adjustment_matrix),
      nobs = NROW(y),
      intercept = TRUE
    )
  }

  binary_fit_args <- list(
    x = x_binary, y = y, family = family, family_string = family_string,
    has_offset = has_offset
  )
  if (!is_cox) {
    binary_fit_args$x_has_intercept <- TRUE
  }
  binary_fit_args <- c(binary_fit_args, list(...))
  fit_binary <- do.call(fit_model, binary_fit_args)
  metrics_binary <- rbind(
    Binary = calculate_model_metrics(fit_binary, n_obs = n_obs)
  )
  params_binary[[xi]]$data_xi <- binary_matrix

  # Assemble one table. The rows are deliberately stable: null, binary, then
  # positive-only and positive-plus-binary forms interleaved by complexity.
  # For ordinary FP degree d this gives exactly 2*d + 4 candidates and works
  # unchanged for FP3, FP4, or any larger degree supported by the package.
  p_names <- setdiff(rownames(fit_p$metrics), "null")
  bp_names <- setdiff(rownames(fit_bp$metrics), "null")
  if (length(p_names) != length(bp_names)) {
    stop("Internal error: SAZ IC branches produced unequal candidate sets.",
         call. = FALSE)
  }

  interleaved_names <- as.vector(rbind(p_names, bp_names))
  candidate_names <- c("null", "Binary", interleaved_names)
  metric_rows <- c(
    list(fit_bp$metrics["null", , drop = FALSE], metrics_binary),
    lapply(interleaved_names, function(model_name) {
      source_fit <- if (model_name %in% bp_names) fit_bp else fit_p
      source_fit$metrics[model_name, , drop = FALSE]
    })
  )
  metrics <- do.call(rbind, metric_rows)
  rownames(metrics) <- candidate_names

  n_power <- max(NCOL(fit_p$powers), NCOL(fit_bp$powers), 1L)
  pad_power_row <- function(fit, model_name) {
    ensure_length(fit$powers[model_name, , drop = FALSE], n_power)
  }
  power_rows <- c(
    list(
      ensure_length(fit_bp$powers["null", , drop = FALSE], n_power),
      rep(NA_real_, n_power)
    ),
    lapply(interleaved_names, function(model_name) {
      source_fit <- if (model_name %in% bp_names) fit_bp else fit_p
      pad_power_row(source_fit, model_name)
    })
  )
  powers_joint <- do.call(rbind, power_rows)
  rownames(powers_joint) <- candidate_names

  # Minimize the requested IC globally. Exact ties are resolved explicitly by
  # smaller adjusted model df, simpler positive-component form, and then fewer
  # retained components. This prevents candidate row order from acting as an
  # implicit statistical preference between equally scoring models.
  eligible <- seq_len(nrow(metrics))
  if (xi %in% keep) {
    eligible <- eligible[rownames(metrics)[eligible] != "null"]
  }
  model_best <- select_saz_ic_winner(
    metrics = metrics,
    criterion = criterion_lower,
    eligible = eligible,
    acd = isTRUE(acdx[[xi]])
  )
  selected_name <- rownames(metrics)[model_best]

  if (selected_name == "null") {
    selected_params <- fit_bp$current_adj_params
    selected_decision <- saz_decision_codes[["continuous_only"]]
  } else if (selected_name == "Binary") {
    selected_params <- params_binary
    selected_decision <- saz_decision_codes[["binary_only"]]
  } else if (selected_name %in% bp_names) {
    selected_params <- fit_bp$current_adj_params
    selected_decision <- saz_decision_codes[["cont_binary"]]
  } else {
    selected_params <- fit_p$current_adj_params
    selected_decision <- saz_decision_codes[["continuous_only"]]
  }
  spike_decision[[xi]] <- selected_decision

  list(
    keep = xi %in% keep,
    acd = isTRUE(acdx[[xi]]),
    powers = powers_joint,
    # A selector returns the selected powers as a numeric vector. Keep the
    # complete candidate collection as a matrix in `powers`, but drop the row
    # dimension for the single winning candidate. This also matches the
    # documented selector contract and the output of select_linear().
    power_best = powers_joint[model_best, , drop = TRUE],
    metrics = metrics,
    model_best = model_best,
    statistic = NA,
    pvalue = NA,
    zero = zero[xi],
    catzero = TRUE,
    spike = TRUE,
    selection_mode = "joint_ic",
    spike_decision = spike_decision,
    current_adj_params = selected_params,
    transform_cache = fit_p$transform_cache
  )
}

#' Fit Reduced Models for P-value SAZ Stage 2
#'
#' Reuses the Stage-1 transformed focal and adjustment matrices, fits the
#' continuous-only and binary-only reductions, and returns them with the
#' already fitted full model. This helper is not used by joint AIC/BIC SAZ
#' selection.
#'
#' @param stage1_selection P-value Stage-1 selector result containing
#'   `current_adj_params[[xi]]$data_xi` and `$data_adj`.
#' @param xi Character scalar naming the focal term.
#' @param y,weights,offset,family,family_string,method,strata,nocenter,control,
#'   rownames,has_offset Passed to `fit_model()`.
#' @param calculate_gaussian_deviance Logical; compute Gaussian deviance for
#'   Stage-2 F-tests.
#' @param fitter Fitting backend.
#'
#' @return The full Stage-1 fit, two reduced fits, and reused design matrices.
#' @keywords internal
#' @noRd
fit_saz_reduced_models <- function(stage1_selection,
                                   xi,
                                   y,
                                   weights,
                                   offset,
                                   family,
                                   family_string,
                                   method,
                                   strata,
                                   nocenter,
                                   control,
                                   rownames,
                                   has_offset,
                                   calculate_gaussian_deviance = FALSE,
                                   fitter = "base") {
  # Step 1: Recover xi's already-transformed stage-1 design from the cache --
  # find_best_fpm_step() (called during stage 1) stores the winning
  # transformed xi matrix under current_adj_params[[xi]]$data_xi precisely so
  # that stage 2 can reuse it here instead of re-deriving the FP/ACD
  # transformation and re-shifting/scaling xi from scratch.
  params_xi <- stage1_selection$current_adj_params[[xi]]

  data_xi <- params_xi$data_xi
  data_xi_colnames <- colnames(data_xi)

  # This should be unreachable in normal operation: evaluate_saz_stage2() is
  # only invoked (from find_best_fp_step()) when xi is both a spike variable
  # and was selected in stage 1, which guarantees a "catzero" column was
  # generated for it. A missing column here would indicate stage 1 and stage
  # 2 have gotten out of sync internally.
  if (!"catzero" %in% data_xi_colnames) {
    stop(
      "Internal error: SAZ stage 2 requires a 'catzero' column in the ",
      "selected xi design matrix.",
      call. = FALSE
    )
  }

  # Record the two stage-2 column groups, but do not materialize both
  # component matrices yet. Model 2 is fitted before Model 3, so extracting
  # the binary component only after Model 2 has finished avoids keeping an
  # otherwise-unused n-row copy alive during the first reduced-model fit.
  continuous_cols <- data_xi_colnames != "catzero"

  # Step 2: Recover the (also already-transformed) adjustment matrix --------
  # Reused from stage 1 for the same reason as data_xi: adjustment-variable
  # transformations don't depend on which SAZ component of xi is retained, so
  # there is no need to recompute them for Model 2/Model 3.
  adjustment_matrix <- params_xi$data_adj

  # Normalize the "no adjustment variables" case: build_adjustment_step()'s
  # return contract allows either NULL or an n x 0 matrix depending on
  # whether there were adjustment variables at all; treat both as NULL here
  # so the cbind() logic below only has one case to handle.
  if (is.null(adjustment_matrix) || ncol(adjustment_matrix) == 0L) {
    adjustment_matrix <- NULL
  }

  # Step 3: Prepare the shared fitting arguments --------------------------
  fit_args <- list(
    y = y,
    family = family,
    family_string = family_string,
    fitter = fitter,
    weights = weights,
    offset = offset,
    method = method,
    strata = strata,
    nocenter = nocenter,
    has_offset = has_offset,
    calculate_gaussian_deviance = calculate_gaussian_deviance,
    control = control,
    rownames = rownames
  )

  # Step 4: Assemble and fit the two reduced models sequentially ----------
  # For GLMs, fit_glm() requires an intercept. The historical path first
  # materialized [SAZ component | adjustment] here and fit_glm() immediately
  # copied that whole matrix again merely to prepend the intercept. Assemble
  # the final [Intercept | SAZ component | adjustment] matrix directly and
  # tell fit_model() that the intercept is already present. This removes one
  # full design-matrix allocation whenever adjustment variables are present.
  #
  # Cox models are intentionally left on the historical no-intercept path:
  # coxph.fit() must not receive an ordinary intercept, and when there are no
  # adjustment variables the component matrix can be passed through without
  # an additional assembly allocation.
  is_cox <- identical(family_string, "cox")
  if (!is_cox) {
    fit_args$x_has_intercept <- TRUE
  }

  xi_continuous <- data_xi[, continuous_cols, drop = FALSE]
  x_fit2 <- if (is_cox) {
    if (is.null(adjustment_matrix)) {
      xi_continuous
    } else {
      cbind(xi_continuous, adjustment_matrix)
    }
  } else {
    assemble_design_matrix(
      blocks = list(xi_continuous, adjustment_matrix),
      nobs = NROW(y),
      intercept = TRUE
    )
  }
  # x_fit2 now owns (or references, for the Cox/no-adjustment case) everything
  # required for Model 2; the extracted continuous-component copy is no longer
  # needed as a separate local binding.
  rm(xi_continuous, envir = environment())

  fit2 <- do.call(fit_model, c(list(x = x_fit2), fit_args))
  # Do not retain Model 2's temporary design while constructing Model 3.
  rm(x_fit2, envir = environment())

  # Materialize the binary component only now, after Model 2 and its temporary
  # design have been released, to keep the stage-2 live working set small.
  xi_binary <- data_xi[, "catzero", drop = FALSE]
  x_fit3 <- if (is_cox) {
    if (is.null(adjustment_matrix)) {
      xi_binary
    } else {
      cbind(xi_binary, adjustment_matrix)
    }
  } else {
    assemble_design_matrix(
      blocks = list(xi_binary, adjustment_matrix),
      nobs = NROW(y),
      intercept = TRUE
    )
  }
  rm(xi_binary, envir = environment())

  fit3 <- do.call(fit_model, c(list(x = x_fit3), fit_args))
  # Neither reduced-model design is consumed after fitting. Release Model 3's
  # matrix before constructing the compact return value below.
  rm(x_fit3, envir = environment())

  # Step 5: Return only objects required by stage-2 scoring/decision-making.
  # In particular, do not retain the two temporary reduced-model design
  # matrices: they can be much larger than the fit statistics and are never
  # consumed by production code after fit2/fit3 have been obtained.
  list(
    fit1 = stage1_selection,
    fit2 = fit2,
    fit3 = fit3,
    data_xi = data_xi,
    adjustment_matrix = adjustment_matrix
  )
}
#' Compute Model Metrics for Candidate Spike-at-zero Models
#'
#' This function computes fit statistics for the three candidate models
#' used in the spike-at-zero (SAZ) algorithm. Model 1 metrics are extracted
#' directly from the previously fitted model, while metrics for Model 2
#' (FPm/linear only) and Model 3 (binary-only) are computed using
#' `calculate_model_metrics`.
#'
#' @param fit1 Fitted object for Model 1 (complex model from stage 1 of SAZ).
#' @param fit2 Fitted object for Model 2 (FPm/linear only plus adjusted covariates).
#' @param fit3 Fitted object for Model 3 (binary-only plus adjusted covariates).
#' @param n_obs Number of observations in the dataset.
#' @param power_best Numeric vector of selected powers for the best FP terms
#' from stage 1 of SAZ algorithm.
#'
#' @details
#' The function determines the degree of the fractional polynomial based on `power_best`.
#' Model 1 metrics are retrieved from the best-fit row of the `fit1` object.
#' Metrics for Models 2 and 3 are calculated using `calculate_model_metrics`.
#'
#' @return A list with three elements:
#'   * `metrics1`: Fit statistics for Model 1.
#'   * `metrics2`: Fit statistics for Model 2.
#'   * `metrics3`: Fit statistics for Model 3.
#'
#' @keywords internal
#' @noRd
compute_saz_stage2_metrics <- function(fit1, fit2, fit3, n_obs, power_best) {

  # power_best is stage 1's selected power vector for xi, which can contain
  # NA components for inactive ACD terms (e.g. c(NA, 1) for "linear only in
  # A(x)", where the x-component was dropped). Dropping NAs and counting what
  # remains gives the number of powers actually estimated for the surviving
  # continuous component - this becomes the AIC/BIC df_additional penalty for
  # Model 2 below (see calculate_model_metrics()'s df_additional: each
  # estimated FP power beyond an ordinary linear coefficient costs one extra
  # df). ACD can produce NA like c(NA,1) so degree will reduce to 1 and in
  # this case additional parameters = 0 since the power = 1, see mfpa paper
  power_best <- power_best[!is.na(power_best)]
  degree <- length(power_best)

  # Extract Model 1 metrics safely
  if (is.null(fit1$metrics) || is.null(fit1$model_best)) {
    stop("fit1 must contain 'metrics' and 'model_best' elements.")
  }

  if (fit1$model_best > nrow(fit1$metrics) || fit1$model_best < 1) {
    stop("fit1$model_best is out of bounds for fit1$metrics.")
  }

  # Model 1's metrics were already computed during stage 1 model selection;
  # just pick out the row for whichever candidate stage 1 actually selected
  # (fit1$model_best), rather than recomputing.
  metrics1 <- fit1$metrics[fit1$model_best, ]

  # Compute Model 2 metrics with degree adjustment: deg2 = 0 only in the
  # special case where the surviving continuous component is exactly linear
  # (a single power equal to 1) - a linear term is an ordinary regression
  # coefficient with no extra estimated power, so it earns no df_additional
  # penalty. Any other surviving power vector (nonlinear FP, or an ACD power
  # that isn't exactly 1) means `degree` estimated powers were spent and each
  # contributes one extra df, exactly as for ordinary (non-spike) FP terms.
  deg2 <- if (length(power_best) == 1L && power_best == 1) {
    0L
  } else {
    degree
  }

  metrics2 <- tryCatch(
    calculate_model_metrics(fit2, n_obs, deg2),
    error = function(e) stop("Failed to compute metrics for fit2: ", e$message)
  )

  # Compute Model 3 metrics: no df_additional argument is passed here (it
  # defaults to 0 in calculate_model_metrics()), because Model 3 is just the
  # binary zero-indicator - a single ordinary coefficient, with no FP power to
  # estimate and therefore no extra df to account for.
  metrics3 <- tryCatch(
    calculate_model_metrics(fit3, n_obs),
    error = function(e) stop("Failed to compute metrics for fit3: ", e$message)
  )
  return(list(metrics1 = metrics1, metrics2 = metrics2, metrics3 = metrics3))
}

#' Compute Stage 2 Spike-at-Zero Model-Selection Decision
#'
#' This internal helper compares competing regression models according to a
#' specified selection criterion and returns the stage 2 spike-at-zero (SAZ)
#' decision. Stage 2 decides whether to retain both the continuous component
#' and the binary zero-indicator component, or whether one of these components
#' can be removed.
#'
#' @param metrics A list containing model fit statistics for three candidate
#'   models:
#'   \itemize{
#'     \item \code{metrics1}: Model 1, the full SAZ model selected in stage 1,
#'       usually containing both the continuous FP/linear/ACD component and the
#'       binary zero-indicator component.
#'     \item \code{metrics2}: Model 2, the continuous FP/linear/ACD component
#'       only.
#'     \item \code{metrics3}: Model 3, the binary zero-indicator component only.
#'   }
#'   Each element must contain named values for:
#'   \itemize{
#'     \item \code{logl}: log-likelihood;
#'     \item \code{df}: model degrees of freedom;
#'     \item \code{aic}: Akaike information criterion;
#'     \item \code{bic}: Bayesian information criterion;
#'     \item \code{deviance_gaussian}: Gaussian deviance, required when
#'       \code{ftest = TRUE};
#'     \item \code{df_resid}: residual degrees of freedom, required when
#'       \code{ftest = TRUE}.
#'   }
#' @param criterion Character string specifying the selection criterion. Must be
#'   one of \code{"pvalue"}, \code{"aic"}, or \code{"bic"}.
#' @param alpha Numeric significance threshold used for the Stage-2
#'   component-removal tests when \code{criterion = "pvalue"}. A component is
#'   removed only when its comparison p-value is strictly greater than
#'   \code{alpha}; therefore \code{alpha = 1} also preserves both SAZ
#'   components at the exact \code{p = 1} boundary. This is the same
#'   functional-form significance level used in Stage 1, not the
#'   variable-inclusion threshold \code{select}.
#' @param n_obs Integer. Number of observations, used for F-tests.
#' @param ftest Logical. If \code{TRUE}, use F-tests instead of likelihood-ratio
#'   tests when \code{criterion = "pvalue"}.
#'
#' @details
#' When \code{criterion = "pvalue"}, Model 1 is treated as the full SAZ model,
#' while Models 2 and 3 are treated as reduced alternatives. Two nested
#' comparisons are performed.
#'
#' The first comparison, Model 2 versus Model 1, tests whether the binary
#' zero-indicator component adds information beyond the continuous FP/linear/ACD
#' component. A small p-value means that dropping the binary component makes the
#' model significantly worse.
#'
#' The second comparison, Model 3 versus Model 1, tests whether the continuous
#' FP/linear/ACD component adds information beyond the binary zero-indicator
#' component. A small p-value means that dropping the continuous component makes
#' the model significantly worse.
#'
#' Let \code{p_drop_binary} denote the p-value for Model 2 versus Model 1, and
#' let \code{p_drop_continuous} denote the p-value for Model 3 versus Model 1.
#'
#' The p-value decision rule is:
#'
#' \tabular{llll}{
#'   \strong{Condition} \tab \strong{Interpretation} \tab
#'   \strong{Retained component(s)} \tab \strong{Selected model} \cr
#'   \code{p_drop_binary <= alpha}, \code{p_drop_continuous <= alpha} \tab
#'   Both reductions are significantly worse than Model 1 \tab
#'   Continuous FP/linear/ACD + binary zero indicator \tab
#'   Model 1 \cr
#'   \code{p_drop_binary <= alpha}, \code{p_drop_continuous > alpha} \tab
#'   Removing the binary component is harmful; removing the continuous component is acceptable \tab
#'   Binary zero indicator only \tab
#'   Model 3 \cr
#'   \code{p_drop_binary > alpha}, \code{p_drop_continuous <= alpha} \tab
#'   Removing the continuous component is harmful; removing the binary component is acceptable \tab
#'   Continuous FP/linear/ACD only \tab
#'   Model 2 \cr
#'   \code{p_drop_binary > alpha}, \code{p_drop_continuous > alpha} \tab
#'   Neither reduction is significantly worse than Model 1 \tab
#'   Reduced component selected by a complexity-adjusted comparison \tab
#'   Model 2 if \code{BIC(Model 2) < BIC(Model 3)}, otherwise Model 3 \cr
#' }
#'
#' If \code{ftest = FALSE}, the nested comparisons use likelihood-ratio tests.
#' If \code{ftest = TRUE}, the nested comparisons use F-tests. In the case
#' where both reductions are non-significant, Models 2 and 3 are
#' non-nested and can have different complexity: Model 3 has a one-df binary
#' component, whereas Model 2 can carry a one-df linear component or a
#' \eqn{2m}-df FP\emph{m} component selected in Stage 1, for any supported
#' degree \emph{m}. Their raw deviances or log-likelihoods are
#' therefore not used as a neutral tie-break. Instead, the model with the
#' smaller BIC is selected so that the additional continuous-component
#' complexity is penalized. An exact BIC tie selects Model 3, the binary-only
#' representation.
#'
#' @return A list with two elements:
#'   \itemize{
#'     \item \code{decision}: Integer indicating the selected model:
#'       \code{1} = both components, \code{2} = continuous FP/linear/ACD only,
#'       and \code{3} = binary zero-indicator only.
#'     \item \code{pvalue}: Named numeric vector containing
#'       \code{p_drop_binary} and \code{p_drop_continuous} when
#'       \code{criterion = "pvalue"}.
#'   }
#'
#' @keywords internal
#' @noRd
compute_saz_stage2_decision <- function(metrics,
                                        criterion,
                                        alpha,
                                        n_obs,
                                        ftest = FALSE) {
  criterion <- tolower(criterion)

  if (!identical(criterion, "pvalue")) {
    stop(
      "compute_saz_stage2_decision() is only for p-value SAZ Stage 2.",
      call. = FALSE
    )
  }

  if (ftest) {
      # Test whether the binary zero-indicator component is needed.
      # Reduced model: Model 2 = continuous FP/linear/ACD only.
      # Full model:    Model 1 = continuous FP/linear/ACD + binary.
      stats1 <- calculate_f_test(
        deviances = c(
          metrics$metrics2["deviance_gaussian"],
          metrics$metrics1["deviance_gaussian"]
        ),
        dfs_resid = c(
          metrics$metrics2["df_resid"],
          metrics$metrics1["df_resid"]
        ),
        n_obs = n_obs
      )

      # Test whether the continuous FP/linear/ACD component is needed.
      # Reduced model: Model 3 = binary zero-indicator only.
      # Full model:    Model 1 = continuous FP/linear/ACD + binary.
      stats2 <- calculate_f_test(
        deviances = c(
          metrics$metrics3["deviance_gaussian"],
          metrics$metrics1["deviance_gaussian"]
        ),
        dfs_resid = c(
          metrics$metrics3["df_resid"],
          metrics$metrics1["df_resid"]
        ),
        n_obs = n_obs
      )
    } else {
      # Test whether the binary zero-indicator component is needed.
      # Reduced model: Model 2 = continuous FP/linear/ACD only.
      # Full model:    Model 1 = continuous FP/linear/ACD + binary.
      stats1 <- calculate_lr_test(
        logl = c(metrics$metrics2["logl"], metrics$metrics1["logl"]),
        dfs = c(metrics$metrics2["df"], metrics$metrics1["df"])
      )

      # Test whether the continuous FP/linear/ACD component is needed.
      # Reduced model: Model 3 = binary zero-indicator only.
      # Full model:    Model 1 = continuous FP/linear/ACD + binary.
      stats2 <- calculate_lr_test(
        logl = c(metrics$metrics3["logl"], metrics$metrics1["logl"]),
        dfs = c(metrics$metrics3["df"], metrics$metrics1["df"])
      )
    }

    p_drop_binary <- stats1$pvalue
    p_drop_continuous <- stats2$pvalue

    # This mirrors the decision table in the function documentation above.
    # Reduced components are accepted only when their p-value strictly exceeds
    # alpha. This deliberately follows the MFP endpoint convention used by the
    # Stage-1 closed tests: alpha = 1 is a forcing value, so an exact p = 1
    # retains the corresponding full-model component by convention rather than
    # triggering simplification. So:
    #   both p-values small  -> neither component can be safely dropped -> Model 1
    #   only p_drop_binary small -> the binary component is needed but the
    #     continuous one isn't -> keep binary only -> Model 3
    #   only p_drop_continuous small -> the continuous component is needed but
    #     the binary one isn't -> keep continuous only -> Model 2
    #   neither small -> both components are individually dispensable. Models
    #     2 and 3 are non-nested and may have different df, so raw deviance
    #     would systematically favor the more flexible continuous component.
    #     Compare their already-computed BIC values instead. The strict
    #     inequality deliberately sends an exact BIC tie to Model 3, whose
    #     binary-only representation is no more complex than Model 2.
    decision <- if (p_drop_binary <= alpha && p_drop_continuous <= alpha) {
      saz_decision_codes[["cont_binary"]]  # Model 1: both components
    } else if (p_drop_binary <= alpha && p_drop_continuous > alpha) {
      saz_decision_codes[["binary_only"]]  # Model 3: binary zero-indicator only
    } else if (p_drop_binary > alpha && p_drop_continuous <= alpha) {
      saz_decision_codes[["continuous_only"]]  # Model 2: continuous FP/linear/ACD only
    } else {
      # If neither reduced model is significantly worse than the full model,
      # use BIC to compare the non-nested reduced models while accounting for
      # the Stage-1 continuous component's one-df linear or 2m-df FPm
      # complexity, for any supported FP degree m.
      if (metrics$metrics2["bic"] < metrics$metrics3["bic"]) {
        saz_decision_codes[["continuous_only"]]
      } else {
        saz_decision_codes[["binary_only"]]
      }
    }

  list(
    decision = decision,
    pvalue = c(
      p_drop_binary = p_drop_binary,
      p_drop_continuous = p_drop_continuous
    )
  )
}


#' Evaluate Stage 2 of the Spike-at-Zero (SAZ) Algorithm for One Variable
#'
#' Internal helper used by \code{find_best_fp_step()}.
#'
#' Stage 2 is specific to `criterion = "pvalue"` and is evaluated only when
#' stage 1 selected a
#' non-null functional form for a spike-at-zero variable \code{xi}. This
#' function fits the two reduced candidate models (continuous-only and
#' binary-only), compares them against the full continuous + binary model
#' already selected in stage 1, updates \code{spike_decision[[xi]]}
#' accordingly, and builds the printable stage-2 metrics table.
#'
#' @param fit1 Stage 1 selection object returned by a p-value
#'   \code{select_*()} function (for example, \code{select_ra2()}).
#'   Must contain \code{metrics} (a matrix with a row for the selected model)
#'   and \code{model_best} (the row index of that model). Also used by
#'   \code{fit_saz_reduced_models()} to reuse the stage-1 transformed design.
#' @param xi Character scalar; focal variable name.
#' @param power_best Named numeric vector of best powers selected in stage 1
#'   for \code{xi}. Used only to determine the FP degree for stage 2 metric
#'   computation.
#' @param y,weights,offset,family,family_string,method,strata,nocenter,control,
#'   rownames,has_offset Passed through to \code{fit_saz_reduced_models()}.
#' @param n_obs Numeric; number of observations (or events, for Cox models).
#' @param criterion,alpha,ftest Passed through to
#'   \code{compute_saz_stage2_decision()}. `criterion` must be `"pvalue"`;
#'   AIC/BIC SAZ selection is joint and never calls this helper. Stage 2 uses \code{alpha} for
#'   component-removal tests; \code{select} is used only for Stage-1 variable
#'   inclusion.
#' @param spike_decision Named numeric vector of current spike decisions. This
#'   function only updates \code{spike_decision[[xi]]}; all other entries are
#'   returned unchanged.
#' @param verbose Logical. If \code{TRUE}, prints the stage-2 selection table
#'   via \code{print_mfp_step(..., stage2 = TRUE)}.
#'
#' @return A list with:
#'   * \code{spike_decision}: the input vector with \code{spike_decision[[xi]]}
#'     updated to the stage-2 decision (\code{1}, \code{2}, or \code{3}).
#'   * \code{spike_metrics}: a list with \code{metrics} (matrix of stage-2
#'     model metrics, row-labeled), \code{spike_decision} (the scalar decision
#'     for \code{xi}), and \code{pvalue}. This mirrors the object historically
#'     attached to \code{fit1$spike_metrics} and is returned directly so
#'     callers and tests do not need to inspect \code{fit1} to retrieve it.
#'
#' @keywords internal
#' @noRd
evaluate_saz_stage2 <- function(fit1,
                                xi,
                                power_best,
                                y,
                                weights,
                                offset,
                                family,
                                family_string,
                                method,
                                strata,
                                nocenter,
                                control,
                                rownames,
                                has_offset,
                                n_obs,
                                criterion,
                                alpha,
                                ftest,
                                spike_decision,
                                verbose,
                                fitter = "base") {

  if (!identical(tolower(criterion), "pvalue")) {
    stop(
      "Internal error: SAZ Stage 2 is available only for p-value selection.",
      call. = FALSE
    )
  }

  # evaluate_saz_stage2() is the single entry point find_best_fp_step() calls
  # for SAZ stage 2 (see its own documentation): it orchestrates fitting the
  # two reduced models, scoring all three, deciding which to keep, and
  # building the printable summary - stitching together the four smaller
  # functions above into the one operation callers actually need.

  # Step 1: Fit the two reduced SAZ candidate models. Model 1 (full continuous +
  # binary) is already available via fit1 and is not refitted here.
  models <- fit_saz_reduced_models(
    stage1_selection = fit1,
    xi = xi,
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
    calculate_gaussian_deviance = isTRUE(ftest)
  )

  # Step 2: Score all three models on a comparable basis (log-likelihood,
  # AIC, BIC, etc., with the correct df_additional penalty for each).
  metrics <- compute_saz_stage2_metrics(
    fit1 = models$fit1,
    fit2 = models$fit2,
    fit3 = models$fit3,
    n_obs = n_obs,
    power_best = power_best
  )

  # Step 3: Decide whether both components, continuous-only, or binary-only
  # is needed, and record that decision for xi (all other variables' entries
  # in spike_decision are left untouched).
  decision <- compute_saz_stage2_decision(metrics, criterion, alpha, n_obs, ftest)
  spike_decision[xi] <- decision$decision

  # Step 4: Build the printable stage-2 metrics table, matching model labels
  # from stage 1, so the table reads consistently with what stage 1 already
  # printed for xi (e.g. "FP2 + Binary" rather than an unrelated label).
  f_names <- rownames(fit1$metrics)[fit1$model_best]
  spike_metrics_mat <- do.call(rbind, metrics)

  # Stage 1's row label for a SAZ variable's selected model is built as
  # "<continuous label> + Binary" (see select_ra2()/select_ic() and friends,
  # where has_binary_xi appends " + Binary" to model names). Splitting on
  # " + " recovers the two component names stage 2 needs to label its own
  # continuous-only and binary-only rows with matching terminology, e.g.
  # "FP2 + Binary" splits into "FP2" and "Binary".
  stage2_names <- strsplit(f_names, " \\+ ")[[1]]

  # For df = 1, stage 1 can be a linear SAZ model whose label has only one
  # component name. Ensure the full-model label always shows an explicit
  # binary component, so the three stage-2 rows read as: full model,
  # continuous/linear component only, binary component only.
  if (length(stage2_names) == 1L) {
    f_names <- paste0(stage2_names, " + Binary")
    stage2_names <- c(stage2_names, "Binary")
  }

  # Row order matches compute_saz_stage2_metrics()'s list order
  # (metrics1, metrics2, metrics3): full model first, then continuous-only,
  # then binary-only.
  rownames(spike_metrics_mat) <- c(f_names, stage2_names)

  spike_metrics <- list(
    metrics = spike_metrics_mat,
    spike_decision = spike_decision[xi],
    pvalue = decision$pvalue
  )

  # Printing (if requested) is done here, immediately after the decision is
  # finalized, rather than deferred to the caller, because the caller
  # (find_best_fp_step()) doesn't otherwise have a natural point at which to
  # print stage-2-specific output distinct from the stage-1 printout it
  # already produced.
  if (verbose) {
    fit1$spike_metrics <- spike_metrics
    print_mfp_step(xi = xi, criterion = criterion, fit = fit1, stage2 = TRUE)
  }

  list(
    spike_decision = spike_decision,
    spike_metrics = spike_metrics
  )
}
