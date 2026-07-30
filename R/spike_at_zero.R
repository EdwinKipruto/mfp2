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
#' Maps SAZ decision codes to existing print, plot, and verbose-output labels.
#' This preserves current user-facing text while keeping the mapping in one
#' place.
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
    out[decision_int == saz_decision_codes[["cont_binary"]]] <- "cont + binary"
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
      return("spike: zero indicator only")
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
        "spike: both components ",
        continuous_label,
        " + zero indicator"
      ))
    }

    return(full_label)
  }

  # continuous_only (decision 2): only the continuous FP/linear/ACD component
  # survives; the binary zero indicator is dropped, so binary_label is never
  # referenced in this branch.
  if (decision_int == saz_decision_codes[["continuous_only"]]) {
    if (style == "plot_title") {
      return(paste0("spike: positive-part only ", continuous_label))
    }

    return(continuous_label)
  }

  # decision_int matched none of the three known codes (e.g. corrupted data);
  # fall back to the generic placeholder rather than erroring, since this
  # function is used in display/labeling contexts, not validation contexts.
  unknown
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
#'   in \code{spike} are examined. For variables undergoing spike-at-zero
#'   eligibility checks, nonpositive values should already have been temporarily
#'   recoded to zero by \code{fit_mfp()}.
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
#' @param min_saz_component_prop Numeric in \eqn{(0, 0.5)}. Minimum required
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
#'     a standard continuous predictor -- no binary indicator, no zero recoding.}
#'   \item{User specified \code{spike} and \code{zero = TRUE}}{After reset:
#'     \code{spike = FALSE}, \code{catzero = FALSE}, \code{zero = TRUE}.
#'     Nonpositive values are still recoded to zero before FP transformation,
#'     but no binary indicator is added and the SAZ algorithm is not run.}
#'   \item{User specified \code{spike} and \code{catzero = TRUE}}{After reset:
#'     \code{spike = FALSE}, \code{catzero = TRUE}, \code{zero = TRUE}. The
#'     binary structural-zero indicator \code{I(x == 0)} on the recoded scale,
#'     equivalent to \code{I(original x <= 0)}, is still added as explicitly
#'     requested via \code{catzero}, but the SAZ selection algorithm is not run.}
#' }
#'
#' @details
#' Two conditions cause a variable's spike indicator to be reset:
#' \enumerate{
#'   \item \strong{Insufficient component representation} -- either the
#'     structural-zero proportion or the positive-observation proportion is
#'     below \code{min_saz_component_prop}.
#'   \item \strong{Binary variable} -- exactly two unique finite values. The
#'     positive part would contain a single unique value, making FP
#'     transformation degenerate.
#' }
#' A warning is issued for each reset reason, identifying the affected variables.
#'
#' @keywords internal
#' @noRd
reset_spike <- function(x, spike, user_catzero, user_zero,
                        min_saz_component_prop = 0.10) {

  # Early exit: no spike variables to evaluate. In this case there was no
  # spike-implied cascade to preserve, so return the user's original zero and
  # catzero settings unchanged.
  if (!any(spike)) {
    return(list(spike = spike, catzero = user_catzero, zero = user_zero))
  }

  names_spike <- names(spike)[spike]

  # Component proportions for each requested spike-at-zero variable.
  # At this point fit_mfp() has already recoded nonpositive values to zero for
  # zero/catzero/spike variables, so x == 0 represents the structural-zero group
  # and x > 0 represents the positive continuous component.
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

  # Binary variables are not eligible for spike-at-zero modelling. A binary
  # covariate already represents a two-level effect, so adding a separate zero
  # spike indicator would be redundant or non-identifiable.
  is_binary <- vapply(
    names_spike,
    function(v) {
      length(unique(x[is.finite(x[, v]), v])) == 2L
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
      prop_zero < min_saz_component_prop |
      prop_positive < min_saz_component_prop
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
            "all observed values are positive; zero proportion is 0 < `min_saz_component_prop = %s`",
            min_saz_component_prop
          ))
        }

        if (np == 0L && nz > 0L) {
          return(sprintf(
            "all observed values are structural zeros; positive observation proportion is 0 < `min_saz_component_prop = %s`",
            min_saz_component_prop
          ))
        }

        if (is.na(p0) || p0 < min_saz_component_prop) {
          return(sprintf(
            "zero proportion is %.4g < `min_saz_component_prop = %s`",
            p0,
            min_saz_component_prop
          ))
        }

        if (is.na(pp) || pp < min_saz_component_prop) {
          return(sprintf(
            "positive observation proportion is %.4g < `min_saz_component_prop = %s`",
            pp,
            min_saz_component_prop
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
#' @param min_saz_component_prop Numeric in `(0, 0.5)`. Minimum required
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
                                    min_saz_component_prop = 0.10) {
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

  # Apply a temporary cascade only for eligibility checking: spike implies
  # catzero implies zero, exactly mirroring the cascade fit_mfp() applies
  # later, so that the component proportions checked below are computed under
  # the same zero-recoding a variable would actually receive if judged
  # eligible.
  catzero_for_spike <- catzero
  zero_for_spike <- zero

  catzero_for_spike[spike] <- TRUE
  zero_for_spike[catzero_for_spike] <- TRUE

  # Build temporary data for the eligibility check. This does not modify the
  # original x used for fitting or preprocessing.
  x_for_spike <- x

  zero_vars <- intersect(
    names(zero_for_spike)[zero_for_spike],
    colnames(x_for_spike)
  )

  # Temporarily recode nonpositive values to zero for every zero-cascade
  # variable, matching what fit_mfp() will eventually do to the real x. This
  # is what lets reset_spike() (called next) count the zero vs. positive
  # component using x == 0 / x > 0 directly, as documented there.
  for (v in zero_vars) {
    x_for_spike[x_for_spike[, v] <= 0, v] <- 0
  }

  # Delegate the actual eligibility rule (component proportions, binary
  # check, cascade restoration for ineligible variables) to reset_spike();
  # this function's job is only to prepare the right temporary inputs for it.
  reset_spike(
    x = x_for_spike,
    spike = spike,
    user_catzero = user_catzero,
    user_zero = user_zero,
    min_saz_component_prop = min_saz_component_prop
  )
}

#' Cap FP Degrees of Freedom for Spike-at-Zero Positive Components
#'
#' Internal helper used by \code{fit_mfp()} after \code{reset_spike()}.
#'
#' For ordinary variables, \code{assign_df()} caps the maximum FP degrees of
#' freedom according to the number of distinct values in the full variable. For
#' spike-at-zero variables, the continuous FP part is fitted to the positive
#' component, while zero/nonpositive observations are represented structurally
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
#' @param x Numeric matrix or data frame after temporary zero recoding for spike
#'   eligibility. For spike variables, nonpositive values should already be
#'   represented as zero.
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


#' Fit Reduced Models for SAZ Stage 2
#'
#' Fits the two reduced candidate models required for stage 2 of the
#' spike-at-zero (SAZ) algorithm.
#'
#' Stage 2 compares three nested models for the current variable \code{xi}:
#'
#' \describe{
#'   \item{Model 1}{
#'     FP/ACD(\code{xi}) + binary zero indicator(\code{xi}) + adjustment
#'     variables. This model is already fitted in stage 1.
#'   }
#'   \item{Model 2}{
#'     FP/ACD(\code{xi}) only + adjustment variables.
#'   }
#'   \item{Model 3}{
#'     Binary zero indicator(\code{xi}) only + adjustment variables.
#'   }
#' }
#'
#' This function fits only Models 2 and 3. Model 1 is represented by the
#' already-selected stage-1 object passed via \code{stage1_selection}; it is not
#' refitted here.
#'
#' The function reuses the already transformed stage-1 \code{xi} design matrix
#' and adjustment matrix stored in
#' \code{stage1_selection$current_adj_params[[xi]]}. The \code{xi} design is
#' split by the \code{"catzero"} column: non-\code{"catzero"} columns form the
#' continuous FP/ACD component, and \code{"catzero"} forms the binary
#' structural-zero component.
#'
#' @param stage1_selection Stage 1 selection object, as returned by one of the
#'   \code{select_*()} functions (e.g. \code{select_ra2()}, \code{select_ic()}).
#'   Used both as the Model 1 fit (returned unchanged as \code{fit1}) and to
#'   reuse the already-transformed \code{xi} design/adjustment matrices stored
#'   in \code{stage1_selection$current_adj_params[[xi]]}, avoiding a full
#'   retransformation of \code{xi} and its adjustment set.
#' @param xi Character scalar; focal variable name.
#' @param y,weights,offset,family,family_string,method,strata,nocenter,control,
#'   rownames,has_offset Passed through to \code{fit_model()} for fitting
#'   Model 2 and Model 3.
#' @param calculate_gaussian_deviance Logical. If `TRUE`, compute the scalar
#'   Gaussian deviance required for stage-2 F-tests.
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{fit1}: alias for \code{stage1_selection} (Model 1; not
#'       refitted here, included for a uniform \code{fit1}/\code{fit2}/
#'       \code{fit3} naming convention downstream).
#'     \item \code{fit2}: fitted Model 2 (continuous FP/ACD component only,
#'       plus adjustment variables).
#'     \item \code{fit3}: fitted Model 3 (binary zero-indicator only, plus
#'       adjustment variables).
#'     \item \code{x}: a list with \code{model2} and \code{model3}, the design
#'       matrices used to fit \code{fit2} and \code{fit3} respectively.
#'     \item \code{data_xi}: the reused stage-1 transformed design matrix for
#'       \code{xi} (continuous FP/ACD column(s) plus the \code{"catzero"}
#'       column).
#'     \item \code{adjustment_matrix}: the reused stage-1 adjustment matrix, or
#'       \code{NULL} if there were no adjustment variables or the stage-1
#'       adjustment matrix had zero columns.
#'   }
#'
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

  # Split xi's stage-1 design into its two components: the continuous
  # FP/linear/ACD column(s), and the single "catzero" binary indicator
  # column. Model 2 uses only the former, Model 3 uses only the latter (see
  # the function documentation for the three-model comparison).
  xi_continuous <- data_xi[, data_xi_colnames != "catzero", drop = FALSE]
  xi_binary <- data_xi[, "catzero", drop = FALSE]

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

  # Step 3: Assemble the design matrices for Model 2 and Model 3 -----------
  if (is.null(adjustment_matrix)) {
    x_fit2 <- xi_continuous
    x_fit3 <- xi_binary
  } else {
    x_fit2 <- cbind(xi_continuous, adjustment_matrix)
    x_fit3 <- cbind(xi_binary, adjustment_matrix)
  }

  # Step 4: Fit both reduced models with the shared fitting arguments -------
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

  fit2 <- do.call(fit_model, c(list(x = x_fit2), fit_args))
  fit3 <- do.call(fit_model, c(list(x = x_fit3), fit_args))

  # Step 5: Return everything compute_saz_stage2_metrics()/decision-making
  # downstream will need: the three fits (fit1 is just an alias for the
  # unmodified stage1_selection, included for a uniform fit1/fit2/fit3
  # naming convention), the two design matrices actually used, and the
  # reused stage-1 pieces (data_xi, adjustment_matrix) in case a caller needs
  # to inspect them directly.
  list(
    fit1 = stage1_selection,
    fit2 = fit2,
    fit3 = fit3,
    x = list(
      model2 = x_fit2,
      model3 = x_fit3
    ),
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
#' @param select Numeric significance threshold used when
#'   \code{criterion = "pvalue"}.
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
#'   \code{p_drop_binary < select}, \code{p_drop_continuous < select} \tab
#'   Both reductions are significantly worse than Model 1 \tab
#'   Continuous FP/linear/ACD + binary zero indicator \tab
#'   Model 1 \cr
#'   \code{p_drop_binary < select}, \code{p_drop_continuous >= select} \tab
#'   Removing the binary component is harmful; removing the continuous component is acceptable \tab
#'   Binary zero indicator only \tab
#'   Model 3 \cr
#'   \code{p_drop_binary >= select}, \code{p_drop_continuous < select} \tab
#'   Removing the continuous component is harmful; removing the binary component is acceptable \tab
#'   Continuous FP/linear/ACD only \tab
#'   Model 2 \cr
#'   \code{p_drop_binary >= select}, \code{p_drop_continuous >= select} \tab
#'   Neither reduction is significantly worse than Model 1 \tab
#'   Better-fitting reduced component \tab
#'   Model 2 if \code{logLik(Model 2) > logLik(Model 3)}, otherwise Model 3 \cr
#' }
#'
#' If \code{ftest = FALSE}, the nested comparisons use likelihood-ratio tests.
#' If \code{ftest = TRUE}, the nested comparisons use F-tests. For
#' information-criterion based selection, no hypothesis tests are used; the
#' model with the smallest requested information criterion is selected.
#'
#' @return A list with two elements:
#'   \itemize{
#'     \item \code{decision}: Integer indicating the selected model:
#'       \code{1} = both components, \code{2} = continuous FP/linear/ACD only,
#'       and \code{3} = binary zero-indicator only.
#'     \item \code{pvalue}: Named numeric vector containing
#'       \code{p_drop_binary} and \code{p_drop_continuous} when
#'       \code{criterion = "pvalue"}; otherwise \code{NA}.
#'   }
#'
#' @keywords internal
#' @noRd
compute_saz_stage2_decision <- function(metrics,
                                        criterion,
                                        select,
                                        n_obs,
                                        ftest = FALSE) {
  criterion <- tolower(criterion)

  if (!criterion %in% c("pvalue", "aic", "bic")) {
    stop(
      "! criterion must be one of 'pvalue', 'aic', or 'bic'.",
      call. = FALSE
    )
  }

  if (criterion == "pvalue") {
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

    # This mirrors the decision table in the function documentation above:
    # a small p-value means dropping that component makes the model
    # significantly worse, i.e. that component must be *kept*. So:
    #   both p-values small  -> neither component can be safely dropped -> Model 1
    #   only p_drop_binary small -> the binary component is needed but the
    #     continuous one isn't -> keep binary only -> Model 3
    #   only p_drop_continuous small -> the continuous component is needed but
    #     the binary one isn't -> keep continuous only -> Model 2
    #   neither small -> both components are individually dispensable; break
    #     the tie by whichever reduced model actually fits better (higher
    #     log-likelihood), rather than defaulting to either one arbitrarily.
    decision <- if (p_drop_binary < select && p_drop_continuous < select) {
      saz_decision_codes[["cont_binary"]]  # Model 1: both components
    } else if (p_drop_binary < select && p_drop_continuous >= select) {
      saz_decision_codes[["binary_only"]]  # Model 3: binary zero-indicator only
    } else if (p_drop_binary >= select && p_drop_continuous < select) {
      saz_decision_codes[["continuous_only"]]  # Model 2: continuous FP/linear/ACD only
    } else {
      # If neither reduced model is significantly worse than the full model,
      # choose the better-fitting reduced model.
      if (metrics$metrics2["logl"] > metrics$metrics3["logl"]) {
        saz_decision_codes[["continuous_only"]]
      } else {
        saz_decision_codes[["binary_only"]]
      }
    }

    return(list(
      decision = decision,
      pvalue = c(
        p_drop_binary = p_drop_binary,
        p_drop_continuous = p_drop_continuous
      )
    ))
  }

  # For AIC/BIC there is no hypothesis test at all: simply compute the
  # criterion for all three models and take whichever is smallest. Naming
  # the values by saz_decision_codes' names (in the same order the values are
  # assembled: metrics1/2/3 <-> cont_binary/continuous_only/binary_only) lets
  # which.min()'s returned name be looked straight back up in
  # saz_decision_codes to get the corresponding integer decision code.
  if (criterion == "aic") {
    aic_values <- c(
      metrics$metrics1["aic"],
      metrics$metrics2["aic"],
      metrics$metrics3["aic"]
    )
    names(aic_values) <- names(saz_decision_codes)

    decision <- saz_decision_codes[[names(which.min(aic_values))]]

    return(list(decision = decision, pvalue = NA_real_))
  }

  # criterion == "bic": identical logic to the AIC branch above, just using
  # the BIC column instead.
  bic_values <- c(
    metrics$metrics1["bic"],
    metrics$metrics2["bic"],
    metrics$metrics3["bic"]
  )
  names(bic_values) <- names(saz_decision_codes)

  decision <- saz_decision_codes[[names(which.min(bic_values))]]

  list(decision = decision, pvalue = NA_real_)
}


#' Evaluate Stage 2 of the Spike-at-Zero (SAZ) Algorithm for One Variable
#'
#' Internal helper used by \code{find_best_fp_step()}.
#'
#' Stage 2 of the SAZ algorithm is evaluated only when stage 1 selected a
#' non-null functional form for a spike-at-zero variable \code{xi}. This
#' function fits the two reduced candidate models (continuous-only and
#' binary-only), compares them against the full continuous + binary model
#' already selected in stage 1, updates \code{spike_decision[[xi]]}
#' accordingly, and builds the printable stage-2 metrics table.
#'
#' @param fit1 Stage 1 selection object, as returned by one of the
#'   \code{select_*()} functions (e.g. \code{select_ra2()}, \code{select_ic()}).
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
#' @param criterion,select,ftest Passed through to
#'   \code{compute_saz_stage2_decision()}.
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
                                select,
                                ftest,
                                spike_decision,
                                verbose,
                                fitter = "base") {

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
  decision <- compute_saz_stage2_decision(metrics, criterion, select, n_obs, ftest)
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