#' Reset Spike-at-Zero Indicators and Undo Cascade for Ineligible Variables
#'
#' Evaluates whether each variable flagged as a spike-at-zero (\code{spike})
#' is eligible for the SAZ algorithm based on its zero proportion and
#' cardinality. For ineligible variables, \code{spike} is reset to
#' \code{FALSE} and the \code{catzero} and \code{zero} flags are restored to
#' the user's original pre-cascade values. This correctly undoes the implicit
#' cascade (\code{catzero[spike] <- TRUE}, \code{zero[catzero] <- TRUE}) that
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
#' \code{zero} remain \code{TRUE} for
#' ineligible variables even though the user never asked for them. This
#' function corrects that by accepting the pre-cascade user values and
#' restoring them for any variable where \code{spike} is reset.
#'
#' @param x A numeric matrix or data frame with column names. Only columns
#'   named in \code{spike} are examined. Assumes that the values for variables
#'   to undergo spike have zeros and positives values only.
#' @param spike A named logical vector. \code{TRUE} indicates the column was
#'   flagged as a spike-at-zero variable. Should be the post-cascade version
#'   (after \code{catzero[spike] <- TRUE} has been applied in
#'   \code{fit_mfp()}).
#' @param user_catzero A named logical vector of the user's original
#'   \code{catzero} specification \strong{before} the cascade
#'   \code{catzero[spike] <- TRUE} was applied. Must have the same names as
#'   \code{spike}.
#' @param user_zero A named logical vector of the user's original \code{zero}
#'   specification \strong{before} the cascade \code{zero[catzero] <- TRUE}
#'   was applied. Must have the same names as \code{spike}.
#' @param min_prop Numeric in \eqn{(0, 0.5)}. Minimum proportion of zeros
#'   required to retain the spike indicator. If the observed zero proportion
#'   is below \code{min_prop}, the variable has too few zeros for a meaningful
#'   spike model. Default \code{0.05}.
#' @param max_prop Numeric in \eqn{(0.5, 1)}. Maximum proportion of zeros
#'   allowed to retain the spike indicator. If the observed zero proportion
#'   exceeds \code{max_prop}, the positive part is too sparse for reliable FP
#'   fitting. Default \code{0.95}.
#'
#' @return A named list with three elements, each a named logical vector of
#'   the same length as \code{spike}:
#' \describe{
#'   \item{\code{spike}}{Updated spike vector. Entries for ineligible
#'     variables are set to \code{FALSE}; all other entries are unchanged.}
#'   \item{\code{catzero}}{Updated catzero vector. For ineligible variables
#'     (where \code{spike} was just reset), restored to \code{user_catzero}.
#'     All other entries retain their post-cascade values, i.e. \code{TRUE}
#'     for variables that remain eligible spike variables.}
#'   \item{\code{zero}}{Updated zero vector. For ineligible variables,
#'     restored to \code{user_zero}. All other entries retain their
#'     post-cascade values.}
#' }
#'
#' @section Cascade restoration examples:
#' For a variable whose \code{spike} is reset, the outcome depends on what
#' the user originally specified:
#' \describe{
#'   \item{User specified \code{spike} only (not \code{catzero} or
#'     \code{zero})}{After reset: \code{spike = FALSE},
#'     \code{catzero = FALSE}, \code{zero = FALSE}. The variable is treated
#'     as a standard continuous predictor -- no binary indicator, no zero
#'     recoding.}
#'   \item{User specified \code{spike} and \code{zero = TRUE}}{After reset:
#'     \code{spike = FALSE}, \code{catzero = FALSE}, \code{zero = TRUE}.
#'     Non-positive values are still recoded to zero before FP transformation,
#'     but no binary indicator is added and the SAZ algorithm is not run.}
#'   \item{User specified \code{spike} and \code{catzero = TRUE}}{After
#'     reset: \code{spike = FALSE}, \code{catzero = TRUE},
#'     \code{zero = TRUE}. The binary structural-zero indicator
#'     \code{I(x == 0)} on the recoded scale, equivalent to
#'     \code{I(original x <= 0)}, is still added (as explicitly requested 
#'     via \code{catzero}), but the SAZ
#'     selection algorithm is not run.}
#' }
#'
#' @details
#' Three conditions cause a variable's spike indicator to be reset:
#' \enumerate{
#'   \item \strong{Too few zeros} -- zero proportion below \code{min_prop}.
#'     There are not enough zero observations to reliably model a spike at
#'     zero or to estimate its contribution separately from the FP function.
#'   \item \strong{Too many zeros} -- zero proportion above \code{max_prop}.
#'     The positive part has too few observations for reliable FP power
#'     selection and fitting.
#'   \item \strong{Binary variable} -- exactly two unique finite values.
#'     The positive part would contain a single unique value, making FP
#'     transformation degenerate.
#' }
#' A warning is issued for each reset, identifying the affected variables and
#' the reason. Variables reset for both proportion and binary reasons receive
#' two separate warnings.
#'
#' @keywords internal
#' @noRd
reset_spike <- function(x, spike, user_catzero, user_zero,
                        min_prop = 0.05, max_prop = 0.95) {
  
  # Early exit: no spike variables to evaluate. In this case there was no
  # spike-implied cascade to preserve, so return the user's original zero and
  # catzero settings unchanged.
  if (!any(spike)) {
    return(list(spike = spike, catzero = user_catzero, zero = user_zero))
  }
  
  names_spike <- names(spike)[spike]
  
  # Proportion of structural zeros for each requested spike-at-zero variable.
  # At this point fit_mfp() has already recoded nonpositive values to zero for
  # zero/catzero/spike variables, so x == 0 represents the zero/nonpositive
  # group used by the spike-at-zero machinery.
  prop_zero <- colMeans(x[, names_spike, drop = FALSE] == 0, na.rm = TRUE)
  
  # Binary variables are not eligible for spike-at-zero modelling. A binary
  # covariate already represents a two-level effect, so adding a separate zero
  # spike indicator would be redundant or non-identifiable.
  is_binary <- vapply(names_spike, function(v) {
    length(unique(x[is.finite(x[, v]), v])) == 2L
  }, logical(1L))
  
  # Identify variables whose spike option must be reset by each reason.
  to_reset_prop   <- names_spike[
    prop_zero <= 0 | prop_zero >= 1 |
      prop_zero < min_prop | prop_zero > max_prop
  ]
  
  to_reset_binary <- names_spike[is_binary]
  to_reset        <- union(to_reset_prop, to_reset_binary)
  
  # Warn about proportion-based resets.
  if (length(to_reset_prop) > 0L) {
    warning(
      "The spike-at-zero option has been reset for the following variable(s) ",
      "because the zero proportion is outside [", min_prop, ", ", max_prop, "]: ",
      paste(to_reset_prop, collapse = ", "), ".",
      call. = FALSE
    )
  }
  
  # Warn about binary-variable resets.
  if (length(to_reset_binary) > 0L) {
    warning(
      "The spike-at-zero option has been reset for the following variable(s) ",
      "because they are binary (exactly two unique finite values): ",
      paste(to_reset_binary, collapse = ", "), ".",
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


#' Fit Reduced Models for SAZ Stage 2
#'
#' Fits the two reduced candidate models required for stage 2 of the
#' spike-at-zero (SAZ) algorithm. Stage 2 compares three nested models for
#' the current variable \code{xi}:
#'
#' \describe{
#'   \item{Model 1}{FP/ACD(\code{xi}) + binary zero indicator(\code{xi})
#'     + adjustment variables. Already fitted in stage 1 by \code{select_fct()}
#'     via \code{find_best_fpm_step()} or \code{fit_linear_step()}.}
#'   \item{Model 2}{FP/ACD(\code{xi}) only + adjustment variables.
#'     Model 1 minus the binary indicator column.}
#'   \item{Model 3}{Binary zero indicator(\code{xi}) only + adjustment variables.
#'     Model 1 minus the FP/ACD columns.}
#' }
#'
#' This function fits only Models 2 and 3. It reuses the already transformed
#' stage-1 design data stored in \code{current_adj_params[[xi]]$data_xi}; it
#' does not call FP or ACD transformation helpers. The transformed \code{xi}
#' matrix is split by the \code{"catzero"} column, where non-\code{"catzero"}
#' columns form the continuous FP/ACD component and the \code{"catzero"}
#' column forms the binary structural-zero component.
#'
#' @section Adjustment matrix invariant:
#' The adjustment matrix must be \strong{identical} across Models 1, 2, and 3.
#' This function guarantees that by reusing the adjustment matrix that
#' \code{transform_data_step()} built during the stage-1 fit, passed in as
#' \code{current_adj_params}. It does not reconstruct adjustment variables.
#' This prevents loss of catzero binary indicators, spike-decision routing,
#' or ACD structure for adjustment variables.
#'
#' @param xi Character string, the name of the current variable of interest.
#' @param y Response vector or \code{Surv} object.
#' @param weights Observation weight vector.
#' @param offset Offset vector.
#' @param family Family specification (character or function).
#' @param method Tie-handling method for Cox models.
#' @param strata Stratification factor for Cox models.
#' @param nocenter Numeric vector for Cox centering control.
#' @param control List of fitting control parameters.
#' @param rownames Row-name vector for Cox models.
#' @param has_offset Logical; whether an offset was originally specified.
#' @param current_adj_params Named list of stage-1 design data and adjustment
#'   parameters from the winning stage-1 model. The element
#'   \code{current_adj_params[[xi]]$data_xi} provides the already transformed
#'   \code{xi} matrix, and \code{current_adj_params[[xi]]$data_adj} provides
#'   the adjustment matrix reused across all three models.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{\code{fit2}}{Fitted model object for Model 2 (continuous-only
#'       \code{xi} + adjustment), as returned by \code{fit_model()}.}
#'     \item{\code{fit3}}{Fitted model object for Model 3 (binary-only
#'       \code{xi} + adjustment), as returned by \code{fit_model()}.}
#'   }
#'
#' @seealso \code{compute_metrics()}, \code{compute_decision()},
#'   \code{find_best_fp_step()}, \code{transform_data_step()}
#'
#' @keywords internal
#' @noRd
fit_saz_reduced_models <- function(xi, y, weights, offset, family,
                                   method, strata, nocenter, control,
                                   rownames, has_offset,
                                   current_adj_params) {
  data_xi <- current_adj_params[[xi]]$data_xi
  if (is.null(data_xi)) {
    stop(
      "Internal error: current_adj_params[[", xi,
      "]]$data_xi is NULL in fit_saz_reduced_models(). ",
      "SAZ stage 2 requires the transformed xi data from stage 1.",
      call. = FALSE
    )
  }
  data_xi <- as.matrix(data_xi)
  if (!"catzero" %in% colnames(data_xi)) {
    stop(
      "Internal error: current_adj_params[[", xi,
      "]]$data_xi must contain a 'catzero' column in fit_saz_reduced_models().",
      call. = FALSE
    )
  }
  
  xi_continuous <- data_xi[, colnames(data_xi) != "catzero", drop = FALSE]
  xi_binary <- data_xi[, "catzero", drop = FALSE]
  
  adjustment_matrix <- current_adj_params[[xi]]$data_adj
  if (!is.null(adjustment_matrix) && ncol(adjustment_matrix) == 0L) {
    adjustment_matrix <- NULL
  }
  
  x_fit2 <- cbind(xi_continuous, adjustment_matrix)
  x_fit3 <- cbind(xi_binary, adjustment_matrix)
  
  fit_args <- list(
    y          = y,
    family     = family,
    weights    = weights,
    offset     = offset,
    method     = method,
    strata     = strata,
    nocenter   = nocenter,
    has_offset = has_offset,
    control    = control,
    rownames   = rownames
  )
  
  fit2 <- do.call(fit_model, c(list(x = x_fit2), fit_args))
  fit3 <- do.call(fit_model, c(list(x = x_fit3), fit_args))
  
  list(fit2 = fit2, fit3 = fit3)
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
compute_metrics <- function(fit1, fit2, fit3, n_obs, power_best) {
  
  # ACD can produce NA like c(NA,1) so degree will reduce to 1 and in this case
  # additional parameters = 0 since the power = 1, see mfpa paper
  power_best <- power_best[!is.na(power_best)]
  degree <- length(power_best[!is.na(power_best)])
  
  # Extract Model 1 metrics safely
  if (is.null(fit1$metrics) || is.null(fit1$model_best)) {
    stop("fit1 must contain 'metrics' and 'model_best' elements.")
  }
  
  if (fit1$model_best > nrow(fit1$metrics) || fit1$model_best < 1) {
    stop("fit1$model_best is out of bounds for fit1$metrics.")
  }
  
  metrics1 <- fit1$metrics[fit1$model_best, ]
  
  # Compute Model 2 metrics with degree adjustment
  deg2 <- if (degree == 1 && power_best == 1) 0 else degree
  metrics2 <- tryCatch(
    calculate_model_metrics(fit2, n_obs, deg2),
    error = function(e) stop("Failed to compute metrics for fit2: ", e$message)
  )
  
  # Compute Model 3 metrics
  metrics3 <- tryCatch(
    calculate_model_metrics(fit3, n_obs),
    error = function(e) stop("Failed to compute metrics for fit3: ", e$message)
  )
  return(list(metrics1 = metrics1, metrics2 = metrics2, metrics3 = metrics3))
}

#' Compute Model Selection Decision
#'
#' This function compares competing regression models based on a specified 
#' selection criterion (`pvalue`, `aic`, or `bic`) and returns a decision 
#' about which model to retain. It is used in stage 2 of the 
#' spike-at-zero (SAZ) algorithm, where a decision is made between alternative 
#' reduced or complex models after stage 1 selection.
#'
#' @param metrics A list containing model fit statistics for three candidate models:
#'   * `metrics1`: statistics for the most complex model selected in stage 1 of the SAZ algorithm, typically including both FPm/linear and binary terms.
#'   * `metrics2`: statistics for a reduced model, e.g., FPm/linear only.
#'   * `metrics3`: statistics for an alternative reduced model, e.g., binary only.
#' Each element must provide named values for:
#'   * `logl`: log-likelihood,
#'   * `df`: model degrees of freedom,
#'   * `aic`: Akaike Information Criterion,
#'   * `bic`: Bayesian Information Criterion,
#'   * `deviance_gaussian`: Gaussian deviance (required if `ftest = TRUE`),
#'   * `df_resid`: residual degrees of freedom (required if `ftest = TRUE`).
#' @param criterion A string specifying the selection criterion. 
#'   One of `"pvalue"`, `"aic"`, or `"bic"`.
#' @param select Numeric threshold for significance testing (used only when 
#'   `criterion = "pvalue"`).
#' @param n_obs Integer. Number of observations in the data (used for F-tests).
#' @param ftest Logical. If `TRUE`, use F-tests instead of likelihood ratio 
#'   tests when `criterion = "pvalue"`.
#'
#' @details
#' - When `criterion = "pvalue"`, the function compares:
#'   * Model 1 (both FPm/linear and binary terms),
#'   * Model 2 (FPm/linear only),
#'   * Model 3 (binary only).
#'   Significance is determined by comparing complex vs. reduced models 
#'   using either likelihood ratio tests or F-tests.
#'
#' - When `criterion = "aic"` or `"bic"`, the function selects the model 
#'   with the lowest information criterion value.
#'
#' @return A list with two elements:
#'   - `decision`: An integer indicating the selected model 
#'     (1 = both terms, 2 = FPm/linear only, 3 = binary only).
#'   - `pvalue`: A numeric vector of p-values from the pairwise tests, 
#'     or `NA` when `criterion` is `"aic"` or `"bic"`.
#'
#' @keywords internal
compute_decision <- function(metrics, criterion, select, n_obs, ftest = FALSE) {
  criterion <- tolower(criterion)
  
  if (criterion == "pvalue") {
    if (ftest) { 
      # Complex vs reduced model (FPm + Binary vs FPm)
      stats1  <- calculate_f_test(
        deviances = c(metrics$metrics2["deviance_gaussian"],
                      metrics$metrics1["deviance_gaussian"]),
        dfs_resid = c(metrics$metrics2["df_resid"],
                      metrics$metrics1["df_resid"]),
        n_obs = n_obs
      )
      # Complex vs reduced model (FPm + Binary vs Binary)
      stats2  <- calculate_f_test(
        deviances = c(metrics$metrics3["deviance_gaussian"],
                      metrics$metrics1["deviance_gaussian"]),
        dfs_resid = c(metrics$metrics3["df_resid"],
                      metrics$metrics1["df_resid"]),
        n_obs = n_obs
      )
      
      
    } else {
      
      # reduced model (FPm/Linear-mode 2) vs complex model (FPm/linear + binary-model 1)
      stats1 <- calculate_lr_test(
        c(metrics$metrics2["logl"], metrics$metrics1["logl"]),
        c(metrics$metrics2["df"], metrics$metrics1["df"])
      )
      
      # reduced model (binary only-mode 3) vs complex model (FPm/linear + binary-model 1)
      stats2 <- calculate_lr_test(
        c(metrics$metrics3["logl"], metrics$metrics1["logl"]),
        c(metrics$metrics3["df"], metrics$metrics1["df"])
      )
    }
    
    p12 <- stats1$pvalue
    p13 <- stats2$pvalue
    
    # Decision logic
    decision <- if (p12 <= select & p13 <= select) {
      1  # Model 1: both terms
    } else if (p12 <= select & p13 > select) {
      3  # Model 3: binary only
    } else if (p12 > select & p13 <= select) {
      2  # Model 2: FPm/linear only
    } else {
      # If both are not significant, use log-likelihood
      if (metrics$metrics2["logl"] > metrics$metrics3["logl"]) 2 else 3
    }
    
    # Return both decision and p-values
    return(list(
      decision = decision,
      pvalue = c(p12, p13)
    ))
    
  } else if (criterion == "aic") {
    decision <- which.min(c(metrics$metrics1["aic"], 
                            metrics$metrics2["aic"], 
                            metrics$metrics3["aic"]))
    return(list(decision = decision, pvalue = NA))
    
  } else { # criterion == "bic"
    decision <- which.min(c(metrics$metrics1["bic"], 
                            metrics$metrics2["bic"], 
                            metrics$metrics3["bic"]))
    return(list(decision = decision, pvalue = NA))
  }
}