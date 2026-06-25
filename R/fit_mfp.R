#' Function for fitting a model using the MFP, MFPA or spike-at-zero algorithm
#'
#' This internal function implements the Multivariable Fractional Polynomial (MFP),
#' MFP with Approximate Cumulative Distribution (MFPA), and spike-at-zero (SAZ)
#' algorithms. It is not exported and is intended to be called from \code{mfp2()}.
#' While most parameters are documented in \code{mfp2()}, their form may differ here.
#' The function does not perform argument checks and expects that all inputs
#' have been properly prepared by \code{mfp2()}.
#'
#' @param x A numeric matrix of dimensions nobs x nvars. Does not contain an
#'   intercept, but columns are already expanded into dummy variables as
#'   necessary. Data are assumed to be shifted and scaled.
#' @param y A vector for the response variable or a \code{Surv} object.
#' @param weights A vector of observation weights of length nobs.
#' @param offset A vector of length nobs of offsets.
#' @param cycles An integer representing the maximum number of iteration cycles
#'   during which FP powers for all predictors are updated.
#' @param scale A named numeric vector of length nvars of scaling factors. Not
#'   applied during power selection cycles, but used in step 4 to backscale
#'   \code{x} before the final model fit so that coefficients are on the
#'   \eqn{\phi(x + \text{shift})} scale.
#' @param shift A named numeric vector of length nvars of shift terms. Not
#'   applied (shifting was done upstream), but re-ordered to conform to
#'   \code{xorder} and stored in the returned object.
#' @param df A numeric vector of length nvars of degrees of freedom.
#' @param center A logical vector of length nvars indicating if variables are
#'   to be centred.
#' @param family Either a character string specifying the model family
#'   (e.g., \code{"gaussian"}, \code{"binomial"}, \code{"poisson"},
#'   \code{"cox"}) or a function that returns a GLM family object such as
#'   \code{stats::gaussian(link = "identity")}. For Cox models, only the
#'   character string \code{"cox"} is allowed.
#' @param family_string A character string representing the selected family,
#'   e.g., \code{"gaussian"}.
#' @param criterion A character string defining the criterion used to select
#'   variables and FP models of different degrees.
#' @param select A numeric vector of length nvars indicating significance levels
#'   used during MFP backfitting to decide whether each predictor is retained.
#' @param alpha A numeric vector of length nvars indicating significance levels
#'   for tests between FP models of different degrees.
#' @param keep A character vector with names of variables to be kept in the
#'   model regardless of selection criteria.
#' @param xorder A string determining the order of entry of the covariates into
#'   the model-selection algorithm.
#' @param powers A named list of numeric values setting the permitted FP powers
#'   for each covariate.
#' @param method A character string specifying the method for tie handling in
#'   Cox regression.
#' @param strata A factor of all possible combinations of stratification
#'   variables. Returned from \code{survival::strata()}.
#' @param nocenter A numeric vector with a list of values for fitting Cox
#'   models. See \code{survival::coxph()} for details.
#' @param acdx A logical vector of length nvars indicating which continuous
#'   variables should undergo the approximate cumulative distribution (ACD)
#'   transformation.
#' @param ftest Logical. If \code{TRUE} and \code{family = "gaussian"}, use an
#'   F-test rather than a chi-square likelihood-ratio test.
#' @param control A list with parameters for model fit. See
#'   \code{survival::coxph()} or \code{stats::glm()} for details.
#' @param zero A logical vector indicating which columns of \code{x} should
#'   treat non-positive values as zero before FP transformation. Must match
#'   the length and order of the columns of \code{x}.
#' @param catzero A logical vector indicating which columns of \code{x} should
#'   treat non-positive values as zero AND have a binary structural-zero
#'   indicator automatically created and included in the model. Internally,
#'   values \code{x <= 0} are first recoded to zero; the indicator is then
#'   computed as \code{I(x == 0)} on the recoded scale, equivalent to
#'   \code{I(original x <= 0)}. Must
#'   match the length and order of the columns of \code{x}.
#' @param spike A logical vector indicating which columns of \code{x} contain
#'   a spike at zero and should be assessed using the SAZ algorithm. Must
#'   match the length and order of the columns of \code{x}.
#' @param min_prop Numeric in \eqn{(0, 0.5)}. Minimum proportion of zeros
#'   required for the SAZ algorithm to be applied. Default \code{0.05}.
#' @param max_prop Numeric in \eqn{(0.5, 1)}. Maximum proportion of zeros
#'   allowed for the SAZ algorithm to be applied. Default \code{0.95}.
#' @param force_max_fp A logical vector of length nvars. If \code{TRUE} for a
#'   variable, forces selection of the most complex functional form at the
#'   degree specified by \code{df}, bypassing AIC/BIC comparison against
#'   simpler forms. Has no effect when \code{criterion = "pvalue"}.
#' @param has_offset logical indicating whether an offset was specified before
#' missing offsets were replaced by zeros internally.
#' @param verbose Logical. If \code{TRUE}, progress information is printed
#'   during fitting. Default \code{FALSE}.
#'
#' @section Algorithm:
#' \enumerate{
#'   \item \strong{Variable ordering.} Variables are ordered according to
#'     \code{xorder}. This may involve fitting a preliminary regression model.
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
#' In step 2, the following operations are applied in order:
#' \enumerate{
#'   \item The pre-cascade user values of \code{catzero} and \code{zero} are
#'     saved.
#'   \item The cascade is applied: \code{catzero[spike] <- TRUE} and
#'     \code{zero[catzero] <- TRUE}.
#'   \item \code{reset_spike()} is called for variables where spike eligibility
#'     criteria are not met (zero proportion outside \code{[min_prop, max_prop]}
#'     or binary variables). For reset variables, \code{spike} is set to
#'     \code{FALSE} and \code{catzero}/\code{zero} are restored to the
#'     user-specified pre-cascade values -- so only what the user explicitly
#'     requested is preserved.
#' }
#'
#' @return See \code{mfp2()} for details on the returned \code{mfp2} object.
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
                    min_prop,
                    max_prop,
                    force_max_fp,
                    has_offset,
                    verbose) {
  
  variables_x <- colnames(x)
  
  # Resolve GLM family objects once for all repeated internal model fits.
  # Public mfp2.default() already does this, but keeping it here makes direct
  # internal calls to fit_mfp() avoid repeated stats::gaussian()/binomial()/
  # poisson() construction as well. Cox remains the character string "cox".
  family_fit <- resolve_fit_model_family(family)
  
  if (verbose) {
    cat("\ni Initial degrees of freedom:\n")
    print(matrix(df, nrow = 1, dimnames = list("df", variables_x)),
          quote = FALSE)
  }
  
  # step 1: order variables ----------------------------------------------------
  variables_ordered <- variables_x
  
  if (length(variables_x) > 1) {
    variables_ordered <- order_variables(
      xorder        = xorder, x = x, y = y, family = family_fit,
      family_string = family_string, weights = weights, offset = offset,
      strata        = strata, method = method, control = control,
      nocenter      = nocenter
    )
  }
  
  if (verbose) {
    cat(sprintf("\ni Visiting order: %s\n",
                paste0(variables_ordered, collapse = ", ")))
  }
  
  # step 2: pre-process input --------------------------------------------------
  powers_current <- setNames(as.list(rep(1, ncol(x))), variables_ordered)
  
  # Re-order all per-variable vectors to match variables_ordered
  alpha        <- setNames(alpha,   variables_x)[variables_ordered]
  select       <- setNames(select,  variables_x)[variables_ordered]
  df           <- setNames(df,      variables_x)[variables_ordered]
  center       <- setNames(center,  variables_x)[variables_ordered]
  shift        <- setNames(shift,   variables_x)[variables_ordered]
  scale        <- setNames(scale,   variables_x)[variables_ordered]
  acdx         <- setNames(acdx,    variables_x)[variables_ordered]
  zero         <- setNames(zero,    variables_x)[variables_ordered]
  catzero      <- setNames(catzero, variables_x)[variables_ordered]
  spike        <- setNames(spike,   variables_x)[variables_ordered]
  force_max_fp <- setNames(force_max_fp, variables_x)[variables_ordered]
  powers       <- powers[variables_ordered]
  
  # Reorder x once to match the fixed variable order used by powers_current.
  # Hot-path transformation code can then access columns by name without
  # defensively copying/reordering the full matrix in every step.
  if (!identical(colnames(x), variables_ordered)) {
    x <- x[, variables_ordered, drop = FALSE]
  }
  if (!identical(colnames(x), variables_ordered)) {
    stop("Internal error: x column order does not match variables_ordered.",
         call. = FALSE)
  }
  
  # Assert repeated powers of 1 are not supported
  diff_one <- vapply(
    powers,
    function(v) {
      v <- v[!is.na(v)]
      length(v) > 1L && all(v == 1)
    },
    logical(1L)
  )
  
  if (any(diff_one)) {
    dfx <- df[diff_one]
    vars_invalid <- names(dfx)[dfx > 1]
    
    if (length(vars_invalid) > 0L) {
      stop(
        paste(
          "The powers of some variables are repeated and all equal to 1.",
          "Repeated powers equal to 1 are not supported.",
          sprintf(
            "i This applies to: %s.",
            paste0(vars_invalid, collapse = ", ")
          ),
          sep = "\n"
        ),
        call. = FALSE
      )
    }
  }
  
  # Force variables into the model by setting p-value threshold to 1
  if (!is.null(keep)) {
    select[which(names(select) %in% keep)] <- 1
  }
  
  # ACD transformation setup ---------------------------------------------------
  if (any(acdx)) {
    acdx           <- reset_acd(x, acdx)
    variables_acd  <- names(acdx)[acdx]
    powers_current <- modifyList(
      powers_current,
      sapply(variables_acd, function(v) c(1, NA), simplify = FALSE)
    )
    df[variables_ordered %in% variables_acd] <- 4
  }
  
  # Spike-at-zero handling -----------------------------------------------------
  # Order of operations is critical:
  #
  # 1. Save user-specified catzero and zero BEFORE the cascade. These are the
  #    values the user explicitly requested -- not what was implied by spike.
  #
  # 2. Apply the cascade: spike implies catzero implies zero. This is done here
  #    (as well as in mfp2()) in case fit_mfp() is called directly.
  #
  # 3. Call reset_spike() for variables that do not meet eligibility criteria.
  #    The new reset_spike() accepts the pre-cascade user values and restores
  #    catzero and zero to those values for reset variables -- so a user who
  #    only specified spike (not catzero or zero) gets a clean revert to
  #    standard continuous predictor treatment when spike is reset.
  
  user_catzero <- catzero   # pre-cascade user intent
  user_zero    <- zero      # pre-cascade user intent
  
  catzero[spike] <- TRUE    # spike implies catzero
  zero[catzero]  <- TRUE    # catzero implies zero
  
  # Temporary recoded data used only for spike eligibility checks.
  # Do not mutate the real x here.
  x_for_spike <- x
  
  zero_aligned_for_spike <- zero[colnames(x_for_spike)]
  cols_to_zero_for_spike <- which(zero_aligned_for_spike)
  
  if (length(cols_to_zero_for_spike) > 0L) {
    for (j in cols_to_zero_for_spike) {
      x_for_spike[x_for_spike[, j] <= 0, j] <- 0
    }
  }
  
  # Reset ineligible spike variables.
  # reset_spike() receives the temporarily recoded data, so its x == 0 check is
  # equivalent to checking original x <= 0 for effective zero variables.
  #
  # For reset variables, catzero and zero are restored to the original
  # user-specified values saved above. This is why user_catzero and user_zero
  # must be saved before the cascade.
  if (any(spike)) {
    result  <- reset_spike(
      x            = x_for_spike,
      spike        = spike,
      user_catzero = user_catzero,
      user_zero    = user_zero,
      min_prop     = min_prop,
      max_prop     = max_prop
    )
    
    spike   <- result$spike
    catzero <- result$catzero
    zero    <- result$zero
  }
  
  # Spike decision initialisation.
  # 2 means standard FP algorithm by default.
  spike_decision        <- rep(2L, length(variables_ordered))
  names(spike_decision) <- variables_ordered
  
  # Final zero recoding of the real x ------------------------------------------
  # Now that reset_spike() has produced the final zero/catzero/spike vectors, we
  # can safely mutate the actual x used by the MFP cycles.
  #
  # Only variables with final zero == TRUE are recoded. Therefore, variables whose
  # spike request was rejected and whose user-specified zero/catzero status was
  # FALSE remain untouched.
  
  zero_x       <- zero
  zero_aligned <- zero[colnames(x)]
  cols_to_zero <- which(zero_aligned)
  
  if (length(cols_to_zero) > 0L) {
    for (j in cols_to_zero) {
      x[x[, j] <= 0, j] <- 0
    }
    
    # Once x <= 0 has been physically recoded to 0, downstream transformation
    # functions no longer need to perform zero recoding again during the MFP
    # cycles. This avoids double handling of zero variables.
    zero_x <- setNames(rep(FALSE, length(variables_ordered)), variables_ordered)
    #zero_x[names(cols_to_zero)] <- FALSE
  }
  
  # Binary structural-zero indicators for catzero variables.
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
  
  # ACD parameters (cached for speed)
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
  
  # Number of events
  n_obs <- ifelse(family_string == "cox", sum(y[, 2]), nrow(x))
  
  # step 3: MFP backfitting cycles --------------------------------------------
  j         <- 1L
  converged <- FALSE
  
  prev_adj_params       <- vector("list", length = length(variables_ordered))
  names(prev_adj_params) <- variables_ordered
  
  while (j <= cycles) {
    if (verbose) {
      cat("\n---------------------")
      cat(sprintf("\ni Running MFP Cycle %d\n", j))
      cat("---------------------\n")
    }
    
    fit_best_cycle <- find_best_fp_cycle(
      x               = x,
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
      zero            = zero_x,        # all FALSE after non-positives set to 0
      catzero         = catzero_mat_list,  # named list of binary indicators
      spike           = spike,
      spike_decision  = spike_decision,
      acd_parameter   = acd_parameter,
      prev_adj_params = prev_adj_params,
      force_max_fp    = force_max_fp,
      has_offset      = has_offset,
      n_obs           = n_obs,
      verbose        = verbose
    )
    
    powers_updated        <- fit_best_cycle$powers_current
    spike_decision_updated <- fit_best_cycle$spike_decision
    prev_adj_params       <- fit_best_cycle$prev_adj_params
    
    
    #powers_same <- identical(powers_current, powers_updated)
    powers_same <- identical(
      normalize_powers_for_convergence(powers_current, spike_decision),
      normalize_powers_for_convergence(powers_updated, spike_decision_updated)
    )
    
    spike_same <- identical(spike_decision, spike_decision_updated)
    
    
    if (powers_same && spike_same) {
      converged <- TRUE
      if (verbose)
        cat(sprintf(
          "\ni Fractional polynomial fitting algorithm converged after %d cycle(s).\n", j
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
  
  # step 4: final transformation ----------------------------------------------
  # Backscale x before final FP transformation so that coefficients are on
  # the phi(x + shift) scale, matching what a user would expect from mfp2().
  # Scaling was applied upstream for numerical stability during power selection
  # and has no effect on power selection itself.
  if (any(scale != 1)) {
    x <- backscale_matrix(x, scale)
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
  
  data_transformed <- transform_matrix(
    x                  = x,
    power_list         = powers_current,
    center             = center,
    acdx               = acdx,
    acd_parameter_list = acd_parameter_final,
    zero               = zero,   # logical needed for centering even after zeroing
    catzero            = catzero,
    spike              = spike,
    reset_zero         = FALSE,
    spike_decision     = spike_decision
  )
  
  # Update catzero for final metadata.
  # Final-time `catzero` is logical metadata. It should reflect whether a
  # *_bin column is present in the final design matrix.
  catzero_effective <- catzero
  
  for (v in names(catzero_effective)) {
    if (isTRUE(spike[[v]])) {
      if (as.integer(spike_decision[[v]]) == 2L) {
        # Spike continuous-only or null: no spike *_bin column.
        catzero_effective[[v]] <- FALSE
      } else if (as.integer(spike_decision[[v]]) == 3L) {
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
  
  for (v in names(spike_decision)) {
    if (spike[v] && spike_decision[v] == 3L) {
      powers_current[[v]] <- NA
    }
  }
  
  # step 5: fit final model ---------------------------------------------------
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
  
  # Build and return mfp2 object ----------------------------------------------
  fit <- modifyList(
    modelfit$fit,
    list(
      centers         = data_transformed$centers,
      acd_parameter   = data_transformed$acd_parameter,
      convergence_mfp = converged,
      x_original = x[, names(powers_current[
        !sapply(powers_current, function(p) all(is.na(p))) |
          (spike & spike_decision == 3L)
      ]), drop = FALSE],
      y               = y,
      fp_terms        = create_fp_terms(powers_current, acdx, df, select, alpha,
                                        criterion, zero, catzero_effective, spike,
                                        spike_decision),
      transformations = data.frame(shift = shift, scale = scale, center = center),
      fp_powers       = powers_current,
      acd             = acdx,
      zero            = zero,
      catzero         = catzero_effective,
      catzero_list = catzero_mat_list,
      spike              = spike,
      spike_dec       = spike_decision
    )
  )
  
  class(fit) <- c("mfp2", class(fit))
  fit
}


#' Resolve a model family once for repeated internal fits
#'
#' @param family Character family name, family function, family object, or "cox".
#' @return A resolved GLM family object, or the character string "cox".
#' @keywords internal
#' @noRd
resolve_fit_model_family <- function(family) {
  if (is.character(family)) {
    if (identical(family, "cox")) {
      return(family)
    }
    return(switch(
      family,
      gaussian = stats::gaussian(),
      binomial = stats::binomial(),
      poisson  = stats::poisson()
    ))
  }
  
  if (is.function(family)) {
    return(family())
  }
  
  family
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
                               n_obs
                               ) {
  
  # order of names of powers does not change
  names_x <- names(powers_current)
  
  for (xi in names_x) {
    # iterate through all predictors xi and update xi's best FP power
    # in terms of loglikelihood
    # the result can be NA (variable not significant), linear, FP1, FP2, ...
    # note that the adjustment set and powers are given by powers_current
    # which is updated in each step
    fit_best_fp_step <- find_best_fp_step(
      x = x, # the order of columns does not matter since internal codes uses column names
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
    # Update parameters
    powers_current[[xi]] <- fit_best_fp_step$power_best
    spike_decision <- fit_best_fp_step$spike_decision
    # Store the returned adjustments keyed by xi
    prev_adj_params[[xi]] <- fit_best_fp_step$current_adj_params[[xi]]
  }
  
  list(powers_current = powers_current, spike_decision = spike_decision, 
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
  out <- Map(
    function(power, decision) {
      if (identical(as.integer(decision), 3L)) {
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
#'   * `1` – include both the continuous transformed term(s) and the binary
#'     spike-at-zero indicator.
#'   * `2` – include only the continuous transformed term(s).
#'   * `3` – include only the binary spike-at-zero indicator.
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
  
  if (length(spike_decision) != 1L ||
      is.na(spike_decision) ||
      !(as.integer(spike_decision) %in% c(1L, 2L, 3L))) {
    stop("`spike_decision` must be a single integer 1, 2, or 3.")
  }
  
  if (length(catzero) != 1L || is.na(catzero) || !is.logical(catzero)) {
    stop("`catzero` must be a single TRUE/FALSE value.")
  }
  
  spike_decision <- as.integer(spike_decision)
  
  # Binary-only spike: exactly one binary indicator column.
  if (spike_decision == 3L) {
    return(1L)
  }
  
  powers <- as.numeric(powers)
  
  # Unselected variable.
  if (all(is.na(powers))) {
    return(0L)
  }
  
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
#' To be used in \code{fit_mfp()}.
#' 
#' @param power_list list of powers created in `fit_mfp()`.
#' 
#' @return 
#' a matrix.
#' @keywords internal
#' @noRd
convert_powers_list_to_matrix <- function(power_list) {
  # Check the maximum number of powers i.e  FP2 has 2 while FP1 has 1
  psize <- sapply(power_list, length)
  maxp <- max(psize)
  
  # Create a new nested list of same length. This means that if FP1 was choosen
  # for x then the second power should be NA
  new_list_powers <- vector(mode = "list", length = length(power_list))
  for (i in 1:maxp) {
    new_list_powers[[i]] <- sapply(power_list, function(x) x[i])
  }
  # combine the powers and rename.
  matp <- do.call(cbind, new_list_powers)
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
#' @param fp_powers powers of the created FP terms.
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
                            spike_decision) {
                              
  vars <- names(fp_powers)
  
  acdx <- acdx[vars]
  df <- df[vars]
  select <- select[vars]
  alpha <- alpha[vars]
  zero <- zero[vars]
  catzero <- catzero[vars]
  spike <- spike[vars]
  spike_decision <- spike_decision[vars]
  
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
    
    #selected = sapply(fp_powers, function(p) ifelse(all(is.na(p)), FALSE, TRUE)),
    selected = mapply(function(p, sd) {
      if (sd == 3L) TRUE  # binary-only spike is still selected
      else !all(is.na(p))
    }, fp_powers, spike_decision),
    # final degrees of freedom
    df_final = mapply(
      calculate_df,
      fp_powers,
      spike_decision,
      catzero,
      SIMPLIFY = TRUE
    ), 
    convert_powers_list_to_matrix(fp_powers)
  )

    rownames(fp_terms) <- names(fp_powers)
  
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
  
  # Multiply each column by the corresponding scalar
  unscale_x <- sweep(x, 2, scalex[vnames], FUN = "*")
  return(unscale_x)
}

