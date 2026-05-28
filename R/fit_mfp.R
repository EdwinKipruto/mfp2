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
#'   treat non-positive values as zero AND have a binary indicator
#'   \eqn{I(x = 0)} automatically created and included in the model. Must
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
                    verbose) {
  
  variables_x <- colnames(x)
  
  if (verbose) {
    cat("\ni Initial degrees of freedom:\n")
    print(matrix(df, nrow = 1, dimnames = list("df", variables_x)),
          quote = FALSE)
  }
  
  # step 1: order variables ----------------------------------------------------
  variables_ordered <- variables_x
  
  if (length(variables_x) > 1) {
    variables_ordered <- order_variables(
      xorder        = xorder, x = x, y = y, family = family,
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
  force_max_fp <- force_max_fp[variables_ordered]
  powers       <- powers[variables_ordered]
  
  # Assert repeated powers of 1 is not supported
  diff_one <- unlist(lapply(powers, function(v) all(v %in% 1) && length(v) != 1))
  if (any(diff_one)) {
    dfx <- df[diff_one]
    if (any(dfx > 1))
      stop(
        "The powers of some variables are repeated and all equal to 1. ",
        "Set df for those variables to 1 or change the powers.\n",
        sprintf("i This applies to: %s.", paste0(names(dfx > 1), collapse = ", ")),
        call. = FALSE
      )
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
  
  if (any(spike)) {
    result  <- reset_spike(x, spike, user_catzero, user_zero, min_prop, max_prop)
    spike   <- result$spike
    catzero <- result$catzero
    zero    <- result$zero
  }
  
  # Spike decision initialisation (2 = standard FP algorithm, default)
  spike_decision       <- rep(2, length(variables_ordered))
  names(spike_decision) <- variables_ordered
  
  # Set non-positive values to 0 for zero/catzero variables ------------------
  zero_x        <- zero
  zero_aligned  <- zero[colnames(x)]
  cols_to_zero  <- which(zero_aligned)
  
  if (length(cols_to_zero) > 0L) {
    for (j in cols_to_zero) {
      x[x[, j] <= 0, j] <- 0
    }
    zero_x <- setNames(rep(FALSE, ncol(x)), variables_ordered)
  }
  
  # Binary indicators Z for catzero variables (Z=1 if x=0, Z=0 if x>0)
  catzero_list       <- lapply(names(catzero), function(v) {
    if (catzero[[v]]) as.integer(x[, v] == 0) else NULL
  })
  names(catzero_list) <- names(catzero)
  
  # ACD parameters (cached for speed)
  acd_parameter       <- lapply(seq_len(ncol(x)), function(i) {
    if (acdx[i]) fit_acd(x[, i], powers = powers[[i]], zero = zero[i]) else NULL
  })
  names(acd_parameter) <- colnames(x)
  
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
      x              = x,
      y              = y,
      powers_current = powers_current,
      df             = df,
      weights        = weights,
      offset         = offset,
      family         = family,
      family_string  = family_string,
      criterion      = criterion,
      select         = select,
      alpha          = alpha,
      keep           = keep,
      powers         = powers,
      ftest          = ftest,
      control        = control,
      rownames       = rownames(x),
      strata         = strata,
      nocenter       = nocenter,
      method         = method,
      acdx           = acdx,
      zero           = zero_x,        # all FALSE after non-positives set to 0
      catzero        = catzero_list,  # named list of binary indicators
      spike          = spike,
      spike_decision = spike_decision,
      acd_parameter  = acd_parameter,
      prev_adj_params = prev_adj_params,
      force_max_fp   = force_max_fp,
      verbose        = verbose
    )
    
    powers_updated        <- fit_best_cycle$powers_current
    spike_decision_updated <- fit_best_cycle$spike_decision
    prev_adj_params       <- fit_best_cycle$prev_adj_params
    
    if (identical(powers_current, powers_updated)) {
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
  
  data_transformed <- transform_matrix(
    x              = x,
    power_list     = powers_current,
    center         = center,
    acdx           = acdx,
    zero           = zero,   # logical needed for centering even after zeroing
    catzero        = catzero,
    spike_decision = spike_decision
  )
  
  # step 5: fit final model ---------------------------------------------------
  modelfit <- fit_model(
    x        = data_transformed$x_transformed,
    y        = y,
    family   = family,
    weights  = weights,
    offset   = offset,
    method   = method,
    strata   = strata,
    control  = control,
    rownames = rownames(data_transformed$x_transformed),
    nocenter = nocenter,
    fast     = FALSE
  )
  
  # Build and return mfp2 object ----------------------------------------------
  fit <- modifyList(
    modelfit$fit,
    list(
      centers         = data_transformed$centers,
      acd_parameter   = data_transformed$acd_parameter,
      convergence_mfp = converged,
      x_original      = x[, names(powers_current[
        !sapply(powers_current, function(p) all(is.na(p)))
      ]), drop = FALSE],
      y               = y,
      fp_terms        = create_fp_terms(powers_current, acdx, df, select, alpha,
                                        criterion, zero, catzero, spike,
                                        spike_decision),
      transformations = data.frame(shift = shift, scale = scale, center = center),
      fp_powers       = powers_current,
      acd             = acdx,
      zero            = zero,
      catzero         = catzero,
      catzero_list    = catzero_list,
      spike_dec       = spike_decision
    )
  )
  
  class(fit) <- c("mfp2", class(fit))
  fit
}

#' Helper to order variables for mfp2 algorithm
#' 
#' To be used in \code{fit_mfp()}.
#' 
#' @param xorder a string determining the order of entry of the covariates
#' into the model-selection algorithm. The default is `ascending`, which enters
#' them by ascending p-values, or decreasing order of significance in a
#' multiple regression (i.e. most significant first).
#' `descending` places them in reverse significance order, whereas 
#' `original` respects the original order in `x`.
#' @param x a design matrix of dimension n * p where n is the number of
#' observations and p the number of predictors including intercept for glms,
#' or excluding intercept for Cox models. 
#' @param y  a vector of responses for glms, or a `Surv` object generated using
#' the [survival::Surv()] function for Cox models. 
#' @param family a character string naming a family function supported by
#' `glm()` or "cox" for Cox models.
#' @param family_string A character string representing the selected family, 
#'   e.g., "gaussian".
#' @param weights,offset parameters for both glm and Cox models, see either
#' [stats::glm()] or [survival::coxph()] depending on family. 
#' @param strata,method,control,nocenter Cox model specific parameters, see
#' [survival::coxph()].
#' @param ... passed to `order_variables_by_significance`.
#' 
#' @return 
#' A vector of the variable names in `x`, ordered according to `xorder`.
#' 
#' @import utils
order_variables <- function(xorder = "ascending",
                            x = NULL, 
                            ...) {
  names_ordered <- colnames(x)
  
  if (xorder != "original") {
    names_ordered <- order_variables_by_significance(xorder = xorder, x = x, ...)
  }
  
  names_ordered
}

#' @describeIn order_variables Order by significance in regression model. The 
#' number of columns of `x` should be greater than 1 for Cox models.
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
                                            nocenter) {
  
  # Convert factors to dummy variables if it exists...take it to the main function
  # x = model.matrix(as.formula(paste("~", paste(colnames(x), collapse="+"))),
  #                  data = as.data.frame(x))
  # save the family in factor form to be used later for gaussian
  fam <- family
  # number of rows of x or observations
  n <- dim(x)[1]
  
  if (family_string != "cox") {
    # glm.fit requires a function as a family not a character name
    if (is.character(family)) {
      family <- get(family, mode = "function", envir = parent.frame())
    }
    if (is.function(family)) family <- family()
    # family <- family()
    fit.full <- glm.fit(
      x = cbind(rep(1, n), x), y = y, weights = weights, offset = offset,
      family = family
    )
    
    # Rank of the glm full model
    p1 <- fit.full$rank
    
    # There is an additional scale parameter to estimate in OLS regression. The
    # binomial and Poisson regression models have no scale parameter.
    if (family_string == "gaussian") p1 <- p1 + 1
    
    # loglikelihood of the full model calculated based on  aic = -2logL + 2k
    logl.full <- p1 - fit.full$aic / 2 
    
    # Deviance of the full model
    dev.full <- -2 * logl.full
    
    # we need to calculate p-values for each variable using likelihood ratio test
    varnames <- colnames(x)
    ns <- length(varnames)
    p.value <- loglikx <- dev <- df.reduced <- numeric(ns)
    names(p.value) <- names(dev) <- names(df.reduced) <- varnames
    for (i in 1:ns) {
      # remove one variable at a time and fit the reduced model. 
      # only works if you have more than one variable due to (-i)
      fit.reduced <- glm.fit(
        x = cbind(rep(1, n), x[, -i, drop = FALSE]),
        y = y,
        weights = weights, 
        offset = offset, 
        family = family
      )
      # calculate the deviance of the reduced model-model without the ith variable
      p2 <- fit.reduced$rank
      
      if (family_string == "gaussian") p2 <- p2 + 1
      # loglikelihood of the reduced model
      logl.reduced <- p2 - fit.reduced$aic / 2 
      
      # Deviance of the reduced model
      dev[i] <- -2 * logl.reduced
      
      # degrees of freedom of the reduced model
      df.reduced[i] <- p2
      
      # loglik difference: -2(logL.reduced - logL.full)
      teststatic <- -2 * logl.reduced + 2 * logl.full
      
      # calculate the Chi-square p.value
      p.value[i] <- pchisq(teststatic, df = p1 - p2, lower.tail = FALSE)
    }
  } else { # cox model
    fit.full <- fit_cox(
      x = x, 
      y = y, 
      strata = strata,
      weights = weights,
      offset = offset,
      control = control, 
      method = method,
      rownames = rownames(x),
      nocenter = nocenter
    ) 
    
    # Degrees of freedom for the full cox model
    p1 <- fit.full$df
    
    # loglikelihood of the full cox model
    logl.full <- fit.full$logl
    
    # Deviance of the full cox model
    dev.full <- -2 * logl.full
    
    # we need to calculate p-values for each variable using likelihood ratio test
    varnames <- colnames(x)
    ns <- length(varnames)
    p.value <- loglikx <- dev <- df.reduced <- numeric(ns)
    names(p.value) <- names(dev) <- names(df.reduced) <- varnames
    for (i in 1:ns) {
      # remove one variable at a time and fit the reduced model
      fit.reduced <- fit_cox(
        x = x[, -i, drop = FALSE],
        y = y, strata = strata,
        weights = weights, 
        offset = offset, 
        control = control,
        method = method,
        rownames = rownames(x),
        nocenter = nocenter
      )
      # Degrees of freedom for the reduced model
      p2 <- fit.reduced$df
      
      # loglikelihood of the reduced model
      logl.reduced <- fit.reduced$logl
      
      # Deviance of the reduced model
      dev[i] <- -2 * logl.reduced
      
      # degrees of freedom of the reduced model
      df.reduced[i] <- p2
      
      # loglik difference: -2(logL.reduced - logL.full)
      teststatic <- -2 * logl.reduced + 2 * logl.full
      
      # calculate the p.value
      p.value[i] <- pchisq(teststatic, df = p1 - p2, lower.tail = FALSE)
    }
  }
  
  # Order the p-values based on xorder
  pvalues <- switch(xorder,
                    "descending" = sort(p.value, decreasing = TRUE),
                    "ascending"  = sort(p.value, decreasing = FALSE),
                     "original"  = p.value 
  )
  
  names(pvalues)
}

#' Helper to reset acd transformation for variables with few values
#' 
#' To be used in \code{fit_mfp()}.
#' This function resets the `acdx` parameter (logical vector) of variables with
#' less than 5 distinct values to `FALSE`.
#' 
#' @param x a design matrix of dimension nobs x nvars where nvars is the number 
#' of predictors excluding an intercept.  
#' @param acdx a named logical vector of length nvars indicating which continuous
#' variables should undergo the approximate cumulative distribution (ACD) 
#' transformation. May be ordered differently than the columns of `x`.
#' 
#' @return 
#' Logical vector of same length as `acdx`.
reset_acd <- function(x, acdx) {
  # exit early if all acdx values are FALSE
  if (!any(acdx)) return(acdx)
  
  names_acd <- names(acdx)[which(acdx == TRUE)]
  
  # number of unique values of each column in acdx
  n_unique <- apply(x[, names_acd, drop = FALSE], 2, 
                    function(col) length(unique(col)))
  
  ind_reset <- which(n_unique < 5)
  
  if (length(ind_reset) > 0) {
    acdx[names_acd][ind_reset] <- FALSE
    warning("i For any variable with fewer than 5 unique values no acd transformation can be performed.\n", 
            sprintf("i The requested acd transform has been reset to FALSE for the following variables: %s.", 
                    paste0(names_acd[ind_reset], collapse = ", ")))
  }
  
  acdx
}
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
#' If \code{reset_spike()} only resets \code{spike} (as in the original
#' implementation), \code{catzero} and \code{zero} remain \code{TRUE} for
#' ineligible variables even though the user never asked for them. This
#' function corrects that by accepting the pre-cascade user values and
#' restoring them for any variable where \code{spike} is reset.
#'
#' @param x A numeric matrix or data frame with column names. Only columns
#'   named in \code{spike} are examined.
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
#'     \code{zero = TRUE}. The binary indicator \eqn{I(x = 0)} is still
#'     added (as explicitly requested via \code{catzero}), but the SAZ
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
reset_spike <- function(x, spike, user_catzero, user_zero,
                        min_prop = 0.05, max_prop = 0.95) {
  
  # Early exit: no spike variables to evaluate
  if (!any(spike)) {
    return(list(spike = spike, catzero = user_catzero, zero = user_zero))
  }
  
  names_spike <- names(spike)[spike]
  
  # Proportion of zeros for each spike variable
  prop_zero <- colMeans(x[, names_spike, drop = FALSE] == 0, na.rm = TRUE)
  
  # Binary: exactly two unique finite values in the column
  is_binary <- vapply(names_spike, function(v) {
    length(unique(x[is.finite(x[, v]), v])) == 2L
  }, logical(1L))
  
  # Identify variables to reset by each reason
  to_reset_prop   <- names_spike[prop_zero < min_prop | prop_zero > max_prop]
  to_reset_binary <- names_spike[is_binary]
  to_reset        <- union(to_reset_prop, to_reset_binary)
  
  # Warn about proportion resets
  if (length(to_reset_prop) > 0L) {
    warning(
      "The spike-at-zero option has been reset for the following variable(s) ",
      "because the zero proportion is outside [", min_prop, ", ", max_prop, "]: ",
      paste(to_reset_prop, collapse = ", "), ".",
      call. = FALSE
    )
  }
  
  # Warn about binary resets
  if (length(to_reset_binary) > 0L) {
    warning(
      "The spike-at-zero option has been reset for the following variable(s) ",
      "because they are binary (exactly two unique values): ",
      paste(to_reset_binary, collapse = ", "), ".",
      call. = FALSE
    )
  }
  
  # Initialise output vectors from user-supplied pre-cascade originals
  catzero <- user_catzero
  zero    <- user_zero
  
  if (length(to_reset) > 0L) {
    # Reset spike
    spike[to_reset] <- FALSE
    
    # Restore catzero and zero to the user's original pre-cascade values for
    # the reset variables. This undoes the implicit cascade that was applied
    # in fit_mfp() before this function was called:
    #   catzero[spike] <- TRUE   (spike implies catzero)
    #   zero[catzero]  <- TRUE   (catzero implies zero)
    # After restoration, only what the user explicitly requested remains.
    # See the @section Cascade restoration examples for the three scenarios.
    catzero[to_reset] <- user_catzero[to_reset]
    zero[to_reset]    <- user_zero[to_reset]
  }
  
  list(spike = spike, catzero = catzero, zero = zero)
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
#' 
#' @return 
#' A list with updated components `powers_current` (current FP powers for all 
#' variables), `spike_decision` (updated spike-at-zero decisions), and
#' `prev_adj_params` (adjustment variable transformations to be used in the next
#' cycle).
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
                               force_max_fp
                               ) {
  
  # order of names of powers does not change
  names_x <- names(powers_current)
  
  for (i in 1:ncol(x)) {
    # iterate through all predictors xi and update xi's best FP power
    # in terms of loglikelihood
    # the result can be NA (variable not significant), linear, FP1, FP2, ...
    # note that the adjustment set and powers are given by powers_current
    # which is updated in each step
    xi <- names_x[i]
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
      verbose = verbose
    )
    # Update parameters
    powers_current[[i]] <- fit_best_fp_step$power_best
    spike_decision <- fit_best_fp_step$spike_decision
    # Store the returned adjustments keyed by xi
    prev_adj_params[[names_x[i]]] <- fit_best_fp_step$current_adj_params[[names_x[i]]]
  }
  
  list(powers_current = powers_current, spike_decision = spike_decision, 
       prev_adj_params = prev_adj_params)
}

#' Calculate degrees of freedom for a transformed variable
#' 
#' Helper function used in `fit_mfp()` to determine the final number of 
#' degrees of freedom (df) contributed by a variable, depending on its powers 
#' and the spike-at-zero decision.
#' 
#' @param powers Numeric vector of selected powers for the variable. Can 
#' contain `NA` values (e.g., for ACD terms). If all values are `NA`, 
#' the variable is considered unselected.
#' @param spike_decision Integer scalar (1, 2, or 3) specifying spike-at-zero 
#'   handling:
#'   * `1` – include both the continuous transformed term(s) and the binary 
#'     spike-at-zero indicator.
#'   * `2` – include only the continuous transformed term(s).
#'   * `3` – include only the binary spike-at-zero indicator.
#
#' @details 
#' * If all entries in `powers` are `NA`, the df is `0` regardless of 
#'   `spike_decision`.
#' * If `spike_decision = 3`, the df is `1` (only the binary indicator).
#' * If the variable is modeled linearly (exactly `powers = 1`), then df = 1.
#' * Otherwise, for fractional polynomials of degree *m* (number of 
#'   non-`NA` powers), df = 2 * m.
#' * If `spike_decision = 1`, one additional df is added for the binary 
#'   indicator.
#' An example calculation: if p is the power(s) and p = c(1,2), then df = 4 
#' but if p = NA then df = 0.
#' 
#' @return 
#' Integer scalar giving the degrees of freedom for the variable.
#' @examples
#'\dontrun{
#' calculate_df(c(1, 2), 2)   # df = 4 (two powers, no spike)
#' calculate_df(1, 1)         # df = 2 (linear + binary spike)
#' calculate_df(c(NA, NA), 1) # df = 0 (unselected variable)
#' calculate_df(2, 3)        # df = 1 (binary spike only)
#' }
#' @keywords internal
calculate_df <- function(powers, spike_decision) {
  
  #------------------
  # Input checks
  #------------------
  if (length(spike_decision) != 1L || !(spike_decision %in% c(1L, 2L, 3L))) {
    stop("`spike_decision` must be a single integer 1, 2, or 3.")
  }
  
  # Convert to numeric in case powers is logical NA
  powers <- as.numeric(powers)
  
  # unselected variable: no df regardless of spike_decision
  if (all(is.na(powers))) {
    return(0L)
  }
  
  # case: binary only
  if (spike_decision == 3) {
    return(1L)
  }
  # base df from powers
  # Remove NAs in powers (2,NA) or (NA,1) etc. Happens because of acd
  # make sure to drop names by using as.numeric
  p <- as.numeric(powers[!is.na(powers)])
    # df of linear function
  if (identical(p, 1)) {
    df <- 1L
    # df of fpm. Note that df = 2m where m is the degree. if length = 1 then
    # degree = 1 and df = 2 etc
  } else {
    df <- 2L * length(p)
  }
  
  # add spike binary if spike_decision = 1
  if (spike_decision == 1) {
    df <- df + 1L
  }
  
  return(df)
}

#' Helper to convert a nested list with same or different length into a matrix
#' 
#' To be used in \code{fit_mfp()}.
#' 
#' @param power_list list of powers created in `fit_mfp()`.
#' 
#' @return 
#' a matrix.
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
    # presence / absence in final model encoded by NAs in fp_powers
    selected = sapply(fp_powers, function(p) ifelse(all(is.na(p)), FALSE, TRUE)),
    # final degrees of freedom
    df_final = mapply(calculate_df, fp_powers, spike_decision), 
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
#'
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

