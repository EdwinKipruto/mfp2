
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
#' @keywords internal
#' @noRd
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
#' @keywords internal
#' @noRd
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
  # If there is only one predictor, there is no ordering problem to solve.
  # Return the current column name unchanged.
  if (ncol(x) <= 1L) {
    return(colnames(x))
  }
  
  # Store predictor names once. These names are used both to label p-values
  # and to return the final visiting order.
  predictor_names <- colnames(x)
  
  # Number of candidate predictors to rank.
  n_predictors <- length(predictor_names)
  
  # Store one likelihood-ratio-test p-value per predictor.
  # A smaller p-value means the model fit worsens more when that predictor is
  # removed, so the predictor is treated as more important.
  p_values <- numeric(n_predictors)
  names(p_values) <- predictor_names
  
  if (family_string != "cox") {
    # GLM ordering ------------------------------------------------------------
    
    # Number of observations. Used to create the intercept column once.
    n_obs <- nrow(x)
    
    # glm.fit() requires a family object, not a character string.
    # For example, "gaussian" must become gaussian().
    if (is.character(family)) {
      family <- get(
        family,
        mode = "function",
        envir = parent.frame()
      )
    }
    
    # If the user supplied a family function, evaluate it to get the family
    # object expected by glm.fit().
    if (is.function(family)) {
      family <- family()
    }
    
    # Build the full GLM design matrix once.
    # x does not contain an intercept at this point, so column 1 is added here.
    x_full <- cbind(
      "(Intercept)" = rep.int(1, n_obs),
      x
    )
    
    # Fit the full model containing the intercept and all predictors.
    full_fit <- glm.fit(
      x = x_full,
      y = y,
      weights = weights,
      offset = offset,
      family = family
    )
    
    # Effective parameter count for the full model.
    # For non-Gaussian GLMs, the rank is the number of estimated regression
    # coefficients. For Gaussian GLMs, the residual scale/dispersion is also
    # estimated, so one extra parameter is counted for the likelihood/AIC
    # relationship used below.
    full_df <- full_fit$rank
    
    if (family_string == "gaussian") {
      full_df <- full_df + 1L
    }
    
    # glm.fit() stores AIC = -2 * logLik + 2 * k.
    # Rearranging gives logLik = k - AIC / 2, where k is full_df.
    full_loglik <- full_df - full_fit$aic / 2
    
    for (predictor_index in seq_len(n_predictors)) {
      # Remove one predictor at a time from the full design matrix.
      #
      # Important indexing detail:
      #   x_full column 1 is the intercept.
      #   x column 1 is predictor 1.
      #   Therefore predictor_index in x corresponds to column
      #   predictor_index + 1 in x_full.
      #
      # The +1 is only needed in the GLM branch because we manually added an
      # intercept column to x_full.
      reduced_x <- x_full[, -(predictor_index + 1L), drop = FALSE]
      
      # Fit the reduced model without the current predictor.
      reduced_fit <- glm.fit(
        x = reduced_x,
        y = y,
        weights = weights, 
        offset = offset, 
        family = family
      )
      
      # Effective parameter count for the reduced model.
      reduced_df <- reduced_fit$rank
      
      # Same Gaussian adjustment as for the full model: count the estimated
      # residual scale/dispersion parameter.
      if (family_string == "gaussian") {
        reduced_df <- reduced_df + 1L
      }
      
      # Recover reduced-model log-likelihood from AIC.
      reduced_loglik <- reduced_df - reduced_fit$aic / 2
      
      # Likelihood-ratio statistic:
      #   -2 * (logLik_reduced - logLik_full)
      #
      # This is equivalent to:
      #   2 * (logLik_full - logLik_reduced)
      #
      # Larger values indicate that removing the predictor worsens the model.
      lrt_statistic <- -2 * reduced_loglik + 2 * full_loglik
      
      # Difference in effective degrees of freedom between full and reduced
      # models. This is usually 1 for a single numeric predictor, but can be
      # different if the model matrix is rank-deficient.
      lrt_df <- full_df - reduced_df
      
      # Convert the likelihood-ratio statistic to a p-value.
      p_values[predictor_index] <- pchisq(
        lrt_statistic,
        df = lrt_df,
        lower.tail = FALSE
      )
    }
    
  } else {
    # Cox ordering ------------------------------------------------------------
    
    # Preserve row names once and reuse them in all Cox fits.
    row_names <- rownames(x)
    
    # Fit the full Cox model containing all predictors.
    full_fit <- fit_cox(
      x = x, 
      y = y, 
      strata = strata,
      weights = weights,
      offset = offset,
      control = control, 
      method = method,
      rownames = row_names,
      nocenter = nocenter
    ) 
    
    # Effective degrees of freedom and log-likelihood for the full Cox model.
    full_df <- full_fit$df
    full_loglik <- full_fit$logl
    
    for (predictor_index in seq_len(n_predictors)) {
      # Remove one predictor at a time.
      #
      # No +1 is needed here because Cox design matrix x does not have an
      # added intercept column. Predictor predictor_index is column
      # predictor_index in x.
      reduced_x <- x[, -predictor_index, drop = FALSE]
      
      # Fit the reduced Cox model without the current predictor.
      reduced_fit <- fit_cox(
        x = reduced_x,
        y = y,
        strata = strata,
        weights = weights, 
        offset = offset, 
        control = control,
        method = method,
        rownames = row_names,
        nocenter = nocenter
      )
      
      # Effective degrees of freedom and log-likelihood for the reduced model.
      reduced_df <- reduced_fit$df
      reduced_loglik <- reduced_fit$logl
      
      # Likelihood-ratio statistic:
      #   -2 * (logLik_reduced - logLik_full)
      #
      # Larger values indicate that removing the predictor worsens the model.
      lrt_statistic <- -2 * reduced_loglik + 2 * full_loglik
      
      # Difference in degrees of freedom between the full and reduced Cox
      # models.
      lrt_df <- full_df - reduced_df
      
      # Convert the likelihood-ratio statistic to a p-value.
      p_values[predictor_index] <- pchisq(
        lrt_statistic,
        df = lrt_df,
        lower.tail = FALSE
      )
    }
  }
  
  # Return variable names ordered by likelihood-ratio-test p-value.
  #
  # ascending:
  #   Most significant predictors first. This is the default MFP visiting order.
  #
  # descending:
  #   Least significant predictors first.
  #
  # original:
  #   Preserve the input column order.
  if (xorder == "ascending") {
    return(names(sort(p_values, decreasing = FALSE)))
  }
  
  if (xorder == "descending") {
    return(names(sort(p_values, decreasing = TRUE)))
  }
  
  predictor_names
}