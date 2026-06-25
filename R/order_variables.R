
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
  if (ncol(x) <= 1L) {
    return(colnames(x))
  }
  
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