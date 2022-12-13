
#' Fits glm and Cox proportional hazards regression models and returns the
#' p-values of variables ordered or not ordered as defined in option xorder.
#'
#' @param x A design matrix of dimension n * p+1 where n is the number of
#' observations and p the number of predictors including intercept for glm but
#' n*p for cox
#' @param y  A vector of response variable. For family = "cox", y is a
#' Surv object from the survival package, generated using Surv() function
#' @param family A character string naming a family function
#' @param xorder Determines the order of entry of the covariates into the
#' model-selection algorithm. The default is "ascending", which enters them in
#' decreasing order of significance in a multiple linear
#'  regression (most significant first). "descending" places them in reverse
#'  significance order, whereas "original" respects the original order in x
#' @param method A character string specifying the method for tie handling in cox
#' @param strata,weights,offset,control,rownames,nocenter parameters for glm or cox model. see glm() or coxph()

#' @return p-values for each predictor variable
#' @export
#'
order_variables <- function(x, y, weights, offset, family, method, strata, control,
                     rownames, nocenter, xorder) {
  # Convert factors to dummy variables if it exists...take it to the main function
  # x = model.matrix(as.formula(paste("~", paste(colnames(x), collapse="+"))),
  #                  data = as.data.frame(x))
  # save the family in factor form to be used later for gaussian
  fam <- family
  # number of rows of x or observations
  n <- dim(x)[1]
  # ============================Gaussian and Glm Models==========================
  if (family != "cox") {
    # glm.fit requires a function as a family not a character name
    if (is.character(family)) {
      family <- get(family, mode = "function", envir = parent.frame())
    }
    if (is.function(family)) family <- family()
    # family <- family()
    fit.full <- glm.fit(x = cbind(rep(1, n), x), y = y, weights = weights, offset = offset, family = family) # full model
    # Rank of the glm full model
    p1 <- fit.full$rank
    # There is an additional scale parameter to estimate in OLS regression. The
    # binomial and Poisson regression models have no scale parameter.
    if (fam == "gaussian") p1 <- p1 + 1
    # loglikelihood of the full model
    logl.full <- p1 - fit.full$aic / 2 # not that aic = -2logL + 2k
    # Deviance of the full model
    dev.full <- -2 * logl.full
    # we need to calculate p-values for each variable using likelihood ratio test
    varnames <- colnames(x)
    ns <- length(varnames)
    p.value <- loglikx <- dev <- df.reduced <- numeric(ns)
    names(p.value) <- names(dev) <- names(df.reduced) <- varnames
    for (i in 1:ns) {
      # remove one variable at a time and fit the reduced model
      fit.reduced <- glm.fit(
        x = cbind(rep(1, n), x[, -i, drop = FALSE]), y = y,
        weights = weights, offset = offset, family = family
      )
      # calculate the deviance of the reduced model
      p2 <- fit.reduced$rank
      if (fam == "gaussian") p2 <- p2 + 1
      # loglikelihood of the reduced model
      logl.reduced <- p2 - fit.reduced$aic / 2 # not that aic = -2logL + 2k
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
      x = x, y = y, strata = strata, weights = weights,
      offset = offset, control = control, method = method,
      rownames = rownames, nocenter = nocenter
    ) # full model
    # Deviance of the full cox model
    p1 <- fit.full$df
    # loglikelihood of the full model
    logl.full <- fit.full$logl
    # Deviance of the full model
    dev.full <- -2 * logl.full
    # we need to calculate p-values for each variable using likelihood ratio test
    varnames <- colnames(x)
    ns <- length(varnames)
    p.value <- loglikx <- dev <- df.reduced <- numeric(ns)
    names(p.value) <- names(dev) <- names(df.reduced) <- varnames
    for (i in 1:ns) {
      # remove one variable at a time and fit the reduced model
      fit.reduced <- fit_cox(
        x = x[, -i, drop = FALSE], y = y, strata = strata,
        weights = weights, offset = offset, control = control,
        method = method, rownames = rownames,
        nocenter = nocenter
      )
      # calculate the deviance of the reduced model
      p2 <- fit.reduced$df
      # loglikelihood of the reduced model
      logl.reduced <- fit.reduced$logl
      # Deviance of the reduced model
      dev[i] <- -2 * logl.reduced
      # degrees of freedom of the reduced model
      df.reduced[i] <- fit.reduced$df
      # loglik difference: -2(logL.reduced - logL.full)
      teststatic <- -2 * logl.reduced + 2 * logl.full
      # calculate the p.value
      p.value[i] <- pchisq(teststatic, df = p1 - p2, lower.tail = FALSE)
    }
  }
  # Order the p-values based on xorder
  pvalues <- switch(xorder,
    "ascending" = sort(p.value, decreasing = F),
    "descending" = sort(p.value, decreasing = F),
    "original" = p.value
  )
  return(list(pvalues = pvalues, dev = dev, df.reduced = df.reduced))
}
