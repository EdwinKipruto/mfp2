#' Function to find the best FP1 function based on deviance, aic, bic or sse
#' 
#' @param x The original untransformed matrix of predictors. Each continuous
#' predictor is assumed to have  shifted and scaled.
#' @param xi the name of the continuous predictor for which the FP function will 
#' be estimated. There are no binary or two-level variables allowed.
#' @param allpowers a named list of FP powers of all variables of interest, 
#' including xi. Not that these powers are updated during backfitting or MFP 
#' cycles.
#' @param powers a set of FP powers.
#' @param y the response variable.
#' @param method method for handling ties in cox model. see coxph() for 
#' explanation.
#' @param weights see glm or coxph for explanation.
#' @param offset see glm or coxph for explanation.
#' @param family A character name specifying the family name i,e "gaussian",
#' "binomial", "poisson" or "cox".
find_best_fp1_step <- function(y, 
                          x, 
                          xi, 
                          allpowers, 
                          powers,
                          family, 
                          method,
                          weights,
                          offset, 
                          strata, 
                          control,
                          rownames,
                          nocenter,
                          acdx) {
  # Generate FP1 data for x of interest (xi). A list with 8 new variables are
  # generated if default FP set is used
  df1 <- extract_adjustment_data(
    x = x, xi = xi, allpowers = allpowers, df = 2, acdx = acdx,
    powers = powers
  )
  # Matrix of adjustment variables
  adjdata <- df1$adjdata
  # List of FP1 data.
  fpdata <- df1$fpdata
  # N = number of observation and log(n) for bic calculation
  N <- nrow(x)
  logn <- log(N)
  # =============================================================================
  # Fit a null model-model without x of interest
  # =============================================================================
  fitnull <- fit_model(
    x = adjdata, y = y, family = family, method = method,
    weights = weights, offset = offset, strata = strata,
    control = control, rownames = rownames, nocenter = nocenter
  )
  # Total number of FP powers in the adjustment model
  # tFP <- calculate_number_fp_powers(df1$adjustpowers)
  # Total number of parameters including estimated FP powers in the null model
  # dfnull <- fitnull$df + tFP
  dfnull <- fitnull$df
  # Deviance, AIC and BIC of the null model
  devnull <- -2 * fitnull$logl
  # NULL model: aic = dev1 + 2(k1 + d) = dev1 + 2k1 + 2d
  # Lin. model: aic = dev2 + 2(k2 + d) = dev2 + 2k2 + 2d
  # FP1. model: aic = dev3 + 2(k3 + d) = dev3 + 2k3 + 2d
  # where k is the total number of parameters in the model (all betas including
  # the adjustment variables plus the powers of variable of interest and scale
  # parameter in gaussian). d is the number of powers in the adjustment model.
  # In each model, this results in a 2d constant quantity. As a result, we will
  # only use betas in calculating AIC and BIC of models instead of FP powers in
  # the adjustment model because they have no bearing on model selection.
  aic.null <- devnull + 2 * dfnull
  bic.null <- devnull + logn * dfnull
  # Sum of squares errors relevant for Gaussian model- for F test
  sse.null <- fitnull$SSE
  dev.roy.null <- deviance_gaussian(RSS = sse.null, weights = weights, n = N)
  df.sse.null <- N - (dfnull - 1) # subtract scale parameter
  # df.sse.null2 <- N-fitnull$df
  # =============================================================================
  # Fit 8 linear models for x of interest while adjusting for other variables.
  # =============================================================================
  nv <- length(fpdata)
  devs <- dev.roy <- sse <- aic <- bic <- dfx <- dfp1 <- numeric(nv)
  # Fit linear models for each FP1
  for (i in seq_len(nv)) {
    # combine each FP1 variable for x of interest with adjustment variables
    xout <- cbind(fpdata[[i]], adjdata)
    colnames(xout) <- c("newx", colnames(adjdata))
    fit1 <- fit_model(
      x = xout, y = y, family = family, method = method,
      weights = weights, offset = offset, strata = strata,
      control = control, rownames = rownames, nocenter = nocenter
    )
    # Deviance of the fitted model
    devs[i] <- -2 * fit1$logl
    # save the number of parameters of the fitted model. This does not take
    # into account the estimated FP powers.
    dfx[i] <- fit1$df
    # Calculate AIC and BIC for the fitted model. 1 is added because of FP1
    # power term for x of interest.
    # dfp1[i] = fit1$df + 1 + tFP
    # aic[i] = devs[i]  + 2*(dfp1[i])
    # bic[i] = devs[i]  + logn*(dfp1[i])
    dfp1[i] <- (fit1$df + 1) - 1 # we subtract scale parameter
    aic[i] <- devs[i] + 2 * (dfp1[i])
    bic[i] <- devs[i] + logn * (dfp1[i])
    # SSE relevant for a Gaussian model
    sse[i] <- fit1$SSE
    dev.roy[i] <- deviance_gaussian(RSS = sse[i], weights = weights, n = N)
  }
  # position of a linear function in set s: Fails when powers does not include 1
  lin.pos <- match(1, powers)
  # Deviance, AIC, BIC and SSE of a linear function
  dev.linear <- devs[lin.pos]
  dev.roy.linear <- dev.roy[lin.pos]
  # aic.linear = dev.linear + 2*(dfx[lin.pos]+tFP)
  # bic.linear = dev.linear + logn*(dfx[lin.pos]+tFP)
  aic.linear <- dev.linear + 2 * (dfx[lin.pos])
  bic.linear <- dev.linear + logn * (dfx[lin.pos])
  sse.linear <- sse[lin.pos]
  # Degrees of freedom of sse assuming linearity
  # df.sse.linear <- N-(dfx[lin.pos]+ tFP)
  df.sse.linear <- N - (dfx[lin.pos] - 1)
  # df.sse.linear2 <- N-(dfx[lin.pos]-1) # TO REMOVE

  # Degrees of freedom of sse assuming non-linearity, subtract 1 because dfp1 includes scale parameter
  df.sse.bestfp1 <- ifelse(which.min(sse) == lin.pos, df.sse.linear, N - (dfp1[which.min(sse)]))
  # df.sse.bestfp12 <- ifelse(which.min(sse)==lin.pos, df.sse.linear2, N-dfp1[which.min(sse)]) # TO REMOVE
  # combine deviance/AIC/BIC of null, linear and best FP1
  dev.all <- setNames(c(devnull, dev.linear, devs[which.min(devs)]), c("Null", "Linear", "FP1"))
  # deviance of gaussian calculated using royston formula instead of loglikelihood for gaussian
  dev.roy.all <- setNames(c(dev.roy.null, dev.roy.linear, dev.roy[which.min(dev.roy)]), c("Null", "Linear", "FP1"))
  aic.all <- setNames(c(aic.null, aic.linear, aic[which.min(aic)]), c("Null", "Linear", "FP1"))
  bic.all <- setNames(c(bic.null, bic.linear, bic[which.min(bic)]), c("Null", "Linear", "FP1"))
  sse.all <- setNames(c(sse.null, sse.linear, sse[which.min(devs)]), c("Null", "Linear", "FP1"))
  # # deviance of gaussian calculated using royston formula instead of loglikelihood for gaussian
  # dev.gaus.royston <- unlist(lapply(sse.all, function(x) deviance_gaussian(RSS = x, weights = weights, n = N)))
  # names(dev.gaus.royston) <- c("Null","Linear","FP1")
  # combine degrees of freedom for null, linear and best fp1
  df.all <- setNames(c(df.sse.null, df.sse.linear, df.sse.bestfp1), c("Null", "Linear", "FP1"))
  # The best FP1 power based on deviance, aic, bic, and sse
  s1 <- powers
  aic[lin.pos] <- aic.linear # correct aic and bic for the linear position
  bic[lin.pos] <- bic.linear
  fn.bestfp1 <- setNames(
    c(
      s1[which.min(devs)],
      s1[which.min(aic)],
      s1[which.min(bic)],
      s1[which.min(sse)],
      s1[which.min(dev.roy)]
    ),
    c("dev", "aic", "bic", "sse", "dev.roy")
  )
  # return
  outx <- list(
    dev.all = dev.all, aic.all = aic.all, bic.all = bic.all,
    sse.all = sse.all, df.all = df.all, fn.bestfp1 = fn.bestfp1,
    dev.roy.all = dev.roy.all
  )
  return(outx)
}
