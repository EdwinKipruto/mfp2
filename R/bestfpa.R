# A function that fits 64 FP1 models for x of interest and acd(x) and returns best
# FP1(P1,P2) function and corresponding deviance. Best FP function is the one
# that has the smallest deviance or aic or bic or RSS.Linear model is fitted for
# each pair of powers as we adjust for the other variables in the model, the
# adjustment variables can also include acd variables.

# @param x The original untransformed matrix of predictors. Each continuous
# predictor is assumed to have  shifted and scaled.
# xi = the name of the continuous predictor for which the FP function will be
# estimated. There are no binary or two-level variables allowed.
# allpowers = a named list of FP powers of all variables of interest, including
# xi. Not that these powers are updated during backfitting or MFP cycles.
# @param powers The set of FP powers.
# @param y The response variable

# @param method method for handling ties in cox model. see coxph() for explanation
# @param weights see glm or coxph for explanation
# @param offset see glm or coxph for explanation
# @param family A character name specifying the family name i,e "gaussian",
# "binomial", "poisson" or "cox"
# acdx = an indicator of acd variables. Either TRUE or FALSE
bestfpa <- function(y, x, xi, allpowers, powers, family, method, weights,
                    offset, strata, control, rownames, nocenter, acdx) {
  # Generate FPa data for x of interest (xi). If the default FP power set is
  # used, 64 pairs of new variables are created.
  df1 <- adjustment.data(
    x = x, xi = xi, allpowers = allpowers, df = 4,
    powers = powers, acdx = acdx
  )
  # Matrix of adjustment variables
  adjdata <- df1$adjdata
  # print(head(adjdata))
  # FPa data for x of interest
  fpdata <- df1$fpdata
  nv <- length(fpdata)
  # log(n) for bic calculation
  N <- nrow(x)
  logn <- log(N)
  # =============================================================================
  # Fit a model without xi and axi--the null model. Model M6 in R&S 2016
  # =============================================================================
  fitnull <- model.fit(
    x = adjdata, y = y, family = family, method = method,
    weights = weights, offset = offset, strata = strata,
    control = control, rownames = rownames, nocenter = nocenter
  )
  # Total number of parameters including estimated FP powers in the null model
  dfnull <- fitnull$df
  # Deviance, AIC and BIC of the null model
  devnull <- -2 * fitnull$logl
  aic.null <- devnull + 2 * dfnull
  bic.null <- devnull + logn * dfnull
  sse.null <- fitnull$SSE
  dev.roy.null <- dev.gaussian(RSS = sse.null, weights = weights, n = N)
  df.null <- N - (dfnull - 1)

  # =============================================================================
  # Fit a model with linear in xi--Model M4 in R&S 2016
  # =============================================================================
  xk <- cbind(x[, xi], adjdata)
  colnames(xk) <- c(xi, colnames(adjdata))
  fit.lin.xi <- model.fit(
    x = xk, y = y, family = family, method = method,
    weights = weights, offset = offset, strata = strata,
    control = control, rownames = rownames,
    nocenter = nocenter
  )
  # Total number of parameters in the fitted model.
  dflinxi <- fit.lin.xi$df
  # Deviance, AIC, BIC and SSE of the null model
  devlinxi <- -2 * fit.lin.xi$logl
  aic.linxi <- devlinxi + 2 * dflinxi
  bic.linxi <- devlinxi + logn * dflinxi
  sse.linxi <- fit.lin.xi$SSE
  dev.roy.linxi <- dev.gaussian(RSS = sse.linxi, weights = weights, n = N)
  df.linxi <- N - (dflinxi - 1)
  # =============================================================================
  # Fit a model with linear in axi = acd(xi)--Model M5 in R&S 2016
  # =============================================================================
  axi <- acd(x = x[, xi], s = powers, shift = 0, scale = 1)$acd
  xkk <- cbind(axi, adjdata)
  colnames(xkk) <- c(xi, colnames(adjdata))
  fit.lin.axi <- model.fit(
    x = xkk, y = y, family = family, method = method,
    weights = weights, offset = offset, strata = strata,
    control = control, rownames = rownames,
    nocenter = nocenter
  )
  # Total number of parameters in the fitted model.
  dflinaxi <- fit.lin.axi$df
  # Deviance, AIC, BIC and SSE of the null model
  devlinaxi <- -2 * fit.lin.axi$logl
  aic.linaxi <- devlinaxi + 2 * dflinaxi
  bic.linaxi <- devlinaxi + logn * dflinaxi
  sse.linaxi <- fit.lin.axi$SSE
  dev.roy.linaxi <- dev.gaussian(RSS = sse.linaxi, weights = weights, n = N)
  # subtract 1 = scale parameter and 1 = FP used for acd calculation
  df.linaxi <- N - (dflinaxi - 1)
  # =============================================================================
  # Fit best FP1 model to xi--Model M2 in R&S 2016
  # =============================================================================
  # change acdx for xi to temporarily false so that we  fit bestFP1 for xi
  # while adjusting other variables
  acdxi <- replace(acdx, which(names(acdx) %in% xi), F)
  fit.fp1.xi <- bestfp1(
    y = y, x = x, xi = xi, allpowers = allpowers,
    powers = powers, family = family, method = method,
    weights = weights, offset = offset, strata = strata,
    control = control, rownames = rownames,
    nocenter = nocenter, acdx = acdxi
  )
  # Deviance, AIC, BIC and SSE
  devfp1xi <- fit.fp1.xi$dev.all[3]
  dev.roy.fp1xi <- fit.fp1.xi$dev.roy.all[3]
  aic.fp1xi <- fit.fp1.xi$aic.all[3]
  bic.fp1xi <- fit.fp1.xi$bic.all[3]
  sse.fp1xi <- fit.fp1.xi$sse.all[3]
  df.fp1xi <- fit.fp1.xi$df.all[3]
  # best FP1 function (dev,aic, bic, sse)
  bestfpxi <- fit.fp1.xi$fn.bestfp1
  # =============================================================================
  # Fit best FP1 model to axi = acd(xi)--Model M3 in R&S 2016
  # =============================================================================
  # replace the column of xi with acd(xi) and estimate the best fp for new xi
  # set acdx = F so that 8 fp variables will be generated for the new xi
  xx <- x
  xx[, which(colnames(x) == xi)] <- axi
  fit.fp1.axi <- bestfp1(
    y = y, x = xx, xi = xi, allpowers = allpowers,
    powers = powers, family = family, method = method,
    weights = weights, offset = offset, strata = strata,
    control = control, rownames = rownames,
    nocenter = nocenter, acdx = acdxi
  )

  # Deviance, AIC, BIC and SSE of the best FP1 model
  devfp1axi <- fit.fp1.axi$dev.all[3] # best fp1 dev in position 3
  dev.roy.fp1axi <- fit.fp1.axi$dev.roy.all[3]
  aic.fp1axi <- fit.fp1.axi$aic.all[3]
  bic.fp1axi <- fit.fp1.axi$bic.all[3]
  sse.fp1axi <- fit.fp1.axi$sse.all[3]

  df.fp1axi <- fit.fp1.axi$df.all[3]
  # best FP1 function (dev,aic, bic, sse)
  bestfpaxi <- fit.fp1.axi$fn.bestfp1
  # =============================================================================
  # Fit 64 linear models for xi and axi of interest while adjusting for other
  # variables. Model M1 in R&S 2016
  # =============================================================================
  devs <- dev.roy <- sse <- aic <- bic <- dfp1 <- numeric(nv)
  for (i in seq_len(nv)) {
    # combine each FP1(p1,p2) variable for x of interest with adjustment variables
    xout <- cbind(fpdata[[i]], adjdata)
    colnames(xout) <- c("newx1", "newx2", colnames(adjdata))
    # Fit the  model which can be a glm or cox depending on the family chosen
    fit1 <- model.fit(
      x = xout, y = y, family = family, method = method,
      weights = weights, offset = offset, strata = strata,
      control = control, rownames = rownames,
      nocenter = nocenter
    )
    # save degrees of freedom from a model fit. regression coefficients and
    # total estimated FP powers. 2 is add because of FP power of xi and acd(xi)
    dfp1[i] <- (fit1$df + 2) - 1 # 1 is scaled parameter in df
    # Deviance, AIC, BIC and SSE.
    devs[i] <- -2 * fit1$logl
    aic[i] <- devs[i] + 2 * (dfp1[i])
    bic[i] <- devs[i] + logn * (dfp1[i])
    sse[i] <- fit1$SSE
    dev.roy[i] <- dev.gaussian(RSS = sse[i], weights = weights, n = N)
  }
  # Best FP1(p1,p2) function based on dev, aic, bic and sse
  s <- df1$powers
  fn.bestfp1 <- list(
    dev = s[which.min(devs), ], aic = s[which.min(aic), ],
    bic = s[which.min(bic), ], sse = s[which.min(sse), ],
    dev.roy = s[which.min(dev.roy), ]
  )
  # =============================================================================
  # combine deviances, aic, bic, sse etc
  # order: M6=NULL, M4=linear(xi), M2=FP1(xi), M3=FP1(xia), M1=FP1(xi, xia),
  #        M5=linear(xia)
  # =============================================================================
  dev.all <- c(devnull, devlinxi, devfp1xi, devfp1axi, devs[which.min(devs)], devlinaxi)
  dev.roy.all <- c(dev.roy.null, dev.roy.linxi, dev.roy.fp1xi, dev.roy.fp1axi, dev.roy[which.min(dev.roy)], dev.roy.linaxi)
  aic.all <- c(aic.null, aic.linxi, aic.fp1xi, aic.fp1axi, aic[which.min(aic)], aic.linaxi)
  bic.all <- c(bic.null, bic.linxi, bic.fp1xi, bic.fp1axi, bic[which.min(bic)], bic.linaxi)
  sse.all <- c(sse.null, sse.linxi, sse.fp1xi, sse.fp1axi, sse[which.min(sse)], sse.linaxi)
  df.all <- c(df.null, df.linxi, df.fp1xi, df.fp1axi, N - dfp1[which.min(sse)], df.linaxi) # subtract scale parameter
  names(dev.all) <- names(aic.all) <- names(bic.all) <- names(sse.all) <- names(df.all) <- c("NULL", "Linearxi", "FP1xi", "FP1xa", "FP1xixia", "Linearxia")
  # The functions are in pairs because of (xi,xia). If xi is eliminated and xia is
  # selected then we have (NA,p2), if both are eliminated we have (NA,NA) etc
  all.funs <- vector(mode = "list", length = 5)
  names(all.funs) <- c("dev", "aic", "bic", "sse", "dev.roy")
  for (i in 1:5) {
    all.funs[[i]] <- list(
      null = c(NA, NA), linxi = c(1, NA), FP1xi = c(bestfpxi[i], NA),
      FP1axi = c(NA, bestfpaxi[1]), FP1xiaxi = fn.bestfp1[[i]], linaxi = c(NA, 1)
    )
  }
  outs <- list(
    dev.all = dev.all, aic.all = aic.all, bic.all = bic.all,
    sse.all = sse.all, df.all = df.all, all.funs = all.funs, dev.roy.all = dev.roy.all
  )
  return(outs)
}
