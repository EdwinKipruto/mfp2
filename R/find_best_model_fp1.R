# A function that select the best model between a null and a linear model for
#  x of interest.Only suitable for variables with df = 1
find_best_model_fp1 <- function(x, y, xi, allpowers, df, weights, offset, control,
                                family, criterion, select, alpha, keep, powers, method,
                                strata, ftest, rownames, nocenter, verbose, acdx) {
  #  Set df = 1 in the extract_adjustment_data() because linearity is assumed and fpdata
  # would be original x of interest (untransformed)
  xadjv <- extract_adjustment_data(
    x = x, xi = xi, allpowers = allpowers, df = 1,
    powers = powers, acdx = acdx
  ) # degree, s and scale does not play any role here
  # Number of observations and a quantity for BIC calculation
  N <- dim(x)[1L]
  logn <- log(N)
  # Fit a null model, which is a model that has no x of interest but only adjustment variables.
  modelnull <- fit_model(
    x = xadjv$adjdata, y = y, family = family, weights = weights,
    offset = offset, method = method, strata = strata,
    control = control, rownames = rownames,
    nocenter = nocenter
  )
  # The total number of parameters in the model. The scale parameter is included
  # in Gaussian. FP powers in the models is not considered since its not estimated
  dfnull <- modelnull$df
  # Deviance, AIC, BIC and SSE for a null model
  devnull <- -2 * modelnull$logl
  aic.null <- devnull + 2 * dfnull
  bic.null <- devnull + logn * dfnull
  sse.null <- modelnull$SSE
  dev.roy.null <- deviance_gaussian(RSS = sse.null, weights = weights, n = N)
  # Residual degrees of freedom. It includes scale parameter. for consistency
  # with stata we get rid of it
  df.null.res <- N - (dfnull - 1)
  # Fit a model based on the assumption that x is linear while accounting for other variables.
  xx <- cbind(xadjv$fpdata, xadjv$adjdata)
  colnames(xx) <- c("xnew", colnames(xadjv$adjdata))
  modellinear <- fit_model(
    x = xx, y = y, family = family, weights = weights,
    offset = offset, method = method, strata = strata,
    control = control, rownames = rownames,
    nocenter = nocenter
  )
  # Deviance, AIC, BIC and SSE for a linear model
  devlin <- -2 * modellinear$logl
  dflin <- modellinear$df # includes scale parameter, which is fine since glm uses it in aic or bic calculation
  aic.lin <- devlin + 2 * dflin
  bic.lin <- devlin + logn * dflin
  sse.linear <- modellinear$SSE
  dev.roy.lin <- deviance_gaussian(RSS = sse.linear, weights = weights, n = N)
  # Residual degrees of freedom. It includes scale parameter. for consistency
  # with stata mfpa results we get rid of it
  df.linear.res <- N - (dflin - 1)
  # Combine the deviance, AIC,BIC and SSE for a null and linear model
  dev.all <- c(devnull, devlin)
  dev.roy.all <- c(dev.roy.null, dev.roy.lin)
  aic.all <- c(aic.null, aic.lin)
  bic.all <- c(bic.null, bic.lin)
  sse.all <- c(sse.null, sse.linear)
  df.all <- c(df.null.res, df.linear.res)
  names(dev.all) <- names(dev.roy.all) <- names(aic.all) <- names(bic.all) <- names(sse.all) <- names(df.all) <- c("null", "linear")
  # Deviance difference
  dev.diff <- as.numeric(devnull - devlin)
  # Model selection based on AIC, BIC or P-values
  if (criterion == "AIC") {
    index.bestmodel <- which.min(aic.all) # 1 = variable removed, 2 = linear
    # keep xi in the model
    if (xi %in% keep) {
      index.bestmodel <- which.min(aic.all[-1]) + 1
    }
  } else if (criterion == "BIC") {
    index.bestmodel <- which.min(bic.all)
    if (xi %in% keep) {
      index.bestmodel <- which.min(bic.all[-1]) + 1
    }
  } else { # criterion=="pvalue"
    if (ftest) {
      stats <- calculcate_f_statistic_stata(
        dev.reduced = dev.roy.all[1], dev.full = dev.roy.all[2], d1 = 1,
        d2 = df.all[2], n = N
      )
      pvalue <- stats$pval
      fstatistic <- stats$fstatistic
      dev.diff <- stats$dev.diff
    } else {
      pvalue <- pchisq(q = dev.diff, df = dflin - dfnull, lower.tail = F)
    }
    names(pvalue) <- names(dev.diff) <- c("Null vs Linear")
    index.bestmodel <- ifelse(pvalue > select, 1, 2)
  }
  if (verbose) {
    if (criterion == "pvalue") {
      if (ftest) {
          print_mfp_summary_2(
          namex = xi, dev.all = dev.roy.all, df.res = df.all, df.den = df.all, dev.diff = dev.diff, f = fstatistic,
          pvalues = pvalue, best.function = list(1), index.bestmodel = index.bestmodel, acd = acdx[xi]
        ) # acdx[xi]
      } else {
          print_mfp_summary_1(
          namex = xi, dev.all = dev.all, dev.diff = dev.diff,
          pvalues = pvalue,
          index.bestmodel = index.bestmodel,
          best.function = list(1), acd = acdx[xi]
        ) # acdx[xi]
      }
    } else {
      switch(criterion,
        "AIC" = print_mfp_summary_3(xi, gic = aic.all, keep = keep, best.function = list(1), acd = F),
        "BIC" = print_mfp_summary_3(xi, gic = bic.all, keep = keep, best.function = list(1), acd = F)
      )
    }
  }

  fit <- list(
    overall.best.fn = as.numeric(ifelse(index.bestmodel == 1, NA, 1)),
    dev.all = if (ftest) {
      dev.roy.all
    } else {
      dev.all
    }, aic.all = aic.all,
    bic.all = bic.all, sse.all = sse.all,
    index.bestmodel = as.numeric(index.bestmodel),
    dev.diff = dev.diff,
    pvalues = if (criterion == "pvalue") {
      pvalue
    } else {
      NA
    },
    fstatistic = ifelse(ftest && criterion == "pvalue", fstatistic, NA),
    df.all = df.all
  )
  return(fit)
}
