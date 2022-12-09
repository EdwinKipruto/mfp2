# A function that estimates the best FP functions for x of interest. The
# function selection procedure (FSP) is used if the p-value criterion is chosen,
# whereas the criteria AIC and BIC select the model with the smallest AIC and BIC,
# respectively.
overall.best.fp <- function(x, y, xi, allpowers, df, weights, offset, family,
                            criterion, select, init, alpha, keep, powers, method, strata,
                            ftest, control, rownames, nocenter, verbose, acdx) {
  N <- dim(x)[1L]
  # If df = 1 then we just fit usual linear models and test: NULL vs Linear
  if (df == 1) {
    fits <- overall.best.linear(
      x = x, y = y, xi = xi, allpowers = allpowers,
      weights = weights, offset = offset, family = family,
      criterion = criterion, select = select, alpha = alpha,
      keep = keep, powers = powers, method = method,
      strata = strata, ftest = ftest, control = control,
      rownames = rownames, nocenter = nocenter,
      verbose = verbose,
      acdx = acdx
    )
    # The best power estimated. Is NA if Null model was selected and 1 if linear
    best.fp.power <- fits$overall.best.fn
    # Deviance, aic and bic for null and linear model
    dev.all <- fits$dev.all
    aic.all <- fits$aic.all
    bic.all <- fits$bic.all
    # index of the best model. 1 = Null model while 2 = linear model
    index.bestmodel <- fits$index.bestmodel
    # Deviance difference and corresponding pvalue between NULL and linear model
    dev.diff <- fits$dev.diff
    pvalue <- fits$pvalues
    fit <- list(
      overall.best.fn = as.numeric(best.fp.power), dev.all = dev.all, aic.all = aic.all,
      bic.all = bic.all, index.bestmodel = index.bestmodel, dev.diff = dev.diff, pvalues = pvalue
    )
  } else {
    if (acdx[xi]) {
      # compute deviances, aic, bic and sse for model M1-M6
      bfpa <- bestfpa(
        y = y, x = x, xi = xi, allpowers = allpowers, powers = powers, family = family,
        method = method, weights = weights, offset = offset,
        strata = strata, control = control, rownames = rownames,
        nocenter = nocenter, acdx = acdx
      )
      # Deviance, aic and bic for c(null, lin(xi), fp1(xi), fp1(acd(xi)), fp1(xi,acd(xi)), lin(acd(xi))) in that order
      dev.all <- bfpa$dev.all
      dev.roy.all <- bfpa$dev.roy.all
      aic.all <- bfpa$aic.all
      bic.all <- bfpa$bic.all
      df.all <- bfpa$df.all
      # choose the best function based on criterion
      if (criterion == "AIC") {
        # if index.bestmodel =
        # 1 then both xi and acd(xi) were removed,
        # 2 = linear(xi), meaning acd(xi) was removed
        # 3 = fp1(xi), meaning acd(xi) was removed
        # 4 = fp1(acd(xi)),meaning xi was removed
        # 5 = fp1(xi,acd(xi)), both xi and acd(xi) selected
        # 6 = lin(acd(xi))), meaning xi was removed but acd(xi) is linear
        index.bestmodel <- which.min(aic.all)
        # keep xi in the model
        if (xi %in% keep) {
          index.bestmodel <- which.min(aic.all[-1]) + 1 # we add 1 so that linear  = 2 and best fp1 = 3
        }
      } else if (criterion == "BIC") {
        index.bestmodel <- which.min(bic.all)
        if (xi %in% keep) {
          index.bestmodel <- which.min(bic.all[-1]) + 1 # we add 1 so that linear  = 2 and best fp1 = 3
        }
      } else {
        # Calculate p-values for 1). M1 vs Null  2). M1 vs M4  3). M1 vs M2
        # 4) M1 vs M3 and 5). M3 vs M5
        if (ftest) {
          mstats <- ftest.royston(dev = dev.roy.all, resid.df = df.all, n = N, acd = T)
          pvalue <- mstats$pvalues
          dev.diff <- mstats$dev.diff
          fstatistic <- mstats$fstatistic
        } else {
          mstats <- model.pvalues.chi(dev = dev.all, acd = T) # can we use dev.all or dev.roy.all? for gaussian
          pvalue <- mstats$pvalues
          dev.diff <- mstats$dev.diff
        }
        # Functions are ordered like: 1 = null, 2 = Lin(xi), 3 = FP1(xi),
        # 4 = FP1(axi), 5 = FP1(xi, axi) and 6 = lin(axi)
        index.bestmodel <- best.model.index.acd(pvalue = pvalue, select = select, alpha = alpha)
      }
      # all best functions are ordered in list. Each element of a list contains:
      # c(NULL, lin(xi),fp1(xi), fp1(acd(xi)), fp1(xi,acd(xi)), lin(acd(xi)))
      # where the element of a nested list is dev, aic, bic and sse in that order
      bestfuns <- bfpa$all.funs
      if (criterion == "pvalue") {
        if (ftest) {
          best.fp.power <- bestfuns[[5]][[index.bestmodel]]
        } else {
          best.fp.power <- bestfuns[[1]][[index.bestmodel]]
        }
      } else {
        best.fp.power <- switch(criterion,
          "AIC" = bestfuns[[2]][[index.bestmodel]],
          "BIC" = bestfuns[[3]][[index.bestmodel]]
        )
      }
      # Printing on the screen
      if (verbose) {
        if (criterion == "pvalue") {
          if (ftest) {
            fprint2(
              namex = xi, dev.all = dev.roy.all, df.res = df.all, dev.diff = dev.diff, f = fstatistic,
              df.den = df.all, pvalues = pvalue, best.function = bestfuns[[5]],
              index.bestmodel = index.bestmodel, acd = T
            )
          } else {
            fprint1(
              namex = xi,
              dev.all = dev.all,
              dev.diff = dev.diff,
              pvalues = pvalue,
              index.bestmodel = index.bestmodel,
              best.function = bestfuns[[1]], acd = T
            )
          }
        } else {
          switch(criterion,
            "AIC" = fprint3(xi, gic = aic.all, keep = keep, best.function = bestfuns[[2]], acd = T),
            "BIC" = fprint3(xi, gic = bic.all, keep = keep, best.function = bestfuns[[3]], acd = T)
          )
        }
      }
      #----------------------------usual mfp procedure------------------------
    } else {
      # this part is needed for FPm if degree>2
      bfp1 <- bestfp1(
        y = y, x = x, xi = xi, allpowers = allpowers, powers = powers, family = family,
        method = method, weights = weights, offset = offset,
        strata = strata, control = control, rownames = rownames,
        nocenter = nocenter, acdx = acdx
      )
      # Deviance, aic, bic and sse of a null model, linear and best FP1 function
      dev.all <- bfp1$dev.all # based on loglikelihood
      dev.roy.all <- bfp1$dev.roy.all # based on royston formula
      aic.all <- bfp1$aic.all
      bic.all <- bfp1$bic.all
      sse.all <- bfp1$sse.all
      # degrees of freedom of sse-important for f statistic
      df.all <- bfp1$df.all
      # A vector of best FP1 function selected based on deviance, aic, bic, sse and dev.roy
      # in that order i.e bestfp1 = c(dev.fun =,...,sse.fun=, dev.roy.fun =  )
      bestfp1x <- bfp1$fn.bestfp1
      # calculate the maximum permitted degree. we know that df = 2m where m is
      # the degree e.g df = 2 = 2(1) is fp1,  df = 4 = 2(2) is fp2, df = 6 = 2(3) is fp3 etc.
      degree <- df / 2
      # If degree = 1, calculate p-values for null vs. linear and linear vs. FP1
      if (degree == 1) {
        # Choose best model between null, linear and best FP1 based on AIC or BIC.
        if (criterion == "AIC") {
          # if index.bestmodel = 1 variable was removed, 2 = linear and 3 is best fp1
          index.bestmodel <- which.min(aic.all)
          # If the user wants to force some variables into the model then we should
          # compare only the aic for linear and other fp1 and update index of best model
          if (xi %in% keep) {
            index.bestmodel <- which.min(aic.all[-1]) + 1 # we add 1 so that linear  = 2 and best fp1 = 3
          }
        } else if (criterion == "BIC") {
          index.bestmodel <- which.min(bic.all)
          if (xi %in% keep) {
            index.bestmodel <- which.min(bic.all[-1]) + 1 # we add 1 so that linear  = 2 and best fp1 = 3
          }
        } else { # p-values
          if (ftest) {
            modelparms <- ftest.royston(dev = dev.roy.all, resid.df = df.all, n = N, acd = F)
            fstatistic <- modelparms$fstatistic
          } else {
            modelparms <- model.pvalues.chi(dev.all, acd = F)
          }
          # P-values for each test-we assign names
          pvalue <- modelparms$pvalues
          dev.diff <- modelparms$dev.diff
          names(pvalue) <- names(dev.diff) <- c("Null", "Linear") # FP1 vs NULL and FP1 vs Linear
          if (pvalue[1] > select) {
            index.bestmodel <- 1
          } else {
            index.bestmodel <- ifelse(pvalue[2] > alpha, 2, 3)
          }
        }
      } else {
        # Here we fit other fpm models where m can be 2, 3, and so on
        fpx <- bestfpm.all(
          y = y, x = x, xi = xi, allpowers = allpowers, powers = powers, family = family,
          method = method, weights = weights, offset = offset,
          strata = strata, control = control, rownames = rownames,
          nocenter = nocenter, degree = degree, acdx = acdx
        )
        # best fp functions based on dev, aic, bic, sse and dev.roy. It's a list with the
        # first element being fp1, then fp2, fp3 and so on.
        bestfp.dev <- append(fpx$fun$dev, list(bestfp1x[1]), 0) # c(FP1, FP2,...)
        bestfp.dev.roy <- append(fpx$fun$dev.roy, list(bestfp1x[5]), 0) # c(FP1, FP2,...)
        bestfp.aic <- append(fpx$fun$aic, list(bestfp1x[2]), 0)
        bestfp.bic <- append(fpx$fun$bic, list(bestfp1x[3]), 0)
        bestfp.sse <- append(fpx$fun$sse, list(bestfp1x[4]), 0)
        # A vector of Deviance, AIC and BIC for all models
        dev.all <- c(dev.all, fpx$dev) # NULL, linear, FP1,....
        dev.roy.all <- c(dev.roy.all, fpx$dev.roy) # NULL, linear, FP1,....
        aic.all <- c(aic.all, fpx$aic)
        bic.all <- c(bic.all, fpx$bic)
        df.all <- c(df.all, fpx$df.best.fpm.sse)
        # select the best model based on AIC or BIC or Chi-square P-values
        if (criterion == "AIC") {
          # if index.bestmodel = 1 then the variable was removed, 2 = linear, 3= best fp1 and so on
          index.bestmodel <- which.min(aic.all)
          if (xi %in% keep) {
            index.bestmodel <- which.min(aic.all[-1]) + 1 # we add 1 so that linear  = 2 and best fp1 = 3
          }
        } else if (criterion == "BIC") {
          index.bestmodel <- which.min(bic.all)
          if (xi %in% keep) {
            index.bestmodel <- which.min(bic.all[-1]) + 1 # we add 1 so that linear  = 2 and best fp1 = 3
          }
        } else { # p-values
          if (ftest) {
            modelparms <- ftest.royston(dev = dev.roy.all, resid.df = df.all, n = N, acd = F)
            fstatistic <- modelparms$fstatistic
          } else {
            modelparms <- model.pvalues.chi(dev.all, acd = F)
          }
          pvalue <- modelparms$pvalues
          dev.diff <- modelparms$dev.diff
          dd <- length(dev.all) - 2 - 1 # subtract Null, linear and max permitted. If 0 then we have pvalues for NULL vs FP1 and Lin vs FP1
          names(pvalue) <- names(dev.diff) <- c("Null", "Linear", if (dd != 0) {
            paste0("FP", seq_len(dd))
          })
          # Variable selection
          if (pvalue[1] > select) {
            index.bestmodel <- 1
            # function selection
          } else {
            # compare pvalues for linear, fp1,...fpm with alpha
            pos <- which(pvalue[-1] > alpha) # note pvalue[1] belongs to FPM vs NULL
            # if alpha = 1 then pvalue[-1]>alpha might all be FALSE.
            index.bestmodel <- ifelse(length(pos) == 0, length(dev.all), pos[1] + 1) # We add 1 because position 1 is null model
          }
        }
      }
      # Select the best overall FP power
      if (index.bestmodel == 1) {
        best.fp.power <- NA # Variable eliminated
      } else if (index.bestmodel == 2) {
        best.fp.power <- 1 # Linear function
      } else if (index.bestmodel == 3) { # best FP1 function
        best.fp.power <- switch(criterion,
          "pvalue" = if (ftest) {
            bestfp1x[5]
          } else {
            bestfp1x[1]
          },
          "AIC" = bestfp1x[2],
          "BIC" = bestfp1x[3]
        )
      } else { # best FPm function
        best.fp.power <- switch(criterion,
          "pvalue" = if (ftest) {
            bestfp.dev.roy[[index.bestmodel - 2]]
          } else {
            bestfp.dev[[index.bestmodel - 2]]
          }, # bestfp.dev = c(FP1, FP2,...) SO subtract 2 to get correct index in bestfp.dev
          "AIC" = bestfp.aic[[index.bestmodel - 2]], # e.g FP2 = index.model = 4 but it is in position 2 in bestfp.dev
          "BIC" = bestfp.bic[[index.bestmodel - 2]]
        )
      }
      # Print results if verbose = T
      if (verbose) {
        if (criterion == "pvalue") {
          best.function1 <- if (degree == 1) {
            list(lin = 1, fpm = if (ftest) {
              bestfp1x[5]
            } else {
              bestfp1x[1]
            })
          } else {
            append(if (ftest) {
              bestfp.dev.roy
            } else {
              bestfp.dev
            }, list(1), 0)
          }
          if (ftest) {
            fprint2(
              namex = xi, dev.all = dev.roy.all, df.res = df.all, dev.diff = dev.diff, f = fstatistic,
              df.den = df.all, pvalues = pvalue, best.function = best.function1,
              index.bestmodel = index.bestmodel, acd = F
            )
          } else {
            fprint1(
              namex = xi, dev.all = dev.all, dev.diff = dev.diff,
              pvalues = pvalue, index.bestmodel = index.bestmodel,
              best.function = best.function1, acd = F
            )
          }
          # AIC and BIC display
        } else {
          best.function.aic <- if (degree == 1) {
            list(lin = 1, fpm = bestfp1x[2])
          } else {
            append(bestfp.aic, list(1), 0)
          }
          best.function.bic <- if (degree == 1) {
            list(lin = 1, fpm = bestfp1x[3])
          } else {
            append(bestfp.bic, list(1), 0)
          }
          switch(criterion,
            "AIC" = fprint3(xi, gic = aic.all, keep = keep, best.function = best.function.aic, acd = F),
            "BIC" = fprint3(xi, gic = bic.all, keep = keep, best.function = best.function.bic, acd = F)
          )
        }
      }
    }
    fit <- list(
      overall.best.fn = as.numeric(best.fp.power), dev.all = if (ftest) {
        dev.roy.all
      } else {
        dev.all
      },
      aic.all = aic.all,
      bic.all = bic.all, index.bestmodel = index.bestmodel,
      dev.diff = if (criterion == "pvalue") {
        dev.diff
      } else {
        NA
      },
      pvalues = if (criterion == "pvalue") {
        pvalue
      } else {
        NA
      }
    )
  }
  return(fit)
}
