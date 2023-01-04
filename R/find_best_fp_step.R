#' Function to estimate the best FP functions for a variable
#' 
#' See [mfpa()] for a brief summary on the notation used here. 
#' 
#' @details 
#' The function selection procedure (FSP) is used if the p-value criterion is 
#' chosen, whereas the criteria AIC and BIC select the model with the smallest 
#' AIC and BIC, respectively.
find_best_fp_step <- function(x, y, xi, allpowers, df, weights, offset, family,
                            criterion, select, init, alpha, keep, powers, method, strata,
                            ftest, control, rownames, nocenter, verbose, acdx) {
  N <- dim(x)[1L]
  # If df = 1 then we just fit usual linear models and test: NULL vs Linear
  if (df == 1) {
    fits <- find_best_fp_step1(
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
      bfpa <- find_best_acd_fp1(
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
          mstats <- calculate_f_test_royston(dev = dev.roy.all, resid.df = df.all, n = N, acd = T)
          pvalue <- mstats$pvalues
          dev.diff <- mstats$dev.diff
          fstatistic <- mstats$fstatistic
        } else {
          mstats <- calculate_chisquare_test(dev = dev.all, acd = T) # can we use dev.all or dev.roy.all? for gaussian
          pvalue <- mstats$pvalues
          dev.diff <- mstats$dev.diff
        }
        # Functions are ordered like: 1 = null, 2 = Lin(xi), 3 = FP1(xi),
        # 4 = FP1(axi), 5 = FP1(xi, axi) and 6 = lin(axi)
        index.bestmodel <- find_index_best_model_acd(pvalue = pvalue, 
                                                     select = select, 
                                                     alpha = alpha)
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
              print_mfp_summary_2(
              namex = xi, dev.all = dev.roy.all, df.res = df.all, dev.diff = dev.diff, f = fstatistic,
              df.den = df.all, pvalues = pvalue, best.function = bestfuns[[5]],
              index.bestmodel = index.bestmodel, acd = T
            )
          } else {
              print_mfp_summary_1(
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
            "AIC" = print_mfp_summary_3(xi, gic = aic.all, keep = keep, best.function = bestfuns[[2]], acd = T),
            "BIC" = print_mfp_summary_3(xi, gic = bic.all, keep = keep, best.function = bestfuns[[3]], acd = T)
          )
        }
      }
      #----------------------------usual mfp procedure------------------------
    } else {
      # this part is needed for FPm if degree>2
      bfp1 <- find_best_fp1(
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
            modelparms <- calculate_f_test_royston(dev = dev.roy.all, resid.df = df.all, n = N, acd = F)
            fstatistic <- modelparms$fstatistic
          } else {
            modelparms <- calculate_chisquare_test(dev.all, acd = F)
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
        fpx <- calculate_metrics_fpm(
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
            modelparms <- calculate_f_test_royston(dev = dev.roy.all, resid.df = df.all, n = N, acd = F)
            fstatistic <- modelparms$fstatistic
          } else {
            modelparms <- calculate_chisquare_test(dev.all, acd = F)
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
              print_mfp_summary_2(
              namex = xi, dev.all = dev.roy.all, df.res = df.all, dev.diff = dev.diff, f = fstatistic,
              df.den = df.all, pvalues = pvalue, best.function = best.function1,
              index.bestmodel = index.bestmodel, acd = F
            )
          } else {
              print_mfp_summary_1(
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
            "AIC" = print_mfp_summary_3(xi, gic = aic.all, keep = keep, best.function = best.function.aic, acd = F),
            "BIC" = print_mfp_summary_3(xi, gic = bic.all, keep = keep, best.function = best.function.bic, acd = F)
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

#' Helper to find best model when acd transformation is desired
#' 
#' To be used in [find_best_fp_step()].
#' 
#' @param pvalue vector of pvalues of: Null vs M1, lin(xi) vs M1, fp1(xi) vs M1, 
#' FP1(acd(xi)) vs M1 and lin(acd(xi)) vs FP1(acd(xi)) in that order. 
#' Note M1 = FP1(xi, acd(xi)).
find_index_best_model_acd <- function(pvalue, 
                                      select, 
                                      alpha) {
  if (pvalue[1] > select) {
    index.bestmodel <- 1 # Null vs M1: if not sig then NULL model is chosen
  } else {
    # Linear(xi) vs M1: if not sig then choose linear xi is chosen
    if (pvalue[2] > alpha) {
      index.bestmodel <- 2
    } else {
      # FP1(xi) vs M1: if not sig then choose FP1(xi)
      if (pvalue[3] > alpha) {
        index.bestmodel <- 3
      } else {
        # FP1(acd(xi)) vs M1: if not sig then choose FP1(acd(xi))
        if (pvalue[4] > alpha) {
          # since the FP1(axi) vs M1 is not sig, we can check whether linear(axi)
          # is better than FP1(axi)
          if (pvalue[5] > alpha) {
            index.bestmodel <- 6
          } else {
            index.bestmodel <- 4
          }
        } else {
          index.bestmodel <- 5
        }
      }
    }
  }
  
  index.bestmodel
}

#' Helper to calculate metrics for models
#' 
#' To be used in [find_best_fp_step()].
#' 
#' @details 
#' If the maximum allowed degree is 5, this function will calculate metrics for
#' FP2 through FP5, which will then be combined with metrics for FP1 estimated 
#' by [find_best_fp1()]. The linear function makes the formula more complicated, 
#' so there are separate functions for FP1 and FPm.
calculate_metrics_fpm <- function(y, 
                                  x, 
                                  xi, 
                                  allpowers, 
                                  powers, 
                                  family, 
                                  weights, 
                                  offset, 
                                  strata,
                                  control,
                                  method,
                                  rownames, 
                                  nocenter, 
                                  degree, 
                                  acdx) {
  # we output best fpm parameters from degree 2 to m for degree 1 see bestfp1
  mm <- seq(2, degree)
  out <- vector(mode = "list", length = length(mm))
  for (k in seq_along(mm)) {
    out[[k]] <- find_best_fpm(
      y = y, x = x, xi = xi, allpowers = allpowers,
      powers = powers, family = family, weights = weights,
      offset = offset, strata = strata, control = control,
      method = method, rownames = rownames,
      nocenter = nocenter, degree = mm[k], acdx = acdx
    )
  }
  # save AIC, BIC, DEV and SSE for m = 2,3,...
  dev <- unlist(lapply(out, `[[`, 1), use.names = F) # deviances are in position 1 of nested list
  aic <- unlist(lapply(out, `[[`, 2), use.names = F) # aic are in position 2 of nested list
  bic <- unlist(lapply(out, `[[`, 3), use.names = F) # bic are in position 3 of nested list
  sse <- unlist(lapply(out, `[[`, 4), use.names = F) # sse are in position 4 of nested list
  dev.roy <- unlist(lapply(out, `[[`, 11), use.names = F) # royston deviances are in position 11 of nested list
  
  # degree of freedom: regression coefficients plus estimated FP powers
  df.best.fpm.sse <- unlist(lapply(out, `[[`, 9), use.names = F) # df are in position 9 of nested list
  names(dev) <- names(dev.roy) <- names(aic) <- names(bic) <- names(sse) <- names(df.best.fpm.sse) <- paste0("FP", mm)
  # save corresponding fp powers
  fn.bestfpm <- lapply(out, `[[`, 10) # functions are in position 10 of nested list
  fun.dev <- lapply(fn.bestfpm, `[[`, 1)
  fun.aic <- lapply(fn.bestfpm, `[[`, 2)
  fun.bic <- lapply(fn.bestfpm, `[[`, 3)
  fun.sse <- lapply(fn.bestfpm, `[[`, 4)
  fun.dev.roy <- lapply(fn.bestfpm, `[[`, 5)
  
  list(
    dev = dev, aic = aic, bic = bic, sse = sse, df.best.fpm.sse = df.best.fpm.sse,
    fun = list(dev = fun.dev, aic = fun.aic, bic = fun.bic, sse = fun.sse, dev.roy = fun.dev.roy),
    dev.roy = dev.roy
  )
}
