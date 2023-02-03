#' Function to estimate the best FP functions for a single variable
#' 
#' See [mfpa()] for a brief summary on the notation used here and 
#' [fit_mfp()] for an overview of the fitting procedure.  
#' 
#' @param x an input matrix of dimensions nobs x nvars. Does not contain 
#' intercept, but columns are already expanded into dummy variables as 
#' necessary. Data are assumed to be shifted and scaled. 
#' @param y a vector for the response variable or a `Surv` object.
#' @param xi a character string indicating the name of the current variable 
#' of interest, for which the best fractional polynomial transformation is
#' to be found in the current step. 
#' @param weights a vector of observation weights of length nobs. 
#' @param offset a vector of length nobs of offsets.
#' @param df a numeric vector indicating the maximum degrees of freedom for the 
#' variable of interest `xi`.
#' @param powers_current a list of length equal to the number of variables, 
#' indicating the fp powers to be used in the current step for all variables 
#' (except `xi`). 
#' @param family a character string representing a family object.
#' @param criterion a character string defining the criterion used to select 
#' variables and FP models of different degrees.
#' @param select a numeric value indicating the significance level
#' for backward elimination of `xi`.
#' @param alpha a numeric value indicating the significance level
#' for tests between FP models of different degrees for `xi`. 
#' @param keep a character vector that with names of variables to be kept 
#' in the model. 
#' @param powers a numeric vector that sets the permitted FP powers for all 
#' covariates.
#' @param method a character string specifying the method for tie handling in 
#' Cox regression.
#' @param strata a factor of all possible combinations of stratification 
#' variables. Returned from [survival::strata()]. 
#' @param nocenter a numeric vector with a list of values for fitting Cox 
#' models. See [survival::coxph()] for details.
#' @param acdx a logical vector of length nvars indicating continuous variables 
#' to undergo the approximate cumulative distribution (ACD) transformation.
#' @param ftest a logical indicating the use of the F-test for Gaussian models.
#' @param control a list with parameters for model fit.
#' @param rownames a parameter for Cox models.
#' @param verbose a logical; run in verbose mode.
#' 
#' @details 
#' The function selection procedure (FSP) is used if the p-value criterion is 
#' chosen, whereas the criteria AIC and BIC select the model with the smallest 
#' AIC and BIC, respectively.
#' 
#' It uses transformations for all other variables to assess the FP form of 
#' the current variable of interest. This function covers three main use cases: 
#' 
#' * the linear case (`df = 1`) to test between null and linear models (see
#' [find_best_linear_step()]).
#' * the case that an acd transformation is requested (`acdx` is `TRUE` 
#' for `xi`) for the variable of interest (see [find_best_acd_step()]).
#' * the (usual) case of the normal mfp algorithm to assess non-linear 
#' functional forms (see [find_best_fp1_step()] and [find_best_fpm_step()]). 
#' 
#' Note that these cases do not encompass the setting that a variable is not
#' selected, because the evaluation is done for each variable in each cycle.
#' A variable which was de-selected in earlier cycles may be added to the 
#' working model again. Also see [find_best_fp_cycle()].
#' 
#' The adjustment in each step uses the current fp powers given in 
#' `powers_current` for all other variables to determine the adjustment set 
#' and transformations in the  working model.
#' 
#' Note that the algorithm starts by setting all `df = 1`, and higher fps
#' are evaluated in turn starting from the first step in the first cycle.
#' 
#' @return 
#' A numeric vector indicating the best powers for `xi`. Entries can be 
#' `NA` if variable is to be removed from the working model. 
find_best_fp_step <- function(x,
                              y, 
                              xi,
                              weights, 
                              offset, 
                              df, 
                              powers_current, 
                              family,
                              criterion, 
                              select, 
                              alpha,
                              keep,
                              powers, 
                              method, 
                              strata,
                              nocenter, 
                              acdx, 
                              ftest, 
                              control,
                              rownames, 
                              verbose) {
  N <- dim(x)[1L]
  
  if (df == 1) {
    
    # linear case --------------------------------------------------------------
    # if df = 1 then we just fit usual linear models and test: NULL vs Linear
    fit <- find_best_linear_step(
      x = x, y = y, xi = xi, powers_current = powers_current,
      weights = weights, offset = offset, family = family,
      criterion = criterion, select = select, alpha = alpha,
      keep = keep, powers = powers, method = method,
      strata = strata, ftest = ftest, control = control,
      rownames = rownames, nocenter = nocenter,
      verbose = verbose,
      acdx = acdx
    )

    power_best = as.numeric(fit$power_best)
  } else if (acdx[xi]) {
    # acd case ---------------------------------------------------------------
    # compute deviances, aic, bic and sse for model M1-M6
    bfpa <- find_best_acd_step(
      y = y, x = x, xi = xi, powers_current = powers_current, powers = powers, family = family,
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
    
    power_best <- as.numeric(best.fp.power)
    
  } else {
    # usual mfp case ---------------------------------------------------------
    # this part is needed for FPm if degree>2
    bfp1 <- find_best_fp1_step(
      y = y, x = x, xi = xi, powers_current = powers_current, powers = powers, family = family,
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
        y = y, x = x, xi = xi, powers_current = powers_current, powers = powers, family = family,
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
    
    power_best <- as.numeric(best.fp.power)
  }
  
  power_best
}

#' Functions to find the best FP functions for a single variable
#' 
#' Handles the FP1 (`find_best_fp1_step`) and the higher order FP 
#' (`find_best_fpm_step`) cases. For parameter definitions, see
#' [find_best_fp_step()].
#' 
#' @return 
#' A list with several components giving the best power found and 
#' performance indices.
find_best_fp1_step <- function(y, 
                               x, 
                               xi, 
                               powers_current, 
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
  df1 <- transform_data_step(
    x = x, xi = xi, powers_current = powers_current, df = 2, acdx = acdx,
    powers = powers
  )
  # Matrix of adjustment variables
  adjdata <- df1$data_adj
  # List of FP1 data.
  fpdata <- df1$data_fp
  # N = number of observation and log(n) for bic calculation
  N <- nrow(x)
  logn <- log(N)
  
  # Fit a null model-model without x of interest
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
  sse.null <- fitnull$sse
  dev.roy.null <- deviance_stata(rss = sse.null, weights = weights, n = N)
  df.sse.null <- N - (dfnull - 1) # subtract scale parameter
  # df.sse.null2 <- N-fitnull$df
  
  # Fit 8 linear models for x of interest while adjusting for other variables.
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
    # sse relevant for a Gaussian model
    sse[i] <- fit1$sse
    dev.roy[i] <- deviance_stata(rss = sse[i], weights = weights, n = N)
  }
  # position of a linear function in set s: Fails when powers does not include 1
  lin.pos <- match(1, powers)
  # Deviance, AIC, BIC and sse of a linear function
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
  # dev.gaus.royston <- unlist(lapply(sse.all, function(x) deviance_stata(rss = x, weights = weights, n = N)))
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
  
  outx
}

#' @describeIn find_best_fp1_step Find higher order FP functions.
find_best_fpm_step <- function(y, 
                               x, 
                               xi, 
                               powers_current, 
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
  
  if (degree < 2) stop("Degree must be >= 2. Here we are interest on best FPm where
                    m>=2. For m = 1 see find_best_fp1_step()")
  # Generate FP data for x of interest (xi) and adjustment variables
  m <- degree
  df1 <- transform_data_step(
    x = x, xi = xi, powers_current = powers_current,
    df = 2 * m, powers = powers, acdx = acdx
  )
  # Matrix of adjustment data
  adjdata <- df1$data_adj
  # List of FP data for xi (continuous) of interest.
  fpdata <- df1$data_fp
  # Total FP powers different from 1 estimated in adjustment model
  # tFP <- calculate_number_fp_powers(df1$adjustpowers)
  # The length of generated FP variables for x of interest.
  nv <- length(fpdata)
  # Generate all possible FPm powers
  fpmpowers <- generate_powers_fp(degree = m, powers = powers)
  # log(n) for bic calculation
  n <- nrow(x)
  logn <- log(n)
  devs <- devs.royston <- sse <- aic <- bic <- dfpm <- dfx <- numeric(nv)
  xnames <- paste0("newx", seq_along(1:m))
  for (i in seq_len(nv)) {
    # combine FP variables for x of interest with adjustment variables
    xout <- cbind(fpdata[[i]], adjdata)
    colnames(xout)[1:m] <- xnames
    fit1 <- fit_model(
      x = xout, y = y, family = family, method = method,
      weights = weights, offset = offset, strata = strata,
      control = control, rownames = rownames,
      nocenter = nocenter
    )
    # Deviance of the fitted model
    devs[i] <- -2 * fit1$logl
    # number of regression coefficients (plus variance if gaussian)
    dfx[i] <- fit1$df
    # save degrees of freedom of sse for fpm i.e n-#parameters in fpm model
    # useful for F test. we add m because of m fp powers and tFP because of
    # estimated FP powers in adjustment model
    # dfpm[i] = n-(fit1$df + m + tFP)
    dfpm[i] <- n - ((fit1$df - 1) + m) # subtract scale parameter
    # AIC and BIC of the fitted model. Add m because of the m FPm powers.
    # aic[i] = devs[i] + 2*(fit1$df + m + tFP) #
    # bic[i] = devs[i] + logn*(fit1$df + m + tFP)
    aic[i] <- devs[i] + 2 * (fit1$df + m) #
    bic[i] <- devs[i] + logn * (fit1$df + m)
    # sse and deviance for gaussian family.
    sse[i] <- fit1$sse
    devs.royston[i] <- deviance_stata(rss = sse[i], weights = weights, n = dim(x)[1L])
  }
  # Best FPm function based on dev (calculated using loglik), aic, bic, sse and dev(calculated using sse)
  fn.bestfpm <- list(
    dev = fpmpowers[which.min(devs), ],
    aic = fpmpowers[which.min(aic), ],
    bic = fpmpowers[which.min(bic), ],
    sse = fpmpowers[which.min(sse), ],
    dev.r = fpmpowers[which.min(devs.royston), ]
  )
  
  mt <- list(
    dev.best.fpm = devs[which.min(devs)], # Deviance of the best FPm function
    aic.best.fpm = aic[which.min(aic)], # aic of the best FPm function
    bic.best.fpm = bic[which.min(bic)],
    sse.best.fpm = sse[which.min(sse)],
    dev.all.FPm = devs,
    aic.all.FPm = aic,
    bic.all.FPm = bic,
    sse.all.FPm = sse,
    df.best.fpm.sse = dfpm[which.min(sse)],
    fn.bestfpm = fn.bestfpm,
    dev.roy.best.fpm = devs.royston[which.min(devs.royston)]
  ) # Deviance of the best FPm function based on royston formula
  
  mt
}

#' Helper to select between null and linear term for a single variable
#' 
#' To be used in [find_best_fp_step()]. Only used if `df = 1` for a variable.
#' For parameter explanations, see [find_best_fp_step()]. All parameters 
#' captured by `...` are passed on to [fit_model()].
#' 
#' @details 
#' This function assesses a single variable of interest `xi` regarding its
#' functional form in the current working model as indicated by
#' `powers_current`, with the choice between a excluding `xi` ("null model") and
#' including a linear term ("linear fp") for `xi`.
#' 
#' @return 
#' A list with several components giving the best power found (`power_best`) and 
#' performance indices. The returned best power may be `NA`, indicating the
#' variable has been removed from the model.
find_best_linear_step <- function(x, 
                                  xi, 
                                  y, 
                                  powers_current, 
                                  criterion, 
                                  select, 
                                  alpha, 
                                  keep, 
                                  powers, 
                                  ftest, 
                                  verbose, 
                                  acdx, 
                                  ...) {
  
  n_obs <- dim(x)[1L]

  # transform all data as given by current working model
  # set variable of interest to linear term only
  x_transformed <- transform_data_step(
    x = x, xi = xi, df = 1,
    powers_current = powers_current, acdx = acdx, powers = powers
  ) 
  
  # fit null model
  # i.e. a model that does not contain xi but only adjustment variables
  model_null <- fit_model(x = x_transformed$data_adj, y = y, ...)
  
  # fit a model based on the assumption that xi is linear 
  model_linear <- fit_model(
    x = cbind(x_transformed$data_fp, x_transformed$data_adj), y = y, ...
  )
  
  # model metrics as matrix
  metrics <- rbind(
    null = calculate_model_metrics(model_null, n_obs), 
    linear = calculate_model_metrics(model_linear, n_obs)
  )  
  dev_diff <- metrics["null", "deviance_rs"] - metrics["linear", "deviance_rs"]
  
  # model selection based on AIC, BIC or P-values, or keep xi if indicated
  # 1 = variable removed, 2 = linear
  criterion = tolower(criterion)
  if (xi %in% keep) {
    model_best = 2
  } else if (grepl("aic|bic", criterion)) {
    model_best <- which.min(metrics[, criterion, drop = TRUE])
  } else { 
    # criterion == "pvalue"
     
    if (ftest) {
      stats <- calculate_f_test(
        deviances = metrics[, "deviance_stata"], 
        dfs_resid = metrics[, "df_resid"],
        n_obs = n_obs
      )
      pvalue <- stats$p_value
      fstatistic <- stats$statistic
      dev_diff <- stats$dev_diff
    } else {
      pvalue <- calculate_lr_test(metrics[, "logl"], metrics[, "df"])$pvalue
    }
    
    names(pvalue) <- c("Null vs Linear")
    names(dev_diff) <- names(pvalue)
    model_best <- ifelse(pvalue > select, 1, 2)
  }
  
  if (verbose) {
    if (criterion == "pvalue") {
      if (ftest) {
        print_mfp_summary_2(
          namex = xi, 
          dev.all = metrics[, "deviance_stata"], 
          df.res = metrics[, "df_resid"], 
          df.den = metrics[, "df_resid"], 
          dev.diff = dev_diff,
          f = fstatistic,
          pvalues = pvalue,
          best.function = list(1), 
          index.bestmodel = model_best, 
          acd = acdx[xi]
        ) 
      } else {
        print_mfp_summary_1(
          namex = xi, 
          dev.all = metrics[, "deviance_rs"], 
          dev.diff = dev_diff,
          pvalues = pvalue,
          index.bestmodel = model_best,
          best.function = list(1), 
          acd = acdx[xi]
        ) 
      }
    } else {
      print_mfp_summary_3(
        xi, 
        gic = metrics[, criterion], 
        keep = keep, 
        best.function = list(1),
        acd = FALSE
      )
    }
  }
  
  list(
    power_best = as.numeric(ifelse(model_best == 1, NA, 1)),
    dev = ifelse(ftest, metrics[, "deviance_stata"], 
                      metrics[, "deviance_rs"]), 
    aic = metrics[, "aic"],
    bic = metrics[, "bic"], 
    sse = metrics[, "sse"],
    df_resid = metrics[, "df_resid"],
    dev_diff = dev_diff,
    pvalues = ifelse(criterion == "pvalue", pvalue, NA),
    fstatistic = ifelse(ftest && criterion == "pvalue", fstatistic, NA),
    index_model_best = as.numeric(model_best)
  )
}

#' Helper to find best model involving acd transformation for a variable
#' 
#' A function that fits 64 FP1 models for x and acd(x) and returns best
#' FP1(P1,P2) function and corresponding deviance. For parameter explanations
#' see [find_best_fp_step()].
#' 
#' @return 
#' A list with several components giving the best power found and 
#' performance indices.
find_best_acd_step <- function(y, x, xi, powers_current, powers, family, method, weights,
                               offset, strata, control, rownames, nocenter, acdx) {
  # Generate FPa data for x of interest (xi). If the default FP power set is
  # used, 64 pairs of new variables are created.
  df1 <- transform_data_step(
    x = x, xi = xi, powers_current = powers_current, df = 4,
    powers = powers, acdx = acdx
  )
  # Matrix of adjustment variables
  adjdata <- df1$data_adj
  # print(head(adjdata))
  # FPa data for x of interest
  fpdata <- df1$data_fp
  nv <- length(fpdata)
  # log(n) for bic calculation
  N <- nrow(x)
  logn <- log(N)
  
  # Fit a model without xi and axi--the null model. Model M6 in R&S 2016
  fitnull <- fit_model(
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
  sse.null <- fitnull$sse
  dev.roy.null <- deviance_stata(rss = sse.null, weights = weights, n = N)
  df.null <- N - (dfnull - 1)
  
  # Fit a model with linear in xi--Model M4 in R&S 2016
  xk <- cbind(x[, xi], adjdata)
  colnames(xk) <- c(xi, colnames(adjdata))
  fit.lin.xi <- fit_model(
    x = xk, y = y, family = family, method = method,
    weights = weights, offset = offset, strata = strata,
    control = control, rownames = rownames,
    nocenter = nocenter
  )
  # Total number of parameters in the fitted model.
  dflinxi <- fit.lin.xi$df
  # Deviance, AIC, BIC and sse of the null model
  devlinxi <- -2 * fit.lin.xi$logl
  aic.linxi <- devlinxi + 2 * dflinxi
  bic.linxi <- devlinxi + logn * dflinxi
  sse.linxi <- fit.lin.xi$sse
  dev.roy.linxi <- deviance_stata(rss = sse.linxi, weights = weights, n = N)
  df.linxi <- N - (dflinxi - 1)
  
  # Fit a model with linear in axi = acd(xi)--Model M5 in R&S 2016
  axi <- fit_acd(x = x[, xi], powers = powers)$acd
  xkk <- cbind(axi, adjdata)
  colnames(xkk) <- c(xi, colnames(adjdata))
  fit.lin.axi <- fit_model(
    x = xkk, y = y, family = family, method = method,
    weights = weights, offset = offset, strata = strata,
    control = control, rownames = rownames,
    nocenter = nocenter
  )
  # Total number of parameters in the fitted model.
  dflinaxi <- fit.lin.axi$df
  # Deviance, AIC, BIC and sse of the null model
  devlinaxi <- -2 * fit.lin.axi$logl
  aic.linaxi <- devlinaxi + 2 * dflinaxi
  bic.linaxi <- devlinaxi + logn * dflinaxi
  sse.linaxi <- fit.lin.axi$sse
  dev.roy.linaxi <- deviance_stata(rss = sse.linaxi, weights = weights, n = N)
  # subtract 1 = scale parameter and 1 = FP used for acd calculation
  df.linaxi <- N - (dflinaxi - 1)
  
  # Fit best FP1 model to xi--Model M2 in R&S 2016
  # change acdx for xi to temporarily false so that we  fit bestFP1 for xi
  # while adjusting other variables
  acdxi <- replace(acdx, which(names(acdx) %in% xi), F)
  fit.fp1.xi <- find_best_fp1_step(
    y = y, x = x, xi = xi, powers_current = powers_current,
    powers = powers, family = family, method = method,
    weights = weights, offset = offset, strata = strata,
    control = control, rownames = rownames,
    nocenter = nocenter, acdx = acdxi
  )
  # Deviance, AIC, BIC and sse
  devfp1xi <- fit.fp1.xi$dev.all[3]
  dev.roy.fp1xi <- fit.fp1.xi$dev.roy.all[3]
  aic.fp1xi <- fit.fp1.xi$aic.all[3]
  bic.fp1xi <- fit.fp1.xi$bic.all[3]
  sse.fp1xi <- fit.fp1.xi$sse.all[3]
  df.fp1xi <- fit.fp1.xi$df.all[3]
  # best FP1 function (dev,aic, bic, sse)
  bestfpxi <- fit.fp1.xi$fn.bestfp1
  
  # Fit best FP1 model to axi = acd(xi)--Model M3 in R&S 2016
  # replace the column of xi with acd(xi) and estimate the best fp for new xi
  # set acdx = F so that 8 fp variables will be generated for the new xi
  xx <- x
  xx[, which(colnames(x) == xi)] <- axi
  fit.fp1.axi <- find_best_fp1_step(
    y = y, x = xx, xi = xi, powers_current = powers_current,
    powers = powers, family = family, method = method,
    weights = weights, offset = offset, strata = strata,
    control = control, rownames = rownames,
    nocenter = nocenter, acdx = acdxi
  )
  
  # Deviance, AIC, BIC and sse of the best FP1 model
  devfp1axi <- fit.fp1.axi$dev.all[3] # best fp1 dev in position 3
  dev.roy.fp1axi <- fit.fp1.axi$dev.roy.all[3]
  aic.fp1axi <- fit.fp1.axi$aic.all[3]
  bic.fp1axi <- fit.fp1.axi$bic.all[3]
  sse.fp1axi <- fit.fp1.axi$sse.all[3]
  
  df.fp1axi <- fit.fp1.axi$df.all[3]
  # best FP1 function (dev,aic, bic, sse)
  bestfpaxi <- fit.fp1.axi$fn.bestfp1
  
  # Fit 64 linear models for xi and axi of interest while adjusting for other
  # variables. Model M1 in R&S 2016
  devs <- dev.roy <- sse <- aic <- bic <- dfp1 <- numeric(nv)
  for (i in seq_len(nv)) {
    # combine each FP1(p1,p2) variable for x of interest with adjustment variables
    xout <- cbind(fpdata[[i]], adjdata)
    colnames(xout) <- c("newx1", "newx2", colnames(adjdata))
    # Fit the  model which can be a glm or cox depending on the family chosen
    fit1 <- fit_model(
      x = xout, y = y, family = family, method = method,
      weights = weights, offset = offset, strata = strata,
      control = control, rownames = rownames,
      nocenter = nocenter
    )
    # save degrees of freedom from a model fit. regression coefficients and
    # total estimated FP powers. 2 is add because of FP power of xi and acd(xi)
    dfp1[i] <- (fit1$df + 2) - 1 # 1 is scaled parameter in df
    # Deviance, AIC, BIC and sse.
    devs[i] <- -2 * fit1$logl
    aic[i] <- devs[i] + 2 * (dfp1[i])
    bic[i] <- devs[i] + logn * (dfp1[i])
    sse[i] <- fit1$sse
    dev.roy[i] <- deviance_stata(rss = sse[i], weights = weights, n = N)
  }
  # Best FP1(p1,p2) function based on dev, aic, bic and sse
  s <- df1$powers_fp
  fn.bestfp1 <- list(
    dev = s[which.min(devs), ], aic = s[which.min(aic), ],
    bic = s[which.min(bic), ], sse = s[which.min(sse), ],
    dev.roy = s[which.min(dev.roy), ]
  )
  
  # combine deviances, aic, bic, sse etc
  # order: M6=NULL, M4=linear(xi), M2=FP1(xi), M3=FP1(xia), M1=FP1(xi, xia),
  #        M5=linear(xia)
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
  
  outs
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
#' by [find_best_fp1_step()]. The linear function makes the formula more complicated, 
#' so there are separate functions for FP1 and FPm.
calculate_metrics_fpm <- function(y, 
                                  x, 
                                  xi, 
                                  powers_current, 
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
    out[[k]] <- find_best_fpm_step(
      y = y, x = x, xi = xi, powers_current = powers_current,
      powers = powers, family = family, weights = weights,
      offset = offset, strata = strata, control = control,
      method = method, rownames = rownames,
      nocenter = nocenter, degree = mm[k], acdx = acdx
    )
  }
  # save AIC, BIC, DEV and sse for m = 2,3,...
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

#' Function to extract and transform adjustment variables
#' 
#' @param x a matrix of predictors that includes the variable of interest `xi`. 
#' It is assumed that continuous variables have already been shifted and scaled.
#' @param xi name of the continuous predictor for which the FP function will be
#' estimated. There are no binary or two-level variables allowed. All variables
#' except `xi` are referred to as "adjustment variables".
#' @param powers_current a named list of FP powers of all variables of interest, 
#' including `xi`. Note that these powers are updated during backfitting or MFP 
#' cycles.
#' @param df a numeric vector of degrees of freedom for `xi`.
#' @param powers a set of allowed FP powers.
#' 
#' @details
#' After extracting the adjustment variables this function, using their
#' corresponding FP powers stored in `powers_current`, transforms them. 
#' This is necessary When evaluating x of interest, as we must account for other
#' variables, which can be transformed or untransformed, depending on the
#' individual powers. It's worth noting that some powers can be NA, indicating
#' that the variable has been left out of the adjustment variables. It also
#' returns the FP data, which is dependent on the degrees of freedom. For example,
#' `df = 2` is equivalent to FP degree one, resulting in the generation of 8 
#' variables. If `acdx` is set to `TRUE`, however, 64 variables are generated.
#' 
#' @return 
#' A list containing the following elements: 
#' 
#' * `powers_fp`: fp powers used for `data_fp`.
#' * `data_fp`: all possible fp transformations for `xi`, see the `data` 
#' component of the output of [generate_transformations_fp()] and 
#' [generate_transformations_acd()]. 
#' * `powers_adj`: fp powers for adjustment variables in `data_adj`.
#' * `data_adj`: adjustment data, i.e. transformed input data for adjustment
#' variables.
transform_data_step <- function(x,
                                xi,
                                powers_current,
                                df, 
                                powers,
                                acdx) {
  
  # sort x based on the names of the powers
  names_powers_current <- names(powers_current)
  x <- x[, names_powers_current, drop = FALSE]
  
  # extract the matrix of adjustment variables
  vars_adj <- names_powers_current[!(names_powers_current %in% xi)]
  x_adj <- x[, vars_adj, drop = FALSE]
  
  # transformations of adjustment variables
  powers_adj <- powers_current[vars_adj]
  acdx_adj <- unname(acdx[vars_adj])
  
  # generate adjustment data 
  # check whether all adjustment powers = NA
  if (all(is.na(unlist(powers_adj, use.names = FALSE)))) {
    # all adjustment variables were eliminated in MFP backfitting process
    data_adj <- NULL
    powers_adj <- NULL
  } else {
    data_adj <- vector(mode = "list", length = ncol(x_adj))
    
    for (i in 1:ncol(x_adj)) {
      if (acdx_adj[i]) {
        data_adj[[i]] <- transform_vector_acd(
          x = x_adj[, i, drop = TRUE], power = powers_adj[[i]],
          powers = powers
        )
      } else {
        data_adj[[i]] <- transform_vector_fp(
          x = x_adj[, i, drop = TRUE], power = powers_adj[[i]]
        )
      }
    }
    
    # combine into data.frame
    # note some variables may have been extended to more than one column
    # in the loop above due to transformation
    data_adj <- do.call(cbind, data_adj)
    
    # assign arbitrary names to adjustment matrix 
    # the names of adjustment variables at this stage are not relevant
    colnames(data_adj) <- paste0("var_adj", 1:ncol(data_adj))
  }
  
  # generate fp data
  data_xi <- x[, xi, drop = TRUE]
  
  if (length(unique(data_xi)) <= 3 || df == 1) {
    # if a variable has less than 4 levels we do not generate FP data
    # when df = 1 we do not generate FP data because we assume linearity
    data_fp <- data_xi
    powers_fp <- 1
  } else {
    if (acdx[xi]) {
      fpd <- generate_transformations_acd(data_xi, powers = powers)
    } else {
      # note that degree is df / 2
      fpd <- generate_transformations_fp(data_xi, degree = df / 2, powers = powers)
    }
    data_fp <- fpd$data
    powers_fp <- fpd$powers
  }
  
  list(
    powers_fp = powers_fp, 
    data_fp = data_fp, 
    powers_adj = powers_adj,
    data_adj = data_adj
  )
}
