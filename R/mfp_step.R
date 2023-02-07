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
#' [find_best_linear_step()]). This step differs from the mfp case because
#' linear models only use 1 df, while estimation of (every) fp power adds 
#' another df. 
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
#' @section Functional form selection:
#' There are 3 criteria to decide for the current best functional form of a 
#' continuous variable. These work as follows. 
#' 
#' In case `criterion = "pvalue"` the function selection procedure as outlined 
#' in Chapters 4 and 6 of Royston and Sauerbrei (2008) is used. Briefly, the
#' first step is to test the best FPm function against a null model at level
#' `select` with 2m df. If not significant, the variable is excluded. 
#' Otherwise, the fps are tested in order against the highest fp. 
#' That is, the next step is to test the best FPm versus linear at level `alpha` 
#' with 2m - 1 df. If not significant, use a linear model. 
#' Otherwise the next step is to test the best FPm versus the best FP1 at 
#' level `alpha` with 2m - 2 df. If not significant, use the best FP1 model. 
#' And so on, until FPm-1, which is tested at level `alpha` with 2 df.
#' If the final test is not significant, use a FPm-1 model, otherwise use FPm. 
#' 
#' Note that the "best" FPx model used in each step is given by the model using
#' a FPx transformation for the variable of interest and having the highest 
#' likelihood of all such models given the current powers for all other
#' variables, as outlined in Section 4.8 of Royston and Sauerbrei (2008).
#' These best FPx models are computed in [find_best_fpm_step()].
#' 
#' For the other criteria `aic` and `bic` all FP models up to the desired degree
#' are fitted and the model with the lowest value for the information criteria 
#' is chosen as the final one. 
#'  
#' @return 
#' A numeric vector indicating the best powers for `xi`. Entries can be 
#' `NA` if variable is to be removed from the working model. 
#' 
#' @references 
#' Royston, P. and Sauerbrei, W., 2008. \emph{Multivariable Model - Building: 
#' A Pragmatic Approach to Regression Anaylsis based on Fractional Polynomials 
#' for Modelling Continuous Variables. John Wiley & Sons.}\cr
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
    
    # model selection based on AIC, BIC or P-values, or keep xi if indicated
    if (xi %in% keep) {
      power_best = 1
      print(sprintf("Variable %s kept.", xi))
    } else {
      # if df = 1 then we just fit usual linear models and test: NULL vs Linear
      fit <- find_best_linear_step(
        x = x, y = y, xi = xi, family = family,
        powers_current = powers_current, powers = powers, acdx = acdx,
        select = select, ftest = ftest, 
        weights = weights, offset = offset, strata = strata, 
        method = method, control = control, rownames = rownames, 
        nocenter = nocenter
      )
      
      power_best = fit$power_best[[tolower(criterion)]]
      
      # TODO: simplify this part
      if (verbose) {
        if (criterion == "pvalue") {
          if (ftest) {
            print_mfp_summary_2(
              namex = xi, 
              dev.all = fit$metrics[, "deviance_stata"], 
              df.res = fit$metrics[, "df_resid"], 
              df.den = fit$metrics[, "df_resid"], 
              dev.diff = diff(fit$metrics[, "deviance_stata"]),
              f = ifelse(ftest, fit$statistic, NA),
              pvalues = fit$pvalue,
              best.function = list(1), 
              index.bestmodel = fit$model_best[criterion], 
              acd = acdx[xi]
            ) 
          } else {
            print_mfp_summary_1(
              namex = xi, 
              dev.all = fit$metrics[, "deviance_rs"], 
              dev.diff = diff(fit$metrics[, "deviance_rs"]),
              pvalues = fit$pvalue,
              index.bestmodel = fit$model_best[criterion],
              best.function = list(1), 
              acd = acdx[xi]
            ) 
          }
        } else {
          print_mfp_summary_3(
            xi, 
            gic = fit$metrics[, criterion], 
            keep = keep, 
            best.function = list(1),
            acd = FALSE
          )
        }
      }
    }
  } else if (acdx[xi]) {
    # acd case -----------------------------------------------------------------
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
    # usual mfp case -----------------------------------------------------------
    # in this case, we seek the best fp for a variable which may use more 
    # than only linear effects and does not use an acd transformation
    
    degree <- as.numeric(df / 2)
    
    fit <- select_ra2(
      y = y, x = x, xi = xi, degree = degree, keep = keep, 
      powers_current = powers_current, acdx = acdx, powers = powers, 
      select = select, alpha = alpha, ftest = ftest, 
      family = family, method = method, weights = weights, offset = offset,
      strata = strata, control = control, rownames = rownames,
      nocenter = nocenter 
    ) 
    
    # TODO AIC BIC
    
    # TODO: printing
     
    # if (verbose) {
    #   if (criterion == "pvalue") {
    #     best.function1 <- if (degree == 1) {
    #       list(lin = 1, fpm = fit$power_best)
    #     } else {
    #       append(if (ftest) {
    #         NA
    #       } else {
    #         NA
    #       }, list(1), 0)
    #     }
    #     
    #     if (ftest) {
    #       print_mfp_summary_2(
    #         namex = xi, 
    #         dev.all = fit$metrics[, "deviance_rs"], 
    #         df.res = fit$metrics[, "df_resid"],
    #         dev.diff = NA, 
    #         f = fit$statistic,
    #         df.den = fit$metrics[, "df_resid"],
    #         pvalues = fit$pvalue, 
    #         best.function = best.function1,
    #         index.bestmodel = fit$model_best, 
    #         acd = F
    #       )
    #     } else {
    #       print_mfp_summary_1(
    #         namex = xi, 
    #         dev.all = fit$metrics[, "deviance_rs"], 
    #         dev.diff = NA,
    #         pvalues = fit$pvalue, 
    #         index.bestmodel = fit$model_best,
    #         best.function = best.function1, 
    #         acd = F
    #       )
    #     }
    #     # AIC and BIC display
    #   } else {
    #     best.function.aic <- if (degree == 1) {
    #       list(lin = 1, fpm = bestfp1x[2])
    #     } else {
    #       append(bestfp.aic, list(1), 0)
    #     }
    #     best.function.bic <- if (degree == 1) {
    #       list(lin = 1, fpm = bestfp1x[3])
    #     } else {
    #       append(bestfp.bic, list(1), 0)
    #     }
    #     switch(criterion,
    #            "AIC" = print_mfp_summary_3(xi, gic = aic.all, keep = keep, best.function = best.function.aic, acd = F),
    #            "BIC" = print_mfp_summary_3(xi, gic = bic.all, keep = keep, best.function = best.function.bic, acd = F)
    #     )
    #   }
    # }
    
    power_best <- as.numeric(fit$power_best)
  }
  
  power_best
}


#' Function to find the best FP functions of given degree for a single variable
#' 
#' Handles the FP1 and the higher order FP cases. For parameter definitions, see
#' [find_best_fp_step()].
#' 
#' @details 
#' The "best" model is determined by the highest likelihood (or smallest 
#' deviance by our definition as minus twice the log-likelihood). This is also 
#' the case for the use of information criteria, as all models investigated in 
#' this function have the same df, so the penalization term is equal for all
#' models and only their likelihoods differ.
#' 
#' Note that the estimation of each fp power adds a degree of freedom. Thus, 
#' all fp1s have 2 df, all fp2s have 4 df and so on.
#' 
#' In the case that `degree = 1`, the linear model (fp power of 1) is NOT 
#' returned, as it is not considered to be a fractional polynomial in this
#' algorithm (as a linear model has only one df, whereas the same function 
#' regarded as fp would have 2 fp).
#' 
#' @return 
#' A list with several components giving the best power found (`power_best`) and 
#' performance indices.
find_best_fpm_step <- function(x, 
                               xi,
                               degree,
                               y, 
                               powers_current,
                               powers,  
                               acdx, 
                               ...) {
  
  n_obs <- dim(x)[1L]
  
  if (degree == 1) {
    # remove linear model
    powers = setdiff(powers, c(1))
  }
    
  
  # generate FP data for x of interest (xi) and adjustment variables
  x_transformed <- transform_data_step(
    x = x, xi = xi, df = 2 * degree,
    powers_current = powers_current, powers = powers, acdx = acdx
  )
  
  metrics = list()
  
  for (i in seq_along(x_transformed$data_fp)) {
    # combine FP variables for x of interest with adjustment variables
    fit <- fit_model(
      x = cbind(x_transformed$data_fp[[i]], x_transformed$data_adj), y = y, ...
    )
    
    # use degree many additional degrees of freedom
    # TODO: here we should likely add additional df for each fp term in the model
    # using calculate_number_fp_powers
    # note: this doesn't change WHICH model is the best, since all use the 
    # same additional df, but it may affect further comparison between 
    # different fp models
    p = x_transformed$powers_fp[i, , drop = TRUE]
    metrics[[paste(p, collapse = " ")]] = calculate_model_metrics(
      fit, n_obs, degree
    )
  }
  
  metrics = do.call(rbind, metrics) 
  model_best = as.numeric(which.max(metrics[, "logl"]))
   
  list(
    powers = x_transformed$powers_fp, 
    power_best = x_transformed$powers_fp[model_best, , drop = TRUE], 
    metrics = metrics, 
    model_best = model_best
  ) 
}

fit_null_linear_step <- function(x, 
                                 xi, 
                                 y, 
                                 powers_current,
                                 powers,
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
  
  list(
    powers = NA,
    metrics = rbind(null = calculate_model_metrics(model_null, n_obs))  
  )
}

fit_linear_step <- function(x, 
                            xi, 
                            y, 
                            powers_current,
                            powers,
                            acdx, 
                            ...) {
  
  n_obs <- dim(x)[1L]
  
  # transform all data as given by current working model
  # set variable of interest to linear term only
  x_transformed <- transform_data_step(
    x = x, xi = xi, df = 1,
    powers_current = powers_current, acdx = acdx, powers = powers
  ) 
  
  # fit a model based on the assumption that xi is linear 
  model_linear <- fit_model(
    x = cbind(x_transformed$data_fp[[1]], x_transformed$data_adj), y = y, ...
  )
  
  # respect acd
  metrics <- rbind(linear = calculate_model_metrics(model_linear, n_obs))  
  if (acdx[xi])
    rownames(metrics) <- "A(linear)"
  
  list(
    powers = x_transformed$powers_fp,
    metrics = metrics
  )
}

#' Helper to assess null and linear term for a single variable
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
#' A list with several components giving the best power found (`power_best`) as
#' a numeric vector according to different criteria ("aic", "bic" or "pvalue"), 
#' and the associated performance indices. 
#' The returned best power may be `NA`, indicating the variable has been 
#' removed from the model.
find_best_linear_step <- function(x, 
                                  xi, 
                                  y, 
                                  powers_current,
                                  powers, 
                                  ftest,  
                                  acdx, 
                                  select,
                                  ...) {
  
  n_obs <- dim(x)[1L]

  fit <- fit_null_linear_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx, 
    ...
  )
  
  # compute best model according to different criteria
  if (ftest) {
    stats <- calculate_f_test(
      deviances = fit$metrics[, "deviance_stata"], 
      dfs_resid = fit$metrics[, "df_resid"],
      n_obs = n_obs
    )
  } else {
    stats <- calculate_lr_test(fit$metrics[, "logl"], fit$metrics[, "df"])
  }
  pvalue <- stats$pvalue
  names(pvalue) <- c("null vs linear")
  statistic <- stats$statistic 
  names(statistic) <- c("null vs linear")
  
  model_best <- c(
    ifelse(pvalue > select, 1, 2),
    which.min(fit$metrics[, "aic", drop = TRUE]), 
    which.min(fit$metrics[, "bic", drop = TRUE])
  )
  # make sure names are correct and do not carry over
  names(model_best) = c("pvalue", "aic", "bic")
  
  list(
    powers = fit$powers,
    power_best = ifelse(model_best == 1, NA, 1),
    metrics = fit$metrics,
    model_best = model_best,
    pvalue = pvalue,
    statistic = statistic
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

#' Function selection procedure 
select_ra2 <- function(x, 
                       xi,
                       degree,
                       y, 
                       powers_current, 
                       select, 
                       alpha, 
                       keep, 
                       powers, 
                       ftest,  
                       acdx, 
                       ...) {
  
  if (degree < 1)
    return(NULL)
  
  # simplify testing by defining test helper function
  if (ftest) {
    calculate_test <- function(metrics, n_obs) {
      calculate_f_test(
        deviances = metrics[, "deviance_rs", drop = TRUE],
        dfs_resid = metrics[, "df_resid", drop = TRUE],
        n_obs = n_obs
      )
    }
  } else {
    calculate_test <- function(metrics, ...) {
      calculate_lr_test(
        logl = metrics[, "logl", drop = TRUE], 
        dfs = metrics[, "df", drop = TRUE] 
      )
    }
  }

  n_obs = nrow(x)
  fpmax = paste0("FP", degree)
  
  # output list
  res <- list(
    power_best = NULL, 
    metrics = NULL, 
    model_best = NULL, 
    statistic = NULL, 
    pvalue = NULL
  )
  
  # fit highest fp and null / linear model for initial steps
  fit_fpmax <- find_best_fpm_step(
    x = x, xi = xi, degree = degree, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx, ...
  )
  fit_lin <- fit_null_linear_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx, ...
  )
  res$metrics <- rbind(
    fit_fpmax$metrics[fit_fpmax$model_best, ],
    fit_lin$metrics
  )
  rownames(res$metrics) <- c(fpmax, "null", "linear")
  
  # selection procedure
  
  # test for overall significance
  # df for tests are degree * 2
  stats <- calculate_test(res$metrics[c("null", fpmax), ], n_obs)
  res$statistic <- stats$statistic
  names(res$statistic) <- sprintf("%s vs null", fpmax)
  res$pvalue <- stats$pvalue
  names(res$pvalue) <- names(res$statistic)
  
  if (stats$pvalue >= select && !(xi %in% keep)) {
    # not selected and not forced into model
    res$power_best = NA
    res$model_best = 2
    return(res)
  }
  
  # test for non-linearity
  # df for tests are degree * 2 - 1
  stats <- calculate_test(res$metrics[c("linear", fpmax), ], n_obs)
  old_names <- names(res$statistic)
  res$statistic <- c(res$statistic, stats$statistic)
  names(res$statistic) <- c(old_names, sprintf("%s vs linear", fpmax))
  res$pvalue <- c(res$pvalue, stats$pvalue)
  names(res$pvalue) <- names(res$statistic)
  
  if (stats$pvalue >= alpha) {
    # no non-linearity detected
    res$power_best = 1
    res$model_best = 3
    return(res)
  }
  
  # tests for functional form - do this for all fps with lower degrees
  # dfs for tests are decreasing
  if (degree > 1) {
    
    for (current_degree in 1:(degree - 1)) {
      fpm = paste0("FP", current_degree)
      
      fit_fpm <- find_best_fpm_step(
        x = x, xi = xi, degree = current_degree, y = y, 
        powers_current = powers_current, powers = powers, acdx = acdx, ...
      )
      
      old_names = rownames(res$metrics)
      res$metrics <- rbind(
        res$metrics, 
        fit_fpm$metrics[fit_fpm$model_best, ]
      )
      rownames(res$metrics) <- c(old_names, fpm)
      
      stats <- calculate_test(res$metrics[c(fpm, fpmax), ], n_obs)
      old_names <- names(res$statistic)
      res$statistic <- c(res$statistic, stats$statistic)
      names(res$statistic) <- c(old_names, sprintf("%s vs %s", fpmax, fpm))
      res$pvalue <- c(res$pvalue, stats$pvalue)
      names(res$pvalue) <- names(res$statistic)
      
      if (stats$pvalue >= alpha) {
        # non-linearity detected, but lower than maximum degree
        res$power_best = fit_fpm$powers[fit_fpm$model_best, , drop = FALSE]
        res$model_best = nrow(res$metrics)
        return(res)
      }
    }
    
  }
  
  # return highest power
  res$power_best = fit_fpmax$powers[fit_fpmax$model_best, , drop = FALSE]
  res$model_best = 1
  
  res
}

#' Function selection procedure 
select_ra2_acd <- function(x, 
                           xi,
                           degree,
                           y, 
                           powers_current, 
                           select, 
                           alpha, 
                           keep, 
                           powers, 
                           ftest,  
                           acdx, 
                           ...) {
  
  
  # TODO: simplify, using local testing function to prevent if else all the
  # time
  
  if (degree < 1)
    return(NULL)
  
  n_obs = nrow(x)
  fpmax = paste0("FP", degree)
  
  res <- list(
    power_best = NULL, 
    metrics = NULL, 
    model_best = NULL, 
    statistic = NULL, 
    pvalue = NULL
  )
  
  # fit highest fp and null / linear model for initial steps
  fit_fpmax <- find_best_fpm_step(
    x = x, xi = xi, degree = degree, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx, ...
  )
  fit_lin <- fit_null_linear_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx, ...
  )
  res$metrics <- rbind(
    fit_fpmax$metrics[fit_fpmax$model_best, ],
    fit_lin$metrics
  )
  rownames(res$metrics) <- c(fpmax, "null", "linear")
  
  # selection procedure
  
  # test for overall significance
  if (ftest) {
    stats <- calculate_f_test(
      deviances = res$metrics[c("null", fpmax), "deviance_rs", drop = TRUE],
      dfs_resid = res$metrics[c("null", fpmax), "df_resid", drop = TRUE],
      n_obs = n_obs
    )
  } else {
    stats <- calculate_lr_test(
      logl = res$metrics[c("null", fpmax), "logl", drop = TRUE], 
      dfs = res$metrics[c("null", fpmax), "df", drop = TRUE] 
    )  
  }
  res$statistic <- stats$statistic
  names(res$statistic) <- sprintf("%s vs null", fpmax)
  res$pvalue <- stats$pvalue
  names(res$pvalue) <- names(res$statistic)
  
  if (stats$pvalue >= select && !(xi %in% keep)) {
    # not selected and not forced into model
    res$power_best = NA
    res$model_best = 2
    return(res)
  }
  
  # test for non-linearity
  if (ftest) {
    stats <- calculate_f_test(
      deviances = res$metrics[c("linear", fpmax), "deviance_rs", drop = TRUE],
      dfs_resid = res$metrics[c("linear", fpmax), "df_resid", drop = TRUE],
      n_obs = n_obs
    )
  } else {
    stats <- calculate_lr_test(
      logl = res$metrics[c("linear", fpmax), "logl", drop = TRUE], 
      dfs = res$metrics[c("linear", fpmax), "df", drop = TRUE] 
    )  
  }
  old_names <- names(res$statistic)
  res$statistic <- c(res$statistic, stats$statistic)
  names(res$statistic) <- c(old_names, sprintf("%s vs linear", fpmax))
  res$pvalue <- c(res$pvalue, stats$pvalue)
  names(res$pvalue) <- names(res$statistic)
  
  if (stats$pvalue >= alpha) {
    # no non-linearity detected
    res$power_best = 1
    res$model_best = 3
    return(res)
  }
  
  # tests for functional form - do this for all fps with lower degrees
  if (degree > 1) {
    
    for (current_degree in 1:(degree - 1)) {
      fpm = paste0("FP", current_degree)
      
      fit_fpm <- find_best_fpm_step(
        x = x, xi = xi, degree = current_degree, y = y, 
        powers_current = powers_current, powers = powers, acdx = acdx, ...
      )
      
      old_names = rownames(res$metrics)
      res$metrics <- rbind(
        res$metrics, 
        fit_fpm$metrics[fit_fpm$model_best, ]
      )
      rownames(res$metrics) <- c(old_names, fpm)
      
      if (ftest) {
        stats <- calculate_f_test(
          deviances = res$metrics[c(fpm, fpmax), "deviance_rs", drop = TRUE],
          dfs_resid = res$metrics[c(fpm, fpmax), "df_resid", drop = TRUE],
          n_obs = n_obs
        )
      } else {
        stats <- calculate_lr_test(
          logl = res$metrics[c(fpm, fpmax), "logl", drop = TRUE], 
          dfs = res$metrics[c(fpm, fpmax), "df", drop = TRUE] 
        )  
      }
      old_names <- names(res$statistic)
      res$statistic <- c(res$statistic, stats$statistic)
      names(res$statistic) <- c(old_names, sprintf("%s vs %s", fpmax, fpm))
      res$pvalue <- c(res$pvalue, stats$pvalue)
      names(res$pvalue) <- names(res$statistic)
      
      if (stats$pvalue >= alpha) {
        # non-linearity detected, but lower than maximum degree
        res$power_best = fit_fpm$powers[fit_fpm$model_best, , drop = FALSE]
        res$model_best = nrow(res$metrics)
        return(res)
      }
    }
    
  }
  
  # return highest power
  res$power_best = fit_fpmax$powers[fit_fpmax$model_best, , drop = FALSE]
  res$model_best = 1
  
  res
}

select_ic <- function() {
  
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
#' @param acdx a logical vector indicating the use of acd transformation.
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
#' variables. If `acdx` for the current variables of interest is set to `TRUE`, 
#' however, 64 variables are generated.
#' 
#' When `df = 1`, this function returns data unchanged, i.e. a "linear" 
#' transformation with power equal to 1. In case `acdx[xi] = TRUE`, the 
#' acd transformation is applied. 
#' 
#' @return 
#' A list containing the following elements: 
#' 
#' * `powers_fp`: fp powers used for `data_fp`.
#' * `data_fp`: a list with all possible fp transformations for `xi`, see the 
#' `data` component of the output of [generate_transformations_fp()] and 
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
    data_adj <- matrix(nrow = 0, ncol = 0)
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
  
  if (length(unique(data_xi)) <= 3) {
    # if a variable has less than 4 levels we do not generate FP data
    data_fp <- list(data_xi)
    powers_fp <- 1
  } else {
    if (acdx[xi]) {
      # when df = 1 -> degree = 0
      # and we return the data unchanged, i.e. with power = 1
      fpd <- generate_transformations_acd(data_xi, degree = floor(df / 2),
                                          powers = powers)
    } else {
      # note that degree is df / 2
      fpd <- generate_transformations_fp(data_xi, degree = floor(df / 2),
                                         powers = powers)
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
