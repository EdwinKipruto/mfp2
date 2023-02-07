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
#' [select_linear()]). This step differs from the mfp case because
#' linear models only use 1 df, while estimation of (every) fp power adds 
#' another df. 
#' * the case that an acd transformation is requested (`acdx` is `TRUE` 
#' for `xi`) for the variable of interest (see [find_best_fpm_step()]).
#' * the (usual) case of the normal mfp algorithm to assess non-linear 
#' functional forms (see [find_best_fpm_step()]). 
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
#' continuous variable. 
#' 
#' The first option for `criterion = "pvalue"` is the function selection 
#' procedure as outlined in e.g. Chapters 4 and 6 of Royston and 
#' Sauerbrei (2008), also abbreviated as "RA2".
#' It is a closed testing procedure and is implemented in [select_ra2()] and
#' extended for ACD transformation in [select_ra2_acd()] according to 
#' Royston and Sauerbrei (2016). 
#' 
#' For the other criteria `aic` and `bic` all FP models up to the desired degree
#' are fitted and the model with the lowest value for the information criteria 
#' is chosen as the final one. This is implemented in [select_ic()].
#'  
#' @return 
#' A numeric vector indicating the best powers for `xi`. Entries can be 
#' `NA` if variable is to be removed from the working model. 
#' 
#' @references 
#' Royston, P. and Sauerbrei, W., 2008. \emph{Multivariable Model - Building: 
#' A Pragmatic Approach to Regression Anaylsis based on Fractional Polynomials 
#' for Modelling Continuous Variables. John Wiley & Sons.}\cr
#' Royston, P. and Sauerbrei, W., 2016. \emph{mfpa: Extension of mfp using the
#' ACD covariate transformation for enhanced parametric multivariable modeling. 
#' The Stata Journal, 16(1), pp.72-87.}
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
      fit <- select_linear(
        x = x, y = y, xi = xi, family = family,
        powers_current = powers_current, powers = powers, acdx = acdx,
        select = select, ftest = ftest, 
        weights = weights, offset = offset, strata = strata, 
        method = method, control = control, rownames = rownames, 
        nocenter = nocenter, criterion = criterion
      )
      
      power_best = fit$power_best
      
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

    fit <- select_ra2_acd(
      y = y, x = x, xi = xi, keep = keep, 
      powers_current = powers_current, acdx = acdx, powers = powers, 
      select = select, alpha = alpha, ftest = ftest, 
      family = family, method = method, weights = weights, offset = offset,
      strata = strata, control = control, rownames = rownames,
      nocenter = nocenter 
    ) 
    
    # TODO: AIC / BIC
    
    # TODO: Printing on the screen
    # if (verbose) {
    #   if (criterion == "pvalue") {
    #     if (ftest) {
    #       print_mfp_summary_2(
    #         namex = xi, dev.all = dev.roy.all, df.res = df.all, dev.diff = dev.diff, f = fstatistic,
    #         df.den = df.all, pvalues = pvalue, best.function = bestfuns[[5]],
    #         index.bestmodel = index.bestmodel, acd = T
    #       )
    #     } else {
    #       print_mfp_summary_1(
    #         namex = xi,
    #         dev.all = dev.all,
    #         dev.diff = dev.diff,
    #         pvalues = pvalue,
    #         index.bestmodel = index.bestmodel,
    #         best.function = bestfuns[[1]], acd = T
    #       )
    #     }
    #   } else {
    #     switch(criterion,
    #            "AIC" = print_mfp_summary_3(xi, gic = aic.all, keep = keep, best.function = bestfuns[[2]], acd = T),
    #            "BIC" = print_mfp_summary_3(xi, gic = bic.all, keep = keep, best.function = bestfuns[[3]], acd = T)
    #     )
    #   }
    # }
    
    power_best <- as.numeric(fit$power_best)
    
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
#' algorithm. 
#' A linear model has only one df, whereas the same function regarded as fp 
#' would have 2 fp.
#' 
#' @section ACD transformation:
#' This function also handles the case of ACD transformations if `acdx` is set
#' to `TRUE` for `xi`. In this case, if `degree = 1`, then 7 models are
#' assessed (like for the non-acd case it excludes the linear case), 
#' and if `degree = 2`, then 64 models are assessed (unlike the 36 models 
#' for non-acd transformation). Other settings for `degree` are currently not
#' supported when used with ACD transformations.
#' 
#' @return 
#' A list with several components: 
#' 
#' * `acd`: logical indicating if an ACD transformation was applied for `xi`.
#' * `powers`: fp powers investigated in step. 
#' * `power_best`: the best power found. `power_best` will always be a 
#' two-column matrix when an ACD transformation is used, otherwise the number 
#' of columns will depend on `degree`. 
#' * `metrics`: a matrix with performance indices for all models investigated. 
#' Same number of rows as, and indexed by, `powers`.
#' * `model_best`: row index of best model in `metrics`.
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
    # remove linear model for normal data, but keep for acd transformation
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
    p = sprintf("%g", x_transformed$powers_fp[i, , drop = TRUE])
    # respect acd
    if (acdx[xi])
      p[length(p)] = sprintf("A(%s)", p[length(p)])
    
    metrics[[paste(p, collapse = " ")]] = calculate_model_metrics(
      fit, n_obs, degree
    )
  }
  
  metrics = do.call(rbind, metrics) 
  model_best = as.numeric(which.max(metrics[, "logl"]))
   
  list(
    acd = acdx[xi],
    powers = x_transformed$powers_fp, 
    power_best = x_transformed$powers_fp[model_best, , drop = TRUE], 
    metrics = metrics, 
    model_best = model_best
  ) 
}

#' Function to fit null model excluding variable of interest
#' 
#' "Null" model here refers to a model which does not include the variable 
#' of interest `xi`. 
#' For parameter definitions, see [find_best_fp_step()]. All parameters 
#' captured by `...` are passed on to [fit_model()].
#' 
#' @return 
#' A list with two entries: 
#' 
#' * `powers`: fp power(s) of `xi` in fitted model - in this case `NA`.
#' * `metrics`: a matrix with performance indices for fitted model.
fit_null_step <- function(x, 
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

#' Function to fit linear model for variable of interest
#' 
#' "Linear" model here refers to a model which includes the variable 
#' of interest `xi` with a fp power of 1. Note that `xi` may be ACD transformed
#' if indicated by `acdx[xi]`.
#' For parameter definitions, see [find_best_fp_step()]. All parameters 
#' captured by `...` are passed on to [fit_model()].
#' 
#' @return 
#' A list with two entries: 
#' 
#' * `powers`: fp power(s) of `xi` (or its ACD transformation) in fitted model.
#' * `metrics`: a matrix with performance indices for fitted model.
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
    rownames(metrics) <- "linear(., A(x))"
  
  list(
    powers = x_transformed$powers_fp,
    metrics = metrics
  )
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
#' Note that this function handles an ACD transformation for `xi` as well. 
#' 
#' @return 
#' A list with several components:
#' 
#' * `acd`: logical indicating if an ACD transformation was applied for `xi`.
#' * `powers`: fp powers investigated in step, indexing `metrics`. 
#' * `power_best`: a numeric vector with the best power found. The returned 
#' best power may be `NA`, indicating the variable has been removed from the 
#' model.
#' * `metrics`: a matrix with performance indices for all models investigated. 
#' Same number of rows as, and indexed by, `powers`.
#' * `model_best`: row index of best model in `metrics`.
#' * `pvalue`: p-value for comparison of linear and null model.
#' * `statistic`: test statistic used, depends on `ftest`.
select_linear <- function(x, 
                          xi, 
                          y, 
                          powers_current,
                          powers, 
                          criterion,
                          ftest,  
                          acdx, 
                          select,
                          ...) {
  
  n_obs <- dim(x)[1L]
  
  fit_null <- fit_null_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx, 
    ...
  )
  fit_linear <- fit_linear_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx, 
    ...
  )
  powers = rbind(fit_null$powers, fit_linear$powers)
  metrics = rbind(fit_null$metrics, fit_linear$metrics)
  
  # compute best model according to different criteria
  if (ftest) {
    stats <- calculate_f_test(
      deviances = metrics[, "deviance_stata"], 
      dfs_resid = metrics[, "df_resid"],
      n_obs = n_obs
    )
  } else {
    stats <- calculate_lr_test(metrics[, "logl"], metrics[, "df"])
  }
  pvalue <- stats$pvalue
  names(pvalue) <- c("null vs linear")
  statistic <- stats$statistic 
  names(statistic) <- c("null vs linear")
  
  model_best <- switch(
    tolower(criterion), 
    "pvalue" = ifelse(pvalue > select, 1, 2), 
    "aic" = which.min(metrics[, "aic", drop = TRUE]), 
    "bic" = which.min(metrics[, "bic", drop = TRUE])
  ) 

  list(
    acd = acdx[xi],
    powers = powers,
    power_best = powers[model_best, ],
    metrics = metrics,
    model_best = model_best,
    pvalue = pvalue,
    statistic = statistic
  )
}

#' Function selection procedure based on closed testing procedure
#' 
#' Used in [find_best_fp_step()] when `criterion = "pvalue"`.
#' For parameter explanations, see [find_best_fp_step()]. All parameters 
#' captured by `...` are passed on to [fit_model()].
#' 
#' @details  
#' In case `criterion = "pvalue"` the function selection procedure as outlined 
#' in Chapters 4 and 6 of Royston and Sauerbrei (2008) is used. 
#' 
#' * \emph{Step 1}: test the best FPm function against a null model at level
#' `select` with 2m df. If not significant, the variable is excluded. 
#' Otherwise continue with step 2.
#' * \emph{Step 2}: test the best FPm versus a linear model at level `alpha` 
#' with 2m - 1 df. If not significant, use a linear model. 
#' Otherwise continue with step 3.
#' * \emph{Step 3}: test the best FPm versus the best FP1 at 
#' level `alpha` with 2m - 2 df. If not significant, use the best FP1 model. 
#' Otherwise, repeat this step for all remaining higher order FPs until 
#' FPm-1, which is tested at level `alpha` with 2 df against FPm.
#' If the final test is not significant, use a FPm-1 model, otherwise use FPm. 
#' 
#' Note that the "best" FPx model used in each step is given by the model using
#' a FPx transformation for the variable of interest and having the highest 
#' likelihood of all such models given the current powers for all other
#' variables, as outlined in Section 4.8 of Royston and Sauerbrei (2008).
#' These best FPx models are computed in [find_best_fpm_step()].
#' 
#' @return 
#' A list with several components:
#' 
#' * `acd`: logical indicating if an ACD transformation was applied for `xi`, 
#' i.e. `FALSE` in this case.
#' * `powers`: (best) fp powers investigated in step, indexing `metrics`. 
#' * `power_best`: a numeric vector with the best power found. The returned 
#' best power may be `NA`, indicating the variable has been removed from the 
#' model.
#' * `metrics`: a matrix with performance indices for all models investigated. 
#' Same number of rows as, and indexed by, `powers`.
#' * `model_best`: row index of best model in `metrics`.
#' * `pvalue`: p-value for comparison of linear and null model.
#' * `statistic`: test statistic used, depends on `ftest`.
#' 
#' @references 
#' Royston, P. and Sauerbrei, W., 2008. \emph{Multivariable Model - Building: 
#' A Pragmatic Approach to Regression Anaylsis based on Fractional Polynomials 
#' for Modelling Continuous Variables. John Wiley & Sons.}
#' 
#' @seealso 
#' [select_ra2_acd()]
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
    acd = FALSE, 
    powers = NULL, 
    power_best = NULL, 
    metrics = NULL, 
    model_best = NULL, 
    statistic = NULL, 
    pvalue = NULL
  )
  
  # fit highest fp and null model for initial step
  fit_fpmax <- find_best_fpm_step(
    x = x, xi = xi, degree = degree, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx, ...
  )
  fit_null <- fit_null_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx, ...
  )
  res$metrics <- rbind(
    fit_fpmax$metrics[fit_fpmax$model_best, ],
    fit_null$metrics
  )
  rownames(res$metrics) <- c(fpmax, "null")
  res$powers <- rbind(fit_fpmax$power_best, NA)
  
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
  fit_lin <- fit_linear_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx, ...
  )
  
  old_names = rownames(res$metrics)
  res$metrics <- rbind(
    res$metrics, 
    fit_lin$metrics
  )
  rownames(res$metrics) <- c(old_names, "linear")
  res$powers = rbind(res$powers, 
                     ensure_length(fit_lin$powers, ncol(res$powers)))
  
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
      res$powers <- rbind(res$powers, 
                          ensure_length(fit_fpm$power_best, ncol(res$powers)))
      
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

#' Function selection procedure for ACD based on closed testing procedure
#' 
#' Used in [find_best_fp_step()] when `criterion = "pvalue"` and an 
#' ACD transformation is requested for `xi`.
#' For parameter explanations, see [find_best_fp_step()]. All parameters 
#' captured by `...` are passed on to [fit_model()].
#' 
#' @details  
#' This function extends the algorithm used in [select_ra2()] to allow the 
#' usage of ACD transformations. The implementation follows the description 
#' in Royston and Sauerbrei (2016). The procedure is outlined in detail in 
#' the corresponding section in the documentation of [mfpa()].
#' 
#' @return 
#' A list with several components:
#' 
#' * `acd`: logical indicating if an ACD transformation was applied for `xi`, 
#' i.e. `FALSE` in this case.
#' * `powers`: (best) fp powers investigated in step, indexing `metrics`. 
#' * `power_best`: a numeric vector with the best power found. The returned 
#' best power may be `NA`, indicating the variable has been removed from the 
#' model.
#' * `metrics`: a matrix with performance indices for all models investigated. 
#' Same number of rows as, and indexed by, `powers`.
#' * `model_best`: row index of best model in `metrics`.
#' * `pvalue`: p-value for comparison of linear and null model.
#' * `statistic`: test statistic used, depends on `ftest`.
#' 
#' @references 
#' Royston, P. and Sauerbrei, W., 2016. \emph{mfpa: Extension of mfp using the
#' ACD covariate transformation for enhanced parametric multivariable modeling. 
#' The Stata Journal, 16(1), pp.72-87.}
#' 
#' @seealso 
#' [select_ra2()]
select_ra2_acd <- function(x, 
                           xi,
                           y, 
                           powers_current, 
                           select, 
                           alpha, 
                           keep, 
                           powers, 
                           ftest,  
                           acdx, 
                           ...) {
  
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
  fpmax = "FP1(x, A(x))"
  acdx_reset_xi = acdx
  acdx_reset_xi[xi] = FALSE
  
  # output list
  res <- list(
    acd = TRUE,
    powers = NULL, 
    power_best = NULL, 
    metrics = NULL, 
    model_best = NULL, 
    statistic = NULL, 
    pvalue = NULL
  )
  
  # fit highest fp and null model for initial step
  fit_fpmax <- find_best_fpm_step(
    x = x, xi = xi, degree = 2, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx, ...
  )
  fit_null <- fit_null_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx, ...
  )
  res$metrics <- rbind(
    fit_fpmax$metrics[fit_fpmax$model_best, ],
    fit_null$metrics
  )
  rownames(res$metrics) <- c(fpmax, "null")
  res$powers <- rbind(fit_fpmax$power_best, fit_null$powers)
  
  # test for overall significance
  # df for tests are degree * 2 = 4
  stats <- calculate_test(res$metrics[c("null", fpmax), ], n_obs)
  res$statistic <- stats$statistic
  names(res$statistic) <- sprintf("%s vs null", fpmax)
  res$pvalue <- stats$pvalue
  names(res$pvalue) <- names(res$statistic)
  
  if (stats$pvalue >= select && !(xi %in% keep)) {
    # not selected and not forced into model
    res$power_best = matrix(c(NA, NA), ncol = 2)
    res$model_best = 2
    return(res)
  }
  
  # test for non-linearity in x
  # df for tests are degree * 2 - 1 = 3
  fit_lin <- fit_linear_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx_reset_xi, ...
  )
  
  old_names = rownames(res$metrics)
  res$metrics <- rbind(
    res$metrics, 
    fit_lin$metrics
  )
  rownames(res$metrics) <- c(old_names, "linear")
  res$powers <- rbind(res$powers, c(1, NA))
  
  stats <- calculate_test(res$metrics[c("linear", fpmax), ], n_obs)
  old_names <- names(res$statistic)
  res$statistic <- c(res$statistic, stats$statistic)
  names(res$statistic) <- c(old_names, sprintf("%s vs linear", fpmax))
  res$pvalue <- c(res$pvalue, stats$pvalue)
  names(res$pvalue) <- names(res$statistic)
  
  if (stats$pvalue >= alpha) {
    # no non-linearity detected
    res$power_best = matrix(c(1, NA), ncol = 2)
    res$model_best = 3
    return(res)
  }
  
  # test for functional form, comparison with FP1(x, .)
  fit <- find_best_fpm_step(
    x = x, xi = xi, degree = 1, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx_reset_xi, ...
  )
  
  old_names = rownames(res$metrics)
  res$metrics <- rbind(
    res$metrics, 
    fit$metrics[fit$model_best, ]
  )
  rownames(res$metrics) <- c(old_names, "FP1(x, .)")
  res$powers <- rbind(res$powers, c(fit$power_best, NA))
  
  stats <- calculate_test(res$metrics[c("FP1(x, .)", fpmax), ], n_obs)
  old_names <- names(res$statistic)
  res$statistic <- c(res$statistic, stats$statistic)
  names(res$statistic) <- c(old_names, sprintf("%s vs FP1(x, .)", fpmax))
  res$pvalue <- c(res$pvalue, stats$pvalue)
  names(res$pvalue) <- names(res$statistic)
  
  if (stats$pvalue >= alpha) {
    # FP1(x, .) is good enough
    res$power_best = matrix(c(fit$power_best, NA), ncol = 2)
    res$model_best = 4
    return(res)
  }
  
  # test for functional form, comparison with FP1(., A(x))
  fit_fp1a <- find_best_fpm_step(
    x = x, xi = xi, degree = 1, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx, ...
  )
  
  old_names = rownames(res$metrics)
  res$metrics <- rbind(
    res$metrics, 
    fit_fp1a$metrics[fit_fp1a$model_best, ]
  )
  rownames(res$metrics) <- c(old_names, "FP1(., A(x))")
  res$powers <- rbind(res$powers, fit_fp1a$power_best)
  
  stats <- calculate_test(res$metrics[c("FP1(., A(x))", fpmax), ], n_obs)
  old_names <- names(res$statistic)
  res$statistic <- c(res$statistic, stats$statistic)
  names(res$statistic) <- c(old_names, sprintf("%s vs FP1(., A(x))", fpmax))
  res$pvalue <- c(res$pvalue, stats$pvalue)
  names(res$pvalue) <- names(res$statistic)
  
  if (stats$pvalue < alpha) {
    # FP1(x, A(x)) is the best
    res$power_best = fit_fpmax$power_best
    res$model_best = 1
    return(res)
  }
  
  # return best model between FP1(., A(x)) and linear(., A(x))
  fit_lineara <- fit_linear_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx, ...
  )
  
  old_names = rownames(res$metrics)
  res$metrics <- rbind(
    res$metrics, 
    fit_lineara$metrics
  )
  rownames(res$metrics) <- c(old_names, "linear(., A(x))")
  res$powers <- rbind(res$powers, fit_lineara$powers)
  
  stats <- calculate_test(res$metrics[c("linear(., A(x))", "FP1(., A(x))"), ], n_obs)
  old_names <- names(res$statistic)
  res$statistic <- c(res$statistic, stats$statistic)
  names(res$statistic) <- c(old_names, sprintf("%s vs linear(., A(x))", "FP1(., A(x))"))
  res$pvalue <- c(res$pvalue, stats$pvalue)
  names(res$pvalue) <- names(res$statistic)
  
  if (stats$pvalue < alpha) {
    # use FP1(., A(x))
    res$power_best = fit_fp1a$power_best
    res$model_best = 5
    return(res)
  }
  
  # use linear(., A(x))
  res$power_best = matrix(c(NA, 1), ncol = 2)
  res$model_best = 6
  
  res
}

#' Function selection procedure based on information criteria
#' 
#' Used in [find_best_fp_step()] when `criterion = "aic"` or `"bic"`.
#' For parameter explanations, see [find_best_fp_step()]. All parameters 
#' captured by `...` are passed on to [fit_model()].
#' 
#' @details  
#' In case an information criterion is used to select the best model the 
#' selection procedure simply fits all relevant models and selects the best
#' one according to the given criterion. 
#' 
#' "Relevant" models for a given degree are the null model excluding the 
#' variable of interest, the linear model and all best FP models up to the 
#' specified degree. 
#' 
#' Note that the "best" FPx model used in this function are given by the models
#' using a FPx transformation for the variable of interest and having the 
#' highest likelihood of all such models given the current powers for all other
#' variables, as outlined in Section 4.8 of Royston and Sauerbrei (2008).
#' These best FPx models are computed in [find_best_fpm_step()].
#' Keep in mind that for a fixed number of degrees of freedom (i.e. fixed m),
#' the model with the highest likelihood is the same as the model with the best
#' information criterion of any kind since all the models share the same 
#' penalty term. 
#' 
#' When a variable is forced into the model by including it in `keep`, then 
#' this function will not exclude it from the model (by setting its power to 
#' `NA`), but will only choose its functional form. 
#' 
#' @return 
#' A list with several components:
#' 
#' * `acd`: logical indicating if an ACD transformation was applied for `xi`, 
#' i.e. `FALSE` in this case.
#' * `powers`: (best) fp powers investigated in step, indexing `metrics`. 
#' * `power_best`: a numeric vector with the best power found. The returned 
#' best power may be `NA`, indicating the variable has been removed from the 
#' model.
#' * `metrics`: a matrix with performance indices for all best models 
#' investigated. Same number of rows as, and indexed by, `powers`.
#' * `model_best`: row index of best model in `metrics`.
#' * `pvalue`: p-value for comparison of linear and null model, `NA` in this
#' case..
#' * `statistic`: test statistic used, depends on `ftest`, `NA` in this 
#' case.
#' 
#' @seealso 
#' [select_ra2()]
select_ic <- function(x, 
                      xi,
                      degree,
                      y, 
                      powers_current,
                      keep, 
                      criterion,
                      powers,   
                      acdx, 
                      ...) {
  
  if (degree < 1)
    return(NULL)
  
  fpmax = paste0("FP", degree)
  
  # output list
  res <- list(
    keep = xi %in% keep,
    acd = FALSE, 
    powers = NULL, 
    power_best = NULL, 
    metrics = NULL, 
    model_best = NULL, 
    statistic = NA, 
    pvalue = NA
  )
  
  # fit all relevant models
  fit_null <- fit_null_step(
      x = x, xi = xi, y = y, 
      powers_current = powers_current, powers = powers, acdx = acdx, ...
    )
  fit_lin <- fit_linear_step(
      x = x, xi = xi, y = y, 
      powers_current = powers_current, powers = powers, acdx = acdx, ...
    )
  
  fits_fpm = list()
  for (m in 1:degree) {
    fits_fpm[[sprintf("FP%g", m)]] <- find_best_fpm_step(
      x = x, xi = xi, degree = m, y = y, 
      powers_current = powers_current, powers = powers, acdx = acdx, ...
    )
  }
  
  # output summary - only output best fpm models
  len_max <- ncol(fits_fpm[[fpmax]]$powers)
  
  res$powers <- lapply(fits_fpm, function(x) {
    ensure_length(x$powers[x$model_best, , drop = FALSE], len_max)
  })
  res$powers <- do.call(
    rbind, 
    c(list(null = ensure_length(fit_null$powers, len_max), 
           linear = ensure_length(fit_lin$powers, len_max)), 
      res$powers)
  )
  
  res$metrics <- lapply(fits_fpm, function(x) {
    x$metrics[x$model_best, , drop = FALSE] 
  })
  res$metrics <- do.call(
    rbind, 
    c(list(null = fit_null$metrics, linear = fit_lin$metrics), res$metrics)
  )
  
  ind_select = 1:nrow(res$metrics)
  if (xi %in% keep) 
    # prevent selection of null model
    ind_select = 2:nrow(res$metrics)
  
  res$model_best = which.min(res$metrics[ind_select, tolower(criterion), 
                                         drop = TRUE])
  res$power_best = res$powers[res$model_best, , drop = FALSE]
  
  res
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

#' Helper function to ensure vectors have a specified length
#' 
#' Used to make sure dimensions of matrix rows match.
ensure_length <- function(x, 
                          size, 
                          fill = NA) {
  if (length(x) == size)
    return(x)
    
  x_new = rep(NA, size)
  x_new[1:length(x)] = x
  
  x_new
}
