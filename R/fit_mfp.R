#' Function to fit a model using the mfpa algorithm
#' 
#' Not exported. To be called from the [`mfpa()`] function. Most parameters
#' are explained in the documentation of `mfpa()`, but their form may differ
#' in this function. This function does not check its arguments and expects that 
#' its input is prepared in `mfpa()`.
#' 
#' @param x an input matrix of dimensions nobs x nvars. Does not contain 
#' intercept, but columns are already expanded into dummy variables as 
#' necessary. Data are assumed to be shifted and scaled. 
#' @param y a vector for the response variable or a `Surv` object.
#' @param weights a vector of observation weights of length nobs. 
#' @param offset a vector of length nobs of offsets.
#' @param cycles an integer, maximum number of iteration cycles (i.e. update
#' powers for all predictors). 
#' @param scale a numeric vector of length nvars of scaling factors. Not applied,
#' but re-ordered to conform to `xorder`.
#' @param shift a numeric vector of length nvars of shifts. Not applied, 
#' but re-ordered to conform to `xorder`.
#' @param df a numeric vector of length nvars of degrees of freedom.
#' @param center a logical vector of length nvars indicating if variables are 
#' to be centered.
#' @param family a character string representing a family object.
#' @param criterion a character string defining the criterion used to select 
#' variables and FP models of different degrees.
#' @param select a numeric vector of length nvars indicating significance levels
#' for backward elimination.
#' @param alpha a numeric vector of length nvars indicating significance levels 
#' for tests between FP models of different degrees. 
#' @param keep a character vector that with names of variables to be kept 
#' in the model. 
#' @param xorder a string determining the order of entry of the covariates
#' into the model-selection algorithm. 
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
#' @param verbose a logical; run in verbose mode.
#' 
#' @section Algorithm: 
#' 
#' * Step 1: order variables according to `xorder`. This step may involve 
#' fitting a regression model to determine order of significance. 
#' * Step 2: input data pre-processing. Setting initial powers for fractional 
#' polynomial terms, checking if acd transformation is required and allowed.
#' Note that the initial powers of all variables are always set to 1, and higher
#' fps are only evaluated in turn for each variables in the first cycle of the 
#' algorithm. See e.g. Sauerbrei and Royston (1999).
#' * Step 3: run mfp algorithm cycles. See [find_best_fp_cycle()] for more 
#' details.
#' * Step 4: fit final model using estimated powers.
#' 
#' @return 
#' See [mfpa()] for details on the returned object.
#' 
#' @references 
#' Sauerbrei, W. and Royston, P., 1999. \emph{Building multivariable prognostic 
#' and diagnostic models: transformation of the predictors by using fractional 
#' polynomials. J Roy Stat Soc a Sta, 162:71-94.}
#' 
#' @seealso 
#' [mfpa()], [find_best_fp_cycle()]
fit_mfp <- function(x, 
                    y, 
                    weights,
                    offset,
                    cycles, 
                    scale, 
                    shift,
                    df, 
                    center, 
                    family, 
                    criterion,
                    select, 
                    alpha, 
                    keep, 
                    xorder,
                    powers, 
                    method, 
                    strata, 
                    nocenter,
                    acdx, 
                    ftest,
                    control,
                    verbose) {
  
  variables_x <- colnames(x) 

  # step 1: order variables ----------------------------------------------------
  variables_ordered = order_variables(
    xorder = xorder, 
    x = x, y = y, family = family,  weights = weights, offset = offset, 
    strata = strata, method = method, control = control, nocenter = nocenter
  ) 

  # step 2: pre-process input --------------------------------------------------
  # named list of initial fp powers set to 1 ordered by xorder
  fp_powers <- setNames(as.list(rep(1, ncol(x))), variables_ordered)
  
  # name and reorder input vectors by xorder
  alpha <- setNames(alpha, variables_x)[variables_ordered]
  select <- setNames(select, variables_x)[variables_ordered]
  df <- setNames(df, variables_x)[variables_ordered]
  center <- setNames(center, variables_x)[variables_ordered]
  shift <- setNames(shift, variables_x)[variables_ordered]
  scale <- setNames(scale, variables_x)[variables_ordered]
  acdx <- setNames(acdx, variables_x)[variables_ordered]
  
  # force variables into the model by setting p-value to 1
  if (!is.null(keep)) {
    select[which(names(select) %in% keep)] <- 1
  }

  # acd transformation setup
  if (any(acdx == TRUE)) {
    # reset acdx of variables with less than 5 distinct values to False
    acdx <- reset_acd(x, acdx)
    
    # assign two powers to acd variables (1, NA): 
    # the first is for xi, and the second is for acd(xi). 
    # Initially NA is assigned to acd(xi), to be updated in step 3.
    variables_acd <- names(acdx)[acdx == TRUE]
    fp_powers_acd <- sapply(variables_acd, function(v) c(1, NA), 
                            simplify = FALSE, USE.NAMES = TRUE)
    # update initial powers 
    fp_powers <- modifyList(x = fp_powers, val = fp_powers_acd)
    # override df of acd variables by setting them to 4
    df[which(variables_ordered %in% variables_acd)] <- 4
  }

  # step 3: mfp cycles ---------------------------------------------------------
  # initialize cycle counter 
  j <- 1
  converged <- FALSE
  
  # run cycles and update the powers in each step
  while (j <= cycles) {
    if (verbose) {
      cat(sprintf("\ni Running MFP Cycle %d", j))
    }
    
    # estimated powers for the j-th cycle
    fp_powers_updated <- find_best_fp_cycle(
      x = x,
      y = y,
      fp_powers = fp_powers,
      df = df,
      weights = weights,
      offset = offset,
      family = family,
      criterion = criterion,
      select = select,
      alpha = alpha,
      keep = keep,
      powers = powers,
      ftest = ftest,
      control = control,
      rownames = rownames,
      strata = strata,
      nocenter = nocenter,
      method = method,
      acdx = acdx, 
      verbose = verbose
    )

    # check for convergence (i.e. no change in powers in model)
    if (identical(fp_powers, fp_powers_updated)) {
      converged <- TRUE
      cat(
        sprintf(
          "\ni Fractional polynomial fitting algorithm converged after %d cycles.\n", 
          j)
      )   
      break
    } else {
      # update the powers of the variables at the end of each cycle
      fp_powers <- fp_powers_updated
      j <- j + 1
    }
  }

  if (!converged) {
    warning(sprintf("i No convergence after %d cycles.", cycles), 
            "i Results of the last iteration reported.")
  }
  
  # step 4: fit final model with estimated functional forms --------------------
  # transform x using the final FP powers selected. 
  # x has already been shifted and scaled.
  X <- transform_matrix(
    x = x, power_list = fp_powers, center = center, acdx = acdx
  )

  modelfit <- fit_model(
    x = X, y = y, family = family, weights = weights, offset = offset,
    method = method, strata = strata, control = control,
    rownames = rownames, nocenter = nocenter
  )
  
  # create mfpa object ---------------------------------------------------------
  
  # common components for glms and cox
  fit <- list(
    coefficients = modelfit$fit$coefficients,
    residuals = modelfit$fit$residuals,
    linear.predictors = modelfit$fit$linear.predictors,
    weights = modelfit$fit$weights,
    prior.weights = modelfit$fit$prior.weights, 
    df.residual = modelfit$fit$df.residual,
    df.null = modelfit$fit$df.null,
    X = X, 
    # untransformed and scaled x for selected variables
    # selected means that not all powers are NA
    x = x[, names(fp_powers[!sapply(fp_powers, function(x) all(is.na(x)))]), 
          drop = F],
    y = y, 
    fp_terms = create_fp_terms(fp_powers, acdx,
                               df, select, alpha, criterion),
    transformations = data.frame(shift = shift, 
                                 scale = scale, 
                                 center = center),
    fp_powers = fp_powers,
    acd = acdx
  )
  
  # expand list to conform to glm or coxph objects
  if (family != "cox") {
    fit <- c(
      fit, 
      list(
        fitted.values = modelfit$fit$fitted.values,
        effects = modelfit$fit$effects,
        R = modelfit$fit$R, 
        rank = modelfit$fit$rank,
        qr = modelfit$fit$qr, 
        family = modelfit$fit$family,
        na.action = NULL,
        deviance = modelfit$fit$deviance,
        aic = modelfit$fit$aic,
        null.deviance = modelfit$fit$null.deviance, 
        prior.weights = modelfit$fit$prior.weights, 
        df.residual = modelfit$fit$df.residual,
        df.null = modelfit$fit$df.null
      )
    )
    class(fit) <- c("mfpa", "glm", "lm")
  } else {
    fit <- c(
      fit, 
      list( 
        family = "cox",
        var = modelfit$fit$var, 
        loglik = modelfit$fit$loglik,
        score = modelfit$fit$score,
        means = modelfit$fit$means,
        method = modelfit$fit$method,
        n = nrow(y),
        nevent = sum(y[, ncol(y)]),
        # na.action = options()$na.action, # default in coxph. 
        fail = if (is.character(modelfit$fit)) {
          "fail"
        } # to work on this later-might not be correct
      )
    )
    # add wald test in order to use summary.coxph()
    # this calculation follows coxph() source code
    nabeta <- !is.na(fit$coefficients)
    fit$wald.test <- survival::coxph.wtest(
      fit$var[nabeta, nabeta], fit$coefficients[nabeta],
      control$toler.chol
    )$test
    class(fit) <- c("mfpa", "coxph")
  }

  fit
}

#' Helper to order variables for mfpa algorithm
#' 
#' To be used in [fit_mfp()].
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
#' @param weights,offset parameters for both glm and Cox models, see either
#' [stats::glm()] or [survival::coxph()] depending on family. 
#' @param strata,method,control,nocenter Cox model specific parameters, see
#' [survival::coxph()].
#' 
#' @return 
#' A vector of the variable names in `x`, ordered according to `xorder`.
#' 
#' @import utils
order_variables <- function(xorder = "ascending",
                            x = NULL, 
                            ...) {
  names_ordered = colnames(x)
  
  if (xorder != "original") {
    names_ordered = order_variables_by_significance(xorder = xorder, x = x, ...)
  }
  
  names_ordered
}

#' @describeIn order_variables Order by significance in regression model.
order_variables_by_significance <- function(xorder, 
                                            x, 
                                            y,
                                            family,
                                            weights, 
                                            offset, 
                                            strata, 
                                            method, 
                                            control,
                                            nocenter) {
  
  # Convert factors to dummy variables if it exists...take it to the main function
  # x = model.matrix(as.formula(paste("~", paste(colnames(x), collapse="+"))),
  #                  data = as.data.frame(x))
  # save the family in factor form to be used later for gaussian
  fam <- family
  # number of rows of x or observations
  n <- dim(x)[1]
  
  if (family != "cox") {
    # glm.fit requires a function as a family not a character name
    if (is.character(family)) {
      family <- get(family, mode = "function", envir = parent.frame())
    }
    if (is.function(family)) family <- family()
    # family <- family()
    fit.full <- glm.fit(
      x = cbind(rep(1, n), x), y = y, weights = weights, offset = offset,
      family = family
    ) # full model
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
      rownames = rownames(x), nocenter = nocenter
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
        method = method, rownames = rownames(x),
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
                    "descending" = sort(p.value, decreasing = TRUE), 
                    sort(p.value, decreasing = FALSE) # default ascending
  )
  
  names(pvalues)
}

#' Helper to reset acd transformation for variables with few values
#' 
#' To be used in [fit_mfp()].
#' This function resets the `acdx` parameter (logical vector) of variables with
#' less than 5 distinct values to `FALSE`.
#' 
#' @param x a design matrix of dimension nobs x nvars where nvars is the number 
#' of predictors excluding an intercept.  
#' @param acdx a named logical vector of length nvars indicating continuous
#' variables to undergo the approximate cumulative distribution (ACD) 
#' transformation. May be ordered differently than the columns of `x`.
#' 
#' @return 
#' Logical vector of same length as `acdx`.
reset_acd <- function(x, 
                      acdx) {
  
  names_acd <- names(acdx)[which(acdx == TRUE)]
  
  # number of unique values of each column in acdx
  n_unique <- apply(x[, names_acd, drop = FALSE], 2, 
                    function(col) length(unique(col)))
  ind_reset <- which(n_unique < 5)
  
  if (length(ind_reset) > 0) {
    acdx[names_acd][ind_reset] <- FALSE
    warning("i For any variable with fewer than 5 unique values no acd transformation can be performed.\n", 
            sprintf("i The requested acd transform has been reset to FALSE for the following variables: %s.", 
                    paste0(names_acd[ind_reset], collapse = ", ")))
  }
  
  acdx
}

#' Helper to run cycles of the mfp algorithm 
#' 
#' This function estimates the best FP functions for all predictors in the 
#' current cycle. To be used in [fit_mfp()].
#' 
#' @details 
#' A cycle is defined as a pass through all of the predictors in the input
#' matrix `x`. A step is defined as the assessment of a single predictor. 
#' This algorithm is described in Sauerbrei et al (2006) and given in detail
#' in Royston and Sauerbrei (2008), in particular chapter 6.
#' 
#' Briefly, a cycle works as follows. It has as input the data matrix and a set 
#' of current best fp powers for each variable. In each step, the fp powers of a 
#' single covariate are assessed. 
#' To do this, all of the other variables are transformed as defined by the 
#' current powers (this is done in [transform_data_step()]) and the 
#' fp powers of the variable of interest are tested using the closed test 
#' procedure (done in [find_best_fp_step()]). Some of the adjustment variables 
#' may have their fp power set to `NA`, which means they were deselected from 
#' the working model and are not used in that step.
#' The results from all steps are collected and returned, completing a cycle.
#' 
#' Note that each cycle goes through all variables. In particular that means
#' even if a variable was eliminated in an earlier cycle, it will re-enter
#' each new cycle at the beginning to be evaluated to be added to the 
#' working model or to be eliminated again. 
#' 
#' The current adjustment set is always given through the current fp powers, 
#' which are updated in each step. 
#' 
#' @references 
#' Royston, P. and Sauerbrei, W., 2008. \emph{Multivariable Model - Building: 
#' A Pragmatic Approach to Regression Anaylsis based on Fractional Polynomials 
#' for Modelling Continuous Variables. John Wiley & Sons.}\cr
#' Sauerbrei, W., Meier-Hirmer, C., Benner, A. and Royston, P., 2006. 
#' \emph{Multivariable regression model building by using fractional 
#' polynomials: Description of SAS, STATA and R programs. 
#' Comput Stat Data Anal, 50(12): 3464-85.}
#' Sauerbrei, W. and Royston, P., 1999. \emph{Building multivariable prognostic 
#' and diagnostic models: transformation of the predictors by using fractional 
#' polynomials. J Roy Stat Soc a Sta, 162:71-94.}
find_best_fp_cycle <- function(x, 
                               y, 
                               fp_powers, 
                               df, 
                               weights, 
                               offset, 
                               family, 
                               criterion,
                               select, 
                               alpha, 
                               keep, 
                               powers, 
                               method, 
                               strata, 
                               verbose, 
                               ftest, 
                               control,
                               rownames,
                               nocenter, 
                               acdx) {
  
  # for printing distinguish between p-value and information criteria
  if (verbose) {
    print_mfp_summary(criterion, ftest = ftest)
  }
  
  names_x <- names(fp_powers)
  for (i in 1:ncol(x)) {
    # iterate through all predictors xi and update xi's best FP power
    # in terms of loglikelihood
    # the result can be NA (variable not significant), linear, FP1, FP2, ...
    # note that the adjustment set and powers are given by fp_powers
    # which is updated in each step
    fp_powers[[i]] <- find_best_fp_step(
      x = x, 
      y = y,
      xi = names_x[i],
      fp_powers = fp_powers,
      weights = weights,
      offset = offset,
      df = df[i], 
      select = select[i], 
      alpha = alpha[i],
      keep = keep,
      family = family,
      criterion = criterion,
      powers = powers,
      method = method,
      strata = strata,
      ftest = ftest,
      control = control,
      rownames = rownames,
      nocenter = nocenter,
      acdx = acdx,
      verbose = verbose
    )
  }
  
  fp_powers
}

#' Helper to calculates the final degrees of freedom for the selected model
#' 
#' To be used in [fit_mfp()].
#' 
#' @param p power of a variable.
#' 
#' @details 
#' An example calculation: if p is the power(s) and p = c(1,2), then df = 4 
#' but if x = NA then df = 0.
calculate_df <- function(p) {
  if (all(is.na(p))) {
    # df for unselected variable
    df <- 0
  } else {
    # Remove NAs in powers (2,NA) or (NA,1) etc. Happens because of acd
    p <- p[!is.na(p)]
    # df of linear function
    if (identical(p, 1)) {
      df <- 1
      # df of fpm. Note that df = 2m where m is the degree. if length = 1 then
      # degree = 1 and df = 2 etc
    } else {
      df <- 2 * length(p)
    }
  }
  
  df
}

#' Helper to convert a nested list with same or different length into a matrix
#' 
#' To be used in [fit_mfp()].
convert_powers_list_to_matrix <- function(power_list) {
  # Check the maximum number of powers i.e  FP2 has 2 while FP1 has 1
  psize <- sapply(power_list, length)
  maxp <- max(psize)
  # Create a new nested list of same length. This means that if FP1 was choosen
  # for x then the second power should be NA
  new_list_powers <- vector(mode = "list", length = length(power_list))
  for (i in 1:maxp) {
    new_list_powers[[i]] <- sapply(power_list, function(x) x[i])
  }
  # combine the powers and rename.
  matp <- do.call(cbind, new_list_powers)
  colnames(matp) <- paste0("power", 1:maxp)
  
  matp
}

#' Helper to create overview table of fp terms
#' 
#' To be used in [fit_mfp()].
#' 
#' @return 
#' Data.frame with overview of all fp terms. Each row represents a variable, 
#' with rownames giving the name of the variable. Variables with acd 
#' transformation are denoted by (A). The data.frame comprises the following 
#' columns: 
#' 
#' * `df_initial`: initial degrees of freedom. 
#' * `select`: significance level for backward elimination.
#' * `alpha`: significance level for fractional polyomial terms.
#' * `status`: encodes presence in the model as 1, and absence as 0.
#' * `df_final`: final estimated degrees of freedom.
#' * `powerN`: one or more columns with the final estimated fp powers (numbered
#' 1 to N).
create_fp_terms <- function(fp_powers, 
                            acdx, 
                            df,
                            select, 
                            alpha, 
                            criterion) {
  
  fp_terms <- data.frame(
    # initial degrees of freedom
    df_initial = df, 
    select = select, 
    alpha = alpha, 
    # presence / absence in final model encoded by NAs in fp_powers
    status = sapply(fp_powers, function(p) ifelse(all(is.na(p)), 0, 1)),
    # final degrees of freedom
    df_final = sapply(fp_powers, calculate_df), 
    convert_powers_list_to_matrix(fp_powers)
  )
  rownames(fp_terms) <- sapply(names(fp_powers), 
                               function(n) ifelse(acdx[n], paste0("(A)", n), n)) 
  
  if (criterion != "pvalue") {
    fp_terms$select <- criterion
  }
  
  fp_terms
}
