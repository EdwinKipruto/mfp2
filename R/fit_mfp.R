#' Function for fitting a model using the MFP or MFPA algorithm
#' 
#' This function is not exported and is intended to be called from 
#' the [`mfp2()`] function. While most parameters are explained in 
#' the documentation of `mfp2()`, their form may differ in this 
#' function. Note that this function does not check its arguments 
#' and expects that its input has been prepared in `mfp2()` function.
#' 
#' @param x an input matrix of dimensions nobs x nvars. Does not contain 
#' intercept, but columns are already expanded into dummy variables as 
#' necessary. Data are assumed to be shifted and scaled. 
#' @param y a vector for the response variable or a `Surv` object.
#' @param weights a vector of observation weights of length nobs. 
#' @param offset a vector of length nobs of offsets.
#' @param cycles an integer representing the maximum number of 
#' iteration cycles during which FP powers for all predictor are updated. 
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
#' @param keep a character vector with names of variables to be kept 
#' in the model. 
#' @param xorder a string determining the order of entry of the covariates
#' into the model-selection algorithm. 
#' @param powers a named list of numeric values that sets the permitted FP 
#' powers for each covariate.
#' @param method a character string specifying the method for tie handling in 
#' Cox regression model.
#' @param strata a factor of all possible combinations of stratification 
#' variables. Returned from [survival::strata()]. 
#' @param nocenter a numeric vector with a list of values for fitting Cox 
#' models. See [survival::coxph()] for details.
#' @param acdx a logical vector of length nvars indicating which continuous
#' variables should undergo the approximate cumulative distribution (ACD) 
#' transformation.
#' @param ftest a logical indicating the use of the F-test for Gaussian models.
#' @param control a list with parameters for model fit. See [survival::coxph()]
#' or [stats::glm()] for details. 
#' @param zero A logical vector indicating, by position, which columns of 
#' \code{x} should treat nonpositive values (zero or negative) as zero before 
#' transformation. Must be the same length and in the same order as the columns 
#' of \code{x}.
#' @param catzero A logical vector similar to \code{zero}, indicating which 
#' columns of \code{x} should treat nonpositive values as zero and also have a binary 
#' indicator automatically created and included in the model. Must match the length 
#' and order of the columns in \code{x}. See \code{\link{mfp2}} for details.
#' @param spike A logical vector indicating which columns of \code{x} contain
#' a spike at zero. The length and order of \code{spike} must match those of
#' the columns in \code{x}.
#' @param verbose Logical; if \code{TRUE}, additional information will be printed
#' during model fitting steps. Useful for understanding internal processing. 
#' Default is \code{FALSE}.
#' 
#' @section Algorithm: 
#' 
#' * Step 1: order variables according to `xorder`. This step may involve 
#' fitting a regression model to determine order of significance. 
#' * Step 2: input data pre-processing. Setting initial powers for fractional 
#' polynomial terms, checking if acd transformation is required and allowed.
#' Note that the initial powers of all variables are always set to 1, and higher
#' FPs are only evaluated in turn for each variables in the first cycle of the 
#' algorithm. See e.g. Sauerbrei and Royston (1999).
#' * Step 3: run mfp algorithm cycles. See [find_best_fp_cycle()] for more 
#' details.
#' * Step 4: fit final model using estimated powers.
#' 
#' @return 
#' See [mfp2()] for details on the returned object.
#' 
#' @references 
#' Sauerbrei, W. and Royston, P., 1999. \emph{Building multivariable prognostic 
#' and diagnostic models: transformation of the predictors by using fractional 
#' polynomials. J Roy Stat Soc a Sta, 162:71-94.}
#' 
#' @seealso 
#' [mfp2()], [find_best_fp_cycle()]
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
                    zero,
                    catzero,
                    spike,
                    verbose) {
  
  variables_x <- colnames(x)
  
  if (verbose) {
    cat("\ni Initial degrees of freedom:\n")
    print(matrix(df, nrow = 1, dimnames = list("df", variables_x)), 
          quote = FALSE)
  }
  
  # step 1: order variables ----------------------------------------------------
  variables_ordered <- variables_x
  
  # order only variables if they are 2 or more
  if (length(variables_x) > 1) {
  variables_ordered <- order_variables(
    xorder = xorder, 
    x = x, y = y, family = family,  weights = weights, offset = offset, 
    strata = strata, method = method, control = control, nocenter = nocenter
  ) 
  }
  
  if (verbose) {
    cat(sprintf("\ni Visiting order: %s\n", 
                paste0(variables_ordered, collapse = ", ")))
  }
  
  # step 2: pre-process input --------------------------------------------------
  # named list of initial fp powers set to 1 ordered by xorder
  powers_current <- setNames(as.list(rep(1, ncol(x))), variables_ordered)
  
  # name and reorder input vectors by xorder
  alpha <- setNames(alpha, variables_x)[variables_ordered]
  select <- setNames(select, variables_x)[variables_ordered]
  df <- setNames(df, variables_x)[variables_ordered]
  center <- setNames(center, variables_x)[variables_ordered]
  shift <- setNames(shift, variables_x)[variables_ordered]
  scale <- setNames(scale, variables_x)[variables_ordered]
  acdx <- setNames(acdx, variables_x)[variables_ordered]
  zero <- setNames(zero, variables_x)[variables_ordered]
  catzero <- setNames(catzero, variables_x)[variables_ordered]
  spike <- setNames(spike, variables_x)[variables_ordered]
  
  # powers is already named. so we need to sort it based on variables_ordered
  powers <- powers[variables_ordered]
  
  # Spike variables must have catzero = TRUE to force the binary variable into the
  # model, which is later handled by the Spike algorithm. Already done in mfp2()
  #  but repeated here in case a user calls this function directly
  catzero[spike] <- TRUE
  
  # Ensure that any variable marked as 'catzero' is also set to 'zero'
  # already done in mfp2() but repeated here in case a user calls this function

  zero[catzero] <- TRUE
  
  # Assert repeated powers of 1 is not supported. The program will fail when
  # degree is 1 in find_best_fpm_step() due to setdiff(v, c(1))
  diff_one <- unlist(lapply(powers, function(v) all(c(v %in% 1)) && length(v) != 1))
  if (any(diff_one)) {
    dfx <- df[diff_one]
    if (any(dfx > 1))
      stop(" The powers of some variables are repeated and all equal to 1. Set df for those 
           variables to 1 or change the powers.\n",
           sprintf("i This applies to the following variables: %s.", 
                   paste0(names(dfx > 1), collapse = ", ")), call. = FALSE)
  }
  
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
    powers_current_acd <- sapply(variables_acd, function(v) c(1, NA), 
                            simplify = FALSE, USE.NAMES = TRUE)
    # update initial powers 
    powers_current <- modifyList(x = powers_current, val = powers_current_acd)
    
    # override df of acd variables by setting them to 4
    df[which(variables_ordered %in% variables_acd)] <- 4
  }
  
  # Reset spike indicators for variables without zeros
  if (any(spike == TRUE)) {
    spike <- reset_spike(x, spike)
  }
  
  # Set default for spike decision: 
  # 1 = FPM/linear + binary
  # 2 = FPM/linear-usual FP algorithm---default
  # 3 = binary ony. 
  spike_decision <- rep(2, length(variables_ordered))
  names(spike_decision) <- variables_ordered

  #--- Create binary variables for catzero and convert all nonpositive to zero
  # to avoid repetition. The functions will handle infinite values caused
  # by zeros correctly
  
  # store original zero for use later
  zero_x <- zero
  
  # Reorder zero to match column order of x
  zero_aligned <- zero[colnames(x)]
  
  # Identify columns to apply the zeroing
  cols_to_zero <- which(zero_aligned)
  
  if (length(cols_to_zero) > 0) {
    # Set nonpositive values to 0
    for (j in cols_to_zero) {
      x[, j][x[, j] <= 0] <- 0
    }
    
    # reset zero to FALSE
    zero_x <- setNames(rep(FALSE, ncol(x)), variables_ordered)

  }
  
  # Create binary variables for catzero variables
  catzero_list <- lapply(names(catzero), function(v) {
    if (catzero[[v]]) {
      as.integer(x[, v] > 0)
    } else {
      NULL
    }
  })
  
  # Assign names to the list
  names(catzero_list) <- names(catzero)
  
  # TODO: like catzero_list, create acd_parameters to speed up computation
  acd_parameter <- lapply(seq_len(ncol(x)), function(i) {
    if (acdx[i]) {
      fit_acd(x[, i], powers = powers[[i]], zero = zero[i])
    } else {
      NULL
    }
  })
  names(acd_parameter) <- colnames(x)
  # step 3: mfp cycles ---------------------------------------------------------
  # initialize cycle counter 
  j <- 1
  converged <- FALSE
  
  # run cycles and update the powers in each step
  while (j <= cycles) {
    if (verbose) {
      cat("\n---------------------")
      cat(sprintf("\ni Running MFP Cycle %d\n", j))
      cat("---------------------\n")
    }
    
    # estimated powers for the j-th cycle
    fit_best_cycle <- find_best_fp_cycle(
      x = x,
      y = y,
      powers_current = powers_current,
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
      rownames = rownames(x),
      strata = strata,
      nocenter = nocenter,
      method = method,
      acdx = acdx, 
      zero = zero_x,# All are FALSE since the variables have been converted
      catzero = catzero_list, # A named list of binary variables
      spike = spike,
      spike_decision = spike_decision,
      acd_parameter = acd_parameter,
      verbose = verbose
    )
    
    powers_updated <- fit_best_cycle$powers_current
    spike_decision <- fit_best_cycle$spike_decision
    # check for convergence (i.e. no change in powers and variables in model)
    if (identical(powers_current, powers_updated)) {
      converged <- TRUE
      if (verbose) {
      cat(sprintf(
          "\ni Fractional polynomial fitting algorithm converged after %d cycles.\n", 
          j)) 
        }  
      break
    } else {
      # update the powers of the variables at the end of each cycle
      powers_current <- powers_updated
      j <- j + 1
    }
  }

  if (!converged) {
    warning(sprintf("i No convergence after %d cycles.", cycles), 
            "i Results of the last iteration reported.")
  }
  
  # step 4: fit final model with estimated functional forms --------------------
  # transform x using the final FP powers selected. 
  # x has already been shifted and scaled if necessary.
  
  # Apply backscaling only if at least one scaling factor is not equal to 1
  if (any(scale != 1)) {
  x <- backscale_matrix(x,scale)
  }
  
  data_transformed <- transform_matrix(
    x = x,  power_list = powers_current, center = center, acdx = acdx, 
    zero = zero, catzero = catzero
  )

  # fit model, and return full glm or coxph object
  modelfit <- fit_model(
    x = data_transformed$x_transformed,
    y = y, 
    family = family,
    weights = weights,
    offset = offset,
    method = method,
    strata = strata,
    control = control,
    rownames = rownames(data_transformed$x_transformed), 
    nocenter = nocenter,
    fast = FALSE
  )
  
  # create mfp2 object ---------------------------------------------------------
  
  # add components to fitted model object
  fit <- modifyList(
    modelfit$fit,
    list(
      centers = data_transformed$centers,
      acd_parameter = data_transformed$acd_parameter,
      convergence_mfp = converged,
      # untransformed and scaled x for selected variables
      # selected means that not all powers are NA
      x_original = x[, names(powers_current[!sapply(powers_current, function(x) all(is.na(x)))]), 
            drop = F],
      y = y, 
      fp_terms = create_fp_terms(powers_current, acdx, df, select, alpha, 
                                 criterion, zero, catzero),
      transformations = data.frame(shift = shift, 
                                   scale = scale, 
                                   center = center),
      fp_powers = powers_current,
      acd = acdx,
      zero = zero
    )
  )

  class(fit) <- c("mfp2", class(fit))

  fit
}

#' Helper to order variables for mfp2 algorithm
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
#' @param ... passed to `order_variables_by_significance`.
#' 
#' @return 
#' A vector of the variable names in `x`, ordered according to `xorder`.
#' 
#' @import utils
order_variables <- function(xorder = "ascending",
                            x = NULL, 
                            ...) {
  names_ordered <- colnames(x)
  
  if (xorder != "original") {
    names_ordered <- order_variables_by_significance(xorder = xorder, x = x, ...)
  }
  
  names_ordered
}

#' @describeIn order_variables Order by significance in regression model. The 
#' number of columns of `x` should be greater than 1 for Cox models.
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
    )
    
    # Rank of the glm full model
    p1 <- fit.full$rank
    
    # There is an additional scale parameter to estimate in OLS regression. The
    # binomial and Poisson regression models have no scale parameter.
    if (fam == "gaussian") p1 <- p1 + 1
    
    # loglikelihood of the full model calculated based on  aic = -2logL + 2k
    logl.full <- p1 - fit.full$aic / 2 
    
    # Deviance of the full model
    dev.full <- -2 * logl.full
    
    # we need to calculate p-values for each variable using likelihood ratio test
    varnames <- colnames(x)
    ns <- length(varnames)
    p.value <- loglikx <- dev <- df.reduced <- numeric(ns)
    names(p.value) <- names(dev) <- names(df.reduced) <- varnames
    for (i in 1:ns) {
      # remove one variable at a time and fit the reduced model. 
      # only works if you have more than one variable due to (-i)
      fit.reduced <- glm.fit(
        x = cbind(rep(1, n), x[, -i, drop = FALSE]),
        y = y,
        weights = weights, 
        offset = offset, 
        family = family
      )
      # calculate the deviance of the reduced model-model without the ith variable
      p2 <- fit.reduced$rank
      
      if (fam == "gaussian") p2 <- p2 + 1
      # loglikelihood of the reduced model
      logl.reduced <- p2 - fit.reduced$aic / 2 
      
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
      x = x, 
      y = y, 
      strata = strata,
      weights = weights,
      offset = offset,
      control = control, 
      method = method,
      rownames = rownames(x),
      nocenter = nocenter
    ) 
    
    # Degrees of freedom for the full cox model
    p1 <- fit.full$df
    
    # loglikelihood of the full cox model
    logl.full <- fit.full$logl
    
    # Deviance of the full cox model
    dev.full <- -2 * logl.full
    
    # we need to calculate p-values for each variable using likelihood ratio test
    varnames <- colnames(x)
    ns <- length(varnames)
    p.value <- loglikx <- dev <- df.reduced <- numeric(ns)
    names(p.value) <- names(dev) <- names(df.reduced) <- varnames
    for (i in 1:ns) {
      # remove one variable at a time and fit the reduced model
      fit.reduced <- fit_cox(
        x = x[, -i, drop = FALSE],
        y = y, strata = strata,
        weights = weights, 
        offset = offset, 
        control = control,
        method = method,
        rownames = rownames(x),
        nocenter = nocenter
      )
      # Degrees of freedom for the reduced model
      p2 <- fit.reduced$df
      
      # loglikelihood of the reduced model
      logl.reduced <- fit.reduced$logl
      
      # Deviance of the reduced model
      dev[i] <- -2 * logl.reduced
      
      # degrees of freedom of the reduced model
      df.reduced[i] <- p2
      
      # loglik difference: -2(logL.reduced - logL.full)
      teststatic <- -2 * logl.reduced + 2 * logl.full
      
      # calculate the p.value
      p.value[i] <- pchisq(teststatic, df = p1 - p2, lower.tail = FALSE)
    }
  }
  
  # Order the p-values based on xorder
  pvalues <- switch(xorder,
                    "descending" = sort(p.value, decreasing = TRUE),
                    "ascending"  = sort(p.value, decreasing = FALSE),
                     "original"  = p.value 
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
#' @param acdx a named logical vector of length nvars indicating which continuous
#' variables should undergo the approximate cumulative distribution (ACD) 
#' transformation. May be ordered differently than the columns of `x`.
#' 
#' @return 
#' Logical vector of same length as `acdx`.
reset_acd <- function(x, acdx) {
  # exit early if all acdx values are FALSE
  if (!any(acdx)) return(acdx)
  
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

#' Reset spike indicators for variables without zeros
#'
#' This function resets elements of a logical vector \code{spike} to \code{FALSE}
#' if the corresponding variables in \code{x} contain no zero values.
#'
#' @param x A data frame or matrix.
#' @param spike A logical vector indicating which columns of \code{x} are
#'   expected to contain a spike at zero. Its length and order must match
#'   the columns of \code{x}.
#'
#' @return A logical vector of the same length as \code{spike}, with entries
#'   reset to \code{FALSE} where the corresponding variable in \code{x}
#'   contains no zeros.
#'
#' @details
#' A warning is issued listing the variables for which the spike option
#' has been reset.
#'
#' @keywords internal
reset_spike <- function(x, spike) {
  # exit early if all spike values are FALSE
  if (!any(spike)) return(spike)
  
  names_spike <- names(spike)[spike]
  
  # identify columns without zeros
  has_zero <- apply(x[, names_spike, drop = FALSE], 2, function(col) any(col == 0))
  
  ind_reset <- which(!has_zero)
  
  if (length(ind_reset) > 0) {
    spike[names_spike][ind_reset] <- FALSE
    warning(sprintf(
      "The spike option has been reset to FALSE for the following variables because they contain no zero values: %s",
      paste(names_spike[ind_reset], collapse = ", ")
    ))
  }
  
  spike
}


#' Helper to run cycles of the mfp algorithm 
#' 
#' This function estimates the best FP functions for all predictors in the 
#' current cycle. To be used in [fit_mfp()].
#' 
#' @details 
#' A cycle is defined as a complete pass through all the predictors in the input
#' matrix `x`, while a step is defined as the assessment of a single predictor. 
#' This algorithm is described in Sauerbrei et al. (2006) and given in detail
#' in Royston and Sauerbrei (2008), in particular chapter 6.
#' 
#' Briefly, a cycle works as follows: it takes as input the data matrix along with
#' a set of current best fp powers for each variable. In each step, the fp
#' powers of a single covariate are assessed, while adjusting for other
#' covariates. Adjustment variables are transformed using their current
#' fp powers (this is done in [transform_data_step()]) and the fp powers 
#' of the variable of interest are tested using the closed test procedure
#' (conducted in [find_best_fp_step()]).
#' Some of the adjustment variables may have their fp power set to `NA`, 
#' which means they were not selected from the working model and are not used
#' in that step. The results from all steps are returned, completing a cycle.
#' 
#' Note that in each cycle every variable is evaluated.This includes variables
#' that may have been eliminated in previous cycles. They will re-enter each
#' new cycle for potential inclusion in the working model or to be re-evaluated
#' for elimination.
#' 
#' The current adjustment set is always given through the current fp powers, 
#' which are updated in each step (denoted as `powers_current`). 
#'
#' If \code{catzero} variables are supplied, the algorithm will automatically create 
#' the corresponding binary variables and include them in the model. Additionally, 
#' each binary variable and its associated continuous variable will be treated as 
#' one predictor, and they will be tested jointly for inclusion in the model.
#'  
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
#' 
#' @inheritParams fit_mfp
#' @param powers_current a list of length equal to the number of variables, 
#' indicating the fp powers to be used in the current step for all variables 
#' (except `xi`). 
#' @param catzero A named list of binary indicator variables of length \code{ncol(x)} 
#' for nonpositive values, created when specific variables are passed to the 
#' \code{catzero} argument of \code{fit_mfp}. If an element of the list is 
#' \code{NULL}, it indicates that the corresponding variable was not specified by
#' the user in the \code{catzero} argument of \code{fit_mfp}. Here, \code{catzero}
#' is a list of binary variables, not a named logical vector as in \code{fit_mfp}.
#' @param rownames passed to [survival::coxph.fit()].
#' 
#' @return 
#' A list of current FP powers
find_best_fp_cycle <- function(x, 
                               y, 
                               powers_current, 
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
                               zero,
                               catzero, 
                               spike,
                               spike_decision,
                               acd_parameter,
                               acdx) {
  
  names_x <- names(powers_current)
  
  for (i in 1:ncol(x)) {
    # iterate through all predictors xi and update xi's best FP power
    # in terms of loglikelihood
    # the result can be NA (variable not significant), linear, FP1, FP2, ...
    # note that the adjustment set and powers are given by powers_current
    # which is updated in each step
    fit_best_fp_step <- find_best_fp_step(
      x = x, 
      y = y,
      xi = names_x[i],
      powers_current = powers_current,
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
      zero = zero,
      catzero = catzero,
      spike = spike,
      spike_decision = spike_decision,
      acd_parameter = acd_parameter,
      verbose = verbose
    )
    # Update parameters
    powers_current[[i]] <- fit_best_fp_step$power_best
    spike_decision <- fit_best_fp_step$spike_decision
      
  }
  
  list(powers_current = powers_current, spike_decision = spike_decision)
}

#' Helper to calculates the final degrees of freedom for the selected model
#' 
#' To be used in [fit_mfp()].
#' 
#' @param p power of a variable.
#' 
#' @details 
#' An example calculation: if p is the power(s) and p = c(1,2), then df = 4 
#' but if p = NA then df = 0.
#' 
#' @return 
#' returns numeric value denoting the number of degrees of freedom (df).
calculate_df <- function(p) {
  if (all(is.na(p))) {
    # df for unselected variable
    df <- 0
  } else {
    # Remove NAs in powers (2,NA) or (NA,1) etc. Happens because of acd
    # make sure to drop names by using as.numeric
    p <- as.numeric(p[!is.na(p)])
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
#' 
#' @param power_list list of powers created in `fit_mfp()`.
#' 
#' @return 
#' a matrix.
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
#' Dataframe with overview of all fp terms. Each row represents a variable, 
#' with rownames giving the name of the variable. Variables with acd 
#' transformation are prefixed by `A_` by the `print` and `summary` methods. 
#' The dataframe comprises the following columns: 
#' 
#' * `df_initial`: initial degrees of freedom. 
#' * `select`: significance level for backward elimination.
#' * `alpha`: significance level for fractional polyomial terms.
#' * `selected`: logical value encoding presence in the model.
#' * `df_final`: final estimated degrees of freedom.
#' * `acd`: logical value encoding use of ACD transformation.
#' * `powerN`: one or more columns with the final estimated fp powers (numbered
#' 1 to N).
#' 
#' @inheritParams fit_mfp
#' @param fp_powers powers of the created FP terms.
create_fp_terms <- function(fp_powers, 
                            acdx, 
                            df,
                            select, 
                            alpha, 
                            criterion,
                            zero,
                            catzero) {
  
  fp_terms <- data.frame(
    # initial degrees of freedom
    df_initial = df, 
    select = select, 
    alpha = alpha, 
    acd = acdx, 
    zero = zero,
    catzero = catzero,
    # presence / absence in final model encoded by NAs in fp_powers
    selected = sapply(fp_powers, function(p) ifelse(all(is.na(p)), FALSE, TRUE)),
    # final degrees of freedom
    df_final = sapply(fp_powers, calculate_df), 
    convert_powers_list_to_matrix(fp_powers)
  )
  rownames(fp_terms) <- names(fp_powers)
  
  
  if (criterion != "pvalue") {
    fp_terms$select <- toupper(criterion)
    fp_terms$alpha <- toupper(criterion)
  }
  
  fp_terms
}


#' Backscale Columns of a Matrix (Internal)
#'
#' Multiplies each column of a numeric matrix by a corresponding scalar value 
#' from a named vector. Typically used to reverse prior scaling (i.e., backscaling).
#' This is an internal helper function and not intended for direct use by package 
#' users.
#'
#' @param x A numeric matrix with column names, or `NULL`.
#' @param scalex A named numeric vector. Each name must match a column name of `x`.
#'
#' @return A matrix with backscaled columns, or `NULL` if `x` is `NULL`.
#'
backscale_matrix <- function(x, scalex) {
  # If x is NULL, return NULL
  if (is.null(x)) {
    return(NULL)
  }
  
  # Check: x must be a matrix
  if (!is.matrix(x)) {
    stop("`x` must be a matrix.")
  }
  
  # Check: x must be numeric
  if (!is.numeric(x)) {
    stop("`x` must be a numeric matrix.")
  }
  
  # Check: column names must be present
  vnames <- colnames(x)
  if (is.null(vnames)) {
    stop("`x` must have column names.")
  }
  
  # Check: scalex must be a named numeric vector
  if (!is.numeric(scalex) || is.null(names(scalex))) {
    stop("`scalex` must be a named numeric vector.")
  }
  
  # Check: all columns in x must have matching names in scalex
  missing_cols <- setdiff(vnames, names(scalex))
  if (length(missing_cols) > 0) {
    stop("Missing scaling values for column(s): ", paste(missing_cols, collapse = ", "))
  }
  
  # Multiply each column by the corresponding scalar
  unscale_x <- sweep(x, 2, scalex[vnames], FUN = "*")
  return(unscale_x)
}

