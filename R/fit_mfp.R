#' @title Function to fit a model using the mfpa algorithm
#' 
#' @details 
#' Not exported. To be called from the [`mfpa()`] function. Most parameters
#' are explained in the documentation of `mfpa()`, but their form may differ
#' in this version. This function does not check its arguments and expects that 
#' its input is prepared in `mfpa()`.
#' 
#' @param x an input matrix of dimensions nobs x nvars. Does not contain 
#' intercept, but columns are already expanded into dummy variables as 
#' necessary.
#' @param y a vector for the response variable or a `Surv` object.
#' @param weights a vector of observation weights of length nobs. 
#' @param offset a vector of length nobs of offsets.
#' @param cycles an integer, maximum number of iteration cycles. 
#' @param scale a numeric vector of length nvars of scaling factors.
#' @param shift a numeric vector of length nvars of shifts.
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
#' * Step 3: run mfp algorithm cycles. 
#' * Step 4: fit final model using estimated powers.
#' 
#' @return 
#' See [mfpa()] for details on the returned object.
#' 
#' @seealso 
#' [mfpa()]
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
    
    # assign two powers to acd variables (1, NA). The first is for xi, and the
    # second is for acd(xi). NA has been assigned to acd(xi), which will be
    # updated in the MFP cycles.
    variables_acd <- names(acdx)[acdx == TRUE]
    fp_powers_acd <- sapply(variables_acd, function(v) c(1, NA), 
                            simplify = FALSE, USE.NAMES = TRUE)
    # update initial powers with acd having two powers
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
    fp_powers_updated <- iterate_find_best_model_fp(
      x = x, y = y,
      allpowers = fp_powers,
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

    # check for convergence
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
  X <- transform_x_fp(
    x = x, power.list = fp_powers, center = center, acdx = acdx
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

#' Function to create overview table of fp terms
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