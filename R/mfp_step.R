#' Function to estimate the best FP functions for a single variable
#' 
#' See [mfp2()] for a brief summary on the notation used here and 
#' [fit_mfp()] for an overview of the fitting procedure.  
#' 
#' @param x an input matrix of dimensions nobs x nvars. Does not contain 
#' intercept, but columns are already expanded into dummy variables as 
#' necessary. Data are assumed to be shifted and scaled. 
#' @param y a vector for the response variable or a `Surv` object.
#' @param xi a character string indicating the name of the current variable 
#' of interest, for which the best fractional polynomial transformation is
#' to be estimated in the current step. 
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
#' @param keep a character vector with names of variables to be kept 
#' in the model. 
#' @param powers a named list of numeric values that sets the permitted FP 
#' powers for each covariate.
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
#' @param zero a named logical vector
#' @param catzero a named logical vector
#' @param zero A named logical vector indicating, which columns of 
#' \code{x} should treat nonpositive values (zero or negative) as zero before 
#' transformation. Must be the same length as the columns of \code{x}.
#' @param catzero A named list of binary indicator variables of length \code{ncol(x)} 
#' for nonpositive values, created when specific variables are passed to the 
#' \code{catzero} argument of \code{fit_mfp}. If an element of the list is 
#' \code{NULL}, it indicates that the corresponding variable was not specified by
#' the user in the \code{catzero} argument of \code{fit_mfp}. Here, \code{catzero}
#' is a list of binary variables, not a named logical vector as in \code{fit_mfp}.
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
#' another df. This is also the case applied for categorical variables for 
#' which `df` are set to 1.
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
#' `NA` if variable is to be removed from the working model. Note that this 
#' vector may include up to two `NA` entries when ACD transformation is 
#' requested, but otherwise is either a vector with all numeric entries, or a 
#' single `NA`.
#' 
#' @references 
#' Royston, P. and Sauerbrei, W., 2008. \emph{Multivariable Model - Building: 
#' A Pragmatic Approach to Regression Anaylsis based on Fractional Polynomials 
#' for Modelling Continuous Variables. John Wiley & Sons.}\cr
#' 
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
                              zero,
                              catzero,# a named list of binary variables
                              verbose) {

  degree <- as.numeric(df / 2)
  
  # choose appropriate selection function
  if (df == 1) {
    # linear case 
    select_fct <- select_linear
  } else if (acdx[xi]) {
    # acd case 
    if (tolower(criterion) == "pvalue") {
      select_fct <- select_ra2_acd
    } else select_fct <- select_ic_acd
  } else {
    # usual mfp case 
    if (tolower(criterion) == "pvalue") {
      select_fct <- select_ra2
    } else select_fct <- select_ic
  }
  
  fit <- select_fct(
    x = x, xi = xi, keep = keep, degree = degree, acdx = acdx, 
    y = y, family = family, weights = weights, offset = offset, 
    powers_current = powers_current, powers = powers,  
    criterion = criterion, ftest = ftest, select = select, alpha = alpha,
    method = method, strata = strata, nocenter = nocenter, 
    control = control, rownames = rownames, zero = zero, catzero = catzero
  )
  
  if (verbose) {
    print_mfp_step(xi = xi, criterion = criterion, fit = fit)
  }
    
  # prepare power to return
  # remove trailing NAs, unless ACD is used
  power_best <- as.numeric(fit$power_best)
  
  if (!fit$acd) {
    power_best <- na.omit(power_best)
    if (length(power_best) == 0)
      power_best <- NA 
  } 
  
  # create names for powers
  names(power_best) <- name_transformed_variables(xi, 
                                                  length(power_best),
                                                  acd = acdx[xi])
  
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
#' * `zero`: Logical indicating whether a zero transformation was applied to \code{xi}. 
#'   In this case, nonpositive values of \code{xi} were set to zero before transformation, 
#'   and only positive values were transformed.
#' * `catzero`: Logical indicating whether a combination of a zero transformation 
#'   and a binary indicator variable was applied to \code{xi}. This means that 
#'   nonpositive values of \code{xi} were set to zero, only positive values were 
#'   transformed, and an additional binary variable was created to indicate 
#'   whether \code{xi} was positive or nonpositive.
#' 
#' @inheritParams find_best_fp_step
#' @param degree degrees of freedom for fp transformation of `xi`.
#' @param ... parameters passed to `fit_model()`.
find_best_fpm_step <- function(x, 
                               xi,
                               degree,
                               y, 
                               powers_current,
                               powers,  
                               acdx, 
                               family,
                               zero,
                               catzero, # a list of binary variables or null
                               ...) {
  # Number of events
  n_obs <- ifelse(family == "cox", sum(y[, 2]), nrow(x))
  
  if (degree == 1) {
    # remove linear model for normal data, but keep for acd transformation
    # powers = numeric(0) if powers = c(1,1)
    powers <- lapply(powers, function(v) setdiff(v, c(1)))
    
  }
    
  # generate FP data for x of interest (xi) and adjustment variables
  # Takes into account variables that should not be shifted thru 'zero'
  x_transformed <- transform_data_step(
    x = x, xi = xi, df = 2 * degree, powers_current = powers_current, 
    powers = powers, acdx = acdx, zero = zero, catzero = catzero
  )
  
  metrics <- list()
  for (i in seq_along(x_transformed$data_fp)) {
    # combine FP variables for x of interest with adjustment variables
    fit <- fit_model(
      x = cbind(x_transformed$data_fp[[i]], x_transformed$data_adj), y = y, # catzero plays a role here thru data_fp and data_adj
      family = family, ...
    )
    
    # use degree many additional degrees of freedom
    # TODO: here we should likely add additional df for each fp term in the model
    # using calculate_number_fp_powers
    # note: this doesn't change WHICH model is the best, since all use the 
    # same additional df, but it may affect further comparison between 
    # different fp models
    p <- sprintf("%g", x_transformed$powers_fp[i, , drop = TRUE])
    
    # respect acd
    if (acdx[xi])
      p[length(p)] <- sprintf("A(%s)", p[length(p)])
    
    metrics[[paste(p, collapse = " ")]] <- calculate_model_metrics(
      fit, n_obs, degree # fit will include df of binary when catzero is used so this part does not change
    )
  }
  
  metrics <- do.call(rbind, metrics) 
  model_best <- as.numeric(which.max(metrics[, "logl"]))
   
  list(
    acd = acdx[xi],
    powers = x_transformed$powers_fp, 
    power_best = x_transformed$powers_fp[model_best, , drop = TRUE], 
    metrics = metrics, 
    model_best = model_best,
    zero = zero[xi],
    catzero = ifelse(!is.null(catzero[[xi]]), TRUE, FALSE)
  ) 
}

#' Function to fit a null model excluding variable of interest
#' 
#' "Null" model here refers to a model which does not include the variable 
#' of interest `xi`. 
#' For parameter definitions, see [find_best_fp_step()]. All parameters 
#' captured by `...` are passed on to [fit_model()].
#' 
#' @return 
#' A list with two entries: 
#' 
#' * `powers`: FP power(s) of `xi` in fitted model - in this case `NA`.
#' * `metrics`: A matrix with performance indices for fitted model.
#' 
#' @inheritParams find_best_fp_step
#' @param ... Parameters passed to `fit_model()`.
fit_null_step <- function(x, 
                          xi, 
                          y, 
                          powers_current,
                          powers,
                          acdx, 
                          family,
                          zero,
                          catzero,
                          ...) {
  
  # Number of events
  n_obs <- ifelse(family == "cox", sum(y[, 2]), nrow(x))
  
  # transform all data as given by current working model
  # set variable of interest to linear term only
  x_transformed <- transform_data_step(
    x = x, xi = xi, df = 1, powers_current = powers_current, acdx = acdx, 
    powers = powers, zero = zero, catzero = catzero
  ) 
  
  # fit null model
  # i.e. a model that does not contain xi but only adjustment variables
  # In addition, adjustment model can be NULL, so we have intercept only
  model_null <- fit_model(x = x_transformed$data_adj, y = y,
                          family = family, ...) 
  
  list(
    powers = NA,
    metrics = rbind(null = calculate_model_metrics(model_null, n_obs))  
  )
}

#' Function to fit linear model for variable of interest
#' 
#' "Linear" model here refers to a model that includes the variable 
#' of interest \code{xi} with an FP (fractional polynomial) power of 1. 
#' Note that \code{xi} may be ACD-transformed if indicated by \code{acdx[xi]}. 
#' If the variable was passed through the \code{catzero} argument in \code{mfp2()}, 
#' both the continuous variable and its corresponding binary indicator 
#' will be included in the model as linear terms.
#' For parameter definitions, see \code{\link{find_best_fp_step}}. 
#' All parameters captured by \code{...} are passed to \code{\link{fit_model}}.
#' 
#' @return 
#' A list with two entries: 
#' 
#' * `powers`: FP power(s) of `xi` (or its ACD transformation) in fitted model.
#' * `metrics`: A matrix with performance indices for fitted model.
#' 
#' @inheritParams find_best_fp_step
#' @param ... Parameters passed to `fit_model()`.
fit_linear_step <- function(x, 
                            xi, 
                            y, 
                            powers_current,
                            powers,
                            acdx, 
                            family,
                            zero,
                            catzero, 
                            ...) {
  # Number of events in survival models or observation in GLM 
  n_obs <- ifelse(family == "cox", sum(y[, 2]), nrow(x))
  
  # transform all data as given by current working model
  # set variable of interest to linear term only
  x_transformed <- transform_data_step(
    x = x, xi = xi, df = 1, powers_current = powers_current, acdx = acdx,
    powers = powers, zero = zero, catzero = catzero
  ) 
  
  # fit a model based on the assumption that xi is linear. If catzero, its
  # corresponding binary variable will be included in the model
  model_linear <- fit_model(
    x = cbind(x_transformed$data_fp[[1]], x_transformed$data_adj), y = y,
    family = family, ...
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


#' Helper function to select between null and linear term for a single variable
#' 
#' To be used in [find_best_fp_step()]. Only used if `df = 1` for a variable.
#' Handles all criteria for selection.
#' For parameter explanations, see [find_best_fp_step()]. All parameters 
#' captured by `...` are passed to [fit_model()].
#' 
#' @details 
#' This function assesses a single variable of interest `xi` regarding its
#' functional form in the current working model as indicated by
#' `powers_current`, with the choice between excluding `xi` ("null model") and
#' including a linear term ("linear fp") for `xi`.
#' 
#' Note that this function handles an ACD transformation for `xi` as well. 
#' 
#' When a variable is forced into the model by including it in `keep`, then 
#' this function will not exclude it from the model (by setting its power to 
#' `NA`), but will only choose the linear model. 
#' 
#' @return 
#' A list with several components:
#' 
#' * `keep`: logical indicating if `xi` is forced into model.
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
#' * `zero`: Logical indicating whether a zero transformation was applied to \code{xi}. 
#'   In this case, nonpositive values of \code{xi} were set to zero before transformation, 
#'   and only positive values were transformed.
#' * `catzero`: Logical indicating whether a combination of a zero transformation 
#'   and a binary indicator variable was applied to \code{xi}. This means that 
#'   nonpositive values of \code{xi} were set to zero, only positive values were 
#'   transformed, and an additional binary variable was created to indicate 
#'   whether \code{xi} was positive or nonpositive.
#' @param degree not used.
#' @param ... passed to fitting functions. 
#' @inheritParams find_best_fp_step 
select_linear <- function(x, 
                          xi,
                          keep, 
                          degree,
                          acdx, 
                          y, 
                          powers_current,
                          powers,  
                          criterion,
                          ftest, 
                          select, 
                          alpha,
                          family,
                          zero,
                          catzero,
                          ...) {
  
  n_obs <- ifelse(family == "cox", sum(y[, 2]), nrow(x))
  
  # Model 1: Null model
  fit_null <- fit_null_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx,
    family = family, zero = zero, catzero = catzero, ...
  )
  
  # Model 2: Linear model
  fit_linear <- fit_linear_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx, 
    family = family, zero = zero, catzero = catzero, ...
  )
  
  # Extract powers and metrics of interest
  powers <- rbind(fit_null$powers, fit_linear$powers)
  metrics <- rbind(fit_null$metrics, fit_linear$metrics)
  
  # compute F or Chi-square statistic between the two models
  if (ftest) {
    # note that ftest is only TRUE if model is gaussian
    stats <- calculate_f_test(
      deviances = metrics[, "deviance_gaussian"], 
      dfs_resid = metrics[, "df_resid"],
      n_obs = n_obs
    )
  } else {
    stats <- calculate_lr_test(metrics[, "logl"], metrics[, "df"])
  }
  
  # Compute the corresponding p-value
  pvalue <- stats$pvalue
  names(pvalue) <- c("null vs linear")
  statistic <- stats$statistic 
  names(statistic) <- c("null vs linear")
  
  # Check whether the variable should be forced into the model; index 1 denotes
  # a null, while 2 denotes a linear model
  if (xi %in% keep) {
    model_best <- 2
  } else {
    model_best <- switch(
    tolower(criterion), 
    "pvalue" = ifelse(pvalue > select, 1, 2), 
    "aic" = which.min(metrics[, "aic", drop = TRUE]), 
    "bic" = which.min(metrics[, "bic", drop = TRUE])
  ) 
  }
  
  list(
    keep = xi %in% keep, 
    acd = acdx[xi],
    powers = powers,
    power_best = powers[model_best, ],
    metrics = metrics,
    model_best = model_best,
    pvalue = pvalue,
    statistic = statistic,
    zero = zero[xi],
    catzero = ifelse(!is.null(catzero[[xi]]), TRUE, FALSE)
  )
}

#' Function selection procedure based on closed testing procedure
#' 
#' Used in [find_best_fp_step()] when `criterion = "pvalue"`.
#' For parameter explanations, see [find_best_fp_step()]. All parameters 
#' captured by `...` are passed to [fit_model()].
#' 
#' @details  
#' In case `criterion = "pvalue"` the function selection procedure as outlined 
#' in Chapters 4 and 6 of Royston and Sauerbrei (2008) is used. 
#' 
#' * \emph{Step 1}: Test the best FP\emph{m} function against a null model at the 
#' significance level specified by \code{select}, using 2\emph{m} degrees of freedom. 
#' If the test is not significant, the variable is excluded. Otherwise, proceed
#' to Step 2.
#' * \emph{Step 2}: Test the best FP\emph{m} function against a linear model at 
#' the significance level specified by \code{alpha}, using 2\emph{m}-1 degrees 
#' of freedom. If the test is not significant, select the linear model. 
#' Otherwise, proceed to Step 3.
#' * \emph{Step 3}: Test the best FP\emph{m} function against the best FP1 model 
#' at the significance level specified by \code{alpha}, using 2\emph{m}-2 degrees
#' of freedom. If the test is not significant, retain the best FP1 model. Otherwise,
#' repeat this step by comparing FP\emph{m} to all remaining lower-order FP models, 
#' down to FP\emph{m}–1, which is tested with 2 degrees of freedom. If the final 
#' test is not significant, retain the best FP\emph{m}–1 model; otherwise, 
#' retain the best FP\emph{m} model.
#' 
#' Note that the "best" FP\emph{x} model used in each step refers to the model 
#' that applies an FP\emph{x} transformation to the variable of interest and 
#' achieves the highest likelihood among all such models, given the current 
#' power transformations for all other variables. This procedure is described 
#' in Section 4.8 of Royston and Sauerbrei (2008). The best FP\emph{x} models 
#' are computed by \code{\link{find_best_fpm_step}}.
#' 
#' When a variable is forced into the model by including it in the \code{keep} 
#' argument of \code{mfp2()}, this function will not exclude it (i.e., will not 
#' set its power to \code{NA}), but will instead select its functional form.
#' 
#' @return 
#' A list with several components:
#' 
#' * `keep`: logical indicating if `xi` is forced into model.
#' * `acd`: logical indicating if an ACD transformation was applied for `xi`, 
#' i.e. `FALSE` in this case.
#' * `powers`: (best) fp powers investigated in step, indexing `metrics`. 
#' Always starts with highest power, then null, then linear, then FP in 
#' increasing degree (e.g. FP2, null, linear, FP1).
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
#'  
#' @param degree integer > 0 giving the degree for the FP transformation. 
#' @param ... passed to fitting functions [fit_model()]. 
#' @inheritParams find_best_fp_step
select_ra2 <- function(x, 
                       xi,
                       keep, 
                       degree,
                       acdx, 
                       y, 
                       powers_current,
                       powers,  
                       criterion,
                       ftest, 
                       select, 
                       alpha, 
                       family,
                       zero,
                       catzero,
                       ...) {
  
  if (degree < 1) {
    return(NULL)
  }
  
  # Number of events in survival models or observations in GLM
  n_obs <- ifelse(family == "cox", sum(y[, 2]), nrow(x))
  
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
    calculate_test <- function(metrics, n_obs) {
      calculate_lr_test(
        logl = metrics[, "logl", drop = TRUE], 
        dfs = metrics[, "df", drop = TRUE] 
      )
    }
  }

  fpmax <- paste0("FP", degree)
  
  # output list
  res <- list(
    keep = xi %in% keep, 
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
    powers_current = powers_current, powers = powers, acdx = acdx, 
    family = family, zero = zero, catzero = catzero, ...
  )
  fit_null <- fit_null_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx, 
    family = family, zero = zero, catzero = catzero, ...
  )
  res$metrics <- rbind(
    fit_fpmax$metrics[fit_fpmax$model_best, ],
    fit_null$metrics
  )
  rownames(res$metrics) <- c(fpmax, "null")
  res$powers <- rbind(fit_fpmax$power_best, NA)
  
  # Test 1: test for overall significance (null vs best FPm)
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
  
  # Test 2: test for non-linearity (linear vs best FPm)
  # df for tests are degree * 2 - 1
  fit_lin <- fit_linear_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx, 
    family = family, zero = zero, catzero = catzero, ...
  )
  
  old_names <- rownames(res$metrics)
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
  
  # Test 3: test for complexity of the functions (best FP1 vs best FPm)
  # do this for all fps with lower degrees
  # dfs for tests are decreasing
  if (degree > 1) {
    
    for (current_degree in 1:(degree - 1)) {
      fpm = paste0("FP", current_degree)
      
      fit_fpm <- find_best_fpm_step(
        x = x, xi = xi, degree = current_degree, y = y, 
        powers_current = powers_current, powers = powers, acdx = acdx,
        family = family, zero = zero, catzero = catzero, ...
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
  res$power_best <- fit_fpmax$powers[fit_fpmax$model_best, , drop = FALSE]
  res$model_best <- 1
  
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
#' the corresponding section in the documentation of [mfp2()].
#' 
#' When a variable is forced into the model by including it in `keep`, then 
#' this function will not exclude it from the model (by setting its power to 
#' `NA`), but will only choose its functional form. 
#' 
#' @return 
#' A list with several components:
#' 
#' * `keep`: logical indicating if `xi` is forced into model.
#' * `acd`: logical indicating if an ACD transformation was applied for `xi`, 
#' i.e. `FALSE` in this case.
#' * `powers`: (best) fp powers investigated in step, indexing `metrics`. 
#' Ordering: FP1(x, A(x)), null, linear, FP1(x, .), linear(., A(x)), 
#' FP1(., A(x)).
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
#'
#' @param degree integer > 0 giving the degree for the FP transformation. 
#' @param ... passed to fitting functions. 
#' @inheritParams find_best_fp_step
select_ra2_acd <- function(x, 
                           xi,
                           keep, 
                           degree,
                           acdx, 
                           y, 
                           powers_current,
                           powers,  
                           criterion,
                           ftest, 
                           select, 
                           alpha, 
                           family,
                           zero,
                           catzero,
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
    calculate_test <- function(metrics, n_obs) {
      calculate_lr_test(
        logl = metrics[, "logl", drop = TRUE], 
        dfs = metrics[, "df", drop = TRUE] 
      )
    }
  }
  
  #n_obs <- nrow(x)
  n_obs <- ifelse(family == "cox", sum(y[, 2]), nrow(x))
  
  fpmax <- "FP1(x, A(x))"
  acdx_reset_xi <- acdx
  acdx_reset_xi[xi] = FALSE
  
  # output list
  res <- list(
    keep = xi %in% keep,
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
    powers_current = powers_current, powers = powers, acdx = acdx, 
    family = family,zero = zero, catzero = catzero, ...
  )
  fit_null <- fit_null_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx,
    family = family, zero = zero, catzero = catzero, ...
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
    powers_current = powers_current, powers = powers, acdx = acdx_reset_xi,
    family = family,zero = zero, catzero = catzero, ...
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
    powers_current = powers_current, powers = powers, acdx = acdx_reset_xi,
    family = family,zero = zero, catzero = catzero, ...
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
    powers_current = powers_current, powers = powers, acdx = acdx, 
    family = family, zero = zero, catzero = catzero,...
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
    powers_current = powers_current, powers = powers, acdx = acdx,
    family = family, zero = zero, catzero = catzero, ...
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
#' In case an ACD transformation is requested, then the models assessed 
#' are the null model, the linear model in x and A(x), the best FP1 models in 
#' x and A(x), and the best FP1(x, A(x)) model.
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
#' * `keep`: logical indicating if `xi` is forced into model.
#' * `acd`: logical indicating if an ACD transformation was applied for `xi`, 
#' i.e. `FALSE` in this case.
#' * `powers`: (best) fp powers investigated in step, indexing `metrics`. 
#' Ordered by increasing complexity, i.e. null, linear, FP1, FP2 and so on.
#' For ACD transformation, it is null, linear, linear(., A(x)), FP1(x, .),
#' FP1(., A(x)) and FP1(x, A(x)).
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
#' 
#' @param degree integer > 0 giving the degree for the FP transformation. 
#' @param ... passed to fitting functions. 
#' @inheritParams find_best_fp_step
select_ic <- function(x, 
                      xi,
                      keep, 
                      degree,
                      acdx, 
                      y, 
                      powers_current,
                      powers,  
                      criterion,
                      ftest, 
                      select, 
                      alpha, 
                      family,
                      zero,
                      catzero,
                      ...) {
  
  if (degree < 1)
    return(NULL)
  
  fpmax <- paste0("FP", degree)
  
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
  # Null Model
  fit_null <- fit_null_step(
      x = x, xi = xi, y = y, 
      powers_current = powers_current, powers = powers, acdx = acdx,
      family = family,zero = zero, catzero = catzero, ...
    )
  # Linear model
  fit_lin <- fit_linear_step(
      x = x, xi = xi, y = y, 
      powers_current = powers_current, powers = powers, acdx = acdx,
      family = family, zero = zero, catzero = catzero, ...
    )
  
  # All FPm models
  fits_fpm <- list()
  
  for (m in 1:degree) {
    fits_fpm[[sprintf("FP%g", m)]] <- find_best_fpm_step(
      x = x, xi = xi, degree = m, y = y, 
      powers_current = powers_current, powers = powers, acdx = acdx,
      family = family,zero = zero, catzero = catzero, ...
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
  rownames(res$metrics) <- c("null", "linear", names(fits_fpm))
  
  if (xi %in% keep) { 
    # prevent selection of null model
    ind_select <- 2:nrow(res$metrics)
    res$model_best <- which.min(res$metrics[ind_select, tolower(criterion), 
                                           drop = TRUE])
    # shift the positions by 1 since null was omitted when finding the index
    res$model_best <- res$model_best + 1
  }else{
    ind_select <- 1:nrow(res$metrics)
    res$model_best <- which.min(res$metrics[ind_select, tolower(criterion), 
                                            drop = TRUE])
  }

  res$power_best <- res$powers[res$model_best, , drop = FALSE]
  
  res
}

#' @describeIn select_ic Function to select ACD based transformation.
select_ic_acd <- function(x, 
                          xi,
                          keep, 
                          degree,
                          acdx, 
                          y, 
                          powers_current,
                          powers,  
                          criterion,
                          ftest, 
                          select, 
                          alpha, 
                          family,
                          zero,
                          catzero,
                          ...) {
  
  acdx_reset_xi <- acdx
  acdx_reset_xi[xi] <- FALSE
  
  # output list
  res <- list(
    keep = xi %in% keep,
    acd = TRUE, 
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
    powers_current = powers_current, powers = powers, acdx = acdx,
    family = family, zero = zero, catzero = catzero, ...
  )
  fit_lin <- fit_linear_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx_reset_xi,
    family = family, zero = zero, catzero = catzero, ...
  )
  fit_lina <- fit_linear_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx,
    family = family, zero = zero, catzero = catzero, ...
  )
  
  fits <- list(
    "FP1(x, .)" = find_best_fpm_step(
      x = x, xi = xi, degree = 1, y = y, 
      powers_current = powers_current, powers = powers, acdx = acdx_reset_xi,
      family = family, zero = zero, catzero = catzero, ...
    ), 
    "FP1(., A(x))" = find_best_fpm_step(
      x = x, xi = xi, degree = 1, y = y, 
      powers_current = powers_current, powers = powers, acdx = acdx,
      family = family, ...
    ), 
    "FP1(x, A(x))" = find_best_fpm_step(
      x = x, xi = xi, degree = 2, y = y, 
      powers_current = powers_current, powers = powers, acdx = acdx,
      family = family,zero = zero, catzero = catzero, ...
    )
  )

  
  # output summary - only output best fpm models
  len_max <- 2
  
  res$powers <- lapply(fits, function(x) {
    ensure_length(x$powers[x$model_best, , drop = FALSE], len_max)
  })
  res$powers <- do.call(
    rbind, 
    c(list("null" = ensure_length(fit_null$powers, len_max), 
           "linear" = ensure_length(fit_lin$powers, len_max), 
           "linear(., A(x))" = ensure_length(fit_lina$powers, len_max)), 
      res$powers)
  )
  
  res$metrics <- lapply(fits, function(x) {
    x$metrics[x$model_best, , drop = FALSE] 
  })
  res$metrics <- do.call(
    rbind, 
    c(list(
      "null" = fit_null$metrics, 
      "linear" = fit_lin$metrics, 
      "linear(., A(x))" = fit_lina$metrics), 
      res$metrics)
  )
  rownames(res$metrics) <- c("null", "linear", "linear(., A(x))", names(fits))
  
  ind_select = 1:nrow(res$metrics)
  if (xi %in% keep) {
    # prevent selection of null model
    ind_select <- 2:nrow(res$metrics)
    res$model_best <- which.min(res$metrics[ind_select, tolower(criterion), 
                                           drop = TRUE])
    res$model_best <- res$model_best + 1
  }else{
    res$model_best = which.min(res$metrics[ind_select, tolower(criterion), 
                                           drop = TRUE])
  }
  

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
#' @param zero named logical vector of length ncol(x)
#' @param catzero A named list of binary indicator variables of length \code{ncol(x)} 
#' for nonpositive values, created when specific variables are passed to the 
#' \code{catzero} argument of \code{fit_mfp}. If an element of the list is 
#' \code{NULL}, it indicates that the corresponding variable was not specified by
#' the user in the \code{catzero} argument of \code{fit_mfp}. Here, \code{catzero}
#' is a list of binary variables, not a named logical vector as in \code{fit_mfp}.
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
                                acdx,
                                zero,
                                catzero) { 
  
  # sort x based on the names of the powers
  names_powers_current <- names(powers_current)
  x <- x[, names_powers_current, drop = FALSE]
  
  # extract the matrix of adjustment variables
  vars_adj <- names_powers_current[!(names_powers_current %in% xi)]
  x_adj <- x[, vars_adj, drop = FALSE]
  
  # transformations of adjustment variables
  powers_adj <- powers_current[vars_adj]
  acdx_adj <- unname(acdx[vars_adj])
  zero_adj <- unname(zero[vars_adj])
  
  # generate adjustment data 
  # check whether all adjustment powers = NA or length(vars_adj)=0 for 1 variable
  # in the model (univariable fp) 
  if (all(is.na(unlist(powers_adj, use.names = FALSE))) || length(vars_adj) == 0) {
    # all adjustment variables were eliminated in MFP backfitting process
    data_adj <- NULL
    powers_adj <- NULL
  } else {
    data_adj <- vector(mode = "list", length = ncol(x_adj))
    
    for (i in 1:ncol(x_adj)) {
      if (acdx_adj[i]) {
        data_adj[[i]] <- transform_vector_acd(
          x = x_adj[, i, drop = TRUE], 
          power = powers_adj[[i]],
          zero = zero_adj[i],
          powers = powers[[vars_adj[i]]] # powers needed by acd() function
          
        )$acd
      } else {
        data_adj[[i]] <- transform_vector_fp(
          x = x_adj[, i, drop = TRUE], power = powers_adj[[i]], 
          zero = zero_adj[i]
        )
      }
    }
    
    # combine into data.frame
    # note some variables may have been extended to more than one column
    # in the loop above due to transformation
    data_adj <- do.call(cbind, data_adj) 
    
    # Add catzero binary indicators for adjustment variables
    catzero_adj <- catzero[vars_adj]
    catzero_adj <- catzero_adj[!sapply(catzero_adj, is.null)]  # keep only non-null
    
    if (length(catzero_adj) > 0) {
      data_adj <- cbind(data_adj, do.call(cbind, catzero_adj))
    }
    
    # assign arbitrary names to adjustment matrix 
    # the names of adjustment variables at this stage are not relevant
    colnames(data_adj) <- paste0("var_adj", 1:ncol(data_adj))
  }
  
  # generate fp data
  data_xi <- x[, xi, drop = TRUE]
  
  if (length(unique(data_xi)) <= 3) {
    # if a variable has less than 4 levels we do not generate FP data
    # but set nonpositive values to zero, when zero argument is used
    if (zero[xi]) {
      data_xi[data_xi <= 0] <- 0 
    }
    
    # add catzero variable if not null; if null cbind will remove it
    data_xi <- cbind(catzero[[xi]], data_xi)
    
    data_fp <- list(data_xi)
    powers_fp <- 1
  } else {
    if (acdx[xi]) {
      # when df = 1 -> degree = 0
      # and we return the data unchanged, i.e. with power = 1
      fpd <- generate_transformations_acd(data_xi, degree = floor(df / 2),
                                          powers = powers[[xi]], 
                                          zero = zero[xi], 
                                          catzero = catzero[[xi]])
    } else {
      # note that degree is df / 2
      fpd <- generate_transformations_fp(data_xi, degree = floor(df / 2),
                                         powers = powers[[xi]], 
                                         zero = zero[xi],
                                         catzero = catzero[[xi]])
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
#' 
#' @param x input vector or matrix. 
#' @param size length or size of `x` which is desired.
#' @param fill value to fill in if `x` is not of desired length or size.
ensure_length <- function(x, size, fill = NA) {
  if (length(x) == size)
    return(x)
    
  x_new = rep(NA, size)
  x_new[1:length(x)] = x
  
  x_new
}
