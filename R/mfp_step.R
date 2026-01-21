#' Function to estimate the best FP functions for a single variable
#' 
#' See \code{mfp2()} for a brief summary on the notation used here and 
#' \code{fit_mfp()} for an overview of the fitting procedure.  
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
#' @param family Either a character string naming the family (e.g., "gaussian", "binomial", "cox") 
#'   or a function that returns a GLM family object (e.g., stats::gaussian). 
#'   For Cox models, only a character string "cox" is allowed.
#' @param family_string A character string representing the selected family, 
#'   e.g., "gaussian".
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
#' @param spike A logical vector indicating which columns of \code{x} contain
#' a spike at zero. The length and order of \code{spike} must match those of
#' the columns in \code{x}.
#' @param acd_parameter Named list of ACD parameters produced by \code{fit_acd()}, 
#' with length equal to \code{ncol(x)}. Each list element corresponds to a variable; 
#' if an element is \code{NULL}, the variable was not specified in the 
#' \code{acdx} argument of \code{fit_mfp}.
#' @param spike_decision Named vector indicating how spike-at-zero (SAZ) 
#' variables are handled. Each element corresponds to a variable and encodes 
#' the selected strategy: `1` = include FP for positive values plus binary SAZ, 
#' `2` = treat as continuous FP only, `3` = include binary SAZ only.
#' @param prev_adj_params Named list storing adjustment variable transformations 
#' from previous steps for each variable.
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
#' \code{select_linear()}). This step differs from the mfp case because
#' linear models only use 1 df, while estimation of (every) fp power adds 
#' another df. This is also the case applied for categorical variables for 
#' which `df` are set to 1.
#' * the case that an acd transformation is requested (`acdx` is `TRUE` 
#' for `xi`) for the variable of interest (see \code{find_best_fpm_step()}).
#' * the (usual) case of the normal mfp algorithm to assess non-linear 
#' functional forms (see \code{find_best_fpm_step()}). 
#' 
#' Note that these cases do not encompass the setting that a variable is not
#' selected, because the evaluation is done for each variable in each cycle.
#' A variable which was de-selected in earlier cycles may be added to the 
#' working model again. Also see \code{find_best_fp_cycle()}.
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
#' It is a closed testing procedure and is implemented in \code{select_ra2()} and
#' extended for ACD transformation in \code{select_ra2_acd()} according to 
#' Royston and Sauerbrei (2016). 
#' 
#' For the other criteria `aic` and `bic` all FP models up to the desired degree
#' are fitted and the model with the lowest value for the information criteria 
#' is chosen as the final one. This is implemented in \code{select_ic()}.
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
                              family_string,
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
                              spike,
                              spike_decision,
                              acd_parameter,
                              prev_adj_params,
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
  
  fit1 <- select_fct(
    x = x, xi = xi, keep = keep, degree = degree, acdx = acdx, 
    y = y, family = family, family_string = family_string, 
    weights = weights, offset = offset, 
    powers_current = powers_current, powers = powers,  
    criterion = criterion, ftest = ftest, select = select, alpha = alpha,
    method = method, strata = strata, nocenter = nocenter, 
    control = control, rownames = rownames, zero = zero, catzero = catzero,
    spike = spike, spike_decision = spike_decision,
    acd_parameter = acd_parameter, prev_adj_params = prev_adj_params
  )
  
  if (verbose) {
    print_mfp_step(xi = xi, criterion = criterion, fit = fit1)
  }
    
  # prepare power to return. Remove trailing NAs, unless ACD is used
  power_best <- as.numeric(fit1$power_best)
  
  if (!fit1$acd) {
    #power_best <- na.omit(power_best) 
    power_best <- power_best[!is.na(power_best)]
    if (length(power_best) == 0) {
      power_best <- NA 
    }
  } 
  
  # create names for powers
  names(power_best) <- name_transformed_variables(
    xi,length(power_best),acd = acdx[xi]
    )
  
  # Return best power, spike decision and current adjustment parameters
  if (all(is.na(power_best)) || !spike[xi]) { 
    return(list(power_best = power_best, spike_decision = spike_decision,
                current_adj_params = fit1$current_adj_params))
  }
  
  # ----------------------------------------------------------------------------
  # Evaluate spike at zero (SAZ) variables to update spike_decision 
  # This is stage 2 of SAZ algorithm. Computed only when variable is selected
  # ----------------------------------------------------------------------------
  # fit candidate models to evaluate Spike at zero variables
  models <- fit_spike_zero_models(x = x, xi = xi, y = y, weights = weights, 
            offset = offset, family = family, method = method, strata = strata,
            nocenter = nocenter, control = control, rownames = rownames, 
            powers_current = powers_current, powers = powers, acdx = acdx, 
            zero = zero, catzero = catzero, power_best = power_best,
            spike = spike, spike_decision = spike_decision,
            acd_parameter = acd_parameter)
  
  # compute metrics for the three models to decide on SAZ models
  n_obs <- ifelse(family_string == "cox", sum(y[, 2]), nrow(x))
  metrics <- compute_metrics(fit1, models$fit2, models$fit3, n_obs,
            power_best)
  
  # Update spike_decision. That is decision on whether both components,  
  # FPm/linear only, or binary only is needed
  decision <- compute_decision(metrics, criterion, select, n_obs, ftest)
  spike_decision[xi] <- decision$decision
  
  # add spike metrics to fit1 for printing stage 2
  f_names <- rownames(fit1$metrics)[fit1$model_best]
  spike_metrics <- do.call(rbind, metrics)
  rownames(spike_metrics) <- c(f_names, strsplit(f_names, " \\+ ")[[1]])
  fit1$spike_metrics <- list(metrics = spike_metrics, 
                             spike_decision = spike_decision[xi],
                             pvalue = decision$pvalue)
  if (verbose) {
    print_mfp_step(xi = xi, criterion = criterion, fit = fit1, stage2 = TRUE)
  }
  
  return(list(power_best = power_best, spike_decision = spike_decision,
              current_adj_params = fit1$current_adj_params))

}

#' Fit Candidate Models for a single variable for Stage 2 of SAZ
#'
#' This function fits the two candidate regression models required for stage 2 
#' of the spike-at-zero (SAZ) algorithm. Stage 2 ideally considers three models 
#' (Model 1: complex, Model 2: FPm/linear only, Model 3: binary only), but 
#' Model 1 is already fitted in stage 1. This function therefore fits Model 2, 
#' which includes fractional polynomial/linear terms plus adjusted covariates, 
#' and Model 3, which includes binary-only terms plus adjusted covariates.
#'
#' @details
#' The function performs the following steps: First, it transforms the predictors
#' according to selected FP powers and other adjustments from stage 1. Second, 
#' it fits models 2 and 3 for comparison. Then returns both model fits (`fit2` 
#' and `fit3`) and the transformed predictor matrix
#' 
#' For the variable of interest `xi`:
#' * In **fit2**, the function sets `spike_decision[xi] = 2` and
#'   `catzero[xi] = FALSE` before calling `transform_matrix()`. This guarantees
#'   that the binary indicator of `xi` is excluded, so only its FP/ACD 
#'   transformation is used.  
#' * In **fit3**, the binary indicator for `xi` is taken directly from
#'   `catzero[[xi]]` and combined with adjustment variables (the other
#'   transformed predictors).  
#'
#' For variables other than `xi`, binary creation and spike-at-zero handling
#' follow the rules in `transform_matrix()`. In particular:
#' * spike_decision = 1 or 2 → FP/ACD transformation included; binary depends on `catzero`.
#' * spike_decision = 3 → FP/ACD removed, only binary retained.
#'
#' @return A list with the following elements:
#'   * `fit2`: Model 2 fit (FPm/linear only + adjusted covariates).
#'   * `fit3`: Model 3 fit (binary only + adjusted covariates).
#'   * `x_transformed`: List containing transformed predictor matrices.
#' @inheritParams find_best_fp_step
#' @keywords internal
fit_spike_zero_models <- function(x, xi, y, weights, offset, family, method, strata,
                                 nocenter, control, rownames, powers_current, powers,
                                 acdx, zero, catzero, power_best,
                                 spike, spike_decision, acd_parameter) {
  
  # Update the estimated FP power for variable xi
  powers_current[[xi]] <- as.vector(power_best)
  
  # Prepare centering vector: no variables are centered at this stage
  # should be centered after the algorithm converges in fit_mfp()
  center_vector <- setNames(rep(FALSE, ncol(x)), colnames(x))
  
  # set spike_decision for xi to 2 to get correct transformation
  spike_decision[xi] <- 2
  
  # transform_matrix() requires logical catzero. so we need to convert list to 
  # logical. Any variable with spike_decision = 2 will be assigned catzero = FALSE by
  # the function which is desired for our variable of interest. if spike_decision = 3,
  # binary indicator will only be returned and by default catzero will be TRUE
  # for this variable  
  catzero_logical <- sapply(catzero, function(x) !is.null(x))
  
  # prevent binary creation for xi inside transform_matrix()
  catzero_logical[xi] <- FALSE
  
  # Transform x using the updated powers
  x_transformed <- transform_matrix(
    x = x,  
    power_list = powers_current,
    center = center_vector,
    acdx = acdx,
    keep_x_order = TRUE, 
    check_binary = TRUE,
    acd_parameter_list = acd_parameter,
    zero = NULL, # already handled upstream in fit_mfp() 
    catzero = catzero_logical,
    spike = spike,
    spike_decision = spike_decision
    )
  
  #---------------------------------------------------
  # Fit2: continuous-only (xi binary never created here)
  #---------------------------------------------------
  fit2 <- fit_model(
    x = x_transformed$x_transformed,
    y = y, family = family, weights = weights, offset = offset, 
    method = method, strata = strata, nocenter = nocenter, 
    control = control, rownames = rownames
  )
  
  #---------------------------------------------------
  # Fit3: binary-only + adjustment variables
  #---------------------------------------------------
  
  # adjustment variables = all except xi
  adjustment_variable_names  <- setdiff(names(powers_current), xi)
  
  # Only bind adjustment variables if there are any
  adjustment_matrix  <- if (length(adjustment_variable_names) > 0) {
    do.call(cbind, x_transformed$x_trafo[adjustment_variable_names])
    } else {
      NULL
    }
  
  # Build binary column with explicit name
  binary_col <- as.matrix(catzero[[xi]])
  colnames(binary_col) <- paste0(xi, "_bin")
  
  # Combine binary + adjustments
  if (!is.null(adjustment_matrix)) {
    binary_plus_adjustment <- cbind(binary_col, adjustment_matrix)
  } else {
    binary_plus_adjustment <- binary_col
  }

  fit3 <- fit_model(
    x = binary_plus_adjustment,
    y = y, family = family, weights = weights, offset = offset, 
    method = method, strata = strata, nocenter = nocenter, 
    control = control, rownames = rownames
  )
  
  return(list(fit2 = fit2, fit3 = fit3, x_transformed = x_transformed))
}

#' Compute Model Metrics for Candidate Spike-at-zero Models
#'
#' This function computes fit statistics for the three candidate models 
#' used in the spike-at-zero (SAZ) algorithm. Model 1 metrics are extracted 
#' directly from the previously fitted model, while metrics for Model 2 
#' (FPm/linear only) and Model 3 (binary-only) are computed using 
#' `calculate_model_metrics`.
#'
#' @param fit1 Fitted object for Model 1 (complex model from stage 1 of SAZ).
#' @param fit2 Fitted object for Model 2 (FPm/linear only plus adjusted covariates).
#' @param fit3 Fitted object for Model 3 (binary-only plus adjusted covariates).
#' @param n_obs Number of observations in the dataset.
#' @param power_best Numeric vector of selected powers for the best FP terms
#' from stage 1 of SAZ algorithm.
#'
#' @details
#' The function determines the degree of the fractional polynomial based on `power_best`. 
#' Model 1 metrics are retrieved from the best-fit row of the `fit1` object. 
#' Metrics for Models 2 and 3 are calculated using `calculate_model_metrics`.
#'
#' @return A list with three elements:
#'   * `metrics1`: Fit statistics for Model 1.
#'   * `metrics2`: Fit statistics for Model 2.
#'   * `metrics3`: Fit statistics for Model 3.
#'
#' @keywords internal
compute_metrics <- function(fit1, fit2, fit3, n_obs, power_best) {
  
  # ACD can produce NA like c(NA,1) so degree will reduce to 1 and in this case
  # additional parameters = 0 since the power = 1, see mfpa paper
  power_best <- power_best[!is.na(power_best)]
  degree <- length(power_best[!is.na(power_best)])
  
  # Extract Model 1 metrics safely
  if (is.null(fit1$metrics) || is.null(fit1$model_best)) {
    stop("fit1 must contain 'metrics' and 'model_best' elements.")
  }
  
  if (fit1$model_best > nrow(fit1$metrics) || fit1$model_best < 1) {
    stop("fit1$model_best is out of bounds for fit1$metrics.")
  }
  
  metrics1 <- fit1$metrics[fit1$model_best, ]
  
  # Compute Model 2 metrics with degree adjustment
  deg2 <- if (degree == 1 && power_best == 1) 0 else degree
  metrics2 <- tryCatch(
    calculate_model_metrics(fit2, n_obs, deg2),
    error = function(e) stop("Failed to compute metrics for fit2: ", e$message)
  )
  
  # Compute Model 3 metrics
  metrics3 <- tryCatch(
    calculate_model_metrics(fit3, n_obs),
    error = function(e) stop("Failed to compute metrics for fit3: ", e$message)
  )
  return(list(metrics1 = metrics1, metrics2 = metrics2, metrics3 = metrics3))
}

#' Compute Model Selection Decision
#'
#' This function compares competing regression models based on a specified 
#' selection criterion (`pvalue`, `aic`, or `bic`) and returns a decision 
#' about which model to retain. It is used in stage 2 of the 
#' spike-at-zero (SAZ) algorithm, where a decision is made between alternative 
#' reduced or complex models after stage 1 selection.
#'
#' @param metrics A list containing model fit statistics for three candidate models:
#'   * `metrics1`: statistics for the most complex model selected in stage 1 of the SAZ algorithm, typically including both FPm/linear and binary terms.
#'   * `metrics2`: statistics for a reduced model, e.g., FPm/linear only.
#'   * `metrics3`: statistics for an alternative reduced model, e.g., binary only.
#' Each element must provide named values for:
#'   * `logl`: log-likelihood,
#'   * `df`: model degrees of freedom,
#'   * `aic`: Akaike Information Criterion,
#'   * `bic`: Bayesian Information Criterion,
#'   * `deviance_gaussian`: Gaussian deviance (required if `ftest = TRUE`),
#'   * `df_resid`: residual degrees of freedom (required if `ftest = TRUE`).
#' @param criterion A string specifying the selection criterion. 
#'   One of `"pvalue"`, `"aic"`, or `"bic"`.
#' @param select Numeric threshold for significance testing (used only when 
#'   `criterion = "pvalue"`).
#' @param n_obs Integer. Number of observations in the data (used for F-tests).
#' @param ftest Logical. If `TRUE`, use F-tests instead of likelihood ratio 
#'   tests when `criterion = "pvalue"`.
#'
#' @details
#' - When `criterion = "pvalue"`, the function compares:
#'   * Model 1 (both FPm/linear and binary terms),
#'   * Model 2 (FPm/linear only),
#'   * Model 3 (binary only).
#'   Significance is determined by comparing complex vs. reduced models 
#'   using either likelihood ratio tests or F-tests.
#'
#' - When `criterion = "aic"` or `"bic"`, the function selects the model 
#'   with the lowest information criterion value.
#'
#' @return A list with two elements:
#'   - `decision`: An integer indicating the selected model 
#'     (1 = both terms, 2 = FPm/linear only, 3 = binary only).
#'   - `pvalue`: A numeric vector of p-values from the pairwise tests, 
#'     or `NA` when `criterion` is `"aic"` or `"bic"`.
#'
#' @keywords internal
compute_decision <- function(metrics, criterion, select, n_obs, ftest = FALSE) {
  criterion <- tolower(criterion)
  
  if (criterion == "pvalue") {
    if (ftest) { 
      # Complex vs reduced model (FPm + Binary vs FPm)
      stats1  <- calculate_f_test(
        deviances = c(metrics$metrics2["deviance_gaussian"],
                      metrics$metrics1["deviance_gaussian"]),
        dfs_resid = c(metrics$metrics2["df_resid"],
                      metrics$metrics1["df_resid"]),
        n_obs = n_obs
      )
      # Complex vs reduced model (FPm + Binary vs Binary)
      stats2  <- calculate_f_test(
        deviances = c(metrics$metrics3["deviance_gaussian"],
                      metrics$metrics1["deviance_gaussian"]),
        dfs_resid = c(metrics$metrics3["df_resid"],
                      metrics$metrics1["df_resid"]),
        n_obs = n_obs
      )
      
      
    } else {
    
    # reduced model (FPm/Linear-mode 2) vs complex model (FPm/linear + binary-model 1)
    stats1 <- calculate_lr_test(
      c(metrics$metrics2["logl"], metrics$metrics1["logl"]),
      c(metrics$metrics2["df"], metrics$metrics1["df"])
    )
    
    # reduced model (binary only-mode 3) vs complex model (FPm/linear + binary-model 1)
    stats2 <- calculate_lr_test(
      c(metrics$metrics3["logl"], metrics$metrics1["logl"]),
      c(metrics$metrics3["df"], metrics$metrics1["df"])
    )
    }
    
    p12 <- stats1$pvalue
    p13 <- stats2$pvalue
    
    # Decision logic
    decision <- if (p12 <= select & p13 <= select) {
      1  # Model 1: both terms
    } else if (p12 <= select & p13 > select) {
      3  # Model 3: binary only
    } else if (p12 > select & p13 <= select) {
      2  # Model 2: FPm/linear only
    } else {
      # If both are not significant, use log-likelihood
      if (metrics$metrics2["logl"] > metrics$metrics3["logl"]) 2 else 3
    }
    
    # Return both decision and p-values
    return(list(
      decision = decision,
      pvalue = c(p12, p13)
    ))
    
  } else if (criterion == "aic") {
    decision <- which.min(c(metrics$metrics1["aic"], 
                            metrics$metrics2["aic"], 
                            metrics$metrics3["aic"]))
    return(list(decision = decision, pvalue = NA))
    
  } else { # criterion == "bic"
    decision <- which.min(c(metrics$metrics1["bic"], 
                            metrics$metrics2["bic"], 
                            metrics$metrics3["bic"]))
    return(list(decision = decision, pvalue = NA))
  }
}


#' Function to find the best FP functions of given degree for a single variable
#' 
#' Handles the FP1 and the higher order FP cases. For parameter definitions, see
#' \code{find_best_fp_step()}.
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
                               family_string,
                               zero,
                               catzero, # a list of binary variables or null
                               spike_decision, # a numeric vector
                               acd_parameter,
                               prev_adj_params,
                               ...) {
  # Number of events
  n_obs <- ifelse(family_string == "cox", sum(y[, 2]), nrow(x))
  
  if (degree == 1) {
    # remove linear model for normal data, but keep for acd transformation
    # powers = numeric(0) if powers = c(1,1)
    powers <- lapply(powers, function(v) setdiff(v, c(1)))
    
  }
    
  # generate FP data for x of interest (xi) and adjustment variables
  # Takes into account variables that should not be shifted thru 'zero'
  x_transformed <- transform_data_step(
    x = x, xi = xi, df = 2 * degree, powers_current = powers_current, 
    powers = powers, acdx = acdx, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter,
    prev_adj_params = prev_adj_params
  )
  
  # if (parallel) { # VERY SLOW. NO NEED FOR PARALLEL
  #   metrics_list <- foreach(i = seq_along(x_transformed$data_fp),
  #                          .packages = "mfp2",
  #                          .export = c("fit_model", "calculate_model_metrics")) %dopar% {
  #   # combine FP variables for x of interest with adjustment variables
  #   fit <- fit_model(
  #     x = cbind(x_transformed$data_fp[[i]], x_transformed$data_adj),
  #     y = y,
  #     family = family, ...
  #   )
  # 
  #   # use degree many additional degrees of freedom
  #   p <- sprintf("%g", x_transformed$powers_fp[i, , drop = TRUE])
  # 
  #   # respect acd
  #   if (acdx[xi]) {
  #     p[length(p)] <- sprintf("A(%s)", p[length(p)])
  #   }
  #   
  #   # compute metrics
  #   metric_val <- calculate_model_metrics(fit, n_obs, degree)
  # 
  #   # return named list element
  #   list(name = paste(p, collapse = " "), value = metric_val)
  # 
  # }
  # 
  # # Convert list to named data.table
  # metrics <- do.call(rbind, lapply(metrics_list, `[[`, "value"))
  # rownames(metrics) <- sapply(metrics_list, `[[`, "name")
  # } else {
  # Sequential approach
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
    if (acdx[xi]) {
      p[length(p)] <- sprintf("A(%s)", p[length(p)])
    }
    metrics[[paste(p, collapse = " ")]] <- calculate_model_metrics(
      fit, n_obs, degree # fit will include df of binary when catzero is used so this part does not change
    )
  }
  
  metrics <- do.call(rbind, metrics) 
  #}
  model_best <- as.numeric(which.max(metrics[, "logl"]))
  
  list(
    acd = acdx[xi],
    powers = x_transformed$powers_fp, 
    power_best = x_transformed$powers_fp[model_best, , drop = TRUE], 
    metrics = metrics, 
    model_best = model_best,
    zero = zero[xi],
    catzero = ifelse(!is.null(catzero[[xi]]), TRUE, FALSE),
    current_adj_params = x_transformed$current_params
  ) 
}

#' Function to fit a null model excluding variable of interest
#' 
#' "Null" model here refers to a model which does not include the variable 
#' of interest `xi`. 
#' For parameter definitions, see \code{find_best_fp_step()}. All parameters 
#' captured by `...` are passed on to \code{fit_model()}.
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
                          family_string,
                          zero,
                          catzero,
                          spike_decision,
                          acd_parameter,
                          prev_adj_params,
                          ...) {
  
  # Number of events
  n_obs <- ifelse(family_string == "cox", sum(y[, 2]), nrow(x))
  
  # transform all data as given by current working model
  # set variable of interest to linear term only
  x_transformed <- transform_data_step(
    x = x, xi = xi, df = 1, powers_current = powers_current, acdx = acdx, 
    powers = powers, zero = zero, catzero = catzero, 
    spike_decision = spike_decision, acd_parameter = acd_parameter,
    prev_adj_params = prev_adj_params
  ) 
  
  # An empty matrix can be returned which creates problem for cox model
  # convert it to NULL 
  x_tran <- x_transformed$data_adj
  if (is.null(x_tran) || ncol(x_tran) == 0L) {
    x_tran <- NULL
  }
  
  # fit null model
  # i.e. a model that does not contain xi but only adjustment variables
  # In addition, adjustment model can be NULL, so we have intercept only
  model_null <- fit_model(x = x_tran, y = y,
                          family = family, ...) 
  
  list(
    powers = NA,
    metrics = rbind(null = calculate_model_metrics(model_null, n_obs)),
    current_adj_params = x_transformed$current_params
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
#' For parameter definitions, see \code{find_best_fp_step}. 
#' All parameters captured by \code{...} are passed to \code{fit_model}.
#' 
#' @return 
#' A list with two entries: 
#' 
#' * `powers`: FP power(s) of `xi` (or its ACD transformation) in fitted model.
#' * `metrics`: A matrix with performance indices for fitted model.
#' * `prev_adj_params`: Previously adjusted parameters.
#' @inheritParams find_best_fp_step
#' @param ... Parameters passed to `fit_model()`.
fit_linear_step <- function(x, 
                            xi, 
                            y, 
                            powers_current,
                            powers,
                            acdx, 
                            family,
                            family_string,
                            zero,
                            catzero, 
                            spike_decision,
                            acd_parameter,
                            prev_adj_params,
                            ...) {
  # Number of events in survival models or observation in GLM 
  n_obs <- ifelse(family_string == "cox", sum(y[, 2]), nrow(x))
  
  # transform all data as given by current working model
  # set variable of interest to linear term only
  x_transformed <- transform_data_step(
    x = x, xi = xi, df = 1, powers_current = powers_current, acdx = acdx,
    powers = powers, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter,
    prev_adj_params = prev_adj_params
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
    metrics = metrics,
    current_adj_params = x_transformed$current_params
  )
}

#' Helper function to select between null and linear term for a single variable
#' 
#' To be used in \code{find_best_fp_step()}. Only used if `df = 1` for a variable.
#' Handles all criteria for selection.
#' For parameter explanations, see \code{find_best_fp_step()}. All parameters 
#' captured by `...` are passed to \code{fit_model()}.
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
#' * `spike_decision`: Spike decision flag for xi.
#' * `prev_adj_params`: Previously adjusted parameters.
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
                          family_string,
                          zero,
                          catzero,
                          spike,
                          spike_decision,
                          acd_parameter,
                          prev_adj_params,
                          ...) {
  
  n_obs <- ifelse(family_string == "cox", sum(y[, 2]), nrow(x))
  
  # Model 1: Null model
  fit_null <- fit_null_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx,
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision,
    acd_parameter = acd_parameter, prev_adj_params = prev_adj_params, ...
  )
  
  # Model 2: Linear model
  fit_linear <- fit_linear_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx, 
    family = family, family_string = family_string, zero = zero, catzero = catzero, spike_decision = spike_decision,
    acd_parameter = acd_parameter, prev_adj_params = prev_adj_params, ...
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
    catzero = ifelse(!is.null(catzero[[xi]]), TRUE, FALSE),
    spike = spike[xi],
    current_adj_params = fit_null$current_adj_params
  )
}

#' Function selection procedure based on closed testing procedure
#' 
#' Used in \code{find_best_fp_step()} when `criterion = "pvalue"`.
#' For parameter explanations, see \code{find_best_fp_step()}. All parameters 
#' captured by `...` are passed to \code{fit_model()}.
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
#' are computed by \code{find_best_fpm_step}.
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
#' * `spike_decision`: Spike decision flag for xi.
#' * `prev_adj_params`: Previously adjusted parameters.
#' @references 
#' Royston, P. and Sauerbrei, W., 2008. \emph{Multivariable Model - Building: 
#' A Pragmatic Approach to Regression Anaylsis based on Fractional Polynomials 
#' for Modelling Continuous Variables. John Wiley & Sons.}
#' 
#' @seealso 
#' \code{select_ra2_acd()}
#'  
#' @param degree integer > 0 giving the degree for the FP transformation. 
#' @param ... passed to fitting functions \code{fit_model()}. 
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
                       family_string,
                       zero,
                       catzero,
                       spike,
                       spike_decision,
                       acd_parameter,
                       prev_adj_params,
                       ...) {
  
  if (degree < 1) {
    return(NULL)
  }
  
  # Number of events in survival models or observations in GLM
  n_obs <- ifelse(family_string == "cox", sum(y[, 2]), nrow(x))
  
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

  #fpmax <- paste0("FP", degree)
  fpmax <- ifelse(spike[xi] || !is.null(catzero[[xi]]),  paste0("FP", degree, " + Binary"), paste0("FP", degree))
  
  # output list
  res <- list(
    keep = xi %in% keep, 
    acd = FALSE, 
    powers = NULL, 
    power_best = NULL, 
    metrics = NULL, 
    model_best = NULL, 
    statistic = NULL, 
    pvalue = NULL,
    spike = spike[xi],
    current_adj_params = NULL
  )
  
  # fit highest fp and null model for initial step
  fit_fpmax <- find_best_fpm_step(
    x = x, xi = xi, degree = degree, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx, 
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter,
    prev_adj_params = prev_adj_params, ...
  )
  fit_null <- fit_null_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx, 
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter,
    prev_adj_params = prev_adj_params, ...
  )
  res$metrics <- rbind(
    fit_fpmax$metrics[fit_fpmax$model_best, ],
    fit_null$metrics
  )
  rownames(res$metrics) <- c(fpmax, "null")
  res$powers <- rbind(fit_fpmax$power_best, NA)
  
  # return also the adjustment parameters
  res$current_adj_params <- fit_fpmax$current_adj_params
  
  # Test 1: test for overall significance (null vs best FPm)
  # df for tests are degree * 2
  stats <- calculate_test(res$metrics[c("null", fpmax), ], n_obs)
  res$statistic <- stats$statistic
  names(res$statistic) <- sprintf("%s vs null", fpmax)
  res$pvalue <- stats$pvalue
  names(res$pvalue) <- names(res$statistic)
  
  # Ensure keep is not NULL
  current_keep <- if (is.null(keep)) character(0) else keep
  
  if (stats$pvalue >= select && !(xi %in% current_keep)) {
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
    family = family,family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter,
    prev_adj_params = prev_adj_params, ...
  )
  
  old_names <- rownames(res$metrics)
  res$metrics <- rbind(
    res$metrics, 
    fit_lin$metrics
  )
  lin_names <- ifelse(spike[xi] || !is.null(catzero[[xi]]),"linear + Binary", "linear")
  rownames(res$metrics) <- c(old_names, lin_names)
  res$powers = rbind(res$powers, 
                     ensure_length(fit_lin$powers, ncol(res$powers)))
  
  stats <- calculate_test(res$metrics[c(lin_names, fpmax), ], n_obs)
  old_names <- names(res$statistic)
  res$statistic <- c(res$statistic, stats$statistic)
  names(res$statistic) <- c(old_names, sprintf("%s vs %s", fpmax, lin_names))
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
    # Vector of lower degrees
    lower_degrees <- 1:(degree - 1)
    
    # Construct FP names once
    fp_names <- paste0("FP", lower_degrees, if (spike[xi] || !is.null(catzero[[xi]])) " + Binary" else "")
    
    for (i in seq_along(lower_degrees)) {
      current_degree <- lower_degrees[i]
      fpm <- fp_names[i]
      
      fit_fpm <- find_best_fpm_step(
        x = x, xi = xi, degree = current_degree, y = y, 
        powers_current = powers_current, powers = powers, acdx = acdx,
        family = family, family_string = family_string, zero = zero, catzero = catzero,
        spike_decision = spike_decision, acd_parameter = acd_parameter,
        prev_adj_params = prev_adj_params, ...
      )
      # Append metrics
      old_names = rownames(res$metrics)
      res$metrics <- rbind(
        res$metrics, 
        fit_fpm$metrics[fit_fpm$model_best, ]
      )
      rownames(res$metrics) <- c(old_names, fpm)
      
      # Append powers
      res$powers <- rbind(res$powers, 
                          ensure_length(fit_fpm$power_best, ncol(res$powers)))
      
      # Calculate test statistics
      stats <- calculate_test(res$metrics[c(fpm, fpmax), ], n_obs)
      old_names <- names(res$statistic)
      res$statistic <- c(res$statistic, stats$statistic)
      names(res$statistic) <- c(old_names, sprintf("%s vs %s", fpmax, fpm))
      
      # Append p-values
      res$pvalue <- c(res$pvalue, stats$pvalue)
      names(res$pvalue) <- names(res$statistic)
      
      # Stop early if non-linearity detected
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
#' Used in \code{find_best_fp_step()} when `criterion = "pvalue"` and an 
#' ACD transformation is requested for `xi`.
#' For parameter explanations, see \code{find_best_fp_step()}. All parameters 
#' captured by `...` are passed on to \code{fit_model()}.
#' 
#' @details  
#' This function extends the algorithm used in \code{select_ra2()} to allow the 
#' usage of ACD transformations. The implementation follows the description 
#' in Royston and Sauerbrei (2016). The procedure is outlined in detail in 
#' the corresponding section in the documentation of \code{mfp2()}.
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
#' * `spike_decision`: Spike decision flag for xi.
#' * `prev_adj_params`: Previously adjusted parameters.
#' 
#' @references 
#' Royston, P. and Sauerbrei, W., 2016. \emph{mfpa: Extension of mfp using the
#' ACD covariate transformation for enhanced parametric multivariable modeling. 
#' The Stata Journal, 16(1), pp.72-87.}
#' 
#' @seealso 
#' \code{select_ra2()}
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
                           family_string,
                           zero,
                           catzero,
                           spike,
                           spike_decision,
                           acd_parameter,
                           prev_adj_params,
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
  n_obs <- ifelse(family_string == "cox", sum(y[, 2]), nrow(x))
  
  #fpmax <- "FP1(x, A(x))"
  fpmax <- ifelse(spike[xi] || !is.null(catzero[[xi]]),  "FP1(x, A(x)) +  Binary", "FP1(x, A(x))")
  
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
    pvalue = NULL,
    spike = spike[xi],
    current_adj_params = NULL
  )
  
  # fit highest fp and null model for initial step
  fit_fpmax <- find_best_fpm_step(
    x = x, xi = xi, degree = 2, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx, 
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter,
    prev_adj_params = prev_adj_params, ...
  )
  fit_null <- fit_null_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx,
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter,
    prev_adj_params = prev_adj_params, ...
  )
  res$metrics <- rbind(
    fit_fpmax$metrics[fit_fpmax$model_best, ],
    fit_null$metrics
  )
  rownames(res$metrics) <- c(fpmax, "null")
  res$powers <- rbind(fit_fpmax$power_best, fit_null$powers)
  
  res$current_adj_params <- fit_fpmax$current_adj_params
  
  # test for overall significance
  # df for tests are degree * 2 = 4
  stats <- calculate_test(res$metrics[c("null", fpmax), ], n_obs)
  res$statistic <- stats$statistic
  names(res$statistic) <- sprintf("%s vs null", fpmax)
  res$pvalue <- stats$pvalue
  names(res$pvalue) <- names(res$statistic)
  
  # Ensure keep is not NULL
  current_keep <- if (is.null(keep)) character(0) else keep
  
  if (stats$pvalue >= select && !(xi %in% current_keep)) {
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
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter,
    prev_adj_params = prev_adj_params, ...
  )
  
  old_names = rownames(res$metrics)
  res$metrics <- rbind(
    res$metrics, 
    fit_lin$metrics
  )
  
  lin_names <- ifelse(spike[xi] || !is.null(catzero[[xi]]),"linear + Binary", "linear")
  
  rownames(res$metrics) <- c(old_names, lin_names)
  res$powers <- rbind(res$powers, c(1, NA))
  
  stats <- calculate_test(res$metrics[c(lin_names, fpmax), ], n_obs)
  old_names <- names(res$statistic)
  res$statistic <- c(res$statistic, stats$statistic)
  names(res$statistic) <- c(old_names, sprintf("%s vs %s", fpmax, lin_names))
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
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter,
    prev_adj_params = prev_adj_params, ...
  )
  
  old_names = rownames(res$metrics)
  res$metrics <- rbind(
    res$metrics, 
    fit$metrics[fit$model_best, ]
  )
  FP_names <- ifelse(spike[xi] || !is.null(catzero[[xi]]),"FP1(x, .) + Binary", "FP1(x, .)")
  rownames(res$metrics) <- c(old_names, FP_names)
  res$powers <- rbind(res$powers, c(fit$power_best, NA))
  
  stats <- calculate_test(res$metrics[c(FP_names, fpmax), ], n_obs)
  old_names <- names(res$statistic)
  res$statistic <- c(res$statistic, stats$statistic)
  names(res$statistic) <- c(old_names, sprintf("%s vs %s", fpmax, FP_names))
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
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter,
    prev_adj_params = prev_adj_params, ...
  )
  
  old_names = rownames(res$metrics)
  res$metrics <- rbind(
    res$metrics, 
    fit_fp1a$metrics[fit_fp1a$model_best, ]
  )
  
  FP1_names <- ifelse(spike[xi] || !is.null(catzero[[xi]]),"FP1(., A(x)) + Binary", "FP1(., A(x))")
  
  rownames(res$metrics) <- c(old_names, FP1_names)
  res$powers <- rbind(res$powers, fit_fp1a$power_best)
  
  stats <- calculate_test(res$metrics[c(FP1_names, fpmax), ], n_obs)
  old_names <- names(res$statistic)
  res$statistic <- c(res$statistic, stats$statistic)
  names(res$statistic) <- c(old_names, sprintf("%s vs %s", fpmax, FP1_names))
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
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter,
    prev_adj_params = prev_adj_params, ...
  )
  
  old_names = rownames(res$metrics)
  res$metrics <- rbind(
    res$metrics, 
    fit_lineara$metrics
  )
  
  linx_names <- ifelse(spike[xi] || !is.null(catzero[[xi]]),"linear(., A(x)) + Binary", "linear(., A(x))")
  
  rownames(res$metrics) <- c(old_names, linx_names)
  res$powers <- rbind(res$powers, fit_lineara$powers)
  
  stats <- calculate_test(res$metrics[c(linx_names, FP1_names), ], n_obs)
  old_names <- names(res$statistic)
  res$statistic <- c(res$statistic, stats$statistic)
  names(res$statistic) <- c(old_names, sprintf("%s vs %s", FP1_names, linx_names))
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
#' Used in \code{find_best_fp_step()} when `criterion = "aic"` or `"bic"`.
#' For parameter explanations, see \code{find_best_fp_step()}. All parameters 
#' captured by `...` are passed on to \code{fit_model()}.
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
#' These best FPx models are computed in \code{find_best_fpm_step()}.
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
#' * `spike_decision`: Spike decision flag for xi.
#' * `prev_adj_params`: Previously adjusted parameters.
#' 
#' @seealso 
#' \code{select_ra2()}
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
                      family_string,
                      zero,
                      catzero,
                      spike,
                      spike_decision,
                      acd_parameter,
                      prev_adj_params,
                      ...) {
  
  if (degree < 1) {
    return(NULL)
  }
  
  #fpmax <- paste0("FP", degree)
  fpmax <- ifelse(spike[xi] || !is.null(catzero[[xi]]),  paste0("FP", degree, " + Binary"), paste0("FP", degree))
  
  # output list
  res <- list(
    keep = xi %in% keep,
    acd = FALSE, 
    powers = NULL, 
    power_best = NULL, 
    metrics = NULL, 
    model_best = NULL, 
    statistic = NA, 
    pvalue = NA,
    spike = spike[xi],
    current_adj_params = NULL
  )
  
  # fit all relevant models
  # Null Model
  fit_null <- fit_null_step(
      x = x, xi = xi, y = y, 
      powers_current = powers_current, powers = powers, acdx = acdx,
      family = family,family_string = family_string, zero = zero, catzero = catzero,
      spike_decision = spike_decision, acd_parameter = acd_parameter,
      prev_adj_params = prev_adj_params, ...
    )
  # Linear model
  fit_lin <- fit_linear_step(
      x = x, xi = xi, y = y, 
      powers_current = powers_current, powers = powers, acdx = acdx,
      family = family, family_string = family_string, zero = zero, catzero = catzero, 
      spike_decision = spike_decision, acd_parameter = acd_parameter,
      prev_adj_params = prev_adj_params, ...
    )
  
  res$current_adj_params <- fit_null$current_adj_params
  
  
  # All FPm models
  fits_fpm <- list()
  
  for (m in seq_len(degree)) {
    fits_fpm[[m]] <- find_best_fpm_step(
      x = x, xi = xi, degree = m, y = y, 
      powers_current = powers_current, powers = powers, acdx = acdx,
      family = family, family_string = family_string, zero = zero, catzero = catzero,
      spike_decision = spike_decision,acd_parameter = acd_parameter,
      prev_adj_params = prev_adj_params, ...
    )
  }

  if (spike[xi] || !is.null(catzero[[xi]])) {
    names(fits_fpm) <- sprintf("FP%g + Binary", seq_len(degree))
  } else {
    names(fits_fpm) <- sprintf("FP%g", seq_len(degree))
  }
  
  # output summary - only output best fpm models
  len_max <- ncol(fits_fpm[[fpmax]]$powers)
  
  res$powers <- lapply(fits_fpm, function(x) {
    ensure_length(x$powers[x$model_best, , drop = FALSE], len_max)
  })
  # Save the original FP names
  fp_names <- names(res$powers)
  
  res$powers <- do.call(
    rbind, 
    c(list(ensure_length(fit_null$powers, len_max), 
           ensure_length(fit_lin$powers, len_max)), 
      res$powers)
  )
  
  # Assign row names: null, linear (with + Binary if spike), then original FP names
  rownames(res$powers) <- c(
    "null",
    if (spike[xi] || !is.null(catzero[[xi]])) "linear + Binary" else "linear",
    fp_names
  )

  res$metrics <- lapply(fits_fpm, function(x) {
    x$metrics[x$model_best, , drop = FALSE] 
  })
  res$metrics <- do.call(
    rbind, 
    c(list(null = fit_null$metrics, linear = fit_lin$metrics), res$metrics)
  )
  rownames(res$metrics) <- c("null", ifelse(spike[xi] || !is.null(catzero[[xi]]),"linear + Binary", "linear"), 
                             names(fits_fpm))
  
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
                          family_string,
                          zero,
                          catzero,
                          spike,
                          spike_decision,
                          acd_parameter,
                          prev_adj_params,
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
    pvalue = NA,
    spike = spike[xi],
    current_adj_params = NULL
  )
  
  # fit all relevant models
  fit_null <- fit_null_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx,
    family = family, family_string = family_string, zero = zero, catzero = catzero, 
    spike_decision = spike_decision, acd_parameter = acd_parameter,
    prev_adj_params = prev_adj_params, ...
  )
  fit_lin <- fit_linear_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx_reset_xi,
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter,
    prev_adj_params = prev_adj_params, ...
  )
  fit_lina <- fit_linear_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx,
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter,
    prev_adj_params = prev_adj_params, ...
  )
  
  res$current_adj_params <- fit_null$current_adj_params
  
  suffix <- if (spike[xi] || !is.null(catzero[[xi]])) " + Binary" else ""
  fits <- setNames(
    list(
     find_best_fpm_step(
      x = x, xi = xi, degree = 1, y = y, 
      powers_current = powers_current, powers = powers, acdx = acdx_reset_xi,
      family = family, family_string = family_string, zero = zero, catzero = catzero,
      spike_decision = spike_decision,acd_parameter = acd_parameter,
      prev_adj_params = prev_adj_params, ...
    ), 
     find_best_fpm_step(
      x = x, xi = xi, degree = 1, y = y, 
      powers_current = powers_current, powers = powers, acdx = acdx,
      family = family, family_string = family_string, zero = zero, catzero = catzero, 
      spike_decision = spike_decision, acd_parameter = acd_parameter,
      prev_adj_params = prev_adj_params, ...
    ), 
     find_best_fpm_step(
      x = x, xi = xi, degree = 2, y = y, 
      powers_current = powers_current, powers = powers, acdx = acdx,
      family = family, family_string = family_string, zero = zero, catzero = catzero,
      spike_decision = spike_decision,acd_parameter = acd_parameter,
      prev_adj_params = prev_adj_params, ...
    )
  ),
  c(
  paste0("FP1(x, .)", suffix),
  paste0("FP1(., A(x))", suffix), 
  paste0("FP1(x, A(x))", suffix)
  )
  )

  # output summary - only output best fpm models
  len_max <- 2
  
  res$powers <- lapply(fits, function(x) {
    ensure_length(x$powers[x$model_best, , drop = FALSE], len_max)
  })
  
  # DO WE NEED TO CHANGE THE NAMES DUE TO SPIKE?
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
  rownames(res$metrics) <- c(
    "null", 
    paste0("linear", suffix), 
    paste0("linear(., A(x))", suffix), 
    names(fits)
    )
  
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
#' This function prepares transformed data for a focal predictor `xi` and its
#' adjustment variables. Adjustment variables are transformed using either
#' fractional polynomials or acd transformations, depending on their assigned
#' powers and parameters. Spike-at-zero effects can be incorporated using
#' binary indicators. Previously computed adjustment variables can be reused
#' if their parameters have not changed.
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
#' @param acd_parameter Named list of ACD parameters produced by \code{fit_acd()}, 
#' with length equal to \code{ncol(x)}.
#' @param spike_decision a named numeric vector with the same names as the 
#' `x` matrix. Each element controls how the corresponding adjustment 
#' variable contributes to the adjustment matrix:
#' 
#' * `1`: combine both the transformed adjustment variable and its binary 
#'   indicator (`catzero`).
#' * `2`: include only the transformed adjustment variable.
#' * `3`: include only the binary indicator (`catzero`).
#' This applies only to adjustment variables, not the focal predictor `xi`.
#' @param prev_adj_params Named list containing results from a previous call to 
#' this function. Used to avoid recomputing adjustment variables if both powers 
#' and `spike_decision` remain unchanged for a certain variable.
#' @details
#' This function extracts the adjustment variables and applies the corresponding
#' FP or ACD transformations based on `powers_current`. When evaluating the variable
#' of interest `xi`, it is necessary to account for other variables in the model, 
#' which may be transformed or untransformed depending on their individual powers. 
#' Some powers may be `NA`, indicating that the corresponding variable has been 
#' excluded from the adjustment set. 
#' 
#' To improve efficiency, the function avoids recomputing adjustment variables if 
#' their parameters (`powers` and `spike_decision`) have not changed from a 
#' previous step, as provided in `prev_adj_params`.
#' 
#' The role of `spike_decision` is to determine how each adjustment variable is
#' represented in the presence of potential spike-at-zero effects. For every 
#' adjustment variable, the function can include the transformed variable, the 
#' binary indicator for nonpositive values (`catzero`), or both. This makes it 
#' possible to model spike-at-zero behavior directly in the adjustment matrix, 
#' while preserving flexibility in how variables are included.
#'  
#' The function also returns the FP data for predictor `xi` of interest, which 
#' depends on the specified degrees of freedom. For example,
#' `df = 2` is equivalent to FP degree one, resulting in the generation of 8 
#' variables. If `acdx` for the current variables of interest is set to `TRUE`, 
#' however, 64 variables are generated.
#' 
#' When `df = 1`, this function returns data unchanged, i.e. a "linear" 
#' transformation with power equal to 1. In case `acdx[xi] = TRUE`, the 
#' acd transformation is applied. 
#' 
#' @return 
#' A list containing `power_best` (numeric vector of best FP powers for `xi`, 
#' possibly including `NA`), `spike_decision` (updated spike-at-zero strategy), 
#' and `current_adj_params` (adjustment variable transformations used in this step).
#' \item{powers_fp}{Numeric vector of FP powers used for `xi`.}
#' \item{data_fp}{List of transformed data for `xi`.}
#' \item{powers_adj}{Named list of FP powers used for adjustment variables.}
#' \item{data_adj}{Matrix of transformed adjustment variables, or `NULL` if none.}
#' \item{current_params}{Named list containing adjustment variable parameters 
#'  for reuse in later steps, including \code{powers_adj}, 
#'   \code{spike_decision_adj}, \code{data_adj_list}, and \code{data_adj}.}
#' @keywords internal
transform_data_step <- function(x,
                                xi,
                                powers_current,
                                df, 
                                powers,
                                acdx,
                                zero,
                                catzero,
                                spike_decision,
                                acd_parameter,
                                prev_adj_params
                                ) { 
  # sort x based on the names of the powers
  names_powers_current <- names(powers_current)
  x <- x[, names_powers_current, drop = FALSE]
  
  # adjustment variables = all except xi
  vars_adj <- setdiff(names_powers_current, xi)
  x_adj <- if (length(vars_adj) > 0) x[, vars_adj, drop = FALSE] else NULL
  
  # Initialize powers and spike decision for adjustment variables in case they
  # do not exist
  powers_adj <- NULL
  spike_decision_adj <- NULL
  # Prepare a container for adjusted data
  data_adj_list <- list()
  
  if (length(vars_adj) > 0) {
    # transformations of adjustment variables
    powers_adj <- powers_current[vars_adj]
    acdx_adj <- acdx[vars_adj]
    zero_adj <- zero[vars_adj]
    acd_parameter_adj <- acd_parameter[vars_adj]
    spike_decision_adj <- spike_decision[vars_adj]

    # For each adjustment variable, check if recomputation is needed. This is
    # only done if FP powers and spike decision change
    for (varname in vars_adj) {
      xvec <- x_adj[, varname, drop = TRUE]
      power_current <- as.numeric(powers_adj[[varname]])
      spike_current <- as.numeric(spike_decision[varname])
      # Extract previous adjustments, powers, and spike_decision
      prev_for_var <- NULL
      if (!is.null(prev_adj_params[[xi]])) {
        prev_xi <- prev_adj_params[[xi]]
        
        prev_data_adj <- prev_xi$data_adj_list[[varname]]
        prev_power <- as.numeric(prev_xi$powers_adj[[varname]])
        prev_spike <- as.numeric(prev_xi$spike_decision_adj[[varname]])
        
        if (!is.null(prev_data_adj)) {
          prev_for_var <- list(
            data_adj = prev_data_adj,
            powers = prev_power,
            spike_decision = prev_spike
          )
        }
      }
      
      # Check if powers or spike_decision changed
      # Determine if recomputation is needed
      recompute <- TRUE
      if (!is.null(prev_for_var)) {
        powers_same <- identical(prev_for_var$powers, power_current)
        spike_same <- identical(prev_for_var$spike_decision, spike_current)
        recompute <- !(powers_same && spike_same)
      }
      
      # Transform if needed, otherwise reuse previous data
      # # If power is NA or no transformation needed, create empty matrix
      if (all(is.na(power_current))) {
      adj_mat <- matrix(nrow = nrow(x), ncol = 0)
      } else if (!recompute) {
        adj_mat <- prev_for_var$data_adj
        if (is.vector(adj_mat)) adj_mat <- matrix(adj_mat, ncol = 1)
      } else {
      # Transform variable if recomputation is needed
      if (acdx_adj[varname]) {
        transformed <- transform_vector_acd(
          x = xvec, 
          power = power_current,
          zero = zero_adj[varname],
          acd_parameter = acd_parameter_adj[[varname]],
          powers = powers[[varname]] # powers needed by acd() function
          
        )$acd
      } else {
        transformed <- transform_vector_fp(
          x = xvec,
          power = power_current, 
          zero = zero_adj[varname]
        )
      }
        # normalize to matrix
        if (is.null(transformed)) {
          adj_mat <- matrix(nrow = nrow(x), ncol = 0)
        } else if (is.null(dim(transformed))) {
          adj_mat <- matrix(transformed, ncol = 1)
        } else {
          adj_mat <- as.matrix(transformed)
        }
      # catzero for this variable (may be NULL) if catzero = FALSE
      # which is not possibe if spike = TRUE
        # catzero (may be NULL)
        cz <- catzero[[varname]]
        if (!is.null(cz)) {
          cz_mat <- as.matrix(cz)
        } else {
          cz_mat <- NULL
        }
      
      # apply spike_decision rule
      spike_val <- spike_decision[[varname]]
      if (spike_val == 1L) {
        adj_mat  <- cbind(cz_mat, transformed)
      } else if (spike_val == 2L) {
        adj_mat  <- transformed
      } else if (spike_val == 3L) {
        adj_mat  <- cz_mat
      }
      }
      
      # assign informative column names that preserve mapping to the variable
      # assign names only if nonempty
    
      if (!is.null(adj_mat) && ncol(adj_mat) > 0) {
        colnames(adj_mat) <- paste0(varname, "_adj", seq_len(ncol(adj_mat)))
      }
      
      # collect results per variable keyed by variable name
      data_adj_list[[varname]] <- adj_mat

    }
    }
 
  # combine all adjustment columns into a single matrix (or NULL)
  data_adj <- if (length(data_adj_list) > 0) do.call(cbind, data_adj_list) else NULL
  
  # generate fp data for xi (unchanged, spike_decision not applied here)
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
                                          catzero = catzero[[xi]],
                                          acd_parameter = acd_parameter[[xi]])
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
  
  #Store everything under xi in current_params
  current_params <- list()
  current_params[[xi]] <- list(
    powers_adj = powers_adj,
    spike_decision_adj = spike_decision_adj,
    data_adj_list = data_adj_list,
    data_adj = data_adj
    )
  
  # Return results and current parameters for next step
  list(
    powers_fp = powers_fp, 
    data_fp = data_fp, 
    powers_adj = powers_adj,
    data_adj = data_adj,
    current_params = current_params
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
