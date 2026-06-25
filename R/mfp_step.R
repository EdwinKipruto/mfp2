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
#' @param catzero A named list of structural-zero indicators. Each element is
#' either `NULL` or an n x 1 integer/numeric matrix. Non-NULL elements indicate
#' variables for which a binary structural-zero column is available.
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
#' @param force_max_fp A logical vector of length \code{nvars}, named by
#'  variable name. If \code{TRUE} for variable \code{xi}, forces
#' \code{select_ic()} and \code{select_ic_acd()} to select the most complex
#' functional form available — the highest-likelihood FP model at the degree
#'   specified by \code{df} for non-ACD variables, or \code{FP1(x, A(x))} for
#'   ACD variables — without competing it against simpler forms (null, linear,
#'   or lower-degree FP) under AIC/BIC. The best power combination within the
#'   selected form is still determined by \code{find_best_fpm_step()} via
#'   deviance minimisation, which is equivalent to AIC/BIC minimisation at
#'   fixed degrees of freedom. Has no effect when \code{criterion = "pvalue"}
#'   since \code{select_ra2()} is used in that case and \code{alpha = 1}
#'   already guarantees acceptance of the most complex form. 
#' @param has_offset logical indicating whether an offset was specified before
#' missing offsets were replaced by zeros internally.
#' @param precomputed_adj Optional internal adjustment object returned by
#' \code{build_adjustment_step()}. When supplied, adjustment-variable
#' transformations are reused and only current-variable candidate
#' transformations are generated.
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
#' @keywords internal
#' @noRd
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
                              force_max_fp,
                              has_offset,
                              n_obs,
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
    weights = weights, offset = offset, force_max_fp = force_max_fp,
    powers_current = powers_current, powers = powers,  
    criterion = criterion, ftest = ftest, select = select, alpha = alpha,
    method = method, strata = strata, nocenter = nocenter, n_obs = n_obs,
    control = control, rownames = rownames, zero = zero, catzero = catzero,
    spike = spike, spike_decision = spike_decision, has_offset = has_offset,
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
  
  # Stage 1 selected null.
  # Nothing more to do.
  # For spike variables, clear stale binary-only decision.
  if (all(is.na(power_best))) {
    # Stage 2 is conditional on stage 1 selecting the full
    # continuous + binary spike model. If stage 1 selects null,
    # any previous binary-only decision is stale.
    if (isTRUE(spike[[xi]])) {
      # Here 2L is just the neutral non-binary-only state and its appropiate
      # since the variable will be eliminated
      spike_decision[[xi]] <- 2L
    }
    
    return(list(
      power_best = power_best,
      spike_decision = spike_decision,
      current_adj_params = fit1$current_adj_params
    ))
  }
  
  if (!isTRUE(spike[[xi]])) {
    # Stage 1 selected a non-null model,
    # but xi is not a spike variable.
    # Nothing more to do.
    return(list(
      power_best = power_best,
      spike_decision = spike_decision,
      current_adj_params = fit1$current_adj_params
    ))
  }
  
  # If we get here:
  # power_best is not NA
  # AND xi is a spike variable
  # Therefore run SAZ stage 2.
  # ----------------------------------------------------------------------------
  # Evaluate spike at zero (SAZ) variables to update spike_decision 
  # This is stage 2 of SAZ algorithm. Computed only when variable is selected
  # ----------------------------------------------------------------------------
  # fit candidate models to evaluate Spike at zero variables
  models <- fit_saz_reduced_models(
    stage1_selection = fit1,
    xi = xi,
    y = y,
    weights = weights,
    offset = offset,
    family = family,
    family_string = family_string,
    method = method,
    strata = strata,
    nocenter = nocenter,
    control = control,
    rownames = rownames,
    has_offset = has_offset
  )
  
  # compute metrics for the three models to decide on SAZ models
  #n_obs <- ifelse(family_string == "cox", sum(y[, 2]), nrow(x))
  
  metrics <- compute_saz_stage2_metrics(
    fit1 = models$fit1,
    fit2 = models$fit2,
    fit3 = models$fit3,
    n_obs = n_obs,
    power_best = power_best
  )
  
  # Update spike_decision. That is decision on whether both components,  
  # FPm/linear only, or binary only is needed
  decision <- compute_saz_stage2_decision(metrics, criterion, select, n_obs, ftest)
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
#' * `catzero`: Logical in
#' dicating whether a combination of a zero transformation 
#'   and a binary indicator variable was applied to \code{xi}. This means that 
#'   nonpositive values of \code{xi} were set to zero, only positive values were 
#'   transformed, and an additional binary variable was created to indicate 
#'   whether \code{xi} was positive or nonpositive.
#' @inheritParams find_best_fp_step
#' @param degree degrees of freedom for fp transformation of `xi`.
#' @param ... parameters passed to `fit_model()`.
#' @keywords internal
#' @noRd
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
                               spike,
                               spike_decision, # a numeric vector
                               acd_parameter,
                               prev_adj_params,
                               has_offset,
                               precomputed_adj = NULL,
                               n_obs,
                               ...
                               ) {
  # Number of events
  #n_obs <- ifelse(family_string == "cox", sum(y[, 2]), nrow(x))
  
  if (degree == 1) {
    # Degree-1 candidate generation is for the current variable xi only.
    # Exclude power 1 for xi because the corresponding linear model is fitted
    # separately by fit_linear_step().
    #
    # Do not remove power 1 from the full powers list: transform_data_step()
    # also receives powers and may pass powers[[v]] to ACD adjustment-variable
    # transformations. Adjustment variables must keep their original allowed
    # power sets.
    #
    # In the ACD case, generate_powers_acd(degree = 1) produces candidates
    # of the form c(NA, p). Thus p = 1 corresponds to the ACD-linear candidate
    # already fitted separately and should not be duplicated for xi.
    powers[[xi]] <- setdiff(powers[[xi]], 1)
  }
  
  # generate FP data for x of interest (xi) and adjustment variables
  # Takes into account variables that should not be shifted thru 'zero'
  x_transformed <- transform_data_step(
    x = x, xi = xi, df = 2 * degree, powers_current = powers_current, 
    powers = powers, acdx = acdx, zero = zero, catzero = catzero,
    spike = spike, spike_decision = spike_decision, 
    acd_parameter = acd_parameter,
    prev_adj_params = prev_adj_params,
    precomputed_adj = precomputed_adj
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

  data_adj <- x_transformed$data_adj
  data_fp <- x_transformed$data_fp
  has_adj <- !is.null(data_adj) && NCOL(data_adj) > 0L
  use_glm_intercept_template <- !identical(family_string, "cox")
  
  first_xi <- data_fp[[1L]]
  n_xi_cols <- NCOL(first_xi)
  design_mat <- NULL
  x_has_intercept <- FALSE
  
  if (use_glm_intercept_template) {
    # For GLM candidate fits, build the full model matrix including the
    # intercept once. fit_glm() is told that the intercept is already present,
    # so it does not cbind another intercept for every candidate.
    intercept_col <- matrix(
      rep.int(1, nrow(first_xi)),
      ncol = 1L,
      dimnames = list(NULL, "(Intercept)")
    )
    design_mat <- if (has_adj) {
      cbind(intercept_col, first_xi, data_adj)
    } else {
      cbind(intercept_col, first_xi)
    }
    xi_cols <- seq_len(n_xi_cols) + 1L
    x_has_intercept <- TRUE
  } else if (has_adj) {
    # Cox models must not include an intercept. Copy the adjustment block once
    # and overwrite only the current-variable block per candidate.
    design_mat <- cbind(first_xi, data_adj)
    xi_cols <- seq_len(n_xi_cols)
  }
  metrics <- vector("list", length(data_fp))
  for (i in seq_along(data_fp)) {
    # combine FP variables for x of interest with adjustment variables.
    # The adjustment block is identical for all candidates, so copy it into the
    # working design matrix once and overwrite only the xi block per candidate.
    # For GLMs, the reusable matrix also includes the intercept column.
    data_xi <- data_fp[[i]]
    if (!is.null(design_mat) && NCOL(data_xi) == n_xi_cols) {
      design_mat[, xi_cols] <- data_xi
      fit <- fit_model(
        x               = design_mat, # catzero plays a role here thru data_fp and data_adj
        y               = y, 
        family          = family,
        family_string   = family_string,
        has_offset      = has_offset,
        x_has_intercept = x_has_intercept,
        ...
      )
    } else {
      if (has_adj) {
        # Defensive fallback: candidate generation should produce a fixed
        # number of xi columns within one degree, but preserve behaviour if it
        # ever does not.
        x_fit <- cbind(data_xi, data_adj)
      } else {
        x_fit <- data_xi
      }
      if (use_glm_intercept_template) {
        x_fit <- cbind("(Intercept)" = rep.int(1, nrow(data_xi)), x_fit)
      }
      fit <- fit_model(
        x               = x_fit,
        y               = y, 
        family          = family,
        family_string   = family_string,
        has_offset      = has_offset,
        x_has_intercept = use_glm_intercept_template,
        ...
      )
    }
    
    metrics[[i]] <- calculate_model_metrics(
      fit, # fit will include df of binary when catzero is used so this part does not change
      n_obs,
      degree
    )
  
  }
  
  metrics <- do.call(rbind, metrics) 
  
  model_best <- as.numeric(which.max(metrics[, "logl"]))
  x_transformed$current_params[[xi]]$data_xi <- x_transformed$data_fp[[model_best]]
  
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
#' @keywords internal
#' @noRd
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
                          spike,
                          spike_decision,
                          acd_parameter,
                          prev_adj_params,
                          has_offset,
                          precomputed_adj = NULL,
                          n_obs,
                          ...
                          ) {
  
  # Number of events
  #n_obs <- ifelse(family_string == "cox", sum(y[, 2]), nrow(x))
  
  # Null model uses adjustment variables only. Do not call
  # transform_data_step(df = 1) here: that would generate the current-variable
  # linear candidate and then discard it.
  adj <- if (is.null(precomputed_adj)) {
    build_adjustment_step(
      x = x,
      xi = xi,
      powers_current = powers_current,
      powers = powers,
      acdx = acdx,
      zero = zero,
      catzero = catzero,
      spike = spike,
      spike_decision = spike_decision,
      acd_parameter = acd_parameter,
      prev_adj_params = prev_adj_params
    )
  } else {
    precomputed_adj
  }
  
  current_params <- list()
  current_params[[xi]] <- list(
    powers_adj = adj$powers_adj,
    spike_decision_adj = adj$spike_decision_adj,
    data_adj_list = adj$data_adj_list,
    data_adj = adj$data_adj
  )
  
  # An empty matrix can be returned which creates problem for cox model
  # convert it to NULL 
  x_tran <- adj$data_adj
  if (is.null(x_tran) || ncol(x_tran) == 0L) {
    x_tran <- NULL
  }
  
  # fit null model
  # i.e. a model that does not contain xi but only adjustment variables
  # In addition, adjustment model can be NULL, so we have intercept only
  model_null <- fit_model(x = x_tran, 
                          y = y,
                          family = family,
                          family_string   = family_string,
                          has_offset = has_offset,
                          ...
                          ) 
  
  list(
    powers = NA,
    metrics = rbind(null = calculate_model_metrics(model_null, n_obs)),
    current_adj_params = current_params
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
#' @keywords internal
#' @noRd
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
                            spike,
                            spike_decision,
                            acd_parameter,
                            prev_adj_params,
                            has_offset,
                            n_obs,
                            ...,
                            precomputed_adj = NULL) {
  # Number of events in survival models or observation in GLM 
  #n_obs <- ifelse(family_string == "cox", sum(y[, 2]), nrow(x))
  
  # transform all data as given by current working model
  # set variable of interest to linear term only
  x_transformed <- transform_data_step(
    x = x, xi = xi, df = 1, powers_current = powers_current, acdx = acdx,
    powers = powers, zero = zero, catzero = catzero, spike = spike,
    spike_decision = spike_decision, acd_parameter = acd_parameter,
    prev_adj_params = prev_adj_params,
    precomputed_adj = precomputed_adj
  ) 
  x_transformed$current_params[[xi]]$data_xi <- x_transformed$data_fp[[1]]
  
  # fit a model based on the assumption that xi is linear. If catzero, its
  # corresponding binary variable will be included in the model
  model_linear <- fit_model(
    x = cbind(x_transformed$data_fp[[1]], x_transformed$data_adj),
    y = y,
    family = family, 
    family_string = family_string, 
    has_offset = has_offset,
    ...
  )
  
  # respect acd
  metrics <- rbind(
    linear = calculate_model_metrics(
      obj = model_linear,
      n_obs = n_obs,
      df_additional = 0
    )
  )
  
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
#' @param force_max_fp not used
#' @param ... passed to fitting functions. 
#' @inheritParams find_best_fp_step 
#' @keywords internal
#' @noRd
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
                          force_max_fp,
                          has_offset,
                          n_obs,
                          ...) {
  
  #n_obs <- ifelse(family_string == "cox", sum(y[, 2]), nrow(x))
  
  precomputed_adj <- build_adjustment_step(
    x = x,
    xi = xi,
    powers_current = powers_current,
    powers = powers,
    acdx = acdx,
    zero = zero,
    catzero = catzero,
    spike = spike,
    spike_decision = spike_decision,
    acd_parameter = acd_parameter,
    prev_adj_params = prev_adj_params
  )
  
  # Model 1: Null model
  fit_null <- fit_null_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx,
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, has_offset = has_offset, spike = spike,
    acd_parameter = acd_parameter, prev_adj_params = prev_adj_params,
    precomputed_adj = precomputed_adj, n_obs = n_obs, ...
  )
  
  # Model 2: Linear model
  fit_linear <- fit_linear_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx, 
    family = family, family_string = family_string, zero = zero, 
    catzero = catzero, spike_decision = spike_decision, spike = spike,
    has_offset = has_offset, n_obs = n_obs,
    acd_parameter = acd_parameter, prev_adj_params = prev_adj_params,
    precomputed_adj = precomputed_adj, ...
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
    current_adj_params = if (model_best == 1) fit_null$current_adj_params else fit_linear$current_adj_params
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
#' @keywords internal
#' @noRd
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
                       force_max_fp,
                       has_offset,
                       n_obs,
                       ...) {
  
  if (degree < 1) {
    return(NULL)
  }
  
  
  # Number of events in survival models or observations in GLM
  #n_obs <- ifelse(family_string == "cox", sum(y[, 2]), nrow(x))
  
  # simplify testing by defining test helper function
  if (ftest) {
    calculate_test <- function(metrics, n_obs) {
      calculate_f_test(
        deviances = metrics[, "deviance_gaussian", drop = TRUE],
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
  
  # Model labels get a binary suffix when the focal variable contributes a
  # structural-zero indicator, either through spike handling or catzero.
  # Compute this once to keep row names, test labels, and candidate names
  # consistent throughout the selection function.
  has_binary_xi <- isTRUE(spike[[xi]]) || !is.null(catzero[[xi]])
  binary_suffix <- if (has_binary_xi) " + Binary" else ""
  
  #fpmax <- ifelse(spike[xi] || !is.null(catzero[[xi]]),  paste0("FP", degree, " + Binary"), paste0("FP", degree))
  fpmax <- paste0("FP", degree, binary_suffix)
  
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
  
  
  precomputed_adj <- build_adjustment_step(
    x = x,
    xi = xi,
    powers_current = powers_current,
    powers = powers,
    acdx = acdx,
    zero = zero,
    catzero = catzero,
    spike = spike,
    spike_decision = spike_decision,
    acd_parameter = acd_parameter,
    prev_adj_params = prev_adj_params
  )
  # fit highest fp and null model for initial step
  fit_fpmax <- find_best_fpm_step(
    x = x, xi = xi, degree = degree, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx, 
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter, spike = spike,
    prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
    precomputed_adj = precomputed_adj, ...
  )
  fit_null <- fit_null_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx, 
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter, spike = spike,
    prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
    precomputed_adj = precomputed_adj, ...
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
    res$current_adj_params <- fit_null$current_adj_params
    return(res)
  }
  
  # Test 2: test for non-linearity (linear vs best FPm)
  # df for tests are degree * 2 - 1
  fit_lin <- fit_linear_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx, 
    family = family,family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter, spike = spike,
    prev_adj_params = prev_adj_params, has_offset = has_offset,n_obs = n_obs,
    precomputed_adj = precomputed_adj, ...
  )
  
  old_names <- rownames(res$metrics)
  res$metrics <- rbind(
    res$metrics, 
    fit_lin$metrics
  )
  
  lin_names <- paste0("linear", binary_suffix)
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
    res$current_adj_params <- fit_lin$current_adj_params
    return(res)
  }
  
  # Test 3: test for complexity of the functions (best FP1 vs best FPm)
  # do this for all fps with lower degrees
  # dfs for tests are decreasing
  if (degree > 1) {
    # Vector of lower degrees
    lower_degrees <- 1:(degree - 1)
    
    # Construct FP names once
    fp_names <- paste0("FP", lower_degrees, binary_suffix)
    
    for (i in seq_along(lower_degrees)) {
      current_degree <- lower_degrees[i]
      fpm <- fp_names[i]
      
      fit_fpm <- find_best_fpm_step(
        x = x, xi = xi, degree = current_degree, y = y, 
        powers_current = powers_current, powers = powers, acdx = acdx,
        family = family, family_string = family_string, zero = zero, catzero = catzero,
        spike_decision = spike_decision, acd_parameter = acd_parameter, spike = spike,
        prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
        precomputed_adj = precomputed_adj, ...
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
        res$current_adj_params <- fit_fpm$current_adj_params
        return(res)
      }
    }
    
  }
  
  # return highest power
  res$power_best <- fit_fpmax$powers[fit_fpmax$model_best, , drop = FALSE]
  res$model_best <- 1
  res$current_adj_params <- fit_fpmax$current_adj_params
  
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
#' @keywords internal
#' @noRd
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
                           force_max_fp,
                           has_offset,
                           n_obs,
                           ...) {
  
  # simplify testing by defining test helper function
  if (ftest) {
    calculate_test <- function(metrics, n_obs) {
      calculate_f_test(
        deviances = metrics[, "deviance_gaussian", drop = TRUE],
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
  
  # Model labels get a binary suffix when the focal variable contributes a
  # structural-zero indicator, either through spike handling or catzero.
  # Compute this once to keep row names, test labels, and candidate names
  # consistent throughout the selection function.
  has_binary_xi <- isTRUE(spike[[xi]]) || !is.null(catzero[[xi]])
  binary_suffix <- if (has_binary_xi) " + Binary" else ""
  fpmax <- paste0("FP1(x, A(x))", binary_suffix)

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
  
  
  precomputed_adj <- build_adjustment_step(
    x = x,
    xi = xi,
    powers_current = powers_current,
    powers = powers,
    acdx = acdx,
    zero = zero,
    catzero = catzero,
    spike = spike,
    spike_decision = spike_decision,
    acd_parameter = acd_parameter,
    prev_adj_params = prev_adj_params
  )
  # fit highest fp and null model for initial step
  fit_fpmax <- find_best_fpm_step(
    x = x, xi = xi, degree = 2, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx, 
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter, spike = spike,
    prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
    precomputed_adj = precomputed_adj, ...
  )
  fit_null <- fit_null_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx,
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter, spike = spike,
    prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
    precomputed_adj = precomputed_adj, ...
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
    res$current_adj_params <- fit_null$current_adj_params
    return(res)
  }
  
  # test for non-linearity in x
  # df for tests are degree * 2 - 1 = 3
  fit_lin <- fit_linear_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx_reset_xi,
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter, spike = spike,
    prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
    precomputed_adj = precomputed_adj, ...
  )
  
  old_names = rownames(res$metrics)
  res$metrics <- rbind(
    res$metrics, 
    fit_lin$metrics
  )
  
  lin_names <- paste0("linear", binary_suffix)
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
    res$current_adj_params <- fit_lin$current_adj_params
    return(res)
  }
  
  # test for functional form, comparison with FP1(x, .)
  fit <- find_best_fpm_step(
    x = x, xi = xi, degree = 1, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx_reset_xi,
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter, spike = spike,
    prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
    precomputed_adj = precomputed_adj, ...
  )
  
  old_names = rownames(res$metrics)
  res$metrics <- rbind(
    res$metrics, 
    fit$metrics[fit$model_best, ]
  )
  
  FP_names <- paste0("FP1(x, .)", binary_suffix)
  
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
    res$current_adj_params <- fit$current_adj_params
    return(res)
  }
  
  # test for functional form, comparison with FP1(., A(x))
  fit_fp1a <- find_best_fpm_step(
    x = x, xi = xi, degree = 1, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx, 
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter, spike = spike,
    prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
    precomputed_adj = precomputed_adj, ...
  )
  
  old_names = rownames(res$metrics)
  res$metrics <- rbind(
    res$metrics, 
    fit_fp1a$metrics[fit_fp1a$model_best, ]
  )
  
  FP1_names <- paste0("FP1(., A(x))", binary_suffix)
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
    res$current_adj_params <- fit_fpmax$current_adj_params
    return(res)
  }
  
  # return best model between FP1(., A(x)) and linear(., A(x))
  fit_lineara <- fit_linear_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx,
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter,spike = spike,
    prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
    precomputed_adj = precomputed_adj, ...
  )
  
  old_names = rownames(res$metrics)
  res$metrics <- rbind(
    res$metrics, 
    fit_lineara$metrics
  )
  
  linx_names <- paste0("linear(., A(x))", binary_suffix)
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
    res$current_adj_params <- fit_fp1a$current_adj_params
    return(res)
  }
  
  # use linear(., A(x))
  res$power_best = matrix(c(NA, 1), ncol = 2)
  res$model_best = 6
  res$current_adj_params <- fit_lineara$current_adj_params
  
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
#' @keywords internal
#' @noRd
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
                      force_max_fp,
                      has_offset,
                      n_obs,
                      ...) {
  
  if (degree < 1) {
    return(NULL)
  }
  
  # Model labels get a binary suffix when the focal variable contributes a
  # structural-zero indicator, either through spike handling or catzero.
  # Compute this once to keep row names, test labels, and candidate names
  # consistent throughout the selection function.
  has_binary_xi <- isTRUE(spike[[xi]]) || !is.null(catzero[[xi]])
  binary_suffix <- if (has_binary_xi) " + Binary" else ""
  fpmax <- paste0("FP", degree, binary_suffix)
  
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
  
  
  precomputed_adj <- build_adjustment_step(
    x = x,
    xi = xi,
    powers_current = powers_current,
    powers = powers,
    acdx = acdx,
    zero = zero,
    catzero = catzero,
    spike = spike,
    spike_decision = spike_decision,
    acd_parameter = acd_parameter,
    prev_adj_params = prev_adj_params
  )
  # Fit all relevant models ----------------------------------------------------
  # Null Model
  fit_null <- fit_null_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx,
    family = family,family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter, spike = spike,
    prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
    precomputed_adj = precomputed_adj, ...
  )
  # Linear model
  fit_lin <- fit_linear_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx,
    family = family, family_string = family_string, zero = zero, catzero = catzero, 
    spike_decision = spike_decision, acd_parameter = acd_parameter, spike = spike,
    prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
    precomputed_adj = precomputed_adj, ...
  )
  
  res$current_adj_params <- fit_null$current_adj_params
  
  
  # All FPm models
  fits_fpm <- list()
  
  for (m in seq_len(degree)) {
    fits_fpm[[m]] <- find_best_fpm_step(
      x = x, xi = xi, degree = m, y = y, 
      powers_current = powers_current, powers = powers, acdx = acdx,
      family = family, family_string = family_string, zero = zero, catzero = catzero,
      spike_decision = spike_decision,acd_parameter = acd_parameter,spike = spike,
      prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
      precomputed_adj = precomputed_adj, ...
    )
  }
  
  names(fits_fpm) <- paste0("FP", seq_len(degree), binary_suffix)
  
  # Assemble output summary ----------------------------------------------------
  # output summary - only output best fpm models
  len_max <- ncol(fits_fpm[[fpmax]]$powers)
  
  candidate_names <- c(
    "null",
    paste0("linear", binary_suffix),
    names(fits_fpm)
  )
  
  res$powers <- lapply(fits_fpm, function(x) {
    ensure_length(x$powers[x$model_best, , drop = FALSE], len_max)
  })
  
  res$powers <- do.call(
    rbind,
    c(
      list(
        ensure_length(fit_null$powers, len_max),
        ensure_length(fit_lin$powers, len_max)
      ),
      res$powers
    )
  )
  rownames(res$powers) <- candidate_names
  
  res$metrics <- lapply(fits_fpm, function(x) {
    x$metrics[x$model_best, , drop = FALSE]
  })
  
  res$metrics <- do.call(
    rbind,
    c(
      list(
        fit_null$metrics,
        fit_lin$metrics
      ),
      res$metrics
    )
  )
  rownames(res$metrics) <- candidate_names
  
  # Select best model ----------------------------------------------------------
  if (isTRUE(force_max_fp[xi])) {
    # Skip functional form competition: always select the most complex FP model
    # at the requested degree (last row of res$metrics = FPm or FPm + Binary).
    # The best power combination within that degree was already found by
    # find_best_fpm_step() via deviance minimisation, which is equivalent to
    # AIC/BIC minimisation at fixed df (same penalty for all power combinations).
    res$model_best <- nrow(res$metrics)
  } else if (xi %in% keep) { 
    # Prevent selection of null model; choose best among linear through FPm.
    ind_select <- 2:nrow(res$metrics)
    res$model_best <- which.min(
      res$metrics[ind_select, tolower(criterion), drop = TRUE]
    )
    # shift by 1 since null row was excluded
    res$model_best <- res$model_best + 1
  }else{
    # Unrestricted: choose best among null through FPm.
    ind_select <- 1:nrow(res$metrics)
    res$model_best <- which.min(
      res$metrics[ind_select, tolower(criterion), drop = TRUE]
    )
  }
  
  res$power_best <- res$powers[res$model_best, , drop = FALSE]
  res$current_adj_params <- if (res$model_best == 1L) {
    fit_null$current_adj_params
  } else if (res$model_best == 2L) {
    fit_lin$current_adj_params
  } else {
    fits_fpm[[res$model_best - 2L]]$current_adj_params
  }
  
  res
}

#' @describeIn select_ic Function to select ACD based transformation.
#' @keywords internal
#' @noRd
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
                          force_max_fp,
                          has_offset,
                          n_obs,
                          ...) {
  
  acdx_reset_xi <- acdx
  acdx_reset_xi[xi] <- FALSE
  # Model labels get a binary suffix when the focal variable contributes a
  # structural-zero indicator, either through spike handling or catzero.
  # Compute this once to keep row names, test labels, and candidate names
  # consistent throughout the selection function.
  has_binary_xi <- isTRUE(spike[[xi]]) || !is.null(catzero[[xi]])
  binary_suffix <- if (has_binary_xi) " + Binary" else ""
  
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
  
  
  precomputed_adj <- build_adjustment_step(
    x = x,
    xi = xi,
    powers_current = powers_current,
    powers = powers,
    acdx = acdx,
    zero = zero,
    catzero = catzero,
    spike = spike,
    spike_decision = spike_decision,
    acd_parameter = acd_parameter,
    prev_adj_params = prev_adj_params
  )
  # fit all relevant models
  fit_null <- fit_null_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx,
    family = family, family_string = family_string, zero = zero, catzero = catzero, 
    spike_decision = spike_decision, acd_parameter = acd_parameter, spike = spike,
    prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
    precomputed_adj = precomputed_adj, ...
  )
  fit_lin <- fit_linear_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx_reset_xi,
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter, spike = spike,
    prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs, 
    precomputed_adj = precomputed_adj, ...
  )
  fit_lina <- fit_linear_step(
    x = x, xi = xi, y = y, 
    powers_current = powers_current, powers = powers, acdx = acdx,
    family = family, family_string = family_string, zero = zero, catzero = catzero,
    spike_decision = spike_decision, acd_parameter = acd_parameter,spike = spike,
    prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
    precomputed_adj = precomputed_adj, ...
  )
  
  res$current_adj_params <- fit_null$current_adj_params
  

  fits <- setNames(
    list(
      find_best_fpm_step(
        x = x, xi = xi, degree = 1, y = y, 
        powers_current = powers_current, powers = powers, acdx = acdx_reset_xi,
        family = family, family_string = family_string, zero = zero, catzero = catzero,
        spike_decision = spike_decision,acd_parameter = acd_parameter,spike = spike,
        prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
        precomputed_adj = precomputed_adj, ...
      ), 
      find_best_fpm_step(
        x = x, xi = xi, degree = 1, y = y, 
        powers_current = powers_current, powers = powers, acdx = acdx,
        family = family, family_string = family_string, zero = zero, catzero = catzero, 
        spike_decision = spike_decision, acd_parameter = acd_parameter,spike = spike,
        prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
        precomputed_adj = precomputed_adj, ...
      ), 
      find_best_fpm_step(
        x = x, xi = xi, degree = 2, y = y, 
        powers_current = powers_current, powers = powers, acdx = acdx,
        family = family, family_string = family_string, zero = zero, catzero = catzero,
        spike_decision = spike_decision,acd_parameter = acd_parameter,spike = spike,
        prev_adj_params = prev_adj_params, has_offset = has_offset, n_obs = n_obs,
        precomputed_adj = precomputed_adj, ...
      )
    ),
    c(
      paste0("FP1(x, .)", binary_suffix),
      paste0("FP1(., A(x))", binary_suffix), 
      paste0("FP1(x, A(x))", binary_suffix)
    )
  )
  
  # Assemble output summary ----------------------------------------------------
  # output summary - only output best fpm models
  len_max <- 2L
  
  candidate_names <- c(
    "null",
    paste0("linear", binary_suffix),
    paste0("linear(., A(x))", binary_suffix),
    names(fits)
  )
  
  res$powers <- lapply(fits, function(x) {
    ensure_length(x$powers[x$model_best, , drop = FALSE], len_max)
  })
  
  res$powers <- do.call(
    rbind,
    c(
      list(
        ensure_length(fit_null$powers, len_max),
        ensure_length(fit_lin$powers, len_max),
        ensure_length(fit_lina$powers, len_max)
      ),
      res$powers
    )
  )
  rownames(res$powers) <- candidate_names
  
  res$metrics <- lapply(fits, function(x) {
    x$metrics[x$model_best, , drop = FALSE]
  })
  
  res$metrics <- do.call(
    rbind,
    c(
      list(
        fit_null$metrics,
        fit_lin$metrics,
        fit_lina$metrics
      ),
      res$metrics
    )
  )
  rownames(res$metrics) <- candidate_names
  # Select best model ---------------------------------------------------------- 
  if (isTRUE(force_max_fp[xi])) {
    # Skip functional form competition: always select the most complex ACD
    # model — FP1(x, A(x)) — which is the last row of res$metrics.
    # The best power combination within it was already found by
    # find_best_fpm_step() via deviance minimisation.
    res$model_best <- nrow(res$metrics)
  } else if (xi %in% keep) {
    # Prevent selection of null model; choose best among linear through FP1(x, A(x)).
    ind_select <- 2:nrow(res$metrics)
    res$model_best <- which.min(
      res$metrics[ind_select, tolower(criterion), drop = TRUE]
    )
    res$model_best <- res$model_best + 1L
  }else{
    # Unrestricted: choose best among null through FP1(x, A(x))
    ind_select = 1:nrow(res$metrics)
    res$model_best <- which.min(
      res$metrics[ind_select, tolower(criterion), drop = TRUE]
    )
  }
  
  res$power_best = res$powers[res$model_best, , drop = FALSE]
  res$current_adj_params <- if (res$model_best == 1L) {
    fit_null$current_adj_params
  } else if (res$model_best == 2L) {
    fit_lin$current_adj_params
  } else if (res$model_best == 3L) {
    fit_lina$current_adj_params
  } else {
    fits[[res$model_best - 3L]]$current_adj_params
  }
  
  res
}

#' Build transformed adjustment data for one MFP step
#'
#' Internal helper used by transform_data_step() and selection functions to
#' avoid recomputing the adjustment matrix repeatedly for the same focal
#' variable, powers_current, and spike_decision state.
#' @inheritParams transform_data_step
#' @keywords internal
#' @noRd
build_adjustment_step <- function(x,
                                  xi,
                                  powers_current,
                                  powers,
                                  acdx,
                                  zero,
                                  catzero,
                                  spike,
                                  spike_decision,
                                  acd_parameter,
                                  prev_adj_params) {
  # ---------------------------------------------------------------------------
  # Purpose
  # ---------------------------------------------------------------------------
  # Build the adjustment matrix for the current focal variable `xi`.
  #
  # In the MFP cycle, each variable is updated conditional on all other variables.
  # For focal variable `xi`, the adjustment variables are all variables except
  # `xi`. Their transformed columns are collected into `data_adj`.
  #
  # This function also maintains a small per-variable cache:
  # if an adjustment variable has the same FP powers and the same spike decision
  # as in the previous call for this focal variable, its transformed matrix is
  # reused instead of recomputed.
  #
  # Important internal invariant:
  #   catzero[[varname]] is either NULL or an n x 1 numeric/integer matrix.
  #
  # It is NOT a vector inside this inner fitting code.
  
  # Use powers_current names for the canonical variable order.
  # Accessing x columns by name avoids copying/reordering the full x matrix.
  names_powers_current <- names(powers_current)
  
  # Adjustment variables are all variables except the current focal variable.
  vars_adj <- setdiff(names_powers_current, xi)
  
  # These are returned even when there are no adjustment variables.
  powers_adj <- NULL
  spike_decision_adj <- NULL
  
  # Per-variable transformed adjustment matrices.
  # Each element is keyed by the adjustment variable name.
  data_adj_list <- list()
  
  if (length(vars_adj) > 0L) {
    # Reusable empty matrix for variables that contribute no adjustment columns,
    # for example eliminated variables with all-NA powers.
    empty_adj_mat <- matrix(nrow = nrow(x), ncol = 0L)
    
    # Previous cached adjustment information for this focal variable.
    # prev_adj_params is indexed by xi.
    prev_xi <- prev_adj_params[[xi]]
    has_prev <- !is.null(prev_xi)
    
    # Slice all adjustment-variable metadata once.
    # This avoids repeated indexing into the full objects inside the loop.
    powers_adj <- powers_current[vars_adj]
    
    # A variable is considered eliminated when all stored powers are NA.
    # Eliminated variables contribute zero adjustment columns unless they are
    # binary-only spike variables, which are handled before this branch.
    eliminated <- vapply(
      powers_adj,
      function(p) all(is.na(as.numeric(p))),
      logical(1L)
    )
    
    acdx_adj <- acdx[vars_adj]
    zero_adj <- zero[vars_adj]
    acd_parameter_adj <- acd_parameter[vars_adj]
    spike_adj <- spike[vars_adj]
    spike_decision_adj <- spike_decision[vars_adj]
    
    # Convert spike decisions once to integer.
    # This is important because identical(3, 3L) is FALSE in R.
    # Cache comparisons below use identical(), so current and previous spike
    # decisions must have the same type.
    spike_decision_int_adj <- as.integer(spike_decision_adj)
    names(spike_decision_int_adj) <- vars_adj
    
    # -------------------------------------------------------------------------
    # Precompute binary-only spike flags
    # -------------------------------------------------------------------------
    # spike_decision == 3 means the continuous FP part is ignored and the
    # variable is represented only by the structural-zero indicator.
    #
    # This branch must be checked before the eliminated/all-NA-powers branch,
    # because for binary-only spike variables the FP powers are intentionally
    # irrelevant.
    spike_binary_only_flags <- vapply(
      vars_adj,
      function(v) {
        isTRUE(spike_adj[[v]]) && identical(spike_decision_int_adj[[v]], 3L)
      },
      logical(1L)
    )
    names(spike_binary_only_flags) <- vars_adj
    
    # -------------------------------------------------------------------------
    # Precompute normalized power keys for cache comparison
    # -------------------------------------------------------------------------
    # For ordinary variables, the cache key is the current FP powers.
    #
    # For binary-only spike variables, the FP powers are irrelevant. The helper
    # normalize_powers_for_convergence() converts their power key to NA_real_
    # whenever spike_decision == 3. That prevents unnecessary recomputation if
    # only the ignored continuous powers changed.
    current_power_keys_adj <- normalize_powers_for_convergence(
      powers_adj,
      spike_decision_int_adj
    )
    
    # Previous normalized power keys are needed only if a previous cache exists.
    prev_power_keys_adj <- NULL
    
    if (has_prev) {
      # Keep previous spike decisions as integer for consistent comparison.
      prev_spike_decision_int_adj <- as.integer(
        prev_xi$spike_decision_adj[vars_adj]
      )
      names(prev_spike_decision_int_adj) <- vars_adj
      
      prev_power_keys_adj <- normalize_powers_for_convergence(
        prev_xi$powers_adj[vars_adj],
        prev_spike_decision_int_adj
      )
    }
    
    # -------------------------------------------------------------------------
    # Build adjustment columns variable by variable
    # -------------------------------------------------------------------------
    for (varname in vars_adj) {
      
      # Current spike decision for this adjustment variable.
      # Needed for cache comparison and catzero/spike handling.
      spike_current <- spike_decision_int_adj[[varname]]
      
      # -----------------------------------------------------------------------
      # Extract previous cache entry for this adjustment variable
      # -----------------------------------------------------------------------
      # We only need two per-variable previous values here:
      #   1. prev_data_adj: the matrix to reuse on a cache hit
      #   2. prev_spike_decision: the previous spike decision for comparison
      #
      # Previous powers are compared through the precomputed prev_power_keys_adj,
      # so no per-variable previous-power object is needed here.
      prev_data_adj <- NULL
      prev_spike_decision <- NULL
      
      if (has_prev) {
        prev_data_adj <- prev_xi$data_adj_list[[varname]]
        
        # Must be integer to match spike_current.
        prev_spike_decision <- prev_spike_decision_int_adj[[varname]]
      }
      
      # -----------------------------------------------------------------------
      # Cache invalidation
      # -----------------------------------------------------------------------
      # Recompute only if either:
      #   1. the normalized power key changed, or
      #   2. the spike decision changed.
      #
      # Other transformation-relevant inputs such as acdx, zero, acd_parameter,
      # spike, and catzero are fixed after preprocessing within one fit_mfp()
      # run, so they are not part of this cache key.
      recompute <- TRUE
      
      if (!is.null(prev_data_adj)) {
        prev_power_key <- prev_power_keys_adj[[varname]]
        current_power_key <- current_power_keys_adj[[varname]]
        
        powers_same <- identical(prev_power_key, current_power_key)
        spike_decision_same <- identical(
          prev_spike_decision,
          spike_current
        )
        
        recompute <- !(powers_same && spike_decision_same)
      }
      
      # -----------------------------------------------------------------------
      # Build or reuse adjustment matrix for this variable
      # -----------------------------------------------------------------------
      
      if (spike_binary_only_flags[[varname]]) {
        # Binary-only spike:
        # use only the structural-zero indicator matrix.
        #
        # Under the internal invariant, catzero[[varname]] is already an n x 1
        # matrix, so do not wrap it in matrix() or as.matrix().
        cz <- catzero[[varname]]
        
        if (is.null(cz)) {
          stop(
            "Internal error: binary-only spike variable '", varname,
            "' has no catzero indicator.",
            call. = FALSE
          )
        }
        
        adj_mat <- cz
        
      } else if (eliminated[[varname]]) {
        # Eliminated non-binary-only variable:
        # contributes no adjustment columns.
        adj_mat <- empty_adj_mat
        
      } else if (!recompute) {
        # Cache hit:
        # reuse the previous transformed matrix for this variable.
        adj_mat <- prev_data_adj
        
      } else {
        xvec <- x[, varname, drop = TRUE]
        power_current <- as.numeric(powers_adj[[varname]])
        
        # Cache miss:
        # recompute the transformed continuous/FP part.
        if (isTRUE(acdx_adj[[varname]])) {
          transformed <- transform_vector_acd(
            x = xvec,
            power = power_current,
            zero = zero_adj[[varname]],
            acd_parameter = acd_parameter_adj[[varname]],
            
            # Full candidate powers are needed by acd() when acd_parameter = NULL.
            powers = powers[[varname]]
          )$acd
        } else {
          transformed <- transform_vector_fp(
            x = xvec,
            power = power_current,
            zero = zero_adj[[varname]]
          )
        }
        
        # Normalize transformed output to matrix form.
        # Downstream code expects every adjustment contribution to be a matrix,
        # including one-column transformations.
        if (is.null(transformed)) {
          adj_mat <- empty_adj_mat
        } else if (is.null(dim(transformed))) {
          adj_mat <- matrix(transformed, ncol = 1L)
        } else {
          adj_mat <- as.matrix(transformed)
        }
        
        # ---------------------------------------------------------------------
        # Add structural-zero indicator if needed
        # ---------------------------------------------------------------------
        # catzero handling depends on spike_decision:
        #
        #   spike_decision == 1:
        #     include both catzero and continuous FP part.
        #
        #   spike_decision == 2:
        #     include continuous FP part only.
        #
        #   spike_decision == 3:
        #     include catzero only.
        #
        # Non-spike variables with catzero available get catzero prepended to
        # their transformed FP columns.
        cz_mat <- catzero[[varname]]
        
        if (!is.null(cz_mat)) {
          if (isTRUE(spike_adj[[varname]])) {
            spike_val <- spike_current
            
            if (identical(spike_val, 1L)) {
              adj_mat <- cbind(cz_mat, adj_mat)
            } else if (identical(spike_val, 2L)) {
              adj_mat <- adj_mat
            } else if (identical(spike_val, 3L)) {
              adj_mat <- cz_mat
            }
          } else {
            adj_mat <- cbind(cz_mat, adj_mat)
          }
        }
      }
      
      # Assign informative column names.
      #
      # The final adjustment matrix is constructed by cbind-ing all entries of
      # data_adj_list. Prefixing with varname preserves the mapping from columns
      # back to the original adjustment variable.
      if (!is.null(adj_mat) && ncol(adj_mat) > 0L) {
        colnames(adj_mat) <- paste0(varname, "_adj", seq_len(ncol(adj_mat)))
      }
      
      # Store per-variable matrix for cache reuse in later iterations.
      data_adj_list[[varname]] <- adj_mat
    }
  }
  
  # Combine all adjustment-variable matrices into one adjustment matrix.
  #
  # If there are no adjustment variables, return NULL. If some variables are
  # eliminated, their entries are 0-column matrices and do not add columns.
  data_adj <- if (length(data_adj_list) > 0L) {
    do.call(cbind, data_adj_list)
  } else {
    NULL
  }
  
  list(
    powers_adj = powers_adj,
    spike_decision_adj = spike_decision_adj,
    data_adj_list = data_adj_list,
    data_adj = data_adj
  )
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
#' @param precomputed_adj Optional internal adjustment object returned by
#' \code{build_adjustment_step()}. When supplied, adjustment-variable
#' transformations are reused and only current-variable candidate
#' transformations are generated.
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
#' @noRd
transform_data_step <- function(x,
                                xi,
                                powers_current,
                                df,
                                powers,
                                acdx,
                                zero,
                                catzero,
                                spike,
                                spike_decision,
                                acd_parameter,
                                prev_adj_params,
                                precomputed_adj = NULL
) {
  # Access x columns by name; avoid copying/reordering x for every candidate
  # transformation call.
  if (is.null(precomputed_adj)) {
    adj <- build_adjustment_step(
      x = x,
      xi = xi,
      powers_current = powers_current,
      powers = powers,
      acdx = acdx,
      zero = zero,
      catzero = catzero,
      spike = spike,
      spike_decision = spike_decision,
      acd_parameter = acd_parameter,
      prev_adj_params = prev_adj_params
    )
  } else {
    adj <- precomputed_adj
  }
  
  # generate fp data for xi (unchanged, spike_decision not applied here)
  data_xi <- x[, xi, drop = TRUE]
  
  if (length(unique(data_xi)) <= 3) {
    # if a variable has less than 4 levels we do not generate FP data
    # but set nonpositive values to zero, when zero argument is used
    if (zero[xi]) {
      data_xi[data_xi <= 0] <- 0
    }
    
    # add catzero variable if not null; if null cbind will remove it
    data_xi_mat <- matrix(data_xi, ncol = 1L)
    
    if (!is.null(catzero[[xi]])) {
      data_xi_mat <- cbind(catzero[[xi]], data_xi_mat)
      colnames(data_xi_mat)[1L] <- "catzero"
    }
    
    data_fp <- list(data_xi_mat)
    powers_fp <- matrix(1, nrow = 1L, ncol = 1L)
    
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
    powers_adj = adj$powers_adj,
    spike_decision_adj = adj$spike_decision_adj,
    data_adj_list = adj$data_adj_list,
    data_adj = adj$data_adj
  )
  
  # Return results and current parameters for next step
  list(
    powers_fp = powers_fp,
    data_fp = data_fp,
    powers_adj = adj$powers_adj,
    data_adj = adj$data_adj,
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
#' @keywords internal
#' @noRd
ensure_length <- function(x, size, fill = NA) {
  if (length(x) == size)
    return(x)
  
  x_new = rep(NA, size)
  x_new[1:length(x)] = x
  
  x_new
}
