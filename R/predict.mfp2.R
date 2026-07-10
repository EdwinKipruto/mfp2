#' Predict Method for `mfp2`
#'
#' Obtain predictions from a fitted `mfp2` object.
#'
#' @description
#' This method predicts from the final fitted model stored in an `mfp2` object.
#' For ordinary predictions it delegates to the underlying `predict.glm()` or
#' `predict.coxph()` method after reconstructing the transformed design matrix
#' used during fitting. For `type = "terms"` and `type = "contrasts"`, it returns
#' variable-specific partial predictors with conditional standard errors.
#'
#' @details
#' If `newdata` is supplied, the variables are prepared using the transformation
#' metadata saved in `object`: fitted shifts, scales, selected fractional
#' polynomial powers, centering constants, zero-handling indicators, catzero
#' indicators, and spike-at-zero decisions. This ensures that prediction uses the
#' same transformed design scale as model fitting.
#'
#' Prediction stops when shifted `newdata` is outside the domain required by the
#' final fitted transformation. This check is applied only to variables whose
#' selected powers require strictly positive input, such as logarithmic powers,
#' repeated powers, negative powers, or non-integer powers. Variables whose final
#' fitted zero-handled or binary-only SAZ representation permits non-positive
#' values are exempt from this positivity stop.
#'
#' For full-model predictions, standard errors are requested with `se.fit = TRUE`
#' and are computed by the underlying `predict.glm()` or `predict.coxph()` method.
#' For Cox models, this method always calls `predict.coxph()` with
#' `reference = "zero"`. The transformed design matrix in `mfp2` has already been
#' centered where required, so asking `predict.coxph()` to subtract an additional
#' sample or strata reference would shift the linear predictor incorrectly.
#'
#' For `type = "terms"`, standard errors are computed from the fitted covariance
#' matrix for the columns belonging to each selected variable. If
#' `add_intercept = TRUE` and the model is not Cox, the intercept contribution is
#' included in both the term value and its standard error. If
#' `add_intercept = FALSE`, the intercept contribution is excluded from both.
#'
#' For `type = "contrasts"`, the returned value is the difference between the
#' partial predictor at each value and the partial predictor at the reference
#' value. Intercepts do not contribute to contrasts, and the standard errors are
#' computed from the transformed difference. These standard errors are
#' conditional on the final selected model and do not include model-selection
#' uncertainty.
#'
#' @section Terms prediction:
#' If `type = "terms"`, this function computes partial linear predictors for
#' selected variables in the final model. A single original variable may be
#' represented by multiple fitted columns, for example FP2 terms, zero-handled
#' positive-part terms, or a catzero binary indicator. The method collects the
#' relevant fitted columns for each variable and multiplies them by their fitted
#' coefficients.
#'
#' @section Contrasts:
#' If `type = "contrasts"`, this function computes variable-specific contrasts
#' relative to reference values. Reference values supplied through `ref` must be
#' on the original data scale; they are shifted, transformed, and centered using
#' the fitted transformation metadata. If `ref = NULL` for a variable, the method
#' uses the mean shifted value for continuous variables and the minimum shifted
#' value for binary-like variables.
#'
#' @param object A fitted object of class `mfp2`.
#' @param newdata Optional matrix or data frame containing variables for
#'   prediction. Column names must identify the original predictors used by the
#'   fit. Formula-created predictors are reconstructed when possible.
#' @param type Prediction type. The default is `"link"` for GLM models and
#'   `"lp"` for Cox models. Use `"terms"` for variable-specific partial
#'   predictors and `"contrasts"` for variable-specific contrasts.
#' @param se.fit Logical scalar. For full-model predictions, if `TRUE`, request
#'   standard errors from `predict.glm()` or `predict.coxph()`. For
#'   `type = "terms"` and `type = "contrasts"`, the returned data frames always
#'   contain an `se` column, and `se.fit` is ignored.
#' @param terms Character vector of original variable names for term or contrast
#'   prediction. Only variables selected in the final fitted model are used. If
#'   `NULL`, all selected variables are used.
#' @param terms_seq Character scalar controlling values used for term prediction.
#'   `"equidistant"` generates `nseq` equally spaced values over the observed
#'   range. `"data"` uses observed values directly.
#' @param alpha Significance level used for confidence intervals in term and
#'   contrast predictions.
#' @param ref Named list of reference values for `type = "contrasts"`. Values
#'   must be supplied on the original variable scale.
#' @param strata Optional stratum values used when predicting from Cox models
#'   with supplied \code{newdata}. This is needed for default/matrix-interface
#'   Cox fits that used the fit-time \code{strata} argument, because the final
#'   Cox model contains an internal \code{strata(strata_)} term. Supply one
#'   value per prediction row for a single stratification factor, or a matrix or
#'   data frame with one row per prediction row for multiple stratification
#'   factors. Ordinary vector or factor strata are passed through as raw
#'   high-level values so that \code{predict.coxph()} can evaluate the stored
#'   \code{strata(strata_)} term itself. Formula-interface strata are
#'   reconstructed automatically from \code{newdata} when the original strata
#'   variable(s) are present.
#' @param newoffset Optional numeric vector of offsets for prediction when the
#'   fitted model used an offset and `newdata` is supplied.
#' @param nseq Positive integer giving the number of equally spaced values used
#'   when `terms_seq = "equidistant"`.
#' @param add_intercept Logical scalar. For `type = "terms"`, controls whether
#'   the model intercept is included in GLM term values and term standard errors.
#'   It has no effect for contrasts and does not apply to Cox models.
#' @param ... Further arguments passed to `predict.glm()` or `predict.coxph()`
#'   for full-model predictions.
#'
#' @return
#' For full-model predictions, the return value follows `predict.glm()` or
#' `predict.coxph()`. In particular, if `se.fit = TRUE`, the result may be a list
#' containing fitted values and standard errors according to the underlying
#' method.
#'
#' For `type = "terms"` or `type = "contrasts"`, the result is a named list with
#' one data frame per selected variable. Each data frame contains:
#' \itemize{
#'   \item `variable`: values on the original scale before fitted shifting.
#'   \item `variable_pre`: values after fitted shifting, before or after binary
#'     SAZ coding as appropriate.
#'   \item `value`: partial linear predictor or contrast.
#'   \item `se`: conditional standard error.
#'   \item `lower`: lower confidence limit.
#'   \item `upper`: upper confidence limit.
#' }
#'
#' @examples
#' data("prostate")
#' x <- as.matrix(prostate[, 2:8])
#' y <- as.numeric(prostate$lpsa)
#'
#' fit <- mfp2(x, y, verbose = FALSE)
#' predict(fit)
#' predict(fit, se.fit = TRUE)
#' predict(fit, type = "terms")
#'
#' @seealso
#' \code{mfp2()}, [stats::predict.glm()], [survival::predict.coxph()]
#'
#' @method predict mfp2
#' @export
predict.mfp2 <- function(object, 
                         newdata = NULL, 
                         type = NULL,
                         se.fit = FALSE,
                         terms = NULL,
                         terms_seq = c("equidistant", "data"),
                         alpha = 0.05,
                         ref = NULL, 
                         strata = NULL, 
                         newoffset = NULL, 
                         nseq = 100,
                         add_intercept = TRUE,
                         ...) {
  
  terms_seq <- match.arg(terms_seq)
  
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be between 0 and 1.", call. = FALSE)
  }
  
  if (!is.numeric(nseq) || nseq <= 0) {
    stop("'nseq' must be a positive integer.", call. = FALSE)
  }
  
  if (!is.logical(se.fit) || length(se.fit) != 1L || is.na(se.fit)) {
    stop("'se.fit' must be a single TRUE or FALSE value.", call. = FALSE)
  }
  
  # assert that the object must be mfp2
  if (!inherits(object, "mfp2")) { 
    stop("The object is not an mfp2 object.", call. = FALSE)
  }
  
  # set defaults and match arguments
  if (is.null(type)) {
    type <- ifelse(object$family_string == "cox", "lp", "link")
  }
  
  if (is.null(terms)) {
    terms <- get_selected_variable_names(object)
  }
  
  if (is.null(ref)) {
    ref <- setNames(lapply(terms, function(v) NULL), terms)
  }
  
  # checks for newdata
  if (!is.null(newdata) && is.null(colnames(newdata))) {
    stop("Newdata must have column names", call. = FALSE)
  }
  
  newdata_raw <- newdata
  
  if (!is.null(newdata)) {
    newdata <- reconstruct_formula_newdata(object, newdata)
  }
  
  if (!is.null(newdata) && anyNA(newdata)) {
    stop("! newdata must not contain any NA (missing data).\n", 
         "i Please remove any missing data before passing newdata to this function.",
         call. = FALSE)
  }
  
  # TODO: add checks for missing strata and offset in case they were used in fit
  if (type == "contrasts" && length(ref) != sum(names(ref) != "", na.rm = TRUE)) {
    warning(
      paste0(
        "i The supplied reference values (ref) must all be named.\n",
        "i predict() continues but uses means (if variables are continuous) ",
        "or min (if binary) instead of the reference values."
      ),
      call. = FALSE
    )
  }
  
  if (type %in% c("terms", "contrasts")) {
    n_term1 <- length(terms)
    terms <- intersect(terms, get_selected_variable_names(object))
    # length of terms after intersections
    n_term2 <- length(terms)
    if (n_term2 == 0) {
      warning("i All the terms supplied are not in the final model.\n", 
              "i predict() continues but returns an empty list.", call. = FALSE) 
    } else if (n_term2 < n_term1)
      warning(
        paste0(
          "i Some terms supplied are not in the final model.\n",
          "i predict() continues but returns an empty list for those terms ",
          "not in the model."
        ),
        call. = FALSE
      )
    # return warning if the names(ref) != names(terms)
    if (!all(sapply(ref, is.null)) && any(!names(ref) %in% terms))
      warning("i Some of names of reference values are not in terms.\n", 
              "i predict() continues but does not consider them.", call. = FALSE)
    
    # Extract the fitted intercept for GLM term predictions. The intercept is
    # included only when add_intercept = TRUE and type = "terms". It is never
    # included for Cox partial predictors or contrasts.
    cf <- coef(object)
    if ("(Intercept)" %in% names(cf)) {
      intercept <- cf["(Intercept)"]
    } else {
      intercept <- 0
    }
    
    # Remove intercept if NA or not requested
    if (is.na(intercept) || !add_intercept) {
      intercept <- 0
    }
    
    # Intercept does not apply for contrasts
    if (type == "contrasts") {
      intercept <- 0
    }
    res_list <- list()
    for (t in terms) {
      # define sequence of variable data as named list
      if (terms_seq == "equidistant") {
        if (!is.null(newdata)) {
          x_range <- range(newdata[, t]) + object$transformations[t,"shift"]
        }else {
          x_range <- range(object$x_original[, t]) # already shifted
        }
        
        x_seq  <- matrix(
          seq(x_range[1], x_range[2], length.out = nseq),
          ncol = 1
        )
        colnames(x_seq) <- t
        
        # no need to apply pretransformation (shift), already done
        x_trafo <- prepare_newdata_for_predict(object, 
                                               x_seq, 
                                               apply_pre = FALSE,
                                               allow_missing_predictors = TRUE)
        x_trafo <- as.matrix(x_trafo)
      } else {
        # no equidistant
        if (!is.null(newdata)) {
          x_seq <- newdata[, t, drop = FALSE] + object$transformations[t,"shift"]
        } else {
          x_seq <- object$x_original[, t, drop = FALSE] # already shifted
        }
        
        #x_names <- object$fp_powers[[t]]
        # in acd we might have a power and NA so we need to remove NA
        #x_trafo <- object$x[, names(x_names[!is.na(x_names)]), drop = FALSE]
        x_trafo <- as.matrix(prepare_newdata_for_predict(
          object, 
          x_seq, 
          apply_pre = FALSE,
          allow_missing_predictors = TRUE))
      }
      # Use the prepared design column names to select coefficients. This is
      # essential for FPm, zero-handled variables, catzero indicators, and SAZ
      # decisions because one original variable may map to multiple fitted
      # columns.
      term_coef <- coef(object)[colnames(x_trafo)]
      
      # Replace NA coefficients with 0 (rank-deficient or aliased terms),
      # matching the convention in predict.coxph().
      term_coef[is.na(term_coef)] <- 0
      
      # create output data.frame
      # backtransform variable to original scale
      variable = (as.numeric(x_seq)) - object$transformations[t,"shift"]
      variable_pre = as.numeric(x_seq)
      
      if (object$spike_dec[t] == saz_decision_codes[["binary_only"]]) {
        #variable <- object$catzero_list[[t]]
        #variable_pre <- variable
        variable <- as.integer(x_seq[, t] <= 0)
        variable_pre <- variable
      }
      
      res <- data.frame(
        variable = variable,
        variable_pre = variable_pre,
        value = x_trafo %*% term_coef + intercept
      )
      
      # For term predictions there is no reference value. For contrasts,
      # x_ref_trafo is the fitted-design representation of the reference point.
      x_ref_trafo <- NULL
      if (type == "contrasts") {
        # compute transformations for reference level
        # note that intercepts do not play a role here and that
        # (f(x) - f(x_ref)) * coef == f(x) * coef - f(x_ref) * coef
        
        x_ref <- ref[[t]]
        if (is.null(x_ref)) {
          if (!is.null(newdata)) {
            v <- newdata[, t] + object$transformations[t,"shift"]
          } else {
            v <- object$x_original[, t] # already shifted
          }
          x_ref <- ifelse(length(unique(na.omit(v))) == 2,
                          min(v, na.rm = TRUE),
                          mean(v, na.rm = TRUE))
          
        } else {
          # pretransform given reference level
          x_ref <- (x_ref + object$transformations[t,"shift"]) 
        }
        # make sure it is a named matrix
        x_ref <- matrix(x_ref, nrow = 1, ncol = 1)
        colnames(x_ref) <- t
        
        # transform x_ref using estimated FP powers
        x_ref_trafo <- as.matrix(prepare_newdata_for_predict(
          object, x_ref, apply_pre = FALSE, check_binary = FALSE,
          reset_zero = FALSE, allow_missing_predictors = TRUE))
        
        # compute contrasts, no intercepts necessary
        res$value <- res$value - as.numeric(x_ref_trafo %*% term_coef)
      }
      
      res$se <- calculate_standard_error(
        object,
        x_trafo,
        x_ref_trafo,
        include_intercept = add_intercept && type == "terms"
      )
      mult <- qnorm(1 - (alpha / 2))
      res$lower <- res$value - mult * res$se
      res$upper <- res$value + mult * res$se
      
      res_list[[t]] <- res
    }
    names(res_list) <- terms
    
    return(res_list)
  } 
  
  # Full-model predictions. At this point newdata, if supplied, is converted
  # to the same transformed design scale used when fitting the final model.
  
  if (!is.null(newdata)) {
    
    if (
      is.null(newoffset) &&
      !is.null(newdata_raw) &&
      !is.null(object$formula_offset_terms)
    ) {
      newoffset <- reconstruct_formula_offset_newdata(object, newdata_raw)
    }
    
    # check whether offset was used in the model
    has_offset <- isTRUE(object$has_offset)
    if (has_offset) {
      if (is.null(newoffset)) {
        stop(
          "No newoffset provided for prediction, yet offset was used in mfp2",
          call. = FALSE
        )
      }
      
      if (!is.numeric(newoffset)) {
        stop("! newoffset must be numeric.", call. = FALSE)
      }
      
      if (length(newoffset) != nrow(newdata)) {
        stop(
          "The length of newoffset must be equal to the number of rows of newdata",
          call. = FALSE
        )
      }
      
    } else {
      newoffset <- NULL
    }
    
    if (
      object$family_string == "cox" &&
      is.null(strata) &&
      !is.null(newdata_raw) &&
      !is.null(object$formula_strata_terms)
    ) {
      strata <- reconstruct_formula_strata_newdata(object, newdata_raw)
    }
    
    newdata <- prepare_newdata_for_predict(
      object,
      newdata,
      strata = strata,
      offset = newoffset,
      check_binary = FALSE
    )
    
    # Strip "mfp2" from the class vector so stats::predict() dispatches to
    # predict.coxph() or predict.glm() instead of recursing into predict.mfp2().
    obj_base <- object
    class(obj_base) <- setdiff(class(obj_base), "mfp2")
    
    if (object$family_string == "cox") {
      # Cox reference convention: mfp2 stores and predicts on an already
      # prepared design scale. Use reference = "zero" so predict.coxph() does
      # not subtract its own sample or stratum means a second time.
      pred <- stats::predict(
        obj_base,
        newdata = newdata,
        type = type,
        se.fit = se.fit,
        reference = "zero",
        ...
      )
    } else {
      pred <- stats::predict(
        obj_base,
        newdata = newdata,
        type = type,
        se.fit = se.fit,
        ...
      )
    }
    
    return(pred)
  }
  
  # no newdata supplied
  obj_base <- object
  class(obj_base) <- setdiff(class(obj_base), "mfp2")
  if (object$family_string == "cox") {
    return(
      stats::predict(
        obj_base,
        type = type,
        se.fit = se.fit,
        reference = "zero",
        ...
      )
    )
  } else {
    return(
      stats::predict(
        obj_base,
        type = type,
        se.fit = se.fit,
        ...
      )
    )
  }
}

#' Transform Linear Predictor to Response or Risk
#' 
#' Converts linear predictors (`nfit`) from a model to the appropriate scale
#' for interpretation or prediction, depending on the model family, link function,
#' and type of prediction. This is an internal helper function for GLMs, survival models,
#' and other regression frameworks.
#' 
#' @param nfit Numeric vector of linear predictors \eqn{X\beta}.
#' @param family Character string specifying the model family. Supported families include:
#'   - `"gaussian"`: linear regression
#'   - `"binomial"`: binary outcomes or proportions
#'   - `"poisson"`: count data
#'   - `"cox"`: proportional hazards model
#'   - other families fallback to returning `nfit`
#' @param link Character string specifying the link function. Defaults:
#'   - Gaussian: `"identity"`  
#'   - Binomial: `"logit"`  
#'   - Poisson: `"log"`  
#'   Supported links:
#'   - Gaussian: `"identity"`, `"log"`, `"inverse"`  
#'   - Binomial: `"logit"`, `"probit"`, `"cloglog"`, `"cauchit"`, `"identity"`, `"log"`  
#'   - Poisson: `"log"`, `"identity"`, `"sqrt"`
#' @param type Character string specifying the type of prediction:
#'   - `"response"` or `"risk"`: returns predictions on the response/probability scale
#'   - `NULL` (default): returns the linear predictor itself
#' 
#' @return Numeric vector of predictions on the requested scale.
#' 
#' @examples
#' \dontrun{
#' # Binomial example
#' lp_bin <- c(-1, 0, 1)
#' transform_linear_predictor(lp_bin, family = "binomial", link = "logit", type = "response")
#' 
#' # Poisson example
#' lp_pois <- c(0.5, 1, 1.5)
#' transform_linear_predictor(lp_pois, family = "poisson", link = "log", type = "response")
#' 
#' # Gaussian example
#' lp_gauss <- c(1, 2, 3)
#' transform_linear_predictor(lp_gauss, family = "gaussian", link = "log", type = "response")
#' }
#' @keywords internal
#' @noRd
transform_linear_predictor <- function(nfit, family, link = NULL, type = NULL) {
  
  # Set default links if not provided
  if (is.null(link)) {
    link <- switch(family,
                   gaussian = "identity",
                   binomial = "logit",
                   poisson  = "log",
                   NULL)
  }
  
  switch(family,
         
         # Gaussian family with multiple links
         gaussian = {
           if (!is.null(type) && type == "response") {
             switch(link,
                    identity = nfit,
                    log      = exp(nfit),
                    inverse  = 1 / nfit,
                    stop("Unknown Gaussian link"))
           } else nfit
         },
         
         # Binomial family with multiple links
         binomial = {
           if (!is.null(type) && type == "response") {
             switch(link,
                    logit   = 1 / (1 + exp(-nfit)),
                    probit  = pnorm(nfit),
                    cloglog = 1 - exp(-exp(nfit)),
                    cauchit = pcauchy(nfit),
                    identity = nfit,
                    log     = exp(nfit),
                    stop("Unknown binomial link"))
           } else nfit
         },
         
         # Poisson family with multiple links
         poisson = {
           if (!is.null(type) && type == "response") {
             switch(link,
                    log      = exp(nfit),
                    identity = nfit,
                    sqrt     = nfit^2,
                    stop("Unknown Poisson link"))
           } else nfit
         },
         
         # Cox proportional hazards: exponentiate for risk
         cox = {
           if (!is.null(type) && type %in% c("response", "risk")) exp(nfit) else nfit
         },
         
         # Default fallback: return linear predictor
         nfit
  )
}

#' Does a fitted term require positive raw prediction input?
#'
#' Internal helper used by `predict.mfp2()`.
#'
#' The prediction method applies the shift values learned during fitting before
#' constructing transformed variables. This helper determines whether the
#' shifted raw covariate for variable `v` must be strictly positive, based on
#' the final fitted model metadata.
#'
#' For ordinary FP terms, all selected powers in `object$fp_powers[[v]]` apply
#' directly to the raw covariate.
#'
#' For ACD terms, `object$fp_powers[[v]][1]` applies to the raw covariate and
#' `object$fp_powers[[v]][2]` applies to the ACD component A(x). Therefore only
#' the first power is checked against the raw covariate here.
#'
#' Variables fitted with zero-handling are exempt because nonpositive values are
#' represented structurally. Variables with `spike_decision = 3` are also exempt
#' because the continuous component is dropped and only the spike/catzero
#' indicator is retained.
#'
#' @param object A fitted `mfp2` object.
#' @param v Character scalar giving the variable name.
#'
#' @return A single logical value.
#'
#' @keywords internal
#' @noRd
requires_positive_raw_input <- function(object, v) {
  powers <- object$fp_powers[[v]]
  
  if (is.null(powers) || all(is.na(powers))) {
    return(FALSE)
  }
  
  if (!is.null(object$zero) && isTRUE(object$zero[[v]])) {
    return(FALSE)
  }
  
  if (!is.null(object$spike_dec) &&
      v %in% names(object$spike_dec) &&
      !is.na(object$spike_dec[[v]]) &&
      as.integer(object$spike_dec[[v]]) == saz_decision_codes[["binary_only"]]) {
    return(FALSE)
  }
  
  is_acd <- !is.null(object$acd) &&
    v %in% names(object$acd) &&
    isTRUE(object$acd[[v]])
  
  raw_powers <- if (is_acd) {
    # For ACD terms:
    #   powers[1] applies to raw x
    #   powers[2] applies to A(x)
    #
    # Only the raw-x power determines whether shifted raw newdata must be
    # strictly positive at this stage.
    powers[1L]
  } else {
    powers
  }
  
  fp_power_requires_positive_input(raw_powers)
}

#' Rebuild formula-interface newdata with model.matrix()
#'
#' Formula-fitted mfp2 objects are fitted on the expanded numeric design matrix
#' returned by model.matrix(). This helper lets users pass ordinary newdata with
#' the original formula variables, including factors, and reconstructs the same
#' expanded columns used during fitting. If newdata already contains all fitted
#' design-matrix columns, it is returned unchanged for backwards compatibility.
#'
#' @keywords internal
#' @noRd
reconstruct_formula_newdata <- function(object, newdata) {
  if (!isTRUE(object$formula_interface)) {
    return(newdata)
  }
  
  expected <- rownames(object$transformations)
  if (!is.null(expected) && length(expected) > 0L &&
      all(expected %in% colnames(newdata))) {
    return(newdata)
  }
  
  if (is.null(object$formula_terms)) {
    stop(
      "! This object was fitted using the formula interface, but formula terms ",
      "needed for prediction were not stored.",
      call. = FALSE
    )
  }
  
  newdata_df <- as.data.frame(newdata)
  
  mf <- stats::model.frame(
    object$formula_terms,
    data = newdata_df,
    na.action = stats::na.pass,
    xlev = object$formula_xlevels
  )
  
  mm <- stats::model.matrix(
    object$formula_terms,
    data = mf,
    contrasts.arg = object$formula_contrasts
  )
  
  keep_cols <- colnames(mm) != "(Intercept)"
  mm <- mm[, keep_cols, drop = FALSE]
  
  column_map <- object$formula_column_map
  if (!is.null(column_map)) {
    mapped <- unname(column_map[colnames(mm)])
    colnames(mm) <- ifelse(is.na(mapped), colnames(mm), mapped)
  }
  
  if (!is.null(object$formula_model_matrix_columns)) {
    missing <- setdiff(object$formula_model_matrix_columns, colnames(mm))
    if (length(missing) > 0L) {
      stop(
        "! Could not reconstruct required model-matrix column(s) from newdata: ",
        paste(missing, collapse = ", "),
        call. = FALSE
      )
    }
    mm <- mm[, object$formula_model_matrix_columns, drop = FALSE]
  }
  
  mm
}

#' Reconstruct formula-level Cox strata from prediction data
#'
#' Reconstructs the original formula-level Cox stratification term from
#' `newdata` for `mfp2` objects fitted through the formula interface.
#'
#' This helper deliberately mimics the high-level `survival::coxph()` strata
#' handling:
#'
#' - for a single strata variable, it returns the evaluated column from the
#'   model frame;
#' - for multiple strata variables, it combines them with
#'   `survival::strata(..., shortlabel = TRUE)`;
#' - it does **not** convert the result to integer codes.
#'
#' Integer conversion belongs only in the low-level Cox fitting path immediately
#' before calling `survival::coxph.fit()`, analogous to the `istrat` object used
#' internally by `survival::coxph()`.
#'
#' @param object A fitted `mfp2` object.
#' @param newdata A user-supplied prediction data frame before formula
#'   reconstruction has been applied.
#'
#' @return Either `NULL`, if no formula-level strata metadata are available, or
#'   a Cox strata object suitable for passing as the `strata` argument to
#'   `prepare_newdata_for_predict()`.
#'
#' @keywords internal
#' @noRd
reconstruct_formula_strata_newdata <- function(object, newdata) {
  if (!isTRUE(object$formula_interface) ||
      is.null(object$formula_strata_terms)) {
    return(NULL)
  }
  
  newdata_df <- as.data.frame(newdata)
  
  mf <- tryCatch(
    stats::model.frame(
      object$formula_strata_terms,
      data = newdata_df,
      na.action = stats::na.pass,
      xlev = object$formula_strata_xlevels
    ),
    error = function(e) {
      stop(
        "! This `mfp2` object was fitted with formula-level Cox strata, ",
        "but the strata term could not be reconstructed from `newdata`.\n",
        "i Include the original strata variable(s) in `newdata`, or supply ",
        "`strata = ...` explicitly to `predict()`.\n",
        "i Original error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  if (NROW(mf) != NROW(newdata_df)) {
    stop(
      "! Formula-level Cox strata reconstruction returned the wrong number of rows.",
      call. = FALSE
    )
  }
  
  if (anyNA(mf)) {
    stop(
      "! Reconstructed Cox strata contain missing values.\n",
      "i Please remove missing strata values from `newdata` before prediction.",
      call. = FALSE
    )
  }
  
  if (NCOL(mf) == 1L) {
    mf[[1L]]
  } else {
    do.call(
      survival::strata,
      c(as.list(mf), list(shortlabel = TRUE))
    )
  }
}

#' Rebuild Formula-Level Offset from Prediction Newdata
#'
#' Internal helper used by predict.mfp2(). If an mfp2 object was fitted with
#' offset() in the formula interface, this reconstructs the offset vector from
#' ordinary newdata.
#'
#' @keywords internal
#' @noRd
reconstruct_formula_offset_newdata <- function(object, newdata) {
  if (!isTRUE(object$formula_interface) ||
      is.null(object$formula_offset_terms)) {
    return(NULL)
  }
  
  newdata_df <- as.data.frame(newdata)
  
  mf <- tryCatch(
    stats::model.frame(
      object$formula_offset_terms,
      data = newdata_df,
      na.action = stats::na.pass,
      xlev = object$formula_offset_xlevels
    ),
    error = function(e) {
      stop(
        "! This `mfp2` object was fitted with a formula-level offset, ",
        "but the offset could not be reconstructed from `newdata`.\n",
        "i Include the original offset variable(s) in `newdata`, or supply ",
        "`newoffset = ...` explicitly to `predict()`.\n",
        "i Original error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  out <- stats::model.offset(mf)
  
  if (is.null(out)) {
    stop(
      "! Formula-level offset could not be evaluated from `newdata`.",
      call. = FALSE
    )
  }
  
  out <- as.vector(out)
  
  if (length(out) != NROW(newdata_df)) {
    stop(
      "! Formula-level offset reconstruction returned the wrong number of rows.",
      call. = FALSE
    )
  }
  
  if (anyNA(out) || any(!is.finite(out))) {
    stop(
      "! Reconstructed offset contains missing or non-finite values.",
      call. = FALSE
    )
  }
  
  out
}

#' Helper function to prepare newdata for predict function
#' 
#' To be used in \code{predict.mfp2()}.
#' 
#' @param object fitted `mfp2` model object.
#' @param newdata dataset to be prepared for predictions. Its columns can be
#' a subset of the columns used for fitting the model. 
#' @param strata,offset passed from \code{predict.mfp2()}. For Cox
#'   prediction, \code{strata} must be kept as a high-level vector/factor or
#'   combined multi-column strata object; integer conversion is not performed
#'   here because the stored Cox formula evaluates \code{strata(strata_)}.
#' @param apply_pre logical indicating whether the fitted pre-transformation
#' is applied or not.
#' @param apply_center logical indicating whether the fitted centers are applied
#' after transformation or not.
#' @param check_binary passed to \code{transform_vector_fp()}.
#' @param reset_zero Logical. If `TRUE`, variables marked as `zero = TRUE`
#' but containing only positive values are reset to `FALSE` before transformation.
#' The prediction helper defaults to `FALSE` because prediction must preserve
#' the zero/catzero/spike structure learned during model fitting. Adaptive
#' resetting based on the contents of `newdata` can create a design matrix that
#' is inconsistent with the fitted coefficients. Parameter of
#' \code{transform_matrix()}.
#' @param allow_missing_predictors Logical. If `FALSE`, the default, `newdata`
#'   must contain all original predictors recorded in `object$transformations`.
#'   Missing predictors trigger an error before the data are subset or
#'   transformed. If `TRUE`, missing predictors are allowed and only the
#'   predictors present in both `newdata` and the fitted object are transformed.
#'   This is intended for internal partial-prediction or term-specific
#'   prediction paths, not for ordinary full-model prediction.
#' @return A dataframe of transformed newdata
#' @keywords internal
#' @noRd
prepare_newdata_for_predict <- function(object, 
                                        newdata, 
                                        strata = NULL, 
                                        offset = NULL, 
                                        apply_pre = TRUE, 
                                        apply_center = TRUE,
                                        check_binary = TRUE,
                                        reset_zero = FALSE,
                                        allow_missing_predictors = FALSE
) {
  newdata <- as.matrix(newdata)
  n_newdata <- nrow(newdata)
  
  # step 0: check that newdata has column names
  if (is.null(colnames(newdata)) || any(colnames(newdata) == "")) {
    stop("! newdata must have non-empty column names.", call. = FALSE)
  }
  
  # step 1: check expected predictors before subsetting
  # expected predictors are the original variables for which transformations
  # were stored during model fitting.
  expected <- rownames(object$transformations)
  
  missing <- setdiff(expected, colnames(newdata))
  if (length(missing) > 0L && !isTRUE(allow_missing_predictors)) {
    stop(
      "! Missing required predictor(s) in newdata: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  
  # step 2: subset as appropriate
  # keep only predictors known to the fitted object. Extra columns in newdata
  # are ignored. The order is determined by the intersection and then used
  # consistently in all downstream metadata lookups through vnames.
  vnames <- intersect(colnames(newdata), expected)
  
  if (length(vnames) == 0L) {
    stop("! No usable predictors were found in newdata.", call. = FALSE)
  }
  
  # sorting is not relevant as we always pass vnames
  newdata <- newdata[, vnames, drop = FALSE]
  
  if (apply_pre) {
    # step 3: shift data using shifting factors from the training data.
    # The final model coefficients are already on the original scale after
    # shifting, so scaling is intentionally not applied here.
    newdata <- sweep(newdata, 2, object$transformations[vnames, "shift"], "+")
    
    # Step 4: validate positive-domain requirements after applying fitted shifts.
    #
    # This check uses the final fitted model metadata rather than the values in
    # `newdata` to decide which variables require strictly positive input. This is
    # important because a prediction batch may contain only one or two unique values
    # for a continuous FP variable.
    positive_raw_vars <- vnames[
      vapply(
        vnames,
        function(v) requires_positive_raw_input(object, v),
        logical(1L)
      )
    ]
    
    if (length(positive_raw_vars) > 0L) {
      bad_vars <- positive_raw_vars[
        vapply(
          positive_raw_vars,
          function(v) {
            any(!is.na(newdata[, v]) & newdata[, v] <= 0)
          },
          logical(1L)
        )
      ]
      
      if (length(bad_vars) > 0L) {
        bad_counts <- vapply(
          bad_vars,
          function(v) {
            sum(!is.na(newdata[, v]) & newdata[, v] <= 0)
          },
          integer(1L)
        )
        
        bad_summary <- paste0(
          bad_vars,
          " (",
          bad_counts,
          ifelse(bad_counts == 1L, " row", " rows"),
          ")"
        )
        
        stop(
          "After applying the shift values learned during fitting, some values in ",
          "`newdata` remain non-positive for variables whose fitted transformation ",
          "requires strictly positive input.\n",
          "i Problematic variable(s): ",
          paste(bad_summary, collapse = ", "),
          ".\n",
          "i Prediction is undefined for these values. Refit with larger shift ",
          "values, restrict `newdata` to the fitted domain, or use zero-handling ",
          "where appropriate.",
          call. = FALSE
        )
      }
    }
    
    # already the coefficients are in original scale after shifting
    # newdata <- sweep(newdata, 2, object$transformations[vnames, "scale"], "/")
  }
  
  # Prediction must reuse the ACD parameters estimated on the training data
  # whenever an ACD-transformed continuous component is active. Otherwise,
  # transform_vector_acd() could refit the ACD transformation on newdata,
  # making predictions depend on the prediction batch and hiding upstream
  # naming/order bugs in object$acd_parameter.
  #
  # ACD variables whose powers are all NA are skipped here. Such variables were
  # either eliminated during model selection or, for spike_decision = 3, retained
  # only through the binary spike indicator. In both cases, transform_vector_acd()
  # returns NULL before applying ACD parameters, so no stored parameter object is
  # needed for prediction.
  acd_vars <- vnames[
    vapply(
      vnames,
      function(v) {
        isTRUE(object$fp_terms[v, "acd"]) &&
          !all(is.na(object$fp_powers[[v]]))
      },
      logical(1)
    )
  ]
  
  missing_acd <- acd_vars[
    vapply(
      acd_vars,
      function(v) {
        is.null(object$acd_parameter) ||
          !v %in% names(object$acd_parameter) ||
          is.null(object$acd_parameter[[v]])
      },
      logical(1)
    )
  ]
  
  if (length(missing_acd) > 0L) {
    stop(
      "Missing stored ACD parameters for prediction variable(s): ",
      paste(missing_acd, collapse = ", "),
      call. = FALSE
    )
  }
  
  # Step 5: transform shifted data using the fitted transformation metadata.
  # Prediction must preserve the fitted zero/catzero/spike structure. Therefore
  # reset_zero defaults to FALSE in this helper; otherwise a newdata batch that
  # contains only positive values could drop a fitted catzero indicator and create
  # a prediction design matrix inconsistent with the fitted coefficients.
  # do not center in this step
  x_trans <- transform_matrix(
    newdata,
    power_list = object$fp_powers[vnames], 
    center = setNames(rep(FALSE, length(vnames)), vnames),
    keep_x_order = TRUE,
    acdx = setNames(object$fp_terms[vnames, "acd"], vnames),
    acd_parameter_list = object$acd_parameter[vnames],
    check_binary = check_binary,
    zero = object$zero[vnames],
    catzero = object$catzero[vnames],
    spike = if (!is.null(object$spike)) {
      object$spike[vnames]
    } else {
      object$fp_terms[vnames, "spike"]
    },
    spike_decision = object$spike_dec[vnames],
    reset_zero = reset_zero
  )
  
  newdata <- x_trans$x_transformed
  
  # step 6: center the transformed data
  if (apply_center && !is.null(object$centers)) {
    newdata <- center_matrix(
      newdata, 
      centers = object$centers[colnames(newdata)],
      zero = x_trans$zero_expanded
    ) 
  }
  
  # step 7: convert to data frame for downstream predict methods.
  # Intercept-only final models may have no transformed predictor columns, but
  # prediction must still preserve the number of rows in newdata.
  if (NCOL(newdata) == 0L) {
    newdata <- data.frame(row.names = seq_len(n_newdata))
  } else {
    newdata <- data.frame(newdata)
  }
  
  # step 8: add Cox strata column if required
  if (object$family_string == "cox") {
    if (!is.null(strata)) {
      strata_n <- if (is.vector(strata) || is.factor(strata)) {
        length(strata)
      } else {
        NROW(strata)
      }
      
      if (strata_n != nrow(newdata)) {
        stop(
          "! `strata` must have one value or row per prediction row.",
          call. = FALSE
        )
      }
      
      if (anyNA(strata)) {
        stop("! `strata` must not contain missing values.", call. = FALSE)
      }
      
      # The stored Cox model contains the formula term `strata(strata_)`.
      # Therefore prediction newdata must contain the same high-level `strata_`
      # variable that was used at fitting. Do not pre-wrap ordinary vector or
      # factor strata with survival::strata(), because predict.coxph() will
      # evaluate strata(strata_) itself. Pre-wrapping would create levels such
      # as "1" instead of fit-time levels such as "strata_=1" and trigger
      # false "new levels" errors.
      #
      # Multiple matrix/data-frame strata columns are the only case where this
      # helper combines values first: the default/matrix interface stores one
      # synthetic `strata_` variable, so prediction must recreate that one
      # combined high-level variable before predict.coxph() evaluates
      # strata(strata_).
      newdata$strata_ <- if (is.matrix(strata) || is.data.frame(strata)) {
        do.call(
          survival::strata,
          c(as.list(as.data.frame(strata)), list(shortlabel = TRUE))
        )
      } else {
        strata
      }
    }
  }
  
  # step 9: add offset column if supplied
  if (!is.null(offset)) {
    if (!is.numeric(offset)) {
      stop("! offset must be numeric.", call. = FALSE)
    }
    
    if (length(offset) != nrow(newdata)) {
      stop("! offset must have one value per observation.", call. = FALSE)
    }
    
    newdata$offset_ <- offset
  }
  
  # step 10: return prepared prediction data
  newdata
}
#' Helper function to compute standard error of a partial predictor
#' 
#' To be used in \code{predict.mfp2()}.
#' 
#' @param model fitted `mfp2` object.
#' @param X transformed input matrix with variables of interest for partial predictor.
#' @param xref transformed reference value for variable of interest. Default is
#'  `NULL`, in which case this function computes standard errors without reference 
#' values.
#' @param include_intercept logical indicating whether intercept variance and
#' covariance should be included when no reference value is supplied. This must
#' match the `add_intercept` choice used for term predictions.
#' 
#' @details 
#' See pages 91-92 and following in the book by Royston and Sauerbrei 2008
#' for the formulas and mathematical details.
#' 
#' @return 
#' Numeric vector of conditional standard errors.
#' 
#' @references
#' Royston, P. and Sauerbrei, W., 2008. \emph{Multivariable Model - Building: 
#' A Pragmatic Approach to Regression Analysis based on Fractional Polynomials 
#' for Modelling Continuous Variables. John Wiley & Sons.}\cr
#' @keywords internal
#' @noRd
calculate_standard_error <- function(model, 
                                     X, 
                                     xref = NULL,
                                     include_intercept = TRUE) { 
  
  vcovx <- vcov(object = model)
  
  # this might happen if subset is used with very few observations
  if (any(is.nan(vcovx)))
    warning("i NaN detected in the covariance matrix of the model.",
            "i Standard errors for calculation of confidence intervals may not exist")
  
  # get rid of variance and covariance of intercept if any
  xnames <- colnames(X)
  ind <- match(xnames, colnames(vcovx))
  
  if (anyNA(ind)) {
    stop("Cannot compute SE; covariance matrix lacks columns: ",
         paste(xnames[is.na(ind)], collapse = ", "))
  }
  
  if (!is.null(xref)) {
    # Subtract the reference value: f(x)-f(xref)
    X <- sweep(X, 2, xref, "-")
  }
  
  # Augment X by the intercept only when the reported term value includes
  # the intercept. Contrasts and Cox partial predictors do not include intercept
  # variance.
  if (isTRUE(include_intercept) && model$family_string != "cox" && is.null(xref)) {
    intercept_ind <- match("(Intercept)", colnames(vcovx))
    if (is.na(intercept_ind)) {
      stop(
        "Cannot compute SE with intercept; covariance matrix lacks `(Intercept)`.",
        call. = FALSE
      )
    }
    X <- cbind("(Intercept)" = 1, X)
    ind <- c(intercept_ind, ind)
  }
  
  # the following computation is equivalent to the formula in the book but
  # uses matrix multiplications for efficiency
  vcovx <- vcovx[ind, ind, drop = FALSE]
  # similar to v = diag(x%*%vcovx%*%t(x))
  
  # Vectorised form of diag(X %*% vcovx %*% t(X)):
  v <- rowSums((X %*% vcovx) * X)
  
  sqrt(v)
}
