#' Predict Method for `mfp2` Fits
#' 
#' Obtains predictions from an `mfp2` object.
#' 
#' @details 
#' To prepare the `newdata` for prediction, this function applies any 
#' necessary shifting and scaling based on the factors obtained from the
#' training data. 
#' It is important to note that if the shifting factors are not sufficiently 
#' large as estimated from the training data, variables in `newdata` may end up
#' with negative values, which can cause prediction errors if non-linear 
#' functional forms are used. A warning is given in this case by the function.
#' The next step involves transforming the data using the selected
#' fractional polynomial (FP) power. If necessary, centering of variables is 
#' conducted. Once the transformation (and centering) is complete, the 
#' transformed data is passed to either `predict.glm()` or `predict.coxph()`, 
#' depending on the chosen family of models and when type is not 
#' `terms` and `contrasts`.
#'
#' @section Terms prediction:
#' This function allows to compute the partial linear predictors
#' for each variable selected into the final model if `type = "terms"`. Note 
#' that the results returned from this function are different from those of
#' `predict.glm()` and `predict.coxph()` since these functions do not take
#' into account that a single variable can be represented by multiple terms.
#' This functionality is useful to assess model fit, since it also allows to 
#' draw data points based on residuals. 
#' 
#' @section Contrasts: 
#' This functions allows to compute contrasts with reference to a specified
#' variable value if `type = "contrasts"`. In this case, the fitted partial
#' predictors will be centered at the reference value (i.e. 0), and also 
#' confidence intervals will have width 0 at that point. 
#' 
#' @param object a fitted object of class `mfp2`.
#' @param newdata optionally, a matrix with column names in which to look for 
#' variables with which to predict. See [mfp2()] for details.
#' @param type the type of prediction required.  The default is on the
#' scale of the linear predictors. See `predict.glm()` or `predict.coxph()` for
#' details. In case `type = "terms"`, see the Section on `Terms prediction`.
#' In case `type = "contrasts"`, see the Section on `Contrasts`.
#' @param terms a character vector of variable names specifying for which 
#' variables term or contrast predictions are desired.
#' Only used in case `type = "terms"` or `type = "contrasts"`.
#' If `NULL` (the default) then all selected variables in the final model will 
#' be used. In any case, only variables used in the final model are used, even
#' if more variable names are passed.
#' @param terms_seq a character string specifying how the range of variable 
#' values for term predictions are handled. The default `equidistant` computes
#' the range of the data range and generates an equidistant sequence of
#' 100 points from the minimum to the maximum values to properly show the functional 
#' form estimated in the final model. 
#' The option `data` uses the observed data values directly, but these may not 
#' adequately reflect the functional form of the data, especially when extreme
#' values or influential points are present.
#' @param alpha significance level used for computing confidence
#' intervals in terms prediction.
#' @param ref a named list of reference values used when `type = "contrasts"`.
#' Note that any variable requested in `terms`, but not having an entry in this
#' list (or if the entry is `NULL`) then the mean value (or minimum for binary
#' variables) will be used as reference. Values are specified on the original 
#' scale of the variable since the program will internally scale it using the
#' scaling factors obtained from [find_scale_factor()]. By default, this function
#' uses the means (for continuous variables) and minima (for binary variables) as
#' reference values. 
#' @param ... further arguments passed to `predict.glm()` or `predict.coxph()`.
#' 
#' @return 
#' For any `type` other than `"terms"` the output conforms to the output
#' of `predict.glm()` or `predict.coxph()`.
#' 
#' If `type = "terms"` or `type = "contrasts"`, then a named list with entries
#' for each variable requested in `terms` (excluding those not present in the
#' final model).
#' Each entry is a `data.frame` with the following columns:
#' 
#' * `variable`: variable values on original scale.
#' * `variable_pre`: variable with pre-transformation applied, i.e. shifted, 
#' scaled and centered as required.
#' * `value`: partial linear predictor or contrast (depending on `type`).
#' * `se`: standard error of partial linear predictor or contrast.
#' * `lower`: lower limit of confidence interval.
#' * `upper`: upper limit of confidence interval.
#' 
#' @seealso 
#' [mfp2()], [stats::predict.glm()], [survival::predict.coxph()]
#' 
#' @method predict mfp2
#' @export 
predict.mfp2 <- function(object, 
                         newdata = NULL, 
                         type = NULL,
                         terms = NULL,
                         terms_seq = c("equidistant", "data"),
                         alpha = 0.05,
                         ref = NULL, 
                         strata = NULL, 
                         offset = NULL, 
                         ...) {
  # check offset and strata
 # if(is.null(offset)){
 #  }
  
  # assert that the object must be mfp2
  if (!inherits(object, "mfp2")) 
    stop("The object is not an mfp2 object.", call. = FALSE)
  
  # set defaults and match arguments
  if (is.null(type)) 
    type <- ifelse(object$family_string == "cox", "lp", "link")
  
  if (is.null(terms))
    terms <- get_selected_variable_names(object)
  # TODO return error when terms given are not in the final model
  if (is.null(ref))  
    ref <- setNames(lapply(terms, function(v) NULL), terms)
  
  terms_seq <- match.arg(terms_seq)
  
  # TODO: add checks for missing strata and offset in case they were used in fit
  # TODO: add checks for correct specification of ref
  if (type == "contrasts" && length(ref) != sum(names(ref) != "", na.rm = TRUE))
    warning("i The supplied reference values (ref) must all be named.\n", 
            "i predict() continues but uses means (if variables are continous) or
            min (if binary) instead of the reference values.", call. = FALSE)
  
   if (type %in% c("terms", "contrasts")) {
    n_term1 <- length(terms)
    terms <- intersect(terms, get_selected_variable_names(object))
    # length of terms after intersections
    n_term2 <- length(terms)
    if (n_term2 == 0) {
      warning("i All the terms supplied are not in the final model.\n", 
              "i predict() continues but returns an empty list.", call. = FALSE) 
    } else if (n_term2 < n_term1)
      warning("i Some terms supplied are not in the final model.\n", 
              "i predict() continues but returns an empty list for those terms 
                 not in the model.", call. = FALSE)
    # TODO: return warning if the names(ref) != names(terms)
    if (!all(sapply(ref, is.null)) && any(!names(ref) %in% terms))
      warning("i Some of names of reference values are not in terms.\n", 
              "i predict() continues but does not consider them.", call. = FALSE)
 
    res_list <- list()
    for (t in terms) {
      
      # define sequence of variable data as named list
      # TODO: definition of sequence should be more flexible
      if (terms_seq == "equidistant") {
        x_range <- range(object$x_original[, t])
        x_seq  <- matrix(
          seq(x_range[1], x_range[2], length.out = 100),
          ncol = 1
        )
        colnames(x_seq) <- t
        
        # no need to apply pretransformation (shift and scaling), already done 
        # in x_original
        x_trafo <- as.matrix(prepare_newdata_for_predict(object, 
                                                         x_seq, 
                                                         apply_pre = FALSE))
      } else {
        # use data
        x_seq <- object$x_original[, t, drop = FALSE]
        x_names <- object$fp_powers[[t]]
        # in acd we might have a power and NA so we need to remove NA
        x_trafo <- object$x[, names(x_names[!is.na(x_names)]), drop = FALSE]
      }
      
      term_coef <- coef(object)[colnames(x_trafo)]
      
      # create output data.frame
    
      # intercepts do not play a role for Cox models, or contrasts
      intercept = coef(object)["(Intercept)"]
      if (is.na(intercept) || type == "contrasts") 
        intercept = 0
      
      res <- data.frame(
        # backtransform variable to original scale
        variable = (as.numeric(x_seq) * object$transformations[t,"scale"]) - 
          object$transformations[t,"shift"],
        variable_pre = as.numeric(x_seq),
        value = x_trafo %*% term_coef + intercept
      )
      
      x_ref_trafo <- NULL
      if (type == "contrasts") {
        # compute transformations for reference level
        # note that intercepts do not play a role here and that
        # (f(x) - f(x_ref)) * coef == f(x) * coef - f(x_ref) * coef
        
        x_ref <- ref[[t]]
        if (is.null(x_ref)) {
          v <- object$x_original[, t]
          if (length(unique(v)) == 2) {
            x_ref <- min(v, na.rm = TRUE)
          } else x_ref <- mean(v, na.rm = TRUE) 
        } else {
          # pretransform given reference level
          x_ref <- (x_ref + object$transformations[t,"shift"]) / object$transformations[t,"scale"]
        }
        # make sure it is a named matrix
        x_ref <- matrix(x_ref, nrow = 1, ncol = 1)
        colnames(x_ref) <- t
        
        # transform x_ref 
        x_ref_trafo <- as.matrix(prepare_newdata_for_predict(
          object, x_ref, apply_pre = FALSE, check_binary = FALSE))
        
        # compute contrasts, no intercepts necessary
        res$value <- res$value - as.numeric(x_ref_trafo %*% term_coef)
      }
      
      res$se <- calculate_standard_error(object, x_trafo, x_ref_trafo)
      mult <- qnorm(1 - (alpha / 2))
      res$lower <- res$value - mult * res$se
      res$upper <- res$value + mult * res$se
      
      res_list[[t]] <- res
    }
    names(res_list) <- terms
    
    return(res_list)
  } 
  
  # predict values
  # transform newdata using the FP powers from the training model
  if (!is.null(newdata)) {
    # assert that x has no missing data
    if (anyNA(newdata)) stop("! newdata must not contain any NA (missing data).\n", 
                       "i Please remove any missing data before passing newdata to this function.",
                       call. = FALSE)
    
    newdata <- prepare_newdata_for_predict(object, newdata)
    
    if (object$family_string == "cox") {
      return(getFromNamespace("predict.coxph", "survival")(
        object = object, newdata = newdata, type = type, ...
      ))
    } else {
      # TODO: offsets for glm?
      return(stats::predict.glm(
        object = object, newdata = newdata, type = type, ...
      ))
    } 
  }
  
  # no newdata supplied
  if (object$family_string == "cox") {
    getFromNamespace("predict.coxph", "survival")(
      object = object, type = type, ...
    )
  } else {
    stats::predict.glm(object = object, type = type, ...)
  }
}

#' Helper function to prepare newdata for predict function
#' 
#' To be used in [predict.mfp2()].
#' 
#' @param object fitted `mfp2` model object.
#' @param newdata dataset to be prepared for predictions. Its columns can be
#' a subset of the columns used for fitting the model. 
#' @param strata,offset passed from [predict.mfp2()].
#' @param apply_pre logical indicating wether the fitted pre-transformation
#' is applied or not.
#' @param apply_center logical indicating whether the fitted centers are applied
#' after transformation or not.
#' @param check_binary passed to [transform_vector_fp()].
prepare_newdata_for_predict <- function(object, 
                                        newdata, 
                                        strata = NULL, 
                                        offset = NULL, 
                                        apply_pre = TRUE, 
                                        apply_center = TRUE,
                                        check_binary = TRUE) {
  newdata <- as.matrix(newdata)
  
  # subset as appropriate
  vnames <- intersect(colnames(newdata), rownames(object$transformations))
  # sorting is not relevant as we always pass vnames
  newdata <- newdata[, vnames, drop = FALSE]
  
  if (apply_pre) {
    # step 1: shift and scale data using using shifting and scaling factors 
    # from the training data 
    # sort columns of newdata based on the names of shift/scale
    newdata <- sweep(newdata, 2, object$transformations[vnames, "shift"], "+")
    
    if (!all(newdata > 0)) 
      warning("i After shifting using training data some values in newdata remain negative.",
              "i Predictions for such observations may not be available in case of non-linear transformations.")
    
    newdata <- sweep(newdata, 2, object$transformations[vnames, "scale"], "/")
  }
  
  # step 2: transform the shifted and scaled data
  # do not center in this step
  newdata <- transform_matrix(
    newdata,
    power_list = object$fp_powers[vnames], 
    center = setNames(rep(FALSE, length(vnames)), vnames),
    keep_x_order = TRUE,
    acdx = setNames(object$fp_terms[vnames, "acd"], vnames),
    acd_parameter_list = object$acd_parameter[vnames],
    check_binary = check_binary
  )$x_transformed
  
  # step 3: center the transformed data
  if (apply_center && !is.null(object$centers)) {
    newdata <- center_matrix(newdata, object$centers[colnames(newdata)])
  }
  
  newdata <- data.frame(newdata)
  
  if (object$family_string == "cox") {
    # use of NextMethod here does not work as expected
    
    # add strata and offset as required
    if (!is.null(strata))
      newdata$strata_ <- survival::strata(strata, shortlabel = TRUE)
    
    if (!is.null(offset))
      newdata$offset_ <- offset
  } 
  
  newdata
}

#' Helper function to compute standard error of a partial predictor
#' 
#' To be used in [predict.mfp2()].
#' 
#' @param model fitted `mfp2` object.
#' @param X transformed input matrix with variables of interest for partial predictor.
#' @param xref transformed reference value for variable of interest. Default is
#'  `NULL`, in which case this function computes standard errors without reference 
#' values.
#' 
#' @details 
#' See pages 91-92 and following in the book by Royston and Sauerbrei 2008
#' for the formulas and mathematical details.
#' 
#' @return 
#' Standard error.
#' 
#' @references
#' Royston, P. and Sauerbrei, W., 2008. \emph{Multivariable Model - Building: 
#' A Pragmatic Approach to Regression Anaylsis based on Fractional Polynomials 
#' for Modelling Continuous Variables. John Wiley & Sons.}\cr
calculate_standard_error <- function(model, 
                                     X, 
                                     xref = NULL) { 

  vcovx <- vcov(object = model)
  
  # this might happen if subset is used with very few observations
  if (any(is.nan(vcovx)))
    warning("i NaN detected in the covariance matrix of the model.",
            "i Standard errors for calculation of confident intervals may not exist")
  
  # get rid of variance and covariance of intercept if any
  xnames <- colnames(X)
  ind <- match(xnames, colnames(vcovx))
  
  if (!is.null(xref)) {
    # Subtract the reference value: f(x)-f(xref)
    X <- sweep(X, 2, xref, "-")
  }

  # augment X by intercept if necessary
  # i.e. when a cox model is used or a reference value is given we don't use
  # the variance and covariances of the intercept
  if (model$family_string != "cox" && is.null(xref)) {
    X <- cbind(1, X)
    ind <- c(1, ind)
  }

  # the following computation is equivalent to the formula in the book but
  # uses matrix multiplications for efficiency
  vcovx <- vcovx[ind, ind, drop = FALSE]
  # similar to v = diag(x%*%vcovx%*%t(x))
  v <- sapply(
    1:nrow(X), 
    function(i, x, vcovx) x[i, , drop = FALSE] %*% vcovx %*% t(x[i, , drop = FALSE]),
    x = X, vcovx = vcovx
  )
  
  sqrt(v)
}
