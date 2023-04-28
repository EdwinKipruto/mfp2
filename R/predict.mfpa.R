#' Predict Method for `mfpa` Fits
#' 
#' Obtains predictions from an `mfpa` object.
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
#' depending on the chosen family of models.
#'
#' @section Terms prediction:
#' This function allows to compute the partial linear predictors, or contrasts,
#' for each variable selected into the final model if `type = "terms"`. Note 
#' that the results returned from this function are different from those of
#' `predict.glm()` and `predict.coxph()` since these functions do not take
#' into account that a single variable can be represented by multiple terms.
#' 
#' @param object a fitted object of class `mfpa`.
#' @param newdata optionally, a matrix with column names in which to look for 
#' variables with which to predict. See [mfpa()] for details.
#' @param type the type of prediction required.  The default is on the
#' scale of the linear predictors. See `predict.glm()` or `predict.coxph()` for
#' details. In case `type = "terms"`, see the Section on `Terms prediction`.
#' @param terms a character vector of variable names specifying for which 
#' variables term predictions are desired. Only used in case `type = "terms"`.
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
#' @param terms_alpha significance level used for computing confidence
#' intervals in terms prediction.
#' @param ... further arguments passed to `predict.glm()` or `predict.coxph()`.
#' 
#' @return 
#' For any `type` other than `"terms"` the output conforms to the output
#' of `predict.glm()` or `predict.coxph()`.
#' 
#' If `type = "terms"`, then a named list with entries for each variable
#' requested in `terms` (excluding those not present in the final model).
#' Each entry is a `data.frame` with the following columns:
#' 
#' * `variable`: variable values (shifted, scaled and centered as required).
#' * `contrast`: partial linear predictor.
#' * `se`: standard error of partial linear predictor.
#' * `lower`: lower limit of confidence interval.
#' * `upper`: upper limit of confidence interval.
#' 
#' @seealso 
#' [mfpa()], [stats::predict.glm()], [survival::predict.coxph()]
#' 
#' @method predict mfpa
#' @export 
predict.mfpa <- function(object, 
                         newdata=NULL, 
                         type = NULL,
                         terms = NULL,
                         terms_seq = c("equidistant", "data"),
                         terms_alpha = 0.05,
                         strata = NULL, 
                         offset = NULL, 
                         ...) {
  
  # set defaults and match arguments
  if (is.null(type)) 
    type <- ifelse(object$family_string == "cox", "lp", "link")
  
  if (is.null(terms))
    terms <- get_selected_variable_names(object)
  
  terms_seq <- match.arg(terms_seq)
  
  # TODO: add checks for missing strata and offset in case they were used in fit
  
  if (type == "terms") {
    terms <- intersect(terms, get_selected_variable_names(object))
    
    res_terms <- list()
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
      } else {
        # use data
        x_seq <- object$x_original[, t, drop = FALSE]
      }
      
      # transform and center if required
      if (terms_seq != "data") {
        # center after transformation using fitted centers
        x_trafo <- transform_matrix(
          x_seq, 
          power_list = object$fp_powers[t], 
          acdx = setNames(object$fp_terms[t, "acd"], t), 
          center = setNames(c(FALSE), t)
        )$x_transformed
        
        if (!is.null(object$centers))
          x_trafo <- center_matrix(x_trafo,
                                   object$centers[colnames(x_trafo)])
      } else {
        # use data directly
        x_trafo <- object$x[, names(object$fp_powers[[t]]), drop = FALSE]
      }
      
      term_coef <- coef(object)[colnames(x_trafo)]
      
      # create output data.frame
      intercept = coef(object)["(Intercept)"]
      if (is.na(intercept)) 
        intercept = 0
      
      res <- data.frame(
        variable = as.numeric(x_seq),
        contrast = x_trafo %*% term_coef + intercept
      )
      res$se <- calculate_standard_error(object, x_trafo)
      mult <- qnorm(1 - (terms_alpha / 2))
      res$lower <- res$contrast - mult * res$se
      res$upper <- res$contrast + mult * res$se
      
      res_terms[[t]] <- res
    }
    names(res_terms) <- terms
    
    return(res_terms)
  }
  
  # transform newdata using the FP powers from the training model
  if (!is.null(newdata)) {
    # step 1: shift and scale data using using shifting and scaling factors 
    # from the training data
    # 
    # subset and sort columns of newdata based on the names of shift/scale
    newdata <- newdata[, rownames(object$transformations), drop = FALSE]
    newdata <- sweep(newdata, 2, object$transformations[,"shift"], "+")
    
    if (!all(newdata > 0)) 
      warning("i After shifting using training data some values in newdata remain negative.",
              "i Predictions for such observations may not be available in case of non-linear transformations.")
    
    newdata <- sweep(newdata, 2, object$transformations[,"scale"], "/")
    
    # step 2: transform the shifted and scaled data
    # do not center in this step
    newdata_transformed <- transform_matrix(
      newdata,
      power_list = object$fp_powers, 
      center = setNames(rep(FALSE, nrow(object$transformations)), 
                        rownames(object$transformations)),
      keep_x_order = TRUE,
      acdx = setNames(object$fp_terms[,"acd"], 
                      rownames(object$fp_terms))
    )
    
    # step 3: center the transformed data
    if (!is.null(object$centers)) {
      newdata <- center_matrix(
        newdata_transformed$x_transformed,
        object$centers
      )
    }
    
    newdata <- data.frame(newdata)
    
    if (object$family_string == "cox") {
      # use of NextMethod here does not work as expected
      
      # add strata and offset as required
      if (!is.null(strata))
        newdata$strata_ = survival::strata(strata, shortlabel = TRUE)
      
      if (!is.null(offset))
        newdata$offset_ = offset
      
      # use getFromNamespace here, as predict.coxph is not exported
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
