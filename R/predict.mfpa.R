#' Predict Method for `mfpa` Fits
#' 
#' Obtains predictions from an `mfpa` object
#' 
#'@details 
#'To prepare the `newdata` for prediction, the first step is to apply any 
#'necessary shifting and scaling based on the factors obtained from the training data. 
#'It's important to note that if the shifting factors are not sufficiently large,
#'variables may end up with negative values, which can cause prediction errors.
#'The next step involves transforming the data using the selected
#'fractional polynomial (FP) power. Once the transformation is complete, the 
#'transformed data is passed to either `predict.glm()` or `predict.coxph()`, 
#'depending on the chosen family of models.
#' 
#'@param object a fitted object of class `mfpa` which inherit from `glm`,`lm` or `coxph`.
#'@param newdata optionally, a data frame in which to look for variables with
#' which to predict. See `predict.glm()` or `predict.coxph()` for details
#'@param type the type of prediction required.  The default is on the
#' scale of the linear predictors. See `predict.glm()` or `predict.coxph()` for details 
#'@param se.fit NOT NEEDED, TO DISCUSS WITH MICHAEL
#'@param dispersion NOT NEEDED 
#'@param terms to be discussed 
#'@param na.action function determining what should be done with missing values in newdata. The default is to predict NA. 
#'@param collapse optional vector of subject identifiers. If specified, the output will contain one entry per subject rather than one entry per observation.
#'@param reference reference for centering predictions, see details below
#'@param ... further arguments passed to or from other methods.
#'@method predict mfpa
#'@export 
predict.mfpa <- function(object, 
                         newdata=NULL, 
                         type = NULL,
                         terms = NULL, 
                         strata = NULL, 
                         offset = NULL, 
                         ...) {
  print("hello")
  if (is.null(type)) type <- ifelse(object$family_string == "cox", "lp", "link")
  
  # TODO: add checks for missing strata and offset in case they were used in fit
  
  # Transform newdata using the FP powers from the training model
  if (!is.null(newdata)) {
    ## STEP 1: shift and scale data using using shifting and scaling factors from the training data
     shift <- object$transformations[,"shift"]
     scale <- object$transformations[,"scale"]
    #subset and sort columns of newdata based on the names of shift/scale
    newdata <- newdata[, rownames(object$transformations), drop = FALSE]
    # shift and scale
    newdata_shifted <- sweep(newdata, 2, shift, "+")
    if (!all(newdata_shifted > 0)) 
      warning("i After shifting using training data some values in newdata remain negative.",
              "i Predictions for such observations may not be available in case of non-linear transformations.")
    
    newdata_shifted_scaled <- sweep(newdata_shifted, 2, scale, "/")
    
    ## STEP 2: TRANSFORM THE SHIFTED AND SCALED DATA
    # Set names to center and acd required by transform_matix()
    center_varx <- setNames(rep(FALSE, nrow(object$transformations)), 
                            rownames(object$transformations))
    acd_varx <- setNames(object$fp_terms[,"acd"], 
                         rownames(object$fp_terms))
    
    # Apply transformating to the newdata
    newdata_transformed <- transform_matrix(
        newdata_shifted_scaled,
        power_list = object$fp_powers, 
        center = center_varx,
        keep_x_order = T,
        acdx = acd_varx
    )
    
    ## STEP 3: CENTER THE TRANSFORMED DATA
    if (!is.null(object$centers)) {
      newdata_transformed <- center_matrix(
        newdata_transformed$x_transformed,
        object$centers
      )
    }
    
    newdata_transformed <- data.frame(newdata_transformed)
    
    if (object$family_string == "cox") {
      # use of NextMethod here does not work as expected
      
      # add strata and offset as required
      if (!is.null(strata))
        newdata_transformed$strata_ = survival::strata(strata, shortlabel = TRUE)
      
      if (!is.null(offset))
        newdata_transformed$offset_ = offset
      
      return(getFromNamespace("predict.coxph", "survival")(object = object, newdata = newdata_transformed, type = type, ...))
    } else {
      # TODO: offsets for glm?
      return(stats::predict.glm(object = object, newdata = newdata_transformed, type = type, ...))
    } 
  }
  
  # no newdata supplied
  if (object$family_string == "cox") {
    getFromNamespace("predict.coxph", "survival")(object = object, type = type, ...)
  } else {
    stats::predict.glm(object = object, type = type, ...)
  }
}
