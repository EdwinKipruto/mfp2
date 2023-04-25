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
#'@export predict.mfpa
predict.mfpa = function(object, newdata=NULL, shift = NULL, scale = NULL, type = c("link", "response","terms", "lp", "risk", "expected", "survival"),
                         terms = NULL, se.fit = FALSE, 
                         dispersion = NULL, na.action = na.pass, collapse,
                         reference=c("strata", "sample", "zero"),...) {
  print("hello")
  type <- match.arg(type)
  reference <- match.arg(reference)
  
  # Transform newdata using the FP powers from the training model
  if(!is.null(newdata)){
    ## STEP 1: SHIFT AND SCALE NEW DATA
    # We use shift and scaling factors from the training data
    # and scaling factor
    if(is.null(shift)) shift = object$transformations[,"shift"]
    if(is.null(scale)) scale = object$transformations[,"scale"]
    #sort newdata based on the names of shift/scale
    newdata = newdata[, rownames(object$transformations), drop=FALSE]
    newdata_shifted_scaled = apply(newdata, 2, function(x) apply_shift_scale(x, shift = shift, scale = scale))
    ## STEP 2: TRANSFORM THE SHIFTED AND SCALED DATA
    # Extract the FP powers for all variables
    powers_vars= object$fp_terms
    # The first 6 columns are not important, only on power terms are needed
    powers_mat=powers_vars[,-c(1:6)]
    # transpose the powers and make a list
    powers_mat = t(powers_mat)
    powers_list=split(powers_mat, rep(1:ncol(powers_mat), each = nrow(powers_mat)))
    # Assign names to the powers
    names(powers_list) = colnames(powers_mat)
    # Set names to center and acd required by transform_matix()
    center_varx = setNames(object$transformations[,"center"], rownames(object$transformations))
    acd_varx = setNames(powers_vars[,"acd"],rownames(powers_vars))
    # Apply transformating to the newdata
    newdata_transformed = data.frame(transform_matrix(newdata_shifted_scaled,power_list = powers_list, center= center_varx,
                               keep_x_order = T,acdx = acd_varx))
    #newdata = newdata[, xnames, drop = FALSE]
  }else{
    # Use already transformed data if newdata is not supplied. 
    newdata_transformed = data.frame(object$x)
  }
  #newdata
  # Make predictions
  fam = object$family$family
  if(fam=="cox"){
    predict.coxph(object = object, newdata = newdata_transformed, type = type,se.fit = se.fit,terms = terms,collapse = collapse,
                  reference = reference,na.action = na.action)
  }else{
    predict.glm(object = object, newdata = newdata_transformed, type = type,se.fit = se.fit,dispersion = dispersion, terms = terms, na.action = na.action)
}
}