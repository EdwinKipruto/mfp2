#' Predict Method for `mfpa` Fits
#' 
#' Obtains predictions from an `mfpa` object
#' 
#'@details 
#' The fp transformation generally transforms `x` as follows. For each pi in
#' 
#'@section Data processing: 
#' An important note on data processing. Variables are shifted and scaled 
#' before being transformed by any powers. That is to ensure positive values
#' and reasonable scales. Note that scaling does not change the estimated 
#' powers, see also [find_scale_factor()].
#' 
#' However, they are centered after transformation. 
#' That is to ensure that the correlation between variables stay intact, 
#' as centering before transformation would affect them. This is described
#' in Sauerbrei et al (2006), as well as in the Stata manual of `mfp`.
#' 
#'@param object a fitted object of class inheriting from \code{"glm", "lm", "coxph"}.
#'@param newdata optionally, a data frame in which to look for variables with
#' which to predict. See predict.glm() or predict.coxph() for details
#'@param type the type of prediction required.  The default is on the
#' scale of the linear predictors. See predict.glm() or predict.coxph() for details 
#'@param se.fit NOT NEEDED, TO DISCUSS WITH MICHAEL
#'@param dispersion NOT NEEDED 
#'@param terms to be discussed 
#'@param na.action function determining what should be done with missing values in newdata. The default is to predict NA. 
#'@param collapse optional vector of subject identifiers. If specified, the output will contain one entry per subject rather than one entry per observation.
#'@param reference reference for centering predictions, see details below
#'@param ... further arguments passed to or from other methods.
#'@method predict mfpa
#'@export 
predict.mfpa = function(object, newdata=NULL, shift = NULL, scale = NULL, type = c("link", "response","terms", "lp", "risk", "expected", "survival"),
                         terms = NULL, se.fit = FALSE, 
                         dispersion = NULL, na.action = na.pass, collapse,
                         reference=c("strata", "sample", "zero"),...) {
  print("hello")
  type <- match.arg(type)
  #fam = object$family$family
  #terms = object$terms
  # Names of coefficients in the model used to sort newdata
  # if(fam=="cox"){
  #   xnames = names(object$coefficients)
  # }else{
  #   # Get rid of intercept
  #   xnames = names(object$coefficients)[-1]
  # }
  # Transform newdata using the FP powers from the training model
  if(!is.null(newdata)){
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
    # shift and scale newdata. It is recommended that the user shift their data or supply shift
    # and scaling factor
    xx = apply(newdata, 2, function(x) apply_shift_scale(x, shift = shift, scale = scale))
    # Apply transformating to the newdata
    newdata = data.frame(transform_matrix(xx,power_list = powers_list, center= center_varx,
                               keep_x_order = T,acdx = acd_varx))
    #newdata = newdata[, xnames, drop = FALSE]
  }else{
    # Use already transformed data if newdata is not supplied. 
    newdata = data.frame(object$x)
  }
  #newdata
  # Make predictions
  fam = object$family$family
  if(fam=="cox"){
    predict.coxph(object = object, newdata = newdata, type = type,se.fit = se.fit,terms = terms,collapse = collapse,
                  reference = reference,na.action = na.action)
  }else{
    predict.glm(object = object, newdata = newdata, type = type,se.fit = se.fit,dispersion = dispersion, terms = terms, na.action = na.action)
}
}