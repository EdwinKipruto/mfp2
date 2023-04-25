predict.mfpa <- function(object, newdata=NULL, shift = NULL, scale = NULL, type = c("link", "response","terms", "lp", "risk", "expected", "survival"),
                         terms = NULL, se.fit = FALSE, 
                         dispersion = NULL, na.action = na.pass, collapse,
                         reference=c("strata", "sample", "zero"),...) {
  type <- match.arg(type)
  fam = object$family$family
  terms = object$terms
  # Names of coefficients in the model used to sort newdata
  if(fam=="cox"){
    xnames = names(object$coefficients)
  }else{
    # Get rid of intercept
    xnames = names(object$coefficients)[-1]
  }
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
    newdata = transform_matrix(xx,power_list = powers_list, center= center_varx,
                               keep_x_order = T,acdx = acd_varx)
    newdata = newdata[, xnames, drop = FALSE]
  }else{
    # Use already transformed data if newdata is not supplied. 
    newdata = out$X[,xnames,drop = FALSE]
  }
  newdata
  # Make predictions
  predict(object = object, newdata = newdata, type = type,se.fit = se.fit,dispersion = dispersion, terms = terms, na.action = na.action)
  
}