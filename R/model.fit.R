#' A function that fits generalized linear models and cox proportional hazard
#' models. It uses glm.fit() and coxph.fit() function in stats and survival package
#'
#' @param x matrix of predictors including column of 1s for intercept except for cox
#' @param y  a vector of outcome variable in glm but a matrix of survival object in cox
#' @param method a character string specifying the method for tie handling. See coxph()
#' @param family a character names specifying glm family in addition to cox as a family
#' @param strata,control,weights,offset,rownames,nocenter parameters for cox model. see coxph()
#' @importFrom stats  family
#@export
model.fit <- function(x, y,family,weights, offset, method, strata,control, rownames,
                      nocenter){
  if(family=="cox"){
    # cox needs more work especially on how to handle strata
    fit <- coxfit(x = x, y = y, strata = strata, weights = weights, offset = offset,
                  control = control,method = method, rownames = rownames,
                  nocenter=nocenter)
  }else{
    # glm.fit() does not allow character as a family- only functions such as gaussian() etc.
    if(is.character(family))
      family <- get(family, mode = "function", envir = parent.frame())
    if(is.function(family)) family <- family()
    # add 1s for intercept
    xx <- cbind(rep(1,length(y)),x)
    # Note that x can be NULL when fitting a NULL model (using only adjustment variable)
    # when all FP powers estimated are NA-all variables removed
    colnames(xx) = if(is.null(ncol(x))){"Intercept"}else{c("Intercept", colnames(x))}
    fit<-glmfit(y = y, x = xx, family = family,weights = weights, offset = offset)
  }
  return(fit)
}
