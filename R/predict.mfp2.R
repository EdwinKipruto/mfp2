#' Predict Method for `mfp2`
#' 
#' Obtains predictions from an `mfp2` object.
#' 
#' @details 
#' To prepare the `newdata` for prediction, this function applies any 
#' necessary shifting based on factors obtained from the training data. 
#' It is important to note that if the shifting factors estimated from the 
#' training data are not sufficiently large, variables in `newdata` may end up
#' being non-positive, which can cause prediction errors when non-linear 
#' functional forms such as logarithms are used. In such cases, the function
#' issues a warning. 
#' The next step involves transforming the data using the selected
#' fractional polynomial (FP) powers. After transformation, variables are
#' centered if `center` was set to TRUE in `mfp2()`. Once transformation 
#' (and centering) is complete, the transformed data is passed 
#' to either `predict.glm()` or `predict.coxph()`, depending on the model 
#' family used, provided that `type` is neither `terms` nor `contrasts` (see the
#' section handling `terms` and `contrasts` for details).
#'
#' @section Terms prediction:
#' If `type = "terms"`, this function computes the partial linear predictors 
#' for each variable included in the final model. Unlike `predict.glm()` and 
#' `predict.coxph()`, this function accounts for the fact that a single variable 
#' may be represented by multiple transformed terms.
#'
#' For a variable modeled using a first-degree fractional polynomial (FP1), 
#' the partial predictor is given by \eqn{\hat{\eta}_j = \hat{\beta}_0 + x_j^* \hat{\beta}_j}, where 
#' \eqn{x_j^*} is the transformed variable (centered if `center = TRUE`).
#'
#' For a second-degree fractional polynomial (FP2), the partial predictor 
#' takes the form \eqn{\hat{\eta}_j = \hat{\beta}_0 + x_{j1}^* \hat{\beta}_{j1} + x_{j2}^* \hat{\beta}_{j2}}, 
#' where \eqn{x_{j1}^*} and \eqn{x_{j2}^*} are the two transformed components 
#' of the original variable (again, centered if `center = TRUE`).
#'
#' This functionality is particularly useful for visualizing the functional
#' relationship of a continuous variable, or for assessing model fit when
#' residuals are included. See also `fracplot()`.
#' 
#' @section Contrasts:
#' If `type = "contrasts"`, this function computes contrasts relative to a 
#' specified reference value for the \eqn{j}th variable (e.g., age = 50). Let 
#' \eqn{x_j} denote the values of the \eqn{j}th variable in `newdata`, and 
#' \eqn{x_j^{\text{ref}}} the reference value. The contrast is defined as the 
#' difference between the partial linear predictor evaluated at the transformed 
#' (and centered, if `center = TRUE`) value \eqn{x_j}, and that evaluated at the
#' transformed reference value \eqn{x_j^{(\text{ref}})}, 
#' i.e., \eqn{f(x_j^*) - f(x_j^{*(\text{ref})})}.
#'
#' For a first-degree fractional polynomial (FP1), the partial predictor is:
#' \deqn{\hat{f}(x_j^*) = \hat{\beta}_0 + x_j^* \hat{\beta}_j}
#' and the contrast is:
#' \deqn{\hat{f}(x_j^*) - \hat{f}(x_j^{*(\text{ref})}) = x_j^* \hat{\beta}_j - x_j^{*(\text{ref})} \hat{\beta}_j}
#'
#' For a second-degree fractional polynomial (FP2), the partial predictor is:
#' \deqn{\hat{f}(x_j^*) = \hat{\beta}_0 + x_{j1}^* \hat{\beta}_{j1} + x_{j2}^* \hat{\beta}_{j2}}
#' and the contrast is:
#' \deqn{\hat{f}(x_j^*) - \hat{f}(x_j^{*(\text{ref})}) = x_{j1}^* \hat{\beta}_{j1} + x_{j2}^* \hat{\beta}_{j2} - x_{j1}^{*(\text{ref})} \hat{\beta}_{j1} - x_{j2}^{*(\text{ref})} \hat{\beta}_{j2}}
#'
#' where \eqn{x_j^*}, \eqn{x_{j1}^*}, and \eqn{x_{j2}^*} are the transformed 
#' (and centered, if applicable) components of the \eqn{j}th variable, and 
#' the \eqn{\hat{\beta}} terms are the corresponding model estimates
#'
#' The reference value \eqn{x_j^{(\text{ref})}} is first **shifted** using the 
#' same shifting factor estimated from the training data, then transformed using 
#' the estimated fractional polynomial (FP) powers, and finally **centered** 
#' (if `center = TRUE`) using the **mean of the transformed (and shifted) values 
#' of \eqn{x_j} in the training data**—ensuring full consistency with the fitted model.
#'
#' If `ref = NULL`, the function uses the **mean of the shifted \eqn{x_j}** as the 
#' reference value when \eqn{x_j} is continuous, or the **minimum of \eqn{x_j}** 
#' (typically 0) when \eqn{x_j} is binary. This provides a natural and interpretable 
#' baseline in the absence of a user-specified reference.
#'
#' The fitted partial predictors are centered at the reference point, 
#' meaning the contrast at that point is zero. Correspondingly, confidence intervals 
#' at the reference value have zero width, reflecting no contrast.
#'
#' This functionality is especially useful for comparing the effect of a variable 
#' relative to a meaningful baseline, such as clinically relevant value.
#' 
#' @param object a fitted object of class `mfp2`.
#' @param newdata optionally, a matrix with column names in which to look for 
#' variables with which to predict. If provided, the variables are internally 
#' shifted using the shifting values stored in `object`. See [mfp2()] for 
#' further details.
#' @param type the type of prediction required.  The default is on the scale of
#' the linear predictors. See `predict.glm()` or `predict.coxph()` for details. 
#' In case `type = "terms"`, see the Section on `Terms prediction`. In case 
#' `type = "contrasts"`, see the Section on `Contrasts`.
#' @param terms a character vector of variable names specifying for which 
#' variables term or contrast predictions are desired. Only used in case 
#' `type = "terms"` or `type = "contrasts"`. If `NULL` (the default) then all 
#' selected variables in the final model will be used. In any case, only 
#' variables used in the final model are used, even if more variable names are 
#' passed.
#' @param terms_seq a character string specifying how the range of variable 
#' values for term predictions are handled. The default `equidistant` computes
#' the range of the data range and generates an equidistant sequence of
#' 100 points from the minimum to the maximum values of shifted values to 
#' properly show the functional form estimated in the final model. 
#' The option `data` uses the observed data values directly, but these may not 
#' adequately reflect the functional form of the data, especially when extreme
#' values or influential points are present.
#' @param alpha significance level used for computing confidence intervals in 
#' terms prediction.
#' @param ref a named list of reference values used when `type = "contrasts"`.
#' Note that any variable requested in `terms`, but not having an entry in this
#' list (or if the entry is `NULL`) then the mean value of shifted data 
#' (or minimum for binary variables) will be used as reference. Values should be 
#' specified on the original scale of the variable since the program will 
#' internally scale it using the scaling factors obtained from 
#' [find_scale_factor()]. By default, this function uses the means 
#' (for continuous variables) and minimum (for binary variables) as
#' reference values. 
#' @param strata stratum levels used for predictions. 
#' @param newoffset A vector of offsets used for predictions. This parameter is 
#' important when newdata is supplied. The offsets are directly added to the 
#' linear predictor without any transformations.
#' @param nseq Integer specifying how many values to generate when 
#' `terms_seq = "equidistant"`. Default is 100.
#' @param ... further arguments passed to `predict.glm()` or `predict.coxph()`.
#' 
#' @examples
#'
#' # Gaussian model
#' data("prostate")
#' x = as.matrix(prostate[,2:8])
#' y = as.numeric(prostate$lpsa)
#' # default interface
#' fit1 = mfp2(x, y, verbose = FALSE)
#' predict(fit1) # make predictions
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
#' * `variable`: variable values on original scale (without shifting).
#' * `variable_pre`: variable with pre-transformation applied, i.e. shifted, and centered as required.
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
                         newoffset = NULL, 
                         nseq = 100,
                         ...) {
 
  terms_seq <- match.arg(terms_seq)
  
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
  
  # TODO: add checks for missing strata and offset in case they were used in fit
  if (type == "contrasts" && length(ref) != sum(names(ref) != "", na.rm = TRUE)) {
    warning("i The supplied reference values (ref) must all be named.\n", 
            "i predict() continues but uses means (if variables are continous) or
            min (if binary) instead of the reference values.", call. = FALSE)
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
      warning("i Some terms supplied are not in the final model.\n", 
              "i predict() continues but returns an empty list for those terms 
                 not in the model.", call. = FALSE)
    # return warning if the names(ref) != names(terms)
    if (!all(sapply(ref, is.null)) && any(!names(ref) %in% terms))
      warning("i Some of names of reference values are not in terms.\n", 
              "i predict() continues but does not consider them.", call. = FALSE)
   # 
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
        # in x_original
        x_trafo <- as.matrix(prepare_newdata_for_predict(object, 
                                                         x_seq, 
                                                         apply_pre = FALSE))
      } else {
        # use data
        if (!is.null(newdata)) {
          x_seq <- newdata[, t, drop = FALSE] + object$transformations[t,"shift"]
        } else {
          x_seq <- object$x_original[, t, drop = FALSE] # already shifted
        }
        x_names <- object$fp_powers[[t]]
        # in acd we might have a power and NA so we need to remove NA
        #x_trafo <- object$x[, names(x_names[!is.na(x_names)]), drop = FALSE]
        x_trafo <- as.matrix(prepare_newdata_for_predict(object, 
                                                         x_seq, 
                                                         apply_pre = FALSE))
      }
      # 
      term_coef <- coef(object)[colnames(x_trafo)]
      
      # create output data.frame
    
      # Check if the intercept is actually in the model. Intercepts is zero for 
      # Cox model
      intercept <- ifelse("(Intercept)" %in% names(coef(object)), 
                           coef(object)["(Intercept)"],
                           0)
      # Intercept do not play a role for type = "contrasts"
      if (type == "contrasts") {
        intercept <- 0
      }
      
      # 
      res <- data.frame(
        # backtransform variable to original scale
        variable = (as.numeric(x_seq)) - object$transformations[t,"shift"],
        variable_pre = as.numeric(x_seq),
        value = x_trafo %*% term_coef + intercept
      )
      
      #
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
  
  # predicted values
  # transform newdata using the FP powers from the training model
  if (!is.null(newdata)) {
    # assert that x has no missing data
    if (anyNA(newdata)) {
      stop("! newdata must not contain any NA (missing data).\n", 
                       "i Please remove any missing data before passing newdata to this function.",
                       call. = FALSE)
    }
    
    newdata <- prepare_newdata_for_predict(object, newdata, check_binary = FALSE)
    betas <- object$coefficients
  
    # check whether offset was used in the model
    is_offset <- all(object$offset == 0)
    
    if (object$family_string == "cox") {
      # subset newdata based on names of coefficients in the model
      newdata <- as.matrix(newdata[,names(betas),drop = FALSE])
      nfit <- newdata %*% betas
    } else {
      # remove intercept before subsetting 
      newdata <- as.matrix(newdata[,names(betas)[-1],drop = FALSE])
      nfit <- cbind(1,newdata) %*% betas

    } 
  
    if (!is_offset) {
      if (is.null(newoffset)) {
        stop("No newoffset provided for prediction, yet offset was used in mfp2", call. = FALSE)
      }
       nr <- dim(newdata)[1]
       if (nr != length(newoffset)) {
         stop("The length of newoffset must be equal to the number of rows of newdata", call. = FALSE)
       }
      nfit <- nfit + newoffset
    }
    
    if (object$family_string == "gaussian")
      return(nfit)
    
    # response = risk score in coxph
    return(switch(type,
                   response = exp(nfit),
                   risk = exp(nfit),
                   nfit
    ))
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
#' @param apply_pre logical indicating whether the fitted pre-transformation
#' is applied or not.
#' @param apply_center logical indicating whether the fitted centers are applied
#' after transformation or not.
#' @param check_binary passed to [transform_vector_fp()].
#' 
#' @return A dataframe of transformed newdata
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
    
  # exclude binary variables before returning an error
  index <-  unname(which(apply(newdata, 2, function(x) length(unique(x))) > 2))
    
    if (!all(newdata[,index] > 0)) 
      warning("i After shifting using training data some values in newdata remain negative.",
              "i Predictions for such observations may not be available in case of non-linear transformations.")
    # already the coefficients are in original scale after shifting
    #newdata <- sweep(newdata, 2, object$transformations[vnames, "scale"], "/")
  }
  
  # step 2: transform the shifted data
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
