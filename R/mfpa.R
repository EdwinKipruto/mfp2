#' Multivariable fractional polynomial models with extensions of sigmoid 
#' functions
#'
#' Selects the multivariable fractional polynomial (FP) model that best predicts
#' the outcome variable. It also has the ability to model a sigmoid relationship
#' between x and an outcome variable using the approximate cumulative 
#' distribution (ACD) transformation proposed by Royston (2016).
#'
#' @section Details on `family` option:
#'
#' `mfpa()` supports the family object as used by [stats:glm()]. The built in
#' families are specifed via a character string. `mfpa(..., family="binomial")` 
#' fits a logistic regression model, while `mfpa(..., family="gaussian")`
#' fits a linear regression (ordinary least squares) model.
#'
#' For Cox models, the response should preferably be a `Surv` object,
#' created by the [survival::Surv()] function, and the `family = "cox"`. 
#' Only right-censored data are currently supported. To fit stratified Cox
#' models, the strata option can be used. Currently `mfpa()` only supports
#' a single factor as strata.
#'
#' @section Details on scaling and centering:
#'
#' After adjusting the location of \eqn{x} so that its minimum value is positive,
#' creating \eqn{x'}, automatic scaling will divide each value of \eqn{x'} by 
#' \eqn{10^p} where the exponent \eqn{p} is given by 
#' \deqn{p = sign(k) \times floor(|k|) \quad \text{where} \quad k = log_{10} (max(x')- min(x'))}
#'
#' The FP transformation of \eqn{x'} is centered on the mean of the observed 
#' values of \eqn{x'}. For example, for the FP1 model \eqn{\beta_0 + \beta_1x^p},
#' the actual model fitted by the software would be 
#' \eqn{\beta'_0 + \beta'_1(x'^p-mean(x'^p))}. This approach ensures that
#' the revised constant \eqn{\beta'_0} equals the fitted value of the FP
#' function at the mean of \eqn{x'}.
#'
#' @section Details on  approximate cumulative distribution transformation:
#' The approximate cumulative distribution (ACD) transformation (Royston 2014a) 
#' converts each predictor, \eqn{x}, smoothly to an approximation, \eqn{acd(x)}, 
#' of its empirical cumulative distribution function. 
#' This is done by smoothing a probit transformation of 
#' the scaled ranks of \eqn{x}. \eqn{acd(x)} could be used instead of \eqn{x} 
#' as a covariate. This has the advantage of providing sigmoid curves, something
#' that regular FP functions cannot achieve. 
#' Details of the precise definition and some possible uses of the ACD 
#' transformation in a univariate context are given by Royston (2014a). 
#' Royston (2014b) describes how one could go further and replace FP2
#' functions with a pair of FP1 functions, one in \eqn{x} and the other in 
#' \eqn{acd(x)}.
#' 
#' This alternative class of four-parameter functions seems to provide about
#' the same flexibility as the standard FP2 family, but the ACD component offers
#' the additional possibility of sigmoid functions.
#' Royston (2014b) discusses how the extended class of functions known as
#' \eqn{FP1(p1, p2)}, namely
#' \deqn{FP1(p1, p2) = \beta_1 x^{p1} + \beta_2 acd(x)^{p2}}
#' can be fitted optimally by seeking the best combination of all 64 pairs of
#' powers (p1, p2). The optimisation is invoked by use of the `acdx` parameter.
#' Royston (2014b) also described simplification of the chosen function through
#' model reduction by applying significance testing to six sub-families of
#' functions,M1-M6, giving models M1 (most complex) through M6 (null):
#' \itemize{
#' \item{M1.}{FP1(p1, p2) (no simplification)}
#' \item{M2.}{FP1(p1, .) (regular FP1 function of \eqn{x})}
#' \item{M3.}{FP1(., p2) (regular FP1 function of \eqn{acd(x)})}
#' \item{M4.}{FP1(1, .) (linear function of \eqn{x})}
#' \item{M5.}{FP1(., 1) (linear function of \eqn{acd(x)})}
#' \item{M6.}{Null (\eqn{x} omitted entirely)}
#' }
#' Selection among these six sub-functions is performed by a closed test 
#' procedure known as the FSPA. It maintains the familywise type 1 error 
#' probability for selecting \eqn{x} at the value determined by the 
#' `select` parameter. To obtain a 'final' model, a structured sequence of up 
#' to five tests is carried out, the first at the significance level specified 
#' by the `select` parameter, and the remainder at the significance level 
#' provided by the `alpha` option. 
#' The sequence of tests is as follows:
#' \itemize{
#' \item{Test 1.}{Compare the deviances of models 6 and 1 on 4 d.f. 
#'   If not significant then stop and omit \eqn{x}, otherwise continue to step 2.}
#' \item{Test 2.}{Compare the deviances of models 4 and 1 on 3 d.f. 
#'   If not significant then accept model 4 and stop. Otherwise, continue to step 3.}
#' \item{Test 3.}{Compare the deviance of models 2 and 1 on 2 d.f. 
#'   If not significant then accept model 2 and stop. Otherwise continue to step 4.}
#' \item{Test 4.}{Compare the deviance of models 3 and 1 on 2 d.f. 
#'   If significant then model 1 cannot be simplified; accept model 1 and stop.  
#'   Otherwise continue to step 5.}
#' \item{Test 5.}{Compare the deviances of models 5 and 3 on 1 d.f. 
#'   If significant then model 3 cannot be simplified; accept model 3. 
#'   Otherwise, accept model 5. End of procedure.}
#' }
#' The result is the selection of one of the six models. 
#'
#' @param x input matrix, of dimension nobs x nvars; each row is an observation 
#'   vector.
#' @param y vector for response variable. For `family="binomial"` should be  a
#'   variable with two levels (see [stats::glm()]). 
#'   For `family="cox"` it must be a `Surv` object containing  2 columns.
#' @param weights observation weights. Default is `NULL` which assigns a weight 
#'   of 1 to each observation
#' @param offset a vector of length nobs that is included in the linear
#'   predictor. Useful for the poisson family (e.g. log of exposure time).
#'   Default is `NULL` which assigns an offset  of 0 to each observation.
#'   If supplied, then values must also be supplied to the `predict()` function.
#' @param cycles maximum number of iteration cycles. Default is 5.
#' @param scale If the values of the variable are too large or too small,
#' the reported results of fractional polynomials may be difficult to
#' interpret. Scaling can be done automatically or by directly specifying the
#' scaling values so that the magnitude of the `x` values are not too large.
#' By default scaling factors are estimated by the program (see Details section
#' below). Rather than letting `mfpa()` automatically choose the scaling 
#' factors, you may specify scale factors for each variable manually. 
#' If scaling is not required set `scale = 1` to disable it.
#' @param shift Fractional polynomials are only defined for positive term
#' variables. By default, `mfpa()` will assume that all variable values are 
#' positive and attempt to compute fractional powers of the input variables. 
#' By default, estimation of shifting factors will be performed.
#' If the positive value assumption is incorrect, user-supplied shifting factors
#' may also be used to make all variables positive. 
#' If shifting is not required, set `shift = 0` to disable it.
#' @param df a vector of values (or single value) that sets the degrees of
#' freedom (df) for each predictor. The df (not counting the intercept) are
#' twice the degree of a fractional polynomial (FP). For example, an FP2 has 
#' 4 df, while FP3 has 6 df. A single value `default` can be specified and the
#' df for all predictors is taken to be that specific value. 
#' The program overrides default df based on the number of distinct (unique) 
#' values for a variable as follows: 
#' 2-3 distinct values are assigned `df = 1` (linear), 4-5 distinct values are
#' assigned `df = min(2, default)` and >=6 distinct values are assigned  
#' `df = default`.
#' @param center logical; defines the centering of the variables The default is
#' mean centering, except for binary covariates, where the covariate is centered
#' using the lower of the two distinct values of the covariate. 
#' See Details section below.
#' @param family A character string representing a `glm()` family object as well
#' as Cox models. For more information, see Details section below.
#' @param criterion a criterion used to select variables and FP models of
#' different degrees. Default is to use p-values in which the user can specify
#' the significance level (or use default level of 0.05) for variable and
#' functional form selection (see `select` and `alpha` parameters below).
#' If the user select the BIC or AIC criterion then the program would ignore 
#' the nominal significance levels and select variables and functional forms
#' using the chosen information criterion.
#' @param select sets the nominal significance levels for variable selection by
#' backward elimination. A variable is dropped if its removal causes a 
#' non-significant increase in deviance. The rules for `select` are the same as
#' those for `df` (see above). The default nominal significance level is 0.05 
#' for all variables. Setting the nominal significance level to be 1 for a 
#' given variable forces it into the model, leaving others to be selected. 
#' Using the default selection level of 1 for all variables forces them all 
#' into the model.
#' @param alpha sets the significance levels for testing between FP models of 
#' different degrees. The rules for `alpha` are the same as those for `df` 
#' (see above). The default nominal significance level is 0.05 for all 
#' variables. Example: `alpha = 0.05` specifies that all variables have an FP 
#' selection level of 0.05. 
#' @param keep keep one or more variables in the model. In case that 
#' `criterion = "pvalue"`, this is equivalent to setting the selection level for
#' the variables to keep to 1. However, this option also keeps the specified
#' variables in the model when when using the BIC or AIC criteria. 
#' @param xorder determines the order of entry of the covariates into
#' the model-selection algorithm. The default is `ascending`, which uses them
#' in decreasing order of significance in a multiple linear regression
#' (most significant first). `descending` places them in reverse significance 
#' order, whereas `original` respects the original order of the covariates. 
#' @param powers sets the permitted FP powers for all covariates. Default is
#' `powers = c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)`, where 0 means natural logarithm.
#' @param ties a character string specifying the method for tie handling in 
#' Cox regression. If there are no tied death times all the methods are 
#' equivalent. Default is the Breslow method. This argument is used for Cox 
#' models only and has no effect for other model families. 
#' See [survival::coxph()] for details.
#' @param strata variables to be used for stratification in a Cox model. 
#' A new factor, whose levels are all possible combinations of the variables 
#' supplied as arguments will be created. Default is `NULL` and a Cox model 
#' without stratification would be fitted. See [survival::coxph()] for details.
#' @param nocenter an optional list of values for fitting Cox models.
#' See [survival::coxph()] for details.
#' @param acdx a vector of names of continuous variables that undergo the
#' approximate cumulative distribution (ACD) transformation.
#' It also invokes `FSPA` to determine the best-fitting FP1(p1, p2) model 
#' (see details section). The variable representing the ACD transformation of 
#' `x` is named `A(x)`.
#' @param ftest logical; for normal error models with small samples, critical 
#' points from the F-distribution can be used instead of Chi-Square 
#' distribution. Default uses the later. This argument is used for Gaussian 
#' models only and has no effect for other model families.
#' @param verbose logical; run in verbose mode.
#' @references
#' Royston, P. 2014. \emph{A smooth covariate rank transformation for use in
#' regression models with asigmoid dose-response function. 
#' Stata Journal 14(2): 329-341.}\cr
#' Royston, P. and Sauerbrei, W., 2016. \emph{mfpa: Extension of mfp using the
#' ACD covariate transformation for enhanced parametric multivariable modeling. 
#' The Stata Journal, 16(1), pp.72-87.}
#' @export
mfpa <- function(x, y, weights = NULL, offset = NULL, cycles = 5,
                 scale = NULL, shift = NULL, df = 4, center = T,
                 family = c("gaussian", "poisson", "binomial", "cox"),
                 criterion = c("pvalue", "BIC", "AIC"),
                 select = 0.05, alpha = 0.05,
                 keep = NULL,
                 xorder = c("ascending", "descending", "original"),
                 powers = NULL,
                 ties = c("breslow", "efron"),
                 strata = NULL,
                 nocenter = NULL,
                 acdx = NULL,
                 ftest = F,
                 verbose = T) {
  cl <- match.call()
  
  # match arguments 
  criterion <- match.arg(criterion)
  xorder <- match.arg(xorder)
  family <- match.arg(family)
  ties <- match.arg(ties)
  
  # assertions 
  # assert dimension of x
  np <- dim(x)
  nobs <- as.integer(np[1])
  nvars <- as.integer(np[2])
  # assert that x is a matrix
  if (is.null(np)) {
      stop("! The dimensions of x must not be Null.\n",
           "i Please make sure that x is a matrix or data.frame with at least one row and column.")
  }
  # assert that x has column names
  vnames <- colnames(x)
  if (is.null(vnames)) stop("! The column names of x must not be Null.\n",
                            "i Please set column names for x.")
  # Transform factor variables to dummies if any using model.matrix
  x <- model.matrix(~ . - 1, data.frame(x)) 
  
  # assert that x has no missing data
  if (anyNA(x)) stop("! x must not contain any NA (missing data).\n", 
                     "i Please remove any missing data before passing x to this function.")
  
  # assert that weights are positive and of appropriate dimensions
  if (!is.null(weights)) {
    if (any(weights < 0)) {
        stop("! Weights must not be negative.")
    }
    if (length(weights) != nobs) {
      stop("! The number of observations (rows in x) and weights must match.\n", 
           sprintf("i The number of rows in x is %d, but the number of elements in weights is %d.", 
                   nobs, length(weights)))
    }
  }
  
  # assert dimensions of offset
  if (!is.null(offset)) {
      if (length(offset) != nobs) {
          stop("! The number of observations (rows in x) and offset must match.\n", 
               sprintf("i The number of rows in x is %d, but the number of elements in offset is %d.", 
                       nobs, length(offset)))
      }
  }
  
  # assert alpha between 0 and 1
  if (any(alpha > 1) || any(alpha < 0)) {
      stop("! alpha must not be < 0 or > 1.")
  }
  
  # assert length of alpha 
  if (length(alpha) != 1 && length(alpha) != nvars) {
      stop("! alpha must be a single number, or the number of variables (columns in x) and alpha must match.\n", 
           sprintf("i The number of variables in x is %d, but the number of elements in alpha is %d.", 
                   nvars, length(alpha)))
  }
  
  # assert select between 0 and 1
  if (any(select > 1) || any(select < 0)) {
      stop("! select must not be < 0 or > 1.")
  }
  # assert length of select 
  if (length(select) != 1 && length(select) != nvars) {
      stop("! select must be a single number, or the number of variables (columns in x) and select must match.\n", 
           sprintf("i The number of variables in x is %d, but the number of elements in select is %d.", 
                   nvars, length(select)))
  }
  
  # assert keep is a subset of x
  if (!is.null(keep)) {
      if (!all(keep %in% colnames(x))) {
          warning("i The set of variables named in keep is not a subset of the variables in x.\n", 
                  "i mfpa() continues with the intersection of keep and colnames(x).")
      }
  }
  
  # assert shift vector is of correct dimension
  if (!is.null(shift)) {
      if (length(shift) != 1 && length(shift) != nvars) {
          stop("! shift must either be NULL, a single number, or the number of variables (columns in x) and shift must match.\n", 
               sprintf("i The number of variables in x is %d, but the number of elements in shift is %d.", 
                       nvars, length(shift)))
      }
  }
  
  # assert scale vector is of correct dimension
  if (!is.null(scale)) {
      if (length(scale) != 1 && length(scale) != nvars) {
          stop("! scale must either be NULL, a single number, or the number of variables (columns in x) and scale must match.\n", 
               sprintf("i The number of variables in x is %d, but the number of elements in shift is %d.", 
                       nvars, length(scale)))
      }
  }
  
  # assert center vector is of correct dimension
  if (length(center) != 1) {
      if (length(center) != nvars) {
          stop("! center must either be of length 1, or the number of variables (columns in x) and center must match.\n", 
               sprintf("i The number of variables in x is %d, but the number of elements in center is %d.", 
                       nvars, length(center)))
      }
  }
  
  # assert acdx is a subset of x
  if (!is.null(acdx)) {
      if (!all(acdx %in% colnames(x))) {
          warning("i The set of variables named in acdx is not a subset of the variables in x.\n", 
                  "i mfpa() continues with the intersection of acdx and colnames(x).")
      }
  }
  
  # assert df is positive
  if (any(df) <= 0) {
      stop("! df must not be 0 or negative.\n", 
           sprintf("i All df must be either 1 (linear) or 2m, where m is the degree of FP."))
  }
  
  if (length(df) == 1) {
      # assert df is 1 or even 
      if (df != 1 && df %% 2 != 0) {
          stop("! Any df > 1 must not be odd.\n", 
               sprintf("i df = %d was passed, but df must be either 1 (linear) or 2m, where m is the degree of FP.", 
                       df))
      } 
  } else {
      # assert length of df
      if (length(df) != nvars) {
          stop("! df must be a single number, or the number of variables (columns in x) and df must match.\n", 
               sprintf("i The number of variables in x is %d, but the number of elements in df is %d.", 
                       nvars, length(df)))
      }
      # assert all df are 1 or even 
      if (any(df != 1 & df %% 2 != 0)) {
          stop("! Any df > 1 must not be odd.\n", 
               sprintf("i All df must be either 1 (linear) or 2m, where m is the degree of FP.", 
                       df))
      }
  }
  
  # assert ftest and family are compatible
  if (ftest && family != "gaussian") {
      warning(sprintf("i F-test not suitable for family = %s.\n", family),
              "i mfpa() reverts to use Chi-square instead.")
  }
  
  if (family == "cox") {
      # assert y is a Surv object
      if (!survival::is.Surv(y)) {
          stop("! Response y must be a survival::Surv object.")
      }
      
      # assert dimensions of y 
      if (nrow(y) != nobs) {
          stop("! Number of observations in y and x must match.", 
               sprintf("i The number of observations in y is %d, but the number of observations in x is %d.", 
                       nrow(y), nobs))
      }
      
      # assert right censoring (other censoring types are not implemented yet)
      type <- attr(y, "type")
      if (type != "right") {
          stop(sprintf("! Type of censoring must not be %s.", type), 
               "i Currently only right censoring is supported by mfpa().")
      }
      
      if (!is.null(strata)) {
          # assert stratification factors are in x
          if (!all(strata %in% colnames(x))) {
              warning("i The set of variables named in strata is not a subset of the variables in x.\n", 
                      "i mfpa() continues with the intersection of strata and colnames(x).")
          }
      }
  } else {
      # assert type of y
      if (!is.vector(y)) {
          stop(sprintf("! Outcome y must not be of class %s.", class(y)), 
               "i Please convert y to a vector.")
      }
      if (length(y) != nobs) {
          stop("! Number of observations in y and x must match.", 
               sprintf("i The number of observations in y is %d, but the number of observations in x is %d.", 
                       length(y), nobs))
      }
  }

  # set defaults
  if (is.null(powers)) {
      # default FP powers proposed by Royston and Sauerbrei (2008)
      powers <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
  }
  if (is.null(weights)) {
      weights <- rep.int(1, nobs)
  }
  if (is.null(offset)) {
      offset <- rep.int(0, nobs)
  }
  if (length(select) == 1) {
      select <- rep(select, nvars)
  } 
  if (length(alpha) == 1) {
      alpha <- rep(alpha, nvars)
  } 
  if (is.null(shift)) {
      shift <- apply(x, 2, shift.factor)
  } else if (length(shift == 1)) {
      shift <- rep(shift, nvars)
  }
  if (is.null(scale)) {
      scale <- apply(x, 2, scalefn)
  } else if (length(scale == 1)) {
      scale <- rep(scale, nvars)
  }
  if (length(center) == 1) {
      center <- rep(center, nvars)    
  }
  if (ftest && family != "gaussian") {
      ftest = FALSE
  }
  
  keep <- intersect(keep, colnames(x))
  strata <- intersect(strata, colnames(x))
  # convert acdx to logical vector
  if (is.null(acdx)) {
      acdx <- rep(F, nvars)
  } else {
      # Get rid of duplicates names
      acdx <- unique(acdx)
      acdx <- intersect(acdx, colnames(x))
      # the variables that undergo acd transformation have acdx = TRUE
      acdx <- replace(rep(FALSE, nvars), 
                      which(vnames %in% acdx), rep(TRUE, length(acdx)))
  }
  
  # further variables
  istrata <- strata
  # control is specific to coxph models, plays no role for other models
  control <- survival::coxph.control() 

  # set df 
  if (length(df) == 1) {
    if (df != 1) {
      # we assign variables different df based on number of unique values
      df.list <- df.assign(x = x, df.default = df)
    } else {
      df.list <- rep(df, nvars)
    }
  } else {
    # ensure that variables with <= 3 unique values have df = 1
    nux <- apply(x, 2, function(v) length(unique(v)))
    index <- nux <= 3
    if (any(df[index] != 1)) {
        warning("i For any variable with fewer than 4 unique valies the df are set to 1 (linear) by mfpa().\n", 
                sprintf("i This applies to variables %s.", 
                        colnames(x)[index & df != 1]))
        df[index] = 1
    }
    df.list <- df
  }

  # data preparation: shift, scale
  x <- sweep(x, 2, shift, "+")
  x <- sweep(x, 2, scale, "/")
  
  # data preparation: stratification
  if (family == "cox" && !is.null(strata)) {
      istrata <- survival::strata(x[, strata], shortlabel = TRUE)
      # drop the variable(s) in x used for stratification
      x <- x[, -c(which(colnames(x) %in% strata)), drop = FALSE]
  }
  
  # fit model and make model specific adaptions
  fit <- mfp.fit(
      x = x, y = y, weights = weights, offset = offset, cycles = cycles,
      scale = scale, shift = shift, df = df.list, keep,
      center = center, criterion = criterion, xorder = xorder,
      powers = powers, family = family, method = ties,
      select = select, alpha = alpha, strata = istrata,
      ftest = ftest, verbose = verbose, control = control,
      nocenter = nocenter, rownames = row.names(x), acdx = acdx
  )
  fit$call <- cl
  
  if (family == "cox") {
      class(fit) <- c("mfpa", "coxph")
      # add wald test in order to use summary.coxph()
      # this calculation follows coxph() source code
      nabeta <- !is.na(fit$coefficients)
      fit$wald.test <- survival::coxph.wtest(
          fit$var[nabeta, nabeta], fit$coefficients[nabeta],
          control$toler.chol
      )$test
  } else {
      class(fit) <- c("mfpa", "glm", "lm")
  }
  
  fit
}
