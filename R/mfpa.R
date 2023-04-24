#' Multivariable fractional polynomial models with extensions of sigmoid 
#' functions
#'
#' Selects the multivariable fractional polynomial (FP) model that best predicts
#' the outcome variable. It also has the ability to model a sigmoid relationship
#' between x and an outcome variable using the approximate cumulative 
#' distribution (ACD) transformation proposed by Royston (2016).
#' 
#' @section Brief summary of FPs:
#' 
#' In the following we denote fractional polynomials for a variable \eqn{x} by 
#' increasing complexity as either FP1 or FP2. \eqn{FP2(p1, p2)} is the most 
#' flexible FP tranformation where 
#' \deqn{FP2(p1, p2) = \beta_1 x^{p1} + \beta_2 x^{p2}.}
#' The (fractional) powers \eqn{p1} and \eqn{p2} are taken from a set
#' of allowed powers, usually {-2, -1, -0.5, 0, 0.5, 1, 2, 3} where the power
#' 0 indicates the logarithm. The optimal FP is then found by a closed testing
#' procedure that seeks the best combination from all 36 pairs of
#' powers \eqn{(p1, p2)}. Functions that only involve a single power of
#' the variable are denoted as FP1, i.e. 
#' \deqn{FP1(p1) = \beta_1 x^{p1}.}
#' For details see e.g. Sauerbrei et al (2006). 
#'
#' @section Details on `family` option:
#'
#' `mfpa()` supports the family object as used by [stats::glm()]. The built in
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
#' @section Details on shifting, scaling and centering:
#' 
#' Fractional polynomials are only defined for positive variables due to the 
#' use of logarithms. Thus, `mfpa()` estimates shifts for each variables to 
#' achieve this, or assumes that this is the case when computing fractional 
#' powers of the input variables in case that shifting is disabled manually. 
#' 
#' If the values of the variables are too large or too small, the reported 
#' results of fractional polynomials may be difficult to interpret. 
#' Scaling can be done automatically or by directly specifying the
#' scaling values so that the magnitude of the `x` values are not too large.
#' By default scaling factors are estimated by the program as follows.
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
#' So in brief: shifting is required to make input values positive, scaling
#' helps to bring the values to a reasonable range. Both operations are 
#' applied before applying the FP transformation to an input variable. 
#' Centering, however, is done after applying the FP transformation.
#' Also see [transform_vector_fp()] for some more details.
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
#' This alternative class of four-parameter functions provides about
#' the same flexibility as the standard FP2 family, but the ACD component offers
#' the additional possibility of sigmoid functions.
#' Royston (2014b) discusses how the extended the class of functions known as
#' \eqn{FP1(p1, p2)}, namely
#' \deqn{FP1(p1, p2) = \beta_1 x^{p1} + \beta_2 acd(x)^{p2}}
#' can be fitted optimally by seeking the best combination of all 64 pairs of
#' powers (p1, p2). The optimisation is invoked by use of the `acdx` parameter.
#' Royston (2014b) also described simplification of the chosen function through
#' model reduction by applying significance testing to six sub-families of
#' functions,M1-M6, giving models M1 (most complex) through M6 (null):
#' \itemize{
#' \item{M1: }{FP1(p1, p2) (no simplification)}
#' \item{M2: }{FP1(p1, .) (regular FP1 function of \eqn{x})}
#' \item{M3: }{FP1(., p2) (regular FP1 function of \eqn{acd(x)})}
#' \item{M4: }{FP1(1, .) (linear function of \eqn{x})}
#' \item{M5: }{FP1(., 1) (linear function of \eqn{acd(x)})}
#' \item{M6: }{Null (\eqn{x} omitted entirely)}
#' }
#' Selection among these six sub-functions is performed by a closed test 
#' procedure known as the function-selection pocedure FSPA. 
#' It maintains the familywise type 1 error 
#' probability for selecting \eqn{x} at the value determined by the 
#' `select` parameter. To obtain a 'final' model, a structured sequence of up 
#' to five tests is carried out, the first at the significance level specified 
#' by the `select` parameter, and the remainder at the significance level 
#' provided by the `alpha` option. 
#' The sequence of tests is as follows:
#' \itemize{
#' \item{Test 1: }{Compare the deviances of models 6 and 1 on 4 d.f. 
#' If not significant then stop and omit \eqn{x}, otherwise continue to step 2.}
#' \item{Test 2: }{Compare the deviances of models 4 and 1 on 3 d.f. 
#' If not significant then accept model 4 and stop. Otherwise, continue to step 3.}
#' \item{Test 3: }{Compare the deviance of models 2 and 1 on 2 d.f. 
#' If not significant then accept model 2 and stop. Otherwise continue to step 4.}
#' \item{Test 4: }{Compare the deviance of models 3 and 1 on 2 d.f. 
#' If significant then model 1 cannot be simplified; accept model 1 and stop.  
#' Otherwise continue to step 5.}
#' \item{Test 5: }{Compare the deviances of models 5 and 3 on 1 d.f. 
#' If significant then model 3 cannot be simplified; accept model 3. 
#' Otherwise, accept model 5. End of procedure.}
#' }
#' The result is the selection of one of the six models. 
#'
#' @param x an input matrix of dimension nobs x nvars. Each row is an 
#' observation vector.
#' @param y a vector for the response variable. For `family="binomial"` it 
#' should be  a vector with two levels (see [stats::glm()]). 
#' For `family="cox"` it must be a `Surv` object containing  2 columns.
#' @param weights a vector of observation weights of length nobs. 
#' Default is `NULL` which assigns a weight of 1 to each observation.
#' @param offset a vector of length nobs that is included in the linear
#' predictor. Useful for the poisson family (e.g. log of exposure time).
#' Default is `NULL` which assigns an offset  of 0 to each observation.
#' If supplied, then values must also be supplied to the `predict()` function.
#' @param cycles an integer, maximum number of iteration cycles. Default is 5.
#' @param scale a numeric vector of length nvars or single numeric specifying 
#' scaling factors. If a single numeric, then the value will be replicated as
#' necessary. Default is `NULL` which lets the program estimate the scaling 
#' factors (see Details section).
#' If scaling is not required set `scale = 1` to disable it.
#' @param shift a numeric vector of length nvars or a single numeric specifying
#' shift terms. If a single numeric, then the value will be replicated as
#' necessary. Default is `NULL` which lets the program estimate the shifts
#' (see Details section).
#' If shifting is not required, set `shift = 0` to disable it.
#' @param df a numeric vector of length nvars or a single numeric that sets the 
#' (default) degrees of freedom (df) for each predictor. If a single numeric, 
#' then the value will be replicated as necessary. The df (not counting 
#' the intercept) are twice the degree of a fractional polynomial (FP). 
#' For example, an FP2 has 4 df, while FP3 has 6 df. 
#' The program overrides default df based on the number of distinct (unique) 
#' values for a variable as follows: 
#' 2-3 distinct values are assigned `df = 1` (linear), 4-5 distinct values are
#' assigned `df = min(2, default)` and >= 6 distinct values are assigned  
#' `df = default`.
#' @param center a logical determining whether variables are centered before 
#' model fit. The default `TRUE` implies mean centering, except for binary 
#' covariates, where the covariate is centered using the lower of the two 
#' distinct values of the covariate. See Details section below.
#' @param family a character string representing a `glm()` family object as well
#' as Cox models. For more information, see Details section below.
#' @param criterion a character string defining the criterion used to select 
#' variables and FP models of different degrees. 
#' Default is to use p-values in which case the user can specify
#' the significance level (or use default level of 0.05) for variable and
#' functional form selection (see `select` and `alpha` parameters below).
#' If the user specifies the BIC (`bic`) or AIC (`aic`) criteria the program  
#' ignores the nominal significance levels and selects variables and functional 
#' forms using the chosen information criterion.
#' @param select a numeric vector of length nvars or a single numeric that 
#' sets the nominal significance levels for variable selection by
#' backward elimination. If a single numeric, then the value will be replicated
#' as necessary. The default nominal significance level is 0.05 
#' for all variables. Setting the nominal significance level to be 1 for  
#' certain variables forces them into the model, leaving all other variables
#' to be selected. 
#' @param alpha a numeric vector of length nvars or a single numeric that 
#' sets the significance levels for testing between FP models of 
#' different degrees. If a single numeric, then the value will be replicated
#' as necessary. The default nominal significance level is 0.05 for all 
#' variables. 
#' @param keep a character vector that with names of variables to be kept 
#' in the model. In case that `criterion = "pvalue"`, this is equivalent to
#' setting the selection level for the variables in `keep` to 1. 
#' However, this option also keeps the specified variables in the model when 
#' using the BIC or AIC criteria. 
#' @param xorder a string determining the order of entry of the covariates
#' into the model-selection algorithm. The default is `ascending`, which enters
#' them by ascending p-values, or decreasing order of significance in a
#' multiple regression (i.e. most significant first).
#' `descending` places them in reverse significance order, whereas 
#' `original` respects the original order in `x`.
#' @param powers a numeric vector that sets the permitted FP powers for all 
#' covariates. Default is `powers = c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)`,
#' where 0 means natural logarithm. Duplicates are removed, and powers are 
#' sorted before further processing in the program. 
#' @param ties a character string specifying the method for tie handling in 
#' Cox regression. If there are no tied death times all the methods are 
#' equivalent. Default is the Breslow method. This argument is used for Cox 
#' models only and has no effect for other model families. 
#' See [survival::coxph()] for details.
#' @param strata a character vector specifying the variables to be used for 
#' stratification in a Cox model. A new factor, whose levels are all possible 
#' combinations of the variables supplied as arguments will be created. 
#' Default is `NULL` and a Cox model without stratification would be fitted. 
#' See [survival::coxph()] for details. Currently only a single stratification
#' factor is supported by `mfpa()`.
#' @param nocenter a numeric vector with a list of values for fitting Cox 
#' models. See [survival::coxph()] for details.
#' @param acdx a character vector of names of continuous variables to undergo 
#' the approximate cumulative distribution (ACD) transformation.
#' It also invokes the function-selection procedure to determine the 
#' best-fitting FP1(p1, p2) model (see Details section). 
#' The variable representing the ACD transformation of `x` is named `A(x)`.
#' @param ftest a logical; for normal error models with small samples, critical 
#' points from the F-distribution can be used instead of Chi-Square 
#' distribution. Default `FALSE` uses the latter. This argument is used for 
#' Gaussian models only and has no effect for other model families.
#' @param control a list object with parameters controlling model fit details. 
#' Returned by either [stats::glm.control()] or [survival::coxph.control()]. 
#' Default is `NULL` to use default parameters for the given model class. 
#' @param verbose a logical; run in verbose mode.
#' 
#' @return 
#' `mfpa()` returns an object of class inheriting from `glm` or `copxh`, 
#' depending on the `family` parameter. 
#' 
#' The function `summary()` (i.e. [summary.mfpa()]) can be used to obtain or
#' print a summary of the results. 
#' 
#' The generic accessor function `coef()` can be used to extract the vector of 
#' coefficients from the fitted model object.
#' 
#' An object of class `mfpa` is a list containing at least the following 
#' components:  
#' \itemize{
#' \item{convergence_mfp: }{logical value indicating convergence of mfp algorithm.}
#' \item{coefficients: }{a named vector of coefficients.}
#' \item{residuals: }{working residuals as returned by [stats::glm()] or 
#' martingale residuals as returned by [survival::coxph()], depending on 
#' `family`.}
#' \item{family: }{either the [stats::family()] object used for families supported
#' by [stats::glm()] or "cox" for Cox proportional hazards models.}
#' \item{linear.predictors: }{the vector of linear predictors, i.e. the linear fit
#' on link scale. For `family = "cox"` this vector has been centered, see
#' [survival::predict.coxph()] for more details.}
#' \item{fp_terms: }{a data.frame with information on fractional polynomial 
#' terms.}
#' \item{transformations: }{a data.frame with information on shifting, scaling
#' and centering for all variables.}   
#' \item{fp_powers: }{a list with all powers of fractional polynomial terms.}
#' \item{acd: }{a vector with information for which variables the acd 
#' transformation was applied.}
#' \item{x: }{the scaled and shifted input matrix but without transformations.}
#' \item{y: }{the original outcome variable.}
#' \item{X: }{the final transformed input matrix used to fit the final model.}
#' }
#' The `mfpa` object may contain further information depending on family.
#' 
#' @references
#' Royston, P. and Sauerbrei, W., 2008. \emph{Multivariable Model - Building: 
#' A Pragmatic Approach to Regression Anaylsis based on Fractional Polynomials 
#' for Modelling Continuous Variables. John Wiley & Sons.}\cr
#' Sauerbrei, W., Meier-Hirmer, C., Benner, A. and Royston, P., 2006. 
#' \emph{Multivariable regression model building by using fractional 
#' polynomials: Description of SAS, STATA and R programs. 
#' Comput Stat Data Anal, 50(12): 3464-85.}\cr
#' Royston, P. 2014. \emph{A smooth covariate rank transformation for use in
#' regression models with a sigmoid dose-response function. 
#' Stata Journal 14(2): 329-341.}\cr
#' Royston, P. and Sauerbrei, W., 2016. \emph{mfpa: Extension of mfp using the
#' ACD covariate transformation for enhanced parametric multivariable modeling. 
#' The Stata Journal, 16(1), pp.72-87.}\cr
#' Sauerbrei, W. and Royston, P., 1999. \emph{Building multivariable prognostic 
#' and diagnostic models: transformation of the predictors by using fractional 
#' polynomials. J Roy Stat Soc a Sta, 162:71-94.}
#' 
#' @seealso 
#' [summary.mfpa()], [coef.mfpa()]
#' 
#' @export
mfpa <- function(x, 
                 y, 
                 weights = NULL, 
                 offset = NULL, 
                 cycles = 5,
                 scale = NULL, 
                 shift = NULL, 
                 df = 4, 
                 center = TRUE,
                 family = c("gaussian", "poisson", "binomial", "cox"),
                 criterion = c("pvalue", "aic", "bic"),
                 select = 0.05, 
                 alpha = 0.05,
                 keep = NULL,
                 xorder = c("ascending", "descending", "original"),
                 powers = NULL,
                 ties = c("breslow", "efron"),
                 strata = NULL,
                 nocenter = NULL,
                 acdx = NULL,
                 ftest = FALSE,
                 control = NULL, 
                 verbose = TRUE) {
  
  # this function prepares everything for fitting the actual mfpa model
  
  cl <- match.call()
  
  # match arguments ------------------------------------------------------------
  criterion <- match.arg(criterion)
  xorder <- match.arg(xorder)
  family <- match.arg(family)
  ties <- match.arg(ties)
  
  # assertions -----------------------------------------------------------------
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
  if (any(df <= 0)) {
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
          stop(sprintf("! Outcome y must not be of class %s.", 
                       paste0(class(y), collapse = ", ")), 
               "i Please convert y to a vector.")
      }
      if (length(y) != nobs) {
          stop("! Number of observations in y and x must match.", 
               sprintf("i The number of observations in y is %d, but the number of observations in x is %d.", 
                       length(y), nobs))
      }
  }

  # set defaults ---------------------------------------------------------------
  if (is.null(powers)) {
      # default FP powers proposed by Royston and Sauerbrei (2008)
      powers <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
  }
  powers <- sort(unique(powers))
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
      shift <- apply(x, 2, find_shift_factor)
  } else if (length(shift) == 1) {
      shift <- rep(shift, nvars)
  }
  if (is.null(scale)) {
      scale <- apply(x, 2, find_scale_factor)
  } else if (length(scale) == 1) {
      scale <- rep(scale, nvars)
  }
  if (length(center) == 1) {
      center <- rep(center, nvars)    
  }
  if (ftest && family != "gaussian") {
      ftest <- FALSE
  }
  if (is.null(control)) {
    if (family == "cox") {
      control = survival::coxph.control()
    } else {
      control = stats::glm.control()
    }
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
  
  # set df ---------------------------------------------------------------------
  if (length(df) == 1) {
    if (df != 1) {
      # we assign variables different df based on number of unique values
      df.list <- assign_df(x = x, df_default = df)
    } else {
      df.list <- rep(df, nvars)
    }
  } else {
    # ensure that variables with <= 3 unique values have df = 1
    nux <- apply(x, 2, function(v) length(unique(v)))
    index <- nux <= 3
    if (any(df[index] != 1)) {
        warning("i For any variable with fewer than 4 unique values the df are set to 1 (linear) by mfpa().\n", 
                sprintf("i This applies to the following variables: %s.", 
                        paste0(colnames(x)[index & df != 1], collapse = ", ")))
        df[index] <- 1
    }
    df.list <- df
  }

  # data preparation -----------------------------------------------------------
  # shift, scale 
  x <- sweep(x, 2, shift, "+")
  x <- sweep(x, 2, scale, "/")
  
  # stratification
  istrata <- strata
  if (family == "cox" && !is.null(strata)) {
      istrata <- survival::strata(x[, strata], shortlabel = TRUE)
      # drop the variable(s) in x used for stratification
      x <- x[, -c(which(colnames(x) %in% strata)), drop = FALSE]
  }
  
  # fit model ------------------------------------------------------------------
  fit <- fit_mfp(
      x = x, y = y, 
      weights = weights, offset = offset, cycles = cycles,
      scale = scale, shift = shift, df = df.list, center = center, 
      family = family, criterion = criterion, select = select, alpha = alpha, 
      keep = keep, xorder = xorder, powers = powers, 
      method = ties, strata = istrata, nocenter = nocenter,
      acdx = acdx, ftest = ftest, 
      control = control, 
      verbose = verbose
  )
  fit$call <- cl
  
  fit
}

#' Extract coefficients from object of class `mfpa`
#' 
#' This function is a method for the generic [stats::coef()] function for 
#' objects of class `mfpa`. 
#' 
#' @param object an object of class `mfpa`, usually, a result of a call to
#' [mfpa()].
#' 
#' @return 
#' Named numeric vector of coefficients extracted from the model `object`.
#' 
#' @export
coef.mfpa <- function(object) {
  object$coefficients
}

#' Summarizing `mfpa` model fits
#' 
#' This function is a method for the generic [base::summary()] function for
#' objects of class `mfpa`.
#' 
#' @param object an object of class `mfpa`, usually, a result of a call to
#' [mfpa()].
#' @param ... further arguments passed to the summary functions for `glm()` 
#' ([stats::summary.glm()], i.e. families supported by `glm()`) or `coxph()` 
#' ([survival::summary.coxph()], if `object$family = "cox"`).
#' 
#' @return 
#' An object returned from [stats::summary.glm()] or
#' [survival::summary.coxph()], depending on the family parameter of `object`.
#' 
#' @seealso 
#' [mfpa()], [stats::glm()], [stats::summary.glm()], [survival::coxph()],
#' [survival::summary.coxph()]
#' 
#' @export
summary.mfpa <- function(object, ...) {
  NextMethod("summary", object)
}

#' Print method for objects of class `mfpa`
#' 
#' Enhances printing by information on data processing and fractional 
#' polynomials.
#' 
#' @export
print.mfpa <- function(x,
                       digits = max(3L, getOption("digits") - 3L), 
                       signif.stars = FALSE, 
                       ...) {
  # shift and scaling factors with centering values
  cat("Shifting, Scaling and Centering of covariates", "\n")
  print.data.frame(x$transformations)
  cat("\n")
  
  # Final MFP Powers
  cat("Final Multivariable Fractional Polynomial for y", "\n")
  print.data.frame(x$fp_terms)
  cat("\n")
  
  cat(sprintf("MFP algorithm convergence: %s\n", x$convergence_mfp))
  
  # print model object using underlying print function
  NextMethod("print", x)
}

#' Helper to assign degrees of freedom
#' 
#' Determine the number of unique values in a variable. To be used in [mfpa()].
#' 
#' @details 
#' Variables with fewer than or equal to three unique values, for example,
#' will be assigned df = 1. df = 2 will be assigned to variables with 4-5 
#' unique values, and df = 4 will be assigned to variables with unique values 
#' greater than or equal to 6.
#' 
#' @param x input matrix.
#' @param df_default default df to be used. Default is 4.
#' 
#' @return 
#' Vector of length `ncol(x)` with degrees of freedom for each variable in `x`.
assign_df <- function(x, 
                      df_default = 4) {
  
  # default degrees of freedom for each variable
  df <- rep(df_default, ncol(x))
  # unique values of each column of x
  nu <- apply(x, 2, function(x) length(unique(x)))
  
  index1 <- which(nu <= 3)
  index2 <- which(nu >= 4 & nu <= 5)
  if (length(index1) != 0) {
    df <- replace(df, index1, rep(1, length(index1)))
  }
  if (length(index2) != 0) {
    df <- replace(df, index1, rep(1, length(index1)))
  }
  
  df
}
