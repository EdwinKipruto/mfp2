#' Multivariable fractional polynomial models with extensions of sigmoid functions
#'
#' selects the multivariable fractional polynomial model that best predicts
#' the outcome variable. It also has the ability to model a sigmoid relationship
#' between x and an outcome variable using the approximate cumulative distribution
#' (ACD) transformation proposed by Royston (2016).
#'
#' ## Details on `family` option
#'
#' mfpa supports the family object as used by stats:glm(). The built in
#' families are specifed via a character string. \code{mfpa(...,
#' family="binomial")} fits a logistic regression model, while \code{mfpa(...,
#' family="gaussian")} fits OLS regression model.
#'
#' For Cox models, the response should preferably be a \code{Surv} object,
#' created by the \code{Surv()} function in \pkg{survival} package and the \code{family = cox}. Only
#' right-censored data are currently supported. To fit stratified Cox
#' models, strata option can be used.
#'
#' ## Details on scaling and centering
#'
#' After adjusting the location of \eqn{x} so that its minimum value is positive,
#' creating \eqn{x^*}, automatic scaling will divide each value of \eqn{x^*} by \eqn{10^p} where the exponent \eqn{p} is given by
#' \deqn{ p = sign(k)\times floor(|k|) \quad \text{where} \quad k = log_{10} (max(x^*)- min(x^*))}
#'
#' The FP transformation of \eqn{x^*} is centered on the mean of the observed values of \eqn{x^*}. For example,
#' for the FP1 model \eqn{\beta_0 + \beta_1x^p}, the actual model fitted by the software would
#' be \eqn{\beta^*_0 + \beta^*_1(x^{*p}-mean(x^{*p}))}. This approach ensures that
#' the revised constant \eqn{\beta^*_0} equals the fitted value of the FP function at the mean of \eqn{x^*}.
#'
#' ## Details on  approximate cumulative distribution (ACD) transformation `acdx`
#'
#' The ACD transformation (Royston 2014a) converts each predictor, x, smoothly
#' to an approximation, acd(x), of its empirical cumulative distribution
#' function. This is done by smoothing a probit transformation of the scaled
#' ranks of x on x. acd(x) could be used instead of x as a covariate. This
#' has the advantage of providing sigmoid curves in x, something that regular FP
#' functions cannot achieve. Details of the precise definition and some possible
#' uses of the ACD transformation in a univariate context are given by Royston
#'  (2014a). Royston (2014b) describes how one could go further and replace FP2
#'  functions with a pair of FP1 functions, one in x and the other in ACD(x).
#'  This alternative class of four-parameter functions seems to provide about
#'  the same flexibility as the standard FP2 family. The ACD component offers
#'  the additional possibility of sigmoid functions, as just described.
#'  Royston (2014b) discusses how the extended class of functions known as
#'  FP1(p1, p2), namely
#'  \deqn{FP1(p1, p2) = \beta_1 * x^{p1} + \beta_2 * ACD(x)^{p2}}
#' can be fitted optimally by seeking the best combination of all 64 pairs of
#' powers (p1, p2). The optimisation is invoked by use of the \code{acd()} option.
#' Royston (2014b) also described simplification of the chosen function through
#' model reduction by applying significance testing to six sub-families of
#'  functions,M1-M6, giving models M1 (most complex) through M6 (null, x omitted):
#' \itemize{
#' \item{M1.}{FP1(p1, p2) (no simplification)}
#' \item{M2.} {FP1(p1, .) (regular FP1 function of x)}
#' \item{M3.} {FP1(., p2) (regular FP1 function of ACD(x))}
#' \item{M4.} {FP1(1, .) (linear function of x)}
#' \item{M5.} {FP1(., 1) (linear function of ACD(x))}
#' \item{M6.} {Null (x omitted entirely)}
#' }
#' Selection among these six sub-functions is performed by a closed test procedure
#' known as the FSPA. It maintains the familywise type 1 error probability for
#' selecting x at the value determined by the \code{select} option. To obtain a
#' 'final' model, a structured sequence of up to five tests is carried out, the
#' first at the \code{select} significance level and the remainder at the
#' significance level provided by the \code{alpha} option.
#'  The sequence of tests is as follows:
#' \itemize{
#'  \item{Test 1.}{Compare the deviances of models 6 and 1 on 4 d.f. If
#'   non-significant then stop and omit x, otherwise continue to step 2.}
#'  \item{Test 2.}{Compare the deviances of models 4 and 1 on 3 d.f. If
#'  non-significant then accept model 4 and stop. Otherwise, continue to step 3.}
#' \item{Test 3.}{Compare the deviance of models 2 and 1 on 2 d.f. If non-significant then accept
#' model 2 and stop. Otherwise continue to step 4.}
#' \item{Test 4.}{Compare the deviance of models 3 and 1 on 2 d.f. If significant then model 1 cannot
#' be simplified; accept model 1 and stop.  Otherwise continue to step 5}
#' \item{Test 5.}{Compare the deviances of models 5 and 3 on 1 d.f. If significant then model 3 cannot
#' be simplified; accept model 3. Otherwise, accept model 5.  End of procedure}
#' }
#' The result is selection of one of the six models. The FSPA procedure is
#' automatically invoked by the acd(varlist) option for each member of varlist
#'  in the iterative fitting algorithm.

#'
#' @param x input matrix, of dimension nobs x nvars; each row is an observation vector.
#' @param y response variable. Quantitative for family="gaussian", or
#' family="poisson" (non-negative counts). For family="binomial" should be  a
#' variable with two levels (see glm). For family "cox" it must be a Surv object
#' containing  2 columns.
#' @param weights observation weights. Default is NULL which assigns a weight of 1 to each observation
#' @param offset a vector of length \code{nobs} that is included in the linear
#' predictor. Useful for the \code{"poisson"} family (e.g. log of exposure time).
#' Default is \code{NULL} which assigns an offset  of 0 to each observation.
#' If supplied, then values must also be supplied to the \code{predict} function.
#' @param cycles maximum number of iteration cycles. Default is 5
#' @param scale If the values of the variable are too large or too small,
#' the reported results of fractional polynomial (fp) may be difficult to
#' interpret. Scaling can be done automatically or by directly specifying the
#' scaling values so that the magnitude of the x values are not too large.
#' By default scaling factors are estimated by the program (see Details section below).
#' Rather than letting the mfpa automatically choose the scaling factors,
#' you may specify scale factors for each variable.If scaling is not required
#' set \code{scale = 1}.
#' @param shift fractional polynomials are only defined for positive term
#' variables. By default, fp will assume that the variable x is positive and
#' attempt to compute fractional powers of x. If the positive value assumption is
#' incorrect, user-supplied shifting factors will be used to make x positive,
#' otherwise, estimation of shifting factors will be performed. if shifting is
#' not required, set \code{shift = 0}.
#' @param df a vector of values (or single value) that sets the degrees of
#' freedom for each predictor. The df (not counting the intercept) is twice the
#' degree of the FP. For example, an FP2 has 4 df, while FP3 has 6 df. A single
#' value can be specified and the df for all predictors is taken to be that
#' specific value. The program can overrides default df based on the number of
#' distinct (unique) values as follows: 2-3 distinct values are assigned df = 1 (linear),
#' 4-5 distinct values are assigned df = min(2, default) and >=6 distinct values
#' are assigned default df.
#' @param center logical; defines the centering of the covariates. The default is mean
#' centering, except for binary covariates, where the covariate is centered using the lower
#' of the two distinct values of the covariate. see Details section below.
#' @param family A character string representing a `glm()` family object as well as Cox models. For more
#' information, see Details section below.
#' @param criterion a criterion used to select variables and FP models of
#' different degrees. Default is p-values in which the user can specify
#' the significance level (or use default) for variable and function selection
#' (see select and alpha options below). If the user select the BIC or AIC
#' criterion then the program would ignore the nominal significance levels and
#' select variable and functions using BIC or AIC.
#' @param select sets the nominal significance levels for variable selection by backward
#' elimination. A variable is dropped if its removal causes a non-significant increase in
#' deviance. The rules for select are the same as those for df (see above). The
#' default nominal significance level is 0.05 for all variable. Setting the
#' nominal significance level to be 1 for a given variable forces it into the
#' model, leaving others to be selected or not. Using the default selection level
#'  of 1 for all variables forces them all into the model.
#' @param alpha sets the significance levels for testing between FP models of different degrees.
#' The rules for alpha are the same as those for df option (see above). The
#' default nominal significance level is 0.05 for all variables.Example:
#' alpha =0.01 specifies that all variables have an FP selection level of 0.01.
#' @param keep keep one or more variables in the model. The selection level for
#' these variables will be set to 1 when the criterion is pvalue.Important for
#' BIC or AIC criterion since it is possible to set nominal significance level
#' of 1 in pvalue criterion.
#' @param xorder determines the order of entry of the covariates into
#' the model-selection algorithm. The default is \code{ascending}, which enters them in
#' decreasing order of significance in a multiple linear regression
#' (most significant first). \code{Descending} places them in reverse significance order,
#' whereas \code{original} respects the original order of the covariates.
#' @param powers sets the permitted FP powers for all covariates.Default is
#' powers = (-2, -1, -0.5, 0, 0.5, 1, 2, 3), where 0 means natural logarithm.
#' @param ties a character string specifying the method for tie handling. If
#' there are no tied death times all the methods are equivalent. Default is
#' Breslow method. This argument is used for Cox models only and has no
#' effect for other model families. See 'coxph' for details.
#' @param strata variables to be used for stratification. A new factor,
#' whose levels are all possible combinations of the variables supplied as
#' arguments will be created. Default is NULL and Cox model without stratification would be fitted.
#' See \code{"coxph"} from the \code{"survival"} package
#' @param nocenter an optional list of values. See \code{"coxph"} from the \code{"survival"} package
#' @param acdx a vector of names of continuous variables that undergo
#' approximate cumulative distribution (ACD) transformation.
#' Creates the ACD transformation of each variable of \code{acdx}. It also invokes the FSPA
#' to determine the best-fitting FP1(p1, p2) model (see details section). For a
#' given continuous predictor \code{x}, depending on the values of \code{select}
#'  and \code{alpha}, mfpa simplifies the FP1(p1, p2) model to select one of
#'  the six sub-models described in details section below. The variable
#'  representing the ACD transformation of \code{x} is named \code{A(x)}.
#' @param ftest logical;for normal error models with small samples, critical points
#' from the \code{F distribution} can be used instead of \code{Chi-Square distribution}. Default uses
#' the later. This argument is used for Gaussian models only and has no effect for
#'  other model families.
#' @param verbose logical; run in verbose mode (default FALSE)
#' @references
#' Royston, P. 2014. \emph{A smooth covariate rank transformation for use in
#' regression models with asigmoid dose-response function.  Stata Journal 14(2): 329-341.}\cr
#' Royston, P. and Sauerbrei, W., 2016. \emph{mfpa: Extension of mfp using the
#' ACD covariate transformation for enhanced parametric multivariable modeling. The Stata Journal, 16(1), pp.72-87.}
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
  method <- match.arg(ties)
  # capture extra arguments
  # see survival for extra arguments
  # extraArgs <- list(...)
  # if (length(extraArgs)) {
  #   controlargs <- names(formals(coxph.control)) #legal arg names
  #   indx <- pmatch(names(extraArgs), controlargs, nomatch=0L)
  #   if (any(indx==0L))
  #     stop(gettextf("Argument %s not matched",
  #                   names(extraArgs)[indx==0L]), domain = NA)
  # }
  ### Need to do this first so defaults in call can be satisfied
  # set default FP powers proposed by Royston and Sauerbrei (2008)
  if (is.null(powers)) powers <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
  # check dimension of x
  np <- dim(x)
  nobs <- as.integer(np[1])
  nvars <- as.integer(np[2])
  # Assert that x must be a matrix
  if (is.null(np)) stop("x must be a matrix with 1 or more columns")
  # Assert that x must have column names
  vnames <- colnames(x)
  if (is.null(vnames)) stop("The column names of x is null")
  # Transform factor variables to dummies if any using model.matrix
  # Convert factors to dummy variables if it exists...take it to the main function
  # x = model.matrix(as.formula(paste("~", paste(colnames(x), collapse="+"))),
  #                  data = as.data.frame(x))
  x <- model.matrix(~ . - 1, data.frame(x)) # entered df, alpha, select, etc will not be equal to ncol(x) if factor variables exist if f
  # we can compare ncol of x original and from model.matrix if different then some variables where
  # converted to dummies so we need to replicate the df, scale, alpha etc for those variables
  # Check for missing data in x
  if (anyNA(x)) stop("Missing data detected in x.You must eliminate")
  # check weights
  if (!is.null(weights)) {
    if (any(weights < 0)) stop("negative weights not allowed")
    if (length(weights) != nobs) {
      stop(paste("number of elements in weights (", length(weights), ") not equal
                 to the number of rows of x (", nobs, ")", sep = ""))
    }
  }
  # Check offset-offset is also required in cox see coxph.fit() github
  if (!is.null(offset) && length(offset) != nobs) {
    stop(gettextf("number of offsets is %d should equal %d (number of observations)", length(offset), nobs), domain = NA)
  }
  ## define weights and offset when NULL...confirm if default weights are also equal to 1 in cox
  if (is.null(weights)) weights <- rep.int(1, nobs)
  if (is.null(offset)) offset <- rep.int(0, nobs)
  # alpha: p-values for testing between FP models; must lie between 0 and 1
  if (any(alpha > 1) || any(alpha < 0)) stop("Any alpha must not be  > 1 or < 0 ")
  # If length (alpha) = 1 then we replicate supplied alpha otherwise we check
  # whether the length of alpha is equal to ncol(x)
  if (length(alpha) == 1) {
    alpha.list <- rep(alpha, nvars)
  } else {
    if (length(alpha) != nvars) stop(gettextf("number of alpha is %d should equal to 1 or %d (number of variables)", length(alpha), nvars), domain = NA)
    alpha.list <- alpha
  }
  # Order variables based on likelihood ratio test

  # select: the nominal p-values (significance levels) for variable selection by
  # backward elimination; must lie between 0 and 1
  if (any(select > 1) || any(select < 0)) stop("Any select must not be  > 1 or < 0")
  if (length(select) == 1) {
    select.list <- rep(select, nvars)
  } else {
    if (length(select) != nvars) {
      stop(gettextf("number of select is %d should equal to 1 or %d (number of variables)", length(select), nvars), domain = NA)
    }
    select.list <- select
  }
  # =============================================================================
  # keep option...finish this part
  # =============================================================================
  if (!is.null(keep)) {
    if (!all(keep %in% colnames(x))) {
      stop("The names of variables to force into the model is not a subset of columns of x")
    }
  }
  # df: sets the df for each predictor. df=4: FP model with maximum permitted
  # degree m=2 (default), df=2: FP model with maximum permitted degree m=1:
  # df=1: Linear FP model
  # =============================================================================
  # df section: write a function that summarizes this part
  # =============================================================================
  if (length(df) == 1) {
    # check whether df = 1, 2, 4, 6,... even numbers after 1
    if (df != 1) {
      # apply modulus e.g 2, 4, 6, etc (even numbers) has modulus of 0 i.e no remainder when divided by 2
      if (df %% 2 != 0) stop("The df is ", df, " must be either 1 = linear or  2m = FPm where m is the degree of FP")
      # we assign variables different df based on number of unique values
      df.list <- df.assign(x = x, df.default = df)
    } else {
      df.list <- rep(df, nvars)
    }
  } else {
    # the length of df must be equal to the number of predictors including dummy variables
    if (length(df) != nvars) stop(gettextf("The length of df is %d should equal be to 1 or %d (number of variables)", length(df), nvars), domain = NA)
    # check whether the supplied df are correct
    if (any(df != 1)) {
      index1 <- which(df != 1)
      if (any(df[index1] %% 2 != 0)) stop("The df of ", paste0(df[index1][which(df[index1] %% 2 != 0)], sep = " "), " is incorrect. Must be either 1 = linear or  2m = FPm where m is the degree of FP")
    }
    # Assert that variables with <=3 unique values must have df = 1
    nux <- apply(x, 2, function(x) length(unique(x)))
    indexx <- which(nux <= 3)
    if (any(df[indexx] != 1)) stop("The df of variable ", paste0(vnames[indexx][df[indexx] != 1], sep = " "), " must be 1 because it has fewer than 3 unique values")
    df.list <- df
  }
  # ===============================================================================
  # adjustment factors for shifting x to positive values
  if (is.null(shift)) {
    shift <- apply(x, 2, function(x) shift.factor(x))
    x <- sweep(x, 2, shift, "+") # shift to positive
  } else {
    if (length(shift) != nvars) {
      stop(gettextf("number of shift is %d should equal to %d (number of variables)", length(shift), nvars), domain = NA)
    }
    x <- sweep(x, 2, shift, "+")
  }
  # scale: specifies the predictors to be scaled see scalefn() for details
  if (is.null(scale)) {
    scale.list <- apply(x, 2, scalefn)
    x <- sweep(x, 2, scale.list, "/")
  } else {
    if (length(scale) != nvars) stop(gettextf("number of scale is %d should equal to NULL or %d (number of variables)", length(scale), nvars), domain = NA)
    scale.list <- scale
    x <- sweep(x, 2, scale.list, "/")
  }
  # center variables
  if (length(center) == 1) {
    center.list <- rep(center, nvars)
  } else {
    if (length(center) != nvars) {
      stop(gettextf("number of center is %d should equal to 1 or %d (number of variables)", length(center), nvars), domain = NA)
    }
  }
  # logical indicator specifying whether a variable has acd or not. Default is false
  if (is.null(acdx)) {
    acdx <- rep(F, nvars)
  } else {
    # Get rid of duplicates names
    acdx <- unique(acdx)
    if (any(!(acdx %in% vnames))) stop("The variable ", paste0(setdiff(acdx, vnames), sep = ","), " in acdx is not in x")
    # The variables that undergo acd transformation have acdx = T
    acdx <- replace(rep(F, nvars), which(vnames %in% acdx), rep(T, length(acdx)))
    # Reset acdx to false for acd variables with unique values <5
  }
  #-----------------------------------------------------------------------------
  # Cox specific setup
  #-----------------------------------------------------------------------------
  if (family == "cox") {
    # Make sure y is a survival object
    if (!is.Surv(y)) stop("Response must be a survival object")
    # y is a matrix and its rows must be equal to number of observations
    n <- nrow(y)
    if (n != nobs) stop(paste("number of observations in y (", n, ") not equal to the number of rows of x (", nobs, ")", sep = ""))
    # We consider only right censoring. Other options like left, interval, counting,
    # interval2 and mstate implemented in survival package will be considered in the future
    type <- attr(y, "type")
    if (type != "right") stop(paste("Cox FP model doesn't support \"", type, "\" survival data", sep = ""))
    # set defaults for coxph.fit: 1-control, 2-rownames, 3- init, 4-strata
    control <- coxph.control() # coxph.control(...)
    rownames <- row.names(x)
    # ===============================================================================
    # stratification
    # one strata works well but more than one strata fails. check the code in coxph()
    # ===============================================================================
    istrata <- strata
    if (!is.null(strata)) {
      # check whether variables considered for stratification are in x
      if (!all(strata %in% colnames(x))) {
        stop("The names of variables for stratifications is not a subset of columns of x")
      }
      # More checks needed if this strata is working correctly
      istrata <- survival::strata(x[, strata], shortlabel = T)
      # Drop the variable(s) in x used for stratification
      x <- x[, -c(which(colnames(x) %in% strata)), drop = F]
    }
    # strata(GBSG$htreat,GBSG$menostat)
    # Glm specific setup
  } else {
    # Treat one-column matrix of response as vector...combine glm and surv y
    if (is.matrix(y)) {
      dimy <- ncol(y)
      if (dimy == 1) {
        y <- drop(y)
        # Assert that y is a vector not a matrix with two or more columns.
      } else {
        stop(paste0("y must be a vector or a matrix with one column not ", dimy, "columns", sep = ""))
      }
    }
    nrowy <- length(y)
    if (nrowy != nobs) {
      stop(paste("number of observations in y (", nrowy, ") not equal to the number of rows of x (", nobs, ")", sep = ""))
    }
  }
  # Restrict F test to Gaussian response
  if (ftest && family != "gaussian") {
    stop(paste("F test not suitable for \"", family, "\" family. Set ftest = F to use Chi-square instead.", sep = ""))
  }


  # AIC or BIC
  # if(criterion=="AIC"||criterion=="BIC")

  # run mfp
  if (family == "cox") {
    fit <- mfp.fit(
      x = x, y = y, weights = weights, offset = offset, cycles = cycles,
      scale = scale.list, shift = shift, df = df.list, keep,
      center = center.list, criterion = criterion, xorder = xorder,
      powers = powers, family = "cox", method = method,
      select = select.list, alpha = alpha.list, strata = istrata,
      ftest = ftest, verbose = verbose, control = control,
      nocenter = nocenter, rownames = rownames, acdx = acdx
    )
    # if (is.character(fit$fit)){
    #   fit <- list(fail = fit)
    #   class(fit)<- c("mfpa", "coxph")
    # }else{
    class(fit) <- c("mfpa", "coxph")
    # }
    # add wald test in order to use summary.coxph(). this calculation is
    # obtained from coxph() source code
    nabeta <- !is.na(fit$coefficients)
    temp <- fit$coefficients[nabeta]
    fit$wald.test <- coxph.wtest(
      fit$var[nabeta, nabeta], temp,
      control$toler.chol
    )$test
  } else {
    fit <- mfp.fit(
      x = x, y = y, weights = weights, offset = offset, cycles = cycles,
      scale = scale.list, shift = shift, df = df.list, keep,
      center = center.list, criterion = criterion, xorder = xorder,
      powers = powers, family = family, method = method,
      select = select.list, alpha = alpha.list, strata = istrata,
      ftest = ftest, verbose = verbose, control = control,
      nocenter = nocenter, rownames = rownames, acdx = acdx
    )
    class(fit) <- c("mfpa", "glm", "lm")
  }
  fit$call <- cl
  fit
}
