#' Multivariable fractional polynomial models with extensions 
#'
#' Selects the multivariable fractional polynomial (FP) model that best predicts
#' the outcome variable. It also has the ability to model a sigmoid relationship
#' between `x` and an outcome variable `y` using the approximate cumulative 
#' distribution (ACD) transformation proposed by Royston (2016).
#' This function provides two interfaces for input data: one for inputting 
#' data matrix `x` and  outcome vector `y` directly and the other for using a
#' `formula` object together with a data.frame `data`. Both interfaces are
#' equivalent in terms of functionality.
#' 
#' @section Brief summary of FPs:
#' 
#' In the following we denote fractional polynomials for a variable \eqn{x} by 
#' increasing complexity as either FP1 or FP2. In this example, 
#' \eqn{FP2(p1, p2)} is the most flexible FP transformation where 
#' \deqn{FP2(p1, p2) = \beta_1 x^{p1} + \beta_2 x^{p2}.}
#' The (fractional) powers \eqn{p1} and \eqn{p2} are taken from a set
#' of allowed powers, usually {-2, -1, -0.5, 0, 0.5, 1, 2, 3} where the power
#' 0 indicates the natural logarithm. The best FP2 is then estimated by a 
#' closed testing procedure that seeks the best combination from all 36 pairs of
#' powers \eqn{(p1, p2)}. Functions that only involve a single power of
#' the variable are denoted as FP1, i.e. 
#' \deqn{FP1(p1) = \beta_1 x^{p1}.}
#' For details see e.g. Sauerbrei et al (2006). 
#'
#' @section Details on `family` option:
#'
#' `mfp2()` supports the family object as used by [stats::glm()]. The built in
#' families are specified via a character string. `mfp2(..., family="binomial")` 
#' fits a logistic regression model, while `mfp2(..., family="gaussian")`
#' fits a linear regression (ordinary least squares) model.
#'
#' For Cox models, the response should preferably be a `Surv` object,
#' created by the [survival::Surv()] function, and the `family = "cox"`. 
#' Only right-censored data are currently supported. To fit stratified Cox
#' models, the `strata` option can be used, or alternatively `strata` terms 
#' can be included in the model formula when using the formula interface 
#' `mfp2.formula`. 
#'
#' @section Details on shifting, scaling, centering:
#' 
#' Fractional polynomials are defined only for positive variables due to the 
#' use of logarithms and other powers. Thus, `mfp2()` estimates shifts for 
#' each variables to ensure positivity or assumes that the variables are 
#' already positive when computing fractional powers of the input variables
#' in case that shifting is disabled manually. 
#' 
#' If the values of the variables are too large or too small, it is important to
#' conduct variable scaling to reduce the chances of numerical underflow or 
#' overflow which can lead to inaccuracies and difficulties in estimating the
#' model. Scaling can be done automatically or by directly specifying the
#' scaling values so that the magnitude of the `x` values are not too extreme.
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
#' the revised constant \eqn{\beta'_0} or baseline hazard function in a Cox
#' model retains a meaningful interpretation.
#' 
#' So in brief: shifting is required to make input values positive, scaling
#' helps to bring the values to a reasonable range. Both operations are 
#' conducted before estimating the FP powers for an input variable. 
#' Centering, however, is done after estimating the FP functions each variable.
#' Centering before estimating the FP powers may result in different powers and
#' should be avoided. Also see [transform_vector_fp()] for some more details.
#' 
#' @section Details on the `subset` argument: 
#' Note that subsetting occurs after data pre-processing, but before model
#' selection and fit. In detail, when the option `subset` is used and scale, 
#' shift or centering values are to be estimated, then `mfp2()` first estimates
#' these using the full dataset (no subsetting), then applies subsetting, then 
#' proceeds to do model selection and fit on the subset of the data specified. 
#' 
#' Therefore, subsetting in `mfp2()` is not equivalent to subsetting the data
#' before passing it to `mfp2()`, and thus cannot be used to implement e.g. 
#' cross-validation or to remove `NA`. This should be done by the caller 
#' beforehand. However, it does allow to use the same data pre-processing 
#' for different subsets of the data. An example usecase is when separate
#' models are to be estimated for women and men in the dataset, but a common
#' data pre-processing should be applied. In this case the `subset` option 
#' can be used to restrict model selection to either women or men, but share
#' data processing between the two models. 
#'
#' @section Details on  approximate cumulative distribution transformation:
#' 
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
#' It maintains the family-wise type 1 error 
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
#' @section Details on model specification using a `formula`:
#' `mfp2` supports model specifications using two different interfaces: one
#' which allows passing of the data matrix `x` and outcome vector `y` directly
#' (as done in e.g. [stats::glm.fit()] or `glmnet`) and another which conforms
#' to the formula interface used by many commonly used R modelling functions 
#' such as [stats::glm()] or [survival::coxph()]. 
#' 
#' Both interfaces are equivalent in terms of possible fitted models, only the
#' details of specification differ. In the standard interface all details 
#' regarding FP-transformations are given as vectors. In the formula interface
#' all details are specified using special `fp` terms. These support the 
#' specification of degrees of freedom (`df`), confidence level for 
#' variable selection (`select`), confidence level for functional form 
#' selection (`alpha`), shift values (`shift`), scale values (`scale`), 
#' centering (`center`) and the ACD-transformation (`acd`). Values specified 
#' by these functions override the values specified as defaults and passed to 
#' the `mfp2()` function. 
#' 
#' The formula may also contain `strata` terms to fit stratified Cox models, or
#' an `offset` term to specify a model offset.
#' 
#' Note that for a formula using `.`, such as `y ~ .` the function may not
#' fit a linear model, but may also do selection of variable and functional 
#' forms using FP-transformations, depending on the default settings of `df`, 
#' `select` and `alpha` passed as arguments to `mfp2()`. 
#' For example, using `y ~ .` with default settings for all other arguments
#' means that `mfp2()` will apply FP transformation with 4 df to all 
#' continuous variables and use alpha equal to 0.05 to select functional forms, 
#' and use the selection algorithm with significance level 0.05 for all 
#' variables. 
#' 
#' @section Compatibility with `mfp` package: 
#' `mfp2` is an extension of the `mfp` package and can be used to reproduce
#' the results from a model fitted by `mfp`. Since both packages provide
#' an implementation of the MFP algorithm, both packages use functions of the 
#' same name. Thus, if you load both packages by a call to `library` there will
#' be namespace conflicts and only the functions of the package loaded later
#' will be working properly. 
#'
#' @param x for `mfp2.default`: `x` is an input matrix of dimensions 
#' nobs x nvars. Each row is an observation vector.
#' @param y for `mfp2.default`: `y` is a vector for the response variable.
#' For `family = "binomial"` it should be  a vector with two levels (see 
#' [stats::glm()]). For `family = "cox"` it must be a [survival::Surv()] object
#' containing 2 columns.
#' @param formula for `mfp2.formula`: an object of class `formula`: a symbolic 
#' description of the model to be fitted. Special `fp` terms can be used to 
#' define fp-transformations. The details of model specification are given 
#' under ‘Details’.
#' @param data for `mfp2.formula`: a `data.frame` which contains all variables
#' specified in `formula`.
#' @param weights a vector of observation weights of length nobs. 
#' Default is `NULL` which assigns a weight of 1 to each observation.
#' @param offset a vector of length nobs that is included in the linear
#' predictor. Useful for the poisson family (e.g. log of exposure time).
#' Default is `NULL` which assigns an offset  of 0 to each observation.
#' If supplied, then values must also be supplied to the `predict()` function.
#' @param cycles an integer, maximum number of iteration cycles. Default is 5.
#' @param scale a numeric vector of length nvars or single numeric specifying 
#' scaling factors. If a single numeric, then the value will be replicated as
#' necessary. The formula interface `mfp2.formula` only supports single numeric 
#' input to set a default value, individual values can be set using `fp` terms
#' in the `formula` input. 
#' Default is `NULL` which lets the program estimate the scaling factors 
#' (see Details section). If scaling is not required set `scale = 1` to disable 
#' it.
#' @param shift a numeric vector of length nvars or a single numeric specifying
#' shift terms. If a single numeric, then the value will be replicated as
#' necessary. The formula interface `mfp2.formula` only supports single numeric 
#' input to set a default value, individual values can be set using `fp` terms
#' in the `formula` input.
#' Default is `NULL` which lets the program estimate the shifts
#' (see Details section). If shifting is not required, set `shift = 0` to 
#' disable it.
#' @param df a numeric vector of length nvars or a single numeric that sets the 
#' (default) degrees of freedom (df) for each predictor. If a single numeric, 
#' then the value will be replicated as necessary. The formula interface
#' `mfp2.formula` only supports single numeric input to set a default value, 
#' individual values can be set using `fp` terms in the `formula` input. 
#' The df (not counting the intercept) are twice the degree of a fractional 
#' polynomial (FP). For example, an FP2 has 4 df, while FPm has 2*m df. 
#' The program overrides default df based on the number of distinct (unique) 
#' values for a variable as follows: 
#' 2-3 distinct values are assigned `df = 1` (linear), 4-5 distinct values are
#' assigned `df = min(2, default)` and >= 6 distinct values are assigned  
#' `df = default`. 
#' @param center a logical determining whether variables are centered before 
#' model fit. The default `TRUE` implies mean centering, except for binary 
#' covariates, where the covariate is centered using the lower of the two 
#' distinct values of the covariate. See Details section below.
#' @param subset subset	an optional vector specifying a subset of observations
#' to be used in the fitting process. Default is `NULL` and all observations are 
#' used. See Details below. 
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
#' sets the nominal significance levels for variable selection on each predictor
#' by backward elimination. If a single numeric, then the value will be replicated
#' as necessary. The formula interface `mfp2.formula` only supports single numeric 
#' input to set a default value, individual values can be set using `fp` terms
#' in the `formula` input. The default nominal significance level is 0.05 
#' for all variables. Setting the nominal significance level to be 1 for  
#' certain variables forces them into the model, leaving all other variables
#' to be selected. 
#' @param alpha a numeric vector of length nvars or a single numeric that 
#' sets the significance levels for testing between FP models of 
#' different degrees. If a single numeric, then the value will be replicated
#' as necessary. The formula interface `mfp2.formula` only supports single numeric 
#' input to set a default value, individual values can be set using `fp` terms
#' in the `formula` input. The default nominal significance level is 0.05 for all 
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
#' @param powers a named list of numeric values that sets the permitted FP 
#' powers for each covariate. The default is NULL, and each covariate is assigned
#' `powers = c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)`, where 0 means the natural 
#' logarithm. Powers are sorted before further
#' processing in the program. If some variables are not assigned powers, the 
#' default powers will be assigned. The formula interface offers two options
#' for supplying powers: through the 'powers' argument and the 'fp()' function. 
#' So, if the user supplies powers in both options for a certain variable, the 
#' powers supplied through 'fp()' will be given preference.
#' @param ties a character string specifying the method for tie handling in 
#' Cox regression. If there are no tied death times all the methods are 
#' equivalent. Default is the Breslow method. This argument is used for Cox 
#' models only and has no effect on other model families. 
#' See [survival::coxph()] for details.
#' @param strata a numeric vector or matrix of variables that define strata
#' to be used for stratification in a Cox model. A new factor, whose levels are 
#' all possible combinations of the variables supplied will be created. 
#' Default is `NULL` and a Cox model without stratification would be fitted. 
#' See [survival::coxph()] for details. 
#' @param nocenter a numeric vector with a list of values for fitting Cox 
#' models. See [survival::coxph()] for details.
#' @param acdx a numeric vector of names of continuous variables to undergo 
#' the approximate cumulative distribution (ACD) transformation.
#' It also invokes the function-selection procedure to determine the 
#' best-fitting FP1(p1, p2) model (see Details section). Not present in the 
#' formula interface `mfp2.formula` and to be set using `fp` terms in the 
#' `formula` input.
#' The variable representing the ACD transformation of `x` is named `A(x)`.
#' @param ftest a logical; for normal error models with small samples, critical 
#' points from the F-distribution can be used instead of Chi-Square 
#' distribution. Default `FALSE` uses the latter. This argument is used for 
#' Gaussian models only and has no effect for other model families.
#' @param control a list object with parameters controlling model fit details. 
#' Returned by either [stats::glm.control()] or [survival::coxph.control()]. 
#' Default is `NULL` to use default parameters for the given model class. 
#' @param verbose a logical; run in verbose mode.
#' @param ... not used.
#' @examples
#'
#' # Gaussian model
#' data("prostate")
#' x = as.matrix(prostate[,2:8])
#' y = as.numeric(prostate$lpsa)
#' # default interface
#' fit1 = mfp2(x, y, verbose = FALSE)
#' fit1$fp_terms
#' fracplot(fit1) # generate plots
#' coef(fit1)
#' print(fit1)
#' summary(fit1)
#' # formula interface
#' fit1b = mfp2(lpsa ~ fp(age) + fp(svi, df = 1) + fp(pgg45) + fp(cavol) + fp(weight) +
#' fp(bph) + fp(cp), data = prostate)
#' 
#' # logistic regression model
#' data("pima")
#' xx <- as.matrix(pima[, 2:9])
#' yy <- as.vector(pima$y)
#' fit2 <- mfp2(xx, yy, family = "binomial", verbose = FALSE)
#' fit2$fp_terms
#' fracplot(fit2)
#' 
#' # Cox regression model
#' data("gbsg")
#' # create dummy variable for grade using ordinal coding
#' gbsg <- create_dummy_variables(gbsg, var_ordinal = "grade", drop_variables = TRUE)
#' xd <- as.matrix(gbsg[, -c(1, 6, 10, 11)])
#' yd <- survival::Surv(gbsg$rectime, gbsg$censrec)
#' # fit mfp and keep hormon in the model
#' fit3 <- mfp2(xd, yd, family = "cox", keep = "hormon", verbose = FALSE)
#' fit3$fp_terms
#' fracplot(fit3)
#' 
#' @return 
#' `mfp2()` returns an object of class inheriting from `glm` or `copxh`, 
#' depending on the `family` parameter. 
#' 
#' The function `summary()` (i.e. [summary.mfp2()]) can be used to obtain or
#' print a summary of the results.  
#' The generic accessor function `coef()` can be used to extract the vector of 
#' coefficients from the fitted model object. 
#' The generic `predict()` can be used to obtain predictions from the fitted 
#' model object.
#' 
#' An object of class `mfp2` is a list containing all entries as for `glm`
#' or `coxph`, and in addition the following entries:  
#' \itemize{
#' \item{convergence_mfp: }{logical value indicating convergence of mfp algorithm.}
#' \item{fp_terms: }{a data.frame with information on fractional polynomial 
#' terms.}
#' \item{transformations: }{a data.frame with information on shifting, scaling
#' and centering for all variables.}   
#' \item{fp_powers: }{a list with all powers of fractional polynomial terms. 
#' Each entry of the list is named according to the transformation of the 
#' variable.}
#' \item{acd: }{a vector with information for which variables the acd 
#' transformation was applied.}
#' \item{x_original: }{the scaled and shifted input matrix but without
#' transformations.}
#' \item{y: }{the original outcome variable.}
#' \item{x: }{the final transformed input matrix used to fit the final model.}
#' \item{call_mfp: }{the call to the `mfp2()` function.}
#' \item{family_string: }{the family stored as character string.}
#' }
#' The `mfp2` object may contain further information depending on family.
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
#' [summary.mfp2()], [coef.mfp2()], [predict.mfp2()], [fp()]
#' 
#' @export
mfp2 <- function(x, ...){
  UseMethod("mfp2", x)
}

#' @describeIn mfp2 Default method using input matrix `x` and outcome vector `y`.
#' @export
mfp2.default <- function(x, 
                         y, 
                         weights = NULL, 
                         offset = NULL, 
                         cycles = 5,
                         scale = NULL, 
                         shift = NULL, 
                         df = 4, 
                         center = TRUE,
                         subset = NULL,
                         family = c("gaussian", "poisson", "binomial", "cox"),
                         criterion = c("pvalue", "aic", "bic"),
                         select = 0.05, 
                         alpha = 0.05,
                         keep = NULL,
                         xorder = c("ascending", "descending", "original"),
                         powers = NULL,
                         ties = c("breslow", "efron", "exact"),
                         strata = NULL,
                         nocenter = NULL,
                         acdx = NULL,
                         ftest = FALSE,
                         control = NULL, 
                         verbose = TRUE,
                         ...) {
  
  # this function prepares everything for fitting the actual mfp2 model
  
  cl <- match.call()
  
  # match arguments ------------------------------------------------------------
  criterion <- match.arg(criterion)
  xorder <- match.arg(xorder)
  family <- match.arg(family)
  ties <- match.arg(ties)
  
  # assertions -----------------------------------------------------------------
  # assert that x is a matrix
  if (!is.matrix(x))
    stop("! x must be a matrix", call. = F)
  
  # assert that x must not contain character values
  if (any(is.character(x)))
    stop("! x contains characters values.\n",
         "i Please convert categorical variables to dummy variables.", 
         call. = F)
  
  # check dimension of x
  np <- dim(x)
  nobs <- as.integer(np[1])
  nvars <- as.integer(np[2])
  
  # assert that x is a matrix
  if (is.null(np)) {
      stop("! The dimensions of x must not be missing.\n",
           "i Please make sure that x is a matrix with at least one row and column.", 
           call. = FALSE)
  }
  # assert that x must have column names
  vnames <- colnames(x)
  if (is.null(vnames)) 
    stop("! The column names of x must not be missing.\n",
         "i Please set column names for x.",
         call. = FALSE)

  # assert that x has no missing data
  if (anyNA(x)) 
    stop("! x must not contain any NA (missing data).\n",   
         "i Please remove any missing data before passing x to this function.")
  
 # assert that subset must be a vector and does not contain negative values
  if (!is.null(subset)) {
    if (!is.vector(subset))
      stop(sprintf("! Subset must not be of class %s.\n", 
                   paste0(class(subset), collapse = ", ")), 
           "i Please convert subset to a vector.", 
           call. = FALSE)
    
    if (any(subset < 0))
      stop("! Subset must not contain negative values.", call. = FALSE)
  }  
  
  # assert that weights are positive and of appropriate dimensions
  if (!is.null(weights)) {
    if (any(weights < 0)) {
        stop("! Weights must not be negative.", call. = FALSE)
    }
    if (length(weights) != nobs) {
      stop("! The number of observations (rows in x) and weights must match.\n", 
           sprintf("i The number of rows in x is %d, but the number of elements in weights is %d.", 
                   nobs, length(weights)), 
           call. = FALSE)
    }
  }
  
  # assert that the length of offset must be equal to the number of observations
  if (!is.null(offset)) {
      if (length(offset) != nobs) {
          stop("! The number of observations (rows in x) and offset must match.\n", 
               sprintf("i The number of rows in x is %d, but the number of elements in offset is %d.", 
                       nobs, length(offset)))
      }
  }
  
  # assert that alpha must be between 0 and 1
  if (any(alpha > 1) || any(alpha < 0)) {
      stop("! alpha must not be < 0 or > 1.")
  }
  
  # assert length of alpha 
  if (length(alpha) != 1 && length(alpha) != nvars) {
      stop("! alpha must be a single number, or the number of variables (columns in x) and alpha must match.\n", 
           sprintf("i The number of variables in x is %d, but the number of elements in alpha is %d.", 
                   nvars, length(alpha)))
  }
  
  # assert that select must be between 0 and 1
  if (any(select > 1) || any(select < 0)) {
      stop("! select must not be < 0 or > 1.")
  }
  # assert length of select 
  if (length(select) != 1 && length(select) != nvars) {
      stop("! select must be a single number, or the number of variables (columns in x) and select must match.\n", 
           sprintf("i The number of variables in x is %d, but the number of elements in select is %d.", 
                   nvars, length(select)))
  }
  
  # assert that keep is a subset of x
  if (!is.null(keep)) {
      if (!all(keep %in% vnames)) {
          warning("i The set of variables named in keep is not a subset of the variables in x.\n", 
                  "i mfp2() continues with the intersection of keep and colnames(x).")
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
      if (!all(acdx %in% vnames)) {
          warning("i The set of variables named in acdx is not a subset of the variables in x.\n", 
                  "i mfp2() continues with the intersection of acdx and colnames(x).")
      }
  }
  
  # assert df is positive
  if (any(df <= 0)) {
      stop("! df must not be 0 or negative.\n", 
           "i All df must be either 1 (linear) or 2m, where m is the degree of FP.")
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
              "i mfp2() reverts to use Chi-square instead.")
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
               "i Currently only right censoring is supported by mfp2().")
      }
      
      # assert numeric
      if (is.factor(strata)) {
        strata <- as.numeric(strata)
      }
      
      if (!is.null(strata)) {
          # assert stratification factors are of correct length
          if (is.vector(strata)) {
            strata_len = length(strata)
          } else strata_len = nrow(strata)
          if (strata_len != nrow(x)) {
              stop("! The length of stratification factor(s) and the number of observations in x must match.\n")
          }
      }
  } else {
      # assert type of y
      if (is.matrix(y)||is.data.frame(y)) {
          stop(sprintf("! Outcome y must not be of class %s.", 
                       paste0(class(y), collapse = ", ")), 
               "i Please convert y to a vector.", call. = FALSE)
      }
      if (length(y) != nobs) {
          stop("! Number of observations in y and x must match.", 
               sprintf("i The number of observations in y is %d, but the number of observations in x is %d.", 
                       length(y), nobs))
      }
  }
  
  # set defaults ---------------------------------------------------------------
  # Assert that powers must be NULL or list when provided
  if (!is.null(powers) && !is.list(powers))
    stop("Powers must be a list")
  
  # set defaults ---------------------------------------------------------------
    # default FP powers proposed by Royston and Altman (1994)
    powx <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
    power_list <- setNames(replicate(nvars, powx, simplify = FALSE), vnames)
    
    # deal with user supplied powers
    if (!is.null(powers)) {
      if (length(powers) != sum(names(powers) != "", na.rm = TRUE))
         stop(" All the powers supplied in the argument must have names",
           call. = FALSE)
    
    # check the names of supplied powers
    dd <- which(!names(powers) %in% vnames)
    if (length(dd) !=0)
      stop(" The names of all powers must be in the column names of x.\n",
           sprintf("i This applies to the following powers: %s.", 
                   paste0(names(powers)[dd], collapse = ", ")), call. = FALSE)
    
    if (!all(sapply(powers, is.numeric)))
      stop("All elements of powers must be numeric", call. = FALSE)
    
    # sort powers
    powers <- lapply(powers, function(v) sort(v))
    
    # modify the default powers, some variables assigned default powers
    power_list <- modifyList(power_list, Filter(Negate(is.null), powers))
    
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
  
  keep <- intersect(keep, vnames)
  # convert acdx to logical vector
  if (is.null(acdx)) {
      acdx <- rep(F, nvars)
  } else {
      # Get rid of duplicates names
      acdx <- unique(acdx)
      acdx <- intersect(acdx, vnames)
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
        warning("i For any variable with fewer than 4 unique values the df are set to 1 (linear) by mfp2().\n", 
                sprintf("i This applies to the following variables: %s.", 
                        paste0(vnames[index & df != 1], collapse = ", ")))
        df[index] <- 1
    }
    df.list <- df
  }

  # data preparation -----------------------------------------------------------
  # shift and scale 
  x <- sweep(x, 2, shift, "+")
  x <- sweep(x, 2, scale, "/")
  
  # stratification for cox model
  istrata <- strata
  if (family == "cox" && !is.null(strata)) {
      istrata <- survival::strata(strata, shortlabel = TRUE)
      # convert strata to integers as conducted by coxph
      istrata <- as.integer(istrata)
  }

    # data subsetting ----------------------------------------------------------
  if (!is.null(subset)) {
    # check the sample size if it is too small, likely to occur after subsetting
    if (length(subset) < 5)
      stop("! The length of subset is too small (<5) to fit an mfp model.", 
           sprintf("i The number of observations is %d", length(subset)), 
           call. = FALSE)
    
    x <- x[subset, , drop = FALSE]
    y <- if (family != "cox") {
      y[subset]
    } else y[subset, , drop = FALSE]
    
    weights <- weights[subset]
    offset <- offset[subset]
    istrata <- istrata[subset]
  }
                 
  # fit model ------------------------------------------------------------------
  fit <- fit_mfp(
      x = x, y = y, 
      weights = weights, offset = offset, cycles = cycles,
      scale = scale, shift = shift, df = df.list, center = center, 
      family = family, criterion = criterion, select = select, alpha = alpha, 
      keep = keep, xorder = xorder, powers = power_list, 
      method = ties, strata = istrata, nocenter = nocenter,
      acdx = acdx, ftest = ftest, 
      control = control, 
      verbose = verbose
  )
  
  # add additional information to fitted object
  # original mfp2 call
  fit$call_mfp <- cl
  fit$family_string <- family
  fit$offset <- offset
  
  fit
}

#' @describeIn mfp2 Provides formula interface for `mfp2`.
#' @export
mfp2.formula <- function(formula, 
                         data, 
                         weights = NULL, 
                         offset = NULL, 
                         cycles = 5,
                         scale = NULL, 
                         shift = NULL, 
                         df = 4, 
                         center = TRUE,
                         subset = NULL,
                         family = c("gaussian", "poisson", "binomial", "cox"),
                         criterion = c("pvalue", "aic", "bic"),
                         select = 0.05, 
                         alpha = 0.05,
                         keep = NULL,
                         xorder = c("ascending", "descending", "original"),
                         powers = NULL,
                         ties = c("breslow", "efron", "exact"),
                         strata = NULL,
                         nocenter = NULL,
                         ftest = FALSE,
                         control = NULL,
                         verbose = TRUE,
                         ...) {
  # capture the call
  call <- match.call()
  family <- match.arg(family)
 
  # assert that data must be provided
  if (missing(data))
    stop("! data argument is missing.\n",
         "i An input data.frame is required for the use of mfp2.",
         call. = FALSE)
  
  # assert that data has column names
  if (is.null(colnames(data)))
    stop("! data must have column names.\n",
         "i Please set column names.")
  
  # assert that a formula must be provided
  if (missing(formula)) 
    stop("! formula is missing.", call. = FALSE)
  
  if (!inherits(formula, "formula"))
    stop("method is only for formula objects", call. = FALSE)
  
  # assert length of df, alpha, select, center, shift, scale, acdx equal to one
  if (length(df) != 1)
    stop("! df must be a single numeric.", 
         "i Use the fp() function to set different df values in the input formula.",
         call. = FALSE)
  
  if (length(alpha) != 1)
    stop("! alpha must be a single numeric.", 
         "i Use the fp() function to set different alpha values in the input formula.",
         call. = FALSE)
  
  if (length(select) != 1)
    stop("! select must be a single numeric.", 
         "i Use the fp() function to set different select values in the input formula.",
         call. = FALSE)
  
  if (!is.null(scale) && length(scale) != 1)
    stop("! scale must be a single numeric or NULL.",
         "i Use the fp() function to set different scaling factors in the input formula.", 
         call. = FALSE)
  
  if (length(center) != 1)
    stop("! center must be a single logical value.", 
         "i Use the fp() function to set different center values in the input formula.",
         call. = FALSE)
  
  if (!is.null(shift) && length(shift) != 1)
    stop("! shift must be a single numeric.", 
         "i Use the fp() function to set different shift values in the input formula.",
         call. = FALSE)
  
  if(!is.null(powers) && !is.list(powers))
    stop(" Powers must be a named list or set it to NULL", call. = FALSE)
  
  # model.frame preserves the attributes of the data unlike model.matrix
  mf <- stats::model.frame(formula, data = data, drop.unused.levels = TRUE)
  
  # check whether no predictor exist in the model i.e y~1: 
  labels <-  attr(terms(mf), "term.labels")
  if (length(labels)==0)
    stop("No predictors are provided for model fitting.\n At least one predictor is required", call. = FALSE)
  
  # stratification for Cox models ----------------------------------------------

  # strata not allowed in the formula if the family is not cox
  specials <- "strata"
  terms_formula <- terms(formula, specials = specials, data = data)

  # remember position of strata variables to drop from input data if necessary
  terms_drop <- NULL
  if (!is.null(attr(terms_formula,"specials")$strata)) {
    if (family == "cox") {
      
      # check whether strata is both in the formula and in the argument
      if (!is.null(call$strata))
        warning("i strata appear both in the formula and as an input argument.\n",
                "i The information in the formula is used and the input argument ignored.",
                call. = FALSE)
      
      # untangle the terms for strata as in coxph
      # this function returns the strata names, e.g "strata(x1)" 
      # and its position in the terms when outcome is excluded
      stemp <- survival::untangle.specials(terms_formula,
                                           special = "strata", 
                                           order = 1)
      
      if (length(stemp$vars) == 1) {
        # only one strata exists in the formula   
        strata <- mf[[stemp$vars]]
      } else {
        # more than one strata exists in the formula
        strata <- mf[, stemp$vars]
      }
      
      # extract the position of strata variables in the terms to be dropped
      terms_drop <- stemp$terms
    } else {
      stop("! strata are only allowed for Cox models.\n", 
           "i Please remove any strata terms from the model formula.",
           call. = FALSE)
    }
  }
  
  # drop strata variables if necessary before using model.matrix()
  if (!is.null(terms_drop)) 
    terms_model <- terms_formula[-terms_drop]
  else terms_model <- terms_formula
  
  # offset ---------------------------------------------------------------------
  term_offset <- attr(terms_formula, "offset")
  if (!is.null(term_offset) && length(term_offset) > 1)
    stop("! Only one offset in the formula is allowed.", call. = FALSE)
  
  # check whether offset is both in the formula and as an argument.
  if (!is.null(term_offset) && !is.null(call$offset)) {
    warning("i Offset appears both in the formula and as an input argument.\n", 
            "i The information in the model formula is used and the input argument is ignored.", 
            call. = FALSE)
    offset <- as.vector(model.offset(mf))
  }
  
  # data preparation -----------------------------------------------------------
  
  y <- model.extract(mf, "response")
  if (family != "cox") 
    y <- as.numeric(y)
  
  x <- model.matrix(terms_model, mf)
  # remove intercept if necessary
  # intercept is coded as entry 0 in attribute assigned by model.matrix
  # intercept is always the first column
  if (0 %in% attr(x, "assign")) 
    x <- x[, -1, drop = FALSE]
  
  nx <- ncol(x) 
  names_x <- colnames(x)
  
  # select variables that undergo fp transformation and extract their attributes
  fp_pos <- grep("fp(.*)", colnames(mf))
  
  if (length(fp_pos) > 0) {
    fp_data <- mf[, fp_pos, drop = FALSE]
    
    # extract names of the variables that undergo fp transformation
    fp_vars <- unname(sapply(fp_data, function(v) attr(v, "name")))
    
    # check for variables used more than once in fp() function 
    fp_vars_duplicates <- fp_vars[duplicated(fp_vars)]
    if (length(fp_vars_duplicates) != 0)
      stop("! Variables should be used only once in the fp() within the formula.\n", 
           sprintf("i The following variable(s) are duplicated in fp() function: %s.", 
                   paste0(fp_vars_duplicates, collapse = ", ")), 
           call. = FALSE)
    
    # check for variables used in fp() as well as other parts of the formula
    vars_duplicates <- which(colnames(mf) %in% fp_vars)
    if (length(vars_duplicates) != 0)
      stop("! Variables used in the fp() should not be included in other parts of the formula.\n", 
           sprintf("i This applies to the following variable(s): %s.", 
                   paste0(colnames(mf)[vars_duplicates], collapse = ", ")), 
           call. = FALSE)
    
    # replace names such as fp(x1) by real name "x1" in the x matrix
    names_x <- replace(names_x, grep("fp(.*)", names_x), fp_vars)
    colnames(x) <- names_x
  }
  
  # call default method---------------------------------------------------------
  
  # if fp() is not used in the formula, it reduces to mfp2.default() 
  df_list <-  setNames(as.list(assign_df(x = x, df_default = df)), names_x)
  
  # scaling
  if (is.null(scale)) {
    scale_list <- setNames(as.list(apply(x, 2, find_scale_factor)), names_x)
  } else {
    scale_list <- setNames(rep(list(scale), nx), names_x)
    
  }
  
  # shifting
  if (is.null(shift)) {
    shift_list <- setNames(as.list(apply(x, 2, find_shift_factor)), names_x)
  } else {
    shift_list <- setNames(rep(list(shift), nx), names_x)
  }
  
  # center, alpha, select and acd
  center_list <- setNames(rep(list(center), nx), names_x)
  alpha_list <- setNames(rep(list(alpha), nx), names_x)
  select_list <- setNames(rep(list(select), nx), names_x)
  acdx_list <- setNames(rep(list(FALSE), nx), names_x)
  
  
  # default FP powers proposed by Royston and Altman (1994)
  powx <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
  power_list <- setNames(replicate(nx, powx, simplify = FALSE),names_x)
  
  # deal with user supplied powers in the argument not through fp()
  if (!is.null(powers)) {
    if (length(powers) != sum(names(powers) != "", na.rm = TRUE))
      stop(" All the powers supplied in the argument must have names",
           call. = FALSE)
    
    # check the names of supplied powers
    dd <- which(!names(powers) %in% names_x)
    if (length(dd) !=0)
      stop(" The names of all powers must be in the column names of x.\n",
           sprintf("i This applies to the following powers: %s.", 
                   paste0(names(powers)[dd], collapse = ", ")), call. = FALSE)
    
    # sort powers
    powers <- lapply(powers, function(v) sort(v))
    
    # modify the default powers, some variables assigned default powers
    power_list <- modifyList(power_list, Filter(Negate(is.null), powers))
    
  }
  
  # if fp() is used in the formula
  if (length(fp_pos) != 0) {
    # modify the default parameters based on the user inputs
    df_list <- modifyList(df_list, 
                          setNames(lapply(fp_data, attr, "df"), fp_vars))
    scale_list <- modifyList(scale_list, 
                             Filter(Negate(is.null),setNames(lapply(fp_data, attr,
                                                                    "scale"), fp_vars)))
    shift_list <- modifyList(shift_list, 
                             Filter(Negate(is.null),setNames(lapply(fp_data, 
                                                                    attr, "shift"), fp_vars)))
    center_list <- modifyList(center_list,
                              setNames(lapply(fp_data, attr, "center"), fp_vars))
    alpha_list <- modifyList(alpha_list, 
                             setNames(lapply(fp_data, attr, "alpha"), fp_vars))
    select_list <- modifyList(select_list, 
                              setNames(lapply(fp_data, attr, "select"), fp_vars))
    acdx_list <- modifyList(acdx_list, 
                            setNames(lapply(fp_data, attr, "acd"), fp_vars))
    
    # We give preference to powers supplied in the fp() function over the power argument.
    powerx <- Filter(Negate(is.null),setNames(lapply(fp_data, attr, "powers"), fp_vars))
    
    nax <- intersect(names(powerx), names(powers))
    if (length(nax)!= 0)
      warning("i Powers are specified in both the `fp()` function within\n the formula and as an argument.                 The argument term ignored.\n", 
              sprintf("i This applies to the following variables: %s.", 
                      paste0(nax, collapse = ", ")), call. = FALSE)
    
    power_list <- modifyList(power_list, powerx)
    
  }
  
  # acd requires variable names or NULL 
  acdx_vector <- unlist(acdx_list)
  if (sum(acdx_vector) == 0) 
    acdx_vector <- NULL
  else
    acdx_vector <- names(acdx_vector[acdx_vector])

  mfp2.default(x = x, 
               y = y, 
               weights = weights, 
               offset = offset, 
               cycles = cycles,
               scale = unlist(scale_list), 
               shift = unlist(shift_list), 
               df = unlist(df_list), 
               center = unlist(center_list),
               subset = subset,
               family = family,
               criterion = criterion,
               select = unlist(select_list), 
               alpha = unlist(alpha_list),
               keep = keep,
               xorder = xorder,
               powers = power_list,
               ties = ties,
               strata = strata,
               nocenter = nocenter,
               acdx = acdx_vector,
               ftest = ftest,
               control = control,
               verbose = verbose
  )
}

#' Extract coefficients from object of class `mfp2`
#' 
#' This function is a method for the generic [stats::coef()] function for 
#' objects of class `mfp2`. 
#' 
#' @param object an object of class `mfp2`, usually, a result of a call to
#' [mfp2()].
#' @param ... not used.
#' 
#' @return 
#' Named numeric vector of coefficients extracted from the model `object`.
#' 
#' @export
coef.mfp2 <- function(object, ...) {
  object$coefficients
}

#' Summarizing `mfp2` model fits
#' 
#' This function is a method for the generic [base::summary()] function for
#' objects of class `mfp2`.
#' 
#' @param object an object of class `mfp2`, usually, a result of a call to
#' [mfp2()].
#' @param ... further arguments passed to the summary functions for `glm()` 
#' ([stats::summary.glm()], i.e. families supported by `glm()`) or `coxph()` 
#' ([survival::summary.coxph()], if `object$family = "cox"`).
#' 
#' @return 
#' An object returned from [stats::summary.glm()] or
#' [survival::summary.coxph()], depending on the family parameter of `object`.
#' 
#' @seealso 
#' [mfp2()], [stats::glm()], [stats::summary.glm()], [survival::coxph()],
#' [survival::summary.coxph()]
#' 
#' @export
summary.mfp2 <- function(object, ...) {
  NextMethod("summary", object)
}

#' Print method for objects of class `mfp2`
#' 
#' Enhances printing by information on data processing and fractional 
#' polynomials.
#' 
#' @param x `mfp2` object to be printed.
#' @param ... passed to `print` methods of underlying model class. A useful 
#' option as the `digits` argument, indicating printed digits.
#' 
#' @export
print.mfp2 <- function(x,
                       ...) {
  # shift and scaling factors with centering values
  cat("Shifting, Scaling and Centering of covariates", "\n")
  print.data.frame(x$transformations)
  cat("\n")
  
  # Final MFP Powers
  cat("Final Multivariable Fractional Polynomial for y", "\n")

  # Exclude 'acd' column if all its elements are false
  fp_terms <- x$fp_terms
  if(all(!fp_terms["acd"])) {
    fp_terms <- fp_terms[, !names(fp_terms) %in% "acd"]
  }
  
  print.data.frame(fp_terms)
  cat("\n")
  
  cat(sprintf("MFP algorithm convergence: %s\n", x$convergence_mfp))
  
  # print model object using underlying print function
  NextMethod("print", x)
}

#' Helper to assign attributes to a variable undergoing FP-transformation
#' 
#' Used in formula interface to `mfp2()`.
#' 
#' @param x a vector representing a continuous variable undergoing 
#' fp-transformation.
#' @param df,alpha,select,shift,scale,center,acdx See [mfp2::mfp2()]) for details. 
#' @param powers a vector of powers to be evaluated for `x`. Default is `NULL` 
#' and `powers = c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)` will be used.
#' @param ... used in alias `fp2` to pass arguments.
#' 
#' @examples
#'
#' xr = 1:10
#' fp(xr)
#' fp2(xr)
#' @return 
#' The vector `x` with new attributes relevant for fp-transformation. All 
#' arguments passed to this function will be stored as attributes. 
#' 
#' @export
fp <- function(x, 
               df = 4, 
               alpha = 0.05,
               select = 0.05, 
               shift = NULL, 
               scale = NULL,
               center = TRUE, 
               acdx = FALSE, 
               powers = NULL) {
  
  name <- deparse(substitute(x))
  
  # Assert that a factor variable must not be subjected to fp transformation
  if (is.factor(x))
    stop(name," is a factor variable and should not be passed to the fp() function.")
  
  attr(x, "df") <- df
  attr(x, "alpha") <- alpha
  attr(x, "select") <- select
  attr(x, "shift") <- shift
  attr(x, "scale") <- scale
  attr(x, "center") <- center
  attr(x, "acd") <- acdx
  attr(x, "powers") <- powers
  attr(x, "name") <- name
  
  x
}

#' @describeIn fp Alias for `fp()` - use in formula when both `mfp` and `mfp2` are loaded to avoid name shadowing.
#' @export
fp2 <- function(...) {
  fp(...)
}

#' Helper function to extract selected variables from fitted `mfp2` object
#' 
#' Simply extracts all variables for which not all powers are estimated to 
#' be `NA`. The names refer to the original names in the dataset and do not
#' include transformations.
#' 
#' @param object fitted `mfp2` object.
#' 
#' @examples
#'
#' # Gaussian model
#' data("prostate")
#' x = as.matrix(prostate[,2:8])
#' y = as.numeric(prostate$lpsa)
#' # default interface
#' fit = mfp2(x, y, verbose = FALSE)
#' get_selected_variable_names(fit)
#' 
#' @return 
#' Character vector of names, ordered as defined by `xorder` in [mfp2()].
#' 
#' @export
get_selected_variable_names <- function(object) {
  nms <- rownames(object$fp_terms)
  nms[object$fp_terms[, "selected"]]
}

#' Helper to assign degrees of freedom
#' 
#' Determine the number of unique values in a variable. To be used in [mfp2()].
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
#' @examples
#' x <- matrix(1:100, nrow = 10)
#' assign_df(x)
#'
#' @return 
#' Vector of length `ncol(x)` with degrees of freedom for each variable in `x`.
#' 
#' @export
assign_df <- function(x, 
                      df_default = 4) {
  
  # default degrees of freedom for each variable
  df <- rep(df_default, ncol(x))
  # unique values of each column of x
  nu <- apply(x, 2, function(x) length(unique(x)))
  
  index1 <- which(nu <= 3)
  index2 <- which(nu >= 4 & nu <= 5)
  # set df to 1
  if (length(index1) != 0) {
    df[index1] <- 1
  }
  if (length(index2) != 0) {
    df[index2] <- 2
  }
  
  return(df)
}
