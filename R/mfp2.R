#' Multivariable Fractional Polynomial Models with Extensions 
#'
#' Selects the multivariable fractional polynomial (MFP) model that best predicts
#' the outcome variable. It also supports approximate cumulative 
#' distribution (ACD) transformations for continuous variables, which enable 
#' modeling of sigmoid relationships between predictors `x` and the outcome `y` 
#' (Royston, 2014; Royston and Sauerbrei, 2016). In addition, it implements 
#' spike-at-zero (SAZ) modeling, a method for appropriately handling 
#' semi-continuous predictors with a non-negligible proportion of zero values 
#' (Becher et al., 2012). The function provides two interfaces for input data: 
#' one for supplying the data matrix `x` and the outcome `y` directly, where `y` 
#' may be a numeric vector for continuous, count, or binary outcomes, a two-level
#' factor or a two-column grouped-count matrix for binomial outcomes, or a
#' `Surv` object for Cox proportional hazards models; and another for using a `formula` object
#' together with a dataframe `data`. Both interfaces are equivalent in functionality.
#' 
#' @section Brief summary of fractional polynomials (FPs):
#' 
#' Fractional polynomials (FPs) provide a flexible framework for modeling
#' nonlinear relationships between a continuous predictor \eqn{x} and an outcome.
#' In general, we denote an FP of degree \eqn{m} as \eqn{FPm(p_1, \dots, p_m)}, 
#' where \eqn{p_1, \dots, p_m} are the selected powers and \eqn{m \ge 1}.
#' 
#' The most commonly used cases are:
#' 
#' * **FP1:** a single-term transformation, \eqn{FP1(p_1) = \beta_1 x^{p_1}}, 
#'   representing the simplest FP model.
#' 
#' * **FP2:** a two-term transformation, \eqn{FP2(p_1, p_2) = \beta_1 x^{p_1} + \beta_2 x^{p_2}}, 
#'   where \eqn{p_1 \ne p_2}, providing greater flexibility to capture nonlinear 
#'   effects.
#' 
#' When \eqn{p_1 = p_2} (repeated powers), the FP2 model is defined as
#' \deqn{FP2(p_1, p_2) = \beta_1 x^{p_1} + \beta_2 x^{p_1} \log(x).}
#' 
#' The powers \eqn{p_1} and \eqn{p_2} are usually chosen from a predefined set
#' \eqn{S = \{-2, -1, -0.5, 0, 0.5, 1, 2, 3\}}, where a power of 0 indicates
#' the natural logarithm. The best FP2 model is then selected using a closed
#' testing procedure that evaluates all 36 pairs of powers \eqn{(p_1, p_2)}.
#' 
#' For further details, see Sauerbrei et al. (2006) and Royston and Sauerbrei (2008). 
#' For the effects of influential points on FP functions, see Sauerbrei et al. (2023).
#'
#' @section Details on `family` option:
#'
#' `mfp2()` supports the `family` argument as used by \code{stats::glm()}. 
#' Families can be specified either as a character string or as a function 
#' returning a GLM family object (e.g., `stats::gaussian(link = "identity")` 
#' or `stats::binomial(link = "logit")`).  
#' 
#' Only the following families are supported at the moment: `"gaussian"`, 
#' `"binomial"`, `"poisson"`, and `"cox"`.  
#' 
#' Examples with character strings:  
#' `mfp2(..., family = "binomial")` fits a logistic regression model.  
#' `mfp2(..., family = "gaussian")` fits a linear regression model using ordinary least squares.  
#' 
#' Examples with family functions and custom links:  
#' `mfp2(..., family = gaussian(link = "log"))` fits a linear model with a log 
#' link. `mfp2(..., family = binomial(link = "probit"))` fits a logistic 
#' regression model with a probit link.  
#' 
#' For Cox proportional hazards models, the response should be a `Surv` object 
#' (created with `survival::Surv()`), and the `family` argument should be set to 
#' `"cox"`. Only right-censored data are currently supported.  
#' Stratified Cox models can be specified using the `strata` argument, or by 
#' including `strata` terms in the model formula when using the formula interface 
#' `mfp2.formula`.
#'
#' @section Details on shifting, scaling, centering:
#' 
#' Fractional polynomials are defined only for positive variables due to the 
#' use of logarithms and other powers. Thus, `mfp2()` estimates shifts for 
#' each variable to ensure positivity or assumes that the variables are 
#' already positive when computing fractional powers of the input variables
#' if shifting is disabled manually. 
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
#' After the final FP powers are estimated, the program backscales \eqn{x'} to 
#' the original scale \eqn{x}, ensuring that the final regression coefficients 
#' are expressed in the original scale of the data. The FP transformation of 
#' \eqn{x} is centered on the mean of the observed values of \eqn{x}. For example,
#' for the FP1 model \eqn{\beta_0 + \beta_1x^p},the actual model fitted by the 
#' software would be \eqn{\beta'_0 + \beta'_1(x^p-mean(x^p))}. This approach 
#' ensures that the revised constant \eqn{\beta'_0} or baseline hazard function 
#' in a Cox model retains a meaningful interpretation.
#' 
#' So in brief: shifting is required to make input values positive, scaling
#' helps to bring the values to a reasonable range. Both operations are 
#' conducted before estimating the FP powers for an input variable. 
#' Centering, however, is done after estimating the FP functions for each 
#' variable.
#' 
#' Additionally, any variable marked as `zero`, `catzero`, `spike`, or with 
#' degrees of freedom equal to 1 (`df = 1`) has its shift automatically set to 
#' zero. This ensures that variables for which transformations are unnecessary,
#' such as linear terms, are not artificially shifted, and that the zero
#' component remains at zero for variables with spike-at-zero or catzero
#' handling.
#' 
#' Centering before estimating the FP powers may result in different powers and
#' should be avoided. Also see \code{transform_vector_fp()} for some more details.
#' 
#' @section Details on the `subset` argument:
#' Subsetting in `mfp2()` occurs after shift and scale parameters have been
#' estimated and applied, but before model selection and fitting. Specifically,
#' if `subset` is used and shift or scale parameters need to be estimated,
#' `mfp2()` first estimates these parameters using the full dataset, without
#' applying the subset. The subset is then applied before model selection and
#' fitting are performed on the selected observations.
#'
#' Consequently, using `subset` in `mfp2()` is not equivalent to subsetting the
#' data before calling the function. The `subset` argument should not be used for
#' tasks such as cross-validation or removal of missing values; those should be
#' handled explicitly before fitting the model. The argument is primarily useful
#' when the same shift and scale transformation should be applied across multiple
#' analysis subsets. For example, separate models may be fitted for women and men
#' while estimating the shift and scale parameters from the same full dataset.
#' In that setting, `subset` restricts model selection and fitting to the chosen
#' group while retaining consistent shift and scale parameters across groups.
#' 
#' 
#' @section Handling of Categorical Variables:
#'
#' The \code{mfp2.formula()} method expands categorical variables through
#' \code{\link[stats]{model.matrix}}, ensuring that the model is fitted on
#' a fully numeric design matrix. This affects how predictors are represented
#' internally and how arguments such as \code{keep} and \code{select} are
#' interpreted.
#'
#' ## Unordered Factors
#'
#' Unordered categorical variables (created using \code{factor()}) are encoded
#' using the default treatment contrasts (\code{\link[stats]{contr.treatment}}).
#' The first level of the factor is treated as the reference category, and a
#' dummy variable is created for each of the remaining levels. For example:
#'
#' \preformatted{
#' x <- factor(c("A", "B", "C"))
#' model.matrix(~ x)
#' # (Intercept) xB xC
#' }
#'
#' When the argument \code{keep} includes the name of a categorical variable,
#' all dummy variables derived from that factor (e.g., \code{xB}, \code{xC})
#' are automatically retained in the model. The corresponding entries in
#' \code{select} are internally set to 1 to prevent their removal during
#' variable selection.
#'
#' ## Ordered Factors
#'
#' Ordered categorical variables (created using \code{ordered()}) are, by default,
#' encoded using polynomial contrasts (\code{\link[stats]{contr.poly}}), which
#' produce orthogonal contrasts representing trend effects across ordered levels
#' (e.g., \code{x.L}, \code{x.Q}, \code{x.C}).
#'
#' If the analysis requires dummy coding that reflects ordinal thresholds-
#' for instance, comparing level A versus the rest, or levels A and B versus
#' the rest—the user must define and assign an appropriate contrast matrix
#' before calling \code{mfp2.formula()}. This can be done using the
#' \code{contr.cumulative()} function provided in this package:
#'
#' \preformatted{
#' data$x <- ordered(data$x, levels = c("A", "B", "C", "D"))
#' contrasts(data$x) <- contr.cumulative(levels(data$x))
#' fit <- mfp2(y ~ x, data = data)
#' }
#'
#' This approach generates cumulative dummy variables that preserve the
#' ordinal nature of the variable while allowing for interpretable
#' threshold-type effects. Users should ensure that the contrast specification
#' reflects the intended interpretation prior to fitting the model.
#'
#' ## The \code{mfp2.default()} Method
#'
#' The default method \code{mfp2.default()} does not perform any automatic
#' expansion of factor variables. It assumes that the input design matrix
#' \code{x} consists solely of numeric predictors. Therefore, when calling
#' \code{mfp2.default()} directly, the user must create any required dummy or
#' contrast-coded variables manually before model fitting.
#'
#' @section Details on  approximate cumulative distribution transformation:
#' 
#' The approximate cumulative distribution (ACD) transformation is a method
#' to model continuous covariates flexibly in regression models. Instead of
#' including the raw variable \eqn{X} directly, a smooth function approximating
#' its empirical cumulative distribution function (ecdf) is used.
#'
#' **Method**
#'
#' Let \eqn{x_1, \dots, x_n} be a sample from the distribution of \eqn{X}.
#' The ACD transformation proceeds in three steps:
#'
#' 1. **Inverse normal transformation of ranks:**
#' Compute the rank of each \eqn{x_i} in the sample and transform it using the
#' standard normal inverse CDF (probit):
#' \deqn{z_i = \Phi^{-1} \Big( \frac{\text{rank}(x_i) - 0.5}{n} \Big),}
#' where \eqn{\Phi^{-1}} is the inverse standard normal CDF. This maps the
#' empirical distribution of \eqn{X} to approximately standard normal values.
#'
#' 2. **Power-linear approximation:**
#' Fit a one-term fractional polynomial regression of \eqn{z_i} on a shifted
#' and powered version of \eqn{X}:
#' \deqn{\hat{z}_i = \hat{\beta}_0 + \hat{\beta}_1 (x_i + \text{shift})^p,}
#' where \eqn{p} is the best-fitting power, and \eqn{\text{shift}} ensures all
#' values are positive if necessary. Ordinary least squares is used to estimate
#' \eqn{\hat{\beta}_0} and \eqn{\hat{\beta}_1}. A power of 0 corresponds to a
#' log transformation.
#'
#' 3. **Back-transformation to the (0,1) scale:**
#' The fitted values \eqn{\hat{z}_i} are transformed back to the interval (0,1)
#' using the standard normal CDF:
#' \deqn{\text{ACD}(x_i) = a_i = \Phi(\hat{z}_i) = \Phi(\hat{\beta}_0 + \hat{\beta}_1 (x_i + \text{shift})^p),}
#' producing a smooth approximation of the ecdf.
#'
#' **Interpretation**
#'
#' The ACD transformation maps \eqn{X} to an approximately uniform scale on (0,1). 
#' When a regression model is specified as 
#' \eqn{E(Y) = \beta_0 + \beta_1 \text{ACD}(X)}, the expected value of \eqn{Y} 
#' changes smoothly from \eqn{\beta_0} at \eqn{\text{ACD}(X) \approx 0} 
#' (corresponding to the minimum of \eqn{X}) to \eqn{\beta_0 + \beta_1} at 
#' \eqn{\text{ACD}(X) \approx 1} (corresponding to the maximum of \eqn{X}). 
#' 
#' Because the mapping from \eqn{X} to \eqn{\text{ACD}(X)} is S-shaped - changing
#' slowly at the extremes and more rapidly in the central range - the resulting 
#' relationship between \eqn{Y} and \eqn{X} is typically nonlinear and 
#' sigmoid-shape. This sigmoid shape cannot be achieved by standard fractional 
#' polynomial (FP) functions.
#' 
#' Intuitively, the ACD transformation "stretches" the middle of \eqn{X}'s 
#' distribution while compressing the tails. A linear effect of 
#' \eqn{\text{ACD}(X)} in the regression thus produces slow changes in \eqn{Y} 
#' for extreme values of \eqn{X} and faster changes for central values, 
#' generating a characteristic sigmoid curve.
#' 
#' Details of the precise definition and some possible uses of the ACD 
#' transformation in a univariate context are given by Royston (2014).
#'  
#' **FP1 with ACD component**
#'
#' Royston (2014) and Royston and Sauerbrei (2016) describe extending the FP2
#' family by replacing one FP1 term with an ACD-transformed term:
#' \deqn{FP1(p_1, p_2) = \beta_1 x^{p_1} + \beta_2 \text{ACD}(x)^{p_2},}
#' which allows modeling of sigmoid-like effects while maintaining flexibility
#' similar to standard FP2 functions.
#'
#' The powers \eqn{p_1} and \eqn{p_2} are chosen from a predefined set
#' \eqn{S = \{-2, -1, -0.5, 0, 0.5, 1, 2, 3\}}, with 0 corresponding to a
#' logarithmic transformation. All 64 combinations of \eqn{(p_1, p_2)} are
#' considered during model selection.
#'
#' **Model simplification**
#'
#' After fitting, the chosen FP1+ACD function can be simplified via a function
#' selection procedure (FSPA) that evaluates six nested sub-families:
#' \itemize{
#'   \item M1: FP1(p1, p2) (no simplification)
#'   \item M2: FP1(p1, .) (regular FP1 in \eqn{x})
#'   \item M3: FP1(., p2) (FP1 in \eqn{ACD(x)})
#'   \item M4: FP1(1, .) (linear in \eqn{x})
#'   \item M5: FP1(., 1) (linear in \eqn{ACD(x)})
#'   \item M6: null (\eqn{x} omitted)
#' }
#' 
#' **Closed test procedure to choose final model**
#' 
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
#' \item Test 1: Compare the deviances of models 6 and 1 on 4 d.f. 
#' If not significant then stop and omit \eqn{x}, otherwise continue to step 2.
#' \item Test 2: Compare the deviances of models 4 and 1 on 3 d.f. 
#' If not significant then accept model 4 and stop. Otherwise, continue to step 3.
#' \item Test 3: Compare the deviance of models 2 and 1 on 2 d.f. 
#' If not significant then accept model 2 and stop. Otherwise continue to step 4.
#' \item Test 4: Compare the deviance of models 3 and 1 on 2 d.f. 
#' If significant then model 1 cannot be simplified; accept model 1 and stop.  
#' Otherwise continue to step 5.
#' \item Test 5: Compare the deviances of models 5 and 3 on 1 d.f. 
#' If significant then model 3 cannot be simplified; accept model 3. 
#' Otherwise, accept model 5. End of procedure.
#' }
#' The result is the selection of one of the six models. For details see 
#' Royston and Sauerbrei (2016).
#' 
#' **Information criteria**
#'
#' As an alternative criterion, the fit of models M1-M6 can be evaluated using 
#' AIC or BIC. The model with the smallest AIC or BIC is selected.
#' @section Handling Nonpositive Values:
#'
#' The \code{zero_vars}, \code{catzero_vars}, and \code{spike_vars} options
#' provide related mechanisms for covariates with nonpositive values in
#' fractional polynomial models.
#'
#' For variables in \code{zero_vars}, FP transformations are applied only to
#' positive values. Nonpositive values are counted in the zero component and are
#' recoded to zero before transformation. In the formula interface, use
#' \code{fp(x, zero = TRUE)}.
#'
#' For variables in \code{catzero_vars}, \code{zero_vars} handling is also
#' applied and a binary zero-component indicator is added. The indicator is
#' computed after recoding as \code{I(x == 0)}, equivalently
#' \code{I(original x <= 0)}. In the formula interface, use
#' \code{fp(x, catzero = TRUE)}.
#'
#' For variables in \code{spike_vars}, the SAZ algorithm evaluates whether the
#' zero-component indicator and the FP function for the positive component are
#' both needed. \code{spike_vars} implies \code{catzero_vars} and
#' \code{zero_vars} only while the variable remains SAZ-active after eligibility
#' checks. In the formula interface, use \code{fp(x, spike = TRUE)}.
#' 
#' @section Details on spike-at-zero (SAZ) modelling:
#' The spike-at-zero (SAZ) procedure is intended for semi-continuous covariates
#' with a zero component and a positive component. When SAZ is active for a
#' variable, the model considers both a binary zero-component indicator and a
#' fractional-polynomial function for the positive component.
#'
#' SAZ model selection has two stages. Stage 1 selects the functional form for
#' the positive component while retaining the zero-component indicator. For
#' \code{criterion = "pvalue"}, this follows the usual closed-test logic for FP
#' terms, with the zero-component indicator included. For \code{criterion = "aic"}
#' or \code{criterion = "bic"}, the model with the best information criterion is
#' selected.
#'
#' Stage 2 assesses whether the final model should retain both components, only
#' the positive-component FP function, or only the zero-component indicator. For
#' \code{criterion = "pvalue"}, this is based on tests comparing the selected
#' full SAZ model with reduced models. For \code{criterion = "aic"} or
#' \code{criterion = "bic"}, the corresponding reduced models are compared by
#' information criterion.
#'
#' If the zero-component indicator is removed, the final model contains only the
#' selected FP function for the positive component. If the FP function is
#' removed, the final model contains only the zero-component indicator.
#' 
#' @section Spike-at-zero eligibility:
#' A variable requested through \code{spike_vars}, or through
#' \code{fp(x, spike = TRUE)} in the formula interface, enters the SAZ algorithm
#' only if both components are sufficiently represented. During this eligibility
#' check, nonpositive values are counted in the zero component and positive
#' values are counted in the positive component.
#'
#' Let \eqn{p_0} denote the zero-component proportion and \eqn{p_+} denote the
#' positive-component proportion. A requested SAZ variable remains eligible only
#' if
#' \deqn{p_0 \geq m}
#' and
#' \deqn{p_+ \geq m,}
#' where \eqn{m =} \code{min_saz_component_prop}.
#'
#' With the default \code{min_saz_component_prop = 0.10}, at least 10 percent
#' of observations must be in each component. Variables failing this check, or
#' variables that are binary, have their spike flag reset to \code{FALSE}.
#' After reset, any \code{zero} or \code{catzero} handling is retained only if
#' it was explicitly requested by the user.
#'
#' For retained SAZ variables, \code{mfp2()} may reduce the maximum FP degrees
#' of freedom when the positive component has few distinct values, using the
#' same principle as ordinary MFP modelling.
#' 
#' @section Formula interface:
#' \code{mfp2()} supports both a matrix interface, \code{mfp2(x, y, ...)}, and
#' a formula interface, \code{mfp2(formula, data, ...)}. The formula method
#' constructs a model frame, processes the predictors, and then calls
#' \code{mfp2.default()} internally.
#'
#' Fractional-polynomial settings for individual predictors can be supplied
#' with \code{fp()} terms in the formula. These settings include degrees of
#' freedom \code{df}, selection level \code{select}, functional-form level
#' \code{alpha}, \code{shift}, \code{scale}, \code{center}, ACD handling
#' \code{acdx}, nonpositive-value handling \code{zero}, zero-component
#' indicators \code{catzero}, spike-at-zero assessment \code{spike}, and
#' candidate powers \code{powers}. Settings supplied inside \code{fp()} take
#' precedence over corresponding global defaults. For example,
#' \code{mfp2(y ~ fp(x, df = 2), df = 4, data = dat)} fits \code{x} with
#' \code{df = 2}; the global \code{df = 4} is ignored for that term.
#'
#' The formula interface expands categorical predictors using
#' \code{\link[stats]{model.matrix}}. When \code{keep} names a categorical
#' predictor, all dummy variables generated from that predictor are retained and
#' protected from exclusion during model selection.
#'
#' Variable-name arguments such as \code{keep}, \code{zero_vars},
#' \code{catzero_vars}, \code{spike_vars}, and \code{acdx} must refer to
#' existing column names in \code{x}. Unknown names are treated as errors to
#' avoid silently ignoring misspelled model specifications.
#'
#' The formula may contain \code{strata()} terms for stratified Cox models and
#' \code{offset()} terms for model offsets. For a formula using \code{.}, such
#' as \code{y ~ .}, continuous variables are still subject to the global
#' \code{df}, \code{select}, and \code{alpha} defaults unless overridden inside
#' \code{fp()} terms.

#' @section Compatibility with `mfp` package: 
#' `mfp2` is an extension of the `mfp` package and can be used to reproduce
#' the results from a model fitted by `mfp`. Since both packages implement the 
#' MFP algorithm, they use functions with the same names (e.g `fp()`). Therefore,
#' if you load both packages using a call to `library`, there will
#' be namespace conflicts and only the functions from the package loaded last
#' will work properly. 
#' 
#' @section Convergence and Troubleshooting: 
#' Typically, `mfp2` requires two to five cycles to achieve convergence. Lack of 
#' convergence involves oscillation between two or more models and is extremely
#' rare. If convergence problems occur, consider adjusting the nominal 
#' significance levels for variable selection (`select`) or functional form 
#' selection (`alpha`)
#'
#' @param x for `mfp2.default`: `x` is an input matrix of dimensions 
#' nobs x nvars. Each row is an observation vector.
#' @param y Response variable for `mfp2.default`.
#' For `family = "gaussian"` and `family = "poisson"`, `y` must be a numeric
#' vector with one value per observation. For `family = "binomial"`, `y` may be
#' a numeric 0/1 vector, a two-level factor, or a numeric two-column matrix of
#' grouped binomial counts of the form `cbind(successes, failures)`, as accepted
#' by [stats::glm()] with `family = binomial()`. For `family = "cox"`, `y` must
#' be a [survival::Surv()] object with two columns.
#' @param formula for `mfp2.formula`: an object of class `formula`: a symbolic 
#' description of the model to be fitted. Special `fp` terms can be used to 
#' define fp-transformations. For `family = binomial()`, the left-hand side may
#' be a numeric 0/1 response, a two-level factor response, or a grouped binomial
#' response written as `cbind(successes, failures)`. The details of model
#' specification are given under "Details" section.
#' @param data for `mfp2.formula`: a `data.frame` which contains all variables
#' specified in `formula`.
#' @param weights Optional numeric vector of observation weights. 
#' Default is `NULL` which assigns a weight of 1 to each observation. In the 
#' formula interface, the vector is supplied directly and is not looked up in 
#' `data`.
#' @param offset a vector of length nobs that is included in the linear
#' predictor. Useful for the poisson family (e.g. log of exposure time).
#' Default is `NULL` which assigns an offset  of 0 to each observation.
#' If supplied, then values must also be supplied to the `predict()` function.
#' @param cycles an integer, specifying the maximum number of iteration cycles. 
#' Default is 5.
#' @param scale a numeric vector of length `nvars` or single numeric specifying 
#' scaling factors. If a single numeric, then the value will be replicated as
#' necessary. The formula interface `mfp2.formula` only supports single numeric 
#' input to set a default value, individual values can be set using `fp` terms
#' in the `formula` input. 
#' Default is `NULL` which lets the program estimate the scaling factors 
#' (see Details section). If scaling is not required set `scale = 1` to disable 
#' it. The final regression coefficients are expressed in the original scale of
#' the data.
#' @param shift a numeric vector of length `nvars` or a single numeric specifying
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
#' final model fitting. The default `TRUE` implies mean centering, except for
#' binary covariates, where the covariate is centered using the lower of the two 
#' distinct values of the covariate. See Details section below.
#' @param subset an optional vector specifying a subset of observations
#' to be used in the fitting process. Default is `NULL` and all observations are 
#' used. See Details below. weights/offset/strata must be aligned to original 
#' data rows, even when subset is supplied.
#' @param family Either a character string specifying the model family 
#'   (e.g., "gaussian", "binomial", "poisson", "cox") or a function that 
#'   returns a GLM family object, such as `stats::gaussian(link = "identity")` 
#'   or `stats::binomial(link = "logit")`. For Cox models, only the character 
#'   string `"cox"` is allowed. See \strong{Details on family option}.
#' @param criterion a character string specifying the criterion used to select 
#' variables and FP models of different degrees. 
#' Default is to use p-values in which case the user can specify
#' the nominal significance level (or use default level of 0.05) for variable and
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
#' @param keep a character vector with names of variables to be kept 
#' in the model. In case that `criterion = "pvalue"`, this is equivalent to
#' setting the selection level for the variables in `keep` to 1. 
#' However, this option also keeps the specified variables in the model when 
#' using the BIC or AIC criteria. 
#' @param xorder a string determining the order of entry of the covariates
#' into the model-selection algorithm (backfitting algorithm). The default is 
#' `ascending`, which enters them by ascending p-values, or decreasing order of 
#' significance in a multiple regression (i.e. most significant first).
#' `descending` places them in reverse significance order, whereas 
#' `original` respects the original order in `x`.
#' @param powers a named list of numeric values specifying the set of candidate 
#' FP powers for each covariate. The default is `NULL`, in which case every 
#' covariate is assigned `powers = c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)`, where 
#' `0` denotes the natural logarithm. Powers are sorted before further 
#' processing. If some variables are not explicitly assigned powers, the default 
#' set is used. In the formula interface, powers can be specified either through 
#' the `powers` argument or within the `fp()` function. If both are provided for 
#' the same variable, the specification inside `fp()` takes precedence. Each 
#' variable with `df > 1` must have at least one candidate power other than `1`,
#' because the linear model (power 1) is fitted separately in the closed-test
#' procedure. A single non-unity candidate power is valid; for example,
#' `powers = list(x = 2)` produces the repeated-power FP2 candidate `c(2, 2)`.
#' To restrict a variable to a purely linear effect, set `df = 1`.
#' @param ties a character string specifying the method for tie handling in 
#' Cox regression. If there are no tied death times all the methods are 
#' equivalent. Default is the Breslow method. This argument is used for Cox 
#' models only and has no effect on other model families. 
#' See \code{survival::coxph()} for details.
#' @param strata a numeric vector or matrix of variables that define strata
#' to be used for stratification in a Cox model. A new factor, whose levels are 
#' all possible combinations of the variables supplied will be created. 
#' Default is `NULL` and a Cox model without stratification would be fitted. 
#' See \code{survival::coxph()} for details. 
#' @param nocenter a numeric vector with a list of values for fitting Cox 
#' models. See \code{survival::coxph()} for details.
#' @param acdx a character vector giving the names of continuous variables to 
#' undergo the approximate cumulative distribution (ACD) transformation. Using 
#' this also triggers the function-selection procedure for ACD (FSPA) to 
#' determine the best-fitting FP1(p1, p2) model (see Details). This argument is 
#' not available in the formula interface (`mfp2.formula`), where the user
#' should instead use `fp()`. Within `fp()` terms, the ACD transformation is 
#' specified via the logical argument `acdx` rather than by listing variable 
#' names. The variable representing the ACD transformation of `x` is named `A_x`.
#' @param ftest a logical indicating whether `mfp2` should use critical values 
#' from the F-distribution instead of the Chi-Square distribution for small-sample 
#' normal error models. This affects variable selection, functional form selection, 
#' and spike-at-zero testing when evaluating fractional polynomial terms. Default 
#' is `FALSE`, in which case the Chi-Square distribution is used. This argument 
#' applies only to Gaussian models.
#' @param control a list of parameters controlling the fitting process, as 
#' returned by [stats::glm.control()] or [survival::coxph.control()]. Default 
#' is `NULL`, in which case `mfp2` uses the default settings for the specified 
#' model family.
#' @param zero_vars A character vector specifying variables for which
#'   non-positive values should be treated as the zero component and recoded to
#'   zero before FP transformation. In \code{fp()} terms, use
#'   \code{zero = TRUE}. See \strong{Handling Nonpositive Values}.
#'
#' @param catzero_vars A character vector specifying variables for which a
#'   binary zero-component indicator should be added in addition to zero
#'   handling. \code{catzero_vars} implies \code{zero_vars}. In \code{fp()}
#'   terms, use \code{catzero = TRUE}. See
#'   \strong{Handling Nonpositive Values}.
#'
#' @param spike_vars A character vector specifying variables to be assessed with
#'   the spike-at-zero algorithm. \code{spike_vars} implies
#'   \code{catzero_vars} and \code{zero_vars} only for variables that remain
#'   SAZ-active after eligibility checks. In \code{fp()} terms, use
#'   \code{spike = TRUE}. See \strong{Spike-at-zero eligibility} and
#'   \strong{Handling Nonpositive Values}.
#' @param min_saz_component_prop
#'   Numeric in \eqn{(0, 0.5)}. Minimum required proportion in each component of
#'   a spike-at-zero covariate: the zero component and the positive
#'   continuous component. Default \code{0.10}. A requested spike-at-zero
#'   variable is retained for SAZ modelling only if both the zero
#'   proportion and the positive-observation proportion are at least this value.
#'   Variables failing this eligibility check have their spike flag reset to
#'   \code{FALSE}; the resulting treatment then depends on whether the user also
#'   specified \code{zero_vars} or \code{catzero_vars}.
#' @param force_max_fp For \code{mfp2.default()}, either a single logical value
#'   or a logical vector of length \code{nvars}. A single value is replicated
#'   across all predictors. If \code{TRUE} for a predictor, the most complex
#'   fractional-polynomial function allowed by that predictor's \code{df} is
#'   forced when using \code{criterion = "aic"} or \code{criterion = "bic"}.
#'   It has no effect when \code{criterion = "pvalue"}, where similar behavior
#'   can be obtained by setting \code{select = 1} and \code{alpha = 1}.
#'   The default is \code{FALSE}. For \code{mfp2.formula()}, set this per
#'   variable inside \code{fp()}.
#' @param verbose Logical specifying whether to print progress messages. Default
#'  is TRUE.
#' @param \dots Additional arguments passed to methods. The current
#'   \code{mfp2.default()} and \code{mfp2.formula()} methods do not require
#'   additional arguments; variable-specific formula options should be supplied
#'   inside \code{fp()} or \code{fp2()} terms.
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
#' summary(fit1)
#' # formula interface
#' fit1b = mfp2(lpsa ~ fp(age) + fp(svi, df = 1) + fp(pgg45) + fp(cavol) + fp(weight) +
#' fp(bph) + fp(cp), data = prostate, verbose = FALSE)
#' 
#' # Handling semi-continuous variables in the default interface
#' # Here, "exposure" is treated as a spike-at-zero variable.
#' # fit_saz <- mfp2(x, y, spike_vars = "exposure")
#' # Formula interface: use logical flags inside fp()
#' # fit_saz_f <- mfp2(y ~ fp(exposure, spike = TRUE), data = dat)
#'
#' @return 
#' The returned object inherits from \code{glm} or \code{coxph} and adds MFP
#' metadata such as selected FP terms, transformations, powers, ACD flags,
#' zero/catzero/spike information, and convergence status.
#' 
#' An object of class `mfp2` is a list containing all entries as for `glm`
#' or `coxph`, and in addition the following entries:  
#' \itemize{
#' \item convergence_mfp: logical value indicating convergence of mfp algorithm.
#' \item fp_terms: a data.frame with information on fractional polynomial 
#' terms.
#' \item transformations: a data.frame with information on shifting, scaling
#' and centering for all variables.  
#' \item fp_powers: a list with all powers of fractional polynomial terms. 
#' Each entry of the list is named according to the transformation of the 
#' variable.
#' \item acd: a vector with information for which variables the acd 
#' transformation was applied.
#' \item x_original: the scaled and shifted input matrix but without
#' transformations.
#' \item y: the original outcome variable.
#' \item x: the final transformed input matrix used to fit the final model.
#' \item call_mfp: the call to the `mfp2()` function.
#' \item family: the family stored as either character string or function.
#' \item family_string: the family stored as character string.
#' \item zero: named logical vector indicating, for each variable, whether only
#'  positive values were transformed.  
#' \item catzero: named logical vector indicating which columns in `x` treated 
#' nonpositive values as zero and included an additional binary indicator in 
#' the model.  
#' \item catzero_list: A list of binary variables created when `catzero` is set to
#' TRUE. Returns `NULL` if `catzero` is FALSE.
#' \item spike_dec: named numeric vector encoding the SAZ stage-2 decision
#' for each variable. Value 1 retains both the continuous FP component and the
#' binary zero indicator. Value 2 retains the continuous FP component only
#' (binary indicator dropped). Value 3 retains the binary zero indicator only
#' (continuous FP component dropped). For non-spike variables the default
#' value is 2 (standard FP, no SAZ).
#' \item has_offset: logical value indicating whether an offset was specified
#' before missing offsets were replaced by zeros internally.
#' }
#' The `mfp2` object may contain further information depending on family.
#' 
#' @references
#' Royston, P. and Sauerbrei, W., 2008. \emph{Multivariable Model - Building: 
#' A Pragmatic Approach to Regression Analysis based on Fractional Polynomials 
#' for Modelling Continuous Variables. John Wiley & Sons.}\cr
#' 
#' Sauerbrei, W., Meier-Hirmer, C., Benner, A. and Royston, P., 2006. 
#' \emph{Multivariable regression model building by using fractional 
#' polynomials: Description of SAS, STATA and R programs. 
#' Comput Stat Data Anal, 50(12): 3464-85.}\cr
#' 
#' Royston, P. 2014. \emph{A smooth covariate rank transformation for use in
#' regression models with a sigmoid dose-response function. 
#' Stata Journal 14(2): 329-341.}\cr
#' 
#' Royston, P. and Sauerbrei, W., 2016. \emph{mfpa: Extension of mfp using the
#' ACD covariate transformation for enhanced parametric multivariable modeling. 
#' The Stata Journal, 16(1), pp.72-87.}\cr
#' 
#' Sauerbrei, W. and Royston, P., 1999. \emph{Building multivariable prognostic 
#' and diagnostic models: transformation of the predictors by using fractional 
#' polynomials. J Roy Stat Soc a Sta, 162:71-94.}\cr
#' 
#' Sauerbrei, W., Kipruto, E. and Balmford, J., 2023. \emph{Effects of influential 
#' points and sample size on the selection and replicability of multivariable 
#' fractional polynomial models. Diagnostic and Prognostic Research, 7(1), p.7.}\cr
#' 
#' Becher, H., Lorenz, E., Royston, P. and Sauerbrei, W., 2012. \emph{Analysing 
#' covariates with spike at zero: a modified FP procedure and conceptual issues.
#' Biometrical journal, 54(5), pp.686-700.}\cr
#' 
#' Lorenz, E., Jenkner, C., Sauerbrei, W. and Becher, H., 2019. \emph{Modeling exposures 
#' with a spike at zero: simulation study and practical application to survival data. 
#' Biostatistics & Epidemiology, 3(1), pp.23-37.}
#' 
#' @seealso 
#' \code{summary.mfp2()}, \code{coef.mfp2()}, \code{predict.mfp2()}, \code{fp()}
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
                         family = "gaussian",
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
                         zero_vars = NULL,
                         catzero_vars = NULL,
                         spike_vars = NULL,
                         min_saz_component_prop = 0.10,
                         force_max_fp = FALSE,
                         verbose = TRUE,
                         ...) {
  
  # this function prepares everything for fitting the actual mfp2 model
  
  cl <- match.call()
  
  # Public interface:
  # zero_vars, catzero_vars and spike_vars are character vectors of variable names.
  #
  # Internal representation:
  # The existing internal names zero, catzero and spike are retained because
  # they are converted below to named logical vectors over the columns of x.
  zero <- zero_vars
  catzero <- catzero_vars
  spike <- spike_vars
  
  # match arguments ------------------------------------------------------------
  criterion <- match.arg(criterion)
  xorder    <- match.arg(xorder)
  ties      <- match.arg(ties)
  
  # assertions -----------------------------------------------------------------
  # ----Family
  family_info <- normalize_family_argument(
    family,
    family_arg = deparse(substitute(family))
  )
  
  family        <- family_info$family
  family_string <- family_info$family_string
  # assert that x is a matrix
  if (!is.matrix(x)) {
    stop("! x must be a matrix", call. = FALSE)
  }
  
  # assert that x must not contain character values
  if (any(is.character(x))) {
    stop("! x contains characters values.\n",
         "i Please convert categorical variables to dummy variables.", 
         call. = FALSE)
  }
  
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
  if (is.null(vnames)) {
    stop("! The column names of x must not be missing.\n",
         "i Please set column names for x.",
         call. = FALSE)
  }
  
  if (!is.numeric(x)) {
    stop(
      "! `x` must be a numeric matrix.",
      sprintf("i Current storage mode is: %s.", typeof(x)),
      call. = FALSE
    )
  }
  
  # assert that x has no missing data
  if (anyNA(x)) {
    stop("! x must not contain any NA (missing data).\n",   
         "i Please remove any missing data before passing x to this function.",
         call. = FALSE)
  }
  
  if (any(!is.finite(x))) {
    stop(
      "! `x` must contain only finite, non-missing numeric values.",
      call. = FALSE
    )
  }
  
  # assert that subset must be a vector and does not contain negative values
if (!is.null(subset)) {
  if (is.logical(subset)) {
    if (length(subset) != nobs || anyNA(subset)) {
      stop(
        "! Logical subset must have length equal to the number of observations in x and contain no NA.",
        call. = FALSE
      )
    }
    subset <- which(subset)
  } else if (is.numeric(subset)) {
    if (
      anyNA(subset) ||
      any(!is.finite(subset)) ||
      any(subset != as.integer(subset)) ||
      any(subset < 1L) ||
      any(subset > nobs)
    ) {
      stop(
        "! Numeric subset must contain valid positive row indices within the range of x.",
        call. = FALSE
      )
    }
    subset <- as.integer(subset)
  } else {
    stop(
      "! subset must be either a logical vector or a numeric/integer vector of row indices.",
      call. = FALSE
    )
  }

  if (length(subset) < 5L) {
    stop(
      "! The selected subset is too small (<5) to fit an mfp model.",
      sprintf("i The number of selected observations is %d.", length(subset)),
      call. = FALSE
    )
  }
} 
  
  # Validate weights -----------------------------------------------------------
  if (!is.null(weights)) {
    if (!is.numeric(weights)) {
      stop(
        "! `weights` must be numeric.",
        sprintf("i Current type is: %s.", typeof(weights)),
        call. = FALSE
      )
    }
    
    if (length(weights) != nobs) {
      stop(
        "! The number of observations in x and weights must match.",
        sprintf(
          "i The number of rows in x is %d, but the number of elements in weights is %d.",
          nobs, length(weights)
        ),
        call. = FALSE
      )
    }
    
    if (anyNA(weights) || any(!is.finite(weights))) {
      stop(
        "! `weights` must contain only finite, non-missing values.",
        call. = FALSE
      )
    }
    
    if (any(weights < 0)) {
      stop(
        "! `weights` must not be negative.",
        call. = FALSE
      )
    }
  }
  
  # Validate offset ------------------------------------------------------------
  if (!is.null(offset)) {
    if (!is.numeric(offset)) {
      stop(
        "! `offset` must be numeric.",
        sprintf("i Current type is: %s.", typeof(offset)),
        call. = FALSE
      )
    }
    
    if (length(offset) != nobs) {
      stop(
        "! The number of observations in x and offset must match.",
        sprintf(
          "i The number of rows in x is %d, but the number of elements in offset is %d.",
          nobs, length(offset)
        ),
        call. = FALSE
      )
    }
    
    if (anyNA(offset) || any(!is.finite(offset))) {
      stop(
        "! `offset` must contain only finite, non-missing values.",
        call. = FALSE
      )
    }
  }
  
  # Validate scalar and vector options ----------------------------------------
  
  # cycles controls the maximum number of MFP backfitting cycles.
  # It must be a positive integer-like scalar.
  validate_positive_integer_scalar(cycles, "cycles")
  cycles <- as.integer(cycles)
  
  # verbose, and ftest are scalar logical flags.
  validate_logical_vector(verbose, "verbose", allowed_lengths = 1L)
  validate_logical_vector(ftest, "ftest", allowed_lengths = 1L)

  # alpha and select are probabilities, either scalar or one value per predictor.
  validate_probability_vector(alpha, "alpha", nvars)
  validate_probability_vector(select, "select", nvars)
  
  # center may be scalar or one logical value per predictor.
  validate_logical_vector(center, "center", allowed_lengths = c(1L, nvars))
  
  # force_max_fp may be scalar or one logical value per predictor.
  validate_logical_vector(force_max_fp, "force_max_fp", allowed_lengths = c(1L, nvars))
  
  # shift may be NULL, scalar, or length nvars.
  # NA is intentionally allowed because later code interprets NA as automatic
  # shift estimation for that variable.
  validate_numeric_vector(
    arg = shift,
    name = "shift",
    nvars = nvars,
    allow_null = TRUE,
    allow_na = TRUE,
    strictly_positive = FALSE
  )
  
  # scale may be NULL, scalar, or length nvars.
  # NA is intentionally allowed because later code interprets NA as automatic
  # scale estimation for that variable. Non-missing scale values must be > 0.
  validate_numeric_vector(
    arg = scale,
    name = "scale",
    nvars = nvars,
    allow_null = TRUE,
    allow_na = TRUE,
    strictly_positive = TRUE
  )
  
  # Validate variable-name arguments ------------------------------------------
  # These arguments define the model specification. Unknown names are treated as
  # errors because silently dropping them can hide spelling mistakes, e.g.
  # zero_vars = "expsoure" instead of "exposure".
  
  validate_variable_names(keep, "keep", vnames)
  validate_variable_names(zero, "zero_vars", vnames)
  validate_variable_names(catzero, "catzero_vars", vnames)
  validate_variable_names(acdx, "acdx", vnames)
  validate_variable_names(spike, "spike_vars", vnames)
  
  # Validate df ---------------------------------------------------------------
  if (!is.numeric(df)) {
    stop(
      "! `df` must be numeric.",
      sprintf("i Current type is: %s.", typeof(df)),
      call. = FALSE
    )
  }
  
  if (!length(df) %in% c(1L, nvars)) {
    stop(
      sprintf(
        "! `df` must be a single number or a numeric vector of length %d; got length %d.",
        nvars, length(df)
      ),
      call. = FALSE
    )
  }
  
  if (anyNA(df) || any(!is.finite(df))) {
    stop(
      "! `df` must contain only finite, non-missing values.",
      call. = FALSE
    )
  }
  
  if (any(df != as.integer(df))) {
    stop(
      "! `df` must contain integer values.",
      "i Valid values are 1 for linear terms, or even positive integers 2, 4, 6, ... for fractional polynomials.",
      call. = FALSE
    )
  }
  
  if (any(df <= 0)) {
    stop(
      "! `df` must contain only positive values.",
      "i Valid values are 1 for linear terms, or even positive integers 2, 4, 6, ... for fractional polynomials.",
      call. = FALSE
    )
  }
  
  if (any(df != 1 & df %% 2 != 0)) {
    stop(
      "! Any `df` value greater than 1 must be even.",
      "i Valid values are 1 for linear terms, or even positive integers 2, 4, 6, ... for fractional polynomials.",
      call. = FALSE
    )
  }
  
  # Revert to Chi-square when family is not Gaussian
  if (ftest && family_string != "gaussian") {
    warning(
      sprintf("i F-test not suitable for family = %s.\n", family_string),
      "i mfp2() reverts to use Chi-square instead.",
      call. = FALSE
    )
    ftest <- FALSE
  }
  
  # Validate response y ----------------------------------------------------------
  # Keep all family-specific response validation centralized in family.R.
  # This avoids drift between mfp2.default(), mfpi.default(), and future methods.
  validate_family_response(
    y = y,
    family_string = family_string,
    nobs = nobs
  )
  
  # Validate Cox-specific auxiliary inputs --------------------------------------
  # validate_family_response() checks the Cox response itself, but strata is an
  # auxiliary argument and therefore still needs to be checked here.
  if (family_string == "cox" && !is.null(strata)) {
    if (is.factor(strata)) {
      strata <- as.numeric(strata)
    }
    
    strata_len <- if (is.vector(strata)) length(strata) else NROW(strata)
    
    if (strata_len != nobs) {
      stop(
        "! The length of stratification factor(s) and the number of observations in x must match.",
        sprintf(
          "i The length of strata is %d, but the number of observations in x is %d.",
          strata_len, nobs
        ),
        call. = FALSE
      )
    }
  }
  
  # Validate spike-at-zero component proportion -------------------------------
  if (!is.numeric(min_saz_component_prop) ||
      length(min_saz_component_prop) != 1L ||
      anyNA(min_saz_component_prop) ||
      !is.finite(min_saz_component_prop) ||
      min_saz_component_prop <= 0 ||
      min_saz_component_prop >= 0.5) {
    stop(
      "! `min_saz_component_prop` must be a single finite numeric value in the open interval (0, 0.5).",
      call. = FALSE
    )
  }
  
  # Validate and build powers list ----------------------------------------------
  # validate_fp_power_list() should no longer make df-dependent decisions here. 
  # At this point, df may still later change because of:
  # ACD forcing df = 4 
  # SAZ positive-component df capping
  # low-cardinality df resets
  power_list <- validate_fp_power_list(
    powers = powers,
    vnames = vnames,
    arg_name = "powers"
  )
  
  # Default weights and offset
  if (is.null(weights)) {
    weights <- rep.int(1, nobs)
  }
  
  has_offset <- !is.null(offset)
  
  if (is.null(offset)) {
    offset <- rep.int(0, nobs)
  }
  
  # Default select and alpha
  if (length(select) == 1) {
    select <- rep(select, nvars)
  } 
  
  if (length(alpha) == 1) {
    alpha <- rep(alpha, nvars)
  } 
  
  # Default shift values
  if (is.null(shift)) {
    shift <- rep(NA_real_, nvars)
  } else {
    if (length(shift) == 1) {
      shift <- rep(shift, nvars)
    }
    
    if (length(shift) != nvars) {
      stop(
        "! shift must either be NULL, a single number, or the number of variables ",
        "(columns in x) and shift must match.",
        call. = FALSE
      )
    }
    
    shift <- setNames(shift, vnames)

  }
  
  # Default center
  if (length(center) == 1) {
    center <- rep(center, nvars)    
  }
  
  # Default force_max_fp
  if (length(force_max_fp) == 1L) {
    force_max_fp <- rep(force_max_fp, nvars)
  }
  force_max_fp <- setNames(force_max_fp, vnames)
  
  if (criterion == "pvalue" && any(force_max_fp)) {
    warning(
      "i `force_max_fp` has no effect when criterion = 'pvalue'.",
      "i Use select = 1 and alpha = 1 to force the most complex FP under p-value selection.",
      call. = FALSE
    )
  }
  
  # Set default control parameters for model fitting if not provided
  if (is.null(control)) {
    if (family_string == "cox") {
      control <- survival::coxph.control()
    } else {
      control <- stats::glm.control()
    }
  }
  
  # Convert zero_vars and catzero_vars to logical vectors
  if (is.null(zero)) {
    zero <- setNames(rep(FALSE, nvars), vnames)
  } else {
    zero_input_vars <- zero
    zero <- setNames(rep(FALSE, nvars), vnames)
    zero[vnames %in% zero_input_vars] <- TRUE
  }
  
  if (is.null(catzero)) {
    catzero <- setNames(rep(FALSE, nvars), vnames)
  } else {
    catzero_input_vars <- catzero
    catzero <- setNames(rep(FALSE, nvars), vnames)
    catzero[vnames %in% catzero_input_vars] <- TRUE
  }
  
  # Only check variables where zero is TRUE
  if (any(zero)) {
    vars_to_check <- vnames[zero]
    
    bad_vars <- vars_to_check[
      apply(x[, vars_to_check, drop = FALSE], 2, function(col) {
        all(col > 0, na.rm = TRUE)
      })
    ]
    
    if (length(bad_vars) > 0) {
      warning(
        "The following variables were marked through 'zero_vars' but contain only positive values. ",
        "Setting 'zero' and 'catzero' to FALSE for: ",
        paste(bad_vars, collapse = ", ")
      )
      
      zero[bad_vars] <- FALSE
      catzero[bad_vars] <- FALSE
    }
  }
  
  if (any(catzero)) {
    vars_to_check_cz <- vnames[catzero]
    bad_vars_cz <- vars_to_check_cz[
      apply(x[, vars_to_check_cz, drop = FALSE], 2, function(col) {
        all(col > 0, na.rm = TRUE)
      })
    ]
    if (length(bad_vars_cz) > 0) {
      warning(
        "The following variables were marked through 'catzero_vars' but contain ",
        "only positive values. Setting 'catzero' to FALSE for: ",
        paste(bad_vars_cz, collapse = ", ")
      )
      catzero[bad_vars_cz] <- FALSE
    }
  }
  
  #-------- Deal with acdx -----------------------------------------------------
  # Convert acdx to a named logical vector
  if (is.null(acdx)) {
    acdx <- setNames(rep(FALSE, nvars), vnames)
  } else {
    acdx <- unique(acdx)

    acdx_vec <- setNames(rep(FALSE, nvars), vnames)
    acdx_vec[acdx] <- TRUE
    acdx <- acdx_vec
  }
  
  #-------- Deal with spike_vars ----------------------------------------------
  # Convert spike_vars to a named logical vector
  if (is.null(spike)) {
    spike <- setNames(rep(FALSE, nvars), vnames)
  } else {
    spike <- unique(spike)

    spike_input_vars <- spike
    spike <- setNames(rep(FALSE, nvars), vnames)
    spike[spike_input_vars] <- TRUE
  }
  
  # Identify binary columns, exactly two unique non-missing values.
  # Binary variables are not eligible for zero/catzero/spike handling here:
  #   - zero handling is unnecessary because the variable is already discrete;
  #   - catzero would create a redundant indicator;
  #   - spike implies catzero and would therefore create the same inconsistency.
  binary_vars <- apply(x, 2, function(col) length(unique(col[!is.na(col)])) == 2)
  binary_names <- vnames[binary_vars]
  
  # Reset zero, catzero, and spike for binary variables, with a single warning.
  # This must reset spike as well as zero/catzero. Otherwise a binary variable
  # requested through spike_vars can be left in an inconsistent transient state:
  #   spike = TRUE, catzero = FALSE, zero = FALSE
  # until fit_mfp() attempts to repair or further process it.
  if (length(binary_names) > 0) {
    zero_binary    <- intersect(binary_names, names(zero)[zero])
    catzero_binary <- intersect(binary_names, names(catzero)[catzero])
    spike_binary   <- intersect(binary_names, names(spike)[spike])
    
    vars_to_reset <- Reduce(union, list(zero_binary, catzero_binary, spike_binary))
    
    if (length(vars_to_reset) > 0) {
      warning(
        "The following binary variables were marked through 'zero_vars', ",
        "'catzero_vars', or 'spike_vars' but are binary and will be reset to FALSE: ",
        paste(vars_to_reset, collapse = ", "), ".",
        call. = FALSE
      )
      zero[vars_to_reset]    <- FALSE
      catzero[vars_to_reset] <- FALSE
      spike[vars_to_reset]   <- FALSE
    }
  }
  
  # Resolve spike-at-zero eligibility before df and shift/scale preprocessing ----
  # This must happen before assign_df() and before shift values are chosen.
  # A variable requested as spike-at-zero but found ineligible should revert to
  # the user's explicit zero/catzero choices before ordinary FP preprocessing.
  if (any(spike)) {
    saz_flags <- resolve_saz_eligibility(
      x                      = x,
      spike                  = spike,
      catzero                = catzero,
      zero                   = zero,
      min_saz_component_prop = min_saz_component_prop
    )
    
    spike   <- saz_flags$spike
    catzero <- saz_flags$catzero
    zero    <- saz_flags$zero
  }
  
  # set df ---------------------------------------------------------------------
  if (length(df) == 1) {
    if (df != 1) {
      df.list <- assign_df(x = x, df_default = df)
    } else {
      df.list <- rep(df, nvars)
    }
  } else {
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
  
  # Compute effective zero/catzero by applying the spike -> catzero -> zero
  # cascade. This is needed to correctly set shift = 0 for spike variables,
  # whose structural zeros must remain at zero after shifting. The original
  # (un-cascaded) zero, catzero, and spike vectors are passed unchanged to
  # fit_mfp(), which saves user intent before re-applying the cascade internally.
  effective_catzero <- catzero
  effective_zero <- zero
  effective_catzero[spike] <- TRUE        # spike implies catzero
  effective_zero[effective_catzero] <- TRUE # catzero implies zero
  
  # Reset shift = 0 for variables with effective zero or catzero status
  # (including those implied by spike), since only the positive part is
  # transformed and the zero group must remain at zero.
  # Similar logic applies to variables with df = 1.
  shift_to_zero <- rep(FALSE, length(vnames))
  names(shift_to_zero) <- vnames
  
  shift_to_zero[names(effective_zero)[effective_zero]] <- TRUE
  shift_to_zero[names(effective_catzero)[effective_catzero]] <- TRUE
  shift_to_zero[vnames[df.list == 1]] <- TRUE
  
  shift[shift_to_zero] <- 0
  
  shift_missing <- is.na(shift)
  if (any(shift_missing)) {
    shift[shift_missing] <- apply(
      x[, shift_missing, drop = FALSE],
      2,
      find_shift_factor
    )
  }
  
  # data preparation -----------------------------------------------------------
  # Apply shift transformation to the x matrix
  x <- sweep(x, 2, shift, "+")
  
  # scaling should be done after shifting
  if (is.null(scale)) {
    scale <- apply(x, 2, find_scale_factor)
  } else {
    if (length(scale) == 1) {
      scale <- rep(scale, nvars) 
    }
    
    if (length(scale) != nvars) {
      stop(
        "! scale must either be NULL, a single number, or the number of variables ",
        "(columns in x) and scale must match.",
        call. = FALSE
      )
    }
    
    # NA values mean automatic scaling for those variables.
    # This allows the formula interface to mix global/default automatic scaling
    # with per-variable scale values from fp(..., scale = ...).
    scale_missing <- is.na(scale)
    
    if (any(!scale_missing & scale <= 0)) {
      stop("All non-missing scale values must be positive.", call. = FALSE)
    }
    
    if (any(scale_missing)) {
      scale[scale_missing] <- apply(
        x[, scale_missing, drop = FALSE],
        2,
        find_scale_factor
      )
    }
  }
  
  # Identify variables that require fractional polynomial transformation
  nonlinear_variables <- which(df.list != 1)
  nonlinear_names <- vnames[nonlinear_variables]
  
  # Exclude variables that have effective zero or catzero status (including
  # spike-implied) from the positivity check, as their zeros are intentional
  all_zero_vars <- union(names(effective_zero)[effective_zero],
                         names(effective_catzero)[effective_catzero])
  vars_to_check <- setdiff(nonlinear_names, all_zero_vars)
  
  if (length(vars_to_check) > 0) {
    check_indices <- match(vars_to_check, vnames)
    xd <- x[, check_indices, drop = FALSE]
    
    neg_cols <- which(colSums(xd <= 0, na.rm = TRUE) > 0)
    
    if (length(neg_cols) > 0) {
      cols_with_negatives <- colnames(xd)[neg_cols]
      
      stop(
        "i The shifting factors are insufficient to ensure positive values for fractional polynomial transformation.\n", 
        sprintf("i Problematic variables: %s", 
                paste0(cols_with_negatives, collapse = ", ")), 
        "\ni Consider increasing the shift values for these variables.",
        call. = FALSE
      )
    }
  }
  
  # scale x after shifting
  x <- sweep(x, 2, scale, "/")
  
  # stratification for cox model
  istrata <- strata
  if (family_string == "cox" && !is.null(strata)) {
    istrata <- survival::strata(strata, shortlabel = TRUE)
    istrata <- as.integer(istrata)
  }
  
  # data subsetting ------------------------------------------------------------
  if (!is.null(subset)) {

    x <- x[subset, , drop = FALSE]
    
    if (is.matrix(y)) {
      y <- y[subset, , drop = FALSE]
    } else {
      y <- y[subset]
    }
    
    weights <- weights[subset]
    offset <- offset[subset]
    istrata <- istrata[subset]
  }
  
  # Fail early if the default-interface design matrix is not estimable.
  validate_default_design_rank(
    x = x,
    intercept = family_string != "cox"
  )
  
  # fit model ------------------------------------------------------------------
  fit <- fit_mfp(
    x = x, y = y, 
    weights = weights, offset = offset, cycles = cycles, 
    scale = scale, shift = shift, df = df.list, center = center, 
    family = family, family_string = family_string, criterion = criterion, 
    select = select, alpha = alpha, keep = keep, xorder = xorder, 
    powers = power_list, method = ties, strata = istrata, nocenter = nocenter, 
    acdx = acdx, ftest = ftest, force_max_fp = force_max_fp,
    control = control, zero = zero, catzero = catzero, spike = spike,
    min_saz_component_prop = min_saz_component_prop, 
    saz_pre_resolved = TRUE,
    has_offset = has_offset,
    verbose = verbose
  )
  
  # add additional information to fitted object
  fit$call_mfp <- cl
  fit$family <- family
  fit$family_string <- family_string
  fit$offset <- offset
  fit$has_offset <- has_offset
  
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
                         family = "gaussian",
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
                         min_saz_component_prop = 0.10,
                         verbose = TRUE,
                         ...) {
  # capture the call
  call <- match.call()
  
  # match arguments ------------------------------------------------------------
  criterion <- match.arg(criterion)
  xorder <- match.arg(xorder)
  ties <- match.arg(ties)
  
  dots <- list(...)
  
  if ("acd_vars" %in% names(dots)) {
    stop(
      "`acd_vars` is not supported as an argument to `mfp2.formula()`. ",
      "Use `fp(..., acdx = TRUE)` or `fp2(..., acdx = TRUE)` inside the formula ",
      "to request ACD handling for formula terms.",
      call. = FALSE
    )
  }
  
  # Allowed families

  family_string <- get_family_string_formula(
    family,
    family_arg = deparse(substitute(family))
  )
  
  
  if (!is.null(strata) && family_string != "cox") {
    stop(
      "! `strata` is only allowed for Cox models.\n",
      "i Please use `family = \"cox\"` or remove the `strata` argument.",
      call. = FALSE
    )
  }
  
  # assert that data must be provided
  if (missing(data)) {
    stop("! data argument is missing.\n",
         "i An input data.frame is required for the use of mfp2.",
         call. = FALSE)
  }
  
  # assert that data has column names
  if (is.null(colnames(data))) {
    stop("! data must have column names.\n",
         "i Please set column names.", call. = FALSE)
  }
  
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }
  
  # assert that a formula must be provided
  if (missing(formula)) {
    stop("! formula is missing.", call. = FALSE)
  }
  
  if (!inherits(formula, "formula")) {
    stop("method is only for formula objects", call. = FALSE)
  }
  
  # Validate scalar formula-interface defaults ---------------------------------
  # In the formula interface, df, alpha, select, shift, scale, and center are
  # global scalar defaults. Per-variable values must be supplied inside fp().
  
  validate_formula_scalar <- function(arg,
                                      arg_name,
                                      type,
                                      allow_null = FALSE,
                                      allow_na = FALSE,
                                      hint = NULL) {
    if (is.null(arg)) {
      if (allow_null) return(invisible(TRUE))
      
      msg <- sprintf("! `%s` must not be NULL.", arg_name)
      if (!is.null(hint)) msg <- paste0(msg, "\n", hint)
      
      stop(msg, call. = FALSE)
    }
    
    ok_type <- switch(
      type,
      numeric = is.numeric(arg),
      logical = is.logical(arg),
      stop("Internal error: unsupported validator type.", call. = FALSE)
    )
    
    if (!ok_type || length(arg) != 1L) {
      expected <- switch(
        type,
        logical = "a single TRUE or FALSE value",
        numeric = "a single numeric value",
        type
      )
      
      actual <- paste0(arg, collapse = ", ")
      
      msg <- sprintf(
        "! `%s` must be %s.\ni You supplied: %s.",
        arg_name, expected, actual
      )
      
      if (!is.null(hint)) msg <- paste0(msg, "\n", hint)
      
      stop(msg, call. = FALSE)
    }
    
    if (!allow_na && anyNA(arg)) {
      expected <- switch(
        type,
        logical = "TRUE or FALSE",
        numeric = "a numeric value",
        type
      )
      
      msg <- sprintf(
        "! `%s` must not be NA.\ni Please use %s.",
        arg_name, expected
      )
      
      if (!is.null(hint)) msg <- paste0(msg, "\n", hint)
      
      stop(msg, call. = FALSE)
    }
    
    if (type == "numeric" && !anyNA(arg) && !is.finite(arg)) {
      msg <- sprintf("! `%s` must be finite.", arg_name)
      if (!is.null(hint)) msg <- paste0(msg, "\n", hint)
      
      stop(msg, call. = FALSE)
    }
    
    invisible(TRUE)
  }
  
  validate_formula_probability <- function(arg, arg_name, hint = NULL) {
    validate_formula_scalar(
      arg = arg,
      arg_name = arg_name,
      type = "numeric",
      allow_null = FALSE,
      allow_na = FALSE,
      hint = hint
    )
    
    if (arg < 0 || arg > 1) {
      msg <- sprintf(
        "! `%s` must be between 0 and 1.\ni You supplied: %s.",
        arg_name, arg
      )
      
      if (!is.null(hint)) msg <- paste0(msg, "\n", hint)
      
      stop(msg, call. = FALSE)
    }
    
    invisible(TRUE)
  }
  
  validate_formula_scalar(
    df,
    "df",
    "numeric",
    allow_null = FALSE,
    allow_na = FALSE,
    hint = "i Use fp(..., df = value) to set different df values for individual variables."
  )
  
  if (df != as.integer(df)) {
    stop(
      sprintf(
        "! `df` must be an integer value.\ni You supplied: %s.\ni Use fp(..., df = value) to set different df values for individual variables.",
        df
      ),
      call. = FALSE
    )
  }
  
  if (df <= 0) {
    stop(
      sprintf(
        "! `df` must be positive.\ni You supplied: %s.\ni Use fp(..., df = value) to set different df values for individual variables.",
        df
      ),
      call. = FALSE
    )
  }
  
  if (df != 1 && df %% 2 != 0) {
    stop(
      sprintf(
        "! `df` must be 1 or an even positive integer.\ni You supplied: %s.\ni Use fp(..., df = value) to set different df values for individual variables.",
        df
      ),
      call. = FALSE
    )
  }
  
  validate_formula_probability(
    alpha,
    "alpha",
    hint = "i Use fp(..., alpha = value) to set different alpha values for individual variables."
  )
  
  validate_formula_probability(
    select,
    "select",
    hint = "i Use fp(..., select = value) to set different select values for individual variables."
  )
  
  validate_formula_scalar(
    scale,
    "scale",
    "numeric",
    allow_null = TRUE,
    allow_na = TRUE,
    hint = "i Use fp(..., scale = value) to set different scaling factors for individual variables."
  )
  
  if (!is.null(scale) && !is.na(scale) && scale <= 0) {
    stop(
      sprintf(
        "! `scale` must be positive, NULL, or NA.\ni You supplied: %s.\ni Use fp(..., scale = value) to set different scaling factors for individual variables.",
        scale
      ),
      call. = FALSE
    )
  }
  
  validate_formula_scalar(
    shift,
    "shift",
    "numeric",
    allow_null = TRUE,
    allow_na = TRUE,
    hint = "i Use fp(..., shift = value) to set different shift values for individual variables."
  )
  
  validate_formula_scalar(
    center,
    "center",
    "logical",
    allow_null = FALSE,
    allow_na = FALSE,
    hint = "i Use fp(..., center = TRUE) or fp(..., center = FALSE) to set different center values for individual variables."
  )

  n_data <- nrow(data)
  
  if (!is.null(weights)) {
    if (!is.numeric(weights)) {
      stop("! `weights` must be numeric.", call. = FALSE)
    }
    
    if (length(weights) != n_data) {
      stop(
        sprintf(
          "! `weights` must have one value per row of `data`.\ni `weights` has length %d, but `data` has %d rows.",
          length(weights), n_data
        ),
        call. = FALSE
      )
    }
    
    if (anyNA(weights) || any(!is.finite(weights))) {
      stop(
        "! `weights` must contain only finite, non-missing values.",
        call. = FALSE
      )
    }
    
    if (any(weights < 0)) {
      stop("! `weights` must not contain negative values.", call. = FALSE)
    }
  }
  
  if (!is.null(offset)) {
    if (!is.numeric(offset)) {
      stop("! `offset` must be numeric.", call. = FALSE)
    }
    
    if (length(offset) != n_data) {
      stop(
        sprintf(
          "! `offset` must have one value per row of `data`.\ni `offset` has length %d, but `data` has %d rows.",
          length(offset), n_data
        ),
        call. = FALSE
      )
    }
    
    if (anyNA(offset) || any(!is.finite(offset))) {
      stop(
        "! `offset` must contain only finite, non-missing values.",
        call. = FALSE
      )
    }
  }
  
  if (!is.null(subset)) {
    if (is.logical(subset)) {
      if (length(subset) != n_data || anyNA(subset)) {
        stop(
          "! Logical `subset` must have one TRUE/FALSE value per row of `data` and contain no NA.",
          call. = FALSE
        )
      }
    } else if (is.numeric(subset)) {
      if (
        anyNA(subset) ||
        any(!is.finite(subset)) ||
        any(subset != as.integer(subset)) ||
        any(subset < 1L) ||
        any(subset > n_data)
      ) {
        stop(
          "! Numeric `subset` must contain valid positive row indices within `data`.",
          call. = FALSE
        )
      }
    } else {
      stop(
        "! `subset` must be either a logical vector or a numeric/integer vector of row indices.",
        call. = FALSE
      )
    }
  }
  
  terms_formula <- stats::terms(
    formula,
    specials = "strata",
    data = data
  )
  
  has_formula_strata <- !is.null(attr(terms_formula, "specials")$strata)
  
  if (!is.null(strata) && !has_formula_strata) {
    strata_n <- if (is.vector(strata) || is.factor(strata)) {
      length(strata)
    } else {
      NROW(strata)
    }
    
    if (strata_n != n_data) {
      stop(
        sprintf(
          "! `strata` must have one value or row per row of `data`.\ni `strata` has %d rows/values, but `data` has %d rows.",
          strata_n, n_data
        ),
        call. = FALSE
      )
    }
    
    if (anyNA(strata)) {
      stop("! `strata` must not contain missing values.", call. = FALSE)
    }
  }
  
  # model.frame preserves the attributes of the data unlike model.matrix.
  # Use na.pass here so that mfp2.formula() does not silently drop rows through
  # the global/default na.action option. Missingness is checked explicitly below.
  
  mf <- stats::model.frame(
    terms_formula,
    data = data,
    drop.unused.levels = TRUE,
    na.action = stats::na.pass
  )
  
  
  # Future work:
  # Support ordered factors as grouped formula terms rather than rejecting them.
  # `mfp2.formula()` should preserve the `assign` metadata returned by
  # model.matrix(), mapping each original formula term to all generated design
  # columns. `mfp2.default()` and the MFP step functions should then operate on
  # term blocks, not only on individual matrix columns. For example, an ordered
  # factor `x4` expanded to `x4.L` and `x4.Q` should be tested, selected, kept,
  # or dropped as one block. This will require extending the internal variable
  # metadata from a flat variable-name vector to a structure such as:
  #
  #   term_blocks <- list(
  #     x1 = list(type = "fp", columns = "x1"),
  #     x2 = list(type = "linear", columns = "x2"),
  #     x4 = list(type = "factor_block", columns = c("x4.L", "x4.Q"))
  #   )
  #
  # The default matrix interface can keep its current one-column-per-variable
  # behavior by constructing trivial one-column blocks internally.
  
  # check whether no predictor exists in the model, i.e. y ~ 1
  labels <- attr(terms_formula, "term.labels")
  if (length(labels) == 0) {
    stop(
      "! No predictors are provided for model fitting.\n",
      "i At least one predictor is required.",
      call. = FALSE
    )
  }

  if (anyNA(mf)) {
    na_rows <- which(!stats::complete.cases(mf))
    
    msg <- sprintf(
      "! Missing values are not allowed in variables used by the model formula.\ni Rows with missing values: %s.",
      paste(utils::head(na_rows, 10L), collapse = ", ")
    )
    
    if (length(na_rows) > 10L) {
      msg <- paste0(msg, sprintf("\ni Showing first 10 of %d affected rows.", length(na_rows)))
    }
    
    msg <- paste0(msg, "\ni Please remove or impute missing values before calling mfp2().")
    
    stop(msg, call. = FALSE)
  }
  
 
  # stratification and offset handling -----------------------------------------

  # Terms to remove before constructing the model matrix.
  # offset() and strata() are handled separately and must not enter x.
  terms_drop <- integer(0L)
  
  # stratification for Cox models ---------------------------------------------
  if (!is.null(attr(terms_formula, "specials")$strata)) {
    if (family_string == "cox") {
      
      # check whether strata is both in the formula and in the argument
      if (!is.null(call$strata)) {
        warning("i strata appear both in the formula and as an input argument.\n",
                "i The information in the formula is used and the input argument ignored.",
                call. = FALSE)
      }
      
      # untangle the terms for strata as in coxph
      stemp <- survival::untangle.specials(
        terms_formula,
        special = "strata",
        order = 1
      )
      
      missing_strata_vars <- setdiff(stemp$vars, names(mf))
      if (length(missing_strata_vars) > 0L) {
        stop(
          sprintf(
            "! Could not extract strata variable(s) from the model frame: %s.",
            paste(missing_strata_vars, collapse = ", ")
          ),
          call. = FALSE
        )
      }
      
      if (length(stemp$vars) == 1) {
        strata <- mf[[stemp$vars]]
      } else {
        strata <- mf[, stemp$vars]
      }
      
      terms_drop <- c(terms_drop, stemp$terms)
    } else {
      stop("! strata are only allowed for Cox models.\n", 
           "i Please remove any strata terms from the model formula.",
           call. = FALSE)
    }
  }
  
  # offset ---------------------------------------------------------------------
  term_offset <- attr(terms_formula, "offset")
  
  if (!is.null(term_offset) && length(term_offset) > 1) {
    stop("! Only one offset in the formula is allowed.", call. = FALSE)
  }
  
  if (!is.null(term_offset)) {
    if (!is.null(call$offset)) {
      warning(
        "i Offset appears both in the formula and as an input argument.\n",
        "i The information in the model formula is used and the input argument is ignored.",
        call. = FALSE
      )
    }
    
    offset <- as.vector(stats::model.offset(mf))
    terms_drop <- c(terms_drop, term_offset)
  }
  
  if (!is.null(offset)) {
    if (!is.numeric(offset)) {
      stop("! `offset` must be numeric.", call. = FALSE)
    }
    
    if (length(offset) != nrow(mf)) {
      stop("! `offset` must have one value per observation.", call. = FALSE)
    }
    
    if (anyNA(offset) || any(!is.finite(offset))) {
      stop(
        "! `offset` must contain only finite, non-missing values.",
        call. = FALSE
      )
    }
  }
  
  # Drop strata() and offset() terms before using model.matrix().
  if (length(terms_drop) > 0L) {
    terms_model <- terms_formula[-unique(terms_drop)]
  } else {
    terms_model <- terms_formula
  }
  
  remaining_labels <- attr(terms_model, "term.labels")
  
  if (length(remaining_labels) == 0L) {
    stop(
      "! No predictors are provided for model fitting after removing strata() and offset() terms.\n",
      "i At least one non-strata, non-offset predictor is required.",
      call. = FALSE
    )
  }
  
  # data preparation -----------------------------------------------------------
  y <- model.extract(mf, "response")

  mm <- stats::model.matrix(terms_model, mf)
  # Ordered factors must be rejected before model.matrix().
  #
  # model.matrix() expands ordered factors into polynomial contrast columns
  # such as .L, .Q, and .C. The MFP engine would then treat those contrast
  # columns as separate predictors, which breaks variable-wise selection and
  # degree testing.
  ordered_check_vars <- all.vars(
    stats::delete.response(terms_model)
  )
  
  check_ordered_factor_variables(
    vars = ordered_check_vars,
    data = data,
    interface = "mfp2.formula"
  )
  
  assign <- attr(mm, "assign")
  term_labels <- attr(terms_model, "term.labels")

  # Remove intercept, but keep assign aligned with the remaining columns.
  keep_cols <- colnames(mm) != "(Intercept)"
  x <- mm[, keep_cols, drop = FALSE]
  assign <- assign[keep_cols]
  
  nx <- ncol(x) 
  names_x <- colnames(x)
  
  
  # Map each original formula term to the model-matrix columns it generated.
  # This is used to expand `keep` safely for factor terms without prefix matching.
  term_to_columns <- split(colnames(x), term_labels[assign])
  
  # Remove any empty/NA assignments defensively.
  term_to_columns <- term_to_columns[!is.na(names(term_to_columns))]
  
  # select variables that undergo fp transformation and extract their attributes
  fp_pos <- which(is_fp_term(colnames(mf)))
  
  if (length(fp_pos) > 0) {
    fp_data <- mf[, fp_pos, drop = FALSE]
    
    # extract names of the variables that undergo fp transformation
    fp_vars <- unname(sapply(fp_data, function(v) attr(v, "name")))
    
    # check for variables used more than once in fp() function 
    fp_vars_duplicates <- fp_vars[duplicated(fp_vars)]
    if (length(fp_vars_duplicates) != 0) {
      stop("! Variables should be used only once in the fp() within the formula.\n", 
           sprintf("i The following variable(s) are duplicated in fp() function: %s.", 
                   paste0(fp_vars_duplicates, collapse = ", ")), 
           call. = FALSE)
    }
    
    # check for variables used in fp() as well as other parts of the formula
    vars_duplicates <- which(colnames(mf) %in% fp_vars)
    if (length(vars_duplicates) != 0) {
      stop("! Variables used in the fp() should not be included in other parts of the formula.\n", 
           sprintf("i This applies to the following variable(s): %s.", 
                   paste0(colnames(mf)[vars_duplicates], collapse = ", ")), 
           call. = FALSE)
    }
    
    # Capture fp() term labels before renaming them in the model matrix.
    fp_terms <- names_x[is_fp_term(names_x)]
    
    # replace names such as fp(x1) by real name "x1" in the x matrix
    #names_x <- replace(names_x, grep("^fp\\(.*\\)$", names_x), fp_vars)
    names_x <- replace(names_x, which(is_fp_term(names_x)), fp_vars)
    
    colnames(x) <- names_x
    
    # Update the formula-term -> model-matrix-column map after renaming fp() columns.
    # This prevents keep expansion from using unsafe prefix matching.
    # Both keep = "fp(x1)" and keep = "x1" resolve to the renamed model column "x1".
    if (length(fp_terms) > 0L) {
      for (i in seq_along(fp_terms)) {
        fp_term <- fp_terms[i]
        fp_var <- fp_vars[i]
        
        if (fp_term %in% names(term_to_columns)) {
          old_cols <- term_to_columns[[fp_term]]
          new_cols <- ifelse(old_cols %in% fp_terms, fp_var, old_cols)
          
          term_to_columns[[fp_term]] <- new_cols
          term_to_columns[[fp_var]] <- new_cols
        }
      }
    }
  }
  
  # call default method --------------------------------------------------------
  
  # if fp() is not used in the formula, it reduces to mfp2.default()
  df_list <- setNames(as.list(assign_df(x = x, df_default = df)), names_x)
  
  # scaling
  # NA means automatic scaling. mfp2.default() will estimate these values
  # after applying its full shift logic, including zero/catzero/spike/df resets.
  if (is.null(scale)) {
    scale_list <- setNames(rep(list(NA_real_), nx), names_x)
  } else {
    scale_list <- setNames(rep(list(scale), nx), names_x)
  }
  
  # shifting
  # NA means automatic shift. mfp2.default() will estimate the shift and then
  # apply its full zero/catzero/spike/df reset logic.
  if (is.null(shift)) {
    shift_list <- setNames(rep(list(NA_real_), nx), names_x)
  } else {
    shift_list <- setNames(rep(list(shift), nx), names_x)
  }
  
  # center, alpha, select
  center_list <- setNames(rep(list(center), nx), names_x)
  alpha_list <- setNames(rep(list(alpha), nx), names_x)
  select_list <- setNames(rep(list(select), nx), names_x)
  force_max_fp_list <- setNames(rep(list(FALSE), nx), names_x)
  
 
  # Powers supplied inside fp() are collected separately so they can override
  # top-level `powers` after both have been validated against the final model
  # matrix column names and effective df values.
  powerx <- list()
  
  # if fp() is used in the formula
  acdx <- NULL
  zero <- NULL
  catzero <- NULL
  spike <- NULL
  
  if (length(fp_pos) != 0) {
    df_list <- modifyList(df_list, 
                          setNames(lapply(fp_data, attr, "df"), fp_vars))
    
    scale_list <- modifyList(
      scale_list, 
      Filter(Negate(is.null), setNames(lapply(fp_data, attr, "scale"), fp_vars))
    )
    
    shift_list <- modifyList(
      shift_list, 
      Filter(Negate(is.null), setNames(lapply(fp_data, attr, "shift"), fp_vars))
    )
    
    center_list <- modifyList(center_list,
                              setNames(lapply(fp_data, attr, "center"), fp_vars))
    
    alpha_list <- modifyList(alpha_list, 
                             setNames(lapply(fp_data, attr, "alpha"), fp_vars))
    
    select_list <- modifyList(select_list, 
                              setNames(lapply(fp_data, attr, "select"), fp_vars))
    
    force_max_fp_list <- modifyList(
      force_max_fp_list,
      setNames(lapply(fp_data, attr, "force_max_fp"), fp_vars)
    )
    
    # Preference is given to powers supplied in fp() over powers argument
    powerx <- Filter(
      Negate(is.null),
      setNames(lapply(fp_data, attr, "powers"), fp_vars)
    )
    
    nax <- if (is.null(powers)) {
      character(0L)
    } else {
      intersect(names(powerx), names(powers))
    }
    
    if (length(nax) != 0) {
      warning(
        "i Powers are specified both in `fp()` and in the `powers` argument.",
        "\ni The `fp()`-specific powers take precedence; the corresponding `powers` argument entries are ignored.",
        sprintf(
          "\ni This applies to: %s.",
          paste(nax, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  
    # In the formula interface, zero, catzero and spike are logical attributes
    # attached to individual fp() terms. Before calling mfp2.default(), they are
    # converted to character vectors of variable names and passed as *_vars.
    acdx <- setNames(sapply(fp_data, attr, "acd"), fp_vars)
    zero <- setNames(sapply(fp_data, attr, "zero"), fp_vars)
    catzero <- setNames(sapply(fp_data, attr, "catzero"), fp_vars)
    spike <- setNames(sapply(fp_data, attr, "spike"), fp_vars)
    
    if (sum(acdx) == 0) {
      acdx <- NULL
    } else {
      acdx <- names(acdx[acdx])
    }
    
    if (sum(zero) == 0) {
      zero <- NULL
    } else {
      zero <- names(zero[zero])
    }
    
    if (sum(catzero) == 0) {
      catzero <- NULL
    } else {
      catzero <- names(catzero[catzero])
    }
    
    if (sum(spike) == 0) {
      spike <- NULL
    } else {
      spike <- names(spike[spike])
    }
  
  }
  
  # Validate and build candidate FP power list --------------------------------
  df_vec <- unlist(df_list, use.names = TRUE)
  
  powers_arg <- powers
  
  if (!is.null(powers_arg) && length(powerx) > 0L) {
    powers_arg[names(powerx)] <- NULL
    
    if (length(powers_arg) == 0L) {
      powers_arg <- NULL
    }
  }
  
  # formula-interface candidate-power parsing should not reject df-dependent 
  # cases before final df modifications happen.
  power_list <- validate_fp_power_list(
    powers = powers_arg,
    vnames = names_x,
    arg_name = "powers"
  )
  
  if (length(powerx) > 0L) {
    fp_power_list <- validate_fp_power_list(
      powers = powerx,
      vnames = names_x,
      arg_name = "fp() powers"
    )
    
    power_list[names(powerx)] <- fp_power_list[names(powerx)]
  }
  
  # --- handle factor variables when keep is used ------------------------------
  # Expand keep using exact formula-term and model-matrix-column matching.
  # This avoids unsafe prefix matching, e.g. keep = "age" matching "age_group".
  if (!is.null(keep)) {
    
    expanded_keep <- character(0)
    
    # Valid keep entries are either final model-matrix columns or original formula terms.
    valid_keep <- unique(c(colnames(x), names(term_to_columns)))
    bad_keep <- setdiff(keep, valid_keep)
    
    if (length(bad_keep) > 0L) {
      stop(
        sprintf(
          "! Unknown variable(s) in keep: %s.",
          paste(bad_keep, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    
    for (k in keep) {
      if (k %in% colnames(x)) {
        # Exact final column match, e.g. keep = "age"
        expanded_keep <- c(expanded_keep, k)
      } else if (k %in% names(term_to_columns)) {
        # Exact formula-term match, e.g. keep = "sex" expanding to factor columns
        expanded_keep <- c(expanded_keep, term_to_columns[[k]])
      }
    }
    
    keep <- unique(expanded_keep)
  }
  
  # Store the raw model-matrix column names before fp() columns are renamed.
  # These are needed by predict.mfp2() to rebuild the same expanded design
  # matrix from ordinary formula-style newdata.
  formula_column_map <- setNames(names_x, colnames(mm)[keep_cols])
  
  fit <- mfp2.default(x = x, 
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
               acdx = acdx,
               ftest = ftest,
               control = control,
               zero_vars = zero,
               catzero_vars = catzero,
               spike_vars = spike,
               force_max_fp = unlist(force_max_fp_list),
               min_saz_component_prop = min_saz_component_prop,
               verbose = verbose)
  fit$formula_interface <- TRUE
  fit$formula <- formula
  fit$formula_terms <- stats::delete.response(terms_model)
  fit$formula_contrasts <- attr(mm, "contrasts")
  fit$formula_xlevels <- .getXlevels(terms_model, mf)
  fit$formula_column_map <- formula_column_map
  fit$formula_model_matrix_columns <- names_x
  fit$formula_term_to_columns <- term_to_columns
  
  fit
}

#' Helper function to identify fp terms when fp() or fp2() is used in the formula
#'
#' @param z A character vector.
#'
#' @return A logical vector indicating whether each element matches 
#' an `fp()` or `fp2()` term.
#'
#' @keywords internal
#' @noRd
is_fp_term <- function(z) {
  grepl("^fp2?\\(.*\\)$", z)
}

#' Extract coefficients from object of class `mfp2`
#' 
#' This function is a method for the generic [stats::coef()] function for 
#' objects of class `mfp2`. 
#' 
#' @param object an object of class `mfp2`, usually, a result of a call to
#' \code{mfp2()}.
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
#' \code{mfp2()}.
#' @param ... further arguments passed to the summary functions for `glm()` 
#' ([stats::summary.glm()], i.e. families supported by `glm()`) or `coxph()` 
#' ([survival::summary.coxph()], if `object$family = "cox"`).
#' 
#' @return 
#' An object returned from [stats::summary.glm()] or
#' [survival::summary.coxph()], depending on the family parameter of `object`.
#' 
#' @seealso 
#' \code{mfp2()}, [stats::glm()], [stats::summary.glm()], [survival::coxph()],
#' [survival::summary.coxph()]
#' 
#' @export
summary.mfp2 <- function(object, ...) {
  if (!inherits(object, "mfp2")) {
    stop("The object is not an mfp2 object.", call. = FALSE)
  }
  
  if (identical(object$family_string, "cox")) {
    return(
      getFromNamespace("summary.coxph", "survival")(
        object,
        ...
      )
    )
  }
  
  stats::summary.glm(
    object,
    ...
  )
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
#' @return Two dataframes: the first one contains preprocessing parameters
#'   (shifting, scaling, and centering), and the second one includes additional
#'   parameters such as `df`, `select`, and `alpha` passed through `mfp2`. It
#'   also returns a list of the final model fitted, which can be either a GLM or
#'   Cox model depending on the chosen family.
#' @export
print.mfp2 <- function(x, ...) {
  # Shift and scaling factors with centering values
  cat("Shifting, Scaling and Centering of covariates", "\n")
  print.data.frame(x$transformations)
  cat("\n")
  
  # Final MFP powers
  cat("Final Multivariable Fractional Polynomial for y", "\n")
  
  # Prepare a print-only copy. Do not modify x$fp_terms itself because internal
  # code may rely on the original column names and numeric spike_dec coding.
  fp_terms <- x$fp_terms
  
  # Convert internal numeric spike_dec to a user-facing SAZ decision label,
  # preserving the original column position so power columns are not moved.
  if ("spike_dec" %in% names(fp_terms)) {
    selected <- if ("selected" %in% names(fp_terms)) {
      fp_terms[["selected"]]
    } else {
      rep(TRUE, nrow(fp_terms))
    }
    
    is_saz <- if ("spike" %in% names(fp_terms)) {
      fp_terms[["spike"]]
    } else {
      rep(FALSE, nrow(fp_terms))
    }
    
    dec <- suppressWarnings(as.integer(fp_terms[["spike_dec"]]))
    
    saz_decision <- rep("not SAZ", nrow(fp_terms))
    saz_decision[is_saz & !selected] <- "not selected"
    saz_decision[is_saz & selected & dec == 1L] <- "cont + binary"
    saz_decision[is_saz & selected & dec == 2L] <- "continuous only"
    saz_decision[is_saz & selected & dec == 3L] <- "binary only"
    saz_decision[is_saz & selected & !(dec %in% c(1L, 2L, 3L))] <- "unknown"
    
    fp_terms[["spike_dec"]] <- saz_decision
    names(fp_terms)[names(fp_terms) == "spike_dec"] <- "saz_decision"
  }
  
  # Rename only the printed column, not the internal object field.
  names(fp_terms)[names(fp_terms) == "catzero"] <- "catzero_final"
  
  print.data.frame(fp_terms)
  cat("\n")
  
  # Print compact explanations for columns that are easy to misread.
  if (any(c("zero", "catzero_final", "spike", "saz_decision") %in%
          names(fp_terms))) {
    cat("Column notes:\n")
    
    if ("zero" %in% names(fp_terms)) {
      cat("  zero: TRUE means FP transformations are applied only to the positive part of the variable.\n")
    }
    
    if ("catzero_final" %in% names(fp_terms)) {
      cat("  catzero_final: TRUE means a binary indicator for the non-positive component is included in the final model.\n")
    }
    
    if ("spike" %in% names(fp_terms)) {
      cat("  spike: TRUE means SAZ was requested and not reset by the eligibility checks.\n")
    }
    
    if ("saz_decision" %in% names(fp_terms)) {
      cat("  saz_decision: final SAZ form selected; \"not selected\" means SAZ was eligible but the variable was dropped.\n")
    }
    
    cat("\n")
  }
  
  cat(sprintf("MFP algorithm convergence: %s\n", x$convergence_mfp))
  
  # Print model object using underlying print function.
  NextMethod("print", x)
}

#' Helper to assign attributes to a variable undergoing FP-transformation
#'
#' Used in the formula interface to `mfp2()`.
#'
#' @param x A vector representing a continuous variable undergoing
#'   FP-transformation.
#' @param df,alpha,select,shift,scale,center,acdx See [mfp2::mfp2()] for
#'   details.
#' @param zero,catzero,spike Logical formula-interface versions of
#'   \code{zero_vars}, \code{catzero_vars}, and \code{spike_vars}. See
#'   \code{\link{mfp2}}.
#' @param force_max_fp Logical. If \code{TRUE}, the most complex functional form
#'   allowed by \code{df} is selected for this variable. This argument is only
#'   used when \code{criterion = "aic"} or \code{criterion = "bic"} in
#'   \code{mfp2()}. It has no effect when \code{criterion = "pvalue"}, because
#'   the same behaviour can be obtained by setting \code{select = 1} and
#'   \code{alpha = 1}. Default is \code{FALSE}.
#' @param powers Numeric vector of candidate powers to be evaluated for `x`.
#'   Must contain at least two values when supplied. If `NULL`, the default
#'   powers `c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)` are used.
#' @param ... Used in alias `fp2()` to pass arguments.
#'
#' @return The vector `x` with attributes relevant for FP-transformation. All
#'   arguments passed to this function are stored as attributes.
#'
#' @examples
#' xr <- 1:10
#' fp(xr)
#' fp2(xr)
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
               powers = NULL,
               zero = FALSE,
               catzero = FALSE,
               spike = FALSE,
               force_max_fp = FALSE
               ) {
  
  name <- deparse(substitute(x))
  
  # Assert that a factor variable must not be subjected to fp transformation.
  if (is.factor(x)) {
    stop(
      name, " is a factor variable and should not be passed to the fp() function.",
      call. = FALSE
    )
  }
  
  # fp() is part of the formula interface. It stores user-supplied options as
  # attributes that mfp2.formula() later extracts. We validate only shape/type
  # issues that can corrupt attribute extraction. Full range and semantic
  # validation is still performed by mfp2.default() after formula options have
  # been expanded to the final per-variable vectors.
  
  validate_scalar_fp <- function(arg, arg_name, var_name, type,
                                 allow_null = FALSE, allow_na = FALSE) {
    if (is.null(arg)) {
      if (allow_null) return(invisible(TRUE))
      stop(
        sprintf("! In fp(%s), `%s` must not be NULL.", var_name, arg_name),
        call. = FALSE
      )
    }
    
    ok_type <- switch(
      type,
      numeric = is.numeric(arg),
      logical = is.logical(arg),
      stop("Internal error: unsupported validator type.", call. = FALSE)
    )
    
    if (!ok_type || length(arg) != 1L) {
      expected <- switch(
        type,
        logical = "a single TRUE or FALSE value",
        numeric = "a single numeric value",
        type
      )
      
      actual <- paste0(arg, collapse = ", ")
      
      stop(
        sprintf(
          "! In fp(%s), `%s` must be %s.\ni You supplied: %s.",
          var_name, arg_name, expected, actual
        ),
        call. = FALSE
      )
    }
    
    if (!allow_na && anyNA(arg)) {
      expected <- switch(
        type,
        logical = "TRUE or FALSE",
        numeric = "a numeric value",
        type
      )
      
      stop(
        sprintf(
          "! In fp(%s), `%s` must not be NA.\ni Please use %s.",
          var_name, arg_name, expected
        ),
        call. = FALSE
      )
    }
    
    invisible(TRUE)
  }
  
  # Numeric scalar attributes.
  # shift and scale may be NULL because NULL means "use the global mfp2()
  # behavior". They may also be NA only if supplied as scalar NA, meaning
  # automatic shift/scale estimation for this variable.
  validate_scalar_fp(df, "df", name, "numeric", allow_null = FALSE, allow_na = FALSE)
  validate_scalar_fp(alpha, "alpha", name, "numeric", allow_null = FALSE, allow_na = FALSE)
  validate_scalar_fp(select, "select", name, "numeric", allow_null = FALSE, allow_na = FALSE)
  validate_scalar_fp(shift, "shift", name, "numeric", allow_null = TRUE, allow_na = TRUE)
  validate_scalar_fp(scale, "scale", name, "numeric", allow_null = TRUE, allow_na = TRUE)
  
  # Logical scalar attributes.
  validate_scalar_fp(center, "center", name, "logical", allow_null = FALSE, allow_na = FALSE)
  validate_scalar_fp(acdx, "acdx", name, "logical", allow_null = FALSE, allow_na = FALSE)
  validate_scalar_fp(zero, "zero", name, "logical", allow_null = FALSE, allow_na = FALSE)
  validate_scalar_fp(catzero, "catzero", name, "logical", allow_null = FALSE, allow_na = FALSE)
  validate_scalar_fp(spike, "spike", name, "logical", allow_null = FALSE, allow_na = FALSE)
  validate_scalar_fp(force_max_fp, "force_max_fp", name, "logical", allow_null = FALSE, allow_na = FALSE)
  
  # `powers` is a candidate base-power set. Duplicate values are removed here so
  # formula-level powers use the same semantics as the top-level `powers` list.
  # Repeated selected powers are generated later by replacement.
  if (!is.null(powers)) {
    powers <- normalize_fp_power_vector(
      powers,
      context = sprintf("`powers` in fp(%s)", name)
    )
  }
  
  attr(x, "df") <- df
  attr(x, "alpha") <- alpha
  attr(x, "select") <- select
  attr(x, "shift") <- shift
  attr(x, "scale") <- scale
  attr(x, "center") <- center
  attr(x, "acd") <- acdx
  attr(x, "powers") <- powers
  attr(x, "zero") <- zero
  attr(x, "catzero") <- catzero
  attr(x, "spike") <- spike
  attr(x, "force_max_fp") <- force_max_fp
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
#' Character vector of names, ordered as defined by `xorder` in \code{mfp2()}.
#' 
#' @export
get_selected_variable_names <- function(object) {
  nms <- rownames(object$fp_terms)
  nms[object$fp_terms[, "selected"]]
}

#' Assign Degrees of Freedom Based on Variable Cardinality
#'
#' For each column of a predictor matrix, determines an appropriate degree of
#' freedom (df) for fractional polynomial (FP) modelling based on the number
#' of distinct values in that column. Variables with very few unique values
#' cannot support a high-degree FP transformation, so their df is reduced
#' automatically. Accepts either a single scalar default or a per-variable
#' vector, making it suitable for both the scalar and vector branches of
#' \code{mfpi()}'s df handling.
#'
#' @section Assignment rules:
#' Let \eqn{u} denote the number of distinct non-missing values in a column
#' and \eqn{d} the user-supplied df for that column. The assigned df is:
#'
#' \tabular{lll}{
#'   \strong{Unique values} \tab \strong{Assigned df} \tab \strong{Rationale} \cr
#'   \eqn{u \le 3}          \tab \code{1} (linear)    \tab Too few values to estimate a
#'                                                          curve; effectively binary or
#'                                                          ternary. \cr
#'   \eqn{4 \le u \le 5}    \tab \code{min(2, d)}     \tab Enough variation for FP1 but
#'                                                          not FP2; capped at 2. \cr
#'   \eqn{u \ge 6}          \tab \code{d}             \tab Sufficient variation; retain
#'                                                          the requested default. \cr
#' }
#'
#' @param x A numeric matrix. Each column is a predictor variable.
#' @param df_default Either a single positive integer, or an integer vector of
#'   length \code{ncol(x)}, giving the desired df for each column. Must contain
#'   only \code{1} or positive even numbers (\eqn{2m} for FP degree \eqn{m}).
#'   Default is \code{4} (FP2).
#'
#' @return An integer vector of length \code{ncol(x)} giving the assigned df
#'   for each column, in the same order as the columns of \code{x}. Values may
#'   be lower than the corresponding element of \code{df_default} when
#'   cardinality rules override the user-supplied value.
#'
#' @examples
#' x <- cbind(
#'   binary     = c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1),
#'   few_levels = c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5),
#'   continuous = 1:10
#' )
#'
#' # Scalar default: same starting df for all columns
#' assign_df(x, df_default = 4)
#' # binary -> 1, few_levels -> 2, continuous -> 4
#'
#' # Per-variable vector: different starting df per column
#' assign_df(x, df_default = c(1, 4, 2))
#' # binary -> 1 (cardinality), few_levels -> 2 (cap), continuous -> 2
#'
#' @export
assign_df <- function(x, df_default = 4L) {
  
  if (!is.matrix(x)) {
    stop("`x` must be a matrix.", call. = FALSE)
  }
  
  v_names <- colnames(x)
  
  if (is.null(v_names)) {
    stop("x must have names")
  }
  
  n_vars <- ncol(x)
  
  # Expand scalar to per-variable vector; validate length
  if (length(df_default) == 1L) {
    df_default <- rep(as.integer(df_default), n_vars)
  } else if (length(df_default) == n_vars) {
    df_default <- as.integer(df_default)
  } else {
    stop(paste0("`df_default` must be a single integer or a vector of length ",
                n_vars, " (ncol(x)); got length ", length(df_default), "."),
         call. = FALSE)
  }
  
  if (any(df_default < 1L))
    stop("`df_default` must contain only positive integers.", call. = FALSE)
  
  # Count distinct values per column
  nu <- apply(x, 2L, function(col) length(unique(col)))
  
  df <- df_default
  
  # Rule 1: <= 3 unique values -> force linear (df = 1)
  idx_low <- which(nu <= 3L)
  if (length(idx_low) > 0L)
    df[idx_low] <- 1L
  
  # Rule 2: 4-5 unique values -> cap at min(2, requested df per variable)
  idx_mid <- which(nu >= 4L & nu <= 5L)
  if (length(idx_mid) > 0L) {
    df[idx_mid] <- pmin(2L, df_default[idx_mid])
  }
  names(df) <- colnames(x)

  df
}
