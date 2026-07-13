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
#' Stratified Cox models can be specified using the \code{strata} argument,
#' or by including \code{strata()} terms in the model formula when using
#' \code{mfp2.formula()}. In both interfaces, strata are kept as high-level
#' factor/vector values for the final \code{survival::coxph()} formula fit;
#' they are converted to integer codes only inside the low-level
#' \code{survival::coxph.fit()} path. This matches \code{survival::coxph()}
#' semantics and preserves labels for prediction.
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
#' All dummy columns generated by one unordered factor are treated as a single
#' conceptual term. The factor is ordered, selected, retained, or removed as one
#' block, and significance-based decisions use a joint likelihood-ratio test
#' with one degree of freedom per non-reference dummy column.
#'
#' When the argument \code{keep} includes the name of a categorical variable,
#' all dummy variables derived from that factor (e.g., \code{xB}, \code{xC})
#' are automatically retained in the model.
#'
#' ## Ordered Factors
#'
#' Ordered factors, created using \code{ordered()} or
#' \code{factor(..., ordered = TRUE)}, are supported by the formula interface.
#' They are expanded using the contrasts configured for the factor, which are
#' \code{\link[stats]{contr.poly}} by default. For a factor with \eqn{k}
#' levels, the resulting \eqn{k - 1} contrast columns are treated as one
#' conceptual term and are tested, retained, or removed jointly.
#'
#' The resulting test is an omnibus \eqn{k - 1}-degree-of-freedom test of the
#' complete ordered-factor effect. It is not a one-degree-of-freedom linear trend
#' test and does not impose monotonicity. To fit a prespecified one-degree-of-
#' freedom ordinal trend, supply an explicit numeric score and set \code{df = 1},
#' for example \code{fp(severity_score, df = 1)}.
#'
#' All model-matrix columns generated from unordered or ordered factor terms are
#' fixed linear columns: their effective \code{df} is 1, shift is 0, and scale
#' is 1. Fractional-polynomial, ACD, zero, catzero, and spike-at-zero processing
#' does not apply to categorical contrast blocks.
#'
#' ## Prediction for Factors
#'
#' For formula fits, \code{predict.mfp2()} accepts ordinary \code{newdata}
#' containing the original factor variables. The fitted factor levels and
#' contrasts are reused to reconstruct the same model-matrix columns as at fit
#' time. New or unknown factor levels are rejected by the usual
#' \code{model.frame()} / \code{model.matrix()} checks.
#'
#' Simple wrappers such as \code{factor(x9)} and \code{ordered(x9)} use
#' \code{x9} as the conceptual term name in selection output, \code{keep},
#' \code{term_to_columns}, and term prediction. The generated design columns
#' and final coefficient names retain their standard \code{model.matrix()}
#' names, such as \code{factor(x9)2}.
#'
#' ## The \code{mfp2.default()} Method
#'
#' The default method \code{mfp2.default()} does not perform any automatic
#' expansion of factor variables. It assumes that the input design matrix
#' \code{x} consists solely of numeric predictors. Therefore, when calling
#' \code{mfp2.default()} directly, the user must create any required dummy or
#' contrast-coded variables manually before model fitting. When one or more
#' columns encode one conceptual categorical term, supply their shared mapping
#' through \code{term_groups} so they are selected and tested jointly.
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
#' The formula interface expands unordered and ordered categorical predictors
#' using \code{\link[stats]{model.matrix}}. All columns generated from one
#' factor term are treated as one fixed linear block and tested jointly. When
#' \code{keep} names the original factor term, its complete contrast block is
#' retained and protected from exclusion during model selection.
#'
#' Variable-name arguments such as \code{keep}, \code{zero_vars},
#' \code{catzero_vars}, \code{spike_vars}, and \code{acdx} must refer to
#' existing column names in \code{x}. Unknown names are treated as errors to
#' avoid silently ignoring misspelled model specifications.
#'
#' The formula may contain \code{strata()} or \code{survival::strata()}
#' terms for stratified Cox models and \code{offset()} or
#' \code{stats::offset()} terms for model offsets. Namespace-qualified formula
#' specials are normalized internally so that they retain standard formula
#' special semantics. For a formula using \code{.}, such
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
#' To avoid this conflict, use \code{fp2()} instead of \code{fp()} inside the
#' \code{mfp2()} formula whenever both packages are loaded in the same
#' session, for example:
#'
#' \preformatted{
#' library(mfp)
#' library(mfp2)
#' fit <- mfp2(y ~ fp2(x1) + fp2(x2), data = dat)
#' }
#'
#' \code{fp2()} is a simple alias for \code{fp()} (see \code{\link{fp}}) and
#' accepts exactly the same arguments; it exists solely so that its name does
#' not collide with `mfp`'s own \code{fp()}.
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
#' @param term_groups For \code{mfp2.default()}, an optional named list mapping
#'   conceptual categorical term names to one or more design-matrix columns, for
#'   example \code{list(race = c("raceB", "raceC"))}. A one-column mapping
#'   supports binary factors whose conceptual name differs from their dummy
#'   column name. Columns not listed are treated as identity singleton terms.
#'   Explicitly mapped terms are selected and tested
#'   jointly and must be linear-only (effective \code{df = 1}) with no ACD,
#'   zero, catzero, or spike-at-zero handling. The formula method constructs
#'   this mapping automatically from \code{model.matrix()} and does not expose
#'   a separate user-facing argument.
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
#' @param center a logical determining whether the fractional polynomial (or
#' ACD) transformed function of each variable is centered before final model
#' fitting. The default `TRUE` centers by the mean of the transformed values;
#' the exception is binary covariates, which are not power-transformed and are
#' instead centered by subtracting the lower of their two distinct values.
#' See Details section below.
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
#' @param strata For Cox models, optional stratification values supplied to
#'   the default/matrix interface. May be a vector or factor with one value per
#'   observation, or a matrix/data frame with one row per observation for
#'   multiple stratification variables. Multiple columns are combined into one
#'   high-level strata object. The final Cox formula fit stores this as an
#'   internal \code{strata(strata_)} term; ordinary vector/factor strata are not
#'   converted to integer codes until the low-level \code{coxph.fit()} path.
#'   In the formula interface, use \code{strata()} or
#'   \code{survival::strata()} in the formula. Default is \code{NULL}, giving
#'   an unstratified Cox model. See \code{survival::coxph()} for details.
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
#' \dontrun{
#' # Unordered factor: all treatment-contrast columns are selected jointly
#' set.seed(1)
#' d_factor <- data.frame(
#'   y = rnorm(60),
#'   age = runif(60, 20, 80),
#'   group = factor(rep(c("A", "B", "C"), each = 20))
#' )
#' fit_factor <- mfp2(y ~ age + group, data = d_factor, keep = "group", verbose = FALSE)
#' predict(fit_factor, newdata = d_factor[1:4, ])
#'
#' # Ordered factor: default polynomial contrasts form one joint term
#' d_factor$severity <- ordered(
#'   rep(c("mild", "moderate", "severe"), length.out = nrow(d_factor)),
#'   levels = c("mild", "moderate", "severe")
#' )
#' fit_ordered <- mfp2(y ~ age + severity, data = d_factor, keep = "severity", verbose = FALSE)
#' predict(
#'   fit_ordered,
#'   newdata = data.frame(
#'     age = c(40, 40),
#'     severity = ordered(
#'       c("mild", "severe"),
#'       levels = levels(d_factor$severity)
#'     )
#'   )
#' )
#'
#' # Matrix interface: declare manually created dummy columns as one term
#' x_grouped <- cbind(
#'   age = d_factor$age,
#'   groupB = as.integer(d_factor$group == "B"),
#'   groupC = as.integer(d_factor$group == "C")
#' )
#' fit_grouped <- mfp2(
#'   x_grouped,
#'   d_factor$y,
#'   term_groups = list(group = c("groupB", "groupC")),
#'   keep = "group",
#'   verbose = FALSE
#' )
#'
#' }
#'
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
#' \item term_to_columns: a complete named list mapping every conceptual term
#' to its raw design-matrix column or columns. Simple factor wrappers use the
#' source variable name, while the raw columns retain model.matrix() names.
#' \item formula_factor_info: for formula fits with simple factor main effects,
#' fitted level and design-block metadata used by term and contrast prediction.
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
#' \code{summary.mfp2()}, \code{coef.mfp2()}, \code{predict.mfp2()}, \code{fp()},
#' \code{fp2()}
#'
#' @export
mfp2 <- function(x, ...){
  UseMethod("mfp2", x)
}

#' Identify Formula Terms Containing Factors
#'
#' Uses the terms-object factor incidence matrix and the evaluated model frame to
#' identify formula terms that contain at least one unordered or ordered factor.
#' The returned names are original formula term labels. Callers may map simple
#' factor wrappers to source-variable conceptual names before constructing the
#' term-to-column lookup.
#'
#' @param terms_object A terms object used to build the predictor model matrix.
#' @param model_frame The evaluated model frame.
#'
#' @return Character vector of factor-containing formula term labels.
#'
#' @keywords internal
#' @noRd
identify_formula_factor_terms <- function(terms_object, model_frame) {
  incidence <- attr(terms_object, "factors")
  
  if (is.null(incidence) || nrow(incidence) == 0L || ncol(incidence) == 0L) {
    return(character(0L))
  }
  
  variable_labels <- rownames(incidence)
  is_factor_variable <- vapply(
    variable_labels,
    function(variable) {
      variable %in% names(model_frame) && is.factor(model_frame[[variable]])
    },
    logical(1L)
  )
  
  if (!any(is_factor_variable)) {
    return(character(0L))
  }
  
  factor_incidence <- incidence[is_factor_variable, , drop = FALSE]
  colnames(incidence)[colSums(factor_incidence != 0) > 0L]
}


#' Extract the Source Variable from a Simple Factor Formula Term
#'
#' Recognises a bare variable or a direct call to factor(), ordered(),
#' as.factor(), or as.ordered() whose first argument is a single variable. The
#' returned source name is used as the conceptual MFP term name, while the
#' original formula expression remains available through the fitted terms
#' object for model-matrix reconstruction.
#'
#' @param term_label Character scalar containing a formula term label.
#'
#' @return Character scalar source-variable name, or NULL for terms that cannot
#'   be reduced safely to one source variable.
#'
#' @keywords internal
#' @noRd
formula_factor_source_name <- function(term_label) {
  expr <- tryCatch(str2lang(term_label), error = function(e) NULL)
  
  if (is.null(expr)) {
    return(NULL)
  }
  
  if (is.symbol(expr)) {
    return(as.character(expr))
  }
  
  if (!is.call(expr) || length(expr) < 2L) {
    return(NULL)
  }
  
  call_head <- expr[[1L]]
  function_name <- if (is.symbol(call_head)) {
    as.character(call_head)
  } else if (
    is.call(call_head) && length(call_head) == 3L &&
    as.character(call_head[[1L]]) %in% c("::", ":::") &&
    is.symbol(call_head[[3L]])
  ) {
    as.character(call_head[[3L]])
  } else {
    NULL
  }
  
  if (is.null(function_name) ||
      !function_name %in% c("factor", "ordered", "as.factor", "as.ordered") ||
      !is.symbol(expr[[2L]])) {
    return(NULL)
  }
  
  as.character(expr[[2L]])
}


#' Resolve User-Facing Conceptual Names for Formula Terms
#'
#' Formula term labels are retained for model.frame()/model.matrix()
#' reconstruction, but simple factor wrappers are represented to the MFP engine
#' by their source variable name. For example, factor(x9) is represented by the
#' conceptual term x9 while its design columns remain factor(x9)2 and
#' factor(x9)3.
#'
#' @param terms_object Predictor terms object.
#' @param model_frame Evaluated model frame.
#'
#' @return Named character vector mapping original formula labels to conceptual
#'   term names.
#'
#' @keywords internal
#' @noRd
formula_conceptual_term_names <- function(terms_object, model_frame) {
  term_labels <- attr(terms_object, "term.labels")
  out <- stats::setNames(term_labels, term_labels)
  factor_labels <- identify_formula_factor_terms(terms_object, model_frame)
  
  for (label in intersect(factor_labels, term_labels)) {
    source_name <- formula_factor_source_name(label)
    if (!is.null(source_name)) {
      out[[label]] <- source_name
    }
  }
  
  duplicated_names <- unique(unname(out)[
    duplicated(unname(out)) | duplicated(unname(out), fromLast = TRUE)
  ])
  
  if (length(duplicated_names) > 0L) {
    details <- vapply(
      duplicated_names,
      function(name) {
        labels <- names(out)[unname(out) == name]
        paste0(name, " <- ", paste(labels, collapse = ", "))
      },
      character(1L)
    )
    stop(
      "! Multiple formula terms resolve to the same conceptual variable: ",
      paste(details, collapse = "; "), ".\n",
      "i Use each source variable only once in the model formula.",
      call. = FALSE
    )
  }
  
  out
}


#' Build a Formula Term-to-Column Lookup
#'
#' @param x_columns Design-matrix column names after intercept removal.
#' @param assign Model-matrix assign vector aligned with x_columns.
#' @param term_labels Original formula term labels.
#' @param conceptual_names Named mapping from formula labels to conceptual names.
#'
#' @return Named list mapping conceptual terms to their design columns in formula
#'   order.
#'
#' @keywords internal
#' @noRd
build_formula_term_to_columns <- function(x_columns,
                                          assign,
                                          term_labels,
                                          conceptual_names) {
  out <- lapply(
    seq_along(term_labels),
    function(index) x_columns[assign == index]
  )
  names(out) <- unname(conceptual_names[term_labels])
  out[lengths(out) > 0L]
}


#' Test Whether a Conceptual Term Uses an Explicit Column Mapping
#'
#' A mapped term either spans multiple raw columns or maps to one raw column
#' whose name differs from the conceptual term. The latter is the usual design
#' for a two-level factor such as treatment -> treatmentB.
#'
#' @param term Conceptual term name.
#' @param columns Character vector of raw columns.
#'
#' @return Logical scalar.
#'
#' @keywords internal
#' @noRd
term_uses_column_mapping <- function(term, columns) {
  length(columns) != 1L || !identical(columns[[1L]], term)
}


#' Identify Explicitly Mapped Conceptual Terms
#'
#' @param term_to_columns Complete conceptual-term lookup.
#'
#' @return Named logical vector aligned with term_to_columns.
#'
#' @keywords internal
#' @noRd
mapped_term_flags <- function(term_to_columns) {
  vapply(
    names(term_to_columns),
    function(term) term_uses_column_mapping(term, term_to_columns[[term]]),
    logical(1L)
  )
}


#' Build Formula Factor Prediction Metadata
#'
#' Stores the exact fitted design row associated with each observed level of a
#' simple factor main effect. Metadata are keyed by the conceptual source
#' variable name, while the raw design-column names remain unchanged.
#'
#' @param factor_terms Character vector of original factor-containing formula
#'   labels returned by identify_formula_factor_terms().
#' @param terms_object The predictor terms object.
#' @param model_frame The evaluated model frame.
#' @param term_name_map Named mapping from original formula labels to conceptual
#'   term names.
#' @param term_to_columns Conceptual term-to-design-column lookup.
#' @param x Numeric design matrix after any fp() column renaming.
#'
#' @return Named list of metadata for simple factor main effects. More complex
#'   factor-containing terms, such as interactions, are omitted because their
#'   design rows also depend on the interacting variables.
#'
#' @keywords internal
#' @noRd
build_formula_factor_info <- function(factor_terms,
                                      terms_object,
                                      model_frame,
                                      term_name_map,
                                      term_to_columns,
                                      x) {
  incidence <- attr(terms_object, "factors")
  
  if (length(factor_terms) == 0L || is.null(incidence)) {
    return(list())
  }
  
  out <- list()
  
  for (formula_term in intersect(factor_terms, colnames(incidence))) {
    participating <- rownames(incidence)[incidence[, formula_term] != 0]
    factor_variables <- participating[
      vapply(
        participating,
        function(variable) {
          variable %in% names(model_frame) && is.factor(model_frame[[variable]])
        },
        logical(1L)
      )
    ]
    
    conceptual_term <- unname(term_name_map[[formula_term]])
    if (length(participating) != 1L || length(factor_variables) != 1L ||
        is.null(conceptual_term) || !conceptual_term %in% names(term_to_columns)) {
      next
    }
    
    frame_variable <- factor_variables[[1L]]
    values <- model_frame[[frame_variable]]
    columns <- term_to_columns[[conceptual_term]]
    level_names <- levels(values)
    
    design_by_level <- matrix(
      NA_real_,
      nrow = length(level_names),
      ncol = length(columns),
      dimnames = list(level_names, columns)
    )
    
    for (level in level_names) {
      row_index <- which(as.character(values) == level)[1L]
      if (!is.na(row_index)) {
        design_by_level[level, ] <- x[row_index, columns, drop = TRUE]
      }
    }
    
    out[[conceptual_term]] <- list(
      levels = level_names,
      ordered = is.ordered(values),
      columns = columns,
      design_by_level = design_by_level
    )
  }
  
  out
}


#' Normalize Design-Matrix Columns into Conceptual Terms
#'
#' Validates an optional grouped-term specification and returns a complete map
#' from conceptual term names to raw design-matrix columns. Unmentioned columns
#' are represented as singleton terms. Term order follows the first occurrence
#' of each term in the original design-matrix column order.
#'
#' @param x_columns Character vector of design-matrix column names.
#' @param term_groups Optional named list mapping conceptual terms to raw column
#'   names.
#'
#' @return A named list mapping every conceptual term to one or more raw columns.
#' @keywords internal
#' @noRd
normalize_term_groups <- function(x_columns, term_groups = NULL) {
  if (is.null(term_groups)) {
    return(stats::setNames(as.list(x_columns), x_columns))
  }
  
  if (!is.list(term_groups)) {
    stop("! `term_groups` must be NULL or a named list.", call. = FALSE)
  }
  
  group_names <- names(term_groups)
  if (is.null(group_names) || length(group_names) != length(term_groups) ||
      anyNA(group_names) || any(!nzchar(group_names))) {
    stop("! `term_groups` must have one non-empty name per list element.",
         call. = FALSE)
  }
  
  if (anyDuplicated(group_names)) {
    duplicated_names <- unique(group_names[duplicated(group_names)])
    stop(
      sprintf(
        "! `term_groups` contains duplicated term name(s): %s.",
        paste(duplicated_names, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  invalid_element <- vapply(
    term_groups,
    function(cols) {
      !is.character(cols) || length(cols) < 1L || anyNA(cols) ||
        any(!nzchar(cols))
    },
    logical(1L)
  )
  if (any(invalid_element)) {
    stop(
      sprintf(
        "! Every `term_groups` element must contain at least one non-missing, non-empty column name. Omit identity singleton mappings; they are added automatically. Invalid term(s): %s.",
        paste(names(term_groups)[invalid_element], collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  referenced_columns <- unlist(term_groups, use.names = FALSE)
  missing_columns <- setdiff(referenced_columns, x_columns)
  if (length(missing_columns) > 0L) {
    stop(
      sprintf(
        "! `term_groups` references unknown column(s) in `x`: %s.",
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  duplicated_columns <- unique(referenced_columns[duplicated(referenced_columns)])
  if (length(duplicated_columns) > 0L) {
    stop(
      sprintf(
        "! A column may appear in only one `term_groups` entry. Duplicated column(s): %s.",
        paste(duplicated_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  colliding_names <- intersect(group_names, x_columns)
  invalid_collisions <- colliding_names[!vapply(
    colliding_names,
    function(term) term %in% term_groups[[term]],
    logical(1L)
  )]
  if (length(invalid_collisions) > 0L) {
    stop(
      sprintf(
        "! A grouped term name may match a column name only when that column is a member of the same group. Invalid collision(s): %s.",
        paste(invalid_collisions, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  column_to_term <- stats::setNames(x_columns, x_columns)
  for (term in group_names) {
    column_to_term[term_groups[[term]]] <- term
  }
  
  term_names <- unique(unname(column_to_term[x_columns]))
  out <- lapply(
    term_names,
    function(term) x_columns[unname(column_to_term[x_columns]) == term]
  )
  names(out) <- term_names
  out
}

#' Validate a Setting for Explicitly Mapped Terms
#'
#' Checks a per-column setting for every explicitly mapped term, including a
#' binary factor represented by one non-identity dummy column, and reports the
#' first offending term with the setting name.
#'
#' @param term_to_columns Complete conceptual-term lookup.
#' @param values Named vector indexed by raw design-matrix columns.
#' @param setting Character setting name used in the error message.
#' @param predicate Function returning a logical vector of valid values.
#' @param requirement Character description of the required value.
#'
#' @return Invisibly returns \code{TRUE}.
#' @keywords internal
#' @noRd
validate_grouped_term_setting <- function(term_to_columns,
                                          values,
                                          setting,
                                          predicate,
                                          requirement) {
  grouped_terms <- names(term_to_columns)[mapped_term_flags(term_to_columns)]
  for (term in grouped_terms) {
    cols <- term_to_columns[[term]]
    valid <- predicate(values[cols])
    if (length(valid) != length(cols) || anyNA(valid) || !all(valid)) {
      stop(
        sprintf(
          "! Grouped term '%s' requires `%s` to be %s for every member column: %s.",
          term, setting, requirement, paste(cols, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

#' Collapse a Per-Column Option to One Value per Conceptual Term
#'
#' For one-column terms the corresponding column value is returned. For
#' multi-column terms all member-column values must be identical because the
#' selection engine has one setting per conceptual term.
#'
#' @param values Named atomic vector or list indexed by raw columns.
#' @param term_to_columns Complete conceptual-term lookup.
#' @param setting Character setting name used in error messages.
#'
#' @return An object of the same broad type as \code{values}, named by term.
#' @keywords internal
#' @noRd
collapse_option_to_terms <- function(values, term_to_columns, setting) {
  is_list <- is.list(values)
  collapsed <- lapply(names(term_to_columns), function(term) {
    cols <- term_to_columns[[term]]
    vals <- values[cols]
    if (length(cols) > 1L) {
      reference <- unname(vals[[1L]])
      same <- vapply(
        vals[-1L],
        function(value) identical(unname(value), reference),
        logical(1L)
      )
      if (length(same) > 0L && !all(same)) {
        stop(
          sprintf(
            "! Grouped term '%s' must use the same `%s` value for all member columns: %s.",
            term, setting, paste(cols, collapse = ", ")
          ),
          call. = FALSE
        )
      }
    }
    vals[[1L]]
  })
  names(collapsed) <- names(term_to_columns)
  
  if (is_list) {
    return(collapsed)
  }
  
  out <- unlist(collapsed, use.names = FALSE)
  names(out) <- names(term_to_columns)
  out
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
                         ...,
                         term_groups = NULL) {
  
  # mfp2.default() does not fit the model itself. It validates and normalizes
  # every argument, derives shift/scale/df settings and zero/catzero/spike/acd
  # flags for each predictor, transforms `x` accordingly, and then delegates
  # the actual multivariable FP selection and fitting to fit_mfp().
  
  # Step 1: Capture the call and rename public-facing arguments -----------------
  cl <- match.call()
  
  # Display the public generic name rather than the S3 method name.
  cl[[1L]] <- quote(mfp2)
  
  # Public interface:
  # zero_vars, catzero_vars and spike_vars are character vectors of variable names.
  #
  # Internal representation:
  # The existing internal names zero, catzero and spike are retained because
  # they are converted below to named logical vectors over the columns of x.
  zero <- zero_vars
  catzero <- catzero_vars
  spike <- spike_vars
  
  # Step 2: Resolve multiple-choice arguments to their single selected value ----
  criterion <- match.arg(criterion)
  xorder    <- match.arg(xorder)
  ties      <- match.arg(ties)
  
  # Step 3: Validate family and the input matrix `x` ----------------------------
  family_info <- normalize_family_argument(
    family,
    family_arg = deparse(substitute(family))
  )
  
  family        <- family_info$family
  family_string <- family_info$family_string
  
  # `x` must be a plain numeric matrix: mfp2.default() does not expand factors,
  # so any categorical predictors must already be coded as dummy columns.
  if (!is.matrix(x)) {
    stop("! x must be a matrix", call. = FALSE)
  }
  
  # character values would silently coerce during FP transformation, so reject
  # them explicitly and point the user to dummy-coding instead.
  if (any(is.character(x))) {
    stop("! x contains characters values.\n",
         "i Please convert categorical variables to dummy variables.", 
         call. = FALSE)
  }
  
  # nobs = number of observations (rows), nvars = number of predictors (columns).
  # These two counts are used throughout the rest of the function to validate
  # the length of vector arguments such as df, select, alpha, shift, scale.
  np <- dim(x)
  nobs <- as.integer(np[1])
  nvars <- as.integer(np[2])
  
  # dim() returns NULL for objects without dimensions (e.g. a plain vector),
  # which would otherwise make nobs/nvars silently become NA above.
  if (is.null(np)) {
    stop("! The dimensions of x must not be missing.\n",
         "i Please make sure that x is a matrix with at least one row and column.", 
         call. = FALSE)
  }
  
  # column names are required because variable-name arguments (keep, zero_vars,
  # catzero_vars, spike_vars, acdx, powers) are matched against them downstream.
  vnames <- colnames(x)
  if (is.null(vnames)) {
    stop("! The column names of x must not be missing.\n",
         "i Please set column names for x.",
         call. = FALSE)
  }
  
  if (!is.null(term_groups) && anyDuplicated(vnames)) {
    stop("! The column names of `x` must be unique when `term_groups` is supplied.",
         call. = FALSE)
  }
  
  # Build the complete conceptual-term lookup once. This is a singleton map
  # when term_groups = NULL, preserving the historical one-column-per-variable
  # behavior exactly.
  term_to_columns <- normalize_term_groups(vnames, term_groups)
  column_to_term <- stats::setNames(
    rep(names(term_to_columns), lengths(term_to_columns)),
    unlist(term_to_columns, use.names = FALSE)
  )
  
  if (!is.numeric(x)) {
    stop(
      "! `x` must be a numeric matrix.",
      sprintf("i Current storage mode is: %s.", typeof(x)),
      call. = FALSE
    )
  }
  
  # missing data is not supported: FP transformations and the backfitting
  # algorithm assume a complete design matrix. Users must remove/impute NAs
  # beforehand rather than relying on implicit row deletion.
  if (anyNA(x)) {
    stop("! x must not contain any NA (missing data).\n",   
         "i Please remove any missing data before passing x to this function.",
         call. = FALSE)
  }
  
  # Inf/-Inf would silently break shift/scale estimation and FP power fitting.
  if (any(!is.finite(x))) {
    stop(
      "! `x` must contain only finite, non-missing numeric values.",
      call. = FALSE
    )
  }
  
  # Step 4: Validate and normalize `subset` --------------------------------------
  # `subset` may be supplied as either a logical mask (one value per row of x)
  # or a vector of row indices; both forms are normalized to integer indices.
  # Note: subsetting is applied later, after shift/scale estimation (see the
  # "Details on subset" documentation section), not at this validation step.
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
    
    # A minimum of 5 observations is required so that model fitting and the
    # FP selection tests below have a chance of being numerically well-defined.
    if (length(subset) < 5L) {
      stop(
        "! The selected subset is too small (<5) to fit an mfp model.",
        sprintf("i The number of selected observations is %d.", length(subset)),
        call. = FALSE
      )
    }
  } 
  
  # Step 5: Validate `weights` and `offset` --------------------------------------
  # weights: optional observation weights, one per row of x.
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
  
  # offset: optional known linear-predictor term, one per row of x
  # (e.g. log of exposure time for a Poisson model).
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
  
  # Step 6: Validate scalar and per-predictor vector options --------------------
  
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
  
  # Step 7: Validate variable-name arguments --------------------------------------
  # These arguments define the model specification. Unknown names are treated as
  # errors because silently dropping them can hide spelling mistakes, e.g.
  # zero_vars = "expsoure" instead of "exposure".
  
  if (!is.null(keep)) {
    if (!is.character(keep) || anyNA(keep) || any(!nzchar(keep))) {
      stop("! `keep` must be a character vector without missing or empty names.",
           call. = FALSE)
    }
    unknown_keep <- setdiff(keep, unique(c(vnames, names(term_to_columns))))
    if (length(unknown_keep) > 0L) {
      stop(
        sprintf("! Unknown variable(s) in keep: %s.",
                paste(unknown_keep, collapse = ", ")),
        call. = FALSE
      )
    }
  }
  validate_variable_names(zero, "zero_vars", vnames)
  validate_variable_names(catzero, "catzero_vars", vnames)
  validate_variable_names(acdx, "acdx", vnames)
  validate_variable_names(spike, "spike_vars", vnames)
  
  # Step 8: Validate `df` ---------------------------------------------------------
  # df must translate into a valid FP degree: 1 (linear) or an even positive
  # integer 2, 4, 6, ... (representing FP1, FP2, FP3, ...).
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
  
  # Step 9: Validate the response `y` and family-specific auxiliary inputs -------
  
  # The F-test correction for small-sample normal error models only applies to
  # Gaussian models; silently revert to the Chi-square test for other families.
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
    strata_len <- if (is.vector(strata) || is.factor(strata)) {
      length(strata)
    } else {
      NROW(strata)
    }
    
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
    
    if (anyNA(strata)) {
      stop(
        "! `strata` must not contain missing values.",
        call. = FALSE
      )
    }
  }
  
  # Step 10: Validate spike-at-zero proportion and the candidate power list -----
  # min_saz_component_prop is the minimum share of observations required in
  # both the zero component and the positive component for a variable
  # requested via spike_vars to remain SAZ-eligible (see resolve_saz_eligibility()
  # below).
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
  
  # Step 11: Apply defaults and expand scalars to per-predictor vectors ---------
  
  # Default weights (all observations equally weighted) and offset (no offset).
  # has_offset is recorded before defaulting so the fitted object can report
  # whether the user actually supplied an offset (see @return: has_offset).
  if (is.null(weights)) {
    weights <- rep.int(1, nobs)
  }
  
  has_offset <- !is.null(offset)
  
  if (is.null(offset)) {
    offset <- rep.int(0, nobs)
  }
  
  # Expand scalar select/alpha to one value per predictor.
  if (length(select) == 1) {
    select <- rep(select, nvars)
  } 
  
  if (length(alpha) == 1) {
    alpha <- rep(alpha, nvars)
  } 
  
  # Default shift values: NA marks a variable for automatic shift estimation
  # further below (find_shift_factor()); a supplied value is used as-is.
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
  
  # Expand scalar center to one value per predictor.
  if (length(center) == 1) {
    center <- rep(center, nvars)    
  }
  
  # Expand scalar force_max_fp to one value per predictor.
  if (length(force_max_fp) == 1L) {
    force_max_fp <- rep(force_max_fp, nvars)
  }
  force_max_fp <- setNames(force_max_fp, vnames)
  
  # force_max_fp only affects information-criterion-based selection; warn the
  # user if they set it while still using the default p-value criterion.
  if (criterion == "pvalue" && any(force_max_fp)) {
    warning(
      "i `force_max_fp` has no effect when criterion = 'pvalue'.",
      "i Use select = 1 and alpha = 1 to force the most complex FP under p-value selection.",
      call. = FALSE
    )
  }
  
  # Fall back to the family-appropriate default fitting-control object
  # (glm.control() or coxph.control()) if the user did not supply one.
  if (is.null(control)) {
    if (family_string == "cox") {
      control <- survival::coxph.control()
    } else {
      control <- stats::glm.control()
    }
  }
  
  # Step 12: Resolve zero / catzero / acdx / spike flags for each predictor -----
  # Convert the character-vector user interface (zero_vars, catzero_vars, ...)
  # into named logical vectors over all columns of x, so downstream code can
  # simply index by variable name.
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
  
  # zero_vars only makes sense for variables that actually contain nonpositive
  # values; warn and reset the flag for any variable that is already all-positive.
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
  
  # Same rationale as above, applied to catzero_vars.
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
  
  # acdx: character vector of variable names -> named logical vector over vnames.
  if (is.null(acdx)) {
    acdx <- setNames(rep(FALSE, nvars), vnames)
  } else {
    acdx <- unique(acdx)
    
    acdx_vec <- setNames(rep(FALSE, nvars), vnames)
    acdx_vec[acdx] <- TRUE
    acdx <- acdx_vec
  }
  
  # spike_vars: character vector of variable names -> named logical vector over vnames.
  if (is.null(spike)) {
    spike <- setNames(rep(FALSE, nvars), vnames)
  } else {
    spike <- unique(spike)
    
    spike_input_vars <- spike
    spike <- setNames(rep(FALSE, nvars), vnames)
    spike[spike_input_vars] <- TRUE
  }
  
  # Explicitly mapped terms are categorical design blocks, including binary
  # factors represented by one non-identity dummy column. They cannot use any
  # continuous-variable extension or structural-zero representation.
  validate_grouped_term_setting(
    term_to_columns, acdx, "acdx", function(v) !v, "FALSE"
  )
  validate_grouped_term_setting(
    term_to_columns, zero, "zero_vars", function(v) !v, "FALSE"
  )
  validate_grouped_term_setting(
    term_to_columns, catzero, "catzero_vars", function(v) !v, "FALSE"
  )
  validate_grouped_term_setting(
    term_to_columns, spike, "spike_vars", function(v) !v, "FALSE"
  )
  
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
  
  # Step 13: Resolve spike-at-zero eligibility -----------------------------------
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
  
  # Step 14: Assign effective df per predictor -----------------------------------
  # A scalar df is adjusted per-variable based on cardinality via assign_df();
  # an explicit per-variable df vector is instead only downgraded to 1 (linear)
  # for variables with too few unique values to support a curve.
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
  
  validate_grouped_term_setting(
    term_to_columns,
    stats::setNames(df.list, vnames),
    "df",
    function(v) v == 1,
    "1 (linear-only)"
  )
  
  # Step 15: Cascade spike -> catzero -> zero, then reset shift for affected vars
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
  
  # Step 16: Shift and scale the design matrix -----------------------------------
  # Shift each column of x so that fractional-polynomial powers (which require
  # positive values) can be computed; see "Details on shifting, scaling,
  # centering" in the documentation above.
  x <- sweep(x, 2, shift, "+")
  
  # Scaling must be estimated after shifting, since find_scale_factor() operates
  # on the already-shifted (positive) values.
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
  
  # Step 17: Verify shifted values are strictly positive where FPs are used -----
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
  
  # Apply the (possibly estimated) scale factors, after shifting and after the
  # positivity check above.
  x <- sweep(x, 2, scale, "/")
  
  # Step 18: Build the Cox stratification object --------------------------------
  # Mimic survival::coxph():
  #   - keep the evaluated strata object at the model-frame level
  #   - convert to integer only immediately before coxph.fit()
  strata_keep <- strata
  
  if (family_string == "cox" && !is.null(strata_keep)) {
    if (!isTRUE(attr(strata_keep, "mfp2_strata_keep"))) {
      strata_keep <- if (is.matrix(strata_keep) || is.data.frame(strata_keep)) {
        do.call(
          survival::strata,
          c(as.list(as.data.frame(strata_keep)), list(shortlabel = TRUE))
        )
      } else {
        strata_keep
      }
    }
  }
  
  # Step 19: Apply `subset`, after shift/scale estimation ------------------------
  # This ordering (estimate shift/scale on the full data, then subset) is what
  # makes mfp2()'s `subset` different from subsetting the data beforehand;
  # see "Details on the subset argument" above.
  if (!is.null(subset)) {
    
    x <- x[subset, , drop = FALSE]
    
    if (is.matrix(y)) {
      y <- y[subset, , drop = FALSE]
    } else {
      y <- y[subset]
    }
    
    weights <- weights[subset]
    offset <- offset[subset]
    if (!is.null(strata_keep)) {
      strata_keep <- strata_keep[subset]
    }
  }
  
  # Collapse raw-column options to one value per conceptual term only after all
  # column-level preprocessing has completed. Grouped members must agree on
  # settings that have a single term-level meaning.
  alpha_term <- collapse_option_to_terms(
    stats::setNames(alpha, vnames), term_to_columns, "alpha"
  )
  select_term <- collapse_option_to_terms(
    stats::setNames(select, vnames), term_to_columns, "select"
  )
  df_term <- collapse_option_to_terms(
    stats::setNames(df.list, vnames), term_to_columns, "df"
  )
  center_term <- collapse_option_to_terms(
    stats::setNames(center, vnames), term_to_columns, "center"
  )
  shift_term <- collapse_option_to_terms(
    stats::setNames(shift, vnames), term_to_columns, "shift"
  )
  scale_term <- collapse_option_to_terms(
    stats::setNames(scale, vnames), term_to_columns, "scale"
  )
  acdx_term <- collapse_option_to_terms(acdx, term_to_columns, "acdx")
  zero_term <- collapse_option_to_terms(zero, term_to_columns, "zero_vars")
  catzero_term <- collapse_option_to_terms(catzero, term_to_columns, "catzero_vars")
  spike_term <- collapse_option_to_terms(spike, term_to_columns, "spike_vars")
  force_max_fp_term <- collapse_option_to_terms(
    force_max_fp, term_to_columns, "force_max_fp"
  )
  
  powers_term <- lapply(names(term_to_columns), function(term) {
    cols <- term_to_columns[[term]]
    if (term_uses_column_mapping(term, cols)) 1 else power_list[[cols]]
  })
  names(powers_term) <- names(term_to_columns)
  
  keep_term <- if (is.null(keep)) {
    NULL
  } else {
    unique(vapply(
      keep,
      function(value) {
        if (value %in% names(term_to_columns)) value else unname(column_to_term[[value]])
      },
      character(1L)
    ))
  }
  
  # Step 20: Fit the multivariable FP model --------------------------------------
  # Fail early if the design matrix is rank-deficient, rather than letting
  # glm()/coxph() fail deep inside the backfitting cycles with a less clear error.
  validate_default_design_rank(
    x = x,
    intercept = family_string != "cox"
  )
  
  # Delegate the actual variable selection, FP degree/power selection, and
  # model fitting (including SAZ and ACD handling) to fit_mfp().
  fit <- fit_mfp(
    x = x, y = y, 
    weights = weights, offset = offset, cycles = cycles, 
    scale = scale_term, shift = shift_term, df = df_term, center = center_term, 
    family = family, family_string = family_string, criterion = criterion, 
    select = select_term, alpha = alpha_term, keep = keep_term, xorder = xorder, 
    powers = powers_term, method = ties, strata = strata_keep, nocenter = nocenter, 
    acdx = acdx_term, ftest = ftest, force_max_fp = force_max_fp_term,
    control = control, zero = zero_term, catzero = catzero_term, spike = spike_term,
    min_saz_component_prop = min_saz_component_prop, 
    saz_pre_resolved = TRUE,
    term_to_columns = term_to_columns,
    has_offset = has_offset,
    verbose = verbose
  )
  
  # Step 21: Attach mfp2-specific metadata to the fitted object -----------------
  fit$call_mfp <- cl
  fit$family <- family
  fit$family_string <- family_string
  fit$offset <- offset
  fit$has_offset <- has_offset
  
  fit
}


#' @describeIn mfp2 Provides formula interface for `mfp2`.
#' @importFrom utils modifyList
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
  # mfp2.formula() translates a formula + data.frame specification into the
  # matrix/vector inputs expected by mfp2.default(): it expands categorical
  # predictors via model.matrix(), extracts per-variable fp() settings, and
  # then calls mfp2.default() to perform the actual fitting.
  
  
  # Step 1: Capture the call and resolve multiple-choice arguments --------------
  call <- match.call()
  
  criterion <- match.arg(criterion)
  xorder <- match.arg(xorder)
  ties <- match.arg(ties)
  
  # acdx is not a top-level argument in the formula interface: ACD handling is
  # requested per-variable via fp(..., acdx = TRUE), so reject the matrix-interface
  # spelling (acd_vars) with a clear pointer to the correct syntax.
  dots <- list(...)
  
  if ("acd_vars" %in% names(dots)) {
    stop(
      "`acd_vars` is not supported as an argument to `mfp2.formula()`. ",
      "Use `fp(..., acdx = TRUE)` or `fp2(..., acdx = TRUE)` inside the formula ",
      "to request ACD handling for formula terms.",
      call. = FALSE
    )
  }
  
  # Resolve the family argument to its canonical string form (e.g. "gaussian",
  # "cox") used for family-specific branching throughout this function.
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
  
  # Step 2: Validate that `data` and `formula` were supplied and well-formed ----
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
  
  # Keep the user's formula unchanged for fitted-object metadata.
  formula_user <- formula
  
  # Use this only for internal formula parsing/evaluation. This helps
  # because if a user passes survival::strata() in the formula it can
  # fail similarly stats:offset() as suggested by Terry 
  formula_internal <- normalize_formula_special_namespaces(formula_user)
  
  # Step 3: Validate scalar formula-interface defaults ---------------------------
  # In the formula interface, df, alpha, select, shift, scale, and center are
  # global scalar defaults applied to every predictor unless overridden inside
  # fp(). validate_formula_scalar()/validate_formula_probability() (defined in
  # validation_helpers.R) are thin wrappers around the shared scalar/vector
  # validators (nvars = 1L, since every formula-interface default is a single
  # global value); they only add the "use fp(...)" hint that is specific to
  # this interface.
  
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
  
  # strictly_positive = TRUE folds the old separate "scale must be positive"
  # check into the shared validator, since NA (automatic scaling) is still
  # allowed via allow_na = TRUE.
  validate_formula_scalar(
    scale,
    "scale",
    "numeric",
    allow_null = TRUE,
    allow_na = TRUE,
    strictly_positive = TRUE,
    hint = "i Use fp(..., scale = value) to set different scaling factors for individual variables."
  )
  
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
  
  # Step 4: Validate `weights`, `offset`, and `subset` against `data` -----------
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
  
  # Step 5: Parse the formula and validate `strata` against it -------------------
  # "strata" is registered as a special so that strata() terms inside the
  # formula can be detected and later separated from the predictor terms.
  terms_formula <- stats::terms(
    formula_internal,
    specials = "strata",
    data = data
  )
  
  has_formula_strata <- !is.null(attr(terms_formula, "specials")$strata)
  
  # If strata() is not used inside the formula, the `strata` argument (if any)
  # is validated here against the raw data; it is otherwise extracted from the
  # formula further below.
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
  
  # Step 6: Build the model frame and check for missing predictors/values -------
  # model.frame preserves the attributes of the data unlike model.matrix.
  # Use na.pass here so that mfp2.formula() does not silently drop rows through
  # the global/default na.action option. Missingness is checked explicitly below.
  
  mf <- stats::model.frame(
    terms_formula,
    data = data,
    drop.unused.levels = TRUE,
    na.action = stats::na.pass
  )
  
  
  # Factor columns are expanded by model.matrix() below and then mapped back
  # to their original formula term through the `assign` attribute. Both
  # unordered treatment-contrast blocks and ordered polynomial-contrast blocks
  # are therefore selected jointly by the grouped-term fitting path.
  
  # Reject an intercept-only formula (e.g. y ~ 1): mfp2() requires at least one
  # predictor to perform variable/FP selection on.
  labels <- attr(terms_formula, "term.labels")
  if (length(labels) == 0) {
    stop(
      "! No predictors are provided for model fitting.\n",
      "i At least one predictor is required.",
      call. = FALSE
    )
  }
  
  # Because na.action = na.pass was used above, missing values survive into mf
  # and must be checked for explicitly here, with informative row numbers.
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
  
  
  # Step 7: Extract strata() and offset() terms out of the formula --------------
  # strata() and offset() are model-fitting constructs, not ordinary predictors,
  # so their term positions are recorded here and removed before model.matrix()
  # is used to build the predictor design matrix `x` further below.
  terms_drop <- integer(0L)
  
  # Prediction metadata for formula-level fitting constructs. These are not
  # ordinary predictors and are removed from the model matrix, so they must be
  # stored separately for predict.mfp2(newdata = ...).
  formula_strata_terms <- NULL
  formula_strata_xlevels <- NULL
  formula_offset_terms <- NULL
  formula_offset_xlevels <- NULL
  
  # Stratification for Cox models: strata() is only meaningful for family = "cox".
  if (!is.null(attr(terms_formula, "specials")$strata)) {
    if (family_string == "cox") {
      
      # The formula and the `strata` argument are alternative ways to specify
      # strata; if both are given, the formula's strata() term wins.
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
      
      # for predict
      strata_formula <- stats::reformulate(stemp$vars)
      environment(strata_formula) <- environment(formula_internal)
      
      formula_strata_terms <- stats::terms(
        strata_formula,
        data = data
      )
      
      formula_strata_xlevels <- .getXlevels(formula_strata_terms, mf)
      
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
      
      if (length(stemp$vars) == 1L) {
        strata <- mf[[stemp$vars]]
      } else {
        strata <- do.call(
          survival::strata,
          c(as.list(mf[, stemp$vars, drop = FALSE]), list(shortlabel = TRUE))
        )
      }
      
      attr(strata, "mfp2_strata_keep") <- TRUE
      
      terms_drop <- c(terms_drop, stemp$terms)
    } else {
      stop("! strata are only allowed for Cox models.\n", 
           "i Please remove any strata terms from the model formula.",
           call. = FALSE)
    }
  }
  
  # Offset: an offset() term in the formula takes precedence over the `offset`
  # argument (with a warning), mirroring the strata precedence rule above.
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
    
    # for predict.mfp2
    offset_call <- attr(terms_formula, "variables")[[term_offset + 1L]]
    
    offset_formula <- stats::as.formula(
      call("~", offset_call),
      env = environment(formula_internal)
    )
    
    formula_offset_terms <- stats::terms(
      offset_formula,
      data = data
    )
    
    formula_offset_xlevels <- .getXlevels(formula_offset_terms, mf)
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
  
  # Step 8: Build the model matrix and drop the intercept column ----------------
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
  
  # Extract the response and build the full design matrix (with intercept).
  y <- model.extract(mf, "response")
  
  mm <- stats::model.matrix(terms_model, mf)
  
  # Identify unordered and ordered factor terms before the intercept is
  # removed. Their model-matrix columns will be grouped and forced to fixed
  # linear settings below. Ordered factors therefore retain their configured
  # contrasts (contr.poly by default) but are selected as one conceptual term.
  factor_term_labels <- identify_formula_factor_terms(terms_model, mf)
  conceptual_term_names <- formula_conceptual_term_names(terms_model, mf)
  factor_terms <- unname(conceptual_term_names[factor_term_labels])
  
  assign <- attr(mm, "assign")
  term_labels <- attr(terms_model, "term.labels")
  
  # mfp2 fits models without an intercept term internally (the intercept is
  # handled by fit_mfp()/glm()/coxph()); drop the (Intercept) column here but
  # keep `assign` aligned with the remaining columns of x.
  keep_cols <- colnames(mm) != "(Intercept)"
  x <- mm[, keep_cols, drop = FALSE]
  assign <- assign[keep_cols]
  
  nx <- ncol(x) 
  names_x <- colnames(x)
  
  
  # Map conceptual terms to the exact model-matrix columns they generated.
  # Simple factor wrappers use the source variable name as the conceptual term;
  # their low-level design columns retain model.matrix() names such as
  # factor(x9)2 and factor(x9)3.
  term_to_columns <- build_formula_term_to_columns(
    x_columns = colnames(x),
    assign = assign,
    term_labels = term_labels,
    conceptual_names = conceptual_term_names
  )
  
  # Step 9: Extract fp()/fp2() term settings and rename fp() columns -----------
  # Locate which columns of the model frame come from fp() or fp2() terms, so
  # their per-variable settings (df, scale, shift, powers, acdx, zero, ...) can
  # be pulled from the attributes fp() attached to them.
  fp_pos <- which(is_fp_term(colnames(mf)))
  
  if (length(fp_pos) > 0) {
    fp_data <- mf[, fp_pos, drop = FALSE]
    
    # The real variable name (e.g. "age") is stored as an attribute on the
    # fp()-transformed column (whose model-frame name is literally "fp(age)").
    fp_vars <- unname(sapply(fp_data, function(v) attr(v, "name")))
    
    # A variable must not be wrapped in fp() more than once in the same formula.
    fp_vars_duplicates <- fp_vars[duplicated(fp_vars)]
    if (length(fp_vars_duplicates) != 0) {
      stop("! Variables should be used only once in the fp() within the formula.\n", 
           sprintf("i The following variable(s) are duplicated in fp() function: %s.", 
                   paste0(fp_vars_duplicates, collapse = ", ")), 
           call. = FALSE)
    }
    
    # A variable wrapped in fp() must not also appear elsewhere in the formula
    # (e.g. `y ~ fp(age) + age`), since that would create two representations
    # of the same predictor.
    vars_duplicates <- which(colnames(mf) %in% fp_vars)
    if (length(vars_duplicates) != 0) {
      stop("! Variables used in the fp() should not be included in other parts of the formula.\n", 
           sprintf("i This applies to the following variable(s): %s.", 
                   paste0(colnames(mf)[vars_duplicates], collapse = ", ")), 
           call. = FALSE)
    }
    
    # Capture fp() term labels (e.g. "fp(age)") before renaming them in the model matrix.
    fp_terms <- names_x[is_fp_term(names_x)]
    
    # Rename fp()/fp2() model-matrix columns to the underlying variable name,
    # e.g. "fp(age)" becomes "age", so downstream code (and the fitted object)
    # exposes ordinary variable names rather than formula syntax.
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
  
  # Step 10: Build per-variable default option lists from the global scalars ---
  # Every predictor starts out with the global df/scale/shift/center/alpha/select
  # defaults; fp()-term-specific settings (Step 11) then override individual
  # entries of these lists where the user supplied them.
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
  
  
  # Step 11: Override defaults with per-variable fp() term settings -------------
  # Powers supplied inside fp() are collected separately so they can override
  # top-level `powers` after both have been validated against the final model
  # matrix column names and effective df values.
  powerx <- list()
  
  # acdx/zero/catzero/spike default to NULL (i.e. no variables flagged) unless
  # fp() terms request them below.
  acdx <- NULL
  zero <- NULL
  catzero <- NULL
  spike <- NULL
  
  if (length(fp_pos) != 0) {
    # modifyList() replaces only the named entries supplied by fp() terms,
    # leaving the global default in place for every other variable.
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
  
  # Step 11b: Force factor-generated columns to fixed linear settings ---------
  # A factor term is represented by a contrast block rather than by one
  # continuous covariate. Every generated column is therefore an ordinary
  # linear coefficient. Using shift = 0 and scale = 1 also ensures all members
  # of an ordered-factor contrast block have identical grouped-term metadata.
  factor_columns <- unique(unlist(
    term_to_columns[intersect(factor_terms, names(term_to_columns))],
    use.names = FALSE
  ))
  factor_columns <- intersect(factor_columns, names_x)
  
  if (length(factor_columns) > 0L) {
    df_list[factor_columns] <- rep(list(1L), length(factor_columns))
    shift_list[factor_columns] <- rep(list(0), length(factor_columns))
    scale_list[factor_columns] <- rep(list(1), length(factor_columns))
    force_max_fp_list[factor_columns] <- rep(list(FALSE), length(factor_columns))
  }
  
  # Step 12: Validate and build the final candidate FP power list ---------------
  # df_vec is currently unused further down but is kept available here in case
  # future df-dependent power validation needs it.
  df_vec <- unlist(df_list, use.names = TRUE)
  
  # fp()-specific powers (powerx) take precedence over the top-level `powers`
  # argument, so any variable present in both is dropped from `powers_arg`
  # before validation to avoid a duplicate/conflicting specification.
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
  
  # Step 13: Expand `keep` to full model-matrix column names --------------------
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
  
  # Step 14: Delegate to mfp2.default() and attach formula-specific metadata ---
  # Store the raw model-matrix column names before fp() columns are renamed.
  # These are needed by predict.mfp2() to rebuild the same expanded design
  # matrix from ordinary formula-style newdata.
  formula_column_map <- setNames(names_x, colnames(mm)[keep_cols])
  
  formula_factor_info <- build_formula_factor_info(
    factor_terms = factor_term_labels,
    terms_object = terms_model,
    model_frame = mf,
    term_name_map = conceptual_term_names,
    term_to_columns = term_to_columns,
    x = x
  )
  
  # Pass every factor mapping, including two-level factors that generate one
  # non-identity dummy column. Ordinary identity singleton terms are added by
  # normalize_term_groups() and need not be supplied explicitly.
  grouped_formula_terms <- term_to_columns[
    intersect(factor_terms, names(term_to_columns))
  ]
  if (length(grouped_formula_terms) > 0L) {
    grouped_formula_terms <- grouped_formula_terms[
      mapped_term_flags(grouped_formula_terms)
    ]
  }
  if (length(grouped_formula_terms) == 0L) {
    grouped_formula_terms <- NULL
  }
  
  fit <- mfp2.default(x = x, 
                      y = y,
                      term_groups = grouped_formula_terms,
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
  fit$formula <- formula_user
  fit$formula_terms <- stats::delete.response(terms_model)
  fit$formula_contrasts <- attr(mm, "contrasts")
  fit$formula_xlevels <- .getXlevels(terms_model, mf)
  fit$formula_column_map <- formula_column_map
  fit$formula_model_matrix_columns <- names_x
  fit$formula_term_to_columns <- term_to_columns
  fit$formula_factor_info <- formula_factor_info
  
  fit$formula_strata_terms <- formula_strata_terms
  fit$formula_strata_xlevels <- formula_strata_xlevels
  fit$formula_offset_terms <- formula_offset_terms
  fit$formula_offset_xlevels <- formula_offset_xlevels
  # Replace the internal mfp2.default() call with the original
  # user-facing formula-interface call.
  call[[1L]] <- quote(mfp2)
  
  formula_position <- match("formula", names(call))
  
  if (!is.na(formula_position)) {
    names(call)[formula_position] <- ""
  }
  
  fit$call_mfp <- call
  
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
  # Recognize both ordinary and namespace-qualified FP terms:
  # fp(x), fp2(x), mfp2::fp(x), and mfp2::fp2(x).
  grepl("^(mfp2::)?fp2?\\(.*\\)$", z)
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
  # Ensure that the supplied object was created by mfp2().
  if (!inherits(object, "mfp2")) {
    stop("The object is not an mfp2 object.", call. = FALSE)
  }
  
  # Dispatch to the summary method for the underlying fitted model:
  # summary.glm() for GLMs or summary.coxph() for Cox models.
  result <- NextMethod("summary")
  
  # Replace the internal glm()/coxph() call with the original
  # user-facing mfp2() call.
  if (!is.null(object$call_mfp)) {
    result$call <- object$call_mfp
  }
  
  result
}

#' Print method for objects of class `mfp2`
#'
#' Prints a structured summary of an \code{mfp2} model, including the original
#' model call, selection criterion, convergence status, selected and excluded
#' variables, covariate preprocessing, a readable breakdown of the final
#' functional form chosen for every variable (split into standard MFP, ACD,
#' and spike-at-zero (SAZ) results), the complete raw settings table, final
#' coefficients, and model-fit measures.
#'
#' @details
#' The "Summary of Function Selection" section reports one row per variable,
#' split across three sub-tables so that ACD- and SAZ-specific detail never
#' has to be crammed into the same columns as ordinary variables:
#' \itemize{
#'   \item \strong{Standard MFP}: variables that used neither ACD nor SAZ.
#'   \item \strong{Approximate Cumulative Distribution (ACD), non-spike
#'     variables}: ACD variables not assessed by the SAZ algorithm.
#'   \item \strong{Spike-at-Zero (SAZ)}: all spike-eligible variables,
#'     including ones that were also ACD-transformed (marked via their own
#'     \code{ACD} column), since the SAZ decision is the more consequential
#'     fact for those variables.
#' }
#' Within every table, selected variables are listed before excluded ones.
#'
#' Each table's \code{Function} column is a single human-readable label
#' (e.g. \code{"linear"}, \code{"FP(1, 2)"}, \code{"FP(0.5) (x > 0)"},
#' \code{"FP(-1, -1) (x > 0) + binary"}, \code{"binary indicator only"}, or
#' \code{"out"}), built from the selected powers together with the
#' \code{zero}/\code{catzero} status.
#'
#' The subsequent "Detailed Settings" table reproduces the original,
#' unabridged \code{fp_terms} columns (renamed only for readability, e.g.
#' \code{catzero} -> \code{catzero_final} and \code{spike_dec} ->
#' \code{saz_decision}), also sorted with selected variables first. Notes
#' about \code{select}/\code{alpha} divergence and the
#' \code{catzero}-implies-\code{zero} relationship, plus definitions for any
#' column whose meaning may not be obvious, are printed immediately below it.
#'
#' The reported model-fit values use family-specific definitions.
#'
#' For generalized linear models:
#' \itemize{
#'   \item the null-model value is the null deviance returned by
#'     \code{glm.fit()} or \code{glm()};
#'   \item the full-linear value is the residual deviance of the model containing
#'     all candidate predictors as ordinary linear terms; and
#'   \item the final-MFP value is the residual deviance of the selected MFP
#'     model.
#' }
#'
#' For Cox proportional-hazards models, the reported values are minus twice the
#' corresponding null or fitted partial log-likelihood.
#'
#' The inherited \code{print.glm()} or \code{print.coxph()} method is not used.
#' Those methods print family-specific footer information such as residual
#' degrees of freedom and AIC, but they do not report the full-linear reference
#' model used by the MFP procedure.
#'
#' The selection criterion is read from \code{x$criterion_mfp} when available.
#' For older objects without that field, it is inferred from the
#' \code{select} and \code{alpha} columns of \code{x$fp_terms}.
#'
#' Older serialized \code{mfp2} objects may not contain all three model-fit
#' fields. Missing values are printed as \code{NA} rather than causing the
#' print method to fail.
#'
#' @param x An object of class \code{"mfp2"}.
#' @param detailed_settings Logical. If \code{FALSE}, omits the "Detailed
#'   Settings" table (the complete, unabridged \code{fp_terms} columns)
#'   entirely, along with its associated \code{Note:} lines and column
#'   definitions (which are controlled independently by \code{notes}, but
#'   have nothing left to annotate when this table isn't printed). Default is
#'   \code{TRUE}.
#' @param notes Logical. If \code{FALSE}, suppresses all auto-generated
#'   explanatory text printed after the Detailed Settings table -- both the
#'   \code{Note:} lines (e.g. the \code{catzero}/\code{zero} relationship,
#'   \code{select}/\code{alpha} divergence) and the column-definition
#'   paragraphs (\code{acd}, \code{zero}, \code{spike}, \code{saz_decision}).
#'   Has no effect when \code{detailed_settings = FALSE}. Default is
#'   \code{TRUE}.
#' @param ... Further arguments passed to \code{print.default()} when printing
#'   the coefficient vector. The \code{digits} argument, when supplied, is also
#'   used to format the model-fit values.
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.mfp2 <- function(x, detailed_settings = TRUE, notes = TRUE, ...) {
  
  # ---------------------------------------------------------------------------
  # Step 1: Resolve display settings
  # ---------------------------------------------------------------------------
  
  dots <- list(...)
  digits <- dots$digits
  
  if (is.null(digits)) {
    digits <- max(3L, getOption("digits") - 3L)
  }
  
  output_width <- 78L
  
  make_boundary <- function(character = "-") {
    paste(rep(character, output_width), collapse = "")
  }
  
  # Full-width heading, used for top-level sections.
  print_section_heading <- function(title) {
    cat(make_boundary("-"), "\n", sep = "")
    cat(title, "\n")
    cat(make_boundary("-"), "\n\n", sep = "")
  }
  
  # Light heading, used for the Standard MFP / ACD / SAZ sub-tables nested
  # inside "Summary of Function Selection": an underline matching the title's
  # own width, so it reads as "one topic, several angles" rather than three
  # unrelated top-level sections.
  print_subsection_heading <- function(title) {
    cat(title, "\n", sep = "")
    cat(paste(rep("-", nchar(title)), collapse = ""), "\n", sep = "")
  }
  
  print_variable_list <- function(label, variables) {
    variable_text <- if (length(variables) > 0L) {
      paste(variables, collapse = ", ")
    } else {
      "none"
    }
    
    complete_text <- paste0(label, ": ", variable_text)
    
    wrapped_text <- strwrap(
      complete_text,
      width = output_width,
      exdent = nchar(label) + 2L
    )
    
    cat(paste(wrapped_text, collapse = "\n"), "\n", sep = "")
  }
  
  format_convergence <- function(value) {
    if (length(value) != 1L || is.na(value)) {
      return("unknown")
    }
    
    if (isTRUE(value)) {
      return("yes")
    }
    
    if (identical(value, FALSE)) {
      return("no")
    }
    
    "unknown"
  }
  
  format_criterion <- function(value) {
    if (length(value) != 1L || is.na(value)) {
      return(NULL)
    }
    
    value <- as.character(value)
    criterion_key <- tolower(gsub("[^[:alnum:]]", "", value))
    
    switch(
      criterion_key,
      "pvalue" = "p-value",
      "aic" = "AIC",
      "bic" = "BIC",
      value
    )
  }
  
  # Build a single human-readable functional-form label from a variable's
  # selected powers together with its zero/catzero status. `acd_prefix`
  # should be FALSE whenever the calling table already has dedicated
  # power-on-x / power-on-A(x) columns (the ACD table), and TRUE only when
  # the table has no such columns but the row is nonetheless ACD-transformed
  # (an ACD variable that also went through the SAZ table).
  format_function_label <- function(powers, zero, catzero, acd_prefix = FALSE,
                                    selected = TRUE) {
    # An unselected variable is always "out", regardless of any lingering
    # catzero/zero flag values -- those describe what was *requested*, not
    # what survived selection, and the two can disagree for eliminated
    # spike variables depending on how spike_decision was left set.
    if (!isTRUE(selected)) {
      return("out")
    }
    
    has_continuous <- length(powers) > 0L
    
    if (has_continuous) {
      base_label <- if (length(powers) == 1L && isTRUE(powers == 1)) {
        "linear"
      } else {
        sprintf("FP(%s)", paste(powers, collapse = ", "))
      }
      
      if (isTRUE(acd_prefix)) {
        base_label <- paste0("ACD ", base_label)
      }
      
      if (isTRUE(zero)) {
        base_label <- paste0(base_label, " (x > 0)")
      }
      
      if (isTRUE(catzero)) {
        base_label <- paste0(base_label, " + binary")
      }
      
      return(base_label)
    }
    
    if (isTRUE(catzero)) {
      return("binary indicator only")
    }
    
    "out"
  }
  
  # Extract a variable's non-NA selected powers from a single fp_terms row,
  # given the names of the power1, power2, ... columns.
  extract_powers <- function(row, power_cols) {
    p <- as.numeric(row[power_cols])
    p[!is.na(p)]
  }
  
  # Sort a data.frame so that selected rows (selected_status = TRUE) come
  # first, preserving original relative order within each group.
  sort_selected_first <- function(df, selected_status) {
    df[order(!selected_status), , drop = FALSE]
  }
  
  
  # ---------------------------------------------------------------------------
  # Step 2: Prepare the detailed MFP table and selection indicators
  # ---------------------------------------------------------------------------
  
  # Work on a print-only copy. Do not modify `x$fp_terms`, because downstream
  # package code may depend on its original column names and numeric SAZ codes.
  fp_terms <- x$fp_terms
  
  variable_names <- rownames(fp_terms)
  
  if (is.null(variable_names) &&
      !is.null(x$fp_powers) &&
      length(x$fp_powers) == nrow(fp_terms)) {
    variable_names <- names(x$fp_powers)
  }
  
  if (is.null(variable_names) ||
      length(variable_names) != nrow(fp_terms)) {
    variable_names <- paste0("V", seq_len(nrow(fp_terms)))
  }
  
  rownames(fp_terms) <- variable_names
  
  if ("selected" %in% names(fp_terms)) {
    selected_status <- fp_terms[["selected"]]
    
    if (!is.logical(selected_status)) {
      selected_status <- tolower(as.character(selected_status)) %in%
        c("true", "t", "yes", "y", "1")
    }
  } else if ("df_final" %in% names(fp_terms)) {
    selected_status <- !is.na(fp_terms[["df_final"]]) &
      fp_terms[["df_final"]] > 0
  } else {
    selected_status <- rep(TRUE, nrow(fp_terms))
  }
  
  if (anyNA(selected_status)) {
    selected_fallback <- if ("df_final" %in% names(fp_terms)) {
      !is.na(fp_terms[["df_final"]]) &
        fp_terms[["df_final"]] > 0
    } else {
      rep(FALSE, nrow(fp_terms))
    }
    
    selected_status[is.na(selected_status)] <-
      selected_fallback[is.na(selected_status)]
  }
  
  acd_flag <- if ("acd" %in% names(fp_terms)) {
    as.logical(fp_terms[["acd"]])
  } else {
    rep(FALSE, nrow(fp_terms))
  }
  
  zero_flag <- if ("zero" %in% names(fp_terms)) {
    as.logical(fp_terms[["zero"]])
  } else {
    rep(FALSE, nrow(fp_terms))
  }
  
  catzero_flag <- if ("catzero" %in% names(fp_terms)) {
    as.logical(fp_terms[["catzero"]])
  } else {
    rep(FALSE, nrow(fp_terms))
  }
  
  spike_flag <- if ("spike" %in% names(fp_terms)) {
    sp <- fp_terms[["spike"]]
    if (!is.logical(sp)) {
      tolower(as.character(sp)) %in% c("true", "t", "yes", "y", "1")
    } else {
      sp
    }
  } else {
    rep(FALSE, nrow(fp_terms))
  }
  
  decision_code <- if ("spike_dec" %in% names(fp_terms)) {
    suppressWarnings(as.integer(fp_terms[["spike_dec"]]))
  } else {
    rep(NA_integer_, nrow(fp_terms))
  }
  
  power_cols <- grep("^power[0-9]+$", names(fp_terms), value = TRUE)
  powers_by_row <- lapply(seq_len(nrow(fp_terms)), function(i) {
    extract_powers(fp_terms[i, , drop = FALSE], power_cols)
  })
  
  df_initial_all <- if ("df_initial" %in% names(fp_terms)) fp_terms[["df_initial"]] else rep(NA, nrow(fp_terms))
  df_final_all   <- if ("df_final" %in% names(fp_terms)) fp_terms[["df_final"]] else rep(NA, nrow(fp_terms))
  df_change_all  <- sprintf("%s -> %s", df_initial_all, df_final_all)
  
  selected_variables <- variable_names[selected_status]
  excluded_variables <- variable_names[!selected_status]
  
  # SAZ decision text, one label per variable (only meaningful for
  # spike-eligible variables; used both in the SAZ table and in
  # Detailed Settings' saz_decision column).
  saz_decision_text <- rep("not SAZ", nrow(fp_terms))
  saz_decision_text[spike_flag & !selected_status] <- "not selected"
  saz_decision_text[spike_flag & selected_status] <- saz_decision_label(
    decision_code[spike_flag & selected_status],
    style = "print",
    unknown = "unknown"
  )
  # saz_decision_label(style = "print") returns "cont + binary" for the
  # combined decision; standardize to the fuller wording used elsewhere in
  # this print method.
  saz_decision_text[saz_decision_text == "cont + binary"] <- "continuous + binary"
  
  
  # ---------------------------------------------------------------------------
  # Step 3: Determine the model-selection criterion
  # ---------------------------------------------------------------------------
  
  criterion_label <- format_criterion(x$criterion_mfp)
  
  if (is.null(criterion_label)) {
    criterion_values <- character(0L)
    
    if ("select" %in% names(fp_terms)) {
      criterion_values <- c(criterion_values, as.character(fp_terms[["select"]]))
    }
    
    if ("alpha" %in% names(fp_terms)) {
      criterion_values <- c(criterion_values, as.character(fp_terms[["alpha"]]))
    }
    
    criterion_values <- unique(toupper(trimws(criterion_values)))
    criterion_values <- criterion_values[!is.na(criterion_values) & nzchar(criterion_values)]
    
    if ("AIC" %in% criterion_values) {
      criterion_label <- "AIC"
    } else if ("BIC" %in% criterion_values) {
      criterion_label <- "BIC"
    } else {
      criterion_label <- "p-value"
    }
  }
  
  is_pvalue_criterion <- identical(criterion_label, "p-value")
  
  
  # ---------------------------------------------------------------------------
  # Step 4: Print the main output banner
  # ---------------------------------------------------------------------------
  
  cat(make_boundary("="), "\n", sep = "")
  cat("MFP Model Fit\n")
  cat(make_boundary("="), "\n\n", sep = "")
  
  
  # ---------------------------------------------------------------------------
  # Step 5: Print the original mfp2 call
  # ---------------------------------------------------------------------------
  
  print_section_heading("Model Call")
  
  if (!is.null(x$call_mfp)) {
    print(x$call_mfp)
  } else {
    cat("(original mfp2 call unavailable)\n")
  }
  
  cat("\n")
  
  
  # ---------------------------------------------------------------------------
  # Step 6: Print the model-selection summary
  # ---------------------------------------------------------------------------
  
  print_section_heading("Selection Summary")
  
  cat("Converged: ", format_convergence(x$convergence_mfp), "\n", sep = "")
  cat("Criterion: ", criterion_label, "\n", sep = "")
  
  # A variable counts as "linear" if its sole continuous term is FP1 with
  # power exactly 1 (and it is not ACD-transformed), or if it is a
  # binary-only spike variable (no continuous component at all, just a 0/1
  # indicator). Every other selected variable is "nonlinear": any FP2
  # (including FP(1, 1), which contains a log(x) term), any FP1 with power
  # != 1, and any ACD variable regardless of its specific powers. Whether a
  # zero/catzero indicator is additionally present does not, by itself,
  # change this classification.
  is_plain_linear <- !acd_flag &
    vapply(powers_by_row, function(p) length(p) == 1L && isTRUE(p == 1), logical(1L))
  is_binary_only <- vapply(powers_by_row, length, integer(1L)) == 0L & catzero_flag
  
  linear_status <- selected_status & (is_plain_linear | is_binary_only)
  nonlinear_status <- selected_status & !linear_status
  
  print_variable_list(
    label = "Selected variables (linear)",
    variables = variable_names[linear_status]
  )
  
  print_variable_list(
    label = "Selected variables (nonlinear)",
    variables = variable_names[nonlinear_status]
  )
  
  print_variable_list(
    label = "Excluded variables",
    variables = excluded_variables
  )
  
  cat("\n")
  
  
  # ---------------------------------------------------------------------------
  # Step 7: Print covariate preprocessing information
  # ---------------------------------------------------------------------------
  
  print_section_heading("Covariate Preprocessing")
  
  if (!is.null(x$transformations)) {
    print.data.frame(x$transformations, right = FALSE)
  } else {
    cat("(preprocessing information unavailable)\n")
  }
  
  cat("\n")
  
  
  # ---------------------------------------------------------------------------
  # Step 8: Summary of Function Selection (Standard MFP / ACD / SAZ)
  # ---------------------------------------------------------------------------
  
  print_section_heading("Summary of Function Selection")
  
  function_label_plain <- mapply(
    FUN = format_function_label,
    powers = powers_by_row,
    zero = zero_flag,
    catzero = catzero_flag,
    selected = selected_status,
    MoreArgs = list(acd_prefix = FALSE),
    SIMPLIFY = TRUE
  )
  
  # Partition: every variable belongs to exactly one of the three tables.
  # Spike-eligible variables go to SAZ regardless of ACD status (SAZ is the
  # more consequential decision for them); among the rest, ACD variables go
  # to the ACD table; everyone else goes to Standard MFP.
  in_saz <- spike_flag
  in_acd_only <- acd_flag & !spike_flag
  in_standard <- !acd_flag & !spike_flag
  
  # --- Standard MFP ----------------------------------------------------------
  
  standard_table <- data.frame(
    Variable = variable_names[in_standard],
    Selected = ifelse(selected_status[in_standard], "yes", "no"),
    `df (init -> final)` = df_change_all[in_standard],
    Function = unname(function_label_plain[in_standard]),
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  standard_table <- sort_selected_first(standard_table, selected_status[in_standard])
  
  print_subsection_heading("Standard MFP")
  print.data.frame(standard_table, row.names = FALSE, right = FALSE)
  cat("\n")
  
  # --- ACD (non-spike) --------------------------------------------------------
  
  if (any(in_acd_only)) {
    power1_all <- if ("power1" %in% names(fp_terms)) fp_terms[["power1"]] else rep(NA, nrow(fp_terms))
    power2_all <- if ("power2" %in% names(fp_terms)) fp_terms[["power2"]] else rep(NA, nrow(fp_terms))
    
    format_na_dot <- function(v) ifelse(is.na(v), ".", format(v, trim = TRUE))
    
    acd_table <- data.frame(
      Variable = variable_names[in_acd_only],
      Selected = ifelse(selected_status[in_acd_only], "yes", "no"),
      `df (init -> final)` = df_change_all[in_acd_only],
      `Power on x` = format_na_dot(power1_all[in_acd_only]),
      `Power on A(x)` = format_na_dot(power2_all[in_acd_only]),
      Function = unname(function_label_plain[in_acd_only]),
      row.names = NULL,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    acd_table <- sort_selected_first(acd_table, selected_status[in_acd_only])
    
    print_subsection_heading(
      "Approximate Cumulative Distribution (ACD) -- non-spike variables"
    )
    print.data.frame(acd_table, row.names = FALSE, right = FALSE)
    cat("\n")
  }
  
  # --- Spike-at-Zero (SAZ) -----------------------------------------------------
  
  if (any(in_saz)) {
    function_label_saz <- mapply(
      FUN = format_function_label,
      powers = powers_by_row[in_saz],
      zero = zero_flag[in_saz],
      catzero = catzero_flag[in_saz],
      acd_prefix = acd_flag[in_saz],
      selected = selected_status[in_saz],
      SIMPLIFY = TRUE
    )
    
    saz_table <- data.frame(
      Variable = variable_names[in_saz],
      Selected = ifelse(selected_status[in_saz], "yes", "no"),
      `df (init -> final)` = df_change_all[in_saz],
      ACD = ifelse(acd_flag[in_saz], "yes", "no"),
      Decision = saz_decision_text[in_saz],
      Function = unname(function_label_saz),
      row.names = NULL,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    saz_table <- sort_selected_first(saz_table, selected_status[in_saz])
    
    print_subsection_heading(
      "Spike-at-Zero (SAZ) -- may include ACD variables"
    )
    print.data.frame(saz_table, row.names = FALSE, right = FALSE)
    cat("\n")
  }
  
  
  # ---------------------------------------------------------------------------
  # Step 9: Detailed Settings (complete raw fp_terms table)
  # ---------------------------------------------------------------------------
  
  if (isTRUE(detailed_settings)) {
    
    print_section_heading("Detailed Settings")
    
    detailed_table <- data.frame(
      Selected = ifelse(selected_status, "yes", "no"),
      `df (init -> final)` = df_change_all,
      row.names = variable_names,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    
    if (is_pvalue_criterion) {
      if ("select" %in% names(fp_terms)) detailed_table[["select"]] <- fp_terms[["select"]]
      if ("alpha" %in% names(fp_terms)) detailed_table[["alpha"]] <- fp_terms[["alpha"]]
    }
    detailed_table[["acd"]] <- ifelse(acd_flag, "yes", "no")
    detailed_table[["zero"]] <- ifelse(zero_flag, "yes", "no")
    detailed_table[["catzero_final"]] <- ifelse(catzero_flag, "yes", "no")
    detailed_table[["spike"]] <- ifelse(spike_flag, "yes", "no")
    detailed_table[["saz_decision"]] <- saz_decision_text
    if ("power1" %in% names(fp_terms)) detailed_table[["power1"]] <- fp_terms[["power1"]]
    if ("power2" %in% names(fp_terms)) detailed_table[["power2"]] <- fp_terms[["power2"]]
    
    detailed_table <- sort_selected_first(detailed_table, selected_status)
    
    print.data.frame(detailed_table, right = FALSE)
    cat("\n")
    
    # --- Notes and column definitions ---------------------------------------
    #
    # Everything here is auxiliary explanatory text, not part of the data
    # itself, so it is all controlled by the `notes` argument together.
    
    if (isTRUE(notes)) {
      
      note_lines <- character(0L)
      
      if (any(catzero_flag)) {
        note_lines <- c(
          note_lines,
          "catzero implies zero.",
          paste0(
            "catzero_final = yes means a binary indicator for the ",
            "non-positive component is included in the final model."
          )
        )
      }
      
      if (is_pvalue_criterion && all(c("select", "alpha") %in% names(fp_terms))) {
        select_num <- suppressWarnings(as.numeric(fp_terms[["select"]]))
        alpha_num  <- suppressWarnings(as.numeric(fp_terms[["alpha"]]))
        
        mode_value <- function(v) {
          v <- v[!is.na(v)]
          if (length(v) == 0L) return(NA_real_)
          tab <- table(v)
          as.numeric(names(tab)[which.max(tab)])
        }
        
        select_common <- mode_value(select_num)
        alpha_common  <- mode_value(alpha_num)
        
        diverges <- (!is.na(select_num) & select_num != select_common) |
          (!is.na(alpha_num) & alpha_num != alpha_common)
        
        if (any(diverges)) {
          diverging_vars <- variable_names[diverges]
          note_lines <- c(note_lines, sprintf(
            "%s use%s select = %s, alpha = %s (all others use select = %s, alpha = %s).",
            paste(diverging_vars, collapse = ", "),
            if (length(diverging_vars) == 1L) "s" else "",
            paste(unique(format(select_num[diverges], trim = TRUE)), collapse = "/"),
            paste(unique(format(alpha_num[diverges], trim = TRUE)), collapse = "/"),
            format(select_common, trim = TRUE),
            format(alpha_common, trim = TRUE)
          ))
        }
      }
      
      for (note in note_lines) {
        cat("Note: ", note, "\n", sep = "")
      }
      
      if (length(note_lines) > 0L) cat("\n")
      
      # --- Column definitions ---------------------------------------------
      # catzero_final's meaning is covered in the Notes above (together with
      # the catzero-implies-zero relationship), not repeated here.
      
      if (any(acd_flag)) {
        cat(
          "acd:\n",
          "  yes means the variable underwent an approximate cumulative ",
          "distribution\n  (ACD) transformation. acd and spike are independent ",
          "settings, so a\n  variable can be yes for both; ACD variables that ",
          "were also assessed by\n  the spike-at-zero algorithm are listed in ",
          "the SAZ table above (with an\n  ACD column marking them), not in the ",
          "ACD table.\n\n",
          sep = ""
        )
      }
      
      if (any(zero_flag)) {
        cat(
          "zero:\n",
          "  yes means FP transformations are applied only to the positive ",
          "component\n  of the variable.\n\n",
          sep = ""
        )
      }
      
      if (any(spike_flag)) {
        cat(
          "spike:\n",
          "  yes means SAZ modelling was requested and remained eligible after ",
          "the\n  eligibility checks.\n\n",
          sep = ""
        )
        
        cat(
          "saz_decision:\n",
          "  Final SAZ status for each variable.\n\n",
          sep = ""
        )
      }
      
    }
    
  }
  
  
  # ---------------------------------------------------------------------------
  # Step 10: Print final-model coefficients
  # ---------------------------------------------------------------------------
  
  print_section_heading("Final Model Coefficients")
  
  coefficients <- stats::coef(x)
  
  if (length(coefficients) > 0L) {
    print.default(coefficients, ...)
  } else {
    cat("(none)\n")
  }
  
  cat("\n")
  
  
  # ---------------------------------------------------------------------------
  # Step 11: Collect the three stored model-fit values
  # ---------------------------------------------------------------------------
  
  get_deviance_value <- function(value) {
    if (length(value) != 1L) {
      return(NA_real_)
    }
    
    suppressWarnings(as.numeric(value))
  }
  
  deviance_values <- c(
    "Null model" = get_deviance_value(x$null_deviance),
    "Full linear model" = get_deviance_value(x$linear_deviance),
    "Final MFP model" = get_deviance_value(x$mfp_deviance)
  )
  
  
  # ---------------------------------------------------------------------------
  # Step 12: Format and print the model-fit values
  # ---------------------------------------------------------------------------
  
  formatted_deviance <- vapply(
    deviance_values,
    function(value) {
      if (is.na(value)) {
        return("NA")
      }
      
      if (!is.finite(value)) {
        return(as.character(value))
      }
      
      format(value, digits = digits, trim = TRUE)
    },
    character(1L)
  )
  
  deviance_labels <- paste0(trimws(names(formatted_deviance)), ":")
  label_width <- max(nchar(deviance_labels))
  
  print_section_heading("Model Deviances")
  
  for (i in seq_along(formatted_deviance)) {
    cat(
      sprintf(
        "%-*s  %s\n",
        label_width,
        deviance_labels[[i]],
        unname(formatted_deviance[[i]])
      )
    )
  }
  
  cat("\n")
  cat(make_boundary("="), "\n", sep = "")
  
  invisible(x)
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
  #
  # validate_scalar_fp() (defined in validation_helpers.R) is a thin wrapper
  # around the shared scalar/vector validators (nvars = 1L, since every fp()
  # term attribute is a single value); it attaches a hint identifying which
  # fp() call and which argument raised the error.
  
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
#' @keywords internal
#' @noRd
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
