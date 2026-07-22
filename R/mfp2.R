#' Multivariable Fractional Polynomial Models with Extensions
#'
#' Fits multivariable fractional polynomial (MFP) models with simultaneous selection
#' of variables and functional forms for continuous predictors. In addition to 
#' standard MFP modelling, `mfp2()` supports approximate cumulative 
#' distribution (ACD) functions for sigmoid-shaped associations and
#' spike-at-zero (SAZ) modelling for semi-continuous predictors. Models may be
#' specified with a formula and data frame or with a numeric predictor matrix
#' and response. Both interfaces provide the same core modelling capabilities 
#' but differ in input handling as described below. Supported families are 
#' Gaussian, binomial, and Poisson. Cox proportional hazards models are also 
#' supported, currently limited to right-censored data specified as 
#' `survival::Surv(time, event)`.
#'
#' @section Fractional-polynomial model selection:
#' Fractional polynomials represent continuous predictor effects using powers
#' selected from a predefined set. An FP1 function has one transformed term,
#'
#' \deqn{\beta_1 x^{p_1},}
#'
#' and an FP2 function has two transformed terms,
#'
#' \deqn{\beta_1 x^{p_1} + \beta_2 x^{p_2}.}
#'
#' For repeated powers, \eqn{p_1 = p_2 = p}, the second term is
#' \eqn{x^p \log(x)}. A power of zero denotes \eqn{\log(x)}. The default
#' candidate set is `c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)`.
#'
#' With `criterion = "pvalue"`, variable and functional-form selection use the
#' MFP closed-testing procedure controlled by `select` and `alpha`. With
#' `criterion = "aic"` or `"bic"`, candidate models are compared using the
#' corresponding information criterion. Selection is repeated until the model
#' stabilises or the maximum number of `cycles` is reached.
#'
#' A predictor may be omitted, retained as linear, or represented by an FP1 or
#' FP2 function, subject to its `df`, candidate `powers`, and selection settings.
#' Variables named in `keep` remain in the model, although their functional form
#' may still be selected.
#'
#' See `vignette("MFP_Introduction", package = "mfp2")` for the complete
#' closed-test procedure, the multivariable backfitting algorithm, influential
#' observations, sample-size considerations, and other methodological details.
#'
#' @section Formula and matrix interfaces:
#' The formula interface, `mfp2(formula, data, ...)`, is recommended for most
#' users. Variable-specific settings are supplied with [fp()] or [fp2()] terms;
#' term-specific settings override corresponding global defaults. Factors are
#' expanded automatically, and formula terms may include `offset()` and, for
#' Cox models, `strata()`. Continuous variables included through `.` remain
#' subject to global MFP settings unless overridden in an `fp()` term.
#'
#' The matrix interface, `mfp2(x, y, ...)`, requires a numeric predictor matrix.
#' Categorical predictors must already be represented by numeric indicator or
#' contrast columns. Use `term_groups` when several columns represent one
#' conceptual predictor so that they are selected, retained, and tested jointly.
#'
#' Both interfaces use the same model-selection procedure, but differ in how
#' factors, offsets, strata, and variable-specific options are supplied.
#' Predictor names supplied through `keep`, `powers`, `acdx`, `zero_vars`,
#' `catzero_vars`, or `spike_vars` are validated against the applicable formula
#' terms, matrix columns, or `term_groups`; unknown names cause an error.
#'
#' @section Model families and responses:
#' Supported families are `"gaussian"`, `"binomial"`, `"poisson"`, and
#' `"cox"`. Gaussian, binomial, and Poisson models accept either a character
#' family name or a corresponding family object, including supported alternative
#' links. See [stats::family()] and [stats::glm()] for general GLM behaviour.
#'
#' Gaussian models use a finite numeric response vector. Poisson models use a
#' finite, nonnegative numeric response vector. Binomial models accept a numeric
#' response in `[0, 1]`, a two-level factor, or a two-column numeric matrix of
#' nonnegative grouped counts with a positive row total. Cox models require
#' `family = "cox"` and an ordinary
#' right-censored response created with [survival::Surv()], using two columns for
#' follow-up time and event status.
#'
#' @section Shifting, scaling, and centering:
#' FP transformations involving logarithms, negative powers, or fractional
#' powers require positive inputs. By default, `mfp2()` estimates a shift where
#' needed and a power-of-ten scale factor to avoid numerically extreme values.
#' Set `shift = 0` or `scale = 1` to disable the corresponding automatic step.
#' The preprocessing order is shift, scale, FP or ACD transformation, and then
#' centering.
#'
#' In the matrix interface, an unnamed scalar is a global setting, whereas a
#' named vector is variable-specific and may name only a subset of columns. For
#' example, `shift = c(age = 20)` fixes the shift for `age` at 20 and leaves all
#' other shifts automatic; `scale = c(age = 10)` behaves analogously. Unspecified
#' columns are represented internally by missing-value sentinels and are passed
#' through the existing automatic estimation steps. Users must not supply `NA`
#' as an explicit shift or scale value.
#'
#' An explicit shift is never silently enlarged. For every predictor ultimately
#' fitted with a nonlinear FP transformation, all shifted values must be strictly
#' positive. If a supplied shift is insufficient, fitting stops and reports the
#' affected variable or variables.
#'
#' Linear terms (`df = 1`) and variables using `zero`, `catzero`, or `spike`
#' handling receive a zero shift so that their zero component remains at zero.
#' With `center = TRUE`, ordinary transformed continuous terms are centered on
#' their fitted-sample mean. For zero-handled FP columns, the centering constant
#' is calculated from nonzero transformed values and structural-zero rows remain
#' zero. Binary terms are centered by subtracting their lower observed value.
#' The final fitted basis is shifted but unscaled, so a coefficient may multiply
#' a term such as `phi(x + shift)` rather than the raw predictor `x`. New
#' observations supplied to [predict.mfp2()] are given on the original predictor
#' scale; the fitted shift, scale, and centering values are reapplied
#' automatically.
#'
#' The exact automatic rules and worked calculations are described in
#' `vignette("mfp2_introduction", package = "mfp2")`.
#'
#' @section Categorical predictors and grouped terms:
#' In the formula interface, unordered and ordered factors are expanded with
#' [stats::model.matrix()]. All columns generated from one factor are treated as
#' one fixed linear term and are selected, retained, or removed jointly. Naming
#' a factor in `keep` retains its complete contrast block.
#'
#' Ordered factors use their configured contrasts, which are polynomial
#' contrasts by default. The resulting test is an omnibus test of the complete
#' ordered-factor effect; it is not a one-degree-of-freedom trend test and does
#' not impose monotonicity. For a prespecified ordinal trend, supply an explicit
#' numeric score and restrict it to a linear effect with `df = 1`.
#'
#' FP, ACD, zero, catzero, and SAZ processing do not apply to factor contrast
#' blocks. For formula fits, prediction uses the original factor variables and
#' fitted levels and contrasts. In the matrix interface, manually created
#' indicator or contrast columns should be grouped with `term_groups`.
#'
#' Detailed examples, including factor behaviour after subsetting, are provided
#' in `vignette("mfp2_introduction", package = "mfp2")`.
#'
#' @section ACD modelling:
#' The ACD extension is intended for continuous predictors whose associations
#' may have a smooth sigmoid shape that is difficult to represent adequately
#' with the standard FP1 and FP2 families. The transformation approximates the
#' empirical cumulative distribution by smoothing inverse-normal ranks and maps
#' the predictor to an approximately uniform scale on `(0, 1)`.
#'
#' In the matrix interface, list candidate variables in `acdx`. In the formula
#' interface, use `fp(x, acdx = TRUE)`. Requesting ACD sets the term's maximum
#' `df` to 4 so that the complete FSPA candidate family can be assessed. An ACD
#' request is reset to standard MFP handling when the fitting data contain fewer
#' than five distinct values for that term. For an eligible ACD term, `mfp2()`
#' evaluates an ACD-extended two-component function and simpler alternatives,
#' including standard FP1, ACD-only, linear, and omitted forms. Selection uses
#' the ACD function-selection procedure when `criterion = "pvalue"`, or direct
#' AIC/BIC comparison for information-criterion selection.
#'
#' See `vignette("mfp2_ACD", package = "mfp2")` for the complete construction
#' of the transformation, the six candidate model families, the five-step
#' closed-test sequence, interpretation, and worked examples.
#'
#' @section Nonpositive values and spike-at-zero modelling:
#' The `zero`, `catzero`, and `spike` options address related but distinct
#' modelling goals for covariates containing nonpositive values:
#'
#' \itemize{
#'   \item `zero` applies the continuous FP function to positive values only;
#'     nonpositive values form the zero component and are recoded to zero before
#'     transformation.
#'   \item `catzero` adds a fixed binary indicator for the nonpositive component
#'     while also fitting the positive-component function.
#'   \item `spike` invokes SAZ selection to determine whether the final model
#'     retains both components, only the positive-component function, or only
#'     the binary zero-component indicator.
#' }
#'
#' In the matrix interface, use `zero_vars`, `catzero_vars`, and `spike_vars`.
#' In the formula interface, use `fp(x, zero = TRUE)`,
#' `fp(x, catzero = TRUE)`, or `fp(x, spike = TRUE)`.
#'
#' `catzero` implies `zero`. While SAZ remains eligible, `spike` implies both
#' `catzero` and `zero`. Binary predictors are not eligible for any of these
#' options, and requests for them are reset. Requests for `zero` or `catzero` on
#' an all-positive predictor are also reset because there is no nonpositive
#' component to represent.
#'
#' SAZ selection has two stages. Stage 1 selects the positive-component
#' functional form while retaining the binary component. Stage 2 selects among
#' the three retained representations. A requested SAZ variable is eligible only
#' when both the nonpositive and positive components contain at least
#' `min_saz_component_prop` of the observations and the variable is not binary.
#' With the default value `0.10`, each component must contain at least 10 percent
#' of observations. If eligibility fails, the spike request is reset; explicit
#' `zero` or `catzero` handling is retained only when requested separately.
#'
#' See `vignette("mfp2_spike", package = "mfp2")` for the complete two-stage
#' algorithm, eligibility equations, information-criterion alternatives,
#' component interpretations, and worked examples.
#'
#' @section Subsetting:
#' When `subset` is used, automatically estimated shift and scale values are
#' obtained from the full input data before the subset is applied. Model
#' selection and final fitting then use only the selected observations.
#' Consequently, `subset` is not equivalent to subsetting the data before
#' calling `mfp2()`.
#'
#' In the formula interface, the subset expression is evaluated once with names
#' resolved first from `data` and then from the formula environment. In the
#' matrix interface, supply a logical vector or unique numeric/integer row
#' positions evaluated in the calling environment. Supplied numeric positions
#' retain their order. `weights`, `offset`, and `strata` must remain aligned with
#' the original rows.
#'
#' The argument is intended for analyses that require common shift and scale
#' values across several analysis subsets. It should not be used to remove
#' missing observations or to construct cross-validation folds; prepare those
#' analysis data explicitly before fitting. If subsetting removes factor levels
#' or reduces the estimable dimension of a grouped matrix term, fitting may stop.
#' Detailed factor, contrast, and grouped-term behaviour is documented in
#' `vignette("mfp2_introduction", package = "mfp2")`.
#'
#' @section Cox models:
#' Cox models currently support ordinary right-censored [survival::Surv()]
#' responses of the form `Surv(time, event)`. Event coding is interpreted by
#' [survival::Surv()]; the usual coding is `0` for censoring and `1` for an event.
#' Stratification may be supplied through `strata` in either interface or with
#' `strata()` or `survival::strata()` in a formula. When both are supplied in a
#' formula fit, the formula term takes precedence. Multiple formula strata terms
#' are combined; multiple variables supplied through `strata` may be given as a
#' matrix or data frame. Formula offsets likewise take precedence
#' over an external `offset` argument. The `ties`, `nocenter`, and Cox control
#' settings follow [survival::coxph()]. See [predict.mfp2()] for the additional
#' data required by Cox prediction types such as survival probabilities.
#'
#' @section Compatibility with the `mfp` package:
#' Both `mfp` and `mfp2` export a function named `fp()`. When both packages are
#' attached, use [fp2()] in `mfp2` formulas to avoid namespace ambiguity:
#'
#' \preformatted{
#' fit <- mfp2(y ~ fp2(x1) + fp2(x2), data = dat)
#' }
#'
#' `fp2()` is an alias for [fp()] and accepts the same arguments.
#'
#' @section Convergence and inference:
#' MFP selection typically stabilises within a small number of cycles.
#' Convergence status is available in `convergence_mfp`. If the procedure does
#' not converge, consider increasing `cycles` and reviewing borderline
#' selection decisions, influential observations, and the chosen `select` and
#' `alpha` values.
#'
#' Standard errors, confidence intervals, and p-values from the final fitted
#' model are conditional on the selected model and do not account for uncertainty
#' introduced by variable selection, functional-form selection, ACD selection,
#' or SAZ selection.
#'
#' @param x For `mfp2.default()`, a numeric predictor matrix with one row per
#'   observation and one column per design variable.
#' @param term_groups For `mfp2.default()`, an optional named list mapping a
#'   conceptual term to one or more columns of `x`, for example
#'   `list(race = c("raceB", "raceC"))`. Mapped columns are selected and tested
#'   jointly and must be fixed linear terms (`df = 1`) without ACD, zero,
#'   catzero, or SAZ handling. Unmapped columns are treated as separate terms.
#' @param y Response for `mfp2.default()`. Gaussian models require a finite
#'   numeric vector, and Poisson models require a finite nonnegative numeric
#'   vector. Binomial models accept a numeric binary vector of 0s and 1s, a 
#'   two-level factor, or a two-column matrix of successes and failures. Cox 
#'   models require a two-column right-censored [survival::Surv()] object.
#' @param formula For `mfp2.formula()`, a model formula. Use [fp()] or [fp2()]
#'   to set variable-specific FP, ACD, zero, catzero, or SAZ options.
#' @param data For `mfp2.formula()`, a data frame containing the variables in
#'   `formula`.
#' @param weights Optional numeric observation weights. In the formula interface,
#'   the vector is supplied directly and is not looked up in `data`.
#' @param offset Optional numeric offset with one value per observation. In the
#'   formula interface, an offset may alternatively be included with `offset()`.
#'   Corresponding values must be supplied when predicting from a model that used
#'   an external offset.
#' @param cycles Maximum number of MFP backfitting cycles. Default `5`.
#' @param scale For `mfp2.default()`, `NULL`, a single unnamed positive numeric
#'   value, or a named positive numeric vector for one or more columns of `x`.
#'   An unnamed scalar is applied to every column. Named values are matched to
#'   `colnames(x)`, so their order does not matter; names must be non-empty,
#'   unique, and known columns. A named partial vector fixes only the specified
#'   scales, while scales for unspecified columns are estimated automatically.
#'   Unnamed multi-value vectors are not accepted. `NULL` estimates every scale
#'   automatically, while `scale = 1` disables scaling globally. Binary
#'   predictors are always assigned scale `1`, even when another value is
#'   supplied, because rescaling a two-level predictor is unnecessary. In
#'   `mfp2.formula()`, the top-level value is a
#'   scalar global default; use `fp()` or `fp2()` for variable-specific values.
#'   Final regression coefficients are expressed on the original data scale.
#' @param shift For `mfp2.default()`, `NULL`, a single unnamed finite numeric
#'   value, or a named finite numeric vector for one or more columns of `x`.
#'   An unnamed scalar is applied to every column. Named values are matched to
#'   `colnames(x)`, so their order does not matter; names must be non-empty,
#'   unique, and known columns. A named partial vector fixes only the specified
#'   shifts, while shifts for unspecified columns are estimated automatically.
#'   Unnamed multi-value vectors are not accepted. `NULL` estimates every shift
#'   automatically, while `shift = 0` disables shifting globally. Binary
#'   predictors are always assigned shift `0`, even when another value is
#'   supplied. Each explicit
#'   shift used for a nonlinear FP term must make that predictor strictly
#'   positive; otherwise fitting fails and identifies the affected variable.
#'   In `mfp2.formula()`, the top-level value is a scalar global default; use
#'   `fp()` or `fp2()` for variable-specific values.
#' @param df Maximum degrees of freedom for each predictor's
#'   fractional-polynomial function. Must be `1` (linear) or a positive even
#'   integer (`2` for FP1, `4` for FP2, etc.). Default `4`.
#'   In `mfp2.default()`, accepts a single value (applied to all predictors) or
#'   a numeric vector of length `ncol(x)`. In `mfp2.formula()`, only a single
#'   global default is accepted; override it for individual terms with
#'   `fp(x, df = value)`.
#'   The effective `df` may be reduced automatically when a predictor has few
#'   distinct values: 2--3 distinct values force `df = 1`; 4--5 distinct values
#'   cap `df` at `min(2, requested)`; 6 or more distinct values use the
#'   requested value unchanged.
#' @param center Logical value controlling centering of final transformed terms.
#'   Default `TRUE`. Ordinary continuous transformed terms are centered on their
#'   mean; zero-handled FP columns use the mean of their nonzero transformed
#'   values while zero rows remain zero. Binary terms are centered by
#'   subtracting their lower observed value.
#' @param subset Optional observations used for model selection and fitting.
#'   Formula calls accept an expression evaluated in `data` and the formula
#'   environment. Matrix calls accept a logical vector or unique numeric/integer
#'   row positions. Automatic shift and scale values are estimated before the
#'   subset is applied. See the Subsetting section.
#' @param family Model family. Either a character string (`"gaussian"`,
#'   `"binomial"`, `"poisson"`, or `"cox"`) or a [stats::family()] object such
#'   as `binomial(link = "logit")`. Cox models accept only the character string
#'   `"cox"`. Default `"gaussian"`.
#' @param criterion Selection criterion. One of `"pvalue"` (default), `"aic"`,
#'   or `"bic"`. With `"pvalue"`, variable inclusion and functional-form
#'   complexity are controlled by `select` and `alpha` through nested closed
#'   tests. With `"aic"` or `"bic"`, candidate models are compared directly by
#'   their information-criterion value and `select`/`alpha` are ignored.
#' @param select Nominal significance level for variable selection. Default
#'  `0.05`. In `mfp2.default()`, accepts a single value
#'   (applied to all predictors) or a numeric vector of length `ncol(x)`. In
#'   `mfp2.formula()`, only a single global default is accepted; override it
#'   for individual terms with `fp(x, select = value)`. Setting `select = 1`
#'   for a predictor forces it into the model (equivalent to naming it in
#'   `keep`). Ignored when `criterion` is `"aic"` or `"bic"`.
#' @param alpha Significance level for closed tests between FP functions of
#'   different degrees (e.g. FP2 versus FP1 versus linear). Default `0.05`.
#'   In `mfp2.default()`, accepts a single value (applied to all predictors)
#'   or a numeric vector of length `ncol(x)`. In `mfp2.formula()`, only a
#'   single global default is accepted; override it for individual terms with
#'   `fp(x, alpha = value)`. Ignored when `criterion` is `"aic"` or
#'   `"bic"`.
#' @param keep Character vector of predictor names that must remain in the final
#' model regardless of the selection criterion. Under `criterion = "pvalue"`,
#' this is equivalent to setting `select = 1` for each named variable. The
#' functional form of kept variables may still be simplified by the
#' `alpha`-level closed test or information criterion. Names must match
#' formula term names or, for matrix fits, column names or `term_groups`
#' names.
#' @param xorder Order in which predictors enter the backfitting algorithm.
#' One of `"ascending"` (default), `"descending"`, or `"original"`.
#' `"ascending"` processes predictors from smallest to largest p-value in an
#' initial linear model (most significant first). `"descending"` reverses that
#' order. `"original"` preserves the column order of `x` or the left-to-right
#' order of formula terms.
#' @param powers Optional named list of candidate FP powers for individual predictors.
#' The default set is `c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)`, where 0 stands
#' for the logarithm. Names must match the corresponding formula terms or
#' matrix column names. For predictors with `df > 1`, the candidate set must
#' include at least one power other than 1; the algorithm needs a non-linear 
#' candidate to search over. To restrict a predictor to a linear effect, 
#' set `df = 1` rather than limiting its power set. Per-term settings in
#' [fp()] take precedence over entries in this list.
#' @param ties Method for handling tied event times in Cox models. One of
#' `"breslow"` (default), `"efron"`, or `"exact"`. All methods are equivalent
#' when no ties exist. Ignored for non-Cox families. See
#' [survival::coxph()] for details.
#' @param strata Optional stratification variable(s) for Cox models. Accepts a
#' vector, factor, or a matrix/data frame of multiple stratification
#' variables. In formula fits, `strata()` terms can be used instead and take
#' precedence if both are supplied.
#' @param nocenter Optional numeric vector of predictor indices whose columns
#' should not be internally centered by [survival::coxph()]. Applies to Cox
#' models only. See [survival::coxph()] for details.
#' @param acdx Character vector naming continuous predictors to be assessed with
#' the approximate cumulative distribution (ACD) extension. An ACD request
#' triggers the function-selection procedure for ACD (FSPA), which evaluates
#' extended two-component functions alongside simpler alternatives. See the
#' ACD modelling section. The ACD-transformed version of predictor `x` is
#' named `A_x` in the model output. This argument applies to
#' `mfp2.default()` only; in the formula interface use
#' `fp(x, acdx = TRUE)`.
#' @param ftest Logical. If `TRUE`, use F-distribution critical values instead
#' of chi-squared critical values for the selection tests in Gaussian models.
#' This may improve small-sample behaviour for variable selection,
#' functional-form selection, and spike-at-zero testing. Default `FALSE`.
#' Ignored for non-Gaussian families.
#' @param control Fitting controls from [stats::glm.control()] or
#' [survival::coxph.control()]. `NULL` uses the relevant defaults.
#' @param zero_vars Character vector of continuous predictors whose nonpositive
#' values should be recoded to zero before FP transformation. Only positive
#' values undergo the FP function; nonpositive values contribute zero to the
#' linear predictor. The shift for these variables is forced to zero. Applies
#' to `mfp2.default()` only; in the formula interface use
#' `fp(x, zero = TRUE)`. See the section on nonpositive values and
#' spike-at-zero modelling.
#' @param catzero_vars Character vector of continuous predictors that combine an
#' FP transformation of their positive values with a binary indicator for
#' nonpositive values. The indicator enters the model as a separate fixed
#' covariate. Requesting `catzero` implies `zero` handling for the same
#' variable. Applies to `mfp2.default()` only; in the formula interface use
#' `fp(x, catzero = TRUE)`.
#' @param spike_vars Character vector of continuous predictors to be assessed
#' with the two-stage spike-at-zero (SAZ) algorithm. SAZ selection determines
#' whether the final model retains both the binary zero-indicator and
#' continuous FP component, only the continuous component, or only the binary
#' indicator. Requesting `spike` implies both `catzero` and `zero` handling
#' while SAZ remains eligible. A predictor is eligible for SAZ only when both
#' the nonpositive and positive components meet the `min_saz_component_prop`
#' threshold and the variable is not binary; ineligible spike requests are
#' reset. Applies to `mfp2.default()` only; in the formula interface use
#' `fp(x, spike = TRUE)`.
#' @param min_saz_component_prop Numeric in `(0, 0.5)`. Minimum required
#'   proportion in each component of a spike-at-zero covariate: the zero
#'   component and the positive continuous component. Default `0.10`. A
#'   spike-at-zero candidate is retained for SAZ modelling only if both the
#'   zero proportion and the positive-observation proportion meet this
#'   threshold. Variables that fail this check have their spike flag reset to
#'   `FALSE`. For large samples the default can be lowered, since even a 
#'   small proportion may yield enough observations for reliable estimation.
#' @param force_max_fp_vars For `mfp2.default()`, an optional character
#'   vector naming predictors for which the most complex FP function allowed by
#'   `df` should be forced into the model. For the named predictors,
#'   variable selection and functional-form simplification are bypassed. Under
#'   `criterion = "pvalue"`, this is implemented internally by setting
#'   `select = 1` and `alpha = 1`; the option also applies under
#'   `criterion = "aic"` and `criterion = "bic"`. The default is
#'   `NULL`. For `mfp2.formula()`, set this option for individual
#'   terms using `fp(x, force_max_fp = TRUE)`.
#' @param verbose Logical value indicating whether progress messages are printed.
#' Default `TRUE`.
#' @param \dots Additional arguments passed to methods. Variable-specific formula
#' options should be supplied inside [fp()] or [fp2()].
#' @examples
#'
#' # Gaussian model: formula interface (recommended for most users)
#' data("prostate")
#' fit_formula <- mfp2(
#'   lpsa ~ fp(age) + fp(svi, df = 1) + fp(pgg45) + fp(cavol) +
#'     fp(weight) + fp(bph) + fp(cp),
#'   data = prostate,
#'   verbose = FALSE
#' )
#'
#' fit_formula
#' summary(fit_formula)
#' get_selected_variable_names(fit_formula)
#'
#' plots <- plot(fit_formula)
#' plots[[1]]
#'
#' predict(fit_formula, newdata = prostate[1:5, ])
#'
#' # Matrix interface
#' x <- as.matrix(prostate[, 2:8])
#' y <- prostate$lpsa
#' fit_matrix <- mfp2(x, y, verbose = FALSE)
#' fit_matrix
#'
#' # Named partial settings apply only to the named columns. Here age is fixed,
#' # while shifts and scales for every other matrix column remain automatic.
#' fit_partial <- mfp2(
#'   x,
#'   y,
#'   shift = c(age = 20),
#'   scale = c(age = 10),
#'   verbose = FALSE
#' )
#' fit_partial$transformations[c("age", "weight"), c("shift", "scale")]
#'
#' \donttest{
#' # ACD modelling
#' fit_acd_formula <- mfp2(
#'   lpsa ~ fp(cavol, acdx = TRUE) + fp(age) + svi,
#'   data = prostate,
#'   verbose = FALSE
#' )
#' plot(fit_acd_formula, terms = "cavol")[[1]]
#'
#' # Spike-at-zero modelling
#' set.seed(1)
#' n <- 200
#' exposure <- numeric(n)
#' zero_rows <- sample(seq_len(n), size = round(0.25 * n))
#' exposure[-zero_rows] <- rgamma(
#'   n - length(zero_rows),
#'   shape = 2,
#'   rate = 0.5
#' )
#' age <- runif(n, 20, 80)
#' outcome <- 1.5 * (exposure == 0) +
#'   2 * log(ifelse(exposure > 0, exposure, 1)) +
#'   0.02 * age + rnorm(n)
#'
#' d_spike <- data.frame(outcome, exposure, age)
#'
#' # Formula interface
#' fit_spike_formula <- mfp2(
#'   outcome ~ fp(exposure, spike = TRUE) + fp(age),
#'   data = d_spike,
#'   verbose = FALSE
#' )
#' plot(fit_spike_formula, terms = "exposure")[[1]]
#'
#' # Matrix interface
#' x_spike <- as.matrix(d_spike[, c("exposure", "age")])
#' fit_spike_matrix <- mfp2(
#'   x_spike,
#'   d_spike$outcome,
#'   spike_vars = "exposure",
#'   verbose = FALSE
#' )
#'
#' # Unordered factor: all treatment-contrast columns are selected jointly
#' set.seed(42)
#' d_factor <- data.frame(
#'   y = rnorm(60),
#'   age = runif(60, 20, 80),
#'   group = factor(rep(c("A", "B", "C"), each = 20))
#' )
#' fit_factor <- mfp2(
#'   y ~ age + group,
#'   data = d_factor,
#'   keep = "group",
#'   verbose = FALSE
#' )
#' predict(fit_factor, newdata = d_factor[1:4, ])
#'
#' # Ordered factor: the default polynomial contrasts form one joint term
#' d_factor$severity <- ordered(
#'   rep(c("mild", "moderate", "severe"), length.out = nrow(d_factor)),
#'   levels = c("mild", "moderate", "severe")
#' )
#' fit_ordered <- mfp2(
#'   y ~ age + severity,
#'   data = d_factor,
#'   keep = "severity",
#'   verbose = FALSE
#' )
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
#' # Matrix interface: manually created dummy columns can be grouped as one term
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
#' }
#'
#' @return
#' An object of class `mfp2`. It also inherits from the underlying fitted-model
#' class:
#'
#' \itemize{
#'   \item `glm` for Gaussian, binomial, and Poisson models;
#'   \item `coxph` for Cox proportional hazards models.
#' }
#'
#' Standard fitted-model components, such as coefficients, residuals, fitted
#' values, and the model call, can be accessed using the usual methods for
#' `glm` or `coxph` objects.
#'
#' The following additional components are part of the user-facing `mfp2`
#' result:
#'
#' \describe{
#'   \item{convergence_mfp}{
#'     Logical value indicating whether the MFP selection algorithm converged.
#'   }
#'
#'   \item{fp_terms}{
#'     A data frame with one row per model term. It summarises the initial and
#'     final degrees of freedom, selection settings, selected status,
#'     fractional-polynomial powers, and the effective ACD, zero, catzero, and
#'     spike flags used after validation and eligibility checks.
#'   }
#'
#'   \item{fp_powers}{
#'     A named list containing the selected fractional-polynomial powers for
#'     each model term. Terms excluded from the final model contain missing
#'     powers.
#'   }
#'
#'   \item{transformations}{
#'     A data frame describing the shifting, scaling, and centering applied to
#'     the model terms.
#'   }
#'
#'   \item{acd}{
#'     A named logical vector indicating which terms were assessed with the ACD
#'     extension. The selected powers identify whether the final representation
#'     retained an ACD component.
#'   }
#'
#'   \item{zero}{
#'     A named logical vector indicating which continuous terms were
#'     transformed using only their positive values.
#'   }
#'
#'   \item{catzero}{
#'     A named logical vector indicating which terms include a separate
#'     indicator for nonpositive values.
#'   }
#'
#'   \item{spike_dec}{
#'     A named numeric vector containing the final spike-at-zero decision for
#'     each term. A value of `1` retains both the zero indicator and continuous
#'     component, `2` retains the continuous component only, and `3` retains
#'     the zero indicator only. Non-spike terms use the standard value `2`.
#'   }
#'
#'   \item{call_mfp}{
#'     The call used to fit the model.
#'   }
#' }
#'
#' Use `print()`, `summary()`, `coef()`, `predict()`, and `plot()` for the main
#' model operations. Use [get_selected_variable_names()] to obtain the names of
#' terms retained in the selected model.
#'
#' Other components may be stored to support prediction, plotting, and model
#' reconstruction. Undocumented components should be regarded as internal and
#' may change between package versions.
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
#' [fp()], [fp2()], [summary.mfp2()], [coef.mfp2()], [predict.mfp2()],
#' [plot.mfp2()], [get_selected_variable_names()], [transform_vector_fp()]
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


#' Resolve a Simple fp()/fp2() Formula Term to Its Source Variable
#'
#' @param term_label Original formula term label.
#'
#' @return Source variable name for a simple fp()/fp2() wrapper, otherwise
#'   `NULL`.
#'
#' @keywords internal
#' @noRd
formula_fp_source_name <- function(term_label) {
  expr <- tryCatch(str2lang(term_label), error = function(e) NULL)
  if (is.null(expr) || !is.call(expr) || length(expr) < 2L) {
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
  
  if (is.null(function_name) || !function_name %in% c("fp", "fp2") ||
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
  
  for (label in term_labels) {
    source_name <- if (label %in% factor_labels) {
      formula_factor_source_name(label)
    } else {
      formula_fp_source_name(label)
    }
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


#' Expand a Scalar df Default for Explicitly Mapped Terms
#'
#' A scalar `df` is the default complexity for ordinary continuous predictors.
#' Explicitly mapped terms, including multi-column categorical contrasts and
#' one-column binary-factor mappings, are supplied fixed linear design blocks
#' and therefore receive `df = 1` before cardinality-based df assignment.
#'
#' Per-column df vectors are intentionally returned unchanged so that the
#' existing grouped-term validation can reject explicit nonlinear settings.
#'
#' @param df Numeric scalar or per-column vector.
#' @param vnames Raw design-matrix column names.
#' @param term_to_columns Complete conceptual-term lookup.
#'
#' @return `df` unchanged when it is already per-column; otherwise a named
#'   per-column vector with explicitly mapped columns set to 1.
#' @keywords internal
#' @noRd
expand_scalar_df_for_mapped_terms <- function(df, vnames, term_to_columns) {
  if (length(df) != 1L) {
    return(df)
  }
  
  df_default <- stats::setNames(rep(df, length(vnames)), vnames)
  mapped_terms <- names(term_to_columns)[mapped_term_flags(term_to_columns)]
  
  if (length(mapped_terms) > 0L) {
    mapped_columns <- unlist(
      term_to_columns[mapped_terms],
      use.names = FALSE
    )
    df_default[mapped_columns] <- 1L
  }
  
  df_default
}


#' Apply Automatic Scale Defaults to Explicitly Mapped Terms
#'
#' Explicitly mapped terms, including multi-column categorical contrasts and
#' one-column binary-factor mappings, are supplied fixed model-matrix blocks.
#' Columns whose scale remains automatic therefore receive `scale = 1`; any
#' explicit scalar or named per-column setting remains authoritative.
#'
#' @param scale A normalized named per-column numeric vector.
#' @param vnames Raw design-matrix column names.
#' @param term_to_columns Complete conceptual-term lookup.
#' @param automatic Logical scalar or per-column logical vector indicating
#'   which scale values remain automatic. Named vectors are aligned to `vnames`.
#'
#' @return `scale`, with automatic explicitly mapped columns set to 1.
#' @keywords internal
#' @noRd
expand_scale_for_mapped_terms <- function(scale,
                                          vnames,
                                          term_to_columns,
                                          automatic = FALSE) {
  scale_default <- scale[vnames]
  
  if (length(automatic) == 1L) {
    if (!is.logical(automatic) || anyNA(automatic)) {
      stop("Internal scale-automatic metadata is malformed.", call. = FALSE)
    }
    automatic <- stats::setNames(rep(automatic, length(vnames)), vnames)
  } else {
    if (!is.logical(automatic) || length(automatic) != length(vnames) ||
        anyNA(automatic)) {
      stop("Internal scale-automatic metadata is malformed.", call. = FALSE)
    }
    if (!is.null(names(automatic))) {
      if (anyNA(names(automatic)) || any(!nzchar(names(automatic))) ||
          anyDuplicated(names(automatic)) ||
          !setequal(names(automatic), vnames)) {
        stop("Internal scale-automatic metadata is malformed.", call. = FALSE)
      }
      automatic <- automatic[vnames]
    } else {
      names(automatic) <- vnames
    }
  }
  
  mapped_terms <- names(term_to_columns)[mapped_term_flags(term_to_columns)]
  
  if (length(mapped_terms) > 0L) {
    mapped_columns <- unlist(
      term_to_columns[mapped_terms],
      use.names = FALSE
    )
    automatic_mapped <- mapped_columns[automatic[mapped_columns]]
    scale_default[automatic_mapped] <- 1
  }
  
  scale_default
}


#' Expand Conceptual-Term Metadata to Raw Design Columns
#'
#' Grouped terms store one set of selection and transformation settings for a
#' conceptual term, while model fitting and prediction operate on the raw
#' columns in its design block. This helper performs that expansion in one
#' place. Explicitly mapped terms are always represented as ordinary linear
#' columns and cannot inherit continuous-only transformations.
#'
#' @param term_to_columns Named list mapping conceptual terms to raw columns.
#' @param powers Named list of selected powers by conceptual term.
#' @param raw_columns Raw columns to return, in the required output order.
#' @param center,acdx,zero,catzero,spike Optional named term-level settings.
#' @param spike_decision Optional named integer decision codes by term.
#' @param acd_parameter Optional named list of fitted ACD parameters by term.
#'
#' @return Named list of metadata aligned with \code{raw_columns}.
#'
#' @keywords internal
#' @noRd
expand_term_metadata_to_columns <- function(
    term_to_columns,
    powers,
    raw_columns = unlist(term_to_columns, use.names = FALSE),
    center = NULL,
    acdx = NULL,
    zero = NULL,
    catzero = NULL,
    spike = NULL,
    spike_decision = NULL,
    acd_parameter = NULL) {
  if (!is.list(term_to_columns) || is.null(names(term_to_columns)) ||
      anyDuplicated(names(term_to_columns))) {
    stop("The conceptual-term column lookup is malformed.", call. = FALSE)
  }
  
  invalid_terms <- vapply(
    term_to_columns,
    function(columns) {
      !is.character(columns) || length(columns) == 0L || anyNA(columns) ||
        any(!nzchar(columns))
    },
    logical(1L)
  )
  all_columns <- unlist(term_to_columns, use.names = FALSE)
  if (any(invalid_terms) || anyDuplicated(all_columns)) {
    stop("The conceptual-term column lookup contains invalid raw columns.",
         call. = FALSE)
  }
  
  raw_columns <- as.character(raw_columns)
  unknown <- setdiff(raw_columns, all_columns)
  if (length(unknown) > 0L) {
    stop(
      "Unknown raw column(s): ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }
  
  column_to_term <- stats::setNames(
    rep(names(term_to_columns), lengths(term_to_columns)),
    all_columns
  )
  raw_terms <- unname(column_to_term[raw_columns])
  grouped_by_term <- mapped_term_flags(term_to_columns)
  grouped <- unname(grouped_by_term[raw_terms])
  
  term_value <- function(values, term, default) {
    if (is.null(values) || is.null(names(values)) ||
        !term %in% names(values)) {
      return(default)
    }
    values[[term]]
  }
  
  powers_columns <- lapply(seq_along(raw_columns), function(index) {
    term <- raw_terms[[index]]
    term_powers <- term_value(powers, term, NA_real_)
    
    if (grouped[[index]]) {
      if (length(term_powers) == 0L || all(is.na(term_powers))) {
        NA_real_
      } else {
        1
      }
    } else {
      term_powers
    }
  })
  names(powers_columns) <- raw_columns
  
  center_columns <- acdx_columns <- zero_columns <- catzero_columns <-
    spike_columns <- logical(length(raw_columns))
  spike_decision_columns <- integer(length(raw_columns))
  acd_parameter_columns <- vector("list", length(raw_columns))
  
  for (index in seq_along(raw_columns)) {
    term <- raw_terms[[index]]
    
    center_columns[[index]] <- isTRUE(term_value(center, term, FALSE))
    if (grouped[[index]]) {
      acdx_columns[[index]] <- FALSE
      zero_columns[[index]] <- FALSE
      catzero_columns[[index]] <- FALSE
      spike_columns[[index]] <- FALSE
      spike_decision_columns[[index]] <-
        saz_decision_codes[["continuous_only"]]
      acd_parameter_columns[index] <- list(NULL)
    } else {
      acdx_columns[[index]] <- isTRUE(term_value(acdx, term, FALSE))
      zero_columns[[index]] <- isTRUE(term_value(zero, term, FALSE))
      catzero_columns[[index]] <- isTRUE(term_value(catzero, term, FALSE))
      spike_columns[[index]] <- isTRUE(term_value(spike, term, FALSE))
      spike_decision_columns[[index]] <- as.integer(term_value(
        spike_decision,
        term,
        saz_decision_codes[["continuous_only"]]
      ))
      acd_parameter_columns[index] <- list(
        term_value(acd_parameter, term, NULL)
      )
    }
  }
  
  names(center_columns) <- names(acdx_columns) <- names(zero_columns) <-
    names(catzero_columns) <- names(spike_columns) <-
    names(spike_decision_columns) <- names(acd_parameter_columns) <-
    raw_columns
  
  list(
    terms = stats::setNames(raw_terms, raw_columns),
    powers = powers_columns,
    center = center_columns,
    acdx = acdx_columns,
    zero = zero_columns,
    catzero = catzero_columns,
    spike = spike_columns,
    spike_decision = spike_decision_columns,
    acd_parameter = acd_parameter_columns
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
#' @param reason Optional explanation appended to the error. When omitted, a
#'   setting-specific explanation is supplied for continuous-only options and
#'   fractional-polynomial search controls.
#'
#' @return Invisibly returns \code{TRUE}.
#' @keywords internal
#' @noRd
validate_grouped_term_setting <- function(term_to_columns,
                                          values,
                                          setting,
                                          predicate,
                                          requirement,
                                          reason = NULL) {
  if (is.null(reason)) {
    if (setting %in% c(
      "acdx", "acd_vars", "zero_vars", "catzero_vars", "spike_vars"
    )) {
      reason <- paste(
        "Grouped terms are fixed categorical design blocks;",
        "continuous-variable transformations and structural-zero options",
        "apply only to singleton continuous predictors."
      )
    } else if (setting %in% c("df", "force_max_fp_vars")) {
      reason <- paste(
        "Grouped terms are fixed linear design blocks and cannot enter the",
        "fractional-polynomial transformation search."
      )
    }
  }
  
  reason_suffix <- if (
    is.character(reason) && length(reason) == 1L && !is.na(reason) &&
    nzchar(reason)
  ) {
    paste0(" ", reason)
  } else {
    ""
  }
  
  grouped_terms <- names(term_to_columns)[mapped_term_flags(term_to_columns)]
  for (term in grouped_terms) {
    cols <- term_to_columns[[term]]
    valid <- predicate(values[cols])
    if (length(valid) != length(cols) || anyNA(valid) || !all(valid)) {
      stop(
        sprintf(
          paste0(
            "! Grouped term '%s' requires `%s` to be %s for every member ",
            "column: %s.%s"
          ),
          term,
          setting,
          requirement,
          paste(cols, collapse = ", "),
          reason_suffix
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

#' @describeIn mfp2 Matrix interface: accepts a numeric predictor matrix `x`
#' and response vector `y`. Categorical predictors must be pre-coded as
#' numeric columns; group related columns with `term_groups`.
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
                         force_max_fp_vars = NULL,
                         term_groups = NULL,
                         verbose = TRUE,
                         ...
) {
  
  # mfp2.default() does not fit the model itself. It validates and normalizes
  # every argument, derives shift/scale/df settings and zero/catzero/spike/acd
  # flags for each predictor, transforms `x` accordingly, and then delegates
  # the actual multivariable FP selection and fitting to fit_mfp().
  
  # Step 1: Capture the call and rename public-facing arguments -----------------
  cl <- match.call()
  
  # Display the public generic name rather than the S3 method name.
  cl[[1L]] <- quote(mfp2)
  
  # Public interface:
  # zero_vars, catzero_vars, spike_vars and force_max_fp_vars are character
  # vectors of variable names.
  #
  # Internal representation:
  # The existing internal names are retained because they are converted below
  # to named logical vectors over the columns of x.
  zero <- zero_vars
  catzero <- catzero_vars
  spike <- spike_vars
  force_max_fp <- force_max_fp_vars
  
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
  
  # Formula dispatch can attach a private, full-row preprocessing matrix to
  # the subset-specific fitting design. Extract it immediately and work with
  # two explicit matrices thereafter:
  #   * x: the design used for validation, selection, and fitting;
  #   * preprocess_x: the source for automatic full-data preprocessing choices.
  # Direct matrix calls carry no attribute, so both names refer to the same x.
  formula_generated_settings <- !is.null(
    attr(x, "mfp2_preprocess_x", exact = TRUE)
  )
  preprocessing <- extract_preprocess_matrix(x)
  x <- preprocessing$x
  preprocess_x <- preprocessing$preprocess_x
  x_input <- x
  
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
    
    if (anyDuplicated(subset)) {
      stop(
        "! `subset` must not contain duplicated row indices.",
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
  
  # Direct matrix users may supply NULL, an unnamed global scalar, or named
  # per-column settings. Named shift and scale vectors may be partial;
  # unspecified entries remain NA so the existing preprocessing step estimates
  # them automatically. Formula methods may also use NA in internally
  # constructed setting vectors.
  shift <- normalize_named_numeric_setting(
    value = shift,
    column_names = vnames,
    argument_name = "shift",
    allow_null = TRUE,
    scalar_recycle = TRUE,
    strictly_positive = FALSE,
    allow_na = formula_generated_settings,
    allow_partial_named = TRUE
  )
  scale <- normalize_named_numeric_setting(
    value = scale,
    column_names = vnames,
    argument_name = "scale",
    allow_null = TRUE,
    scalar_recycle = TRUE,
    strictly_positive = TRUE,
    allow_na = formula_generated_settings,
    allow_partial_named = TRUE
  )
  scale <- expand_scale_for_mapped_terms(
    scale = scale,
    vnames = vnames,
    term_to_columns = term_to_columns,
    automatic = is.na(scale)
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
  validate_variable_names(force_max_fp,"force_max_fp_vars", vnames)
  
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
  
  # `shift` was normalized above. NA marks a variable for automatic shift
  # estimation further below; supplied values are already aligned by name.
  
  # Expand scalar center to one value per predictor.
  if (length(center) == 1) {
    center <- rep(center, nvars)    
  }
  
  # Convert the public character-vector interface into the named logical vector
  # expected by the internal MFP fitting functions.
  force_max_fp <- setNames(
    vnames %in% force_max_fp,
    vnames
  )
  
  # Under p-value selection, forcing the maximum FP requires both retention of
  # the variable and acceptance of the most complex permitted FP degree.
  if (criterion == "pvalue" && any(force_max_fp)) {
    select[force_max_fp] <- 1
    alpha[force_max_fp] <- 1
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
      apply(preprocess_x[, vars_to_check, drop = FALSE], 2, function(col) {
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
      apply(preprocess_x[, vars_to_check_cz, drop = FALSE], 2, function(col) {
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
  validate_grouped_term_setting(
    term_to_columns, force_max_fp, "force_max_fp_vars", function(v) !v, "FALSE"
  )
  
  # Identify binary columns, exactly two unique non-missing values.
  # Binary variables are not eligible for zero/catzero/spike handling here:
  #   - zero handling is unnecessary because the variable is already discrete;
  #   - catzero would create a redundant indicator;
  #   - spike implies catzero and would therefore create the same inconsistency.
  binary_vars <- apply(preprocess_x, 2, function(col) length(unique(col[!is.na(col)])) == 2)
  binary_names <- vnames[binary_vars]
  
  # Binary predictors must remain on their original two-level scale. Automatic
  # preprocessing already returns scale = 1 for them via find_scale_factor(),
  # but an explicitly supplied scale would otherwise override that safeguard.
  # Reset it here, before scaling is estimated and applied below.
  if (length(binary_names) > 0L) {
    scale[binary_names] <- 1
  }
  
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
      x                      = preprocess_x,
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
    # A scalar df controls ordinary continuous predictors. Explicitly mapped
    # terms are supplied fixed linear design blocks, so assign df = 1 to their
    # raw columns before applying cardinality rules to the remaining columns.
    df_default <- expand_scalar_df_for_mapped_terms(
      df = df,
      vnames = vnames,
      term_to_columns = term_to_columns
    )
    
    if (any(df_default != 1L)) {
      df.list <- assign_df(x = preprocess_x, df_default = df_default)
    } else {
      df.list <- df_default
    }
  } else {
    nux <- apply(preprocess_x, 2, function(v) length(unique(v)))
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
      preprocess_x[, shift_missing, drop = FALSE],
      2,
      find_shift_factor
    )
  }
  
  # Step 16: Shift and scale the design matrix -----------------------------------
  # Shift each column of x so that fractional-polynomial powers (which require
  # positive values) can be computed; see "Details on shifting, scaling,
  # centering" in the documentation above.
  preprocess_x_shifted <- sweep(preprocess_x, 2, shift, "+")
  x <- sweep(x, 2, shift, "+")
  
  # Scaling is estimated after shifting because find_scale_factor() operates on
  # shifted values. `scale` is already validated, named, and column-aligned.
  scale_missing <- is.na(scale)
  
  if (any(scale_missing)) {
    scale[scale_missing] <- apply(
      preprocess_x_shifted[, scale_missing, drop = FALSE],
      2,
      find_scale_factor
    )
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
    xd <- preprocess_x_shifted[, check_indices, drop = FALSE]
    
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
  
  # Step 19: Apply `subset`, after full-data preprocessing ----------------------
  # Matrix calls reach this branch with their original numeric design. Before
  # rows are removed, verify that every explicitly mapped block keeps the same
  # estimable dimension. Unlike the formula method, this path has no factor-level
  # or contrast recipe from which to rebuild a reduced categorical design.
  if (!is.null(subset)) {
    validate_grouped_subset_rank(
      x_full = x_input,
      x_fit = x_input[subset, , drop = FALSE],
      term_to_columns = term_to_columns
    )
    
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
    shift, term_to_columns, "shift"
  )
  scale_term <- collapse_option_to_terms(
    scale, term_to_columns, "scale"
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


#' @describeIn mfp2 Formula interface: accepts a model formula and data frame.
#' Factors are expanded automatically. Use [fp()] or [fp2()] to set
#' variable-specific options.
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
  
  # Step 4: Evaluate `subset` and validate observation-level inputs ------------
  n_data <- nrow(data)
  
  # Match the data-mask semantics used by standard formula methods. Capture the
  # user's expression before the `subset` promise is forced, evaluate it exactly
  # once with columns of `data` taking precedence, and then reuse only the
  # resolved value. The normalized formula environment inherits the user's
  # formula environment, so caller-scope objects remain available as a fallback.
  subset_expr <- substitute(subset)
  subset <- eval(
    subset_expr,
    envir = data,
    enclos = environment(formula_internal)
  )
  
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
      if (anyDuplicated(subset)) {
        stop(
          "! `subset` must not contain duplicated row indices.",
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
  
  # Step 6: Build full-data and fitted model frames ------------------------------
  # The formula interface deliberately keeps two model frames:
  #   * mf_full supplies the original continuous values used by automatic
  #     preprocessing, preserving the documented pre-subset behavior;
  #   * mf reuses mf_full when subset is NULL; otherwise it is rebuilt on
  #     fit_rows and supplies the response, fitted factor levels, contrasts,
  #     design matrix, and prediction metadata.
  #
  # Use na.pass so the global na.action option cannot silently alter either row
  # set. Missingness is diagnosed explicitly against the full model frame below.
  mf_full <- stats::model.frame(
    terms_formula,
    data = data,
    drop.unused.levels = TRUE,
    na.action = stats::na.pass
  )
  # Keep the ordinary formula path independent of subset-specific helpers.
  # A genuine subset is resolved once; otherwise all evaluated model-frame rows
  # are retained in their original order.
  if (is.null(subset)) {
    fit_rows <- seq_len(nrow(mf_full))
  } else {
    fit_rows <- formula_subset_rows(subset, nrow(mf_full))
  }
  if (length(fit_rows) < 5L) {
    stop(
      sprintf("! The selected subset is too small (<5) to fit an mfp model.\ni The number of selected observations is %d.", length(fit_rows)),
      call. = FALSE
    )
  }
  
  
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
  if (anyNA(mf_full)) {
    na_rows <- which(!stats::complete.cases(mf_full))
    
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
  
  # With no subset, reuse the complete evaluated frame unchanged. This keeps the
  # ordinary formula path simple and avoids unnecessary row copying or factor
  # reconstruction. Only a genuine subset needs a retained-row frame with unused
  # levels removed before model.matrix() rebuilds treatment or ordered contrasts.
  if (is.null(subset)) {
    mf <- mf_full
  } else {
    # Avoid passing the local `fit_rows` symbol through model.frame(subset = ...):
    # model.frame() evaluates that expression in the formula/data environment,
    # where a method-local object is not visible.
    mf <- subset_formula_model_frame(mf_full, fit_rows)
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
  
  if (is.null(term_offset) && !is.null(offset) && !is.null(subset)) {
    offset <- offset[fit_rows]
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
  
  validate_formula_factor_levels(terms_model, mf)
  
  # Extract the fitted response and build both designs independently. `mm` is
  # authoritative for categorical columns, term mappings, fitting, and prediction.
  # `mm_full` is only a source of full-data values for automatic preprocessing of
  # stable, same-named continuous columns; its factor coding is never fitted.
  y <- model.extract(mf, "response")
  
  mm_full <- stats::model.matrix(terms_model, mf_full)
  mm <- stats::model.matrix(terms_model, mf)
  
  # Identify unordered and ordered factor terms before the intercept is
  # removed. Their model-matrix columns will be grouped and forced to fixed
  # linear settings below. Ordered factors therefore retain their configured
  # contrasts (contr.poly by default) but are selected as one conceptual term.
  factor_term_labels <- identify_formula_factor_terms(terms_model, mf)
  conceptual_term_names <- formula_conceptual_term_names(terms_model, mf)
  # Retain the original formula-label -> conceptual-term relationship for
  # prediction. This map is updated below when fp()/fp2() expressions are
  # renamed to their source-variable names.
  formula_prediction_term_names <- conceptual_term_names
  factor_terms <- unname(conceptual_term_names[factor_term_labels])
  
  assign <- attr(mm, "assign")
  term_labels <- attr(terms_model, "term.labels")
  
  # mfp2 fits models without an intercept term internally (the intercept is
  # handled by fit_mfp()/glm()/coxph()); drop the (Intercept) column here but
  # keep `assign` aligned with the remaining columns of x.
  keep_cols <- colnames(mm) != "(Intercept)"
  x <- mm[, keep_cols, drop = FALSE]
  assign <- assign[keep_cols]
  full_keep_cols <- colnames(mm_full) != "(Intercept)"
  x_full <- mm_full[, full_keep_cols, drop = FALSE]
  
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
    
    # Capture the raw fp()/fp2() model-matrix column names before renaming.
    fp_columns <- names_x[is_fp_term(names_x)]
    if (length(fp_columns) != length(fp_vars)) {
      stop(
        "! Internal formula parsing error: the number of fp()/fp2() terms in the model matrix does not match the model frame.",
        call. = FALSE
      )
    }
    
    # Rename fp()/fp2() model-matrix columns to the underlying variable name,
    # e.g. "fp(age)" becomes "age", so downstream code (and the fitted object)
    # exposes ordinary variable names rather than formula syntax.
    names_x <- replace(names_x, which(is_fp_term(names_x)), fp_vars)
    colnames(x) <- names_x
    full_fp_pos <- which(is_fp_term(colnames(x_full)))
    if (length(full_fp_pos) != length(fp_vars)) {
      stop(
        "! Internal formula parsing error: full-data and fitted fp()/fp2() columns do not align.",
        call. = FALSE
      )
    }
    colnames(x_full)[full_fp_pos] <- fp_vars
    
    # The conceptual lookup is keyed by source-variable name. Replace its raw
    # formula-generated column with the renamed column without creating a second
    # conceptual entry for the original fp() label.
    for (i in seq_along(fp_columns)) {
      fp_column <- fp_columns[[i]]
      fp_var <- fp_vars[[i]]
      if (fp_var %in% names(term_to_columns)) {
        term_to_columns[[fp_var]] <- ifelse(
          term_to_columns[[fp_var]] == fp_column,
          fp_var,
          term_to_columns[[fp_var]]
        )
      }
    }
  }
  
  # Step 10: Build per-variable default option lists from the global scalars ---
  # Every predictor starts out with the global df/scale/shift/center/alpha/select
  # defaults; fp()-term-specific settings (Step 11) then override individual
  # entries of these lists where the user supplied them. Cardinality-based df
  # defaults retain the original full-data source for continuous predictors.
  x_preprocess <- align_formula_preprocess_matrix(x, x_full)
  df_list <- setNames(
    as.list(assign_df(x = x_preprocess, df_default = df)),
    names_x
  )
  
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
    
    # Accept final model-matrix columns, conceptual term names, and original
    # formula labels such as fp(x). Formula labels resolve through the explicit
    # fit-time label-to-conceptual-name map.
    valid_keep <- unique(c(
      colnames(x),
      names(term_to_columns),
      names(formula_prediction_term_names)
    ))
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
        # Exact final column match, e.g. keep = "age".
        expanded_keep <- c(expanded_keep, k)
      } else {
        conceptual_name <- if (k %in% names(formula_prediction_term_names)) {
          unname(formula_prediction_term_names[[k]])
        } else {
          k
        }
        if (conceptual_name %in% names(term_to_columns)) {
          # Exact conceptual or original formula-term match, including fp(x).
          expanded_keep <- c(
            expanded_keep,
            term_to_columns[[conceptual_name]]
          )
        }
      }
    }
    
    keep <- unique(expanded_keep)
  }
  
  # Step 14: Delegate using the subset-specific fitting design ------------------
  # The private attribute transports the aligned full-data preprocessing source
  # into mfp2.default(); it is removed immediately on entry there. Responses and
  # other observation-level inputs are subset here, so pass subset = NULL below
  # and avoid applying the row selection twice. Formula metadata is derived only
  # from mf/mm, never from the full-data categorical representation.
  x <- attach_formula_preprocess_matrix(x, x_full)
  if (!is.null(subset)) {
    if (!is.null(weights)) weights <- weights[fit_rows]
    if (!is.null(strata) && !has_formula_strata) {
      strata <- if (is.matrix(strata) || is.data.frame(strata)) {
        strata[fit_rows, , drop = FALSE]
      } else {
        strata[fit_rows]
      }
    }
  }
  
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
  
  # Convert the per-column logical settings collected from fp() terms into the
  # public variable-name interface expected by mfp2.default().
  force_max_fp_values <- unlist(
    force_max_fp_list,
    use.names = TRUE
  )
  
  force_max_fp_vars <- unique(
    names(force_max_fp_values)[force_max_fp_values]
  )
  
  if (length(force_max_fp_vars) == 0L) {
    force_max_fp_vars <- NULL
  }
  
  # Construct explicit per-column numeric vectors from the formula defaults and
  # fp()/fp2() overrides. vapply() preserves the outer variable names while
  # discarding any incidental name attached to an individual scalar value.
  scale_vec <- vapply(scale_list, function(value) {
    as.numeric(value[[1L]])
  }, numeric(1L))
  shift_vec <- vapply(shift_list, function(value) {
    as.numeric(value[[1L]])
  }, numeric(1L))
  
  fit <- mfp2.default(x = x, 
                      y = y,
                      term_groups = grouped_formula_terms,
                      weights = weights, 
                      offset = offset, 
                      cycles = cycles,
                      scale = scale_vec,
                      shift = shift_vec,
                      df = unlist(df_list), 
                      center = unlist(center_list),
                      subset = NULL,
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
                      force_max_fp_vars = force_max_fp_vars,
                      min_saz_component_prop = min_saz_component_prop,
                      verbose = verbose
  )
  fit$formula_interface <- TRUE
  fit$formula <- formula_user
  fit$formula_terms <- stats::delete.response(terms_model)
  fit$formula_contrasts <- attr(mm, "contrasts")
  fit$formula_xlevels <- .getXlevels(terms_model, mf)
  fit$formula_column_map <- formula_column_map
  fit$formula_model_matrix_columns <- names_x
  fit$formula_term_to_columns <- term_to_columns
  fit$formula_prediction_term_names <- formula_prediction_term_names
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

#' Extract Coefficients from an `mfp2` Model
#'
#' Returns the named coefficient vector from a fitted [mfp2()] model. This is a
#' method for the generic [stats::coef()] function.
#'
#' @param object A fitted [mfp2()] object.
#' @param ... Not used.
#'
#' @return
#' Named numeric vector of coefficients from the final selected model.
#'
#' @seealso [mfp2()], [summary.mfp2()], [predict.mfp2()]
#'
#' @export
coef.mfp2 <- function(object, ...) {
  object$coefficients
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
#' The subsequent "Detailed Settings" table reproduces the original,
#' unabridged \code{fp_terms} columns, also sorted with selected variables first.
#' Notes about \code{select}/\code{alpha} divergence and the
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
#' @seealso [mfp2()], [summary.mfp2()], [get_selected_variable_names()]
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
  # Step 5: Print the original mfp2 call and one-line meta block
  # ---------------------------------------------------------------------------
  #
  # A plain "Call:" label (no dashed rule) opens the output, followed by the
  # deparsed mfp2() call, a one-line family / criterion / convergence block,
  # and an observation/event line. This matches the opening produced by
  # print.summary.mfp2(); the two methods therefore present the call and
  # top-level metadata identically. The dashed-rule section headings resume
  # only from Selection Summary onward, where they bracket the more
  # substantial tabular sections.
  
  cat("Call:\n")
  if (!is.null(x$call_mfp)) {
    print(x$call_mfp)
  } else {
    cat("(original mfp2 call unavailable)\n")
  }
  
  # Family / Criterion / Converged: three fields on one line, joined by " | "
  # for compactness. The criterion label was pre-formatted earlier in this
  # method; convergence uses the same helper the Selection Summary previously
  # used, so the yes/no wording is unchanged.
  cat(sprintf(
    "Family: %s | Criterion: %s | Converged: %s\n",
    x$family_string,
    criterion_label,
    format_convergence(x$convergence_mfp)
  ))
  
  # Observations / Events: for Cox models the number of observed events is
  # more informative than the total sample size alone, so both are shown.
  # For non-Cox families only the observation count is meaningful.
  n_obs <- tryCatch(NROW(x$y), error = function(e) NA_integer_)
  if (identical(x$family_string, "cox") && inherits(x$y, "Surv")) {
    status_col <- ncol(x$y)
    n_events <- tryCatch(
      sum(x$y[, status_col] == 1, na.rm = TRUE),
      error = function(e) NA_integer_
    )
    cat(sprintf("Observations: %s | Events: %s\n", n_obs, n_events))
  } else {
    cat(sprintf("Observations: %s\n", n_obs))
  }
  
  cat("\n")
  
  
  # ---------------------------------------------------------------------------
  # Step 6: Print the model-selection summary
  # ---------------------------------------------------------------------------
  #
  # The Selection Summary section previously opened with "Converged:" and
  # "Criterion:" lines. Both fields are now shown in the top meta block above,
  # so they are omitted here to avoid duplication. The section retains its
  # dashed-rule heading and its three "Selected variables" lines.
  
  print_section_heading("Selection Summary")
  
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
  # Step 11: Render the Model Fit block via the shared helper
  # ---------------------------------------------------------------------------
  #
  # The Model Fit block (formerly "Model Deviances") is produced by
  # mfp2_format_model_fit_block(), the same helper called by
  # print.summary.mfp2(). This guarantees the two methods display the same
  # three numbers under the same convention (-2 log L on all rows; df on the
  # MFP-adjusted scale) and the same explanatory note, so there is exactly
  # one place in the package where the block's format lives.
  mfp2_format_model_fit_block(
    values          = mfp2_summary_model_fit_values(x),
    digits          = digits,
    heading_printer = print_section_heading
  )
  
  cat("\n")
  cat(make_boundary("="), "\n", sep = "")
  
  invisible(x)
}

#' Specify Fractional-Polynomial Options in a Formula
#'
#' Marks a continuous variable for fractional-polynomial modelling by
#' [mfp2()] or [mfpi()] and specifies any variable-specific modelling options.
#' `fp()` is intended for use inside a model formula. It does not transform
#' the variable immediately; the transformation is selected and constructed
#' during model fitting.
#' 
#' @details
#' By default, the function allows selection among omission, a linear effect,
#' FP1, and FP2, subject to the model-selection settings supplied to [mfp2()]
#' or [mfpi()].
#'
#' Variable-specific options can be used to restrict the maximum complexity,
#' force the variable into the model, request ACD or spike-at-zero modelling,
#' or control how zero and nonpositive values are handled.
#'
#' `fp()` should normally be applied only to continuous numeric variables.
#' Categorical variables should be included directly in the formula as
#' factors.
#'
#' @param x A vector representing a continuous variable undergoing
#'   FP-transformation.
#' @param df Maximum degrees of freedom for this variable's FP function.
#'   Must be `1` (linear) or a positive even integer. Default `4` (FP2).
#'   See [mfp2()] for the cardinality-based reduction rules.
#' @param alpha Significance level for the closed test between FP functions of
#'   different degrees. Default `0.05`. Overrides the global `alpha` from
#'   [mfp2()] for this variable.
#' @param select Significance level for backward-elimination variable
#'   selection. Default `0.05`. Set to `1` to force this variable into the
#'   model. Overrides the global `select` from [mfp2()].
#' @param shift Numeric shift override for this variable, or `NULL` (default)
#'   to inherit the global setting from [mfp2()]. See [mfp2()] for details.
#' @param scale Numeric scale override for this variable, or `NULL` (default)
#'   to inherit the global setting from [mfp2()]. See [mfp2()] for details.
#' @param center Logical centering override for this variable. Default `TRUE`.
#'   See [mfp2()] for details.
#' @param acdx Logical. If `TRUE`, request ACD modelling for this variable.
#'   Default `FALSE`. See the ACD modelling section in [mfp2()].
#' @param zero Logical. If `TRUE`, only positive values undergo the FP
#'   function; nonpositive values contribute zero to the linear predictor.
#'   Default `FALSE`. Equivalent to listing the variable in `zero_vars` in
#'   `mfp2.default()`.
#' @param catzero Logical. If `TRUE`, a fixed binary indicator for nonpositive
#'   values is added alongside the positive-component FP function. Implies
#'   `zero = TRUE`. Default `FALSE`. Equivalent to listing the variable in
#'   `catzero_vars` in `mfp2.default()`.
#' @param spike Logical. If `TRUE`, request two-stage spike-at-zero (SAZ)
#'   selection for this variable. While SAZ remains eligible, implies both
#'   `catzero` and `zero`. Default `FALSE`. Equivalent to listing the variable
#'   in `spike_vars` in `mfp2.default()`. See the SAZ section in [mfp2()] for
#'   eligibility requirements.
#' @param force_max_fp Logical. If `TRUE`, the most complex functional form
#'   allowed by `df` is forced for this variable, bypassing both variable
#'   selection and functional-form simplification. Under `criterion = "pvalue"`,
#'   this is implemented by internally setting `select = 1` and `alpha = 1`
#'   for the variable; the option also applies under `"aic"` and `"bic"`.
#'   Default `FALSE`. Equivalent to naming the variable in `force_max_fp_vars`
#'   in `mfp2.default()`.
#' @param powers Numeric vector of candidate powers to be evaluated for `x`.
#'   If `NULL` (default), the global candidate set from [mfp2()] is used,
#'   which defaults to `c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)`. When `df > 1`,
#'   the candidate set must include at least one power other than `1` so the
#'   algorithm can search over nonlinear forms. Duplicate values are removed
#'   automatically. Overrides the corresponding entry in the `powers` list
#'   argument of [mfp2()].
#' @param ... Not used directly by `fp()`. Present so that the `fp2()` alias
#'   can forward all arguments.
#'
#' @return
#' The input vector `x` with modelling options attached as attributes. The
#' returned value is interpreted by [mfp2()] or [mfpi()] during model fitting;
#' calling `fp()` outside a formula has no modelling effect.
#'
#' @seealso [mfp2()], [mfpi()], [fp2()]
#'
#' @examples
#' data("prostate")
#'
#' # Allow the default fractional-polynomial search for age.
#' fit1 <- mfp2(
#'   lpsa ~ fp(age) + svi,
#'   data = prostate,
#'   verbose = FALSE
#' )
#'
#' # Restrict age to at most one degree of freedom.
#' fit2 <- mfp2(
#'   lpsa ~ fp(age, df = 1) + svi,
#'   data = prostate,
#'   verbose = FALSE
#' )
#'
#' \donttest{
#' # Request ACD modelling for a continuous predictor.
#' fit_acd <- mfp2(
#'   lpsa ~ fp(cavol, acdx = TRUE) + fp(age) + svi,
#'   data = prostate,
#'   verbose = FALSE
#' )
#' }
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

#' Extract Names of Selected Variables from an `mfp2` Model
#'
#' Returns the names of predictors retained in the final model after MFP
#' selection. A variable is considered selected when at least one of its
#' fractional-polynomial powers is not `NA`. Excluded variables have all
#' powers set to `NA`.
#'
#' @param object A fitted [mfp2()] object.
#'
#' @examples
#' # Gaussian model
#' data("prostate")
#' x <- as.matrix(prostate[, 2:8])
#' y <- as.numeric(prostate$lpsa)
#' fit <- mfp2(x, y, verbose = FALSE)
#' get_selected_variable_names(fit)
#'
#' @return
#' Character vector of selected predictor names, ordered according to the
#' `xorder` used in the [mfp2()] call.
#'
#' @seealso [mfp2()], [print.mfp2()]
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