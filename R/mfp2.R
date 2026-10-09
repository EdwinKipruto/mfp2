#' Multivariable Fractional Polynomial Models with Extensions
#'
#' `mfp2()` selects predictors for a regression model and chooses a linear or
#' fractional-polynomial form for each selected continuous predictor. It also
#' supports approximate cumulative distribution (ACD) transformations for
#' sigmoid-shaped relationships and spike-at-zero modelling for predictors
#' with many exact zeros. Models can be specified with a formula and data frame
#' or with a numeric predictor matrix and an outcome. Both interfaces use the
#' same selection procedure. The formula interface handles categorical predictors
#' automatically; the matrix interface requires them to be coded as numeric
#' columns. Supported regression models include Gaussian, binomial, Poisson, Gamma,
#' inverse Gaussian, negative binomial, multinomial, and ordinal models.
#' For time-to-event outcomes, it supports Cox, parametric survival,
#' and Fine-Gray competing-risks models. It also supports generalized
#' estimating equations (GEE) for clustered data.
#'
#' @section Fractional-polynomial model selection:
#' For a continuous predictor `x`, a linear effect has the form
#' \eqn{\beta_1 x}. A first-degree fractional polynomial (FP1) replaces
#' `x` with a selected power:
#'
#' \deqn{\beta_1 x^{p_1}.}
#'
#' A second-degree fractional polynomial (FP2) uses two powers:
#'
#' \deqn{\beta_1 x^{p_1} + \beta_2 x^{p_2}.}
#'
#' By convention, a power of zero means \eqn{\log(x)}. If the two powers
#' are equal, the second term is \eqn{x^p\log(x)}; for two zero powers,
#' the terms are \eqn{\log(x)} and \eqn{\{\log(x)\}^2}. Here `x` denotes
#' the positive input after any shifting and scaling. The default candidate
#' powers are `c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)`.
#'
#' The `df` setting limits the maximum FP degree: `df = 1` allows only a
#' linear effect, `df = 2` allows up to FP1, and the default `df = 4`
#' allows up to FP2. Higher degrees can be allowed. Predictors with only
#' two or three distinct values are restricted to `df = 1`. With four or
#' five distinct values, a default or global `df` is capped at `2`, but
#' an explicitly set predictor-specific value is retained. Selection may
#' still omit a predictor or choose a simpler form. Predictors in `keep`
#' cannot be omitted, but their functional forms may still be selected.
#' A selected FP1 with power `1` has the same fitted basis as the fixed
#' linear form, but retains the FP power-search degrees of freedom. The
#' `fp_terms$searched_fp` column distinguishes these cases: `TRUE` means
#' the final continuous form came from an FP candidate class; `FALSE` means
#' it came from the fixed linear candidate or no continuous form was retained.
#' The flag identifies the form that won selection, not whether FP candidates
#' were considered. For a single-predictor term with power `1`, this gives
#' `df_final = 2` for FP1 and `df_final = 1` for the fixed linear form.
#'
#' With `criterion = "pvalue"`, a closed-test-style procedure uses `select`
#' for predictor inclusion and `alpha` for functional-form comparisons.
#' With `"aic"` or `"bic"`, candidate models are compared using the chosen
#' information criterion.
#'
#' See `mfp` vignette for methodological details.
#'
#' @section Formula and matrix interfaces:
#' For most users, the formula interface is the simpler choice:
#' `mfp2(formula, data, ...)`. It handles factors automatically. Use [fp()]
#' or [fp2()] to set options for individual continuous predictors; other
#' predictors, including those added with `.`, use the applicable top-level
#' settings. An `fp()` term has its own default `select = 0.05`, so set
#' `select` inside the term if a different value is needed. Formulas may
#' also include `offset()` and, for survival models, `strata()`.
#'
#' The matrix interface, `mfp2(x, y, ...)`, takes a numeric predictor
#' matrix. Categorical predictors must first be coded as numeric columns.
#' If several columns represent one predictor, use `term_groups` to have
#' them selected and tested together.
#'
#' Both interfaces use the same selection procedure. Predictor names must
#' be unique, non-empty, and free of backticks. Names given in options such
#' as `keep` and `powers` must match predictors in the model.
#'
#' @section Model families and responses:
#' Set `family` to choose the model. Standard choices are `"gaussian"` for
#' numeric outcomes, `"binomial"` for binary outcomes, `"poisson"` for counts,
#' and `"Gamma"` or `"inverse.gaussian"` for positive numeric outcomes.
#' Use `"negbin"` for negative-binomial count regression. The first five
#' also accept a GLM family function or object, so you can choose a link
#' supported by [stats::glm()]. Quasi families are not supported.
#'
#' Binomial responses may be values from 0 to 1, a two-level factor, or a
#' two-column matrix of success and failure counts. Poisson and negative-
#' binomial responses must be nonnegative integer counts. Count matrices
#' must contain at least one trial in every row.
#'
#' For outcomes with at least three categories, use [multinomial_family()]
#' for unordered categories or [ordinal_family()] for ordered categories.
#' Multinomial models also accept a matrix of counts by category. See
#' **Multinomial models** and **Ordinal models** for response details.
#'
#' For survival outcomes, use `"cox"` with a right-censored [survival::Surv()]
#' response, [survreg_family()] for parametric survival regression, or
#' [finegray_family()] for competing risks with a multi-state `Surv`
#' response. The character values `"survreg"` and `"finegray"` use the
#' respective default settings.
#'
#' For clustered observations, use [gee_family()] with a Gaussian,
#' binomial, Poisson, or Gamma response family. The response follows the
#' rules for that family. See **Generalized estimating equations**.
#'
#' @section Generalized estimating equations:
#' Specify a GEE model with [gee_family()]. The required `id` argument
#' identifies clusters; the package groups their observations before fitting.
#' The optional `waves` argument specifies visit order within clusters.
#' [gee_family()] sets the response family, working correlation,
#' standard-error method, and scale settings. The final model is returned
#' as a `geeglm` object.
#'
#' With `criterion = "pvalue"`, selection uses approximate comparisons of
#' overall robust Wald statistics, following Stata's `mfp: xtgee`
#' calculation. During verbose fitting, the selection tables display
#' `Wald chi-sq` for candidates and `Wald diff.` for comparisons. The robust
#' covariance uses the finite-cluster correction \eqn{K/(K-1)}, where
#' \eqn{K} is the number of clusters.
#'
#' With `criterion = "aic"` or `"bic"`, selection uses quasi-likelihood
#' scores:
#'
#' \deqn{-2Q/\phi + 2p \quad\text{or}\quad -2Q/\phi + \log(n)p,}
#'
#' respectively. Here \eqn{Q} is the quasi-likelihood, \eqn{\phi} is a single
#' reference dispersion shared by all comparisons, and \eqn{n} is the number
#' of clusters by default (or observations if
#' `gee_family(qbic_penalty = "observations")`).
#' The reference dispersion is estimated once from a full model containing all
#' predictors at their maximum allowed FP degrees, unless supplied through
#' `gee_family(reference_dispersion = ...)`. Reference powers are fixed from
#' the permitted sets before the one reference fit. Geepack continues to
#' estimate scale for each GEE
#' fit; that fitted scale is used for fitting and reporting, not for these
#' information-criterion comparisons. \eqn{p} counts
#' mean-model coefficients, including the intercept, plus an allowance for
#' selecting FP powers. The first is an MFP-adjusted QICu score; the second
#' is a package-specific BIC-like score.
#'
#' For AIC or BIC selection, [print.mfp2()] and [summary.mfp2()] show the
#' full-linear and selected-MFP scores under `Selection Scores`. Smaller
#' scores are preferred. For GEE p-value selection, they omit the `Model Fit`
#' block; Wald comparisons appear during verbose fitting.
#'
#' @section Shifting, scaling, and centering:
#' FP functions involving logarithms or non-integer or negative powers need
#' positive inputs. By default, `mfp2()` adds a shift when needed and chooses
#' a power-of-ten scale for ordinary FP predictors to keep values numerically
#' manageable during selection. ACD predictors may be shifted but are not
#' scaled. Use `shift = 0` or `scale = 1` to turn off the corresponding
#' automatic step for ordinary FP predictors.
#'
#' In `mfp2.default()`, a single unnamed value applies to all predictors.
#' A named value applies only to the specified column: for example,
#' `shift = c(age = 20)` fixes the shift for `age`, while shifts for other
#' columns are chosen automatically. The same rule applies to `scale`.
#' Do not supply `NA` as an explicit shift or scale value.
#'
#' A supplied shift is not increased automatically. If it leaves values
#' nonpositive for a nonlinear FP that requires positive inputs, fitting
#' stops with an error. Linear predictors (`df = 1`) and predictors using
#' `zero`, `catzero`, or `spike` receive a shift of zero.
#'
#' With `center = TRUE`, transformed continuous columns are centered on
#' their mean in the fitted sample. For a predictor with exact-zero handling,
#' its transformed values are set to zero at the original zeros *before*
#' this mean is subtracted. Thus, its centered values at zero need not be
#' zero. Binary columns have their lower observed value subtracted.
#' Formula-generated factor columns follow the appropriate binary or
#' mean-centering rule. Centering changes the intercept but not fitted values.
#'
#' The final fitted FP basis uses shifted, unscaled values; scaling is used
#' during selection. Supply new observations to [predict.mfp2()] on their
#' original predictor scale. Prediction applies the fitted transformation
#' and centering settings automatically.
#'
#' @section Categorical predictors and grouped terms:
#' In the formula interface, factors are expanded into contrast columns
#' automatically. Each factor is selected or removed as one term. To keep a
#' factor, give `keep` its name, not the names of its contrast columns:
#'
#' \preformatted{
#' mfp2(y ~ age + race, data = dat, keep = "race")
#' }
#'
#' Factors are not FP-transformed, and ACD and zero options do not apply
#' to them. Ordered factors use their configured contrasts, which are
#' polynomial by default. They are tested as a whole, rather than with a
#' single trend test.
#'
#' In the matrix interface, create the indicator or contrast columns yourself.
#' Use `term_groups` to select columns representing one categorical predictor
#' together. For example, if `x` has columns `age`, `raceB`, and `raceC`:
#'
#' \preformatted{
#' mfp2(x, y, term_groups = list(race = c("raceB", "raceC")),
#'      keep = "race")
#' }
#'
#' Here `keep = "race"` retains both `raceB` and `raceC`. With centering
#' enabled, formula-generated treatment indicators keep their zero reference;
#' other contrast columns are centered on their fitted-sample means. See
#' **Shifting, scaling, and centering**.
#'
#' @section ACD modelling:
#' The approximate cumulative distribution (ACD) transformation can model
#' a smooth S-shaped relationship between a continuous predictor and the
#' outcome (Royston, 2014; Royston and Sauerbrei, 2016). It uses smoothed
#' ranks to transform predictor values to a scale between 0 and 1.
#'
#' In the matrix interface, specify predictors in `acd_vars`. In a formula,
#' use `fp(x, acd = TRUE)`. Model selection considers ACD-based forms,
#' ordinary fractional polynomials, a linear effect, and exclusion of the
#' predictor. If the predictor has fewer than five distinct values, standard
#' MFP is used instead. ACD predictors may be shifted but are not scaled.
#' Selection uses the ACD function-selection procedure with
#' `criterion = "pvalue"`; with `"aic"` or `"bic"`, candidate models are
#' compared using the chosen criterion.
#'
#' When ACD is combined with `zero = TRUE`, the transformation is fitted
#' using positive values only; at least two distinct positive values are
#' required. Exact zeros receive a transformed value of zero before
#' centering. During prediction and plotting, positive values use the fitted
#' transformation, and exact zeros again receive zero before centering.
#'
#' See `vignette("mfp2_ACD", package = "mfp2")` for the ACD
#' function-selection procedure.
#'
#' @section Zero, catzero, and spike-at-zero modelling:
#' Three options are available for nonnegative predictors that contain
#' exact zeros.
#'
#' \itemize{
#'   \item `zero` applies the continuous transformation to positive values
#'     only. Its value at zero is set to zero before centering.
#'   \item `catzero` transforms positive values in the same way as `zero`
#'     and adds a binary variable that equals 1 when the predictor is zero
#'     and 0 otherwise.
#'   \item `spike` uses the spike-at-zero selection algorithm to decide whether
#'   the final model includes the positive-value function, the zero indicator,
#'   both, or neither.
#' }
#'
#' In the matrix interface, use `zero_vars`, `catzero_vars`, or `spike_vars`.
#' In a formula, use `fp(x, zero = TRUE)`, `fp(x, catzero = TRUE)`, or
#' `fp(x, spike = TRUE)`. `catzero` includes `zero` handling. An eligible
#' `spike` request includes both `catzero` and `zero` handling. If ACD is
#' also requested, its transformation uses positive values only.
#'
#' Negative values cause an error and are never treated as zeros. These
#' options do not apply to binary predictors. `zero` and `catzero` are
#' ignored if the predictor has no exact zeros.
#'
#' Spike-at-zero selection requires at least `min_saz_prop` of observations
#' with zero values and at least that proportion with positive values.
#' The default is `0.10`. If either group is too small, spike-at-zero
#' selection is skipped. Separately requested `zero` or `catzero` handling
#' remains in effect. With `criterion = "pvalue"`, spike-at-zero selection uses
#' a two-stage procedure controlled by `select` and `alpha`. With
#' `criterion = "aic"` or `"bic"`, it uses a one-stage procedure that
#' compares all candidate models directly. Further details on spike-at-zero
#' modelling are given by Kipruto and Sauerbrei (2026).
#'
#'
#' @section Subsetting:
#' When you use `subset`, `mfp2()` determines automatic shifts and scales
#' from the full input data before selecting observations. Model selection
#' and fitting then use only the selected observations. The result may
#' differ from a fit where the data were subset before calling `mfp2()`.
#'
#' In the formula interface, you can select rows using a column of `data`.
#' For example, if `dat` has a column called `group`, this fits the model
#' using only rows where `group` is `"A"`:
#'
#' \preformatted{
#' mfp2(y ~ age, data = dat, subset = group == "A")
#' }
#'
#' In the matrix interface, supply a logical vector or unique row numbers.
#' The package applies the subset to `x`, `y`, `weights`, `offset`, and
#' `strata` together. Row numbers retain the order you specify.
#'
#' Use `subset` when analyses of different groups should use shifts and
#' scales determined from the same full dataset. For cross-validation or
#' removal of missing values, prepare the analysis data before fitting.
#' The subset must leave more observations than predictor columns.
#'
#' @section Survival models:
#' Cox models require a right-censored [survival::Surv()] response, usually
#' `Surv(time, event)` with `0` for censoring and `1` for an event. For
#' stratification, use `strata()` in a formula or the `strata` argument with
#' matrix input. The separate `strata` argument in the formula interface is
#' deprecated. Cox settings such as `ties`, `nocenter`, and `control` follow
#' [survival::coxph()].
#'
#' Parametric survival models use [survreg_family()] and return a `survreg`
#' object. Here, `strata` allows separate scale parameters, as in
#' [survival::survreg()].
#'
#' Fine--Gray models use [finegray_family()] with a multi-state `Surv`
#' response. By default, `strata()` applies to both the censoring
#' distribution and the baseline subdistribution hazard; change this with
#' `strata_action`. Supply `id` for start–stop responses; it is optional
#' when each subject has one row.
#'
#' See [predict.mfp2()] for the data needed to predict survival probabilities.
#'
#' @section Multinomial models:
#' Use `family = "multinomial"` or [multinomial_family()] for an outcome
#' with at least three unordered classes. The response can contain class
#' labels or be a matrix of counts, with one column per class. The first
#' class is the reference; use `multinomial_family(reference = ...)` to
#' choose another.
#'
#' For each continuous predictor, the selected FP powers are shared across
#' classes, but each non-reference class has its own coefficients. Selection
#' degrees of freedom account for both the class-specific coefficients and
#' the search for shared powers. The final fit inherits from `multinom`.
#'
#' @section Ordinal models:
#' Use [ordinal_family()] when the outcome has at least three categories
#' with a meaningful order. An ordered factor specifies that order.
#' Otherwise, numeric values are ordered from smallest to largest, and
#' character values and unordered factor levels are ordered alphabetically.
#' For example, to set the order of named categories:
#'
#' \preformatted{
#' severity <- factor(
#'   c("low", "high", "medium", "low"),
#'   levels = c("low", "medium", "high"),
#'   ordered = TRUE
#' )
#' }
#'
#' The default `link = "logistic"` fits a proportional-odds model. It uses
#' the same fitted function for each predictor when comparing higher outcome
#' categories with lower ones. The final fit is an `orm` object and requires
#' the optional `rms` package. See [ordinal_family()] for other links and
#' how to interpret coefficient signs.
#'
#' @section Compatibility with the `mfp` package:
#' Both `mfp` and `mfp2` provide `fp()`. In formulas passed to `mfp2()`
#' or `mfpi()`, `fp()` always uses the `mfp2` version, even when both
#' packages are loaded. `fp()` and [fp2()] are equivalent in these formulas:
#'
#' \preformatted{
#' fit1 <- mfp2(y ~ fp(x1) + fp(x2), data = dat)
#' fit2 <- mfp2(y ~ fp2(x1) + fp2(x2), data = dat)
#' }
#'
#' @section Convergence and inference:
#' MFP selection usually converges within 2–4 cycles. With `verbose = TRUE`,
#' progress messages report when it converges. Lack of convergence involves
#' oscillation between two or more models and is extremely rare. If selection
#' does not converge within the allowed number of `cycles`, `mfp2()`
#' returns the model obtained after the final cycle. Consider
#' increasing `cycles`, changing `select` or `alpha`, or using
#' `criterion = "aic"` or `"bic"`.
#'
#' Standard errors, confidence intervals, and p-values from the final model
#' do not account for uncertainty from selecting predictors or their
#' functional forms, including ACD and spike-at-zero selection.
#'
#' @param x For `mfp2.default()`, a numeric matrix with one row per observation
#' and one column per predictor. Column names are required.
#' @param y Outcome for `mfp2.default()`. In the formula interface, put the
#' outcome on the left side of the formula. Its required format depends on
#' `family`; see **Model families and responses**.
#' @param term_groups For `mfp2.default()`, an optional named list of columns
#' to select and test together as one term. This is useful for dummy columns
#' representing a categorical predictor. For example,
#' `list(race = c("raceB", "raceC"))` treats both columns as the predictor
#' `race`. Grouped columns must be linear (`df = 1`) and cannot use ACD or
#' zero, catzero, or spike-at-zero modelling. Columns not listed are handled
#' individually.
#' @param formula For `mfp2.formula()`, a model formula. Use [fp()] or [fp2()]
#'   to set variable-specific FP, ACD, zero, catzero, or SAZ options.
#' @param data For `mfp2.formula()`, a data frame containing the variables in
#' `formula`.
#' @param weights Optional numeric observation weights. All weights must
#' be strictly positive. In the formula interface, an expression is evaluated
#' first in `data` and then in the formula environment.
#' @param offset A known term added to the model, such as the log of exposure
#' time. With formula input you can write `offset()` in the formula instead.
#' When predicting for new data, supply the matching offset.
#' @param subset Optional observations used for model selection and fitting.
#' Formula calls accept an expression evaluated in `data` and the formula
#' environment. Matrix calls accept a logical vector or unique numeric/integer
#' row positions. Automatic shift and scale values are estimated before the
#' subset is applied. See the Subsetting section.
#' @param cycles Maximum number of MFP backfitting cycles. Default `5`.
#' @param scale Divisor used when preparing predictors for FP selection.
#' `NULL` chooses scales automatically; `1` turns off automatic scaling.
#' In `mfp2.default()`, use a named value such as `c(age = 10)` to set
#' one column's scale. In formulas, set it inside `fp()` or `fp2()`.
#' Binary and ACD predictors always use scale `1`. See **Shifting,
#' scaling, and centering**.
#' @param shift Value added to predictors before FP transformation.
#' `NULL` chooses shifts automatically; `0` turns off automatic shifting.
#' In `mfp2.default()`, use a named value such as `c(age = 20)` to set
#' one column's shift. In formulas, set it inside `fp()` or `fp2()`.
#' Binary and linear predictors, and those using exact-zero handling,
#' always have shift `0`. See **Shifting, scaling, and centering**.
#' @param df Maximum FP complexity: `1` allows linear only, `2` allows up
#' to FP1, and `4` allows up to FP2. The default is `4`. Selection may
#' choose a simpler form or omit the predictor. In `mfp2.default()`, use
#' a named value such as `c(age = 2)` for one column; in formulas, use
#' `fp(age, df = 2)`. Predictors with few distinct values may be limited
#' to a simpler form. See **Fractional-polynomial model selection**.
#' @param center If `TRUE` (default), center predictors after any FP or ACD
#' transformation; `FALSE` leaves them uncentered. Continuous terms use
#' their mean and binary terms use their lower observed value in the data.
#' Centering does not change fitted values. In `mfp2.default()`,
#' use one logical value for all predictors or a named vector for selected
#' columns; unspecified columns use `TRUE`. In `mfp2.formula()`, use
#' `fp()` or `fp2()` to set this for an individual predictor.
#' @param family Model family; default `"gaussian"`. Supported GLM families
#' are `"gaussian"`, `"binomial"`, `"poisson"`, `"Gamma"`, and
#' `"inverse.gaussian"`, supplied as a name, function, or family object.
#' Quasi families are not supported. Use `"negbin"` for negative-binomial
#' models; `"multinomial"` or [multinomial_family()] for multinomial models;
#' `"ordinal"` or [ordinal_family()] for cumulative-link ordinal models;
#' `"cox"` for Cox models; [survreg_family()] for parametric survival
#' models; [finegray_family()] for Fine--Gray models; or [gee_family()]
#' for generalized estimating equations. GEE supports Gaussian,
#' binomial, Poisson, and Gamma response families.
#' @param fitter Backend for GLM candidate fits during the FP search.
#' `"base"` (default) uses [stats::glm.fit()]; `"fastglm"` uses the
#' optional `fastglm` package. If `fastglm` is unavailable, ordinary GLMs
#' use the base backend with a warning. Negative-binomial models
#' always use `fastglm` and require it to be installed. Survival,
#' multinomial, ordinal, and GEE models use their own fitters and ignore
#' this setting.
#' @param criterion Model-selection criterion: `"pvalue"` (default),
#' `"aic"`, or `"bic"`. With `"pvalue"`, `select` controls predictor
#' inclusion and `alpha` controls FP functional form and spike-at-zero component
#' comparisons. With `"aic"` or `"bic"`, the model with the lowest score is
#' selected. For GEE models, `"aic"` and `"bic"` use quasi likelihood
#' based scores; see **Model selection criteria** for their definitions.
#' @param select Significance level for predictor inclusion when
#' `criterion = "pvalue"`; default `0.05`. In `mfp2.default()`, a single
#' value applies to all predictors. A named vector sets the selection level
#' for the specified predictors; all others use the default of `0.05`. In
#' `mfp2.formula()`, the top-level value applies only to terms outside
#' `fp()` or `fp2()`; those terms each default to `0.05` unless set
#' explicitly, for example `fp(x, select = 1)`. Under p-value selection,
#' `select = 1` forces inclusion. For spike-at-zero terms, `select`
#' controls initial inclusion; `alpha` controls subsequent component
#' removal. `select` does not affect AIC or BIC selection; use `keep`
#' to retain a term under those criteria.
#' @param alpha Significance level for choosing between FP forms when
#' `criterion = "pvalue"`; default `0.05`. A simpler form is chosen only
#' when its comparison p-value is greater than `alpha`, so `alpha = 1`
#' prevents this simplification. For spike-at-zero terms, `alpha` also
#' controls tests that remove either component. In `mfp2.default()`,
#' a single value applies to all predictors. A named vector sets values
#' for specified predictors; all others use `0.05`. In
#' `mfp2.formula()`, the top-level value applies only to terms outside
#' `fp()` or `fp2()`. Set `alpha` within `fp()` or `fp2()` to change it
#' for those terms. `alpha` does not affect AIC or BIC selection.
#' @param keep Character vector naming predictors to retain in the final
#' model under any selection criterion. For a formula interface, use a
#' predictor name  from the formula, such as `keep = "age"`. For a matrix
#' interface, use a column name from `x` or a name in `term_groups`. Keeping a
#' predictor does not prevent selection of a simpler functional form.
#' @param xorder Order in which predictors are examined during backfitting:
#' `"ascending"` (default), `"descending"`, or `"original"`.
#' To rank predictors, each is tested as a whole by removing it from
#' the full starting model. `"ascending"` visits the smallest
#' p-value first, `"descending"` the largest, and `"original"` follows
#' the input term order.
#' @param powers Optional named list of candidate FP powers for individual
#' predictors, for example `list(age = c(0, 1, 2))`. Predictors not named
#' use `c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)`, where `0` denotes a
#' logarithm. Each list element must contain at least one finite numeric
#' power; duplicate powers are removed. Use the predictor's column name
#' as the list name. To allow only a linear effect, set `df = 1`. Powers
#' supplied inside `fp()` or `fp2()` take precedence over this list.
#' @param ties Method for handling tied event times in Cox and Fine-Gray models.
#' Supported values are `"breslow"` (default) and `"efron"`. `"exact"` is not
#' supported. The argument has no effect for other families. See
#' [survival::coxph()] for details.
#' @param strata Optional stratification for Cox, parametric survival, and
#' Fine-Gray models. In `mfp2.default()`, supply one stratum label per
#' observation as a vector or factor. To stratify by several variables,
#' supply a matrix or data frame with one variable per column. In
#' `mfp2.formula()`, put `strata(...)` in the formula. The separate
#' `strata` argument is deprecated for formula fits; if both are supplied,
#' the formula term takes precedence.
#' @param id A subject or cluster identifier, with one value per observation.
#' Use the same ID for observations from the same subject or cluster. IDs
#' cannot be missing. Fine--Gray models require `id` for start--stop
#' multi-state data; it is optional when there is one row per subject. GEE
#' models always require `id` and group observations by cluster before
#' fitting, preserving their order within each cluster. In formula fits,
#' `id` is looked up first in `data`, then in the formula environment. It
#' is not used by other model families.
#' @param waves Optional index giving the order of repeated observations
#' within each `id` cluster in GEE models. Supply one positive integer per
#' observation; values must be unique within each cluster. `waves` does not
#' reorder rows. See [geepack::geeglm()] for further details. In formula
#' fits, `waves` is looked up first in `data`, then in the formula environment.
#' @param nocenter Numeric values that exempt a Cox or Fine--Gray predictor
#' column from centering when every value in that column is among them.
#' The default, `c(-1, 0, 1)`, leaves indicator columns uncentered, as in
#' [survival::coxph()]. Set `NULL` to allow all predictor columns to be
#' centered. Ignored for other model families.
#' @param acd_vars Character vector naming continuous columns of `x` to
#' consider for approximate cumulative distribution (ACD) modelling.
#' ACD functions are compared with simpler alternatives and need not be
#' selected. Applies only to `mfp2.default()`; in formulas, use
#' `fp(x, acd = TRUE)` or `fp2(x, acd = TRUE)`. See **ACD modelling**.
#' @param acdx Deprecated compatibility alias for `acd_vars` in
#' `mfp2.default()`. Use `acd_vars` in matrix interface.
#' @param ftest Logical. If `TRUE`, use F-tests for selection in ordinary
#' Gaussian regression models when `criterion = "pvalue"`. These tests cover
#' variable selection, functional-form selection, and spike-at-zero testing,
#' and may improve small-sample behaviour. Default `FALSE`. Not applicable
#' to GEE or other model families.
#' @param control Optional named list of model-fitting settings. If `NULL`,
#' the package uses its default settings for the selected model family.
#' The available settings depend on the model family: see
#' [stats::glm.control()] for ordinary GLMs, [survival::coxph.control()] for
#' Cox and Fine-Gray models, [survival::survreg.control()] for parametric
#' survival models, and [geepack::geese.control()] for GEE models.
#' Multinomial models accept `maxit`, `reltol`, `abstol`, `trace`, and
#' `MaxNWts`; ordinal models accept `maxit`, `eps`, `tol`, and `trace`.
#' @param zero_vars Names of nonnegative continuous predictors with exact zeros.
#' Only positive values are transformed; zeros remain zero before centering.
#' If a predictor is also in `acd_vars`, ACD is fitted using only its positive
#' values. The shift is fixed at zero. Applies to `mfp2.default()` only; in
#' formulas, use `fp(x, zero = TRUE)`. Negative values are rejected. See zero
#' and spike-at-zero modelling section.
#' @param catzero_vars Names of nonnegative continuous predictors with exact
#' zeros. Positive values are transformed, and a separate indicator for zero
#' values is added to the model. This also enables `zero_vars` behavior for
#' each predictor. Applies to `mfp2.default()` only; in formulas, use
#' `fp(x, catzero = TRUE)`. Negative values are rejected. See the zero
#' and spike-at-zero modelling section.
#' @param spike_vars Names of nonnegative continuous predictors to assess with
#' spike-at-zero selection. The model may retain the transformed positive
#' values, an indicator for zeros, or both. Selection uses two-stage tests
#' when `criterion = "pvalue"` and compares candidate models when using AIC
#' or BIC. Each group must meet `min_saz_prop`; binary predictors are
#' ineligible. Ineligible requests are ignored. Applies to `mfp2.default()`
#' only; in formulas, use `fp(x, spike = TRUE)`.
#' @param min_saz_prop Minimum proportion of observations required in each
#' group (`x = 0` and `x > 0`) for spike-at-zero selection. Must be between
#' 0 and 0.5; default `0.10`. Predictors that do not meet this threshold
#' are not assessed. In large samples, a lower threshold may still leave
#' enough observations in each group.
#' @param force_max_fp_vars Names of predictors whose most complex permitted
#' form should be kept in the model, regardless of the selection criterion.
#' The maximum is set by `df`; when `df = 1`, it is a linear form. For an
#' eligible predictor in `spike_vars`, both the zero indicator and the most
#' complex permitted function of positive values are kept. Default `NULL`.
#' Applies to `mfp2.default()` only; in formulas, use
#' `fp(x, force_max_fp = TRUE)`.
#' @param verbose Logical. If `TRUE`, print progress messages during model
#' fitting. Default `TRUE`.
#' @param \dots Used to pass arguments from `mfp2()` to its fitting method.
#'   Additional arguments supplied to the fitting method are not used.
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
#'   lpsa ~ fp(cavol, acd = TRUE) + fp(age) + svi,
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
#'
#' # Examples for every supported model family -----------------------------
#' set.seed(2026)
#' n_family <- 120
#' family_data <- data.frame(
#'   x1 = runif(n_family, 0.5, 3),
#'   x2 = rnorm(n_family)
#' )
#' eta <- with(family_data, 0.2 + 0.3 * x1 - 0.2 * x2)
#' mu <- exp(eta)
#' simulate_inverse_gaussian <- function(mu, shape) {
#'   squared_normal <- stats::rnorm(length(mu))^2
#'   candidate <- mu + (mu^2 * squared_normal) / (2 * shape) -
#'     (mu / (2 * shape)) * sqrt(
#'       4 * mu * shape * squared_normal + mu^2 * squared_normal^2
#'     )
#'   choose_candidate <- stats::runif(length(mu)) <= mu / (mu + candidate)
#'   ifelse(choose_candidate, candidate, mu^2 / candidate)
#' }
#' family_data$y_gaussian <- eta + rnorm(n_family)
#' family_data$y_binomial <- rbinom(n_family, 1, plogis(eta - 1))
#' family_data$y_poisson <- rpois(n_family, mu)
#' family_data$y_gamma <- rgamma(n_family, shape = 5, scale = mu / 5)
#' family_data$y_inverse_gaussian <- simulate_inverse_gaussian(mu, shape = 8)
#' family_data$y_negbin <- rnbinom(n_family, mu = mu, size = 2)
#' family_data$y_multinomial <- sample.int(
#'   3L, n_family, replace = TRUE, prob = c(0.45, 0.35, 0.20)
#' )
#' family_data$y_ordinal <- ordered(
#'   cut(
#'     eta + stats::rlogis(n_family),
#'     breaks = c(-Inf, -0.5, 0.5, 1.5, Inf),
#'     labels = c("low", "medium", "high", "very high")
#'   )
#' )
#'
#' fit_gaussian <- mfp2(
#'   y_gaussian ~ x1 + x2, data = family_data, family = "gaussian",
#'   df = 1, select = 1, alpha = 1, cycles = 1, verbose = FALSE
#' )
#' fit_binomial <- mfp2(
#'   y_binomial ~ x1 + x2, data = family_data, family = "binomial",
#'   df = 1, select = 1, alpha = 1, cycles = 1, verbose = FALSE
#' )
#' fit_poisson <- mfp2(
#'   y_poisson ~ x1 + x2, data = family_data, family = "poisson",
#'   df = 1, select = 1, alpha = 1, cycles = 1, verbose = FALSE
#' )
#' fit_gamma <- mfp2(
#'   y_gamma ~ x1 + x2, data = family_data,
#'   family = stats::Gamma(link = "log"),
#'   df = 1, select = 1, alpha = 1, cycles = 1, verbose = FALSE
#' )
#' fit_inverse_gaussian <- mfp2(
#'   y_inverse_gaussian ~ x1 + x2, data = family_data,
#'   family = stats::inverse.gaussian(link = "log"),
#'   df = 1, select = 1, alpha = 1, cycles = 1, verbose = FALSE
#' )
#'
#' # Negative binomial uses the optional fastglm backend.
#' if (requireNamespace("fastglm", quietly = TRUE)) {
#'   fit_negbin <- mfp2(
#'     y_negbin ~ x1 + x2, data = family_data,
#'     family = "negbin",
#'     df = 1, select = 1, alpha = 1, cycles = 1, verbose = FALSE
#'   )
#' }
#'
#' # Integer class labels are converted to a multinomial factor response.
#' fit_multinomial <- mfp2(
#'   y_multinomial ~ x1 + x2, data = family_data,
#'   family = multinomial_family(),
#'   df = 1, select = 1, alpha = 1, cycles = 1, verbose = FALSE
#' )
#'
#' # Ordinal regression uses the optional rms package.
#' if (requireNamespace("rms", quietly = TRUE)) {
#'   fit_ordinal <- mfp2(
#'     y_ordinal ~ x1 + x2, data = family_data,
#'     family = ordinal_family(),
#'     df = 1, select = 1, alpha = 1, cycles = 1, verbose = FALSE
#'   )
#' }
#'
#' # Cox proportional hazards and parametric survival models.
#' event_time <- rexp(n_family, rate = exp(eta - 2))
#' censor_time <- rexp(n_family, rate = 0.15)
#' survival_data <- transform(
#'   family_data,
#'   time = pmin(event_time, censor_time),
#'   status = as.integer(event_time <= censor_time)
#' )
#' fit_cox <- mfp2(
#'   survival::Surv(time, status) ~ x1 + x2,
#'   data = survival_data, family = "cox",
#'   df = 1, select = 1, alpha = 1, cycles = 1, verbose = FALSE
#' )
#' fit_survreg <- mfp2(
#'   survival::Surv(time, status) ~ x1 + x2,
#'   data = survival_data, family = survreg_family(dist = "weibull"),
#'   df = 1, select = 1, alpha = 1, cycles = 1, verbose = FALSE
#' )
#'
#' # Fine--Gray subdistribution hazards for a competing-risks response.
#' finegray_data <- transform(
#'   family_data,
#'   ftime = rexp(n_family, rate = 0.12),
#'   event = factor(
#'     sample(
#'       c("censor", "cause1", "cause2"), n_family, replace = TRUE,
#'       prob = c(0.35, 0.40, 0.25)
#'     ),
#'     levels = c("censor", "cause1", "cause2")
#'   )
#' )
#' fit_finegray <- mfp2(
#'   survival::Surv(ftime, event) ~ x1 + x2,
#'   data = finegray_data, family = finegray_family(etype = "cause1"),
#'   df = 1, select = 1, alpha = 1, cycles = 1, verbose = FALSE
#' )
#' }
#'
#' @return
#' An object of class `mfp2`, built from the final fitted model and extended
#' with MFP selection and transformation information. Depending on the model
#' family, it also inherits from:
#'
#' \itemize{
#'   \item `glm` for Gaussian, binomial, Poisson, Gamma, and inverse-Gaussian models;
#'   \item `fastglm_nb` and `fastglm` for negative-binomial models;
#'   \item `geeglm` for GEE models;
#'   \item `multinom` and `nnet` for multinomial logistic models;
#'   \item `orm` for ordinal regression models;
#'   \item `coxph` for Cox and Fine--Gray models;
#'   \item `survreg` for parametric survival models.
#' }
#'
#' Standard fitted-model components, such as coefficients, residuals, fitted
#' values, and the model call, remain available. `mfp2()` adds the following
#' components to describe its selection results:
#'
#' \describe{
#'   \item{convergence_mfp}{
#'     Logical value indicating whether the MFP selection algorithm converged.
#'   }
#'
#'   \item{fp_terms}{
#'     A data frame with one row per model term. It records the MFP complexity
#'     setting (`df_setting`) separately from the initial and final model degrees
#'     of freedom (`df_initial` and `df_final`), together with selection settings,
#'     selected status, searched-FP identity (`searched_fp`),
#'     fractional-polynomial powers, the ACD, zero, catzero, and spike settings
#'     actually applied after validation and eligibility checks, and `prop_zero`
#'     for retained spike-at-zero terms.
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
#'     indicator for exact-zero values.
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
#' Only documented components form part of the stable public interface. Other
#' components may change between package versions.
#'
#' @references
#' Royston, P. and Sauerbrei, W., 2008. \emph{Multivariable Model - Building:
#' A Pragmatic Approach to Regression Analysis based on Fractional Polynomials
#' for Modelling Continuous Variables. John Wiley & Sons.}\cr
#'
#' Liang, K.-Y. and Zeger, S. L., 1986. \emph{Longitudinal data analysis using
#' generalized linear models. Biometrika, 73(1): 13--22.}\cr
#'
#' Pan, W., 2001. \emph{Akaike's information criterion in generalized
#' estimating equations. Biometrics, 57(1): 120--125.}\cr
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
#' [plot.mfp2()], [get_selected_variable_names()], [transform_vector_fp()],
#' [gee_family()], [multinomial_family()], [ordinal_family()],
#' [survreg_family()], [finegray_family()]
#'
#' @export
mfp2 <- function(x, ...) {
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
#' @return Named list of metadata for simple factor main effects. Each entry
#'   records the source factor variable, its levels and ordering, the fitted raw
#'   design columns, and the exact design row for every observed level. More
#'   complex factor-containing terms, such as interactions, are omitted because
#'   their design rows also depend on the interacting variables.
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
      variable = frame_variable,
      levels = level_names,
      ordered = is.ordered(values),
      columns = columns,
      design_by_level = design_by_level
    )
  }

  out
}


#' Determine Centering Methods for Formula-Generated Factor Columns
#'
#' Formula processing knows which numeric design columns came from categorical
#' terms, information that is unavailable once only a numeric matrix remains.
#' Treatment-style 0/1 columns retain minimum centering, while other contrast
#' columns use their fitted-sample means even when a particular contrast happens
#' to contain only two distinct numeric values.
#'
#' @param factor_terms Character vector of conceptual factor-term names.
#' @param term_to_columns Conceptual term-to-design-column lookup.
#' @param x Numeric fitting design matrix after intercept removal and any
#'   \code{fp()} column renaming.
#'
#' @return A named character vector containing \code{"mean"} or
#'   \code{"minimum"} for factor-generated columns. Returns \code{NULL} when
#'   the formula contains no factor design columns.
#' @keywords internal
#' @noRd
formula_factor_center_methods <- function(factor_terms, term_to_columns, x) {
  factor_columns <- unique(unname(unlist(
    term_to_columns[intersect(factor_terms, names(term_to_columns))],
    use.names = FALSE
  )))
  factor_columns <- intersect(factor_columns, colnames(x))

  if (length(factor_columns) == 0L) {
    return(NULL)
  }

  methods <- vapply(factor_columns, function(column) {
    values <- unique(x[, column])
    is_indicator <- length(values) <= 2L && all(values %in% c(0, 1))

    if (is_indicator) "minimum" else "mean"
  }, character(1L))

  stats::setNames(unname(methods), factor_columns)
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
                         ties = c("breslow", "efron"),
                         strata = NULL,
                         nocenter = c(-1, 0, 1),
                         acd_vars = NULL,
                         ftest = FALSE,
                         control = NULL,
                         zero_vars = NULL,
                         catzero_vars = NULL,
                         spike_vars = NULL,
                         min_saz_prop = 0.10,
                         force_max_fp_vars = NULL,
                         term_groups = NULL,
                         verbose = TRUE,
                         fitter = c("base", "fastglm"),
                         id = NULL,
                         waves = NULL,
                         acdx = NULL,
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

  # Normalize the deprecated public spelling at the API boundary. Internal
  # ACD code continues to use `acdx` as its named logical-vector representation.
  acd_vars_supplied <- !missing(acd_vars)
  acdx_supplied <- !missing(acdx)

  if (acd_vars_supplied && acdx_supplied) {
    stop(
      "`acd_vars` and its deprecated alias `acdx` cannot both be supplied.",
      call. = FALSE
    )
  }

  if (acdx_supplied) {
    .Deprecated(
      new = "acd_vars",
      package = "mfp2",
      old = "acdx",
      msg = "`acdx` is deprecated; use `acd_vars` instead."
    )
    acd_vars <- acdx
  }

  acdx <- acd_vars

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
  ties      <- resolve_mfp_ties(ties)
  fitter <- match.arg(fitter)

  # Step 3: Validate family and the input matrix `x` ----------------------------
  family_info <- normalize_family_argument(
    family,
    family_arg = deparse(substitute(family))
  )

  family        <- family_info$family
  family_string <- family_info$family_string
  fitter <- resolve_fitter(fitter, family_string)

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
  center_method <- attr(x, "mfp2_center_method", exact = TRUE)
  attr(x, "mfp2_center_method") <- NULL
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

  # Predictor names are used throughout MFP to match variable-specific
  # settings and to construct the single final formula-based refit. Validate
  # them once here, outside every candidate-fitting loop. Backticks are rejected
  # because the existing formula builder uses backticks to quote non-syntactic
  # user names; spaces, hyphens, and other non-syntactic names remain valid.
  vnames <- colnames(x)
  validate_predictor_names(vnames, object = "`x`")
  center_method <- normalize_center_method(center_method, vnames)

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
  }

  # Step 5: Validate `weights` and `offset` --------------------------------------
  # Use one strict observation-weight contract for every family. Zero weights
  # are rejected even if a fitter would otherwise accept them because they can
  # make likelihood-based MFP comparisons undefined (for example, Gaussian
  # glm() can return an infinite AIC/log-likelihood with zero prior weights).
  validate_model_weights(
    weights = weights,
    nobs = nobs
  )

  # Ordinary families use one offset per row. Multinomial models accept an
  # n x C class-offset matrix or n x (C - 1) reference-logit matrix.
  if (!is.null(offset)) {
    if (!is.numeric(offset)) {
      stop(
        "! `offset` must be numeric.",
        sprintf("i Current type is: %s.", typeof(offset)),
        call. = FALSE
      )
    }

    offset_rows <- if (is.matrix(offset)) nrow(offset) else length(offset)
    if (offset_rows != nobs) {
      stop(
        "! The number of observations in x and offset must match.",
        sprintf(
          "i The number of rows in x is %d, but the offset has %d observation rows.",
          nobs, offset_rows
        ),
        call. = FALSE
      )
    }

    if (identical(family_string, "multinomial")) {
      c_classes <- mfp2_multinomial_n_classes(y)
      if (!is.matrix(offset) || !ncol(offset) %in% c(c_classes - 1L, c_classes)) {
        stop(
          "! Multinomial `offset` must be an n x C class matrix or an n x (C - 1) reference-logit matrix.",
          call. = FALSE
        )
      }
    } else if (is.matrix(offset)) {
      stop("! `offset` must be a numeric vector for this family.", call. = FALSE)
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

  # These switches are scalar logical flags.
  validate_logical_vector(verbose, "verbose", allowed_lengths = 1L)
  validate_logical_vector(ftest, "ftest", allowed_lengths = 1L)

  # alpha and select accept an unnamed global scalar or named partial
  # overrides. Omitted named entries use the public defaults and all matching is
  # by column name rather than by position.
  alpha_setting <- normalize_named_override_setting(
    value = alpha,
    column_names = vnames,
    default = 0.05,
    argument_name = "alpha"
  )
  alpha <- alpha_setting$value

  select_setting <- normalize_named_override_setting(
    value = select,
    column_names = vnames,
    default = 0.05,
    argument_name = "select"
  )
  select <- select_setting$value

  validate_probability_vector(alpha, "alpha", nvars)
  validate_probability_vector(select, "select", nvars)

  # center follows the same matrix-interface convention as df/select/alpha:
  # an unnamed scalar is global, while named partial vectors are matched to
  # colnames(x) and omitted columns retain the public default (TRUE).
  center <- normalize_named_logical_setting(
    value = center,
    column_names = vnames,
    default = TRUE,
    argument_name = "center"
  )

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
  validate_variable_names(acdx, "acd_vars", vnames)
  validate_variable_names(spike, "spike_vars", vnames)
  validate_variable_names(force_max_fp,"force_max_fp_vars", vnames)

  # Enforce the exact-zero domain on the original, unshifted preprocessing
  # values before any option cascade, SAZ eligibility check, or shift/scale
  # calculation. A variable requested through any of the three options must be
  # nonnegative even if a later eligibility rule would otherwise reset it.
  validate_zero_option_covariates(
    x = preprocess_x,
    variables = unique(c(zero, catzero, spike))
  )

  # Step 8: Validate `df` ---------------------------------------------------------
  # df accepts an unnamed global scalar or named partial overrides. Ordinary
  # omitted columns use df = 4; explicitly mapped design blocks use their
  # structural linear default df = 1 unless the caller overrides them.
  df_fallback <- expand_scalar_df_for_mapped_terms(
    df = 4L,
    vnames = vnames,
    term_to_columns = term_to_columns
  )
  df_setting <- normalize_named_override_setting(
    value = df,
    column_names = vnames,
    default = df_fallback,
    argument_name = "df"
  )
  df <- df_setting$value
  df_supplied <- df_setting$supplied
  df_global_scalar <- df_setting$global_scalar
  df_scalar_value <- if (df_global_scalar && length(df) > 0L) {
    unname(df[[1L]])
  } else {
    NULL
  }

  # df must translate into a valid FP degree: 1 (linear) or an even positive
  # integer (2*degree) 2, 4, 6, ... (representing FP1, FP2, FP3, ...).
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
  # GEE validates the response against its inner GLM response family; every
  # other family validates against its own canonical name.
  validate_family_response(
    y = y,
    family_string = if (identical(family_string, "gee")) {
      family$response_family$family
    } else {
      family_string
    },
    nobs = nobs
  )

  # Validate survival-model auxiliary inputs -----------------------------------
  # validate_family_response() checks the response itself, but strata is an
  # auxiliary argument and therefore still needs to be checked here. For Cox it
  # defines baseline-hazard strata, for survreg it defines scale strata, and for
  # Fine--Gray it defines censoring-distribution strata.
  if (mfp2_family_is_survival(family_string) && !is.null(strata)) {
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
  # `id` is a finegray and GEE-only auxiliary parameter
  if (!is.null(id)) {
    if (!family_string %in% c("finegray", "gee")) {
      stop(
        "! `id` is only used with `family = finegray_family()` or `family = gee_family()`.",
        call. = FALSE
      )
    }
    if (!is.atomic(id) || length(id) != nobs || anyNA(id) ||
        is.matrix(id) || !is.null(dim(id))) {
      stop(
        "! `id` must be an atomic vector or factor with one non-missing value per observation.",
        call. = FALSE
      )
    }
  } else if (identical(family_string, "gee")) {
    stop(
      "! `id` is required for `family = gee_family()`: supply a cluster ",
      "identifier with one value per observation.",
      call. = FALSE
    )
  }

  # `waves` is a GEE-only auxiliary aligned per observation. Full structural
  # validation (integer-like, unique within cluster) happens in
  # prepare_gee_family(); here only the shape and family are checked.
  if (!is.null(waves)) {
    if (!identical(family_string, "gee")) {
      stop("! `waves` is only used with `family = gee_family()`.", call. = FALSE)
    }
    if (!is.atomic(waves) || length(waves) != nobs || anyNA(waves) ||
        is.matrix(waves) || !is.null(dim(waves))) {
      stop(
        "! `waves` must be an atomic vector with one non-missing value per observation.",
        call. = FALSE
      )
    }
  }

  # Step 10: Validate spike-at-zero proportion and the candidate power list -----
  # min_saz_prop is the minimum share of observations required in
  # both the zero component and the positive component for a variable
  # requested via spike_vars to remain SAZ-eligible (see resolve_saz_eligibility()
  # below).
  if (!is.numeric(min_saz_prop) ||
      length(min_saz_prop) != 1L ||
      anyNA(min_saz_prop) ||
      !is.finite(min_saz_prop) ||
      min_saz_prop <= 0 ||
      min_saz_prop >= 0.5) {
    stop(
      "! `min_saz_prop` must be a single finite numeric value in the open interval (0, 0.5).",
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

  # select and alpha were normalized above to complete named vectors.

  # `shift` was normalized above. NA marks a variable for automatic shift
  # estimation further below; supplied values are already aligned by name.

  # center was normalized above to a complete named vector in colnames(x) order.

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

  # Resolve defaults, validate partial lists, and check fitter compatibility
  # once before the repeated candidate-fit path begins.
  control <- normalize_fit_control(
    control = control,
    family_string = family_string,
    fitter = fitter
  )

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

  # zero_vars only makes sense for variables that actually contain exact-zero
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

  # acd_vars: character vector of variable names -> named logical vector over vnames.
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
    term_to_columns, acdx, "acd_vars", function(v) !v, "FALSE"
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

  # Resolve ACD eligibility on the actual fitting rows before scale factors are
  # chosen. Matrix-interface subsets are applied later so that automatic
  # preprocessing can use the full data, but ACD eligibility itself must reflect
  # the observations entering the model. Formula-interface calls already supply
  # the subset-specific fitting matrix in `x`.
  if (any(acdx)) {
    acd_eligibility_x <- if (is.null(subset)) {
      x
    } else {
      x[subset, , drop = FALSE]
    }
    acdx <- reset_acd(acd_eligibility_x, acdx)
  }

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

  # ACD uses the shifted covariate directly, without the ordinary MFP scaling
  # step. Reset active ACD variables to scale 1 before missing scales are
  # estimated and before the design matrix is divided below. Requests reset by
  # reset_acd() remain ordinary FP variables and retain their usual scale rule.
  acd_names <- names(acdx)[acdx]
  if (length(acd_names) > 0L) {
    scale[acd_names] <- 1
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
      min_saz_prop = min_saz_prop
    )

    spike   <- saz_flags$spike
    catzero <- saz_flags$catzero
    zero    <- saz_flags$zero
  }

  # Step 14: Assign effective df per predictor -----------------------------------
  # An unnamed scalar retains the historical scalar-default cardinality rules.
  # For named inputs, explicitly supplied values retain the historical vector
  # behavior, while omitted values filled from the package default receive the
  # ordinary scalar-default cardinality reductions.
  if (df_global_scalar) {
    # A scalar df controls ordinary continuous predictors. Explicitly mapped
    # terms are supplied fixed linear design blocks, so assign df = 1 to their
    # raw columns before applying cardinality rules to the remaining columns.
    df_default <- expand_scalar_df_for_mapped_terms(
      df = df_scalar_value,
      vnames = vnames,
      term_to_columns = term_to_columns
    )

    if (any(df_default != 1L)) {
      df.list <- assign_df(x = preprocess_x, df_default = df_default)
    } else {
      df.list <- df_default
    }
  } else {
    defaulted <- !df_supplied
    if (any(defaulted)) {
      df_from_defaults <- assign_df(x = preprocess_x, df_default = df)
      df[defaulted] <- df_from_defaults[defaulted]
    }

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
  if (formula_generated_settings) {
    # Formula calls may carry a full-data preprocessing matrix that differs
    # from the fitting design (for example when `subset` is used). Shift
    # both matrices independently so full-data preprocessing semantics are
    # preserved exactly.
    preprocess_x_shifted <- sweep(preprocess_x, 2, shift, "+")
    x <- sweep(x, 2, shift, "+")
  } else {
    # Direct matrix calls have no separate preprocessing matrix:
    # extract_preprocess_matrix() returns the same input matrix for both
    # `x` and `preprocess_x`. Shift that matrix once and let the
    # shifted-but-unscaled result serve both roles. The later scaling
    # assignment rebinds `x` to a new matrix, so `preprocess_x_shifted`
    # remains available for scale estimation and positivity checks.
    # This removes one redundant O(n * p) sweep and one full temporary
    # matrix allocation without changing the numerical preprocessing.
    x <- sweep(x, 2, shift, "+")
    preprocess_x_shifted <- x
  }

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
  # positivity check above. Active ACD variables have scale 1 and therefore pass
  # through unchanged on their shifted scale.
  x <- sweep(x, 2, scale, "/")

  # Step 18: Build the Cox stratification object --------------------------------
  # Normalize all Cox strata to one factor before any candidate fit. This
  # mirrors survival::coxph(): multiple strata variables are combined first,
  # while a single character/numeric/logical vector is treated as categorical
  # labels. Integer conversion is deferred to the coxph.fit() boundary.
  strata_keep <- strata

  if (mfp2_family_is_survival(family_string) && !is.null(strata_keep)) {
    strata_keep <- normalize_cox_strata(strata_keep, nobs = nobs)
  }

  id_keep <- id
  waves_keep <- waves

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
    offset <- if (is.matrix(offset)) {
      offset[subset, , drop = FALSE]
    } else {
      offset[subset]
    }
    if (!is.null(strata_keep)) {
      # Keep strata aligned with the fitted rows and remove levels no longer
      # represented after subsetting.
      strata_keep <- droplevels(strata_keep[subset])
    }
    if (!is.null(id_keep)) id_keep <- id_keep[subset]
    if (!is.null(waves_keep)) waves_keep <- waves_keep[subset]
  }

  # geepack treats a contiguous run of equal IDs as one cluster. Group the
  # fitted rows once for all candidate models and the retained geeglm fit.
  # This stable order preserves within-cluster visit order when waves is absent
  # and keeps x, response, weights, offset, id and waves aligned. It is relative
  # to the already-subset fitted sample, not to the original full data.
  gee_input_row_order <- NULL
  if (identical(family_string, "gee")) {
    gee_input_row_order <- mfp2_gee_cluster_order(id_keep)
    if (!identical(gee_input_row_order, seq_len(NROW(y)))) {
      x <- x[gee_input_row_order, , drop = FALSE]
      y <- if (is.matrix(y)) {
        y[gee_input_row_order, , drop = FALSE]
      } else {
        y[gee_input_row_order]
      }
      if (!is.null(weights)) weights <- weights[gee_input_row_order]
      if (!is.null(offset)) offset <- offset[gee_input_row_order]
      id_keep <- id_keep[gee_input_row_order]
      if (!is.null(waves_keep)) waves_keep <- waves_keep[gee_input_row_order]
    }
  }

  # Require more fitted observations than predictor columns. For matrix calls,
  # this check runs after `subset` has been applied; formula calls already pass
  # their subset-specific model matrix into this shared path.
  dimx <- dim(x)
  nrow_x <- dimx[1L]
  ncol_x <- dimx[2L]
  if (nrow_x <= ncol_x) {
    stop(
      sprintf(
        "! The number of observations must be greater than the number of predictor columns (n = %d, p = %d).",
        nrow_x,
        ncol_x
      ),
      call. = FALSE
    )
  }

  # Collapse raw-column options to one value per conceptual term only after all
  # column-level preprocessing has completed. Grouped members must agree on
  # settings that have a single term-level meaning.
  alpha_term <- collapse_option_to_terms(
    alpha, term_to_columns, "alpha"
  )
  select_term <- collapse_option_to_terms(
    select, term_to_columns, "select"
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
  acdx_term <- collapse_option_to_terms(acdx, term_to_columns, "acd_vars")
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
  prepared_family <- prepare_family_for_fit(
    family = family,
    family_string = family_string,
    y = y,
    weights = weights,
    strata = strata_keep,
    id = id_keep,
    offset = offset,
    has_offset = has_offset,
    waves = waves_keep,
    control = control
  )
  family <- prepared_family$family
  strata_keep <- prepared_family$strata

  validate_default_design_rank(
    x = x,
    intercept = mfp2_family_has_intercept(family_string)
  )

  # Delegate the actual variable selection, FP degree/power selection, and
  # model fitting (including SAZ and ACD handling) to fit_mfp().
  fit <- fit_mfp(
    x = x, y = y,
    weights = weights, offset = offset, cycles = cycles,
    scale = scale_term, shift = shift_term, df = df_term, center = center_term,
    family = family, family_string = family_string, fitter = fitter,
    criterion = criterion,
    select = select_term, alpha = alpha_term, keep = keep_term, xorder = xorder,
    powers = powers_term, method = ties, strata = strata_keep, nocenter = nocenter,
    acdx = acdx_term, ftest = ftest, force_max_fp = force_max_fp_term,
    control = control, zero = zero_term, catzero = catzero_term, spike = spike_term,
    min_saz_prop = min_saz_prop,
    saz_pre_resolved = TRUE,
    term_to_columns = term_to_columns,
    center_method = center_method,
    has_offset = has_offset,
    verbose = verbose
  )

  # Step 21: Attach mfp2-specific metadata to the fitted object -----------------
  fit$call_mfp <- cl
  if (!identical(family_string, "negbin")) {
    fit$family <- mfp2_strip_prepared_family(family)
  }
  fit$family_string <- family_string
  if (identical(family_string, "finegray")) {
    # Preserve the expanded offset stored by the native coxph fit: its mean is
    # part of predict.coxph()'s linear-predictor origin. Summary refits instead
    # use this original-row vector and expand it through the cached row map.
    fit$mfp2_original_offset <- offset
  } else {
    fit$offset <- offset
  }
  fit$has_offset <- has_offset
  if (identical(family_string, "gee")) {
    fit$gee_input_row_order <- gee_input_row_order
  }

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
                         ties = c("breslow", "efron"),
                         strata = NULL,
                         nocenter = c(-1, 0, 1),
                         ftest = FALSE,
                         control = NULL,
                         min_saz_prop = 0.10,
                         verbose = TRUE,
                         fitter = c("base", "fastglm"),
                         id = NULL,
                         waves = NULL,
                         ...) {
  # mfp2.formula() translates a formula + data.frame specification into the
  # matrix/vector inputs expected by mfp2.default(): it expands categorical
  # predictors via model.matrix(), extracts per-variable fp() settings, and
  # then calls mfp2.default() to perform the actual fitting.


  # Step 1: Capture the call and resolve multiple-choice arguments -------------
  call <- match.call()

  # Match coxph-style formula usage while preserving the released mfp2 API for
  # a transition period. Formula strata remain authoritative when both forms
  # are supplied, and this is the only warning emitted for that combination.
  strata_supplied <- !missing(strata)
  if (strata_supplied) {
    .Deprecated(
      new = "strata() in the formula",
      package = "mfp2",
      old = "strata",
      msg = paste0(
        "`strata` as an argument to `mfp2.formula()` is deprecated; ",
        "include `strata(...)` in the formula instead. If both are supplied, ",
        "the formula term is used."
      )
    )
  }

  # Preserve the unevaluated observation-level arguments before any validation
  # can force their promises. Formula methods conventionally look these names
  # up in `data` first and then in the formula environment. Resolve each one
  # exactly once below, before subsetting, so expressions with side effects are
  # not repeated and complete-vector validation retains its existing semantics.
  weights_expr <- substitute(weights)
  offset_expr <- substitute(offset)
  subset_expr <- substitute(subset)
  strata_expr <- substitute(strata)
  id_expr <- substitute(id)
  waves_expr <- substitute(waves)

  criterion <- match.arg(criterion)
  xorder <- match.arg(xorder)
  ties <- resolve_mfp_ties(ties)
  fitter <- match.arg(fitter)

  # acd_vars is not a top-level argument in the formula interface: ACD handling is
  # requested per-variable via fp(..., acd = TRUE), so reject the matrix-interface
  # spelling (acd_vars) with a clear pointer to the correct syntax.
  dots <- list(...)

  if ("acd_vars" %in% names(dots)) {
    stop(
      "`acd_vars` is not supported as an argument to `mfp2.formula()`. ",
      "Use `fp(..., acd = TRUE)` or `fp2(..., acd = TRUE)` inside the formula ",
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
  # fail similarly stats:offset() as suggested by Terry Therneau
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

  # Step 4: Evaluate and validate observation-level inputs ---------------------
  n_data <- nrow(data)

  # Match the data-mask semantics used by standard formula methods. Columns of
  # `data` take precedence and caller-scope objects remain available through
  # the normalized formula environment. Reuse only these resolved values from
  # this point onward.
  weights <- eval(
    weights_expr,
    envir = data,
    enclos = environment(formula_internal)
  )
  offset <- eval(
    offset_expr,
    envir = data,
    enclos = environment(formula_internal)
  )
  subset <- eval(
    subset_expr,
    envir = data,
    enclos = environment(formula_internal)
  )
  strata <- eval(
    strata_expr,
    envir = data,
    enclos = environment(formula_internal)
  )
  id <- eval(
    id_expr,
    envir = data,
    enclos = environment(formula_internal)
  )
  waves <- eval(
    waves_expr,
    envir = data,
    enclos = environment(formula_internal)
  )

 # strata required only for survival models
  if (!is.null(strata) && !mfp2_family_is_survival(family_string)) {
    stop(
      "! `strata` is only allowed for survival models.\n",
      "i Use `family = \"cox\"`, `survreg_family()`, or `finegray_family()`, or remove `strata`.",
      call. = FALSE
    )
  }

  # id required only for fine-gray and GEE models
  if (!is.null(id)) {
    if (!family_string %in% c("finegray", "gee")) {
      stop(
        "! `id` is only used with `family = finegray_family()` or `family = gee_family()`.",
        call. = FALSE
      )
    }
    if (!is.atomic(id) || length(id) != n_data || anyNA(id) ||
        is.matrix(id) || !is.null(dim(id))) {
      stop(
        "! `id` must be a one-dimensional vector or factor with one value for every row of `data` and no missing values.\n",
        "i Use the same ID for rows from the same subject or cluster.",
        call. = FALSE
      )
    }
  } else if (identical(family_string, "gee")) {
    stop(
      "! `id` is required for `family = gee_family()`: supply a cluster ",
      "identifier evaluated in `data`.",
      call. = FALSE
    )
  }

  # waves required only for GEE models
  if (!is.null(waves)) {
    if (!identical(family_string, "gee")) {
      stop("! `waves` is only used with `family = gee_family()`.", call. = FALSE)
    }
    if (!is.atomic(waves) || length(waves) != n_data || anyNA(waves) ||
        is.matrix(waves) || !is.null(dim(waves))) {
      stop(
        "! `waves` must have one non-missing positive integer for every row of `data`.\n",
        "i Use `waves` to number observations within each `id` cluster; numbers must be unique within a cluster.",
        call. = FALSE
      )
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

  # Apply the same strictly-positive weight contract as the matrix interface.
  # Validation covers the complete user-supplied vector, including rows later
  # excluded by `subset`, so the public API has one unambiguous weight rule.
  validate_model_weights(
    weights = weights,
    nobs = n_data
  )

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
      strata_unit <- if (is.vector(strata) || is.factor(strata)) {
        "values"
      } else {
        "rows"
      }

      stop(
        sprintf(
          "! `strata` must provide one value or row for each row of `data`.\n",
          "i `strata` has %d %s; `data` has %d rows.",
          strata_n, strata_unit, n_data
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

  # Survival-model strata: baseline hazards for Cox, scale for survreg, and the
  # censoring distribution for Fine--Gray.
  if (!is.null(attr(terms_formula, "specials")$strata)) {
    if (mfp2_family_is_survival(family_string)) {

      # untangle the terms for strata as in survival::coxph()
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


      terms_drop <- c(terms_drop, stemp$terms)
    } else {
      stop("! strata are only allowed for survival models.\n",
           "i Please remove any strata terms from the model formula.",
           call. = FALSE)
    }
  }

  # Offset: an offset() term in the formula takes precedence over the `offset`
  # argument, with a warning.
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

    offset <- stats::model.offset(mf)
    if (!identical(family_string, "multinomial")) {
      offset <- as.vector(offset)
    }
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
    offset <- if (is.matrix(offset)) {
      offset[fit_rows, , drop = FALSE]
    } else {
      offset[fit_rows]
    }
  }

  if (!is.null(offset)) {
    if (!is.numeric(offset)) {
      stop("! `offset` must be numeric.", call. = FALSE)
    }

    offset_rows <- if (is.matrix(offset)) nrow(offset) else length(offset)
    if (offset_rows != nrow(mf)) {
      stop("! `offset` must have one value per observation.", call. = FALSE)
    }
    if (identical(family_string, "multinomial")) {
      c_classes <- mfp2_multinomial_n_classes(y)
      if (!is.matrix(offset) || !ncol(offset) %in% c(c_classes - 1L, c_classes)) {
        stop(
          "! Multinomial `offset` must be an n x C class matrix or an n x (C - 1) reference-logit matrix.",
          call. = FALSE
        )
      }
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

    # The model-matrix column was renamed above (for example, "fp(age)" to
    # "age"). Update the existing `age` entry in `term_to_columns` to use
    # the new column name.
    for (i in seq_along(fp_columns)) {
      fp_var <- fp_vars[[i]]
      if (fp_var %in% names(term_to_columns)) {
        columns <- term_to_columns[[fp_var]]
        columns[columns == fp_columns[[i]]] <- fp_var
        term_to_columns[[fp_var]] <- columns
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
  # Keep the normalized df vector aligned with candidate powers for downstream
  # fitting metadata and compatibility with internal validation hooks.
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
    if (!is.null(id)) id <- id[fit_rows]
    if (!is.null(waves)) waves <- waves[fit_rows]
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

  # A numeric matrix alone cannot distinguish a two-valued polynomial contrast
  # from an ordinary binary predictor. Transport the formula-derived distinction
  # privately so the final transformation can mean-center non-indicator factor
  # contrasts without changing direct matrix-interface behavior.
  attr(x, "mfp2_center_method") <- formula_factor_center_methods(
    factor_terms = factor_terms,
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
                      fitter = fitter,
                      criterion = criterion,
                      select = unlist(select_list),
                      alpha = unlist(alpha_list),
                      keep = keep,
                      xorder = xorder,
                      powers = power_list,
                      ties = ties,
                      strata = strata,
                      id = id,
                      waves = waves,
                      nocenter = nocenter,
                      acd_vars = acdx,
                      ftest = ftest,
                      control = control,
                      zero_vars = zero,
                      catzero_vars = catzero,
                      spike_vars = spike,
                      force_max_fp_vars = force_max_fp_vars,
                      min_saz_prop = min_saz_prop,
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

#' Identify fractional polynomial terms
#'
#' Checks whether each string represents an `fp()` or `fp2()` call,
#' optionally qualified with `mfp2::` or `mfp2:::`.
#'
#' @param z A character vector of formula terms.
#' @return A logical vector with one value per element of `z`.
#' @keywords internal
#' @noRd
is_fp_term <- function(z) {
  grepl("^(mfp2:::?)?fp2?\\(.*\\)$", z)
}

#' Extract Coefficients from an `mfp2` Model
#'
#' Returns coefficients from a fitted [mfp2()] model. This is a method for the
#' generic [stats::coef()] function.
#'
#' @param object A fitted [mfp2()] object.
#' @param ... Not used.
#'
#' @return
#' A named numeric vector for scalar-response models. For multinomial models,
#' a matrix with one row per non-reference outcome and one column per fitted
#' term; powers are common across rows but coefficients are outcome specific.
#'
#' @seealso [mfp2()], [summary.mfp2()], [predict.mfp2()]
#'
#' @export
coef.mfp2 <- function(object, ...) {
  object$coefficients
}


#' Confidence Intervals for `mfp2` Model Parameters
#'
#' Wald confidence intervals for the coefficients of a fitted [mfp2()] model,
#' a method for the generic [stats::confint()].
#'
#' The intervals are
#' \eqn{\hat\beta \pm z_{1-\alpha/2}\,\mathrm{SE}(\hat\beta)}, with standard
#' errors taken from [stats::vcov()] of the fitted object. These are the same
#' Wald intervals summarised by [summary.mfp2()] (which displays its 95%
#' interval with the conventional rounded multiplier 1.96, whereas `confint()`
#' uses the exact normal quantile and honours `level`), and they apply uniformly
#' to every scalar-coefficient family. They are used in preference to the default
#' profile-likelihood intervals of [stats::confint.glm()] because those refit
#' the model on `mfp2`'s internally transformed design and fail; for GEE
#' (`gee_family()`) the covariance is the robust sandwich estimator, which is
#' the correct basis for GEE inference (and `geepack` supplies no `confint`
#' method). Multinomial and ordinal models, whose coefficients are not a single
#' named vector, defer to the interval method of the underlying fit.
#'
#' @param object A fitted [mfp2()] object.
#' @param parm Optional character vector of parameter names, or numeric indices,
#'   selecting which coefficients to report. Defaults to all.
#' @param level The confidence level. Default `0.95`.
#' @param ... Passed to the underlying interval method for multinomial/ordinal
#'   models.
#'
#' @return A matrix with a row per selected parameter and two columns giving the
#'   lower and upper confidence limits.
#'
#' @seealso [mfp2()], [summary.mfp2()], [vcov()]
#'
#' @export
confint.mfp2 <- function(object, parm, level = 0.95, ...) {
  # Multinomial and ordinal models carry non-vector coefficient structures
  # (per-logit matrices, threshold intercepts); keep their own interval method.
  if (object$family_string %in% c("multinomial", "ordinal")) {
    return(NextMethod())
  }

  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      level <= 0 || level >= 1) {
    stop("`level` must be a single number in (0, 1).", call. = FALSE)
  }

  a <- (1 - level) / 2
  probs <- c(a, 1 - a)
  pct <- paste(
    format(100 * probs, trim = TRUE, scientific = FALSE, digits = 3),
    "%"
  )

  cf <- object$coefficients

  # A selection run may retain no covariates at all (every variable dropped);
  # such a model has no coefficients to interval, so return an empty matrix
  # rather than delegating to a base method that errors on the empty fit.
  if (is.null(cf) || length(cf) == 0L) {
    return(matrix(numeric(0), nrow = 0L, ncol = 2L,
                  dimnames = list(character(0), pct)))
  }

  nm <- names(cf)
  if (is.null(nm)) nm <- as.character(seq_along(cf))

  # Wald standard errors from the fitted covariance (robust sandwich for GEE).
  # If the covariance is unavailable, the interval limits are NA rather than a
  # hard error, so a usable object is always returned.
  V <- tryCatch(as.matrix(stats::vcov(object)), error = function(e) NULL)
  ses <- stats::setNames(rep(NA_real_, length(nm)), nm)
  if (!is.null(V) && !is.null(rownames(V))) {
    shared <- nm[nm %in% rownames(V)]
    if (length(shared) > 0L) {
      ses[shared] <- sqrt(diag(V)[shared])
    }
  }

  ci <- unname(cf) + outer(unname(ses), stats::qnorm(probs))
  dimnames(ci) <- list(nm, pct)

  if (!missing(parm)) {
    sel <- if (is.character(parm)) parm else nm[parm]
    ci <- ci[sel, , drop = FALSE]
  }
  ci
}

#' Extract Residuals from an `mfp2` Model
#'
#' Returns residuals from a fitted [mfp2()] model. For GEE and multinomial
#' models, `type` can be `"deviance"` (the default), `"pearson"`,
#' `"working"`, or `"response"`.
#'
#' @details
#' For GEE models, residuals are based on the marginal fitted means.
#' Deviance residuals are calculated from the response family's deviance
#' function and the model's prior weights. Pearson, working, and response
#' residuals are obtained from the underlying GEE fit.
#'
#' For multinomial models, residuals are calculated separately for each
#' non-reference category. Each category is treated as a binary outcome,
#' using its fitted class probability. The result has one column per
#' non-reference category.
#'
#' For other models, residual calculation is passed to the underlying fit,
#' including its choice of default residual type.
#'
#' @param object A fitted [mfp2()] model.
#' @param ... For GEE and multinomial models, use `type` to select
#' `"deviance"` (the default), `"pearson"`, `"working"`, or `"response"`
#' residuals. For other models, `type` and any other arguments follow
#' the residual method for that model's family; see, for example,
#' [stats::residuals.glm()] for GLM models.
#'
#' @return A numeric vector for GEE and other scalar-response models, or a
#' numeric matrix with one row per observation and one column per
#' non-reference category for multinomial models.
#'
#' @seealso [mfp2()], [plot()], [stats::residuals()]
#' @export
residuals.mfp2 <- function(object, ...) {
  if (identical(object$family_string, "gee")) {
    dots <- list(...)
    type <- if (!is.null(dots[["type"]])) {
      match.arg(
        dots[["type"]],
        c("deviance", "pearson", "working", "response")
      )
    } else {
      "deviance"
    }

    if (identical(type, "deviance")) {
      # Deviance residuals are not provided by geepack, but they only require
      # the marginal mean and the family deviance function (not the working
      # correlation), so compute them directly. This matches Stata `fracplot`.
      family <- object$family$response_family
      mu <- as.numeric(object$fitted.values)
      y <- as.numeric(object$y)
      wt <- object$prior.weights
      if (is.null(wt) || length(wt) != length(y)) wt <- rep(1, length(y))
      dev <- family$dev.resids(y, mu, wt)
      res <- sign(y - mu) * sqrt(pmax(dev, 0))
      names(res) <- names(object$fitted.values)
      return(res)
    }

    obj <- mfp2_restore_gee_family(object)
    class(obj) <- setdiff(class(obj), "mfp2")
    res <- stats::residuals(obj, type = type)
    # geeglm stores its linear predictor and fitted values as one-column
    # matrices, so residuals.geeglm() returns an n-by-1 matrix. Coerce to a
    # plain (optionally named) vector, which is what callers and the
    # component-plus-residual plots expect.
    if (is.matrix(res) && ncol(res) == 1L) {
      res <- drop(res)
    }
    return(res)
  }

  if (identical(object$family_string, "multinomial")) {
    dots <- list(...)
    type <- if (!is.null(dots[["type"]])) {
      match.arg(
        dots[["type"]],
        c("deviance", "pearson", "working", "response")
      )
    } else {
      "deviance"
    }
    return(mfp2_multinomial_residuals(object, type = type))
  }

  NextMethod()
}


#' Residuals for a Multinomial `mfp2` Model
#'
#' `nnet::multinom` provides no `residuals()` method, so
#' `residuals.default()` is used and returns the raw response residuals
#' \eqn{y_{ik}-\hat\pi_{ik}} for every class (including the reference), with the
#' `type` argument silently ignored. This helper computes proper per-logit
#' residuals instead.
#'
#' Each non-reference logit \eqn{q} is treated as a binary sub-problem with
#' indicator \eqn{y_{iq}=\mathbf{1}\{\text{class}_i=q\}} and fitted probability
#' \eqn{\hat\pi_{iq}} (the marginal class probability from
#' `predict(type = "response")`). Deviance residuals then follow the binomial
#' convention
#' \eqn{\mathrm{sign}(y_{iq}-\hat\pi_{iq})\sqrt{-2[y_{iq}\log\hat\pi_{iq} +
#' (1-y_{iq})\log(1-\hat\pi_{iq})]}}, matching the residual used for the
#' component-plus-residual plot. Pearson, working, and response residuals use
#' the corresponding binary definitions.
#'
#' @param object A fitted multinomial [mfp2()] object.
#' @param type Residual scale, one of `"deviance"` (default), `"pearson"`,
#' `"working"`, or `"response"`.
#'
#' @return A numeric `n` by `Q` matrix, where `Q` is the number of
#' non-reference logits. Columns are named by the non-reference class labels,
#' in the order carried by the fitted coefficient matrix.
#'
#' @keywords internal
#' @noRd
mfp2_multinomial_residuals <- function(object, type = "deviance") {
  probabilities <- stats::predict(object, type = "response")
  probabilities <- as.matrix(probabilities)

  # Non-reference logits, in the order carried by the coefficient matrix, so
  # residual columns align with coef(), vcov() blocks, and term predictions.
  logits <- rownames(object$mfp2_coefficient_matrix)
  if (is.null(logits)) {
    logits <- setdiff(object$class_levels, object$reference_class)
  }

  y_labels <- as.character(object$y)
  eps <- .Machine$double.eps

  out <- matrix(
    NA_real_, nrow(probabilities), length(logits),
    dimnames = list(rownames(probabilities), logits)
  )
  for (q in logits) {
    y_q <- as.numeric(y_labels == q)
    p_q <- pmin(pmax(probabilities[, q], eps), 1 - eps)
    resid_q <- y_q - p_q
    out[, q] <- switch(
      type,
      response = resid_q,
      pearson  = resid_q / sqrt(p_q * (1 - p_q)),
      working  = resid_q / (p_q * (1 - p_q)),
      deviance = sign(resid_q) * sqrt(pmax(
        -2 * (y_q * log(p_q) + (1 - y_q) * log(1 - p_q)), 0
      ))
    )
  }
  out
}


#' Specify a Fractional-Polynomial Term
#'
#' Use `fp()` or `fp2()` inside an [mfp2()] or [mfpi()] formula to mark a
#' continuous predictor for fractional-polynomial (FP) selection. They
#' specify the same model term and record options for fitting; neither
#' transforms the values when the formula is created.
#'
#' `fp2()` was introduced to avoid a name conflict with `mfp::fp()`.
#' That conflict is now resolved: formulas passed to `mfp2()` or `mfpi()`
#' use this package's `fp()` even when `mfp` is attached.
#'
#' Add categorical predictors to the formula without `fp()` or `fp2()`.
#' Convert numeric category codes with `factor()`.
#'
#' @details
#' The default `df = 4` permits selection of a linear effect, a first-degree
#' FP (FP1), or a second-degree FP (FP2). Selection may also exclude the
#' predictor. Use `df = 2` to permit at most FP1, or `df = 1` for a linear
#' effect only. The fitting procedure may reduce the permitted complexity
#' when a predictor has few distinct values.
#'
#' Candidate powers default to `c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)`. A power
#' of zero represents a logarithm. `powers` changes the candidates searched;
#' it does not specify the powers of the final fitted function.
#'
#' Options set inside `fp()` apply to this term. In particular, `alpha` and
#' `select` default to `0.05` here even if different values are supplied to
#' the surrounding `mfp2()` call. Under p-value selection, `select = 1`
#' forces inclusion and `alpha = 1` prevents simplification of the permitted
#' functional form.
#'
#' @param x Continuous numeric predictor to include in the formula. Do not
#'   wrap a factor in `fp()`.
#' @param df Maximum FP complexity: `1` for linear only, `2` for up to FP1,
#'   or `4` (default) for up to FP2. Higher positive even values permit
#'   higher FP degrees.
#' @param alpha Significance level for comparisons between functional forms
#'   under p-value selection. Default `0.05`. For eligible spike-at-zero
#'   terms, also controls tests that remove components.
#' @param select Significance level for variable inclusion under p-value
#'   selection. Default `0.05`; use `1` to force inclusion. It does not
#'   control inclusion under AIC or BIC selection.
#' @param shift Numeric value added before FP transformation. `NULL` (default)
#'   inherits the setting from `mfp2()`; `0` disables shifting for this term.
#' @param scale Positive divisor applied after shifting. `NULL` (default)
#'   inherits the setting from `mfp2()`; `1` disables scaling for this term.
#' @param center Whether to center the fitted terms. Default `TRUE`.
#' @param acd Whether to request approximate cumulative distribution (ACD)
#'   modelling. Default `FALSE`.
#' @param powers Numeric vector of candidate FP powers for this predictor.
#'   `NULL` (default) uses the global candidate set. Duplicate candidates
#'   are removed.
#' @param zero If `TRUE`, transform only values with `x > 0`; exact zeros
#'   contribute zero to the continuous component. Requires nonnegative `x`.
#'   Default `FALSE`.
#' @param catzero If `TRUE`, add a binary indicator for `x = 0` alongside
#'   the positive component. Implies `zero = TRUE`. Default `FALSE`.
#' @param spike If `TRUE`, request spike-at-zero (SAZ) model selection,
#'   subject to eligibility checks. While eligible, it considers a positive
#'   component and a binary indicator for `x = 0`. Requires nonnegative `x`.
#'   Default `FALSE`.
#' @param force_max_fp If `TRUE`, retain the most complex form permitted by
#'   `df` and bypass variable removal and functional-form simplification.
#'   For an eligible SAZ term, retain both components. Default `FALSE`.
#' @param acdx Deprecated alias for `acd`. Use `acd`; supplying both
#'   arguments is an error.
#'
#' @return The input vector with modelling options attached as attributes.
#'   The options take effect when the term is used in an [mfp2()] or
#'   [mfpi()] formula.
#'
#' @examples
#' data("prostate")
#'
#' # Search up to FP1 for age.
#' fit_fp <- mfp2(
#'   lpsa ~ fp(age, df = 2) + svi,
#'   data = prostate,
#'   verbose = FALSE
#' )
#'
#' # fp2() accepts the same options and produces the same model specification.
#' fit_fp2 <- mfp2(
#'   lpsa ~ fp2(age, df = 2) + svi,
#'   data = prostate,
#'   verbose = FALSE
#' )
#'
#' @seealso [fp2()], [mfp2()], [mfpi()]
#' @export
fp <- function(x,
               df = 4,
               alpha = 0.05,
               select = 0.05,
               shift = NULL,
               scale = NULL,
               center = TRUE,
               acd = FALSE,
               powers = NULL,
               zero = FALSE,
               catzero = FALSE,
               spike = FALSE,
               force_max_fp = FALSE,
               acdx = NULL
) {

  acd_supplied <- !missing(acd)
  acdx_supplied <- !missing(acdx)

  if (acd_supplied && acdx_supplied) {
    stop(
      "`acd` and its deprecated alias `acdx` cannot both be supplied.",
      call. = FALSE
    )
  }

  if (acdx_supplied) {
    .Deprecated(
      new = "acd",
      package = "mfp2",
      old = "acdx",
      msg = "`acdx` is deprecated; use `acd` instead."
    )
    acd <- acdx
  }

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
  validate_scalar_fp(acd, "acd", name, "logical", allow_null = FALSE, allow_na = FALSE)
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
  attr(x, "acd") <- acd
  attr(x, "powers") <- powers
  attr(x, "zero") <- zero
  attr(x, "catzero") <- catzero
  attr(x, "spike") <- spike
  attr(x, "force_max_fp") <- force_max_fp
  attr(x, "name") <- name

  x
}
#' @rdname fp
#' @param ... Arguments passed by `fp2()` to `fp()`.
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
#' Assigns a degree-of-freedom (df) setting to each column of a predictor
#' matrix. Columns with few distinct values receive a lower df than requested.
#' The fitting pipeline supplies columns without missing values or constant
#' predictors; this helper does not check those conditions itself.
#'
#' @section Assignment rules:
#' Let \eqn{u} be the number of distinct values in a column and \eqn{d} its
#' requested df. For columns supplied by the fitting pipeline:
#'
#' \tabular{lll}{
#'   \strong{Distinct values} \tab \strong{Assigned df} \tab \strong{Rule} \cr
#'   \eqn{2 \le u \le 3} \tab \code{1}         \tab Use a linear term. \cr
#'   \eqn{4 \le u \le 5} \tab \code{min(2, d)} \tab Cap df at 2. \cr
#'   \eqn{u \ge 6}       \tab \code{d}         \tab Keep the requested df. \cr
#' }
#'
#' The implementation also assigns df 1 when a directly supplied column has
#' fewer than two distinct values. A directly supplied column containing
#' missing values has those values counted by \code{unique()}.
#'
#' @param x A predictor matrix with column names. Each column represents one
#'   variable.
#' @param df_default A single value for all columns or a vector with one
#'   value per column. Values are converted to integers and must be positive
#'   after conversion. Default \code{4}.
#'
#' @return A named integer vector with one assigned df per column, in the
#'   same order as \code{x}. Names come from \code{colnames(x)}.
#'
#' @examples
#' x <- cbind(
#'   binary     = c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1),
#'   few_levels = c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5),
#'   continuous = 1:10
#' )
#'
#' assign_df(x, df_default = 4)
#' # binary = 1, few_levels = 2, continuous = 4
#'
#' assign_df(x, df_default = c(1, 4, 2))
#' # binary = 1, few_levels = 2, continuous = 2
#'
#' @keywords internal
#' @noRd
assign_df <- function(x, df_default = 4L) {
  if (!is.matrix(x)) {
    stop("`x` must be a matrix.", call. = FALSE)
  }

  v_names <- colnames(x)
  if (is.null(v_names)) {
    stop("`x` must have column names.", call. = FALSE)
  }

  # A scalar applies to every column; a vector supplies one value per column.
  n_vars <- length(v_names)
  if (!(length(df_default) %in% c(1L, n_vars))) {
    stop(
      "`df_default` must have length 1 or the number of columns (",
      n_vars, ").",
      call. = FALSE
    )
  }

  df <- rep(as.integer(df_default), length.out = n_vars)
  if (anyNA(df) || any(df < 1L)) {
    stop("`df_default` must contain only positive integers.", call. = FALSE)
  }

  # Count distinct values separately for each predictor.
  nu <- vapply(
    seq_len(n_vars),
    function(j) length(unique(x[, j])),
    integer(1L)
  )

  # Two or three distinct values permit only a linear term.
  df[nu <= 3L] <- 1L

  # Four or five distinct values permit at most FP1 (df = 2).
  mid <- nu >= 4L & nu <= 5L
  df[mid] <- pmin(df[mid], 2L)

  # Align the result with the matrix columns.
  stats::setNames(df, v_names)
}
