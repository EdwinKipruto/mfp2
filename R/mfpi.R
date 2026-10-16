#' Model Interactions with a Categorical Grouping Variable
#'
#' Investigate whether the association between an outcome and a continuous,
#' binary, or categorical variable differs between groups. Continuous
#' associations can be modelled using fractional polynomials (FPs).
#'
#' @details
#' `mfpi()` first selects an adjustment model using [mfp2()]. It then evaluates
#' each variable in `interaction_vars` separately, comparing a model with a
#' common association across groups against a model with group-specific
#' associations.
#'
#' Each comparison includes group main effects and the selected adjustment
#' terms, excluding the variable currently being tested. The function returns
#' a separate interaction model for each tested variable; it does not combine
#' the selected interactions into one model.
#'
#' @param x A predictor matrix or data frame for the default method. Do not
#'   include an intercept column. See *Data requirements*.
#' @param y The response for the default method, with one observation per row
#'   of `x`. Its format must match `family`; see *Supported models*.
#' @param formula A model formula containing the response, grouping variable,
#'   interaction variables, and adjustment variables. See *Formula interface*.
#' @param data A data frame containing the variables used in `formula`.
#' @param group_var A character string naming the categorical grouping
#'   variable, which must have at least two observed groups.
#' @param interaction_vars A non-empty character vector naming the variables
#'   to test for interaction with `group_var`.
#' @param interaction_forms An optional named character vector specifying
#'   `"linear"`, `"fp1"`, or `"fp2"` for each tested variable. Unspecified forms
#'   default to `"fp1"` for continuous variables and `"linear"` for binary and
#'   categorical variables. See *Interaction forms*.
#' @param flex How FP powers are selected: `"flex1"`, `"flex2"`, `"flex3"`
#'   (default), or `"flex4"`. See *Flexibility levels* and
#'   *Statistical interpretation*.
#' @param p_interact The interaction p-value threshold, in `(0, 1]`.
#'   Default: `0.05`. Used only with `criterion = "pvalue"`.
#' @param min_improvement A positive minimum AIC or BIC improvement required
#'   to select an interaction. `NULL` uses 2. See *Interaction selection*.
#' @param include_group_var Whether to force the grouping variable into the
#'   initial adjustment model. Default: `FALSE`. See *Adjustment model*.
#' @param show_models Whether to print summaries of selected interaction
#'   models when `verbose = TRUE`. Default: `FALSE`.
#' @param weights Optional strictly positive observation weights.
#'   See *Weights, offsets, and strata*.
#' @param offset An optional offset for the default method, added to the linear
#'   predictor in both fitting stages. See *Weights, offsets, and strata*.
#' @param cycles The maximum number of MFP fitting cycles. A positive integer;
#'   default: 10.
#' @param scale Positive scaling factors, or `NULL` for automatic scaling.
#'   See *Predictor settings* and *Shifting, scaling, and centring*.
#' @param shift Finite shifting constants, or `NULL` for automatic shifting.
#'   See *Predictor settings* and *Shifting, scaling, and centring*.
#' @param df The maximum FP complexity in the adjustment model: 1 for linear,
#'   2 for FP1, 4 for FP2, or a larger even positive integer. Default: 4.
#'   See *Predictor settings*.
#' @param center Whether to centre transformed predictors. Default: `TRUE`.
#'   See *Predictor settings* and *Shifting, scaling, and centring*.
#' @param subset An optional logical vector or vector of unique positive row
#'   indices selecting observations. The formula method also accepts an
#'   expression using `data`. See *Data requirements*.
#' @param family The regression family. Default: `"gaussian"`.
#'   See *Supported models*.
#' @param fitter The GLM fitting method: `"base"` (default) or `"fastglm"`.
#'   See *Fitting options*.
#' @param criterion The selection criterion for adjustment and interaction
#'   models: `"pvalue"` (default), `"aic"`, or `"bic"`.
#' @param select The variable-selection p-value threshold for the adjustment
#'   model, in `[0, 1]`. Default: `0.05`. See *Adjustment model* and
#'   *Predictor settings*.
#' @param alpha The FP-degree selection threshold for the adjustment model,
#'   in `[0, 1]`. Default: `0.05`. See *Adjustment model* and
#'   *Predictor settings*.
#' @param keep An optional character vector naming adjustment variables or
#'   complete categorical terms that must be retained.
#' @param force_max_fp_vars An optional character vector naming adjustment
#'   variables that must be retained at the maximum FP degree allowed by `df`.
#'   See *Adjustment model*.
#' @param xorder The order of variable assessment in the adjustment model:
#'   `"ascending"` (default), `"descending"`, or `"original"`.
#'   See *Adjustment model*.
#' @param powers An optional named list of candidate FP powers for adjustment
#'   and interaction functions. See *Interaction forms*.
#' @param ties The Cox method for tied event times: `"efron"` (default) or
#'   `"breslow"`. Ignored for other families.
#' @param strata Optional survival-model strata for the default method.
#'   See *Weights, offsets, and strata*.
#' @param nocenter Values identifying Cox predictor columns that should not
#'   be recentred by the Cox fitter. Default: `c(-1, 0, 1)`, as in
#'   [survival::coxph()]. Use `NULL` to allow recentring of all columns.
#' @param acd_vars Adjustment variables to model using an approximate
#'   cumulative distribution (ACD) transformation. Available directly only
#'   in the default method. See *Special transformations*.
#' @param zero_vars Variables for which exact zeros are treated separately
#'   from the positive continuous component. See *Special transformations*.
#' @param catzero_vars Adjustment variables requiring a positive-part function
#'   and a separate indicator for exact zeros. See *Special transformations*.
#' @param spike_vars Adjustment variables to assess using the spike-at-zero
#'   procedure. See *Special transformations*.
#' @param min_saz_prop The minimum proportion required in each component of
#'   a spike-at-zero adjustment variable. A number in `(0, 0.5)`;
#'   default: `0.10`.
#' @param ftest Whether to use F-tests instead of likelihood-ratio chi-square
#'   tests for Gaussian p-value comparisons in both fitting stages.
#'   Default: `FALSE`. Ignored for other families.
#' @param control Control settings for the regression fitter, or `NULL` for
#'   defaults. See *Fitting options*.
#' @param winsorize Whether to Winsorise continuous interaction variables.
#'   Default: `FALSE`. See *Winsorisation*.
#' @param winsorize_probs The lower and upper Winsorisation probabilities,
#'   in increasing order within `[0, 1]`. Default: `c(0.01, 0.99)`.
#' @param center_type The centring reference for interaction functions:
#'   `"grand"` (default) or `"group"`. See
#'   *Shifting, scaling, and centring*.
#' @param p_adjust_method A method accepted by [stats::p.adjust()].
#'   Default: `"none"`. See *Multiple testing*.
#' @param verbose Whether to print progress and selection results.
#'   Default: `TRUE`.
#' @param digits The number of decimal places in printed interaction results,
#'   excluding FP powers and degrees of freedom. A positive integer;
#'   default: 3.
#' @param term_groups An optional named list grouping numeric columns into
#'   complete categorical terms for the default method.
#'   See *Categorical variables*.
#' @param ... Additional arguments are not supported.
#'
#' @section Interaction forms:
#' Specify the interaction form before fitting:
#'
#' - `"linear"` fits a linear continuous effect, a binary effect, or a complete
#'   categorical effect in each group, without searching FP powers.
#' - `"fp1"` selects one FP power.
#' - `"fp2"` selects a two-term FP function, including repeated powers.
#'
#' The function does not select between these three forms. Unspecified forms,
#' including when `interaction_forms = NULL`, default to `"fp1"` for continuous
#' variables and `"linear"` for binary and categorical variables. Explicit
#' requests for `"fp1"` or `"fp2"` for binary or categorical variables are
#' changed to `"linear"` with a warning.
#'
#' The default candidate powers are
#' `c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)`, where zero denotes the natural
#' logarithm. This gives eight FP1 functions and 36 FP2 combinations.
#' Power 1 is included, so an FP1 search can select a linear function.
#' Repeated FP2 powers use the usual logarithmic second term.
#'
#' Supply `powers` as a named list to change the candidate set for individual
#' variables, for example `powers = list(age = c(0, 0.5, 1))`.
#' These candidates apply to both adjustment and interaction functions.
#'
#' @section Flexibility levels:
#' For FP interactions, `flex` determines how powers are selected:
#'
#' - `"flex1"` selects powers in the model without interaction and uses them
#'   in every group.
#' - `"flex2"` selects shared powers in the interaction model and uses those
#'   powers in the model without interaction.
#' - `"flex3"` selects powers separately for the models with and without
#'   interaction. The interaction powers remain common to all groups.
#' - `"flex4"` allows each group to have its own powers and selects the
#'   common function separately.
#'
#' Group-specific coefficients are estimated under every option.
#' `flex` has no effect on interactions specified as `"linear"`.
#'
#' @section Adjustment model:
#' The adjustment model is selected once using the MFP procedure in [mfp2()].
#' Its selected terms and transformations are then held fixed during the
#' separate interaction analyses. Each analysis excludes the tested variable
#' from its adjustment terms because its main effect and interaction are
#' fitted explicitly.
#'
#' By default, the grouping variable is excluded from the initial adjustment
#' model. Set `include_group_var = TRUE` to force it into that model.
#' Group main effects are included in every interaction comparison regardless
#' of this setting.
#'
#' With `criterion = "pvalue"`, `select` controls variable retention and
#' `alpha` controls FP-degree selection. Setting `select = 1` retains a
#' variable; setting `alpha = 1` requires the maximum permitted FP degree.
#' With `"aic"` or `"bic"`, adjustment-model selection uses the corresponding
#' information criterion.
#'
#' Use `keep` to retain specified adjustment terms. Use `force_max_fp_vars`
#' to retain specified variables at their maximum permitted FP degree under
#' any criterion; powers within that degree are still selected from the data.
#' For a spike-at-zero adjustment variable, this retains its full permitted
#' representation.
#'
#' `xorder = "ascending"` processes terms from the smallest to the largest
#' p-value for their contribution to the initial full linear model.
#' Each p-value compares that model with a model omitting the complete term,
#' while retaining the other terms. `"descending"` reverses this order;
#' `"original"` preserves the predictor or formula term order.
#'
#' These settings control adjustment-model selection. They do not determine
#' the interaction form, which is set by `interaction_forms`.
#' If no adjustment terms are selected, the interaction analyses proceed
#' without them.
#'
#' @section Interaction selection:
#' Each requested interaction is selected when:
#'
#' - its raw or adjusted p-value is strictly below `p_interact` for
#'   `criterion = "pvalue"`;
#' - `AIC_main_minus_int` is strictly greater than `min_improvement` for
#'   `criterion = "aic"`;
#' - `BIC_main_minus_int` is strictly greater than `min_improvement` for
#'   `criterion = "bic"`.
#'
#' Positive AIC or BIC differences favour the interaction. The default
#' improvement threshold is 2. `min_improvement` does not affect the
#' adjustment-model selection thresholds.
#'
#' The reported AIC and BIC values use MFPI complexity allowances for searched
#' FP powers and should not be compared directly with [stats::AIC()] or
#' [stats::BIC()] results.
#'
#' @section Statistical interpretation:
#' Interaction p-values are based on likelihood-ratio comparisons, or
#' Gaussian F-tests when requested. They do not account fully for uncertainty
#' from selecting the adjustment model and FP functions.
#'
#' P-values for `"flex3"`, and `"flex4"` are approximate. With
#' `"flex3"` and `"flex4"`, different selected powers can make the compared
#' models non-nested.
#' If the interaction model has a larger deviance than the model without
#' interaction beyond numerical tolerance, its p-value is reported as `NA`;
#' AIC and BIC comparisons remain available.
#'
#' For two groups and a single outcome linear predictor, interaction tests
#' use 1 degree of freedom for FP1 and 2 for FP2 with `"flex1"`, `"flex2"`,
#' and `"flex3"`. With `"flex4"`, the corresponding values are 2 and 4.
#' These follow Royston and Sauerbrei (2014, Web Appendix, Table 2).
#'
#' More flexible functions can be unstable with sparse data. Inspect the fitted
#' functions and the distributions of the tested variables when interpreting an
#' interaction.
#'
#' @section Multiple testing:
#' P-values are unadjusted by default. Set `p_adjust_method` to a method such
#' as `"holm"`, `"bonferroni"`, or `"BH"` to adjust across the requested
#' interaction tests.
#'
#' There is one planned test per variable in `interaction_vars`. Unavailable
#' p-values remain `NA` and are included in the planned number of tests.
#' An interaction with an unavailable p-value cannot be selected by the
#' p-value criterion.
#'
#' This adjustment applies only to `criterion = "pvalue"` and does not
#' correct for model-selection uncertainty.
#'
#' @section Data requirements:
#' The grouping variable and all interaction variables must be included in
#' `x` or `formula`. The grouping variable cannot also be an interaction
#' variable or a member of `term_groups`.
#'
#' For the default method, supply numeric predictors in a matrix or data frame.
#' Data frames also accept logical predictors and a factor or character
#' grouping variable. Predictor names must be unique, non-missing, non-empty,
#' and contain no backticks.
#'
#' Missing and non-finite fitting values are not supported. The response,
#' weights, offset, and strata must correspond to the predictor rows.
#' At least five observations and two groups must remain after subsetting.
#' Each group needs at least two observations, or three if any FP2
#' interaction is requested.
#'
#' These minimum counts do not guarantee that a model can be estimated.
#' The requested effects must have sufficient variation within groups and
#' must be distinguishable from the other model terms. FP searches for
#' continuous variables with five or fewer distinct values produce a warning.
#'
#' The reference group is the first observed factor level, the first retained
#' distinct character value, or the smallest observed logical or numeric
#' value. To choose a reference explicitly, supply a factor with the desired
#' level order.
#'
#' @section Formula interface:
#' Supply original variables in `data`. Use `fp()` or `fp2()` to specify
#' variable-specific FP and preprocessing settings. These settings override
#' the corresponding global defaults.
#'
#' Global `zero_vars`, `catzero_vars`, `spike_vars`, and
#' `force_max_fp_vars` settings are combined with the corresponding options
#' in the formula wrappers.
#'
#' The wrappers do not set the interaction form. For example,
#' `fp(age, df = 4)` permits FP2 in the adjustment model, while
#' `interaction_forms = c(age = "linear")` requests a linear age interaction.
#'
#' Use `fp(variable, acd = TRUE)` for ACD transformations, `offset(...)`
#' for offsets, and `strata(...)` for survival strata. The direct `acd_vars`,
#' `offset`, and `strata` arguments are unavailable in the formula method.
#'
#' @section Categorical variables:
#' In the formula interface, categorical predictors are represented and
#' assessed as complete terms. They can be adjustment variables or interaction
#' variables with form `"linear"`. Do not wrap categorical predictors in
#' `fp()` or `fp2()`.
#'
#' In the default interface, represent a categorical predictor by numeric
#' contrast columns and identify them with a named list, for example
#' `term_groups = list(stage = c("stageII", "stageIII"))`.
#' Each column can belong to only one term.
#'
#' Name the complete term, such as `"stage"`, in `interaction_vars` or `keep`.
#' Individual columns of a grouped term cannot be tested separately.
#' Grouped terms are linear and cannot use special transformations or
#' `force_max_fp_vars`.
#'
#' Binary and categorical interaction variables retain their coding and are
#' not shifted, scaled, centred, Winsorised, or subjected to an FP search.
#' For formula-based prediction, supply categorical variables in their
#' original form.
#'
#' @section Predictor settings:
#' In the default method, `df`, `select`, `alpha`, `center`, `shift`, and
#' `scale` accept either an unnamed scalar or a named vector.
#' A scalar supplies a common setting. A named vector sets only the specified
#' columns of `x`; omitted columns use the argument's default.
#' Names must be unique and match predictor columns. Unnamed vectors of
#' length greater than one are not accepted.
#'
#' In the formula method, these global arguments must be scalars.
#' Use `fp()` or `fp2()` for variable-specific settings.
#'
#' @section Shifting, scaling, and centring:
#' Predictors are shifted and scaled before FP modelling. `NULL` requests
#' automatic values; `shift = 0` disables shifting and `scale = 1` disables
#' scaling. Explicit values must be finite, and scales must be positive.
#' A supplied shift is not increased automatically if it leaves values
#' unsuitable for a nonlinear FP transformation.
#'
#' Linear-only and structural-zero predictors use a zero shift.
#' Binary predictors and the grouping variable use zero shift and unit scale.
#'
#' For interaction functions, `center_type = "grand"` centres transformed
#' values using a common mean across groups; `"group"` uses within-group
#' means. With `"flex4"`, a separate centring value is calculated for each
#' group's selected function. Fitted centring values are reused for prediction.
#'
#' For a structural-zero interaction variable, centring means use positive
#' observations only, and the continuous contribution remains zero at exact
#' zeros. Group centring therefore also changes the reference for the positive
#' component and should be chosen deliberately. Each group must retain at
#' least two positive observations for a linear or FP1 form, or three for FP2.
#'
#' Structural-zero adjustment variables follow the [mfp2()] convention:
#' their uncentred continuous functions are zero at exact zeros, and centring
#' subtracts the fitted mean of the complete column from every observation.
#'
#' @section Special transformations:
#' The following options require nonnegative variables. Negative values must
#' be recoded explicitly if they are intended to represent structural zeros.
#'
#' - `zero_vars` applies the continuous function only to positive values,
#'   with zero continuous contribution at exact zeros before centring.
#' - `catzero_vars` additionally includes an indicator for exact zeros and
#'   implies positive-part modelling.
#' - `spike_vars` assesses alternative representations of the zero and positive
#'   components using the spike-at-zero procedure.
#'
#' For adjustment variables, these can also be specified using
#' `fp(variable, zero = TRUE)`, `catzero = TRUE`, or `spike = TRUE`.
#' ACD transformations are available for adjustment variables through
#' `acd_vars` or `fp(variable, acd = TRUE)`. With structural-zero handling,
#' ACD is estimated from positive values only.
#'
#' Interaction variables support explicit `zero_vars` handling, but not ACD,
#' categorical-zero indicators, or spike-at-zero selection. Unsupported
#' options are disabled with a warning; explicit zero handling is preserved.
#' A zero-handled interaction variable must retain positive observations
#' after subsetting. Binary variables use ordinary binary modelling.
#'
#' Spike-at-zero assessment requires at least `min_saz_prop` of observations
#' in each component. Variables that fail this requirement retain the
#' applicable zero or categorical-zero handling without spike-at-zero selection.
#' If a variable contains no zeros, special zero handling is disabled.
#'
#' @section Winsorisation:
#' With `winsorize = TRUE`, continuous interaction values below or above the
#' specified empirical quantiles are replaced by those limits.
#' No observations are removed.
#'
#' Limits are calculated from the fitting observations after shifting and
#' scaling and are reused during prediction. Binary, categorical, and
#' structural-zero interaction variables are excluded.
#' Inspect influential observations and report any Winsorisation used.
#'
#' @section Supported models:
#' Supported GLM families are `"gaussian"`, `"binomial"`, `"poisson"`,
#' `"Gamma"`, and `"inverse.gaussian"`, supplied as a family name, function,
#' or object. Responses must satisfy the corresponding family requirements.
#' Grouped binomial responses require integer success and failure counts.
#'
#' Additional models are:
#'
#' - `"negbin"` for negative-binomial regression;
#' - `"multinomial"` or [multinomial_family()] for outcomes with at least
#'   three classes, supplied as class labels or a matrix of class counts;
#' - `"ordinal"` or [ordinal_family()] for ordered outcomes with at least
#'   three categories;
#' - `"cox"` for a right-censored [survival::Surv()] response;
#' - [survreg_family()] for parametric survival regression with a
#'   non-counting-process `Surv` response.
#'
#' Quasi-likelihood, GEE, and Fine--Gray models are not supported.
#'
#' Multinomial models share FP powers across outcome logits but estimate
#' separate regression coefficients for each logit. Their interaction degrees
#' of freedom account for these separate coefficients.
#'
#' Ordinal models use common slopes across response cut-points.
#' Supply an ordered factor to define the category order explicitly;
#' otherwise, numeric categories are ordered increasingly and character or
#' unordered-factor categories alphabetically. Ordinal fitting requires
#' the optional `rms` package (version 7.0-0 or later).
#'
#' @section Weights, offsets, and strata:
#' Weights must be finite and strictly positive and are used in both fitting
#' stages. Ordinal models use the supplied case weights without rescaling them.
#'
#' Offsets must be finite and are used in both stages. In the default method,
#' supply the `offset` argument. In the formula interface, include `offset(...)`
#' in the formula; a direct `offset` argument is not available. Weight and
#' formula offset expressions are evaluated in `data` and then in the formula
#' environment.
#'
#' For Cox models, strata define separate baseline hazards. For `survreg`
#' models, they define separate scale parameters. The default method accepts
#' stratum labels or a matrix or data frame of stratification variables.
#' The formula method uses `strata(...)`.
#'
#' @section Fitting options:
#' `fitter = "fastglm"` requires the optional `fastglm` package.
#' Ordinary GLMs warn and retry with the base fitter if it is unavailable
#' or fitting fails. Negative-binomial models require `fastglm` and select
#' it automatically. Survival models are unaffected by `fitter`.
#'
#' Supply `control` using the appropriate settings for the model, usually
#' [stats::glm.control()], [survival::coxph.control()], or
#' [survival::survreg.control()]. For `fastglm`, `epsilon` and `maxit`
#' control regression convergence; `trace = TRUE` is not supported.
#' For negative-binomial models, these controls do not change the separate
#' dispersion-estimation settings.
#'
#' @return
#' An object of class `"mfpi"`. Use `summary()`, `coef()`, `vcov()`,
#' `predict()`, and `plot()` to examine the results.
#'
#' Documented components include:
#'
#' - `all_model_metrics`: interaction results for every evaluated variable,
#'   including p-values, degrees of freedom, information-criterion comparisons,
#'   and selected powers.
#' - `best_model_metrics`: results for interactions selected by `criterion`.
#' - `all_interaction_models`: fitted interaction models for all evaluated
#'   variables.
#' - `best_interaction_model`: fitted models for selected interactions.
#' - `adjustment_model`: the selected [mfp2()] adjustment model, accepting
#'   prediction data on the original covariate scale.
#' - `interaction_forms`: the form used for each tested variable.
#' - `group_level_map`: the correspondence between fitted group values
#'   and original group labels.
#' - `winsorize_limits`: fitted Winsorisation limits after shifting and scaling,
#'   when requested.
#' - `criterion`: the selection criterion.
#' - `p_adjust_method`: the multiple-testing adjustment method.
#' - `call`: the matched function call.
#'
#' @references
#' Royston, P. and Sauerbrei, W. (2004). A new approach to modelling
#' interactions between treatment and continuous covariates in clinical trials
#' by using fractional polynomials. *Statistics in Medicine*, 23, 2509–2525.
#' \doi{10.1002/sim.1815}.
#'
#' Royston, P. and Sauerbrei, W. (2009). Two techniques for investigating
#' interactions between treatment and continuous covariates in clinical trials.
#' *The Stata Journal*, 9, 230–251.
#' \doi{10.1177/1536867X0900900204}.
#'
#' Royston, P. and Sauerbrei, W. (2013). Interaction of treatment with a
#' continuous variable: simulation study of significance level for several
#' methods of analysis. *Statistics in Medicine*, 32, 3788–3803.
#' \doi{10.1002/sim.5813}.
#'
#' Royston, P. and Sauerbrei, W. (2014). Interaction of treatment with a
#' continuous variable: simulation study of power for several methods of
#' analysis. *Statistics in Medicine*, 33, 4695–4708.
#' \doi{10.1002/sim.6308}.
#'
#' @examples
#' data("prostate")
#'
#' # Specify each interaction form separately from the adjustment settings.
#' fit <- mfpi(
#'   lpsa ~ svi + fp(age) + fp(cavol) + fp(weight) + fp(bph) + fp(cp),
#'   data = prostate,
#'   group_var = "svi",
#'   interaction_vars = c("cavol", "age"),
#'   interaction_forms = c(cavol = "fp1", age = "linear"),
#'   flex = "flex1",
#'   verbose = FALSE
#' )
#'
#' print(fit)
#' summary(fit)
#'
#' # Examine an interaction even if it was not selected.
#' pred <- predict(
#'   fit,
#'   terms = "cavol",
#'   model = "all",
#'   type = "both",
#'   grid = TRUE
#' )
#' pred$functions
#' pred$differences
#'
#' \donttest{
#' if (requireNamespace("patchwork", quietly = TRUE)) {
#'   plot(fit, terms = "cavol", model = "all", plot_type = "both")
#' }
#' }
#'
#' @seealso [mfp2()], [fp()], [coef.mfpi()], [vcov.mfpi()],
#'   [predict.mfpi()], [plot.mfpi()], [summary.mfpi()],
#'   [multinomial_family()], [ordinal_family()], [survreg_family()]
#'
#' @export
mfpi <- function(x, ...) {
  UseMethod("mfpi", x)
}

#' Check Whether a Model Family Is Supported by MFPI
#'
#' Rejects Fine--Gray and GEE models. Fine--Gray models use a weighted partial
#' pseudo-likelihood for which the standard likelihood-ratio comparisons used
#' by MFPI are not appropriate. GEE models use estimating equations rather
#' than a full likelihood. The MFPI interface also lacks the cluster arguments
#' required for GEE fitting.
#'
#' Both families are supported by [mfp2()].
#'
#' @param family_string A character string giving the canonical family name.
#'
#' @return Invisibly returns `NULL` unless the family is Fine--Gray or GEE,
#'   in which case it raises an error.
#'
#' @keywords internal
#' @noRd
validate_mfpi_family <- function(family_string) {
  if (identical(family_string, "finegray")) {
    stop(
      "! `mfpi()` does not support Fine--Gray models.",
      call. = FALSE
    )
  }
  if (identical(family_string, "gee")) {
    stop(
      "! `mfpi()` does not support GEE models (`gee_family()`). ",
      "GEE is available in `mfp2()`, which takes the required `id` (and ",
      "optional `waves`) cluster arguments.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Describe MFPI Interaction Terms
#'
#' Identifies continuous, binary, and categorical interaction terms and maps
#' each term to its predictor columns. Continuous and binary terms use one
#' column; categorical terms use the complete contrast block recorded in
#' `term_to_columns`.
#'
#' @param x A numeric predictor matrix.
#' @param interaction_vars Character vector naming the interaction terms.
#' @param term_to_columns A named list mapping terms to predictor column names.
#'
#' @return A named list describing each term's type, predictor columns, and
#'   block width, with observed-level information for discrete terms.
#'
#' @keywords internal
#' @noRd
mfpi_interaction_specs <- function(x, interaction_vars, term_to_columns) {
  specs <- stats::setNames(vector("list", length(interaction_vars)),
                           interaction_vars)

  for (term in interaction_vars) {
    columns <- term_to_columns[[term]]
    if (is.null(columns) && term %in% colnames(x)) columns <- term
    columns <- unique(as.character(columns))

    if (length(columns) == 0L || any(!columns %in% colnames(x))) {
      stop(
        paste0(
          "! Interaction variable `", term,
          "` could not be resolved to fitted design columns."
        ),
        call. = FALSE
      )
    }

    identity_term <- length(columns) == 1L && identical(columns, term)
    observed <- x[, columns, drop = FALSE]
    n_unique <- if (identity_term) length(unique(observed[, 1L])) else NA_integer_
    kind <- if (!identity_term || length(columns) > 1L) {
      "categorical"
    } else if (n_unique == 2L) {
      "binary"
    } else {
      "continuous"
    }

    level_design <- NULL
    level_labels <- NULL
    if (kind == "binary") {
      values <- sort(unique(observed[, 1L]))
      level_design <- matrix(values, ncol = 1L,
                             dimnames = list(NULL, columns))
      level_labels <- as.character(values)
    } else if (kind == "categorical") {
      keep <- !duplicated(as.data.frame(observed))
      level_design <- observed[keep, , drop = FALSE]
      level_labels <- apply(level_design, 1L, paste, collapse = ":")
    }

    specs[[term]] <- list(
      term = term,
      kind = kind,
      discrete = kind != "continuous",
      columns = columns,
      width = length(columns),
      level_labels = level_labels,
      level_design = level_design
    )
  }

  specs
}

#' @describeIn mfpi Default method accepting a numeric matrix or data frame
#'   \code{x} and response vector \code{y}.
#'
#' @method mfpi default
#' @export
mfpi.default <- function(
    x,
    y,
    group_var         = NULL,
    interaction_vars  = NULL,
    interaction_forms = NULL,
    flex              = c("flex3", "flex1", "flex2", "flex4"),
    p_interact        = 0.05,
    min_improvement   = NULL,
    include_group_var = FALSE,
    show_models       = FALSE,
    weights           = NULL,
    offset            = NULL,
    cycles            = 10,
    scale             = NULL,
    shift             = NULL,
    df                = 4,
    center            = TRUE,
    subset            = NULL,
    family            = "gaussian",
    criterion         = c("pvalue", "aic", "bic"),
    select            = 0.05,
    alpha             = 0.05,
    keep              = NULL,
    force_max_fp_vars = NULL,
    xorder            = c("ascending", "descending", "original"),
    powers            = NULL,
    ties              = c("efron", "breslow"),
    strata            = NULL,
    nocenter          = c(-1, 0, 1),
    acd_vars          = NULL,
    zero_vars         = NULL,
    catzero_vars      = NULL,
    spike_vars        = NULL,
    min_saz_prop      = 0.10,
    ftest             = FALSE,
    control           = NULL,
    winsorize         = FALSE,
    winsorize_probs   = c(0.01, 0.99),
    center_type       = c("grand", "group"),
    p_adjust_method   = "none",
    verbose           = TRUE,
    digits            = 3,
    term_groups       = NULL,
    fitter            = c("base", "fastglm"),
    ...
) {
  # mfpi.default() validates and normalizes every argument, resolves the
  # zero/catzero/spike/acd flags and shift/scale/df settings per predictor
  # (mirroring mfp2.default()), then delegates to fit_mfpi() to select the
  # adjustment model and test each interaction_vars interaction in turn.
  cl <- match.call()

  dots <- match.call(expand.dots = FALSE)$...
  mfpi_check_unused_dots(
    dots = dots,
    context = "`mfpi.default()`"
  )

  # Step 1: Resolve multiple-choice arguments and the model family ------------
  # Match enumerated arguments -------------------------------------------------
  xorder      <- match.arg(xorder)
  criterion   <- match.arg(criterion)
  flex        <- match.arg(flex)
  ties        <- resolve_mfp_ties(ties)
  center_type <- match.arg(center_type)

  # Resolve family -------------------------------------------------------------
  family_info <- normalize_family_argument(
    family,
    family_arg = deparse(substitute(family))
  )

  family        <- family_info$family
  family_string <- family_info$family_string
  validate_mfpi_family(family_string)
  fitter <- resolve_fitter(fitter, family_string)

  # Ordinal proportional-odds models via rms::orm are supported. The k-1

  # cut-point intercepts are handled transparently by fit_ordinal(), which
  # returns rank = p (slopes only) and df = n_int + p. The rank validation
  # in validate_mfpi_fitted_rank() accounts for this convention.

  # Step 2: Prepare x and group_var, then run basic sanity checks -------------
  # Formula dispatch can attach a full-row preprocessing design to the
  # subset-specific fitting matrix. Preserve that private attribute while
  # prepare_mfpi_default_x() normalizes x/group_var, then extract it into an
  # explicit preprocess_x matrix. Direct matrix calls use x for both roles.
  preprocess_attr <- attr(x, "mfp2_preprocess_x", exact = TRUE)
  formula_generated_settings <- !is.null(preprocess_attr)
  attr(x, "mfp2_preprocess_x") <- NULL

  # Prepare x
  prepared_x <- prepare_mfpi_default_x(
    x = x,
    group_var = group_var
  )

  x <- prepared_x$x
  group_var <- prepared_x$group_var
  group_input_levels <- prepared_x$group_input_levels
  group_levels_original <- prepared_x$group_levels_original

  if (!is.null(preprocess_attr)) {
    attr(x, "mfp2_preprocess_x") <- preprocess_attr
  }
  preprocessing <- extract_preprocess_matrix(x)
  x <- preprocessing$x
  preprocess_x <- preprocessing$preprocess_x
  x_input <- x

  vnames <- colnames(x)

  # Keep a complete conceptual-term lookup throughout MFPI. Unmentioned
  # columns remain singleton terms; grouped entries represent one categorical
  # adjustment term generated by multiple design-matrix columns.
  term_to_columns <- normalize_term_groups(vnames, term_groups)
  # Expand user-facing conceptual term names to their raw matrix columns.
  # This lets matrix users refer either to a complete grouped term or to an
  # ordinary singleton column in keep/flag arguments.
  expand_term_names <- function(values) {
    if (is.null(values)) return(NULL)
    unique(unlist(lapply(values, function(value) {
      if (value %in% names(term_to_columns)) {
        term_to_columns[[value]]
      } else {
        value
      }
    }), use.names = FALSE))
  }

  if (anyNA(x)) {
    stop(
      "! `x` must not contain missing values (NA).\n",
      "i Remove or impute missing data before calling `mfpi()`.",
      call. = FALSE
    )
  }

  if (!is.numeric(x) || any(!is.finite(x))) {
    stop(
      "! `x` must be numeric after preprocessing and contain finite values.",
      call. = FALSE
    )
  }

  # Enforce the exact-zero domain as soon as the original unshifted matrix and
  # term mapping are available. This precedes group/continuous-variable
  # eligibility checks as well as all option cascading and preprocessing, so a
  # later validation or reset cannot mask negative requested covariates.
  validate_zero_option_covariates(
    x = preprocess_x,
    variables = intersect(
      unique(c(
        expand_term_names(zero_vars),
        expand_term_names(catzero_vars),
        expand_term_names(spike_vars)
      )),
      vnames
    )
  )

  # Basic dimension checks on x ------------------------------------------------
  np    <- dim(x)
  if (is.null(np)) {
    stop(
      "! `x` must be a matrix with at least one row and one column.\n",
      "i `dim(x)` returned NULL.",
      call. = FALSE
    )
  }

  nobs  <- as.integer(np[1L])
  nvars <- as.integer(np[2L])
  # Step 3: Validate group_var, p_adjust_method, and interaction_vars ----------------
  # Validate group_var ---------------------------------------------------------
  if (is.null(group_var)) {
    stop(
      "! `group_var` is NULL.\n",
      "i Provide the name of the categorical grouping variable (e.g. treatment indicator).",
      call. = FALSE
    )
  }

  if (length(group_var) != 1L) {
    stop(
      paste0("! `group_var` must name exactly one variable; ",
             length(group_var), " were supplied."),
      call. = FALSE
    )
  }

  if (!group_var %in% vnames) {
    stop(paste0("! `group_var = '", group_var, "'` is not a column of `x`."),
         call. = FALSE)
  }

  if (group_var %in% interaction_vars) {
    stop(
      paste0("! `group_var = '", group_var,
             "'` must not also appear in `interaction_vars`."),
      call. = FALSE
    )
  }

  nlev <- length(unique(x[, group_var]))
  if (nlev < 2L) {
    stop(
      paste0("! `group_var` must have at least two distinct values; found ", nlev, "."),
      call. = FALSE
    )
  }

  if (!is.character(p_adjust_method) || length(p_adjust_method) != 1L) {
    stop("! `p_adjust_method` must be a single character string.", call. = FALSE)
  }

  if (!p_adjust_method %in% stats::p.adjust.methods) {
    stop(
      "! Invalid `p_adjust_method`: ", p_adjust_method, ".\n",
      "i Supported methods are: ",
      paste(stats::p.adjust.methods, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  # Validate interaction_vars ---------------------------------------------------------
  if (is.null(interaction_vars) || !is.character(interaction_vars) || length(interaction_vars) == 0L) {
    stop(
      paste0(
        "! `interaction_vars` must be a non-empty character vector naming ",
        "continuous, binary, or categorical variables to test for interaction ",
        "with `group_var = '", group_var, "'`."
      ),
      call. = FALSE
    )
  }

  valid_interaction_terms <- unique(c(vnames, names(term_to_columns)))
  missing_cont <- setdiff(interaction_vars, valid_interaction_terms)
  if (length(missing_cont) > 0L) {
    stop(
      paste0("! The following variables in `interaction_vars` are not predictors in `x`: ",
             paste(missing_cont, collapse = ", "), "."),
      call. = FALSE
    )
  }

  # Grouped interaction variables must be requested by their conceptual term
  # name so the complete fitted design block is tested jointly. The grouping
  # variable itself cannot participate in a mapped predictor term.
  grouped_terms <- names(term_to_columns)[mapped_term_flags(term_to_columns)]
  grouped_columns <- unlist(term_to_columns[grouped_terms], use.names = FALSE)

  if (group_var %in% grouped_terms || group_var %in% grouped_columns) {
    stop("! `group_var` must not be included in `term_groups`.", call. = FALSE)
  }

  grouped_members <- interaction_vars[
    interaction_vars %in% grouped_columns & !interaction_vars %in% grouped_terms
  ]
  if (length(grouped_members) > 0L) {
    owners <- vapply(grouped_members, function(column) {
      names(term_to_columns)[vapply(
        term_to_columns,
        function(columns) column %in% columns,
        logical(1L)
      )][1L]
    }, character(1L))
    stop(
      paste0(
        "! Grouped categorical interaction columns must be requested by their ",
        "complete term name: ",
        paste0(grouped_members, " -> ", owners, collapse = ", "), "."
      ),
      call. = FALSE
    )
  }

  # Resolve forms after identifying complete continuous and discrete terms.
  interaction_specs <- mfpi_interaction_specs(
    x = preprocess_x,
    interaction_vars = interaction_vars,
    term_to_columns = term_to_columns
  )
  interaction_forms <- resolve_mfpi_interaction_forms(
    interaction_forms = interaction_forms,
    interaction_specs = interaction_specs
  )
  for (term in names(interaction_specs)) {
    interaction_specs[[term]]$form <- interaction_forms[[term]]
  }
  interaction_columns <- unique(unlist(
    lapply(interaction_specs, `[[`, "columns"),
    use.names = FALSE
  ))
  discrete_interaction_columns <- unique(unlist(
    lapply(Filter(function(spec) isTRUE(spec$discrete), interaction_specs),
           `[[`, "columns"),
    use.names = FALSE
  ))
  continuous_interaction_vars <- names(Filter(
    function(spec) identical(spec$kind, "continuous"),
    interaction_specs
  ))

  # Near-categorical continuous variables remain eligible, but the warning is
  # relevant only to FP searches. Binary and factor terms use the discrete
  # linear path and therefore do not trigger it.
  n_unique_continuous <- vapply(continuous_interaction_vars, function(v) {
    length(unique(preprocess_x[, interaction_specs[[v]]$columns, drop = TRUE]))
  }, integer(1L))
  few_unique <- continuous_interaction_vars[
    n_unique_continuous > 2L & n_unique_continuous <= 5L &
      interaction_forms[continuous_interaction_vars] != "linear"
  ]

  if (length(few_unique) > 0L) {
    warning(
      paste0(
        "! The following variable(s) in `interaction_vars` have 5 or fewer unique ",
        "values: ", paste(few_unique, collapse = ", "), ".\n",
        "  FP transformation may be unreliable for near-categorical variables. ",
        "Consider whether these variables are truly continuous."
      ),
      call. = FALSE
    )
  }

  # Step 4: Validate the response, survival inputs, weights, and offset --------
  # Validate response y --------------------------------------------------------
  validate_family_response(
    y             = y,
    family_string = family_string,
    nobs          = nobs
  )

  # Validate survival-model stratification ------------------------------------
  if (mfp2_family_is_survival(family_string) && !is.null(strata)) {
    strata_len <- if (is.vector(strata) || is.factor(strata)) {
      length(strata)
    } else {
      NROW(strata)
    }

    if (strata_len != nobs) {
      stop(
        paste0(
          "! Length of `strata` (", strata_len, ") must equal the number of ",
          "rows in `x` (", nobs, ")."
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

  # Validate offset -------------------------------------------------------------
  if (!is.null(offset)) {
    offset_rows <- if (is.matrix(offset)) nrow(offset) else length(offset)
    if (!is.numeric(offset) || offset_rows != nobs || anyNA(offset) ||
        any(!is.finite(offset))) {
      stop(
        "! `offset` must be finite numeric data with one row per observation.",
        call. = FALSE
      )
    }
    if (identical(family_string, "multinomial")) {
      c_classes <- mfp2_multinomial_n_classes(y)
      if (!is.matrix(offset) || !ncol(offset) %in% c(c_classes - 1L, c_classes)) {
        stop("! Multinomial `offset` must have C or C - 1 columns.", call. = FALSE)
      }
    } else if (is.matrix(offset)) {
      stop("! `offset` must be a vector for this family.", call. = FALSE)
    }
  }

  # Step 5: Validate alpha/select/df, keep, and force_max_fp_vars -------------
  # Matrix inputs use unnamed global scalars or named partial overrides. All
  # named values are matched to columns rather than assigned positionally.
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
  if (df_setting$global_scalar && length(df) > 0L) {
    df <- expand_scalar_df_for_mapped_terms(
      df = unname(df[[1L]]),
      vnames = vnames,
      term_to_columns = term_to_columns
    )
  }

  validate_probability_vector(alpha, "alpha", nvars)
  validate_probability_vector(select, "select", nvars)

  # Validate keep. A grouped categorical term may be named directly; it is
  # converted to one conceptual keep entry before adjustment-model fitting.
  valid_term_names <- unique(c(vnames, names(term_to_columns)))
  if (!is.null(keep) && !all(keep %in% valid_term_names)) {
    warning(
      "i Some entries in `keep` are not predictor columns or conceptual terms; ",
      "continuing with the valid entries.",
      call. = FALSE
    )
    keep <- intersect(keep, valid_term_names)
  }

  # Validate force_max_fp_vars ------------------------------------------------
  # Convert character vector of variable names to named logical vector
  # (the internal format expected by fit_mfpi and flex functions).
  force_max_fp <- setNames(rep(FALSE, nvars), vnames)
  if (!is.null(force_max_fp_vars)) {
    if (!is.character(force_max_fp_vars)) {
      stop("! `force_max_fp_vars` must be a character vector or NULL.",
           call. = FALSE)
    }
    bad_names <- setdiff(force_max_fp_vars, valid_term_names)
    if (length(bad_names) > 0L) {
      warning("i Some variables in `force_max_fp_vars` are not columns of `x`; ",
              "they will be ignored: ", paste(bad_names, collapse = ", "), ".",
              call. = FALSE)
    }
    force_max_fp_vars <- intersect(force_max_fp_vars, valid_term_names)
    force_max_fp_columns <- expand_term_names(force_max_fp_vars)
    if (length(force_max_fp_columns) > 0L) {
      force_max_fp[force_max_fp_columns] <- TRUE
    }
  }

  # Step 6: Validate zero_vars / catzero_vars / spike_vars / acd_vars ---------
  # Validate zero_vars / catzero_vars / spike_vars -----------------------------
  if (!is.null(zero_vars) && !all(zero_vars %in% valid_term_names)) {
    warning("i Some variables in `zero_vars` are not columns of `x`; they will be ignored.",
            call. = FALSE)
  }

  if (!is.null(catzero_vars) && !all(catzero_vars %in% valid_term_names)) {
    warning("i Some variables in `catzero_vars` are not columns of `x`; they will be ignored.",
            call. = FALSE)
  }

  if (!is.null(spike_vars) && !all(spike_vars %in% valid_term_names)) {
    warning("i Some variables in `spike_vars` are not columns of `x`; they will be ignored.",
            call. = FALSE)
  }

  zero_vars_in    <- intersect(expand_term_names(zero_vars),    vnames)
  catzero_vars_in <- intersect(expand_term_names(catzero_vars), vnames)
  spike_vars_in   <- intersect(expand_term_names(spike_vars),   vnames)

  # For interaction_vars, catzero/spike will be dropped later, so do not treat
  # zero + catzero overlap on interaction_vars as fatal.
  overlap <- intersect(zero_vars_in, setdiff(catzero_vars_in, interaction_columns))

  if (length(overlap) > 0L) {
    stop(
      paste0(
        "! The following variables appear in both `zero_vars` and `catzero_vars`: ",
        paste(overlap, collapse = ", "),
        ". Use only one option per variable."
      ),
      call. = FALSE
    )
  }

  zero_vars    <- zero_vars_in
  catzero_vars <- catzero_vars_in
  spike_vars   <- spike_vars_in

  # Validate acd_vars ----------------------------------------------------------
  if (!is.null(acd_vars) && !all(acd_vars %in% valid_term_names)) {
    warning("i Some variables in `acd_vars` are not columns of `x`; they will be ignored.",
            call. = FALSE)
  }

  # Step 7: Validate shift/scale/center/df, SAZ proportion, and ftest ---------
  # Validate shift and scale ---------------------------------------------------
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

  # Validate and normalize center ----------------------------------------------
  # Developer note: mfpi.formula() may pass a complete named logical vector
  # after reading fp()/fp2() attributes. Direct default-method calls may instead
  # use an unnamed global scalar or a named partial vector. In all named cases,
  # matching is by colnames(x), never by predictor position.
  center <- normalize_named_logical_setting(
    value = center,
    column_names = vnames,
    default = TRUE,
    argument_name = "center"
  )

  # Validate df ----------------------------------------------------------------
  if (any(df <= 0L)) {
    stop("! All values of `df` must be positive (1 for linear, 2m for FP degree m).",
         call. = FALSE)
  }

  invalid_df <- df != 1L & df %% 2L != 0L
  if (any(invalid_df)) {
    stop(
      paste0(
        "! Each element of `df` must be 1 (linear) or an even number 2m. ",
        "Invalid variable(s): ",
        paste(names(df)[invalid_df], collapse = ", "), "."
      ),
      call. = FALSE
    )
  }

  # Validate spike-at-zero component proportion -------------------------------
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

  # Warn if ftest is incompatible with family ------------------------------
  if (ftest && family_string != "gaussian") {
    warning(
      paste0("i `ftest = TRUE` is only applicable to Gaussian models; ",
             "reverting to chi-square test for family = '", family_string, "'."),
      call. = FALSE
    )
    ftest <- FALSE
  }

  # Validate and build powers list ----------------------------------------------
  # `powers` defines candidate base-power sets, not selected FP power vectors.
  # Duplicate candidate values are removed. Repeated selected FP powers are
  # generated later by replacement when fitting FP2 or higher-degree models.
  # Step 8: Validate the powers list, subset, and interaction selection
  # thresholds (p_interact / min_improvement) ---------------------------------
  # (p_adjust_method was already fully validated in Step 3, including the
  # check against stats::p.adjust.methods, so it is not re-validated here.)
  power_list <- validate_fp_power_list(
    powers = powers,
    vnames = vnames,
    arg_name = "powers"
  )

  # Validate subset ------------------------------------------------------------
  if (!is.null(subset)) {
    if (is.logical(subset)) {
      if (length(subset) != nobs || anyNA(subset)) {
        stop(
          "! Logical `subset` must have one TRUE/FALSE value per observation and contain no NA.",
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
          "! Numeric `subset` must contain valid positive row indices within the range of `x`.",
          call. = FALSE
        )
      }

      subset <- as.integer(subset)

    } else {
      stop(
        "! `subset` must be either a logical vector or a numeric/integer vector of row indices.",
        call. = FALSE
      )
    }

    if (anyDuplicated(subset)) {
      stop(
        "! `subset` must not contain duplicated row indices.",
        call. = FALSE
      )
    }

    if (length(subset) < 5L) {
      stop(
        paste0(
          "! After subsetting, only ", length(subset),
          " observations remain; at least 5 are required."
        ),
        call. = FALSE
      )
    }
  }

  # Apply one strictly-positive observation-weight contract to every family.
  # Zero weights are rejected even when a fitter accepts them because they can
  # make likelihood-based MFP/MFPI comparisons undefined for some families.
  validate_model_weights(
    weights = weights,
    nobs = nobs
  )

  # Validate interaction-selection thresholds -----------------------------------
  if (!is.numeric(p_interact) ||
      length(p_interact) != 1L ||
      anyNA(p_interact) ||
      !is.finite(p_interact) ||
      p_interact <= 0 ||
      p_interact > 1) {
    stop(
      "! `p_interact` must be a single finite numeric value in (0, 1].",
      call. = FALSE
    )
  }

  if (!is.null(min_improvement)) {
    if (!is.numeric(min_improvement) ||
        length(min_improvement) != 1L ||
        anyNA(min_improvement) ||
        !is.finite(min_improvement) ||
        min_improvement <= 0) {
      stop(
        "! `min_improvement` must be NULL or a single finite positive numeric value.",
        call. = FALSE
      )
    }
  }

  # Step 9: Apply defaults and expand scalar options to named per-variable
  # vectors, then resolve the force_max_fp/alpha interaction ------------------
  # Set default min_improvement -------------------------------------------------
  if (is.null(min_improvement)) {
    min_improvement <- switch(criterion, pvalue = p_interact, aic = 2, bic = 2)
  }

  # Default weights (all observations equally weighted) and offset (none);
  # has_offset is recorded before defaulting so the fitted object can report
  # whether the user actually supplied an offset.
  if (is.null(weights)) weights <- rep.int(1, nobs)
  has_offset <- !is.null(offset)

  if (is.null(offset))  offset  <- rep.int(0, nobs)

  # select, alpha, and center were normalized above to complete named vectors
  # ordered exactly like colnames(x).

  # For adjustment-model variables under p-value selection, force_max_fp
  # requires both retention of the variable and acceptance of its most complex
  # permitted FP degree.
  if (criterion == "pvalue" && any(force_max_fp)) {
    select[force_max_fp] <- 1
    alpha[force_max_fp] <- 1
  }

  # Automatic shifts are deliberately estimated later, after SAZ eligibility and
  # the zero/catzero/spike cascade have reached their final retained state.
  # Retained zero-handled variables need shift = 0; variables reset from SAZ to
  # ordinary FP need ordinary positivity shifting.
  # `shift` and `scale` were normalized above and are already aligned by name.
  # Missing internal formula sentinels are resolved after the SAZ cascade and
  # after shifting, respectively.

  # Resolve and validate fitting controls once for both MFPI stages and every
  # nested MFP candidate fit.
  control <- normalize_fit_control(
    control = control,
    family_string = family_string,
    fitter = fitter
  )

  # Retain valid raw-column or conceptual-term entries. preprocess_data()
  # converts these to the term names expected by the grouped MFP core.
  keep <- intersect(keep, valid_term_names)

  # Step 10: Build zero/catzero/spike flags, suppress them for interaction_vars, and
  # resolve spike-at-zero eligibility for the remaining adjustment variables -
  # Process zero_vars / catzero_vars / spike_vars to named logical vectors -----
  zero_flag    <- setNames(rep(FALSE, nvars), vnames)
  catzero_flag <- setNames(rep(FALSE, nvars), vnames)
  spike_flag   <- setNames(rep(FALSE, nvars), vnames)

  if (length(zero_vars)    > 0L) zero_flag[zero_vars]       <- TRUE
  if (length(catzero_vars) > 0L) catzero_flag[catzero_vars] <- TRUE
  if (length(spike_vars)   > 0L) spike_flag[spike_vars]     <- TRUE

  # interaction_vars may keep explicit zero handling, but cannot use catzero/spike.
  # Developer note: this cleanup must happen BEFORE the cascade below. Otherwise
  # spike_vars/catzero_vars would imply zero_vars for interaction_vars, which is not the
  # intended MFPI stage-2 contract. Only explicit zero_vars should survive for
  # interaction variables.
  spike_in_cont   <- intersect(names(spike_flag)[spike_flag], interaction_columns)
  catzero_in_cont <- intersect(names(catzero_flag)[catzero_flag], interaction_columns)

  if (length(spike_in_cont) > 0L) {
    warning(
      "The following `interaction_vars` were also marked as `spike_vars`: ",
      paste(spike_in_cont, collapse = ", "),
      ". Spike-at-zero handling is not supported for MFPI interaction variables; ",
      "`spike` and the implied `catzero` indicator are set to FALSE for these variables. ",
      "Explicit `zero_vars` settings are preserved.",
      call. = FALSE
    )

    spike_flag[spike_in_cont]   <- FALSE
    catzero_flag[spike_in_cont] <- FALSE
  }

  if (length(catzero_in_cont) > 0L) {
    warning(
      "The following `interaction_vars` were also marked as `catzero_vars`: ",
      paste(catzero_in_cont, collapse = ", "),
      ". Structural-zero binary indicators are not supported for MFPI interaction variables; ",
      "`catzero` is set to FALSE for these variables. ",
      "Explicit `zero_vars` settings are preserved.",
      call. = FALSE
    )

    catzero_flag[catzero_in_cont] <- FALSE
  }

  # Resolve spike-at-zero eligibility before the cascade and before shift/scale ---
  # This must happen after unsupported cont_var spike/catzero flags have been
  # removed, but before the spike -> catzero -> zero cascade is applied.
  #
  # resolve_saz_eligibility() needs the user's explicit zero/catzero choices.
  # Therefore, do not call it after catzero_flag[spike_flag] <- TRUE or
  # zero_flag[catzero_flag] <- TRUE, because that would lose the distinction
  # between user-specified flags and spike-implied flags.
  if (any(spike_flag)) {
    saz_flags <- resolve_saz_eligibility(
      x                      = preprocess_x,
      spike                  = spike_flag,
      catzero                = catzero_flag,
      zero                   = zero_flag,
      min_saz_prop = min_saz_prop
    )

    spike_flag   <- saz_flags$spike
    catzero_flag <- saz_flags$catzero
    zero_flag    <- saz_flags$zero
  }

  # Enforce the cascade after SAZ eligibility has been resolved.
  # Retained spike variables imply catzero, and retained catzero variables imply
  # zero. Variables whose spike flag was reset keep only the user's explicit
  # zero/catzero choices.
  catzero_flag[spike_flag] <- TRUE
  zero_flag[catzero_flag]  <- TRUE

  # Step 11: Reset zero/catzero/spike for all-positive or binary variables,
  # and build the acd_vars flag (suppressed for interaction_vars) -------------------
  # Reset zero/catzero/spike for variables that contain only positive values -----
  if (any(zero_flag | catzero_flag | spike_flag)) {
    vars_pos <- names(zero_flag)[zero_flag | catzero_flag | spike_flag]

    xpos <- preprocess_x[, vars_pos, drop = FALSE]
    all_positive <- vars_pos[colSums(xpos == 0) == 0L]

    if (length(all_positive) > 0L) {
      warning(
        paste0(
          "i The following variables in `zero_vars`, `catzero_vars`, or ",
          "`spike_vars` contain only positive values and have been reset to ",
          "standard processing: ",
          paste(all_positive, collapse = ", "), "."
        ),
        call. = FALSE
      )

      zero_flag[all_positive]    <- FALSE
      catzero_flag[all_positive] <- FALSE
      spike_flag[all_positive]   <- FALSE
    }
  }

  # Reset zero/catzero for binary variables ------------------------------------
  binary_vars  <- apply(preprocess_x, 2L,
                        function(col) length(unique(col[!is.na(col)])) == 2L)
  binary_names <- vnames[binary_vars]
  if (length(binary_names) > 0L) {
    # Binary predictors, including factor dummy columns, must remain on their
    # original two-level scale. Automatic preprocessing already returns scale = 1
    # for them, but an explicitly supplied scale would otherwise override that
    # safeguard and rescale the columns before the MFPI adjustment fit.
    scale[binary_names] <- 1

    reset_vars <- intersect(
      binary_names,
      Reduce(
        union,
        list(
          names(zero_flag)[zero_flag],
          names(catzero_flag)[catzero_flag],
          names(spike_flag)[spike_flag]
        )
      )
    )
    if (length(reset_vars) > 0L) {
      warning(
        paste0(
          "i The following binary variables were marked in `zero_vars`, ",
          "`catzero_vars`, or `spike_vars` but are binary; resetting to standard ",
          "processing: ",
          paste(reset_vars, collapse = ", "), "."
        ),
        call. = FALSE
      )

      zero_flag[reset_vars]    <- FALSE
      catzero_flag[reset_vars] <- FALSE
      spike_flag[reset_vars]   <- FALSE
    }
  }

  acd_flag <- setNames(rep(FALSE, nvars), vnames)

  if (!is.null(acd_vars)) {
    # Keep only variables that exist in `x`; unknown variables were warned about
    # in the earlier validation block.
    acd_vars <- unique(intersect(expand_term_names(acd_vars), vnames))

    # Developer note:
    # ACD transformations are disabled for MFPI interaction variables. The
    # interaction variable's functional form is handled by the flex0-flex4
    # search, and applying ACD upstream would change the interaction scale and
    # break reconstruction of f_j(x).
    acd_in_cont <- intersect(acd_vars, interaction_columns)

    if (length(acd_in_cont) > 0L) {
      warning(
        "ACD transformation is disabled for interaction_vars: ",
        paste(acd_in_cont, collapse = ", "),
        ".",
        call. = FALSE
      )

      acd_vars <- setdiff(acd_vars, interaction_columns)
    }

    acd_flag[acd_vars] <- TRUE
  }

  # Step 12: Assign per-variable df, then reset shift for zero/catzero/linear
  # variables (their shift is meaningless since only the positive part, or no
  # nonlinear transform at all, is ever fitted) -------------------------------
  # Set degrees of freedom per variable ----------------------------------------
  # df was normalized above to a complete named vector. Omitted named
  # entries already contain ordinary defaults, with explicitly mapped design
  # blocks receiving their structural df = 1 fallback.
  df_default <- df
  df_list <- setNames(assign_df(x = preprocess_x, df_default = df_default), vnames)

  # Discrete interaction terms are ordinary, uncentred design columns. They
  # must never inherit continuous-variable preprocessing or FP settings from
  # global arguments. Keeping the original indicator coding also makes their
  # coefficients and predictions directly comparable with glm(), lm(), and
  # nnet::multinom() treatment-contrast fits.
  if (length(discrete_interaction_columns) > 0L) {
    df_list[discrete_interaction_columns] <- 1L
    center[discrete_interaction_columns] <- FALSE
    shift[discrete_interaction_columns] <- 0
    scale[discrete_interaction_columns] <- 1
    power_list[discrete_interaction_columns] <- rep(
      list(1), length(discrete_interaction_columns)
    )
    zero_flag[discrete_interaction_columns] <- FALSE
    catzero_flag[discrete_interaction_columns] <- FALSE
    spike_flag[discrete_interaction_columns] <- FALSE
    acd_flag[discrete_interaction_columns] <- FALSE
    force_max_fp[discrete_interaction_columns] <- FALSE
  }

  # Multi-column adjustment terms are fixed linear design blocks. They cannot
  # use continuous-variable extensions or structural-zero representations.
  validate_grouped_term_setting(
    term_to_columns, acd_flag, "acd_vars", function(v) !v, "FALSE"
  )
  validate_grouped_term_setting(
    term_to_columns, zero_flag, "zero_vars", function(v) !v, "FALSE"
  )
  validate_grouped_term_setting(
    term_to_columns, catzero_flag, "catzero_vars", function(v) !v, "FALSE"
  )
  validate_grouped_term_setting(
    term_to_columns, spike_flag, "spike_vars", function(v) !v, "FALSE"
  )
  validate_grouped_term_setting(
    term_to_columns, force_max_fp, "force_max_fp_vars", function(v) !v, "FALSE"
  )
  validate_grouped_term_setting(
    term_to_columns, df_list, "df", function(v) v == 1, "1 (linear-only)"
  )

  # Reset shift to 0 for zero/catzero variables and linear terms ---------------
  # Only the positive part of these variables is FP-transformed so no shift
  # is needed (matching the logic in mfp2.default()).
  # Step 13: Shift and scale x (excluding group_var), then verify positivity --
  shift_to_zero <- rep(FALSE, nvars)
  names(shift_to_zero) <- vnames
  shift_to_zero[names(zero_flag)[zero_flag]]       <- TRUE
  shift_to_zero[names(catzero_flag)[catzero_flag]] <- TRUE
  shift_to_zero[vnames[df_list == 1L]]             <- TRUE
  shift[shift_to_zero] <- 0

  # At this point  the code has already handled zero_flag, catzero_flag,
  # spike_flag, etc. so a good place to estimate shifting factors
  # shift = 0 for retained zero/catzero variables and df = 1 variables
  # automatic find_shift_factor() only for remaining NA shifts
  shift_missing <- is.na(shift)
  if (any(shift_missing)) {
    shift[shift_missing] <- apply(
      preprocess_x[, shift_missing, drop = FALSE],
      2L,
      find_shift_factor
    )
  }

  # `group_var` is categorical/index metadata, not an FP predictor.
  # It must remain on its input coding until preprocess_data() remaps it to
  # internal 0, 1, ..., K - 1 levels.
  if (group_var %in% names(shift)) {
    shift[group_var] <- 0
  }

  # Apply shift and scale to x -------------------------------------------------
  preprocess_x_shifted <- sweep(preprocess_x, 2L, shift, "+")
  x <- sweep(x, 2L, shift, "+")

  # `group_var` is categorical/index metadata, not an FP predictor. It must not be
  # automatically scaled because preprocess_data() later remaps it to internal
  # 0, 1, ..., K - 1 group codes. Keeping scale[group_var] = 1 also ensures that
  # stored prediction metadata does not imply a raw continuous-variable scaling
  # operation for the grouping column.
  if (group_var %in% names(scale)) {
    scale[group_var] <- 1
  }

  # Compute scale factors on shifted x, matching standalone mfp2. User-supplied
  # scale values are left untouched; only missing scales are estimated.
  scale_missing <- is.na(scale)

  if (any(scale_missing)) {
    scale[scale_missing] <- apply(
      preprocess_x_shifted[, scale_missing, drop = FALSE],
      2L,
      find_scale_factor
    )
  }

  # Positivity check for variables that undergo nonlinear FP transformation ----
  nonlinear_names <- vnames[df_list != 1L]
  all_zero_names  <- union(names(zero_flag)[zero_flag],
                           names(catzero_flag)[catzero_flag])
  check_names     <- setdiff(nonlinear_names, all_zero_names)

  if (length(check_names) > 0L) {
    xcheck  <- preprocess_x_shifted[, check_names, drop = FALSE]
    neg_idx <- which(colSums(xcheck <= 0, na.rm = TRUE) > 0L)
    if (length(neg_idx) > 0L) {
      stop(
        paste0(
          "! Shifting factors are insufficient to ensure positive values for FP ",
          "transformation.\n",
          "i Problematic variables: ",
          paste(colnames(xcheck)[neg_idx], collapse = ", "),
          ".\ni Consider increasing the shift values for these variables."
        ),
        call. = FALSE
      )
    }
  }

  # Scale the x based on the estimated scaling factors
  x <- sweep(x, 2L, scale, "/")

  # Step 14: Prepare survival stratification, apply `subset`, and re-validate
  # group_var against the (possibly reduced) subsetted data -------------------
  # Prepare survival stratification -------------------------------------------
  # Normalize all strata to one factor before MFP/MFPI candidate fitting.
  # Multiple columns are combined as in survival::coxph(), while a single
  # vector is treated as categorical labels regardless of storage mode.
  # Integer conversion belongs only at the survival::coxph.fit() boundary.
  strata_keep <- strata

  if (mfp2_family_is_survival(family_string) && !is.null(strata_keep)) {
    strata_keep <- normalize_cox_strata(strata_keep, nobs = nobs)
  }

  # Apply subset ---------------------------------------------------------------
  # This is the direct matrix path. Keep the user-supplied numeric columns, but
  # reject mapped blocks that lose within-block rank because MFPI has no factor
  # levels or contrast function from which to regenerate them. Singleton columns
  # are checked separately for complete loss of variation.
  if (!is.null(subset)) {
    validate_grouped_subset_rank(
      x_full = x_input,
      x_fit = x_input[subset, , drop = FALSE],
      term_to_columns = term_to_columns
    )
    x       <- x[subset, , drop = FALSE]

    y <- if (is.matrix(y)) {
      y[subset, , drop = FALSE]
    } else {
      y[subset]
    }
    weights <- weights[subset]
    offset  <- if (is.matrix(offset)) offset[subset, , drop = FALSE] else offset[subset]
    if (!is.null(strata_keep)) {
      # Subsetting can remove complete strata, so drop unused levels after the
      # row restriction while preserving row alignment.
      strata_keep <- droplevels(strata_keep[subset])
    }
  }

  # zero_vars applies the interaction function only to x > 0. Validate this
  # requirement after the optional subset has been applied, because a variable
  # can contain positive values in the full data yet lose all of them in the
  # actual fitting population. This shared boundary gives every flexibility
  # method the same variable-specific error instead of exposing flex2/flex4
  # transformation assertions to users.
  zero_interaction_vars <- intersect(
    interaction_vars,
    names(zero_flag)[zero_flag]
  )

  if (length(zero_interaction_vars) > 0L) {
    has_positive <- vapply(
      zero_interaction_vars,
      function(variable) any(x[, variable] > 0),
      logical(1L)
    )

    if (any(!has_positive)) {
      stop(
        "Zero-handled variables in `interaction_vars` must contain at least one ",
        "positive value in the fitting data after subsetting. Problematic ",
        "variables: ",
        paste(zero_interaction_vars[!has_positive], collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }

  # Run the general post-subset variation check after the zero-specific support
  # check so an all-zero zero-handled continuous variable receives the more
  # informative missing-positive-component error. Formula calls arrive with an
  # already-subset fitting matrix and carry `formula_generated_settings = TRUE`.
  if (!is.null(subset) || formula_generated_settings) {
    validate_subset_predictor_variation(x, exclude = group_var)
  }

  # Refresh group-level metadata after subsetting. The metadata are created
  # before subset is applied, so a subset can remove one or more group levels.
  # Passing stale levels into preprocess_data() can create all-zero group dummies.
  observed_group_levels <- sort(unique(as.numeric(x[, group_var])))

  keep_group_levels <- group_input_levels %in% observed_group_levels

  group_input_levels <- group_input_levels[keep_group_levels]
  group_levels_original <- group_levels_original[keep_group_levels]

  # A dataset with two groups can pass initial validation, then subset can
  # remove all observations from one group. Later internals may fail inside
  # dummy creation or interaction fitting with a less informative error.

  if (length(unique(x[, group_var])) < 2L) {
    stop("`group_var` must contain at least two groups after subsetting.", call. = FALSE)
  }

  group_counts <- table(x[, group_var])

  # A singleton group cannot identify a group effect and a group-specific
  # continuous effect simultaneously. With group centering, its only basis
  # value is also centered to exact zero. This is a structural failure rather
  # than merely a small-sample concern, irrespective of the requested form.
  singleton_groups <- names(group_counts)[group_counts < 2L]
  if (length(singleton_groups) > 0L) {
    stop(
      "Each group must contain at least two observations in the fitting data ",
      "after subsetting. Singleton groups: ",
      paste(singleton_groups, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  # FP2 contributes two group-specific basis columns. Together with the group
  # effect, these require at least three observations in every group. A two-row
  # group can still identify a one-column linear/FP1 interaction, so it retains
  # the existing warning below when no FP2 interaction is requested.
  fp2_requested <- any(interaction_forms == "fp2")
  fp2_small_groups <- names(group_counts)[group_counts < 3L]
  if (fp2_requested && length(fp2_small_groups) > 0L) {
    stop(
      "FP2 interactions require at least three observations in every group ",
      "after subsetting. Problematic groups: ",
      paste(fp2_small_groups, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (identical(center_type, "group") && length(zero_interaction_vars) > 0L) {
    group_values <- x[, group_var]
    fitted_group_levels <- sort(unique(group_values))
    insufficient_positive <- character(0L)

    for (variable in zero_interaction_vars) {
      required_positive <- if (identical(interaction_forms[[variable]], "fp2")) {
        3L
      } else {
        2L
      }

      positive_counts <- vapply(
        fitted_group_levels,
        function(group) {
          as.integer(sum(group_values == group & x[, variable] > 0))
        },
        integer(1L)
      )
      names(positive_counts) <- as.character(fitted_group_levels)

      bad_groups <- names(positive_counts)[positive_counts < required_positive]
      if (length(bad_groups) > 0L) {
        insufficient_positive <- c(
          insufficient_positive,
          paste0(
            variable, " (group ", bad_groups, ": ",
            positive_counts[bad_groups], " positive; need ",
            required_positive, ")"
          )
        )
      }
    }

    if (length(insufficient_positive) > 0L) {
      stop(
        "Group-centered zero-handled interactions lack enough positive ",
        "observations to identify their requested basis: ",
        paste(insufficient_positive, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }

  small_groups <- names(group_counts)[group_counts == 2L]
  if (length(small_groups) > 0L) {
    warning(
      "Some groups have only 2 observations after subsetting: ",
      paste(small_groups, collapse = ", "),
      ". Linear or FP1 interaction fitting may be unreliable.",
      call. = FALSE
    )
  }

  # Step 15: Winsorise interaction_vars (unless excluded by zero handling) ----------
  # Winsorise interaction_vars to reduce influence of extreme values ------------------
  # Variables with final zero handling are excluded from Winsorisation so that
  # structural-zero values and the positive distribution are both preserved.
  winsorize_limits <- NULL
  if (isTRUE(winsorize) && length(continuous_interaction_vars) > 0L) {
    exclude_from_win <- names(zero_flag)[zero_flag]
    win <- winsorize_interaction_vars(
      x         = x,
      interaction_vars = continuous_interaction_vars,
      probs     = winsorize_probs,
      zero_vars = exclude_from_win
    )
    x                <- win$x
    winsorize_limits <- win$limits
    if (verbose) {
      message(sprintf(
        "Winsorising %d cont_var(s) at percentiles (%g, %g).",
        length(continuous_interaction_vars), winsorize_probs[1L], winsorize_probs[2L]
      ))

      excluded <- intersect(exclude_from_win, continuous_interaction_vars)
      if (length(excluded) > 0L) {
        message(sprintf(
          "Excluded from Winsorisation because of zero handling: %s.",
          paste(excluded, collapse = ", ")
        ))
      }

      if (!is.null(winsorize_limits)) {
        limits_text <- utils::capture.output(print(winsorize_limits))
        message(
          "Winsorisation limits:\n",
          paste(limits_text, collapse = "\n")
        )
      }
    }
  }

  # Prepare transformed survival data once before either MFPI stage. Candidate
  # fits then reuse the cached survreg response and distribution metadata.
  prepared_family <- prepare_family_for_fit(
    family = family,
    family_string = family_string,
    y = y,
    weights = weights,
    strata = strata_keep,
    offset = offset,
    has_offset = has_offset
  )
  family <- prepared_family$family
  strata_keep <- prepared_family$strata

  # Step 16: Fit the MFPI model and attach return metadata --------------------
  # Fit the MFPI model ---------------------------------------------------------
  fit <- fit_mfpi(
    x                 = x,
    y                 = y,
    family            = family,         # resolved GLM or survival family
    family_string     = family_string,  # character string for branching
    fitter            = fitter,
    weights           = weights,
    offset            = offset,
    cycles            = cycles,
    center            = center,
    criterion         = criterion,
    select            = select,
    alpha             = alpha,
    df                = df_list,
    keep              = keep,
    force_max_fp      = force_max_fp,
    xorder            = xorder,
    fp_powers         = power_list,
    ties              = ties,
    strata            = strata_keep,
    nocenter          = nocenter,
    acd_vars          = acd_flag,
    zero_vars         = zero_flag,
    catzero_vars      = catzero_flag,
    spike_vars        = spike_flag,
    min_saz_prop      = min_saz_prop,
    use_ftest         = ftest,
    control           = control,
    verbose           = verbose,
    group_var         = group_var,
    include_group_var = include_group_var,
    flex              = flex,
    interaction_vars  = interaction_vars,
    interaction_forms = interaction_forms,
    interaction_specs = interaction_specs,
    p_interact        = p_interact,
    show_models       = show_models,
    min_improvement   = min_improvement,
    digits            = digits,
    center_type       = center_type,
    scale             = scale,
    shift             = shift,
    has_offset        = has_offset,
    p_adjust_method   = p_adjust_method,
    group_input_levels     = group_input_levels,
    group_levels_original  = group_levels_original,
    term_to_columns        = term_to_columns
  )

  # Attach Winsorisation metadata for transparency in the returned object
  fit$winsorize         <- isTRUE(winsorize)
  fit$winsorize_probs   <- if (isTRUE(winsorize)) winsorize_probs else NULL
  fit$winsorize_limits  <- winsorize_limits
  fit$p_adjust_method   <- p_adjust_method
  fit$p_interact        <- p_interact
  fit$min_improvement   <- min_improvement
  fit$interaction_forms <- interaction_forms
  fit$interaction_specs <- interaction_specs

  # `fit_mfpi()` stores `x_train_internal`, the post-preprocessing training
  # matrix used for manual prediction reconstruction. It contains shifted/scaled
  # predictors and internally remapped group_var codes. Do not store or use the
  # pre-remap matrix `x` for prediction.

  fit$family_string  <- family_string
  if (identical(family_string, "multinomial")) {
    # Keep the prepared class ordering on the top-level MFPI object. The
    # public family component below is stripped to avoid duplicating the
    # response, while retained term-specific multinom fits keep the prepared
    # family needed for prediction and coefficient labelling.
    fit$mfp2_family <- family
    fit$class_levels <- family$prepared$levels
    fit$reference_class <- family$prepared$reference
    fit$n_logits <- family$prepared$n_logits
  }
  if (identical(family_string, "ordinal")) {
    # Ordinal metadata: the response levels and link function. These are
    # needed for prediction and coefficient labelling. Unlike multinomial,
    # ordinal has flat (scalar) coefficients, so it uses the standard
    # coefficient/prediction path with no per-logit fan-out.
    fit$ordinal_levels <- family$prepared$levels
    fit$ordinal_link   <- family$prepared$link
    fit$n_intercepts   <- family$prepared$n_intercepts
  }
  if (!identical(family_string, "negbin")) {
    fit$family <- mfp2_strip_prepared_family(family)
  }
  fit$zero_vars      <- zero_flag
  fit$center_type    <- center_type
  fit$interaction_vars      <- interaction_vars
  fit$term_to_columns <- fit$adjustment_term_to_columns
  fit$call <- cl

  # Compute the number of events for survival models so that
  # print.mfpi() can display it alongside n.
  if (identical(family_string, "cox")) {
    fit$nevents <- sum(y[, NCOL(y)] > 0, na.rm = TRUE)
  }

  class(fit) <- "mfpi"
  fit
}

#' @describeIn mfpi Formula method using a model formula and data frame.
#'
#' @method mfpi formula
#' @export
mfpi.formula <- function(formula,
                         data,
                         group_var         = NULL,
                         interaction_vars  = NULL,
                         interaction_forms = NULL,
                         flex              = c("flex3", "flex1", "flex2", "flex4"),
                         p_interact        = 0.05,
                         min_improvement   = NULL,
                         include_group_var = FALSE,
                         show_models       = FALSE,
                         weights           = NULL,
                         cycles            = 10,
                         scale             = NULL,
                         shift             = NULL,
                         df                = 4,
                         center            = TRUE,
                         subset            = NULL,
                         family            = "gaussian",
                         criterion         = c("pvalue", "aic", "bic"),
                         select            = 0.05,
                         alpha             = 0.05,
                         keep              = NULL,
                         force_max_fp_vars = NULL,
                         xorder            = c("ascending", "descending", "original"),
                         powers            = NULL,
                         ties              = c("efron", "breslow"),
                         nocenter          = c(-1, 0, 1),
                         zero_vars         = NULL,
                         catzero_vars      = NULL,
                         spike_vars        = NULL,
                         min_saz_prop      = 0.10,
                         ftest             = FALSE,
                         control           = NULL,
                         winsorize         = FALSE,
                         winsorize_probs   = c(0.01, 0.99),
                         center_type       = c("grand", "group"),
                         p_adjust_method   = "none",
                         verbose           = TRUE,
                         digits            = 3,
                         fitter            = c("base", "fastglm"),
                         ...) {

  # mfpi.formula() translates a formula + data.frame specification into the
  # matrix/vector inputs expected by mfpi.default(): it builds the model
  # frame, expands fp()/fp2() term attributes into per-variable option lists,
  # then calls mfpi.default() to perform validation, preprocessing, and fitting.
  call <- match.call()

  # Capture these promises before validation can force them. They are resolved
  # exactly once in the formula data mask below, matching the standard lookup
  # behavior of model-frame-based formula methods.
  weights_expr <- substitute(weights)
  subset_expr <- substitute(subset)

  # Step 1: Resolve multiple-choice arguments, family, and reject unused `...`
  dots <- match.call(expand.dots = FALSE)$...
  if ("offset" %in% names(dots)) {
    stop(
      "`offset` is not supported as an argument to `mfpi.formula()`. ",
      "Include `offset(...)` in the formula instead.",
      call. = FALSE
    )
  }
  if ("strata" %in% names(dots)) {
    stop(
      "`strata` is not supported as an argument to `mfpi.formula()`. ",
      "Include `strata(...)` in the formula instead.",
      call. = FALSE
    )
  }
  mfpi_check_unused_dots(
    dots = dots,
    context = "`mfpi.formula()`"
  )

  xorder      <- match.arg(xorder)
  criterion   <- match.arg(criterion)
  flex        <- match.arg(flex)
  ties        <- resolve_mfp_ties(ties)
  center_type <- match.arg(center_type)
  fitter      <- match.arg(fitter)

  # Family
  family_info <- normalize_family_argument(
    family,
    family_arg = deparse(substitute(family))
  )
  family_string <- family_info$family_string
  validate_mfpi_family(family_string)

  # Unused or misspelled arguments are already rejected above by
  # mfpi_check_unused_dots(), which reports the offending name(s) directly.
  # `offset` and `strata` receive specific formula-syntax guidance above.

  # Step 2: Validate the formula-specific inputs (data, formula) --------------
  # ---------------------------------------------------------------------------
  # Input validation (formula-specific constraints)
  # ---------------------------------------------------------------------------
  if (missing(data)) {
    stop("! data argument is missing.\n",
         "i An input data.frame is required for the formula interface.",
         call. = FALSE)
  }

  if (is.null(colnames(data))) {
    stop("! data must have column names.", call. = FALSE)
  }

  if (missing(formula)) {
    stop("! formula is missing.", call. = FALSE)
  }

  if (!inherits(formula, "formula")) {
    stop("! method is only for formula objects.", call. = FALSE)
  }

  formula_user <- formula
  formula_internal <- normalize_formula_special_namespaces(formula_user)

  # Step 3: Validate that df/alpha/select/scale/shift are scalar (global)
  # defaults only - per-variable settings belong in fp()/fp2() terms ---------
  # The formula interface only supports scalar defaults; per-variable settings
  # must be supplied via fp()/fp2() terms in the formula.
  if (length(df) != 1L) {
    stop(
      "! `df` must be a single numeric value.\n",
      "i Use `fp()` or `fp2()` in the formula to set per-variable `df` values.",
      call. = FALSE
    )
  }

  if (length(alpha) != 1L) {
    stop(
      "! `alpha` must be a single numeric value.\n",
      "i Use `fp()` or `fp2()` in the formula to set per-variable `alpha` values.",
      call. = FALSE
    )
  }

  if (length(select) != 1L) {
    stop(
      "! `select` must be a single numeric value.\n",
      "i Use `fp()` or `fp2()` in the formula to set per-variable `select` values.",
      call. = FALSE
    )
  }

  if (!is.null(scale) && length(scale) != 1L) {
    stop(
      "! `scale` must be a single numeric value, `NA`, or `NULL`.\n",
      "i Use `fp()` or `fp2()` in the formula to set per-variable `scale` values.",
      call. = FALSE
    )
  }

  if (!is.null(shift) && length(shift) != 1L) {
    stop(
      "! `shift` must be a single numeric value, `NA`, or `NULL`.\n",
      "i Use `fp()` or `fp2()` in the formula to set per-variable `shift` values.",
      call. = FALSE
    )
  }

  validate_numeric_vector(
    arg = df,
    name = "df",
    nvars = 1L,
    allow_null = FALSE,
    allow_na = FALSE,
    strictly_positive = TRUE
  )

  if (df != 1L && df %% 2L != 0L) {
    stop(
      "! `df` must be 1 or an even positive integer.\n",
      "i Use `fp()` or `fp2()` in the formula to set per-variable `df` values.",
      call. = FALSE
    )
  }

  validate_probability_vector(
    arg = alpha,
    name = "alpha",
    nvars = 1L
  )

  validate_probability_vector(
    arg = select,
    name = "select",
    nvars = 1L
  )

  validate_logical_vector(
    arg = center,
    name = "center",
    allowed_lengths = 1L
  )

  validate_numeric_vector(
    arg = scale,
    name = "scale",
    nvars = 1L,
    allow_null = TRUE,
    allow_na = TRUE,
    strictly_positive = TRUE
  )

  validate_numeric_vector(
    arg = shift,
    name = "shift",
    nvars = 1L,
    allow_null = TRUE,
    allow_na = TRUE,
    strictly_positive = FALSE
  )

  if (!is.null(powers) && !is.list(powers)) {
    stop("! powers must be a named list or NULL.", call. = FALSE)
  }


  # Discrete interaction forms are resolved after model.matrix() has
  # created their complete treatment-contrast blocks. At that point the
  # default method can enforce the linear-only rule and compute the correct
  # block width for the interaction degrees of freedom.

  # Evaluate formula observation arguments once, then validate them.
  n_data <- nrow(data)

  # Match standard formula-method lookup rules: data columns are searched first,
  # followed by the normalized formula environment (which inherits the user's
  # formula environment). Rebinding each argument to its resolved value prevents
  # later branches from re-evaluating expressions with side effects.
  weights <- eval(
    weights_expr,
    envir = data,
    enclos = environment(formula_internal)
  )
  subset <- eval(
    subset_expr,
    envir = data,
    enclos = environment(formula_internal)
  )
  strata <- NULL
  if (!is.null(subset)) {
    if (is.logical(subset)) {
      if (length(subset) != n_data || anyNA(subset)) {
        stop("! Logical `subset` must have one TRUE/FALSE value per row of `data` and contain no NA.", call. = FALSE)
      }
    } else if (is.numeric(subset)) {
      if (anyNA(subset) || any(!is.finite(subset)) ||
          any(subset != as.integer(subset)) || any(subset < 1L) ||
          any(subset > n_data)) {
        stop("! Numeric `subset` must contain valid positive integer row indices within `data`.", call. = FALSE)
      }
      if (anyDuplicated(subset)) {
        stop("! `subset` must not contain duplicated row indices.", call. = FALSE)
      }
    } else {
      stop("! `subset` must be either a logical vector or a numeric/integer vector of row indices.", call. = FALSE)
    }
  }

  # Apply the same strictly-positive weight contract as the matrix interface.
  # The complete user-supplied vector is validated even when `subset` is used,
  # keeping the public weight rule simple and family-independent.
  validate_model_weights(
    weights = weights,
    nobs = n_data
  )

  # Keep the ordinary formula path independent of subset-specific helpers.
  # A genuine subset is resolved once; otherwise all data rows are retained in
  # their original order.
  if (is.null(subset)) {
    fit_rows <- seq_len(n_data)
  } else {
    fit_rows <- formula_subset_rows(subset, n_data)
  }
  if (length(fit_rows) < 5L) {
    stop(
      paste0("! After subsetting, only ", length(fit_rows),
             " observations remain; at least 5 are required."),
      call. = FALSE
    )
  }

  # Step 5: Build full-data and fitted model frames ----------------------------
  # mf_full preserves the complete continuous values needed by automatic
  # preprocessing. With subset = NULL, mf reuses mf_full; otherwise it is rebuilt
  # from fit_rows and is authoritative for factor levels, contrast construction,
  # the response, fitting, and prediction metadata.
  # This mirrors standard formula-model behavior after the preprocessing source
  # has been fixed.
  mf_full <- stats::model.frame(
    formula_internal,
    data = data,
    drop.unused.levels = TRUE,
    na.action = stats::na.fail
  )
  # Reuse the complete frame when no subset was requested. For a genuine subset,
  # create a retained-row frame and drop unused factor levels before model.matrix()
  # rebuilds the fitted categorical coding.
  if (is.null(subset)) {
    mf <- mf_full
  } else {
    # Do not pass the local `fit_rows` symbol through model.frame(subset = ...):
    # model.frame() would re-evaluate it in the formula/data environment.
    mf <- subset_formula_model_frame(mf_full, fit_rows)
  }

  labels <- attr(terms(mf), "term.labels")
  if (length(labels) == 0L)
    stop("! No predictors found in formula. At least one predictor is required.",
         call. = FALSE)

  # Step 6: Extract strata() terms for survival models ------------------------
  # ---------------------------------------------------------------------------
  # Handle strata terms for survival models
  # ---------------------------------------------------------------------------
  specials      <- "strata"
  terms_formula <- stats::terms(
    formula_internal,
    specials = specials,
    data = data
  )
  terms_drop <- NULL
  formula_strata_terms <- NULL
  formula_strata_xlevels <- NULL
  formula_offset_terms <- NULL
  formula_offset_xlevels <- NULL

  if (!is.null(attr(terms_formula, "specials")$strata)) {
    if (mfp2_family_is_survival(family_string)) {
      stemp <- survival::untangle.specials(
        terms_formula,
        special = "strata",
        order = 1
      )

      strata_formula <- stats::reformulate(stemp$vars)
      environment(strata_formula) <- environment(formula_internal)
      formula_strata_terms <- stats::terms(
        strata_formula,
        data = data
      )
      formula_strata_xlevels <- .getXlevels(formula_strata_terms, mf)

      strata <- if (length(stemp$vars) == 1L) {
        mf[[stemp$vars]]
      } else {
        do.call(
          survival::strata,
          c(as.list(mf[, stemp$vars, drop = FALSE]), list(shortlabel = TRUE))
        )
      }

      terms_drop <- stemp$terms
    } else {
      stop("! strata are only allowed for survival models.\n",
           "i Please remove strata terms from the formula.", call. = FALSE)
    }
  }

  terms_model <- if (!is.null(terms_drop)) terms_formula[-terms_drop]
  else terms_formula

  validate_formula_factor_levels(terms_model, mf)

  # Step 7: Identify categorical formula terms before model.matrix() expands
  # them. Both unordered and ordered factors are retained using their configured
  # contrasts and are grouped for joint adjustment-model selection.
  factor_term_labels <- identify_formula_factor_terms(terms_model, mf)
  conceptual_term_names <- formula_conceptual_term_names(terms_model, mf)
  # Retain the original formula-label -> conceptual-term relationship for
  # prediction. This map is updated below when fp()/fp2() expressions are
  # renamed to their source-variable names.
  formula_prediction_term_names <- conceptual_term_names
  factor_terms <- unname(conceptual_term_names[factor_term_labels])

  # Step 8: Extract an offset() term from the formula, if present ------------
  # ---------------------------------------------------------------------------
  # Handle offset
  # ---------------------------------------------------------------------------
  offset <- NULL
  term_offset <- attr(terms_formula, "offset")
  if (!is.null(term_offset) && length(term_offset) > 1L)
    stop("! Only one offset term is allowed in the formula.", call. = FALSE)

  if (!is.null(term_offset)) {
    # Always extract offset from model frame when formula contains offset().
    # Also keep the offset expression itself so predict.mfpi() can rebuild
    # expressions such as offset(log(exposure)) from raw newdata.
    offset <- stats::model.offset(mf)
    if (!identical(family_string, "multinomial")) offset <- as.vector(offset)

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

  # Step 9: Extract y, build the model matrix, and drop the intercept --------
  # ---------------------------------------------------------------------------
  # Extract y and x
  # ---------------------------------------------------------------------------
  y <- stats::model.extract(mf, "response")

  if (!mfp2_family_is_survival(family_string) && survival::is.Surv(y))
    stop(sprintf(
      paste0(
        "! Response is a Surv object but family = '%s'. Use family = 'cox' ",
        "or survreg_family() as appropriate."
      ),
      family_string
    ), call. = FALSE)

  # Build the two designs independently. x is the actual fitted design and may
  # have fewer or differently coded factor columns after unused levels are
  # dropped. x_full is used only to recover full-data preprocessing values for
  # stable, same-named continuous columns.
  x_full <- stats::model.matrix(terms_model, mf_full)
  x <- stats::model.matrix(terms_model, mf)

  # Developer note: save formula reconstruction metadata before modifying the
  # model matrix. The subsequent intercept removal, factor-group handling, and
  # fp()/fp2() renaming can otherwise drop or obscure these attributes. These
  # fields allow predict.mfpi() to rebuild the same design from
  # ordinary formula-style newdata.
  x_contrasts <- attr(x, "contrasts")
  x_xlevels   <- .getXlevels(terms_model, mf)

  # Save the assign attribute before subsetting (matrix subsetting drops it)
  x_assign <- attr(x, "assign")

  # Remove intercept column (entry 0 in model.matrix "assign" attribute)
  if (0L %in% x_assign) {
    intercept_col <- which(x_assign == 0L)
    x        <- x[, -intercept_col, drop = FALSE]
    x_assign <- x_assign[-intercept_col]
  }
  full_assign <- attr(x_full, "assign")
  if (0L %in% full_assign) {
    x_full <- x_full[, full_assign != 0L, drop = FALSE]
  }

  # Map conceptual terms to the exact raw design columns they generated.
  # Simple factor wrappers use the source variable name, while model.matrix()
  # column names such as factor(stage)II remain unchanged.
  term_labels <- attr(terms_model, "term.labels")
  term_to_columns <- build_formula_term_to_columns(
    x_columns = colnames(x),
    assign = x_assign,
    term_labels = term_labels,
    conceptual_names = conceptual_term_names
  )

  # Step 10: Convert a categorical group_var into a numeric column ----------
  # ---------------------------------------------------------------------------
  # Handle categorical group_var from the original data.
  # model.matrix() expands factor/character predictors into dummy columns, so
  # the original column name may disappear. mfpi.default() expects group_var as
  # a single numeric column. We therefore remove any dummy columns for group_var
  # and add back numeric codes, while preserving the original labels as a matrix
  # attribute for downstream reporting.
  # ---------------------------------------------------------------------------
  group_levels_original <- NULL

  if (!is.null(group_var)) {
    gv_raw <- mf[[group_var]]

    if (is.factor(gv_raw)) {
      gv_factor <- droplevels(gv_raw)
      group_levels_original <- levels(gv_factor)
      gv_numeric <- as.integer(gv_factor)
    } else if (is.character(gv_raw)) {
      group_levels_original <- unique(gv_raw)
      gv_factor <- factor(gv_raw, levels = group_levels_original)
      gv_numeric <- as.integer(gv_factor)
    } else if (is.logical(gv_raw)) {
      group_levels_original <- as.character(sort(unique(gv_raw)))
      gv_numeric <- as.numeric(gv_raw)
    } else {
      gv_numeric <- as.numeric(gv_raw)
      group_levels_original <- as.character(sort(unique(gv_numeric)))
    }

    if (!group_var %in% colnames(x)) {
      # group_var was expanded by model.matrix. Identify and remove its dummies
      # using the saved assign attribute.
      gv_term_idx <- which(term_labels == group_var)

      if (length(gv_term_idx) == 1L) {
        gv_cols <- which(x_assign == gv_term_idx)

        if (length(gv_cols) > 0L) {
          x        <- x[, -gv_cols, drop = FALSE]
          x_assign <- x_assign[-gv_cols]
        }
      }

      x <- cbind(x, gv_numeric)
      colnames(x)[ncol(x)] <- group_var
    }

    # group_var is handled by MFPI's dedicated categorical-group machinery, not
    # as a grouped adjustment term. Keep it as one singleton conceptual term.
    term_to_columns[[group_var]] <- group_var

    attr(x, "mfpi_group_levels_original") <- group_levels_original
  }

  nx      <- ncol(x)
  names_x <- colnames(x)

  # Step 11: Detect fp()/fp2() terms and rename their columns to plain
  # variable names, then initialize per-variable option lists from the
  # global scalar defaults ----------------------------------------------------
  # ---------------------------------------------------------------------------
  # Detect fp()/fp2() terms in the model frame and extract their attributes
  # ---------------------------------------------------------------------------

  fp_pos <- which(is_fp_term(colnames(mf)))
  # ---------------------------------------------------------------------------
  # Resolve final variable names: fp()/fp2() terms rename to raw variable names in x.
  # We do the rename now so all parameter lists use the correct final names.
  # ---------------------------------------------------------------------------
  if (length(fp_pos) > 0L) {
    fp_data_pre <- mf[, fp_pos, drop = FALSE]
    fp_vars_pre <- unname(vapply(fp_data_pre, function(v) attr(v, "name"), character(1L)))
    fp_x_pos <- which(is_fp_term(names_x))

    if (length(fp_x_pos) != length(fp_vars_pre)) {
      stop(
        "! Internal formula parsing error: the number of fp()/fp2() terms in the model matrix does not match the model frame.",
        call. = FALSE
      )
    }

    old_fp_names <- names_x[fp_x_pos]
    names_x <- replace(names_x, fp_x_pos, fp_vars_pre)
    colnames(x) <- names_x
    full_fp_pos <- which(is_fp_term(colnames(x_full)))
    if (length(full_fp_pos) != length(fp_vars_pre)) {
      stop(
        "! Internal formula parsing error: full-data and fitted fp()/fp2() columns do not align.",
        call. = FALSE
      )
    }
    colnames(x_full)[full_fp_pos] <- fp_vars_pre

    # Keep the conceptual lookup aligned with the renamed fp()/fp2() columns.
    # formula_conceptual_term_names() already maps each original fp() formula
    # label to its source-variable name for prediction metadata.
    for (i in seq_along(old_fp_names)) {
      fp_column <- old_fp_names[[i]]
      fp_var <- fp_vars_pre[[i]]
      if (fp_var %in% names(term_to_columns)) {
        term_to_columns[[fp_var]] <- ifelse(
          term_to_columns[[fp_var]] == fp_column,
          fp_var,
          term_to_columns[[fp_var]]
        )
      }
    }
  }

  # Initialise per-variable parameter lists from global defaults.
  df_list <- setNames(rep(list(df), nx), names_x)
  shift_list <- if (is.null(shift)) {
    setNames(rep(list(NA_real_), nx), names_x)
  } else {
    setNames(rep(list(shift), nx), names_x)
  }

  scale_list <- if (is.null(scale)) {
    setNames(rep(list(NA_real_), nx), names_x)
  } else {
    setNames(rep(list(scale), nx), names_x)
  }

  center_list <- setNames(rep(list(center), nx), names_x)
  alpha_list  <- setNames(rep(list(alpha),  nx), names_x)
  select_list <- setNames(rep(list(select), nx), names_x)

  # acdx, zero, catzero, spike: initialise to FALSE for all variables.
  # Overridden below only for variables wrapped in fp(). Argument-level
  # zero_vars / catzero_vars / spike_vars are merged separately after.
  acdx_list    <- setNames(rep(list(FALSE), nx), names_x)
  zero_list    <- setNames(rep(list(FALSE), nx), names_x)
  catzero_list <- setNames(rep(list(FALSE), nx), names_x)
  spike_list   <- setNames(rep(list(FALSE), nx), names_x)
  force_max_fp_list <- setNames(rep(list(FALSE), nx), names_x)

  # Default FP candidate powers (Royston and Altman, 1994)
  powx        <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
  power_list  <- setNames(replicate(nx, powx, simplify = FALSE), names_x)

  # Override with user-supplied powers argument (before fp() overrides)
  if (!is.null(powers)) {
    if (is.null(names(powers)) ||
        length(powers) != sum(nchar(names(powers)) > 0L, na.rm = TRUE))
      stop("! All elements of powers must have names.", call. = FALSE)
    bad <- setdiff(names(powers), names_x)
    if (length(bad) > 0L)
      stop(paste0("! powers names not found in x: ",
                  paste(bad, collapse = ", "), "."), call. = FALSE)
    fp_powers  <- lapply(powers, sort)
    power_list <- modifyList(power_list, Filter(Negate(is.null), fp_powers))
  }

  # Step 12: Apply fp()/fp2() attribute overrides onto the default lists -----
  # ---------------------------------------------------------------------------
  # Apply fp() term overrides
  # ---------------------------------------------------------------------------
  if (length(fp_pos) > 0L) {
    fp_data <- mf[, fp_pos, drop = FALSE]
    # fp_vars already computed as fp_vars_pre in the rename block above;
    # reuse here under the canonical name fp_vars.
    fp_vars <- fp_vars_pre

    # Guard against duplicated fp()/fp2() usage
    dups <- fp_vars[duplicated(fp_vars)]
    if (length(dups) > 0L)
      stop(paste0("! Variables used more than once in fp()/fp2(): ",
                  paste(dups, collapse = ", "), "."), call. = FALSE)

    # Guard against variables appearing both inside fp()/fp2() and as plain terms.
    # colnames(mf) contains "fp(age)"/"fp2(age)" for FP terms and the response variable.
    # Exclude the response (first column when LHS is present) and FP columns;
    # check remaining predictor names against fp_vars.
    response_pos <- attr(stats::terms(formula), "response")
    response_col <- if (response_pos > 0L) names(mf)[response_pos] else character(0L)

    pred_cols  <- setdiff(colnames(mf), response_col)
    plain_cols <- pred_cols[!is_fp_term(pred_cols)]
    also_plain   <- intersect(plain_cols, fp_vars)

    if (length(also_plain) > 0L)
      stop(paste0("! Variables in fp()/fp2() must not appear elsewhere in the formula: ",
                  paste(also_plain, collapse = ", "), "."), call. = FALSE)

    # Note: fp()/fp2() column names were already renamed above (before list init)
    # so names_x and colnames(x) already use the real variable names.

    # Helper: apply fp() attribute overrides to a named default list.
    # attr_name: the fp() attribute to extract.
    # drop_null: if TRUE, entries where fp() returned NULL are dropped so
    #   the estimated default (e.g. auto shift/scale) is not overwritten.
    apply_fp_override <- function(base_list, attr_name, drop_null = FALSE) {
      overrides <- setNames(lapply(fp_data, attr, attr_name), fp_vars)
      if (drop_null)
        overrides <- Filter(Negate(is.null), overrides)
      modifyList(base_list, overrides)
    }

    df_list      <- apply_fp_override(df_list,      "df")
    scale_list   <- apply_fp_override(scale_list,   "scale",   drop_null = TRUE)
    shift_list   <- apply_fp_override(shift_list,   "shift",   drop_null = TRUE)
    center_list  <- apply_fp_override(center_list,  "center")
    alpha_list   <- apply_fp_override(alpha_list,   "alpha")
    select_list  <- apply_fp_override(select_list,  "select")
    acdx_list    <- apply_fp_override(acdx_list,    "acd")
    zero_list    <- apply_fp_override(zero_list,    "zero")
    catzero_list <- apply_fp_override(catzero_list, "catzero")
    spike_list   <- apply_fp_override(spike_list,   "spike")
    force_max_fp_list <- apply_fp_override(force_max_fp_list, "force_max_fp")

    # fp()/fp2() powers take precedence over the powers argument
    fp_pow_override <- Filter(Negate(is.null),
                              setNames(lapply(fp_data, attr, "powers"), fp_vars))
    conflict <- if (!is.null(powers)) {
      intersect(names(fp_pow_override), names(powers))
    } else {
      character(0L)
    }

    if (length(conflict) > 0L)
      warning(paste0("i Powers specified in both fp()/fp2() and powers argument; ",
                     "fp()/fp2() values take precedence for: ",
                     paste(conflict, collapse = ", "), "."), call. = FALSE)
    power_list <- modifyList(power_list, fp_pow_override)
  }

  # Factor-generated contrast columns are one fixed linear adjustment term.
  # Preserve their fitted contrast coding but disable FP preprocessing.
  factor_columns <- unique(unlist(
    term_to_columns[intersect(factor_terms, names(term_to_columns))],
    use.names = FALSE
  ))
  factor_columns <- intersect(factor_columns, names_x)

  if (length(factor_columns) > 0L) {
    df_list[factor_columns] <- rep(list(1L), length(factor_columns))
    shift_list[factor_columns] <- rep(list(0), length(factor_columns))
    scale_list[factor_columns] <- rep(list(1), length(factor_columns))
    acdx_list[factor_columns] <- rep(list(FALSE), length(factor_columns))
    zero_list[factor_columns] <- rep(list(FALSE), length(factor_columns))
    catzero_list[factor_columns] <- rep(list(FALSE), length(factor_columns))
    spike_list[factor_columns] <- rep(list(FALSE), length(factor_columns))
    force_max_fp_list[factor_columns] <- rep(list(FALSE), length(factor_columns))
    power_list[factor_columns] <- rep(list(1), length(factor_columns))
  }

  # Step 13: Convert list-form per-variable parameters to plain vectors,
  # and merge fp()-derived zero/catzero/spike/acd variable names with the
  # corresponding *_vars arguments ---------------------------------------------
  # ---------------------------------------------------------------------------
  # Convert list-form parameters to vectors / character vectors for mfpi.default
  # ---------------------------------------------------------------------------

  # Per-variable numeric vectors. vapply extracts the first (and only) element
  # of each list slot, producing a clean named numeric vector regardless of
  # whether the slot holds an auto-estimated scalar or a user override.
  scale_vec  <- vapply(scale_list, function(v) as.numeric(v[[1L]]), numeric(1L))
  shift_vec  <- vapply(shift_list, function(v) as.numeric(v[[1L]]), numeric(1L))

  # df: unlist gives a named numeric vector. mfpi.default() will apply
  # assign_df() cardinality overrides on top of these values.
  df_vec     <- unlist(df_list)

  # center: unlist gives a named logical vector.
  center_vec <- unlist(center_list)

  # select and alpha: unlist gives named numeric vectors.
  select_vec <- unlist(select_list)
  alpha_vec  <- unlist(alpha_list)

  # acd_vars: character vector of variable names with acd = TRUE, or NULL
  acdx_vec  <- unlist(acdx_list)
  acd_vars  <- if (any(acdx_vec)) names(acdx_vec[acdx_vec]) else NULL

  # zero_vars from fp(): merge with argument zero_vars
  zero_fp         <- unlist(zero_list)
  zero_from_fp    <- if (any(zero_fp)) names(zero_fp[zero_fp]) else character(0L)
  zero_vars_final <- union(zero_from_fp, zero_vars)
  zero_vars_final <- if (length(zero_vars_final) == 0L) NULL else zero_vars_final

  # catzero_vars from fp(): merge with argument catzero_vars
  catzero_fp         <- unlist(catzero_list)
  catzero_from_fp    <- if (any(catzero_fp)) names(catzero_fp[catzero_fp]) else character(0L)
  catzero_vars_final <- union(catzero_from_fp, catzero_vars)
  catzero_vars_final <- if (length(catzero_vars_final) == 0L) NULL else catzero_vars_final

  # spike_vars from fp(): merge with argument spike_vars
  spike_fp         <- unlist(spike_list)
  spike_from_fp    <- if (any(spike_fp)) names(spike_fp[spike_fp]) else character(0L)
  spike_vars_final <- union(spike_from_fp, spike_vars)
  spike_vars_final <- if (length(spike_vars_final) == 0L) NULL else spike_vars_final

  # force_max_fp_vars from fp(): merge with argument force_max_fp_vars
  fmfp_fp         <- unlist(force_max_fp_list)
  fmfp_from_fp    <- if (any(fmfp_fp)) names(fmfp_fp[fmfp_fp]) else character(0L)
  force_max_fp_vars_final <- union(fmfp_from_fp, force_max_fp_vars)
  force_max_fp_vars_final <- if (length(force_max_fp_vars_final) == 0L) NULL else force_max_fp_vars_final

  # Expand formula-term keep entries using exact term-to-column mapping. This
  # keeps all contrast columns of a categorical term together.
  if (!is.null(keep)) {
    valid_keep <- unique(c(
      colnames(x),
      names(term_to_columns),
      names(formula_prediction_term_names)
    ))
    bad_keep <- setdiff(keep, valid_keep)
    if (length(bad_keep) > 0L) {
      stop(
        paste0("! Unknown variable(s) in keep: ", paste(bad_keep, collapse = ", "), "."),
        call. = FALSE
      )
    }
    keep <- unique(unlist(lapply(keep, function(value) {
      if (value %in% colnames(x)) {
        return(value)
      }
      conceptual_name <- if (value %in% names(formula_prediction_term_names)) {
        unname(formula_prediction_term_names[[value]])
      } else {
        value
      }
      term_to_columns[[conceptual_name]]
    }), use.names = FALSE))
  }

  # Supply all categorical mappings to the grouped MFP adjustment core,
  # including binary factors that generate one non-identity dummy column.
  grouped_formula_terms <- term_to_columns[
    intersect(factor_terms, names(term_to_columns))
  ]
  if (length(grouped_formula_terms) > 0L) {
    grouped_formula_terms <- grouped_formula_terms[
      mapped_term_flags(grouped_formula_terms)
    ]
  }
  if (length(grouped_formula_terms) == 0L) grouped_formula_terms <- NULL

  # Step 14: Delegate using the subset-specific fitting design ----------------
  # The private attribute transports the full-data preprocessing source and
  # tells mfpi.default() that the fitting matrix has already been subset. The
  # default method performs zero-specific support validation before the general
  # singleton-column variation check.
  x <- attach_formula_preprocess_matrix(x, x_full)
  if (!is.null(subset)) {
    if (!is.null(weights)) weights <- weights[fit_rows]
  }

  # ---------------------------------------------------------------------------
  # Delegate to mfpi.default
  # ---------------------------------------------------------------------------
  fit <-  mfpi.default(
    x                 = x,
    y                 = y,
    group_var         = group_var,
    interaction_vars         = interaction_vars,
    interaction_forms    = interaction_forms,
    flex              = flex,
    p_interact        = p_interact,
    min_improvement   = min_improvement,
    include_group_var = include_group_var,
    show_models       = show_models,
    weights           = weights,
    offset            = offset,
    cycles            = cycles,
    scale             = scale_vec,
    shift             = shift_vec,
    df                = df_vec,
    center            = center_vec,
    subset            = NULL,
    family            = family,
    fitter            = fitter,
    criterion         = criterion,
    select            = select_vec,
    alpha             = alpha_vec,
    keep              = keep,
    force_max_fp_vars = force_max_fp_vars_final,
    xorder            = xorder,
    powers            = power_list,
    ties              = ties,
    strata            = strata,
    nocenter          = nocenter,
    acd_vars          = acd_vars,
    zero_vars         = zero_vars_final,
    catzero_vars      = catzero_vars_final,
    spike_vars        = spike_vars_final,
    min_saz_prop = min_saz_prop,
    ftest             = ftest,
    control           = control,
    winsorize         = winsorize,
    winsorize_probs   = winsorize_probs,
    center_type       = center_type,
    p_adjust_method   = p_adjust_method,
    verbose           = verbose,
    digits            = digits,
    term_groups       = grouped_formula_terms
  )
  fit$formula_interface <- TRUE
  fit$formula <- formula_user
  fit$formula_terms <- stats::delete.response(terms_model)
  fit$formula_contrasts <- x_contrasts
  fit$formula_xlevels <- x_xlevels
  # predict.mfpi() needs to know the final fit-time column names after intercept
  # removal, group-variable replacement, and fp() renaming.
  fit$formula_design_columns <- colnames(x)
  fit$formula_term_to_columns <- term_to_columns
  fit$formula_prediction_term_names <- formula_prediction_term_names
  fit$formula_factor_terms <- factor_terms
  fit$formula_offset_terms <- formula_offset_terms
  fit$formula_offset_xlevels <- formula_offset_xlevels
  fit$formula_strata_terms <- formula_strata_terms
  fit$formula_strata_xlevels <- formula_strata_xlevels

  # Replace design-pattern labels for formula-supplied discrete interactions
  # with the original binary/factor levels. The corresponding contrast rows are
  # retained so prediction and plotting can evaluate one point per level
  # without recomputing or guessing the fit-time coding.
  for (term in interaction_vars) {
    spec <- fit$interaction_specs[[term]]
    if (is.null(spec) || !isTRUE(spec$discrete) || !term %in% names(mf)) next

    raw_term <- mf[[term]]
    level_labels <- if (is.factor(raw_term)) {
      levels(droplevels(raw_term))
    } else if (is.character(raw_term)) {
      unique(raw_term)
    } else if (is.logical(raw_term)) {
      as.character(sort(unique(raw_term)))
    } else {
      as.character(sort(unique(raw_term)))
    }

    raw_labels <- as.character(raw_term)
    first_rows <- match(level_labels, raw_labels)
    if (anyNA(first_rows)) next

    spec$level_labels <- level_labels
    spec$level_design <- x[first_rows, spec$columns, drop = FALSE]
    spec$formula_discrete <- TRUE
    fit$interaction_specs[[term]] <- spec
  }
  fit$call <- call

  fit
}

#' Format unused dot arguments for error messages
#'
#' @param dots Pairlist or list of arguments captured from \code{...}.
#'
#' @return Character vector of display labels.
#'
#' @keywords internal
#' @noRd
mfpi_dot_labels <- function(dots) {
  if (length(dots) == 0L) return(character(0L))

  dot_names <- names(dots)
  if (is.null(dot_names)) {
    dot_names <- rep("", length(dots))
  }

  named <- !is.na(dot_names) & nzchar(dot_names)
  labels <- character(length(dots))

  # Named arguments are labeled by name (e.g. `foo`); unnamed ones are
  # labeled by a truncated deparse of their value, since there is no name to
  # show for them.
  labels[named] <- paste0("`", dot_names[named], "`")

  if (any(!named)) {
    labels[!named] <- vapply(
      dots[!named],
      function(z) paste(deparse(z, width.cutoff = 60L), collapse = " "),
      character(1L)
    )
  }

  labels
}


#' Reject unused dot arguments
#'
#' @param dots Pairlist or list of arguments captured from \code{...}.
#' @param context Character scalar naming the interface being validated.
#'
#' @return Invisibly returns \code{TRUE} when no arguments are supplied.
#'
#' @keywords internal
#' @noRd
mfpi_check_unused_dots <- function(dots, context) {
  if (length(dots) == 0L) return(invisible(TRUE))

  stop(
    paste0(
      "Unused argument(s) supplied to ", context, ": ",
      paste(mfpi_dot_labels(dots), collapse = ", "),
      "."
    ),
    call. = FALSE
  )
}

#' Prepare Predictors for the Default MFPI Method
#'
#' Converts a numeric matrix or data frame to a numeric predictor matrix while
#' preserving the grouping variable's original labels.
#'
#' In a data frame, the grouping variable may be factor, character, logical,
#' or numeric. Other predictors must be numeric or logical; logical predictors
#' are converted to numeric values.
#'
#' Group order follows factor levels, first appearance for character values,
#' or increasing order for logical and numeric values.
#'
#' @param x A numeric predictor matrix or data frame.
#' @param group_var A character string naming the grouping variable.
#'   Numeric column indices are not supported.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{x}}{The numeric predictor matrix.}
#'   \item{\code{group_var}}{The validated grouping-variable name.}
#'   \item{\code{group_input_levels}}{Numeric group values in the returned
#'     matrix, in group order.}
#'   \item{\code{group_levels_original}}{Original group labels corresponding
#'     to \code{group_input_levels}.}
#' }
#'
#' @keywords internal
#' @noRd
prepare_mfpi_default_x <- function(x, group_var) {
  # The public API requires group_var to be a column name. Do not support
  # numeric column indices here, because allowing both names and positions makes
  # data-frame input ambiguous and harder to document.
  if (!is.character(group_var) ||
      length(group_var) != 1L ||
      is.na(group_var) ||
      !nzchar(group_var)) {
    stop("`group_var` must be a single character variable name.", call. = FALSE)
  }

  # Step 1: Data-frame input ---------------------------------------------------
  # ---------------------------------------------------------------------------
  # Data-frame input
  # ---------------------------------------------------------------------------
  # Data frames can preserve factor and character labels. We allow categorical
  # input only for group_var, convert it to numeric codes for fitting, and keep
  # the original labels separately for display.
  if (is.data.frame(x)) {
    # Validate names before resolving group_var. This prevents duplicate or
    # malformed names from making x[[group_var]] ambiguous. The check is
    # performed once per public call and is outside all MFPI fitting loops.
    validate_predictor_names(names(x), object = "`x`")

    # group_var must be present as a named column.
    if (!group_var %in% names(x)) {
      stop("`group_var` was not found in `x`.", call. = FALSE)
    }

    # Extract the user-supplied grouping variable before any conversion.
    group_raw <- x[[group_var]]

    # Missing group labels cannot be mapped reliably to internal group levels.
    if (anyNA(group_raw)) {
      stop("`group_var` must not contain missing values.", call. = FALSE)
    }

    # Convert the grouping variable to numeric input codes and preserve
    # user-facing labels. The numeric input codes are not the final internal
    # levels; preprocess_data() later remaps them to 0, 1, ..., K - 1.
    if (is.factor(group_raw)) {
      # Factor input respects the user's factor-level order. The first retained
      # factor level is the reference label.
      group_factor <- droplevels(group_raw)
      group_levels_original <- levels(group_factor)
      group_numeric <- as.integer(group_factor)
      group_input_levels <- seq_along(group_levels_original)
    } else if (is.character(group_raw)) {
      # Character input has no explicit level order, so preserve first-seen
      # order. The first distinct value encountered is the reference label.
      group_levels_original <- unique(group_raw)
      group_factor <- factor(group_raw, levels = group_levels_original)
      group_numeric <- as.integer(group_factor)
      group_input_levels <- seq_along(group_levels_original)
    } else if (is.logical(group_raw)) {
      # Logical input is treated as a two-level categorical variable when both
      # values are observed. Sorting puts FALSE before TRUE.
      group_input_levels <- sort(unique(as.numeric(group_raw)))
      group_levels_original <- as.character(as.logical(group_input_levels))
      group_numeric <- as.numeric(group_raw)
    } else if (is.numeric(group_raw) || is.integer(group_raw)) {
      # Numeric/integer input is already usable for modelling. Sorting makes the
      # smallest observed value the reference level.
      group_input_levels <- sort(unique(as.numeric(group_raw)))
      group_levels_original <- as.character(group_input_levels)
      group_numeric <- as.numeric(group_raw)
    } else {
      stop(
        "`group_var` must be factor, character, logical, integer, or numeric.",
        call. = FALSE
      )
    }

    # Only group_var may be categorical. Other columns are not dummy-coded by
    # mfpi.default(); users should either supply numeric encodings explicitly or
    # use the formula method where appropriate.
    non_group <- setdiff(names(x), group_var)

    bad <- non_group[!vapply(x[non_group], function(z) {
      is.numeric(z) || is.integer(z) || is.logical(z)
    }, logical(1L))]

    if (length(bad) > 0L) {
      stop(
        "Only `group_var` may be categorical in `mfpi.default()`. ",
        "Non-group variables must be numeric, integer, or logical. ",
        "Problem variable(s): ",
        paste(bad, collapse = ", "),
        call. = FALSE
      )
    }

    # Work on a copy so the caller's data frame is not modified.
    x_work <- x

    # Replace the original grouping variable with numeric input codes.
    x_work[[group_var]] <- group_numeric

    # Convert logical non-group variables to numeric values so as.matrix()
    # produces a numeric matrix rather than a mixed-type matrix.
    for (nm in non_group) {
      if (is.logical(x_work[[nm]])) {
        x_work[[nm]] <- as.numeric(x_work[[nm]])
      }
    }

    # Convert the validated data frame to the numeric matrix expected by the
    # fitting internals.
    x_mat <- as.matrix(x_work)
    storage.mode(x_mat) <- "double"

    return(list(
      x = x_mat,
      group_var = group_var,
      group_input_levels = group_input_levels,
      group_levels_original = group_levels_original
    ))
  }

  # Step 2: Matrix input --------------------------------------------------------
  # ---------------------------------------------------------------------------
  # Matrix input
  # ---------------------------------------------------------------------------
  # Matrix input preserves the historical default-method contract: the predictor
  # matrix must already be numeric, and group_var must name one of its columns.
  if (!is.matrix(x)) {
    stop("`x` must be a matrix or data frame.", call. = FALSE)
  }

  # The matrix interface uses the same predictor-name contract as data-frame
  # input. Validate before group_var lookup so all downstream name matching is
  # unambiguous.
  validate_predictor_names(colnames(x), object = "`x`")

  if (!is.numeric(x)) {
    stop(
      "`x` must be a numeric matrix, or a data frame with only `group_var` categorical.",
      call. = FALSE
    )
  }

  if (!group_var %in% colnames(x)) {
    stop("`group_var` was not found in `x`.", call. = FALSE)
  }

  # For numeric matrix input, the observed group values are both the input
  # levels and, unless overridden by an internal formula-method attribute, the
  # user-facing labels.
  group_input_levels <- sort(unique(as.numeric(x[, group_var])))

  # Formula input may arrive here as a numeric matrix with this attribute set.
  # When present, it carries the original factor/character labels from the model
  # frame. Ordinary matrix input will not have this attribute.
  group_levels_original <- attr(
    x,
    "mfpi_group_levels_original",
    exact = TRUE
  )

  if (is.null(group_levels_original)) {
    group_levels_original <- as.character(group_input_levels)
  }

  list(
    x = x,
    group_var = group_var,
    group_input_levels = group_input_levels,
    group_levels_original = as.character(group_levels_original)
  )
}
