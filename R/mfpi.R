#' Model Interactions Between a Categorical Variable and Continuous Variables
#'
#' Use `mfpi()` to investigate whether the relationship between an outcome and
#' one or more continuous variables differs across the levels of a categorical
#' grouping variable.
#'
#' The function supports Gaussian, binomial, Poisson, negative-binomial, and Cox models. Nonlinear
#' relationships can be modelled with fractional polynomials.
#'
#' @details
#' `mfpi()` evaluates each variable in `cont_vars` separately. For every
#' continuous variable, it compares:
#'
#' - a model in which the continuous-variable function is the same across
#'   groups; and
#' - a model in which the function can differ across groups.
#'
#' The result is one term-specific interaction model for each variable in
#' `cont_vars`. The function does not fit all tested interactions together in
#' one final model.
#'
#' Before testing the interactions, `mfpi()` builds an MFP adjustment model from
#' the supplied predictors. Variables selected in this first step are used as
#' adjustment variables in the interaction analyses. The continuous variable
#' currently being tested is excluded from its own adjustment set because it is
#' already represented by its main effect and interaction.
#'
#' @section Basic use:
#' The main arguments are:
#'
#' - `group_var`: the categorical variable defining the groups;
#' - `cont_vars`: the continuous variables to test for interaction with
#'   `group_var`;
#' - `cont_var_forms`: the form used for each tested interaction;
#' - `flex`: how fractional-polynomial powers are chosen across groups;
#' - `criterion`: how the interaction and adjustment models are selected.
#'
#' For example:
#'
#' ```r
#' fit <- mfpi(
#'   outcome ~ treatment + fp(age) + fp(biomarker) + stage,
#'   data = dat,
#'   group_var = "treatment",
#'   cont_vars = c("age", "biomarker"),
#'   cont_var_forms = c(age = "linear", biomarker = "fp1")
#' )
#' ```
#'
#' This fits one interaction analysis for `age` and another for `biomarker`.
#'
#' @section Choosing the interaction form:
#' Use `cont_var_forms` to specify the form tested for each variable in
#' `cont_vars`:
#'
#' - `"linear"`: a linear group-by-variable interaction;
#' - `"fp1"`: a first-degree fractional-polynomial interaction;
#' - `"fp2"`: a second-degree fractional-polynomial interaction.
#'
#' The form is specified before fitting. `mfpi()` does not compare linear, FP1,
#' and FP2 forms and choose among them. A variable not named in
#' `cont_var_forms` uses `"fp1"`.
#'
#' In the formula interface, wrapping a variable in `fp()` or `fp2()` does not
#' set its MFPI interaction form. The wrappers provide variable-specific FP and
#' preprocessing settings. Some of those settings are also reused when the
#' interaction function is fitted, but the interaction degree itself is chosen
#' only through `cont_var_forms`.
#'
#' @section Choosing the flexibility level:
#' The `flex` argument controls how FP powers are chosen for variables specified
#' as `"fp1"` or `"fp2"`. It has no effect on variables specified as
#' `"linear"`.
#'
#' - `"flex1"`: FP powers are chosen from the model without interaction and
#'   then used for every group. This is the most restrictive option.
#' - `"flex2"`: FP powers are chosen using the interaction model, but the same
#'   powers are used for every group and for the model without interaction.
#' - `"flex3"`: powers for the model without interaction and the interaction
#'   model are chosen separately, but all groups share the same interaction
#'   powers. This is the default.
#' - `"flex4"`: each group can have different interaction powers. This is the
#'   most flexible option.
#'
#' More flexible choices use more information from the data and can be less
#' stable in small samples. For `flex2`, power selection is not fully reflected
#' in the interaction p-value. For `flex3` and `flex4`, the compared models can
#' be non-nested when different powers are chosen. Interpret p-values from these
#' options with care. See the MFPI vignette for methodological details.
#'
#' @section Selecting interactions:
#' The `criterion` argument is used both to build the MFP adjustment model and
#' to decide whether each term-specific interaction model is selected.
#'
#' - With `criterion = "pvalue"`, an interaction is selected when its raw or
#'   adjusted p-value is strictly below `p_interact`.
#' - With `criterion = "aic"`, an interaction is selected when
#'   `AIC_main_minus_int` is strictly greater than `min_improvement`.
#' - With `criterion = "bic"`, an interaction is selected when
#'   `BIC_main_minus_int` is strictly greater than `min_improvement`.
#'
#' The default `min_improvement` is 2 for AIC and BIC. `p_interact` is used only
#' for p-value selection.
#'
#' The AIC and BIC columns reported by `mfpi()` are interaction-comparison
#' quantities. Do not compare their individual `AIC_main`, `AIC_interaction`,
#' `BIC_main`, or `BIC_interaction` values with values returned by
#' [stats::AIC()] or [stats::BIC()]. Use `AIC_main_minus_int` or
#' `BIC_main_minus_int` for the MFPI comparison.
#'
#' @section Multiple interaction tests:
#' When several variables are listed in `cont_vars`, the default p-values are
#' not adjusted for multiple testing. Set `p_adjust_method` to a method accepted
#' by [stats::p.adjust()], such as `"holm"`, `"bonferroni"`, or `"BH"`,
#' when adjustment is appropriate.
#'
#' Adjustment is applied to the planned set of tests: one test for each variable
#' in `cont_vars`. A failed test remains `NA` and cannot be selected, but it is
#' still counted in the number of planned tests.
#'
#' `p_adjust_method` affects only `criterion = "pvalue"`.
#'
#' @section Adjustment variables:
#' The adjustment model is fitted once using the MFP procedure from [mfp2()].
#' The selected adjustment terms and their fitted transformations are reused in
#' every term-specific interaction model.
#'
#' The tested continuous variable is removed from its own adjustment set. The
#' group variable is also excluded from the ordinary adjustment matrix because
#' group main effects are included directly in every interaction comparison.
#'
#' Set `include_group_var = TRUE` only when the grouping variable should also be
#' forced into the first-stage MFP adjustment model. This does not control
#' whether group main effects appear in the interaction models; those main
#' effects are always included.
#'
#' The arguments `select`, `alpha`, `df`, `keep`, `force_max_fp_vars`,
#' `xorder`, and related FP settings control the MFP adjustment step. They
#' do not replace `cont_var_forms`, which controls the tested interaction form.
#'
#' If the adjustment model selects no covariates, the interaction models are
#' fitted without adjustment variables.
#'
#' @section Shift and scale settings:
#' In `mfpi.default()`, an unnamed scalar is a global setting. A named vector is
#' matched to `colnames(x)` and may specify only a subset of columns. Thus,
#' `shift = c(age = 20)` fixes only the shift for `age`, and
#' `scale = c(age = 10)` fixes only its scale. Unspecified settings remain
#' automatic and are estimated through the existing preprocessing pipeline.
#' User-supplied values must be finite, scales must be strictly positive, and
#' explicit `NA` values are not accepted.
#'
#' Explicit shifts are not increased automatically. If a supplied shift leaves
#' zero or negative values for a predictor fitted with a nonlinear FP form,
#' fitting stops and identifies the affected variable. The MFPI grouping column
#' and predictors handled as linear or structural-zero terms retain their
#' existing special shift and scale rules.
#'
#' @section Supplying data:
#' With the formula interface, supply the original variables in `data`.
#' Categorical adjustment variables can be factors and are handled as complete
#' terms. Their contrast columns are selected or excluded together. The
#' `group_var` and every variable in `cont_vars` must be included in the model
#' formula.
#'
#' With the default matrix or data-frame interface, all non-group predictors
#' must be numeric. Predictor names must be unique, non-missing, non-empty, and
#' must not contain the backtick character; spaces, hyphens, and other non-syntactic
#' names remain supported. The grouping variable may be factor, character,
#' logical, integer, or numeric when `x` is a data frame. When a categorical
#' adjustment variable is represented by several numeric columns, use
#' `term_groups` to identify those columns as one adjustment term.
#'
#' Variables in `cont_vars` must be single numeric columns. They cannot be
#' binary or members of a grouped categorical term. A variable with five or
#' fewer distinct values produces a warning because FP modelling can be
#' unreliable for a near-categorical variable.
#'
#' The reference group is determined as follows:
#'
#' - for a factor, the first factor level;
#' - for a character variable, the first distinct value in the supplied data
#'   that remains in the analysis;
#' - for logical, integer, or numeric input, the smallest observed value.
#'
#' To choose a specific reference group, convert `group_var` to a factor and set
#' its levels before fitting.
#'
#' Missing and non-finite values are not allowed in variables used for fitting.
#' Use `subset` to fit the model to selected observations. The selected rows are
#' used for model fitting and for calculating Winsorisation limits.
#'
#' Every group must retain at least two fitting observations. FP2 interactions
#' require at least three observations per group because two group-specific FP
#' columns must be estimated in addition to the group effect. Groups containing
#' exactly two observations remain permitted for linear and FP1 interactions,
#' with a warning that fitting may be unreliable.
#'
#' @section Non-positive values and special transformations:
#' Explicit `zero_vars` handling is supported for variables in `cont_vars`.
#' Non-positive values are then treated as structural zero, and the selected
#' linear or FP function is applied to the positive values. At least one
#' positive observation must remain in the fitting data after `subset` is
#' applied; otherwise fitting stops and identifies every affected variable.
#' With `center_type = "group"`, a zero-handled interaction additionally needs
#' at least two positive observations per group for a linear or FP1 form and at
#' least three for FP2. Fewer observations produce zero or rank-deficient
#' centered interaction columns.
#'
#' `catzero_vars`, `spike_vars`, and ACD transformations are supported for
#' adjustment variables but not for variables being tested in `cont_vars`.
#' When one of these settings is requested for a tested variable, it is disabled
#' with a warning. An explicit `zero_vars` setting for that variable is
#' preserved.
#'
#' Spike-at-zero adjustment variables must have enough observations in both the
#' non-positive and positive components. `min_saz_component_prop` sets the
#' minimum required proportion in each component.
#'
#' @section Winsorisation:
#' Set `winsorize = TRUE` to limit the influence of extreme values in
#' `cont_vars`. Values below and above the probabilities in `winsorize_probs`
#' are replaced by the corresponding cutoffs. No observations are removed.
#'
#' The default probabilities are `c(0.01, 0.99)`. The cutoffs are calculated
#' from the observations used for fitting, after the model's shift and scale
#' have been applied. The same cutoffs are applied automatically during
#' prediction.
#'
#' Variables that remain under structural-zero handling are not Winsorised, so
#' their structural-zero and positive components remain unchanged.
#'
#' Inspect continuous variables for influential observations even when
#' Winsorisation is used. Any Winsorisation should be reported with the analysis.
#'
#' @section Centering:
#' `center` controls whether transformed predictors are centered. For the
#' interaction functions, `center_type` determines the centering reference:
#'
#' - `"grand"` uses one mean across all observations;
#' - `"group"` uses a separate mean within each group.
#'
#' The default is `"grand"`. The fitted centering values are reused for
#' prediction and plotting.
#'
#' @section Formula offsets and Cox strata:
#' An offset can be supplied through the `offset` argument or included in the
#' formula with `offset(...)`. A formula offset takes precedence when both are
#' supplied. The offset is used in both the adjustment and interaction models.
#'
#' For Cox models, strata can be supplied through `strata` or included in the
#' formula with `strata(...)`. Formula strata take precedence when both are
#' supplied. Only right-censored [survival::Surv()] responses are supported.
#'
#' Internal response, offset, and strata columns used by the final stored models
#' are named collision-safely. User predictors named `y`, `offset_`, `strata_`,
#' or with a `..mfp2_` prefix remain ordinary predictors and are not overwritten.
#'
#' @param x For `mfpi.default()`, a numeric matrix or data frame containing the
#'   predictors. Column names must be unique, non-missing, non-empty, and must
#'   not contain backticks (`); spaces, hyphens, and other non-syntactic names
#'   are supported. Predictor values must not be missing or non-finite. Do not
#'   add an intercept column. If `x` is a data frame, only `group_var` can be
#'   categorical; all other columns must be numeric.
#'
#' @param y For `mfpi.default()`, the response. Supply a finite numeric vector
#'   for Gaussian models, non-negative counts for Poisson models, non-negative
#'   integer counts for negative-binomial models, a valid binomial response for
#'   binomial models, or a right-censored [survival::Surv()] object for Cox models. It
#'   must have the same number of observations as `x`.
#'
#' @param formula For `mfpi.formula()`, a model formula containing the response,
#'   `group_var`, the variables in `cont_vars`, and any adjustment variables.
#'   Use `fp()` or `fp2()` for variable-specific adjustment-model settings.
#'   Specify the interaction forms separately through `cont_var_forms`.
#'
#' @param data For `mfpi.formula()`, a data frame containing the variables used
#'   in `formula`.
#'
#' @param term_groups For `mfpi.default()`, an optional named list that groups
#'   numeric columns representing one categorical adjustment term. For example,
#'   `list(stage = c("stageII", "stageIII"))`. Each member column can belong
#'   to only one group. Grouped terms are fitted as fixed linear blocks: they
#'   are not eligible for FP transformation searches and cannot be used in
#'   `cont_vars`, `zero_vars`, `catzero_vars`, `spike_vars`, `acd_vars`, or
#'   `force_max_fp_vars`. In `mfpi.formula()`, factor predictors are grouped
#'   automatically; this argument is not needed.
#'
#' @param group_var A single character string naming the categorical grouping
#'   variable. It must have at least two observed groups and must not appear in
#'   `cont_vars` or `term_groups`. After applying `subset`, every group must
#'   contain at least two observations, or at least three when an FP2
#'   interaction is requested.
#'
#' @param cont_vars A non-empty character vector naming the continuous variables
#'   whose interactions with `group_var` are evaluated. Each variable must be a
#'   single numeric column with more than two distinct values.
#'
#' @param cont_var_forms An optional named character vector specifying the
#'   functional form of the interaction between each variable in `cont_vars`
#'   and the grouping variable. Permitted values are `"linear"`, `"fp1"`, and
#'   `"fp2"`. All entries must be named by the corresponding variable in
#'   `cont_vars`; variables not named default to `"fp1"`. Passing `NULL`
#'   (the default) applies `"fp1"` to every tested variable.
#'
#'   - `"linear"`: fixes the interaction at power 1 with no power search.
#'     The interaction term is a single product of the grouping indicator(s)
#'     and the untransformed continuous variable.
#'   - `"fp1"`: searches the supplied FP1 candidate powers (by default
#'     \eqn{\mathcal{P} = \{-2, -1, -0.5, 0, 0.5, 1, 2, 3\}}) and selects
#'     the single power that best fits the data. Because \eqn{p = 1} is a
#'     member of the default candidate set, a linear interaction can still be
#'     selected as the best-fitting FP1 function.
#'   - `"fp2"`: searches all ordered pairs of the supplied candidate powers,
#'     including repeated powers (e.g., \eqn{(3, 3)}), and selects the
#'     two-term FP2 combination that best fits the data. With the default
#'     eight-element power set, this evaluates 36 candidate pairs.
#'
#' @param flex A character string specifying how FP powers are chosen across
#'   groups: `"flex1"`, `"flex2"`, `"flex3"`, or `"flex4"`. The
#'   default is `"flex3"`. See **Choosing the flexibility level**.
#'
#' @param p_interact A finite number in `(0, 1]`. With
#'   `criterion = "pvalue"`, an interaction is selected only when its raw or
#'   adjusted p-value is strictly below this value. The default is `0.05`.
#'
#' @param min_improvement `NULL` or a single positive number giving the minimum
#'   AIC or BIC improvement required to select an interaction. The interaction
#'   is selected only when the relevant `AIC_main_minus_int` or
#'   `BIC_main_minus_int` is strictly greater than this value. When `NULL`
#'   (the default), the threshold is 2 for `criterion = "aic"` or `"bic"`.
#'   This argument is used only with AIC or BIC selection and has no effect
#'   when `criterion = "pvalue"`.
#'
#' @param include_group_var A single logical value. If `TRUE`, the group term is
#'   forced into the first-stage MFP adjustment model. Group main effects are
#'   included in the interaction models regardless of this setting. The default
#'   is `FALSE`.
#'
#' @param show_models A single logical value. If `TRUE` and `verbose = TRUE`,
#'   print the regression summary for each selected interaction model. The
#'   default is `FALSE`.
#'
#' @param weights An optional finite numeric vector with one value per
#'   observation. All weights must be strictly positive; zero and negative
#'   weights are not supported. This restriction keeps likelihood-based model
#'   comparisons well-defined across all supported families. The weights are
#'   used in both the adjustment and interaction models.
#'
#' @param offset An optional finite numeric vector with one value per
#'   observation, added to the linear predictor in both fitting stages. A formula
#'   `offset(...)` term takes precedence over this argument.
#'
#' @param cycles A positive integer giving the maximum number of MFP fitting
#'   cycles. The default is 10.
#'
#' @param scale For `mfpi.default()`, `NULL`, a single unnamed positive numeric
#'   value, or a named positive numeric vector for one or more columns of `x`.
#'   An unnamed scalar is applied to every column. Named values are matched to
#'   `colnames(x)`, so their order does not matter; names must be non-empty,
#'   unique, and known columns. A named partial vector fixes only the specified
#'   scales, while scales for unspecified columns are estimated automatically.
#'   Unnamed multi-value vectors are not accepted. `NULL` estimates every scale
#'   automatically and `scale = 1` disables scaling globally. Binary predictors,
#'   including factor dummy columns, are always assigned scale `1`, even when
#'   another value is supplied. The MFPI grouping column remains unscaled under the
#'   existing categorical-group handling. In `mfpi.formula()`, the top-level
#'   value is a scalar global default; use `fp()` or `fp2()` for
#'   variable-specific values.
#'
#' @param shift For `mfpi.default()`, `NULL`, a single unnamed finite numeric
#'   value, or a named finite numeric vector for one or more columns of `x`.
#'   An unnamed scalar is applied to every column. Named values are matched to
#'   `colnames(x)`, so their order does not matter; names must be non-empty,
#'   unique, and known columns. A named partial vector fixes only the specified
#'   shifts, while shifts for unspecified columns are estimated automatically.
#'   Unnamed multi-value vectors are not accepted. `NULL` estimates every shift
#'   automatically and `shift = 0` disables shifting globally. Binary
#'   predictors, including factor dummy columns, are always assigned shift `0`,
#'   even when another value is supplied. Each explicit
#'   shift used for a nonlinear FP term must make that predictor strictly
#'   positive; otherwise fitting fails and identifies the affected variable.
#'   Linear and explicitly zero-handled variables still use a shift of zero,
#'   as does the MFPI grouping column. In `mfpi.formula()`, the top-level value
#'   is a scalar global default; use `fp()` or `fp2()` for variable-specific
#'   values.
#'
#' @param df Maximum FP complexity considered in the MFP adjustment model.
#'   Valid values are 1 for linear or an even positive integer: 2 for FP1, 4
#'   for FP2, and so on. In `mfpi.default()`, an unnamed scalar is applied to
#'   all predictor columns. A named numeric vector may override one or more
#'   columns of `x`; values are matched by `colnames(x)`, order does not matter,
#'   and omitted columns use the default `df = 4`. Unnamed multi-value vectors
#'   are not accepted. Variables with few distinct values may be restricted
#'   automatically to a simpler form. `df` controls the adjustment model only
#'   and does not set the MFPI interaction form; use `cont_var_forms` for that.
#'   In `mfpi.formula()`, only a scalar is accepted; use
#'   `fp(variable, df = ...)` for per-variable values.
#'
#' @param center Controls whether transformed predictors are centered before
#'   fitting. In `mfpi.default()`, supply either a single unnamed logical value
#'   applied to every predictor column or a named logical vector for one or more
#'   columns of `x`. Named values are matched to `colnames(x)`, so their order
#'   does not matter; names must be non-empty, unique, and known columns.
#'   Omitted columns in a named partial specification use the default
#'   `center = TRUE`. Unnamed multi-value logical vectors are not accepted. In
#'   `mfpi.formula()`, only a scalar is accepted; use
#'   `fp(variable, center = ...)` for per-variable values. See `center_type` for
#'   how centering references are computed in the interaction models.
#'
#' @param subset An optional logical vector or vector of unique positive row
#'   indices selecting observations for fitting. In the formula interface, an
#'   expression using columns of `data` is also accepted. At least five
#'   observations and two groups must remain.
#'
#' @param family The model family. Use `"gaussian"`, `"binomial"`,
#'   `"poisson"`, `"negbin"`, or `"cox"`; a GLM family function; or a GLM
#'   family object. Custom GLM links are supported through family objects.
#'   Negative binomial models require `family = "negbin"` and
#'   `fitter = "fastglm"`. Cox models require `family = "cox"`.
#'
#' @param fitter Fitting backend for repeated GLM candidate fits. `"base"`
#'   (default) uses [stats::glm.fit()]. `"fastglm"` uses the compiled
#'   `fastglm::fastglm.default()` fitter and may be faster during FP candidate
#'   search. The optional `fastglm` package must be installed; ordinary GLM
#'   fits fall back to `"base"` with a warning if it is unavailable or a fit
#'   fails. Cox models are unaffected. Negative binomial models are available
#'   only with `fitter = "fastglm"` and have no base fallback.
#'
#' @param criterion The selection criterion used for both the MFP adjustment
#'   model and the interaction comparisons: `"pvalue"`, `"aic"`, or `"bic"`.
#'   With `"pvalue"`, variable retention and FP degree selection are governed
#'   by `select` and `alpha`; interactions are selected when the interaction
#'   p-value is below `p_interact`. With `"aic"` or `"bic"`, the corresponding
#'   information criterion drives both adjustment-model construction and
#'   interaction selection via `min_improvement`. The default is `"pvalue"`.
#'
#' @param select Variable-selection threshold for the MFP adjustment model.
#'   Values must lie in `[0, 1]`. In `mfpi.default()`, an unnamed scalar is
#'   applied to all predictor columns. A named numeric vector may override one
#'   or more columns of `x`; values are matched by `colnames(x)`, order does not
#'   matter, and omitted columns use the default `select = 0.05`. Unnamed
#'   multi-value vectors are not accepted. With `criterion = "pvalue"`, setting
#'   `select = 1` forces a variable into the model. With AIC or BIC selection,
#'   `select` is not used directly for thresholding but variables named in
#'   `force_max_fp_vars` are still forced to `select = 1`. In `mfpi.formula()`,
#'   only a scalar is accepted; use `fp(variable, select = ...)` for
#'   per-variable values.
#'
#' @param alpha Significance level for FP-degree selection in the MFP
#'   adjustment model. Values must lie in `[0, 1]`. In `mfpi.default()`, an
#'   unnamed scalar is applied to all predictor columns. A named numeric vector
#'   may override one or more columns of `x`; values are matched by
#'   `colnames(x)`, order does not matter, and omitted columns use the default
#'   `alpha = 0.05`. Unnamed multi-value vectors are not accepted. It is used
#'   only when `criterion = "pvalue"`; setting `alpha = 1` forces the maximum
#'   FP degree allowed by `df`. In `mfpi.formula()`, only a scalar is accepted;
#'   use `fp(variable, alpha = ...)` for per-variable values.
#'
#' @param keep An optional character vector naming adjustment variables or
#'   grouped adjustment terms that must be included in the adjustment model.
#'   For `mfpi.default()`, grouped-term names created by `term_groups` are
#'   accepted. For `mfpi.formula()`, factor variable names are accepted and
#'   all contrast columns of that factor are retained together.
#'
#' @param force_max_fp_vars An optional character vector naming adjustment
#'   variables that must use the maximum FP degree allowed by `df`. The best
#'   powers within that degree are still chosen by the data. Under
#'   `criterion = "pvalue"`, this also forces `select = 1` and `alpha = 1` for
#'   the named variables, so they are both retained and kept at full FP
#'   complexity. Under AIC and BIC the same maximum-degree requirement is
#'   enforced directly. If a forced adjustment variable is also an eligible
#'   spike-at-zero term, its complete maximum SAZ representation is retained
#'   and SAZ Stage 2 is skipped. In a formula, the corresponding option is
#'   `fp(variable, force_max_fp = TRUE)` or
#'   `fp2(variable, force_max_fp = TRUE)`.
#'
#' @param xorder The order in which variables are processed by the MFP
#'   adjustment algorithm: `"ascending"` orders by increasing univariate
#'   significance, `"descending"` by decreasing significance, and `"original"`
#'   preserves the column order of `x` or the formula term order. The default
#'   is `"ascending"`.
#'
#' @param powers An optional named list of candidate FP power sets. The default
#'   uses `c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)`, where zero represents the
#'   natural logarithm. These candidates are used when FP powers are selected
#'   for both adjustment and interaction functions.
#'
#' @param ties The method for handling tied event times in Cox models. Supported
#'   values are `"breslow"` (default) and `"efron"`. `"exact"` is not supported
#'   because MFPI relies on MFP-based Cox candidate fits that do not implement
#'   the exact partial likelihood. The argument otherwise has no effect for
#'   non-Cox families.
#'
#' @param strata Optional Cox strata with one value or row per observation. A
#'   vector or factor supplies one categorical stratum label per observation;
#'   character, numeric, integer, and logical values are accepted. A matrix or
#'   data frame may supply multiple stratification variables, one per column,
#'   which are combined before Cox fitting. A [survival::strata()] object is
#'   also accepted. Formula `strata(...)` terms take precedence over this argument.
#'
#' @param nocenter Numeric set of values used to identify Cox design-matrix
#'   columns that should not be internally recentered. A column is left
#'   uncentered when all of its values are contained in `nocenter`. The default
#'   `c(-1, 0, 1)` matches [survival::coxph()] and typically leaves indicator
#'   and dummy columns uncentered. Set `NULL` to allow all eligible columns to
#'   be recentered. It is used only for Cox models.
#'
#' @param acd_vars For `mfpi.default()`, an optional character vector naming
#'   adjustment variables that use the approximate cumulative distribution
#'   (ACD) transformation instead of a standard FP transformation. ACD is
#'   disabled with a warning for variables in `cont_vars` because it would
#'   conflict with the interaction-form search. In `mfpi.formula()`, this
#'   argument is not available; use `fp(variable, acd = TRUE)` or
#'   `fp2(variable, acd = TRUE)` in the formula instead.
#'
#' @param zero_vars An optional character vector naming variables for which
#'   non-positive values are treated as structural zeros. When enabled, the
#'   selected FP or linear function is applied only to positive values, while
#'   non-positive values receive a separate fixed contribution. The shift for
#'   a zero-handled variable is set to 0. This option can be used for
#'   variables in `cont_vars` and for adjustment variables. In the formula
#'   interface, this can also be specified with `fp(variable, zero = TRUE)`.
#'   A zero-handled variable in `cont_vars` must retain at least one positive
#'   observation after applying `subset`.
#'
#' @param catzero_vars An optional character vector naming adjustment variables
#'   that use both a positive-part FP function and a binary structural-zero
#'   indicator as an additional covariate. Setting `catzero_vars` implies
#'   `zero_vars` for the same variables. This option is disabled with a
#'   warning for variables in `cont_vars`; only explicit `zero_vars` settings
#'   are preserved for those. In the formula interface, use
#'   `fp(variable, catzero = TRUE)`.
#'
#' @param spike_vars An optional character vector naming adjustment variables to
#'   assess with the spike-at-zero procedure. Spike-at-zero testing requires
#'   enough observations in both the non-positive and positive components (see
#'   `min_saz_component_prop`). Setting `spike_vars` implies both
#'   `catzero_vars` and `zero_vars` for the same variables. Variables that fail
#'   the component-proportion check are automatically downgraded to `catzero`
#'   or `zero` handling. This option is disabled with a warning for variables
#'   in `cont_vars`. In the formula interface, use
#'   `fp(variable, spike = TRUE)`.
#'
#' @param min_saz_component_prop A finite number in `(0, 0.5)` giving the
#'   minimum proportion required in both components of a spike-at-zero
#'   adjustment variable. The default is `0.10`.
#'
#' @param ftest A single logical value. For Gaussian models with
#'   `criterion = "pvalue"`, use F-tests instead of chi-square tests in the
#'   adjustment and interaction steps. The default is `FALSE`. It is ignored
#'   for non-Gaussian families.
#'
#' @param control Control settings for the underlying fit, usually from
#'   [stats::glm.control()] or [survival::coxph.control()]. `NULL` uses the
#'   relevant defaults. For GLMs, the same controls are used during MFP
#'   selection and the final fit. With `fitter = "fastglm"`, `epsilon` is
#'   mapped to `tol` and `maxit` is passed through; `trace = TRUE` is not
#'   supported by the fastglm backend. For `family = "negbin"`, these settings
#'   apply to the inner IRLS fit; fastglm_nb-specific outer controls retain
#'   their defaults.
#'
#' @param winsorize A single logical value. If `TRUE`, Winsorise variables in
#'   `cont_vars` before fitting. The default is `FALSE`.
#'
#' @param winsorize_probs A numeric vector of length two giving the lower and
#'   upper probabilities used for Winsorisation. The values must be in `[0, 1]`
#'   and in increasing order. The default is `c(0.01, 0.99)`.
#'
#' @param center_type How interaction functions are centered when
#'   `center = TRUE`: `"grand"` or `"group"`. The default is `"grand"`. For a
#'   group-centered zero-handled interaction, each group must retain enough
#'   positive observations to identify the requested interaction basis: two
#'   for a linear or FP1 form and three for FP2.
#'
#' @param p_adjust_method A character string accepted by [stats::p.adjust()].
#'   It controls multiplicity adjustment across the variables in `cont_vars`
#'   when `criterion = "pvalue"`. The default is `"none"`.
#'
#' @param verbose A single logical value indicating whether progress information
#'   is printed. The default is `TRUE`.
#'
#' @param digits A positive integer controlling the number of significant digits
#'   used in printed interaction results. The default is 3.
#'
#' @param ... Additional arguments are not currently supported and produce an
#'   error.
#'
#' @return
#' An object of class `"mfpi"`. Use `print()`, `summary()`, `coef()`, `vcov()`,
#' `predict()`, and `plot()` for the main results.
#'
#' Important components include:
#'
#' - `best_model_metrics`: a data frame of results for interactions selected
#'   by the chosen criterion. Contains columns such as the deviance, degrees
#'   of freedom, p-values, AIC and BIC comparisons, and selected powers;
#' - `all_model_metrics`: the same columns for every variable evaluated in
#'   `cont_vars`, regardless of whether its interaction was selected;
#' - `best_interaction_model`: a named list of fitted model objects for
#'   selected interactions;
#' - `all_interaction_models`: a named list of fitted model objects for every
#'   evaluated continuous variable, including interactions that were not
#'   selected;
#' - `adjustment_model`: the fitted [mfp2()] adjustment model. Its
#'   preprocessing metadata are restored to the MFPI-level shifts after Stage 1,
#'   so printing reports the shifts applied to the raw covariates and direct
#'   `predict()` calls on this nested model correctly preprocess raw `newdata`;
#' - `cont_var_forms`: the interaction form (`"linear"`, `"fp1"`, or `"fp2"`)
#'   specified for each tested variable;
#' - `group_level_map`: the mapping between internal group codes and the
#'   original group labels;
#' - `winsorize_limits`: the Winsorisation limits on the processed scale used
#'   for fitting, when Winsorisation was requested;
#' - `criterion`: the selection criterion used;
#' - `p_adjust_method`: the multiplicity adjustment method applied;
#' - `call`: the matched model call.
#'
#' Additional components are stored for prediction, plotting, and summaries and
#' should not be treated as a stable public interface.
#'
#' @references
#' Royston, P. and Sauerbrei, W. (2004). A new approach to modelling
#' interactions between treatment and continuous covariates in clinical trials
#' by using fractional polynomials. *Statistics in Medicine*, 23, 2509-2525.
#'
#' Royston, P. and Sauerbrei, W. (2013). Interaction of treatment with a
#' continuous variable: simulation study of significance level for several
#' methods of analysis. *Statistics in Medicine*, 32, 3788-3803.
#'
#' Royston, P. and Sauerbrei, W. (2014). Interaction of treatment with a
#' continuous variable: simulation study of power for several methods of
#' analysis. *Statistics in Medicine*, 33, 4695-4708.
#'
#' @examples
#' # --- Formula interface: test interactions with fp() adjustment terms ---
#' data("prostate")
#'
#' fit <- mfpi(
#'   lpsa ~ svi + fp(age) + fp(cavol) + fp(weight) + fp(bph) + fp(cp),
#'   data = prostate,
#'   group_var = "svi",
#'   cont_vars = c("cavol", "age"),
#'   cont_var_forms = c(cavol = "fp1", age = "linear"),
#'   flex = "flex1",
#'   verbose = FALSE
#' )
#'
#' # Interaction results for all evaluated variables.
#' fit$all_model_metrics
#'
#' # Results and fitted models for selected interactions.
#' fit$best_model_metrics
#' summary(fit)
#'
#' # Fitted group functions and their differences.
#' pred <- predict(
#'   fit,
#'   terms = "cavol",
#'   model = "all",
#'   type = "both",
#'   grid = TRUE
#' )
#'
#' pred$functions
#' pred$differences
#'
#' \donttest{
#' if (requireNamespace("patchwork", quietly = TRUE)) {
#'   plot(fit, terms = "cavol", model = "all", plot_type = "both")
#' }
#' }
#'
#' # --- Formula with factor adjustment term ---
#' set.seed(1)
#' dat <- data.frame(
#'   y = rnorm(180),
#'   trt = factor(rep(c("control", "active"), each = 90)),
#'   age = runif(180, 30, 80),
#'   stage = factor(rep(c("I", "II", "III"), length.out = 180))
#' )
#'
#' fit_factor <- mfpi(
#'   y ~ trt + fp(age) + stage,
#'   data = dat,
#'   group_var = "trt",
#'   cont_vars = "age",
#'   cont_var_forms = c(age = "fp1"),
#'   keep = "stage",
#'   verbose = FALSE
#' )
#'
#' summary(fit_factor)
#'
#' # --- Matrix interface with term_groups and explicit shift/scale ---
#' stage_matrix <- stats::model.matrix(~ stage, data = dat)[, -1, drop = FALSE]
#' x <- cbind(
#'   trt = as.integer(dat$trt) - 1L,
#'   age = dat$age,
#'   stage_matrix
#' )
#'
#' fit_matrix <- mfpi(
#'   x = x,
#'   y = dat$y,
#'   group_var = "trt",
#'   cont_vars = "age",
#'   cont_var_forms = c(age = "fp1"),
#'   term_groups = list(stage = colnames(stage_matrix)),
#'   keep = "stage",
#'   shift = c(age = 20),
#'   scale = c(age = 10),
#'   verbose = FALSE
#' )
#'
#' # Only age uses the explicit settings; other adjustment columns remain
#' # automatic and the grouping column retains MFPI's categorical defaults.
#' fit_matrix$shift
#' fit_matrix$scale
#' summary(fit_matrix)
#'
#' @seealso [mfp2()], [fp()], [coef.mfpi()], [vcov.mfpi()], [predict.mfpi()],
#'   [plot.mfpi()], [summary.mfpi()], [print.mfpi()]
#'
#' @export
mfpi <- function(x, ...) {
  UseMethod("mfpi", x)
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
    cont_vars         = NULL,
    cont_var_forms    = NULL,
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
    ties              = c("breslow", "efron"),
    strata            = NULL,
    nocenter          = c(-1, 0, 1),
    acd_vars          = NULL,
    zero_vars         = NULL,
    catzero_vars      = NULL,
    spike_vars        = NULL,
    min_saz_component_prop = 0.10,
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
  # adjustment model and test each cont_vars interaction in turn.
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
  fitter <- resolve_fitter(fitter, family_string)

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
  # Step 3: Validate group_var, p_adjust_method, and cont_vars ----------------
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

  if (group_var %in% cont_vars) {
    stop(
      paste0("! `group_var = '", group_var,
             "'` must not also appear in `cont_vars`."),
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

  # Validate cont_vars ---------------------------------------------------------
  if (is.null(cont_vars) || !is.character(cont_vars) || length(cont_vars) == 0L) {
    stop(
      paste0(
        "! `cont_vars` must be a non-empty character vector naming continuous ",
        "variables to test for interaction with `group_var = '", group_var, "'`."
      ),
      call. = FALSE
    )
  }

  missing_cont <- setdiff(cont_vars, vnames)
  if (length(missing_cont) > 0L) {
    stop(
      paste0("! The following variables in `cont_vars` are not columns of `x`: ",
             paste(missing_cont, collapse = ", "), "."),
      call. = FALSE
    )
  }

  # MFPI interaction variables and the grouping variable must each be one raw
  # identity-mapped column. Explicit mappings are adjustment covariates only,
  # including binary factors represented by one non-identity dummy column.
  grouped_terms <- names(term_to_columns)[mapped_term_flags(term_to_columns)]
  grouped_columns <- unlist(term_to_columns[grouped_terms], use.names = FALSE)

  if (group_var %in% grouped_terms || group_var %in% grouped_columns) {
    stop("! `group_var` must not be included in `term_groups`.", call. = FALSE)
  }

  grouped_cont <- cont_vars[
    cont_vars %in% grouped_terms | cont_vars %in% grouped_columns
  ]
  if (length(grouped_cont) > 0L) {
    stop(
      paste0(
        "! Variables in `cont_vars` must be singleton continuous columns and ",
        "cannot be grouped categorical terms: ",
        paste(grouped_cont, collapse = ", "), "."
      ),
      call. = FALSE
    )
  }

  # Validate and expand cont_var_forms ------------------------------------------
  # cont_var_forms must be either NULL (default: all fp1) or a named
  # character vector whose names are a subset of cont_vars and whose values
  # are each one of "linear", "fp1", "fp2". Missing cont_vars entries are
  # silently filled with "fp1". Unnamed vectors and unnamed scalars are
  # rejected to enforce explicit per-variable specification.
  valid_forms <- c("linear", "fp1", "fp2")

  if (is.null(cont_var_forms)) {
    cont_var_forms <- stats::setNames(
      rep("fp1", length(cont_vars)), cont_vars
    )
  } else {
    if (!is.character(cont_var_forms)) {
      stop(
        "! `cont_var_forms` must be a named character vector or NULL.",
        call. = FALSE
      )
    }
    if (is.null(names(cont_var_forms)) || any(names(cont_var_forms) == "")) {
      stop(
        paste0(
          "! Every entry of `cont_var_forms` must be named with a variable ",
          "from `cont_vars`. Unnamed entries are not allowed.\n",
          "  Example: cont_var_forms = c(age = \"fp2\", bmi = \"fp1\")"
        ),
        call. = FALSE
      )
    }
    unknown_names <- setdiff(names(cont_var_forms), cont_vars)
    if (length(unknown_names) > 0L) {
      stop(
        paste0(
          "! The following name(s) in `cont_var_forms` are not in `cont_vars`: ",
          paste(unknown_names, collapse = ", "), "."
        ),
        call. = FALSE
      )
    }
    bad_values <- cont_var_forms[!cont_var_forms %in% valid_forms]
    if (length(bad_values) > 0L) {
      stop(
        paste0(
          "! Invalid value(s) in `cont_var_forms`: ",
          paste(unique(bad_values), collapse = ", "), ".\n",
          "  Allowed values are: \"linear\", \"fp1\", \"fp2\"."
        ),
        call. = FALSE
      )
    }
    # Fill any cont_vars not mentioned with "fp1"
    missing_vars <- setdiff(cont_vars, names(cont_var_forms))
    if (length(missing_vars) > 0L) {
      cont_var_forms <- c(
        cont_var_forms,
        stats::setNames(rep("fp1", length(missing_vars)), missing_vars)
      )
    }
    # Re-order to match cont_vars order
    cont_var_forms <- cont_var_forms[cont_vars]
  }

  # Check cont_vars are numeric -------------------------------------------------
  non_numeric <- cont_vars[!vapply(cont_vars, function(v)
    is.numeric(x[, v, drop = TRUE]), logical(1L))]
  if (length(non_numeric) > 0L) {
    stop(
      paste0(
        "! The following variable(s) in `cont_vars` are not numeric: ",
        paste(non_numeric, collapse = ", "), ".\n",
        "  FP transformation requires continuous numeric variables. ",
        "Convert categorical variables to numeric dummy columns before ",
        "calling `mfpi()`."
      ),
      call. = FALSE
    )
  }

  # Check cont_vars have sufficient unique values --------------------------------
  # Binary or near-categorical variables produce degenerate FP transformations.
  # `x` has already been checked above to contain no missing or non-finite values.
  # Count unique observed values directly.
  n_unique <- vapply(cont_vars, function(v) {
    length(unique(preprocess_x[, v, drop = TRUE]))
  }, integer(1L))


  binary_cont <- cont_vars[n_unique <= 2L]
  few_unique  <- cont_vars[n_unique > 2L & n_unique <= 5L]

  if (length(binary_cont) > 0L) {
    stop(
      paste0(
        "! The following variable(s) in `cont_vars` are binary (<=2 unique ",
        "values): ", paste(binary_cont, collapse = ", "), ".\n",
        "  Binary variables cannot be FP-transformed. Use `group_var` for ",
        "binary treatment indicators, or omit these variables from `cont_vars`."
      ),
      call. = FALSE
    )
  }

  if (length(few_unique) > 0L) {
    warning(
      paste0(
        "! The following variable(s) in `cont_vars` have 5 or fewer unique ",
        "values: ", paste(few_unique, collapse = ", "), ".\n",
        "  FP transformation may be unreliable for near-categorical variables. ",
        "Consider whether these variables are truly continuous."
      ),
      call. = FALSE
    )
  }

  # Step 4: Validate the response, Cox auxiliary inputs, weights, and offset --
  # Validate response y --------------------------------------------------------
  validate_family_response(
    y             = y,
    family_string = family_string,
    nobs          = nobs
  )

  # Validate Cox-specific auxiliary inputs ------------------------------------
  if (family_string == "cox" && !is.null(strata)) {
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
    if (
      !is.numeric(offset) ||
      length(offset) != nobs ||
      anyNA(offset) ||
      any(!is.finite(offset))
    ) {
      stop(
        "! `offset` must be a finite numeric vector of length nobs.",
        call. = FALSE
      )
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

  # For cont_vars, catzero/spike will be dropped later, so do not treat
  # zero + catzero overlap on cont_vars as fatal.
  overlap <- intersect(zero_vars_in, setdiff(catzero_vars_in, cont_vars))

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

  if (is.null(control)) {
    control <- if (family_string == "cox") survival::coxph.control() else
      stats::glm.control()
  }

  # Retain valid raw-column or conceptual-term entries. preprocess_data()
  # converts these to the term names expected by the grouped MFP core.
  keep <- intersect(keep, valid_term_names)

  # Step 10: Build zero/catzero/spike flags, suppress them for cont_vars, and
  # resolve spike-at-zero eligibility for the remaining adjustment variables -
  # Process zero_vars / catzero_vars / spike_vars to named logical vectors -----
  zero_flag    <- setNames(rep(FALSE, nvars), vnames)
  catzero_flag <- setNames(rep(FALSE, nvars), vnames)
  spike_flag   <- setNames(rep(FALSE, nvars), vnames)

  if (length(zero_vars)    > 0L) zero_flag[zero_vars]       <- TRUE
  if (length(catzero_vars) > 0L) catzero_flag[catzero_vars] <- TRUE
  if (length(spike_vars)   > 0L) spike_flag[spike_vars]     <- TRUE

  # cont_vars may keep explicit zero handling, but cannot use catzero/spike.
  # Developer note: this cleanup must happen BEFORE the cascade below. Otherwise
  # spike_vars/catzero_vars would imply zero_vars for cont_vars, which is not the
  # intended MFPI stage-2 contract. Only explicit zero_vars should survive for
  # interaction variables.
  spike_in_cont   <- intersect(names(spike_flag)[spike_flag], cont_vars)
  catzero_in_cont <- intersect(names(catzero_flag)[catzero_flag], cont_vars)

  if (length(spike_in_cont) > 0L) {
    warning(
      "The following `cont_vars` were also marked as `spike_vars`: ",
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
      "The following `cont_vars` were also marked as `catzero_vars`: ",
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
      min_saz_component_prop = min_saz_component_prop
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
  # and build the acd_vars flag (suppressed for cont_vars) -------------------
  # Reset zero/catzero/spike for variables that contain only positive values -----
  if (any(zero_flag | catzero_flag | spike_flag)) {
    vars_pos <- names(zero_flag)[zero_flag | catzero_flag | spike_flag]

    xpos <- preprocess_x[, vars_pos, drop = FALSE]
    all_positive <- vars_pos[colSums(xpos <= 0) == 0L]

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
    acd_in_cont <- intersect(acd_vars, cont_vars)

    if (length(acd_in_cont) > 0L) {
      warning(
        "ACD transformation is disabled for cont_vars: ",
        paste(acd_in_cont, collapse = ", "),
        ".",
        call. = FALSE
      )

      acd_vars <- setdiff(acd_vars, cont_vars)
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

  # Step 14: Prepare Cox stratification, apply `subset`, and re-validate
  # group_var against the (possibly reduced) subsetted data -------------------
  # Prepare Cox stratification ------------------------------------------------
  # Normalize all strata to one factor before MFP/MFPI candidate fitting.
  # Multiple columns are combined as in survival::coxph(), while a single
  # vector is treated as categorical labels regardless of storage mode.
  # Integer conversion belongs only at the survival::coxph.fit() boundary.
  strata_keep <- strata

  if (family_string == "cox" && !is.null(strata_keep)) {
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

    y <- if (family_string == "cox" || is.matrix(y)) {
      y[subset, , drop = FALSE]
    } else {
      y[subset]
    }
    weights <- weights[subset]
    offset  <- offset[subset]
    if (!is.null(strata_keep)) {
      # Subsetting can remove complete strata, so drop unused levels after the
      # row restriction while preserving row alignment.
      strata_keep <- droplevels(strata_keep[subset])
    }
    validate_subset_predictor_variation(x, exclude = group_var)
  }

  # zero_vars applies the interaction function only to x > 0. Validate this
  # requirement after the optional subset has been applied, because a variable
  # can contain positive values in the full data yet lose all of them in the
  # actual fitting population. This shared boundary gives every flexibility
  # method the same variable-specific error instead of exposing flex2/flex4
  # transformation assertions to users.
  zero_cont_vars <- intersect(
    cont_vars,
    names(zero_flag)[zero_flag]
  )

  if (length(zero_cont_vars) > 0L) {
    has_positive <- vapply(
      zero_cont_vars,
      function(variable) any(x[, variable] > 0),
      logical(1L)
    )

    if (any(!has_positive)) {
      stop(
        "Zero-handled variables in `cont_vars` must contain at least one ",
        "positive value in the fitting data after subsetting. Problematic ",
        "variables: ",
        paste(zero_cont_vars[!has_positive], collapse = ", "),
        ".",
        call. = FALSE
      )
    }
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
  fp2_requested <- any(cont_var_forms == "fp2")
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

  if (identical(center_type, "group") && length(zero_cont_vars) > 0L) {
    group_values <- x[, group_var]
    fitted_group_levels <- sort(unique(group_values))
    insufficient_positive <- character(0L)

    for (variable in zero_cont_vars) {
      required_positive <- if (identical(cont_var_forms[[variable]], "fp2")) {
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

  # Step 15: Winsorise cont_vars (unless excluded by zero handling) ----------
  # Winsorise cont_vars to reduce influence of extreme values ------------------
  # Variables with final zero handling are excluded from Winsorisation so that
  # structural-zero values and the positive distribution are both preserved.
  winsorize_limits <- NULL
  if (isTRUE(winsorize) && length(cont_vars) > 0L) {
    exclude_from_win <- names(zero_flag)[zero_flag]
    win <- winsorize_cont_vars(
      x         = x,
      cont_vars = cont_vars,
      probs     = winsorize_probs,
      zero_vars = exclude_from_win
    )
    x                <- win$x
    winsorize_limits <- win$limits
    if (verbose) {
      message(sprintf(
        "Winsorising %d cont_var(s) at percentiles (%g, %g).",
        length(cont_vars), winsorize_probs[1L], winsorize_probs[2L]
      ))

      excluded <- intersect(exclude_from_win, cont_vars)
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

  # Step 16: Fit the MFPI model and attach return metadata --------------------
  # Fit the MFPI model ---------------------------------------------------------
  fit <- fit_mfpi(
    x                 = x,
    y                 = y,
    family            = family,         # GLM object or "cox"
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
    min_saz_component_prop = min_saz_component_prop,
    use_ftest         = ftest,
    control           = control,
    verbose           = verbose,
    group_var         = group_var,
    include_group_var = include_group_var,
    flex              = flex,
    cont_vars         = cont_vars,
    cont_var_forms    = cont_var_forms,
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
  fit$cont_var_forms    <- cont_var_forms

  # `fit_mfpi()` stores `x_train_internal`, the post-preprocessing training
  # matrix used for manual prediction reconstruction. It contains shifted/scaled
  # predictors and internally remapped group_var codes. Do not store or use the
  # pre-remap matrix `x` for prediction.
  fit$family_string  <- family_string
  fit$zero_vars      <- zero_flag
  fit$center_type    <- center_type
  fit$cont_vars      <- cont_vars
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

#' @describeIn mfpi Formula interface for \code{mfpi()}.
#'
#' Constructs a model frame from the supplied formula and data, processes
#' predictor variables (including any \code{fp()} or \code{fp2()} terms), and fits an MFPI
#' model by calling \code{mfpi.default()} internally.
#'
#' @section Argument precedence (fp()/fp2() vs global):
#' Arguments supplied inside \code{fp()} or \code{fp2()} are treated as
#' variable-specific settings and take precedence over the corresponding global
#' arguments passed to \code{mfpi.formula()}. Variables not wrapped in
#' \code{fp()} or \code{fp2()} use the global argument values. For
#' \code{zero_vars}, \code{catzero_vars}, \code{spike_vars}, and
#' \code{force_max_fp_vars}, global settings are combined with the
#' corresponding \code{fp()} or \code{fp2()} flags and then checked by
#' \code{mfpi.default()} for unsupported combinations. Cardinality-based
#' adjustments to \code{df} may still be applied internally.
#'
#' @section Parameters restricted to scalar values in the formula interface:
#' The following \code{mfpi.default()} parameters must be supplied as scalar
#' arguments only; per-variable specification via \code{fp()} or \code{fp2()} is not
#' supported: \code{group_var}, \code{cont_vars}, \code{cont_var_forms},
#' \code{flex}, \code{p_interact}, \code{min_improvement},
#' \code{include_group_var}, \code{show_models}, \code{cycles},
#' \code{criterion}, \code{keep}, \code{xorder}, \code{ties}, \code{strata},
#' \code{nocenter}, \code{min_saz_component_prop}, \code{ftest},
#' \code{control}, \code{verbose}, \code{digits}.
#'
#' \code{acd_vars}, by contrast, is not merely scalar-restricted but entirely
#' unavailable as a direct argument to \code{mfpi.formula()} (passing it
#' raises an error); use \code{fp(variable, acd = TRUE)} or
#' \code{fp2(variable, acd = TRUE)} in the formula instead.
#'
#' @section Categorical terms in the formula interface:
#' Unordered and ordered factor predictors other than \code{group_var} are
#' expanded with the contrasts recorded by \code{model.matrix()}. All columns
#' generated by one factor are stored as one conceptual adjustment term and are
#' selected jointly. Do not wrap a categorical predictor in \code{fp()} or
#' \code{fp2()}; those wrappers are for continuous variables. The factor term
#' may be named in \code{keep}, but it cannot be named in \code{cont_vars}.
#' The fitted terms, contrasts, levels, and term-to-column mapping are stored so
#' \code{predict.mfpi()} can accept raw factor-valued \code{newdata}.
#' @section Strata and offset in the formula:
#' For Cox models, \code{strata()} terms may be included directly in the formula,
#' for example \code{Surv(t, d) ~ fp(age) + strata(centre)}. Multiple strata
#' variables are combined using \code{survival::strata(shortlabel = TRUE)}. If
#' \code{strata} is also supplied as an argument, the formula value is used and
#' a warning is issued. Similarly, an \code{offset()} term in the
#' formula takes precedence over the \code{offset} argument. Formula-offset
#' metadata are stored on the fitted object so \code{predict.mfpi()} can
#' reconstruct the offset from raw prediction \code{newdata}; if the needed raw
#' offset variables are absent, supply \code{newoffset} to \code{predict()}.
#'
#' @method mfpi formula
#' @export
mfpi.formula <- function(formula,
                         data,
                         group_var         = NULL,
                         cont_vars         = NULL,
                         cont_var_forms    = NULL,
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
                         ties              = c("breslow", "efron"),
                         strata            = NULL,
                         nocenter          = c(-1, 0, 1),
                         zero_vars         = NULL,
                         catzero_vars      = NULL,
                         spike_vars        = NULL,
                         min_saz_component_prop = 0.10,
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

  # Step 1: Resolve multiple-choice arguments, family, and reject unused `...`
  dots <- match.call(expand.dots = FALSE)$...
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

  # Unused/misspelled arguments (including acd_vars, which is not supported
  # in the formula interface - see fp(variable, acd = TRUE) instead) are
  # already rejected above by mfpi_check_unused_dots(), which reports the
  # offending argument name(s) directly.

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


  # Step 4: Validate cont_vars against the *original* data types, before
  # model.matrix() expands factors and would otherwise produce a confusing
  # "variable not found" error instead of a clear "categorical" message -----
  # ---------------------------------------------------------------------------
  # Validate cont_vars against original data types (before model.matrix)
  # ---------------------------------------------------------------------------
  # Check here using original data types because model.matrix dummy-codes
  # factors, which would cause a misleading "variable not found" error later
  # rather than a clear "this is categorical" message.
  if (!is.null(cont_vars)) {
    # Factor or character columns cannot be FP-transformed
    factor_in_cont <- cont_vars[vapply(cont_vars, function(v)
      v %in% names(data) && (is.factor(data[[v]]) || is.character(data[[v]])),
      logical(1L))]
    if (length(factor_in_cont) > 0L) {
      stop(
        "! The following variable(s) in `cont_vars` are categorical ",
        "(factor or character) in `data`: ",
        paste(factor_in_cont, collapse = ", "), ".\n",
        "  FP transformation requires continuous numeric variables. ",
        "Categorical variables cannot be tested for interaction via MFPI.",
        call. = FALSE
      )
    }

    # Numeric binary columns (0/1 or two unique values) in the original data
    binary_in_cont <- cont_vars[vapply(cont_vars, function(v)
      v %in% names(data) &&
        is.numeric(data[[v]]) &&
        length(unique(stats::na.omit(data[[v]]))) <= 2L,
      logical(1L))]
    if (length(binary_in_cont) > 0L) {
      stop(
        "! The following variable(s) in `cont_vars` are binary ",
        "(2 or fewer unique values) in `data`: ",
        paste(binary_in_cont, collapse = ", "), ".\n",
        "  Binary variables cannot be FP-transformed. ",
        "Use `group_var` for binary treatment indicators, or omit these ",
        "variables from `cont_vars`.",
        call. = FALSE
      )
    }
  }

  # Evaluate `subset` once, then validate observation-level inputs.
  n_data <- nrow(data)

  # Match standard formula-method lookup rules: data columns are searched first,
  # followed by the normalized formula environment (which inherits the user's
  # formula environment). Rebinding `subset` to the resolved value prevents any
  # later branch from re-evaluating an expression with side effects.
  subset_expr <- substitute(subset)
  subset <- eval(
    subset_expr,
    envir = data,
    enclos = environment(formula_internal)
  )

  if (!is.null(offset)) {
    if (!is.numeric(offset) || length(offset) != n_data || anyNA(offset) ||
        any(!is.finite(offset))) {
      stop("! `offset` must be a finite numeric vector with one value per row of `data`.", call. = FALSE)
    }
  }
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

  # Step 6: Extract strata() terms for Cox models -----------------------------
  # ---------------------------------------------------------------------------
  # Handle strata terms for Cox models
  # ---------------------------------------------------------------------------
  specials      <- "strata"
  terms_formula <- stats::terms(
    formula_internal,
    specials = specials,
    data = data
  )
  if (!is.null(strata) && is.null(attr(terms_formula, "specials")$strata)) {
    strata_n <- if (is.vector(strata) || is.factor(strata)) length(strata) else NROW(strata)
    if (strata_n != n_data || anyNA(strata)) {
      stop("! `strata` must have one non-missing value or row per row of `data`.", call. = FALSE)
    }
  }

  terms_drop <- NULL
  formula_strata_terms <- NULL
  formula_strata_xlevels <- NULL
  formula_offset_terms <- NULL
  formula_offset_xlevels <- NULL

  if (!is.null(attr(terms_formula, "specials")$strata)) {
    if (family_string == "cox") {
      if (!is.null(call$strata))
        warning("i strata appear in both the formula and as an argument.\n",
                "i Formula strata are used; the argument is ignored.",
                call. = FALSE)

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
      stop("! strata are only allowed for Cox models.\n",
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
  term_offset <- attr(terms_formula, "offset")
  if (!is.null(term_offset) && length(term_offset) > 1L)
    stop("! Only one offset term is allowed in the formula.", call. = FALSE)

  if (!is.null(term_offset)) {
    if (!is.null(call$offset))
      warning("i Offset appears in both the formula and as an argument.\n",
              "i Formula offset is used; the argument is ignored.", call. = FALSE)
    # Always extract offset from model frame when formula contains offset().
    # Also keep the offset expression itself so predict.mfpi() can rebuild
    # expressions such as offset(log(exposure)) from raw newdata.
    offset <- as.vector(stats::model.offset(mf))

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

  if (family_string != "cox" && survival::is.Surv(y))
    stop(sprintf(
      "! Response is a Surv object but family = '%s'. Set family = 'cox'.",
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
  # Validate singleton-column variation before delegation because subset has
  # already been applied. The private attribute transports only the aligned
  # full-data preprocessing source and is removed on entry to mfpi.default().
  # All fitted categorical metadata comes from mf/x, and subset = NULL prevents
  # the default method from applying the row selection a second time.
  if (!is.null(subset)) {
    validate_subset_predictor_variation(x, exclude = group_var)
  }
  x <- attach_formula_preprocess_matrix(x, x_full)
  if (!is.null(subset)) {
    if (!is.null(weights)) weights <- weights[fit_rows]
    if (is.null(term_offset) && !is.null(offset)) offset <- offset[fit_rows]
    if (!is.null(strata) && is.null(attr(terms_formula, "specials")$strata)) {
      strata <- if (is.matrix(strata) || is.data.frame(strata)) {
        strata[fit_rows, , drop = FALSE]
      } else {
        strata[fit_rows]
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Delegate to mfpi.default
  # ---------------------------------------------------------------------------
  fit <-  mfpi.default(
    x                 = x,
    y                 = y,
    group_var         = group_var,
    cont_vars         = cont_vars,
    cont_var_forms    = cont_var_forms,
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
    min_saz_component_prop = min_saz_component_prop,
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

#' Prepare predictor input for \code{mfpi.default()}
#'
#' Converts predictor input supplied to \code{mfpi.default()} into the numeric
#' matrix representation required by the MFPI fitting internals, while
#' preserving user-facing labels for the grouping variable.
#'
#' The default method supports two input forms:
#'
#' \describe{
#'   \item{Numeric matrix input}{
#'     \code{x} is used directly after validation. The grouping variable must be
#'     present as a named numeric column. Its observed numeric values are treated
#'     as the user-facing group labels unless an internal formula-method
#'     attribute, \code{"mfpi_group_levels_original"}, is present.
#'   }
#'   \item{Data-frame input}{
#'     \code{x} may contain a categorical \code{group_var}. The grouping
#'     variable may be factor, character, logical, integer, or numeric. It is
#'     converted to numeric input codes before fitting. Its original labels are
#'     retained in \code{group_levels_original} so that print, summary, plot, and
#'     prediction methods can display meaningful group names such as
#'     \code{"Placebo"} and \code{"Treatment"} rather than internal codes.
#'   }
#' }
#'
#' Only \code{group_var} may be categorical in the default method. All other
#' variables must be numeric, integer, or logical. Logical non-group variables
#' are converted to numeric values before the returned matrix is created.
#' Character and factor variables other than \code{group_var} are rejected
#' because \code{mfpi.default()} does not silently create dummy variables for
#' adjustment covariates or continuous candidate variables.
#'
#' The returned matrix is intended for internal modelling only. Downstream
#' functions such as \code{fit_mfpi()}, \code{preprocess_data()}, and
#' \code{fit_mfp()} should therefore continue to receive a numeric matrix even
#' when the user supplied a data frame to \code{mfpi.default()}.
#'
#' The returned group-level metadata separates two concepts:
#'
#' \describe{
#'   \item{\code{group_input_levels}}{
#'     Numeric values present in the returned matrix column for \code{group_var}.
#'     These values are used by \code{preprocess_data()} to map the grouping
#'     variable to consecutive internal levels \code{0, 1, ..., K - 1}.
#'   }
#'   \item{\code{group_levels_original}}{
#'     User-facing labels corresponding to \code{group_input_levels}. These
#'     labels are used for reporting and display.
#'   }
#' }
#'
#' Reference-level behaviour follows the construction of
#' \code{group_input_levels}. For factor input, the reference label is the first
#' element of \code{levels(droplevels(group_var))}. For character input, the
#' reference label is the first distinct value encountered in the data. For
#' logical, integer, numeric, and matrix input, the reference level is the
#' smallest observed value.
#'
#' @param x Predictor input supplied to \code{mfpi.default()}. Must be either a
#'   numeric matrix or a data frame. If \code{x} is a data frame, only
#'   \code{group_var} may be categorical.
#' @param group_var Character scalar naming the grouping variable. Numeric
#'   column indices are not supported.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{x}}{A numeric matrix suitable for internal MFPI fitting.}
#'   \item{\code{group_var}}{The validated grouping-variable name.}
#'   \item{\code{group_input_levels}}{Numeric input levels of the grouping
#'     variable in the returned matrix.}
#'   \item{\code{group_levels_original}}{Character labels corresponding to
#'     \code{group_input_levels}.}
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
