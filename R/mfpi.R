#' Model Interactions Between a Categorical and Continuous Covariates
#'
#' `mfpi()` investigates interactions between a categorical variable
#' (`group_var`) and one or more continuous covariates, using fractional
#' polynomial (FP) transformations to capture nonlinear effects. The categorical
#' variable is treated as a factor internally; in a typical randomised controlled
#' trial it represents treatment allocation, with the lowest value as the control
#' arm. Confounders and other prognostic variables are adjusted for via FP
#' transformations selected by the MFP algorithm implemented in the
#' [mfp2::mfp2()].
#'
#' @section The MFPI approach:
#' The MFPI procedure assesses whether the association between a continuous
#' covariate and the outcome differs across levels of \code{group_var}. For each
#' variable in \code{cont_vars}, the functional form is specified through
#' \code{cont_var_forms} as \code{"linear"}, \code{"fp1"}, or \code{"fp2"}.
#' Fractional polynomial powers are selected according to the requested
#' \code{flex} level, and the interaction is evaluated using the selection
#' criterion specified by \code{criterion}. For methodological background and
#' examples, see the mfpi vignette.
#'
#' @section Flexibility levels (`flex`):
#' The MFPI procedure allows four levels of flexibility in how FP powers are
#' selected and constrained across levels of the treatment variable. These
#' variants are available for FP1 and FP2 functions; linear interaction models
#' are fitted without FP power selection. The appropriate flexibility level is
#' chosen by the analyst based on subject-matter knowledge and sample-size
#' considerations. The functional form tested for each variable in
#' \code{cont_vars} is the one specified in \code{cont_var_forms}; no
#' data-driven selection among linear, FP1, and FP2 candidates is performed.
#'
#' **`flex1` (default; least flexible).** FP powers are selected for the main
#' effect of the continuous covariate \eqn{z} in a model that excludes the
#' treatment-by-covariate interaction. The same selected powers are then used at
#' each level of the treatment variable \eqn{t}. The interaction is assessed by
#' a likelihood ratio test comparing the main effect model with the interaction
#' model. For FP1, the interaction model has 3 degrees of freedom: one power 
#' and two regression coefficients, \eqn{\beta}. For FP2, it has 6
#' degrees of freedom: two powers and four regression coefficients. The
#' corresponding main-effect models have 2 and 4 degrees of freedom,
#' respectively. Therefore, the interaction test has \eqn{3-2 = 1} degree of 
#' freedom for FP1 and \eqn{6-4 = 2} degrees of freedom for FP2.
#'
#' **`flex2`.** FP powers are selected in a model that includes the interaction,
#' with the powers of \eqn{z} constrained to be the same at each level of the
#' treatment variable \eqn{t}. These same powers are then also used for the 
#' main effect model. Because the selected powers may differ from those obtained 
#' under flex1, the interaction test result may also differ. The degrees of 
#' freedom for the interaction test are the same as for flex1: 1 for FP1 and 2 
#' for FP2.
#'
#' **`flex3`.** Separate FP powers are estimated for the main-effects model and
#' for the within-group functions, though the within-group powers are still
#' constrained to be equal across groups. Because the main-effects and
#' interaction models may use different FP families, they are non-nested. The 
#' interaction test has 1 degree of freedom for FP1 and 2 degrees of freedom for
#' FP2.
#'
#' **`flex4` (most flexible).** FP powers are selected for the main effect and
#' separately at each level of the treatment variable. Unlike `flex3`, the FP
#' powers for \eqn{z} are allowed to differ between treatment groups. The
#' interaction test has 2 degrees of freedom for FP1 and 4 degrees of freedom
#' for FP2.
#'
#' Significance tests for interaction are based on a chi-square distribution
#' with the stated degrees of freedom. For flex3 and flex4, some comparisons
#' may be non-nested when different fractional-polynomial powers are estimated
#' for the main effects and interaction models. In those cases, the resulting
#' P-values may lack a formal theoretical underpinning and the significance
#' levels may be liberal or conservative. As an alternative, information
#' criteria such as AIC or BIC may be used to compare the fitted main-effects
#' and interaction models.
#' 
#' @section Prespecified and exploratory interaction analyses:
#' \code{mfpi()} may be used either to test a pre-specified interaction or to
#' screen several candidate interactions. In a pre-specified analysis, one or
#' more interactions have been chosen in advance, for example from prior
#' evidence or a study protocol, and the reported p-values can be interpreted
#' directly without adjustment.
#'
#' In an exploratory analysis, several variables in \code{cont_vars} are tested
#' for interaction with \code{group_var}. In this setting, the reported p-values
#' are not adjusted for multiplicity unless \code{p_adjust_method} is set, for
#' example to \code{"holm"}. Exploratory findings should therefore be treated as
#' candidate interactions requiring further assessment or validation.
#'
#' When multiplicity adjustment is requested, adjustment is performed over the
#' planned family of tests, one for each variable in \code{cont_vars}. Failed or
#' undefined interaction tests remain \code{NA} and cannot be selected, but they
#' are still counted in the adjustment family.
#'
#' In both settings, the functional form for each variable should be specified
#' in advance using \code{cont_var_forms}. Each variable may be assigned
#' \code{"linear"}, \code{"fp1"}, or \code{"fp2"}; variables not named in
#' \code{cont_var_forms} default to \code{"linear"}. Data-driven selection among
#' these forms is not conducted by \code{mfpi()}, because it can inflate Type I
#' error and introduce selection bias.
#'
#' @section Influential observations or Outliers:
#' FP models are sensitive to outliers or influential observations in the
#' continuous predictor. A small number of influential observations in the tails
#' of \code{cont_vars} can disproportionately determine the chosen functional
#' form and may produce \strong{spuriously significant interactions}. It is
#' therefore important to assess whether extreme values may influence the
#' results before fitting the model. Before running \code{mfpi()}, examine the
#' data for influential points using standard diagnostic approaches
#' (e.g. univariable plots).
#'
#' To reduce the potential influence of outliers, \code{mfpi()} provides the
#' option \code{winsorize = TRUE} to apply
#' \strong{Winsorisation} to each variable in \code{cont_vars}. Values below 
#' the lower percentile and above the upper percentile (specified by 
#' \code{winsorize_probs}, with default \code{c(0.01, 0.99)}) are replaced by 
#' the corresponding percentile cutoff values. Winsorisation preserves the 
#' number of observations (i.e. no rows are removed) and only truncates extreme
#' values; the bulk of the distribution remains unaffected.
#'
#' Variables with final \code{zero_vars} handling are \strong{excluded from
#' Winsorisation} so that structural-zero values and the positive distribution
#' are preserved. For variables in \code{cont_vars}, \code{catzero_vars} and
#' \code{spike_vars} are first suppressed; therefore they exclude a variable
#' from Winsorisation only when explicit \code{zero_vars} handling remains.
#'
#' The Winsorisation cutoffs actually applied are returned in the
#' \code{winsorize_limits} component of the fitted object for transparency.
#' To disable Winsorisation entirely, set \code{winsorize = FALSE}. Any trimming or Winsorisation applied
#' should be reported transparently as part of the initial data analysis.
#'
#' @section Regression families (`family`):
#' `mfpi()` accepts family specifications as character strings, GLM family
#' functions, or GLM family objects. Use `family = "gaussian"` or
#' `family = stats::gaussian()` for linear regression, `family = "binomial"`
#' or `family = stats::binomial()` for logistic regression, and
#' `family = "poisson"` or `family = stats::poisson()` for Poisson regression.
#' Custom links may be supplied through family objects, for example
#' `stats::binomial(link = "probit")`. For Cox proportional-hazards models,
#' set `family = "cox"` and supply a `survival::Surv()` response.
#'
#' @section Adjustment variables:
#' Adjustment variables are selected once using \code{mfp2::mfp2()} on
#' \code{x}. The selected variables and their transformations are then used as
#' covariates in the main-effects and interaction models.
#'
#' For each variable tested in \code{cont_vars}, that variable is removed from
#' its own adjustment set because it is already included in the interaction
#' model. The grouping factor, represented by \code{group_var} and any required
#' internal dummy coding, is also excluded from the adjustment set. Therefore,
#' different variables in \code{cont_vars} may be tested with different
#' adjustment sets. If no adjustment variables are selected, interaction tests
#' are fitted without adjustment.
#'
#' To force adjustment variables into the adjustment model, use the corresponding
#' arguments passed to \code{mfp2::mfp2()}, such as \code{select}, \code{alpha},
#' \code{keep}, and \code{df}.
#'
#' @section Shifting, scaling, and centering:
#' Fractional polynomials require strictly positive input values. `mfpi()`
#' estimates a shift for each variable to ensure positivity, then scales the
#' shifted values to a convenient range, following the same procedure as
#' [mfp2::mfp2()]. Centering is applied after FP powers are estimated.
#' Variables with final \code{zero_vars} or \code{catzero_vars} handling, and
#' linear terms (\code{df = 1}), have their shift automatically set to zero
#' because only the positive values are transformed in these cases, or because
#' no nonlinear transformation is applied. For variables in \code{cont_vars},
#' \code{catzero_vars} and \code{spike_vars} are first suppressed; therefore
#' shift is forced to zero only when explicit \code{zero_vars} handling remains
#' or the variable is treated as linear.
#'
#' @section Handling non-positive values:
#' Non-positive values in continuous predictors can be handled using the
#' \code{zero_vars}, \code{catzero_vars}, and \code{spike_vars} options passed
#' to \code{mfp2::mfp2()}.
#'
#' \describe{
#'   \item{\code{zero_vars}}{FP transformations are applied to the positive
#'     values only; non-positive values are set to zero in the transformed
#'     variable.}
#'   \item{\code{catzero_vars}}{As \code{zero_vars}, with an additional
#'     indicator for zero values.}
#'   \item{\code{spike_vars}}{Uses the spike-at-zero algorithm and implies
#'     \code{catzero_vars}.}
#' }
#'
#' Variables listed in \code{catzero_vars} or \code{spike_vars} are supported
#' during stage 1 adjustment selection. However, they cannot also be listed in
#' \code{cont_vars}, because the additional indicator variables created for
#' \code{catzero_vars} and \code{spike_vars} are not supported in the stage 2
#' of the mfpi algorithm.
#'
#' Therefore, if a variable appears in both \code{cont_vars} and
#' \code{spike_vars} or \code{catzero_vars}, the following happens
#' automatically (with a warning):
#' \enumerate{
#'   \item The \code{spike} and \code{catzero} flags are set to \code{FALSE}
#'     for that variable in both the adjustment model and the interaction
#'     stage.
#'   \item Explicit \code{zero_vars} handling is preserved: if the variable is
#'     also listed in \code{zero_vars}, non-positive values are still recoded to
#'     zero before FP transformation. Zero handling implied only by
#'     \code{catzero_vars} or \code{spike_vars} is not retained for variables in
#'     \code{cont_vars}.
#' }
#'
#' Variables in \code{spike_vars} or \code{catzero_vars} that are
#' \strong{not} in \code{cont_vars} are fully respected in the adjustment
#' model as documented in [mfp2::mfp2()].
#' 
#' @section Information criteria:
#' The individual \code{AIC_main}, \code{AIC_interaction},
#' \code{BIC_main}, and \code{BIC_interaction} values reported here are not
#' the actual AIC/BIC values of the fitted models. They are
#' interaction-specific comparison criteria computed from the model deviances
#' and from the degrees of freedom returned by \code{interaction_model_df()},
#' which count only the interaction-related terms and exclude the intercept
#' and adjustment-model parameters.
#'
#' Therefore, these individual AIC/BIC values should not be compared with
#' AIC/BIC values returned by \code{stats::AIC()}, \code{stats::BIC()}, or
#' other model-fitting functions.
#'
#' The differences \code{AIC_main_minus_int} and
#' \code{BIC_main_minus_int} are correct for comparing the main-effects and
#' interaction candidates because the shared model terms cancel.
#'
#' @param x
#'   For \code{mfpi.default()} only. A numeric matrix of dimension
#'   \eqn{n \times p}, where rows are observations and columns are variables.
#'   Column names are required. The matrix must not contain an intercept
#'   column; binary variables should be coded as 0/1; multi-level categorical
#'   variables (other than \code{group_var}) should be expanded into dummy
#'   variables beforehand. Must be free of missing values.
#' @param y
#'   For \code{mfpi.default()} only. The response object. For Gaussian models,
#'   supply a finite numeric vector of length \eqn{n}. For Poisson models,
#'   supply a finite non-negative numeric vector of length \eqn{n}. For binomial
#'   models, supply a numeric vector with values in \eqn{[0, 1]}, a two-level
#'   factor, or a two-column numeric matrix \code{cbind(successes, failures)}.
#'   For Cox models, supply a right-censored \code{survival::Surv()} object.
#'   Must have the same number of observations as \code{x}.
#' @param formula
#'   For \code{mfpi.formula()} only. A formula object describing the model.
#'   Continuous predictors to be FP-transformed should be wrapped in
#'   \code{fp()} or \code{fp2()}, e.g.
#'   \code{Surv(t, d) ~ fp(age, df = 4) + trt + strata(centre)}.
#'   Variables not wrapped in \code{fp()} or \code{fp2()} receive the global scalar defaults
#'   for \code{df}, \code{alpha}, \code{select}, \code{center}, \code{shift},
#'   and \code{scale}. See the \emph{Argument precedence} section in
#'   \code{mfpi.formula()} for full rules.
#' @param data
#'   For \code{mfpi.formula()} only. A data frame containing all variables
#'   named in \code{formula}. Multi-level categorical variables other than
#'   \code{group_var} must already be expanded to dummy columns, or encoded
#'   as factors which \code{model.matrix()} will expand automatically.
#' @param group_var
#'   A single character string naming the categorical (grouping) variable in
#'   \code{x} (e.g. a treatment or exposure indicator). The analysis tests
#'   whether the relationship between each variable in \code{cont_vars} and
#'   the outcome differs across levels of \code{group_var}. The group variable
#'   may use arbitrary numeric levels. Internally, levels are sorted, the lowest
#'   level is treated as the reference group, and dummy variables are created
#'   for all non-reference levels. Must have at least two distinct non-missing
#'   values.
#' @param cont_vars
#'   A non-empty character vector naming the continuous variables in \code{x}
#'   for which interactions with \code{group_var} are to be investigated. All
#'   named variables must exist as columns of \code{x} and must be numeric.
#'   Binary variables (2 or fewer unique values) are not permitted and will
#'   cause an error. Variables with 5 or fewer unique values will trigger a
#'   warning, as FP transformation may be unreliable for near-categorical
#'   variables.
#' @param cont_var_forms
#'   Optional named character vector specifying the functional form to use for
#'   each variable in \code{cont_vars}. Names must be a subset of
#'   \code{cont_vars}; values must each be one of \code{"linear"},
#'   \code{"fp1"}, or \code{"fp2"}. Variables in \code{cont_vars} that are not
#'   named in \code{cont_var_forms} are assigned \code{"linear"} by default.
#'   \code{NULL} (the default) assigns \code{"linear"} to all variables.
#'
#'   Each variable is tested using exactly the pre-specified functional form;
#'   no data-driven selection among forms is performed. This avoids the
#'   inflated Type I error and selection bias that arise when the functional
#'   form is chosen by comparing linear, FP1, and FP2 candidates on the same
#'   data.
#'
#'   Unnamed vectors and scalars are not accepted; every entry must carry the
#'   name of the target variable to ensure the specification is unambiguous.
#'
#'   Example: to test \code{age} as FP2 and \code{bmi} as FP1 while leaving
#'   all other \code{cont_vars} as linear:
#'   \preformatted{cont_var_forms = c(age = "fp2", bmi = "fp1")}
#' @param flex
#'   A character string controlling how FP powers are estimated and constrained
#'   across groups. One of `"flex1"` (default), `"flex2"`, `"flex3"`, or
#'   `"flex4"`. See the *Flexibility levels* section for details.
#' @param p_interact
#'   Numeric in \eqn{(0, 1]}. Nominal significance level for the interaction
#'   test when \code{criterion = "pvalue"}. The interaction model for a variable
#'   is retained only if its p-value is strictly below \code{p_interact}. Default
#'   is \code{0.05}. Ignored when \code{criterion} is \code{"aic"} or
#'   \code{"bic"}. The functional form tested is the one pre-specified in
#'   \code{cont_var_forms}; no selection among forms takes place.
#' @param min_improvement
#'   Numeric. Minimum improvement in the selection criterion required to retain
#'   an interaction term. Interpretation depends on `criterion`:
#'   \itemize{
#'     \item `"pvalue"`: not used directly (the threshold is `p_interact`);
#'       `min_improvement` defaults to `p_interact` for consistency.
#'     \item `"aic"`: minimum required `AIC_main_minus_int`
#'       (\eqn{= \mathrm{AIC}_\text{main} - \mathrm{AIC}_\text{int}}).
#'       Must be positive. Default is `2`.
#'     \item `"bic"`: same as `"aic"` using `BIC_main_minus_int`. Default is `2`.
#'   }
#'   Among candidates that clear `min_improvement`, the one with the
#'   **largest improvement** (largest `AIC_main_minus_int` or `BIC_main_minus_int`)
#'   is chosen as the final functional form. Small differences in AIC or BIC 
#'   provide little evidence for preferring one model over another; values 
#'   around 2 are often used as a rough guide rather than a strict cutoff.
#' @param include_group_var
#'   Logical. Whether `group_var` should be passed to the MFP algorithm
#'   alongside the other adjustment variables when selecting the adjustment
#'   model in stage one. Default is `FALSE`. When `TRUE`, the nominal
#'   significance level for `group_var` is set to 1, forcing it into the
#'   adjustment model.
#' @param show_models
#'   Logical. If \code{TRUE} and \code{verbose = TRUE}, prints the full
#'   regression coefficient table (via \code{summary()}) for the winning
#'   interaction model of each variable in \code{cont_vars} that cleared the
#'   selection threshold. Has no effect when \code{verbose = FALSE}. When no
#'   variable is selected nothing is printed. Default \code{FALSE}.
#' @param weights
#'   An optional numeric vector of non-negative observation weights of length
#'   \eqn{n}. Applied to both the adjustment model (MFP step) and the
#'   interaction model fitting step. Default \code{NULL} (all weights equal
#'   to 1).
#' @param offset
#'   An optional numeric vector of length \eqn{n} to be added to the linear
#'   predictor. Applied to both the adjustment model (MFP step) and the
#'   interaction model fitting step. Useful for Poisson models (e.g. log of
#'   exposure time). Default \code{NULL} (zero offset for all observations).
#' @param cycles
#'   A positive integer. Maximum number of iteration cycles for the MFP
#'   algorithm. Default is `10`.
#' @param shift
#'   A numeric vector of length \eqn{p} or a single numeric giving shift terms
#'   added to columns of \code{x} before scaling. Default is \code{NULL}, which lets
#'   the program estimate shifting factors automatically via
#'   [mfp2::find_shift_factor()]. Set \code{shift = 0} to disable
#'   shifting. For \code{zero_vars}, \code{catzero_vars}, and linear variables
#'   (\code{df = 1}), shift is forced to 0.
#' @param scale
#'   A numeric vector of length \eqn{p} or a single numeric giving scaling
#'   factors for the columns of \code{x}. Default is \code{NULL}, which lets
#'   the program estimate scaling factors automatically via
#'   [mfp2::find_scale_factor()]. Set \code{scale = 1} to disable scaling.
#'
#'   Internally, \code{x} is divided by \code{scale} before model fitting for
#'   numerical stability during FP power selection. The scale factors are then
#'   passed through the call chain so that both the adjustment model and the
#'   interaction model backscale before their final fit. This ensures all
#'   returned model coefficients are on the \eqn{\phi(x + \text{shift})} scale.
#' @param df
#'   A numeric vector of length \eqn{p} or a single positive integer setting
#'   the default degrees of freedom for each predictor. Degrees of freedom
#'   equal twice the FP degree (e.g. `df = 1` for linear, `df = 2` for FP1, 
#'   `df = 4` for FP2). Regardless of the supplied value, the program overrides
#'   `df` per variable based on the number of distinct values \eqn{u}:
#'   \itemize{
#'     \item \eqn{u \le 3}: `df` is set to `1` (linear). Too few values to
#'       estimate a curve.
#'     \item \eqn{4 \le u \le 5}: `df` is capped at `min(2, df)` (FP1 at
#'       most). Enough variation for a simple curve but not for FP2.
#'     \item \eqn{u \ge 6}: `df` is retained as supplied.
#'   }
#'   Default is `4`.
#' @param center
#'   Logical or a logical vector of length \eqn{p}. Whether to mean-centre
#'   predictors before fitting the final interaction model and adjustment model.
#'   Binary covariates are centered at the lower of their two values rather than
#'   the mean. Default is `TRUE`.
#' @param subset
#'   Optional subset of observations used for model selection and fitting.
#'   May be either a logical vector with one value per observation or a numeric
#'   vector of positive integer row indices. Missing, non-finite, zero,
#'   negative, out-of-range, non-integer, and duplicated indices are not allowed.
#'   Default is `NULL`, meaning that all observations are used. As in
#'   \code{\link[mfp2]{mfp2}}, subsetting is applied after preprocessing
#'   parameters such as shift, scale, and centering constants have been computed
#'   from the full data.
#' @param family
#'   Model family specification. May be one of the character strings
#'   `"gaussian"`, `"binomial"`, `"poisson"`, or `"cox"`; a GLM family function
#'   such as `stats::binomial`; or a GLM family object such as
#'   `stats::binomial(link = "probit")`. Cox models must be specified as
#'   `family = "cox"`.
#' @param criterion
#'   A character string specifying the criterion used in two distinct places:
#'   \enumerate{
#'     \item **Adjustment-variable selection** (Step 1): governs which
#'       predictors and FP degrees survive MFP backfitting.
#'     \item **Interaction test** (Step 2): for each variable in `cont_vars`,
#'       the pre-specified functional form (see `cont_var_forms`) is tested
#'       against its main-effects counterpart using this criterion.
#'   }
#'   One of `"pvalue"` (default), `"aic"`, or `"bic"`.
#' @param select
#'   A numeric vector of length \eqn{p} or a single value in \eqn{[0, 1]}
#'   giving the nominal significance level used during MFP backfitting to
#'   decide whether each predictor is retained in the adjustment model.
#'   At each cycle, a variable is kept only if its contribution is significant
#'   at this level; otherwise it is set to zero (excluded). Setting a
#'   variable's level to \code{1} forces it into the model. Default is 
#'   \code{0.05}.
#' @param alpha
#'   A numeric vector of length \eqn{p} or a single value in \eqn{[0, 1]}
#'   giving the significance level for choosing between FP degrees for each
#'   predictor. Default is `0.05`.
#' @param keep
#'   An optional character vector of variable names to retain in the adjustment
#'   model regardless of selection criteria. When `criterion = "pvalue"`,
#'   equivalent to setting `select = 1` for those variables; also effective
#'   under AIC and BIC criteria.
#' @param force_max_fp_vars
#'   An optional character vector naming variables for which the algorithm
#'   should select the most complex FP functional form at the degree specified
#'   by \code{df}, preventing simplification to a lower degree. The best power
#'   combination within that degree is still selected by the criterion
#'   (equivalently, by deviance minimisation at fixed df).
#'
#'   The mechanism depends on the selection criterion:
#'   \itemize{
#'     \item \code{criterion = "pvalue"}: \code{alpha} is automatically set
#'       to 1 for the named variables, so the significance test always
#'       accepts the most complex form.
#'     \item \code{criterion = "aic"} or \code{"bic"}:
#'       \code{select_ic()} normally competes null, linear, FP1, and
#'       FP2 against each other and may simplify the functional form;
#'       listing a variable here suppresses this simplification.
#'   }
#'   Default \code{NULL} (no variables forced). In the formula interface,
#'   per-variable control is also available via
#'   \code{fp(force_max_fp = TRUE)} or \code{fp2(force_max_fp = TRUE)}.
#' @param xorder
#'   A character string controlling the order in which adjustment covariates
#'   enter the MFP selection algorithm. `"ascending"` (default) enters
#'   variables from most to least significant in a full multiple regression;
#'   `"descending"` reverses this order; `"original"` uses the column order of
#'   `x`.
#' @param powers
#'   A named list of numeric vectors giving the candidate FP powers for each
#'   variable. Default is `NULL`, which uses the standard set proposed by
#'   Royston and Altman (1994): \eqn{\{-2, -1, -0.5, 0, 0.5, 1, 2, 3\}}, where
#'   0 denotes the natural logarithm. Each element must be named after the
#'   corresponding column of `x` and must contain at least two distinct values.
#' @param ties
#'   A character string specifying the method for handling tied event times in
#'   Cox regression. One of `"breslow"` (default), `"efron"`, or `"exact"`.
#'   Ignored for non-Cox families. See [survival::coxph()] for details.
#' @param strata
#'   A numeric vector or matrix defining stratification factors for Cox models.
#'   A single combined factor is created from all supplied variables. Default is
#'   `NULL` (no stratification). Currently only a single stratification factor
#'   is supported.
#' @param nocenter
#'   A numeric vector of values passed to [survival::coxph()]. Cox models only;
#'   ignored otherwise.
#' @param acd_vars
#'   An optional character vector naming continuous variables to be transformed
#'   via the approximate cumulative distribution (ACD) transformation before FP
#'   selection. The transformed variable is named `A(x)`. ACD transformation is
#'   disabled for variables listed in `cont_vars`. See [mfp2::mfp2()] for
#'   details.
#' @param zero_vars
#'   An optional character vector naming variables for which non-positive values
#'   should be treated as structural zero. FP transformations are applied only
#'   to strictly positive values; non-positive values are represented as zero in
#'   the transformed columns. This treatment is compatible with
#'   \code{cont_vars} and is preserved for interaction variables when explicitly
#'   requested through \code{zero_vars}. See the \emph{Handling non-positive
#'   values} section for details.
#' @param catzero_vars
#'   An optional character vector naming variables for which a binary structural
#'   zero indicator \eqn{I(x \le 0)} should be added to the adjustment model
#'   alongside the FP terms for the positive part. For adjustment variables,
#'   \code{catzero_vars} implies \code{zero_vars}. A variable that is not in
#'   \code{cont_vars} may not appear in both \code{zero_vars} and
#'   \code{catzero_vars}. \strong{Note:} if a variable also appears in
#'   \code{cont_vars}, the structural-zero indicator is dropped and
#'   \code{catzero} is set to \code{FALSE} for that variable in both the
#'   adjustment model and the interaction stage, with a warning. Explicit
#'   \code{zero_vars} handling is preserved if the variable was also listed in
#'   \code{zero_vars}; zero handling implied only by \code{catzero_vars} is not
#'   retained for variables in \code{cont_vars}.
#' @param spike_vars
#'   An optional character vector naming variables to be assessed for a spike at
#'   zero using the SAZ algorithm. For adjustment variables, \code{spike_vars}
#'   implies \code{catzero_vars} and therefore \code{zero_vars}. \strong{Note:}
#'   if a variable also appears in \code{cont_vars}, spike-at-zero handling is
#'   dropped and both \code{spike} and the implied \code{catzero} indicator are
#'   set to \code{FALSE} for that variable in both the adjustment model and the
#'   interaction stage, with a warning. Explicit \code{zero_vars} handling is
#'   preserved if the variable was also listed in \code{zero_vars}; zero handling
#'   implied only by \code{spike_vars} is not retained for variables in
#'   \code{cont_vars}.
#' @param min_prop
#'   Numeric in \eqn{(0, 1)}. Minimum proportion of zeros required for the
#'   SAZ algorithm to be applied to a variable in \code{spike_vars}. Default
#'   \code{0.05}. Must be less than \code{max_prop}. If the observed zero
#'   proportion is below \code{min_prop} (too few zeros to model a spike
#'   meaningfully), the spike flag is reset to \code{FALSE} for that variable.
#'   The resulting treatment depends on what the user originally specified
#'   alongside \code{spike_vars}: if \code{catzero_vars} or \code{zero_vars}
#'   were also specified for that variable, those flags are preserved;
#'   otherwise the variable reverts to a standard continuous predictor.
#'   Only affects variables in the adjustment model since \code{spike_vars}
#'   is always suppressed for \code{cont_vars}.
#' @param max_prop
#'   Numeric in \eqn{(0, 1)}. Maximum proportion of zeros allowed for the
#'   SAZ algorithm to be applied to a variable in \code{spike_vars}. Default
#'   \code{0.95}. Must be greater than \code{min_prop}. If the observed zero
#'   proportion exceeds \code{max_prop} (too many zeros; the positive part is
#'   too sparse for reliable FP fitting), the spike flag is reset to
#'   \code{FALSE}. As with \code{min_prop}, the resulting treatment of
#'   \code{catzero} and \code{zero} depends on the user's original
#'   specification. Variables that are binary (exactly two unique values) are
#'   also reset regardless of their zero proportion. Only affects variables
#'   in the adjustment model.
#' @param ftest
#'   Logical. Whether to use an F-test rather than a chi-square
#'   likelihood-ratio test when computing p-values for Gaussian models.
#'   When \code{TRUE} and \code{family = "gaussian"}, the F-test is applied
#'   to both the adjustment-variable selection step and the interaction test 
#'   for each variable in \code{cont_vars}. Recommended when the sample size is
#'   small. Default \code{FALSE}. Has no effect for non-Gaussian families or 
#'   when \code{criterion} is not \code{"pvalue"}. Inactive
#' @param control
#'   A list of control parameters for the underlying fitting routine, as
#'   returned by [stats::glm.control()] (non-Cox families) or
#'   [survival::coxph.control()] (Cox). Default is `NULL`, which uses the
#'   default control parameters for the chosen family.
#' @param winsorize
#'   Logical. Whether to Winsorise the continuous variables in `cont_vars`
#'   before fitting, to reduce the influence of extreme observations on the
#'   selected FP functional form. Default is `FALSE`. See the
#'   \emph{Influential observations} section.
#' @param winsorize_probs
#'   Numeric vector of length 2 giving the lower and upper percentile
#'   cutoffs used for Winsorisation. Default is `c(0.01, 0.99)`, which
#'   truncates approximately the bottom and top 1\% of values for each
#'   `cont_var`. Ignored when `winsorize = FALSE`.
#' @param center_type
#'   Character string controlling how the FP-transformed continuous variables
#'   are centered in the \strong{interaction model} when \code{center = TRUE}.
#'   Has no effect on the adjustment model, which uses standard per-variable
#'   centering via \code{fit_mfp()}. One of:
#'   \describe{
#'     \item{\code{"grand"} (default)}{For each variable in \code{cont_vars},
#'       center each column of the group-specific block matrix by the grand
#'       mean of the FP-transformed variable across all observations.}
#'     \item{\code{"group"}}{Center each group's columns by the within-group
#'       mean, the mean of in-group observations only. Removes
#'       group-specific location effects from the FP-transformed variable
#'       before fitting the interaction model.}
#'   }
#'   The centering constants computed at fit time are stored in
#'   \code{center_vals_list} on the returned object and reused exactly
#'   during evaluation of fitted functions.
#' @param p_adjust_method
#'   Character string specifying the method for adjusting p-values when
#'   \code{criterion = "pvalue"}. Passed to \code{\link[stats]{p.adjust}}.
#'   Accepted values include \code{"none"} (default, no adjustment),
#'   \code{"holm"}, \code{"bonferroni"}, \code{"hochberg"}, \code{"BH"},
#'   and \code{"BY"}. P-values are adjusted over the planned family of
#'   interaction tests, namely one test for each variable in \code{cont_vars}.
#'   Variables for which no valid p-value can be computed remain \code{NA}
#'   and are not selected, but they are still counted in the multiplicity
#'   adjustment. The argument has no effect on AIC/BIC-based selection
#'   decisions.
#' @param verbose
#'   Logical. Whether to print progress information during model fitting.
#'   Default is `TRUE`.
#' @param digits
#'   A positive integer. Minimum number of significant digits displayed when
#'   printing interaction-test results. Default is `3`.
#' @param ...
#'   Currently unused. Reserved for future extensions.
#'
#' @return
#' An object of class \code{"mfpi"}. The object is a list with the following
#' components:
#'
#' \describe{
#'   \item{\code{best_model_metrics}}{A data frame of evaluation metrics for
#'     the retained interaction models. Each row corresponds to one selected
#'     continuous variable. The column \code{type} gives the tested functional
#'     form, one of \code{"linear"}, \code{"fp1"}, or \code{"fp2"}. The
#'     metrics include the main and interaction FP powers and the
#'     criterion-specific interaction evaluation columns. For
#'     \code{criterion = "pvalue"}, the metrics include \code{pvalue} and,
#'     when applicable, \code{p_adjusted}. For \code{criterion = "aic"} or
#'     \code{criterion = "bic"}, p-value columns are omitted and the reported
#'     interaction metrics are based on AIC/BIC only. If no interaction is
#'     retained, this is an empty data frame. The reported AIC/BIC columns are
#'      interaction-specific comparison criteria, not full fitted-model AIC/BIC 
#'      values; see \code{test_interaction()} for details.}
#'   \item{\code{all_model_metrics}}{A data frame with the same columns as
#'     \code{best_model_metrics} but containing metrics for the interaction
#'     model evaluated for each variable in \code{cont_vars}, using the form
#'     specified by \code{cont_var_forms}.}
#'   \item{\code{best_interaction_model}}{A named list of fitted interaction
#'     model objects for the continuous variables whose interactions were
#'     retained by the selection criterion. Each element is named by the
#'     corresponding variable in \code{cont_vars} and contains the fitted
#'     interaction-model object returned by the interaction test, typically
#'     with the underlying regression fit stored in its \code{fit} component.
#'     If no interaction is retained, this is an empty list.}
#'   \item{\code{all_interaction_models}}{A named list with one element per
#'     \code{cont_var}. Each element is itself a named list containing the
#'     fitted interaction model for the functional form specified in
#'     \code{cont_var_forms} for that variable, for example
#'     \code{object$all_interaction_models$age$fp2}.}
#'   \item{\code{center_vals_list}}{A named list (one element per
#'     \code{cont_var}) of centering constants used when fitting the
#'     interaction model. These are the exact values subtracted from the
#'     FP-transformed variables during model fitting and must be reused
#'     (not recomputed) during prediction. \code{NULL} when
#'     \code{center = FALSE}.}
#'   \item{\code{adjust_terms}}{A data frame describing the FP terms
#'     selected for the adjustment model.}
#'   \item{\code{adjustment_model}}{The full adjustment model object
#'     returned by \code{mfp2()}.}
#'   \item{\code{univariable_interactions}}{The full return value of the
#'     internal \code{evaluate_interactions()} call. Contains all
#'     intermediate results and is useful for programmatic access.}
#'   \item{\code{group_var}}{The name of the grouping variable.}
#'   \item{\code{group_levels_new}}{A zero-based integer index
#'     (0, 1, 2, \ldots) corresponding to the sorted levels of
#'     \code{group_var}. This records the level ordering; the observed group
#'     values themselves are not replaced.}
#'   \item{\code{group_levels_original}}{The sorted original levels of
#'     \code{group_var} as they appear in the data. The lowest level is the
#'     reference group.}
#'   \item{\code{flex}}{The flexibility level used (\code{"flex1"} through
#'     \code{"flex4"}).}
#'   \item{\code{family}}{The regression family used for fitting: a GLM
#'     family object for Gaussian, binomial, and Poisson models, or the
#'     character string \code{"cox"} for Cox models.}
#'   \item{\code{nobs}}{Number of observations used in model fitting.}
#'   \item{\code{show_models}}{Value of the \code{show_models} argument.}
#'   \item{\code{winsorize}}{Logical; whether Winsorisation was applied.}
#'   \item{\code{winsorize_probs}}{The probability cutoffs used for
#'     Winsorisation, or \code{NULL} when \code{winsorize = FALSE}.}
#'   \item{\code{winsorize_limits}}{A named list of the actual lower and
#'     upper limits applied to each \code{cont_var} during Winsorisation,
#'     or \code{NULL} when \code{winsorize = FALSE}.}
#'   \item{\code{p_adjust_method}}{Character string: the multiplicity
#'     adjustment method used (e.g. \code{"none"}, \code{"holm"}).}
#'   \item{\code{p_interact}}{Numeric: the significance threshold for
#'     interaction selection.}
#'   \item{\code{min_improvement}}{Numeric: the threshold for AIC/BIC
#'     selection. Only relevant when \code{criterion = "aic"} or
#'     \code{"bic"}.}
#'   \item{\code{var_winners}}{A named list (one element per \code{cont_var})
#'     storing the best candidate for each variable regardless of whether
#'     it was selected. Each element contains \code{fit} (the fitted model),
#'     \code{metric} (evaluation metrics), \code{type} (winning functional
#'     form), \code{score} (the decisive metric value), and
#'     \code{center_vals} (centering constants).}
#'   \item{\code{criterion}}{character specifying criterion used for adjustment
#'   and interaction models}
#' }
#'
#' @references
#' Royston, P. and Sauerbrei, W. (2004). A new approach to modelling
#' interactions between treatment and continuous covariates in clinical trials
#' by using fractional polynomials. \emph{Statistics in Medicine}, 23,
#' 2509--2525.
#'
#' Royston, P. and Sauerbrei, W. (2008). Interactions between treatment and
#' continuous covariates -- a step towards individualising therapy (Editorial).
#' \emph{Journal of Clinical Oncology}, 26, 1397--1399.
#'
#' Royston, P. and Sauerbrei, W. (2013). Interaction of treatment with a
#' continuous variable: simulation study of significance level for several
#' methods of analysis. \emph{Statistics in Medicine}, 32, 3788--3803.
#'
#' Royston, P. and Sauerbrei, W. (2014). Interaction of treatment with a
#' continuous variable: simulation study of power for several methods of
#' analysis. \emph{Statistics in Medicine}, 33, 4695--4708.
#'
#' Royston, P., Sauerbrei, W. and Ritchie, A. (2004). Is treatment with
#' interferon-alpha effective in all patients with metastatic renal carcinoma?
#' A new approach to the investigation of interactions. \emph{British Journal
#' of Cancer}, 90, 794--799.
#'
#' Sauerbrei, W., Royston, P. and Zapien, K. (2007). Detecting an interaction
#' between treatment and a continuous covariate: a comparison of two
#' approaches. \emph{Computational Statistics and Data Analysis}, 51,
#' 4054--4063.
#'
#' Royston, P. and Sauerbrei, W. (2008). \emph{Multivariable Model-Building: A
#' Pragmatic Approach to Regression Analysis Based on Fractional Polynomials
#' for Modelling Continuous Variables}. John Wiley & Sons.
#'
#' @seealso [mfp2::summary.mfpi()], [mfp2::mfp2()], [mfp2::print.mfpi()] and 
#' [mfp2::plot.mfpi()]
#' 
#' @examples
#' data("prostate")
#' # Flexibility level 1: simplest interaction structure
#' flex_1 <- mfpi(
#'   lpsa ~ fp(age) + svi + fp(pgg45) + fp(cavol) + fp(weight) +
#'     fp(bph) + fp(cp),
#'   data = prostate,
#'   cont_vars = c("cavol", "age"),
#'   group_var = "svi",
#'   center = FALSE,
#'   flex = "flex1",
#'   include_group_var = TRUE,
#'   show_models = TRUE,
#'   verbose = TRUE
#' )
#'
#' # Plot both treatment-specific curves and treatment effects
#' plot(flex_1, terms = "cavol", plot_type = "both")
#'
#' # Pre-specify functional forms: cavol as FP2, age as linear (default)
#' flex_1_prespec <- mfpi(
#'   lpsa ~ fp(age) + svi + fp(pgg45) + fp(cavol) + fp(weight) +
#'     fp(bph) + fp(cp),
#'   data = prostate,
#'   cont_vars      = c("cavol", "age"),
#'   cont_var_forms = c(cavol = "fp2"),
#'   group_var      = "svi",
#'   center         = FALSE,
#'   flex           = "flex1",
#'   include_group_var = TRUE,
#'   verbose        = TRUE
#' )
#'
#' # Flexibility level 2:
#' flex_2 <- mfpi(
#'   lpsa ~ fp(age) + svi + fp(pgg45) + fp(cavol) + fp(weight) +
#'     fp(bph) + fp(cp),
#'   data = prostate,
#'   cont_vars = c("cavol", "age"),
#'   group_var = "svi",
#'   center = FALSE,
#'   flex = "flex2",
#'   include_group_var = TRUE,
#'   show_models = TRUE,
#'   verbose = TRUE
#' )
#'
#' plot(flex_2, terms = "cavol", plot_type = "both")
#'
#' # Flexibility level 3: 
#' flex_3 <- mfpi(
#'   lpsa ~ fp(age) + svi + fp(pgg45) + fp(cavol) + fp(weight) +
#'     fp(bph) + fp(cp),
#'   data = prostate,
#'   cont_vars = c("cavol", "age"),
#'   group_var = "svi",
#'   center = FALSE,
#'   flex = "flex3",
#'   include_group_var = TRUE,
#'   show_models = TRUE,
#'   verbose = TRUE
#' )
#'
#' plot(flex_3, terms = "cavol", plot_type = "both")
#'
#' # Flexibility level 4: most flexible specification 
#' flex_4 <- mfpi(
#'   lpsa ~ fp(age) + svi + fp(pgg45) + fp(cavol) + fp(weight) +
#'     fp(bph) + fp(cp),
#'   data = prostate,
#'   cont_vars = c("cavol", "age"),
#'   group_var = "svi",
#'   center = FALSE,
#'   flex = "flex4",
#'   include_group_var = TRUE,
#'   show_models = TRUE,
#'   verbose = TRUE
#' )
#'
#' plot(flex_4, terms = "cavol", plot_type = "both")
#' @export
mfpi <- function(x, ...) {
  UseMethod("mfpi", x)
}

#' @rdname mfpi
#'
#' @details
#' The default method accepts a numeric matrix \code{x} and response vector
#' \code{y}.
#'
#' @export

mfpi.default <- function(
    x,
    y,
    group_var         = NULL,
    cont_vars         = NULL,
    cont_var_forms    = NULL,
    flex              = c("flex1", "flex2", "flex3", "flex4"),
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
    ties              = c("breslow", "efron", "exact"),
    strata            = NULL,
    nocenter          = NULL,
    acd_vars          = NULL,
    zero_vars         = NULL,
    catzero_vars      = NULL,
    spike_vars        = NULL,
    min_prop          = 0.05,
    max_prop          = 0.95,
    ftest             = FALSE,
    control           = NULL,
    winsorize         = FALSE,
    winsorize_probs   = c(0.01, 0.99),
    center_type       = c("grand", "group"),
    p_adjust_method   = "none",
    verbose           = TRUE,
    digits            = 3,
    ...
) {
  cl <- match.call()
  
  # Match enumerated arguments -------------------------------------------------
  xorder      <- match.arg(xorder)
  criterion   <- match.arg(criterion)
  flex        <- match.arg(flex)
  ties        <- match.arg(ties)
  center_type <- match.arg(center_type)
  
  # Resolve family -------------------------------------------------------------
  family_info <- normalize_family_argument(family)
  
  family        <- family_info$family
  family_string <- family_info$family_string
  
  if (!is.matrix(x)) {
    stop("! `x` must be a matrix.", call. = FALSE)
  }
  
  vnames <- colnames(x)
  if (is.null(vnames)) {
    stop(
      "! `x` must have column names.\n",
      "i Use `colnames(x) <- ...` to assign them.",
      call. = FALSE
    )
  }
  
  if (any(is.character(x))) {
    stop(
      "! `x` contains character values.\n",
      "i Convert all categorical variables to numeric dummy columns before calling `mfpi()`.",
      call. = FALSE
    )
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
      "! `x` must be a numeric matrix with finite values.",
      call. = FALSE
    )
  }
  
  # Basic dimension checks on x ------------------------------------------------
  np    <- dim(x)
  nobs  <- as.integer(np[1L])
  nvars <- as.integer(np[2L])
  
  if (is.null(np)) {
    stop(
      "! `x` must be a matrix with at least one row and one column.\n",
      "i `dim(x)` returned NULL.",
      call. = FALSE
    )
  }
  
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
  
  # Validate and expand cont_var_forms ------------------------------------------
  # cont_var_forms must be either NULL (default: all linear) or a named
  # character vector whose names are a subset of cont_vars and whose values
  # are each one of "linear", "fp1", "fp2". Missing cont_vars entries are
  # silently filled with "linear". Unnamed vectors and unnamed scalars are
  # rejected to enforce explicit per-variable specification.
  valid_forms <- c("linear", "fp1", "fp2")
  
  if (is.null(cont_var_forms)) {
    cont_var_forms <- stats::setNames(
      rep("linear", length(cont_vars)), cont_vars
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
    # Fill any cont_vars not mentioned with "linear"
    missing_vars <- setdiff(cont_vars, names(cont_var_forms))
    if (length(missing_vars) > 0L) {
      cont_var_forms <- c(
        cont_var_forms,
        stats::setNames(rep("linear", length(missing_vars)), missing_vars)
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
  n_unique <- vapply(cont_vars, function(v)
    length(unique(x[!is.na(x[, v, drop = TRUE]), v, drop = TRUE])), integer(1L))
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
  
  # Validate response y --------------------------------------------------------
  validate_family_response(
    y             = y,
    family_string = family_string,
    nobs          = nobs
  )
  
  # Validate Cox-specific auxiliary inputs ------------------------------------
  if (family_string == "cox" && !is.null(strata)) {
    if (is.factor(strata)) strata <- as.numeric(strata)
    
    strata_len <- if (is.vector(strata)) length(strata) else NROW(strata)
    
    if (strata_len != nobs) {
      stop(
        paste0(
          "! Length of `strata` (", strata_len, ") must equal the number of ",
          "rows in `x` (", nobs, ")."
        ),
        call. = FALSE
      )
    }
  }
  
  # Validate weights and offset ------------------------------------------------
  if (!is.null(weights)) {
    if (
      !is.numeric(weights) ||
      length(weights) != nobs ||
      anyNA(weights) ||
      any(!is.finite(weights)) ||
      any(weights < 0)
    ) {
      stop(
        paste0(
          "! `weights` must be a finite non-negative numeric vector ",
          "of length nobs."
        ),
        call. = FALSE
      )
    }
  }
  
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
  
  # Validate alpha and select --------------------------------------------------
  if (!is.numeric(alpha) || anyNA(alpha) || any(!is.finite(alpha))) {
    stop("! `alpha` must be numeric, finite, and non-missing.", call. = FALSE)
  }
  
  if (!is.numeric(select) || anyNA(select) || any(!is.finite(select))) {
    stop("! `select` must be numeric, finite, and non-missing.", call. = FALSE)
  }
  
  if (!is.numeric(df) || anyNA(df) || any(!is.finite(df))) {
    stop("! `df` must be numeric, finite, and non-missing.", call. = FALSE)
  }
  
  if (any(alpha < 0) || any(alpha > 1)) {
    stop("! All values of `alpha` must be in [0, 1].", call. = FALSE)
  }
  
  if (length(alpha) != 1L && length(alpha) != nvars) {
    stop(paste0("! `alpha` must be a single number or a vector of length ", nvars,
                "; got length ", length(alpha), "."), call. = FALSE)
  }
  
  if (any(select < 0) || any(select > 1)) {
    stop("! All values of `select` must be in [0, 1].", call. = FALSE)
  }
  
  if (length(select) != 1L && length(select) != nvars) {
    stop(paste0("! `select` must be a single number or a vector of length ", nvars,
                "; got length ", length(select), "."), call. = FALSE)
  }
  
  # Validate keep --------------------------------------------------------
  if (!is.null(keep) && !all(keep %in% vnames)) {
    warning(
      "i Some variables in `keep` are not columns of `x`; ",
      "continuing with the intersection.",
      call. = FALSE
    )
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
    bad_names <- setdiff(force_max_fp_vars, vnames)
    if (length(bad_names) > 0L) {
      warning("i Some variables in `force_max_fp_vars` are not columns of `x`; ",
              "they will be ignored: ", paste(bad_names, collapse = ", "), ".",
              call. = FALSE)
    }
    force_max_fp_vars <- intersect(force_max_fp_vars, vnames)
    if (length(force_max_fp_vars) > 0L) {
      force_max_fp[force_max_fp_vars] <- TRUE
    }
  }
  
  # Validate zero_vars / catzero_vars / spike_vars -----------------------------
  if (!is.null(zero_vars) && !all(zero_vars %in% vnames)) {
    warning("i Some variables in `zero_vars` are not columns of `x`; they will be ignored.",
            call. = FALSE)
  }
  
  if (!is.null(catzero_vars) && !all(catzero_vars %in% vnames)) {
    warning("i Some variables in `catzero_vars` are not columns of `x`; they will be ignored.",
            call. = FALSE)
  }
  
  if (!is.null(spike_vars) && !all(spike_vars %in% vnames)) {
    warning("i Some variables in `spike_vars` are not columns of `x`; they will be ignored.",
            call. = FALSE)
  }
  
  zero_vars_in    <- intersect(zero_vars,    vnames)
  catzero_vars_in <- intersect(catzero_vars, vnames)
  spike_vars_in   <- intersect(spike_vars,   vnames)
  
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
  if (!is.null(acd_vars) && !all(acd_vars %in% vnames)) {
    warning("i Some variables in `acd_vars` are not columns of `x`; they will be ignored.",
            call. = FALSE)
  }
  
  # Validate shift and scale ---------------------------------------------------
  if (!is.null(shift)) {
    if (!is.numeric(shift)) {
      stop("! `shift` must be numeric or NULL.", call. = FALSE)
    }
    
    bad_shift <- !is.na(shift) & !is.finite(shift)
    
    if (any(bad_shift)) {
      stop(
        "! All non-missing values of `shift` must be finite.",
        call. = FALSE
      )
    }
  }
  
  if (!is.null(shift) && length(shift) != 1L && length(shift) != nvars) {
    stop(
      paste0(
        "! `shift` must be a single number or a vector of length ",
        nvars, "; got length ", length(shift), "."
      ),
      call. = FALSE
    )
  }
  
  if (!is.null(scale)) {
    if (!is.numeric(scale)) {
      stop("! `scale` must be numeric or NULL.", call. = FALSE)
    }
    
    bad_scale <- !is.na(scale) & (!is.finite(scale) | scale <= 0)
    
    if (any(bad_scale)) {
      stop(
        "! All non-missing values of `scale` must be finite and positive.",
        call. = FALSE
      )
    }
  }
  
  if (!is.null(scale) && length(scale) != 1L && length(scale) != nvars) {
    stop(
      paste0(
        "! `scale` must be a single number or a vector of length ",
        nvars, "; got length ", length(scale), "."
      ),
      call. = FALSE
    )
  }
  
  # Validate center ------------------------------------------------------------
  # Developer note: mfpi.formula() may pass a per-variable logical vector after
  # reading fp()/fp2() attributes, so the default method is the final validation
  # gate for both direct matrix calls and formula calls.
  if (!is.logical(center) || anyNA(center)) {
    stop(
      "! `center` must be logical and must not contain missing values.",
      call. = FALSE
    )
  }
  
  if (length(center) != 1L && length(center) != nvars) {
    stop(
      paste0(
        "! `center` must be a single logical or a logical vector of length ",
        nvars, "; got length ", length(center), "."
      ),
      call. = FALSE
    )
  }
  
  # Validate df ----------------------------------------------------------------
  if (any(df <= 0L)) {
    stop("! All values of `df` must be positive (1 for linear, 2m for FP degree m).",
         call. = FALSE)
  }
  
  if (length(df) == 1L) {
    if (df != 1L && df %% 2L != 0L)
      stop(paste0("! `df = ", df, "` is invalid. `df` must be 1 (linear) or an even ",
                  "number 2m for FP degree m."), call. = FALSE)
  } else {
    if (length(df) != nvars)
      stop(paste0("! When `df` is a vector it must have length ", nvars,
                  "; got length ", length(df), "."), call. = FALSE)
    invalid_df <- df != 1L & df %% 2L != 0L
    if (any(invalid_df))
      stop(paste0("! Each element of `df` must be 1 (linear) or an even number 2m. ",
                  "Invalid values at positions: ",
                  paste(which(invalid_df), collapse = ", "), "."), call. = FALSE)
  }
  
  # Validate spike proportions -------------------------------------------------
  if (!is.numeric(min_prop) || length(min_prop) != 1L ||
      anyNA(min_prop) || !is.finite(min_prop) ||
      min_prop < 0 || min_prop > 1) {
    stop(
      "! `min_prop` must be a single finite numeric value in [0, 1].",
      call. = FALSE
    )
  }
  
  if (!is.numeric(max_prop) || length(max_prop) != 1L ||
      anyNA(max_prop) || !is.finite(max_prop) ||
      max_prop < 0 || max_prop > 1) {
    stop(
      "! `max_prop` must be a single finite numeric value in [0, 1].",
      call. = FALSE
    )
  }
  
  if (min_prop > max_prop) {
    stop("! `min_prop` cannot be greater than `max_prop`.", call. = FALSE)
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
  
  # Validate and build powers list ------------------------------------------
  if (!is.null(powers) && !is.list(powers)) {
    stop("! `powers` must be a named list.", call. = FALSE)
  }
  
  default_powers <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
  power_list <- setNames(replicate(nvars, default_powers, simplify = FALSE), vnames)
  
  if (!is.null(powers)) {
    if (length(powers) != sum(nchar(names(powers)) > 0L, na.rm = TRUE)) {
      stop("! Every element of `powers` must have a name.", call. = FALSE)
    }
    unknown_names <- setdiff(names(powers), vnames)
    if (length(unknown_names) > 0L) {
      stop(paste0("! The following names in `powers` do not match any column of `x`: ",
                  paste(unknown_names, collapse = ", "), "."), call. = FALSE)
    }
    
    if (!all(vapply(powers, is.numeric, logical(1L)))) {
      stop("! All elements of `powers` must be numeric vectors.", call. = FALSE)
    }
    
    powers <- lapply(powers, sort)
    
    pw_lengths <- vapply(powers, length, integer(1L))

    too_short  <- names(pw_lengths[pw_lengths < 2L])
    if (length(too_short) > 0L) {
      stop(paste0("! Each element of `powers` must contain at least two values. ",
                  "Insufficient values for: ", paste(too_short, collapse = ", "), "."),
           call. = FALSE)
    }
    power_list <- modifyList(power_list, Filter(Negate(is.null), powers))
  }
  
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
  
  # Validate p_adjust_method ----------------------------------------------------
  if (!is.character(p_adjust_method) || length(p_adjust_method) != 1L ||
      is.na(p_adjust_method)) {
    stop("! `p_adjust_method` must be a single character string.", call. = FALSE)
  }
  
  p_adjust_method <- match.arg(
    p_adjust_method,
    choices = c("none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")
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
  
  # Set default min_improvement -------------------------------------------------
  if (is.null(min_improvement)) {
    min_improvement <- switch(criterion, pvalue = p_interact, aic = 2, bic = 2)
  }
  
  # 
  if (is.null(weights)) weights <- rep.int(1, nobs)
  has_offset <- !is.null(offset)
  
  if (is.null(offset))  offset  <- rep.int(0, nobs)
  
  # Expand scalars and guarantee names on all per-variable vectors.
  # fit_mfp() and preprocess_data() rely on named vectors for
  # correct alignment after group_var is removed.
  if (length(select) == 1L) select <- rep(select, nvars)
  if (length(alpha)  == 1L) alpha  <- rep(alpha,  nvars)
  if (length(center) == 1L) center <- rep(center, nvars)
  select <- setNames(select, vnames)
  alpha  <- setNames(alpha,  vnames)
  center <- setNames(center, vnames)
  
  # When criterion = "pvalue", force_max_fp is achieved by setting alpha = 1,
  # which ensures the significance test always accepts the most complex FP form.
  # This gives force_max_fp a consistent interface across all criteria.
  if (criterion == "pvalue" && any(force_max_fp)) {
    alpha[force_max_fp] <- 1
  }
  
  if (is.null(shift)) {
    shift <- rep(NA_real_, nvars)
  } else if (length(shift) == 1L) {
    shift <- rep(shift, nvars)
  }
  
  shift <- setNames(shift, vnames)
  
  shift_missing <- is.na(shift)
  
  if (any(shift_missing)) {
    shift[shift_missing] <- apply(
      x[, shift_missing, drop = FALSE],
      2L,
      find_shift_factor
    )
  }
  
  # NOTE: scale is computed AFTER the shift has been applied to x (see below),
  # because find_scale_factor() must operate on the shifted variable to match
  # standalone mfp2. Here we only normalise a user-supplied scale; automatic
  # computation is deferred until after shifting.
  if (is.null(scale)) {
    scale <- rep(NA_real_, nvars)
  } else if (length(scale) == 1L) {
    scale <- rep(scale, nvars)
  }
  
  scale <- setNames(scale, vnames)
  
  if (is.null(control)) {
    control <- if (family_string == "cox") survival::coxph.control() else
      stats::glm.control()
  }
  
  # Retain only valid keep entries ---------------------------------------
  keep <- intersect(keep, vnames)
  
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
  
  # Cascade only after unsupported cont_var flags have been removed.
  # For adjustment variables, spike -> catzero -> zero is preserved. For
  # cont_vars, spike/catzero were already reset above, so this cascade will not
  # create zero handling unless zero_vars was explicitly supplied.
  catzero_flag[spike_flag] <- TRUE
  zero_flag[catzero_flag]  <- TRUE
  
  # Reset zero/catzero/spike for variables that contain only positive values -----
  if (any(zero_flag | catzero_flag | spike_flag)) {
    vars_pos <- names(zero_flag)[zero_flag | catzero_flag | spike_flag]
    
    all_positive <- vars_pos[
      apply(
        x[, vars_pos, drop = FALSE],
        2L,
        function(col) all(col > 0, na.rm = TRUE)
      )
    ]
    
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
  binary_vars  <- apply(x, 2L,
                        function(col) length(unique(col[!is.na(col)])) == 2L)
  binary_names <- vnames[binary_vars]
  if (length(binary_names) > 0L) {
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
    acd_vars <- unique(intersect(acd_vars, vnames))
    
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
  
  # Set degrees of freedom per variable ----------------------------------------
  # assign_df() handles both scalar and per-variable vector df_default,
  # and applies all three cardinality rules (see ?assign_df).
  df_list <- setNames(assign_df(x = x, df_default = df), vnames)
  
  # Reset shift to 0 for zero/catzero variables and linear terms ---------------
  # Only the positive part of these variables is FP-transformed so no shift
  # is needed (matching the logic in mfp2.default()).
  shift_to_zero <- rep(FALSE, nvars)
  names(shift_to_zero) <- vnames
  shift_to_zero[names(zero_flag)[zero_flag]]       <- TRUE
  shift_to_zero[names(catzero_flag)[catzero_flag]] <- TRUE
  shift_to_zero[vnames[df_list == 1L]]             <- TRUE
  shift[shift_to_zero] <- 0
  
  # Apply shift and scale to x -------------------------------------------------
  x <- sweep(x, 2L, shift, "+")
  
  # Compute scale on the SHIFTED x, matching standalone mfp2. find_scale_factor()
  # must operate on the shifted variable; computing it on the raw (unshifted) x
  # would give a different scaling factor for any variable with a non-zero shift.
  # User-supplied scale (normalised earlier) is left untouched.
  scale_missing <- is.na(scale)
  
  if (any(scale_missing)) {
    scale[scale_missing] <- apply(
      x[, scale_missing, drop = FALSE],
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
    xcheck  <- x[, check_names, drop = FALSE]
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
  
  x <- sweep(x, 2L, scale, "/")
  
  # Prepare stratification for Cox models --------------------------------------
  istrata <- strata
  if (family_string == "cox" && !is.null(strata)) {
    istrata <- as.integer(survival::strata(strata, shortlabel = TRUE))
  }
  
  # Apply subset ---------------------------------------------------------------
  if (!is.null(subset)) {
    x       <- x[subset, , drop = FALSE]
  
    y <- if (family_string == "cox" || is.matrix(y)) {
      y[subset, , drop = FALSE]
    } else {
      y[subset]
    }
    weights <- weights[subset]
    offset  <- offset[subset]
    if (!is.null(istrata)) istrata <- istrata[subset]
  }
  
  # A dataset with two groups can pass initial validation, then subset can 
  # remove all observations from one group. Later internals may fail inside 
  # dummy creation or interaction fitting with a less informative error.
  
  if (length(unique(x[, group_var])) < 2L) {
    stop("`group_var` must contain at least two groups after subsetting.", call. = FALSE)
  }
  
  group_counts <- table(x[, group_var])
  small_groups <- names(group_counts)[group_counts < 3L]
  if (length(small_groups) > 0L) {
    warning(
      "Some groups have fewer than 3 observations after subsetting: ",
      paste(small_groups, collapse = ", "),
      ". FP model fitting may be unreliable.",
      call. = FALSE
    )
  }
  
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
      cat(sprintf(
        "\ni Winsorising %d cont_var(s) at percentiles (%g, %g):\n",
        length(cont_vars), winsorize_probs[1L], winsorize_probs[2L]
      ))
      if (length(exclude_from_win) > 0L) {
        excluded <- intersect(exclude_from_win, cont_vars)
        if (length(excluded) > 0L) {
          cat(sprintf(
            "i Excluded from Winsorisation (zero-handled): %s\n",
            paste(excluded, collapse = ", ")
          ))
        }
      }
      print(winsorize_limits)
    }
  }
  
  # Fit the MFPI model ---------------------------------------------------------
  fit <- fit_mfpi(
    x                 = x,
    y                 = y,
    family            = family,         # GLM object or "cox"
    family_string     = family_string,  # character string for branching
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
    strata            = istrata,
    nocenter          = nocenter,
    acd_vars          = acd_flag,
    zero_vars         = zero_flag,
    catzero_vars      = catzero_flag,
    spike_vars        = spike_flag,
    min_prop          = min_prop,
    max_prop          = max_prop,
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
    p_adjust_method   = p_adjust_method
  )
  
  # Attach Winsorisation metadata for transparency in the returned object
  fit$winsorize         <- isTRUE(winsorize)
  fit$winsorize_probs   <- if (isTRUE(winsorize)) winsorize_probs else NULL
  fit$winsorize_limits  <- winsorize_limits
  fit$p_adjust_method   <- p_adjust_method
  fit$p_interact        <- p_interact
  fit$min_improvement   <- min_improvement
  fit$cont_var_forms    <- cont_var_forms
  
  # Developer note:
  # Store prediction-time reconstruction metadata on the returned mfpi object.
  # `predict.mfpi()` needs the preprocessed training matrix and model settings to
  # rebuild group-specific fitted functions f_j(x) without relying on
  # precomputed fitted functions from the flex functions.
  fit$x_train_scaled <- x
  fit$family_string  <- family_string
  fit$zero_vars      <- zero_flag
  fit$center_type    <- center_type
  fit$cont_vars      <- cont_vars
  
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
#' @section Parameters not available in the formula interface:
#' The following \code{mfpi.default()} parameters must be supplied as scalar
#' arguments only; per-variable specification via \code{fp()} or \code{fp2()} is not
#' supported: \code{group_var}, \code{cont_vars}, \code{cont_var_forms},
#' \code{flex}, \code{p_interact}, \code{min_improvement},
#' \code{include_group_var}, \code{show_models}, \code{cycles},
#' \code{criterion}, \code{keep}, \code{xorder}, \code{ties}, \code{strata},
#' \code{nocenter}, \code{min_prop}, \code{max_prop}, \code{ftest},
#' \code{control}, \code{verbose}, \code{digits}.
#' @section Strata and offset in the formula:
#' For Cox models, \code{strata()} terms may be included directly in the
#' formula (e.g. \code{Surv(t, d) ~ fp(age) + strata(centre)}). If
#' \code{strata} is also supplied as an argument, the formula value is used
#' and a warning is issued. Similarly, an \code{offset()} term in the formula
#' takes precedence over the \code{offset} argument.
#'
#' @seealso [mfp2::print.mfpi()], [mfp2::summary.mfpi()], [mfp2::plot.mfpi()], 
#' [mfp2::mfp2()]
#' @export
mfpi.formula <- function(formula,
                         data,
                         group_var         = NULL,
                         cont_vars         = NULL,
                         cont_var_forms    = NULL,
                         flex              = c("flex1", "flex2", "flex3", "flex4"),
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
                         ties              = c("breslow", "efron", "exact"),
                         strata            = NULL,
                         nocenter          = NULL,
                         zero_vars         = NULL,
                         catzero_vars      = NULL,
                         spike_vars        = NULL,
                         min_prop          = 0.05,
                         max_prop          = 0.95,
                         ftest             = FALSE,
                         control           = NULL,
                         winsorize         = FALSE,
                         winsorize_probs   = c(0.01, 0.99),
                         center_type       = c("grand", "group"),
                         p_adjust_method   = "none",
                         verbose           = TRUE,
                         digits            = 3,
                         ...) {
  
  call <- match.call()
  
  xorder      <- match.arg(xorder)
  criterion   <- match.arg(criterion)
  flex        <- match.arg(flex)
  ties        <- match.arg(ties)
  center_type <- match.arg(center_type)
  
  family_info   <- normalize_family_argument(family)
  family_string <- family_info$family_string
  # ---------------------------------------------------------------------------
  # Input validation (formula-specific constraints)
  # ---------------------------------------------------------------------------
  if (missing(data))
    stop("! data argument is missing.\n",
         "i An input data.frame is required for the formula interface.",
         call. = FALSE)
  
  if (is.null(colnames(data)))
    stop("! data must have column names.", call. = FALSE)
  
  if (missing(formula))
    stop("! formula is missing.", call. = FALSE)
  
  if (!inherits(formula, "formula"))
    stop("! method is only for formula objects.", call. = FALSE)
  
  # The formula interface only supports scalar defaults; per-variable settings
  # must be supplied via fp()/fp2() terms in the formula.
  if (length(df) != 1L)
    stop("! df must be a single numeric.\n",
         "i Use fp() or fp2() in the formula to set per-variable df values.", call. = FALSE)
  
  if (length(alpha) != 1L)
    stop("! alpha must be a single numeric.\n",
         "i Use fp() or fp2() in the formula to set per-variable alpha values.", call. = FALSE)
  
  if (length(select) != 1L)
    stop("! select must be a single numeric.\n",
         "i Use fp() or fp2() in the formula to set per-variable select values.", call. = FALSE)
  
  if (!is.null(scale) && length(scale) != 1L)
    stop("! scale must be a single numeric or NULL.\n",
         "i Use fp() or fp2() in the formula to set per-variable scale values.", call. = FALSE)
  
  if (!is.null(shift) && length(shift) != 1L) {
    stop("! shift must be a single numeric or NULL.\n",
         "i Use fp() or fp2() in the formula to set per-variable shift values.", call. = FALSE)
  }
  
  if (!is.numeric(df) || anyNA(df) || !is.finite(df) || df <= 0) {
    stop(
      "! `df` must be a single positive finite numeric value.\n",
      "i Use `fp()` or `fp2()` in the formula to set per-variable `df` values.",
      call. = FALSE
    )
  }
  
  if (df != 1L && df %% 2L != 0L) {
    stop(
      "! `df` must be 1 or an even positive integer.\n",
      "i Use `fp()` or `fp2()` in the formula to set per-variable `df` values.",
      call. = FALSE
    )
  }
  
  if (!is.numeric(alpha) || anyNA(alpha) || !is.finite(alpha) ||
      alpha < 0 || alpha > 1) {
    stop(
      "! `alpha` must be a single numeric value in [0, 1].\n",
      "i Use `fp()` or `fp2()` in the formula to set per-variable `alpha` values.",
      call. = FALSE
    )
  }
  
  if (!is.numeric(select) || anyNA(select) || !is.finite(select) ||
      select < 0 || select > 1) {
    stop(
      "! `select` must be a single numeric value in [0, 1].\n",
      "i Use `fp()` or `fp2()` in the formula to set per-variable `select` values.",
      call. = FALSE
    )
  }
  
  if (!is.logical(center) || length(center) != 1L || anyNA(center)) {
    stop(
      "! `center` must be a single logical value.\n",
      "i Use `fp()` or `fp2()` in the formula to set per-variable `center` values.",
      call. = FALSE
    )
  }
  
  if (!is.null(scale)) {
    if (!is.numeric(scale) || length(scale) != 1L) {
      stop(
        "! `scale` must be a single numeric value, `NA`, or `NULL`.\n",
        "i Use `fp()` or `fp2()` in the formula to set per-variable `scale` values.",
        call. = FALSE
      )
    }
    if (!is.na(scale) && (!is.finite(scale) || scale <= 0)) {
      stop(
        "! `scale` must be positive, `NA`, or `NULL`.\n",
        "i Use `fp()` or `fp2()` in the formula to set per-variable `scale` values.",
        call. = FALSE
      )
    }
  }
  
  if (!is.null(shift)) {
    if (!is.numeric(shift) || length(shift) != 1L) {
      stop(
        "! `shift` must be a single numeric value, `NA`, or `NULL`.\n",
        "i Use `fp()` or `fp2()` in the formula to set per-variable `shift` values.",
        call. = FALSE
      )
    }
    if (!is.na(shift) && !is.finite(shift)) {
      stop(
        "! `shift` must be finite, `NA`, or `NULL`.\n",
        "i Use `fp()` or `fp2()` in the formula to set per-variable `shift` values.",
        call. = FALSE
      )
    }
  }
  
  if (!is.null(powers) && !is.list(powers)) {
    stop("! powers must be a named list or NULL.", call. = FALSE)
  }
  
  
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
  
  # ---------------------------------------------------------------------------
  # Build model frame and extract y, x
  # ---------------------------------------------------------------------------
  mf <- stats::model.frame(
    formula,
    data = data,
    drop.unused.levels = TRUE,
    na.action = stats::na.fail
  )
  
  labels <- attr(terms(mf), "term.labels")
  if (length(labels) == 0L)
    stop("! No predictors found in formula. At least one predictor is required.",
         call. = FALSE)
  
  # ---------------------------------------------------------------------------
  # Handle strata terms for Cox models
  # ---------------------------------------------------------------------------
  specials      <- "strata"
  terms_formula <- terms(formula, specials = specials, data = data)
  terms_drop    <- NULL
  
  if (!is.null(attr(terms_formula, "specials")$strata)) {
    if (family_string == "cox") {
      if (!is.null(call$strata))
        warning("i strata appear in both the formula and as an argument.\n",
                "i Formula strata are used; the argument is ignored.",
                call. = FALSE)
      
      stemp      <- survival::untangle.specials(terms_formula,
                                                special = "strata", order = 1)
      strata     <- if (length(stemp$vars) == 1L) mf[[stemp$vars]]
      else mf[, stemp$vars]
      terms_drop <- stemp$terms
    } else {
      stop("! strata are only allowed for Cox models.\n",
           "i Please remove strata terms from the formula.", call. = FALSE)
    }
  }
  
  terms_model <- if (!is.null(terms_drop)) terms_formula[-terms_drop]
  else terms_formula
  
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
    # Always extract offset from model frame when formula contains offset()
    offset <- as.vector(stats::model.offset(mf))
  }
  
  # ---------------------------------------------------------------------------
  # Extract y and x
  # ---------------------------------------------------------------------------
  y <- stats::model.extract(mf, "response")
  
  if (family_string != "cox" && survival::is.Surv(y))
    stop(sprintf(
      "! Response is a Surv object but family = '%s'. Set family = 'cox'.",
      family_string
    ), call. = FALSE)
  
  x <- stats::model.matrix(terms_model, mf)
  
  # Developer note: save formula reconstruction metadata before modifying the
  # model matrix. The subsequent intercept removal, factor-group handling, and
  # fp()/fp2() renaming can otherwise drop or obscure these attributes. These
  # fields allow a future predict.mfpi() method to rebuild the same design from
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
  
  # ---------------------------------------------------------------------------
  # Handle group_var that was a factor in the original data.
  # model.matrix() expands factors into dummy columns, so the original column
  # name disappears. mfpi.default() expects group_var as a single numeric
  # column. Fix: remove the dummy columns and add back the raw integer codes.
  # ---------------------------------------------------------------------------
  if (!is.null(group_var) && !group_var %in% colnames(x)) {
    # group_var was expanded by model.matrix. Identify and remove its dummies
    # using the saved assign attribute.
    term_labels <- attr(terms_model, "term.labels")
    gv_term_idx <- which(term_labels == group_var)
    if (length(gv_term_idx) == 1L) {
      gv_cols <- which(x_assign == gv_term_idx)
      if (length(gv_cols) > 0L) {
        x        <- x[, -gv_cols, drop = FALSE]
        x_assign <- x_assign[-gv_cols]
      }
    }
    # Add back the raw group_var values as a single numeric column.
    # For factors, use integer codes (preprocess_data remaps to 0-based).
    gv_raw <- mf[[group_var]]
    if (is.factor(gv_raw)) {
      # Try to recover original numeric levels; fall back to integer codes
      num_levels <- suppressWarnings(as.numeric(levels(gv_raw)))
      if (anyNA(num_levels)) {
        gv_numeric <- as.integer(gv_raw)
      } else {
        gv_numeric <- num_levels[as.integer(gv_raw)]
      }
    } else {
      gv_numeric <- as.numeric(gv_raw)
    }
    x <- cbind(x, gv_numeric)
    colnames(x)[ncol(x)] <- group_var
  }
  
  nx      <- ncol(x)
  names_x <- colnames(x)
  
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
    
    names_x <- replace(names_x, fp_x_pos, fp_vars_pre)
    colnames(x) <- names_x
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
    subset            = subset,
    family            = family,
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
    min_prop          = min_prop,
    max_prop          = max_prop,
    ftest             = ftest,
    control           = control,
    winsorize         = winsorize,
    winsorize_probs   = winsorize_probs,
    center_type       = center_type,
    p_adjust_method   = p_adjust_method,
    verbose           = verbose,
    digits            = digits,
    ...
  )
  fit$formula_interface <- TRUE
  fit$formula <- formula
  fit$formula_terms <- stats::delete.response(terms_model)
  fit$formula_contrasts <- x_contrasts
  fit$formula_xlevels <- x_xlevels
  
  fit
}