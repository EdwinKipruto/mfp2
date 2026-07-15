#' Model Interactions Between a Categorical and Continuous Covariates
#'
#' `mfpi()` investigates interactions between a categorical variable
#' (`group_var`) and one or more continuous covariates, using fractional
#' polynomial (FP) transformations to capture nonlinear effects. The categorical
#' variable is treated as a factor internally; in a typical randomised controlled
#' trial it represents treatment allocation, with the lowest value as the control
#' arm. Confounders and other prognostic variables are adjusted for via FP
#' transformations selected by the MFP algorithm implemented in the
#' \code{mfp2()}.
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
#' **`flex1` (least flexible).** FP powers are selected for the main
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
#' **`flex3` (default).** Separate FP powers are estimated for the main-effects
#' model and for the within-group functions, though the within-group powers are
#' still constrained to be equal across groups. Because the main-effects and
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
#' \code{cont_var_forms} default to \code{"fp1"}. Data-driven selection among
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
#' @section Details on the `subset` argument:
#' In `mfpi.formula()`, the supplied subset expression is evaluated exactly once
#' using standard formula lookup: names are resolved first from columns of `data`
#' and then from the formula environment. Therefore expressions such as
#' `subset = age >= 50`, a logical data column such as `subset = eligible`, and
#' caller-defined row indices are supported; the subset does not itself have to
#' be stored in `data`.
#'
#' In `mfpi.default()`, there is no formula data mask. The supplied logical or
#' numeric/integer vector is evaluated normally in the caller, and matrix column
#' names are not exposed as standalone variables. Use an explicit expression
#' such as `subset = x[, "age"] >= 50` when deriving the subset from a matrix
#' column. In both interfaces, continuous shift and scale calculations retain
#' their full-data source, while model selection and fitting use only the
#' retained observations. Numeric row positions must be unique; unique positions
#' retain the order supplied by the user and are not sorted automatically.
#'
#' When `subset = NULL`, the formula interface reuses the complete evaluated
#' model frame directly. It does not invoke retained-row resolution or rebuild
#' factor coding. Subset-specific row and factor processing occurs only when a
#' non-`NULL` subset is supplied.
#'
#' @section Categorical adjustment terms:
#' Categorical predictors other than \code{group_var} are adjustment terms;
#' they are not candidates for the MFPI interaction test. In the formula
#' interface, unordered and ordered factors are expanded using the fitted
#' \code{model.matrix()} contrasts, and all contrast columns generated by one
#' factor are represented by one conceptual term. The complete block is
#' retained or removed by a joint adjustment-model test. Ordered factors retain
#' their polynomial contrasts unless the formula specifies another contrast.
#' Simple wrappers such as \code{factor(stage)} and \code{ordered(stage)}
#' use \code{stage} as the conceptual adjustment-term name. The generated
#' design columns and fitted coefficient names retain their standard
#' \code{model.matrix()} names.
#'
#' In the default matrix interface, categorical columns must already be encoded
#' numerically and grouped with \code{term_groups}. For example,
#' \code{term_groups = list(stage = c("stageII", "stageIII"))}
#' declares the two columns to be one adjustment term named \code{stage}.
#' With \code{subset}, automatic continuous shift and scale calculations retain
#' their full-data source. The formula interface performs extra factor processing
#' only for a non-`NULL` subset. It preserves factor-specific custom contrasts when
#' all levels remain, and it drops unused levels and regenerates default treatment
#' or ordered polynomial contrasts when no custom contrast was supplied. If a
#' predictor with factor-specific custom contrasts loses a level, fitting stops
#' because a unique reduced contrast basis cannot be inferred safely. The matrix
#' interface has only the supplied numeric columns and cannot reconstruct their
#' factor origin; it rejects a grouped block whose estimable dimension decreases
#' after subsetting. Formula and matrix results may therefore differ when a subset
#' removes factor levels.
#'
#' Columns not mentioned in \code{term_groups} remain singleton terms. Group
#' names and member columns must be unique, and a column cannot belong to more
#' than one group.
#'
#' Grouped categorical terms are fixed linear blocks. Their effective settings
#' are \code{df = 1}, power 1, shift 0, scale 1, no centering, and no ACD,
#' zero, catzero, spike-at-zero, or \code{force_max_fp} processing.
#' \code{keep} may name the conceptual term or one of its member columns; in
#' either case the complete block is retained. A grouped term or member column
#' cannot be listed in \code{cont_vars}, because MFPI interaction candidates
#' must be singleton continuous variables.
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
#'     \item{\code{spike_vars}}{Uses the spike-at-zero algorithm and implies
#'  \code{catzero_vars} while SAZ remains active after eligibility checks.}
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
#' @section Spike-at-zero eligibility:
#' Spike-at-zero handling in \code{mfpi()} follows the same eligibility rule as
#' in \code{mfp2()}. A requested spike-at-zero adjustment variable remains
#' eligible only if both the structural-zero component and the positive
#' continuous component are sufficiently represented. Let \eqn{p_0} be the
#' structural-zero proportion and \eqn{p_+ = 1 - p_0} be the positive-observation
#' proportion. SAZ is applied only if both \eqn{p_0} and \eqn{p_+} are at least
#' \code{min_saz_component_prop}.
#'
#' Variables in \code{cont_vars}, which define the continuous variables tested
#' for interaction with the grouping variable, are not treated as spike-at-zero
#' variables during the MFPI interaction test. Spike-at-zero handling applies to
#' eligible adjustment variables.
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
#'   For \code{mfpi.default()} only. A numeric matrix or data frame of predictor
#'   variables. If \code{x} is a data frame, \code{group_var} may be factor,
#'   character, logical, integer, or numeric. All non-group variables must be
#'   numeric, integer, or logical. Internally, \code{group_var} is remapped to
#'   consecutive integer levels, while its original labels are preserved for
#'   printing, summaries, plotting, and prediction metadata. The matrix must not contain an intercept
#'   column; binary variables should be coded as 0/1. Multi-level categorical
#'   adjustment variables must be expanded into design columns and supplied as
#'   one conceptual block through \code{term_groups}. Must be free of missing
#'   values.
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
#'   named in \code{formula}. Unordered and ordered factors other than
#'   \code{group_var} are expanded by \code{model.matrix()} and their contrast
#'   columns are treated as one fixed linear adjustment term for joint selection.
#' @param term_groups
#'   For \code{mfpi.default()} only. Optional named list defining conceptual
#'   categorical adjustment terms. Each list name is the term name and each
#'   value is a character vector of one or more columns of \code{x}; for
#'   example, \code{list(stage = c("stageII", "stageIII"))}. A one-column
#'   mapping supports a binary factor whose conceptual name differs from its
#'   dummy column name. The columns must exist, must be unique within and across
#'   groups, and may not
#'   include \code{group_var}. Unlisted columns remain singleton terms. Each
#'   grouped block is selected jointly and fitted as a fixed linear effect. It
#'   cannot be listed in \code{cont_vars} or use ACD, zero, catzero,
#'   spike-at-zero, or \code{force_max_fp}. Formula factors are grouped
#'   automatically, so this argument is not available in \code{mfpi.formula()}.
#' @param group_var
#' Character scalar identifying the categorical grouping
#'   variable. In \code{mfpi.default()}, this variable may be factor,
#'   character, logical, integer, or numeric when \code{x} is a data frame. It
#'   is converted to internal consecutive integer levels for fitting. For a
#'   factor \code{group_var}, the first element of \code{levels(group_var)} is
#'   used as the reference level. For a character \code{group_var}, the first
#'   distinct value encountered in the data is used as the reference level. For
#'   logical, integer, or numeric \code{group_var}, the smallest observed value
#'   is used as the reference level. The original labels are preserved for
#'    printing, summaries, plotting, and prediction metadata.
#' @param cont_vars
#'   A non-empty character vector naming singleton continuous columns of
#'   \code{x} whose interactions with \code{group_var} are to be
#'   investigated. Each name must identify an existing numeric column. A
#'   conceptual grouped term, or any raw column belonging to one, is not
#'   permitted because categorical blocks are adjustment terms only. Binary
#'   variables (2 or fewer unique values) are not permitted and cause an error.
#'   Variables with 5 or fewer unique values trigger a warning because FP
#'   transformation may be unreliable for near-categorical variables.
#' @param cont_var_forms
#'   Optional named character vector specifying the functional form to use for
#'   each variable in \code{cont_vars}. Names must be a subset of
#'   \code{cont_vars}; values must each be one of \code{"linear"},
#'   \code{"fp1"}, or \code{"fp2"}. Variables in \code{cont_vars} that are not
#'   named in \code{cont_var_forms} are assigned \code{"fp1"} by default.
#'   \code{NULL} (the default) assigns \code{"fp1"} to all variables.
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
#'   Example: to test \code{age} as FP2 while leaving all other
#'   \code{cont_vars} at the default (\code{"fp1"}):
#'   \preformatted{cont_var_forms = c(age = "fp2")}
#' @param flex
#'   A character string controlling how FP powers are estimated and constrained
#'   across groups. One of `"flex1"`, `"flex2"`, `"flex3"` (default), or
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
#'   Logical. Whether \code{group_var} should enter the stage-one MFP
#'   adjustment model. Default is \code{FALSE}. When \code{TRUE}, the
#'   non-reference group dummy columns are registered as one conceptual
#'   categorical term named by \code{group_var} and forced into the adjustment
#'   model with \code{select = 1}. The block is excluded from the ordinary
#'   stage-two adjustment matrix because the interaction models add the group
#'   main effects separately.
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
#'   Optional subset of observations used for model selection and fitting. In
#'   the formula interface, this may be an expression involving columns of
#'   `data` and objects in the formula environment; it is evaluated exactly once,
#'   with data columns taking precedence. In the matrix interface, it is an
#'   ordinary logical vector with one value per observation or a numeric vector
#'   of positive integer row indices evaluated in the caller. Missing,
#'   non-finite, zero, negative, out-of-range, non-integer, and duplicated
#'   indices are not allowed. Unique numeric indices retain their supplied order.
#'   Default is `NULL`, meaning that all observations
#'   are used. Automatic shift and scale values for continuous predictors are
#'   computed from the full data. The formula method then drops unused factor
#'   levels and rebuilds categorical design columns from the retained
#'   observations. The matrix method cannot reconstruct factor contrasts and
#'   stops if a mapped \code{term_groups} block loses estimable dimension after
#'   subsetting.
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
#'   corresponding column of `x` and must contain at least one finite value.
#'   A single non-unity candidate power is valid (e.g. `powers = list(x = 2)`);
#'   see the `powers` argument of \code{mfp2()} for details.
#' @param ties
#'   A character string specifying the method for handling tied event times in
#'   Cox regression. One of `"breslow"` (default), `"efron"`, or `"exact"`.
#'   Ignored for non-Cox families. See [survival::coxph()] for details.
#' @param strata
#'   Optional stratification factor(s) for Cox models. May be a vector, factor,
#'   matrix, data frame, or \code{survival::strata()} object with one row per
#'   observation. Multiple supplied columns are combined with
#'   \code{survival::strata(shortlabel = TRUE)}. Default \code{NULL} means no
#'   stratification. Ordinary vector/factor strata are kept as high-level
#'   values for formula-based Cox fitting; integer conversion is restricted to
#'   the low-level \code{coxph.fit()} path. In the formula interface,
#'   \code{strata()} and \code{survival::strata()} terms may also be supplied
#'   directly in the formula; formula strata take precedence over this argument.
#' @param nocenter
#'   A numeric vector of values passed to [survival::coxph()]. Cox models only;
#'   ignored otherwise.
#' @param acd_vars
#'   For \code{mfpi.default()} only; not a parameter of \code{mfpi.formula()}
#'   (passing it there raises an error). An optional character vector naming
#'   continuous variables to be transformed
#'   via the approximate cumulative distribution (ACD) transformation before FP
#'   selection. The transformed variable is named `A(x)`. ACD transformation is
#'   disabled for variables listed in `cont_vars`. In the formula interface, use
#'   \code{fp(variable, acd = TRUE)} or \code{fp2(variable, acd = TRUE)} instead.
#'   See [mfp2::mfp2()] for details.
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
#'    \code{spike_vars} implies \code{catzero_vars} and therefore \code{zero_vars}
#'   while the SAZ algorithm remains active after eligibility checks. A variable
#'   that is not in \code{cont_vars} may not appear in both \code{zero_vars} and
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
#' @param min_saz_component_prop
#'   Numeric in \eqn{(0, 0.5)}. Minimum required proportion in each component of
#'   a spike-at-zero covariate: the structural-zero component and the positive
#'   continuous component. Default \code{0.10}. A requested spike-at-zero
#'   adjustment variable is retained for SAZ modelling only if both the
#'   structural-zero proportion and the positive-observation proportion are at
#'   least this value.
#' @param ftest
#'   Logical. Whether to use an F-test rather than a chi-square
#'   likelihood-ratio test when computing p-values for Gaussian models.
#'   When \code{TRUE} and \code{family = "gaussian"}, the F-test is applied
#'   to both the adjustment-variable selection step and the interaction test 
#'   for each variable in \code{cont_vars}. Recommended when the sample size is
#'   small. Default \code{FALSE}. Has no effect for non-Gaussian families or 
#'   when \code{criterion} is not \code{"pvalue"}. 
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
#'   For S3 method compatibility only. Additional arguments are not currently
#'   used; supplying arguments through \code{...} is an error. In the formula
#'   interface, variable-specific FP options should be supplied inside
#'   \code{fp()} or \code{fp2()}.
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
#'   \item{\code{adjustment_model}}{The full grouped adjustment model object
#'     returned by the MFP core. Its \code{fp_terms} component has one row per
#'     conceptual term, including one row for each categorical block.}
#'   \item{\code{term_to_columns}}{Named list mapping conceptual adjustment
#'     terms to the raw predictor columns used by the fitted object. Singleton
#'     entries contain one column; categorical entries contain the complete
#'     contrast block.}
#'   \item{\code{adjustment_term_to_columns}}{The same mapping restricted to
#'     terms available to the stage-one adjustment model after preprocessing.}
#'   \item{\code{univariable_interactions}}{The full return value of the
#'     internal \code{evaluate_interactions()} call. Contains all
#'     intermediate results and is useful for programmatic access.}
#'   \item{\code{group_var}}{The name of the grouping variable.}
#'   \item{\code{group_levels_new}}{A zero-based integer index
#'     (0, 1, 2, \ldots) used internally after \code{group_var} is remapped by
#'     \code{preprocess_data()}.}
#'   \item{\code{group_levels_original}}{User-facing labels of
#'     \code{group_var}, used for printing, summaries, plotting, and prediction
#'     metadata.}
#'   \item{\code{group_level_map}}{A data frame mapping user-facing
#'     \code{group_var} labels to input values and internal integer levels.}
#'   \item{\code{flex}}{The flexibility level used (\code{"flex1"} through
#'     \code{"flex4"}).}
#'   \item{\code{family}}{The regression family used for fitting: a GLM
#'     family object for Gaussian, binomial, and Poisson models, or the
#'     character string \code{"cox"} for Cox models.}
#'   \item{\code{family_string}}{Character scalar giving the normalized family
#'     name used internally for family-specific branching. One of
#'     \code{"gaussian"}, \code{"binomial"}, \code{"poisson"}, or
#'     \code{"cox"}.}
#'   \item{\code{scale}}{Named numeric vector of scale factors used to transform
#'     continuous predictors onto the shifted/scaled FP fitting scale. These
#'     values are reused by \code{predict.mfpi()} and \code{plot.mfpi()} and
#'     should not be recomputed from new data.}
#'   \item{\code{shift}}{Named numeric vector of shifts added to continuous
#'     predictors before scaling. These values define the FP transformation
#'     scale used during fitting and are reused during prediction.}
#'   \item{\code{center_type}}{Character scalar describing the centering
#'     convention used for FP-transformed terms.}
#'   \item{\code{zero_vars}}{Named logical vector indicating which variables
#'     used structural-zero handling during FP transformation.}
#'   \item{\code{cont_vars}}{Character vector of continuous variables tested for
#'     interaction with \code{group_var}.}
#'   \item{\code{cont_var_forms}}{Named character vector giving the requested
#'     interaction functional form for each variable in \code{cont_vars}, for
#'     example \code{"linear"}, \code{"fp1"}, or \code{"fp2"}.}
#'   \item{\code{x_train_internal}}{Numeric matrix containing the
#'     post-preprocessing training predictor matrix used by MFPI. This matrix
#'     contains shifted/scaled predictors and the internally remapped
#'     \code{group_var} codes produced by \code{preprocess_data()}. It is used
#'     by \code{predict.mfpi()} when manual reconstruction of training-data
#'     predictions is required, for example when \code{newdata = NULL} and a
#'     replacement \code{newoffset} is supplied.}
#'   \item{\code{digits}}{Integer giving the number of significant digits used
#'     when printing interaction-test results.}   
#'   \item{\code{nobs}}{Number of observations used in model fitting.}
#'   \item{\code{show_models}}{Value of the \code{show_models} argument.}
#'   \item{\code{winsorize}}{Logical; whether Winsorisation was applied.}
#'   \item{\code{winsorize_probs}}{The probability cutoffs used for
#'     Winsorisation, or \code{NULL} when \code{winsorize = FALSE}.}
#'   \item{\code{winsorize_limits}}{A numeric matrix of Winsorisation limits,
#'     with rows \code{"lower"} and \code{"upper"} and columns named by
#'     \code{cont_vars}. Entries are \code{NA} for variables excluded from
#'     Winsorisation, for example zero-handled variables. \code{NULL} when
#'     \code{winsorize = FALSE}.}
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
#'   \item{\code{call}}{The matched call to \code{mfpi()}.}
#'   \item{\code{formula_interface}}{Logical; \code{TRUE} if the model was
#'     fitted via \code{mfpi.formula()}. Not present when fitted via
#'     \code{mfpi.default()}.}
#'   \item{\code{formula}}{The formula supplied to \code{mfpi.formula()}.
#'     Only present when fitted via the formula interface.}
#'   \item{\code{formula_terms}}{The \code{terms} object for the predictors
#'     (response removed), as used to build the model matrix. Only present
#'     when fitted via the formula interface.}
#'   \item{\code{formula_contrasts}}{The \code{contrasts} attribute of the
#'     model matrix built from \code{formula}/\code{data}, as needed to
#'     reconstruct the same design from new data. Only present when fitted
#'     via the formula interface.}
#'   \item{\code{formula_xlevels}}{Factor levels recorded from the training
#'     model frame (via \code{.getXlevels()}), needed to reconstruct the same
#'     design from new data. Only present when fitted via the formula
#'     interface.}
#'   \item{\code{formula_design_columns}}{Character vector of the final
#'     fit-time column names of \code{x}, after intercept removal,
#'     group-variable replacement, and \code{fp()}/\code{fp2()} renaming. Only
#'     present when fitted via the formula interface.}
#'   \item{\code{formula_term_to_columns}}{Named list mapping each conceptual
#'     predictor term to the design columns it generated. Simple factor wrappers
#'     are named by their source variable; their raw design-column names remain
#'     unchanged. Entries are used to reconstruct complete categorical
#'     adjustment blocks during prediction. Only present for formula fits.}
#'   \item{\code{formula_prediction_term_names}}{Named character vector mapping
#'     original formula term labels to the conceptual names used by MFPI
#'     selection and prediction. Only present for formula fits.}
#'   \item{\code{formula_factor_terms}}{Character vector naming the conceptual
#'     unordered and ordered factor terms treated as fixed-linear adjustment
#'     blocks. Only present for formula fits.}
#'   \item{\code{formula_offset_terms}, \code{formula_offset_xlevels}}{Formula
#'     offset metadata used by \code{predict.mfpi()} to reconstruct a
#'     formula-level offset from raw \code{newdata}. Only present when fitted
#'     via the formula interface; values are \code{NULL} when no formula offset
#'     was used.}
#'   \item{\code{formula_strata_terms}, \code{formula_strata_xlevels}}{Formula
#'     Cox-strata metadata corresponding to any \code{strata()} term supplied
#'     in the formula. Only present when fitted via the formula interface; values
#'     are \code{NULL} when no formula strata were used.}
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
#' @seealso \code{\link{summary.mfpi}}, \code{\link{mfp2}},
#'   \code{\link{print.mfpi}}, and \code{\link{plot.mfpi}}.
#' 
#' @examples
#' \dontrun{
#' data("prostate")
#'
#' # Basic MFPI fit using the simplest interaction structure.
#' fit <- mfpi(
#'   lpsa ~ fp(age) + svi + fp(pgg45) + fp(cavol) + fp(weight) +
#'     fp(bph) + fp(cp),
#'   data = prostate,
#'   cont_vars = c("cavol", "age"),
#'   group_var = "svi",
#'   center = FALSE,
#'   flex = "flex1",
#'   include_group_var = TRUE,
#'   show_models = FALSE,
#'   verbose = FALSE
#' )
#'
#' # Plot group-specific fitted functions and fitted-function differences.
#' # The difference plot is on the model scale, not a probability, odds ratio,
#' # rate ratio, hazard ratio, or full predicted response.
#' if (requireNamespace("patchwork", quietly = TRUE)) {
#'   plot(fit, terms = "cavol", plot_type = "both")
#' }
#'
#' # Pre-specify functional forms: cavol as FP2, age uses the default search.
#' fit_prespecified <- mfpi(
#'   lpsa ~ fp(age) + svi + fp(pgg45) + fp(cavol) + fp(weight) +
#'     fp(bph) + fp(cp),
#'   data = prostate,
#'   cont_vars = c("cavol", "age"),
#'   cont_var_forms = c(cavol = "fp2"),
#'   group_var = "svi",
#'   center = FALSE,
#'   flex = "flex1",
#'   include_group_var = TRUE,
#'   verbose = FALSE
#' )
#'
#' # Formula factors are grouped automatically in the adjustment model.
#' dat <- data.frame(
#'   y = rnorm(180),
#'   trt = factor(rep(c("control", "active"), each = 90)),
#'   age = runif(180, 30, 80),
#'   stage = factor(rep(c("I", "II", "III"), length.out = 180))
#' )
#' fit_factor <- mfpi(
#'   y ~ trt + fp(age) + stage,
#'   data = dat,
#'   group_var = "trt",
#'   cont_vars = "age",
#'   keep = "stage",
#'   verbose = FALSE
#' )
#' fit_factor$term_to_columns$stage
#'
#' # Matrix-interface categorical columns are grouped explicitly.
#' mm <- model.matrix(~ stage, data = dat)[, -1, drop = FALSE]
#' x <- cbind(trt = as.integer(dat$trt) - 1L, age = dat$age, mm)
#' fit_matrix <- mfpi(
#'   x, dat$y,
#'   group_var = "trt",
#'   cont_vars = "age",
#'   term_groups = list(stage = colnames(mm)),
#'   keep = "stage",
#'   verbose = FALSE
#' )
#' }
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
    ties              = c("breslow", "efron", "exact"),
    strata            = NULL,
    nocenter          = NULL,
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
  ties        <- match.arg(ties)
  center_type <- match.arg(center_type)
  
  # Resolve family -------------------------------------------------------------
  family_info <- normalize_family_argument(
    family,
    family_arg = deparse(substitute(family))
  )
  
  family        <- family_info$family
  family_string <- family_info$family_string
  
  # Step 2: Prepare x and group_var, then run basic sanity checks -------------
  # Formula dispatch can attach a full-row preprocessing design to the
  # subset-specific fitting matrix. Preserve that private attribute while
  # prepare_mfpi_default_x() normalizes x/group_var, then extract it into an
  # explicit preprocess_x matrix. Direct matrix calls use x for both roles.
  preprocess_attr <- attr(x, "mfp2_preprocess_x", exact = TRUE)
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
  
  if (anyDuplicated(vnames)) {
    stop("! `x` must have unique column names.", call. = FALSE)
  }
  
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
  
  # Step 5: Validate alpha/select/df, keep, and force_max_fp_vars -------------
  # Validate alpha, select, and df --------------------------------------------
  validate_probability_vector(alpha, "alpha", nvars)
  validate_probability_vector(select, "select", nvars)
  
  validate_numeric_vector(
    arg = df,
    name = "df",
    nvars = nvars,
    allow_null = FALSE,
    allow_na = FALSE,
    strictly_positive = FALSE
  )
  
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
  validate_numeric_vector(
    arg = shift,
    name = "shift",
    nvars = nvars,
    allow_null = TRUE,
    allow_na = TRUE,
    strictly_positive = FALSE
  )
  
  validate_numeric_vector(
    arg = scale,
    name = "scale",
    nvars = nvars,
    allow_null = TRUE,
    allow_na = TRUE,
    strictly_positive = TRUE
  )
  
  # Validate center ------------------------------------------------------------
  # Developer note: mfpi.formula() may pass a per-variable logical vector after
  # reading fp()/fp2() attributes, so the default method is the final validation
  # gate for both direct matrix calls and formula calls.
  validate_logical_vector(
    arg = center,
    name = "center",
    allowed_lengths = c(1L, nvars)
  )
  
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
  
  # Small helper to expand a scalar option to one named value per predictor.
  expand_to_named <- function(val, n, nms) {
    if (length(val) == 1L) val <- rep(val, n)
    setNames(val, nms)
  }
  
  # Expand scalars and guarantee names on all per-variable vectors.
  # fit_mfp() and preprocess_data() rely on named vectors for
  # correct alignment after group_var is removed.
  select <- expand_to_named(select, nvars, vnames)
  alpha  <- expand_to_named(alpha,  nvars, vnames)
  center <- expand_to_named(center, nvars, vnames)
  
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
  
  # Automatic shifts are deliberately estimated later, after SAZ eligibility and
  # the zero/catzero/spike cascade have reached their final retained state.
  # Retained zero-handled variables need shift = 0; variables reset from SAZ to
  # ordinary FP need ordinary positivity shifting.
  shift <- setNames(shift, vnames)
  
  
  # NOTE: scale is computed AFTER the shift has been applied to x (see below),
  # because find_scale_factor() must operate on the shifted variable to match
  # standalone mfp2. NULL or scalar settings control ordinary continuous
  # predictors, while explicitly mapped adjustment blocks retain their supplied
  # model-matrix basis under scale = 1. Per-column vectors remain authoritative.
  scale <- expand_scale_for_mapped_terms(
    scale = scale,
    vnames = vnames,
    term_to_columns = term_to_columns
  )
  scale <- setNames(scale, vnames)
  
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
  # assign_df() applies the cardinality rules after a scalar default has
  # been expanded. Explicitly mapped adjustment terms are supplied fixed linear
  # design blocks, so their columns receive df = 1 under a scalar default.
  # A per-column df vector is left unchanged and is checked below, preserving
  # the existing error for an explicitly nonlinear grouped member.
  df_default <- expand_scalar_df_for_mapped_terms(
    df = df,
    vnames = vnames,
    term_to_columns = term_to_columns
  )
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
  # Keep high-level strata labels here. Integer conversion belongs only in the
  # low-level Cox fitting path immediately before survival::coxph.fit(). This
  # mirrors survival::coxph() and prevents formula-strata labels from being lost.
  strata_keep <- strata
  
  if (family_string == "cox" && !is.null(strata_keep)) {
    if (!isTRUE(attr(strata_keep, "mfp2_strata_keep"))) {
      strata_keep <- if (inherits(strata_keep, "strata")) {
        strata_keep
      } else if (is.matrix(strata_keep) || is.data.frame(strata_keep)) {
        do.call(
          survival::strata,
          c(as.list(as.data.frame(strata_keep)), list(shortlabel = TRUE))
        )
      } else {
        strata_keep
      }
    }
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
    if (!is.null(strata_keep)) strata_keep <- strata_keep[subset]
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
  small_groups <- names(group_counts)[group_counts < 3L]
  if (length(small_groups) > 0L) {
    warning(
      "Some groups have fewer than 3 observations after subsetting: ",
      paste(small_groups, collapse = ", "),
      ". FP model fitting may be unreliable.",
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
#' For Cox models, \code{strata()} or \code{survival::strata()} terms may be
#' included directly in the formula, for example
#' \code{Surv(t, d) ~ fp(age) + strata(centre)}. Multiple strata variables are
#' combined using \code{survival::strata(shortlabel = TRUE)}. If \code{strata}
#' is also supplied as an argument, the formula value is used and a warning is
#' issued. Similarly, an \code{offset()} or \code{stats::offset()} term in the
#' formula takes precedence over the \code{offset} argument. Formula-offset
#' metadata are stored on the fitted object so \code{predict.mfpi()} can
#' reconstruct the offset from raw prediction \code{newdata}; if the needed raw
#' offset variables are absent, supply \code{newoffset} to \code{predict()}.
#'
#' @seealso [mfp2::print.mfpi()], [mfp2::summary.mfpi()], [mfp2::plot.mfpi()], 
#' [mfp2::mfp2()]
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
                         ties              = c("breslow", "efron", "exact"),
                         strata            = NULL,
                         nocenter          = NULL,
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
  ties        <- match.arg(ties)
  center_type <- match.arg(center_type)
  
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
  
  if (!is.null(weights)) {
    if (!is.numeric(weights) || length(weights) != n_data || anyNA(weights) ||
        any(!is.finite(weights)) || any(weights < 0)) {
      stop("! `weights` must be a finite, non-negative numeric vector with one value per row of `data`.", call. = FALSE)
    }
  }
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
      
      attr(strata, "mfp2_strata_keep") <- TRUE
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
    # Require complete, non-empty column names so group_var and downstream
    # variable selections can be resolved unambiguously.
    if (is.null(names(x)) || any(!nzchar(names(x)))) {
      stop("`x` must have valid column names.", call. = FALSE)
    }
    
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
  
  if (is.null(colnames(x))) {
    stop("`x` must have column names.", call. = FALSE)
  }
  
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
