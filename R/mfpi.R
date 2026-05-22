#' Model Interactions Between a Categorical and Continuous Covariates
#'
#' `mfpi()` investigates interactions between a categorical variable
#' (`group_var`) and one or more continuous covariates, using fractional
#' polynomial (FP) transformations to capture nonlinear effects. The categorical
#' variable is treated as a factor internally; in a typical randomised controlled
#' trial it represents treatment allocation, with the lowest value as the control
#' arm. Confounders and other prognostic variables are adjusted for via FP
#' transformations selected by the MFP algorithm implemented in the
#' \pkg{mfp2} package.
#'
#' @section The MFPI approach:
#' Modelling interactions between continuous predictors with nonlinear effects
#' and a grouping variable (e.g. treatment) has traditionally relied on
#' categorising the continuous variable and applying standard interaction tests.
#' Categorisation is statistically inefficient and results can be sensitive to
#' the choice of cut-point, which in data-driven approaches introduces serious
#' bias (Royston, Altman, and Sauerbrei 2006).
#'
#' Multivariable Fractional Polynomials for Interaction (MFPI), proposed by
#' Royston and Sauerbrei (2004), avoids categorisation by modelling the
#' relationship between the continuous predictor and the outcome using FP
#' functions within each level of `group_var`. The algorithm extends the
#' standard MFP variable-selection procedure by simultaneously determining
#' functional forms and testing whether those forms differ across groups. For
#' full details see Sauerbrei, Royston, and Zapien (2007).
#'
#' @section Flexibility levels (`flex`):
#' Four levels of flexibility control how FP powers are estimated and
#' constrained between groups. The appropriate level is chosen by the analyst
#' based on subject-matter knowledge and sample-size considerations.
#'
#' **`flex1` (default — least flexible).** FP powers are estimated once from
#' the pooled main-effects model (ignoring interaction), then used unchanged
#' within each group. A likelihood-ratio test comparing the main-effects and
#' interaction models constitutes the interaction test. Because the models are
#' nested and no extra degrees of freedom are consumed by power estimation,
#' simulation studies confirm that the type-I error rate is close to its nominal
#' value (Royston and Sauerbrei 2014). This is the recommended starting point.
#'
#' **`flex2`.** FP powers are estimated jointly across groups (constrained to be
#' equal), then used for both the main-effects and interaction models. The
#' selected powers may differ from those under `flex1`, potentially altering the
#' test outcome. Degrees of freedom are the same as for `flex1`, but no
#' correction is made for the degrees of freedom consumed during power
#' estimation, so p-values are slightly anti-conservative and should be treated
#' as indicative.
#'
#' **`flex3`.** Separate FP powers are estimated for the main-effects model and
#' for the within-group functions, though the within-group powers are still
#' constrained to be equal across groups. Because the main-effects and
#' interaction models may use different FP families, they are non-nested and the
#' likelihood-ratio p-value is again indicative rather than definitive.
#' Simulation evidence suggests `flex3` with FP1 recovers within-group
#' functional forms more accurately than `flex1` when a true interaction is
#' present.
#'
#' **`flex4` (most flexible).** Like `flex3`, but the within-group powers are
#' allowed to differ between groups. The additional degrees of freedom reduce
#' power to detect interaction compared with `flex1`--`flex3`.
#'
#' **Practical guidance.** FP2 functions are more flexible and reduce the risk
#' of misspecification bias, but they are sensitive to influential extreme
#' values and lose power when the true function is monotone. If subject-matter
#' knowledge supports a non-monotone treatment-effect function, use `flex1` with
#' FP2. Otherwise, `flex3` with FP1 is preferred: it includes the linear
#' function as a special case and tends to produce simple, interpretable, and
#' transferable results (Royston and Sauerbrei 2014).
#'
#' To guard against overfitting, consider examining results within a small
#' number of clinically defined subgroups, and seek independent validation
#' before drawing firm conclusions.
#'
#' @section Influential observations:
#' Unlike categorisation, FP-based methods are sensitive to extreme values of
#' the continuous predictor. Before running `mfpi()`, examine the data for
#' influential points and consider whether a small number of extreme values
#' unduly determine the chosen functional form. Trimming or Winsorising a
#' handful of extreme observations is an acceptable pre-processing step and
#' should be reported transparently as part of the initial data analysis
#' (Royston and Sauerbrei 2014).
#'
#' @section Regression families (`family`):
#' `mfpi()` accepts the same family strings as [stats::glm()]. Use
#' `family = "gaussian"` for linear regression, `family = "binomial"` for
#' logistic regression, and `family = "poisson"` for Poisson regression with an
#' optional offset (e.g. log of exposure time). For Cox proportional-hazards
#' models, set `family = "cox"` and supply a [survival::Surv()] object as the
#' response `y`; only right-censored data are currently supported.
#'
#' @section Adjustment variables:
#' MFPI operates in two stages. In the first stage, the MFP algorithm selects
#' adjustment variables and their functional forms from all columns of `x`
#' except the current continuous variable of interest. If no adjustment
#' variables survive selection, the interaction is tested without adjustment. In
#' the second stage, the selected (and possibly transformed) adjustment
#' variables are included as covariates when fitting the main-effects and
#' interaction models for each continuous variable listed in `cont_vars`.
#'
#' To force adjustment variables into the model without selection or
#' transformation, set `select = 1` and `alpha = 0`, or pass variable names via
#' `force_keep` and set `df = 1` for those variables. See [mfp2::mfp2()]
#' for further details.
#'
#' @section Shifting, scaling, and centering:
#' Fractional polynomials require strictly positive input values. `mfpi()`
#' estimates a shift for each variable to ensure positivity, then scales the
#' shifted values to a convenient range, following the same procedure as
#' [mfp2::mfp2()]. Centering is applied after FP powers are estimated.
#' Variables marked with `zero_vars` or `catzero_vars`, and linear terms
#' (`df = 1`), have their shift automatically set to zero because only the
#' positive values are transformed.
#'
#' @section Handling non-positive values:
#' When a continuous predictor has a mixture of zeros (or negative values) and
#' positive values, standard shifting and FP transformation may be
#' inappropriate. Three options are available, matching those in
#' [mfp2::mfp2()]:
#'
#' * **`zero_vars`** — FP transformations are applied only to strictly positive
#'   values; non-positive values are set to zero in the transformed variable.
#' * **`catzero_vars`** — as `zero_vars`, but a binary indicator
#'   (\eqn{Z = 1} if \eqn{x = 0}, \eqn{Z = 0} if \eqn{x > 0}) is
#'   automatically created and included in the adjustment model alongside the
#'   transformed variable.
#' * **`spike_vars`** — triggers the spike-at-zero (SAZ) algorithm, which
#'   formally evaluates whether the binary indicator and the FP function each
#'   contribute explanatory value beyond the other.
#'
#' Note that `spike_vars` and `catzero_vars` are mutually exclusive per
#' variable. Setting `spike_vars` for a variable implicitly sets
#' `catzero_vars`, which in turn implies `zero_vars`.
#'
#' @param x
#'   A numeric matrix of dimension \eqn{n \times p}, where rows are
#'   observations and columns are variables. Column names are required. The
#'   matrix must not contain an intercept column; binary variables should be
#'   coded as 0/1; multi-level categorical variables (other than `group_var`)
#'   should be expanded into dummy variables beforehand. `x` must be free of
#'   missing values.
#' @param y
#'   The response vector (or matrix for survival data). For Gaussian, binomial,
#'   and Poisson families, supply a numeric or factor vector of length \eqn{n}.
#'   For Cox models, supply a two-column [survival::Surv()] object. Must have
#'   the same number of observations as `x`.
#' @param group_var
#'   A single character string naming the categorical (grouping) variable in
#'   `x`. Interactions between `group_var` and each variable in `cont_vars` are
#'   the primary estimands. Internally, values are remapped to 0, 1, 2, ... in
#'   ascending order; the group with the lowest original value is the reference
#'   category. Must have at least two distinct non-missing values.
#' @param cont_vars
#'   A non-empty character vector naming the continuous variables in `x` for
#'   which interactions with `group_var` are to be investigated. All named
#'   variables must exist as columns of `x`.
#' @param flex
#'   A character string controlling how FP powers are estimated and constrained
#'   across groups. One of `"flex1"` (default), `"flex2"`, `"flex3"`, or
#'   `"flex4"`. See the *Flexibility levels* section for details.
#' @param p_interact
#'   Numeric in \eqn{(0, 1]}. Nominal significance level for the interaction
#'   test when `criterion = "pvalue"`. A candidate interaction model is
#'   retained only if its p-value is strictly below `p_interact`. Default is
#'   `0.05`. Ignored when `criterion` is `"aic"` or `"bic"`.
#'
#'   Among the three candidates (linear, FP1, FP2) that clear the threshold,
#'   the one with the **smallest p-value** is chosen as the final functional
#'   form. When two candidates have identical p-values (which can occur due to
#'   numerical precision or rounding), the tie is broken by selecting the
#'   candidate with the **larger `AIC_main_minus_int`** (i.e. the greater improvement in
#'   AIC over the main-effects model), favouring the more parsimonious fit.
#' @param min_improvement
#'   Numeric. Minimum improvement in the selection criterion required to retain
#'   an interaction term. Interpretation depends on `criterion`:
#'   \itemize{
#'     \item `"pvalue"`: not used directly (the threshold is `p_interact`);
#'       `min_improvement` defaults to `p_interact` for consistency.
#'     \item `"aic"`: minimum required `AIC_main_minus_int`
#'       (\eqn{= \mathrm{AIC}_\text{main} - \mathrm{AIC}_\text{int}}).
#'       Must be positive. Typical values: `1` or `2`. Default is `2`.
#'     \item `"bic"`: same as `"aic"` using `BIC_main_minus_int`. Default is `2`.
#'   }
#'   Among candidates that clear `min_improvement`, the one with the
#'   **largest improvement** (largest `AIC_main_minus_int` or `BIC_main_minus_int`) is chosen as
#'   the final functional form.
#' @param include_group_var
#'   Logical. Whether `group_var` should be passed to the MFP algorithm
#'   alongside the other adjustment variables when selecting the adjustment
#'   model in stage one. Default is `FALSE`. When `TRUE`, the nominal
#'   significance level for `group_var` is set to 1, forcing it into the
#'   adjustment model.
#' @param show_models
#'   Logical. Whether to print regression coefficients and standard errors for
#'   the final interaction models. Default is `FALSE`, which prints only the
#'   selected adjustment variables and model-evaluation metrics.
#' @param weights
#'   An optional numeric vector of non-negative observation weights of length
#'   \eqn{n}. Default is `NULL` (all weights equal to 1).
#' @param offset
#'   An optional numeric vector of length \eqn{n} to be added to the linear
#'   predictor. Useful for Poisson models (e.g. log of exposure time). Default
#'   is `NULL` (zero offset for all observations).
#' @param cycles
#'   A positive integer. Maximum number of iteration cycles for the MFP
#'   algorithm. Default is `10`.
#' @param scale
#'   A numeric vector of length \eqn{p} or a single numeric giving scaling
#'   factors for the columns of `x`. Default is `NULL`, which lets the program
#'   estimate scaling factors automatically. Set `scale = 1` to disable scaling.
#' @param shift
#'   A numeric vector of length \eqn{p} or a single numeric giving shift terms
#'   added to columns of `x` before scaling. Default is `NULL`, which estimates
#'   shift factors automatically. Set `shift = 0` to disable shifting.
#' @param df
#'   A numeric vector of length \eqn{p} or a single positive integer setting
#'   the default degrees of freedom for each predictor. Degrees of freedom
#'   equal twice the FP degree (e.g. `df = 2` for FP1, `df = 4` for FP2).
#'   Regardless of the supplied value, the program overrides `df` per variable
#'   based on the number of distinct values \eqn{u}:
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
#'   predictors before fitting the final interaction model. Binary covariates
#'   are centred at the lower of their two values rather than the mean. Default
#'   is `TRUE`.
#' @param subset
#'   An optional integer vector of positive row indices selecting a subset of
#'   observations. Default is `NULL` (all observations used).
#' @param family
#'   A character string specifying the error distribution. One of
#'   `"gaussian"` (default), `"binomial"`, `"poisson"`, or `"cox"`. See the
#'   *Regression families* section for details.
#' @param criterion
#'   A character string specifying the criterion used in two distinct places:
#'   \enumerate{
#'     \item **Adjustment-variable selection** (Step 2): governs which
#'       predictors and FP degrees survive MFP backfitting.
#'     \item **Functional form selection for the interaction** (Step 3): for
#'       each variable in `cont_vars`, three candidate interaction models are
#'       fitted — linear, FP1, and FP2 — and the criterion selects among them.
#'   }
#'   One of `"pvalue"` (default), `"aic"`, or `"bic"`.
#'   The interaction test always reports p-values, AIC differences, and BIC
#'   differences regardless of this setting; `criterion` controls only which
#'   quantity is used to make the final selection.
#' @param select
#'   A numeric vector of length \eqn{p} or a single value in \eqn{[0, 1]}
#'   giving the nominal significance level for backward elimination of each
#'   predictor in the adjustment model. Setting a variable's level to `1`
#'   forces it into the model. Default is `0.05`.
#' @param alpha
#'   A numeric vector of length \eqn{p} or a single value in \eqn{[0, 1]}
#'   giving the significance level for choosing between FP degrees for each
#'   predictor. Default is `0.05`.
#' @param force_keep
#'   An optional character vector of variable names to retain in the adjustment
#'   model regardless of selection criteria. When `criterion = "pvalue"`,
#'   equivalent to setting `select = 1` for those variables; also effective
#'   under AIC and BIC criteria.
#' @param force_max_fp
#'   A logical scalar or named logical vector of length \eqn{p}. If \code{TRUE}
#'   for a variable, forces \code{select_ic()} to select the most complex
#'   functional form at the degree specified by \code{df} for that variable,
#'   bypassing AIC/BIC comparison against simpler forms. Specifically, when
#'   \code{criterion = "aic"} or \code{"bic"}, \code{select_ic()} would
#'   normally compete null, linear, FP1, and FP2 against each other and may
#'   simplify the functional form; \code{force_max_fp = TRUE} suppresses this
#'   and always selects the most complex FP model at the requested degree.
#'   The best power combination within that degree is still selected by the
#'   criterion (equivalently, by deviance minimisation at fixed df). Has no
#'   effect when \code{criterion = "pvalue"} since \code{alpha = 1} already
#'   guarantees the most complex form is accepted. A scalar \code{FALSE}
#'   (default) applies to all variables. A named vector allows per-variable
#'   control.
#' @param xorder
#'   A character string controlling the order in which adjustment covariates
#'   enter the MFP selection algorithm. `"ascending"` (default) enters
#'   variables from most to least significant in a full multiple regression;
#'   `"descending"` reverses this order; `"original"` uses the column order of
#'   `x`.
#' @param fp_powers
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
#'   A numeric vector of values passed to [survival::coxph()] to suppress
#'   centring for specified predictors. Cox models only; ignored otherwise.
#' @param acd_vars
#'   An optional character vector naming continuous variables to be transformed
#'   via the approximate cumulative distribution (ACD) transformation before FP
#'   selection. The transformed variable is named `A(x)`. ACD transformation is
#'   disabled for variables listed in `cont_vars`. See [mfp2::mfp2()] for
#'   details.
#' @param zero_vars
#'   An optional character vector naming variables for which non-positive values
#'   should be treated as zero. FP transformations are applied only to strictly
#'   positive values; non-positive values are set to zero. See the
#'   *Handling non-positive values* section for details.
#' @param catzero_vars
#'   An optional character vector naming variables for which non-positive values
#'   should be treated as zero and a binary indicator variable automatically
#'   created and included in the adjustment model. Implies `zero_vars` for the
#'   named variables. A variable may not appear in both `zero_vars` and
#'   `catzero_vars`.
#' @param spike_vars
#'   An optional character vector naming variables to be assessed for a spike
#'   at zero using the SAZ algorithm. Implies `catzero_vars` (and therefore
#'   `zero_vars`) for the named variables.
#' @param min_prop
#'   Numeric in \eqn{[0, 1]}. Minimum proportion of zeros required for the
#'   SAZ algorithm to be applied. Default is `0.05`. Variables with fewer zeros
#'   are treated as standard continuous predictors.
#' @param max_prop
#'   Numeric in \eqn{[0, 1]}. Maximum proportion of zeros for SAZ modeling.
#'   Default is `0.95`. Variables with more zeros are also treated as standard
#'   continuous predictors.
#' @param use_ftest
#'   Logical. Whether to use an F-test rather than a chi-square test when
#'   computing p-values for Gaussian models. Recommended when the sample size
#'   is small. Default is `FALSE`. Has no effect for non-Gaussian families or
#'   when `criterion` is not `"pvalue"`.
#' @param control
#'   A list of control parameters for the underlying fitting routine, as
#'   returned by [stats::glm.control()] (non-Cox families) or
#'   [survival::coxph.control()] (Cox). Default is `NULL`, which uses the
#'   default control parameters for the chosen family.
#' @param verbose
#'   Logical. Whether to print progress information during model fitting.
#'   Default is `TRUE`.
#' @param digits
#'   A positive integer. Minimum number of significant digits displayed when
#'   printing interaction-test results. Default is `3`.
#' @param data
#'   For `mfpi.formula()`: a data frame containing all variables named in
#'   `formula`. Multi-level categorical variables other than `group_var` must
#'   already be expanded to dummy columns.
#' @param ...
#'   Currently unused. Reserved for future extensions.
#'
#' @return
#' An object of class `"mfpi"`. Use [mfp2::summary.mfpi()] for a formatted
#' summary. The object is a list with the following components:
#'
#' \describe{
#'   \item{\code{model_evaluation_metrics}}{A data frame of evaluation metrics
#'     for the main-effects and interaction models for each variable in
#'     \code{cont_vars}. Columns include: \code{fp_powers_main} and \code{fp_powers_int}
#'     (selected FP powers); \code{deviance_int} (deviance of the interaction
#'     model); \code{deviance_diff} (deviance of main-effects model minus deviance
#'     of interaction model); \code{pvalue} (likelihood-ratio p-value);
#'     \code{AIC_interaction} and \code{BIC_interaction} (information criteria for the
#'     interaction model); \code{AIC_main_minus_int} and \code{BIC_main_minus_int}
#'     (main-effects minus interaction model for each criterion).}
#'   \item{\code{adjust_terms}}{A data frame describing the FP terms selected
#'     for the adjustment model.}
#'   \item{\code{interaction_models}}{A named list of fitted model objects, one
#'     per variable in \code{cont_vars}.}
#'   \item{\code{fitted_functions}}{A named list of data frames, one per
#'     variable in \code{cont_vars}. Each data frame contains the estimated FP
#'     function \eqn{f_j} and its pointwise standard error \eqn{se(f_j)} at
#'     each level \eqn{j} of \code{group_var}, together with 95\% confidence
#'     intervals (\code{fj_lower}, \code{fj_upper}). Columns \code{fj - f0}
#'     and \code{se(fj - f0)} give the difference relative to the reference
#'     group, with corresponding confidence bounds.}
#'   \item{\code{group_var}}{The name of the grouping variable.}
#'   \item{\code{show_models}}{The value of the \code{show_models} argument.}
#'   \item{\code{flex}}{The chosen flexibility level.}
#'   \item{\code{group_levels_new}}{The remapped levels of \code{group_var}
#'     (0, 1, 2, ...).}
#'   \item{\code{group_levels_original}}{The original levels of
#'     \code{group_var}.}
#'   \item{\code{family}}{The regression family, as a character string.}
#'   \item{\code{nobs}}{The number of observations used in model fitting.}
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
#' Sauerbrei, W. and Royston, P. (2022). Investigating treatment-effect
#' modification by a continuous covariate in IPD meta-analysis: an approach
#' using fractional polynomials. \emph{BMC Medical Research Methodology}, 22,
#' 1--13.
#'
#' Royston, P. and Sauerbrei, W. (2008). \emph{Multivariable Model-Building: A
#' Pragmatic Approach to Regression Analysis Based on Fractional Polynomials
#' for Modelling Continuous Variables}. John Wiley & Sons.
#'
#' @seealso [mfp2::summary.mfpi()], [mfp2::mfp2()]
#'
#' @export
mfpi <- function(x, ...) {
  UseMethod("mfpi", x)
}

#' @describeIn mfpi Default method accepting a numeric matrix \code{x} and
#'   response vector \code{y}.
#' @import mfp2
#' @export
mfpi.default <- function(
    x,
    y,
    group_var         = NULL,
    cont_vars         = NULL,
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
    family            = c("gaussian", "poisson", "binomial", "cox"),
    criterion         = c("pvalue", "aic", "bic"),
    select            = 0.05,
    alpha             = 0.05,
    force_keep        = NULL,
    force_max_fp      = FALSE,
    xorder            = c("ascending", "descending", "original"),
    fp_powers         = NULL,
    ties              = c("breslow", "efron", "exact"),
    strata            = NULL,
    nocenter          = NULL,
    acd_vars          = NULL,
    zero_vars         = NULL,
    catzero_vars      = NULL,
    spike_vars        = NULL,
    min_prop          = 0.05,
    max_prop          = 0.95,
    use_ftest         = FALSE,
    control           = NULL,
    verbose           = TRUE,
    digits            = 3,
    ...
) {
  cl <- match.call()
  
  # Match enumerated arguments -------------------------------------------------
  family    <- match.arg(family)
  xorder    <- match.arg(xorder)
  criterion <- match.arg(criterion)
  flex      <- match.arg(flex)
  ties      <- match.arg(ties)
  
  # Resolve family: convert string to GLM object for non-Cox families ----------
  # mfp2:::fit_mfp() and mfp2:::fit_model() expect a GLM family object for
  # non-Cox families and the character string "cox" for Cox models.
  allowed_families <- c("gaussian", "binomial", "poisson", "cox")
  if (!family %in% allowed_families) {
    stop(
      paste0("! Invalid family: '", family, "'. Allowed values are: ",
             paste(allowed_families, collapse = ", "), "."),
      call. = FALSE
    )
  }
  family_string <- family   # keep string for branching (== "cox", == "gaussian")
  if (family != "cox") {
    family <- tryCatch(
      get(family, mode = "function", envir = parent.frame())(),
      error = function(e) stop(
        paste0("! Cannot construct GLM family object from '", family_string,
               "'. ", conditionMessage(e)),
        call. = FALSE
      )
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
  nlev <- length(unique(x[, group_var]))
  if (nlev < 2L) {
    stop(
      paste0("! `group_var` must have at least two distinct values; found ", nlev, "."),
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
  
  # Validate response y --------------------------------------------------------
  if (family_string == "cox") {
    if (!survival::is.Surv(y)) {
      stop("! For `family = 'cox'`, `y` must be a `survival::Surv` object.",
           call. = FALSE)
    }
    if (nrow(y) != nobs) {
      stop(paste0("! `y` has ", nrow(y), " rows but `x` has ", nobs,
                  " rows; they must match."), call. = FALSE)
    }
    if (attr(y, "type") != "right") {
      stop(paste0("! Only right-censored survival data are currently supported; ",
                  "`y` has censoring type '", attr(y, "type"), "'."),
           call. = FALSE)
    }
    if (is.factor(strata)) strata <- as.numeric(strata)
    if (!is.null(strata)) {
      strata_len <- if (is.vector(strata)) length(strata) else nrow(strata)
      if (strata_len != nobs) {
        stop(
          paste0("! Length of `strata` (", strata_len, ") must equal the number of ",
                 "rows in `x` (", nobs, ")."),
          call. = FALSE
        )
      }
    }
  } else {
    if (is.matrix(y) || is.data.frame(y)) {
      stop(paste0("! For `family = '", family_string, "'`, `y` must be a vector, ",
                  "not an object of class ", paste(class(y), collapse = ", "), "."),
           call. = FALSE)
    }
    if (length(y) != nobs) {
      stop(paste0("! `y` has length ", length(y), " but `x` has ", nobs,
                  " rows; they must match."), call. = FALSE)
    }
  }
  
  # Validate weights and offset ------------------------------------------------
  if (!is.null(weights)) {
    if (any(weights < 0))
      stop("! `weights` must not be negative.", call. = FALSE)
    if (length(weights) != nobs)
      stop(paste0("! Length of `weights` (", length(weights), ") must equal ",
                  "the number of rows in `x` (", nobs, ")."), call. = FALSE)
  }
  if (!is.null(offset) && length(offset) != nobs) {
    stop(paste0("! Length of `offset` (", length(offset), ") must equal ",
                "the number of rows in `x` (", nobs, ")."), call. = FALSE)
  }
  
  # Validate alpha and select --------------------------------------------------
  if (any(alpha < 0) || any(alpha > 1))
    stop("! All values of `alpha` must be in [0, 1].", call. = FALSE)
  if (length(alpha) != 1L && length(alpha) != nvars)
    stop(paste0("! `alpha` must be a single number or a vector of length ", nvars,
                "; got length ", length(alpha), "."), call. = FALSE)
  
  if (any(select < 0) || any(select > 1))
    stop("! All values of `select` must be in [0, 1].", call. = FALSE)
  if (length(select) != 1L && length(select) != nvars)
    stop(paste0("! `select` must be a single number or a vector of length ", nvars,
                "; got length ", length(select), "."), call. = FALSE)
  
  # Validate force_keep --------------------------------------------------------
  if (!is.null(force_keep) && !all(force_keep %in% vnames)) {
    warning(
      "i Some variables in `force_keep` are not columns of `x`; ",
      "continuing with the intersection.",
      call. = FALSE
    )
  }
  
  # Validate force_max_fp ------------------------------------------------------
  if (!is.logical(force_max_fp)) {
    stop("! `force_max_fp` must be logical.", call. = FALSE)
  }
  if (length(force_max_fp) == 1L) {
    force_max_fp <- setNames(rep(force_max_fp, nvars), vnames)
  } else if (length(force_max_fp) != nvars) {
    stop(paste0("! `force_max_fp` must be a single logical or a named logical ",
                "vector of length ", nvars, " (ncol(x)); got length ",
                length(force_max_fp), "."),
         call. = FALSE)
  } else {
    force_max_fp <- setNames(as.logical(force_max_fp), vnames)
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
  if (!is.null(zero_vars) && !is.null(catzero_vars)) {
    overlap <- intersect(zero_vars, catzero_vars)
    if (length(overlap) > 0L) {
      stop(paste0("! The following variables appear in both `zero_vars` and `catzero_vars`: ",
                  paste(overlap, collapse = ", "),
                  ". Use only one option per variable."), call. = FALSE)
    }
  }
  
  # Validate acd_vars ----------------------------------------------------------
  if (!is.null(acd_vars) && !all(acd_vars %in% vnames)) {
    warning("i Some variables in `acd_vars` are not columns of `x`; they will be ignored.",
            call. = FALSE)
  }
  
  # Validate shift and scale ---------------------------------------------------
  if (!is.null(shift) && length(shift) != 1L && length(shift) != nvars) {
    stop(paste0("! `shift` must be NULL, a single number, or a vector of length ",
                nvars, "; got length ", length(shift), "."), call. = FALSE)
  }
  if (!is.null(scale) && length(scale) != 1L && length(scale) != nvars) {
    stop(paste0("! `scale` must be NULL, a single number, or a vector of length ",
                nvars, "; got length ", length(scale), "."), call. = FALSE)
  }
  if (!is.null(scale) && any(scale <= 0 | is.na(scale))) {
    stop("! All values of `scale` must be positive and non-missing.", call. = FALSE)
  }
  
  # Validate center ------------------------------------------------------------
  if (length(center) != 1L && length(center) != nvars) {
    stop(paste0("! `center` must be a single logical or a logical vector of length ",
                nvars, "; got length ", length(center), "."), call. = FALSE)
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
      min_prop < 0 || min_prop > 1)
    stop("! `min_prop` must be a single numeric value in [0, 1].", call. = FALSE)
  if (!is.numeric(max_prop) || length(max_prop) != 1L ||
      max_prop < 0 || max_prop > 1)
    stop("! `max_prop` must be a single numeric value in [0, 1].", call. = FALSE)
  if (min_prop > max_prop)
    stop("! `min_prop` cannot be greater than `max_prop`.", call. = FALSE)
  
  # Warn if use_ftest is incompatible with family ------------------------------
  if (use_ftest && family_string != "gaussian") {
    warning(
      paste0("i `use_ftest = TRUE` is only applicable to Gaussian models; ",
             "reverting to chi-square test for family = '", family_string, "'."),
      call. = FALSE
    )
    use_ftest <- FALSE
  }
  
  # Validate and build fp_powers list ------------------------------------------
  if (!is.null(fp_powers) && !is.list(fp_powers))
    stop("! `fp_powers` must be a named list.", call. = FALSE)
  
  default_powers <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
  power_list <- setNames(replicate(nvars, default_powers, simplify = FALSE), vnames)
  
  if (!is.null(fp_powers)) {
    if (length(fp_powers) != sum(nchar(names(fp_powers)) > 0L, na.rm = TRUE))
      stop("! Every element of `fp_powers` must have a name.", call. = FALSE)
    unknown_names <- setdiff(names(fp_powers), vnames)
    if (length(unknown_names) > 0L)
      stop(paste0("! The following names in `fp_powers` do not match any column of `x`: ",
                  paste(unknown_names, collapse = ", "), "."), call. = FALSE)
    if (!all(sapply(fp_powers, is.numeric)))
      stop("! All elements of `fp_powers` must be numeric vectors.", call. = FALSE)
    fp_powers  <- lapply(fp_powers, sort)
    pw_lengths <- sapply(fp_powers, length)
    too_short  <- names(pw_lengths[pw_lengths < 2L])
    if (length(too_short) > 0L)
      stop(paste0("! Each element of `fp_powers` must contain at least two values. ",
                  "Insufficient values for: ", paste(too_short, collapse = ", "), "."),
           call. = FALSE)
    power_list <- modifyList(power_list, Filter(Negate(is.null), fp_powers))
  }
  
  # Validate subset ------------------------------------------------------------
  if (!is.null(subset)) {
    if (!is.vector(subset))
      stop(paste0("! `subset` must be a vector; got class ",
                  paste(class(subset), collapse = ", "), "."), call. = FALSE)
    if (any(subset < 0L))
      stop("! `subset` must not contain negative indices.", call. = FALSE)
    if (length(subset) < 5L)
      stop(paste0("! After subsetting, only ", length(subset),
                  " observations remain; at least 5 are required."), call. = FALSE)
  }
  
  # Set scalar/vector defaults -------------------------------------------------
  if (is.null(min_improvement)) {
    min_improvement <- switch(criterion, pvalue = p_interact, aic = 2, bic = 2)
  }
  if (is.null(weights)) weights <- rep.int(1, nobs)
  if (is.null(offset))  offset  <- rep.int(0, nobs)
  
  # Expand scalars and guarantee names on all per-variable vectors.
  # mfp2:::fit_mfp() and preprocess_data() rely on named vectors for
  # correct alignment after group_var is removed.
  if (length(select) == 1L) select <- rep(select, nvars)
  if (length(alpha)  == 1L) alpha  <- rep(alpha,  nvars)
  if (length(center) == 1L) center <- rep(center, nvars)
  select <- setNames(select, vnames)
  alpha  <- setNames(alpha,  vnames)
  center <- setNames(center, vnames)
  
  if (is.null(shift)) {
    shift <- apply(x, 2L, mfp2::find_shift_factor)   # already named by apply
  } else if (length(shift) == 1L) {
    shift <- rep(shift, nvars)
  }
  shift <- setNames(shift, vnames)
  
  if (is.null(scale)) {
    scale <- apply(x, 2L, mfp2::find_scale_factor)   # already named by apply
  } else if (length(scale) == 1L) {
    scale <- rep(scale, nvars)
  }
  scale <- setNames(scale, vnames)
  if (is.null(control)) {
    control <- if (family_string == "cox") survival::coxph.control() else
      stats::glm.control()
  }
  
  # Retain only valid force_keep entries ---------------------------------------
  force_keep <- intersect(force_keep, vnames)
  
  # Process zero_vars / catzero_vars / spike_vars to named logical vectors -----
  zero_vars    <- intersect(zero_vars,    vnames)
  catzero_vars <- intersect(catzero_vars, vnames)
  spike_vars   <- intersect(spike_vars,   vnames)
  
  zero_flag    <- setNames(rep(FALSE, nvars), vnames)
  catzero_flag <- setNames(rep(FALSE, nvars), vnames)
  spike_flag   <- setNames(rep(FALSE, nvars), vnames)
  
  if (length(zero_vars)    > 0L) zero_flag[zero_vars]       <- TRUE
  if (length(catzero_vars) > 0L) catzero_flag[catzero_vars] <- TRUE
  if (length(spike_vars)   > 0L) {
    spike_flag[spike_vars] <- TRUE
    catzero_flag[spike_vars] <- TRUE   # spike implies catzero
  }
  catzero_flag[spike_flag]  <- TRUE   # redundant but explicit
  zero_flag[catzero_flag]   <- TRUE   # catzero implies zero
  
  # Reset zero/catzero for variables that contain only positive values ---------
  if (any(zero_flag)) {
    vars_pos <- names(zero_flag)[zero_flag]
    all_positive <- vars_pos[
      apply(x[, vars_pos, drop = FALSE], 2L,
            function(col) all(col > 0, na.rm = TRUE))
    ]
    if (length(all_positive) > 0L) {
      warning(
        paste0("i The following variables in `zero_vars` / `catzero_vars` contain ",
               "only positive values and have been reset to standard processing: ",
               paste(all_positive, collapse = ", "), "."),
        call. = FALSE
      )
      zero_flag[all_positive]    <- FALSE
      catzero_flag[all_positive] <- FALSE
    }
  }
  
  # Reset zero/catzero for binary variables ------------------------------------
  binary_vars  <- apply(x, 2L,
                        function(col) length(unique(col[!is.na(col)])) == 2L)
  binary_names <- vnames[binary_vars]
  if (length(binary_names) > 0L) {
    reset_vars <- intersect(
      binary_names,
      union(names(zero_flag)[zero_flag], names(catzero_flag)[catzero_flag])
    )
    if (length(reset_vars) > 0L) {
      warning(
        paste0("i The following binary variables were marked in `zero_vars` or ",
               "`catzero_vars` but are binary; resetting to standard processing: ",
               paste(reset_vars, collapse = ", "), "."),
        call. = FALSE
      )
      zero_flag[reset_vars]    <- FALSE
      catzero_flag[reset_vars] <- FALSE
    }
  }
  
  # Process acd_vars to a named logical vector ---------------------------------
  acd_flag <- setNames(rep(FALSE, nvars), vnames)
  if (!is.null(acd_vars)) {
    acd_vars <- unique(intersect(acd_vars, vnames))
    # ACD transformation is not supported for cont_vars
    acd_vars <- setdiff(acd_vars, cont_vars)
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
    y       <- if (family_string == "cox") y[subset, , drop = FALSE] else y[subset]
    weights <- weights[subset]
    offset  <- offset[subset]
    if (!is.null(istrata)) istrata <- istrata[subset]
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
    force_keep        = force_keep,
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
    use_ftest         = use_ftest,
    control           = control,
    verbose           = verbose,
    group_var         = group_var,
    include_group_var = include_group_var,
    flex              = flex,
    cont_vars         = cont_vars,
    p_interact        = p_interact,
    show_models       = show_models,
    min_improvement   = min_improvement,
    digits            = digits
  )
  
  class(fit) <- "mfpi"
  fit
}