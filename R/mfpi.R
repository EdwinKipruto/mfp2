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
#'   candidate with the **larger `AIC_diff`** (i.e. the greater improvement in
#'   AIC over the main-effects model), favouring the more parsimonious fit.
#' @param min_improvement
#'   Numeric. Minimum improvement in the selection criterion required to retain
#'   an interaction term. Interpretation depends on `criterion`:
#'   \itemize{
#'     \item `"pvalue"`: not used directly (the threshold is `p_interact`);
#'       `min_improvement` defaults to `p_interact` for consistency.
#'     \item `"aic"`: minimum required `AIC_diff`
#'       (\eqn{= \mathrm{AIC}_\text{main} - \mathrm{AIC}_\text{int}}).
#'       Must be positive. Typical values: `1` or `2`. Default is `2`.
#'     \item `"bic"`: same as `"aic"` using `BIC_diff`. Default is `2`.
#'   }
#'   Among candidates that clear `min_improvement`, the one with the
#'   **largest improvement** (largest `AIC_diff` or `BIC_diff`) is chosen as
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
#'   The program overrides user-supplied values based on the number of unique
#'   values: variables with 2--3 unique values are set to `df = 1` (linear);
#'   4--5 unique values to `df = min(2, default)`; and 6 or more unique values
#'   retain the default. Default is `4`.
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
#' An object of class `"mfpi"`. Use `summary.mfpi()` for a formatted
#' summary. The object is a list with the following components:
#'
#' \describe{
#'   \item{\code{model_evaluation_metrics}}{A data frame of evaluation metrics
#'     for the main-effects and interaction models for each variable in
#'     \code{cont_vars}. Columns include: \code{pow_main} and \code{pow_int}
#'     (selected FP powers); \code{dev_int} (deviance of the interaction
#'     model); \code{dev_diff} (deviance of main-effects model minus deviance
#'     of interaction model); \code{pvalue} (likelihood-ratio p-value);
#'     \code{AIC_int} and \code{BIC_int} (information criteria for the
#'     interaction model); \code{AIC_diff} and \code{BIC_diff}
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
#' @seealso `summary.mfpi()`, [mfp2::mfp2()]
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
  if (length(select) == 1L) select <- rep(select, nvars)
  if (length(alpha)  == 1L) alpha  <- rep(alpha,  nvars)
  if (length(center) == 1L) center <- setNames(rep(center, nvars), vnames)
  
  if (is.null(shift)) {
    shift <- apply(x, 2L, mfp2::find_shift_factor)
  } else if (length(shift) == 1L) {
    shift <- rep(shift, nvars)
  }
  if (is.null(scale)) {
    scale <- apply(x, 2L, mfp2::find_scale_factor)
  } else if (length(scale) == 1L) {
    scale <- rep(scale, nvars)
  }
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
  if (length(df) == 1L) {
    df_list <- if (df != 1L) mfp2::assign_df(x = x, df_default = df) else
      rep(df, nvars)
  } else {
    nux   <- apply(x, 2L, function(v) length(unique(v)))
    small <- nux <= 3L & df != 1L
    if (any(small)) {
      warning(
        paste0("i `df` has been overridden to 1 (linear) for variables with fewer ",
               "than 4 unique values: ",
               paste(vnames[small], collapse = ", "), "."),
        call. = FALSE
      )
      df[small] <- 1L
    }
    df_list <- df
  }
  
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

#' @describeIn mfpi Provides formula interface for `mfpi`.
#' @export
mfpi.formula <- function(formula,
                         data,
                         catvar = NULL,
                         adjcatvar = FALSE,
                         flex = c("flex1", "flex2", "flex3", "flex4"),
                         linearvar = NULL,
                         fp1var = NULL,
                         fp2var = NULL,
                         showmodel = FALSE,
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
                         data_input = c("original", "equidistant"),
                         strata = NULL,
                         nocenter = NULL,
                         ftest = FALSE,
                         control = NULL,
                         verbose = TRUE,
                         digits = 3,
                         ...) {
  # capture the call
  call <- match.call()
  family <- match.arg(family)
  xorder <- match.arg(xorder)
  criterion <- match.arg(criterion)
  flex <- match.arg(flex)
  ties <- match.arg(ties)
  data_input <- match.arg(data_input)
  
  # assert that data must be provided
  if (missing(data))
    stop("! data argument is missing.\n",
         "i An input data.frame is required for the use of mfpi.",
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
  
  # Check if at least one of the following arguments is not NULL: linear, fp1var, or fp2var
  if (is.null(linearvar) && is.null(fp1var) && is.null(fp2var)) {
    stop("i At least one of the following arguments must not be NULL: linear, fp1var, and fp2var.",
         call. = FALSE)
  }
  
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
  
  # stratification for Cox models : TO CHECK THIS PART
  
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
  if (family != "cox" & survival::is.Surv(y)){
    stop(sprintf("The survival object (y) has an unexpected family '%s' instead of the expected 'cox'.",
                 family), call. = FALSE)
  }
  
  if (family != "cox"){
    y <- as.numeric(y)
  }
  
  
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
  df_list <-  setNames(as.list(mfp2::assign_df(x = x, df_default = df)), names_x)
  
  # scaling
  if (is.null(scale)) {
    scale_list <- stats::setNames(as.list(apply(x, 2, mfp2::find_scale_factor)), names_x)
  } else {
    scale_list <- setNames(rep(list(scale), nx), names_x)
    
  }
  
  # shifting
  if (is.null(shift)) {
    shift_list <- setNames(as.list(apply(x, 2, mfp2::find_shift_factor)), names_x)
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
      warning("i Powers are specified in both the `fp()` function within\n the formula and as an argument. The argument term ignored.\n",
              sprintf("i This applies to the following variables: %s.",
                      paste0(nax, collapse = ", ")), call. = FALSE)
    
    power_list <- modifyList(power_list, powerx)
    
  }
  
  # acd requires variable names or NULL
  acdx_vector <- unlist(acdx_list)
  if (sum(acdx_vector) == 0){
    acdx_vector <- NULL
  }else{
    acdx_vector <- names(acdx_vector[acdx_vector])
  }
  
  mfpi.default(x = x, y = y,
               catvar = catvar,
               adjcatvar = adjcatvar,
               flex = flex,
               linearvar = linearvar,
               fp1var = fp1var,
               fp2var = fp2var,
               showmodel = showmodel,
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
               verbose = verbose,
               digits = digits,
               data_input = data_input
  )
}


# -----------------------------------------------------------------------------
# print.mfpi() ----------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Print an \code{mfpi} Object
#'
#' Displays a structured summary of an \code{"mfpi"} object, including the
#' adjustment variables selected by the MFP algorithm and the interaction test
#' results for each continuous variable in \code{cont_vars}. Optionally prints
#' the full regression output for each final interaction model when
#' \code{show_models = TRUE} was set in \code{mfpi()}.
#'
#' @param x An object of class \code{"mfpi"}, as returned by \code{mfpi()}.
#' @param ... Currently unused.
#'
#' @return \code{x} invisibly.
#'
#' @section Output sections:
#' \enumerate{
#'   \item **Adjustment model** — the \code{fp_terms} table from the MFP
#'     algorithm, showing which variables were selected and their estimated FP
#'     powers.
#'   \item **Interaction test** — for each variable in \code{cont_vars}, the
#'     best-selected functional form (linear, FP1, or FP2) and its evaluation
#'     metrics: deviance, degrees of freedom, p-value, AIC, and BIC. Columns:
#'     \itemize{
#'       \item \code{type}: selected functional form (\code{linear}, \code{fp1},
#'         or \code{fp2}).
#'       \item \code{pow_main}: FP powers in the main-effects model.
#'       \item \code{pow_int}: FP powers in the interaction model, one set per
#'         group. For \code{flex4} these may differ across groups.
#'       \item \code{dev_diff}: likelihood-ratio test statistic
#'         \eqn{T = -2\ell_\text{main} - (-2\ell_\text{int})}.
#'       \item \code{df}: degrees of freedom for the LRT
#'         (\eqn{= (K-1)m} for flex1/2/3; \eqn{= 2(K-1)m} for flex4).
#'       \item \code{pvalue}: \eqn{\Pr[\chi^2(df) > T]}.
#'       \item \code{AIC_diff}, \code{BIC_diff}: improvement of interaction
#'         model over main-effects model; positive values favour interaction.
#'     }
#'   \item **Interaction model coefficients** (only when
#'     \code{show_models = TRUE}) — the regression summary for the final
#'     interaction model for each variable.
#' }
#'
#' @method print mfpi
#' @export
print.mfpi <- function(x, ...) {
  
  # ---------------------------------------------------------------------------
  # Section 1: Adjustment model
  # ---------------------------------------------------------------------------
  cat("Adjustment Variables Selected by MFP Algorithm:\n")
  cat(strrep("-", 65), "\n")
  
  fp <- x$adjust_terms
  if (!is.null(fp)) {
    # Show only user-relevant columns; suppress mfp2 bookkeeping columns.
    keep_cols <- intersect(
      c("selected", "df_final", "power1", "power2"),
      colnames(fp)
    )
    # Split into selected and dropped for clarity
    selected_rows <- fp[!is.na(fp$selected) & fp$selected, keep_cols, drop = FALSE]
    dropped_rows  <- fp[!is.na(fp$selected) & !fp$selected, keep_cols, drop = FALSE]
    
    if (nrow(selected_rows) > 0L) {
      cat("Selected:\n")
      print(selected_rows)
    }
    if (nrow(dropped_rows) > 0L) {
      cat(sprintf("\nDropped (%d variables eliminated by MFP):", nrow(dropped_rows)))
      cat("", paste(rownames(dropped_rows), collapse = ", "), "\n")
    }
    if (nrow(selected_rows) == 0L) {
      cat("No adjustment variables were selected.\n")
    }
  } else {
    cat("No adjustment model fitted.\n")
  }
  
  # ---------------------------------------------------------------------------
  # Section 2: Interaction test results
  # ---------------------------------------------------------------------------
  cat("\n", strrep("-", 65), "\n", sep = "")
  cat("Interaction Test\n")
  cat(strrep("-", 65), "\n")
  cat(sprintf(
    "\nInteraction with '%s'  |  %d observations  |  %s strategy\n\n",
    x$group_var, x$nobs, x$flex
  ))
  
  print(x$model_evaluation_metrics)
  
  cat(
    "\nNote: df = LRT degrees of freedom for interaction;",
    "tdf = total df in interaction model.\n"
  )
  
  # For flex4, each group has its own FP powers so pow_int shows K power sets
  if (identical(x$flex, "flex4") && !is.null(x$group_levels_original)) {
    cat(sprintf(
      "\nFor flex4, 'pow_int' shows separate FP powers for each level of '%s': %s.\n",
      x$group_var,
      paste(x$group_levels_original, collapse = ", ")
    ))
  }
  
  # ---------------------------------------------------------------------------
  # Section 3: Interaction model summaries (optional)
  # ---------------------------------------------------------------------------
  if (isTRUE(x$show_models) && length(x$interaction_models) > 0L) {
    
    cat("\n", strrep("-", 65), "\n", sep = "")
    cat("Regression Output for Final Interaction Models:\n")
    cat(strrep("-", 65), "\n")
    
    # model_evaluation_metrics has one row per variable with a 'type' column
    # (the selected functional form) — use it to label the output
    metrics <- x$model_evaluation_metrics
    type_lookup <- if (!is.null(metrics) && "type" %in% names(metrics)) {
      setNames(
        toupper(gsub("fp", "FP", metrics$type)),  # "linear"->"LINEAR" / "fp1"->"FP1"
        metrics$variable
      )
    } else {
      setNames(rep("", length(x$interaction_models)),
               names(x$interaction_models))
    }
    
    for (var_name in names(x$interaction_models)) {
      fit_obj  <- x$interaction_models[[var_name]]
      form_lbl <- if (!is.null(type_lookup[[var_name]])) type_lookup[[var_name]] else ""
      
      cat(sprintf(
        "\n'%s' x '%s'  [%s]\n",
        x$group_var, var_name, form_lbl
      ))
      cat(strrep("-", 45), "\n")
      
      # fit_obj is the object returned by mfp2:::fit_model (fast = FALSE).
      # It carries a $fit component holding the raw glm / coxph object.
      if (!is.null(fit_obj$fit)) {
        print(summary(fit_obj$fit))
      } else {
        # Fallback: print whatever is available
        print(fit_obj)
      }
    }
  }
  
  invisible(x)
}

# S3 summary method for mfpi objects
#
# summary.mfpi() returns a structured list containing the key results from an
# mfpi fit. Unlike print.mfpi(), which formats output for the console,
# summary.mfpi() returns an object that can be inspected programmatically.
# A print method for the returned object formats it for display.
#
# Note on degrees of freedom: the p-values in the regression summary tables
# produced by summary.glm() and summary.coxph() do NOT account for the
# degrees of freedom consumed by estimating FP powers. They should be treated
# as approximate. The correct interaction test p-values, with proper df, are
# in the $model_evaluation_metrics component.


# -----------------------------------------------------------------------------
# summary.mfpi() --------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Summarise an \code{"mfpi"} Object
#'
#' Produces a structured summary of an \code{"mfpi"} model fit, including the
#' adjustment model, interaction test metrics, and the regression output for
#' each final interaction model. The returned object has its own print method
#' for formatted console display.
#'
#' @param object An object of class \code{"mfpi"}, as returned by
#'   \code{mfpi()}.
#' @param ... Currently unused.
#'
#' @return An object of class \code{"summary.mfpi"}, which is a named list
#'   with the following components:
#' \describe{
#'   \item{\code{group_var}}{Name of the grouping variable.}
#'   \item{\code{nobs}}{Number of observations.}
#'   \item{\code{family}}{Regression family (GLM object or \code{"cox"}).}
#'   \item{\code{flex}}{Flexibility level used (\code{"flex1"} to
#'     \code{"flex4"}).}
#'   \item{\code{group_levels_original}}{Original levels of \code{group_var}.}
#'   \item{\code{adjust_terms}}{The full \code{fp_terms} data frame from the
#'     MFP adjustment model, as returned by \code{mfp2:::fit_mfp()}.}
#'   \item{\code{model_evaluation_metrics}}{Tibble of interaction test results
#'     (one row per significant variable): functional form selected, FP powers,
#'     deviance, df, p-value, AIC, and BIC. These p-values correctly account
#'     for the degrees of freedom of the interaction test.}
#'   \item{\code{all_model_metrics}}{Tibble of metrics for all three candidate
#'     models (linear, FP1, FP2) for every variable tested, regardless of
#'     significance.}
#'   \item{\code{model_summaries}}{Named list of regression summaries, one per
#'     variable in \code{cont_vars} that had a significant interaction. Each
#'     element is the output of \code{summary(fit$fit)}, where \code{fit} is
#'     the raw \code{glm} or \code{coxph} object. \strong{Note:} p-values in
#'     these summaries are from the standard \code{glm}/\code{coxph} machinery
#'     and do not account for the degrees of freedom used to estimate FP
#'     powers; use \code{model_evaluation_metrics$pvalue} for the correct
#'     interaction test p-values.}
#' }
#'
#' @seealso \code{print.summary.mfpi()}, \code{mfpi()},
#'   [stats::summary.glm()], [survival::summary.coxph()]
#'
#' @method summary mfpi
#' @export
summary.mfpi <- function(object, ...) {
  
  # Summarise each final interaction model (flat named list by var_name)
  model_summaries <- lapply(object$interaction_models, function(fit_obj) {
    if (!is.null(fit_obj$fit)) {
      summary(fit_obj$fit)
    } else {
      NULL
    }
  })
  
  structure(
    list(
      group_var                = object$group_var,
      nobs                     = object$nobs,
      family                   = object$family,
      flex                     = object$flex,
      group_levels_original    = object$group_levels_original,
      adjust_terms             = object$adjust_terms,
      model_evaluation_metrics = object$model_evaluation_metrics,
      all_model_metrics        = object$all_model_metrics,
      model_summaries          = model_summaries
    ),
    class = "summary.mfpi"
  )
}


# -----------------------------------------------------------------------------
# print.summary.mfpi() --------------------------------------------------------
# -----------------------------------------------------------------------------

#' Print a \code{"summary.mfpi"} Object
#'
#' Formats and prints the structured summary produced by
#' \code{summary.mfpi()}.
#'
#' @param x An object of class \code{"summary.mfpi"}.
#' @param ... Currently unused.
#'
#' @return \code{x} invisibly.
#'
#' @method print summary.mfpi
#' @export
print.summary.mfpi <- function(x, ...) {
  
  cat(strrep("=", 65), "\n")
  cat(sprintf(
    "MFPI Summary  |  group: '%s'  |  n = %d  |  %s\n",
    x$group_var, x$nobs, x$flex
  ))
  cat(strrep("=", 65), "\n")
  
  # --- Adjustment model ------------------------------------------------------
  cat("\nAdjustment Variables (MFP):\n")
  cat(strrep("-", 65), "\n")
  fp <- x$adjust_terms
  if (!is.null(fp)) {
    keep_cols     <- intersect(c("selected", "df_final", "power1", "power2"),
                               colnames(fp))
    selected_rows <- fp[!is.na(fp$selected) &  fp$selected, keep_cols, drop = FALSE]
    dropped_rows  <- fp[!is.na(fp$selected) & !fp$selected, keep_cols, drop = FALSE]
    
    if (nrow(selected_rows) > 0L) {
      print(selected_rows)
    } else {
      cat("  No adjustment variables selected.\n")
    }
    if (nrow(dropped_rows) > 0L) {
      cat(sprintf(
        "\n  Dropped (%d): %s\n",
        nrow(dropped_rows),
        paste(rownames(dropped_rows), collapse = ", ")
      ))
    }
  } else {
    cat("  No adjustment model fitted.\n")
  }
  
  # --- Interaction test results ----------------------------------------------
  cat("\n", strrep("-", 65), "\n", sep = "")
  cat("Interaction Test Results:\n")
  cat(strrep("-", 65), "\n\n")
  print(x$model_evaluation_metrics)
  cat("\nNote: df = LRT df for interaction; tdf = total df in interaction model.\n")
  cat("      p-values correctly account for interaction test degrees of freedom.\n")
  
  if (identical(x$flex, "flex4") && !is.null(x$group_levels_original)) {
    cat(sprintf(
      "\nFor flex4, 'pow_int' shows per-group FP powers (levels: %s).\n",
      paste(x$group_levels_original, collapse = ", ")
    ))
  }
  
  # --- Regression summaries --------------------------------------------------
  if (length(x$model_summaries) > 0L) {
    cat("\n", strrep("-", 65), "\n", sep = "")
    cat("Regression Output for Final Interaction Models:\n")
    cat("(Note: p-values below are from glm/coxph and do not account for\n")
    cat(" FP power estimation df. Use 'model_evaluation_metrics$pvalue'\n")
    cat(" for the correct interaction test p-values.)\n")
    cat(strrep("-", 65), "\n")
    
    metrics <- x$model_evaluation_metrics
    type_lookup <- if (!is.null(metrics) && "type" %in% names(metrics)) {
      setNames(toupper(gsub("fp", "FP", metrics$type)), metrics$variable)
    } else {
      setNames(rep("", length(x$model_summaries)), names(x$model_summaries))
    }
    
    for (var_name in names(x$model_summaries)) {
      sm       <- x$model_summaries[[var_name]]
      form_lbl <- if (!is.null(type_lookup[[var_name]])) type_lookup[[var_name]] else ""
      cat(sprintf("\n'%s' x '%s'  [%s]\n", x$group_var, var_name, form_lbl))
      cat(strrep("-", 45), "\n")
      if (!is.null(sm)) print(sm) else cat("  (model summary unavailable)\n")
    }
  }
  
  invisible(x)
}
