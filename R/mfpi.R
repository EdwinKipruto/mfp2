#' Model Interactions Between a Categorical and Continuous Covariates
#'
#' `mfpi()` investigates interactions between a categorical variable
#' (`group_var`) and one or more continuous covariates, using fractional
#' polynomial (FP) transformations to capture nonlinear effects. The categorical
#' variable is treated as a factor internally; in a typical randomised controlled
#' trial it represents treatment allocation, with the lowest value as the control
#' arm. Confounders and other prognostic variables are adjusted for via FP
#' transformations selected by [mfp2::mfp2()].
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
#' power to detect interaction compared with `flex1`–`flex3`.
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
#' MFPI operates in two stages. In the first stage, [mfp2::mfp2()] selects
#' adjustment variables and their functional forms from all columns of `x`
#' except the current continuous variable of interest. If no adjustment
#' variables survive selection, the interaction is tested without adjustment. In
#' the second stage, the selected (and possibly transformed) adjustment
#' variables are included as covariates when fitting the main-effects and
#' interaction models for each continuous variable listed in `cont_vars`.
#'
#' To force adjustment variables into the model without selection or
#' transformation, set `select = 1` and `alpha = 0`, or pass variable names via
#' `force_keep` and set `df = 1` for those variables. See [mfp2::mfp2()] for
#' further details.
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
#'   the primary estimands. Internally, values are remapped to 0, 1, 2, … in
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
#'   Numeric. The nominal significance level used as the default threshold when
#'   testing for interaction under `criterion = "pvalue"`. Default is `0.05`.
#'   Ignored when `criterion` is `"aic"` or `"bic"`.
#' @param min_improvement
#'   Numeric. The minimum improvement in the selection criterion required to
#'   include an interaction term. When `criterion = "pvalue"` this is a p-value
#'   threshold (e.g. `0.05`); when `criterion = "aic"` or `"bic"` it is the
#'   minimum decrease in the information criterion (must be positive; typical
#'   values are `1` or `2`). Defaults to `p_interact` for p-value selection and
#'   `2` for AIC/BIC selection.
#' @param include_group_var
#'   Logical. Whether `group_var` should be passed to [mfp2::mfp2()] alongside
#'   the other adjustment variables when selecting the adjustment model in stage
#'   one. Default is `FALSE`. When `TRUE`, the nominal significance level for
#'   `group_var` is set to 1, forcing it into the adjustment model.
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
#'   estimate scaling factors automatically (see [mfp2::mfp2()]). Set
#'   `scale = 1` to disable scaling.
#' @param shift
#'   A numeric vector of length \eqn{p} or a single numeric giving shift terms
#'   added to columns of `x` before scaling. Default is `NULL`, which estimates
#'   shift factors automatically (see [mfp2::mfp2()]). Set `shift = 0` to
#'   disable shifting. Shifting is applied to the data used in the final
#'   interaction model.
#' @param df
#'   A numeric vector of length \eqn{p} or a single positive integer setting
#'   the default degrees of freedom for each predictor. Degrees of freedom
#'   equal twice the FP degree (e.g. `df = 2` for FP1, `df = 4` for FP2).
#'   The program overrides user-supplied values based on the number of unique
#'   values: variables with 2–3 unique values are set to `df = 1` (linear);
#'   4–5 unique values to `df = min(2, default)`; and 6 or more unique values
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
#'   A character string specifying the criterion used by [mfp2::mfp2()] to
#'   select adjustment variables and FP degrees. One of `"pvalue"` (default),
#'   `"aic"`, or `"bic"`. Note that this criterion governs adjustment-variable
#'   selection only; the interaction test always reports p-values, AIC
#'   differences, and BIC differences regardless of this setting.
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
#'   `x`. See [mfp2::mfp2()] for details.
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
#'   An optional character vector naming variables that may contain zeros or
#'   negative values, for which the standard positivity check before log or
#'   negative-power transformation is suppressed. Variables in `zero_vars` that
#'   contain only positive values are silently reset to standard processing with
#'   a warning.
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
#' An object of class `"mfpi"`. Use `summary.mfpi()` for a formatted summary.
#' The object is a list with the following components:
#'
#' \describe{
#'   \item{`model_evaluation_metrics`}{A data frame of evaluation metrics for
#'     the main-effects and interaction models for each variable in `cont_vars`.
#'     Columns include: `pow_main` and `pow_int` (selected FP powers);
#'     `dev_int` (deviance of the interaction model); `dev_diff` (deviance of
#'     main-effects model minus deviance of interaction model); `pvalue`
#'     (likelihood-ratio p-value); `AIC_int` and `BIC_int` (information
#'     criteria for the interaction model); `AIC_diff` and `BIC_diff`
#'     (main-effects minus interaction model for each criterion).}
#'   \item{`adjust_terms`}{A data frame describing the FP terms selected for
#'     the adjustment model; see [mfp2::mfp2()] for the column definitions.}
#'   \item{`interaction_models`}{A named list of fitted model objects, one per
#'     variable in `cont_vars`.}
#'   \item{`fitted_functions`}{A named list of data frames, one per variable in
#'     `cont_vars`. Each data frame contains the estimated FP function
#'     \eqn{f_j} and its pointwise standard error \eqn{\text{se}(f_j)} at each
#'     level \eqn{j} of `group_var`, together with 95\% confidence intervals
#'     (`fj_lower`, `fj_upper`). Columns `fj - f0` and `se(fj - f0)` give the
#'     difference relative to the reference group, with corresponding confidence
#'     bounds.}
#'   \item{`group_var`}{The name of the grouping variable, as a character
#'     string.}
#'   \item{`show_models`}{The value of the `show_models` argument.}
#'   \item{`flex`}{The chosen flexibility level, as a character string.}
#'   \item{`group_levels_new`}{The remapped levels of `group_var` (0, 1, 2,
#'     …).}
#'   \item{`group_levels_original`}{The original levels of `group_var`.}
#'   \item{`family`}{The regression family, as a character string.}
#'   \item{`nobs`}{The number of observations used in model fitting.}
#' }
#'
#' @references
#' Royston, P. and Sauerbrei, W. (2004). A new approach to modelling
#' interactions between treatment and continuous covariates in clinical trials
#' by using fractional polynomials. *Statistics in Medicine*, 23, 2509–2525.
#'
#' Royston, P. and Sauerbrei, W. (2008). Interactions between treatment and
#' continuous covariates — a step towards individualising therapy (Editorial).
#' *Journal of Clinical Oncology*, 26, 1397–1399.
#'
#' Royston, P. and Sauerbrei, W. (2013). Interaction of treatment with a
#' continuous variable: simulation study of significance level for several
#' methods of analysis. *Statistics in Medicine*, 32, 3788–3803.
#'
#' Royston, P. and Sauerbrei, W. (2014). Interaction of treatment with a
#' continuous variable: simulation study of power for several methods of
#' analysis. *Statistics in Medicine*, 33, 4695–4708.
#'
#' Royston, P., Sauerbrei, W. and Ritchie, A. (2004). Is treatment with
#' interferon-alpha effective in all patients with metastatic renal carcinoma?
#' A new approach to the investigation of interactions. *British Journal of
#' Cancer*, 90, 794–799.
#'
#' Sauerbrei, W., Royston, P. and Zapien, K. (2007). Detecting an interaction
#' between treatment and a continuous covariate: a comparison of two
#' approaches. *Computational Statistics and Data Analysis*, 51, 4054–4063.
#'
#' Sauerbrei, W. and Royston, P. (2022). Investigating treatment-effect
#' modification by a continuous covariate in IPD meta-analysis: an approach
#' using fractional polynomials. *BMC Medical Research Methodology*, 22, 1–13.
#'
#' Royston, P. and Sauerbrei, W. (2008). *Multivariable Model-Building: A
#' Pragmatic Approach to Regression Analysis Based on Fractional Polynomials
#' for Modelling Continuous Variables*. John Wiley & Sons.
#'
#' @seealso [mfp2::summary.mfpi()], [mfp2::mfp2()]
#'
#' @export
mfpi <- function(x, ...) {
  UseMethod("mfpi", x)
}

#' @describeIn mfpi Default method accepting a numeric matrix `x` and response
#'   vector `y`.
#' @import mfp2
#' @export
mfpi.default <- function(
    x,
    y,
    group_var            = NULL,
    cont_vars            = NULL,
    flex                 = c("flex1", "flex2", "flex3", "flex4"),
    p_interact           = 0.05,
    min_improvement      = NULL,
    include_group_var    = FALSE,
    show_models          = FALSE,
    weights              = NULL,
    offset               = NULL,
    cycles               = 10,
    scale                = NULL,
    shift                = NULL,
    df                   = 4,
    center               = TRUE,
    subset               = NULL,
    family               = c("gaussian", "poisson", "binomial", "cox"),
    criterion            = c("pvalue", "aic", "bic"),
    select               = 0.05,
    alpha                = 0.05,
    force_keep           = NULL,
    xorder               = c("ascending", "descending", "original"),
    fp_powers            = NULL,
    ties                 = c("breslow", "efron", "exact"),
    strata               = NULL,
    nocenter             = NULL,
    acd_vars             = NULL,
    zero_vars            = NULL,
    use_ftest            = FALSE,
    control              = NULL,
    verbose              = TRUE,
    digits               = 3,
    ...
) {
  cl <- match.call()
  
  # Match enumerated arguments -------------------------------------------------
  family           <- match.arg(family)
  xorder           <- match.arg(xorder)
  criterion        <- match.arg(criterion)
  flex             <- match.arg(flex)
  ties             <- match.arg(ties)
  
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
      "i Convert all categorical variables to numeric dummy columns before passing to `mfpi()`.",
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
  
  if (!group_var %in% vnames) {
    stop(
      paste0("! `group_var = '", group_var, "'` is not a column of `x`."),
      call. = FALSE
    )
  }
  
  if (length(group_var) != 1L) {
    stop(
      paste0("! `group_var` must name exactly one variable; ", length(group_var), " were supplied."),
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
      paste0(
        "! The following variables in `cont_vars` are not columns of `x`: ",
        paste(missing_cont, collapse = ", "), "."
      ),
      call. = FALSE
    )
  }
  
  # Validate subset ------------------------------------------------------------
  if (!is.null(subset)) {
    if (!is.vector(subset)) {
      stop(
        paste0("! `subset` must be a vector, not an object of class ",
               paste(class(subset), collapse = ", "), "."),
        call. = FALSE
      )
    }
    if (any(subset < 0L)) {
      stop("! `subset` must not contain negative indices.", call. = FALSE)
    }
    if (length(subset) < 5L) {
      stop(
        paste0("! After subsetting, only ", length(subset), " observations remain; ",
               "at least 5 are required to fit an mfpi model."),
        call. = FALSE
      )
    }
  }
  
  # Validate weights -----------------------------------------------------------
  if (!is.null(weights)) {
    if (any(weights < 0)) {
      stop("! `weights` must not be negative.", call. = FALSE)
    }
    if (length(weights) != nobs) {
      stop(
        paste0("! Length of `weights` (", length(weights), ") must equal ",
               "the number of rows in `x` (", nobs, ")."),
        call. = FALSE
      )
    }
  }
  
  # Validate offset ------------------------------------------------------------
  if (!is.null(offset) && length(offset) != nobs) {
    stop(
      paste0("! Length of `offset` (", length(offset), ") must equal ",
             "the number of rows in `x` (", nobs, ")."),
      call. = FALSE
    )
  }
  
  # Validate alpha and select --------------------------------------------------
  if (any(alpha < 0) || any(alpha > 1)) {
    stop("! All values of `alpha` must be in [0, 1].", call. = FALSE)
  }
  if (length(alpha) != 1L && length(alpha) != nvars) {
    stop(
      paste0("! `alpha` must be a single number or a vector of length ", nvars,
             " (number of columns in `x`); got length ", length(alpha), "."),
      call. = FALSE
    )
  }
  
  if (any(select < 0) || any(select > 1)) {
    stop("! All values of `select` must be in [0, 1].", call. = FALSE)
  }
  if (length(select) != 1L && length(select) != nvars) {
    stop(
      paste0("! `select` must be a single number or a vector of length ", nvars,
             " (number of columns in `x`); got length ", length(select), "."),
      call. = FALSE
    )
  }
  
  # Validate force_keep --------------------------------------------------------
  if (!is.null(force_keep) && !all(force_keep %in% vnames)) {
    stop("! All variables named in `force_keep` must be columns of `x`.", call. = FALSE)
  }
  
  # Validate shift and scale ---------------------------------------------------
  if (!is.null(shift) && length(shift) != 1L && length(shift) != nvars) {
    stop(
      paste0("! `shift` must be NULL, a single number, or a vector of length ", nvars,
             "; got length ", length(shift), "."),
      call. = FALSE
    )
  }
  if (!is.null(scale) && length(scale) != 1L && length(scale) != nvars) {
    stop(
      paste0("! `scale` must be NULL, a single number, or a vector of length ", nvars,
             "; got length ", length(scale), "."),
      call. = FALSE
    )
  }
  
  # Validate center ------------------------------------------------------------
  if (length(center) != 1L && length(center) != nvars) {
    stop(
      paste0("! `center` must be a single logical or a logical vector of length ", nvars,
             "; got length ", length(center), "."),
      call. = FALSE
    )
  }
  
  # Validate acd_vars ----------------------------------------------------------
  if (!is.null(acd_vars) && !all(acd_vars %in% vnames)) {
    stop("! All variables named in `acd_vars` must be columns of `x`.", call. = FALSE)
  }
  
  # Validate df ----------------------------------------------------------------
  if (any(df <= 0L)) {
    stop(
      "! All values of `df` must be positive.\n",
      "i Use `df = 1` for linear terms or `df = 2m` for an FP of degree m.",
      call. = FALSE
    )
  }
  if (length(df) == 1L) {
    if (df != 1L && df %% 2L != 0L) {
      stop(
        paste0("! `df = ", df, "` is invalid. `df` must be 1 (linear) or an even ",
               "number equal to 2m, where m is the FP degree."),
        call. = FALSE
      )
    }
  } else {
    if (length(df) != nvars) {
      stop(
        paste0("! When `df` is a vector it must have length ", nvars,
               " (number of columns in `x`); got length ", length(df), "."),
        call. = FALSE
      )
    }
    invalid_df <- df != 1L & df %% 2L != 0L
    if (any(invalid_df)) {
      stop(
        paste0("! Each element of `df` must be 1 (linear) or an even number (2m for FP degree m). ",
               "Invalid values found at positions: ",
               paste(which(invalid_df), collapse = ", "), "."),
        call. = FALSE
      )
    }
  }
  
  # Warn if use_ftest is incompatible with family ------------------------------
  if (use_ftest && family != "gaussian") {
    warning(
      paste0("i `use_ftest = TRUE` is only applicable to Gaussian models; ",
             "reverting to chi-square test for family = '", family, "'."),
      call. = FALSE
    )
    use_ftest <- FALSE
  }
  
  # Validate response y --------------------------------------------------------
  if (family == "cox") {
    if (!survival::is.Surv(y)) {
      stop("! For `family = 'cox'`, `y` must be a `survival::Surv` object.", call. = FALSE)
    }
    if (nrow(y) != nobs) {
      stop(
        paste0("! `y` has ", nrow(y), " rows but `x` has ", nobs,
               " rows; they must match."),
        call. = FALSE
      )
    }
    if (attr(y, "type") != "right") {
      stop(
        paste0("! Only right-censored survival data are currently supported; ",
               "`y` has censoring type '", attr(y, "type"), "'."),
        call. = FALSE
      )
    }
    if (is.factor(strata)) {
      strata <- as.numeric(strata)
    }
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
      stop(
        paste0("! For family = '", family, "', `y` must be a vector, not an object of class ",
               paste(class(y), collapse = ", "), "."),
        call. = FALSE
      )
    }
    if (length(y) != nobs) {
      stop(
        paste0("! `y` has length ", length(y), " but `x` has ", nobs,
               " rows; they must match."),
        call. = FALSE
      )
    }
  }
  
  # Validate and build fp_powers -----------------------------------------------
  if (!is.null(fp_powers) && !is.list(fp_powers)) {
    stop("! `fp_powers` must be a named list.", call. = FALSE)
  }
  
  default_powers <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
  power_list <- setNames(replicate(nvars, default_powers, simplify = FALSE), vnames)
  
  if (!is.null(fp_powers)) {
    if (length(fp_powers) != sum(nchar(names(fp_powers)) > 0L, na.rm = TRUE)) {
      stop("! Every element of `fp_powers` must have a name.", call. = FALSE)
    }
    unknown_names <- setdiff(names(fp_powers), vnames)
    if (length(unknown_names) > 0L) {
      stop(
        paste0("! The following names in `fp_powers` do not match any column of `x`: ",
               paste(unknown_names, collapse = ", "), "."),
        call. = FALSE
      )
    }
    if (!all(sapply(fp_powers, is.numeric))) {
      stop("! All elements of `fp_powers` must be numeric vectors.", call. = FALSE)
    }
    fp_powers <- lapply(fp_powers, sort)
    pw_lengths <- sapply(fp_powers, length)
    too_short  <- names(pw_lengths[pw_lengths < 2L])
    if (length(too_short) > 0L) {
      stop(
        paste0("! Each element of `fp_powers` must contain at least two values. ",
               "Insufficient values for: ", paste(too_short, collapse = ", "), "."),
        call. = FALSE
      )
    }
    power_list <- modifyList(power_list, Filter(Negate(is.null), fp_powers))
  }
  
  # Warn about zero_vars consistency -------------------------------------------
  if (!is.null(zero_vars) && !all(zero_vars %in% vnames)) {
    warning(
      "i Some variables named in `zero_vars` are not columns of `x`; they will be ignored.",
      call. = FALSE
    )
  }
  
  # Set defaults for scalar/vector arguments -----------------------------------
  if (is.null(min_improvement)) {
    min_improvement <- switch(criterion, pvalue = p_interact, aic = 2, bic = 2)
  }
  if (is.null(weights))  weights <- rep.int(1, nobs)
  if (is.null(offset))   offset  <- rep.int(0, nobs)
  if (length(select) == 1L) select <- rep(select, nvars)
  if (length(alpha)  == 1L) alpha  <- rep(alpha,  nvars)
  
  if (is.null(shift)) {
    shift <- apply(x, 2L, find_shift_factor)
  } else if (length(shift) == 1L) {
    shift <- rep(shift, nvars)
  }
  
  if (is.null(scale)) {
    scale <- apply(x, 2L, find_scale_factor)
  } else if (length(scale) == 1L) {
    scale <- rep(scale, nvars)
  }
  
  if (length(center) == 1L) {
    center <- setNames(rep(center, nvars), vnames)
  }
  
  if (is.null(control)) {
    control <- if (family == "cox") survival::coxph.control() else stats::glm.control()
  }
  
  # Process zero_vars to a logical vector --------------------------------------
  zero_vars <- intersect(zero_vars, vnames)
  zero_flag  <- setNames(rep(FALSE, nvars), vnames)
  
  if (length(zero_vars) > 0L) {
    zero_flag[zero_vars] <- TRUE
    # Reset flag for variables that contain only positive values
    all_positive <- zero_vars[
      apply(x[, zero_vars, drop = FALSE], 2L,
            function(col) all(col > 0, na.rm = TRUE))
    ]
    if (length(all_positive) > 0L) {
      warning(
        paste0("i The following variables in `zero_vars` contain only positive values ",
               "and have been reset to standard processing: ",
               paste(all_positive, collapse = ", "), "."),
        call. = FALSE
      )
      zero_flag[all_positive] <- FALSE
    }
  }
  
  # acd_vars to logical vector -------------------------------------------------
  if (is.null(acd_vars)) {
    acd_flag <- setNames(rep(FALSE, nvars), vnames)
  } else {
    acd_vars <- unique(intersect(acd_vars, vnames))
    acd_flag <- setNames(rep(FALSE, nvars), vnames)
    acd_flag[acd_vars] <- TRUE
  }
  
  # Set degrees of freedom per variable ----------------------------------------
  if (length(df) == 1L) {
    df_list <- if (df != 1L) mfp2::assign_df(x = x, df_default = df) else rep(df, nvars)
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
  
  # Apply shift and scale to x -------------------------------------------------
  x <- sweep(x, 2L, shift, "+")
  x <- sweep(x, 2L, scale, "/")
  
  # Prepare stratification for Cox models --------------------------------------
  istrata <- strata
  if (family == "cox" && !is.null(strata)) {
    istrata <- as.integer(survival::strata(strata, shortlabel = TRUE))
  }
  
  # Apply subset ---------------------------------------------------------------
  if (!is.null(subset)) {
    x       <- x[subset, , drop = FALSE]
    y       <- if (family == "cox") y[subset, , drop = FALSE] else y[subset]
    weights <- weights[subset]
    offset  <- offset[subset]
    istrata <- istrata[subset]
  }
  
  # Fit the MFPI model ---------------------------------------------------------
  fit <- fit_mfpi(
    x                    = x,
    y                    = y,
    family               = family,
    weights              = weights,
    offset               = offset,
    cycles               = cycles,
    center               = center,
    criterion            = criterion,
    select               = select,
    alpha                = alpha,
    df                   = df_list,
    force_keep           = force_keep,
    xorder               = xorder,
    fp_powers            = power_list,
    ties                 = ties,
    strata               = istrata,
    nocenter             = nocenter,
    acd_vars             = acd_flag,
    use_ftest            = use_ftest,
    control              = control,
    verbose              = verbose,
    group_var            = group_var,
    include_group_var    = include_group_var,
    flex                 = flex,
    zero_vars            = zero_flag,
    cont_vars            = cont_vars,
    p_interact           = p_interact,
    show_models          = show_models,
    min_improvement      = min_improvement,
    digits               = digits
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


#' Print mfpi Object
#'
#' This function prints the adjustment and interaction models stored in an mfpi object.
#'
#' @param x An mfpi object.
#' @param ... not used.
#' @return prints the mfpi models
#' @method print mfpi
#' @export
print.mfpi <- function(x, ...) {
  cat("Adjustment Variable Selected Using MFP Algorithm:\n")
  print(x$adj)
  cat("----------------------------------------------------------------\n")
  cat("Interaction Test", "\n")
  cat("----------------------------------------------------------------\n")
  cat(sprintf("\nInteraction with %s (%d Observations). %s Strategy:\n", x$catvar, x$nobs, x$flex))
  
  print(x$model_evaluation_metrics)
  
  cat("\ndf = degrees of freedom for interaction test; tdf = total degrees of freedom for interaction model\n")
  
  if (x$flex == "Flex 4") {
    cat(sprintf("\nThe second column on `pow_int` contains FP powers for each level of %s variable.\n", x$catvar))
  }
  
  if (x$showmodel) {
    cat("-------------------------------------------------------------\n")
    cat("Interaction Models for Each Continuous Variable of Interest:", "\n")
    cat("-------------------------------------------------------------\n")
    
    n <- length(x$interaction_models)
    catx <- c("Linear", "FP1", "FP2")
    
    for (i in 1:n) {
      vat <- x$interaction_models[[i]]
      if (length(vat) != 0) {
        numv <- length(vat)
        xname <- names(vat)
        for (j in 1:numv) {
          cat(sprintf("Interaction between the categorical variable '%s' and '%s' with %s functional form:\n",
                      x$catvar, xname[j], catx[i]))
          print(summary(vat[[j]]$fit))
        }
      }
    }
  }
}


#' Summarizing `mfpi` model fits
#'
#' This function is a method for the generic [base::summary()] function for
#' objects of class `mfpi`.
#'
#' @param object an object of class `mfpi`, usually, a result of a call to
#' \code{mfpi()}.
#' @param ... further arguments passed to the summary functions for `glm()`
#' ([stats::summary.glm()], i.e. families supported by `glm()`) or `coxph()`
#' ([survival::summary.coxph()], if `object$family = "cox"`).
#'
#' @return
#' An object returned from [stats::summary.glm()] or
#' [survival::summary.coxph()], depending on the family parameter of `object`.
#'
#' @seealso
#' [mfp2::mfp2()], [stats::glm()], [stats::summary.glm()], [survival::coxph()],
#' [survival::summary.coxph()]
#'
#' @export
summary.mfpi <- function(object, ...) {
  # Degrees of freedom in glm() or coxph() do not account for fp power terms
  # in calculating pvalues
  lapply(object$interaction_models,
         function(v) lapply(v, function(x) summary(x$fit)))
  
}
