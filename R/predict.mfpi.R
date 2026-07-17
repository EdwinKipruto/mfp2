# -----------------------------------------------------------------------------
# Shared MFPI prediction helpers ----------------------------------------------
# -----------------------------------------------------------------------------

#' Build a Group-Specific Fractional-Polynomial Basis
#'
#' Builds the fractional-polynomial (FP) design block used to evaluate MFPI
#' group-specific functions and ordinary subject-level interaction effects.
#' The function starts from either a one-column continuous variable or an
#' already transformed basis. In the usual prediction path, it duplicates the
#' continuous variable once per group, applies the selected FP powers for each
#' group, checks that the generated block size matches the fitted model
#' metadata, and then assigns the exact coefficient names used in the fitted
#' interaction model.
#'
#' @section Scale convention:
#' When \code{transform_basis = TRUE}, \code{cont_mat} must be on the
#' shifted-but-not-scaled FP scale used by \code{transform_matrix()}, usually
#' \eqn{x + shift}. This is the same scale on which the fitted interaction
#' coefficients were estimated. The caller is responsible for multiplying a
#' stored shifted/scaled variable by its scale factor before calling this helper.
#'
#' @section Column alignment:
#' \code{coefficient_groups} is the authoritative mapping between group levels
#' and fitted model coefficient columns. The function never infers coefficient
#' groups using regular expressions. It validates that the number of columns
#' implied by \code{coefficient_groups} matches the number of FP powers in
#' \code{group_fp_powers} for each group. This is important for FP2 terms,
#' repeated powers, and \code{flex4}, where groups can have different FP power
#' sets.
#'
#' @param cont_mat Numeric matrix. If \code{transform_basis = TRUE}, this must
#'   be a one-column matrix containing the continuous variable on the
#'   shifted-but-not-scaled FP scale. If \code{transform_basis = FALSE}, this
#'   must already be the transformed FP basis with one column per fitted
#'   interaction coefficient.
#' @param group_fp_powers List of numeric FP power vectors, one element per
#'   group. Names may be absent, internal group labels, or the temporary
#'   transformation names generated inside this helper.
#' @param coefficient_groups List mapping each group to the exact fitted
#'   interaction coefficient names for that group's FP block. Each element must
#'   be a non-empty character vector.
#' @param group_levels Vector of group labels in fitted-model order. Values are
#'   coerced to character labels for naming and ordering.
#' @param cont_name Character scalar naming the continuous variable. This is used
#'   only to build temporary input names for \code{transform_matrix()}.
#' @param transform_basis Logical scalar. If \code{TRUE}, duplicate and
#'   transform \code{cont_mat}. If \code{FALSE}, validate \code{cont_mat} as
#'   an already transformed basis and rename its columns to the fitted
#'   coefficient names.
#' @param zero_var Logical scalar. Whether non-positive values should be treated
#'   as structural zeros during FP transformation. Centering and final
#'   structural-zero restoration are handled by higher-level prediction helpers.
#' @param check_binary Logical scalar passed to \code{transform_matrix()}.
#'
#' @return A numeric matrix whose columns are named exactly as
#'   \code{unlist(coefficient_groups, use.names = FALSE)}. The returned matrix
#'   has attributes \code{coefficient_groups}, \code{group_fp_powers}, and
#'   \code{group_levels}, all normalized to fitted-model group order.
#'
#' @keywords internal
#' @noRd
build_group_fp_basis <- function(cont_mat,
                                 group_fp_powers,
                                 coefficient_groups,
                                 group_levels,
                                 cont_name,
                                 transform_basis = TRUE,
                                 zero_var = FALSE,
                                 check_binary = TRUE) {
  # This helper has one job: rebuild the exact FP columns that the fitted
  # interaction model estimated. It therefore validates names, group order,
  # and per-group block sizes before doing any transformation.
  
  if (!is.character(cont_name) || length(cont_name) != 1L ||
      is.na(cont_name) || !nzchar(cont_name)) {
    stop("`cont_name` must be a single non-empty character string.",
         call. = FALSE)
  }
  
  if (!is.logical(transform_basis) || length(transform_basis) != 1L ||
      anyNA(transform_basis)) {
    stop("`transform_basis` must be a single non-missing logical value.",
         call. = FALSE)
  }
  
  if (!is.logical(zero_var) || length(zero_var) != 1L || anyNA(zero_var)) {
    stop("`zero_var` must be a single non-missing logical value.",
         call. = FALSE)
  }
  
  if (!is.logical(check_binary) || length(check_binary) != 1L ||
      anyNA(check_binary)) {
    stop("`check_binary` must be a single non-missing logical value.",
         call. = FALSE)
  }
  
  if (missing(group_levels) || length(group_levels) == 0L ||
      anyNA(group_levels)) {
    stop("`group_levels` must be a non-empty vector without missing values.",
         call. = FALSE)
  }
  
  group_labels <- as.character(group_levels)
  n_groups <- length(group_labels)
  temp_names <- paste0(cont_name, group_labels, "1")
  
  # Coefficient groups are the authoritative mapping from group to fitted model
  # columns. They prevent regex-based name parsing and are required for FP2 and
  # flex4, where groups can have different basis sizes.
  if (!is.list(coefficient_groups) || length(coefficient_groups) != n_groups) {
    stop("`coefficient_groups` must be a list with one element per group.",
         call. = FALSE)
  }
  
  if (!is.null(names(coefficient_groups)) &&
      anyDuplicated(names(coefficient_groups))) {
    stop("`coefficient_groups` must not contain duplicated names.",
         call. = FALSE)
  }
  
  if (is.null(names(coefficient_groups))) {
    names(coefficient_groups) <- group_labels
  } else {
    cg_names <- names(coefficient_groups)
    if (setequal(cg_names, group_labels)) {
      coefficient_groups <- coefficient_groups[group_labels]
    } else if (setequal(cg_names, temp_names)) {
      coefficient_groups <- coefficient_groups[temp_names]
      names(coefficient_groups) <- group_labels
    } else {
      stop(
        paste0(
          "`coefficient_groups` names must match `group_levels` or temporary ",
          "transform names: ", paste(temp_names, collapse = ", "), "."
        ),
        call. = FALSE
      )
    }
  }
  
  if (!is.list(group_fp_powers) || length(group_fp_powers) != n_groups) {
    stop("`group_fp_powers` must be a list with one element per group.",
         call. = FALSE)
  }
  
  if (!is.null(names(group_fp_powers)) &&
      anyDuplicated(names(group_fp_powers))) {
    stop("`group_fp_powers` must not contain duplicated names.",
         call. = FALSE)
  }
  
  if (!is.null(names(group_fp_powers))) {
    gp_names <- names(group_fp_powers)
    if (setequal(gp_names, group_labels)) {
      group_fp_powers <- group_fp_powers[group_labels]
    } else if (setequal(gp_names, temp_names)) {
      group_fp_powers <- group_fp_powers[temp_names]
    } else {
      stop(
        paste0(
          "`group_fp_powers` names must match either `group_levels` or ",
          "temporary transform names: ", paste(temp_names, collapse = ", "), "."
        ),
        call. = FALSE
      )
    }
  }
  
  group_fp_powers_by_group <- stats::setNames(group_fp_powers, group_labels)
  
  expected_group_sizes <- vapply(coefficient_groups, length, integer(1L))
  power_group_sizes <- vapply(group_fp_powers_by_group, length, integer(1L))
  
  if (any(expected_group_sizes <= 0L)) {
    stop("Each element of `coefficient_groups` must be non-empty.",
         call. = FALSE)
  }
  
  if (any(power_group_sizes <= 0L)) {
    stop("Each element of `group_fp_powers` must contain at least one power.",
         call. = FALSE)
  }
  
  bad_power <- vapply(
    group_fp_powers_by_group,
    function(p) !is.numeric(p) || length(p) == 0L || anyNA(p) || any(!is.finite(p)),
    logical(1L)
  )
  
  if (any(bad_power)) {
    stop("`group_fp_powers` must contain finite numeric FP power vectors.",
         call. = FALSE)
  }
  
  if (!identical(as.integer(expected_group_sizes),
                 as.integer(power_group_sizes))) {
    stop(
      paste0(
        "Per-group coefficient block sizes do not match the selected ",
        "group-specific FP powers."
      ),
      call. = FALSE
    )
  }
  
  x_split_names <- unlist(coefficient_groups, use.names = FALSE)
  
  if (!is.character(x_split_names) || anyNA(x_split_names) ||
      any(!nzchar(x_split_names))) {
    stop("`coefficient_groups` must contain non-empty coefficient names.",
         call. = FALSE)
  }
  
  if (anyDuplicated(x_split_names)) {
    stop("`coefficient_groups` must not contain duplicated coefficient names.",
         call. = FALSE)
  }
  
  if (!transform_basis) {
    if (!is.matrix(cont_mat) || !is.numeric(cont_mat)) {
      stop("`cont_mat` must be a numeric matrix when `transform_basis = FALSE`.",
           call. = FALSE)
    }
    
    if (ncol(cont_mat) != length(x_split_names)) {
      stop(
        paste0(
          "`transform_basis = FALSE` requires `cont_mat` to have ",
          length(x_split_names), " columns."
        ),
        call. = FALSE
      )
    }
    
    if (anyNA(cont_mat) || any(!is.finite(cont_mat))) {
      stop("`cont_mat` contains non-finite values.", call. = FALSE)
    }
    
    out <- as.matrix(cont_mat)
    storage.mode(out) <- "double"
    colnames(out) <- x_split_names
    attr(out, "coefficient_groups") <- coefficient_groups
    attr(out, "group_fp_powers") <- group_fp_powers_by_group
    attr(out, "group_levels") <- group_labels
    return(out)
  }
  
  if (!is.matrix(cont_mat) || ncol(cont_mat) != 1L || !is.numeric(cont_mat)) {
    stop("`cont_mat` must be a one-column numeric matrix.", call. = FALSE)
  }
  
  if (is.null(colnames(cont_mat)) || length(colnames(cont_mat)) != 1L) {
    stop("`cont_mat` must have exactly one column name.", call. = FALSE)
  }
  
  if (anyNA(cont_mat) || any(!is.finite(cont_mat))) {
    stop("`cont_mat` must contain only finite numeric values.", call. = FALSE)
  }
  
  # transform_matrix() requires names(power_list) to match colnames(x). The
  # temporary input names are discarded after transformation; final coefficient
  # names come from coefficient_groups.
  x_in <- cont_mat[, rep(1L, n_groups), drop = FALSE]
  colnames(x_in) <- temp_names
  
  group_fp_powers_tm <- group_fp_powers_by_group
  names(group_fp_powers_tm) <- temp_names
  
  center_map <- stats::setNames(rep(FALSE, n_groups), temp_names)
  acd_map <- stats::setNames(rep(FALSE, n_groups), temp_names)
  zero_map <- stats::setNames(rep(isTRUE(zero_var), n_groups), temp_names)
  
  x_out <- transform_matrix(
    x = x_in,
    power_list = group_fp_powers_tm,
    center = center_map,
    acdx = acd_map,
    zero = zero_map,
    check_binary = check_binary
  )$x_transformed
  
  if (is.null(x_out)) {
    stop("FP basis construction returned no transformed columns.",
         call. = FALSE)
  }
  
  x_out <- as.matrix(x_out)
  storage.mode(x_out) <- "double"
  
  if (ncol(x_out) != sum(power_group_sizes)) {
    stop(
      paste0(
        "Transformed fitted-function basis has ", ncol(x_out),
        " columns, but expected ", sum(power_group_sizes), "."
      ),
      call. = FALSE
    )
  }
  
  if (length(x_split_names) != ncol(x_out)) {
    stop("Fitted-function basis columns do not match coefficient groups.",
         call. = FALSE)
  }
  
  if (anyNA(x_out) || any(!is.finite(x_out))) {
    stop("Non-finite values were produced in the fitted-function basis.",
         call. = FALSE)
  }
  
  colnames(x_out) <- x_split_names
  attr(x_out, "coefficient_groups") <- coefficient_groups
  attr(x_out, "group_fp_powers") <- group_fp_powers_by_group
  attr(x_out, "group_levels") <- group_labels
  x_out
}


#' Predict from an MFPI Model
#'
#' Obtain group-specific fitted curves, differences between group curves, or
#' predictions for individual observations from a model fitted with [mfpi()].
#'
#' @details
#' MFPI fits a separate interaction model for each continuous variable examined
#' through `cont_vars`. Predictions are therefore produced separately for each
#' requested variable. Predictions for two different terms come from two
#' different term-specific models; they are not parts of one combined model.
#'
#' Supply continuous predictors on their original scale. Stored shifts, scales,
#' fractional-polynomial functions, centering values, zero handling, and
#' Winsorisation limits are applied automatically.
#'
#' The default `type = NULL` is equivalent to `type = "both"`. Supplying
#' `newdata` by itself does not change this default. To obtain one prediction per
#' row, explicitly use `type = "link"`, `"response"`, or one of the Cox
#' prediction types.
#'
#' @section Common calls:
#' - `predict(fit, terms = "age")` returns the fitted group curves for age and
#'   their differences.
#' - `predict(fit, terms = "age", type = "function", grid = TRUE)` returns
#'   group-specific fitted curves on an equally spaced grid.
#' - `predict(fit, terms = "age", type = "difference", grid = TRUE)` returns
#'   differences between each comparison-group curve and the reference-group
#'   curve.
#' - `predict(fit, terms = "age", type = "response", newdata = new_data)`
#'   returns one response-scale prediction per row for a GLM.
#' - `predict(fit, terms = "age", type = "risk", newdata = new_data)` returns
#'   one relative-risk prediction per row for a Cox model.
#'
#' @section Selecting terms and models:
#' Use `terms` to select continuous variables that were examined by [mfpi()].
#' These are normally names supplied through `cont_vars`, not adjustment-variable
#' names.
#'
#' MFPI stores one term-specific interaction model for each evaluated continuous
#' variable. The `model` argument determines which of these models can be used:
#'
#' - `model = "best"` uses only variables whose interaction was retained by the
#'   MFPI selection criterion. This is the default.
#' - `model = "all"` uses the stored model for every evaluated continuous
#'   variable for which a model is available, including interactions that were
#'   not retained.
#'
#' `model = "all"` does not return every candidate model considered during
#' fitting. It returns the stored term-specific model for each evaluated
#' continuous variable.
#'
#' When `terms = NULL`, all terms available under the selected `model` setting
#' are used. If no interactions were retained, `model = "best"` has no models
#' available; use `model = "all"` to inspect the evaluated terms.
#'
#' @section Choosing the prediction type:
#' Use the following types for fitted curves:
#'
#' - `type = "function"` returns the fitted curve for each group.
#' - `type = "difference"` returns each comparison-group curve minus the
#'   reference-group curve at the same continuous-variable value.
#' - `type = "both"` returns both fitted curves and their differences. This is
#'   the default when `type = NULL`.
#'
#' Use the following types for predictions for individual observations:
#'
#' - For Gaussian, binomial, and Poisson models, `type = "link"` returns the
#'   linear predictor and `type = "response"` applies the inverse link function.
#'   `"lp"` is accepted as another name for `"link"`.
#' - For Cox models, `type = "lp"` returns the relative log-hazard,
#'   `type = "risk"` returns relative hazard, `type = "expected"` returns the
#'   predicted cumulative hazard up to the supplied follow-up time, and
#'   `type = "survival"` returns the estimated survival probability at that
#'   time. `"link"` is accepted as another name for `"lp"`.
#'
#' For Cox models, `type = "risk"` is not the same as a GLM response
#' prediction.
#'
#' @section Fitted curves and group differences:
#' For `type = "function"`, `"difference"`, or `"both"`, only the requested
#' continuous variable is required in `newdata`. The grouping variable and
#' adjustment variables are not required because the function is evaluated for
#' every fitted group at each supplied value.
#'
#' A group-specific fitted curve contains the part of the term-specific model
#' defined by the continuous variable and group. For a GLM, it also includes the
#' model intercept and the group's main effect. Adjustment variables and offsets
#' are not included in fitted-curve output.
#'
#' A fitted difference is the comparison-group curve minus the reference-group
#' curve at the same value. It therefore includes both the group main-effect
#' difference and any difference in the fitted continuous-variable functions.
#'
#' The reference group is the reference level used when [mfpi()] was fitted.
#' For more than two groups, one difference is returned for every non-reference
#' group. Original group labels are used in the returned tables when available.
#'
#' Fitted curves and differences are always reported on the model's
#' linear-predictor scale. With the standard links, this is:
#'
#' - the outcome scale for a Gaussian identity-link model;
#' - the log-odds scale for a binomial logit model;
#' - the log-mean scale for a Poisson log-link model;
#' - a partial log-hazard scale for a Cox model.
#'
#' For a Cox model, these fitted curves do not include the baseline hazard and
#' are not survival probabilities.
#'
#' @section Evaluation values and grids:
#' When `grid = FALSE`, fitted curves are evaluated at the requested variable
#' values in `newdata`, or at the values used to fit the model when
#' `newdata = NULL`.
#'
#' When `grid = TRUE`, an equally spaced grid is created between the smallest
#' and largest available values. The range is based on `newdata` when supplied
#' and on the fitting data otherwise. `n_grid` controls the number of grid
#' points.
#'
#' For a variable fitted with structural-zero handling, the grid is constructed
#' over the positive values. A zero point is also included when structural zeros
#' are present, so the result can contain `n_grid + 1` distinct values.
#'
#' If Winsorisation was used during fitting, the same limits are applied to new
#' values before prediction. Values outside the fitted limits are replaced by
#' the relevant limit. The returned `x` column shows the value actually used for
#' prediction after Winsorisation.
#'
#' `grid` is used only for fitted curves. It is ignored with a warning for
#' predictions for individual observations.
#'
#' @section Predictions for individual observations:
#' For `type = "link"`, `"response"`, `"lp"`, `"risk"`, `"expected"`, or
#' `"survival"`, one prediction is returned for each row.
#'
#' These predictions use the complete term-specific interaction model. Supply:
#'
#' - the continuous variable named in `terms`;
#' - the grouping variable used in the MFPI analysis;
#' - every adjustment variable retained in that term-specific model;
#' - any offset, strata, or follow-up information required for the selected
#'   prediction type.
#'
#' Different terms can require different prediction columns because the tested
#' continuous variable is removed from its own adjustment set.
#'
#' For a model fitted with a formula, normally supply the original variables.
#' Factors and retained formula expressions are reconstructed automatically.
#'
#' For a model fitted from a matrix, supply the required columns with the same
#' names and coding used during fitting. If adjustment columns were grouped with
#' `term_groups`, every member column of a required group must be supplied.
#'
#' New group levels and new factor levels are not supported. Required prediction
#' values must not be missing and numeric values must be finite. Extra columns
#' not used by the selected term-specific model are ignored.
#'
#' Stored transformations and Winsorisation limits are applied automatically.
#' Do not create transformed or group-interaction columns yourself.
#'
#' When `newdata = NULL`, predictions are returned for the observations used to
#' fit each term-specific model.
#'
#' @section Offsets:
#' Offsets are used only for predictions for individual observations. They are
#' not included in fitted curves or fitted differences.
#'
#' When `newdata = NULL`, the fitted offset is used unless `newoffset` is
#' supplied to replace it. A replacement offset must contain one value for each
#' fitting observation.
#'
#' When `newdata` is supplied:
#'
#' - If the fitted formula contains `offset(...)`, include the source variable
#'   or variables needed by that expression in `newdata`. The offset is
#'   evaluated automatically.
#' - If the offset was supplied through the separate `offset` argument, supply
#'   the offset for the new rows through `newoffset`.
#' - A supplied `newoffset` takes precedence over a formula offset.
#'
#' Supply one finite `newoffset` value for each prediction row. `newoffset`
#' cannot be used when the relevant interaction model was fitted without an
#' offset.
#'
#' @section Cox predictions:
#' Cox predictions with `type = "lp"` and `type = "risk"` are relative
#' quantities. They are not absolute hazards or survival probabilities.
#'
#' Use `cox_reference` to choose the covariate reference for these predictions:
#'
#' - `"zero"` uses zero on the fitted predictor scale and is the default.
#' - `"sample"` uses the fitted-sample predictor means.
#' - `"strata"` uses fitted-sample means within each stratum.
#'
#' Changing `cox_reference` changes the displayed linear predictors and risks,
#' but comparisons made using the same reference remain unchanged.
#'
#' `cox_reference` does not change the group used by
#' `type = "difference"`. It applies only to observation-level Cox predictions
#' with `type = "lp"` or `type = "risk"`.
#'
#' Cox predictions with `type = "expected"` or `type = "survival"` also need
#' follow-up information for every prediction row.
#'
#' For a model fitted with a formula, include the variables used in the fitted
#' `Surv()` response. Alternatively, include exactly one column that is a
#' right-censored `Surv` object. Wrap it in `I()` when constructing a data frame,
#' for example:
#'
#' `data.frame(age = age, group = group, followup = I(survival::Surv(time, event)))`.
#'
#' Do not use `cox_reference` with `type = "expected"` or
#' `type = "survival"`.
#'
#' For a stratified Cox model, each prediction row must have valid stratum
#' information. If `strata()` was included in the fitted formula, include its
#' original variable or variables in `newdata`. Otherwise, supply new stratum
#' values through `strata`.
#'
#' @section Standard errors and confidence intervals:
#' For fitted curves and fitted differences, `se.fit = TRUE` returns pointwise
#' standard errors and confidence intervals with confidence level `level`.
#'
#' For predictions for individual observations, `se.fit = TRUE` requests the
#' standard errors supported by the underlying GLM or Cox prediction method.
#' The `level` argument is not used for these predictions.
#'
#' All standard errors and confidence intervals are conditional on the selected
#' term-specific model. They do not include uncertainty from adjustment-variable
#' selection, function selection, interaction selection, or selection among the
#' evaluated terms.
#'
#' @param object An object of class `"mfpi"` returned by [mfpi()].
#'
#' @param newdata An optional data frame or matrix containing values at which to
#'   predict.
#'
#'   For fitted curves, only the requested continuous variables are required.
#'   For predictions for individual observations, also supply the grouping
#'   variable, retained adjustment variables, and any required offset, strata,
#'   or Cox follow-up variables. See **Predictions for individual
#'   observations**.
#'
#' @param terms An optional character vector naming continuous variables examined
#'   by [mfpi()]. If `NULL`, all terms available under the selected `model`
#'   setting are used.
#'
#' @param model A character value selecting the available term-specific models.
#'   `"best"` uses only interactions retained by the MFPI selection criterion.
#'   `"all"` uses the stored model for every evaluated continuous variable.
#'   The default is `"best"`.
#'
#' @param type The prediction type. Use `"function"`, `"difference"`, or
#'   `"both"` for fitted curves. For predictions for individual observations,
#'   Gaussian, binomial, and Poisson models support `"link"` and `"response"`,
#'   while Cox models support `"lp"`, `"risk"`, `"expected"`, and
#'   `"survival"`. If `NULL`, `"both"` is used.
#'
#' @param se.fit A single `TRUE` or `FALSE` value indicating whether standard
#'   errors should be returned. The default is `TRUE`.
#'
#' @param level A number between 0 and 1 giving the confidence level for
#'   pointwise fitted-curve and fitted-difference intervals. The default is
#'   0.95. This argument is not used for predictions for individual
#'   observations.
#'
#' @param grid A single `TRUE` or `FALSE` value. If `TRUE`, fitted curves are
#'   evaluated on an equally spaced grid. If `FALSE`, supplied or fitted data
#'   values are used. This argument is ignored for predictions for individual
#'   observations.
#'
#' @param n_grid An integer of at least 2 giving the number of grid points when
#'   `grid = TRUE`. The default is 200. A structural-zero point can be added in
#'   addition to these grid points.
#'
#' @param strata Optional stratum information for predictions for individual
#'   observations from a stratified Cox model. For one stratification variable,
#'   supply one value per prediction row. For several variables, supply a matrix
#'   or data frame with one row per prediction row.
#'
#'   This argument is not used for non-Cox models, unstratified Cox models, or
#'   fitted-curve predictions.
#'
#' @param cox_reference A character value selecting the covariate reference for
#'   Cox observation-level predictions with `type = "lp"` or `type = "risk"`.
#'   Use `"zero"`, `"sample"`, or `"strata"`. The default is `"zero"`.
#'
#' @param newoffset An optional finite numeric vector containing one offset value
#'   per prediction row. It can be used only for predictions for individual
#'   observations from a term-specific model fitted with an offset.
#'
#' @param ... Reserved for future extensions. Currently, supplied arguments are
#'   ignored with a warning.
#'
#' @return
#' When one term is requested, an object of class `"mfpi_prediction"` is
#' returned. When several terms are requested, the result is a named list of
#' `"mfpi_prediction"` objects with class `"mfpi_prediction_list"`.
#'
#' Each result comes from the term-specific interaction model named by `term`.
#'
#' For `type = "function"` or `type = "both"`, the `functions` component is a
#' data frame containing:
#'
#' - `term`: the continuous-variable name;
#' - `x`: the value used to evaluate the fitted curve;
#' - `group`: the group label;
#' - `fit`: the fitted curve value;
#' - `se.fit`, `lower`, and `upper`: the pointwise standard error and confidence
#'   limits when `se.fit = TRUE`.
#'
#' For `type = "difference"` or `type = "both"`, the `differences` component is
#' a data frame containing:
#'
#' - `term`: the continuous-variable name;
#' - `x`: the value used to evaluate the curves;
#' - `contrast`: a label identifying the comparison;
#' - `group`: the comparison group;
#' - `reference`: the reference group;
#' - `fit`: the comparison-group curve minus the reference-group curve;
#' - `se.fit`, `lower`, and `upper`: the pointwise standard error and confidence
#'   limits when `se.fit = TRUE`.
#'
#' For predictions for individual observations, the `predictions` component is a
#' data frame containing `term` and `fit`, with `se.fit` added when requested and
#' available.
#'
#' Other components contain information used by plotting and internal
#' prediction checks and should not be treated as a stable public interface.
#'
#' @examples
#' \dontrun{
#' data("prostate")
#'
#' fit <- mfpi(
#'   lpsa ~ fp(age) + svi + fp(pgg45) + fp(cavol) + fp(weight) +
#'     fp(bph) + fp(cp),
#'   data = prostate,
#'   group_var = "svi",
#'   cont_vars = c("cavol", "age"),
#'   flex = "flex1",
#'   include_group_var = TRUE,
#'   center = FALSE,
#'   verbose = FALSE
#' )
#'
#' # By default, return fitted curves and differences for retained interactions.
#' predict(fit)
#'
#' # Inspect the stored model for an evaluated term even if it was not retained.
#' curves <- predict(
#'   fit,
#'   terms = "cavol",
#'   model = "all",
#'   type = "both",
#'   grid = TRUE
#' )
#'
#' curves$functions
#' curves$differences
#'
#' # Evaluate fitted curves at selected values. Only cavol is required here.
#' selected_values <- predict(
#'   fit,
#'   newdata = data.frame(cavol = c(0.2, 1, 3)),
#'   terms = "cavol",
#'   model = "all",
#'   type = "function",
#'   grid = FALSE
#' )
#'
#' selected_values$functions
#'
#' # Obtain one prediction per observation from the cavol interaction model.
#' row_predictions <- predict(
#'   fit,
#'   newdata = prostate[1:10, ],
#'   terms = "cavol",
#'   model = "all",
#'   type = "response"
#' )
#'
#' row_predictions$predictions
#'
#' # Several terms return a named list. Each element uses a different
#' # term-specific interaction model.
#' several_curves <- predict(
#'   fit,
#'   terms = c("cavol", "age"),
#'   model = "all",
#'   type = "function",
#'   grid = TRUE
#' )
#'
#' several_curves$cavol$functions
#' several_curves$age$functions
#'
#' # A formula offset is evaluated automatically from newdata.
#' set.seed(1)
#' d <- data.frame(
#'   y = rpois(100, 2),
#'   group = factor(sample(c("control", "treated"), 100, replace = TRUE)),
#'   age = runif(100, 30, 75),
#'   exposure = runif(100, 0.5, 3)
#' )
#'
#' fit_offset <- mfpi(
#'   y ~ age + group + offset(log(exposure)),
#'   data = d,
#'   family = "poisson",
#'   group_var = "group",
#'   cont_vars = "age",
#'   cont_var_forms = c(age = "linear"),
#'   flex = "flex1",
#'   p_interact = 1,
#'   verbose = FALSE
#' )
#'
#' predict(
#'   fit_offset,
#'   newdata = d[1:5, c("age", "group", "exposure")],
#'   terms = "age",
#'   model = "all",
#'   type = "response"
#' )
#'
#' # Cox predictions.
#' lung <- survival::lung
#' lung$status <- as.integer(lung$status == 2)
#' lung <- lung[complete.cases(lung[, c("time", "status", "age", "sex")]), ]
#' lung$sex <- factor(lung$sex)
#'
#' fit_cox <- mfpi(
#'   survival::Surv(time, status) ~ age + sex,
#'   data = lung,
#'   family = "cox",
#'   group_var = "sex",
#'   cont_vars = "age",
#'   cont_var_forms = c(age = "linear"),
#'   flex = "flex1",
#'   p_interact = 1,
#'   verbose = FALSE
#' )
#'
#' relative_data <- lung[1:5, c("age", "sex")]
#' predict(
#'   fit_cox,
#'   newdata = relative_data,
#'   terms = "age",
#'   model = "all",
#'   type = "risk"
#' )
#'
#' # Survival predictions also require the fitted follow-up variables.
#' absolute_data <- lung[1:5, c("time", "status", "age", "sex")]
#' predict(
#'   fit_cox,
#'   newdata = absolute_data,
#'   terms = "age",
#'   model = "all",
#'   type = "survival"
#' )
#' }
#'
#' @seealso [mfpi()], [plot.mfpi()]
#'
#' @method predict mfpi
#' @export
predict.mfpi <- function(object,
                         newdata = NULL,
                         terms = NULL,
                         model = c("best", "all"),
                         type = NULL,
                         se.fit = TRUE,
                         level = 0.95,
                         grid = FALSE,
                         n_grid = 200L,
                         strata = NULL,
                         cox_reference = NULL,
                         newoffset = NULL,
                         ...) {
  # Public S3 entry point. After argument validation, prediction is routed to
  # either the fitted-function path or the ordinary subject-level path.
  
  if (!inherits(object, "mfpi")) {
    stop("`object` must be an object of class \"mfpi\".", call. = FALSE)
  }
  
  model <- match.arg(model)
  family_string <- mfpi_family_string(object)
  type <- mfpi_match_prediction_type(type, family_string)
  
  function_types <- c("both", "function", "difference")
  ordinary_types <- if (identical(family_string, "cox")) {
    c("lp", "risk", "expected", "survival")
  } else {
    c("link", "response")
  }
  is_function_prediction <- type %in% function_types
  is_ordinary_prediction <- type %in% ordinary_types
  
  if (!is.logical(se.fit) || length(se.fit) != 1L || anyNA(se.fit)) {
    stop("`se.fit` must be a single non-missing logical value.",
         call. = FALSE)
  }
  
  if (!is.numeric(level) || length(level) != 1L || anyNA(level) ||
      !is.finite(level) || level <= 0 || level >= 1) {
    stop("`level` must be a single numeric value in (0, 1).",
         call. = FALSE)
  }
  
  if (!is.logical(grid) || length(grid) != 1L || anyNA(grid)) {
    stop("`grid` must be a single non-missing logical value.",
         call. = FALSE)
  }
  
  if (!is.numeric(n_grid) || length(n_grid) != 1L || anyNA(n_grid) ||
      !is.finite(n_grid) || n_grid < 2 || n_grid != as.integer(n_grid)) {
    stop("`n_grid` must be an integer >= 2.", call. = FALSE)
  }
  n_grid <- as.integer(n_grid)
  
  dots <- list(...)
  if (length(dots) > 0L) {
    dot_names <- names(dots)
    if (is.null(dot_names)) dot_names <- rep("<unnamed>", length(dots))
    dot_names[!nzchar(dot_names)] <- "<unnamed>"
    warning(
      "Unused arguments in `...`: ", paste(dot_names, collapse = ", "), ".",
      call. = FALSE
    )
  }
  
  if (is_ordinary_prediction && isTRUE(grid)) {
    warning("`grid` is ignored for ordinary subject-level prediction.",
            call. = FALSE)
  }
  
  if (is_function_prediction) {
    if (!is.null(strata)) {
      stop(
        "`strata` is not used for MFPI fitted-function predictions.",
        call. = FALSE
      )
    }
    if (!is.null(cox_reference)) {
      stop(
        "`cox_reference` is not used for MFPI fitted-function predictions.",
        call. = FALSE
      )
    }
    if (!is.null(newoffset)) {
      stop(
        "`newoffset` is not used for MFPI fitted-function predictions.",
        call. = FALSE
      )
    }
  }
  
  if (!identical(family_string, "cox")) {
    if (!is.null(strata)) {
      stop("`strata` is available only for Cox predictions.", call. = FALSE)
    }
    if (!is.null(cox_reference)) {
      stop("`cox_reference` is available only for Cox predictions.",
           call. = FALSE)
    }
  }
  
  if (identical(family_string, "cox")) {
    if (type %in% c("lp", "risk")) {
      cox_reference <- match_cox_reference(
        value = cox_reference,
        default = "zero",
        argument = "cox_reference"
      )
    } else if (!is.null(cox_reference)) {
      stop(
        "`cox_reference` applies only to Cox predictions with `type = \"lp\"` ",
        "or `type = \"risk\"`.",
        call. = FALSE
      )
    }
  }
  
  fits <- mfpi_get_prediction_fits(object = object, terms = terms, model = model)
  
  out <- stats::setNames(
    lapply(names(fits), function(term) {
      fit_result <- fits[[term]]
      
      if (is_function_prediction) {
        # Fitted-function targets evaluate every group-specific function at
        # each requested x value and return long-format function/contrast data.
        pred_data <- mfpi_prepare_prediction_data(
          object = object,
          term = term,
          newdata = newdata,
          grid = grid,
          n_grid = n_grid,
          purpose = "function"
        )
        
        basis <- mfpi_build_function_basis(
          object = object,
          term = term,
          fit_result = fit_result,
          cont_var_scaled = pred_data$cont_var_scaled,
          x_display = pred_data$x_display
        )
        
        return(mfpi_compute_function_prediction(
          object = object,
          term = term,
          fit_result = fit_result,
          basis = basis,
          type = type,
          se.fit = se.fit,
          level = level,
          reference = NULL
        ))
      }
      
      interaction_model <- fit_result$test_results$interaction_model
      fit_obj <- interaction_model$fit
      fit_is_stratified <- mfpi_fit_has_strata(fit_obj)
      
      if (!is.null(strata) && !fit_is_stratified) {
        stop(
          "`strata` was supplied, but the retained Cox interaction model for ",
          "term `", term, "` is not stratified.",
          call. = FALSE
        )
      }
      
      # Absolute Cox predictions need the response before raw newdata are
      # reduced to the predictors required by the term-specific design.
      cox_prediction_response <- if (
        !is.null(newdata) &&
        identical(family_string, "cox") &&
        type %in% c("expected", "survival")
      ) {
        reconstruct_cox_prediction_response(
          object = object,
          fit_obj = fit_obj,
          newdata = newdata
        )
      } else {
        NULL
      }
      
      # Ordinary prediction is subject-level: each row belongs to one group,
      # so the reconstructed interaction block is masked row by row. Formula-
      # level Cox strata are reconstructed from raw newdata unless explicitly
      # supplied.
      ordinary_strata <- strata
      if (is.null(ordinary_strata) && !is.null(newdata) &&
          !is.null(object$formula_strata_terms)) {
        ordinary_strata <- reconstruct_formula_strata_newdata(object, newdata)
      }
      
      design <- mfpi_build_ordinary_design(
        object = object,
        term = term,
        fit_result = fit_result,
        newdata = newdata,
        strata = ordinary_strata,
        newoffset = newoffset
      )
      
      if (identical(family_string, "cox") &&
          type %in% c("expected", "survival") &&
          !is.null(design$model_newdata)) {
        response <- if (is.null(newdata)) fit_obj$y else cox_prediction_response
        design$model_newdata <- attach_cox_prediction_response(
          fit_obj = fit_obj,
          newdata = design$model_newdata,
          response = response
        )
      }
      
      pred <- mfpi_predict_ordinary(
        object = object,
        term = term,
        fit_result = fit_result,
        model_newdata = design$model_newdata,
        type = type,
        se.fit = se.fit,
        cox_reference = cox_reference
      )
      
      pred_df <- data.frame(term = term, fit = pred$fit, stringsAsFactors = FALSE)
      if (se.fit && !is.null(pred$se.fit)) pred_df$se.fit <- pred$se.fit
      
      structure(
        list(
          term = term,
          type = type,
          predictions = pred_df,
          design = if (isTRUE(design$reconstructed)) design$X else NULL,
          metadata = list(
            ordinary_prediction = TRUE,
            newdata = isTRUE(design$newdata),
            design_matrix_reconstructed = isTRUE(design$reconstructed),
            group_levels = design$group_levels,
            group_internal = design$group_internal,
            offset_used = design$has_offset,
            response_scale = type %in% c("response", "risk"),
            absolute_cox_prediction = type %in% c("expected", "survival"),
            cox_reference = if (type %in% c("lp", "risk")) {
              cox_reference
            } else {
              NULL
            },
            used_model_predict = isTRUE(pred$used_model_predict),
            model_newdata_columns = if (!is.null(design$model_newdata)) {
              colnames(design$model_newdata)
            } else {
              NULL
            }
          )
        ),
        class = c("mfpi_prediction", "list")
      )
    }),
    names(fits)
  )
  
  if (length(out) == 1L) return(out[[1L]])
  structure(out, class = c("mfpi_prediction_list", "list"))
}

# -----------------------------------------------------------------------------
# Fit selection and prediction-data preparation --------------------------------
# -----------------------------------------------------------------------------

#' Select Term-Specific MFPI Fits for Prediction
#'
#' Resolves the terms requested by \code{predict.mfpi()} and returns the stored
#' flex fit result for each term. This is the only helper that decides whether
#' prediction is restricted to selected models (\code{model = "best"}) or may
#' use all stored term-specific winner models (\code{model = "all"}).
#'
#' @param object Object of class \code{"mfpi"} containing \code{var_winners},
#'   \code{best_interaction_model}, and term-specific fit results.
#' @param terms Character vector of requested continuous variables, or
#'   \code{NULL}. When \code{NULL}, all available terms in the requested
#'   \code{model} scope are selected.
#' @param model Character scalar, either \code{"best"} or \code{"all"}.
#'
#' @return A named list. Names are continuous-variable names. Each element is
#'   the stored term-specific flex result, usually \code{object$var_winners[[term]]$fit}.
#'
#' @keywords internal
#' @noRd
mfpi_get_prediction_fits <- function(object, terms, model) {
  # Decide which term-specific fitted interaction models are eligible for
  # prediction before any prediction matrices are constructed.
  
  if (model == "best") {
    available <- names(object$best_interaction_model)
  } else {
    winners <- object$var_winners
    available <- names(winners)[vapply(winners, function(w) {
      !is.null(w) && !is.null(w$fit) &&
        !is.null(w$fit$test_results$interaction_model)
    }, logical(1L))]
  }
  
  available <- available[!is.na(available) & nzchar(available)]
  
  if (length(available) == 0L) {
    stop(paste0("No ", model, " MFPI interaction models are available."),
         call. = FALSE)
  }
  
  if (!is.null(terms)) {
    if (!is.character(terms) || length(terms) == 0L || anyNA(terms) ||
        any(!nzchar(terms))) {
      stop("`terms` must be a non-empty character vector or NULL.",
           call. = FALSE)
    }
    
    missing_terms <- setdiff(terms, available)
    if (length(missing_terms) > 0L) {
      stop(
        paste0(
          "The following requested terms do not have ", model,
          " MFPI interaction models: ", paste(missing_terms, collapse = ", "), "."
        ),
        call. = FALSE
      )
    }
    available <- unique(terms)
  }
  
  out <- stats::setNames(vector("list", length(available)), available)
  for (term in available) {
    w <- object$var_winners[[term]]
    if (is.null(w) || is.null(w$fit)) {
      stop(paste0("No stored MFPI winner is available for term `", term, "`."),
           call. = FALSE)
    }
    out[[term]] <- w$fit
  }
  
  out
}


#' Prepare Fitted-Function Prediction Values
#'
#' Prepares the one-column continuous-variable matrix used to reconstruct
#' fitted functions. For supplied \code{newdata}, raw values are transformed
#' using the same shift, scale, and winsorisation metadata used at model fit
#' time. For \code{newdata = NULL}, the helper reads the stored
#' post-preprocessing training matrix from \code{object$x_train_internal}.
#'
#' @section Scale convention:
#' The returned \code{cont_var_scaled} is on the stored shifted/scaled scale,
#' \eqn{(x + shift) / scale}. The returned \code{x_display} is on the raw
#' display scale, \eqn{x}. Higher-level basis reconstruction multiplies
#' \code{cont_var_scaled} by \code{scale} before applying FP transformations.
#'
#' @section Grid construction:
#' If \code{grid = TRUE}, values are converted to the \eqn{x + shift} scale,
#' an equally spaced grid is built, and the grid is converted back to the stored
#' shifted/scaled scale. For zero-handled variables, the grid is based on the
#' positive part only, and a structural zero point is added when structural
#' zeros were present.
#'
#' @param object Object of class \code{"mfpi"} containing shift, scale,
#'   winsorisation, zero-handling, and stored training-matrix metadata.
#' @param term Character scalar naming the continuous variable to prepare.
#' @param newdata Optional matrix or data frame. If supplied, it must contain a
#'   column named \code{term} on the raw scale.
#' @param grid Logical scalar. Whether to replace observed or supplied values by
#'   a fitted-function evaluation grid.
#' @param n_grid Integer scalar. Number of positive-range grid points when
#'   \code{grid = TRUE}.
#' @param purpose Character scalar. Currently only \code{"function"}; retained
#'   to document that this helper is for fitted-function prediction, not full
#'   ordinary prediction.
#'
#' @return A list with \code{cont_var_scaled}, \code{x_display},
#'   \code{scale_var}, \code{shift_var}, and \code{zero_var}.
#'
#' @keywords internal
#' @noRd
mfpi_prepare_prediction_data <- function(object, term, newdata, grid, n_grid,
                                         purpose = "function") {
  # Fitted-function prediction needs only one continuous variable at a time.
  # This helper returns that variable on both the stored scaled scale and the
  # raw display scale used in returned tables and plots.
  
  if (!identical(purpose, "function")) {
    stop("Internal error: unsupported prediction-data preparation purpose.",
         call. = FALSE)
  }
  
  scale_var <- mfpi_named_scalar(object$scale, term, default = 1)
  shift_var <- mfpi_named_scalar(object$shift, term, default = 0)
  zero_var <- mfpi_get_zero_var(object, term)
  
  apply_winsor <- function(x_scaled_one, var) {
    lim <- object$winsorize_limits
    if (is.null(lim) || is.null(colnames(lim)) || !(var %in% colnames(lim))) {
      return(x_scaled_one)
    }
    lo <- lim["lower", var]
    hi <- lim["upper", var]
    if (is.finite(lo)) x_scaled_one[x_scaled_one < lo, 1L] <- lo
    if (is.finite(hi)) x_scaled_one[x_scaled_one > hi, 1L] <- hi
    x_scaled_one
  }
  
  as_numeric_matrix <- function(x, arg) {
    if (is.data.frame(x)) x <- data.matrix(x)
    if (is.vector(x) && is.null(dim(x))) x <- matrix(x, ncol = 1L)
    if (!is.matrix(x) || !is.numeric(x)) {
      stop(paste0("`", arg, "` must be a numeric matrix or data frame."),
           call. = FALSE)
    }
    if (is.null(colnames(x))) {
      stop(paste0("`", arg, "` must have column names."), call. = FALSE)
    }
    storage.mode(x) <- "double"
    x
  }
  
  if (is.null(newdata)) {
    x_train <- object$x_train_internal
    
    if (is.null(x_train)) {
      stop(
        paste0(
          "`newdata` is NULL, but the MFPI object does not contain ",
          "`x_train_internal`, the post-preprocessing training matrix. ",
          "Refit the model with the current version of `mfpi()` or supply ",
          "`newdata` containing `", term, "`."
        ),
        call. = FALSE
      )
    }
    
    if (!is.matrix(x_train) || is.null(colnames(x_train)) ||
        !(term %in% colnames(x_train))) {
      stop(paste0("Stored training matrix does not contain term `", term, "`."),
           call. = FALSE)
    }
    
    cont_scaled <- x_train[, term, drop = FALSE]
  } else {
    nd <- as_numeric_matrix(newdata, "newdata")
    if (!(term %in% colnames(nd))) {
      stop(paste0("`newdata` must contain a column named `", term, "`."),
           call. = FALSE)
    }
    x_raw <- nd[, term, drop = FALSE]
    if (anyNA(x_raw) || any(!is.finite(x_raw))) {
      stop(paste0("`newdata[[\"", term, "\"]]` must be finite."),
           call. = FALSE)
    }
    cont_scaled <- (x_raw + shift_var) / scale_var
    colnames(cont_scaled) <- term
    cont_scaled <- apply_winsor(cont_scaled, term)
  }
  
  cont_scaled <- as.matrix(cont_scaled)
  storage.mode(cont_scaled) <- "double"
  colnames(cont_scaled) <- term
  
  cont_plus <- cont_scaled * scale_var
  
  if (grid) {
    cont_vec <- as.vector(cont_plus)
    if (isTRUE(zero_var)) {
      zero_present <- any(cont_vec <= 0, na.rm = TRUE)
      grid_source <- cont_vec[cont_vec > 0 & is.finite(cont_vec)]
    } else {
      zero_present <- FALSE
      grid_source <- cont_vec[is.finite(cont_vec)]
    }
    
    if (length(grid_source) == 0L) {
      msg <- if (isTRUE(zero_var)) {
        "Cannot construct MFPI prediction grid: no finite positive values are available."
      } else {
        "Cannot construct MFPI prediction grid: no finite values are available."
      }
      stop(msg, call. = FALSE)
    }
    
    rg <- range(grid_source, na.rm = TRUE)
    grid_plus <- seq(rg[1L], rg[2L], length.out = n_grid)
    if (isTRUE(zero_var) && zero_present) grid_plus <- c(0, grid_plus)
    
    cont_scaled <- matrix(grid_plus / scale_var, ncol = 1L)
    colnames(cont_scaled) <- term
    cont_plus <- cont_scaled * scale_var
  }
  
  list(
    cont_var_scaled = cont_scaled,
    x_display = as.vector(cont_plus - shift_var),
    scale_var = scale_var,
    shift_var = shift_var,
    zero_var = zero_var
  )
}


# -----------------------------------------------------------------------------
# Fitted-function prediction ---------------------------------------------------
# -----------------------------------------------------------------------------

#' Reconstruct the Group-Specific FP Basis for One Term
#'
#' Builds the matrix used to compute MFPI fitted functions for one continuous
#' variable. The function converts stored shifted/scaled prediction values back
#' to the \eqn{x + shift} FP scale, delegates raw FP expansion and coefficient
#' name alignment to \code{build_group_fp_basis()}, applies fit-time centering
#' constants, and restores structural-zero rows to zero after centering.
#'
#' @section Centering contract:
#' Centering constants are read from the stored flex result. They are never
#' recomputed from prediction data. This is required because recomputing centers
#' would create a different design matrix than the one used to estimate the
#' fitted coefficients.
#'
#' @param object Object of class \code{"mfpi"} containing scale, shift, and
#'   zero-handling metadata.
#' @param term Character scalar naming the continuous variable being predicted.
#' @param fit_result Stored term-specific flex result. It must contain
#'   \code{bestfp_interaction}, \code{center_vals}, \code{coefficient_groups}
#'   or equivalent \code{xinteraction} attributes, and
#'   \code{test_results$interaction_model}.
#' @param cont_var_scaled One-column numeric matrix on the shifted/scaled scale,
#'   usually produced by \code{mfpi_prepare_prediction_data()}.
#' @param x_display Numeric vector of raw-scale x values to be stored in the
#'   prediction output.
#'
#' @return A list with \code{x} (centered group-specific FP basis),
#'   \code{x_display}, \code{group_fp_powers}, \code{coefficient_groups},
#'   \code{center_vals}, \code{centered}, \code{scale_var},
#'   \code{shift_var}, and \code{zero_var}.
#'
#' @keywords internal
#' @noRd
mfpi_build_function_basis <- function(object, term, fit_result,
                                      cont_var_scaled, x_display) {
  # Rebuild the group-specific FP block using only fit-time metadata. Centers
  # are reused exactly; they are never recomputed from prediction data.
  
  interaction_model <- fit_result$test_results$interaction_model
  if (is.null(interaction_model)) {
    stop(paste0("No interaction model is stored for term `", term, "`."),
         call. = FALSE)
  }
  
  group_fp_powers <- fit_result$bestfp_interaction
  if (!is.list(group_fp_powers) || length(group_fp_powers) == 0L) {
    stop(paste0("Term `", term, "` does not contain `bestfp_interaction`."),
         call. = FALSE)
  }
  
  coefficient_groups <- fit_result$coefficient_groups
  if (is.null(coefficient_groups)) {
    coefficient_groups <- attr(fit_result$xinteraction, "column_groups")
  }
  mfpi_validate_coefficient_groups(coefficient_groups, term = term)
  
  group_labels <- names(coefficient_groups)
  if (is.null(group_labels) || anyNA(group_labels) || any(!nzchar(group_labels))) {
    group_labels <- as.character(seq_along(coefficient_groups) - 1L)
    names(coefficient_groups) <- group_labels
  }
  
  if (length(group_fp_powers) != length(coefficient_groups)) {
    stop(
      paste0(
        "For term `", term, "`, `bestfp_interaction` length does not match ",
        "`coefficient_groups` length."
      ),
      call. = FALSE
    )
  }
  
  if (!is.matrix(cont_var_scaled) || ncol(cont_var_scaled) != 1L ||
      !is.numeric(cont_var_scaled) || anyNA(cont_var_scaled) ||
      any(!is.finite(cont_var_scaled))) {
    stop("`cont_var_scaled` must be a finite one-column numeric matrix.",
         call. = FALSE)
  }
  
  scale_var <- mfpi_named_scalar(object$scale, term, default = 1)
  shift_var <- mfpi_named_scalar(object$shift, term, default = 0)
  zero_var <- mfpi_get_zero_var(object, term)
  
  cont_plus <- cont_var_scaled * scale_var
  colnames(cont_plus) <- term
  
  transformed <- build_group_fp_basis(
    cont_mat = cont_plus,
    group_fp_powers = group_fp_powers,
    coefficient_groups = coefficient_groups,
    group_levels = group_labels,
    cont_name = term,
    transform_basis = TRUE,
    zero_var = zero_var,
    check_binary = FALSE
  )
  
  centers <- fit_result$center_vals
  centered <- !is.null(centers)
  
  if (centered) {
    if (!is.numeric(centers) || length(centers) == 0L) {
      stop(paste0("For term `", term, "`, `center_vals` must be numeric."),
           call. = FALSE)
    }
    
    if (is.null(names(centers))) {
      if (length(centers) != ncol(transformed)) {
        stop(
          paste0(
            "For term `", term, "`, unnamed `center_vals` length does not ",
            "match the reconstructed basis."
          ),
          call. = FALSE
        )
      }
      names(centers) <- colnames(transformed)
    } else {
      missing_centers <- setdiff(colnames(transformed), names(centers))
      if (length(missing_centers) > 0L) {
        stop(
          paste0(
            "For term `", term, "`, `center_vals` is missing constants for: ",
            paste(missing_centers, collapse = ", "), "."
          ),
          call. = FALSE
        )
      }
      centers <- centers[colnames(transformed)]
    }
    
    if (anyNA(centers) || any(!is.finite(centers))) {
      stop("`center_vals` must contain finite values.", call. = FALSE)
    }
    
    transformed <- sweep(transformed, 2L, centers, "-", check.margin = FALSE)
  } else {
    centers <- NULL
  }
  
  # Structural zeros must be restored after centering; otherwise zero rows would
  # incorrectly receive minus the centering constant as their FP contribution.
  zero_rows <- if (isTRUE(zero_var)) {
    as.vector(cont_plus) <= 0
  } else {
    rep(FALSE, nrow(cont_plus))
  }
  if (isTRUE(zero_var) && any(zero_rows)) transformed[zero_rows, ] <- 0
  
  if (anyNA(transformed) || any(!is.finite(transformed))) {
    stop(paste0("Non-finite values remain in the basis for term `", term, "`."),
         call. = FALSE)
  }
  
  list(
    x = transformed,
    x_display = x_display,
    group_fp_powers = attr(transformed, "group_fp_powers"),
    coefficient_groups = attr(transformed, "coefficient_groups"),
    center_vals = centers,
    centered = centered,
    scale_var = scale_var,
    shift_var = shift_var,
    zero_var = zero_var
  )
}


#' Map Internal MFPI Group Labels to User-Facing Labels
#'
#' Converts the internal group labels used for coefficient names, for example
#' \code{"0"} and \code{"1"}, to the original group labels supplied by the
#' user, for example \code{"Placebo"} and \code{"Active"}. The internal
#' labels must remain available for coefficient lookup, but prediction output
#' should use the original user-facing labels whenever possible.
#'
#' @param object Object of class \code{"mfpi"}.
#' @param group_labels Character vector of internal group labels in fitted-model
#'   order.
#'
#' @return Character vector of display labels with the same length and order as
#'   \code{group_labels}. If no stored mapping is available, \code{group_labels}
#'   is returned unchanged.
#'
#' @keywords internal
#' @noRd
mfpi_prediction_group_display_labels <- function(object, group_labels) {
  # Keep coefficient lookup on internal labels, but return user-facing labels
  # for long-format prediction tables whenever a stored mapping is available.
  
  group_labels <- as.character(group_labels)
  display_labels <- group_labels
  if (length(group_labels) == 0L) return(display_labels)
  
  level_map <- object$group_level_map
  if (is.data.frame(level_map) &&
      all(c("original", "internal") %in% names(level_map))) {
    internal <- as.character(level_map$internal)
    original <- as.character(level_map$original)
    hit <- match(group_labels, internal)
    mapped <- !is.na(hit)
    display_labels[mapped] <- original[hit[mapped]]
    return(display_labels)
  }
  
  original <- object$group_levels_original
  internal <- object$group_levels_new
  if (!is.null(original) && !is.null(internal) &&
      length(original) == length(internal)) {
    hit <- match(group_labels, as.character(internal))
    mapped <- !is.na(hit)
    display_labels[mapped] <- as.character(original)[hit[mapped]]
  }
  
  display_labels
}


#' Compute MFPI Fitted Functions and Differences from a Basis
#'
#' Computes group-specific fitted functions \eqn{\hat f_j(x)}, optional
#' pointwise standard errors and confidence intervals, and fitted-function
#' differences relative to a reference group.
#'
#' @section Statistical calculation:
#' For each group, the fitted function is computed by multiplying that group's
#' FP basis columns by their fitted coefficients and adding the model intercept
#' and group dummy offset when present. Pointwise fitted-function variances use
#' the corresponding submatrix of \code{vcov(interaction_model$fit)}.
#' Difference standard errors are computed analytically from the derivative of
#' \eqn{f_j(x) - f_r(x)} with respect to the fitted coefficients.
#'
#' @param object Object of class \code{"mfpi"}.
#' @param term Character scalar naming the continuous variable being predicted.
#' @param fit_result Stored term-specific flex fit result.
#' @param basis Output from \code{mfpi_build_function_basis()}.
#' @param type Character scalar: \code{"function"}, \code{"difference"}, or
#'   \code{"both"}.
#' @param se.fit Logical scalar. Whether to compute standard errors and
#'   pointwise intervals.
#' @param level Numeric scalar in \code{(0, 1)} giving the pointwise confidence
#'   level.
#' @param reference Optional reference group. May be an internal group label or
#'   an original group label stored in the MFPI object.
#'
#' @return An object of class \code{"mfpi_prediction"} containing
#'   long-format \code{functions}, long-format \code{differences}, and
#'   \code{metadata} components according to \code{type}.
#'
#' @keywords internal
#' @noRd
mfpi_compute_function_prediction <- function(object, term, fit_result, basis,
                                             type, se.fit, level, reference) {
  # Convert a reconstructed FP basis into long-format fitted functions and
  # long-format group contrasts. No historical wide matrix is built here.
  
  interaction_model <- fit_result$test_results$interaction_model
  coef_vec <- interaction_model$coefficients
  if (!is.numeric(coef_vec) || is.null(names(coef_vec))) {
    stop("Interaction-model coefficients must be a named numeric vector.",
         call. = FALSE)
  }
  coef_vec <- stats::setNames(as.numeric(coef_vec), names(coef_vec))
  vcov_mat <- stats::vcov(interaction_model$fit)
  
  family_string <- mfpi_family_string(object)
  group_labels <- names(basis$coefficient_groups)
  group_display_labels <- mfpi_prediction_group_display_labels(
    object = object,
    group_labels = group_labels
  )
  k <- length(group_labels)
  n <- nrow(basis$x)
  crit <- stats::qnorm((1 + level) / 2)
  
  ref_pos <- mfpi_resolve_reference_group(
    reference = reference,
    group_labels = group_labels,
    object = object
  )
  
  has_intercept <- !identical(family_string, "cox") &&
    "(Intercept)" %in% names(coef_vec)
  intercept <- if (has_intercept) unname(coef_vec["(Intercept)"]) else 0
  
  group_name <- object$group_var
  group_dummy_names <- if (k > 1L) paste0(group_name, group_labels[-1L]) else character(0L)
  
  required <- unique(c(
    if (has_intercept) "(Intercept)",
    unlist(basis$coefficient_groups, use.names = FALSE),
    group_dummy_names
  ))
  mfpi_validate_required_coefficients(
    coef_vec = coef_vec,
    vcov_mat = vcov_mat,
    required_names = required,
    context = paste0("MFPI fitted-function prediction for term `", term, "`")
  )
  
  fit_mat <- matrix(NA_real_, nrow = n, ncol = k)
  se_mat <- if (se.fit) matrix(NA_real_, nrow = n, ncol = k) else NULL
  colnames(fit_mat) <- paste0("f", group_labels)
  if (se.fit) colnames(se_mat) <- paste0("se(f", group_labels, ")")
  
  for (g in seq_len(k)) {
    cols_g <- basis$coefficient_groups[[g]]
    x_g <- basis$x[, cols_g, drop = FALSE]
    beta_g <- coef_vec[cols_g]
    
    x_var <- x_g
    x_var_names <- cols_g
    dummy_offset <- 0
    
    if (has_intercept) {
      x_var <- cbind("(Intercept)" = 1, x_var)
      x_var_names <- c("(Intercept)", x_var_names)
    }
    
    if (g > 1L) {
      dummy_name <- group_dummy_names[g - 1L]
      dummy_offset <- unname(coef_vec[dummy_name])
      x_var <- cbind(x_var, 1)
      colnames(x_var)[ncol(x_var)] <- dummy_name
      x_var_names <- c(x_var_names, dummy_name)
    }
    
    fit_mat[, g] <- intercept + as.vector(x_g %*% beta_g) + dummy_offset
    
    if (se.fit) {
      v_sub <- vcov_mat[x_var_names, x_var_names, drop = FALSE]
      var_g <- rowSums((x_var %*% v_sub) * x_var)
      se_mat[, g] <- sqrt(mfpi_sanitize_variance(
        var_g,
        context = paste0("fitted-function SE for term `", term,
                         "`, group `", group_labels[g], "`")
      ))
    }
  }
  
  # Long fitted-function output.
  function_rows <- vector("list", k)
  for (g in seq_len(k)) {
    df_g <- data.frame(
      term = term,
      x = basis$x_display,
      group = group_display_labels[g],
      fit = fit_mat[, g],
      stringsAsFactors = FALSE
    )
    if (se.fit) {
      df_g$se.fit <- se_mat[, g]
      df_g$lower <- df_g$fit - crit * df_g$se.fit
      df_g$upper <- df_g$fit + crit * df_g$se.fit
    }
    function_rows[[g]] <- df_g
  }
  functions_df <- do.call(rbind, function_rows)
  
  # Long fitted-function-difference output. The returned table is
  # intentionally long-format rather than matrix-shaped: group and reference
  # labels remain explicit data columns instead of being encoded in column names
  # such as f1-f0.
  compare_pos <- setdiff(seq_len(k), ref_pos)
  diff_rows <- vector("list", length(compare_pos))
  
  if (length(compare_pos) > 0L) {
    contrast_labels <- paste0(group_labels[compare_pos], "-", group_labels[ref_pos])
    contrast_display_labels <- paste0(
      group_display_labels[compare_pos],
      "-",
      group_display_labels[ref_pos]
    )
    
    for (ii in seq_along(compare_pos)) {
      g <- compare_pos[ii]
      diff_vec <- fit_mat[, g] - fit_mat[, ref_pos]
      
      df_d <- data.frame(
        term = term,
        x = basis$x_display,
        contrast = contrast_display_labels[ii],
        group = group_display_labels[g],
        reference = group_display_labels[ref_pos],
        fit = diff_vec,
        stringsAsFactors = FALSE
      )
      
      if (se.fit) {
        ref_cols <- basis$coefficient_groups[[ref_pos]]
        grp_cols <- basis$coefficient_groups[[g]]
        
        # Derivative of f_g(x) - f_ref(x) with respect to the FP coefficients:
        # -X_ref(x) for the reference group and +X_g(x) for the comparison group.
        d_mat <- cbind(
          -basis$x[, ref_cols, drop = FALSE],
          basis$x[, grp_cols, drop = FALSE]
        )
        d_names <- c(ref_cols, grp_cols)
        
        # Add derivatives for group-dummy offsets. These terms are needed when
        # the comparison group or chosen reference group is not the model's
        # baseline group.
        dummy_terms <- character(0L)
        dummy_values <- numeric(0L)
        if (g > 1L) {
          dummy_terms <- c(dummy_terms, group_dummy_names[g - 1L])
          dummy_values <- c(dummy_values, 1)
        }
        if (ref_pos > 1L) {
          dummy_terms <- c(dummy_terms, group_dummy_names[ref_pos - 1L])
          dummy_values <- c(dummy_values, -1)
        }
        if (length(dummy_terms) > 0L) {
          for (j in seq_along(dummy_terms)) {
            d_mat <- cbind(d_mat, rep(dummy_values[j], nrow(d_mat)))
            d_names <- c(d_names, dummy_terms[j])
          }
        }
        colnames(d_mat) <- d_names
        
        mfpi_validate_required_coefficients(
          coef_vec = coef_vec,
          vcov_mat = vcov_mat,
          required_names = colnames(d_mat),
          context = paste0("MFPI difference prediction for term `", term,
                           "`, contrast `", contrast_labels[ii], "`")
        )
        
        v_sub <- vcov_mat[colnames(d_mat), colnames(d_mat), drop = FALSE]
        var_d <- rowSums((d_mat %*% v_sub) * d_mat)
        se_d <- sqrt(mfpi_sanitize_variance(
          var_d,
          context = paste0("difference SE for term `", term,
                           "`, contrast `", contrast_labels[ii], "`")
        ))
        df_d$se.fit <- se_d
        df_d$lower <- df_d$fit - crit * df_d$se.fit
        df_d$upper <- df_d$fit + crit * df_d$se.fit
      }
      
      diff_rows[[ii]] <- df_d
    }
  }
  
  differences_df <- if (length(diff_rows) > 0L) do.call(rbind, diff_rows) else data.frame()
  
  structure(
    list(
      term = term,
      type = type,
      x = basis$x_display,
      # Canonical fitted-function output. Long format is deliberate because it
      # scales to multiple groups and keeps group labels as data, not column names.
      functions = if (type == "difference") NULL else functions_df,
      # Canonical fitted-difference output. Long format keeps the comparison
      # group and reference group explicit for plotting and downstream use.
      differences = if (type == "function") NULL else differences_df,
      metadata = list(
        group_fp_powers = basis$group_fp_powers,
        center_vals = basis$center_vals,
        coefficient_groups = basis$coefficient_groups,
        group_display_labels = stats::setNames(group_display_labels, group_labels),
        reference = group_labels[ref_pos],
        reference_label = group_display_labels[ref_pos],
        scale_var = basis$scale_var,
        shift_var = basis$shift_var,
        zero_var = basis$zero_var,
        level = level,
        se.fit = se.fit
      )
    ),
    class = c("mfpi_prediction", "list")
  )
}


# -----------------------------------------------------------------------------
# Ordinary subject-level prediction -------------------------------------------
# -----------------------------------------------------------------------------

#' Resolve the Adjustment Design Required by One MFPI Interaction Model
#'
#' @keywords internal
#' @noRd
mfpi_prediction_adjustment_spec <- function(object, term, fit_result) {
  coefficient_groups <- fit_result$coefficient_groups
  if (is.null(coefficient_groups)) {
    coefficient_groups <- attr(fit_result$xinteraction, "column_groups")
  }
  mfpi_validate_coefficient_groups(coefficient_groups, term = term)
  
  group_labels <- names(coefficient_groups)
  if (is.null(group_labels) || anyNA(group_labels) || any(!nzchar(group_labels))) {
    group_labels <- as.character(seq_along(coefficient_groups) - 1L)
    names(coefficient_groups) <- group_labels
  }
  
  adj_model <- object$adjustment_model
  adj_terms <- character(0L)
  if (!is.null(adj_model)) {
    adj_terms <- get_selected_variables(adj_model)
  }
  
  adjustment_lookup <- if (!is.null(adj_model$term_to_columns)) {
    adj_model$term_to_columns
  } else if (!is.null(object$adjustment_term_to_columns)) {
    object$adjustment_term_to_columns
  } else {
    stats::setNames(as.list(adj_terms), adj_terms)
  }
  
  group_var <- object$group_var
  group_dummy_names <- paste0(group_var, group_labels[-1L])
  adj_terms <- setdiff(adj_terms, c(term, group_var, group_dummy_names))
  adj_terms <- unique(adj_terms[!is.na(adj_terms) & nzchar(adj_terms)])
  
  missing_terms <- setdiff(adj_terms, names(adjustment_lookup))
  if (length(missing_terms) > 0L) {
    stop(
      paste0(
        "The stored adjustment term-to-column lookup is missing selected term(s): ",
        paste(missing_terms, collapse = ", "), "."
      ),
      call. = FALSE
    )
  }
  
  raw_adj_cols <- unique(unlist(
    adjustment_lookup[adj_terms],
    use.names = FALSE
  ))
  
  list(
    coefficient_groups = coefficient_groups,
    group_labels = group_labels,
    adjustment_model = adj_model,
    adjustment_lookup = adjustment_lookup,
    adjustment_terms = adj_terms,
    raw_adjustment_columns = raw_adj_cols
  )
}

#' Build the Ordinary Interaction-Model Design Matrix
#'
#' Reconstructs the numeric design matrix required for ordinary GLM or Cox
#' subject-level prediction for one MFPI term. This helper is used only for subject-level prediction. It differs from
#' fitted-function prediction because each new subject belongs to exactly one
#' group, so out-of-group interaction blocks must be set to zero row by row.
#'
#' @section Training-data path:
#' If \code{newdata = NULL}, \code{newoffset = NULL}, and no replacement
#' \code{strata} are supplied, no matrix is reconstructed. The stored model's
#' own \code{predict()} method is called directly, preserving native GLM and
#' Cox prediction behavior. If \code{newdata = NULL} but a
#' replacement \code{newoffset} is supplied, the stored shifted/scaled training
#' matrix is used to reconstruct the design matrix so the supplied new offset can
#' be applied explicitly.
#'
#' @section New-data path:
#' For supplied \code{newdata}, the helper first reconstructs only formula
#' terms and factor contrasts required by the stored interaction model, then
#' maps the grouping variable to
#' internal levels, creates group dummies, rebuilds the group-specific FP
#' interaction block, and masks out-of-group blocks row by row. Selected
#' adjustment terms are expanded through the stored term-to-column lookup so
#' categorical contrast blocks remain complete. The helper then applies the
#' stored adjustment transformations, resolves offsets, and orders columns to
#' match the fitted coefficient vector.
#'
#' @param object Object of class \code{"mfpi"} containing group coding,
#'   preprocessing, adjustment-model, and formula metadata.
#' @param term Character scalar naming the continuous interaction variable.
#' @param fit_result Stored term-specific flex fit result.
#' @param newdata Optional prediction data. If \code{NULL}, the fitted model's
#'   own prediction method is used. If supplied, it must contain the interaction
#'   term, grouping variable, and selected adjustment terms. Formula fits accept
#'   original factor columns; matrix fits require every raw column in each
#'   selected grouped term.
#' @param newoffset Optional numeric offset vector for ordinary prediction. For
#'   supplied \code{newdata}, it must contain one value per new-data row. For
#'   \code{newdata = NULL}, it must contain one value per training row and
#'   triggers reconstruction of the stored training design. Required when the
#'   fitted interaction model used an offset and reconstructed prediction data
#'   are needed.
#'
#' @return A list with \code{X} (design matrix without intercept, or
#'   \code{NULL} when prediction delegates to the fitted model),
#'   \code{offset}, \code{has_offset}, \code{group_internal},
#'   \code{group_levels}, \code{newdata}, and \code{reconstructed}.
#'
#' @keywords internal
#' @noRd
mfpi_build_ordinary_design <- function(object, term, fit_result, newdata,
                                       strata = NULL, newoffset = NULL) {
  # Build formula-compatible subject-level prediction data only when new rows,
  # replacement strata, or a replacement offset must be supplied to the stored
  # fitted model.
  
  interaction_model <- fit_result$test_results$interaction_model
  if (is.null(interaction_model) || is.null(interaction_model$fit)) {
    stop(paste0("No fitted interaction model is stored for term `", term, "`."),
         call. = FALSE)
  }
  
  coef_vec <- interaction_model$coefficients
  if (!is.numeric(coef_vec) || is.null(names(coef_vec))) {
    stop("Interaction-model coefficients must be a named numeric vector.",
         call. = FALSE)
  }
  
  # Training-data ordinary prediction usually delegates to the fitted model.
  # This keeps Cox and GLM fitted-value conventions identical to the original
  # model object when no replacement newoffset is requested.
  #
  # If a new `newoffset` or replacement `strata` are supplied with
  # `newdata = NULL`, reconstruct formula-compatible prediction data from
  # `x_train_internal` so the stored model can consume those replacements.
  newdata_is_training_internal <- FALSE
  if (is.null(newdata)) {
    if (is.null(newoffset) && is.null(strata)) {
      return(list(
        X = NULL,
        offset = NULL,
        model_newdata = NULL,
        has_offset = mfpi_fit_has_offset(interaction_model$fit),
        group_internal = NULL,
        group_levels = NULL,
        newdata = FALSE,
        reconstructed = FALSE
      ))
    }
    
    newdata <- object$x_train_internal
    
    if (is.null(newdata)) {
      stop(
        paste0(
          "Cannot reconstruct training prediction data for replacement ",
          "`newoffset` or `strata` because the MFPI object does not store ",
          "`x_train_internal`, the post-preprocessing training matrix. Refit ",
          "the model with the current version of `mfpi()` or supply explicit ",
          "`newdata`."
        ),
        call. = FALSE
      )
    }
    
    newdata_is_training_internal <- TRUE
  }
  
  coef_names <- names(coef_vec)
  target_cols <- setdiff(coef_names, "(Intercept)")
  
  # Resolve the selected adjustment recipe before formula reconstruction.
  # Prediction should evaluate only formula expressions needed by this stored
  # term-specific interaction model, not every predictor offered at fit time.
  prediction_spec <- mfpi_prediction_adjustment_spec(
    object = object,
    term = term,
    fit_result = fit_result
  )
  coefficient_groups <- prediction_spec$coefficient_groups
  group_labels <- prediction_spec$group_labels
  adj_model <- prediction_spec$adjustment_model
  adjustment_lookup <- prediction_spec$adjustment_lookup
  adj_terms <- prediction_spec$adjustment_terms
  raw_adj_cols <- prediction_spec$raw_adjustment_columns
  
  n_raw <- if (is.data.frame(newdata)) nrow(newdata) else NROW(newdata)
  if (n_raw == 0L) {
    stop("`newdata` must contain at least one row.", call. = FALSE)
  }
  
  # Keep the raw user-supplied data frame for formula-level offset
  # reconstruction. Formula design reconstruction may add/replace columns, but
  # offset expressions such as offset(log(exposure)) must be evaluated against
  # the original newdata variables.
  newdata_raw <- newdata
  
  # Formula-interface prediction should replay the fit-time formula design
  # recipe before numeric coercion. This rebuilds formula-derived columns such
  # as factor dummies using stored terms, contrasts, and xlevels, while preserving
  # raw group and continuous variables needed for MFPI-specific reconstruction.
  required_input_columns <- unique(c(term, object$group_var, raw_adj_cols))
  
  if (!isTRUE(newdata_is_training_internal) &&
      isTRUE(object$formula_interface) &&
      is.data.frame(newdata)) {
    newdata <- mfpi_prepare_formula_newdata(
      object = object,
      newdata = newdata,
      terms = unique(c(term, object$group_var, adj_terms)),
      required_columns = required_input_columns
    )
  }
  
  # Ignore extra supplied columns before numeric validation. In particular,
  # missing values or unseen levels in formula predictors that are absent from
  # this selected interaction model must not affect prediction.
  if (is.data.frame(newdata)) {
    keep_input <- intersect(required_input_columns, names(newdata))
    newdata <- newdata[, keep_input, drop = FALSE]
  } else if (is.matrix(newdata) && !is.null(colnames(newdata))) {
    keep_input <- intersect(required_input_columns, colnames(newdata))
    newdata <- newdata[, keep_input, drop = FALSE]
  }
  
  # Local conversion helper. At this point, formula-interface data frames have
  # already been expanded with fit-time formula metadata. Plain matrix-interface
  # input remains unchanged.
  as_prediction_matrix <- function(data) {
    if (is.data.frame(data)) data <- data.matrix(data)
    if (is.vector(data) && is.null(dim(data))) data <- matrix(data, ncol = 1L)
    
    if (!is.matrix(data) || !is.numeric(data)) {
      stop("`newdata` must be a numeric matrix or coercible data frame.",
           call. = FALSE)
    }
    
    if (is.null(colnames(data))) {
      stop("`newdata` must have column names.", call. = FALSE)
    }
    
    if (anyNA(data) || any(!is.finite(data))) {
      stop("`newdata` must contain only finite numeric values.", call. = FALSE)
    }
    
    storage.mode(data) <- "double"
    data
  }
  
  nd <- as_prediction_matrix(newdata)
  n <- nrow(nd)
  
  if (is.null(newoffset) &&
      !isTRUE(newdata_is_training_internal) &&
      isTRUE(object$formula_interface) &&
      is.data.frame(newdata_raw) &&
      !is.null(object$formula_offset_terms)) {
    newoffset <- reconstruct_formula_offset_newdata(
      object = object,
      newdata = newdata_raw
    )
  }
  
  # Local scaling helper applies the same raw -> shifted/scaled transformation
  # used at fit time and then applies stored winsorisation limits if present.
  scale_vars <- function(vars) {
    vars <- unique(vars[!is.na(vars) & nzchar(vars)])
    if (length(vars) == 0L) return(matrix(nrow = n, ncol = 0L))
    
    missing_vars <- setdiff(vars, colnames(nd))
    if (length(missing_vars) > 0L) {
      stop(
        paste0("`newdata` is missing required predictor column(s): ",
               paste(missing_vars, collapse = ", "), "."),
        call. = FALSE
      )
    }
    
    # When a replacement newoffset is supplied with `newdata = NULL`, `nd` is
    # `x_train_internal`. It has already been shifted, scaled, winsorised, and
    # internally group-remapped at fit time. Do not transform it a second time.
    if (isTRUE(newdata_is_training_internal)) {
      out <- nd[, vars, drop = FALSE]
      storage.mode(out) <- "double"
      return(out)
    }
    
    out <- nd[, vars, drop = FALSE]
    for (v in vars) {
      shift_v <- mfpi_named_scalar(object$shift, v, default = 0)
      scale_v <- mfpi_named_scalar(object$scale, v, default = 1)
      if (scale_v <= 0) {
        stop(paste0("Stored scale for predictor `", v, "` must be positive."),
             call. = FALSE)
      }
      out[, v] <- (out[, v] + shift_v) / scale_v
      
      lim <- object$winsorize_limits
      if (!is.null(lim) && !is.null(colnames(lim)) && v %in% colnames(lim)) {
        lo <- lim["lower", v]
        hi <- lim["upper", v]
        if (is.finite(lo)) out[out[, v] < lo, v] <- lo
        if (is.finite(hi)) out[out[, v] > hi, v] <- hi
      }
    }
    storage.mode(out) <- "double"
    out
  }
  
  # Group values are read from the original data frame when available to avoid
  # losing factor labels during numeric matrix conversion.
  group_var <- object$group_var
  if (is.null(group_var) || length(group_var) != 1L || !nzchar(group_var)) {
    stop("The MFPI object does not store a valid `group_var`.", call. = FALSE)
  }
  
  if (is.data.frame(newdata) && group_var %in% names(newdata)) {
    group_raw <- as.vector(newdata[[group_var]])
  } else {
    if (!(group_var %in% colnames(nd))) {
      stop(paste0("`newdata` must contain the grouping variable `", group_var, "`."),
           call. = FALSE)
    }
    group_raw <- as.vector(nd[, group_var])
  }
  
  if (isTRUE(newdata_is_training_internal)) {
    # `x_train_internal` already contains MFPI's internal numeric group codes.
    # Do not match these values against original labels such as "Placebo" and
    # "Active".
    group_internal <- suppressWarnings(as.numeric(group_raw))
    if (anyNA(group_internal)) {
      stop(
        "Stored training group values must be numeric internal group codes.",
        call. = FALSE
      )
    }
  } else {
    orig <- object$group_levels_original
    new <- object$group_levels_new
    if (!is.null(orig) && !is.null(new) && length(orig) == length(new)) {
      hit <- match(as.character(group_raw), as.character(orig))
      if (anyNA(hit)) {
        bad <- unique(group_raw[is.na(hit)])
        stop(
          paste0(
            "`newdata` contains group level(s) not seen at fit time: ",
            paste(bad, collapse = ", "),
            ". Ordinary prediction uses the group coding and reference level stored ",
            "at model fitting; no coefficients exist for new group levels."
          ),
          call. = FALSE
        )
      }
      group_internal <- as.numeric(new[hit])
    } else {
      group_internal <- suppressWarnings(as.numeric(group_raw))
      if (anyNA(group_internal)) {
        stop("`newdata` group values must match stored group levels or be numeric internal levels.",
             call. = FALSE)
      }
    }
  }
  
  group_levels_numeric <- suppressWarnings(as.numeric(group_labels))
  if (anyNA(group_levels_numeric)) {
    stop(
      paste0("Internal group labels for term `", term,
             "` must be numeric-like for dummy reconstruction."),
      call. = FALSE
    )
  }
  
  # Numeric fallback only proves that `newdata` group values can be interpreted
  # as internal numeric group codes. It does not prove that those codes were
  # observed when the model was fitted.
  #
  # New group levels must be rejected because the fitted MFPI model has no
  # group-specific interaction coefficients for unseen groups.
  unknown_group <- !group_internal %in% group_levels_numeric
  
  if (any(unknown_group)) {
    unknown_values <- unique(group_internal[unknown_group])
    
    stop(
      paste0(
        "`newdata` contains group level(s) not seen at fit time: ",
        paste(unknown_values, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }
  
  group_mat <- matrix(group_internal, ncol = 1L)
  colnames(group_mat) <- group_var
  group_dummies <- create_group_dummies(
    group_mat, 
    levels = group_levels_numeric,
    quiet  = TRUE
  )
  
  # Interaction block: first evaluate all group-specific functions at each row's
  # x value, then zero out FP blocks that do not correspond to the row's group.
  if (!(term %in% colnames(nd))) {
    stop(paste0("`newdata` must contain the interaction term `", term, "`."),
         call. = FALSE)
  }
  
  term_scaled <- scale_vars(term)
  scale_var <- mfpi_named_scalar(object$scale, term, default = 1)
  shift_var <- mfpi_named_scalar(object$shift, term, default = 0)
  x_display <- as.vector(term_scaled[, term, drop = FALSE] * scale_var - shift_var)
  
  basis <- mfpi_build_function_basis(
    object = object,
    term = term,
    fit_result = fit_result,
    cont_var_scaled = term_scaled[, term, drop = FALSE],
    x_display = x_display
  )
  
  z_interaction <- basis$x
  group_chr <- as.character(group_internal)
  for (g in seq_along(basis$coefficient_groups)) {
    rows_outside_g <- group_chr != group_labels[g]
    if (any(rows_outside_g)) {
      z_interaction[rows_outside_g, basis$coefficient_groups[[g]]] <- 0
    }
  }
  
  # Adjustment block. The adjustment model stores one row per conceptual term,
  # while prediction data contain the raw design columns. Expand selected terms
  # through term_to_columns so categorical contrast blocks are reconstructed and
  # centered together exactly as they were during MFPI fitting.
  # Formula newdata has already been expanded with the fit-time terms,
  # contrasts, and factor levels. Every member column of each selected grouped
  # adjustment term must therefore be present before transformation.
  if (!isTRUE(newdata_is_training_internal) &&
      isTRUE(object$formula_interface) &&
      length(raw_adj_cols) > 0L) {
    missing_adj <- setdiff(raw_adj_cols, colnames(nd))
    if (length(missing_adj) > 0L) {
      stop(
        paste0(
          "Formula-interface prediction could not reconstruct selected ",
          "adjustment column(s): ",
          paste(missing_adj, collapse = ", "), "."
        ),
        call. = FALSE
      )
    }
  }
  
  x_adjustment <- NULL
  if (length(adj_terms) > 0L) {
    if (is.null(adj_model)) {
      stop("Selected adjustment variables require `object$adjustment_model`.",
           call. = FALSE)
    }
    
    x_scaled <- scale_vars(raw_adj_cols)
    x_plus <- x_scaled
    for (col in raw_adj_cols) {
      x_plus[, col] <- x_plus[, col] *
        mfpi_named_scalar(object$scale, col, default = 1)
    }
    
    powers_term <- if (!is.null(adj_model$fp_powers)) {
      adj_model$fp_powers[adj_terms]
    } else if (!is.null(adj_model$fp_terms)) {
      get_fp_powers(adj_terms, adj_model$fp_terms)
    } else {
      NULL
    }
    
    if (is.null(powers_term) || length(powers_term) != length(adj_terms) ||
        any(vapply(powers_term, length, integer(1L)) == 0L)) {
      stop("Adjustment model is missing FP powers required for prediction.",
           call. = FALSE)
    }
    
    expanded_adjustment <- expand_term_metadata_to_columns(
      term_to_columns = adjustment_lookup[adj_terms],
      powers = powers_term,
      raw_columns = raw_adj_cols,
      acdx = adj_model$acd,
      zero = adj_model$zero,
      catzero = adj_model$catzero,
      spike = adj_model$spike,
      spike_decision = adj_model$spike_dec,
      acd_parameter = adj_model$acd_parameter
    )
    power_cols <- expanded_adjustment$powers
    acd_cols <- expanded_adjustment$acdx
    zero_cols <- expanded_adjustment$zero
    catzero_cols <- expanded_adjustment$catzero
    spike_cols <- expanded_adjustment$spike
    spike_decision_cols <- expanded_adjustment$spike_decision
    acd_parameter_cols <- expanded_adjustment$acd_parameter
    
    x_trans <- transform_matrix(
      x = x_plus,
      power_list = power_cols,
      center = stats::setNames(rep(FALSE, length(raw_adj_cols)), raw_adj_cols),
      acdx = acd_cols,
      acd_parameter_list = acd_parameter_cols,
      zero = zero_cols,
      catzero = catzero_cols,
      spike = spike_cols,
      spike_decision = spike_decision_cols,
      keep_x_order = FALSE,
      reset_zero = FALSE,
      check_binary = TRUE
    )
    
    if (!is.null(x_trans) && !is.null(x_trans$x_transformed) &&
        ncol(x_trans$x_transformed) > 0L) {
      x_adjustment <- as.matrix(x_trans$x_transformed)
      storage.mode(x_adjustment) <- "double"
      
      centers <- adj_model$centers
      if (!is.null(centers)) {
        missing_centers <- setdiff(colnames(x_adjustment), names(centers))
        if (length(missing_centers) > 0L) {
          stop(
            paste0(
              "Adjustment-model centering constants are missing for: ",
              paste(missing_centers, collapse = ", "), "."
            ),
            call. = FALSE
          )
        }
        
        zero_expanded <- x_trans$zero_expanded[colnames(x_adjustment)]
        if (is.null(zero_expanded) || anyNA(zero_expanded)) {
          zero_expanded <- stats::setNames(
            rep(FALSE, ncol(x_adjustment)),
            colnames(x_adjustment)
          )
        }
        
        x_adjustment <- center_matrix(
          mat = x_adjustment,
          centers = centers[colnames(x_adjustment)],
          zero = zero_expanded
        )
      }
      
      if (anyNA(x_adjustment) || any(!is.finite(x_adjustment))) {
        stop("Non-finite values were produced in the adjustment prediction matrix.",
             call. = FALSE)
      }
    }
  }
  
  pieces <- list(group_dummies, z_interaction)
  if (!is.null(x_adjustment) && ncol(x_adjustment) > 0L) {
    pieces <- c(pieces, list(x_adjustment))
  }
  
  X <- do.call(cbind, pieces)
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  
  if (is.null(colnames(X)) || anyNA(colnames(X)) || any(!nzchar(colnames(X)))) {
    stop("Reconstructed prediction design has missing column names.",
         call. = FALSE)
  }
  
  # Resolve source design columns to the exact names used by the fitted model.
  # Formula-based fits may quote non-syntactic coefficient names, whereas
  # predict.glm()/predict.coxph() still require raw source-column names in
  # newdata. Keep both representations rather than rewriting names heuristically.
  column_map <- interaction_model$transformed_to_model_columns
  if (is.null(column_map) || !is.character(column_map) ||
      is.null(names(column_map)) || anyNA(column_map) ||
      any(!nzchar(column_map)) || anyDuplicated(names(column_map)) ||
      anyDuplicated(unname(column_map))) {
    stop(
      "Internal error: interaction model lacks a valid transformed-to-model column mapping.",
      call. = FALSE
    )
  }
  
  missing_model_cols <- setdiff(target_cols, unname(column_map))
  if (length(missing_model_cols) > 0L) {
    stop(
      paste0(
        "Cannot reconstruct ordinary MFPI prediction matrix for term `", term,
        "`: mapping is missing fitted-model column(s): ",
        paste(missing_model_cols, collapse = ", "), "."
      ),
      call. = FALSE
    )
  }
  
  source_cols <- names(column_map)[match(target_cols, unname(column_map))]
  missing_source_cols <- setdiff(source_cols, colnames(X))
  if (length(missing_source_cols) > 0L) {
    stop(
      paste0(
        "Cannot reconstruct ordinary MFPI prediction matrix for term `", term,
        "`: missing source design column(s): ",
        paste(missing_source_cols, collapse = ", "), "."
      ),
      call. = FALSE
    )
  }
  
  X_source <- X[, source_cols, drop = FALSE]
  X <- X_source
  colnames(X) <- target_cols
  
  if (anyNA(X) || any(!is.finite(X))) {
    stop(paste0("Non-finite values were produced in the prediction matrix for term `",
                term, "`."), call. = FALSE)
  }
  
  fit_obj <- interaction_model$fit
  has_offset <- mfpi_fit_has_offset(fit_obj)
  
  if (!has_offset && !is.null(newoffset)) {
    stop(
      "`newoffset` can be supplied only when the retained interaction model ",
      "was fitted with an offset.",
      call. = FALSE
    )
  }
  
  # When training rows are reconstructed solely to replace strata, preserve the
  # fitted offset rather than requiring the caller to repeat it. For supplied
  # newdata, formula offsets are reconstructed earlier; matrix-interface fits
  # still require an explicit replacement because the raw offset is unavailable.
  if (has_offset && is.null(newoffset) && isTRUE(newdata_is_training_internal)) {
    stored_offset <- fit_obj$offset
    if (!is.null(stored_offset) && inherits(fit_obj, "coxph")) {
      offset_origin <- object$cox_offset_reference
      if (!is.numeric(offset_origin) || length(offset_origin) != 1L ||
          is.na(offset_origin) || !is.finite(offset_origin)) {
        stop(
          "The fitted MFPI object lacks valid Cox offset-reference metadata. ",
          "Refit the model with the current package version.",
          call. = FALSE
        )
      }
      # coxph stores its fitted offset after subtracting the mean training
      # offset. New prediction data must contain the original offset scale,
      # because predict.coxph() performs that centering itself.
      stored_offset <- stored_offset + unname(offset_origin)
    }
    if (is.null(stored_offset)) {
      fitted_frame <- try(stats::model.frame(fit_obj), silent = TRUE)
      if (!inherits(fitted_frame, "try-error")) {
        stored_offset <- stats::model.offset(fitted_frame)
      }
    }
    if (is.null(stored_offset)) {
      stop(
        "The retained interaction model used an offset, but its fitted offset ",
        "could not be recovered. Supply `newoffset` explicitly.",
        call. = FALSE
      )
    }
    newoffset <- stored_offset
  }
  
  if (is.null(newoffset)) {
    if (has_offset) {
      stop(
        paste0(
          "The fitted interaction model used an offset. Supply `newoffset` with ",
          "one value per row of `newdata`, or include the original formula ",
          "offset variable(s) in `newdata`."
        ),
        call. = FALSE
      )
    }
    pred_offset <- rep(0, n)
  } else {
    if (!is.numeric(newoffset) || length(newoffset) != n || anyNA(newoffset) ||
        any(!is.finite(newoffset))) {
      stop("`newoffset` must be a finite numeric vector with one value per prediction row.",
           call. = FALSE)
    }
    pred_offset <- as.numeric(newoffset)
  }
  
  # Build formula-compatible newdata using raw source-column names. The
  # returned diagnostic matrix X uses fitted coefficient names, while native
  # model prediction must receive source-column names from the fit-time frame.
  model_newdata <- as.data.frame(X_source, check.names = FALSE)
  
  expects_offset <- has_offset
  expects_strata <- mfpi_fit_has_strata(fit_obj)
  
  if (expects_offset) {
    model_newdata$offset_ <- pred_offset
  }
  
  if (expects_strata) {
    # Reconstruction of training rows can reuse the fitted stratum membership.
    # This allows a caller to replace only the offset without redundantly
    # supplying the original strata.
    if (is.null(strata) && isTRUE(newdata_is_training_internal)) {
      strata <- fit_obj$strata
    }
    
    if (is.null(strata)) {
      stop(
        paste0(
          "The fitted Cox interaction model is stratified. Supply `strata` ",
          "with one value or row per prediction row, or include the original ",
          "formula strata variable(s) in `newdata`."
        ),
        call. = FALSE
      )
    }
    
    strata_n <- if (is.vector(strata) || is.factor(strata)) length(strata) else NROW(strata)
    if (strata_n != n) {
      stop("`strata` must have one value or row per prediction row.", call. = FALSE)
    }
    if (anyNA(strata)) {
      stop("`strata` must not contain missing values.", call. = FALSE)
    }
    strata_frame <- as.data.frame(strata, check.names = FALSE)
    bad_numeric_strata <- names(strata_frame)[vapply(
      strata_frame,
      function(column) is.numeric(column) && any(!is.finite(column)),
      logical(1L)
    )]
    if (length(bad_numeric_strata) > 0L) {
      stop("Numeric `strata` values must be finite.", call. = FALSE)
    }
    
    model_newdata$strata_ <- if (is.matrix(strata) || is.data.frame(strata)) {
      do.call(
        survival::strata,
        c(as.list(as.data.frame(strata)), list(shortlabel = TRUE))
      )
    } else {
      strata
    }
  }
  
  list(
    X = X,
    offset = pred_offset,
    model_newdata = model_newdata,
    has_offset = has_offset,
    group_internal = group_internal,
    group_levels = group_labels,
    newdata = !isTRUE(newdata_is_training_internal),
    reconstructed = TRUE
  )
}


#' Delegate Ordinary MFPI Prediction to the Stored Model
#'
#' Ordinary subject-level predictions are calculated exclusively by the stored
#' `glm` or `coxph` object. MFPI reconstructs formula-compatible `newdata`, but
#' it does not reproduce native prediction, centering, offset, baseline-hazard,
#' or standard-error calculations manually.
#'
#' @param object Object of class `"mfpi"`.
#' @param term Character scalar naming the continuous variable whose
#'   term-specific model is being used.
#' @param fit_result Stored term-specific flex fit result.
#' @param model_newdata Reconstructed formula-compatible prediction data, or
#'   `NULL` to predict on the stored fitting data.
#' @param type Validated family-native prediction type.
#' @param se.fit Logical scalar. Whether to request standard errors.
#' @param cox_reference Cox covariate reference for `"lp"` or `"risk"`, or
#'   `NULL` for other prediction types and GLMs.
#'
#' @return List with numeric `fit`, optional numeric `se.fit`, and
#'   `used_model_predict = TRUE`.
#'
#' @keywords internal
#' @noRd
mfpi_predict_ordinary <- function(object,
                                  term,
                                  fit_result,
                                  model_newdata = NULL,
                                  type,
                                  se.fit,
                                  cox_reference = NULL) {
  interaction_model <- fit_result$test_results$interaction_model
  fit_obj <- interaction_model$fit
  family_string <- mfpi_family_string(object)
  
  if (is.null(fit_obj)) {
    stop(
      "No stored fitted interaction model is available for term `", term, "`.",
      call. = FALSE
    )
  }
  
  predict_args <- list(
    object = fit_obj,
    type = type,
    se.fit = se.fit
  )
  
  if (!is.null(model_newdata)) {
    predict_args$newdata <- model_newdata
  }
  
  if (identical(family_string, "cox") && type %in% c("lp", "risk")) {
    predict_args$reference <- cox_reference
  }
  
  pred <- do.call(stats::predict, predict_args)
  
  if (is.list(pred)) {
    return(list(
      fit = as.numeric(pred$fit),
      se.fit = if (is.null(pred$se.fit)) NULL else as.numeric(pred$se.fit),
      used_model_predict = TRUE
    ))
  }
  
  list(
    fit = as.numeric(pred),
    se.fit = NULL,
    used_model_predict = TRUE
  )
}

# -----------------------------------------------------------------------------
# Validation and scalar helpers ------------------------------------------------
# -----------------------------------------------------------------------------

#' Extract a Finite Named Scalar from MFPI Metadata
#'
#' Retrieves a scalar metadata value, such as a shift or scale value, from a
#' named vector. If \code{x} is \code{NULL}, or if \code{x} is a longer vector
#' without a matching name, the supplied \code{default} is returned.
#'
#' @param x Numeric vector or \code{NULL}.
#' @param name Character scalar naming the value to retrieve.
#' @param default Numeric scalar used when no matching value is stored.
#'
#' @return Finite numeric scalar.
#'
#' @keywords internal
#' @noRd
mfpi_named_scalar <- function(x, name, default) {
  # Metadata such as shift and scale may be stored as named vectors, length-one
  # vectors, or absent values. Normalize those cases to one finite scalar.
  
  if (is.null(x)) return(default)
  if (!is.null(names(x)) && name %in% names(x)) {
    val <- x[[name]]
  } else if (length(x) == 1L) {
    val <- x[[1L]]
  } else {
    return(default)
  }
  
  if (!is.numeric(val) || length(val) != 1L || anyNA(val) || !is.finite(val)) {
    stop(paste0("Stored scalar for `", name, "` is not finite."),
         call. = FALSE)
  }
  val
}


#' Determine the Fitted Family String
#'
#' Extracts the normalized family name from an MFPI object. The helper accepts
#' either an explicitly stored \code{family_string}, a character \code{family},
#' or a family object/list with a \code{family} element.
#'
#' @param object Object of class \code{"mfpi"}.
#'
#' @return Character scalar family name, such as \code{"gaussian"},
#'   \code{"binomial"}, \code{"poisson"}, or \code{"cox"}.
#'
#' @keywords internal
#' @noRd
mfpi_family_string <- function(object) {
  # Prefer the normalized family string stored by mfpi(); fall back only for
  # older or manually constructed objects.
  
  if (!is.null(object$family_string)) {
    return(unname(as.character(object$family_string)[1L]))
  }
  fam <- object$family
  if (is.character(fam)) return(unname(as.character(fam)[1L]))
  if (is.list(fam) && !is.null(fam$family)) {
    return(unname(as.character(fam$family)[1L]))
  }
  stop("Cannot determine the model family from the MFPI object.", call. = FALSE)
}


#' Resolve an MFPI Prediction Type for the Fitted Family
#'
#' @param type Character scalar or `NULL`. `NULL` selects `"both"`.
#' @param family_string Character scalar identifying the fitted family.
#'
#' @return A validated, family-native prediction type.
#'
#' @keywords internal
#' @noRd
mfpi_match_prediction_type <- function(type, family_string) {
  if (is.null(type)) {
    return("both")
  }
  
  type <- normalize_prediction_type(type, family_string)
  
  choices <- if (identical(family_string, "cox")) {
    c("both", "function", "difference", "lp", "risk", "expected", "survival")
  } else {
    c("both", "function", "difference", "link", "response")
  }
  
  tryCatch(
    match.arg(type, choices),
    error = function(e) {
      alias_note <- if (identical(family_string, "cox")) {
        " The alias 'link' is also accepted for 'lp'."
      } else {
        " The alias 'lp' is also accepted for 'link'."
      }
      stop(
        "For ", if (identical(family_string, "cox")) "Cox" else "GLM",
        " MFPI models, `type` must be one of: ",
        paste(shQuote(choices), collapse = ", "),
        ".", alias_note,
        call. = FALSE
      )
    }
  )
}


#' Detect Stratification in a Stored Cox Model
#'
#' @param fit_obj Fitted model object.
#'
#' @return Logical scalar.
#'
#' @keywords internal
#' @noRd
mfpi_fit_has_strata <- function(fit_obj) {
  prediction_fit_has_strata(fit_obj)
}


#' Retrieve the Zero-Handling Flag for One Term
#'
#' Determines whether a continuous variable was fitted with structural-zero
#' handling. The helper checks prediction metadata stored directly on the MFPI
#' object first, then falls back to adjustment-model metadata when available.
#'
#' @param object Object of class \code{"mfpi"}.
#' @param term Character scalar naming the variable.
#'
#' @return Logical scalar. \code{TRUE} means non-positive values are treated as
#'   structural zeros for FP prediction.
#'
#' @keywords internal
#' @noRd
mfpi_get_zero_var <- function(object, term) {
  # Structural-zero flags may live in newer prediction metadata or in older
  # adjustment-model metadata. Check both locations for backward compatibility.
  
  z <- NULL
  if (!is.null(object$zero_vars)) z <- object$zero_vars
  if (is.null(z) && !is.null(object$adjustment_model$zero)) z <- object$adjustment_model$zero
  if (is.null(z) && !is.null(object$adjustment_model$fp_terms$zero)) {
    z <- object$adjustment_model$fp_terms$zero
    names(z) <- rownames(object$adjustment_model$fp_terms)
  }
  if (!is.null(z) && !is.null(names(z)) && term %in% names(z)) return(isTRUE(z[[term]]))
  FALSE
}


#' Validate Group-to-Coefficient Metadata
#'
#' Checks the \code{coefficient_groups} object required to reconstruct
#' group-specific FP blocks. The object must be a list with one non-empty
#' character vector per group and no duplicated coefficient names across groups.
#'
#' @param groups List mapping group labels to fitted coefficient names.
#' @param term Optional character scalar used only to make error messages more
#'   informative.
#'
#' @return Invisibly \code{TRUE}.
#'
#' @keywords internal
#' @noRd
mfpi_validate_coefficient_groups <- function(groups, term = NULL) {
  # A valid coefficient map is essential: every group must have at least one
  # fitted FP coefficient and no coefficient may belong to two groups.
  
  label <- if (is.null(term)) "`coefficient_groups`" else paste0("`coefficient_groups` for term `", term, "`")
  if (!is.list(groups) || length(groups) < 2L) {
    stop(paste0(label, " must be a list with at least two groups."),
         call. = FALSE)
  }
  lens <- vapply(groups, length, integer(1L))
  if (any(lens == 0L)) stop(paste0(label, " must have non-empty elements."), call. = FALSE)
  cols <- unlist(groups, use.names = FALSE)
  if (!is.character(cols) || anyNA(cols) || any(!nzchar(cols))) {
    stop(paste0(label, " must contain non-empty coefficient names."),
         call. = FALSE)
  }
  if (anyDuplicated(cols)) {
    stop(paste0(label, " must not contain duplicated coefficient names."),
         call. = FALSE)
  }
  invisible(TRUE)
}


#' Validate Required Coefficients and Covariance Entries
#'
#' Ensures that every coefficient needed for prediction exists, is finite, and
#' has corresponding finite rows and columns in the covariance matrix. This
#' catches rank-deficient or non-estimable models before invalid fitted values or
#' standard errors are produced.
#'
#' @param coef_vec Named numeric coefficient vector.
#' @param vcov_mat Named square numeric covariance matrix, usually
#'   \code{vcov(interaction_model$fit)}.
#' @param required_names Character vector of coefficient names required by the
#'   current prediction calculation.
#' @param context Character string describing the calculation, used in error
#'   messages.
#'
#' @return Invisibly \code{TRUE}.
#'
#' @keywords internal
#' @noRd
mfpi_validate_required_coefficients <- function(coef_vec, vcov_mat,
                                                required_names, context) {
  # Prediction is name-driven. Fail early if any required coefficient or
  # covariance entry is missing, non-finite, or non-estimable.
  
  required_names <- unique(required_names)
  if (!is.numeric(coef_vec) || is.null(names(coef_vec))) {
    stop("Model coefficients must be a named numeric vector.", call. = FALSE)
  }
  
  missing_coef <- setdiff(required_names, names(coef_vec))
  if (length(missing_coef) > 0L) {
    stop(paste0("Cannot perform ", context, ": missing coefficients: ",
                paste(missing_coef, collapse = ", "), "."), call. = FALSE)
  }
  
  bad_coef <- required_names[is.na(coef_vec[required_names]) |
                               !is.finite(coef_vec[required_names])]
  if (length(bad_coef) > 0L) {
    stop(
      paste0(
        "Cannot perform ", context, ": non-finite coefficients indicate ",
        "rank deficiency or non-estimable terms: ",
        paste(bad_coef, collapse = ", "), "."
      ),
      call. = FALSE
    )
  }
  
  if (!is.matrix(vcov_mat) || !is.numeric(vcov_mat) ||
      nrow(vcov_mat) != ncol(vcov_mat) || is.null(rownames(vcov_mat)) ||
      is.null(colnames(vcov_mat))) {
    stop("Model covariance matrix must be a named square numeric matrix.",
         call. = FALSE)
  }
  
  missing_cov <- union(setdiff(required_names, rownames(vcov_mat)),
                       setdiff(required_names, colnames(vcov_mat)))
  if (length(missing_cov) > 0L) {
    stop(paste0("Cannot perform ", context,
                ": covariance matrix is missing entries for: ",
                paste(missing_cov, collapse = ", "), "."), call. = FALSE)
  }
  
  v_sub <- vcov_mat[required_names, required_names, drop = FALSE]
  if (anyNA(v_sub) || any(!is.finite(v_sub))) {
    stop(
      paste0(
        "Cannot perform ", context,
        ": covariance matrix contains non-finite entries for required coefficients."
      ),
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}


#' Sanitize Prediction Variances
#'
#' Checks prediction variances before square-rooting them. Non-finite variances
#' and materially negative variances produce errors. Very small negative values
#' within numerical tolerance are truncated to zero.
#'
#' @param v Numeric vector of variances.
#' @param context Character string describing the calculation.
#' @param tol Numeric tolerance used to distinguish numerical roundoff from a
#'   genuinely negative variance.
#'
#' @return Numeric vector with small negative values replaced by zero.
#'
#' @keywords internal
#' @noRd
mfpi_sanitize_variance <- function(v, context,
                                   tol = sqrt(.Machine$double.eps)) {
  # Variances should be non-negative. Tiny negative values can arise from
  # floating-point roundoff; materially negative values indicate a real problem.
  
  if (anyNA(v) || any(!is.finite(v))) {
    stop(paste0("Non-finite variances produced during ", context, "."),
         call. = FALSE)
  }
  if (any(v < -tol)) {
    stop(
      paste0(
        "Negative variances produced during ", context,
        ". This suggests a non-positive-semidefinite covariance matrix or ",
        "invalid reconstruction."
      ),
      call. = FALSE
    )
  }
  pmax(v, 0)
}


#' Resolve the Reference Group for Fitted-Function Differences
#'
#' Converts a user-supplied reference group into the internal group position
#' used by the prediction code. The user may supply either an internal group
#' label or an original group label when original-to-internal group-level
#' metadata are stored in the MFPI object.
#'
#' @param reference Optional reference group label. If \code{NULL}, the first
#'   internal group level is used.
#' @param group_labels Character vector of internal group labels in fitted-model
#'   order.
#' @param object Object of class \code{"mfpi"} containing optional
#'   \code{group_levels_original} and \code{group_levels_new} metadata.
#'
#' @return Integer position of the reference group within \code{group_labels}.
#'
#' @keywords internal
#' @noRd
mfpi_resolve_reference_group <- function(reference, group_labels, object) {
  # Users may specify the reference using internal codes or original labels;
  # prediction calculations need the corresponding internal group position.
  
  if (is.null(reference)) return(1L)
  
  ref_chr <- as.character(reference)
  hit <- match(ref_chr, group_labels)
  if (!is.na(hit)) return(hit)
  
  orig <- object$group_levels_original
  new <- object$group_levels_new
  if (!is.null(orig) && !is.null(new) && length(orig) == length(new)) {
    hit_orig <- match(ref_chr, as.character(orig))
    if (!is.na(hit_orig)) {
      hit_new <- match(as.character(new[hit_orig]), group_labels)
      if (!is.na(hit_new)) return(hit_new)
    }
  }
  
  stop(
    paste0(
      "`reference` does not match an internal or original group level. ",
      "Available internal levels: ", paste(group_labels, collapse = ", "), "."
    ),
    call. = FALSE
  )
}


#' Detect Whether a Fitted Model Used an Offset
#'
#' Inspects a fitted GLM or Cox model object for stored offset information. The
#' result is used to decide whether ordinary new-data prediction must require a
#' new offset vector.
#'
#' @param fit_obj Fitted model object.
#'
#' @return Logical scalar.
#'
#' @keywords internal
#' @noRd
mfpi_fit_has_offset <- function(fit_obj) {
  # Detect offsets both from fitted-object storage and from the model terms.
  # This determines whether reconstructed prediction data require an offset.
  
  if (!is.null(fit_obj$offset)) return(TRUE)
  trm <- try(stats::terms(fit_obj), silent = TRUE)
  if (!inherits(trm, "try-error")) {
    off <- attr(trm, "offset")
    return(length(off) > 0L)
  }
  FALSE
}


#' Prepare New Data for Formula-Interface MFPI Prediction
#'
#' Replays the fit-time formula design recipe for formula-interface
#' \code{mfpi()} objects. The helper rebuilds formula-derived columns using the
#' stored terms object, factor levels, and contrasts, while preserving the raw
#' user-supplied columns needed by MFPI-specific prediction logic.
#'
#' This function does not rerun variable selection, refit models, or recompute
#' fractional-polynomial powers. It only rebuilds the formula-side design
#' columns that \code{model.matrix()} created at fit time.
#'
#' @param object Object of class \code{"mfpi"} containing formula metadata.
#' @param newdata Raw data frame supplied to \code{predict.mfpi()}.
#' @param terms Optional conceptual terms required by the selected interaction
#'   model. When supplied, unrelated original formula expressions are not
#'   evaluated.
#' @param required_columns Optional raw or model-matrix columns that must be
#'   present after reconstruction.
#'
#' @return A data frame containing the original raw columns plus any formula-
#'   derived design columns that were not already present.
#'
#' @keywords internal
#' @noRd
mfpi_prepare_formula_newdata <- function(object, newdata, terms = NULL,
                                         required_columns = NULL) {
  if (!is.data.frame(newdata)) {
    return(newdata)
  }
  
  terms_obj <- object$formula_terms
  if (is.null(terms_obj)) {
    stop(
      "Formula-interface prediction metadata are missing: `formula_terms` is NULL.",
      call. = FALSE
    )
  }
  
  xlevels <- object$formula_xlevels
  if (is.null(xlevels)) {
    xlevels <- list()
  }
  
  xcontrasts <- object$formula_contrasts
  
  if (is.null(terms)) {
    active_terms <- terms_obj
  } else {
    terms <- unique(as.character(terms))
    terms <- terms[!is.na(terms) & nzchar(terms)]
    if (length(terms) == 0L) {
      return(newdata)
    }
    active_terms <- prediction_formula_terms(object, terms)
  }
  
  active_factors <- attr(active_terms, "factors")
  active_frame_variables <- if (is.null(active_factors)) {
    character(0L)
  } else {
    rownames(active_factors)[rowSums(active_factors != 0) > 0L]
  }
  active_xlevels <- xlevels[intersect(names(xlevels), active_frame_variables)]
  
  expected_cols <- required_columns
  if (is.null(expected_cols)) {
    if (is.null(terms)) {
      expected_cols <- object$formula_design_columns
    } else {
      lookup <- object$formula_term_to_columns
      if (is.null(lookup)) lookup <- object$adjustment_term_to_columns
      if (!is.null(lookup)) {
        expected_cols <- unique(unlist(
          lookup[intersect(terms, names(lookup))],
          use.names = FALSE
        ))
      }
    }
  }
  expected_cols <- unique(as.character(expected_cols))
  expected_cols <- expected_cols[!is.na(expected_cols) & nzchar(expected_cols)]
  if (length(expected_cols) > 0L &&
      all(expected_cols %in% names(newdata))) {
    return(newdata)
  }
  
  mf <- tryCatch(
    stats::model.frame(
      active_terms,
      data = newdata,
      na.action = stats::na.pass,
      xlev = active_xlevels,
      drop.unused.levels = FALSE
    ),
    error = function(e) {
      stop(
        paste0(
          "Could not reconstruct the formula model frame for `newdata`: ",
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )
  
  active_contrasts <- xcontrasts[
    intersect(names(xcontrasts), names(mf))
  ]
  
  mm <- tryCatch(
    {
      if (length(active_contrasts) == 0L) {
        stats::model.matrix(active_terms, data = mf)
      } else {
        stats::model.matrix(
          active_terms,
          data = mf,
          contrasts.arg = active_contrasts
        )
      }
    },
    error = function(e) {
      stop(
        paste0(
          "Could not reconstruct the formula model matrix for `newdata`: ",
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )
  
  if ("(Intercept)" %in% colnames(mm)) {
    mm <- mm[, setdiff(colnames(mm), "(Intercept)"), drop = FALSE]
  }
  
  out <- newdata
  
  if (ncol(mm) > 0L) {
    mm_df <- as.data.frame(mm, optional = TRUE)
    names(mm_df) <- colnames(mm)
    
    add_cols <- setdiff(names(mm_df), names(out))
    if (length(add_cols) > 0L) {
      out <- cbind(out, mm_df[, add_cols, drop = FALSE])
    }
  }
  
  if (length(expected_cols) > 0L) {
    missing_cols <- setdiff(expected_cols, names(out))
    if (length(missing_cols) > 0L) {
      stop(
        paste0(
          "`newdata` could not be expanded to the selected interaction-model design. ",
          "Missing required column(s): ",
          paste(missing_cols, collapse = ", "),
          "."
        ),
        call. = FALSE
      )
    }
  }
  
  out
}
