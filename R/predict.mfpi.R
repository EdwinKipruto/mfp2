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
#' Computes predictions from an object of class \code{"mfpi"}.The main prediction
#'  target is the group-specific MFPI fitted function \eqn{f_j(x)}. For each 
#'  requested continuous variable, the method can return
#' fitted functions, pointwise standard errors, confidence intervals, and
#' contrasts between fitted functions. Ordinary link- and response-scale
#' prediction from the term-specific interaction model is also supported.
#'
#' @section Prediction targets:
#' \describe{
#'   \item{\code{type = "function"}}{Returns the group-specific fitted
#'   functions \eqn{f_j(x)}. The \code{functions} component contains one row
#'   for every evaluation \eqn{x} value and every group. When
#'   \code{se.fit = TRUE}, pointwise standard errors and confidence intervals
#'   are included.}
#'   \item{\code{type = "difference"}}{Returns fitted-function contrasts
#'   \eqn{f_j(x) - f_r(x)} relative to a reference group \eqn{r}. The
#'   \code{differences} component contains one row for every evaluation
#'   \eqn{x} value and every non-reference contrast.}
#'   \item{\code{type = "both"}}{Returns both group-specific fitted functions
#'   and fitted-function differences. This is the default fitted-function
#'   prediction target.}
#'   \item{\code{type = "link"}}{Returns ordinary subject-level linear
#'   predictor values from the term-specific interaction model. This target
#'   returns one prediction per row of the fitting data or \code{newdata}.}
#'   \item{\code{type = "response"}}{Returns ordinary subject-level
#'   response-scale predictions. For Gaussian models this is the fitted mean,
#'   for binomial models the fitted probability, for Poisson models the fitted
#'   mean count or rate, and for Cox models the relative risk score
#'   \eqn{\exp(\eta)}. Cox response prediction is not an absolute survival
#'   probability.}
#' }
#'
#' @section Term-specific prediction:
#' MFPI fits one interaction model per tested continuous variable. Predictions
#' are therefore term-specific. If one term is requested, the method returns one
#' object of class \code{"mfpi_prediction"}. If multiple terms are requested,
#' it returns a named list of such objects with class
#' \code{"mfpi_prediction_list"}.
#'
#' @section Choosing terms and model scope:
#' Each tested continuous variable has a single fitted interaction model,
#' compared against its main-effects-only counterpart using the model's
#' selection \code{criterion} (AIC, BIC, or p-value) at fit time.
#' \code{model = "best"} restricts prediction to variables whose interaction
#' was retained by that criterion in the final selected model. \code{model =
#' "all"} widens this to every tested variable's fitted interaction model,
#' including variables that were not retained, as long as a fitted model was
#' stored for that variable.
#'
#' \code{terms} and \code{model} work together: \code{terms} picks which
#' variables to predict, while \code{model} determines which pool of fitted
#' models those variables are drawn from. If \code{terms = NULL}, all
#' variables available within the requested \code{model} scope are used. If
#' \code{terms} is supplied, every requested variable must have a fitted model
#' within that scope, or the call fails with the names of the missing
#' variable(s); this is most likely to happen when a variable's interaction
#' was not retained and \code{model = "best"} is used.
#'
#' @section Fitted-function prediction:
#' Fitted-function prediction reconstructs only the group-specific FP basis
#' needed for \eqn{f_j(x)} and combines that basis with the fitted
#' interaction-model coefficients and covariance matrix. This target is for
#' plotting and interpretation, not subject-level risk or mean prediction.
#'
#' Each group-specific function is evaluated at every supplied \eqn{x} value.
#' Therefore, if there are \eqn{n} evaluation values and \eqn{k} groups, the
#' \code{functions} component contains \eqn{n \times k} rows. With
#' \code{grid = FALSE}, evaluation values are the observed or supplied values
#' in their original row order; duplicates are retained. With
#' \code{grid = TRUE}, values are replaced by an equally spaced grid. For
#' zero-handled variables, the grid is built over the positive part and includes
#' a structural zero point if structural zeros are present.
#'
#' Fitted-function differences are computed relative to a reference group.
#' Unless \code{reference} is supplied, the reference group used at model fitting
#' is reused. Supplying \code{reference} only changes the contrast baseline for
#' fitted-function differences returned by \code{type = "difference"} or
#' \code{type = "both"}; it does not refit the model or change the fit-time
#' group coding used by ordinary prediction.
#'
#' @section Ordinary link and response prediction:
#' Ordinary \code{type = "link"} and \code{type = "response"} prediction are
#' subject-level targets. For \code{newdata = NULL} and \code{newoffset = NULL},
#' the stored fitted model's own prediction method is used. If
#' \code{newdata = NULL} but a replacement \code{newoffset} is supplied, the
#' stored training design is reconstructed so that the supplied new offset is used
#' explicitly. For supplied \code{newdata}, the method reconstructs the full
#' term-specific interaction-model matrix, including group dummies, the
#' row-specific group interaction block, selected adjustment covariates, and
#' offsets when applicable. Formula-level offsets are reconstructed from
#' \code{newdata} when possible; an explicitly supplied \code{newoffset}
#' takes precedence. The reconstructed matrix is then multiplied by the fitted
#' coefficient vector.
#'
#' Ordinary prediction uses the group coding, reference level, contrasts, and
#' centering conventions stored at model fitting. New observations must belong
#' to group levels observed during fitting. New group levels are rejected
#' because no fitted coefficients exist for groups outside the fit-time design.
#'
#' @section Categorical adjustment variables in prediction:
#' For formula fits, \code{newdata} should contain the original factor columns,
#' not manually created dummy columns. The stored predictor terms, factor
#' levels, and contrasts are used to recreate the fit-time design. All contrast
#' columns belonging to a selected categorical adjustment term are then
#' transformed and inserted together. Unseen factor levels are rejected.
#'
#' For default-interface fits using \code{term_groups}, \code{newdata} must
#' contain every raw member column of each selected grouped term. Supplying only
#' part of a block is an error. Categorical adjustment columns are reused as
#' fixed linear columns; no FP, ACD, zero, catzero, or spike transformation is
#' applied at prediction time.
#'
#' @section Cox models:
#' For Cox models, \code{type = "link"} returns the zero-reference linear
#' predictor and \code{type = "response"} returns the relative risk score
#' \eqn{\exp(\eta)}. Baseline survival, cumulative-hazard, expected-event, and
#' absolute-risk prediction are not currently implemented. Formula-level Cox
#' strata used during fitting affect the partial likelihood and
#' stratum-specific baseline hazards. Prediction-time strata for new
#' observations are supplied via the \code{strata} argument; see below.
#'
#' @section Standard errors and confidence intervals:
#' \code{se.fit} and \code{level} apply to both prediction targets, but the
#' quantity they describe differs between them. For fitted functions and
#' fitted-function differences, \code{se.fit = TRUE} adds pointwise standard
#' errors on the model scale, together with \code{lower}/\code{upper}
#' confidence limits at the requested \code{level}; these describe the
#' precision of \eqn{f_j(x)} or \eqn{f_j(x) - f_r(x)} at each evaluation
#' point, not a simultaneous confidence band. For ordinary \code{"link"} and
#' \code{"response"} prediction, \code{se.fit = TRUE} adds a single
#' \code{se.fit} column to \code{predictions}; \code{level} is not used here,
#' since no interval is returned for ordinary prediction. Response-scale
#' standard errors follow the same inverse-link delta-method convention as
#' \code{predict.glm()}, and Cox response-scale standard errors follow the
#' zero-reference relative-risk convention described under \emph{Cox models}.
#'
#' @section Required object metadata:
#' New-data prediction assumes that the \code{"mfpi"} object stores prediction
#' metadata created at fit time. In particular, fitted-function prediction with
#' \code{newdata = NULL} and ordinary prediction with \code{newdata = NULL}
#' plus a replacement \code{newoffset} require \code{object$x_train_internal},
#' the post-preprocessing training matrix containing shifted/scaled predictors
#' and internally remapped \code{group_var} codes. New-data ordinary prediction
#' additionally relies on stored shift, scale, winsorisation limits, group-level
#' mapping, the conceptual-term-to-column lookup, selected adjustment-model
#' powers, centering constants, formula factor levels and contrasts when
#' applicable, and each flex result's \code{coefficient_groups}.
#'
#' @param object An object of class \code{"mfpi"}.
#' @param newdata Optional matrix or data frame. For fitted-function prediction
#'   (\code{type = "function"}, \code{"difference"}, or \code{"both"}), it
#'   must contain the requested continuous-variable columns on their original
#'   raw scale. The grouping variable is not required for fitted-function
#'   prediction because all group-specific functions are evaluated at each
#'   supplied \eqn{x} value. For ordinary prediction (\code{type = "link"} or
#'   \code{"response"}), \code{newdata} must also contain the grouping variable
#'   and all selected adjustment variables required by the term-specific
#'   interaction model. Other predictors from the original formula may be
#'   omitted. For formula fits, supply original factor-valued columns; for
#'   default-interface grouped terms, supply every raw member column. Group
#'   values and factor values used by selected terms must use levels observed
#'   during fitting.
#' @param terms Character vector of continuous variables to predict. If
#'   \code{NULL}, terms are selected from the requested \code{model} scope. See
#'   \emph{Choosing terms and model scope}.
#' @param model Character scalar, either \code{"best"} (retained interaction
#'   models only) or \code{"all"} (every tested variable's fitted interaction
#'   model). See \emph{Choosing terms and model scope}.
#' @param type Character scalar. One of \code{"both"}, \code{"function"},
#'   \code{"difference"}, \code{"link"}, or \code{"response"}.
#' @param se.fit Logical scalar. Whether to compute pointwise standard errors.
#'   For fitted-function targets, these are standard errors of fitted functions
#'   or function differences on the model scale. For ordinary GLM prediction,
#'   response-scale standard errors use the same inverse-link derivative
#'   convention as \code{predict.glm()}. For Cox response prediction,
#'   \code{fit} is the relative risk score \eqn{\exp(\eta)} and
#'   \code{se.fit} follows the zero-reference relative-risk convention.
#' @param level Numeric scalar in \code{(0, 1)}. Confidence level for
#'   pointwise fitted-function and fitted-difference intervals.
#' @param grid Logical scalar. If \code{TRUE}, evaluate fitted functions on a
#'   grid. Ignored for ordinary \code{"link"} and \code{"response"}
#'   prediction.
#' @param n_grid Integer scalar. Number of positive-part grid points when
#'   \code{grid = TRUE}. If a structural zero is added, the returned number of
#'   \code{x} values can be \code{n_grid + 1}.
#' @param reference Optional reference group for fitted-function differences.
#'   May be an internal group label or an original group label if group-level
#'   metadata are stored in \code{object}. If \code{NULL}, the reference group
#'   used at model fitting is reused. This argument affects only fitted-function
#'   differences; it does not change the fit-time reference level used by the
#'   stored interaction model.
#' @param strata Optional prediction-time Cox strata for ordinary \code{"link"}
#'   or \code{"response"} prediction. Ignored for fitted-function prediction.
#'   Supply one value per row, or a matrix/data frame with one row per prediction
#'   row for multiple strata variables.
#' @param newoffset Optional numeric vector for ordinary \code{"link"} or
#'   \code{"response"} prediction. If \code{newdata} is supplied,
#'   \code{newoffset} must contain one value per row of \code{newdata}. If
#'   \code{newdata = NULL}, a supplied \code{newoffset} must contain one value
#'   per training row; the stored training design is reconstructed and the
#'   supplied new offset is used instead of delegating to the fitted model's own
#'   \code{predict()} method. For formula-interface fits with a formula-level
#'   offset, \code{newoffset} may be omitted when \code{newdata} contains the
#'   original offset variable(s), because the offset is reconstructed from the
#'   stored formula metadata. Ignored for fitted-function prediction.
#' @param ... Reserved for future extensions. Currently unused arguments produce
#'   a warning. Cox baseline-survival arguments such as prediction times are
#'   not currently supported.
#'
#' @return If one term is requested, an object of class
#' \code{"mfpi_prediction"}. If multiple terms are requested, a named list of
#' \code{"mfpi_prediction"} objects with class \code{"mfpi_prediction_list"}.
#'
#' For \code{type = "function"}, \code{"difference"}, or \code{"both"}, each
#' prediction object contains:
#' \describe{
#'   \item{\code{term}}{Continuous variable for which prediction was computed.}
#'   \item{\code{type}}{Requested prediction type.}
#'   \item{\code{x}}{Evaluation values on the original raw scale. With
#'   \code{grid = FALSE}, these are observed or supplied values in row order,
#'   including duplicates. With \code{grid = TRUE}, these are grid values.}
#'   \item{\code{functions}}{Data frame of fitted functions, present for
#'   \code{type = "function"} and \code{type = "both"}. Columns are
#'   \code{term}, \code{x}, \code{group}, \code{fit}, and, when
#'   \code{se.fit = TRUE}, \code{se.fit}, \code{lower}, and \code{upper}.}
#'   \item{\code{differences}}{Data frame of fitted-function differences,
#'   present for \code{type = "difference"} and \code{type = "both"}. Columns
#'   are \code{term}, \code{x}, \code{contrast}, \code{group},
#'   \code{reference}, \code{fit}, and, when \code{se.fit = TRUE},
#'   \code{se.fit}, \code{lower}, and \code{upper}.}
#'   \item{\code{metadata}}{List containing group-specific FP powers,
#'   centering constants, coefficient groups, reference group, scale and shift
#'   values, zero-handling flag, confidence level, and \code{se.fit}.}
#' }
#'
#' For \code{type = "link"} or \code{"response"}, each prediction object
#' contains:
#' \describe{
#'   \item{\code{term}}{Continuous variable whose term-specific interaction
#'   model was used.}
#'   \item{\code{type}}{Either \code{"link"} or \code{"response"}.}
#'   \item{\code{predictions}}{Data frame with one row per subject. It
#'   contains \code{term}, \code{fit}, and, when requested and available,
#'   \code{se.fit}.}
#'   \item{\code{design}}{Reconstructed numeric interaction-model matrix.
#'   Present when ordinary prediction uses the manual reconstructed-design path,
#'   either because \code{newdata} was supplied or because a replacement
#'   \code{newoffset} was supplied for training-data prediction. \code{NULL}
#'   when prediction delegates to the stored fitted model's own \code{predict()}
#'   method.}
#'   \item{\code{metadata}}{List describing whether ordinary prediction was
#'   used, whether \code{newdata} was supplied, whether the design matrix was
#'   reconstructed, group metadata, offset usage, and response-scale status.}
#' }
#'
#' @details Fitted-function prediction reuses fit-time centering constants and
#' never recomputes centers from \code{newdata}. Recomputing centers during
#' prediction would define a different FP basis and make the reconstructed
#' fitted functions incompatible with the fitted coefficients.
#'
#' This method does not call or depend on the legacy fitted-function generator.
#' Its shared dependency is \code{build_group_fp_basis()}, which must remain
#' available after \code{gen_fitted_values_per_group()} is removed.
#'
#' @examples
#' \dontrun{
#' # Fit an MFPI model. The exact arguments depend on your analysis.
#' fit <- mfpi(x, y, group_var = "trt", cont_vars = c("age", "bmi"))
#'
#' # Default: selected term(s), fitted functions and differences, observed x.
#' p1 <- predict(fit)
#' names(p1)
#' head(p1$functions)
#' head(p1$differences)
#'
#' # Smooth fitted-function curves for plotting.
#' p2 <- predict(fit, terms = "age", type = "both", grid = TRUE)
#' age_fun <- subset(p2$functions, group == unique(p2$functions$group)[1L])
#' plot(age_fun$x, age_fun$fit, type = "l")
#'
#' # Fitted functions only.
#' p3 <- predict(fit, terms = "age", type = "function")
#'
#' # Formula-fitted categorical adjustment variables are supplied on their
#' # original factor scale; predict.mfpi() recreates the stored contrasts.
#' nd <- training_data[1:20, c("trt", "age", "stage"), drop = FALSE]
#' p_factor <- predict(fit_with_stage, newdata = nd, terms = "age",
#'                     model = "all", type = "link")
#'
#' # Differences relative to the fit-time reference group.
#' p4 <- predict(fit, terms = "age", type = "difference")
#'
#' # Differences relative to an explicitly chosen observed group.
#' p5 <- predict(fit, terms = "age", type = "difference", reference = 0)
#'
#' # Ordinary subject-level response prediction for new observations.
#' new_x <- x[1:10, , drop = FALSE]
#' p6 <- predict(fit, terms = "age", type = "response", newdata = new_x)
#' head(p6$predictions)
#'
#' # Multiple terms return a named prediction list.
#' p7 <- predict(fit, terms = c("age", "bmi"), model = "all", grid = TRUE)
#' names(p7)
#' }
#'
#' @method predict mfpi
#' @export
predict.mfpi <- function(object,
                         newdata = NULL,
                         terms = NULL,
                         model = c("best", "all"),
                         type = c("both", "function", "difference", "link", "response"),
                         se.fit = TRUE,
                         level = 0.95,
                         grid = FALSE,
                         n_grid = 200L,
                         reference = NULL,
                         strata = NULL,
                         newoffset = NULL,
                         ...) {
  # Public S3 entry point. After argument validation, prediction is routed to
  # either the fitted-function path or the ordinary subject-level path.
  
  if (!inherits(object, "mfpi")) {
    stop("`object` must be an object of class \"mfpi\".", call. = FALSE)
  }
  
  model <- match.arg(model)
  type <- match.arg(type)
  
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
  
  n_grid <- as.integer(n_grid)
  if (length(n_grid) != 1L || is.na(n_grid) || n_grid < 2L) {
    stop("`n_grid` must be an integer >= 2.", call. = FALSE)
  }
  
  dots <- list(...)
  if (length(dots) > 0L) {
    warning(
      "Unused arguments in `...`: ", paste(names(dots), collapse = ", "), ".",
      call. = FALSE
    )
  }
  
  if (type %in% c("link", "response") && isTRUE(grid)) {
    warning("`grid` is ignored for `type = \"link\"` and `type = \"response\"`.",
            call. = FALSE)
  }
  
  fits <- mfpi_get_prediction_fits(object = object, terms = terms, model = model)
  
  out <- stats::setNames(
    lapply(names(fits), function(term) {
      fit_result <- fits[[term]]
      
      if (type %in% c("function", "difference", "both")) {
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
          reference = reference
        ))
      }
      
      # Ordinary prediction is subject-level: each row belongs to one group,
      # so the reconstructed interaction block is masked row by row. Formula-level
      # Cox strata are reconstructed from raw newdata unless explicitly supplied.
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
      
      pred <- mfpi_predict_from_design(
        object = object,
        term = term,
        fit_result = fit_result,
        X_new = design$X,
        offset = design$offset,
        model_newdata = design$model_newdata,
        type = type,
        se.fit = se.fit,
        use_model_predict = design$use_model_predict
      )
      
      pred_df <- data.frame(term = term, fit = pred$fit, stringsAsFactors = FALSE)
      if (se.fit && !is.null(pred$se.fit)) pred_df$se.fit <- pred$se.fit
      
      structure(
        list(
          term = term,
          type = type,
          predictions = pred_df,
          design = if (!isTRUE(design$use_model_predict)) design$X else NULL,
          metadata = list(
            ordinary_prediction = TRUE,
            newdata = isTRUE(design$newdata),
            design_matrix_reconstructed = !isTRUE(design$use_model_predict),
            group_levels = design$group_levels,
            group_internal = design$group_internal,
            offset_used = design$has_offset,
            response_scale = type == "response",
            used_model_predict = isTRUE(pred$used_model_predict),
            model_newdata_columns = if (!is.null(design$model_newdata)) colnames(design$model_newdata) else NULL
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
# Ordinary link/response prediction -------------------------------------------
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
#' Reconstructs the numeric design matrix required for ordinary
#' \code{type = "link"} or \code{type = "response"} prediction for one MFPI
#' term. This helper is used only for subject-level prediction. It differs from
#' fitted-function prediction because each new subject belongs to exactly one
#' group, so out-of-group interaction blocks must be set to zero row by row.
#'
#' @section Training-data path:
#' If \code{newdata = NULL} and \code{newoffset = NULL}, no matrix is
#' reconstructed. The returned object tells \code{mfpi_predict_from_design()}
#' to call the fitted model's own \code{predict()} method. This preserves the
#' original GLM or Cox fitted-value behavior. If \code{newdata = NULL} but a
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
#'   fitted interaction model used an offset and manual prediction is used.
#'
#' @return A list with \code{X} (design matrix without intercept, or
#'   \code{NULL} when prediction delegates to the fitted model),
#'   \code{offset}, \code{has_offset}, \code{group_internal},
#'   \code{group_levels}, \code{newdata}, and \code{use_model_predict}.
#'
#' @keywords internal
#' @noRd
mfpi_build_ordinary_design <- function(object, term, fit_result, newdata,
                                       strata = NULL, newoffset = NULL) {
  # Build the subject-level model matrix only when delegation to the fitted
  # model is not possible or not appropriate, such as supplied newdata.
  
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
  # If a new `newoffset` is supplied with `newdata = NULL`, delegation would ignore
  # that new offset and use the offset stored at fit time. To honor the supplied
  # new offset, switch to the same manual reconstructed-design path used for
  # supplied `newdata`, but start from `x_train_internal`, the
  # post-preprocessing training matrix.
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
        use_model_predict = TRUE
      ))
    }
    
    newdata <- object$x_train_internal
    
    if (is.null(newdata)) {
      stop(
        paste0(
          "Cannot apply a supplied `newoffset` with `newdata = NULL` because ",
          "the MFPI object does not store `x_train_internal`, the ",
          "post-preprocessing training matrix. Refit the model with the current ",
          "version of `mfpi()` or supply explicit `newdata`."
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
  
  has_offset <- mfpi_fit_has_offset(interaction_model$fit)
  if (is.null(newoffset)) {
    if (has_offset) {
      stop(
        paste0(
          "The fitted interaction model used an offset. Supply `newoffset` with ",
          "one value per row of `newdata`."
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
    has_offset <- TRUE
  }
  
  # Build formula-compatible newdata using the raw source-column names. The
  # manual matrix X uses fitted coefficient names, while model prediction must
  # receive the names that appeared in the fit-time data frame.
  model_newdata <- as.data.frame(X_source, check.names = FALSE)
  
  fit_formula <- stats::formula(interaction_model$fit)
  formula_text <- paste(deparse(fit_formula), collapse = " ")
  expects_offset <- grepl("offset\\(offset_\\)", formula_text)
  expects_strata <- grepl("strata\\(strata_\\)", formula_text)
  
  if (expects_offset) {
    model_newdata$offset_ <- pred_offset
  }
  
  if (expects_strata) {
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
    use_model_predict = FALSE
  )
}


#' Predict Link or Response Values from a Design Matrix
#'
#' Computes ordinary subject-level predictions for one term-specific MFPI
#' interaction model. When \code{use_model_predict = TRUE}, the helper
#' delegates to the stored fitted model's \code{predict()} method. Otherwise,
#' it multiplies the reconstructed design matrix by the fitted coefficient
#' vector, adds the intercept and offset when applicable, and transforms to
#' the response scale if requested. GLM offsets are added directly. Cox offsets
#' are first expressed relative to the mean training offset, matching
#' `predict.coxph(reference = "zero")`.
#'
#' @section Standard errors:
#' For manual prediction, the helper first computes the standard error of the
#' linear predictor, \eqn{se_\eta}, as a row-wise quadratic form using the
#' fitted covariance matrix. For GLM response-scale prediction, it then matches
#' \code{predict.glm()} by multiplying \eqn{se_\eta} by
#' \code{abs(family$mu.eta(eta))}. For Cox response prediction,
#' \code{fit} is the relative risk score \eqn{\exp(\eta)} and
#' \code{se.fit} follows \code{predict.coxph(type = "risk",
#' reference = "zero")}: \eqn{se_\eta \sqrt{\exp(\eta)}}. This Cox
#' convention differs from the ordinary inverse-link delta-method value
#' \eqn{se_\eta \exp(\eta)}.
#'
#' @param object Object of class \code{"mfpi"}.
#' @param term Character scalar naming the continuous variable whose
#'   term-specific model is being used.
#' @param fit_result Stored term-specific flex fit result containing the fitted
#'   interaction model.
#' @param X_new Reconstructed design matrix without an intercept column, or
#'   \code{NULL} when \code{use_model_predict = TRUE}.
#' @param offset Numeric offset vector for manual prediction, or \code{NULL}
#'   when the fitted model's own prediction method is used. For manually
#'   reconstructed training-data prediction, this can be a replacement newoffset
#'   supplied with \code{newdata = NULL}.
#' @param type Character scalar, either \code{"link"} or \code{"response"}.
#' @param se.fit Logical scalar. Whether to compute standard errors.
#' @param use_model_predict Logical scalar. If \code{TRUE}, call the stored
#'   fitted model's \code{predict()} method instead of manual multiplication.
#'
#' @return List with \code{fit}, optional \code{se.fit}, and \code{link} when
#'   the linear predictor is computed manually.
#'
#' @keywords internal
#' @noRd
mfpi_predict_from_design <- function(object, term, fit_result, X_new, offset,
                                     model_newdata = NULL, type, se.fit,
                                     use_model_predict = FALSE) {
  # Final ordinary-prediction step. Either delegate to the fitted model's
  # predict() method or multiply the reconstructed design by the coefficients.
  
  interaction_model <- fit_result$test_results$interaction_model
  fit_obj <- interaction_model$fit
  family_string <- mfpi_family_string(object)
  
  if (isTRUE(use_model_predict) || !is.null(model_newdata)) {
    # Training-data prediction without a replacement newoffset delegates to the
    # fitted model. Cox uses reference = "zero" so delegated predictions
    # agree with the manual covariate scale; both paths also use the same
    # mean-training-offset origin.
    if (identical(family_string, "cox")) {
      pred_type <- if (type == "link") "lp" else "risk"
      
      predict_args <- list(
        object = fit_obj,
        type = pred_type,
        se.fit = se.fit,
        reference = "zero"
      )
      if (!is.null(model_newdata)) predict_args$newdata <- model_newdata
      pred <- do.call(stats::predict, predict_args)
      
      if (is.list(pred)) {
        fit <- as.numeric(pred$fit)
        
        return(list(
          fit    = fit,
          se.fit = if (!is.null(pred$se.fit)) as.numeric(pred$se.fit) else NULL,
          link   = if (type == "link") fit else NULL,
          used_model_predict = TRUE
        ))
      }
      
      fit <- as.numeric(pred)
      
      return(list(
        fit    = fit,
        se.fit = NULL,
        link   = if (type == "link") fit else NULL,
        used_model_predict = TRUE
      ))
    }
    
    predict_args <- list(object = fit_obj, type = type, se.fit = se.fit)
    if (!is.null(model_newdata)) predict_args$newdata <- model_newdata
    pred <- do.call(stats::predict, predict_args)
    if (is.list(pred)) {
      return(list(
        fit = as.numeric(pred$fit),
        se.fit = if (!is.null(pred$se.fit)) as.numeric(pred$se.fit) else NULL,
        link = NULL,
        used_model_predict = TRUE
      ))
    }
    return(list(fit = as.numeric(pred), se.fit = NULL, link = NULL,
                used_model_predict = TRUE))
  }
  
  coef_vec <- interaction_model$coefficients
  coef_vec <- stats::setNames(as.numeric(coef_vec), names(coef_vec))
  vcov_mat <- stats::vcov(fit_obj)
  
  has_intercept <- !identical(family_string, "cox") &&
    "(Intercept)" %in% names(coef_vec)
  beta_names <- setdiff(names(coef_vec), "(Intercept)")
  
  mfpi_validate_required_coefficients(
    coef_vec = coef_vec,
    vcov_mat = vcov_mat,
    required_names = unique(c(if (has_intercept) "(Intercept)", beta_names)),
    context = paste0("ordinary MFPI prediction for term `", term, "`")
  )
  
  if (!identical(colnames(X_new), beta_names)) {
    missing <- setdiff(beta_names, colnames(X_new))
    if (length(missing) > 0L) {
      stop(paste0("Prediction matrix is missing coefficient column(s): ",
                  paste(missing, collapse = ", "), "."), call. = FALSE)
    }
    X_new <- X_new[, beta_names, drop = FALSE]
  }
  
  # Manual path: construct the linear predictor from the reconstructed
  # design matrix. Cox offsets follow survival::predict.coxph(): the supplied
  # prediction offset is expressed relative to the mean offset in the fitting
  # data, independently of reference = "zero" for the covariates. GLM offsets
  # retain their ordinary uncentred interpretation.
  prediction_offset <- offset
  if (identical(family_string, "cox")) {
    prediction_offset <- offset - mfpi_cox_offset_reference(object)
  }
  
  eta <- as.vector(X_new %*% coef_vec[beta_names]) + prediction_offset
  if (has_intercept) eta <- eta + unname(coef_vec["(Intercept)"])
  
  fam <- object$family
  if (is.function(fam)) fam <- fam()
  
  # Convert the linear predictor to the requested scale. For Cox,
  # response-scale prediction is the relative risk score exp(eta).
  if (type == "link") {
    fit <- eta
  } else if (identical(family_string, "cox")) {
    fit <- exp(eta)
  } else if (is.list(fam) && is.function(fam$linkinv)) {
    fit <- fam$linkinv(eta)
  } else {
    fit <- switch(
      family_string,
      gaussian = eta,
      binomial = stats::binomial()$linkinv(eta),
      poisson = stats::poisson()$linkinv(eta),
      stop("Unsupported family for response-scale prediction.", call. = FALSE)
    )
  }
  
  se_out <- NULL
  if (se.fit) {
    # First compute SE(eta). Response-scale SEs are transformations of this
    # quantity; the Cox risk transformation intentionally follows coxph.
    X_vcov <- X_new
    if (has_intercept) X_vcov <- cbind("(Intercept)" = 1, X_vcov)
    v_sub <- vcov_mat[colnames(X_vcov), colnames(X_vcov), drop = FALSE]
    var_eta <- rowSums((X_vcov %*% v_sub) * X_vcov)
    se_eta <- sqrt(mfpi_sanitize_variance(
      var_eta,
      context = paste0("ordinary MFPI prediction for term `", term, "`")
    ))
    
    if (type == "link") {
      se_out <- se_eta
    } else if (identical(family_string, "cox")) {
      # Match survival::predict.coxph(type = "risk", se.fit = TRUE,
      # reference = "zero"). Cox response prediction returns exp(eta), but
      # coxph reports the risk-score SE as se_eta * sqrt(exp(eta)), not as
      # the ordinary inverse-link delta-method value se_eta * exp(eta).
      se_out <- sqrt(fit) * se_eta
    } else if (is.list(fam) && is.function(fam$mu.eta)) {
      se_out <- abs(fam$mu.eta(eta)) * se_eta
    } else {
      deriv <- switch(
        family_string,
        gaussian = rep(1, length(eta)),
        binomial = stats::binomial()$mu.eta(eta),
        poisson = stats::poisson()$mu.eta(eta),
        stop("Unsupported family for response-scale standard errors.",
             call. = FALSE)
      )
      se_out <- abs(deriv) * se_eta
    }
  }
  
  # Match the delegated predict.glm()/predict.coxph() branches: ordinary
  # prediction outputs are plain unnamed numeric vectors. Matrix row names can
  # otherwise propagate through rowSums() into se_eta and se_out.
  if (!is.null(se_out)) se_out <- as.numeric(se_out)
  
  list(
    fit = as.numeric(fit),
    se.fit = se_out,
    link = as.numeric(eta),
    used_model_predict = FALSE
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
  
  if (!is.null(object$family_string)) return(as.character(object$family_string)[1L])
  fam <- object$family
  if (is.character(fam)) return(fam[1L])
  if (is.list(fam) && !is.null(fam$family)) return(as.character(fam$family)[1L])
  stop("Cannot determine the model family from the MFPI object.", call. = FALSE)
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


#' Retrieve the Cox Offset Reference Stored at Fit Time
#'
#' `survival::predict.coxph()` subtracts the mean training offset from every
#' prediction offset, including when `reference = "zero"`. All Cox models
#' The final adjustment fit calculates that scalar once and the top-level
#' MFPI object stores it as `cox_offset_reference` for every interaction model.
#'
#' @param object A fitted `mfpi` object.
#'
#' @return Finite numeric scalar.
#'
#' @keywords internal
#' @noRd
mfpi_cox_offset_reference <- function(object) {
  reference <- object$cox_offset_reference
  
  if (!is.numeric(reference) || length(reference) != 1L ||
      is.na(reference) || !is.finite(reference)) {
    stop(
      paste0(
        "The fitted MFPI object lacks valid Cox offset-reference metadata. ",
        "Refit the MFPI object with the current package version."
      ),
      call. = FALSE
    )
  }
  
  unname(reference)
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
  # This determines whether manual newdata prediction must require an offset.
  
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