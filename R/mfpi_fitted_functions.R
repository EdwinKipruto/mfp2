# -----------------------------------------------------------------------------
# build_group_fp_basis() -------------------------------------------------------
# -----------------------------------------------------------------------------

#' Build Group-Specific FP Basis for MFPI Fitted-Function Prediction
#'
#' Builds the group-specific fractional-polynomial basis used to evaluate
#' fitted MFPI functions \eqn{f_j(x)}. The helper duplicates a one-column
#' continuous variable once per group, applies the selected group-specific FP
#' powers, and finally renames the transformed columns to the exact coefficient
#' names used by the fitted interaction model.
#'
#' This helper is intentionally separate from
#' \code{gen_fitted_values_per_group()} because the same reconstruction is needed
#' by \code{predict.mfpi()}.
#'
#' @param cont_mat Numeric matrix. If \code{transform_basis = TRUE}, this must be
#'   a one-column matrix containing the continuous variable on the backscaled
#'   scale used for FP transformation, usually \eqn{x + shift}. If
#'   \code{transform_basis = FALSE}, this must already be the transformed basis
#'   with one column per fitted FP coefficient.
#' @param group_fp_powers List of FP power vectors, one element per group. Names
#'   may be absent, group labels, or the temporary transform names generated
#'   internally.
#' @param coefficient_groups List mapping each group to the exact fitted
#'   coefficient names for that group's FP block.
#' @param group_levels Vector of group levels in fitted-model order.
#' @param cont_name Character scalar giving the continuous-variable name.
#' @param transform_basis Logical. If \code{TRUE}, transform \code{cont_mat}
#'   using \code{transform_matrix()}. If \code{FALSE}, treat \code{cont_mat} as
#'   already transformed.
#' @param zero_var Logical scalar. Whether non-positive values should be treated
#'   as structural zeros during FP transformation.
#' @param check_binary Logical scalar passed to \code{transform_matrix()}.
#'
#' @return Numeric matrix with columns named exactly as
#'   \code{unlist(coefficient_groups, use.names = FALSE)}. Attributes
#'   \code{coefficient_groups}, \code{group_fp_powers}, and
#'   \code{group_levels} are attached for downstream checks.
#'
#' @keywords internal
#' @noRd
build_group_fp_basis2 <- function(cont_mat,
                                 group_fp_powers,
                                 coefficient_groups,
                                 group_levels,
                                 cont_name,
                                 transform_basis = TRUE,
                                 zero_var = FALSE,
                                 check_binary = TRUE) {
  
  # Basic scalar validation ----------------------------------------------------
  if (!is.character(cont_name) || length(cont_name) != 1L ||
      is.na(cont_name) || !nzchar(cont_name)) {
    stop(
      "`cont_name` must be a single non-empty character string.",
      call. = FALSE
    )
  }
  
  if (!is.logical(transform_basis) || length(transform_basis) != 1L ||
      anyNA(transform_basis)) {
    stop(
      "`transform_basis` must be a single non-missing logical value.",
      call. = FALSE
    )
  }
  
  if (!is.logical(zero_var) || length(zero_var) != 1L || anyNA(zero_var)) {
    stop(
      "`zero_var` must be a single non-missing logical value.",
      call. = FALSE
    )
  }
  
  if (!is.logical(check_binary) || length(check_binary) != 1L ||
      anyNA(check_binary)) {
    stop(
      "`check_binary` must be a single non-missing logical value.",
      call. = FALSE
    )
  }
  
  if (missing(group_levels) || length(group_levels) == 0L ||
      anyNA(group_levels)) {
    stop(
      "`group_levels` must be a non-empty vector with no missing values.",
      call. = FALSE
    )
  }
  
  group_labels <- as.character(group_levels)
  n_groups     <- length(group_labels)
  temp_names   <- paste0(cont_name, group_labels, "1")
  
  # Validate and normalize coefficient groups ---------------------------------
  if (!is.list(coefficient_groups) || length(coefficient_groups) != n_groups) {
    stop(
      "`coefficient_groups` must be a list with one element per group.",
      call. = FALSE
    )
  }
  
  if (!is.null(names(coefficient_groups)) &&
      anyDuplicated(names(coefficient_groups))) {
    stop(
      "`coefficient_groups` must not contain duplicated names.",
      call. = FALSE
    )
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
          "transform names: ",
          paste(temp_names, collapse = ", "),
          "."
        ),
        call. = FALSE
      )
    }
  }
  
  # Validate and normalize group-specific FP powers ----------------------------
  if (!is.list(group_fp_powers) || length(group_fp_powers) != n_groups) {
    stop(
      "`group_fp_powers` must be a list with one element per group.",
      call. = FALSE
    )
  }
  
  if (!is.null(names(group_fp_powers)) &&
      anyDuplicated(names(group_fp_powers))) {
    stop(
      "`group_fp_powers` must not contain duplicated names.",
      call. = FALSE
    )
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
          "`group_fp_powers` names must match either `group_levels` or the ",
          "temporary transform names: ",
          paste(temp_names, collapse = ", "),
          "."
        ),
        call. = FALSE
      )
    }
  }
  
  # Keep a group-labelled copy for attributes. `transform_matrix()` itself needs
  # temporary names matching the duplicated input columns.
  group_fp_powers_by_group <- stats::setNames(group_fp_powers, group_labels)
  
  expected_group_sizes <- vapply(coefficient_groups, length, integer(1L))
  power_group_sizes    <- vapply(group_fp_powers_by_group, length, integer(1L))
  
  if (any(expected_group_sizes <= 0L)) {
    stop(
      "Each element of `coefficient_groups` must contain at least one coefficient name.",
      call. = FALSE
    )
  }
  
  if (any(power_group_sizes <= 0L)) {
    stop(
      "Each element of `group_fp_powers` must contain at least one FP power.",
      call. = FALSE
    )
  }
  
  bad_power <- vapply(
    group_fp_powers_by_group,
    function(p) {
      !is.numeric(p) || length(p) == 0L || anyNA(p) || any(!is.finite(p))
    },
    logical(1L)
  )
  
  if (any(bad_power)) {
    stop(
      "`group_fp_powers` must contain finite numeric FP power vectors.",
      call. = FALSE
    )
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
    stop(
      "`coefficient_groups` must contain non-empty character coefficient names.",
      call. = FALSE
    )
  }
  
  if (anyDuplicated(x_split_names)) {
    stop(
      "`coefficient_groups` must not contain duplicated coefficient names.",
      call. = FALSE
    )
  }
  
  # Already-transformed path ---------------------------------------------------
  # Developer note:
  # This path is for future prediction code that has already built the exact
  # transformed basis. Do not duplicate the continuous variable in this branch.
  if (!transform_basis) {
    if (!is.matrix(cont_mat) || !is.numeric(cont_mat)) {
      stop(
        "`cont_mat` must be a numeric matrix when `transform_basis = FALSE`.",
        call. = FALSE
      )
    }
    
    if (ncol(cont_mat) != length(x_split_names)) {
      stop(
        paste0(
          "`transform_basis = FALSE` requires `cont_mat` to have ",
          length(x_split_names), " columns, matching `coefficient_groups`."
        ),
        call. = FALSE
      )
    }
    
    if (anyNA(cont_mat) || any(!is.finite(cont_mat))) {
      stop(
        "`cont_mat` contains non-finite values in the supplied FP basis.",
        call. = FALSE
      )
    }
    
    out <- as.matrix(cont_mat)
    storage.mode(out) <- "double"
    colnames(out) <- x_split_names
    
    attr(out, "coefficient_groups") <- coefficient_groups
    attr(out, "group_fp_powers")    <- group_fp_powers_by_group
    attr(out, "group_levels")       <- group_levels
    
    return(out)
  }
  
  # Raw one-column path --------------------------------------------------------
  if (!is.matrix(cont_mat) || ncol(cont_mat) != 1L || !is.numeric(cont_mat)) {
    stop(
      "`cont_mat` must be a one-column numeric matrix when `transform_basis = TRUE`.",
      call. = FALSE
    )
  }
  
  if (is.null(colnames(cont_mat)) || length(colnames(cont_mat)) != 1L) {
    stop(
      "`cont_mat` must have exactly one column name.",
      call. = FALSE
    )
  }
  
  if (anyNA(cont_mat) || any(!is.finite(cont_mat))) {
    stop(
      "`cont_mat` must contain only finite numeric values.",
      call. = FALSE
    )
  }
  
  # Duplicate the continuous variable once per group. The temporary names must
  # match names(power_list) for transform_matrix(). Final coefficient names are
  # assigned only after transformation and validation.
  x_in <- cont_mat[, rep(1L, n_groups), drop = FALSE]
  colnames(x_in) <- temp_names
  
  group_fp_powers_tm <- group_fp_powers_by_group
  names(group_fp_powers_tm) <- temp_names
  
  center_map <- stats::setNames(rep(FALSE, n_groups), temp_names)
  acd_map    <- stats::setNames(rep(FALSE, n_groups), temp_names)
  zero_map   <- stats::setNames(rep(isTRUE(zero_var), n_groups), temp_names)
  
  x_out <- transform_matrix(
    x            = x_in,
    power_list   = group_fp_powers_tm,
    center       = center_map,
    acdx         = acd_map,
    zero         = zero_map,
    check_binary = check_binary
  )$x_transformed
  
  if (is.null(x_out)) {
    stop(
      "FP basis construction returned no transformed columns.",
      call. = FALSE
    )
  }
  
  x_out <- as.matrix(x_out)
  storage.mode(x_out) <- "double"
  
  expected_n_cols <- sum(power_group_sizes)
  
  if (ncol(x_out) != expected_n_cols) {
    stop(
      paste0(
        "Transformed fitted-function basis has ",
        ncol(x_out), " columns, but expected ",
        expected_n_cols, " from group-specific FP powers."
      ),
      call. = FALSE
    )
  }
  
  if (length(x_split_names) != ncol(x_out)) {
    stop(
      "Fitted-function basis columns do not match interaction-model coefficients.",
      call. = FALSE
    )
  }
  
  if (anyNA(x_out) || any(!is.finite(x_out))) {
    stop(
      "Non-finite values were produced in the fitted-function basis.",
      call. = FALSE
    )
  }
  
  # Developer note:
  # transform_matrix() is expected to return transformed columns grouped in the
  # same order as power_list. The checks above verify total column count and
  # per-group block sizes before assigning fitted coefficient names. This makes
  # FP2 and repeated-power alignment explicit.
  colnames(x_out) <- x_split_names
  
  attr(x_out, "coefficient_groups") <- coefficient_groups
  attr(x_out, "group_fp_powers")    <- group_fp_powers_by_group
  attr(x_out, "group_levels")       <- group_levels
  
  x_out
}


# -----------------------------------------------------------------------------
# gen_fitted_values_per_group() -----------------------------------------------
# -----------------------------------------------------------------------------

#' Compute Group-Specific Fitted FP Functions and Confidence Intervals
#'
#' Given a fitted MFPI interaction model and the estimated FP powers for each
#' group, computes the group-specific fitted functions \eqn{\hat f_j(x)}, their
#' pointwise standard errors, confidence intervals, and pairwise differences
#' relative to the reference group.
#'
#' This function is retained for compatibility with the current flex0-flex4
#' fitting path. Long term, this computation should move to \code{predict.mfpi()},
#' and flex functions should return only model metadata rather than precomputed
#' fitted curves.
#'
#' @section Standard errors:
#' Fitted-function standard errors are computed from the relevant submatrix of
#' \code{vcov(interaction_model$fit)}. Difference standard errors are delegated
#' to \code{compute_diff_standard_errors()}, which applies the delta method.
#'
#' @param cont_var One-column numeric matrix of the continuous variable. It must
#'   have a column name and should be on the same scaled working scale passed to
#'   the flex functions.
#' @param group_fp_powers List of FP power vectors, one element per group. For
#'   flex0-flex3, all groups usually share the same powers. For flex4, powers
#'   may differ by group.
#' @param interaction_model Fitted interaction model object returned by
#'   \code{test_interaction()}. Must contain named \code{$coefficients} and a
#'   fitted object in \code{$fit} supporting \code{vcov()}.
#' @param group_var One-column numeric matrix of the grouping variable. It must
#'   have a column name.
#' @param family Regression family object. Retained for interface consistency.
#' @param family_string Character scalar identifying the model family. Cox
#'   models are treated as having no intercept.
#' @param transform Logical. If \code{TRUE}, transform \code{cont_var} using
#'   \code{group_fp_powers}. If \code{FALSE}, \code{cont_var} is assumed to
#'   already be the transformed basis.
#' @param center Logical. Whether to subtract centering constants before fitted
#'   functions are computed.
#' @param center_vals Optional named numeric vector of fit-time centering
#'   constants. This is the preferred path. If \code{NULL} and
#'   \code{center = TRUE}, constants are recomputed from the observed sample.
#' @param center_type Character scalar, \code{"grand"} or \code{"group"}, used
#'   only for fallback centering when \code{center_vals = NULL}.
#' @param use_grid Logical. If \code{TRUE}, evaluate fitted curves on a smooth
#'   200-point grid. This argument is provisional and should move to
#'   \code{predict.mfpi()}.
#' @param scale_var Positive numeric scalar. Multiplies \code{cont_var} before FP
#'   transformation to restore the \eqn{x + shift} scale used by the fitted
#'   coefficients.
#' @param shift_var Numeric scalar. Subtracted from the displayed x-coordinate so
#'   the returned first column is on the raw x scale.
#' @param zero_var Logical scalar. Whether non-positive values are structural
#'   zeros for the interaction variable.
#' @param coefficient_groups Optional list mapping each group to exact fitted
#'   interaction coefficient columns.  \code{var_group()} metadata-based 
#'   coefficient lookup.
#'
#' @return Numeric matrix with columns for raw x, fitted functions, fitted
#'   standard errors, fitted confidence intervals, differences from the reference
#'   group, difference standard errors, and difference confidence intervals.
#'
#' @seealso \code{build_group_fp_basis()},
#'   \code{compute_diff_standard_errors()}, \code{var_group()}
#'
#' @keywords internal
#' @noRd
gen_fitted_values_per_group <- function(cont_var,
                                        group_fp_powers,
                                        interaction_model,
                                        group_var,
                                        family,
                                        family_string,
                                        transform   = TRUE,
                                        center      = FALSE,
                                        center_vals = NULL,
                                        center_type = c("grand", "group"),
                                        use_grid    = FALSE,
                                        scale_var   = 1,
                                        shift_var   = 0,
                                        zero_var    = FALSE,
                                        coefficient_groups = NULL) {
  
  center_type <- match.arg(center_type)
  
  # Input validation -----------------------------------------------------------
  if (!is.logical(transform) || length(transform) != 1L || anyNA(transform)) {
    stop("`transform` must be a single non-missing logical value.",
         call. = FALSE)
  }
  
  if (!is.logical(center) || length(center) != 1L || anyNA(center)) {
    stop("`center` must be a single non-missing logical value.",
         call. = FALSE)
  }
  
  if (!is.logical(use_grid) || length(use_grid) != 1L || anyNA(use_grid)) {
    stop("`use_grid` must be a single non-missing logical value.",
         call. = FALSE)
  }
  
  if (!is.logical(zero_var) || length(zero_var) != 1L || anyNA(zero_var)) {
    stop("`zero_var` must be a single non-missing logical value.",
         call. = FALSE)
  }
  
  if (!is.numeric(scale_var) || length(scale_var) != 1L ||
      anyNA(scale_var) || !is.finite(scale_var) || scale_var <= 0) {
    stop("`scale_var` must be a single finite positive numeric value.",
         call. = FALSE)
  }
  
  if (!is.numeric(shift_var) || length(shift_var) != 1L ||
      anyNA(shift_var) || !is.finite(shift_var)) {
    stop("`shift_var` must be a single finite numeric value.",
         call. = FALSE)
  }
  
  if (!is.matrix(cont_var) || ncol(cont_var) != 1L || !is.numeric(cont_var)) {
    stop("`cont_var` must be a one-column numeric matrix.", call. = FALSE)
  }
  
  cont_name <- colnames(cont_var)
  
  if (is.null(cont_name) || length(cont_name) != 1L) {
    stop("`cont_var` must have exactly one column name.", call. = FALSE)
  }
  
  cont_name <- cont_name[1L]
  
  if (anyNA(cont_var) || any(!is.finite(cont_var))) {
    stop("`cont_var` must contain only finite numeric values.", call. = FALSE)
  }
  
  if (!is.matrix(group_var) || ncol(group_var) != 1L || !is.numeric(group_var)) {
    stop("`group_var` must be a one-column numeric matrix.", call. = FALSE)
  }
  
  group_name <- colnames(group_var)
  
  if (is.null(group_name) || length(group_name) != 1L) {
    stop("`group_var` must have exactly one column name.", call. = FALSE)
  }
  
  group_name <- group_name[1L]
  
  if (anyNA(group_var) || any(!is.finite(group_var))) {
    stop("`group_var` must contain only finite numeric values.", call. = FALSE)
  }
  
  if (nrow(cont_var) != nrow(group_var)) {
    stop(
      "`cont_var` and `group_var` must have the same number of rows.",
      call. = FALSE
    )
  }
  
  if (!is.list(group_fp_powers) || length(group_fp_powers) == 0L) {
    stop(
      "`group_fp_powers` must be a non-empty list.",
      call. = FALSE
    )
  }
  
  if (!center && !is.null(center_vals)) {
    warning(
      "`center_vals` was supplied but `center = FALSE`; centering constants will be ignored.",
      call. = FALSE
    )
  }
  
  # `family` is retained because flex functions pass it consistently. This
  # reconstruction only needs `family_string` for intercept handling.
  invisible(family)
  
  if (!is.character(family_string) || length(family_string) != 1L ||
      is.na(family_string)) {
    stop("`family_string` must be a single character string.", call. = FALSE)
  }
  
  # Extract model coefficients -------------------------------------------------
  coef_vec <- interaction_model$coefficients
  
  if (!is.numeric(coef_vec) || is.null(names(coef_vec))) {
    stop(
      "`interaction_model$coefficients` must be a named numeric vector.",
      call. = FALSE
    )
  }
  
  coef_names <- names(coef_vec)
  
  if (anyDuplicated(coef_names)) {
    stop(
      "`interaction_model$coefficients` must not contain duplicated names.",
      call. = FALSE
    )
  }
  
  # Cox models have no intercept. For other families, include the intercept only
  # if it exists in the fitted coefficient vector.
  has_intercept <- family_string != "cox" && "(Intercept)" %in% coef_names
  intercept     <- if (has_intercept) unname(coef_vec["(Intercept)"]) else 0
  
  if (!is.finite(intercept)) {
    stop(
      "The fitted model intercept is non-finite.",
      call. = FALSE
    )
  }
  
  # Identify the internal group codes used in the fitted interaction model.
  # `group_var` has already been converted to MFPI's internal numeric coding
  # before fitted-function reconstruction is reached. These codes are therefore
  # the correct metadata to use when reconstructing generated FP coefficient names.
  group_vec  <- as.vector(group_var)
  grp_levels <- sort(unique(group_vec))
  n_groups   <- length(grp_levels)
  
  # Identify coefficient blocks by group --------------------------------------
  groups <- if (!is.null(coefficient_groups)) {
    coefficient_groups
  } else {
    # Legacy fallback. Newer code should pass `coefficient_groups` explicitly.
    # The fallback now uses known internal group codes instead of guessing groups
    # from coefficient names with a regular expression.
    var_group(
      var_prefix = cont_name,
      var_names = coef_names,
      group_codes = grp_levels
    )
  }
  
  if (!is.list(groups) || length(groups) == 0L) {
    stop(
      "! Internal error: no coefficient groups were found for fitted-function reconstruction.",
      call. = FALSE
    )
  }
  
  if (length(groups) != n_groups) {
    stop(
      paste0(
        "Number of coefficient groups (", length(groups), ") does not match ",
        "the number of levels in `group_var` (", n_groups, ")."
      ),
      call. = FALSE
    )
  }
  
  if (length(group_fp_powers) != n_groups) {
    stop(
      paste0(
        "`group_fp_powers` has length ", length(group_fp_powers),
        ", but expected one entry per group: ", n_groups, "."
      ),
      call. = FALSE
    )
  }
  
  group_cols <- unname(unlist(groups, use.names = FALSE))
  
  if (!is.character(group_cols) || length(group_cols) == 0L ||
      anyNA(group_cols) || any(!nzchar(group_cols))) {
    stop(
      "! Internal error: coefficient groups must contain non-empty coefficient names.",
      call. = FALSE
    )
  }
  
  if (anyDuplicated(group_cols)) {
    stop(
      "! Internal error: coefficient groups contain duplicated coefficient names.",
      call. = FALSE
    )
  }
  
  missing_group_cols <- setdiff(group_cols, coef_names)
  
  if (length(missing_group_cols) > 0L) {
    stop(
      paste0(
        "! Internal error: interaction model is missing expected FP coefficients: ",
        paste(missing_group_cols, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }
  
  # Backscale cont_var ---------------------------------------------------------
  # cont_var arrives from the flex functions as (x + shift) / scale_var. The
  # fitted interaction coefficients correspond to FP transformations of
  # x + shift, so multiply by scale_var before FP transformation.
  if (scale_var != 1) {
    cont_var <- cont_var * scale_var
  }
  
  # Save the observed sample before optional grid replacement. Fallback centering
  # must be computed from observed rows, not artificial grid rows.
  cont_var_center  <- cont_var
  group_var_center <- group_var
  
  # Optional grid --------------------------------------------------------------
  if (use_grid) {
    cont_vec <- as.vector(cont_var)
    
    if (isTRUE(zero_var)) {
      zero_present <- any(cont_vec <= 0, na.rm = TRUE)
      grid_source  <- cont_vec[cont_vec > 0]
    } else {
      zero_present <- FALSE
      grid_source  <- cont_vec
    }
    
    grid_source <- grid_source[is.finite(grid_source)]
    
    if (length(grid_source) == 0L) {
      msg <- if (isTRUE(zero_var)) {
        "Cannot construct fitted-function grid: no finite positive values are available."
      } else {
        "Cannot construct fitted-function grid: no finite values are available."
      }
      
      stop(msg, call. = FALSE)
    }
    
    grid_range <- range(grid_source, na.rm = TRUE)
    
    grid_values <- seq(
      from       = grid_range[1L],
      to         = grid_range[2L],
      length.out = 200L
    )
    
    if (isTRUE(zero_var) && zero_present) {
      grid_values <- c(0, grid_values)
    }
    
    cont_var <- matrix(grid_values, ncol = 1L)
    colnames(cont_var) <- cont_name
  }
  
  zero_rows_eval <- if (isTRUE(zero_var)) {
    as.vector(cont_var) <= 0
  } else {
    rep(FALSE, nrow(cont_var))
  }
  
  zero_rows_center <- if (isTRUE(zero_var)) {
    as.vector(cont_var_center) <= 0
  } else {
    rep(FALSE, nrow(cont_var_center))
  }
  
  n_obs <- nrow(cont_var)
  
  # Evaluation basis -----------------------------------------------------------
  x_split <- build_group_fp_basis2(
    cont_mat           = cont_var,
    group_fp_powers    = group_fp_powers,
    coefficient_groups = groups,
    group_levels       = grp_levels,
    cont_name          = cont_name,
    transform_basis    = transform,
    zero_var           = zero_var,
    check_binary       = TRUE
  )
  
  missing_x_cols <- setdiff(group_cols, colnames(x_split))
  
  if (length(missing_x_cols) > 0L) {
    stop(
      paste0(
        "! Internal error: fitted-function basis is missing required columns: ",
        paste(missing_x_cols, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }
  
  # Centering ------------------------------------------------------------------
  col_means <- NULL
  
  if (center) {
    if (!is.null(center_vals)) {
      # Preferred path: reuse fit-time centers exactly.
      if (!is.numeric(center_vals) || length(center_vals) == 0L) {
        stop(
          "`center_vals` must be a non-empty numeric vector when supplied.",
          call. = FALSE
        )
      }
      
      if (is.null(names(center_vals))) {
        if (length(center_vals) != ncol(x_split)) {
          stop(
            "`center_vals` must have one value per transformed fitted-function column.",
            call. = FALSE
          )
        }
        
        col_means <- center_vals
        names(col_means) <- colnames(x_split)
      } else {
        missing_centers <- setdiff(colnames(x_split), names(center_vals))
        
        if (length(missing_centers) > 0L) {
          stop(
            paste0(
              "`center_vals` is missing centering constants for: ",
              paste(missing_centers, collapse = ", "),
              "."
            ),
            call. = FALSE
          )
        }
        
        col_means <- center_vals[colnames(x_split)]
      }
      
      if (anyNA(col_means) || any(!is.finite(col_means))) {
        stop(
          "`center_vals` must contain finite, non-missing values.",
          call. = FALSE
        )
      }
    } else {
      # Fallback path: compute centers from observed sample, not grid.
      x_center_split <- build_group_fp_basis(
        cont_mat           = cont_var_center,
        group_fp_powers    = group_fp_powers,
        coefficient_groups = groups,
        group_levels       = grp_levels,
        cont_name          = cont_name,
        transform_basis    = transform,
        zero_var           = zero_var,
        check_binary       = TRUE
      )
      
      if (ncol(x_center_split) != ncol(x_split) ||
          !identical(colnames(x_center_split), colnames(x_split))) {
        stop(
          "! Internal error: fallback-centering basis is not aligned with evaluation basis.",
          call. = FALSE
        )
      }
      
      if (ncol(x_center_split) %% n_groups != 0L) {
        stop(
          "! Internal error: transformed fitted-function columns are not divisible by the number of groups.",
          call. = FALSE
        )
      }
      
      n_fp      <- ncol(x_center_split) %/% n_groups
      col_means <- numeric(ncol(x_center_split))
      names(col_means) <- colnames(x_center_split)
      grp_vec <- as.vector(group_var_center)
      
      for (gi in seq_len(n_groups)) {
        grp_cols <- seq.int((gi - 1L) * n_fp + 1L, gi * n_fp)
        
        if (center_type == "grand") {
          rows_for_centering <- if (isTRUE(zero_var)) {
            !zero_rows_center
          } else {
            rep(TRUE, nrow(x_center_split))
          }
        } else {
          in_grp <- grp_vec == grp_levels[gi]
          
          rows_for_centering <- if (isTRUE(zero_var)) {
            in_grp & !zero_rows_center
          } else {
            in_grp
          }
        }
        
        if (!any(rows_for_centering)) {
          stop(
            "! No valid rows are available for fallback centering in gen_fitted_values_per_group().",
            call. = FALSE
          )
        }
        
        col_means[grp_cols] <- colMeans(
          x_center_split[rows_for_centering, grp_cols, drop = FALSE],
          na.rm = TRUE
        )
      }
      
      if (anyNA(col_means) || any(!is.finite(col_means))) {
        stop(
          "! Could not compute finite fallback centering constants in gen_fitted_values_per_group().",
          call. = FALSE
        )
      }
    }
    
    x_split <- sweep(x_split, 2L, col_means, "-", check.margin = FALSE)
    
    # Centering subtracts non-zero constants from every row. For zero-handled
    # variables, structural-zero observations must remain zero contributions.
    if (isTRUE(zero_var) && any(zero_rows_eval)) {
      x_split[zero_rows_eval, ] <- 0
    }
    
    if (anyNA(x_split) || any(!is.finite(x_split))) {
      stop(
        "! Non-finite values remain in the centered fitted-function basis.",
        call. = FALSE
      )
    }
  }
  
  # Group dummy coefficients ---------------------------------------------------
  group_dummy_names <- character(0L)
  
  if (n_groups > 1L) {
    group_dummy_names <- paste0(group_name, grp_levels[-1L])
    missing_dummies   <- setdiff(group_dummy_names, coef_names)
    
    if (length(missing_dummies) > 0L) {
      stop(
        paste0(
          "! Internal error: interaction model is missing expected group dummy coefficients: ",
          paste(missing_dummies, collapse = ", "),
          "."
        ),
        call. = FALSE
      )
    }
  }
  
  # Variance-covariance matrix and rank-deficiency guards ----------------------
  vcov_mat <- stats::vcov(interaction_model$fit)
  
  if (!is.matrix(vcov_mat) || !is.numeric(vcov_mat) ||
      nrow(vcov_mat) != ncol(vcov_mat)) {
    stop(
      "`vcov(interaction_model$fit)` must be a square numeric matrix.",
      call. = FALSE
    )
  }
  
  if (is.null(rownames(vcov_mat)) || is.null(colnames(vcov_mat))) {
    stop(
      "`vcov(interaction_model$fit)` must have row and column names.",
      call. = FALSE
    )
  }
  
  required_coef_names <- unique(c(
    if (has_intercept) "(Intercept)",
    group_cols,
    group_dummy_names
  ))
  
  missing_coef <- setdiff(required_coef_names, coef_names)
  
  if (length(missing_coef) > 0L) {
    stop(
      paste0(
        "! Internal error: missing required fitted-function coefficients: ",
        paste(missing_coef, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }
  
  missing_vcov <- union(
    setdiff(required_coef_names, rownames(vcov_mat)),
    setdiff(required_coef_names, colnames(vcov_mat))
  )
  
  if (length(missing_vcov) > 0L) {
    stop(
      paste0(
        "! Internal error: covariance matrix is missing required coefficients: ",
        paste(missing_vcov, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }
  
  if (anyNA(coef_vec[required_coef_names]) ||
      any(!is.finite(coef_vec[required_coef_names]))) {
    bad <- required_coef_names[
      is.na(coef_vec[required_coef_names]) |
        !is.finite(coef_vec[required_coef_names])
    ]
    
    stop(
      paste0(
        "! Non-finite coefficients found during fitted-function reconstruction. ",
        "This usually indicates rank deficiency or non-estimable terms: ",
        paste(bad, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }
  
  vcov_required <- vcov_mat[required_coef_names, required_coef_names,
                            drop = FALSE]
  
  if (anyNA(vcov_required) || any(!is.finite(vcov_required))) {
    stop(
      paste0(
        "! Non-finite covariance entries found during fitted-function ",
        "reconstruction. This usually indicates rank deficiency or ",
        "non-estimable terms."
      ),
      call. = FALSE
    )
  }
  
  # Allocate output matrices ---------------------------------------------------
  make_mat <- function(prefix, suffix = "") {
    m <- matrix(0, nrow = n_obs, ncol = n_groups)
    colnames(m) <- paste0(prefix, grp_levels, suffix)
    m
  }
  
  fitted_vals  <- make_mat("f")
  fitted_var   <- make_mat("f", "_var")
  
  crit_val <- stats::qnorm(0.975)
  
  # Compute fitted values and variances per group ------------------------------
  for (i in seq_len(n_groups)) {
    grp_coef_names <- groups[[i]]
    grp_coef       <- coef_vec[grp_coef_names]
    x_grp          <- x_split[, grp_coef_names, drop = FALSE]
    
    x_var       <- x_grp
    x_var_names <- grp_coef_names
    
    if (has_intercept) {
      x_var       <- cbind("(Intercept)" = 1, x_var)
      x_var_names <- c("(Intercept)", x_var_names)
    }
    
    group_coef_offset <- 0
    
    if (i > 1L) {
      group_coef_offset <- unname(coef_vec[group_dummy_names[i - 1L]])
      x_var             <- cbind(x_var, 1)
      x_var_names       <- c(x_var_names, group_dummy_names[i - 1L])
      colnames(x_var)[ncol(x_var)] <- group_dummy_names[i - 1L]
    }
    
    fitted_vals[, i] <- intercept +
      as.vector(x_grp %*% grp_coef) +
      group_coef_offset
    
    vcov_sub <- vcov_mat[x_var_names, x_var_names, drop = FALSE]
    
    fitted_var[, i] <- rowSums((x_var %*% vcov_sub) * x_var)
  }
  
  # Guard against covariance/rank-deficiency artifacts -------------------------
  if (anyNA(fitted_var) || any(!is.finite(fitted_var))) {
    stop(
      "! Non-finite fitted-function variances were produced.",
      call. = FALSE
    )
  }
  
  variance_tol <- sqrt(.Machine$double.eps)
  
  if (any(fitted_var < -variance_tol)) {
    stop(
      paste0(
        "! Negative fitted-function variances were produced. This suggests ",
        "a non-positive-semidefinite covariance matrix or invalid reconstruction."
      ),
      call. = FALSE
    )
  }
  
  fitted_var <- pmax(fitted_var, 0)
  
  # Standard errors and confidence intervals for fitted values -----------------
  fitted_se <- sqrt(fitted_var)
  colnames(fitted_se) <- paste0("se(f", grp_levels, ")")
  
  fitted_lower <- fitted_vals - crit_val * fitted_se
  colnames(fitted_lower) <- paste0("f", grp_levels, "_lower")
  
  fitted_upper <- fitted_vals + crit_val * fitted_se
  colnames(fitted_upper) <- paste0("f", grp_levels, "_upper")
  
  # Function differences relative to reference group ---------------------------
  ref_level <- grp_levels[1L]
  
  if (n_groups > 1L) {
    ref_col <- fitted_vals[, 1L, drop = FALSE]
    
    diff_cols <- fitted_vals[, -1L, drop = FALSE] -
      ref_col[, rep(1L, n_groups - 1L), drop = FALSE]
    
    colnames(diff_cols) <- paste0(
      "f", grp_levels[-1L], "-f", ref_level
    )
    
    diff_se <- compute_diff_standard_errors(
      coefx             = coef_vec,
      cov_betas         = vcov_mat,
      groups            = groups,
      group_name        = group_name,
      group_dummy_names = group_dummy_names,
      xtransformed      = x_split,
      group_fp_powers   = group_fp_powers
    )
    
    if (!is.matrix(diff_se)) {
      diff_se <- as.matrix(diff_se)
    }
    
    if (nrow(diff_se) != n_obs || ncol(diff_se) != n_groups - 1L) {
      stop(
        "! Internal error: difference standard-error matrix has incompatible dimensions.",
        call. = FALSE
      )
    }
    
    if (anyNA(diff_se) || any(!is.finite(diff_se))) {
      stop(
        "! Non-finite difference standard errors were produced.",
        call. = FALSE
      )
    }
    
    colnames(diff_se) <- paste0(
      "se(f", grp_levels[-1L], "-f", ref_level, ")"
    )
    
    diff_lower <- diff_cols - crit_val * diff_se
    colnames(diff_lower) <- paste0(
      "(f", grp_levels[-1L], "-f", ref_level, ")_lower"
    )
    
    diff_upper <- diff_cols + crit_val * diff_se
    colnames(diff_upper) <- paste0(
      "(f", grp_levels[-1L], "-f", ref_level, ")_upper"
    )
  } else {
    diff_cols  <- matrix(nrow = n_obs, ncol = 0L)
    diff_se    <- matrix(nrow = n_obs, ncol = 0L)
    diff_lower <- matrix(nrow = n_obs, ncol = 0L)
    diff_upper <- matrix(nrow = n_obs, ncol = 0L)
  }
  
  # Assemble and return --------------------------------------------------------
  # Fitted values are evaluated on the x + shift scale. The public x-coordinate
  # is returned on the original raw scale.
  cont_var_display <- cont_var - shift_var
  colnames(cont_var_display) <- cont_name
  
  fitted <- cbind(
    cont_var_display,
    fitted_vals,
    fitted_se,
    fitted_lower,
    fitted_upper,
    diff_cols,
    diff_se,
    diff_lower,
    diff_upper
  )
  
  attr(fitted, "fp_centers")         <- if (center) col_means else NULL
  attr(fitted, "center_type")        <- if (center) center_type else NULL
  attr(fitted, "group_fp_powers")    <- group_fp_powers
  attr(fitted, "coefficient_groups") <- groups
  attr(fitted, "scale_var")          <- scale_var
  attr(fitted, "shift_var")          <- shift_var
  attr(fitted, "x_scale")            <- "raw"
  
  fitted
}