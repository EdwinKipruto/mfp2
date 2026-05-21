# Utility functions for MFPI
#
# These functions handle dummy variable creation, group-specific FP variable
# construction, and related data-preparation tasks. `create_dummies()`,
# `create_z_variables()`, and `transform_z_variables()` are exported because
# they may be useful to callers building custom interaction models.
# `var_group()` and `adjust_reference_category()` are small helpers that are
# also exported for convenience.
#
# Naming conventions (shared with flex_functions.R and fit_mfpi.R):
#   group_var  — the categorical grouping variable (column name or matrix)
#   cont_var   — the continuous variable being transformed
#   fp_cand    — candidate FP powers for a single variable
#   fp_degree  — degree of the FP (1 = FP1, 2 = FP2)
#   na_replace — how to handle non-finite values from structural zeros


# -----------------------------------------------------------------------------
# create_dummies() ------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Create Dummy Variables from a Single-Column Numeric Matrix
#'
#' Encodes a categorical variable held in a one-column numeric matrix as a set
#' of binary dummy columns. The lowest numeric level is the reference category
#' and is omitted from the output, consistent with the default treatment
#' contrasts used in R's [stats::model.matrix()].
#'
#' @param x A numeric matrix with exactly one column and a column name. Each
#'   distinct numeric value is treated as a category.
#' @param levels Optional numeric vector specifying all expected levels of `x`.
#'   Useful when some levels may be absent from `x` (e.g. in a subset) but
#'   dummy columns are still required for model-matrix consistency. When
#'   supplied, the values in `x` must be a subset of `levels`; any extra values
#'   in `x` not listed in `levels` raise an error.
#'
#' @return A numeric matrix with `nrow(x)` rows and \eqn{K - 1} columns, where
#'   \eqn{K} is the number of unique levels (from `x` or from `levels`). Column
#'   names follow the pattern `<varname><level>`, e.g. `"trt2"`, `"trt3"` when
#'   the input column is named `"trt"` and the non-reference levels are 2 and 3.
#'
#' @examples
#' x <- matrix(c(3, 1, 2, 1, 3, 2), ncol = 1)
#' colnames(x) <- "group"
#' create_dummies(x)
#'
#' # Enforce a fixed set of levels (e.g. for prediction on a subset)
#' x2 <- matrix(c(1, 1, 1), ncol = 1)
#' colnames(x2) <- "trt"
#' create_dummies(x2, levels = 1:3)   # trt2 and trt3 columns are all-zero
#'
#' @export
create_dummies <- function(x, levels = NULL) {
  
  if (!is.matrix(x))         stop("`x` must be a matrix.",                  call. = FALSE)
  if (ncol(x) != 1L)         stop("`x` must have exactly one column.",      call. = FALSE)
  if (!is.numeric(x))        stop("`x` must be numeric.",                   call. = FALSE)
  if (is.null(colnames(x)))  stop("`x` must have a column name.",           call. = FALSE)
  if (nrow(x) == 0L)         stop("`x` must not be empty.",                 call. = FALSE)
  
  xname <- colnames(x)[1L]
  xvec  <- drop(x)
  
  observed_levels <- sort(unique(xvec))
  
  if (is.null(levels)) {
    all_levels <- observed_levels
  } else {
    if (!is.numeric(levels) || anyNA(levels) || !all(is.finite(levels))) {
      stop("`levels` must be a finite numeric vector with no missing values.",
           call. = FALSE)
    }
    all_levels <- sort(unique(levels))
    
    extra_in_data   <- setdiff(observed_levels, all_levels)
    missing_in_data <- setdiff(all_levels, observed_levels)
    
    msgs <- character(0L)
    if (length(extra_in_data) > 0L)
      msgs <- c(msgs, paste0("Values in `x` not listed in `levels`: ",
                             paste(extra_in_data, collapse = ", "), "."))
    if (length(missing_in_data) > 0L)
      msgs <- c(msgs, paste0("Levels not present in `x` (will produce all-zero columns): ",
                             paste(missing_in_data, collapse = ", "), "."))
    if (length(extra_in_data) > 0L)
      stop(paste(msgs, collapse = "\n"), call. = FALSE)
    if (length(missing_in_data) > 0L)
      message(paste(msgs, collapse = "\n"))
  }
  
  if (length(all_levels) < 2L) {
    stop("At least two distinct levels are required to create dummy variables.",
         call. = FALSE)
  }
  
  f   <- factor(xvec, levels = all_levels)
  mat <- model.matrix(~ f)[, -1L, drop = FALSE]
  colnames(mat) <- paste0(xname, all_levels[-1L])
  mat
}


# -----------------------------------------------------------------------------
# create_z_variables() --------------------------------------------------------
# -----------------------------------------------------------------------------

#' Build Group-Specific FP-Transformed Variables for Interaction Modeling
#'
#' For each level of `group_var`, constructs a copy of `cont_var` that retains
#' the (optionally transformed) covariate values only for observations in that
#' group and is zero elsewhere. The resulting block matrix \eqn{Z \in
#' \mathbb{R}^{n \times KJ}} is used as the interaction design matrix in MFPI,
#' where \eqn{K} is the number of groups and \eqn{J} is the number of FP
#' columns implied by `power`.
#'
#' @section Mathematical definition:
#' Let \eqn{f_j(x_i)} denote the \eqn{j}-th FP transformation of observation
#' \eqn{i} under the supplied `power` vector. For group \eqn{k}, the
#' interaction variable is
#' \deqn{z_{ijk} = \begin{cases} f_j(x_i) & \text{if } g_i = k \\ 0 &
#' \text{otherwise,} \end{cases}}
#' producing \eqn{K \times J} columns in total.
#'
#' @section Column naming:
#' Columns are named `<varname><group><power_index>`, where `<group>` is the
#' numeric group level and `<power_index>` is 1-based. For example, if
#' `cont_var = "age"`, groups are 0 and 1, and `power = c(1, 0.5)`, the
#' columns are `age01`, `age02`, `age11`, `age12`.
#'
#' @param cont_var A one-column numeric matrix representing the continuous
#'   covariate. Must have a column name and no missing or non-finite values.
#' @param group_var A one-column numeric matrix representing group membership.
#'   Must have at least two distinct values and the same number of rows as
#'   `cont_var`.
#' @param power Numeric vector of FP powers applied to `cont_var`. Default is
#'   `1` (linear, no transformation beyond optional shift/scale). A length-2
#'   vector (e.g. `c(0.5, 2)`) produces two columns per group.
#' @param shift Numeric scalar. Added to `cont_var` before transformation. If
#'   `NULL`, estimated automatically by [mfp2::find_shift_factor()].
#' @param scale Numeric scalar. Divides `cont_var` after shifting. If `NULL`,
#'   estimated automatically by [mfp2::find_scale_factor()].
#' @param center Logical. If `TRUE`, each transformed column is mean-centred
#'   (using the column mean of the transformed values across all observations).
#'
#' @return A list with two elements:
#' \describe{
#'   \item{`z`}{Numeric matrix of group-specific transformed variables,
#'     \eqn{n \times KJ}.}
#'   \item{`xtransformed`}{Numeric matrix of the FP-transformed (and optionally
#'     centred) version of `cont_var`, \eqn{n \times J}.}
#' }
#'
#' @seealso [mfp2::transform_vector_fp()], [mfp2::transform_z_variables()]
#'
#' @examples
#' cont <- matrix(1:6, ncol = 1); colnames(cont) <- "age"
#' grp  <- matrix(c(0, 1, 0, 1, 0, 1), ncol = 1); colnames(grp) <- "trt"
#' res  <- create_z_variables(cont, grp, power = c(1, 0.5))
#' head(res$z)
#'
#' @export
create_z_variables <- function(cont_var, group_var, power = 1,
                               shift = NULL, scale = NULL, center = FALSE) {
  
  # Input validation -----------------------------------------------------------
  if (!is.matrix(cont_var) || ncol(cont_var) != 1L)
    stop("`cont_var` must be a one-column numeric matrix.", call. = FALSE)
  xname <- colnames(cont_var)
  if (is.null(xname))
    stop("`cont_var` must have a column name.", call. = FALSE)
  if (!is.matrix(group_var) || ncol(group_var) != 1L || !is.numeric(group_var))
    stop("`group_var` must be a one-column numeric matrix.", call. = FALSE)
  if (nrow(cont_var) != nrow(group_var))
    stop("`cont_var` and `group_var` must have the same number of rows.", call. = FALSE)
  if (anyNA(cont_var) || !all(is.finite(cont_var)))
    stop("`cont_var` must not contain NA, NaN, or infinite values.", call. = FALSE)
  if (anyNA(group_var))
    stop("`group_var` must not contain NA values.", call. = FALSE)
  if (!is.numeric(power) || length(power) == 0L)
    stop("`power` must be a non-empty numeric vector.", call. = FALSE)
  if (var(as.vector(cont_var)) == 0)
    warning("`cont_var` has zero variance; all values are identical.", call. = FALSE)
  if (any(power <= 0) && any(cont_var <= 0))
    stop(
      "`cont_var` must be strictly positive when `power` includes 0 or negative values.",
      call. = FALSE
    )
  
  group_vec    <- as.vector(group_var)
  group_levels <- sort(unique(group_vec))
  k            <- length(group_levels)
  
  if (k < 2L)
    stop("`group_var` must have at least two distinct levels.", call. = FALSE)
  
  # Apply FP transformation to cont_var ----------------------------------------
  if (!identical(power, 1) && !identical(power, 1L)) {
    x_fp <- transform_vector_fp(
      cont_var,
      power        = power,
      shift        = shift,
      scale        = scale,
      check_binary = TRUE,
      name         = xname
    )
  } else {
    # Linear: apply shift and scale manually to stay consistent with mfp2
    s    <- if (!is.null(shift)) shift else find_shift_factor(cont_var)
    sc   <- if (!is.null(scale)) scale else find_scale_factor(cont_var + s)
    x_fp <- (cont_var + s) / sc
  }
  
  if (center) {
    x_fp <- sweep(x_fp, 2L, colMeans(x_fp, na.rm = TRUE), "-")
  }
  
  # Build group-specific block matrix ------------------------------------------
  n   <- nrow(x_fp)
  j   <- ncol(x_fp)
  z   <- matrix(0, nrow = n, ncol = j * k)
  
  for (i in seq_along(group_levels)) {
    rows <- which(group_vec == group_levels[i])
    cols <- seq.int((i - 1L) * j + 1L, i * j)
    z[rows, cols] <- x_fp[rows, , drop = FALSE]
  }
  
  colnames(z) <- sprintf(
    "%s%d%d",
    xname,
    rep(group_levels, each = j),
    rep(seq_len(j),   times = k)
  )
  
  list(z = z, xtransformed = x_fp)
}


# -----------------------------------------------------------------------------
# transform_z_variables() -----------------------------------------------------
# -----------------------------------------------------------------------------

#' Enumerate All FP Power Combinations for Group-Specific Interaction Variables
#'
#' For each combination of FP powers generated by [mfp2::generate_powers_fp()],
#' constructs the full group-specific interaction matrix via
#' [mfp2::transform_matrix()] and collects the results in a list. This is the
#' engine behind `flex2` and `flex4`, where the best power combination is
#' selected by comparing model deviances.
#'
#' @section Structural zeros and non-finite values:
#' Group-specific variables contain structural zeros for observations outside
#' the relevant group. Applying log or negative-power FP transforms to zero
#' produces `-Inf` or `Inf`. These are not data errors; they arise purely from
#' the construction of the group indicator. The `na_replace` argument controls
#' how they are handled: `"NA"` (default) marks them as missing so they can be
#' distinguished from valid zeros such as `log(1) = 0`; `"zero"` replaces them
#' with 0 for fitting routines that cannot tolerate `NA`.
#'
#' @param cont_var A one-column numeric matrix of the continuous covariate.
#'   Must have a column name; no missing or non-finite values permitted.
#' @param group_var A one-column numeric matrix of group membership with at
#'   least two distinct values and the same row count as `cont_var`.
#' @param shift Optional numeric scalar. Shift applied to `cont_var` before
#'   transformation. `NULL` triggers automatic estimation.
#' @param scale Optional numeric scalar. Scale divisor applied after shifting.
#'   `NULL` triggers automatic estimation.
#' @param center Logical. Whether to mean-centre each transformed column.
#' @param acdx Logical. Whether `cont_var` undergoes the ACD transformation.
#'   Default `FALSE`.
#' @param fp_cand Numeric vector of candidate FP powers. Default is the
#'   standard Royston-Altman set `c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)`.
#' @param fp_degree Positive integer. FP degree (1 = FP1, 2 = FP2). Default
#'   `2`.
#' @param na_replace Character string; `"NA"` (default) or `"zero"`. Controls
#'   replacement of non-finite values arising from structural zeros.
#'
#' @return A named list with four elements:
#' \describe{
#'   \item{`z_transformed`}{List of numeric matrices, one per power
#'     combination. Each matrix has \eqn{n} rows and \eqn{K \times J} columns,
#'     where \eqn{K} is the number of groups and \eqn{J} is `fp_degree`.}
#'   \item{`z_untransformed`}{Numeric matrix of group-specific linear (power =
#'     1) variables before FP transformation, used for downstream power
#'     application.}
#'   \item{`powers_matrix`}{Numeric matrix where each row is one FP power
#'     combination tried.}
#'   \item{`znames`}{Character vector of column names for the group-specific
#'     variables.}
#' }
#'
#' @seealso [mfp2::create_z_variables()], [mfp2::generate_powers_fp()],
#'   [mfp2::transform_matrix()]
#'
#' @examples
#' cont <- matrix(1:6, ncol = 1); colnames(cont) <- "age"
#' grp  <- matrix(c(0, 1, 0, 1, 0, 1), ncol = 1); colnames(grp) <- "trt"
#' res  <- transform_z_variables(cont, grp, fp_cand = c(1, 0.5), fp_degree = 2)
#' str(res$z_transformed[[1]])
#'
#' @importFrom mfp2 transform_matrix generate_powers_fp
#' @export
transform_z_variables <- function(cont_var, group_var,
                                  shift = NULL, scale = NULL,
                                  center = FALSE, acdx = FALSE,
                                  fp_cand = NULL, fp_degree = 2L,
                                  na_replace = c("NA", "zero")) {
  
  na_replace <- match.arg(na_replace)
  
  # Input validation -----------------------------------------------------------
  if (!is.matrix(cont_var) || ncol(cont_var) != 1L)
    stop("`cont_var` must be a one-column numeric matrix.", call. = FALSE)
  if (!is.matrix(group_var) || ncol(group_var) != 1L)
    stop("`group_var` must be a one-column numeric matrix.", call. = FALSE)
  if (nrow(cont_var) != nrow(group_var))
    stop("`cont_var` and `group_var` must have the same number of rows.", call. = FALSE)
  if (anyNA(cont_var) || !all(is.finite(cont_var)))
    stop("`cont_var` must not contain NA, NaN, or infinite values.", call. = FALSE)
  if (anyNA(group_var))
    stop("`group_var` must not contain NA values.", call. = FALSE)
  if (length(unique(as.vector(group_var))) < 2L)
    stop("`group_var` must have at least two distinct values.", call. = FALSE)
  
  if (is.null(fp_cand)) {
    fp_cand <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
  }
  if (!is.numeric(fp_cand) || anyNA(fp_cand))
    stop("`fp_cand` must be a numeric vector with no missing values.", call. = FALSE)
  if (!is.numeric(fp_degree) || length(fp_degree) != 1L || fp_degree < 1L)
    stop("`fp_degree` must be a positive integer (1 or 2).", call. = FALSE)
  
  # Build untransformed group-specific variables (linear placeholder) ----------
  z_untrans <- create_z_variables(
    cont_var  = cont_var,
    group_var = group_var,
    power     = 1L,
    scale     = scale,
    shift     = shift,
    center    = FALSE   # centering applied below, after FP transform
  )$z
  
  znames   <- colnames(z_untrans)
  n_groups <- length(znames)
  
  # Generate all power combinations for this degree ---------------------------
  powers_matrix    <- mfp2::generate_powers_fp(degree = fp_degree, powers = fp_cand)
  n_combinations   <- nrow(powers_matrix)
  
  # Per-variable transform control vectors ------------------------------------
  acd_map    <- setNames(rep(acdx,   n_groups), znames)
  center_map <- setNames(rep(center, n_groups), znames)
  
  # Enumerate transformations --------------------------------------------------
  z_trans_list <- vector("list", length = n_combinations)
  
  for (i in seq_len(n_combinations)) {
    power_set <- setNames(
      replicate(n_groups, powers_matrix[i, ], simplify = FALSE),
      znames
    )
    
    transformed <- transform_matrix(
      x          = z_untrans,
      power_list = power_set,
      center     = center_map,
      acdx       = acd_map
    )$x_transformed
    
    # Handle non-finite values from structural zeros -------------------------
    non_finite <- !is.finite(transformed)
    transformed[non_finite] <- if (na_replace == "NA") NA_real_ else 0
    
    z_trans_list[[i]] <- transformed
  }
  
  list(
    z_transformed   = z_trans_list,
    z_untransformed = z_untrans,
    powers_matrix   = powers_matrix,
    znames          = znames
  )
}


# -----------------------------------------------------------------------------
# var_group() -----------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Partition Variable Names into Group-Specific Subsets
#'
#' Given a common variable prefix (e.g. `"age"`) and a character vector of
#' column names, returns a list where each element contains the column names
#' belonging to one group. Group identity is read from the single digit
#' immediately following the prefix.
#'
#' This is a lookup helper used after [mfp2::create_z_variables()] to recover which
#' columns belong to each level of `group_var`.
#'
#' @param var_prefix Character string. The common prefix of the variable names
#'   to partition (e.g. `"age"` matches `"age01"`, `"age02"`, `"age11"`, ...).
#' @param var_names Character vector of all column names to search.
#'
#' @return A list of character vectors, one per unique group digit found
#'   immediately after `var_prefix`.
#'
#' @examples
#' var_group("age", c("age01", "age02", "age11", "age12", "bmi"))
#' # Returns: list(c("age01", "age02"), c("age11", "age12"))
#'
#' @export
var_group <- function(var_prefix, var_names) {
  if (!is.character(var_prefix) || length(var_prefix) != 1L)
    stop("`var_prefix` must be a single character string.", call. = FALSE)
  if (!is.character(var_names) || length(var_names) == 0L)
    stop("`var_names` must be a non-empty character vector.", call. = FALSE)
  
  matched <- grep(paste0("^", var_prefix, "\\d+"), var_names, value = TRUE)
  if (length(matched) == 0L) return(list())
  
  # The group digit sits immediately after the prefix
  prefix_len  <- nchar(var_prefix)
  group_digit <- substr(matched, prefix_len + 1L, prefix_len + 1L)
  unique_grps <- unique(group_digit)
  
  lapply(unique_grps, function(g) {
    grep(paste0("^", var_prefix, g), matched, value = TRUE)
  })
}


# -----------------------------------------------------------------------------
# adjust_reference_category() -------------------------------------------------
# -----------------------------------------------------------------------------

#' Reverse or Relevel the Reference Category of a Factor
#'
#' Reverses the level order of a factor variable and, optionally, sets the
#' reference to the top (highest) or bottom (lowest) level. This is a
#' convenience wrapper around `relevel()` for situations where the
#' default alphabetical or numerical level ordering is not appropriate.
#'
#' When `x` is not a factor, it is returned unchanged.
#'
#' @param x A factor vector, or any other vector (returned as-is if not a
#'   factor).
#' @param use_top_as_ref Logical. If `TRUE`, the highest level (after
#'   reversal) is used as the reference category. If `FALSE` (default), the
#'   lowest level is the reference.
#'
#' @return A factor with the same values as `x` but with the level order
#'   reversed and reference category updated, or the original `x` unchanged if
#'   it is not a factor.
#'
#' @examples
#' x <- factor(c("Low", "Medium", "High", "Low", "High"))
#' adjust_reference_category(x)                       # reference: "High" (first after reversal)
#' adjust_reference_category(x, use_top_as_ref = TRUE)  # reference: "Low" (last after reversal)
#'
#' @importFrom utils tail
#' @export
adjust_reference_category <- function(x, use_top_as_ref = FALSE) {
  if (!is.factor(x)) return(x)
  
  levels(x) <- rev(levels(x))
  ref <- if (use_top_as_ref) tail(levels(x), 1L) else levels(x)[1L]
  relevel(x, ref = ref)
}
