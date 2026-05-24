# Utility functions for MFPI
#
# These functions handle dummy variable creation, group-specific FP variable
# construction, and related data-preparation tasks. `create_group_dummies()`,
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
# create_group_dummies() ------------------------------------------------------------
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
#' create_group_dummies(x)
#'
#' # Enforce a fixed set of levels (e.g. for prediction on a subset)
#' x2 <- matrix(c(1, 1, 1), ncol = 1)
#' colnames(x2) <- "trt"
#' create_group_dummies(x2, levels = 1:3)   # trt2 and trt3 columns are all-zero
#'
#' @export
create_group_dummies <- function(x, levels = NULL) {
  
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

# Utility functions for MFPI interaction design matrices
#
# create_z_variables()    — group-specific FP-transformed block matrix
# transform_z_variables() — enumerate all FP power combinations for flex2/flex4


# -----------------------------------------------------------------------------
# create_z_variables() --------------------------------------------------------
# -----------------------------------------------------------------------------

#' Build Group-Specific FP-Transformed Variables for Interaction Modelling
#'
#' For each level of \code{group_var}, constructs a copy of \code{cont_var}
#' that retains the (optionally FP-transformed) covariate values only for
#' observations in that group and is zero elsewhere. The resulting block matrix
#' \eqn{Z \in \mathbb{R}^{n \times KJ}} is used as the interaction design
#' matrix in MFPI, where \eqn{K} is the number of groups and \eqn{J} is the
#' number of FP columns implied by \code{power}.
#'
#' @section Mathematical definition:
#' Let \eqn{f_j(x_i)} denote the \eqn{j}-th FP transformation of observation
#' \eqn{i} under the supplied \code{power} vector. For group \eqn{k}:
#' \deqn{z_{ijk} = \begin{cases} f_j(x_i) & \text{if } g_i = k \\ 0 &
#' \text{otherwise.} \end{cases}}
#'
#' @section Column naming:
#' Columns are named \code{<varname><group><power_index>}, e.g. for
#' \code{cont_var = "age"}, groups 0/1, and \code{power = c(1, 0.5)}: \code{age01},
#' \code{age02}, \code{age11}, \code{age12}.
#'
#' @param cont_var A one-column numeric matrix of the continuous covariate,
#'   with a column name and no missing or non-finite values.
#' @param group_var A one-column numeric matrix of group membership with at
#'   least two distinct values and the same number of rows as \code{cont_var}.
#' @param power Numeric vector of FP powers. Default \code{1} (linear).
#' @param shift Numeric scalar added to \code{cont_var} before transformation.
#'   If \code{NULL}, estimated by \code{mfp2::find_shift_factor()}.
#' @param scale Numeric scalar dividing \code{cont_var} after shifting.
#'   If \code{NULL}, estimated by \code{mfp2::find_scale_factor()}.
#' @param center Logical. If \code{TRUE}, centres both return matrices using
#'   different strategies — see the \emph{Centering} section. Default
#'   \code{FALSE}.
#' @param zero Logical. If \code{TRUE}, non-positive values in \code{cont_var}
#'   are treated as structural zeros: excluded from mean computation when
#'   \code{center = TRUE} and reset to zero after centering. Default
#'   \code{FALSE}.
#'
#' @section Centering:
#' Two distinct strategies are used depending on the return component:
#'
#' \describe{
#'   \item{\code{z} (group-specific interaction matrix)}{Each column belongs
#'     to one group and contains structural zeros for all out-of-group
#'     observations. When \code{center = TRUE}, each column is centred by its
#'     \strong{within-group mean} — the mean of the non-zero (in-group) values
#'     only, via \code{mfp2::center_matrix(zero = TRUE)}.
#'     This ensures group 0\'s column is centred around group 0\'s mean and
#'     group 1\'s column around group 1\'s mean, which is the correct reference
#'     for each group-specific FP curve. Using the pooled mean would be wrong
#'     because the two groups may have different covariate distributions.}
#'   \item{\code{xtransformed} (pooled FP-transformed \code{cont_var})}{Used in
#'     the main-effects model to represent the overall trend. When
#'     \code{center = TRUE}, centred by the \strong{grand mean} across all
#'     transformed values. No structural-zero special-casing is applied here
#'     because \code{xtransformed} is a standard continuous column — any zeros
#'     that may exist after FP transformation are valid transformed values,
#'     not indicators of group absence.}
#' }
#'
#' @return A list with two elements:
#' \describe{
#'   \item{\code{z}}{Numeric matrix of group-specific FP-transformed variables,
#'     \eqn{n \times KJ}. Each column is non-zero only for its own group.
#'     When \code{center = TRUE}, each column is centred by its within-group
#'     mean.}
#'   \item{\code{xtransformed}}{FP-transformed \code{cont_var} pooled across
#'     all observations, \eqn{n \times J}. When \code{center = TRUE}, centred
#'     by the grand mean of all transformed values (standard centering, no
#'     structural-zero adjustment).}
#' }
#'
#' @seealso \code{mfp2::transform_vector_fp()}, \code{transform_z_variables()}
#'
#' @examples
#' cont <- matrix(1:6, ncol = 1); colnames(cont) <- "age"
#' grp  <- matrix(c(0, 1, 0, 1, 0, 1), ncol = 1); colnames(grp) <- "trt"
#' res  <- create_z_variables(cont, grp, power = c(1, 0.5))
#' head(res$z)
#'
#' @export
create_z_variables <- function(cont_var, group_var, power = 1,
                               shift = NULL, scale = NULL,
                               center = FALSE, zero = FALSE) {
  
  # Input validation -----------------------------------------------------------
  if (!is.matrix(cont_var) || ncol(cont_var) != 1L)
    stop("`cont_var` must be a one-column numeric matrix.", call. = FALSE)
  xname <- colnames(cont_var)
  if (is.null(xname))
    stop("`cont_var` must have a column name.", call. = FALSE)
  if (!is.matrix(group_var) || ncol(group_var) != 1L || !is.numeric(group_var))
    stop("`group_var` must be a one-column numeric matrix.", call. = FALSE)
  if (nrow(cont_var) != nrow(group_var))
    stop("`cont_var` and `group_var` must have the same number of rows.",
         call. = FALSE)
  if (!all(is.finite(cont_var)))
    stop("`cont_var` must not contain NA, NaN, or infinite values.", call. = FALSE)
  if (anyNA(group_var))
    stop("`group_var` must not contain NA values.", call. = FALSE)
  if (!is.numeric(power) || length(power) == 0L)
    stop("`power` must be a non-empty numeric vector.", call. = FALSE)
  if (diff(range(cont_var)) == 0)
    warning("`cont_var` has zero variance; all values are identical.", call. = FALSE)
  
  group_vec    <- as.vector(group_var)
  group_levels <- sort(unique(group_vec))
  k            <- length(group_levels)
  
  if (k < 2L)
    stop("`group_var` must have at least two distinct levels.", call. = FALSE)
  
  # FP-transform cont_var (handles linear and non-linear uniformly) -----------
  # zero = zero ensures non-positive values are set to 0 before transformation,
  # consistent with how mfp2:::fit_mfp() estimated the FP powers when zero = TRUE.
  x_fp <- transform_vector_fp(
    x            = cont_var,
    power        = power,
    shift        = shift,
    scale        = scale,
    zero         = zero,
    check_binary = TRUE,
    name         = xname
  )
  
  # Build group-specific block matrix via vectorised index assignment ---------
  # Always built from uncentred x_fp so each column of z can be centred by
  # its own within-group mean (not the pooled mean across all groups).
  n <- nrow(x_fp)
  j <- ncol(x_fp)
  z <- matrix(0, nrow = n, ncol = j * k)
  
  # group_idx maps each observation to its group position (1-based)
  group_idx <- match(group_vec, group_levels)
  
  # Fill all groups at once: for group position i, columns (i-1)*j+1 : i*j
  for (i in seq_len(k)) {
    rows <- group_idx == i
    cols <- seq.int((i - 1L) * j + 1L, i * j)
    z[rows, cols] <- x_fp[rows, , drop = FALSE]
  }
  
  colnames(z) <- sprintf(
    "%s%d%d",
    xname,
    rep(group_levels, each = j),
    rep(seq_len(j),   times = k)
  )
  
  # Centre each column of z by its within-group mean (mean of non-zero values).
  # zero = TRUE tells center_matrix to compute the mean only over non-zero rows,
  # leaving structural zeros (out-of-group observations) at zero.
  if (center) {
    z <- center_matrix(
      mat     = z,
      centers = NULL,
      zero    = setNames(rep(TRUE, ncol(z)), colnames(z))
    )
  }
  
  # xtransformed: pooled FP-transformed cont_var used in the main-effects model.
  # Centred by the grand mean of all transformed values — no structural zero
  # special-casing needed here since this is the pooled continuous term.
  xtransformed <- x_fp
  if (center) {
    fp_colnames <- colnames(x_fp)
    if (is.null(fp_colnames))
      fp_colnames <- paste0(xname, seq_len(ncol(x_fp)))
    xtransformed <- center_matrix(
      mat     = x_fp,
      centers = NULL,
      zero    = setNames(rep(FALSE, ncol(x_fp)), fp_colnames)
    )
  }
  
  list(z = z, xtransformed = xtransformed)
}


# -----------------------------------------------------------------------------
# transform_z_variables() -----------------------------------------------------
# -----------------------------------------------------------------------------

#' Enumerate All FP Power Combinations for Group-Specific Interaction Variables
#'
#' For each combination of FP powers generated by
#' \code{mfp2::generate_powers_fp()}, constructs the full group-specific
#' interaction matrix and collects the results in a list. Used by
#' \code{flex2()} to find the best power combination by deviance minimisation.
#'
#' @section Structural zeros:
#' Group-specific variables contain structural zeros for observations outside
#' each group. Applying log or negative-power FP transforms to zero produces
#' \code{-Inf}/\code{Inf}. The \code{na_replace} argument controls handling:
#' \code{"NA"} (default) marks them missing; \code{"zero"} replaces with 0.
#'
#' @param cont_var One-column numeric matrix of the continuous covariate with
#'   column name; no missing or non-finite values.
#' @param group_var One-column numeric matrix of group membership with at
#'   least two distinct values and the same row count as \code{cont_var}.
#' @param shift,scale Optional numeric scalars for shifting and scaling
#'   \code{cont_var}. \code{NULL} triggers automatic estimation.
#' @param center Logical. Mean-centre transformed columns. Default
#'   \code{FALSE}.
#' @param acdx Logical. Whether to apply the ACD transformation. Default
#'   \code{FALSE}.
#' @param fp_cand Numeric vector of candidate FP powers. Default is the
#'   standard Royston-Altman set \code{c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)}.
#' @param fp_degree Positive integer. FP degree (1 = FP1, 2 = FP2). Default
#'   \code{2}.
#' @param na_replace Character; \code{"zero"} (default) or \code{"NA"}.
#'   Controls replacement of non-finite values arising from structural zeros.
#'   \code{"zero"} produces a matrix ready for model fitting. \code{"NA"}
#'   preserves the distinction between structural zeros and valid zeros such
#'   as \code{log(1) = 0}.
#'
#' @return A named list with four elements:
#' \describe{
#'   \item{\code{z_transformed}}{List of matrices, one per power combination,
#'     each \eqn{n \times KJ}.}
#'   \item{\code{z_untransformed}}{Numeric matrix of group-specific linear
#'     (power = 1) variables before FP transformation. Useful for inspecting
#'     the pre-transformation group indicator columns.}
#'   \item{\code{powers_matrix}}{Numeric matrix; each row is one power
#'     combination.}
#'   \item{\code{znames}}{Column names for the group-specific variables.}
#' }
#'
#' @seealso \code{create_z_variables()}, \code{mfp2::generate_powers_fp()},
#'   \code{mfp2::transform_matrix()}
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
                                  na_replace = c("zero", "NA")) {
  
  na_replace <- match.arg(na_replace)
  
  # Input validation -----------------------------------------------------------
  if (!is.matrix(cont_var) || ncol(cont_var) != 1L)
    stop("`cont_var` must be a one-column numeric matrix.", call. = FALSE)
  if (!is.matrix(group_var) || ncol(group_var) != 1L)
    stop("`group_var` must be a one-column numeric matrix.", call. = FALSE)
  if (nrow(cont_var) != nrow(group_var))
    stop("`cont_var` and `group_var` must have the same number of rows.", call. = FALSE)
  if (!all(is.finite(cont_var)))
    stop("`cont_var` must not contain NA, NaN, or infinite values.", call. = FALSE)
  if (anyNA(group_var))
    stop("`group_var` must not contain NA values.", call. = FALSE)
  if (length(unique(as.vector(group_var))) < 2L)
    stop("`group_var` must have at least two distinct values.", call. = FALSE)
  if (is.null(fp_cand))
    fp_cand <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
  if (!is.numeric(fp_cand) || anyNA(fp_cand))
    stop("`fp_cand` must be a numeric vector with no missing values.", call. = FALSE)
  fp_degree <- as.integer(fp_degree)
  if (length(fp_degree) != 1L || fp_degree < 1L)
    stop("`fp_degree` must be a positive integer (1 or 2).", call. = FALSE)
  
  # Build untransformed group-specific linear placeholder --------------------
  z_untrans <- create_z_variables(
    cont_var  = cont_var,
    group_var = group_var,
    power     = 1L,
    scale     = scale,
    shift     = shift,
    center    = FALSE   # centering applied after FP transform
  )$z
  
  znames   <- colnames(z_untrans)
  n_groups <- length(znames)
  
  # Generate all power combinations -------------------------------------------
  powers_matrix  <- generate_powers_fp(degree = fp_degree, powers = fp_cand)
  n_combinations <- nrow(powers_matrix)
  
  # Per-variable control vectors ----------------------------------------------
  acd_map    <- setNames(rep(acdx,   n_groups), znames)
  center_map <- setNames(rep(center, n_groups), znames)
  zero_map   <- setNames(rep(TRUE,   n_groups), znames)  # structural zeros
  
  # Enumerate transformations with lapply (cleaner than for-loop) -------------
  z_transformed <- lapply(seq_len(n_combinations), function(i) {
    power_set <- setNames(
      replicate(n_groups, powers_matrix[i, ], simplify = FALSE),
      znames
    )
    # Pass zero = TRUE so transform_matrix handles structural zeros correctly,
    # avoiding manual non-finite replacement after the fact.
    transformed <- transform_matrix(
      x          = z_untrans,
      power_list = power_set,
      center     = center_map,
      acdx       = acd_map,
      zero       = zero_map
    )$x_transformed
    
    if (na_replace == "NA") {
      transformed[!is.finite(transformed)] <- NA_real_
    } else {
      transformed[!is.finite(transformed)] <- 0
    }
    
    transformed
  })
  
  list(
    z_transformed   = z_transformed,
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
#' This is a lookup helper used after \code{create_z_variables()} to recover which
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