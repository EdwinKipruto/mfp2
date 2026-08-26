# -----------------------------------------------------------------------------
# create_z_variables() ---------------------------------------------------------
# -----------------------------------------------------------------------------


#' Precompute Group Row Indices for MFPI Design Construction
#'
#' Builds reusable row-index lists for group-specific MFPI design matrices.
#' This avoids repeatedly constructing logical masks such as `group_idx == g`
#' and `group_idx == g & zero_rows` inside hot transformation loops.
#'
#' @param group_idx Integer vector giving the group index for each row.
#' @param n_groups Integer scalar; number of groups.
#' @param zero_rows Logical vector identifying structural-zero rows.
#' @param valid_rows Logical vector identifying rows used for FP centering.
#' @param zero Logical scalar; whether structural-zero handling is active.
#'
#' @return A list with `rows_by_group`, `valid_by_group`, and `zero_by_group`.
#'
#' @keywords internal
#' @noRd
mfpi_group_row_cache <- function(group_idx,
                                 n_groups,
                                 zero_rows,
                                 valid_rows,
                                 zero) {
  rows_by_group <- unname(split(
    seq_along(group_idx),
    factor(group_idx, levels = seq_len(n_groups))
  ))

  valid_by_group <- lapply(
    rows_by_group,
    function(rows) {
      rows[valid_rows[rows]]
    }
  )

  zero_by_group <- if (isTRUE(zero)) {
    lapply(
      rows_by_group,
      function(rows) {
        rows[zero_rows[rows]]
      }
    )
  } else {
    vector("list", n_groups)
  }

  list(
    rows_by_group  = rows_by_group,
    valid_by_group = valid_by_group,
    zero_by_group  = zero_by_group
  )
}

#' Build Group-Specific Fractional-Polynomial Interaction Variables
#'
#' Constructs a group-specific fractional-polynomial design matrix for modelling
#' interactions between a continuous covariate and a categorical or treatment
#' grouping variable.
#'
#' The function first transforms the full continuous variable using the supplied
#' fractional-polynomial powers. It then inserts the transformed values into
#' group-specific blocks of the design matrix. Values outside the corresponding
#' group are kept as structural zeros.
#'
#' This implementation deliberately avoids transforming a zero-padded
#' group-specific matrix. Instead, the full continuous variable is transformed
#' before group masking. This is important for log and negative FP powers because
#' artificial out-of-group zeros would otherwise generate invalid values such as
#' \code{Inf}, \code{-Inf}, or \code{NaN}.
#'
#' @section Mathematical definition:
#'
#' Let \eqn{x_i} denote the continuous covariate for observation \eqn{i}, and let
#' \eqn{G_i \in \{g_1, \ldots, g_K\}} denote its group membership.
#'
#' Let \eqn{B(x_i)} be the FP basis generated from \eqn{x_i}. For FP1,
#' \eqn{B(x_i)} has one column. For FP2, it has two columns. For example, for
#' powers \eqn{p_1} and \eqn{p_2},
#'
#' \deqn{
#' B(x_i) =
#' \begin{cases}
#' \left(x_i^{p_1}, x_i^{p_2}\right), & p_1 \ne p_2, \\
#' \left(x_i^{p_1}, x_i^{p_1}\log(x_i)\right), & p_1 = p_2 \ne 0, \\
#' \left(\log(x_i), \log(x_i)^2\right), & p_1 = p_2 = 0.
#' \end{cases}
#' }
#'
#' The group-specific design block for group \eqn{g} is
#'
#' \deqn{
#' Z_g(x_i) = \{B(x_i) - c_g\} I(G_i = g),
#' }
#'
#' where \eqn{c_g} is the centering vector for group \eqn{g}. The final design
#' matrix is obtained by concatenating all group-specific blocks:
#'
#' \deqn{
#' Z = \left[Z_{g_1}, Z_{g_2}, \ldots, Z_{g_K}\right].
#' }
#'
#' Rows outside a given group remain structural zeros in that group's block.
#'
#' @section FP transformations:
#'
#' FP transformations are delegated to \code{transform_vector_fp()}.
#' The power convention is:
#'
#' \describe{
#'   \item{\code{p = 1}}{Linear term, \eqn{x}.}
#'   \item{\code{p = 0}}{Log term, \eqn{\log(x)}.}
#'   \item{\code{p != 0}}{Power term, \eqn{x^p}.}
#'   \item{Repeated powers}{For repeated powers, the second and later terms are
#'   multiplied by powers of \eqn{\log(x)} according to the standard FP rule.}
#' }
#'
#' If the continuous variable contains non-positive values, it should usually be
#' shifted before calling this function. With \code{zero = TRUE}, the variable
#' must be nonnegative and only exact zeros are structural zeros.
#'
#' The continuous variable must contain at least two distinct values. A
#' constant variable cannot define an estimable group-specific FP effect:
#' centering turns its FP basis into zero columns, while an uncentered constant
#' basis is collinear with the corresponding group indicators.
#'
#' @section Centering:
#'
#' If \code{center = FALSE}, no centering is applied and \code{center_vals} is
#' returned as \code{NULL}.
#'
#' If \code{center = TRUE}, centering is applied to the FP basis before insertion
#' into the active group rows. The inactive rows of each group-specific block
#' remain structural zeros.
#'
#' The centering strategy is controlled by \code{center_type}:
#'
#' \describe{
#'   \item{\code{"grand"}}{
#'   Each FP basis column is centered by its empirical mean over all valid
#'   observations:
#'
#'   \deqn{
#'   c_g = n^{-1}\sum_i B(x_i).
#'   }
#'
#'   The same centering vector is used for every group. This is useful when all
#'   group-specific functions should be centered against the same empirical
#'   covariate distribution.
#'   }
#'
#'   \item{\code{"group"}}{
#'   Each group-specific FP basis is centered by its empirical mean within that
#'   group:
#'
#'   \deqn{
#'   c_g = n_g^{-1}\sum_{i:G_i=g} B(x_i).
#'   }
#'
#'   This makes each group-specific function centered relative to its own group
#'   covariate distribution.
#'   }
#' }
#'
#' When \code{zero = TRUE}, exact-zero values of \code{cont_var} are excluded
#' from the centering means and reset to zero after centering. Negative values
#' are invalid and must be recoded explicitly if they should represent zero.
#'
#' @section Fitting matrix versus evaluation matrix:
#'
#' The returned \code{z} matrix is a model-fitting design matrix. It contains
#' structural zeros outside each group:
#'
#' \deqn{
#' Z_g(x_i) = 0 \quad \text{if } G_i \ne g.
#' }
#'
#' If \code{return_eval = TRUE}, the function also returns \code{x_eval}. This
#' matrix evaluates every group-specific function for every observed value of
#' \code{cont_var}:
#'
#' \deqn{
#' X_{\mathrm{eval},g}(x_i) = B(x_i) - c_g.
#' }
#'
#' Unlike \code{z}, \code{x_eval} does not contain out-of-group structural zeros.
#' It is therefore the appropriate matrix for fitted-function evaluation,
#' plotting, and prediction over a common set of \eqn{x} values.
#'
#' @section Column naming:
#'
#' Columns are named using the convention:
#'
#' \deqn{
#' \texttt{<varname><group><power_index>}.
#' }
#'
#' For example, if \code{cont_var} is named \code{"age"}, the group levels are
#' \code{0} and \code{1}, and \code{power = c(1, 0.5)}, the generated column
#' names are:
#'
#' \preformatted{
#' age01 age02 age11 age12
#' }
#'
#' @section Stored metadata:
#'
#' The returned design matrices contain attributes used by downstream fitted
#' function and prediction code:
#'
#' \describe{
#'   \item{\code{fp_centers}}{The centering constants used for each
#'   group-specific FP column.}
#'   \item{\code{center}}{Whether centering was applied.}
#'   \item{\code{center_type}}{The centering strategy.}
#'   \item{\code{group_levels}}{The observed group levels.}
#'   \item{\code{group_fp_powers}}{A named list of powers, one element per
#'   group.}
#'   \item{\code{power}}{The FP power vector supplied to the function.}
#' }
#'
#' These constants should be reused during prediction. They should not be
#' recomputed on new prediction data.
#'
#' @param cont_var A one-column numeric matrix containing the continuous
#'   covariate. It must have a column name, contain only finite values, and
#'   contain at least two distinct values.
#' @param group_var A one-column numeric matrix containing group membership. It
#'   must have the same number of rows as \code{cont_var}, contain no missing
#'   values, and have at least two distinct values.
#' @param power Numeric vector of FP powers. Length one gives an FP1 basis;
#'   length two gives an FP2 basis. Repeated powers are handled by
#'   \code{transform_vector_fp()}.
#' @param shift Optional numeric shift applied before FP transformation. Passed
#'   to \code{transform_vector_fp()}.
#' @param scale Optional numeric scale applied before FP transformation. Passed
#'   to \code{transform_vector_fp()}.
#' @param center Logical. If \code{TRUE}, center the FP basis columns before
#'   constructing the group-specific design matrix. Default is \code{FALSE}.
#' @param zero Logical. If \code{TRUE}, exact-zero values of nonnegative
#'   \code{cont_var} are treated as structural zeros: they are excluded from
#'   centering means and reset to zero after centering. Callers must supply
#'   previously validated nonnegative values. Default is \code{FALSE}.
#' @param center_type Character string controlling centering when
#'   \code{center = TRUE}. Either \code{"grand"} or \code{"group"}. Default is
#'   \code{"grand"}.
#' @param return_eval Logical. If \code{TRUE}, also return an evaluation matrix
#'   in which every group-specific function is evaluated for every observed
#'   \code{x} value. Default is \code{FALSE}.
#' @param scale_var Numeric scalar. The scale factor for \code{cont_var},
#'   as computed in \code{mfpi.default()}. Before any FP transformation,
#'   \code{cont_var} is multiplied by \code{scale_var} to restore the
#'   variable before scaling was applied. For standard variables this gives
#'   \eqn{x + \text{shift}}; for \code{zero_var = TRUE} variables (where
#'   shift is forced to 0) this gives the original \eqn{x}, and the
#'   \code{zero} argument then correctly handles exact-zero values.
#'   This ensures model coefficients match the \eqn{\phi(x + \text{shift})}
#'   scale of the adjustment model and standalone \pkg{mfp2}. Default
#'   \code{1} (no backscaling).
#'
#' @return A list with the following elements:
#'
#' \describe{
#'   \item{\code{z}}{Group-specific fitting design matrix. Rows outside the
#'   corresponding group are structural zeros.}
#'   \item{\code{xtransformed}}{Pooled FP-transformed version of
#'   \code{cont_var}. If \code{center = TRUE}, this is centered using grand
#'   means.}
#'   \item{\code{x_eval}}{Group-specific evaluation matrix if
#'   \code{return_eval = TRUE}; otherwise \code{NULL}.}
#'   \item{\code{center_vals}}{Named numeric vector of centering constants, or
#'   \code{NULL} when \code{center = FALSE}.}
#'   \item{\code{x_fp}}{Full FP-transformed continuous variable before group
#'   masking and before centering.}
#' }
#'
#' @section Usage with raw vs pre-shifted/scaled data:
#' This function supports two calling conventions depending on whether
#' \code{cont_var} has been pre-processed upstream:
#'
#' \describe{
#'   \item{Raw data (direct user call)}{Pass the original unscaled,
#'     unshifted \code{cont_var} and let \code{shift}, \code{scale}, and
#'     \code{zero} control the FP transformation:
#'     \preformatted{
#' create_z_variables(
#'   cont_var  = x_raw,
#'   group_var = grp,
#'   power     = c(0.5, 1),
#'   shift     = NULL,   # estimated automatically
#'   scale     = NULL,   # estimated automatically
#'   scale_var = 1       # default: no backscaling
#' )}}
#'   \item{Pre-shifted and scaled data (internal mfpi use)}{When called
#'     from the flex functions inside \code{mfpi()}, \code{cont_var} has
#'     already been shifted and scaled in \code{mfpi.default()}:
#'     \eqn{x_{\text{scaled}} = (x + \text{shift}) / \text{scale\_var}}.
#'     Pass \code{shift = 0}, \code{scale = 1}, and the real
#'     \code{scale_var} so that the function backscales before FP
#'     transformation and coefficients are on the
#'     \eqn{\phi(x + \text{shift})} scale:
#'     \preformatted{
#' create_z_variables(
#'   cont_var  = x_scaled,   # (x + shift) / scale_var
#'   group_var = grp,
#'   power     = bestfp,
#'   shift     = 0,
#'   scale     = 1,
#'   scale_var = scale_var   # real scale factor: x_scaled * scale_var = x + shift
#' )}}
#'   \item{\code{column_groups}}{Named list mapping each group level to the
#'   exact generated FP coefficient columns for that group. Used downstream to
#'   avoid reparsing coefficient names with regular expressions.}
#' }
#'
#' @seealso \code{transform_vector_fp}
#'
#' @examples
#' cont <- matrix(1:6, ncol = 1)
#' colnames(cont) <- "age"
#'
#' grp <- matrix(c(0, 1, 0, 1, 0, 1), ncol = 1)
#' colnames(grp) <- "trt"
#'
#' # FP2 group-specific design, grand centered
#' create_z_variables(
#'   cont_var = cont,
#'   group_var = grp,
#'   power = c(1, 0.5),
#'   center = TRUE,
#'   center_type = "grand"
#' )
#'
#' # Also return an evaluation matrix for fitted-function plotting
#' create_z_variables(
#'   cont_var = cont,
#'   group_var = grp,
#'   power = c(1, 0.5),
#'   center = TRUE,
#'   return_eval = TRUE
#' )
#'
#' @importFrom stats setNames
#' @keywords internal
#' @noRd
create_z_variables <- function(cont_var,
                               group_var,
                               power = 1,
                               shift = NULL,
                               scale = NULL,
                               center = FALSE,
                               zero = FALSE,
                               center_type = c("grand", "group"),
                               return_eval = FALSE,
                               scale_var = 1) {

  center_type <- match.arg(center_type)

  # Input validation -----------------------------------------------------------

  if (!is.matrix(cont_var) || ncol(cont_var) != 1L || !is.numeric(cont_var)) {
    stop("cont_var must be a one-column numeric matrix.", call. = FALSE)
  }

  xname <- colnames(cont_var)

  if (is.null(xname) || length(xname) != 1L) {
    stop("cont_var must have exactly one column name.", call. = FALSE)
  }

  if (!is.matrix(group_var) || ncol(group_var) != 1L || !is.numeric(group_var)) {
    stop("group_var must be a one-column numeric matrix.", call. = FALSE)
  }

  if (nrow(cont_var) != nrow(group_var)) {
    stop("cont_var and group_var must have the same number of rows.",
         call. = FALSE)
  }

  if (!all(is.finite(cont_var))) {
    stop("cont_var must not contain NA, NaN, or infinite values.", call. = FALSE)
  }

  if (!all(is.finite(group_var))) {
    stop("group_var must not contain NA, NaN, or infinite values.", call. = FALSE)
  }

  if (!is.numeric(power) || length(power) == 0L || anyNA(power)) {
    stop("power must be a non-empty numeric vector with no missing values.",
         call. = FALSE)
  }

  group_vec    <- as.vector(group_var)
  group_levels <- sort(unique(group_vec))
  group_idx    <- match(group_vec, group_levels)
  n_groups     <- length(group_levels)

  if (n_groups < 2L) {
    stop("group_var must have at least two distinct values.", call. = FALSE)
  }

  n <- nrow(cont_var)

  if (diff(range(as.vector(cont_var))) == 0) {
    # A constant covariate cannot identify an MFPI effect. With centering its
    # transformed basis becomes exactly zero; without centering each group-
    # specific column is only a scaled copy of a group indicator. Stop before
    # either rank-deficient design can reach the downstream fitting engine.
    stop(
      "cont_var must contain at least two distinct values; a constant ",
      "continuous variable cannot define an estimable MFPI effect.",
      call. = FALSE
    )
  }

  # Backscale cont_var to match mfp2 convention --------------------------------
  # x was divided by scale_var upstream in mfpi.default for numerical stability
  # during power selection. Multiplying by scale_var here restores the
  # scaled-but-not-shifted variable: cont_var_scaled * scale_var = x + shift.
  #
  # This function is always called from the flex functions with shift = 0 and
  # scale = 1 (no further shift or scale applied by transform_vector_fp below),
  # so the backscaled cont_var is passed directly to transform_vector_fp and
  # coefficients are on the phi(x + shift) scale, matching standalone mfp2.
  #
  # For zero_var variables (shift forced to 0 in mfpi.default), backscaling
  # restores the original x. The zero argument then correctly handles
  # exact-zero values in that restored x.
  #
  # Applied AFTER validation so checks run on the pre-backscaled (scaled) values.
  if (!is.null(scale_var) && any(scale_var != 1)) {
    cont_var <- cont_var * scale_var
  }

  # Structural-zero rows from the original continuous variable.
  zero_rows  <- if (zero) as.vector(cont_var) == 0 else rep(FALSE, n)
  valid_rows <- !zero_rows

  if (!any(valid_rows)) {
    stop("No valid rows are available for FP transformation/centering.",
         call. = FALSE)
  }

  # Precompute group-specific row indices once. These are reused when filling
  # the block-sparse interaction design and when resetting structural-zero rows.
  row_cache <- mfpi_group_row_cache(
    group_idx  = group_idx,
    n_groups   = n_groups,
    zero_rows  = zero_rows,
    valid_rows = valid_rows,
    zero       = zero
  )

  rows_by_group  <- row_cache$rows_by_group
  valid_by_group <- row_cache$valid_by_group
  zero_by_group  <- row_cache$zero_by_group

  # Transform full continuous variable first ----------------------------------
  # This avoids applying FP transformations to artificial out-of-group zeros.
  x_fp <- transform_vector_fp(
    x            = cont_var,
    power        = power,
    shift        = shift,
    scale        = scale,
    zero         = zero,
    check_binary = TRUE,
    name         = xname
  )

  x_fp <- as.matrix(x_fp)
  storage.mode(x_fp) <- "double"

  if (any(!is.finite(x_fp[valid_rows, , drop = FALSE]))) {
    stop(
      "FP transformation produced non-finite values among valid rows. ",
      "Check whether cont_var needs shifting/scaling.",
      call. = FALSE
    )
  }

  n_terms <- ncol(x_fp)

  z_names <- paste0(
    xname,
    rep(as.character(group_levels), each = n_terms),
    rep(seq_len(n_terms), times = n_groups)
  )

  column_groups <- stats::setNames(
    lapply(seq_len(n_groups), function(g) {
      z_names[seq.int((g - 1L) * n_terms + 1L, g * n_terms)]
    }),
    as.character(group_levels)
  )

  z <- matrix(0, nrow = n, ncol = n_groups * n_terms)
  colnames(z) <- z_names

  center_vals_work <- numeric(n_groups * n_terms)
  names(center_vals_work) <- z_names

  if (center && center_type == "grand") {
    grand_means <- colMeans(x_fp[valid_rows, , drop = FALSE], na.rm = TRUE)

    if (any(!is.finite(grand_means))) {
      stop("Could not compute finite grand centering constants.",
           call. = FALSE)
    }
  }

  # Fill group-specific fitting design matrix ---------------------------------
  # The fitting design is block-sparse: each group owns one FP block, while rows
  # from other groups remain structural zeros. We still assign group blocks one
  # at a time, but precomputed row indices avoid rebuilding logical masks inside
  # the loop.
  x_fp_shared_centered <- NULL

  if (!center) {
    # No centering: every group receives the same uncentered FP basis.
    x_fp_shared_centered <- x_fp
  } else if (center_type == "grand") {
    # Grand centering: every group receives the same centered FP basis.
    x_fp_shared_centered <- sweep(
      x_fp,
      2L,
      grand_means,
      "-",
      check.margin = FALSE
    )
  }

  for (g in seq_len(n_groups)) {
    rows_g_all <- rows_by_group[[g]]
    cols_g     <- seq.int((g - 1L) * n_terms + 1L, g * n_terms)

    if (center) {
      if (center_type == "grand") {
        centers_g <- grand_means
      } else {
        rows_g_valid <- valid_by_group[[g]]

        if (length(rows_g_valid) == 0L) {
          stop(
            paste0(
              "Cannot compute within-group centering constants for group ",
              group_levels[g], ": no valid rows."
            ),
            call. = FALSE
          )
        }

        centers_g <- colMeans(x_fp[rows_g_valid, , drop = FALSE], na.rm = TRUE)

        if (any(!is.finite(centers_g))) {
          stop(
            paste0(
              "Could not compute finite centering constants for group ",
              group_levels[g], "."
            ),
            call. = FALSE
          )
        }
      }
    } else {
      centers_g <- rep(0, n_terms)
    }

    center_vals_work[cols_g] <- centers_g

    # z_g(x_i) = [B_g(x_i) - center_g] * I(group_i == g).
    # Out-of-group rows remain zero because `z` was initialized as zero.
    if (!center || center_type == "grand") {
      z[rows_g_all, cols_g] <- x_fp_shared_centered[rows_g_all, , drop = FALSE]
    } else {
      z[rows_g_all, cols_g] <- sweep(
        x_fp[rows_g_all, , drop = FALSE],
        2L,
        centers_g,
        "-",
        check.margin = FALSE
      )
    }

    # Restore structural-zero rows after centering. This prevents zero rows
    # from becoming `-center_g` in the active group block.
    if (zero && length(zero_by_group[[g]]) > 0L) {
      z[zero_by_group[[g]], cols_g] <- 0
    }
  }

  # Pooled transformed variable for possible main-effect FP terms --------------
  xtransformed <- x_fp

  if (center) {
    pooled_centers <- colMeans(x_fp[valid_rows, , drop = FALSE], na.rm = TRUE)

    xtransformed <- sweep(
      x_fp,
      2L,
      pooled_centers,
      "-",
      check.margin = FALSE
    )

    if (zero && any(zero_rows)) {
      xtransformed[zero_rows, ] <- 0
    }
  }

  # Optional evaluation matrix -------------------------------------------------
  # Unlike `z`, this evaluates every group-specific function at every observed
  # x value. It does not contain out-of-group structural zeros.
  x_eval <- NULL

  if (return_eval) {
    x_eval <- matrix(0, nrow = n, ncol = n_groups * n_terms)
    colnames(x_eval) <- z_names

    if (!center || center_type == "grand") {
      # With no centering, every group subtracts zero. With grand centering,
      # every group subtracts the same constants. Therefore all group blocks are
      # identical and can be constructed once, then repeated.
      base_cols <- seq_len(n_terms)

      x_eval_base <- sweep(
        x_fp,
        2L,
        center_vals_work[base_cols],
        "-",
        check.margin = FALSE
      )

      x_eval[, ] <- x_eval_base[, rep(seq_len(n_terms), times = n_groups),
                                drop = FALSE]
    } else {
      # Group centering uses different constants by group, so each group block
      # must still be centered separately.
      for (g in seq_len(n_groups)) {
        cols_g <- seq.int((g - 1L) * n_terms + 1L, g * n_terms)

        x_eval[, cols_g] <- sweep(
          x_fp,
          2L,
          center_vals_work[cols_g],
          "-",
          check.margin = FALSE
        )
      }
    }

    # Structural-zero rows are zero in every evaluated group-specific function.
    # Reset once across all columns instead of once per group block.
    if (zero && any(zero_rows)) {
      x_eval[zero_rows, ] <- 0
    }
  }

  # Restore structural-zero rows explicitly, then fail on any remaining
  # non-finite values. Non-finite values outside structural-zero rows indicate an
  # internal transformation or centering problem and must not be silently repaired.
  if (zero && any(zero_rows)) {
    xtransformed[zero_rows, ] <- 0

    if (!is.null(x_eval)) {
      x_eval[zero_rows, ] <- 0
    }
  }

  if (any(!is.finite(z))) {
    stop(
      "Non-finite values remain in the group-specific interaction design.",
      call. = FALSE
    )
  }

  if (any(!is.finite(xtransformed))) {
    stop(
      "Non-finite values remain in the pooled transformed continuous variable.",
      call. = FALSE
    )
  }

  if (!is.null(x_eval) && any(!is.finite(x_eval))) {
    stop(
      "Non-finite values remain in the fitted-function evaluation matrix.",
      call. = FALSE
    )
  }

  center_vals <- if (center) center_vals_work else NULL

  group_fp_powers <- stats::setNames(
    rep(list(power), n_groups),
    as.character(group_levels)
  )

  attr(z, "fp_centers")      <- center_vals
  attr(z, "center")          <- center
  attr(z, "center_type")     <- if (center) center_type else "none"
  attr(z, "group_levels")    <- group_levels
  attr(z, "group_fp_powers") <- group_fp_powers
  attr(z, "power")           <- power
  attr(z, "column_groups") <- column_groups

  if (!is.null(x_eval)) {
    attr(x_eval, "fp_centers")      <- center_vals
    attr(x_eval, "center")          <- center
    attr(x_eval, "center_type")     <- if (center) center_type else "none"
    attr(x_eval, "group_levels")    <- group_levels
    attr(x_eval, "group_fp_powers") <- group_fp_powers
    attr(x_eval, "power")           <- power
    attr(x_eval, "column_groups")   <- column_groups
  }

  out <- list(
    z              = z,
    xtransformed   = xtransformed,
    x_eval         = x_eval,
    center_vals    = center_vals,
    x_fp           = x_fp,
    column_groups  = column_groups
  )

  return(out)
}

# Utility functions for MFPI
#
# These functions handle dummy variable creation, group-specific FP variable
# construction, and related data-preparation tasks used by MFPI. Candidate FP
# enumeration for flex2 is handled directly by the compact-basis workflow in
# `mfpi_flexibilities.R`, so this file contains only the reusable construction
# helpers required by the package.
#
# Naming conventions (shared with flex_functions.R and fit_mfpi.R):
#   group_var  - the categorical grouping variable (column name or matrix)
#   cont_var   - the continuous variable being transformed
#   fp_cand    - candidate FP powers for a single variable
#   fp_degree  - degree of the FP (1 = FP1, 2 = FP2)
#   na_replace - how to handle non-finite values from structural zeros


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
#' @param quiet Logical. If \code{TRUE}, suppress informational messages about
#'   levels supplied in \code{levels} that are not present in \code{x}. This is
#'   used during prediction, where new data may validly contain only a subset of
#'   the fitted group levels. Unseen levels in \code{x} that are not listed in
#'   \code{levels} still produce an error.
#'
#' @return A numeric matrix with `nrow(x)` rows and \eqn{K - 1} columns, where
#'   \eqn{K} is the number of unique levels (from `x` or from `levels`). Column
#'   names follow the pattern `<varname><level>`, e.g. `"trt2"`, `"trt3"` when
#'   the input column is named `"trt"` and the non-reference levels are 2 and 3.
#'
#' @examples
#' \dontrun{
#' x <- matrix(c(3, 1, 2, 1, 3, 2), ncol = 1)
#' colnames(x) <- "group"
#' create_group_dummies(x)
#'
#' # Enforce a fixed set of levels (e.g. for prediction on a subset)
#' x2 <- matrix(c(1, 1, 1), ncol = 1)
#' colnames(x2) <- "trt"
#' create_group_dummies(x2, levels = 1:3)   # trt2 and trt3 columns are all-zero
#'}
#' @keywords internal
#' @noRd
create_group_dummies <- function(x, levels = NULL, quiet = FALSE) {

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

    if (length(extra_in_data) > 0L) {
      stop(
        paste0(
          "Values in `x` not listed in `levels`: ",
          paste(extra_in_data, collapse = ", "),
          "."
        ),
        call. = FALSE
      )
    }

    if (length(missing_in_data) > 0L && !isTRUE(quiet)) {
      message(
        paste0(
          "Levels not present in `x` (will produce all-zero columns): ",
          paste(missing_in_data, collapse = ", "),
          "."
        )
      )
    }
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
# var_group() -----------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Partition Generated FP Variable Names into Group-Specific Subsets
#'
#' Given a continuous-variable name, a set of column names, and the known
#' internal group codes used by MFPI, returns a list where each element contains
#' the generated FP variable names belonging to one group.
#'
#' This helper does not infer group membership from a regular expression.
#' Instead, it constructs the expected generated names directly from known model
#' metadata:
#'
#' \preformatted{
#'   <var_prefix><group_code><power_index>
#' }
#'
#' For example, if \code{var_prefix = "age"}, \code{group_codes = c(0, 1, 10)},
#' and \code{power_indices = 1:2}, the expected generated names are:
#'
#' \preformatted{
#'   age01, age02,
#'   age11, age12,
#'   age101, age102
#' }
#'
#' This design avoids fragile string-pattern matching. In particular, valid
#' variable names containing characters such as \code{.}, \code{+}, \code{(},
#' or \code{)} are treated literally because the function compares complete
#' expected names with \code{var_names}; it never interprets variable names as
#' regular expressions.
#'
#' The final digit in each generated name is treated as the FP power-column
#' index. The preceding suffix characters after \code{var_prefix} identify the
#' internal group code. Multi-digit group codes such as \code{"10"} and
#' \code{"11"} are therefore handled explicitly through \code{group_codes}.
#'
#' This is a lookup helper used after \code{create_z_variables()} to recover
#' which generated columns belong to each level of \code{group_var}.
#'
#' @param var_prefix Character scalar. The original continuous-variable name
#'   used as the prefix of generated FP variables.
#' @param var_names Character vector of column names to search.
#' @param group_codes Vector of known internal group codes. These are usually
#'   \code{0}, \code{1}, ..., or their character equivalents. They should be the
#'   same internal codes used when the group-specific FP columns were generated.
#' @param power_indices Vector of allowed FP power-column indices. The default
#'   \code{1:9} preserves the previous helper's ability to recognise any
#'   single-digit power index while still matching names exactly.
#'
#' @return A list of character vectors, one per group code with at least one
#'   generated column present in \code{var_names}. Groups are returned in the
#'   order supplied by \code{group_codes}. Within each group, column names are
#'   returned in the order in which they appear in \code{var_names}.
#'
#' @examples
#' var_group(
#'   var_prefix = "age",
#'   var_names = c("age01", "age02", "age11", "age12", "age101", "age102", "bmi"),
#'   group_codes = c(0, 1, 10),
#'   power_indices = 1:2
#' )
#' # Returns:
#' # list(
#' #   c("age01", "age02"),
#' #   c("age11", "age12"),
#' #   c("age101", "age102")
#' # )
#'
#' @keywords internal
#' @noRd
var_group <- function(var_prefix,
                      var_names,
                      group_codes,
                      power_indices = 1:9) {
  if (!is.character(var_prefix) || length(var_prefix) != 1L ||
      is.na(var_prefix) || !nzchar(var_prefix)) {
    stop(
      "`var_prefix` must be a single non-empty character string.",
      call. = FALSE
    )
  }

  if (!is.character(var_names) || length(var_names) == 0L ||
      anyNA(var_names)) {
    stop(
      "`var_names` must be a non-empty character vector without missing values.",
      call. = FALSE
    )
  }

  if (missing(group_codes) || length(group_codes) == 0L ||
      anyNA(group_codes) || !is.atomic(group_codes)) {
    stop(
      "`group_codes` must be a non-empty atomic vector without missing values.",
      call. = FALSE
    )
  }

  if (length(power_indices) == 0L || anyNA(power_indices) ||
      !is.atomic(power_indices)) {
    stop(
      "`power_indices` must be a non-empty atomic vector without missing values.",
      call. = FALSE
    )
  }

  group_codes <- unique(as.character(group_codes))
  power_indices <- unique(as.character(power_indices))

  # The naming convention stores the FP power-column index as the final single
  # digit. Restricting `power_indices` to one digit avoids ambiguous names such
  # as paste0("age", "1", "10"), where the boundary between group code and
  # power index would no longer be clear.
  if (any(nchar(power_indices) != 1L) ||
      any(!(power_indices %in% as.character(0:9)))) {
    stop(
      "`power_indices` must contain only single-digit values from 0 to 9.",
      call. = FALSE
    )
  }

  groups <- lapply(group_codes, function(g) {
    # Construct the exact names that are valid for this variable and group.
    # This replaces the old regex-based search:
    #
    #   grep(paste0("^", var_prefix, "\\d{2,}$"), var_names, value = TRUE)
    #
    # Exact matching avoids accidental matches when `var_prefix` contains regex
    # metacharacters such as ".", "+", "(", or ")".
    expected <- paste0(var_prefix, g, power_indices)

    # Preserve the existing column order from `var_names`, because downstream
    # model-matrix code may rely on the order in which columns appear.
    var_names[var_names %in% expected]
  })

  groups[lengths(groups) > 0L]
}
