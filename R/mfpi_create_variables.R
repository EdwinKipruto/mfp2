# -----------------------------------------------------------------------------
# create_z_variables() ---------------------------------------------------------
# -----------------------------------------------------------------------------

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
#' shifted before calling this function unless \code{zero = TRUE} is intended.
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
#' When \code{zero = TRUE}, non-positive values of \code{cont_var} are excluded
#' from the centering means and reset to zero after centering.
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
#'   covariate. It must have a column name and contain only finite values.
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
#' @param zero Logical. If \code{TRUE}, non-positive values of \code{cont_var}
#'   are treated as structural zeros: they are excluded from centering means and
#'   reset to zero after centering. Default is \code{FALSE}.
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
#'   \code{zero} argument then correctly handles non-positive values.
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
#' @seealso \code{\link{transform_z_variables}},
#'   \code{transform_vector_fp}
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
    warning("cont_var has zero variance; all values are identical.",
            call. = FALSE)
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
  # non-positive values in that restored x.
  #
  # Applied AFTER validation so checks run on the pre-backscaled (scaled) values.
  if (!is.null(scale_var) && any(scale_var != 1)) {
    cont_var <- cont_var * scale_var
  }
  
  # Structural-zero rows from the original continuous variable.
  zero_rows  <- if (zero) as.vector(cont_var) <= 0 else rep(FALSE, n)
  valid_rows <- !zero_rows
  
  if (!any(valid_rows)) {
    stop("No valid rows are available for FP transformation/centering.",
         call. = FALSE)
  }
  
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
  for (g in seq_len(n_groups)) {
    in_group <- group_idx == g
    cols_g   <- seq.int((g - 1L) * n_terms + 1L, g * n_terms)
    
    if (center) {
      if (center_type == "grand") {
        centers_g <- grand_means
      } else {
        rows_g <- in_group & valid_rows
        
        if (!any(rows_g)) {
          stop(
            paste0(
              "Cannot compute within-group centering constants for group ",
              group_levels[g], ": no valid rows."
            ),
            call. = FALSE
          )
        }
        
        centers_g <- colMeans(x_fp[rows_g, , drop = FALSE], na.rm = TRUE)
        
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
    
    # z_g(x_i) = [B_g(x_i) - center_g] * I(group_i == g)
    # Out-of-group rows remain structural zeros.
    z[in_group, cols_g] <- sweep(
      x_fp[in_group, , drop = FALSE],
      2L,
      centers_g,
      "-",
      check.margin = FALSE
    )
    
    if (zero && any(in_group & zero_rows)) {
      z[in_group & zero_rows, cols_g] <- 0
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
  # Unlike z, this evaluates every group-specific function at every x value.
  x_eval <- NULL
  
  if (return_eval) {
    x_eval <- matrix(0, nrow = n, ncol = n_groups * n_terms)
    colnames(x_eval) <- z_names
    
    for (g in seq_len(n_groups)) {
      cols_g <- seq.int((g - 1L) * n_terms + 1L, g * n_terms)
      
      x_eval[, cols_g] <- sweep(
        x_fp,
        2L,
        center_vals_work[cols_g],
        "-",
        check.margin = FALSE
      )
      
      if (zero && any(zero_rows)) {
        x_eval[zero_rows, cols_g] <- 0
      }
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
}# -----------------------------------------------------------------------------
# transform_z_variables() ------------------------------------------------------
# -----------------------------------------------------------------------------

#' Enumerate Group-Specific Fractional-Polynomial Interaction Designs
#'
#' Generates group-specific fractional-polynomial interaction design matrices
#' for all candidate FP power combinations. This function is used when searching
#' over possible FP1 or FP2 interaction forms.
#'
#' For each candidate power combination, the full continuous variable is first
#' FP-transformed. The transformed values are then copied into group-specific
#' design blocks. Rows outside each group remain structural zeros.
#'
#' This implementation deliberately avoids transforming zero-padded
#' group-specific variables. Instead, it transforms the full continuous variable
#' first. This is safer for log and negative powers because artificial
#' out-of-group zeros are not treated as real covariate values.
#'
#' @section Mathematical definition:
#'
#' Let \eqn{x_i} be the continuous covariate and \eqn{G_i} the group membership.
#' For a candidate FP power vector \eqn{p}, let \eqn{B_p(x_i)} denote the FP
#' basis generated from \eqn{x_i}.
#'
#' For group \eqn{g}, the fitting design block is
#'
#' \deqn{
#' Z_{g,p}(x_i) = \{B_p(x_i) - c_{g,p}\} I(G_i = g),
#' }
#'
#' where \eqn{c_{g,p}} is the centering vector for group \eqn{g} under power
#' combination \eqn{p}. The full design matrix for that power combination is
#'
#' \deqn{
#' Z_p = [Z_{g_1,p}, Z_{g_2,p}, \ldots, Z_{g_K,p}].
#' }
#'
#' @section Candidate FP powers:
#'
#' Candidate powers are generated by \code{mfp2::generate_powers_fp()} using
#' \code{fp_cand} and \code{fp_degree}.
#'
#' For FP1, each candidate power vector has length one. For FP2, each candidate
#' power vector has length two. Repeated powers are interpreted using the
#' standard fractional-polynomial rule:
#'
#' \describe{
#'   \item{Non-repeated powers}{For powers \eqn{p_1 \ne p_2}, the basis is
#'   \eqn{(x^{p_1}, x^{p_2})}, with \eqn{p = 0} interpreted as \eqn{\log(x)}.}
#'   \item{Repeated nonzero powers}{For powers \eqn{p_1 = p_2 = p \ne 0}, the
#'   basis is \eqn{(x^p, x^p \log(x))}.}
#'   \item{Repeated zero powers}{For powers \eqn{p_1 = p_2 = 0}, the basis is
#'   \eqn{(\log(x), \log(x)^2)}.}
#' }
#'
#' @section Centering:
#'
#' If \code{center = FALSE}, no centering is applied and
#' \code{center_vals_list} is returned as \code{NULL}.
#'
#' If \code{center = TRUE}, centering constants are computed separately for each
#' candidate power combination and stored in \code{center_vals_list}. The
#' centering strategy is controlled by \code{center_type}:
#'
#' \describe{
#'   \item{\code{"grand"}}{
#'   For each candidate FP basis, the empirical column mean over all valid
#'   observations is subtracted. The same centering vector is used for every
#'   group:
#'
#'   \deqn{
#'   c_{g,p} = n^{-1}\sum_i B_p(x_i).
#'   }
#'   }
#'
#'   \item{\code{"group"}}{
#'   For each group, the empirical column mean within that group is subtracted:
#'
#'   \deqn{
#'   c_{g,p} = n_g^{-1}\sum_{i:G_i=g} B_p(x_i).
#'   }
#'   }
#' }
#'
#' When \code{zero = TRUE}, non-positive values of \code{cont_var} are excluded
#' from centering means and reset to zero after centering.
#'
#' @section Structural zeros:
#'
#' There are two distinct notions of zero:
#'
#' \describe{
#'   \item{Out-of-group structural zeros}{Rows outside a group are set to zero
#'   in that group's fitting-design block. These zeros are not transformed.
#'   They are created only after the full continuous variable has been
#'   FP-transformed.}
#'   \item{Non-positive covariate values}{If \code{zero = TRUE}, non-positive
#'   values of \code{cont_var} are treated as structural zeros of the original
#'   covariate. They are excluded from centering and reset to zero after
#'   transformation and centering.}
#' }
#'
#' This separation prevents invalid operations such as \eqn{\log(0)} or
#' \eqn{0^p} for negative powers on artificial out-of-group zeros.
#'
#' @section Fitting matrix versus evaluation matrix:
#'
#' The matrices in \code{z_transformed} are fitting design matrices. For each
#' group \eqn{g}, rows outside that group remain zero:
#'
#' \deqn{
#' Z_{g,p}(x_i) = 0 \quad \text{if } G_i \ne g.
#' }
#'
#' If \code{return_eval = TRUE}, the function also returns
#' \code{x_eval_transformed}. These matrices evaluate every group-specific
#' function for every observed value of \code{cont_var}:
#'
#' \deqn{
#' X_{\mathrm{eval},g,p}(x_i) = B_p(x_i) - c_{g,p}.
#' }
#'
#' Evaluation matrices are useful for computing fitted functions, fitted
#' function differences, standard errors, and plotting curves over a common
#' \eqn{x} grid.
#'
#' @section Column naming:
#'
#' Each generated matrix uses the convention:
#'
#' \deqn{
#' \texttt{<varname><group><power_index>}.
#' }
#'
#' For example, if the continuous variable is named \code{"cavol"}, group levels
#' are \code{0} and \code{1}, and the FP degree is 2, columns are named:
#'
#' \preformatted{
#' cavol01 cavol02 cavol11 cavol12
#' }
#'
#' @section Stored metadata:
#'
#' Each matrix in \code{z_transformed}, and each matrix in
#' \code{x_eval_transformed} when requested, contains attributes:
#'
#' \describe{
#'   \item{\code{fp_centers}}{Centering constants for the corresponding power
#'   combination.}
#'   \item{\code{center}}{Whether centering was applied.}
#'   \item{\code{center_type}}{The centering strategy.}
#'   \item{\code{group_levels}}{The observed group levels.}
#'   \item{\code{group_fp_powers}}{A named list of powers, one element per
#'   group.}
#'   \item{\code{power}}{The candidate power vector used to generate that
#'   matrix.}
#' }
#'
#' These metadata are intended for downstream model fitting, fitted-function
#' generation, and prediction. In particular, centering constants should be
#' reused and not recomputed on prediction data.
#'
#' @param cont_var A one-column numeric matrix containing the continuous
#'   covariate. It must have a column name and contain only finite values.
#' @param group_var A one-column numeric matrix containing group membership. It
#'   must have the same number of rows as \code{cont_var}, contain no missing
#'   values, and have at least two distinct values.
#' @param shift Optional numeric shift applied before FP transformation. Passed
#'   to \code{transform_vector_fp()}.
#' @param scale Optional numeric scale applied before FP transformation. Passed
#'   to \code{mfp2::transform_vector_fp()}.
#' @param center Logical. If \code{TRUE}, center each candidate FP basis before
#'   constructing the group-specific design matrices. Default is \code{FALSE}.
#' @param acdx Logical. Whether to apply the ACD transformation. Currently,
#'   \code{TRUE} is not supported in this structural-zero-safe implementation
#'   and will produce an error.
#' @param fp_cand Numeric vector of candidate FP powers. If \code{NULL}, the
#'   standard Royston-Altman set is used:
#'   \code{c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)}.
#' @param fp_degree Positive integer specifying the FP degree. Use \code{1L} for
#'   FP1 and \code{2L} for FP2.
#' @param center_type Character string controlling centering when
#'   \code{center = TRUE}. Either \code{"grand"} or \code{"group"}. Default is
#'   \code{"grand"}.
#' @param zero Logical. If \code{TRUE}, non-positive values of \code{cont_var}
#'   are treated as structural zeros, excluded from centering means, and reset to
#'   zero after centering.
#' @param return_eval Logical. If \code{TRUE}, also return evaluation matrices
#'   in \code{x_eval_transformed}. Default is \code{FALSE}.
#' @param scale_var Numeric scalar. The scale factor for \code{cont_var},
#'   as computed in \code{mfpi.default()}. Before any FP transformation,
#'   \code{cont_var} is multiplied by \code{scale_var} to restore the
#'   variable before scaling was applied. For standard variables this gives
#'   \eqn{x + \text{shift}}; for \code{zero_var = TRUE} variables (where
#'   shift is forced to 0) this gives the original \eqn{x}, and the
#'   \code{zero} argument then correctly handles non-positive values.
#'   This ensures model coefficients match the \eqn{\phi(x + \text{shift})}
#'   scale of the adjustment model and standalone \pkg{mfp2}. Default
#'   \code{1} (no backscaling).
#'
#' @return A list with the following elements:
#'
#' \describe{
#'   \item{\code{z_transformed}}{Named list of group-specific fitting design
#'   matrices, one per candidate FP power combination.}
#'   \item{\code{z_untransformed}}{Group-specific linear design matrix generated
#'   with \code{power = 1}. This is retained for compatibility and diagnostics.}
#'   \item{\code{x_eval_transformed}}{Named list of group-specific evaluation
#'   matrices if \code{return_eval = TRUE}; otherwise \code{NULL}.}
#'   \item{\code{powers_matrix}}{Numeric matrix of candidate FP power
#'   combinations. Each row corresponds to one element of
#'   \code{z_transformed}.}
#'   \item{\code{znames}}{Column names of \code{z_untransformed}.}
#'   \item{\code{center_vals_list}}{Named list of centering-constant vectors,
#'   one per candidate FP power combination, or \code{NULL} when
#'   \code{center = FALSE}.}
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
#' transform_z_variables(
#'   cont_var  = x_raw,
#'   group_var = grp,
#'   fp_cand   = c(-2, -1, -0.5, 0, 0.5, 1, 2, 3),
#'   fp_degree = 2,
#'   shift     = NULL,   # estimated automatically
#'   scale     = NULL,   # estimated automatically
#'   scale_var = 1       # default: no backscaling
#' )}}
#'   \item{Pre-shifted and scaled data (internal mfpi use)}{When called
#'     from flex2 inside \code{mfpi()}, \code{cont_var} has already been
#'     shifted and scaled:
#'     \eqn{x_{\text{scaled}} = (x + \text{shift}) / \text{scale\_var}}.
#'     Pass \code{shift = 0}, \code{scale = 1}, and the real
#'     \code{scale_var} so that the function backscales before FP
#'     transformation and coefficients are on the
#'     \eqn{\phi(x + \text{shift})} scale:
#'     \preformatted{
#' transform_z_variables(
#'   cont_var  = x_scaled,   # (x + shift) / scale_var
#'   group_var = grp,
#'   fp_cand   = c(-2, -1, -0.5, 0, 0.5, 1, 2, 3),
#'   fp_degree = 2,
#'   shift     = 0,
#'   scale     = 1,
#'   scale_var = scale_var   # real scale factor: x_scaled * scale_var = x + shift
#' )}}
#' \item{\code{column_groups_list}}{Named list, parallel to
#'   \code{z_transformed}, containing group-to-column mappings for each
#'   candidate FP power combination.}
#' }
#'
#' @seealso \code{\link{create_z_variables}},
#'   \code{mfp2::generate_powers_fp},
#'   \code{mfp2::transform_vector_fp}
#'
#' @examples
#' cont <- matrix(1:6, ncol = 1)
#' colnames(cont) <- "age"
#'
#' grp <- matrix(c(0, 1, 0, 1, 0, 1), ncol = 1)
#' colnames(grp) <- "trt"
#'
#' # Enumerate FP2 designs using two candidate powers
#' transform_z_variables(
#'   cont_var = cont,
#'   group_var = grp,
#'   fp_cand = c(0, 1),
#'   fp_degree = 2,
#'   center = TRUE,
#'   center_type = "grand"
#' )
#'
#' # Also return evaluation matrices for fitted-function computation
#' transform_z_variables(
#'   cont_var = cont,
#'   group_var = grp,
#'   fp_cand = c(0, 1),
#'   fp_degree = 2,
#'   center = TRUE,
#'   return_eval = TRUE
#' )
#'
#' @importFrom stats setNames
#' @importFrom mfp2 generate_powers_fp transform_vector_fp
#' @keywords internal
#' @noRd
transform_z_variables <- function(cont_var,
                                  group_var,
                                  shift = NULL,
                                  scale = NULL,
                                  center = FALSE,
                                  acdx = FALSE,
                                  fp_cand = NULL,
                                  fp_degree = 2L,
                                  center_type = c("grand", "group"),
                                  zero = FALSE,
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
    stop("cont_var must not contain NA, NaN, or infinite values.",
         call. = FALSE)
  }
  
  if (!all(is.finite(group_var))) {
    stop("group_var must not contain NA, NaN, or infinite values.", call. = FALSE)
  }
  
  group_vec    <- as.vector(group_var)
  group_levels <- sort(unique(group_vec))
  group_idx    <- match(group_vec, group_levels)
  n_groups     <- length(group_levels)
  
  if (n_groups < 2L) {
    stop("group_var must have at least two distinct values.", call. = FALSE)
  }
  
  if (isTRUE(any(acdx))) {
    stop(
      "acdx = TRUE is not supported in this structural-zero-safe implementation.",
      call. = FALSE
    )
  }
  
  # Backscale cont_var to match mfp2 convention --------------------------------
  # x was divided by scale_var upstream in mfpi.default for numerical stability
  # during power selection. Multiplying by scale_var here restores the
  # shifted-but-not-scaled variable: cont_var_scaled * scale_var = x + shift.
  #
  # This function is always called from flex2 with shift = 0 and scale = 1
  # (no further shift or scale applied by transform_vector_fp below), so the
  # backscaled cont_var is passed directly to transform_vector_fp and
  # coefficients are on the phi(x + shift) scale, matching standalone mfp2.
  #
  # For zero_var variables (shift forced to 0 in mfpi.default), backscaling
  # restores the original x. The zero argument then correctly handles
  # non-positive values in that restored x.
  #
  # Applied AFTER validation so checks run on the pre-backscaled (scaled) values.
  if (!is.null(scale_var) && any(scale_var != 1)) {
    cont_var <- cont_var * scale_var
  }
  
  if (is.null(fp_cand)) {
    fp_cand <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
  }
  
  if (!is.numeric(fp_cand) || length(fp_cand) == 0L ||
      anyNA(fp_cand) || any(!is.finite(fp_cand))) {
    stop(
      "fp_cand must be a non-empty finite numeric vector with no missing values.",
      call. = FALSE
    )
  }
  
  fp_degree <- as.integer(fp_degree)
  
  if (length(fp_degree) != 1L || fp_degree < 1L) {
    stop("fp_degree must be a positive integer.", call. = FALSE)
  }
  
  n <- nrow(cont_var)
  
  zero_rows  <- if (zero) as.vector(cont_var) <= 0 else rep(FALSE, n)
  valid_rows <- !zero_rows
  
  if (!any(valid_rows)) {
    stop("No valid rows are available for FP transformation/centering.",
         call. = FALSE)
  }
  
  # Candidate powers -----------------------------------------------------------
  powers_matrix <- generate_powers_fp(
    degree = fp_degree,
    powers = fp_cand
  )
  
  powers_matrix <- as.matrix(powers_matrix)
  storage.mode(powers_matrix) <- "double"
  
  n_combinations <- nrow(powers_matrix)
  
  z_transformed      <- vector("list", n_combinations)
  center_vals_list   <- if (center) vector("list", n_combinations) else NULL
  x_eval_transformed <- if (return_eval) vector("list", n_combinations) else NULL
  column_groups_list  <- vector("list", n_combinations)
  
  power_names <- apply(
    powers_matrix,
    1L,
    function(p) paste(p, collapse = ",")
  )
  
  # Build untransformed linear design for compatibility -----------------------
  z_untransformed <- create_z_variables(
    cont_var    = cont_var,
    group_var   = group_var,
    power       = 1,
    shift       = shift,
    scale       = scale,
    center      = FALSE,
    zero        = zero,
    center_type = center_type,
    return_eval = FALSE
  )$z
  
  znames <- colnames(z_untransformed)
  
  # Enumerate FP transformations ----------------------------------------------
  for (i in seq_len(n_combinations)) {
    
    power_i <- as.numeric(powers_matrix[i, ])
    
    # Transform full continuous variable first.
    # Do not transform zero-padded group variables.
    x_fp <- transform_vector_fp(
      x            = cont_var,
      power        = power_i,
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
        cols_g <- seq.int((g - 1L) * n_terms + 1L, g * n_terms)
        z_names[cols_g]
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
    
    for (g in seq_len(n_groups)) {
      in_group <- group_idx == g
      cols_g   <- seq.int((g - 1L) * n_terms + 1L, g * n_terms)
      
      if (center) {
        if (center_type == "grand") {
          centers_g <- grand_means
        } else {
          rows_g <- in_group & valid_rows
          
          if (!any(rows_g)) {
            stop(
              paste0(
                "Cannot compute within-group centering constants for group ",
                group_levels[g], ": no valid rows."
              ),
              call. = FALSE
            )
          }
          
          centers_g <- colMeans(x_fp[rows_g, , drop = FALSE], na.rm = TRUE)
          
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
      
      # Fitting design:
      # active rows get centered FP values; inactive rows stay zero.
      z[in_group, cols_g] <- sweep(
        x_fp[in_group, , drop = FALSE],
        2L,
        centers_g,
        "-",
        check.margin = FALSE
      )
      
      if (zero && any(in_group & zero_rows)) {
        z[in_group & zero_rows, cols_g] <- 0
      }
    }
    
    # Optional evaluation matrix ----------------------------------------------
    x_eval <- NULL
    
    if (return_eval) {
      x_eval <- matrix(0, nrow = n, ncol = n_groups * n_terms)
      colnames(x_eval) <- z_names
      
      for (g in seq_len(n_groups)) {
        cols_g <- seq.int((g - 1L) * n_terms + 1L, g * n_terms)
        
        x_eval[, cols_g] <- sweep(
          x_fp,
          2L,
          center_vals_work[cols_g],
          "-",
          check.margin = FALSE
        )
        
        if (zero && any(zero_rows)) {
          x_eval[zero_rows, cols_g] <- 0
        }
      }
    }
    
    # Restore structural-zero rows explicitly. Inactive group blocks were
    # initialized to zero; active structural-zero rows should also carry zero FP
    # contribution. Non-finite values elsewhere are errors.
    if (zero && any(zero_rows)) {
      for (g in seq_len(n_groups)) {
        cols_g <- seq.int((g - 1L) * n_terms + 1L, g * n_terms)
        z[zero_rows, cols_g] <- 0
      }
      
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
    
    if (!is.null(x_eval) && any(!is.finite(x_eval))) {
      stop(
        "Non-finite values remain in the fitted-function evaluation matrix.",
        call. = FALSE
      )
    }
    
    center_vals <- if (center) center_vals_work else NULL
    
    group_fp_powers <- stats::setNames(
      rep(list(power_i), n_groups),
      as.character(group_levels)
    )
    
    attr(z, "fp_centers")      <- center_vals
    attr(z, "center")          <- center
    attr(z, "center_type")     <- if (center) center_type else "none"
    attr(z, "group_levels")    <- group_levels
    attr(z, "group_fp_powers") <- group_fp_powers
    attr(z, "power")           <- power_i
    attr(z, "column_groups") <- column_groups
    
    if (!is.null(x_eval)) {
      attr(x_eval, "fp_centers")      <- center_vals
      attr(x_eval, "center")          <- center
      attr(x_eval, "center_type")     <- if (center) center_type else "none"
      attr(x_eval, "group_levels")    <- group_levels
      attr(x_eval, "group_fp_powers") <- group_fp_powers
      attr(x_eval, "power")           <- power_i
      attr(x_eval, "column_groups") <- column_groups
    }
    
    z_transformed[[i]] <- z
    column_groups_list[[i]] <- column_groups
    
    if (center) {
      center_vals_list[[i]] <- center_vals
    }
    
    if (return_eval) {
      x_eval_transformed[[i]] <- x_eval
    }
  }
  
  names(z_transformed) <- power_names
  names(column_groups_list) <- power_names
  if (center) {
    names(center_vals_list) <- power_names
  }
  
  if (return_eval) {
    names(x_eval_transformed) <- power_names
  }
  
  out <- list(
    z_transformed      = z_transformed,
    z_untransformed    = z_untransformed,
    x_eval_transformed = x_eval_transformed,
    powers_matrix      = powers_matrix,
    znames             = znames,
    center_vals_list   = center_vals_list,
    column_groups_list = column_groups_list
  )
  
  return(out)
}

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
#' Given a common variable prefix and a character vector of column names, returns
#' a list where each element contains the generated variable names belonging to
#' one group.
#'
#' Variable names are assumed to follow the convention
#' `<var_prefix><group_code><power_index>`, where `group_code` identifies the
#' level of `group_var` and `power_index` is the final single digit identifying
#' the fractional-polynomial power column. For example, in `"age105"`, the
#' prefix is `"age"`, the group code is `"10"`, and the power index is `"5"`.
#'
#' This parsing rule supports multi-digit group codes, so groups such as `"1"`,
#' `"10"`, and `"11"` are kept separate.
#'
#' This is a lookup helper used after \code{create_z_variables()} to recover
#' which generated columns belong to each level of `group_var`.
#'
#' @param var_prefix Character string. The common prefix of the generated
#'   variable names to partition. For example, `"age"` matches names such as
#'   `"age01"`, `"age02"`, `"age101"`, and `"age105"`.
#' @param var_names Character vector of column names to search.
#'
#' @return A list of character vectors, one per unique group code found between
#'   `var_prefix` and the final single-digit power index. Groups are returned in
#'   the order in which their first matching column appears in `var_names`.
#'
#' @examples
#' var_group(
#'   "age",
#'   c("age01", "age02", "age11", "age12", "age101", "age102", "bmi")
#' )
#' # Returns:
#' # list(
#' #   c("age01", "age02"),
#' #   c("age11", "age12"),
#' #   c("age101", "age102")
#' # )
#'
#' @examples
#' # FP5 example with groups 0, 1, and 10
#' var_group(
#'   "age",
#'   c(
#'     "age01", "age02", "age03", "age04", "age05",
#'     "age11", "age12", "age13", "age14", "age15",
#'     "age101", "age102", "age103", "age104", "age105"
#'   )
#' )
#'
#' @keywords internal
#' @noRd
var_group <- function(var_prefix, var_names) {
  if (!is.character(var_prefix) || length(var_prefix) != 1L)
    stop("`var_prefix` must be a single character string.", call. = FALSE)
  
  if (!is.character(var_names) || length(var_names) == 0L)
    stop("`var_names` must be a non-empty character vector.", call. = FALSE)
  
  matched <- grep(paste0("^", var_prefix, "\\d{2,}$"), var_names, value = TRUE)
  if (length(matched) == 0L) return(list())
  
  prefix_len <- nchar(var_prefix)
  
  suffix <- substr(matched, prefix_len + 1L, nchar(matched))
  
  # Naming convention: <varname><group><power_index>
  # The final digit is the power index; everything before it is the group code.
  group_code <- substr(suffix, 1L, nchar(suffix) - 1L)
  
  unique_grps <- unique(group_code)
  
  lapply(unique_grps, function(g) {
    matched[group_code == g]
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
#' @keywords internal
#' @noRd
adjust_reference_category <- function(x, use_top_as_ref = FALSE) {
  if (!is.factor(x)) return(x)
  
  levels(x) <- rev(levels(x))
  ref <- if (use_top_as_ref) tail(levels(x), 1L) else levels(x)[1L]
  relevel(x, ref = ref)
}