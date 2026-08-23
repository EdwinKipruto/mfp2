# -----------------------------------------------------------------------------
# winsorize_cont_vars() -------------------------------------------------------
# -----------------------------------------------------------------------------

#' Winsorise Continuous Variables to Reduce the Influence of Extreme
#' Observations
#'
#' Applies symmetric Winsorisation to one or more continuous columns of a
#' predictor matrix. For each variable, values below the lower percentile
#' cutoff are replaced by the lower cutoff value, and values above the upper
#' percentile cutoff are replaced by the upper cutoff value. The number of
#' distinct values is preserved (that is, no rounding is performed); only
#' extreme observations are truncated. Influential observations can 
#' disproportionately determine the selected FP functional form and may produce 
#' spuriously significant interactions. Winsorisation provides a transparent 
#' preprocessing step that reduces the influence of extreme values without 
#' removing observations from the analysis.
#'
#' @section Default percentiles:
#' By default, Winsorisation is applied at the 1st and 99th percentiles,
#' truncating approximately the lower and upper 1 percent of observations for each
#' variable.
#'
#' Variables listed in \code{zero_vars} are \strong{excluded entirely} from
#' Winsorisation because the user has indicated that both the spike at zero
#' and the positive component of the distribution carry substantive modelling
#' meaning. For transparency, their corresponding entries in
#' \code{limits} are recorded as \code{NA}.
#'
#' Missing values are excluded when calculating the cutoffs and remain missing
#' in the returned matrix. A requested variable with no observed values is
#' rejected. If the lower and upper empirical quantiles coincide, the function
#' also stops before modifying the data because applying equal cutoffs would
#' collapse every observed value to a constant.
#'
#' @param x A numeric predictor matrix with column names.
#'
#' @param cont_vars A character vector specifying the continuous variables to
#'   Winsorise. Variables not listed are returned unchanged.
#'
#' @param probs A numeric vector of length 2 giving the lower and upper
#'   percentile cutoffs in \code{[0, 1]}. Default is \code{c(0.01, 0.99)}.
#'   For example, \code{c(0.025, 0.975)} corresponds to a 5 percent Winsorisation
#'   (2.5 percent in each tail).
#'
#' @param zero_vars Optional character vector of variable names that should be
#'   excluded from Winsorisation. Use this argument for spike-at-zero variables
#'   whose full distribution should be preserved.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{\code{x}}{
#'     A numeric matrix with the specified variables Winsorised.
#'   }
#'   \item{\code{limits}}{
#'     A named numeric matrix containing the lower and upper cutoff values
#'     used for each variable. The matrix has two rows, \code{lower} and
#'     \code{upper}, and one column per variable.
#'   }
#' }
#' @examples
#' \dontrun{
#' set.seed(1)
#'
#' x <- cbind(
#'   age = c(rnorm(98, 60, 10), 5, 120),
#'   wt  = c(rnorm(98, 70, 15), 200, 30)
#' )
#'
#' out <- winsorize_cont_vars(
#'   x,
#'   cont_vars = c("age", "wt")
#' )
#'
#' out$limits
#' range(out$x[, "age"])
#' }
#' @keywords internal
#' @noRd
winsorize_cont_vars <- function(x,
                                cont_vars,
                                probs     = c(0.01, 0.99),
                                zero_vars = NULL) {
  
  if (!is.matrix(x))
    stop("! `x` must be a matrix.", call. = FALSE)
  if (is.null(colnames(x)))
    stop("! `x` must have column names.", call. = FALSE)
  if (!is.character(cont_vars) || length(cont_vars) == 0L)
    stop("! `cont_vars` must be a non-empty character vector.", call. = FALSE)
  if (!is.numeric(probs) || length(probs) != 2L ||
      any(probs < 0) || any(probs > 1) || probs[1L] >= probs[2L])
    stop("! `probs` must be a length-2 numeric vector in [0, 1] with ",
         "`probs[1] < probs[2]`.", call. = FALSE)
  
  missing_vars <- setdiff(cont_vars, colnames(x))
  if (length(missing_vars) > 0L)
    stop(paste0("! Variables not found in `x`: ",
                paste(missing_vars, collapse = ", "), "."), call. = FALSE)
  
  limits <- matrix(NA_real_, nrow = 2L, ncol = length(cont_vars),
                   dimnames = list(c("lower", "upper"), cont_vars))
  
  for (v in cont_vars) {
    col <- x[, v]
    
    # Skip variables flagged as zero / spike: the user has indicated that
    # both zeros and the positive distribution carry substantive meaning,
    # so Winsorising the positive tail could distort the modelled
    # relationship. Limits are recorded as NA for transparency.
    if (v %in% zero_vars) {
      limits[, v] <- c(NA_real_, NA_real_)
      next
    }
    
    # Count usable observations rather than physical rows. In particular, an
    # all-NA column must not reach quantile(), where it would produce a backend-
    # dependent low-level error or unusable limits.
    observed <- col[!is.na(col)]

    if (length(observed) == 0L) {
      stop(
        "Cannot Winsorise variable `", v,
        "`: it contains no non-missing values.",
        call. = FALSE
      )
    }

    # Retain the existing behavior for a single usable observation: there is
    # no empirical distribution to Winsorise, so leave the column unchanged
    # and record unavailable limits.
    if (length(observed) < 2L) {
      limits[, v] <- c(NA_real_, NA_real_)
      next
    }

    q <- stats::quantile(observed, probs = probs, names = FALSE)

    if (any(!is.finite(q))) {
      stop(
        "Cannot Winsorise variable `", v,
        "`: non-finite quantile limits were produced.",
        call. = FALSE
      )
    }

    if (q[1L] >= q[2L]) {
      # Equal cutoffs map every value below them upward and every value above
      # them downward, making the entire observed column constant. Detect the
      # cause here instead of passing a rank-deficient variable to MFPI.
      stop(
        "Cannot Winsorise variable `", v,
        "`: the lower and upper quantile limits are identical. The requested ",
        "probabilities would collapse the variable to a single value; choose ",
        "wider cutoffs or exclude this variable from Winsorisation.",
        call. = FALSE
      )
    }

    limits["lower", v] <- q[1L]
    limits["upper", v] <- q[2L]
    
    col[!is.na(col) & col < q[1L]] <- q[1L]
    col[!is.na(col) & col > q[2L]] <- q[2L]
    
    x[, v] <- col
  }
  
  list(x = x, limits = limits)
}
