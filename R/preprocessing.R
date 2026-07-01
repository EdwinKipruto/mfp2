#' Function that calculates an integer used to scale predictor
#' 
#' @details
#' For details on why scaling is useful see the corresponding section in the
#' documentation of \code{mfp2()}.
#' 
#' The determination of the scaling factor is independent (i.e. not affected 
#' by) shifts in the input data, as it only depends on the range of the 
#' input data.
#' 
#' Note that the estimation of powers is unaffected by scaling, the same powers 
#' are found for scaled input data. In extreme cases scaling is necessary to 
#' preserve accuracy, see Royston and Sauerbrei (2008).
#' This formula uses the scaling formula from Section 4.11.1 of 
#' Royston and Sauerbrei (2008). Further information can also be found in the 
#' Stata manual for mfp at https://www.stata.com/manuals/rfp.pdf.
#'
#' @param x a numeric vector already shifted to positive values (see 
#' \code{find_shift_factor()}). This function requires at least 2 distinct values to 
#' work.
#' 
#' #' @examples
#' x = 1:1000
#' find_scale_factor(x)
#' 
#' @return 
#' An integer that can be used to scale `x` to a reasonable range. For binary
#' variables 1 is returned.
#' 
#' @references 
#' Royston, P., and Sauerbrei, W., 2008. \emph{Multivariable Model - Building: 
#' A Pragmatic Approach to Regression Anaylsis based on Fractional Polynomials 
#' for Modelling Continuous Variables. John Wiley & Sons.}
#' 
#' @export
find_scale_factor <- function(x) {

  n_unique <- length(unique(x))
  
  if (n_unique == 1)
    stop("! Input data must not be constant.", 
         "i All values of x are identical, hence log(max(x)-min(x)) = log(0) is not defined.")
  
  if (n_unique == 2)
    return(1)
    
  lrange <- log10(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  
  10^(sign(lrange) * floor(abs(lrange)))
}

#' Function that calculates a value used to shift predictor
#' 
#' @details
#' For details on why shifting is necessary see the corresponding section in the
#' documentation of \code{mfp2()}.
#' 
#' This function implements the formula in Section 4.7 of Royston and 
#' Sauerbrei (2008).
#' 
#' @param x a numeric vector.
#' 
#' @examples
#' x = 1:1000
#' find_shift_factor(x)
#' 
#' @return 
#' A numeric value that can be used to shift `x` to positive values. 
#' If all values are positive, or if `x` is binary then 0 is returned. 
#' 
#' @references 
#' Royston, P., and Sauerbrei, W., 2008. \emph{Multivariable Model - Building: 
#' A Pragmatic Approach to Regression Anaylsis based on Fractional Polynomials 
#' for Modelling Continuous Variables. John Wiley & Sons.}
#' 
#' @export
find_shift_factor <- function(x) {
  
  n_unique <- length(unique(x))
  
  if (all(x > 0) || n_unique <= 2) 
    return(0)
  
  difx <- diff(sort(x))
  eps <- min(difx[difx != 0], na.rm = TRUE)
  
  eps - min(x, na.rm = TRUE)
}

#' Shift and scale vector x
#' 
#' A function that is used to shift x values to positive values if it contains
#' negative or zero values.If all values of x are positive then the original
#' values of x is returned without shifting but scaled if the scaling factor is
#' not equal to 1. If x has already been shifted and scaled then the function does
#' nothing.
#' 
#' @param x A vector of predictor variable
#' @param scale scaling factors for x of interest. Must be positive integers.
#' Default is NULL and  scaling factors are automatically estimated using
#' find_scale_factor() function else it uses user supplied scaling factors. If no scaling
#' is needed just use scale = 1
#' @param shift adjustment factors required for shifting x to positive
#' values. Default is NULL and adjustment factors are estimated automatically
#' using find_shift_factor() function
#' @examples
#' x = 1:1000
#' apply_shift_scale(x)
#' 
#' @returns 
#' A numeric value that has been shifted and scaled.
#'  
#' @export
apply_shift_scale <- function(x, scale = NULL, shift = NULL) {
  # restrict x to be a vector not matrix
  if (is.matrix(x))
    stop("x must be a vector not a matrix", call. = FALSE)
  
  N <- length(x)
  # If adjustment factors are NULL we use R&S formula to shift x to positive values
  if (is.null(shift)) {
    # Check whether all x values are positive. If true-No need of shifting
    if (all(x > 0)) {
      x <- x
    } else {
      # estimate adjustment factors
      shift <- find_shift_factor(x)
      # Shift x to positive
      x <- x + shift
    }
    # use the adjustment factors supplied by the user to shift x to positive values
  } else {
    # term + a
    x <- x + shift
    # check whether all x are now positive
    if (!all(x > 0))
      stop(
        "The minimum value of x after shifting x is ",
        min(x, na.rm = TRUE),
        " which is not > 0. Check your adjustment factors",
        call. = FALSE
      )
    }
  
  # if scale is NULL then scale x for computational stability using R&S formula
  if (is.null(scale)) { # No scaling
    x <- x / find_scale_factor(x)
  } else {
    x <- x / scale
  }
  return(x)
}


#' Validate rank and estimability of a default-interface design matrix
#'
#' Checks that a numeric design matrix supplied to the default matrix interface
#' has complete column names, contains no constant or non-informative columns,
#' and is full column rank after adding an intercept when appropriate.
#'
#' This is an internal defensive check for the matrix/default interface. Unlike
#' the formula interface, the default interface receives only a numeric design
#' matrix and cannot know whether columns came from raw variables, dummy coding,
#' ordered-factor polynomial contrasts, spline bases, or user-defined features.
#' Therefore this helper validates the actual estimability of the supplied
#' columns rather than trying.
#' Therefore this helper validates the actual estimability of the supplied
#' columns rather than trying to infer their origin from column names.
#'
#' For Cox models, `intercept` should be `FALSE` because Cox partial-likelihood
#' models do not include an intercept. For ordinary Gaussian, binomial, poisson,
#' and other GLM-style fits, `intercept` should usually be `TRUE`.
#'
#' @param x A numeric design matrix, or an object coercible to a matrix. Columns
#'   are candidate variables supplied to the default interface.
#' @param intercept Logical scalar. If `TRUE`, an intercept column is included
#'   when checking matrix rank. If `FALSE`, rank is checked on `x` alone.
#'
#' @return Invisibly returns `TRUE` if the design passes validation.
#'
#' @details
#' A column is treated as non-informative when, after removing non-finite values,
#' it has at most one unique value. Rank deficiency is detected using `qr()`.
#' If the design is rank-deficient, the helper reports columns identified by the
#' QR pivot as aliased or not estimable.
#'
#' This check is intentionally conservative. It is meant to fail early for input
#' matrices that cannot support valid model-comparison tests in the MFP
#' selection cycle. Step-level rank checks are still needed because a candidate
#' can become non-estimable in a particular adjustment model even if the initial
#' design passes this validation.
#'
#' @keywords internal
#' @noRd
validate_default_design_rank <- function(x, intercept = TRUE) {
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  
  if (is.null(colnames(x)) || anyNA(colnames(x)) || any(colnames(x) == "")) {
    stop("`x` must have complete, non-empty column names.", call. = FALSE)
  }
  
  bad_constant <- vapply(
    seq_len(ncol(x)),
    function(j) {
      z <- x[, j]
      z <- z[is.finite(z)]
      length(unique(z)) <= 1L
    },
    logical(1L)
  )
  
  if (any(bad_constant)) {
    stop(
      sprintf(
        "The following columns in `x` are constant or non-informative: %s.",
        paste(colnames(x)[bad_constant], collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  X <- if (intercept) cbind("(Intercept)" = 1, x) else x
  
  qr_x <- qr(X)
  full_rank <- qr_x$rank == ncol(X)
  
  if (!full_rank) {
    pivot <- qr_x$pivot
    aliased_positions <- pivot[(qr_x$rank + 1L):ncol(X)]
    aliased_names <- colnames(X)[aliased_positions]
    aliased_names <- setdiff(aliased_names, "(Intercept)")
    
    stop(
      sprintf(
        paste0(
          "The supplied design matrix is rank-deficient. ",
          "The following column(s) are aliased or not estimable in the initial design: %s. ",
          "Remove redundant columns or provide a design in which each candidate variable ",
          "adds an estimable degree of freedom."
        ),
        paste(aliased_names, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}
