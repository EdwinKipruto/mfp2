#' Functions to transform a variable using fractional polynomial powers or acd
#' 
#' These functions generate fractional polynomials for a variable similar to
#' `fracgen` in Stata. `transform_vector_acd` generates the acd transformation
#' for a variable.
#' 
#' @details 
#' The fp transformation generally transforms `x` as follows. For each pi in
#' `power` = (p1, p2, ..., pn) it creates a variable x^pi and returns the
#' collection of variables as a matrix. It may process the data using 
#' shifting, scaling and centering as desired. 
#' 
#' A special case are repeated powers, i.e. when some pi = pj. In this case, 
#' the fp transformations are given by x^pi and x^pi * log(x). In case
#' more than 2 powers are repeated they are repeatedly multiplied with 
#' log(x) terms, e.g. pi = pj = pk leads to x^pi, x^pi * log(x) and 
#' x^pi * log(x)^2. 
#' 
#' Binary variables are not transformed, but may only be centered. 
#' 
#' @section Data processing: 
#' An important note on data processing. Variables are shifted and scaled 
#' before being transformed by any powers. That is to ensure positive values
#' and reasonable scales. Note that scaling does not change the estimated 
#' powers, see also [find_scale_factor()].
#' 
#' However, they are centered after transformation. 
#' That is to ensure that the correlation between variables stay intact, 
#' as centering before transformation would affect them. This is described
#' in Sauerbrei et al (2006), as well as in the Stata manual of `mfp`.
#' 
#' @param x a vector of a predictor variable.
#' @param power a numeric vector indicating the FP power. Default is 1 (linear). 
#' Must be a vector of length 2 for acd transformation.
#' @param scale scaling factor for x of interest. Must be a positive integer
#' or `NULL`. Default is 1, meaning no scaling is applied. 
#' If `NULL`, then scaling factors are automatically estimated by the
#' program. 
#' @param shift shift required for shifting x to positive values. Default is 0, 
#' meaning no shift is applied. If `NULL` then the shift is estimated 
#' automatically using the Royston and Sauerbrei formula iff any `x` <= 0.
#' @param powers passed to [fit_acd()].
#' @param acd_parameter a list usually returned by [fit_acd()]. In particular, 
#' it must have components that define `beta0`, `beta1`, `power`, `shift` and 
#' `scale` which are to be applied when using the acd transformation in 
#' new data.
#' @param center Specification of centering for variable using
#' the mean i.e. `f(x) - mean(f(x))` for continuous variables and 
#' `x - min(x)` for binary variables. Default is no centering.
#' 
#' @details 
#' This function does not handle column names, as this function may be
#' called with an unnamed vector as well. Thus, the returned matrix is
#' always unnamed. 
#' 
#' @return 
#' Returns a matrix of transformed variable(s). The number of columns
#' depends on the number of powers provided, the number of rows is equal to the
#' length of `x`. If all powers are `NA`, then this function returns `NULL`.
#' In case an acd transformation is applied, the acd term is returned
#' as the second column of the matrix. 
#' 
#' @references 
#' Sauerbrei, W., Meier-Hirmer, C., Benner, A. and Royston, P., 2006. 
#' \emph{Multivariable regression model building by using fractional 
#' polynomials: Description of SAS, STATA and R programs. 
#' Comput Stat Data Anal, 50(12): 3464-85.}
#' 
#' @export
transform_vector_fp <- function(x, 
                                power = 1,
                                scale = 1, 
                                center = FALSE, 
                                shift = 0) {
  
  if (all(is.na(power))) { 
    # variable omitted
    return(NULL)
  }
  
  # do not transform x when it is a two level variable but center if necessary
  if (length(unique(x)) <= 2) {
    x <- as.matrix(x)
    # center using the lower(minimum) of the two distinct values of the covariates
    # as also done in stata mfp by Patrick Royston
    if (center) {
      return(x - min(x))
    } else return(x)
  } 
  
  # process input data by shifting and scaling
  if (is.null(shift)) {
    shift <- find_shift_factor(x)
  }
  
  if (is.null(scale)) {
    scale <- find_scale_factor(x)
  }
  
  x <- x + shift
  x <- x / scale
  
  # sort the powers
  # see Royston and Altman 1994 near equation 6 for explanation
  # note that sort removes NAs 
  power <- sort(power)

  # transform data
  x_trafo <- matrix(NA, nrow = length(x), ncol = length(power))
  # transform using first power
  x_trafo[, 1] <- transform_vector_power(x, power[1])
  # transform other powers via loop if necessary
  for (j in seq_len(length(power) - 1)) { 
    k <- j + 1
    if (power[k] == power[j]) {
      # the subsequent power is repeated, e.g. 1, 2, 2
      # repeatedly multiply with log(x), thereby creating powers of log(x)
      x_trafo[, k] <- x_trafo[, j] * log(x)
    } else {
      # the subsequent power is not repeated, e.g. 1, 2, 3
      x_trafo[, k] <- transform_vector_power(x, power[k])
    }
  }

  if (center) {
    x_trafo <- scale(x_trafo, scale = FALSE)
  }
  
  x_trafo
}

#' @describeIn transform_vector_fp Function to generate acd transformation.
#' @export
transform_vector_acd <- function(x, 
                                 power = c(1, 1),
                                 shift = 0, 
                                 powers = NULL, 
                                 scale = 1, 
                                 center = FALSE, 
                                 acd_parameter = NULL) {
  
  if (length(power) != 2) 
    stop("! power must be of length two.", 
         sprintf("i The length of powers supplied is %d.", length(power)))
  
  if (all(is.na(power))) {
    return(NULL)
  } 
  
  if (is.null(acd_parameter)) {
    # estimate acd(x)
    x_acd <- fit_acd(x, powers = powers, shift = shift, scale = scale)$acd
  } else x_acd <- do.call(apply_acd, modifyList(acd_parameter, list(x = x)))
  
  # apply fp transform on x and acd(x)
  x_acd <- transform_vector_fp(x = x_acd, power = power[2],
                               scale = scale, shift = shift, 
                               center = center)
  x_fp <- transform_vector_fp(x = x, power = power[1],
                              scale = scale, shift = shift,
                              center = center)

  cbind(x_fp, x_acd)
}

#' Function to transform each column of matrix using final FP powers or acd
#' 
#' @param x a matrix with all continuous variables shifted and scaled.
#' @param power_list a named list of FP powers to be applied to the columns of
#' `x`. Only variables named in this list are transformed.
#' @param center a named logical vector specifying whether the columns in `x`
#' should be centered.
#' @param acdx a named logical vector specifying the use of acd transformation.
#' @param keep_x_order a logical indicating whether the order of columns
#' should be kept as in the input matrix `x`, of if the columns should be 
#' ordered according to `power_list`. The default is `FALSE`, since 
#' the ordering by `power_list` reflects the `xorder` argument in [mfpa()].
#' @param acd_parameter_list a named list. Only required when transformation
#' are to be applied to new data. Entries must correspond to the entries where
#' `acdx` is set to `TRUE`. Each components is to be passed to 
#' [transform_vector_acd()]. The default value `NULL` indicates that the
#' parameters for the acd transformations are to be estimated.
#' 
#' @details 
#' For details on the transformations see [transform_vector_fp()] and
#' [transform_vector_acd()].
#' 
#' @section Column names: 
#' Generally the original variable names are suffixed with ".i", where
#' i enumerates the powers for a given variable in `power_list`. If a term
#' uses an acd transformation, then the variable is prefixed with "(A)".
#' 
#' @return 
#' If all elements of `power_list` are `NA` then this function returns `NULL`.
#' Otherwise a matrix is returned with transformed variables as named in 
#' `power_list`. The number of columns may possibly be different to the 
#' input matrix due to higher order FP transformations.
#' 
#' @export
transform_matrix <- function(x,
                             power_list, 
                             center, 
                             acdx, 
                             keep_x_order = FALSE, 
                             acd_parameter_list = NULL) {

  if (all(is.na(unlist(power_list)))) {
    # all variables were eliminated
    return(NULL)
  }

  if (keep_x_order)
    # reorder power_list to be in same order as columns in x
    power_list <- power_list[order(match(names(power_list), colnames(x)))]

  # only consider variables in power_list
  names_vars <- names(power_list)
  x <- x[, names_vars, drop = FALSE]
  center <- center[names_vars]
  acdx <- acdx[names_vars]

  x_trafo = list()
  for (name in names_vars) {
    if (acdx[name]) {
      # apply acd transformation
      x_trafo[[name]] <- transform_vector_acd(
        x[, name], power = power_list[[name]], center = center[name], 
        acd_parameter = acd_parameter_list[[name]]
      )
      if (!is.null(x_trafo[[name]]))
        # note that acd components are always in the second column
        colnames(x_trafo[[name]]) <- c(
          paste0(name, ".1"), paste0("(A)", name, ".1")
        )
    } else {
      # apply fp transform
      x_trafo[[name]] <- transform_vector_fp(
        x[, name], power = power_list[[name]], center = center[name]
      )
      if (!is.null(x_trafo[[name]]))
        # note that acd components are always in the second column
        colnames(x_trafo[[name]]) <- paste0(
          name, ".", seq_along(power_list[[name]])
        )
    }
  }
  
  do.call(cbind, x_trafo)
}

#' Simple function to transform vector by a single power
#' 
#' @param x a vector of a predictor variable.
#' @param power single power. 
transform_vector_power <- function(x,
                                   power = 1) {
  if (power == 0)
    return(log(x))
  
  x^power
}
