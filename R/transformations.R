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
#' @param s passed to [acd()].
#' @param center Specification of centering for variable using
#' the mean i.e. `f(x) - mean(f(x))` for continuous variables and 
#' `x - min(x)` for binary variables. Default is no centering.
#' 
#' @return 
#' Returns a matrix of transformed variable(s). The number of columns
#' depends on the number of powers provided, the number of rows is equal to the
#' length of `x`. If all powers are `NA`, then this function returns `NULL`.
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
                                 s = NULL, 
                                 scale = 1, 
                                 center = FALSE) {
  
  if (length(power) != 2) 
    stop("! power must be of length two.", 
         sprintf("i The length of powers supplied is %d.", length(power)))
  
  if (all(is.na(power))) {
    return(NULL)
  } 
  
  # transform x and acdx using the supplied powers
  x_acd <- acd(x, power = NULL, shift = shift, s = s, scale = scale)$acd
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
#' @param power_list a named list of FP powers. 
#' @param center a named logical vector specifying whether the columns in `x`
#' should be centered.
#' @param acdx a named logical vector specifying the use of acd transformation.
#' 
#' @return 
#' If all elements of `power_list` are `NA` then this function returns `NULL`.
#' Otherwise a matrix is returned with transformed data, possibly with more 
#' columns than the input matrix due to higher order FP transformations.
transform_matrix <- function(x,
                             power_list, 
                             center, 
                             acdx) {
  # check whether all power_list are equal to NA-all variables were eliminated
  if (all(is.na(unlist(power_list)))) {
    # return NULL
    xtr.out <- NULL
  } else {
    namx <- names(power_list)
    x <- x[, namx, drop = F]
    # subset scale, center and shift using names of powers selected-
    # scale <- scale[namx] # we assume x has been shifted and scaled
    center <- center[namx]
    acdx <- acdx[namx]
    # transform acd variables if any
    if (any(acdx)) {
      # select FP powers for acd variables
      pow.acd <- power_list[acdx]
      # check whether the acd variables were selected
      if (all(is.na(unlist(pow.acd)))) {
        xtr1 <- NULL
      } else {
        xnames1 <- names(pow.acd)
        x.acdx <- x[, xnames1, drop = F]
        xtrans.acd <- vector(mode = "list", length = length(pow.acd))
        for (i in seq_along(xnames1)) {
          xtrans.acd[[i]] <- transform_vector_acd(
            x = x.acdx[, i, drop = T],
            power = pow.acd[[i]], scale = 1, # scale = scale[i],
            center = center[xnames1][i], shift = NULL
          ) # shift = shift[i])
        }
        # cbind acd transformed variables
        xtr1 <- do.call(cbind, xtrans.acd)
        # rename acd transformed variables. they start with letter "A"
        snames <- unlist(sapply(xnames1, function(x) paste0(c(x, paste0("(A)", x)), ".", c(1, 1)),
                                simplify = F, USE.NAMES = F
        ))
        colnames(xtr1) <- snames[which(!is.na(unlist(pow.acd)))]
      }
      # Variables without acd transformations
      pow <- power_list[!acdx]
      if (all(is.na(unlist(pow)))) {
        xtr2 <- NULL
      } else {
        # get rid of unselected variables
        # transform_vector_fp returns NULL when power = NA
        pow2 <- pow[!is.na(pow)]
        xnames2 <- names(pow2)
        xx <- x[, xnames2, drop = F]
        xtrans <- vector(mode = "list", length = length(pow2))
        for (i in seq_along(xnames2)) {
          # x is already scaled and shifted so transform_vector_fp
          # shifting is irrelevant
          xtrans[[i]] <- transform_vector_fp(
            x = xx[, i, drop = T], power = pow2[[i]], scale = 1, # scale = scale[i],
            center = center[xnames2][i], shift = NULL
          ) # shift = shift[i])
        }
        # cbind non-acd transformed variables-
        xtr2 <- do.call(cbind, xtrans)
        # rename xtr such that if x1 is fp2 then we have x1.1, x1.2 whereas fp1 is x1.1
        # we can as well use make.names(rep(names(fpp), lapply(fpp, length)), sep = ".")
        # but has undesirable names like x1, x1.1 instead of x1.1, x1.2
        colnames(xtr2) <- unlist(sapply(xnames2, function(x) paste0(x, ".", seq_along(pow2[[x]])),
                                        simplify = F, USE.NAMES = F
        ))
      }
      # combine xtr1 and xtr2
      xtr.out <- cbind(xtr1, xtr2)
      # Usual MFP without acd variables
    } else {
      # # Get rid of unselected variables denoted by NA in power_list
      fpp <- power_list[!is.na(power_list)]
      namxx <- names(fpp)
      # subset x
      x <- x[, namxx, drop = F]
      xtransx <- vector(mode = "list", length = length(fpp))
      for (i in seq_along(namxx)) {
        # x is already scaled and shifted so transform_vector_fp shifting is irrelevant
        xtransx[[i]] <- transform_vector_fp(
          x = x[, i, drop = T], power = fpp[[i]], scale = 1, # scale = scale[i],
          center = center[namxx][i], shift = 0
        ) # shift = shift[i])
      }
      xtr.out <- do.call(cbind, xtransx)
      colnames(xtr.out) <- unlist(sapply(namxx, function(x) paste0(x, ".", seq_along(fpp[[x]])),
                                         simplify = F, USE.NAMES = F
      ))
    }
  }
  
  xtr.out
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
