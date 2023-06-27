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
#' shifting and scaling as desired. Centering has to be done after the 
#' data is transformed using these functions, if desired. 
#' 
#' A special case are repeated powers, i.e. when some pi = pj. In this case, 
#' the fp transformations are given by x^pi and x^pi * log(x). In case
#' more than 2 powers are repeated they are repeatedly multiplied with 
#' log(x) terms, e.g. pi = pj = pk leads to x^pi, x^pi * log(x) and 
#' x^pi * log(x)^2. 
#' 
#' Note that the powers pi are assumed to be sorted. That is, this function 
#' sorts them, then proceeds to compute the transformation. For example, 
#' the output will be the same for `power = c(1, 1, 2)` and
#' `power = c(1, 2, 1)`. This is done to make sense of repeated powers and 
#' to uniquely define FPs. In case an ACD transformation is used, there is a 
#' specific order in which powers are processed, which is always the same (but 
#' not necessarily sorted). 
#' Thus, throughout the whole package powers will always be given and processed
#' in either sorted, or ACD specific order and the columns of the matrix 
#' returned by this function will always align with the powers used
#' throughout this package.
#' 
#' Binary variables are not transformed, unless `check_binary` is set to
#' `FALSE`. This is usually not necessary, the only special case to set it to 
#' `FALSE` is when a single value is to be transformed during prediction (e.g.
#' to transform a reference value). When this is done, binary variables are
#' still returned unchanged, but a single value from a continuous variable will
#' be transformed as desired by the fitted transformations. For model fit, 
#' `check_binary` should always be at its default value.
#' 
#' @section Data processing: 
#' An important note on data processing. Variables are shifted and scaled 
#' before being transformed by any powers. That is to ensure positive values
#' and reasonable scales. Note that scaling does not change the estimated 
#' powers, see also [find_scale_factor()].
#' 
#' However, they may be centered after transformation. This is not done by
#' these functions.
#' That is to ensure that the correlation between variables stay intact, 
#' as centering before transformation would affect them. This is described
#' in Sauerbrei et al (2006), as well as in the Stata manual of `mfp`.
#' Also, centering is not recommended, and should only be done for the final
#' model if desired.
#' 
#' @param x a vector of a predictor variable.
#' @param power a numeric vector indicating the FP power. Default is 1 (linear). 
#' Must be a vector of length 2 for acd transformation. Ignores `NA`, unless
#' an ACD transformation is applied in which case power must be a numeric 
#' vector of length 2, and `NA` indicated which parts are used for the final 
#' FP.
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
#' @param name character used to define names for the output matrix. Default
#' is `NULL`, meaning the output will have unnamed columns.
#' @param check_binary a logical indicating whether or not input `x` is checked
#' if it is a binary variable (i.e. has only two distinct values). The default
#' `TRUE` usually only needs to changed when this function is to be used to 
#' transform data for predictions. See Details.
#' 
#' @return 
#' Returns a matrix of transformed variable(s). The number of columns
#' depends on the number of powers provided, the number of rows is equal to the
#' length of `x`. The columns are sorted by increased power.
#' If all powers are `NA`, then this function returns `NULL`.
#' In case an acd transformation is applied, the output is a list with two 
#' entries. The first `acd` is the matrix of transformed variables, the acd 
#' term is returned as the last column of the matrix (i.e. in case that the 
#' power for the normal data is `NA`, then it is the only column in the matrix).
#' The second entry `acd_parameter` returns a list of estimated parameters
#' for the ACD transformation, or simply the input `acd_parameter` if it was
#' not `NULL`.
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
                                shift = 0, 
                                name = NULL, 
                                check_binary = TRUE) {
  
  if (all(is.na(power))) { 
    # variable omitted
    return(NULL)
  }
  
  # do not transform x when it is a two level variable 
  if (check_binary && length(unique(x)) <= 2) {
    x <- as.matrix(x)
    if (!is.null(name))
      colnames(x) <- name_transformed_variables(name, 1)
    
    return(x)
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
  
  if (!is.null(name)) {
    colnames(x_trafo) <- name_transformed_variables(name, ncol(x_trafo))
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
                                 acd_parameter = NULL, 
                                 name = NULL) {
  
  if (length(power) != 2) 
    stop("! power must have length two.", 
         sprintf("i The length of powers supplied is %d.", length(power)))
  
  if (all(is.na(power))) {
    return(NULL)
  } 
  
  if (is.null(acd_parameter)) {
    # estimate acd(x)
    acd_parameter <- fit_acd(x, powers = powers, shift = shift, scale = scale)
    x_acd <- acd_parameter$acd
    # no need to store acd further
    acd_parameter$acd <- NULL
  } else x_acd <- do.call(apply_acd, modifyList(acd_parameter, list(x = x)))
  
  name_acd <- NULL
  if (!is.null(name))
    name_acd <- paste0("A_", name)
  
  # apply fp transform on x (if required) and acd(x)
  # if any of these is NA, transform_vector_fp returns NULL and thus the 
  # component is not used in the final result, as desired
  x_acd <- transform_vector_fp(x = x_acd, power = power[2],
                               scale = scale, shift = shift, name = name_acd)
  x_fp <- transform_vector_fp(x = x, power = power[1],
                              scale = scale, shift = shift, name = name)

  list(
    acd = cbind(x_fp, x_acd),
    acd_parameter = acd_parameter
  )
}

#' Function to transform each column of matrix using final FP powers or acd
#' 
#' @param x a matrix with all continuous variables shifted and scaled.
#' @param power_list a named list of FP powers to be applied to the columns of
#' `x`. Only variables named in this list are transformed.
#' @param center a named logical vector specifying whether the columns in `x`
#' should be centered. Centering will occur after transformations and will be
#' done separately for each individual column of the transformed data matrix. 
#' @param acdx a named logical vector specifying the use of acd transformation.
#' @param keep_x_order a logical indicating whether the order of columns
#' should be kept as in the input matrix `x`, of if the columns should be 
#' ordered according to `power_list`. The default is `FALSE`, since 
#' the ordering by `power_list` reflects the `xorder` argument in [mfp2()].
#' @param acd_parameter_list a named list. Only required when transformation
#' are to be applied to new data. Entries must correspond to the entries where
#' `acdx` is set to `TRUE`. Each components is to be passed to 
#' [transform_vector_acd()]. The default value `NULL` indicates that the
#' parameters for the acd transformations are to be estimated.
#' @param check_binary passed to [transform_vector_fp()].
#' 
#' @details 
#' For details on the transformations see [transform_vector_fp()] and
#' [transform_vector_acd()].
#' 
#' @section Column names: 
#' Generally the original variable names are suffixed with ".i", where
#' i enumerates the powers for a given variable in `power_list`. If a term
#' uses an acd transformation, then the variable is prefixed with `A_`.
#' 
#' @return 
#' If all elements of `power_list` are `NA` then this function returns `NULL`.
#' Otherwise a list with three entries: the first `x_transformed` is a matrix 
#' with transformed variables as named in `power_list`. 
#' The number of columns may possibly be different to the 
#' input matrix due to higher order FP transformations.
#' The second entry `centers` stores the values used to center the variables if
#' for any variable `center = TRUE` (note that usually all variables are 
#' centered, or none of them).
#' The third entry `acd_parameter` stores a named list of estimated 
#' `acd_parameters`. May be empty if no ACD transformation is applied.
#' 
#' @export
transform_matrix <- function(x,
                             power_list, 
                             center, 
                             acdx, 
                             keep_x_order = FALSE, 
                             acd_parameter_list = NULL,
                             check_binary = TRUE) {

  if (all(is.na(unlist(power_list)))) {
    # all variables were eliminated
    return(NULL)
  }
  
  # power_list, center and acdx must have names
  if (is.null(names(power_list))) 
    stop("! List power_list must have names.")
  
  if (is.null(names(center))) 
    stop("! Vector center must have names.")
  
  if (is.null(names(acdx)))
    stop("! Vector acdx must have names.")

  if (keep_x_order)
    # reorder power_list to be in same order as columns in x
    power_list <- power_list[order(match(names(power_list), colnames(x)))]

  # only consider variables in power_list
  names_vars <- names(power_list)
  x <- x[, names_vars, drop = FALSE]
  center <- center[names_vars]
  acdx <- acdx[names_vars]

  x_trafo = list()
  acd_parameter = list()
  for (name in names_vars) {
    if (acdx[name]) {
      # apply acd transformation
      acd <- transform_vector_acd(
        x[, name], power = power_list[[name]],  
        acd_parameter = acd_parameter_list[[name]], name = name
      )
      x_trafo[[name]] <- acd$acd
      acd_parameter[[name]] <- acd$acd_parameter
    } else {
      # apply fp transform
      x_trafo[[name]] <- transform_vector_fp(
        x[, name], power = power_list[[name]], name = name, 
        check_binary = check_binary
      )
    }
  }
  
  # create output matrix
  x_transformed <- do.call(cbind, x_trafo)
  
  # also center values if required
  centers = NULL
  if (any(center)) {
    x_transformed <- center_matrix(x_transformed)
    centers <- attr(x_transformed, "scaled:center")
  }
  
  list(
    x_transformed = x_transformed,
    centers = centers,
    acd_parameter = acd_parameter
  )
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

#' Simple function to center data
#' 
#' @param mat a transformed data matrix. 
#' @param centers a vector of centering values. Length must be equal to the 
#' number of columns in `mat`. If `NULL` (default) then 
#' centering values are determined by the function (see Details).
#' 
#' @details 
#' Centering is done by means for continuous variables (i.e. more than 2
#' distinct values), and the minimum for binary variables. 
#' 
#' It is assumed all categorical variables in the data are represented by 
#' binary dummy variables. 
#' 
#' @return 
#' Transformed data matrix. Has an attribute `scaled:center` that stores 
#' values used for centering.
#' 
#' @export
center_matrix <- function(mat, centers = NULL) {
  
  if (is.null(centers)) {
    centers <- colMeans(mat) 
    
    # identify the column with at most 2 unique values
    index_binary <- which(apply(mat, 2, function(x) length(unique(x))) <= 2)
    
    # replace the means of binary variables with the minimum 
    if (any(index_binary)) {
      centers[index_binary] <- apply(mat[,index_binary, drop = FALSE], 2, min)
    }
  }
  
  scale(mat, center = centers, scale = FALSE)
}

#' Helper function to name transformed variables
#' 
#' @param name character with name of variable being transformed.
#' @param n_powers number of resulting variables from FP-transformation.
#' @param acd logical indicating the use of ACD-transformation
#' 
#' @return 
#' Character vector of names of length `n_powers`.
name_transformed_variables <- function(name, n_powers, acd = FALSE) {
  if (!acd) {
    paste0(name, ".", 1:n_powers)
  } else {
    c(paste0(name, ".1"), paste0("A_", name, ".1"))
  }
}

#' Simple function to create dummy variables for ordinal variables
#' 
#' @param data a dataframe containing the ordinal variable. 
#' @param var_names names of ordinal variables in the data.
#' 
#' @details 
#' This function creates based on ordinal coding described in Royston and Sauerbrei (2008) book (chapter 3, Table 3.5)
#' @return 
#' dataframe with new dummy variables.
#' 
#' @export
create_dummy_ordinal <- function(data, var_names) {
  # assert that data must be provided
  if (missing(data))
    stop("! data argument is missing.\n",
         "i An input data.frame is required for the use of create_dummy_ordinal.",
         call. = FALSE)
  if (!is.data.frame(data))
    stop ("The data must be a data.frame")
  
  for (col in var_names) {
    if (is.null(data[[col]]) || length(data[[col]]) == 0) {
      stop(paste("variable name ", col, "is empty or does not exist."))
    }
    
    unique_levels <- unique(data[[col]])
    if (length(unique_levels) == 1) {
      warning(paste("var_names ", col, "has only one unique level. Skipping dummy variable creation."))
      next
    }
    
    levels_list <- lapply(seq_along(unique_levels)[-length(unique_levels)], function(i) unique_levels[seq(i)])
    
    for (i in seq_along(levels_list)) {
      level <- levels_list[[i]]
      var_name <- paste0(col, "_", paste(level, collapse = "_"))
      data[[var_name]] <- as.integer(!data[[col]] %in% level)
    }
  }
  
  return(data)
}

