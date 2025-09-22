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
#' If a variable is specified in the \code{zero} or \code{catzero} arguments, 
#' nonpositive values (zero or negative) are not shifted. Instead, they are replaced 
#' with zero, and transformation is applied only to the positive values. This approach 
#' is useful in cases where nonpositive values have a qualitatively different interpretation 
#' (e.g., nonsmokers in smoking data) and should not be transformed in the same way 
#' as positive values.
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
#' @param zero Logical indicating whether only positive values of the variable 
#' should be transformed, with nonpositive values (zero or negative) set to zero. 
#' If \code{TRUE}, transformation is applied only to positive values; nonpositive values 
#' are replaced with zero before transformation. If \code{FALSE} (default), all values 
#' are shifted (if needed) to ensure positivity before transformation.
#' @param check_binary a logical indicating whether or not input `x` is checked
#' if it is a binary variable (i.e. has only two distinct values). The default
#' `TRUE` usually only needs to changed when this function is to be used to 
#' transform data for predictions. See Details.
#' 
#' @examples
#' z = 1:10
#' transform_vector_fp(z)
#' transform_vector_acd(z)
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
                                zero = FALSE,
                                check_binary = TRUE) {
  
  if (!is.logical(zero) || length(zero) != 1 || is.na(zero)) {
    stop("`zero` must be a single logical value (TRUE or FALSE).", call. = FALSE)
  }
  
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
  
  if (zero) {
    # Replace all zero or negative values with 0
    x[x <= 0] <- 0
    shift <- 0
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
    # Clean up invalid entries; can arise due to zero argument
    x_trafo[!is.finite(x_trafo[, k]), k] <- 0
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
                                 name = NULL,
                                 zero = FALSE) {
  
  if (length(power) != 2) 
    stop("! power must have length two.", 
         sprintf("i The length of powers supplied is %d.", length(power)))
  
  if (all(is.na(power))) {
    return(NULL)
  } 
  
  if (is.null(acd_parameter)) {
    # estimate acd(x)
    acd_parameter <- fit_acd(x, powers = powers, shift = shift, scale = scale, zero = zero)
    x_acd <- acd_parameter$acd
    # no need to store acd further
    acd_parameter$acd <- NULL
  } else {
    if (zero) {
      x[x <= 0] <- 0
    }
    x_acd <- do.call(apply_acd, modifyList(acd_parameter, list(x = x, zero = FALSE)))
  }
  
  name_acd <- NULL
  if (!is.null(name)) {
    name_acd <- paste0("A_", name)
  }
  
  # apply fp transform on x (if required) and acd(x)
  # if any of these is NA, transform_vector_fp returns NULL and thus the 
  # component is not used in the final result, as desired
  x_acd <- transform_vector_fp(x = x_acd, power = power[2], scale = scale, 
                               shift = shift, name = name_acd, zero = zero)
  x_fp <- transform_vector_fp(x = x, power = power[1], scale = scale, 
                              shift = shift, name = name, zero = zero)

  list(
    acd = cbind(x_fp, x_acd),
    acd_parameter = acd_parameter
  )
}

#' Transform each column of matrix using final FP powers or ACD transformation
#' 
#' This function applies FP and/or ACD transformations to the columns of a matrix,
#' optionally centers the transformed variables, and can generate binary indicators
#' for variables with nonpositive values (catzero). The spike-at-zero variables 
#' is supported via the `spike_decision` argument.
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
#' @param zero A named logical vector specifying, for each variable, whether only 
#' positive values should be transformed. If \code{TRUE}, the transformation is applied 
#' only to positive values; nonpositive values (zero or negative) are replaced with zero. 
#' If \code{FALSE}, all values are used (after shifting, if necessary). 
#' Variable names must match the corresponding column names in \code{x}. The default 
#' is \code{NULL}, meaning no variables are treated as zero-transformed.
#' @param catzero A named logical vector similar to \code{zero}, indicating which 
#' columns in \code{x} should treat nonpositive values as zero and additionally have a 
#' binary indicator created and included in the model. The vector must have names matching 
#' the column names in \code{x}. The default is \code{NULL}, meaning no categorical-zero 
#' variables are used.
#' @param spike_decision Named numeric vector with values 1, 2, or 3 specifying 
#' spike-at-zero handling for each variable. Default is `NULL`, meaning no 
#' spike-at-zero processing is applied. Value 1 includes both the transformed 
#' variable and binary indicator, 2 disables the spike/binary indicator, and 3 
#' keeps only the binary indicator.
#' @details 
#' For details on the transformations see [transform_vector_fp()] and
#' [transform_vector_acd()].
#' @details 
#' Transformations are handled by [transform_vector_fp()] and [transform_vector_acd()]. 
#' The `x_transformed` matrix contains transformed variables, optionally centered. 
#' Binary indicators are appended if specified in `catzero`. 
#' The `spike_decision` argument controls spike-at-zero (SAZ) processing: for variables 
#' with value 2, the binary indicator is disabled reducing to usual transformations;
#' for variables with value 3, only the binary indicator of SAZ is retained and
#' any FP/linear transformations are removed. This ensures correct handling 
#' of variables with SAZ.
#'
#' @section Column names:
#' Transformed variable names are based on the original names. FP terms are 
#' suffixed with ".i" to indicate the power index for that variable. Variables 
#' transformed using ACD are prefixed with "A_". Binary indicators created 
#' for nonpositive values (catzero) are suffixed with "_bin".
#'
#' @examples
#' x = matrix(1:100, nrow = 10)
#' colnames(x) = paste0("x", 1:ncol(x))
#' powx = setNames(replicate(ncol(x), c(1,2), simplify = FALSE), colnames(x))
#' center = setNames(rep(FALSE, ncol(x)), colnames(x))
#' acdx = setNames(rep(FALSE, ncol(x)), colnames(x))
#' transform_matrix(x, powx, center, acdx)
#' 
#' @return 
#' If all elements of `power_list` are `NA`, this function returns `NULL`.
#' Otherwise, it returns a list with four entries. The entry `x_transformed` 
#' is a matrix of transformed variables as specified in `power_list`, possibly 
#' centered. The number of columns may differ from the input matrix due to 
#' higher-order FP transformations. Binary variables are also added if 
#' `catzero` is not `NULL`. The entry `centers` contains the values used to 
#' center variables if `center = TRUE` (typically all variables are centered, 
#' or none of them). The entry `acd_parameter` is a named list of estimated ACD 
#' parameters, which may be empty if no ACD transformation is applied. The 
#' entry `x_trafo` is a named list of transformed `x` without centering or
#' `catzero` variables and ;  this last component is 
#' important for the spike-at-zero algorithm.
#'   
#' @export
transform_matrix <- function(x,
                             power_list, 
                             center, 
                             acdx, 
                             keep_x_order = FALSE, 
                             acd_parameter_list = NULL,
                             check_binary = TRUE,
                             zero = NULL,
                             catzero = NULL,
                             spike_decision = NULL) {
  #-------------------------
  # Input checks
  #-------------------------
  if (!is.matrix(x)) {
    stop("! 'x' must be a matrix.")
  }
  
  #-Ensure x has column names
  x_cols <- colnames(x)
  if (is.null(x_cols)) {
    stop("! Input data 'x' must have column names.")
  }
  
  if (!is.list(power_list)) {
    stop("! 'power_list' must be a list.")
    }
  
  if (all(is.na(unlist(power_list)))) {
    # All variables were eliminated, nothing to transform
    return(NULL)
  }
  
  if (is.null(names(power_list))) { 
    stop("! power_list must have names.")
  }
  
  # Ensure power_list variables are present in x
  pl_names <- names(power_list)
  missing_in_x <- setdiff(pl_names, x_cols)
  if (length(missing_in_x) > 0) {
    stop("! The following variables in 'power_list' are not in 'x': ",
         paste(missing_in_x, collapse = ", "))
  }
  
  # Center
  if (!is.logical(center)) {
    stop("! 'center' must be a logical vector.")
  }
  
  if (is.null(names(center))) {
    stop("! 'center' must have names.")
  }
  
  if (length(center) != length(power_list)) {
    stop("! 'center' must have same length as 'power_list'.")
  }
  
  # acdx
  if (!is.logical(acdx)) {
    stop("! 'acdx' must be a logical vector.")
  }
  
  if (is.null(names(acdx))) {
    stop("! 'acdx' must have names.")
  }
  
  if (length(acdx) != length(power_list)) {
    stop("! 'acdx' must have same length as 'power_list'.")
  }
  
  # zero
  if (is.null(zero)) {
   zero <- setNames(rep(FALSE, length(power_list)), names(power_list))
  } else {
    if (!is.logical(zero)) {
      stop("! 'zero' must be a logical vector.")
    }
    
    if (is.null(names(zero))) {
      stop("! 'zero' must have names.")
    } 
    zero <- zero[names(power_list)]
  }
  
  # catzero
  if (is.null(catzero)) {
    catzero <- setNames(rep(FALSE, length(power_list)), names(power_list))
  } else {
    if (!is.logical(catzero)) {
      stop("! 'catzero' must be a logical vector.")
    }
    
    if (is.null(names(catzero))) {
      stop("! 'catzero' must have names.")
    } 
    catzero <- catzero[names(power_list)]
  }
  
  # spike_decision
  if (!is.null(spike_decision)) {
    if (!is.numeric(spike_decision)) {
      stop("! 'spike_decision' must be numeric (values 1, 2, or 3).")
    }
    
    if (is.null(names(spike_decision))) {
      stop("! 'spike_decision' must have names.")
    }
    
    bad_vals <- setdiff(unique(spike_decision), c(1, 2, 3))
    if (length(bad_vals) > 0) {
      stop("! 'spike_decision' must only contain values 1, 2, or 3.")
    }
    
    # # Automatically set zero = TRUE for variables with spike_decision = 1
    # # this is important for centering. 
    # spike1_vars <- intersect(names(spike_decision)[spike_decision == 1], names(zero))
    # if (length(spike1_vars) > 0) zero[spike1_vars] <- TRUE
  }
  
  
  #-------------------------
  # Check zero variables
  #-------------------------
  if (any(zero)) {
    vars_to_check <- names(zero)[zero]
    bad_vars <- vars_to_check[
      apply(x[, vars_to_check, drop = FALSE], 2, function(col) all(col > 0, na.rm = TRUE))
    ]
    if (length(bad_vars) > 0) {
      warning("These variables were marked as 'zero = TRUE' but contain only positive values. ",
              "Resetting 'zero' and 'catzero' to FALSE for: ",
              paste(bad_vars, collapse = ", "))
      zero[bad_vars] <- FALSE
      catzero[bad_vars] <- FALSE
    }
  }
  
  #----------------------
  # Reorder to match x
  #----------------------
  if (keep_x_order) {
    # reorder power_list to be in same order as columns in x
    power_list <- power_list[order(match(names(power_list), x_cols))]
  }
  
  # only consider variables in power_list
  names_vars <- names(power_list)
  x <- x[, names_vars, drop = FALSE]
  center <- center[names_vars]
  acdx <- acdx[names_vars]
  zero <- zero[names_vars]
  catzero <- catzero[names_vars]
  
  #----------------------
  # Transformations
  #----------------------
  x_trafo <- list()
  acd_parameter <- list()
  for (name in names_vars) {
    if (acdx[name]) {
      # apply acd transformation
      acd <- transform_vector_acd(
        x[, name], power = power_list[[name]],  
        acd_parameter = acd_parameter_list[[name]], name = name, 
        zero = zero[name] 
      )
      x_trafo[[name]] <- acd$acd
      acd_parameter[[name]] <- acd$acd_parameter
    } else {
      # apply fp transform
      x_trafo[[name]] <- transform_vector_fp(
        x[, name], power = power_list[[name]], name = name, 
        check_binary = check_binary, zero = zero[name]
      )
    }
  }
  
  #----------------------
  # Apply spike_decision logic AFTER x_trafo is built
  #----------------------
  if (!is.null(spike_decision)) {
    # Variables where spike is not needed
    vars_no_spike <- names(spike_decision)[spike_decision == 2]
    intersect_vars <- intersect(vars_no_spike, names(catzero))
    if (length(intersect_vars) > 0) {
      catzero[intersect_vars] <- FALSE
    }
    
    # Variables where only spike binary is needed (remove FP/linear transform)
    vars_only_spike <- names(spike_decision)[spike_decision == 3]
    intersect_vars2 <- intersect(vars_only_spike, names(x_trafo))
    if (length(intersect_vars2) > 0) {
      x_trafo[intersect_vars2] <- NULL
    }
  }
  
  #-------------------------
  # Build output matrix
  #-------------------------
  x_transformed <- do.call(cbind, x_trafo)
  
  #----------------------
  # Handle catzero binary vars
  #----------------------
  cat_vars <- names(catzero)[catzero]
  
  if (length(cat_vars) > 0) {
    # Build list of binary variables (skip NA-only power entries)
    catzero_list <- lapply(cat_vars, function(v) {
      if (all(is.na(power_list[[v]]))) {
        return(NULL)  # skip
      } else {
        #as.integer(x[, v] > 0)
        as.integer(x[, v] == 0)  # 1 if zero, 0 if positive
      }
    })
    
    # Remove NULLs
    valid_idx <- !sapply(catzero_list, is.null)
    catzero_list <- catzero_list[valid_idx]
    cat_vars_valid <- cat_vars[valid_idx]
    names(catzero_list) <- cat_vars_valid
    
    # Only bind if there’s something left
    if (length(catzero_list) > 0) {
      catzero_matrix <- do.call(cbind, catzero_list)
      colnames(catzero_matrix) <- paste0(cat_vars_valid, "_bin")
      x_transformed <- cbind(x_transformed, catzero_matrix)
    } 
  } 
  
  #----------------------
  # Centering
  #----------------------
  centers <- NULL
  if (any(center)) {
    # Initialize expanded zero with FALSE for every transformed column
    zero_expanded <- setNames(rep(FALSE, ncol(x_transformed)), colnames(x_transformed))
    
    # For each original variable flagged zero=TRUE, mark its transformed
    # columns (FP or ACD) as zero-handled, but exclude *_bin columns.
    orig_zero_vars <- names(zero)[which(zero)]
    if (length(orig_zero_vars) > 0) {
      for (nm in orig_zero_vars) {
        # Match either "nm" or "A_nm" followed by dot or end-of-string.
        # This captures: nm, nm.1, nm.2, A_nm, A_nm.1, ...
        pat <- paste0("(^", nm, "(\\.|$))|(^A_", nm, "(\\.|$))")
        matches <- grep(pat, colnames(x_transformed), perl = TRUE, value = TRUE)
        if (length(matches) > 0) {
          # Exclude binary indicator columns (those ending with "_bin")
          matches <- matches[!grepl("_bin$", matches)]
          if (length(matches) > 0) zero_expanded[matches] <- TRUE
        }
      }
    }
    
    # Ensure any explicit *_bin columns are FALSE (defensive)
    bin_cols <- grep("_bin$", colnames(x_transformed), value = TRUE)
    if (length(bin_cols) > 0) zero_expanded[bin_cols] <- FALSE
    
    x_transformed <- center_matrix(mat = x_transformed, centers = NULL, zero = zero_expanded) # check centering of zero variables, do we need the positive part?
    centers <- attr(x_transformed, "scaled:center")
  }
  
  list(
    x_transformed = x_transformed,
    centers = centers,
    acd_parameter = acd_parameter,
    x_trafo = x_trafo
  )
}

#' Simple function to transform vector by a single power
#' 
#' @param x a vector of a predictor variable.
#' @param power single power.
#' @param zero Logical indicating whether only positive values of the variable 
#' should be transformed, with nonpositive values (zero or negative) set to zero. 
#' If \code{TRUE}, transformation is applied only to positive values; nonpositive values 
#' are replaced with zero before transformation. 
#' @return A vector of transformed values if power is not equal to 1
transform_vector_power <- function(x, power = 1, zero = FALSE) {
  
  if (zero) {
    x[x <= 0] <- 0
  }
  if (power == 0) {
    xt <- log(x)
  } else {
  xt <- x ^ power
  }
  xt[!is.finite(xt)] <- 0
  
  return(xt)
}

#' Simple function to center data
#' 
#' @param mat a transformed data matrix. 
#' @param centers a vector of centering values. Length must be equal to the 
#' number of columns in `mat`. If `NULL` (default) then 
#' centering values are determined by the function (see Details).
#' @param zero Optional named logical vector indicating which columns treat 
#' zero values specially. Names must match `mat` columns. Default `NULL` means 
#' no zero-specific handling.
#'
#' @details 
#' Centering is done by column means for continuous variables (more than 2 
#' distinct values) and by the minimum for binary variables. For variables 
#' with `zero = TRUE`, the mean is computed only over the non-zero values, 
#' while zero values remain at zero.
#' 
#' It is assumed all categorical variables in the data are represented by 
#' binary dummy variables. 
#' @examples
#' mat <- matrix(1:100, nrow = 10)
#' colnames(mat) <- paste0("x", 1:ncol(mat))
#' zero <- setNames(rep(FALSE, ncol(mat)), colnames(mat))
#' center_matrix(mat, zero = zero)
#' 
#' @return 
#' Transformed data matrix. Has an attribute `scaled:center` that stores 
#' values used for centering.
#' 
#' @export
center_matrix <- function(mat, centers = NULL, zero = NULL) {
  
  if (!is.matrix(mat)) {
    stop("! 'mat' must be a matrix.")
  }
  if (is.null(colnames(mat))) {
    stop("! 'mat' must have column names.")
  }
  
  # Validate zero
  if (!is.null(zero)) {
    if (!is.logical(zero)) {
      stop("! 'zero' must be a logical vector.")
    }
    if (is.null(names(zero))) {
      stop("! 'zero' must have names.")
    }
    if (!setequal(names(zero), colnames(mat))) {
      stop("! 'zero' names must match column names of 'mat'.")
    }
    zero <- zero[colnames(mat)]  # reorder to match mat
  } else {
    zero <- setNames(rep(FALSE, ncol(mat)), colnames(mat))
  }
  
  if (is.null(centers)) {
    centers <- numeric(ncol(mat))
    
    for (j in seq_len(ncol(mat))) {
      x <- mat[, j]
      is_binary <- length(unique(x)) <= 2
      
      if (zero[j]) {
        spike_mask <- x == 0
        centers[j] <- if (all(spike_mask)) 0 else mean(x[!spike_mask], na.rm = TRUE)
      } else if (is_binary) {
        # replace the means of binary variables with the minimum 
        centers[j] <- min(x, na.rm = TRUE)
      } else {
        centers[j] <- mean(x, na.rm = TRUE)
      }
    }
  }
  
  mat_centered <- scale(mat, center = centers, scale = FALSE)
  attr(mat_centered, "scaled:center") <- centers
  mat_centered
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

#' Simple function to create dummy variables for ordinal and nominal variables
#' 
#' @param data A dataframe containing the ordinal variable. 
#' @param var_ordinal Names of ordinal variables in the data for which dummy variables should be created.
#' @param var_nominal Names of nominal variables in the data for which dummy variables should be created.
#' @param drop_variables Specifies whether to drop the original variables after dummy variables have
#'  been created. The default value is FALSE, and the original variables are kept in the data.
#' @details 
#' This function creates dummy variables based on ordinal and categorical coding described in
#' the Royston and Sauerbrei (2008) book (Chapter 3, Table 3.5). It uses the levels of
#' the categorical variable if they exist; otherwise, it will extract the unique values of the
#' variable, sort them, and use them as levels. We recommend that the user sets the levels of
#' categorical variables and specifies their reference group. You can use the factor() function in
#' base R. If the levels are 1, 2, and 3, then 1 will be the reference group. On the other hand,
#' if the levels are 3, 2, and 1, then 3 will be the reference group. In brief, the first
#' level will be taken as the reference group.
#' 
#' @examples
#' data("gbsg")
#' # create dummy variable for grade using ordinal coding
#' gbsg <- create_dummy_variables(gbsg, var_ordinal = "grade", drop_variables = TRUE)
#' head(gbsg)
#' 
#' @return 
#' A dataframe with new dummy variables.
#' 
#' @export
create_dummy_variables <-  function(data, var_ordinal = NULL, var_nominal = NULL, drop_variables = FALSE) {
    # assert that data must be provided
    if (missing(data))
      stop("! data argument is missing.\n",
        "i An input data.frame is required for the use of create_dummy_ordinal.",
        call. = FALSE
      )
  
    if (!is.data.frame(data))
      stop("The data must be a data.frame")
  
    # colnames of data
    xnames <- colnames(data)
    
    if (is.null(xnames))
      stop("the column names of the provided data are empty")
    
    # assert that either var_ordinal or var_nominal must be provided
    if (is.null(var_nominal) && is.null(var_ordinal))
      stop("Either var_nominal or var_ordinal must be provided.")
    
    # Deal with ordinal variables when provided
    if (!is.null(var_ordinal)) {
      if (!is.character(var_ordinal))
        stop("var_ordinal must be a character name.")
      
      index1 <- which(!var_ordinal %in% xnames)
      
      if (length(index1) != 0)
        stop(paste0("Variable ", var_ordinal[index1]," is not in column name of data",
          collapse = ", "
        ))
      
      if (!is.null(var_nominal)) {
        
        if (!is.character(var_nominal))
          stop("var_nominal must be a character name.")
        
        index2 <- which(!var_nominal %in% xnames)
        
        if (length(index2) != 0)
          stop(paste0(
            "Variable ",
            var_nominal[index2],
            " is not in column name of data",
            collapse = ", "
          ))
      }
      
      for (col in var_ordinal) {
        # levels of the variable if it exist: important for reference
        unique_levels <- levels(data[[col]])
        
        # if levels does not exist use unique values
        if (is.null(unique_levels)) {
          unique_levels <- sort(unique(data[[col]]))
        }
        
        if (length(unique_levels) == 1) {
          warning(paste("var_names ", col,
              "has only one unique level. Skipping dummy variable creation."
            )
          )
          next
        }
    
          levels_list <- lapply(seq_along(unique_levels)[-length(unique_levels)], function(i)
            unique_levels[seq(i)])
        
        for (i in seq_along(levels_list)) {
          level <- levels_list[[i]]
         
          var_name <- paste0(col, "_", paste(i, collapse = "_"))
          data[[var_name]] <- as.integer(!data[[col]] %in% level)
        }
      }
    }
    
    # Deal with nominal variables
    if (!is.null(var_nominal)) {
      dummies <- model.matrix( ~ ., data = data[, var_nominal, drop = FALSE])[, -1]
      data <- cbind(data, dummies)
    }
    
    # drop original variables
    if (drop_variables) {
    vars_to_drop <- c(var_nominal, var_ordinal)
    data <- data[, !(colnames(data) %in% vars_to_drop)]
    }
    return(data)
  }

