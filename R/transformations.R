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
#' powers, see also \code{find_scale_factor()}.
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
#' @param powers passed to \code{fit_acd()}.
#' @param acd_parameter a list usually returned by \code{fit_acd()}. In particular, 
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
  
  # process zero variables before binary detection
  if (zero) {
    # Replace all zero or negative values with 0
    x[x <= 0] <- 0
    shift <- 0
  }
  
  # do not transform true binary variables, but not zero-handled continuous variables
  if (!zero && check_binary && length(unique(x)) <= 2) {
    x <- as.matrix(x)
    if (!is.null(name)) {
      colnames(x) <- name_transformed_variables(name, 1)
    }
    
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
  x_acd <- transform_vector_fp(x = x_acd, power = power[2], scale = 1, 
                               shift = 0, name = name_acd, zero = FALSE)
  x_fp <- transform_vector_fp(x = x, power = power[1], scale = scale, 
                              shift = shift, name = name, zero = zero)

  list(
    acd = cbind(x_fp, x_acd),
    acd_parameter = acd_parameter
  )
}


#' Transform each column of a matrix using final FP powers or ACD transformation
#'
#' This function applies fractional-polynomial (FP) and/or approximate cumulative
#' distribution (ACD) transformations to selected columns of a matrix. It can
#' also generate structural-zero binary indicators for semi-continuous variables
#' (`catzero`) and supports spike-at-zero decisions through `spike` and
#' `spike_decision`.
#'
#' @param x A matrix with continuous variables, usually already shifted and
#'   scaled before transformation. The matrix must have column names.
#' @param power_list A named list of FP or ACD powers to be applied to columns
#'   of `x`. Only variables named in this list are transformed.
#' @param center Named logical vector indicating, for each original variable,
#'   whether centering should be applied to the final transformed columns derived
#'   from that variable. If all entries are `FALSE`, no centering is applied and
#'   the returned `centers` component is `NULL`. If all entries are `TRUE`, the
#'   full transformed design matrix is centered in one call to
#'   `center_matrix()`. If values are mixed, only transformed columns derived
#'   from variables with `center = TRUE` are centered. For variables expanded
#'   into multiple FP or ACD columns, each derived column is centered separately.
#' @param acdx Named logical vector specifying whether each variable should use
#'   the ACD transformation.
#' @param keep_x_order Logical. If `TRUE`, transformed variables are ordered
#'   according to their order in the input matrix `x`. If `FALSE`, columns are
#'   ordered according to `power_list`. The default is `FALSE`, because ordering
#'   by `power_list` reflects the `xorder` argument in `mfp2()`.
#' @param acd_parameter_list A named list of ACD parameters. This is only
#'   required when applying previously estimated transformations to new data,
#'   for example during prediction. Entries should correspond to variables where
#'   `acdx` is `TRUE` and are passed to `transform_vector_acd()`. The default
#'   `NULL` means ACD parameters are estimated during transformation.
#' @param check_binary Passed to `transform_vector_fp()`.
#' @param zero Named logical vector specifying, for each variable, whether
#'   nonpositive values should be treated as structural zero before FP or ACD
#'   transformation. If `TRUE`, transformations are applied to the positive part
#'   and nonpositive values are represented as zero in the transformed columns.
#'   If `FALSE`, all values are transformed directly. If `NULL`, no variables
#'   are treated as zero-handled.
#' @param catzero Named logical vector specifying, for each variable, whether a
#'   structural-zero binary indicator should be added. For these variables, the
#'   indicator is computed as `I(x <= 0)` using the input matrix supplied to this
#'   function. This is consistent with zero handling, where nonpositive values
#'   are treated as structural zero before FP or ACD transformation. If `NULL`,
#'   no categorical-zero indicators are added.
#' @param spike Named logical vector indicating which variables are subject to
#'   spike-at-zero handling. If `NULL`, this is equivalent to setting all entries
#'   to `FALSE`.
#' @param spike_decision Optional named numeric vector with values `1`, `2`, or
#'   `3`, specifying the final spike-at-zero representation for each variable.
#'   Value `1` keeps both the transformed FP/ACD component and the structural-zero
#'   binary indicator. Value `2` keeps only the transformed FP/ACD component and
#'   suppresses the binary indicator. Value `3` keeps only the binary indicator
#'   and removes the FP/ACD component. If `NULL`, no spike-at-zero
#'   post-processing is applied.
#' @param reset_zero Logical. If `TRUE`, variables marked as `zero = TRUE` but
#'   containing only positive values in the current data are reset to `FALSE`,
#'   and their `catzero` indicators are also suppressed. The default is `FALSE`
#'   because final fitting and prediction should preserve the zero, catzero, and
#'   spike structure learned from the fitting process rather than adapt it to a
#'   particular data slice.
#'
#' @details
#' Transformations are applied variable by variable. Variables with
#' `acdx = FALSE` are transformed by `transform_vector_fp()`. Variables with
#' `acdx = TRUE` are transformed by `transform_vector_acd()`. The resulting
#' transformed columns are combined into `x_transformed`.
#'
#' Structural-zero binary indicators are appended after the FP/ACD transformed
#' columns. Binary indicator columns are named with the suffix `"_bin"` and are
#' computed as `I(x <= 0)`.
#'
#' The interaction of `spike`, `spike_decision`, and `catzero` determines whether
#' binary indicators are included:
#'
#' * If `spike = TRUE`:
#'   - `spike_decision = 1`: keep both the FP/ACD component and the binary
#'     indicator.
#'   - `spike_decision = 2`: keep the FP/ACD component only and suppress the
#'     binary indicator, regardless of `catzero`.
#'   - `spike_decision = 3`: keep the binary indicator only and remove the
#'     FP/ACD component.
#'
#' * If `spike = FALSE`:
#'   - binary inclusion depends only on `catzero`.
#'   - the FP/ACD component is included whenever powers are specified.
#'
#' If all entries in `power_list` are `NA`, the function normally returns
#' `NULL`. The exception is a binary-only spike decision
#' (`spike_decision = 3`), where the final selected representation is the
#' structural-zero binary indicator.
#'
#' @section Centering:
#' Centering is controlled at the original-variable level through `center`, but
#' is applied at the final transformed-column level. Internally, the function
#' builds a column-to-variable map from the actual columns returned by
#' `transform_vector_fp()` and `transform_vector_acd()`, and from any appended
#' `"_bin"` columns.
#'
#' The centering rule is:
#'
#' * `center[v] = TRUE`: center all final columns derived from original variable
#'   `v`.
#' * `center[v] = FALSE`: leave all final columns derived from original variable
#'   `v` unchanged.
#'
#' If all variables have `center = TRUE`, the whole transformed design matrix is
#' centered in one call to `center_matrix()`. If only some variables have
#' `center = TRUE`, only the corresponding transformed columns are passed to
#' `center_matrix()`, and all other columns are left unchanged.
#'
#' The actual centering constants are determined by `center_matrix()`:
#'
#' * ordinary continuous transformed columns are centered by their mean;
#' * binary columns are centered by their minimum, usually `0`;
#' * zero-handled FP columns are centered using only their nonzero transformed
#'   values, and transformed zero rows are reset to `0` after centering.
#'
#' For variables transformed with ACD, `transform_vector_acd()` returns the FP
#' component first and the ACD component afterwards. Only the FP component may
#' inherit zero-aware centering from the parent variable. ACD component columns
#' represent transformed cumulative probabilities, not structural-zero FP values;
#' therefore exact numeric zeros in ACD columns are centered as ordinary
#' continuous values rather than excluded from the mean or reset to zero after
#' centering.
#'
#' Binary `"_bin"` columns follow the centering decision of their parent
#' variable. If they are passed to `center_matrix()`, they are treated as binary
#' columns and centered by their minimum. Thus ordinary `0/1` indicators usually
#' remain unchanged.
#'
#' @section Zero-expanded centering:
#' The input `zero` is defined per original variable, while `center_matrix()`
#' needs one zero flag per final transformed column. Therefore this function
#' returns `zero_expanded`, a named logical vector aligned with
#' `colnames(x_transformed)`.
#'
#' For ordinary FP columns derived from variables with `zero = TRUE`, and for
#' the first FP component of an ACD transformation when that component is
#' present, `zero_expanded` is `TRUE`. This causes `center_matrix()` to compute
#' the centering constant from nonzero transformed values and then preserve
#' transformed zero rows as exactly zero.
#'
#' ACD component columns, structurally identified during the ACD transformation
#' branch, are forced to `zero_expanded = FALSE` even when their parent variable
#' has `zero = TRUE`. Structural-zero binary indicator columns ending in
#' `"_bin"` are also always `FALSE`, because these columns are binary indicators
#' rather than zero-handled FP columns.
#'
#' @section Column names:
#' Transformed column names are based on the original variable names. FP terms
#' are suffixed with `".i"` to indicate the transformed-power index. ACD-derived
#' columns are prefixed with `"A_"`. Structural-zero binary indicators are
#' suffixed with `"_bin"`.
#'
#' @examples
#' x <- matrix(1:100, nrow = 10)
#' colnames(x) <- paste0("x", seq_len(ncol(x)))
#'
#' powx <- setNames(
#'   replicate(ncol(x), c(1, 2), simplify = FALSE),
#'   colnames(x)
#' )
#'
#' center <- setNames(rep(FALSE, ncol(x)), colnames(x))
#' acdx <- setNames(rep(FALSE, ncol(x)), colnames(x))
#'
#' transform_matrix(x, powx, center, acdx)
#'
#' @return
#' If all elements of `power_list` are `NA` and there is no binary-only spike
#' decision, returns `NULL`. Otherwise, returns a list with the following
#' components:
#'
#' * `x_transformed`: matrix of transformed variables, optionally centered.
#'   The number of columns may differ from the input matrix because FP and ACD
#'   transformations can create multiple columns per original variable.
#'   Structural-zero binary indicators are appended when selected through
#'   `catzero` or `spike_decision`.
#' * `centers`: named numeric vector of centering constants used for the final
#'   transformed columns. If no centering is requested, this is `NULL`. Under
#'   mixed centering, columns that were not centered have centering constant `0`.
#' * `acd_parameter`: named list of estimated or reused ACD parameters. This may
#'   be empty if no ACD transformation is applied.
#' * `x_trafo`: list of transformed FP/ACD components before adding
#'   structural-zero binary indicators and before centering. For binary-only
#'   spike decisions, the corresponding FP/ACD component is removed from this
#'   list.
#' * `zero_expanded`: named logical vector aligned with the columns of
#'   `x_transformed`, indicating which transformed FP/ACD columns require
#'   zero-specific centering behavior in `center_matrix()`.
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
                             spike = NULL,
                             spike_decision = NULL,
                             reset_zero = FALSE) {
  #-------------------------
  # Input checks
  #-------------------------
  if (!is.matrix(x)) {
    stop("! 'x' must be a matrix.", call. = FALSE)
  }
  
  x_cols <- colnames(x)
  if (is.null(x_cols)) {
    stop("! Input data 'x' must have column names.", call. = FALSE)
  }
  
  if (!is.list(power_list)) {
    stop("! 'power_list' must be a list.", call. = FALSE)
  }
  
  if (is.null(names(power_list))) {
    stop("! 'power_list' must have names.", call. = FALSE)
  }
  
  if (anyDuplicated(names(power_list))) {
    stop("! 'power_list' must not contain duplicated names.", call. = FALSE)
  }
  
  pl_names <- names(power_list)
  
  missing_in_x <- setdiff(pl_names, x_cols)
  if (length(missing_in_x) > 0L) {
    stop(
      "! The following variables in 'power_list' are not in 'x': ",
      paste(missing_in_x, collapse = ", "),
      call. = FALSE
    )
  }
  
  # All per-variable control vectors must refer to exactly the same variables
  # as power_list. This prevents silent NA creation after name-based subsetting.
  check_same_names <- function(x, arg, ref_names) {
    x_names <- names(x)
    
    if (is.null(x_names)) {
      stop(sprintf("! '%s' must have names.", arg), call. = FALSE)
    }
    
    if (anyDuplicated(x_names)) {
      stop(sprintf("! '%s' must not contain duplicated names.", arg), call. = FALSE)
    }
    
    missing_names <- setdiff(ref_names, x_names)
    extra_names <- setdiff(x_names, ref_names)
    
    if (length(missing_names) > 0L || length(extra_names) > 0L) {
      msg <- sprintf("! Names of '%s' must match names of 'power_list'.", arg)
      
      if (length(missing_names) > 0L) {
        msg <- paste0(
          msg,
          "\n",
          sprintf(
            "i Missing from '%s': %s.",
            arg,
            paste0(missing_names, collapse = ", ")
          )
        )
      }
      
      if (length(extra_names) > 0L) {
        msg <- paste0(
          msg,
          "\n",
          sprintf(
            "i Not present in 'power_list': %s.",
            paste0(extra_names, collapse = ", ")
          )
        )
      }
      
      stop(msg, call. = FALSE)
    }
    
    x[ref_names]
  }
  
  # center: per-original-variable centering request.
  if (!is.logical(center)) {
    stop("! 'center' must be a logical vector.", call. = FALSE)
  }
  
  if (anyNA(center)) {
    stop("! 'center' must not contain missing values.", call. = FALSE)
  }
  
  center <- check_same_names(center, "center", pl_names)
  
  # acdx: per-original-variable flag for ACD transformation.
  if (!is.logical(acdx)) {
    stop("! 'acdx' must be a logical vector.", call. = FALSE)
  }
  
  if (anyNA(acdx)) {
    stop("! 'acdx' must not contain missing values.", call. = FALSE)
  }
  
  acdx <- check_same_names(acdx, "acdx", pl_names)
  
  # zero: variables where nonpositive values are structurally handled as zero
  # before FP/ACD transformation.
  if (is.null(zero)) {
    zero <- setNames(rep(FALSE, length(pl_names)), pl_names)
  } else {
    if (!is.logical(zero)) {
      stop("! 'zero' must be a logical vector.", call. = FALSE)
    }
    
    if (anyNA(zero)) {
      stop("! 'zero' must not contain missing values.", call. = FALSE)
    }
    
    zero <- check_same_names(zero, "zero", pl_names)
  }
  
  # catzero: variables requiring an additional structural-zero binary indicator.
  if (is.null(catzero)) {
    catzero <- setNames(rep(FALSE, length(pl_names)), pl_names)
  } else {
    if (!is.logical(catzero)) {
      stop("! 'catzero' must be a logical vector.", call. = FALSE)
    }
    
    if (anyNA(catzero)) {
      stop("! 'catzero' must not contain missing values.", call. = FALSE)
    }
    
    catzero <- check_same_names(catzero, "catzero", pl_names)
  }
  
  # spike: variables assessed by the spike-at-zero selection logic.
  if (is.null(spike)) {
    spike <- setNames(rep(FALSE, length(pl_names)), pl_names)
  } else {
    if (!is.logical(spike)) {
      stop("! 'spike' must be a logical vector.", call. = FALSE)
    }
    
    if (anyNA(spike)) {
      stop("! 'spike' must not contain missing values.", call. = FALSE)
    }
    
    spike <- check_same_names(spike, "spike", pl_names)
  }
  
  # spike_decision encodes the selected spike representation:
  #   1 = FP/ACD part + structural-zero binary
  #   2 = FP/ACD part only
  #   3 = structural-zero binary only
  if (!is.null(spike_decision)) {
    if (!is.numeric(spike_decision)) {
      stop(
        "! 'spike_decision' must be numeric, with values 1, 2, or 3.",
        call. = FALSE
      )
    }
    
    if (anyNA(spike_decision)) {
      stop("! 'spike_decision' must not contain missing values.", call. = FALSE)
    }
    
    bad_vals <- setdiff(unique(spike_decision), c(1, 2, 3))
    if (length(bad_vals) > 0L) {
      stop(
        "! 'spike_decision' must only contain values 1, 2, or 3.",
        call. = FALSE
      )
    }
    
    spike_decision <- check_same_names(
      spike_decision,
      "spike_decision",
      pl_names
    )
  }
  
  # If every FP/ACD power is NA, there is usually no transformed continuous
  # component to return. The exception is a binary-only spike decision, where
  # the final selected representation is the *_bin column.
  all_powers <- unlist(power_list, use.names = FALSE)
  all_powers_na <- length(all_powers) == 0L || all(is.na(all_powers))
  
  has_binary_only_spike <- !is.null(spike_decision) &&
    any(spike & spike_decision == 3, na.rm = TRUE)
  
  if (all_powers_na && !has_binary_only_spike) {
    return(NULL)
  }
  
  #-------------------------
  # Optional zero reset
  #-------------------------
  # If requested, variables marked zero=TRUE but containing only positive values
  # are reset because there is no structural-zero/nonpositive group to preserve.
  if (reset_zero && any(zero)) {
    vars_to_check <- names(zero)[zero]
    
    bad_vars <- vars_to_check[
      apply(
        x[, vars_to_check, drop = FALSE],
        2,
        function(col) all(col > 0, na.rm = TRUE)
      )
    ]
    
    if (length(bad_vars) > 0L) {
      warning(
        "These variables were marked as 'zero = TRUE' but contain only positive values. ",
        "Resetting 'zero' and 'catzero' to FALSE for: ",
        paste(bad_vars, collapse = ", "),
        call. = FALSE
      )
      
      zero[bad_vars] <- FALSE
      catzero[bad_vars] <- FALSE
    }
  }
  
  #----------------------
  # Reorder variables
  #----------------------
  # keep_x_order is used when the output design matrix should follow the column
  # order in x rather than the order in power_list.
  if (keep_x_order) {
    power_list <- power_list[order(match(names(power_list), x_cols))]
  }
  
  names_vars <- names(power_list)
  
  # Reindex all per-variable controls after possible reordering of power_list.
  x <- x[, names_vars, drop = FALSE]
  center <- center[names_vars]
  acdx <- acdx[names_vars]
  zero <- zero[names_vars]
  catzero <- catzero[names_vars]
  spike <- spike[names_vars]
  
  if (!is.null(spike_decision)) {
    spike_decision <- spike_decision[names_vars]
  }
  
  #----------------------
  # FP / ACD transformations
  #----------------------
  x_trafo <- list()
  acd_parameter <- list()
  acd_component_cols <- character(0L)
  
  for (name in names_vars) {
    if (isTRUE(acdx[[name]])) {
      # During fitting, acd_parameter_list is usually NULL and parameters are
      # estimated. During prediction, fitted ACD parameters are supplied here.
      acd_parameter_name <- NULL
      
      if (!is.null(acd_parameter_list)) {
        acd_parameter_name <- acd_parameter_list[[name]]
      }
      
      acd <- transform_vector_acd(
        x[, name],
        power = power_list[[name]],
        acd_parameter = acd_parameter_name,
        name = name,
        zero = zero[[name]]
      )
      
      x_trafo[[name]] <- acd$acd
      acd_parameter[[name]] <- acd$acd_parameter
      
      if (!is.null(acd$acd)) {
        cols_acd <- colnames(acd$acd)
        
        if (is.null(cols_acd)) {
          cols_acd <- name
        }
        
        # transform_vector_acd() returns cbind(x_fp, x_acd): the FP part
        # comes first, followed by ACD component columns. ACD component values
        # are cumulative-probability quantities, not structural zeros, so they
        # must not inherit zero-aware centering from the parent variable.
        n_fp_cols <- if (all(is.na(power_list[[name]][1L]))) 0L else 1L
        
        if (length(cols_acd) > n_fp_cols) {
          acd_component_cols <- c(
            acd_component_cols,
            cols_acd[seq.int(n_fp_cols + 1L, length(cols_acd))]
          )
        }
      }
    } else {
      x_trafo[[name]] <- transform_vector_fp(
        x[, name],
        power = power_list[[name]],
        name = name,
        check_binary = check_binary,
        zero = zero[[name]]
      )
    }
  }
  
  #------------------------------
  # Spike decision post-processing
  #------------------------------
  # The continuous FP/ACD part is built first. Then spike_decision controls
  # whether the structural-zero binary is included and whether the continuous
  # component is retained.
  if (!is.null(spike_decision)) {
    for (v in names(spike_decision)) {
      dec <- spike_decision[[v]]
      
      if (isTRUE(spike[[v]])) {
        if (dec == 2L) {
          # FP/ACD only: suppress *_bin while keeping the transformed component.
          catzero[[v]] <- FALSE
        } else if (dec == 3L) {
          # Binary only: remove the transformed component and force *_bin.
          x_trafo[[v]] <- NULL
          catzero[[v]] <- TRUE
        }
        # dec == 1: FP/ACD + *_bin; no change.
      }
    }
  }
  
  #-------------------------
  # Build transformed matrix
  #-------------------------
  # x_trafo can be empty after binary-only spike decisions. In that case, the
  # matrix may still be created below from catzero/spike binary indicators.
  if (length(x_trafo) > 0L) {
    x_transformed <- do.call(cbind, x_trafo)
  } else {
    x_transformed <- NULL
  }
  
  #----------------------
  # Add catzero / spike binary indicators
  #----------------------
  cat_vars <- names(catzero)[catzero]
  bin_col_to_var <- setNames(character(0L), character(0L))
  
  if (length(cat_vars) > 0L) {
    catzero_list <- lapply(cat_vars, function(v) {
      binary_only_spike <- isTRUE(spike[[v]]) &&
        !is.null(spike_decision) &&
        v %in% names(spike_decision) &&
        isTRUE(spike_decision[[v]] == 3L)
      
      # all powers NA usually means the variable was eliminated. The exception
      # is spike_decision == 3, where the binary indicator is the selected term.
      if (all(is.na(power_list[[v]])) && !binary_only_spike) {
        return(NULL)
      }
      
      # Structural-zero indicator is I(x <= 0), consistent with
      # transform_vector_fp(..., zero = TRUE), which treats nonpositive values
      # as structural zero before FP transformation.
      as.integer(x[, v] <= 0)
    })
    
    valid_idx <- !vapply(catzero_list, is.null, logical(1L))
    catzero_list <- catzero_list[valid_idx]
    cat_vars_valid <- cat_vars[valid_idx]
    names(catzero_list) <- cat_vars_valid
    
    if (length(catzero_list) > 0L) {
      catzero_matrix <- do.call(cbind, catzero_list)
      colnames(catzero_matrix) <- paste0(cat_vars_valid, "_bin")
      bin_col_to_var <- c(
        bin_col_to_var,
        setNames(cat_vars_valid, colnames(catzero_matrix))
      )
      
      if (is.null(x_transformed)) {
        x_transformed <- catzero_matrix
      } else {
        x_transformed <- cbind(x_transformed, catzero_matrix)
      }
    }
  }
  
  # No FP/ACD columns and no binary indicators remain.
  if (is.null(x_transformed) || ncol(x_transformed) == 0L) {
    return(NULL)
  }
  
  #----------------------
  # Column-to-variable map
  #----------------------
  # After transformation, original variables may expand into several columns:
  #   x      -> x.1, x.2, ...
  #   acd(x) -> x.1 and/or A_x.1
  #   catzero/spike -> x_bin
  #
  # col_to_var maps every final transformed column back to its original variable.
  # This map is used for both zero_expanded and mixed per-variable centering.
  col_to_var <- setNames(
    rep(NA_character_, ncol(x_transformed)),
    colnames(x_transformed)
  )
  
  for (v in names(x_trafo)) {
    if (!is.null(x_trafo[[v]])) {
      cols_v <- colnames(x_trafo[[v]])
      
      if (is.null(cols_v)) {
        cols_v <- v
      }
      
      cols_v <- intersect(cols_v, colnames(x_transformed))
      
      if (length(cols_v) > 0L) {
        col_to_var[cols_v] <- v
      }
    }
  }
  
  bin_cols <- intersect(names(bin_col_to_var), colnames(x_transformed))
  
  if (length(bin_cols) > 0L) {
    col_to_var[bin_cols] <- bin_col_to_var[bin_cols]
  }
  
  #----------------------
  # Expanded zero vector
  #----------------------
  # zero is defined per original variable, but center_matrix() needs one zero
  # flag per transformed column. zero_expanded is therefore the transformed
  # column-level version of zero.
  #
  # For zero-handled FP columns, including the FP component of an ACD
  # transformation, zero_expanded = TRUE makes center_matrix() compute the
  # center using only nonzero transformed values and then reset transformed
  # zero rows back to 0 after centering.
  #
  # *_bin columns are not FP/ACD transformed zero columns. They are binary
  # indicators, so they must remain zero_expanded = FALSE and are handled by
  # center_matrix()'s binary-minimum rule.
  zero_expanded <- setNames(
    rep(FALSE, ncol(x_transformed)),
    colnames(x_transformed)
  )
  
  mapped_zero_cols <- !is.na(col_to_var) & col_to_var %in% names(zero)
  
  if (any(mapped_zero_cols)) {
    zero_expanded[mapped_zero_cols] <- zero[col_to_var[mapped_zero_cols]]
  }
  
  # ACD component columns are not structural-zero FP columns. Even when the
  # parent variable has zero=TRUE, exact numeric zeros in ACD columns can be
  # legitimate transformed cumulative probabilities and must be centered as
  # ordinary continuous values.
  acd_component_cols <- intersect(acd_component_cols, names(zero_expanded))
  
  if (length(acd_component_cols) > 0L) {
    zero_expanded[acd_component_cols] <- FALSE
  }
  
  if (length(bin_cols) > 0L) {
    zero_expanded[bin_cols] <- FALSE
  }
  
  #----------------------
  # Centering
  #----------------------
  # center is defined per original variable. Mixed centering means:
  #   center[v] = TRUE  -> center all final columns derived from v
  #   center[v] = FALSE -> leave all final columns derived from v unchanged
  #
  # Binary *_bin columns follow their original variable. If centered, they are
  # passed to center_matrix(); center_matrix() centers binary columns by their
  # minimum, so usual 0/1 indicators remain unchanged.
  centers <- NULL
  
  if (any(center)) {
    centers <- setNames(rep(0, ncol(x_transformed)), colnames(x_transformed))
    
    if (all(center)) {
      # Fast path: all variables request centering, so center the whole design
      # matrix in one call.
      x_transformed <- center_matrix(
        mat = x_transformed,
        centers = NULL,
        zero = zero_expanded
      )
      
      centers <- attr(x_transformed, "scaled:center")
    } else {
      # Mixed path: center only columns whose original variable has center=TRUE.
      mapped_center_cols <- !is.na(col_to_var) & col_to_var %in% names(center)
      
      center_flags <- setNames(
        rep(FALSE, length(col_to_var)),
        names(col_to_var)
      )
      
      if (any(mapped_center_cols)) {
        center_flags[mapped_center_cols] <- center[col_to_var[mapped_center_cols]]
      }
      
      cols_to_center <- names(center_flags)[center_flags]
      cols_to_center <- intersect(cols_to_center, colnames(x_transformed))
      
      if (length(cols_to_center) > 0L) {
        centered_part <- center_matrix(
          mat = x_transformed[, cols_to_center, drop = FALSE],
          centers = NULL,
          zero = zero_expanded[cols_to_center]
        )
        
        x_transformed[, cols_to_center] <- centered_part
        centers[cols_to_center] <- attr(centered_part, "scaled:center")
      }
    }
  }
  
  list(
    x_transformed = x_transformed,
    centers = centers,
    acd_parameter = acd_parameter,
    x_trafo = x_trafo,
    zero_expanded = zero_expanded
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
#' @keywords internal
#' @noRd
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
  
  # Validate centers (if provided)
  if (!is.null(centers)) {
    if (!is.numeric(centers)) {
      stop("! 'centers' must be numeric.")
    }
    if (is.null(names(centers))) {
      stop("! 'centers' must have names.")
    }
    if (!setequal(names(centers), colnames(mat))) {
      stop("! 'centers' names must match column names of 'mat'.")
    }
    centers <- centers[colnames(mat)]  # reorder to match mat
  }
  
  # Compute centers if not provided
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
    
    centers <- setNames(centers, colnames(mat))
  }
  
  # Apply centering
  mat_centered <- scale(mat, center = centers, scale = FALSE)
  
  # Reset zero values for variables flagged with zero=TRUE
  for (j in seq_len(ncol(mat))) {
    if (zero[j]) {
      mat_centered[mat[, j] == 0, j] <- 0
    }
  }
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
#' @keywords internal
#' @noRd
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
create_dummy_variables <- function(data, 
                                   var_ordinal = NULL, 
                                   var_nominal = NULL, 
                                   drop_variables = FALSE) {
  
  # assert that data must be provided
  if (missing(data)) {
    stop(
      "! data argument is missing.\n",
      "i An input data.frame is required for the use of create_dummy_variables.",
      call. = FALSE
    )
  }
  
  if (!is.data.frame(data)) {
    stop("The data must be a data.frame.", call. = FALSE)
  }
  
  # colnames of data
  xnames <- colnames(data)
  
  if (is.null(xnames)) {
    stop("The column names of the provided data are empty.", call. = FALSE)
  }
  
  # assert that either var_ordinal or var_nominal must be provided
  if (is.null(var_nominal) && is.null(var_ordinal)) {
    stop("Either var_nominal or var_ordinal must be provided.", call. = FALSE)
  }
  
  # validate ordinal variable names
  if (!is.null(var_ordinal)) {
    if (!is.character(var_ordinal)) {
      stop("var_ordinal must be a character vector.", call. = FALSE)
    }
    
    index1 <- which(!var_ordinal %in% xnames)
    
    if (length(index1) != 0L) {
      stop(
        paste0(
          "Variable ",
          var_ordinal[index1],
          " is not in column names of data.",
          collapse = ", "
        ),
        call. = FALSE
      )
    }
  }
  
  # validate nominal variable names
  if (!is.null(var_nominal)) {
    if (!is.character(var_nominal)) {
      stop("var_nominal must be a character vector.", call. = FALSE)
    }
    
    index2 <- which(!var_nominal %in% xnames)
    
    if (length(index2) != 0L) {
      stop(
        paste0(
          "Variable ",
          var_nominal[index2],
          " is not in column names of data.",
          collapse = ", "
        ),
        call. = FALSE
      )
    }
  }
  
  # A variable cannot be encoded both as ordinal and nominal.
  # These encodings are mutually exclusive and would create duplicated or
  # contradictory dummy variables for the same source variable.
  if (!is.null(var_ordinal) && !is.null(var_nominal)) {
    overlap <- intersect(var_ordinal, var_nominal)
    
    if (length(overlap) > 0L) {
      stop(
        sprintf(
          "Variables cannot be both ordinal and nominal: %s.",
          paste(overlap, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }
  
  # Deal with ordinal variables when provided
  if (!is.null(var_ordinal)) {
    for (col in var_ordinal) {
      
      # levels of the variable if it exists; important for reference
      unique_levels <- levels(data[[col]])
      
      # if levels do not exist, use sorted unique values
      if (is.null(unique_levels)) {
        unique_levels <- sort(unique(data[[col]]))
      }
      
      if (length(unique_levels) == 1L) {
        warning(
          paste(
            "var_names",
            col,
            "has only one unique level. Skipping dummy variable creation."
          ),
          call. = FALSE
        )
        next
      }
      
      levels_list <- lapply(
        seq_along(unique_levels)[-length(unique_levels)],
        function(i) unique_levels[seq_len(i)]
      )
      
      for (i in seq_along(levels_list)) {
        level <- levels_list[[i]]
        
        var_name <- paste0(col, "_", paste(i, collapse = "_"))
        data[[var_name]] <- as.integer(!data[[col]] %in% level)
      }
    }
  }
  
  # Deal with nominal variables when provided
  if (!is.null(var_nominal)) {
    dummies <- stats::model.matrix(
      ~ .,
      data = data[, var_nominal, drop = FALSE]
    )[, -1, drop = FALSE]
    
    data <- cbind(data, dummies)
  }
  
  # drop original variables
  if (drop_variables) {
    vars_to_drop <- c(var_nominal, var_ordinal)
    data <- data[, !(colnames(data) %in% vars_to_drop), drop = FALSE]
  }
  
  data
}
#' Cumulative (Threshold) Contrast Coding for Ordered Factors
#'
#' @description
#' Generates a contrast matrix for ordered factors in which each column
#' represents a cumulative comparison: level k versus all lower levels.
#' This is sometimes called cumulative or sequential dummy coding.
#'
#' @param n Integer or character vector. If an integer, specifies the number
#'   of levels. If a character vector, its elements are used as factor levels.
#'
#' @return A numeric matrix with \code{length(n)} rows and \code{length(n)-1} columns.
#'   Each column is a dummy variable encoding the cumulative comparison of
#'   higher levels against lower levels. Row names correspond to factor levels.
#'
#' @details
#' For an ordered factor with levels A < B < C < D, the resulting matrix
#' produces three columns:
#' Column 1 compares B, C, D versus A; 
#' Column 2 compares C, D versus A and B; 
#' Column 3 compares D versus A, B, and C.
#' Each element of the matrix is 0 or 1, with 1 indicating that the observation
#' belongs to the "higher" category for that threshold.
#'
#' **Column Names:** Column names are formatted as \code{varname_1}, \code{varname_2}, etc.
#' These names are syntactically valid, unique, and correspond to thresholds in order:
#' \code{varname_1} represents the first threshold (level 2 vs level 1),
#' \code{varname_2} represents the second threshold (level 3 vs levels 1 and 2), etc.
#'
#' This contrast matrix can be assigned to an ordered factor via
#' \code{contrasts()} before fitting a model, e.g., in \code{mfp2.formula()}.
#' This approach preserves ordinal information while allowing threshold-type
#' interpretation of regression coefficients.
#'
#' @examples
#' # Create a data frame
#' data <- data.frame(grade = c("A", "B", "C", "D", "A"))
#' # Convert the column to an ordered factor
#' data$grade <- factor(data$grade, levels = c("A", "B", "C", "D"), ordered = TRUE)
#' # Assign the cumulative contrasts to the ordered factor
#' contrasts(data$grade) <- contr.cumulative(levels(data$grade))
#'
#' @export
contr.cumulative <- function(n) {
  if (is.numeric(n)) nlev <- n else nlev <- length(n)
  mat <- matrix(0, nrow = nlev, ncol = nlev - 1)
  for (j in 1:(nlev - 1)) {
    mat[(j + 1):nlev, j] <- 1
  }
  rownames(mat) <- if (is.numeric(n)) as.character(seq_len(n)) else n
  # Safe column names without special characters
  colnames(mat) <- paste0("_", seq_len(nlev - 1))
  mat
}

