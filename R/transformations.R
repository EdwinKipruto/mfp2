#' Transform a Variable Using Fractional-Polynomial Powers
#'
#' Constructs fractional-polynomial design columns from specified powers,
#' analogous to `fracgen` in Stata.
#'
#' @details
#' This is a low-level utility for constructing fractional-polynomial design
#' columns from specified powers. It does not perform model or power selection.
#' Most users should specify continuous variables with [fp()] inside an
#' [mfp2()] formula instead.
#'
#' The FP transformation constructs design columns as follows. For each power
#' pi in `power` = (p1, p2, ..., pn) it computes x^pi and returns the
#' collection as a matrix. The input may be shifted and scaled before
#' transformation. Centering, if desired, should be applied separately
#' after transformation (see the Data processing section).
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
#' Variables are shifted and then scaled before any power transformation is
#' applied. Shifting ensures positive values; scaling avoids numerically
#' extreme magnitudes. Scaling does not change the selected powers.
#'
#' Centering is not performed by this function. When centering is desired, it
#' should be applied after transformation, because centering before
#' transformation alters the correlation structure among predictors. In
#' [mfp2()], centering is applied only to the final selected model when
#' `center = TRUE`. See Sauerbrei et al (2006) for discussion.
#'
#' If `zero = TRUE`, `x` must be nonnegative. Exact-zero values are not shifted
#' and their transformed components remain zero; transformation is applied only
#' to values with `x > 0`. For performance, this low-level helper assumes that
#' callers have already enforced the nonnegative-input contract. The public
#' fitting and prediction interfaces perform that validation once, before any
#' transformation work.
#'
#' @section Domain restrictions:
#' Fractional-polynomial transformations are applied after shifting and scaling.
#' Users should ensure that the shifted and scaled values are in the domain
#' required by the specified powers. In particular, powers involving logarithms
#' require strictly positive input. This includes `power = 0`, which represents
#' `log(x)`, and repeated powers such as `c(2, 2)`, where the second transformed
#' variable is proportional to `x^2 * log(x)`.
#'
#' Negative powers are undefined at zero, and non-integer powers are not
#' real-valued for negative input. Therefore, when using powers such as `0`,
#' `-1`, `-0.5`, `0.5`, or repeated powers, the values of `(x + shift) / scale`
#' should be strictly positive unless `zero = TRUE`.
#'
#' Positive integer powers such as `1`, `2`, and `3` can be evaluated at zero
#' and negative values, provided they are not part of a repeated-power FP basis.
#'
#' These low-level transformation functions assume that the supplied `x`,
#' `shift`, `scale`, and `power` values are compatible. If incompatible values
#' are supplied, the transformation may produce non-finite values or fail with
#' an error from the underlying numerical code. In normal use through `mfp2()`,
#' shift and scale parameters are estimated and checked during model fitting.
#' During prediction, `predict.mfp2()` checks the fitted powers and stops with
#' an informative error if new data fall outside the required transformation
#' domain.
#'
#' If `zero = TRUE`, exact-zero values are treated structurally: they are not
#' shifted, and their transformed continuous components are set to zero.
#'
#' @param x a vector of a predictor variable.
#' @param power a numeric vector indicating the FP power. Default is 1 (linear).
#' Must be a vector of length 2 for acd transformation. Ignores `NA`, unless
#' an ACD transformation is applied in which case power must be a numeric
#' vector of length 2, and `NA` indicated which parts are used for the final
#' FP.
#' @param scale Positive numeric scaling factor applied to `x` after shifting.
#'   Default `1` (no scaling). If `NULL`, the scale factor is estimated
#'   automatically. Scaling does not change the selected powers.
#' @param shift shift required for shifting x to positive values. Default is 0,
#' meaning no shift is applied. If `NULL` then the shift is estimated
#' automatically using the Royston and Sauerbrei formula iff any `x` <= 0.
#' @param name character used to define names for the output matrix. Default
#' is `NULL`, meaning the output will have unnamed columns.
#' @param zero Logical indicating whether only values with `x > 0` should be
#' transformed, with exact-zero values retained as structural zeros. If
#' `TRUE`, `x` must be nonnegative and must have been validated by the caller;
#' negative values must be explicitly recoded if they should represent the zero
#' group. If `FALSE`
#' (default), all values are shifted (if needed) to ensure positivity.
#' @param check_binary a logical indicating whether or not input `x` is checked
#' if it is a binary variable (i.e. has only two distinct values). The default
#' `TRUE` usually only needs to changed when this function is to be used to
#' transform data for predictions. See Details.
#'
#' @examples
#' z = 1:10
#' transform_vector_fp(z)
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

  if (!is.logical(zero) || length(zero) != 1L || is.na(zero)) {
    stop("`zero` must be a single logical value (TRUE or FALSE).", call. = FALSE)
  }

  if (!is.logical(check_binary) ||
      length(check_binary) != 1L ||
      is.na(check_binary)) {
    stop("`check_binary` must be a single logical value (TRUE or FALSE).",
         call. = FALSE)
  }

  if (all(is.na(power))) {
    return(NULL)
  }

  if (is.null(shift)) {
    shift <- find_shift_factor(x)
  }

  if (is.null(scale)) {
    scale <- find_scale_factor(x)
  }

  if (zero) {
    shift <- 0
  }

  power <- sort(power)

  if (!zero && check_binary && length(unique(x)) <= 2L) {
    x_out <- matrix(x, ncol = 1L)
    if (!is.null(name)) {
      colnames(x_out) <- name_transformed_variables(name, 1L)
    }
    return(x_out)
  }

  x_trafo <- transform_fp_core(
    x_raw     = as.numeric(x),
    power     = as.numeric(power),
    shift_val = as.numeric(shift),
    scale_val = as.numeric(scale),
    zero      = zero
  )

  if (!is.null(name)) {
    colnames(x_trafo) <- name_transformed_variables(name, ncol(x_trafo))
  }

  x_trafo
}

#' @describeIn transform_vector_fp Function to generate acd transformation.
#' @param acd_parameter a list usually returned by \code{fit_acd()}. In particular,
#' it must have components that define `beta0`, `beta1`, `power`, `shift` and
#' `scale` which are to be applied when using the acd transformation in
#' new data.
#' @param powers passed to \code{fit_acd()}.
#' @keywords internal
#' @noRd
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
    x_acd <- do.call(
      apply_acd,
      utils::modifyList(acd_parameter, list(x = x, zero = zero))
    )
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
#' @param zero Named logical vector specifying, for each nonnegative variable,
#'   whether exact-zero values should be treated as structural zero before FP or
#'   ACD transformation. If `TRUE`, transformations are applied to the positive
#'   part (`x > 0`) and exact-zero values remain zero in transformed columns.
#'   Callers must supply values already checked against the nonnegative-input
#'   contract and explicitly recode negative values if appropriate.
#'   If `FALSE`, all values are transformed directly. If `NULL`, no variables
#'   are treated as zero-handled.
#' @param catzero Named logical vector specifying, for each variable, whether a
#'   structural-zero binary indicator should be added. For these variables, the
#'   indicator is computed as `I(x == 0)` using the input matrix supplied to this
#'   function. This is consistent with exact-zero handling. If `NULL`,
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
#' computed as `I(x == 0)`.
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
#' * `transformed_column_to_source`: named character vector mapping each final
#'   transformed column to the source column from which it was constructed.
#' * `transformed_column_component`: named character vector identifying each
#'   final transformed column as an FP basis, ACD basis, structural-zero
#'   indicator, or unchanged binary column.
#' * `transformed_column_zero_handled`: named logical vector indicating which
#'   final columns use positive-part centering, where exact-zero rows remain 0.
#' * `transformed_column_centered`: named logical vector indicating whether
#'   centering was requested for the source variable of each final column.
#' @keywords internal
#' @noRd
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

  # ---------------------------------------------------------------------------
  # Input checks
  # ---------------------------------------------------------------------------
  validated <- validate_transform_matrix_args(
    x = x, power_list = power_list, center = center, acdx = acdx,
    zero = zero, catzero = catzero, spike = spike,
    spike_decision = spike_decision
  )
  pl_names       <- validated$pl_names
  center         <- validated$center
  acdx           <- validated$acdx
  zero           <- validated$zero
  catzero        <- validated$catzero
  spike          <- validated$spike
  spike_decision <- validated$spike_decision

  # If every FP/ACD power is NA, there is usually no transformed continuous
  # component to return. The exception is a binary-only spike decision, where
  # the final selected representation is the *_bin column.
  all_powers <- unlist(power_list, use.names = FALSE)
  all_powers_na <- length(all_powers) == 0L || all(is.na(all_powers))

  has_binary_only_spike <- !is.null(spike_decision) &&
    any(spike & spike_decision == saz_decision_codes[["binary_only"]],
        na.rm = TRUE)

  if (all_powers_na && !has_binary_only_spike) {
    return(NULL)
  }

  # ---------------------------------------------------------------------------
  # Optional zero reset
  # ---------------------------------------------------------------------------
  if (reset_zero && any(zero)) {
    reset <- reset_ineligible_zero_flags(x = x, zero = zero, catzero = catzero)
    zero    <- reset$zero
    catzero <- reset$catzero
  }

  # ---------------------------------------------------------------------------
  # Reorder variables
  # ---------------------------------------------------------------------------
  reordered <- reorder_transform_inputs(
    x = x, power_list = power_list, center = center, acdx = acdx,
    zero = zero, catzero = catzero, spike = spike,
    spike_decision = spike_decision, keep_x_order = keep_x_order
  )
  x              <- reordered$x
  power_list     <- reordered$power_list
  center         <- reordered$center
  acdx           <- reordered$acdx
  zero           <- reordered$zero
  catzero        <- reordered$catzero
  spike          <- reordered$spike
  spike_decision <- reordered$spike_decision
  names_vars     <- reordered$names_vars

  # ---------------------------------------------------------------------------
  # Resolve binary-only SAZ variables before constructing FP/ACD components
  # ---------------------------------------------------------------------------
  # This must happen before the transformation loop. Binary-only SAZ variables
  # must not have their positive-value FP/ACD function constructed and
  # discarded afterwards, because that can evaluate log(0), negative powers, or
  # repeated-power log terms at structural-zero rows.
  resolved_saz <- identify_binary_only_spike_vars(
    catzero = catzero, spike = spike, spike_decision = spike_decision
  )
  catzero                 <- resolved_saz$catzero
  binary_only_spike_vars  <- resolved_saz$binary_only_spike_vars

  # ---------------------------------------------------------------------------
  # FP / ACD transformations
  # ---------------------------------------------------------------------------
  fp_acd <- build_fp_acd_columns(
    x = x, names_vars = names_vars, power_list = power_list, acdx = acdx,
    zero = zero, acd_parameter_list = acd_parameter_list,
    check_binary = check_binary, binary_only_spike_vars = binary_only_spike_vars
  )
  x_trafo            <- fp_acd$x_trafo
  acd_parameter      <- fp_acd$acd_parameter
  acd_component_cols <- fp_acd$acd_component_cols

  # ---------------------------------------------------------------------------
  # Spike decision post-processing
  # ---------------------------------------------------------------------------
  spike_applied <- apply_spike_decision_to_columns(
    x_trafo = x_trafo, catzero = catzero, spike = spike,
    spike_decision = spike_decision
  )
  x_trafo <- spike_applied$x_trafo
  catzero <- spike_applied$catzero

  # ---------------------------------------------------------------------------
  # Build transformed matrix, then add catzero / spike binary indicators
  # ---------------------------------------------------------------------------
  bin <- build_binary_indicator_columns(
    x = x, x_trafo = x_trafo, catzero = catzero, spike = spike,
    spike_decision = spike_decision, power_list = power_list
  )
  x_transformed  <- bin$x_transformed
  bin_col_to_var <- bin$bin_col_to_var

  # No FP/ACD columns and no binary indicators remain.
  if (is.null(x_transformed) || ncol(x_transformed) == 0L) {
    return(NULL)
  }

  # ---------------------------------------------------------------------------
  # Column-to-variable map
  # ---------------------------------------------------------------------------
  col_map    <- build_col_to_var_map(
    x_trafo = x_trafo, x_transformed = x_transformed,
    bin_col_to_var = bin_col_to_var
  )
  col_to_var <- col_map$col_to_var
  bin_cols   <- col_map$bin_cols

  # ---------------------------------------------------------------------------
  # Expanded zero vector
  # ---------------------------------------------------------------------------
  zero_expanded <- expand_zero_flags_to_columns(
    x_transformed = x_transformed, col_to_var = col_to_var, zero = zero,
    acd_component_cols = acd_component_cols, bin_cols = bin_cols
  )

  # ---------------------------------------------------------------------------
  # Centering
  # ---------------------------------------------------------------------------
  centered <- apply_mixed_centering(
    x_transformed = x_transformed, center = center, col_to_var = col_to_var,
    zero_expanded = zero_expanded
  )

  # Record the exact structural metadata used to assemble the final design.
  # Display code can then describe each fitted column without parsing generated
  # names or reconstructing how the column was created.
  column_component <- stats::setNames(
    rep("fp_basis", ncol(x_transformed)),
    colnames(x_transformed)
  )
  column_component[intersect(acd_component_cols, names(column_component))] <-
    "acd_basis"
  column_component[intersect(bin_cols, names(column_component))] <-
    "zero_indicator"

  if (isTRUE(check_binary)) {
    binary_sources <- names_vars[vapply(names_vars, function(v) {
      !isTRUE(zero[[v]]) &&
        !isTRUE(acdx[[v]]) &&
        length(unique(x[, v])) <= 2L
    }, logical(1L))]

    identity_columns <- names(col_to_var)[
      !is.na(col_to_var) &
        col_to_var %in% binary_sources &
        column_component == "fp_basis"
    ]
    column_component[identity_columns] <- "identity_binary"
  }

  column_centered <- stats::setNames(
    rep(FALSE, ncol(x_transformed)),
    colnames(x_transformed)
  )
  mapped_center <- !is.na(col_to_var) & col_to_var %in% names(center)
  if (any(mapped_center)) {
    column_centered[mapped_center] <- center[col_to_var[mapped_center]]
  }

  list(
    x_transformed                    = centered$x_transformed,
    centers                          = centered$centers,
    acd_parameter                    = acd_parameter,
    x_trafo                          = x_trafo,
    zero_expanded                    = zero_expanded,
    transformed_column_to_source     = col_to_var,
    transformed_column_component     = column_component,
    transformed_column_zero_handled  = zero_expanded,
    transformed_column_centered      = column_centered
  )
}

# -----------------------------------------------------------------------------
# transform_matrix() helpers --------------------------------------------------
# -----------------------------------------------------------------------------

#' Check That a Named Vector's Names Match a Reference Set
#'
#' Internal helper used by \code{validate_transform_matrix_args()}. Validates
#' that \code{x} has names, that its names are not duplicated, and that its
#' names are exactly \code{ref_names} (no missing, no extra). Returns
#' \code{x} reindexed to \code{ref_names} order.
#'
#' Promoted from a closure previously defined locally inside
#' \code{transform_matrix()} (it never actually captured anything from that
#' enclosing scope), so it can be reused and tested independently of
#' \code{transform_matrix()}.
#'
#' @param x Named vector to check.
#' @param arg Character scalar; the argument name, used only in error
#'   messages.
#' @param ref_names Character vector of expected names, in the order the
#'   result should follow.
#'
#' @return \code{x[ref_names]}.
#' @keywords internal
#' @noRd
check_same_names_as <- function(x, arg, ref_names) {
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

#' Validate and Normalize transform_matrix() Input Arguments
#'
#' Internal helper used by \code{transform_matrix()}. Validates \code{x} and
#' \code{power_list}, then validates and normalizes the six per-variable
#' control vectors (\code{center}, \code{acdx}, \code{zero}, \code{catzero},
#' \code{spike}, \code{spike_decision}) so that every non-\code{NULL} one is a
#' complete, non-missing, \code{power_list}-ordered vector. \code{NULL}
#' optional arguments (\code{zero}, \code{catzero}, \code{spike}) are filled
#' with all-\code{FALSE} defaults; \code{spike_decision} stays \code{NULL} if
#' not supplied.
#'
#' This was previously the first ~190 lines of \code{transform_matrix()}'s
#' body. Extracting it isolates argument validation from the actual
#' transformation logic, so validation edge cases (missing names, duplicated
#' names, out-of-range \code{spike_decision} values, etc.) can be tested
#' directly without constructing a full transformation scenario.
#'
#' @param x Numeric matrix with column names.
#' @param power_list Named list of FP/ACD powers, one element per variable.
#' @param center,acdx Named logical vectors, one entry per \code{power_list}
#'   variable.
#' @param zero,catzero,spike Optional named logical vectors, one entry per
#'   \code{power_list} variable. \code{NULL} is treated as all-\code{FALSE}.
#' @param spike_decision Optional named numeric vector with values in
#'   \code{unname(saz_decision_codes)}, one entry per \code{power_list}
#'   variable.
#'
#' @return A list with \code{pl_names} (character vector, the validated
#'   variable order from \code{power_list}) and the normalized
#'   \code{center}, \code{acdx}, \code{zero}, \code{catzero}, \code{spike},
#'   \code{spike_decision}, each reindexed to \code{pl_names} order
#'   (\code{spike_decision} is \code{NULL} if not supplied).
#'
#' @keywords internal
#' @noRd
validate_transform_matrix_args <- function(x, power_list, center, acdx,
                                           zero, catzero, spike,
                                           spike_decision) {
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

  # center: per-original-variable centering request.
  if (!is.logical(center)) {
    stop("! 'center' must be a logical vector.", call. = FALSE)
  }

  if (anyNA(center)) {
    stop("! 'center' must not contain missing values.", call. = FALSE)
  }

  center <- check_same_names_as(center, "center", pl_names)

  # acdx: per-original-variable flag for ACD transformation.
  if (!is.logical(acdx)) {
    stop("! 'acdx' must be a logical vector.", call. = FALSE)
  }

  if (anyNA(acdx)) {
    stop("! 'acdx' must not contain missing values.", call. = FALSE)
  }

  acdx <- check_same_names_as(acdx, "acdx", pl_names)

  # zero: variables where exact-zero values are structurally handled as zero
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

    zero <- check_same_names_as(zero, "zero", pl_names)
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

    catzero <- check_same_names_as(catzero, "catzero", pl_names)
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

    spike <- check_same_names_as(spike, "spike", pl_names)
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

    bad_vals <- setdiff(unique(spike_decision), unname(saz_decision_codes))
    if (length(bad_vals) > 0L) {
      stop(
        "! 'spike_decision' must only contain values 1L, 2L, or 3L.",
        call. = FALSE
      )
    }

    spike_decision <- check_same_names_as(
      spike_decision,
      "spike_decision",
      pl_names
    )
  }

  list(
    pl_names       = pl_names,
    center         = center,
    acdx           = acdx,
    zero           = zero,
    catzero        = catzero,
    spike          = spike,
    spike_decision = spike_decision
  )
}

#' Reset Zero Flags for Variables With No Structural-Zero Mass
#'
#' Internal helper used by \code{transform_matrix()} when
#' \code{reset_zero = TRUE}. Variables marked \code{zero = TRUE} but whose
#' raw column in \code{x} contains only strictly positive values have no
#' structural exact-zero group to preserve, so both \code{zero} and
#' \code{catzero} are reset to \code{FALSE} for them (with a warning).
#'
#' @param x Numeric matrix with column names.
#' @param zero,catzero Named logical vectors, one entry per variable in
#'   \code{x}.
#'
#' @return A list with the (possibly updated) \code{zero} and \code{catzero}.
#' @keywords internal
#' @noRd
reset_ineligible_zero_flags <- function(x, zero, catzero) {
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

  list(zero = zero, catzero = catzero)
}

#' Reorder x and Per-Variable Control Vectors to a Common Variable Order
#'
#' Internal helper used by \code{transform_matrix()}. When
#' \code{keep_x_order = TRUE}, reorders \code{power_list} to follow the
#' column order of \code{x} rather than its own original order. Either way,
#' reindexes \code{x} and every per-variable control vector
#' (\code{center}, \code{acdx}, \code{zero}, \code{catzero}, \code{spike},
#' \code{spike_decision}) to the resulting variable order, so downstream code
#' can rely on all of them being aligned.
#'
#' @param x Numeric matrix with column names.
#' @param power_list Named list of FP/ACD powers, one element per variable.
#' @param center,acdx,zero,catzero,spike Named vectors aligned with
#'   \code{power_list}.
#' @param spike_decision Named numeric vector aligned with \code{power_list},
#'   or \code{NULL}.
#' @param keep_x_order Logical. If \code{TRUE}, follow \code{colnames(x)}
#'   order instead of \code{power_list}'s own order.
#'
#' @return A list with the reordered \code{x}, \code{power_list}, \code{center},
#'   \code{acdx}, \code{zero}, \code{catzero}, \code{spike},
#'   \code{spike_decision}, and \code{names_vars} (the resulting variable
#'   order, character vector).
#' @keywords internal
#' @noRd
reorder_transform_inputs <- function(x, power_list, center, acdx, zero,
                                     catzero, spike, spike_decision,
                                     keep_x_order) {
  # keep_x_order is used when the output design matrix should follow the column
  # order in x rather than the order in power_list.
  if (keep_x_order) {
    power_list <- power_list[order(match(names(power_list), colnames(x)))]
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

  list(
    x              = x,
    power_list     = power_list,
    center         = center,
    acdx           = acdx,
    zero           = zero,
    catzero        = catzero,
    spike          = spike,
    spike_decision = spike_decision,
    names_vars     = names_vars
  )
}

#' Identify Binary-Only SAZ Variables Before Building FP/ACD Components
#'
#' Internal helper used by \code{transform_matrix()}. Must run before the
#' FP/ACD transformation loop. For every spike variable whose
#' \code{spike_decision} is \code{"continuous_only"}, suppresses the
#' \code{catzero} binary indicator (the continuous component is kept). For
#' every spike variable whose \code{spike_decision} is \code{"binary_only"},
#' forces \code{catzero} on and records the variable so that
#' \code{build_fp_acd_columns()} can skip constructing its positive-value
#' FP/ACD function entirely.
#'
#' This early resolution matters for correctness, not just efficiency:
#' binary-only SAZ variables must not have their positive-value FP/ACD
#' function constructed and discarded afterwards, because doing so can
#' evaluate \code{log(0)}, negative powers, or repeated-power log terms at
#' structural-zero rows.
#'
#' @param catzero Named logical vector, one entry per variable.
#' @param spike Named logical vector, one entry per variable.
#' @param spike_decision Named numeric vector with values in
#'   \code{unname(saz_decision_codes)}, or \code{NULL} (in which case this
#'   function is a no-op).
#'
#' @return A list with the (possibly updated) \code{catzero} and
#'   \code{binary_only_spike_vars} (character vector of variable names whose
#'   positive-value component must not be constructed).
#' @keywords internal
#' @noRd
identify_binary_only_spike_vars <- function(catzero, spike, spike_decision) {
  binary_only_spike_vars <- character(0L)

  if (!is.null(spike_decision)) {
    for (v in names(spike_decision)) {
      if (isTRUE(spike[[v]])) {
        dec <- as.integer(spike_decision[[v]])

        if (dec == saz_decision_codes[["continuous_only"]]) {
          catzero[[v]] <- FALSE
        } else if (dec == saz_decision_codes[["binary_only"]]) {
          catzero[[v]] <- TRUE
          binary_only_spike_vars <- c(binary_only_spike_vars, v)
        }
      }
    }
  }

  list(catzero = catzero, binary_only_spike_vars = binary_only_spike_vars)
}

#' Build the Per-Variable FP or ACD Transformed Columns
#'
#' Internal helper used by \code{transform_matrix()}. For each variable in
#' \code{names_vars}, applies either the ACD transformation
#' (\code{transform_vector_acd()}) or the plain FP transformation
#' (\code{transform_vector_fp()}), depending on \code{acdx}, and tracks which
#' resulting columns are ACD-component columns (as opposed to the FP part),
#' since those must not inherit zero-aware centering later. Variables listed
#' in \code{binary_only_spike_vars} are skipped entirely (their positive-value
#' component must never be constructed; see
#' \code{identify_binary_only_spike_vars()}).
#'
#' @param x Numeric matrix, already reordered to \code{names_vars} order.
#' @param names_vars Character vector of variable names, in processing order.
#' @param power_list Named list of FP/ACD powers, one element per variable.
#' @param acdx Named logical vector: \code{TRUE} to use the ACD
#'   transformation for that variable.
#' @param zero Named logical vector: \code{TRUE} to treat exact-zero values of
#'   a nonnegative variable as structural zero before transformation.
#' @param acd_parameter_list Optional named list of previously-fitted ACD
#'   parameters (supplied during prediction; \code{NULL} during fitting, in
#'   which case parameters are estimated).
#' @param check_binary Passed through to \code{transform_vector_fp()}.
#' @param binary_only_spike_vars Character vector of variable names to skip,
#'   as returned by \code{identify_binary_only_spike_vars()}.
#'
#' @return A list with:
#'   * \code{x_trafo}: named list of transformed column matrices/vectors, one
#'     per variable (\code{NULL} entries for all-NA powers and for
#'     \code{binary_only_spike_vars}).
#'   * \code{acd_parameter}: named list of ACD parameters, one per ACD
#'     variable.
#'   * \code{acd_component_cols}: character vector of column names that are
#'     ACD-component columns (not the FP part of an ACD transformation).
#' @keywords internal
#' @noRd
build_fp_acd_columns <- function(x, names_vars, power_list, acdx, zero,
                                 acd_parameter_list, check_binary,
                                 binary_only_spike_vars) {
  x_trafo <- list()
  acd_parameter <- list()
  acd_component_cols <- character(0L)

  for (name in names_vars) {
    if (name %in% binary_only_spike_vars) {
      x_trafo[[name]] <- NULL
      next
    }

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

  list(
    x_trafo            = x_trafo,
    acd_parameter       = acd_parameter,
    acd_component_cols = acd_component_cols
  )
}

#' Apply the Spike Decision to the Transformed FP/ACD Columns (Defensive Pass)
#'
#' Internal helper used by \code{transform_matrix()}. The SAZ decisions that
#' affect component construction were already applied before the
#' transformation loop (see \code{identify_binary_only_spike_vars()}). This
#' defensive pass keeps the final component list in sync with
#' \code{spike_decision} in case future callers modify the earlier handling:
#' \code{spike_decision == "continuous_only"} suppresses the binary indicator
#' (continuous only); \code{spike_decision == "binary_only"} removes the
#' continuous component and forces the binary indicator (binary only);
#' \code{spike_decision == "continuous_and_binary"} leaves both as-is.
#'
#' @param x_trafo Named list of transformed columns, as returned by
#'   \code{build_fp_acd_columns()}.
#' @param catzero Named logical vector, one entry per variable.
#' @param spike Named logical vector, one entry per variable.
#' @param spike_decision Named numeric vector with values in
#'   \code{unname(saz_decision_codes)}, or \code{NULL} (in which case this
#'   function is a no-op).
#'
#' @return A list with the (possibly updated) \code{x_trafo} and
#'   \code{catzero}.
#' @keywords internal
#' @noRd
apply_spike_decision_to_columns <- function(x_trafo, catzero, spike,
                                            spike_decision) {
  if (!is.null(spike_decision)) {
    for (v in names(spike_decision)) {
      dec <- as.integer(spike_decision[[v]])

      if (isTRUE(spike[[v]])) {
        if (dec == saz_decision_codes[["continuous_only"]]) {
          catzero[[v]] <- FALSE
        } else if (dec == saz_decision_codes[["binary_only"]]) {
          x_trafo[[v]] <- NULL
          catzero[[v]] <- TRUE
        }
        # dec == 1: FP/ACD + *_bin; no change.
      }
    }
  }

  list(x_trafo = x_trafo, catzero = catzero)
}

#' Assemble the Transformed Matrix and Add Binary Indicator Columns
#'
#' Internal helper used by \code{transform_matrix()}. First combines
#' \code{x_trafo} into a single matrix (or \code{NULL} if empty, which can
#' happen after a binary-only spike decision). Then, for every variable
#' flagged in \code{catzero}, builds the structural-zero binary indicator
#' column \code{I(x == 0)} and appends it, consistent with exact-zero handling
#' in \code{transform_vector_fp(..., zero = TRUE)}.
#'
#' @param x Numeric matrix, already reordered to \code{names(catzero)} order.
#' @param x_trafo Named list of transformed columns, as returned by
#'   \code{apply_spike_decision_to_columns()}.
#' @param catzero Named logical vector, one entry per variable.
#' @param spike,spike_decision Used only to detect binary-only spike
#'   variables, whose all-NA power entry in \code{power_list} should not
#'   cause the binary indicator to be skipped.
#' @param power_list Named list of FP/ACD powers, one element per variable.
#'
#' @return A list with:
#'   * \code{x_transformed}: the combined matrix (FP/ACD columns plus any
#'     binary indicator columns), or \code{NULL} if there are none.
#'   * \code{bin_col_to_var}: named character vector mapping each \code{*_bin}
#'     column name to its source variable.
#' @keywords internal
#' @noRd
build_binary_indicator_columns <- function(x, x_trafo, catzero, spike,
                                           spike_decision, power_list) {
  # x_trafo can be empty after binary-only spike decisions. In that case, the
  # matrix may still be created below from catzero/spike binary indicators.
  if (length(x_trafo) > 0L) {
    x_transformed <- do.call(cbind, x_trafo)
  } else {
    x_transformed <- NULL
  }

  cat_vars <- names(catzero)[catzero]
  bin_col_to_var <- setNames(character(0L), character(0L))

  if (length(cat_vars) > 0L) {
    catzero_list <- lapply(cat_vars, function(v) {
      binary_only_spike <- isTRUE(spike[[v]]) &&
        !is.null(spike_decision) &&
        v %in% names(spike_decision) &&
        isTRUE(spike_decision[[v]] == saz_decision_codes[["binary_only"]])

      # all powers NA usually means the variable was eliminated. The exception
      # is spike_decision == 3, where the binary indicator is the selected term.
      if (all(is.na(power_list[[v]])) && !binary_only_spike) {
        return(NULL)
      }

      # Structural-zero indicator is I(x == 0), consistent with
      # transform_vector_fp(..., zero = TRUE).
      as.integer(x[, v] == 0)
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

  list(x_transformed = x_transformed, bin_col_to_var = bin_col_to_var)
}

#' Map Every Transformed Column Back to Its Original Variable
#'
#' Internal helper used by \code{transform_matrix()}. After transformation,
#' an original variable may expand into several final columns (e.g.
#' \code{x -> x.1, x.2, ...}; \code{acd(x) -> x.1} and/or \code{A_x.1};
#' \code{catzero}/\code{spike -> x_bin}). This builds the map from every
#' final column name back to its source variable, which
#' \code{expand_zero_flags_to_columns()} and \code{apply_mixed_centering()}
#' both need for per-variable-derived centering and zero handling.
#'
#' @param x_trafo Named list of transformed columns, as returned by
#'   \code{apply_spike_decision_to_columns()}.
#' @param x_transformed The combined matrix, as returned by
#'   \code{build_binary_indicator_columns()}.
#' @param bin_col_to_var Named character vector mapping \code{*_bin} column
#'   names to their source variable, as returned by
#'   \code{build_binary_indicator_columns()}.
#'
#' @return A list with:
#'   * \code{col_to_var}: named character vector, one entry per column of
#'     \code{x_transformed}, giving its source variable.
#'   * \code{bin_cols}: character vector of the \code{*_bin} column names
#'     actually present in \code{x_transformed}.
#' @keywords internal
#' @noRd
build_col_to_var_map <- function(x_trafo, x_transformed, bin_col_to_var) {
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

  list(col_to_var = col_to_var, bin_cols = bin_cols)
}

#' Expand the Per-Variable Zero Flag to Transformed Columns
#'
#' Internal helper used by \code{transform_matrix()}. \code{zero} is defined
#' per original variable, but \code{center_matrix()} needs one zero flag per
#' transformed column. For zero-handled FP columns (including the FP
#' component of an ACD transformation), \code{zero_expanded = TRUE} makes
#' \code{center_matrix()} compute the center using only nonzero transformed
#' values and reset transformed zero rows back to 0 after centering.
#' ACD-component columns and \code{*_bin} columns are excluded, since ACD
#' components are cumulative-probability quantities (not structural zeros)
#' and binary indicators are handled by \code{center_matrix()}'s
#' binary-minimum rule instead.
#'
#' @param x_transformed The combined transformed matrix.
#' @param col_to_var Named character vector mapping each transformed column
#'   to its source variable, as returned by \code{build_col_to_var_map()}.
#' @param zero Named logical vector, one entry per original variable.
#' @param acd_component_cols Character vector of ACD-component column names.
#' @param bin_cols Character vector of \code{*_bin} column names.
#'
#' @return Named logical vector, one entry per column of \code{x_transformed}.
#' @keywords internal
#' @noRd
expand_zero_flags_to_columns <- function(x_transformed, col_to_var, zero,
                                         acd_component_cols, bin_cols) {
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

  zero_expanded
}

#' Apply Per-Variable Centering to the Transformed Matrix
#'
#' Internal helper used by \code{transform_matrix()}. \code{center} is
#' defined per original variable: \code{center[v] = TRUE} centers all final
#' columns derived from \code{v}; \code{center[v] = FALSE} leaves them
#' unchanged. Binary \code{*_bin} columns follow their original variable; if
#' centered, \code{center_matrix()} centers them by their minimum, so usual
#' 0/1 indicators remain unchanged. Uses a fast path (center the whole matrix
#' in one call) when every variable requests centering, and a mixed path
#' (center only the relevant columns) otherwise.
#'
#' @param x_transformed The combined transformed matrix.
#' @param center Named logical vector, one entry per original variable.
#' @param col_to_var Named character vector mapping each transformed column
#'   to its source variable, as returned by \code{build_col_to_var_map()}.
#' @param zero_expanded Named logical vector, one entry per column of
#'   \code{x_transformed}, as returned by
#'   \code{expand_zero_flags_to_columns()}.
#'
#' @return A list with the (possibly centered) \code{x_transformed} and
#'   \code{centers} (named numeric vector of centering constants actually
#'   applied, or \code{NULL} if no variable requested centering).
#' @keywords internal
#' @noRd
apply_mixed_centering <- function(x_transformed, center, col_to_var,
                                  zero_expanded) {
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

  list(x_transformed = x_transformed, centers = centers)
}



#' Transform a Vector by One Fractional-Polynomial Power
#'
#' Applies a single fractional-polynomial power transformation to a numeric
#' vector. This helper is a scalar-power wrapper around the shared C++ FP
#' transformation core used by \code{transform_vector_fp()}.
#'
#' @details
#' The transformation rule is:
#' \itemize{
#'   \item \code{power = 0}: return \eqn{\log(x)}.
#'   \item \code{power != 0}: return \eqn{x^\code{power}}.
#' }
#'
#' This helper is intentionally limited to one power. Full FP bases with
#' multiple powers, such as \code{c(p1, p2)}, should be constructed with
#' \code{transform_vector_fp()} instead. The scalar-power check prevents a
#' silent error where the C++ core would return multiple columns and this helper
#' would otherwise keep only the first one.
#'
#' When \code{zero = TRUE}, exact-zero values are treated as structural zeros:
#' rows with \code{x == 0} are returned as zero and are not evaluated by
#' \code{log()} or power operations. This prevents invalid evaluations such as
#' \code{log(0)} or \code{0^(-1)}. Callers must supply previously validated
#' nonnegative values when zero handling is active. Missing and
#' non-finite values are propagated
#' by the shared C++ transformation core.
#'
#' The function passes \code{shift = 0} and \code{scale = 1} to the C++ core
#' because its contract is to transform the vector exactly as supplied, not to
#' apply additional preprocessing.
#'
#' @param x Numeric vector to transform and must have positive values.
#' @param power Numeric scalar. Fractional-polynomial power to apply. A value
#'   of \code{0} represents the logarithmic transformation.
#' @param zero Logical scalar. If \code{TRUE}, transform only positive values
#'   of a nonnegative vector and return zero for exact-zero values. The caller
#'   is responsible for enforcing nonnegative input. If \code{FALSE}, transform all
#'   values as supplied.
#'
#' @return
#' A numeric vector of transformed values with the same length as \code{x}.
#'
#' @seealso
#' \code{\link{transform_vector_fp}}
#'
#' @keywords internal
#' @noRd

transform_vector_single_power <- function(x, power = 1, zero = FALSE) {
  # This helper is intentionally a single-power wrapper. It is used where one
  # scalar FP power is expected, not where a full FP basis such as c(p1, p2) is
  # requested. Guarding here prevents a silent bug if a vector of powers is
  # accidentally supplied: transform_fp_core() would return one column per
  # power, and taking [, 1L] would silently discard the remaining columns.
  if (!is.numeric(power) || length(power) != 1L || is.na(power) || !is.finite(power)) {
    stop(
      "'power' must be a single finite, non-missing numeric value.",
      call. = FALSE
    )
  }

  if (!is.logical(zero) || length(zero) != 1L || is.na(zero)) {
    stop("'zero' must be a single non-missing logical value.", call. = FALSE)
  }

  # Delegate to the same C++ kernel used by transform_vector_fp(). This keeps
  # the single-power helper numerically consistent with the package-wide FP
  # implementation:
  #   power = 0   -> log(x)
  #   power != 0  -> x^power
  #   zero = TRUE -> rows with x == 0 are structural zeros and are not evaluated
  #                  by log() or pow(); they remain zero in the returned vector.
  #
  # The helper passes shift = 0 and scale = 1 because its original contract was
  # to transform the vector exactly as supplied, not to apply extra preprocessing.
  transform_fp_core(
    x_raw     = as.numeric(x),
    power     = as.numeric(power),
    shift_val = 0,
    scale_val = 1,
    zero      = zero
  )[, 1L]
}

#' Does an FP power specification require strictly positive input?
#'
#' Internal helper used by transformation and prediction code.
#'
#' Fractional-polynomial terms involving logarithms require strictly positive
#' input. This includes power 0, which represents log(x), and repeated-power
#' terms such as c(2, 2), where the second term is x^2 * log(x).
#'
#' Positive integer powers such as 1, 2, and 3 do not require strictly positive
#' input, because they are well-defined for zero and negative values.
#'
#' Negative powers and non-integer powers are treated as requiring strictly
#' positive input in this package. This is conservative, but consistent with the
#' usual MFP convention that shifted/scaled covariates should be positive for
#' non-linear FP transformations.
#'
#' @param power Numeric vector of selected FP powers.
#'
#' @return A single logical value.
#'
#' @keywords internal
#' @noRd
fp_power_requires_positive_input <- function(power) {
  if (is.null(power) || length(power) == 0L || all(is.na(power))) {
    return(FALSE)
  }

  p <- as.numeric(power[!is.na(power)])

  if (length(p) == 0L) {
    return(FALSE)
  }

  # Repeated powers require log(x), e.g. c(2, 2) gives
  # x^2 and x^2 * log(x).
  if (length(p) > length(unique(p))) {
    return(TRUE)
  }

  # Power 0 represents log(x).
  if (any(p == 0)) {
    return(TRUE)
  }

  # Negative powers are undefined at zero and unstable near zero. For package
  # consistency, require strictly positive input rather than allowing negative
  # bases for negative integer powers.
  if (any(p < 0)) {
    return(TRUE)
  }

  # Non-integer powers are not real-valued for negative inputs.
  if (any(p != floor(p))) {
    return(TRUE)
  }

  # Remaining case: distinct positive integer powers, e.g. 1, 2, 3.
  FALSE
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
#' @keywords internal
#' @noRd
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
    paste0(name, ".", seq_len(n_powers))
  } else {
    c(paste0(name, ".1"), paste0("A_", name, ".1"))
  }
}

#' Create Indicator Variables for Categorical Predictors
#'
#' Creates numeric indicator columns for ordinal or nominal variables.
#'
#' Ordinal variables are encoded using cumulative threshold indicators.
#' Nominal variables are encoded using treatment contrasts, with the first
#' factor level used as the reference category.
#'
#' @details
#' This function is primarily useful when preparing predictors for the matrix
#' interfaces of [mfp2()] or [mfpi()].
#'
#' When using the formula interface, categorical variables can usually be
#' supplied directly as factors. The formula interface constructs the required
#' contrast columns and treats all columns belonging to the factor as one
#' conceptual model term.
#'
#' For an ordinal variable with ordered levels `A < B < C < D`, the function
#' creates cumulative indicators corresponding to:
#'
#' \itemize{
#'   \item `B`, `C`, or `D` versus `A`;
#'   \item `C` or `D` versus `A` or `B`;
#'   \item `D` versus `A`, `B`, or `C`.
#' }
#'
#' The level order determines the comparisons. Users should therefore convert
#' categorical variables to factors with explicitly defined levels before
#' calling this function.
#'
#' A variable cannot be specified in both `var_ordinal` and `var_nominal`.
#'
#' When the resulting columns are passed to the matrix interface, use
#' `term_groups` to identify all indicator columns derived from the same
#' categorical variable. This ensures that they are selected or retained
#' jointly.
#'
#' @param data A data frame containing the variables to encode.
#'
#' @param var_ordinal Optional character vector naming ordinal variables.
#'   Variables are encoded using cumulative threshold indicators according to
#'   their factor-level order. If a variable is not a factor, its sorted unique
#'   values determine the order.
#'
#' @param var_nominal Optional character vector naming nominal variables.
#'   Variables are encoded using treatment contrasts, with the first factor
#'   level as the reference category.
#'
#' @param drop_variables Logical scalar. If `TRUE`, the original categorical
#'   variables are removed after the indicator columns are created. The default
#'   is `FALSE`.
#'
#' @return
#' A data frame containing the original data and the newly created indicator
#' columns. When `drop_variables = TRUE`, the encoded source variables are
#' omitted.
#'
#' @examples
#' data("gbsg")
#'
#' # Define the intended ordering explicitly.
#' gbsg$grade <- ordered(
#'   gbsg$grade,
#'   levels = sort(unique(gbsg$grade))
#' )
#'
#' # Create cumulative indicators for tumour grade.
#' gbsg_encoded <- create_dummy_variables(
#'   data = gbsg,
#'   var_ordinal = "grade",
#'   drop_variables = TRUE
#' )
#'
#' head(gbsg_encoded)
#'
#' \dontrun{
#' # Formula-interface users can ordinarily supply the factor directly.
#' fit_formula <- mfp2(
#'   survival::Surv(rectime, censrec) ~ fp(age) + grade,
#'   data = gbsg,
#'   family = "cox",
#'   keep = "grade",
#'   verbose = FALSE
#' )
#'
#' # For a matrix interface, group the generated grade indicators as one term.
#' grade_columns <- grep("^grade_", names(gbsg_encoded), value = TRUE)
#'
#' x <- as.matrix(
#'   gbsg_encoded[, c("age", "nodes", grade_columns), drop = FALSE]
#' )
#' y <- survival::Surv(gbsg_encoded$rectime, gbsg_encoded$censrec)
#'
#' fit_matrix <- mfp2(
#'   x = x,
#'   y = y,
#'   family = "cox",
#'   term_groups = list(grade = grade_columns),
#'   keep = "grade",
#'   verbose = FALSE
#' )
#' }
#'
#' @keywords internal
#' @noRd
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
#' @return A numeric matrix with \code{nlev} rows and \code{nlev - 1}
#'   columns, where \code{nlev} is either \code{n} when \code{n} is numeric
#'   or \code{length(n)} when level labels are supplied.
#' @keywords internal
#' @noRd
contr.cumulative <- function(n) {
  # Numeric input represents the number of ordered levels. Keep this contract
  # explicit so invalid scalar values do not reach matrix dimensions or create
  # unsafe colon sequences such as 1:0. Non-numeric input supplies the labels.
  if (is.numeric(n)) {
    if (length(n) != 1L || is.na(n) || !is.finite(n) ||
        n < 2 || n != floor(n)) {
      stop(
        "`n` must specify at least two ordered levels.",
        call. = FALSE
      )
    }
    nlev <- as.integer(n)
  } else {
    nlev <- length(n)
    if (nlev < 2L) {
      stop(
        "`n` must specify at least two ordered levels.",
        call. = FALSE
      )
    }
  }

  mat <- matrix(0, nrow = nlev, ncol = nlev - 1L)
  for (j in seq_len(nlev - 1L)) {
    mat[seq.int(j + 1L, nlev), j] <- 1
  }
  rownames(mat) <- if (is.numeric(n)) as.character(seq_len(nlev)) else n
  # Safe column names without special characters
  colnames(mat) <- paste0("_", seq_len(nlev - 1L))
  mat
}
