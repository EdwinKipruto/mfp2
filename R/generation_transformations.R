#' Function to generate all requested FP transformations for a single variable
#'
#' @param x A numeric vector of length \code{nobs}, assumed to have been shifted
#' (except for variables with \code{zero} or \code{catzero} transformations) and
#' scaled.
#' @param degree A numeric value indicating the degree of fractional polynomials (FPs).
#' For ACD transformation, this is assumed to be 2.
#' @param powers A numeric vector specifying the set of allowed fractional
#' polynomial (FP) powers to be used in the transformation.
#' @param zero Logical indicating whether only positive values of the variable
#' should be transformed, with nonpositive values (zero or negative) set to zero.
#' If \code{TRUE}, transformation is applied only to positive values; nonpositive
#' values are replaced with zero before transformation. If \code{FALSE} (default),
#' all values are shifted (if needed) to ensure positivity before transformation.
#' @param catzero An optional n x 1 numeric/integer matrix containing the
#' structural-zero binary indicator for the current variable.
#' @details
#' Any fractional polynomial (FP) transformation is defined by a vector of powers,
#' such as \code{(p1, p2)} for degree 2. These correspond to the terms \code{x^p1}
#' and \code{x^p2}. Therefore, all combinations of the values in \code{powers}
#' are considered (see \code{generate_powers_fp}).
#'
#' A special case arises when powers are repeated, i.e., \code{p1 = p2}. In such cases,
#' the second term is multiplied by \code{log(x)}, following the standard FP convention
#' (see \code{transform_vector_fp}).
#'
#' When the ACD transformation is requested, all pairs of powers of length 2 are evaluated,
#' resulting in 64 unique combinations (see \code{generate_powers_acd}).
#'
#' If `degree = 0` then these functions return the data unchanged for fp,
#' or simply the acd transformation of the input variable, i.e. in both cases
#' the power is set to 1 (linear).
#'
#' If \code{degree = 0}, the function returns the data unchanged for an FP transformation,
#' or applies only the ACD transformation to the input variable. In both cases,
#' the power is set to 1 (linear).
#'
#' When \code{catzero} is used, the transformed (or untransformed) continuous variable is
#' combined with its corresponding binary indicator, representing whether the original
#' value was positive or nonpositive.
#'
#' @return
#' A list with two components:
#'
#' * \code{data}: A list with length equal to the number of possible fractional
#'  polynomial (FP) transformations for the variable of interest. Each entry is
#'  a matrix with \code{nobs} rows. The number of columns equals the
#'  FP \code{degree}, unless \code{catzero = TRUE}, in which case an additional
#'   column is included for the binary indicator variable. For example, with
#'   \code{degree = 2}, \code{catzero = TRUE}, and \code{nobs = 10}, each entry
#'   is a \eqn{10 \times 3} matrix.
#'   The FP-transformed values are not centered. If \code{degree = 0}, the list contains a single
#'   entry with one column (or two columns if \code{catzero = TRUE}), representing the linear
#'   transformation (and binary indicator, if applicable).
#' *  \code{data}: the associated FP powers for each entry in data.
#' * \code{powers}: A matrix of FP powers corresponding to each entry in \code{data}.
#'   Each row contains the powers used for the associated transformation (e.g.,
#'  two columns for \code{degree = 2}, one for \code{degree = 1}, and one for
#'   \code{degree = 0}).
#' @keywords internal
#' @noRd
generate_transformations_fp <- function(x,
                                        degree,
                                        powers,
                                        zero,
                                        catzero = NULL) {

  # Validate catzero, if provided
  if (!is.null(catzero)) {
    if (!is.matrix(catzero) || ncol(catzero) != 1L) {
      stop("`catzero` must be an n x 1 matrix.")
    }

    if (nrow(catzero) != length(x)) {
      stop("`catzero` must have one row per observation in `x`.")
    }

    if (!is.numeric(catzero)) {
      stop("`catzero` must be a numeric/integer matrix.")
    }
  }

  # All possible combinations of powers given degree.
  # Keep this in R because it defines the package's FP candidate semantics.
  combs <- generate_powers_fp(
    degree = degree,
    powers = powers
  )

  # Generate all candidate FP transformations in one C++ call.
  fpdt <- generate_transformations_fp_cpp(
    x = as.numeric(x),
    powers = combs,
    zero = isTRUE(zero),
    catzero = catzero
  )

  list(
    data = fpdt,
    powers = combs
  )
}


#' Generate a compact reusable basis for ordinary FP candidate fitting
#'
#' This is an mfp2 fitting-path optimization. Unlike
#' \code{generate_transformations_fp()}, it does not materialize a complete
#' transformed matrix for every candidate power combination. Instead it
#' computes each distinct FP basis term once and returns a small integer map
#' that identifies which basis columns form each candidate.
#'
#' The structural-zero indicator supplied through \code{catzero} is retained
#' separately because it is identical for all FP candidates and therefore
#' should not be copied into every candidate matrix.
#'
#' @inheritParams generate_transformations_fp
#' @return A list with components \code{basis}, \code{candidate_map},
#'   \code{powers}, and \code{catzero}. \code{basis} is an n x b matrix of
#'   unique FP basis columns. \code{candidate_map} is an integer matrix whose
#'   rows correspond to \code{powers} and whose entries index columns of
#'   \code{basis}.
#' @keywords internal
#' @noRd
generate_transformations_fp_basis <- function(x,
                                              degree,
                                              powers,
                                              zero,
                                              catzero = NULL) {

  if (!is.null(catzero)) {
    if (!is.matrix(catzero) || ncol(catzero) != 1L) {
      stop("`catzero` must be an n x 1 matrix.")
    }

    if (nrow(catzero) != length(x)) {
      stop("`catzero` must have one row per observation in `x`.")
    }

    if (!is.numeric(catzero)) {
      stop("`catzero` must be a numeric/integer matrix.")
    }
  }

  combs <- generate_powers_fp(
    degree = degree,
    powers = powers
  )

  compact <- generate_transformations_fp_basis_cpp(
    x = as.numeric(x),
    powers = combs,
    zero = isTRUE(zero)
  )

  list(
    basis = compact$basis,
    candidate_map = compact$candidate_map,
    powers = combs,
    catzero = catzero
  )
}


#' Materialize one ordinary FP candidate from a compact basis
#'
#' Used only when a complete candidate matrix is genuinely required, notably
#' for the winning stage-1 SAZ transformation that is reused in stage 2.
#' Candidate fitting itself writes basis columns directly into a reusable
#' design matrix and does not call this helper for every candidate.
#'
#' @param fp_basis A compact basis object returned by
#'   \code{generate_transformations_fp_basis()}.
#' @param candidate Integer row index of the requested candidate.
#' @return An n x degree matrix, or n x (degree + 1) when \code{catzero} is
#'   present.
#' @keywords internal
#' @noRd
materialize_fp_basis_candidate <- function(fp_basis, candidate) {
  candidate <- as.integer(candidate)

  if (length(candidate) != 1L || is.na(candidate) ||
      candidate < 1L || candidate > nrow(fp_basis$candidate_map)) {
    stop("Internal error: FP candidate index is out of range.", call. = FALSE)
  }

  source_cols <- fp_basis$candidate_map[candidate, , drop = TRUE]
  fp_data <- fp_basis$basis[, source_cols, drop = FALSE]

  if (is.null(fp_basis$catzero)) {
    return(fp_data)
  }

  out <- matrix(
    0,
    nrow = nrow(fp_data),
    ncol = ncol(fp_data) + 1L
  )
  out[, 1L] <- fp_basis$catzero[, 1L]
  out[, -1L] <- fp_data

  colnames(out) <- c(
    "catzero",
    paste0("V", seq_len(ncol(fp_data)))
  )

  out
}


#' Resolve the base ACD values used during candidate generation
#'
#' During fit_mfp(), fit_acd() is run once per ACD variable and its fitted
#' training A(x) vector is retained in acd_parameter$acd. The model-search hot
#' path can pass that vector explicitly via acd_training_values and avoid
#' applying the same fitted ACD transformation to the same training x again.
#'
#' The cache is deliberately explicit rather than inferred from
#' acd_parameter$acd. This preserves the existing apply-to-new-x semantics for
#' direct/internal callers that supply stored ACD parameters with a different x.
#' If no training cache is supplied, the historical fit/apply behaviour is used.
#'
#' @param x Numeric covariate vector used for candidate generation.
#' @param powers Candidate powers used if ACD parameters must be fitted here.
#' @param zero Logical indicating zero-component handling.
#' @param acd_parameter Optional stored ACD parameter list.
#' @param acd_training_values Optional cached A(x) values for exactly the same
#'   training observations, in the same order as x.
#' @return Numeric vector containing the base A(x) values.
#' @keywords internal
#' @noRd
resolve_acd_base_values <- function(x,
                                    powers,
                                    zero,
                                    acd_parameter = NULL,
                                    acd_training_values = NULL) {
  # Fast training-data path: reuse fit_acd()$acd instead of recalculating the
  # FP term, linear predictor, and pnorm() for every focal ACD search.
  if (!is.null(acd_training_values)) {
    # Keep validation O(1): scanning the complete cached vector here would add
    # an n-length pass to the hot path and partly defeat the optimization.
    if (!is.numeric(acd_training_values) ||
        !is.null(dim(acd_training_values)) ||
        length(acd_training_values) != length(x)) {
      stop(
        paste0(
          "Internal error: `acd_training_values` must be a numeric vector ",
          "with one value per observation in `x`."
        ),
        call. = FALSE
      )
    }

    return(acd_training_values)
  }

  # Historical fallback for standalone/internal calls that have no fit-time
  # training cache. If parameters are absent, fit ACD here as before.
  if (is.null(acd_parameter)) {
    acd_parameter_work <- fit_acd(
      x = x,
      powers = powers,
      shift = 0,
      scale = 1,
      zero = zero
    )

    return(acd_parameter_work$acd)
  }

  # Apply stored parameters to the supplied x. Remove a possibly stored $acd
  # component so direct callers do not accidentally reuse values belonging to
  # another data vector with the same length.
  acd_parameter_apply <- acd_parameter
  acd_parameter_apply$acd <- NULL

  do.call(
    apply_acd,
    utils::modifyList(
      acd_parameter_apply,
      list(
        x = x,
        zero = zero
      )
    )
  )
}


#' Generate a compact shared basis for ACD candidates
#'
#' This is the ACD analogue of \code{generate_transformations_fp_basis()}.
#' Instead of materializing one complete n x degree matrix for every ACD power
#' pair, it stores each unique transformed column once in \code{basis} and
#' represents each candidate by one row of the small integer
#' \code{candidate_map}. The structural-zero indicator, when present, is kept
#' separately because it is invariant across all candidates.
#'
#' For the default eight powers, ACD degree 2 therefore stores 16 transformed
#' columns (8 for x and 8 for A(x)) plus a 64 x 2 integer map, rather than 64
#' separate n x 2 candidate matrices.
#'
#' @inheritParams generate_transformations_acd
#' @return A list with \code{basis}, \code{candidate_map}, \code{powers}, and
#'   \code{catzero}.
#' @keywords internal
#' @noRd
generate_transformations_acd_basis <- function(x,
                                               degree,
                                               powers,
                                               zero,
                                               catzero = NULL,
                                               acd_parameter = NULL,
                                               acd_training_values = NULL) {
  # Match the validation contract of generate_transformations_acd(). Keeping
  # both interfaces aligned allows the compact path to be substituted only in
  # the model-search hot loop without changing public/internal callers that
  # still expect a materialized list of candidates.
  if (!is.null(catzero)) {
    if (!is.matrix(catzero) || ncol(catzero) != 1L) {
      stop("`catzero` must be an n x 1 matrix.")
    }

    if (nrow(catzero) != length(x)) {
      stop("`catzero` must have one row per observation in `x`.")
    }

    if (!is.numeric(catzero)) {
      stop("`catzero` must be a numeric/integer matrix.")
    }
  }

  combs <- generate_powers_acd(
    degree = degree,
    powers = powers
  )

  # ACD candidates have one active transformed term for degrees 0/1 and two
  # active terms for degree 2. Validate that this width is constant so the
  # reusable design matrix can reserve one fixed focal block for the search.
  n_terms_per_candidate <- rowSums(!is.na(combs))
  if (length(unique(n_terms_per_candidate)) != 1L) {
    stop(
      "Internal error: ACD candidates have inconsistent numbers of terms.",
      call. = FALSE
    )
  }
  n_terms <- n_terms_per_candidate[1L]

  # Resolve A(x) once for this candidate basis. During normal mfp2 training,
  # transform_data_step() supplies the A(x) vector already computed by fit_acd(),
  # so no repeated apply_acd() work is needed for the same training observations.
  # Direct/internal calls without that explicit cache retain the historical
  # fit/apply behaviour.
  x_acd_base <- resolve_acd_base_values(
    x = x,
    powers = powers,
    zero = zero,
    acd_parameter = acd_parameter,
    acd_training_values = acd_training_values
  )

  # Store each unique power transform once for each ACD component. The first
  # component is FP(x, p1); the second is FP(A(x), p2). Degree 1 has no active
  # first component, so x_powers is empty and the basis contains A(x) terms only.
  x_powers <- unique(combs[!is.na(combs[, 1L]), 1L])
  acd_powers <- unique(combs[!is.na(combs[, 2L]), 2L])

  n_x_basis <- length(x_powers)
  n_acd_basis <- length(acd_powers)
  basis <- matrix(
    NA_real_,
    nrow = length(x),
    ncol = n_x_basis + n_acd_basis
  )

  if (n_x_basis > 0L) {
    for (j in seq_along(x_powers)) {
      basis[, j] <- transform_vector_fp(
        x = x,
        power = x_powers[j],
        scale = 1,
        shift = 0,
        zero = zero
      )
    }
  }

  if (n_acd_basis > 0L) {
    for (j in seq_along(acd_powers)) {
      basis[, n_x_basis + j] <- transform_vector_fp(
        x = x_acd_base,
        power = acd_powers[j],
        scale = 1,
        shift = 0,
        zero = FALSE
      )
    }
  }

  # candidate_map uses 1-based basis indices because it is consumed directly
  # by copy_fp_basis_candidate_cpp(). No n-length candidate object is created.
  candidate_map <- matrix(
    NA_integer_,
    nrow = nrow(combs),
    ncol = n_terms
  )

  for (i in seq_len(nrow(combs))) {
    k <- 1L

    if (!is.na(combs[i, 1L])) {
      candidate_map[i, k] <- match(combs[i, 1L], x_powers)
      k <- k + 1L
    }

    if (!is.na(combs[i, 2L])) {
      candidate_map[i, k] <- n_x_basis + match(combs[i, 2L], acd_powers)
    }
  }

  list(
    basis = basis,
    candidate_map = candidate_map,
    powers = combs,
    catzero = catzero
  )
}


#' Materialize one ACD candidate from a compact basis
#'
#' Candidate fitting itself uses \code{candidate_map} to copy shared basis
#' columns directly into a reusable design matrix. This helper is called only
#' after model selection when downstream SAZ/current-parameter logic genuinely
#' needs the complete winning focal matrix.
#'
#' @param acd_basis A compact basis object returned by
#'   \code{generate_transformations_acd_basis()}.
#' @param candidate Integer row index of the requested candidate.
#' @return An n x degree matrix, or n x (degree + 1) matrix when
#'   \code{catzero} is present.
#' @keywords internal
#' @noRd
materialize_acd_basis_candidate <- function(acd_basis, candidate) {
  candidate <- as.integer(candidate)

  if (length(candidate) != 1L || is.na(candidate) ||
      candidate < 1L || candidate > nrow(acd_basis$candidate_map)) {
    stop("Internal error: ACD candidate index is out of range.", call. = FALSE)
  }

  source_cols <- acd_basis$candidate_map[candidate, , drop = TRUE]
  acd_data <- acd_basis$basis[, source_cols, drop = FALSE]

  if (is.null(acd_basis$catzero)) {
    return(acd_data)
  }

  out <- matrix(
    0,
    nrow = nrow(acd_data),
    ncol = ncol(acd_data) + 1L
  )
  out[, 1L] <- acd_basis$catzero[, 1L]
  out[, -1L] <- acd_data

  colnames(out) <- c(
    "catzero",
    paste0("V", seq_len(ncol(acd_data)))
  )

  out
}


#' @describeIn generate_transformations_fp Function to generate acd transformations.
#' @param acd_parameter Optional named list of ACD parameters, generated by
#'   \code{fit_acd()}. If `NULL`, the function generates it.
#' @param acd_training_values Optional cached A(x) vector for the same training
#'   observations in `x`. This is a fit-time optimization only; leave `NULL`
#'   when applying stored ACD parameters to different/new data.
#' @importFrom utils modifyList
#' @keywords internal
#' @noRd
generate_transformations_acd <- function(x,
                                         degree,
                                         powers,
                                         zero,
                                         catzero = NULL,
                                         acd_parameter = NULL,
                                         acd_training_values = NULL) {

  # Validate catzero, if provided
  if (!is.null(catzero)) {
    if (!is.matrix(catzero) || ncol(catzero) != 1L) {
      stop("`catzero` must be an n x 1 matrix.")
    }

    if (nrow(catzero) != length(x)) {
      stop("`catzero` must have one row per observation in `x`.")
    }

    if (!is.numeric(catzero)) {
      stop("`catzero` must be a numeric/integer matrix.")
    }
  }

  # All possible ACD power combinations. For ACD, each row has two entries:
  # the first power is applied to x, and the second power is applied to acd(x).
  combs <- generate_powers_acd(
    degree = degree,
    powers = powers
  )

  nfp <- nrow(combs)
  use_catzero <- !is.null(catzero)

  # Resolve the base ACD transformation once. During normal mfp2 training,
  # transform_data_step() passes fit_acd()$acd explicitly, which avoids
  # recalculating A(x) for the same training observations. Calls that do not
  # supply that cache keep the historical fit/apply semantics.
  x_acd_base <- resolve_acd_base_values(
    x = x,
    powers = powers,
    zero = zero,
    acd_parameter = acd_parameter,
    acd_training_values = acd_training_values
  )

  # Cache each unique FP transformation of x and acd(x).
  #
  # For default ACD degree 2, combs has 64 rows, but only 8 unique powers are
  # applied to x and 8 unique powers are applied to acd(x). Caching avoids
  # repeatedly calling transform_vector_fp() for the same power.
  make_power_key <- function(p) {
    if (is.na(p)) {
      return("<NA>")
    }

    as.character(p)
  }

  power_x_values <- unique(combs[, 1L])
  power_acd_values <- unique(combs[, 2L])

  x_fp_cache <- vector("list", length(power_x_values))
  names(x_fp_cache) <- vapply(power_x_values, make_power_key, character(1L))

  for (j in seq_along(power_x_values)) {
    p <- power_x_values[j]

    if (is.na(p)) {
      x_fp_cache[[j]] <- NULL
    } else {
      x_fp_cache[[j]] <- transform_vector_fp(
        x = x,
        power = p,
        scale = 1,
        shift = 0,
        zero = zero
      )
    }
  }

  x_acd_cache <- vector("list", length(power_acd_values))
  names(x_acd_cache) <- vapply(power_acd_values, make_power_key, character(1L))

  for (j in seq_along(power_acd_values)) {
    p <- power_acd_values[j]

    if (is.na(p)) {
      x_acd_cache[[j]] <- NULL
    } else {
      x_acd_cache[[j]] <- transform_vector_fp(
        x = x_acd_base,
        power = p,
        scale = 1,
        shift = 0,
        zero = FALSE
      )
    }
  }

  # Assemble one candidate matrix per ACD power row.
  fpdt <- vector("list", nfp)

  for (i in seq_len(nfp)) {
    x_fp <- x_fp_cache[[make_power_key(combs[i, 1L])]]
    x_acd <- x_acd_cache[[make_power_key(combs[i, 2L])]]

    mat <- cbind(x_fp, x_acd)

    if (use_catzero) {
      mat <- cbind(catzero, mat)
      colnames(mat) <- c("catzero", paste0("V", seq_len(ncol(mat) - 1L)))
    }

    fpdt[[i]] <- mat
  }

  list(
    data = fpdt,
    powers = combs
  )
}
