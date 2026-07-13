#' Predict Method for `mfp2`
#'
#' Obtain predictions from a fitted `mfp2` object.
#'
#' @description
#' This method predicts from the final fitted model stored in an `mfp2` object.
#' For ordinary predictions it delegates to the underlying `predict.glm()` or
#' `predict.coxph()` method after reconstructing the transformed design matrix
#' used during fitting. For `type = "terms"` and `type = "contrasts"`, it returns
#' variable-specific partial predictors with conditional standard errors.
#'
#' @details
#' If `newdata` is supplied, the variables are prepared using the transformation
#' metadata saved in `object`: fitted shifts, scales, selected fractional
#' polynomial powers, centering constants, zero-handling indicators, catzero
#' indicators, and spike-at-zero decisions. This ensures that prediction uses the
#' same transformed design scale as model fitting.
#'
#' Prediction stops when shifted `newdata` is outside the domain required by the
#' final fitted transformation. This check is applied only to variables whose
#' selected powers require strictly positive input, such as logarithmic powers,
#' repeated powers, negative powers, or non-integer powers. Variables whose final
#' fitted zero-handled or binary-only SAZ representation permits non-positive
#' values are exempt from this positivity stop.
#'
#' For full-model predictions, standard errors are requested with `se.fit = TRUE`
#' and are computed by the underlying `predict.glm()` or `predict.coxph()` method.
#' For Cox models, this method always calls `predict.coxph()` with
#' `reference = "zero"`. The transformed design matrix in `mfp2` has already been
#' centered where required, so asking `predict.coxph()` to subtract an additional
#' sample or strata reference would shift the covariate contribution incorrectly.
#' Cox offsets retain the native `predict.coxph()` convention: the prediction
#' offset is expressed relative to the mean offset in the fitting data.
#'
#' For `type = "terms"`, standard errors are computed from the fitted covariance
#' matrix for the columns belonging to each selected variable. If
#' `add_intercept = TRUE` and the model is not Cox, the intercept contribution is
#' included in both the term value and its standard error. If
#' `add_intercept = FALSE`, the intercept contribution is excluded from both.
#'
#' For `type = "contrasts"`, the returned value is the difference between the
#' partial predictor at each value and the partial predictor at the reference
#' value. Intercepts do not contribute to contrasts, and the standard errors are
#' computed from the transformed difference. These standard errors are
#' conditional on the final selected model and do not include model-selection
#' uncertainty.
#'
#' @section Terms prediction:
#' If `type = "terms"`, this function computes partial linear predictors for
#' selected terms in the final model. A conceptual term may be represented by
#' multiple fitted columns, for example an FP2 basis, a zero-handled positive
#' component with an indicator, or a categorical contrast block. The method
#' collects the relevant fitted columns for each term and multiplies them by their
#' coefficients.
#'
#' @section Grouped categorical terms:
#' Formula-created unordered and ordered factors are reconstructed with the
#' fitted factor levels and contrasts. Ordinary predictions therefore accept
#' the original factor variables in \code{newdata}; users do not need to create
#' dummy or polynomial-contrast columns manually.
#'
#' For \code{type = "terms"}, a grouped categorical term is evaluated at the
#' observed rows of the supplied \code{newdata}, or at the fitted rows when
#' \code{newdata = NULL}. The returned data frame includes the raw design-block
#' columns and, for simple formula factor main effects, the corresponding factor
#' level in \code{variable}. The \code{terms_seq} argument applies only to
#' singleton continuous terms.
#'
#' For \code{type = "contrasts"}, a formula factor reference may be supplied as
#' a level label, for example \code{ref = list(race = "A")}. If omitted, the
#' first fitted factor level is used. For a manually declared matrix-interface
#' group, the reference may be a numeric vector with one value per raw group
#' column; if omitted, the first observed design row is used.
#'
#' @section Contrasts:
#' If `type = "contrasts"`, this function computes variable-specific contrasts
#' relative to reference values. Reference values supplied through `ref` must be
#' on the original data scale; they are shifted, transformed, and centered using
#' the fitted transformation metadata. If `ref = NULL` for a variable, the method
#' uses the mean shifted value for continuous variables and the minimum shifted
#' value for binary-like variables.
#'
#' @param object A fitted object of class `mfp2`.
#' @param newdata Optional matrix or data frame containing variables for
#'   prediction. Column names must identify the original predictors used by the
#'   fit. Formula-created predictors are reconstructed when possible.
#' @param type Prediction type. The default is `"link"` for GLM models and
#'   `"lp"` for Cox models. Use `"terms"` for variable-specific partial
#'   predictors and `"contrasts"` for variable-specific contrasts.
#' @param se.fit Logical scalar. For full-model predictions, if `TRUE`, request
#'   standard errors from `predict.glm()` or `predict.coxph()`. For
#'   `type = "terms"` and `type = "contrasts"`, the returned data frames always
#'   contain an `se` column, and `se.fit` is ignored.
#' @param terms Character vector of conceptual term names for term or contrast
#'   prediction. For a simple wrapper such as \code{factor(x9)}, use the source
#'   variable name \code{x9}, not the wrapper expression or an individual
#'   generated contrast column. Only terms selected in the final fitted model
#'   are used. If `NULL`, all selected terms are used.
#' @param terms_seq Character scalar controlling values used for singleton
#'   continuous-term prediction. `"equidistant"` generates `nseq` equally
#'   spaced values over the observed range and `"data"` uses observed values.
#'   Grouped categorical terms are always evaluated at observed prediction rows.
#' @param alpha Significance level used for confidence intervals in term and
#'   contrast predictions.
#' @param ref Named list of reference values for `type = "contrasts"`.
#'   Continuous-term values are supplied on the original variable scale. A
#'   formula factor reference may be a fitted level label. A manually grouped
#'   matrix term may use a numeric vector with one value per raw member column.
#' @param strata Optional stratum values used when predicting from Cox models
#'   with supplied \code{newdata}. This is needed for default/matrix-interface
#'   Cox fits that used the fit-time \code{strata} argument, because the final
#'   Cox model contains an internal \code{strata(strata_)} term. Supply one
#'   value per prediction row for a single stratification factor, or a matrix or
#'   data frame with one row per prediction row for multiple stratification
#'   factors. Ordinary vector or factor strata are passed through as raw
#'   high-level values so that \code{predict.coxph()} can evaluate the stored
#'   \code{strata(strata_)} term itself. Formula-interface strata are
#'   reconstructed automatically from \code{newdata} when the original strata
#'   variable(s) are present.
#' @param newoffset Optional numeric vector of offsets for prediction when the
#'   fitted model used an offset and `newdata` is supplied.
#' @param nseq Positive integer giving the number of equally spaced values used
#'   when `terms_seq = "equidistant"`.
#' @param add_intercept Logical scalar. For `type = "terms"`, controls whether
#'   the model intercept is included in GLM term values and term standard errors.
#'   It has no effect for contrasts and does not apply to Cox models.
#' @param ... Further arguments passed to `predict.glm()` or `predict.coxph()`
#'   for full-model predictions.
#'
#' @return
#' For full-model predictions, the return value follows `predict.glm()` or
#' `predict.coxph()`. In particular, if `se.fit = TRUE`, the result may be a list
#' containing fitted values and standard errors according to the underlying
#' method.
#'
#' For `type = "terms"` or `type = "contrasts"`, the result is a named list with
#' one data frame per selected conceptual term. Each data frame contains:
#' \itemize{
#'   \item `variable`: values on the original scale before fitted shifting.
#'     For simple formula factor terms, this contains factor-level labels.
#'   \item `variable_pre`: values after fitted shifting, before or after binary
#'     SAZ coding as appropriate. For factor terms, this repeats the level label.
#'   \item Raw grouped-term columns: categorical or manually grouped terms also
#'     include the underlying model-matrix columns used to compute the term.
#'   \item `value`: partial linear predictor or contrast.
#'   \item `se`: conditional standard error.
#'   \item `lower`: lower confidence limit.
#'   \item `upper`: upper confidence limit.
#' }
#'
#' @examples
#' data("prostate")
#' x <- as.matrix(prostate[, 2:8])
#' y <- as.numeric(prostate$lpsa)
#'
#' fit <- mfp2(x, y, verbose = FALSE)
#' predict(fit)
#' predict(fit, se.fit = TRUE)
#' predict(fit, type = "terms")
#'
#' \dontrun{
#' # Formula factors are supplied in their original form for prediction.
#' set.seed(1)
#' d <- data.frame(
#'   y = rnorm(45),
#'   age = runif(45, 30, 75),
#'   race = factor(rep(c("A", "B", "C"), each = 15)),
#'   severity = ordered(
#'     rep(c("mild", "moderate", "severe"), length.out = 45),
#'     levels = c("mild", "moderate", "severe")
#'   )
#' )
#' fit_factor <- mfp2(
#'   y ~ age + race + severity, data = d,
#'   keep = c("race", "severity"), verbose = FALSE
#' )
#' predict(fit_factor, newdata = d[1:4, ])
#' predict(fit_factor, newdata = d[1:4, ], type = "terms", terms = "race")
#' predict(
#'   fit_factor,
#'   newdata = d[1:4, ],
#'   type = "contrasts",
#'   terms = "race",
#'   ref = list(race = "A")
#' )
#'
#' }
#'
#' @seealso
#' \code{mfp2()}, [stats::predict.glm()], [survival::predict.coxph()]
#'
#' @method predict mfp2
#' @export
predict.mfp2 <- function(object, 
                         newdata = NULL, 
                         type = NULL,
                         se.fit = FALSE,
                         terms = NULL,
                         terms_seq = c("equidistant", "data"),
                         alpha = 0.05,
                         ref = NULL, 
                         strata = NULL, 
                         newoffset = NULL, 
                         nseq = 100,
                         add_intercept = TRUE,
                         ...) {
  
  terms_seq <- match.arg(terms_seq)
  
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be between 0 and 1.", call. = FALSE)
  }
  
  if (!is.numeric(nseq) || nseq <= 0) {
    stop("'nseq' must be a positive integer.", call. = FALSE)
  }
  
  if (!is.logical(se.fit) || length(se.fit) != 1L || is.na(se.fit)) {
    stop("'se.fit' must be a single TRUE or FALSE value.", call. = FALSE)
  }
  
  # Prediction metadata and fitted-model methods are defined only for mfp2 objects.
  if (!inherits(object, "mfp2")) { 
    stop("The object is not an mfp2 object.", call. = FALSE)
  }
  
  # Resolve the prediction scale. Formula fits store the same conceptual term
  # names used by selection and printing, including source-variable names for
  # simple wrappers such as factor(x9).
  if (is.null(type)) {
    type <- ifelse(object$family_string == "cox", "lp", "link")
  }
  
  selected_internal_terms <- get_selected_variable_names(object)
  
  if (is.null(terms)) {
    terms <- selected_internal_terms
  }
  
  if (is.null(ref)) {
    ref <- stats::setNames(lapply(terms, function(v) NULL), terms)
  }
  
  # Preserve the user's original formula-style data for factor labels, offsets,
  # and strata. Formula fits are then rebuilt on their fitted model-matrix scale.
  if (!is.null(newdata) && is.null(colnames(newdata))) {
    stop("Newdata must have column names", call. = FALSE)
  }
  
  newdata_raw <- newdata
  
  if (!is.null(newdata)) {
    newdata <- reconstruct_formula_newdata(object, newdata)
  }
  
  if (!is.null(newdata) && anyNA(newdata)) {
    stop("! newdata must not contain any NA (missing data).\n", 
         "i Please remove any missing data before passing newdata to this function.",
         call. = FALSE)
  }
  
  if (type == "contrasts" && length(ref) != sum(names(ref) != "", na.rm = TRUE)) {
    warning(
      paste0(
        "i The supplied reference values (ref) must all be named.\n",
        "i predict() continues with the default reference: the mean for ",
        "continuous terms, the minimum for binary terms, or the first fitted ",
        "level for formula factors."
      ),
      call. = FALSE
    )
  }
  
  # Term and contrast prediction is handled term by term. A conceptual term
  # may correspond to one transformed continuous predictor or to a block of
  # dummy/contrast columns for a categorical term.
  if (type %in% c("terms", "contrasts")) {
    requested_terms <- terms
    internal_terms <- requested_terms
    
    selected <- internal_terms %in% selected_internal_terms
    if (!any(selected)) {
      warning(
        "i All the terms supplied are not in the final model.\n",
        "i predict() continues but returns an empty list.",
        call. = FALSE
      )
    } else if (!all(selected)) {
      warning(
        paste0(
          "i Some terms supplied are not in the final model.\n",
          "i predict() continues but omits terms not in the model."
        ),
        call. = FALSE
      )
    }
    
    requested_terms <- requested_terms[selected]
    internal_terms <- internal_terms[selected]
    
    valid_reference_names <- requested_terms
    if (!all(vapply(ref, is.null, logical(1L))) &&
        any(!names(ref) %in% valid_reference_names)) {
      warning(
        "i Some reference names are not selected prediction terms and are ignored.",
        call. = FALSE
      )
    }
    
    # The intercept contributes to GLM term predictions only when requested.
    # It cancels from contrasts and is not present in Cox partial predictors.
    cf <- coef(object)
    intercept <- if ("(Intercept)" %in% names(cf)) cf[["(Intercept)"]] else 0
    if (is.na(intercept) || !add_intercept || type == "contrasts") {
      intercept <- 0
    }
    
    # Resolve conceptual prediction terms once. Each element names the raw
    # design-matrix columns that must be transformed and combined for that term.
    lookup <- prediction_term_to_columns(object)
    res_list <- list()
    
    for (term_index in seq_along(internal_terms)) {
      t <- internal_terms[[term_index]]
      display_term <- requested_terms[[term_index]]
      term_columns <- lookup[[t]]
      factor_info <- if (!is.null(object$formula_factor_info)) {
        object$formula_factor_info[[t]]
      } else {
        NULL
      }
      block_term <- term_uses_column_mapping(t, term_columns) || !is.null(factor_info)
      
      if (block_term) {
        # Grouped terms are predicted as complete design blocks. This includes
        # multi-column factors and one-column formula factors carrying factor
        # metadata for user-facing level labels.
        # Reconstructed `newdata` is on the raw model-matrix scale and still
        # needs the fitted shift. In contrast, `object$x_original` is stored on
        # the shifted-but-not-scaled final-fitting scale, so its shift must not
        # be applied a second time.
        source_needs_preprocessing <- !is.null(newdata)
        
        source_block <- if (!is.null(newdata)) {
          missing_columns <- setdiff(term_columns, colnames(newdata))
          if (length(missing_columns) > 0L) {
            stop(
              "Prediction data are missing grouped-term column(s): ",
              paste(missing_columns, collapse = ", "),
              call. = FALSE
            )
          }
          newdata[, term_columns, drop = FALSE]
        } else {
          missing_columns <- setdiff(term_columns, colnames(object$x_original))
          if (length(missing_columns) > 0L) {
            stop(
              "Fitted data do not contain grouped-term column(s): ",
              paste(missing_columns, collapse = ", "),
              call. = FALSE
            )
          }
          object$x_original[, term_columns, drop = FALSE]
        }
        
        source_block <- as.matrix(source_block)
        x_trafo <- as.matrix(prepare_newdata_for_predict(
          object,
          source_block,
          apply_pre = source_needs_preprocessing,
          allow_missing_predictors = TRUE,
          check_binary = FALSE
        ))
        
        model_columns <- prediction_model_column_names(
          object,
          colnames(x_trafo)
        )
        term_coef <- stats::coef(object)[model_columns]
        term_coef[is.na(term_coef)] <- 0
        
        labels <- decode_grouped_term_values(
          object,
          term = t,
          block = source_block
        )
        
        res <- data.frame(
          variable = labels,
          variable_pre = labels,
          source_block,
          value = as.numeric(x_trafo %*% term_coef + intercept),
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
        
        x_ref_trafo <- NULL
        if (type == "contrasts") {
          reference_block <- resolve_grouped_reference(
            object,
            term = t,
            reference = if (display_term %in% names(ref)) ref[[display_term]] else ref[[t]],
            block = source_block
          )
          x_ref_trafo <- as.matrix(prepare_newdata_for_predict(
            object,
            reference_block,
            apply_pre = TRUE,
            allow_missing_predictors = TRUE,
            check_binary = FALSE
          ))
          res$value <- res$value - as.numeric(x_ref_trafo %*% term_coef)
        }
      } else {
        # Singleton numeric terms retain the ordinary FP prediction path.
        # Evaluation values are represented on the shifted, unscaled scale used
        # by the final fitted transformation.
        if (terms_seq == "equidistant") {
          if (!is.null(newdata)) {
            x_range <- range(newdata[, t]) + object$transformations[t, "shift"]
          } else {
            x_range <- range(object$x_original[, t])
          }
          
          x_seq <- matrix(
            seq(x_range[1L], x_range[2L], length.out = nseq),
            ncol = 1L,
            dimnames = list(NULL, t)
          )
          
          x_trafo <- as.matrix(prepare_newdata_for_predict(
            object,
            x_seq,
            apply_pre = FALSE,
            allow_missing_predictors = TRUE
          ))
        } else {
          x_seq <- if (!is.null(newdata)) {
            newdata[, t, drop = FALSE] + object$transformations[t, "shift"]
          } else {
            object$x_original[, t, drop = FALSE]
          }
          
          x_trafo <- as.matrix(prepare_newdata_for_predict(
            object,
            x_seq,
            apply_pre = FALSE,
            allow_missing_predictors = TRUE
          ))
        }
        
        model_columns <- prediction_model_column_names(
          object,
          colnames(x_trafo)
        )
        term_coef <- stats::coef(object)[model_columns]
        term_coef[is.na(term_coef)] <- 0
        
        variable <- as.numeric(x_seq) - object$transformations[t, "shift"]
        variable_pre <- as.numeric(x_seq)
        
        if (object$spike_dec[t] == saz_decision_codes[["binary_only"]]) {
          variable <- as.integer(x_seq[, t] <= 0)
          variable_pre <- variable
        }
        
        res <- data.frame(
          variable = variable,
          variable_pre = variable_pre,
          value = as.numeric(x_trafo %*% term_coef + intercept)
        )
        
        x_ref_trafo <- NULL
        if (type == "contrasts") {
          x_ref <- if (display_term %in% names(ref)) ref[[display_term]] else ref[[t]]
          if (is.null(x_ref)) {
            v <- if (!is.null(newdata)) {
              newdata[, t] + object$transformations[t, "shift"]
            } else {
              object$x_original[, t]
            }
            x_ref <- if (length(unique(stats::na.omit(v))) == 2L) {
              min(v, na.rm = TRUE)
            } else {
              mean(v, na.rm = TRUE)
            }
          } else {
            x_ref <- x_ref + object$transformations[t, "shift"]
          }
          
          x_ref <- matrix(x_ref, nrow = 1L, ncol = 1L, dimnames = list(NULL, t))
          x_ref_trafo <- as.matrix(prepare_newdata_for_predict(
            object,
            x_ref,
            apply_pre = FALSE,
            check_binary = FALSE,
            reset_zero = FALSE,
            allow_missing_predictors = TRUE
          ))
          res$value <- res$value - as.numeric(x_ref_trafo %*% term_coef)
        }
      }
      
      # Standard errors use the coefficient covariance block corresponding to
      # the transformed columns for this conceptual term.
      res$se <- calculate_standard_error(
        object,
        x_trafo,
        x_ref_trafo,
        include_intercept = add_intercept && type == "terms"
      )
      mult <- stats::qnorm(1 - alpha / 2)
      res$lower <- res$value - mult * res$se
      res$upper <- res$value + mult * res$se
      res_list[[display_term]] <- res
    }
    
    return(res_list)
  }
  
  # Full-model prediction reconstructs every required raw design column,
  # applies the fitted transformations, and then delegates to the underlying
  # glm or coxph prediction method.
  
  if (!is.null(newdata)) {
    
    if (
      is.null(newoffset) &&
      !is.null(newdata_raw) &&
      !is.null(object$formula_offset_terms)
    ) {
      newoffset <- reconstruct_formula_offset_newdata(object, newdata_raw)
    }
    
    # Formula offsets are reconstructed automatically when possible; otherwise
    # an offset must be supplied explicitly for models fitted with one.
    has_offset <- isTRUE(object$has_offset)
    if (has_offset) {
      if (is.null(newoffset)) {
        stop(
          "No newoffset provided for prediction, yet offset was used in mfp2",
          call. = FALSE
        )
      }
      
      if (!is.numeric(newoffset)) {
        stop("! newoffset must be numeric.", call. = FALSE)
      }
      
      if (length(newoffset) != nrow(newdata)) {
        stop(
          "The length of newoffset must be equal to the number of rows of newdata",
          call. = FALSE
        )
      }
      
    } else {
      newoffset <- NULL
    }
    
    if (
      object$family_string == "cox" &&
      is.null(strata) &&
      !is.null(newdata_raw) &&
      !is.null(object$formula_strata_terms)
    ) {
      strata <- reconstruct_formula_strata_newdata(object, newdata_raw)
    }
    
    newdata <- prepare_newdata_for_predict(
      object,
      newdata,
      strata = strata,
      offset = newoffset,
      check_binary = FALSE
    )
    
    # Strip "mfp2" from the class vector so stats::predict() dispatches to
    # predict.coxph() or predict.glm() instead of recursing into predict.mfp2().
    obj_base <- object
    class(obj_base) <- setdiff(class(obj_base), "mfp2")
    
    if (object$family_string == "cox") {
      # Cox reference convention: mfp2 stores and predicts on an already
      # prepared design scale. Use reference = "zero" so predict.coxph() does
      # not subtract its own sample or stratum means a second time. Raw
      # prediction offsets are passed through unchanged: predict.coxph() itself
      # subtracts the stored mean training offset. Validate the fit-time origin
      # whenever an offset was structurally included.
      if (isTRUE(object$has_offset)) {
        mfp2_cox_offset_reference(object)
      }
      
      pred <- stats::predict(
        obj_base,
        newdata = newdata,
        type = type,
        se.fit = se.fit,
        reference = "zero",
        ...
      )
    } else {
      pred <- stats::predict(
        obj_base,
        newdata = newdata,
        type = type,
        se.fit = se.fit,
        ...
      )
    }
    
    return(pred)
  }
  
  # With no newdata, delegate directly to the fitted base model.
  obj_base <- object
  class(obj_base) <- setdiff(class(obj_base), "mfp2")
  if (object$family_string == "cox") {
    if (isTRUE(object$has_offset)) {
      mfp2_cox_offset_reference(object)
    }
    
    return(
      stats::predict(
        obj_base,
        type = type,
        se.fit = se.fit,
        reference = "zero",
        ...
      )
    )
  } else {
    return(
      stats::predict(
        obj_base,
        type = type,
        se.fit = se.fit,
        ...
      )
    )
  }
}

#' Resolve Conceptual Terms to Raw Prediction Columns
#'
#' Returns the complete term lookup stored by grouped-term fits. The names are
#' conceptual terms and each value contains the corresponding raw design-matrix
#' columns in fitted term order.
#'
#' @param object A fitted \code{mfp2} object.
#'
#' @return Named list mapping conceptual terms to raw input columns.
#'
#' @keywords internal
#' @noRd
prediction_term_to_columns <- function(object) {
  terms <- rownames(object$transformations)
  
  if (is.null(terms) || length(terms) == 0L) {
    stop("The fitted object does not contain transformation term names.",
         call. = FALSE)
  }
  
  lookup <- object$term_to_columns
  
  if (is.null(lookup)) {
    stop("The fitted object does not contain `term_to_columns`.", call. = FALSE)
  }
  
  if (!is.list(lookup) || is.null(names(lookup)) || anyDuplicated(names(lookup))) {
    stop("The fitted grouped-term lookup is malformed.", call. = FALSE)
  }
  
  missing_terms <- setdiff(terms, names(lookup))
  if (length(missing_terms) > 0L) {
    stop(
      "The fitted grouped-term lookup is missing term(s): ",
      paste(missing_terms, collapse = ", "),
      call. = FALSE
    )
  }
  
  lookup <- lookup[terms]
  
  invalid <- vapply(
    lookup,
    function(columns) {
      !is.character(columns) || length(columns) == 0L || anyNA(columns) ||
        any(!nzchar(columns))
    },
    logical(1L)
  )
  
  if (any(invalid) || anyDuplicated(unlist(lookup, use.names = FALSE))) {
    stop("The fitted grouped-term lookup contains invalid raw columns.",
         call. = FALSE)
  }
  
  lookup
}

# Resolve transformed design columns to the exact coefficient names stored by
# the final formula-based fit. The map is created once in fit_mfp() from the
# transformed matrix and fitted coefficient order; prediction therefore does
# not need to quote, unquote, sanitize, or otherwise guess model column names.
prediction_model_column_names <- function(object, transformed_columns) {
  model_columns <- unname(
    object$transformed_to_model_columns[transformed_columns]
  )
  
  if (is.null(object$transformed_to_model_columns) || anyNA(model_columns)) {
    stop(
      "Internal error: fitted object lacks the transformed-to-model column mapping.",
      call. = FALSE
    )
  }
  
  model_columns
}



#' Expand Term-Level Metadata for Raw Prediction Columns
#'
#' Grouped fits store selection and transformation metadata once per conceptual
#' term, while model fitting and prediction use the underlying raw design
#' columns. This helper recreates the raw-column metadata needed by
#' \code{transform_matrix()}.
#'
#' @param object A fitted \code{mfp2} object.
#' @param raw_columns Character vector of raw prediction columns.
#'
#' @return List of raw-column metadata aligned with \code{raw_columns}.
#'
#' @keywords internal
#' @noRd
expand_prediction_metadata <- function(object, raw_columns) {
  # Resolve the fitted term structure and build its reverse mapping so each raw
  # design column can inherit the metadata stored for its conceptual term.
  lookup <- prediction_term_to_columns(object)
  column_to_term <- stats::setNames(
    rep(names(lookup), lengths(lookup)),
    unlist(lookup, use.names = FALSE)
  )
  
  unknown <- setdiff(raw_columns, names(column_to_term))
  if (length(unknown) > 0L) {
    stop(
      "Unknown raw prediction column(s): ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }
  
  raw_terms <- unname(column_to_term[raw_columns])
  
  # Determine grouped status through a named term-level vector. This avoids
  # repeatedly subsetting the lookup with duplicated term names when several
  # raw columns belong to the same conceptual term.
  grouped_by_term <- mapped_term_flags(lookup)
  grouped <- unname(grouped_by_term[raw_terms])
  
  # Grouped categorical columns are fixed linear columns. A selected group
  # therefore expands to power 1 for every member; an omitted group expands to NA.
  powers <- lapply(seq_along(raw_columns), function(i) {
    term <- raw_terms[[i]]
    term_powers <- object$fp_powers[[term]]
    
    if (grouped[[i]]) {
      if (is.null(term_powers) || all(is.na(term_powers))) NA_real_ else 1
    } else {
      term_powers
    }
  })
  names(powers) <- raw_columns
  
  term_value <- function(values, term, default = FALSE) {
    if (is.null(values) || !term %in% names(values)) default else values[[term]]
  }
  
  acdx <- zero <- catzero <- spike <- logical(length(raw_columns))
  spike_decision <- integer(length(raw_columns))
  acd_parameter <- vector("list", length(raw_columns))
  
  # Continuous-only options are disabled for grouped categorical blocks.
  # Singleton terms inherit their stored ACD, zero, catzero, and spike settings.
  for (i in seq_along(raw_columns)) {
    term <- raw_terms[[i]]
    
    if (grouped[[i]]) {
      acdx[[i]] <- FALSE
      zero[[i]] <- FALSE
      catzero[[i]] <- FALSE
      spike[[i]] <- FALSE
      spike_decision[[i]] <- saz_decision_codes[["continuous_only"]]
      acd_parameter[i] <- list(NULL)
    } else {
      acdx[[i]] <- isTRUE(term_value(object$acd, term, FALSE))
      zero[[i]] <- isTRUE(term_value(object$zero, term, FALSE))
      catzero[[i]] <- isTRUE(term_value(object$catzero, term, FALSE))
      spike_source <- if (!is.null(object$spike)) object$spike else object$fp_terms[, "spike"]
      spike[[i]] <- isTRUE(term_value(spike_source, term, FALSE))
      spike_decision[[i]] <- as.integer(
        term_value(object$spike_dec, term, saz_decision_codes[["continuous_only"]])
      )
      # Use single-bracket list assignment so a missing ACD parameter is
      # stored as a NULL element without deleting the element from the list.
      acd_parameter[i] <- list(
        if (!is.null(object$acd_parameter) &&
            term %in% names(object$acd_parameter)) {
          object$acd_parameter[[term]]
        } else {
          NULL
        }
      )
    }
  }
  
  names(acdx) <- names(zero) <- names(catzero) <- names(spike) <-
    names(spike_decision) <- names(acd_parameter) <- raw_columns
  
  shifts <- stats::setNames(
    vapply(
      raw_terms,
      function(term) as.numeric(object$transformations[term, "shift"]),
      numeric(1L)
    ),
    raw_columns
  )
  
  list(
    terms = stats::setNames(raw_terms, raw_columns),
    powers = powers,
    shift = shifts,
    acdx = acdx,
    zero = zero,
    catzero = catzero,
    spike = spike,
    spike_decision = spike_decision,
    acd_parameter = acd_parameter
  )
}


#' Decode Formula Factor Levels from a Raw Design Block
#'
#' @param object A fitted \code{mfp2} object.
#' @param term Conceptual factor term name.
#' @param block Raw design block for the term.
#' @return Character vector of level labels or design signatures.
#'
#' @keywords internal
#' @noRd
decode_grouped_term_values <- function(object, term, block) {
  # Match each reconstructed design row to the fitted level-to-design mapping.
  # This remains correct when an inline factor() call changes labels, level
  # order, or contrast coding relative to the raw source variable.
  info <- if (!is.null(object$formula_factor_info)) {
    object$formula_factor_info[[term]]
  } else {
    NULL
  }
  
  if (!is.null(info) && !is.null(info$design_by_level)) {
    design <- as.matrix(info$design_by_level)
    out <- rep(NA_character_, nrow(block))
    
    for (i in seq_len(nrow(block))) {
      matches <- which(vapply(
        seq_len(nrow(design)),
        function(j) {
          all(is.finite(design[j, ])) &&
            isTRUE(all.equal(
              as.numeric(block[i, ]),
              as.numeric(design[j, ]),
              tolerance = 1e-10,
              check.attributes = FALSE
            ))
        },
        logical(1L)
      ))
      
      if (length(matches) > 0L) {
        out[[i]] <- rownames(design)[matches[[1L]]]
      }
    }
    
    if (all(!is.na(out))) {
      return(out)
    }
  }
  
  # Non-factor grouped terms have no level labels; expose their design row as
  # a compact signature so the returned value remains identifiable.
  apply(
    block,
    1L,
    function(row) paste(format(row, trim = TRUE), collapse = ",")
  )
}


#' Resolve a Reference Design Row for a Grouped Term
#'
#' @param object A fitted \code{mfp2} object.
#' @param term Conceptual term name.
#' @param reference User-supplied reference value or \code{NULL}.
#' @param block Observed raw design block used as a fallback reference source.
#'
#' @return One-row numeric matrix with columns matching \code{block}.
#'
#' @keywords internal
#' @noRd
resolve_grouped_reference <- function(object, term, reference, block) {
  # Formula factors accept level labels. Manually grouped matrix terms accept a
  # complete numeric design row with one value per raw member column.
  columns <- colnames(block)
  info <- if (!is.null(object$formula_factor_info)) {
    object$formula_factor_info[[term]]
  } else {
    NULL
  }
  
  if (!is.null(info)) {
    level <- if (is.null(reference)) info$levels[[1L]] else as.character(reference)
    
    if (length(level) == 1L && level %in% rownames(info$design_by_level)) {
      values <- info$design_by_level[level, columns, drop = FALSE]
      if (anyNA(values)) {
        stop(
          "No fitted design row is available for reference level '", level,
          "' of term '", term, "'.",
          call. = FALSE
        )
      }
      return(values)
    }
    
    if (!is.numeric(reference)) {
      stop(
        "Reference for factor term '", term, "' must be one of: ",
        paste(info$levels, collapse = ", "),
        call. = FALSE
      )
    }
  }
  
  # For manually grouped terms, default to the first observed complete design
  # row when no explicit reference is supplied.
  if (is.null(reference)) {
    reference <- block[1L, , drop = FALSE]
  }
  
  if (is.data.frame(reference) || is.matrix(reference)) {
    reference <- as.matrix(reference)
    if (nrow(reference) != 1L) {
      stop("Grouped-term reference must contain exactly one row.", call. = FALSE)
    }
    if (!is.null(colnames(reference))) {
      missing <- setdiff(columns, colnames(reference))
      if (length(missing) > 0L) {
        stop(
          "Grouped-term reference is missing column(s): ",
          paste(missing, collapse = ", "),
          call. = FALSE
        )
      }
      reference <- reference[, columns, drop = FALSE]
    }
  } else {
    reference_names <- names(reference)
    reference <- as.numeric(reference)
    
    if (!is.null(reference_names) && all(nzchar(reference_names))) {
      missing <- setdiff(columns, reference_names)
      if (length(missing) > 0L) {
        stop(
          "Grouped-term reference is missing column(s): ",
          paste(missing, collapse = ", "),
          call. = FALSE
        )
      }
      names(reference) <- reference_names
      reference <- reference[columns]
    }
    
    if (length(reference) != length(columns)) {
      stop(
        "Reference for grouped term '", term, "' must contain ",
        length(columns), " value(s), one per raw column.",
        call. = FALSE
      )
    }
    reference <- matrix(reference, nrow = 1L, dimnames = list(NULL, columns))
  }
  
  storage.mode(reference) <- "double"
  reference
}


#' Transform Linear Predictor to Response or Risk
#' 
#' Converts linear predictors (`nfit`) from a model to the appropriate scale
#' for interpretation or prediction, depending on the model family, link function,
#' and type of prediction. This is an internal helper function for GLMs, survival models,
#' used by mfp2 prediction helpers.
#' 
#' @param nfit Numeric vector of linear predictors \eqn{X\beta}.
#' @param family Character string specifying the model family. Supported families include:
#'   - `"gaussian"`: linear regression
#'   - `"binomial"`: binary outcomes or proportions
#'   - `"poisson"`: count data
#'   - `"cox"`: proportional hazards model
#'   - other families fallback to returning `nfit`
#' @param link Character string specifying the link function. Defaults:
#'   - Gaussian: `"identity"`  
#'   - Binomial: `"logit"`  
#'   - Poisson: `"log"`  
#'   Supported links:
#'   - Gaussian: `"identity"`, `"log"`, `"inverse"`  
#'   - Binomial: `"logit"`, `"probit"`, `"cloglog"`, `"cauchit"`, `"identity"`, `"log"`  
#'   - Poisson: `"log"`, `"identity"`, `"sqrt"`
#' @param type Character string specifying the type of prediction:
#'   - `"response"` or `"risk"`: returns predictions on the response/probability scale
#'   - `NULL` (default): returns the linear predictor itself
#' 
#' @return Numeric vector of predictions on the requested scale.
#' 
#' @examples
#' \dontrun{
#' # Binomial example
#' lp_bin <- c(-1, 0, 1)
#' transform_linear_predictor(lp_bin, family = "binomial", link = "logit", type = "response")
#' 
#' # Poisson example
#' lp_pois <- c(0.5, 1, 1.5)
#' transform_linear_predictor(lp_pois, family = "poisson", link = "log", type = "response")
#' 
#' # Gaussian example
#' lp_gauss <- c(1, 2, 3)
#' transform_linear_predictor(lp_gauss, family = "gaussian", link = "log", type = "response")
#' }
#' @keywords internal
#' @noRd
transform_linear_predictor <- function(nfit, family, link = NULL, type = NULL) {
  
  # Use the canonical link for supported families when the fitted link is absent.
  if (is.null(link)) {
    link <- switch(family,
                   gaussian = "identity",
                   binomial = "logit",
                   poisson  = "log",
                   NULL)
  }
  
  switch(family,
         
         # Gaussian family with multiple links
         gaussian = {
           if (!is.null(type) && type == "response") {
             switch(link,
                    identity = nfit,
                    log      = exp(nfit),
                    inverse  = 1 / nfit,
                    stop("Unknown Gaussian link"))
           } else nfit
         },
         
         # Binomial family with multiple links
         binomial = {
           if (!is.null(type) && type == "response") {
             switch(link,
                    logit   = 1 / (1 + exp(-nfit)),
                    probit  = pnorm(nfit),
                    cloglog = 1 - exp(-exp(nfit)),
                    cauchit = pcauchy(nfit),
                    identity = nfit,
                    log     = exp(nfit),
                    stop("Unknown binomial link"))
           } else nfit
         },
         
         # Poisson family with multiple links
         poisson = {
           if (!is.null(type) && type == "response") {
             switch(link,
                    log      = exp(nfit),
                    identity = nfit,
                    sqrt     = nfit^2,
                    stop("Unknown Poisson link"))
           } else nfit
         },
         
         # Cox proportional hazards: exponentiate for risk
         cox = {
           if (!is.null(type) && type %in% c("response", "risk")) exp(nfit) else nfit
         },
         
         # Unknown families remain on the linear-predictor scale.
         nfit
  )
}

#' Does a fitted term require positive raw prediction input?
#'
#' Internal helper used by `predict.mfp2()`.
#'
#' The prediction method applies the shift values learned during fitting before
#' constructing transformed variables. This helper determines whether the
#' shifted raw covariate for variable `v` must be strictly positive, based on
#' the final fitted model metadata.
#'
#' For ordinary FP terms, all selected powers in `object$fp_powers[[v]]` apply
#' directly to the raw covariate.
#'
#' For ACD terms, `object$fp_powers[[v]][1]` applies to the raw covariate and
#' `object$fp_powers[[v]][2]` applies to the ACD component A(x). Therefore only
#' the first power is checked against the raw covariate here.
#'
#' Variables fitted with zero-handling are exempt because nonpositive values are
#' represented structurally. Variables with `spike_decision = 3` are also exempt
#' because the continuous component is dropped and only the spike/catzero
#' indicator is retained.
#'
#' @param object A fitted `mfp2` object.
#' @param v Character scalar giving the conceptual term name.
#'
#' @return A single logical value.
#'
#' @keywords internal
#' @noRd
requires_positive_raw_input <- function(object, v) {
  powers <- object$fp_powers[[v]]
  
  if (is.null(powers) || all(is.na(powers))) {
    return(FALSE)
  }
  
  if (!is.null(object$zero) && isTRUE(object$zero[[v]])) {
    return(FALSE)
  }
  
  if (!is.null(object$spike_dec) &&
      v %in% names(object$spike_dec) &&
      !is.na(object$spike_dec[[v]]) &&
      as.integer(object$spike_dec[[v]]) == saz_decision_codes[["binary_only"]]) {
    return(FALSE)
  }
  
  is_acd <- !is.null(object$acd) &&
    v %in% names(object$acd) &&
    isTRUE(object$acd[[v]])
  
  raw_powers <- if (is_acd) {
    # For ACD terms:
    #   powers[1] applies to raw x
    #   powers[2] applies to A(x)
    #
    # Only the raw-x power determines whether shifted raw newdata must be
    # strictly positive at this stage.
    powers[1L]
  } else {
    powers
  }
  
  fp_power_requires_positive_input(raw_powers)
}

#' Rebuild formula-interface newdata with model.matrix()
#'
#' Formula-fitted mfp2 objects are fitted on the expanded numeric design matrix
#' returned by model.matrix(). This helper lets users pass ordinary newdata with
#' the original formula variables, including factors, and reconstructs the same
#' expanded columns used during fitting. If all fitted design columns are already
#' present, no formula reconstruction is necessary.
#'
#' @keywords internal
#' @noRd
reconstruct_formula_newdata <- function(object, newdata) {
  # Matrix-interface fits already receive raw design columns and need no formula
  # evaluation. Formula fits are reconstructed with stored levels and contrasts.
  if (!isTRUE(object$formula_interface)) {
    return(newdata)
  }
  
  expected <- unlist(prediction_term_to_columns(object), use.names = FALSE)
  # Accept callers that already supplied the exact fitted model-matrix columns.
  # Return only those columns, in fitted order. Formula newdata may also contain
  # factor-valued strata or other non-predictor variables. Keeping those extra
  # columns until prepare_newdata_for_predict() converts the data to a matrix
  # would coerce the complete predictor block to character. The original raw
  # newdata is retained separately by predict.mfp2() for reconstructing strata
  # and offsets, so it is correct to isolate the fitted numeric design here.
  if (length(expected) > 0L && all(expected %in% colnames(newdata))) {
    return(newdata[, expected, drop = FALSE])
  }
  
  if (is.null(object$formula_terms)) {
    stop(
      "! This object was fitted using the formula interface, but formula terms ",
      "needed for prediction were not stored.",
      call. = FALSE
    )
  }
  
  newdata_df <- as.data.frame(newdata)
  
  # model.frame() evaluates formula expressions and enforces the fitted factor
  # levels; model.matrix() then recreates the original contrast coding.
  mf <- stats::model.frame(
    object$formula_terms,
    data = newdata_df,
    na.action = stats::na.pass,
    xlev = object$formula_xlevels
  )
  
  mm <- stats::model.matrix(
    object$formula_terms,
    data = mf,
    contrasts.arg = object$formula_contrasts
  )
  
  keep_cols <- colnames(mm) != "(Intercept)"
  mm <- mm[, keep_cols, drop = FALSE]
  
  # Restore the exact column names used by the fitted object after any formula
  # sanitisation or internal renaming performed during fitting.
  column_map <- object$formula_column_map
  if (!is.null(column_map)) {
    mapped <- unname(column_map[colnames(mm)])
    colnames(mm) <- ifelse(is.na(mapped), colnames(mm), mapped)
  }
  
  # Keep and order only the predictor columns that were present in the fitted
  # formula design matrix.
  if (!is.null(object$formula_model_matrix_columns)) {
    missing <- setdiff(object$formula_model_matrix_columns, colnames(mm))
    if (length(missing) > 0L) {
      stop(
        "! Could not reconstruct required model-matrix column(s) from newdata: ",
        paste(missing, collapse = ", "),
        call. = FALSE
      )
    }
    mm <- mm[, object$formula_model_matrix_columns, drop = FALSE]
  }
  
  mm
}

#' Reconstruct formula-level Cox strata from prediction data
#'
#' Reconstructs the original formula-level Cox stratification term from
#' `newdata` for `mfp2` objects fitted through the formula interface.
#'
#' This helper deliberately mimics the high-level `survival::coxph()` strata
#' handling:
#'
#' - for a single strata variable, it returns the evaluated column from the
#'   model frame;
#' - for multiple strata variables, it combines them with
#'   `survival::strata(..., shortlabel = TRUE)`;
#' - it does **not** convert the result to integer codes.
#'
#' Integer conversion belongs only in the low-level Cox fitting path immediately
#' before calling `survival::coxph.fit()`, analogous to the `istrat` object used
#' internally by `survival::coxph()`.
#'
#' @param object A fitted `mfp2` object.
#' @param newdata A user-supplied prediction data frame before formula
#'   reconstruction has been applied.
#'
#' @return Either `NULL`, if no formula-level strata metadata are available, or
#'   a Cox strata object suitable for passing as the `strata` argument to
#'   `prepare_newdata_for_predict()`.
#'
#' @keywords internal
#' @noRd
reconstruct_formula_strata_newdata <- function(object, newdata) {
  if (!isTRUE(object$formula_interface) ||
      is.null(object$formula_strata_terms)) {
    return(NULL)
  }
  
  newdata_df <- as.data.frame(newdata)
  
  mf <- tryCatch(
    stats::model.frame(
      object$formula_strata_terms,
      data = newdata_df,
      na.action = stats::na.pass,
      xlev = object$formula_strata_xlevels
    ),
    error = function(e) {
      stop(
        "! This `mfp2` object was fitted with formula-level Cox strata, ",
        "but the strata term could not be reconstructed from `newdata`.\n",
        "i Include the original strata variable(s) in `newdata`, or supply ",
        "`strata = ...` explicitly to `predict()`.\n",
        "i Original error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  if (NROW(mf) != NROW(newdata_df)) {
    stop(
      "! Formula-level Cox strata reconstruction returned the wrong number of rows.",
      call. = FALSE
    )
  }
  
  if (anyNA(mf)) {
    stop(
      "! Reconstructed Cox strata contain missing values.\n",
      "i Please remove missing strata values from `newdata` before prediction.",
      call. = FALSE
    )
  }
  
  if (NCOL(mf) == 1L) {
    mf[[1L]]
  } else {
    do.call(
      survival::strata,
      c(as.list(mf), list(shortlabel = TRUE))
    )
  }
}

#' Retrieve the Cox Offset Reference Stored at Fit Time
#'
#' `survival::predict.coxph()` subtracts the mean training offset from each
#' prediction offset independently of the covariate reference. The final `fit_mfp()` result
#' stores that scalar as `cox_offset_reference`; prediction validates the field
#' before delegating raw offsets to `predict.coxph()`.
#'
#' @param object A fitted `mfp2` Cox object.
#'
#' @return Finite numeric scalar.
#'
#' @keywords internal
#' @noRd
mfp2_cox_offset_reference <- function(object) {
  reference <- object$cox_offset_reference
  
  if (!is.numeric(reference) || length(reference) != 1L ||
      is.na(reference) || !is.finite(reference)) {
    stop(
      paste0(
        "The fitted Cox model lacks valid offset-reference metadata. ",
        "Refit the mfp2 object with the current package version."
      ),
      call. = FALSE
    )
  }
  
  unname(reference)
}


#' Rebuild Formula-Level Offset from Prediction Newdata
#'
#' Internal helper used by predict.mfp2(). If an mfp2 object was fitted with
#' offset() in the formula interface, this reconstructs the offset vector from
#' ordinary newdata.
#'
#' @keywords internal
#' @noRd
reconstruct_formula_offset_newdata <- function(object, newdata) {
  if (!isTRUE(object$formula_interface) ||
      is.null(object$formula_offset_terms)) {
    return(NULL)
  }
  
  newdata_df <- as.data.frame(newdata)
  
  mf <- tryCatch(
    stats::model.frame(
      object$formula_offset_terms,
      data = newdata_df,
      na.action = stats::na.pass,
      xlev = object$formula_offset_xlevels
    ),
    error = function(e) {
      stop(
        "! This `mfp2` object was fitted with a formula-level offset, ",
        "but the offset could not be reconstructed from `newdata`.\n",
        "i Include the original offset variable(s) in `newdata`, or supply ",
        "`newoffset = ...` explicitly to `predict()`.\n",
        "i Original error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  out <- stats::model.offset(mf)
  
  if (is.null(out)) {
    stop(
      "! Formula-level offset could not be evaluated from `newdata`.",
      call. = FALSE
    )
  }
  
  out <- as.vector(out)
  
  if (length(out) != NROW(newdata_df)) {
    stop(
      "! Formula-level offset reconstruction returned the wrong number of rows.",
      call. = FALSE
    )
  }
  
  if (anyNA(out) || any(!is.finite(out))) {
    stop(
      "! Reconstructed offset contains missing or non-finite values.",
      call. = FALSE
    )
  }
  
  out
}

#' Helper function to prepare newdata for predict function
#' 
#' To be used in \code{predict.mfp2()}.
#' 
#' @param object fitted `mfp2` model object.
#' @param newdata Numeric design data on the raw model-matrix scale. Full-model
#'   prediction requires every fitted term. Internal term prediction may supply
#'   a subset consisting of complete conceptual terms.
#' @param strata,offset passed from \code{predict.mfp2()}. For Cox
#'   prediction, \code{strata} must be kept as a high-level vector/factor or
#'   combined multi-column strata object; integer conversion is not performed
#'   here because the stored Cox formula evaluates \code{strata(strata_)}.
#' @param apply_pre logical indicating whether the fitted pre-transformation
#' is applied or not.
#' @param apply_center logical indicating whether the fitted centers are applied
#' after transformation or not.
#' @param check_binary passed to \code{transform_vector_fp()}.
#' @param reset_zero Logical. If `TRUE`, variables marked as `zero = TRUE`
#' but containing only positive values are reset to `FALSE` before transformation.
#' The prediction helper defaults to `FALSE` because prediction must preserve
#' the zero/catzero/spike structure learned during model fitting. Adaptive
#' resetting based on the contents of `newdata` can create a design matrix that
#' is inconsistent with the fitted coefficients. Passed to
#' \code{transform_matrix()}.
#' @param allow_missing_predictors Logical. If `FALSE`, the default, `newdata`
#'   must contain every raw input column required by the fitted term lookup.
#'   For a grouped term, all member columns are required together. If `TRUE`,
#'   entire terms may be omitted and only complete terms present in `newdata`
#'   are transformed. This is intended for internal term-specific prediction,
#'   not for ordinary full-model prediction.
#' @return A data frame containing the transformed prediction design and any
#'   required Cox strata or offset columns.
#' @keywords internal
#' @noRd
prepare_newdata_for_predict <- function(object,
                                        newdata,
                                        strata = NULL,
                                        offset = NULL,
                                        apply_pre = TRUE,
                                        apply_center = TRUE,
                                        check_binary = TRUE,
                                        reset_zero = FALSE,
                                        allow_missing_predictors = FALSE) {
  # Keep formula-level strata and offset variables out of the numeric predictor
  # matrix. A data frame containing a factor stratum would otherwise be coerced
  # wholesale to a character matrix before the fitted predictor columns are
  # selected. Predictor isolation therefore precedes matrix conversion.
  newdata <- as.data.frame(newdata, check.names = FALSE)
  n_newdata <- nrow(newdata)
  
  if (is.null(colnames(newdata)) || any(colnames(newdata) == "")) {
    stop("! newdata must have non-empty column names.", call. = FALSE)
  }
  
  # Resolve the conceptual-term lookup once. Full prediction requires every raw
  # member column; term-specific prediction may include complete term subsets.
  lookup <- prediction_term_to_columns(object)
  expected <- unlist(lookup, use.names = FALSE)
  missing <- setdiff(expected, colnames(newdata))
  
  if (length(missing) > 0L && !isTRUE(allow_missing_predictors)) {
    stop(
      "! Missing required predictor column(s) in newdata: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  
  # Count supplied raw columns per conceptual term. Partial categorical blocks
  # are invalid because they would change the fitted parameterisation.
  present_by_term <- vapply(
    lookup,
    function(columns) sum(columns %in% colnames(newdata)),
    integer(1L)
  )
  partial_terms <- names(lookup)[
    present_by_term > 0L & present_by_term < lengths(lookup)
  ]
  
  if (length(partial_terms) > 0L) {
    details <- vapply(
      partial_terms,
      function(term) {
        absent <- setdiff(lookup[[term]], colnames(newdata))
        paste0(term, " (missing: ", paste(absent, collapse = ", "), ")")
      },
      character(1L)
    )
    stop(
      "! Grouped terms must be supplied with all member columns. Incomplete term(s): ",
      paste(details, collapse = "; "),
      call. = FALSE
    )
  }
  
  active_terms <- names(lookup)[present_by_term == lengths(lookup)]
  
  # Expand complete conceptual terms in lookup order. This keeps every grouped
  # block contiguous and gives transform_matrix() raw-column metadata in the
  # same order as the supplied design block.
  raw_columns <- unlist(lookup[active_terms], use.names = FALSE)
  
  if (length(raw_columns) == 0L) {
    stop("! No usable predictors were found in newdata.", call. = FALSE)
  }
  
  newdata <- newdata[, raw_columns, drop = FALSE]
  
  non_numeric <- raw_columns[
    !vapply(newdata, is.numeric, logical(1L))
  ]
  if (length(non_numeric) > 0L) {
    stop(
      "! Reconstructed predictor column(s) must be numeric: ",
      paste(non_numeric, collapse = ", "),
      call. = FALSE
    )
  }
  
  newdata <- as.matrix(newdata)
  storage.mode(newdata) <- "double"
  metadata <- expand_prediction_metadata(object, raw_columns)
  
  # Apply the fitted shift only to raw user data. Stored x_original data and
  # internally generated term grids are already on the shifted prediction scale.
  if (apply_pre) {
    newdata <- sweep(newdata, 2L, metadata$shift[raw_columns], "+")
    
    positive_terms <- active_terms[
      vapply(
        active_terms,
        function(term) {
          length(lookup[[term]]) == 1L && requires_positive_raw_input(object, term)
        },
        logical(1L)
      )
    ]
    
    if (length(positive_terms) > 0L) {
      bad_terms <- positive_terms[
        vapply(
          positive_terms,
          function(term) {
            column <- lookup[[term]][[1L]]
            any(!is.na(newdata[, column]) & newdata[, column] <= 0)
          },
          logical(1L)
        )
      ]
      
      if (length(bad_terms) > 0L) {
        bad_summary <- vapply(
          bad_terms,
          function(term) {
            column <- lookup[[term]][[1L]]
            count <- sum(!is.na(newdata[, column]) & newdata[, column] <= 0)
            paste0(term, " (", count, if (count == 1L) " row)" else " rows)")
          },
          character(1L)
        )
        
        stop(
          "After applying the shift values learned during fitting, some values in ",
          "`newdata` remain non-positive for terms whose fitted transformation ",
          "requires strictly positive input.\n",
          "i Problematic term(s): ", paste(bad_summary, collapse = ", "), ".",
          call. = FALSE
        )
      }
    }
  }
  
  # ACD prediction needs the parameters estimated at fit time for each active
  # singleton ACD term.
  active_acd <- raw_columns[
    metadata$acdx &
      !vapply(metadata$powers, function(p) is.null(p) || all(is.na(p)), logical(1L))
  ]
  missing_acd <- active_acd[
    vapply(
      active_acd,
      function(column) is.null(metadata$acd_parameter[[column]]),
      logical(1L)
    )
  ]
  
  if (length(missing_acd) > 0L) {
    stop(
      "Missing stored ACD parameters for prediction column(s): ",
      paste(missing_acd, collapse = ", "),
      call. = FALSE
    )
  }
  
  # Recreate the exact uncentred transformed columns used by the final model.
  x_trans <- transform_matrix(
    newdata,
    power_list = metadata$powers,
    center = stats::setNames(rep(FALSE, length(raw_columns)), raw_columns),
    keep_x_order = TRUE,
    acdx = metadata$acdx,
    acd_parameter_list = metadata$acd_parameter,
    check_binary = check_binary,
    zero = metadata$zero,
    catzero = metadata$catzero,
    spike = metadata$spike,
    spike_decision = metadata$spike_decision,
    reset_zero = reset_zero
  )
  
  if (is.null(x_trans)) {
    newdata <- matrix(numeric(0L), nrow = n_newdata, ncol = 0L)
  } else {
    newdata <- x_trans$x_transformed
    
    # Apply the stored final-model centres after transformation.
    if (apply_center && !is.null(object$centers)) {
      newdata <- center_matrix(
        newdata,
        centers = object$centers[colnames(newdata)],
        zero = x_trans$zero_expanded
      )
    }
  }
  
  if (NCOL(newdata) == 0L) {
    newdata <- data.frame(row.names = seq_len(n_newdata))
  } else {
    newdata <- data.frame(newdata, check.names = FALSE)
  }
  
  # Attach high-level strata and offset variables expected by the stored model
  # formula after the predictor design has been reconstructed.
  if (object$family_string == "cox" && !is.null(strata)) {
    strata_n <- if (is.vector(strata) || is.factor(strata)) length(strata) else NROW(strata)
    if (strata_n != nrow(newdata)) {
      stop("! `strata` must have one value or row per prediction row.", call. = FALSE)
    }
    if (anyNA(strata)) {
      stop("! `strata` must not contain missing values.", call. = FALSE)
    }
    
    newdata$strata_ <- if (is.matrix(strata) || is.data.frame(strata)) {
      do.call(
        survival::strata,
        c(as.list(as.data.frame(strata)), list(shortlabel = TRUE))
      )
    } else {
      strata
    }
  }
  
  if (!is.null(offset)) {
    if (!is.numeric(offset)) {
      stop("! offset must be numeric.", call. = FALSE)
    }
    if (length(offset) != nrow(newdata)) {
      stop("! offset must have one value per observation.", call. = FALSE)
    }
    newdata$offset_ <- offset
  }
  
  newdata
}
#' Helper function to compute standard error of a partial predictor
#' 
#' To be used in \code{predict.mfp2()}.
#' 
#' @param model fitted `mfp2` object.
#' @param X transformed input matrix with variables of interest for partial predictor.
#' @param xref transformed reference value for variable of interest. Default is
#'  `NULL`, in which case this function computes standard errors without reference 
#' values.
#' @param include_intercept logical indicating whether intercept variance and
#' covariance should be included when no reference value is supplied. This must
#' match the `add_intercept` choice used for term predictions.
#' 
#' @details 
#' See pages 91-92 and following in the book by Royston and Sauerbrei 2008
#' for the formulas and mathematical details.
#' 
#' @return 
#' Numeric vector of conditional standard errors.
#' 
#' @references
#' Royston, P. and Sauerbrei, W., 2008. \emph{Multivariable Model - Building: 
#' A Pragmatic Approach to Regression Analysis based on Fractional Polynomials 
#' for Modelling Continuous Variables. John Wiley & Sons.}\cr
#' @keywords internal
#' @noRd
calculate_standard_error <- function(model, 
                                     X, 
                                     xref = NULL,
                                     include_intercept = TRUE) { 
  
  vcovx <- vcov(object = model)
  
  # Small or rank-deficient fitted samples can yield undefined covariance entries.
  if (any(is.nan(vcovx)))
    warning("i NaN detected in the covariance matrix of the model.",
            "i Standard errors for calculation of confidence intervals may not exist")
  
  # Use the exact fit-time mapping between transformed matrix columns and
  # fitted coefficient names. This handles non-syntactic formula terms without
  # modifying either the transformed matrix names or the model covariance names.
  xnames <- colnames(X)
  model_columns <- prediction_model_column_names(model, xnames)
  ind <- match(model_columns, colnames(vcovx))
  
  if (anyNA(ind)) {
    stop(
      "Cannot compute SE; covariance matrix lacks mapped columns: ",
      paste(model_columns[is.na(ind)], collapse = ", "),
      call. = FALSE
    )
  }
  
  if (!is.null(xref)) {
    # Subtract the reference value: f(x)-f(xref)
    X <- sweep(X, 2, xref, "-")
  }
  
  # Augment X by the intercept only when the reported term value includes
  # the intercept. Contrasts and Cox partial predictors do not include intercept
  # variance.
  if (isTRUE(include_intercept) && model$family_string != "cox" && is.null(xref)) {
    intercept_ind <- match("(Intercept)", colnames(vcovx))
    if (is.na(intercept_ind)) {
      stop(
        "Cannot compute SE with intercept; covariance matrix lacks `(Intercept)`.",
        call. = FALSE
      )
    }
    X <- cbind("(Intercept)" = 1, X)
    ind <- c(intercept_ind, ind)
  }
  
  # Compute diag(X V X') without constructing the full n-by-n product.
  vcovx <- vcovx[ind, ind, drop = FALSE]
  v <- rowSums((X %*% vcovx) * X)
  
  sqrt(v)
}
