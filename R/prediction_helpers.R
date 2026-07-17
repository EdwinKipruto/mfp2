# Shared prediction helpers ---------------------------------------------------

#' Normalize Linear-Predictor Type Names Across Model Families
#'
#' GLM and Cox prediction methods use different names for the same
#' linear-predictor scale: `predict.glm()` uses `"link"`, whereas
#' `predict.coxph()` uses `"lp"`. This helper accepts either exact spelling and
#' returns the family-native name. Response-scale names are not translated,
#' because a GLM response and a Cox risk score are not generally equivalent.
#'
#' @param type Character scalar.
#' @param family_string Character scalar identifying the fitted family.
#'
#' @return Family-native prediction type.
#'
#' @keywords internal
#' @noRd
normalize_prediction_type <- function(type, family_string) {
  if (!is.character(type) || length(type) != 1L || is.na(type) || !nzchar(type)) {
    stop("`type` must be a single non-missing character value.", call. = FALSE)
  }
  
  # Prediction types occasionally arrive as named scalar character vectors
  # (for example when selected from a named lookup in user code or tests).
  # Names are irrelevant to dispatch but make identical(type, "link") false,
  # so normalize both inputs to unclassed, unnamed scalar strings first.
  type <- unname(as.character(type)[1L])
  family_string <- unname(as.character(family_string)[1L])
  
  if (identical(family_string, "cox") && identical(type, "link")) {
    return("lp")
  }
  
  if (!identical(family_string, "cox") && identical(type, "lp")) {
    return("link")
  }
  
  type
}


#' Match a Cox Covariate Reference
#'
#' Validates the covariate-centering reference understood by
#' `predict.coxph()` for linear-predictor and risk-score predictions.
#'
#' @param value Character scalar or `NULL`.
#' @param default Character scalar returned when `value` is `NULL`.
#' @param argument Character scalar used in error messages.
#'
#' @return One of `"zero"`, `"sample"`, or `"strata"`.
#'
#' @keywords internal
#' @noRd
match_cox_reference <- function(value,
                                default = "zero",
                                argument = "reference") {
  choices <- c("zero", "sample", "strata")
  
  if (is.null(value)) {
    value <- default
  }
  
  choice_text <- "'zero', 'sample', or 'strata'"
  
  if (!is.character(value) || length(value) != 1L ||
      is.na(value) || !nzchar(value)) {
    stop(
      "`", argument, "` must be one of ", choice_text, ".",
      call. = FALSE
    )
  }
  
  value <- unname(as.character(value)[1L])
  
  tryCatch(
    match.arg(value, choices),
    error = function(e) {
      stop(
        "`", argument, "` must be one of ", choice_text, ".",
        call. = FALSE
      )
    }
  )
}


#' Recover the Response Symbol Used by a Stored Cox Model
#'
#' @param fit_obj Fitted `coxph` object.
#'
#' @return Character scalar naming the response variable expected in
#'   prediction data.
#'
#' @keywords internal
#' @noRd
cox_internal_response_name <- function(fit_obj) {
  terms_object <- fit_obj$terms
  response_index <- attr(terms_object, "response")
  variables <- attr(terms_object, "variables")
  
  if (is.null(terms_object) ||
      is.null(response_index) ||
      length(response_index) != 1L ||
      response_index < 1L ||
      is.null(variables) ||
      length(variables) < response_index + 1L) {
    stop(
      "The fitted Cox model lacks response-term metadata needed for ",
      "absolute prediction. Refit the model with the current package version.",
      call. = FALSE
    )
  }
  
  response_expression <- variables[[response_index + 1L]]
  if (!is.symbol(response_expression)) {
    stop(
      "The stored Cox response is not a simple internal variable. ",
      "Refit the model with the current package version.",
      call. = FALSE
    )
  }
  
  as.character(response_expression)
}


#' Derive and Validate a Cox Prediction Response
#'
#' Absolute Cox predictions require a `Surv` response in the prediction model
#' frame. The response is taken from exactly one `Surv` column in `newdata`, or
#' reconstructed by evaluating the left-hand side of the original formula for
#' formula-interface fits.
#'
#' @param object Fitted package object containing formula-interface metadata.
#' @param fit_obj Stored fitted `coxph` object whose response type is required.
#' @param newdata Raw prediction data, before predictor reconstruction.
#'
#' @return Validated `Surv` object aligned with `newdata`.
#'
#' @keywords internal
#' @noRd
reconstruct_cox_prediction_response <- function(object, fit_obj, newdata) {
  newdata_df <- as.data.frame(newdata, check.names = FALSE)
  
  surv_columns <- names(newdata_df)[vapply(
    newdata_df,
    function(column) inherits(column, "Surv"),
    logical(1L)
  )]
  
  if (length(surv_columns) > 1L) {
    stop(
      "Cox predictions of type 'expected' or 'survival' found more than ",
      "one Surv column in `newdata`: ", paste(surv_columns, collapse = ", "),
      ". Supply exactly one prediction response.",
      call. = FALSE
    )
  }
  
  response <- if (length(surv_columns) == 1L) {
    newdata_df[[surv_columns]]
  } else {
    NULL
  }
  
  if (is.null(response) && isTRUE(object$formula_interface)) {
    if (!inherits(object$formula, "formula")) {
      stop(
        "The fitted formula-interface Cox model lacks its original formula ",
        "and cannot derive the prediction response from `newdata`. Refit the ",
        "model with the current package version.",
        call. = FALSE
      )
    }
    
    response_formula <- object$formula
    response_formula[[3L]] <- 1
    
    formula_environment <- environment(response_formula)
    if (is.null(formula_environment)) {
      formula_environment <- parent.frame()
    }
    if (!exists("Surv", envir = formula_environment, inherits = TRUE)) {
      response_environment <- new.env(parent = formula_environment)
      response_environment$Surv <- survival::Surv
      environment(response_formula) <- response_environment
    }
    
    response <- tryCatch(
      stats::model.response(
        stats::model.frame(
          response_formula,
          data = newdata_df,
          na.action = stats::na.pass
        )
      ),
      error = function(e) {
        stop(
          "Cox predictions of type 'expected' or 'survival' require the ",
          "follow-up response information in `newdata`. Include the variables ",
          "used on the left-hand side of the fitted formula, or include exactly one ",
          "Surv column in `newdata`.\nOriginal error: ", conditionMessage(e),
          call. = FALSE
        )
      }
    )
  }
  
  if (is.null(response)) {
    stop(
      "Cox predictions of type 'expected' or 'survival' require the follow-up ",
      "response information in `newdata`. For a formula fit, include the ",
      "original response variables. Otherwise include exactly one Surv column, ",
      "for example `data.frame(x = xnew, y = I(Surv(time, status)))`.",
      call. = FALSE
    )
  }
  
  if (!inherits(response, "Surv")) {
    stop(
      "The Cox prediction response derived from `newdata` is not a ",
      "survival::Surv object.",
      call. = FALSE
    )
  }
  
  if (NROW(response) != NROW(newdata_df)) {
    stop(
      "The Cox prediction response must have one row per row of `newdata`.",
      call. = FALSE
    )
  }
  
  response_matrix <- unclass(response)
  if (anyNA(response_matrix) || any(!is.finite(response_matrix))) {
    stop(
      "The Cox prediction response in `newdata` must contain only finite, ",
      "non-missing values.",
      call. = FALSE
    )
  }
  
  fitted_response <- fit_obj$y
  if (is.null(fitted_response)) {
    fitted_response <- object$y_original
  }
  
  if (!inherits(fitted_response, "Surv")) {
    stop(
      "The fitted Cox model lacks a valid stored Surv response. Refit the ",
      "model with the current package version.",
      call. = FALSE
    )
  }
  
  if (!identical(attr(response, "type"), attr(fitted_response, "type")) ||
      NCOL(response) != NCOL(fitted_response)) {
    stop(
      "The Cox prediction response in `newdata` has a different survival ",
      "type from the fitted model.",
      call. = FALSE
    )
  }
  
  response
}


#' Attach a Cox Prediction Response to Reconstructed New Data
#'
#' @param fit_obj Stored fitted `coxph` object.
#' @param newdata Reconstructed formula-compatible prediction data.
#' @param response Validated `Surv` response.
#'
#' @return `newdata` with the internal Cox response column appended.
#'
#' @keywords internal
#' @noRd
attach_cox_prediction_response <- function(fit_obj, newdata, response) {
  newdata <- as.data.frame(newdata, check.names = FALSE)
  response_name <- cox_internal_response_name(fit_obj)
  
  if (response_name %in% names(newdata)) {
    stop(
      "The reconstructed prediction data already contain the internal Cox ",
      "response column '", response_name, "'.",
      call. = FALSE
    )
  }
  
  if (!inherits(response, "Surv") || NROW(response) != nrow(newdata)) {
    stop("Internal error: invalid or misaligned Cox prediction response.",
         call. = FALSE)
  }
  
  newdata[[response_name]] <- I(response)
  newdata
}

#' Detect Stratification in a Stored Cox Model
#'
#' Centralizes the fitted-model metadata check used by both prediction methods.
#' A model is treated as stratified when the retained `coxph` object records
#' strata directly or its terms object contains a `strata()` special.
#'
#' @param fit_obj Fitted model object.
#'
#' @return Logical scalar.
#'
#' @keywords internal
#' @noRd
prediction_fit_has_strata <- function(fit_obj) {
  if (!inherits(fit_obj, "coxph")) {
    return(FALSE)
  }
  
  if (!is.null(fit_obj$strata) && length(fit_obj$strata) > 0L) {
    return(TRUE)
  }
  
  trm <- fit_obj$terms
  if (is.null(trm)) {
    trm <- try(stats::terms(fit_obj), silent = TRUE)
    if (inherits(trm, "try-error")) {
      return(FALSE)
    }
  }
  
  specials <- attr(trm, "specials")
  if (!is.null(specials$strata) && length(specials$strata) > 0L) {
    return(TRUE)
  }
  
  labels <- attr(trm, "term.labels")
  if (is.null(labels)) {
    return(FALSE)
  }
  any(grepl("(^|:)strata\\(", labels))
}

