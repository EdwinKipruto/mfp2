#' Predict from an `mfp2` Model
#'
#' Obtain predictions from a model fitted with [mfp2()]. The function can
#' return predictions from the complete model, term-specific fitted values,
#' or differences from chosen reference values.
#'
#' @details
#' Supply predictors using their original values. The fitted shifts, FP or ACD
#' functions, centering, zero handling, and spike-at-zero decisions are applied
#' automatically.
#'
#' Do not create fractional-polynomial basis columns, centered columns, or
#' structural-zero indicator columns yourself.
#'
#' Standard errors and confidence intervals are conditional on the final fitted
#' model. They do not include uncertainty from variable selection or from
#' choosing the fitted function.
#'
#' Term and contrast predictions require enough valid covariance information to
#' calculate their standard errors and confidence intervals. Negligible negative
#' variances caused by numerical rounding are treated as zero. The same checks
#' apply to full-model GEE, fastglm, ordinal, and multinomial standard errors.
#' Prediction stops
#' if a variance is non-finite or meaningfully negative and a valid interval
#' cannot be calculated.
#'
#' @section Supplying prediction data:
#' Supply `newdata` with one row for each observation to predict.
#'
#' If you fitted the model using the formula interface, supply the variables
#' named in the formula. For example, for a model fitted with
#' `y ~ age + race`, prediction data could be:
#'
#' ```r
#' newdata <- data.frame(
#'   age = c(50, 60, 70),
#'   race = c("A", "B", "C")
#' )
#' ```
#'
#' Supply `race` itself; you do not need to create separate columns for its
#' categories. The race categories `"A"`, `"B"`, and `"C"` must have been
#' available during fitting.
#'
#' If you fitted the model using the matrix interface, supply the same
#' column names and use the same coding as in the fitting matrix.
#' For example, if race had three categories during fitting, with race
#' category A as the reference category, the fitting matrix could contain
#' `age`, `raceB`, and `raceC`. The same prediction data would then be:
#'
#' ```r
#' newdata <- cbind(
#'   age   = c(50, 60, 70),
#'   raceB = c(0, 1, 0),
#'   raceC = c(0, 0, 1)
#' )
#' ```
#'
#' The first row represents an observation aged 50 in race category A;
#' the second, an observation aged 60 in race category B; and the third,
#' an observation aged 70 in race category C.
#' Race category A has zero in both `raceB` and `raceC`.
#'
#' If you grouped these columns during fitting with
#' `term_groups = list(race = c("raceB", "raceC"))`, include both `raceB`
#' and `raceC` when predicting the group. In `terms`, use the group name
#' `"race"`, rather than the column names. For example:
#'
#' ```r
#' predict(fit, newdata = newdata, type = "terms", terms = "race")
#' ```
#'
#' For predictions from the whole model, include all variables kept in the
#' final model. For term or contrast predictions, include only the variables
#' needed for the requested terms.
#'
#' Extra columns may be needed for offsets or survival groups (`strata`)
#' specified in the formula. For Cox `type = "expected"` or `"survival"`,
#' also include the time and event-status variables in `newdata`, using the
#' names from the fitting formula. The time values specify when predictions
#' are evaluated.
#'
#' Values needed for prediction must not be missing or contain `Inf` or
#' `-Inf`. Extra columns that are not used are ignored.
#'
#' Supply values on their original scale; fitted shifts and transformations
#' are applied automatically. Prediction stops if a transformation cannot
#' use a supplied value. For example, a logarithm requires a positive value
#' after the fitted shift is added.
#'
#' @section Choosing the prediction type:
#' The default is `"lp"` for Cox and Fine--Gray models and `"link"` for
#' all other models. `"link"` and `"lp"` are interchangeable.
#'
#' For GLM (including negative-binomial) and GEE models:
#'
#' - `type = "link"` returns the linear predictor.
#' - `type = "response"` returns the predicted mean on the response scale.
#'
#' For GEE models, standard errors use the covariance estimator chosen
#' during fitting. Response-scale standard errors are approximated from
#' those on the link scale.
#'
#' For parametric `survreg` models:
#'
#' - `type = "link"` returns the fitted location (linear predictor).
#'   `"linear"` is also accepted.
#' - `type = "response"` returns the fitted location transformed back to
#'   the response scale. This is not generally the mean survival time.
#' - `type = "quantile"` returns predicted quantiles on the response scale.
#' - `type = "uquantile"` returns predicted quantiles on the linear-predictor
#'   scale.
#'
#' For `"quantile"` and `"uquantile"`, supply probabilities using `p`,
#' for example, `p = 0.5` for the median. See [survival::predict.survreg()].
#'
#' For Cox models:
#'
#' - `type = "lp"` returns the relative log-hazard.
#' - `type = "risk"` returns the relative hazard, calculated as `exp(lp)`.
#' - `type = "expected"` returns the predicted cumulative hazard up to
#'   the supplied follow-up time.
#' - `type = "survival"` returns the survival probability at that time,
#'   calculated as `exp(-expected)`.
#'
#' For Fine--Gray models:
#'
#' - `type = "lp"` returns the relative log subdistribution hazard.
#' - `type = "risk"` returns the relative subdistribution hazard,
#'   calculated as `exp(lp)`.
#' - `type = "response"` returns cumulative-incidence probabilities at
#'   the times supplied in `times`.
#'
#' For ordinal models:
#'
#' - `type = "link"` returns the linear predictor including the cutoff
#'   intercept selected by `rms::orm()`. Set `add_intercept = FALSE` to
#'   return only the predictor effects and any offset.
#' - `type = "response"` returns category probabilities, with one row per
#'   observation and one column per category.
#' - `type = "mean"` returns the expected category value. Numeric category
#'   labels are used when possible; otherwise, categories are scored
#'   `1, ..., K`, where `K` is the number of categories.
#'
#' For multinomial models:
#'
#' - `type = "link"` returns log-odds relative to the reference category,
#'   with one row per observation and one column per non-reference category.
#' - `type = "response"` returns category probabilities, with one row per
#'   observation and one column per category.
#' - `type = "class"` returns the category with the highest predicted
#'   probability.
#'
#' Use `se.fit = TRUE` to also obtain standard errors for complete-model
#' numeric predictions. Ordinal probability and mean standard errors and
#' multinomial probability standard errors use a first-order approximation
#' (the delta method). Multinomial class labels have no standard error;
#' use `type = "response"` to obtain uncertainty for their probabilities.
#'
#' For all supported models:
#'
#' - `type = "terms"` returns each requested term's contribution on the
#'   linear-predictor scale.
#' - `type = "contrasts"` returns the difference between each term's
#'   contribution and its contribution at a reference value.
#'
#' Term predictions include the intercept by default for GLM, GEE,
#' `survreg`, and multinomial models. Set `add_intercept = FALSE` to omit it.
#' Cox and Fine--Gray terms have no intercept. Ordinal terms omit cutoff
#' intercepts unless `ordinal_cutoff` is supplied. Intercepts cancel from
#' contrast predictions.
#'
#' Term and contrast predictions include standard errors and confidence
#' intervals. For multinomial models, each term's result contains separate
#' predictions for each non-reference category, identified by the `logit`
#' column.
#'
#' @section Term predictions:
#' A single term result can represent an FP2 function, a factor, a zero-handled
#' variable, or several matrix columns grouped with `term_groups`. Contributions
#' from all parts of the term are combined in that result.
#'
#' Use the original variable or group name in `terms`. Use
#' [get_selected_variable_names()] to see the terms included in the final model.
#' A requested term that is not in the final model is omitted with a warning.
#'
#' For a single numeric term:
#'
#' - `terms_seq = "equidistant"` uses `nseq` equally spaced values between the
#'   smallest and largest available values.
#' - `terms_seq = "data"` uses the values in `newdata`, or the fitted values
#'   when `newdata = NULL`.
#'
#' Factors and other grouped terms are evaluated row by row.
#'
#' For GLM (including negative-binomial), GEE, and `survreg` models,
#' `add_intercept = TRUE` includes the fitted intercept in each term result;
#' for multinomial models it includes each logit's own intercept. Set `add_intercept = FALSE` to return only the
#' contribution from the requested term. This argument does not affect Cox,
#' Fine--Gray, or contrast predictions.
#'
#' Ordinal term predictions without `ordinal_cutoff` return the common term
#' contribution and calculate its standard error from the covariance of that
#' term's coefficients. Selecting a cutoff adds that cutoff's fitted intercept
#' and includes its variance and covariance with the term coefficients in the
#' standard error. Ordinal contrasts cancel the intercept at every cutoff.
#'
#' @section Contrasts and reference values:
#' Use `type = "contrasts"` to compare a term's contribution at different
#' values with its contribution at a chosen reference value.
#' Each contrast is calculated as:
#'
#' `term contribution at the prediction value - term contribution at the reference`.
#'
#' The contrast is zero at the reference value. A positive contrast means
#' a larger contribution than at the reference; a negative contrast means
#' a smaller contribution. Differences are reported on the model's
#' linear-predictor scale, such as log-odds for logistic regression.
#'
#' **Example using the formula interface**
#'
#' For a model that includes `age` and the factor `race`, compare ages
#' 50, 60, and 70 with age 50, and race categories A, B, and C with
#' race category A:
#'
#' ```r
#' newdata <- data.frame(
#'   age = c(50, 60, 70),
#'   race = c("A", "B", "C")
#' )
#'
#' predict(
#'   fit_formula,
#'   newdata = newdata,
#'   type = "contrasts",
#'   terms = c("age", "race"),
#'   terms_seq = "data",
#'   ref = list(age = 50, race = "A")
#' )
#' ```
#'
#' The result contains separate age and race comparisons. The age contrasts
#' compare each age with 50; the race contrasts compare each category with A.
#' The contrast for age 50 is zero, as is the contrast for race category A.
#' For example, if the fitted age contributions at ages 60 and 50 are
#' 0.8 and 0.3, the contrast at age 60 is `0.8 - 0.3 = 0.5`.
#'
#' Supply references as a named list: each name identifies a term, and its
#' value specifies the reference. For formula factors, use a category label
#' that was available during fitting, such as `"A"`.
#'
#' **Numeric references**
#'
#' Supply numeric references on the original variable scale. Any fitted
#' shift, FP or ACD transformation, and centering are applied automatically.
#' For example, if the fitted function uses `log(age + 10)`, supply
#' `ref = list(age = 50)`, not 60 or `log(60)`.
#'
#' **References for grouped matrix columns**
#'
#' If race was represented by `raceB` and `raceC`, with race category A
#' as the reference category, and grouped during fitting using
#' `term_groups = list(race = c("raceB", "raceC"))`, compare with
#' race category A by supplying both column values as zero:
#'
#' ```r
#' newdata <- cbind(
#'   raceB = c(0, 1, 0),
#'   raceC = c(0, 0, 1)
#' )
#'
#' predict(
#'   fit_matrix,
#'   newdata = newdata,
#'   type = "contrasts",
#'   terms = "race",
#'   ref = list(race = c(raceB = 0, raceC = 0))
#' )
#' ```
#'
#' The rows compare race categories A, B, and C with race category A.
#' For comparison with race category B instead, use
#' `ref = list(race = c(raceB = 1, raceC = 0))`.
#'
#' For other grouped terms, including formula interactions, supply one
#' numeric reference value for every member column. Use the column names
#' supplied in `term_groups`, or those shown in the term prediction result.
#'
#' **Default references**
#'
#' If a reference is not supplied for a term, the function uses:
#'
#' - the smaller value for a numeric term with two distinct values;
#' - the mean for another single numeric term;
#' - the first fitted category for a formula factor;
#' - the first available row for another grouped term.
#'
#' Numeric and grouped-row defaults use `newdata` when supplied, otherwise
#' the fitting data. Formula factors always use their first fitted category.
#' Supply `ref` explicitly to keep the same comparison when `newdata` changes.
#'
#' The model intercept cancels from the comparison and is not included
#' in contrasts.
#'
#' @section Zero and spike-at-zero variables:
#' Supply zero-handled and spike-at-zero variables on their original scale.
#' Do not create a separate zero-indicator column.
#'
#' The zero handling selected during fitting is reused during prediction. These
#' options require nonnegative covariates: `x = 0` belongs to the zero component
#' and `x > 0` belongs to the positive component. Negative values in `newdata`
#' are rejected; recode them explicitly before prediction if they should
#' scientifically represent the zero group.
#'
#' If the final spike-at-zero model includes both a continuous component and a
#' zero indicator, both are used automatically. If it includes only one
#' component, only that component is used.
#'
#' The same rule is applied to a supplied reference. For a term that includes
#' only the zero indicator, an exact-zero reference is treated as indicator 1
#' and a positive reference as indicator 0. Negative references are rejected.
#'
#' @section Offsets:
#' Offsets are included in predictions from the whole model. They are not
#' included in term or contrast predictions.
#'
#' When `newdata = NULL`, the offsets used during fitting are reused.
#' For new observations, provide offsets as described below.
#'
#' **Offset in the formula**
#'
#' If the fitting formula contains `offset(...)`, include the variables
#' needed to calculate it in `newdata`. The offset is calculated automatically.
#' For example, for a model fitted with
#' `count ~ age + offset(log(exposure))`:
#'
#' ```r
#' newdata <- data.frame(
#'   age = c(50, 60),
#'   exposure = c(1, 2)
#' )
#'
#' predict(fit_formula, newdata = newdata, type = "response")
#' ```
#'
#' Here, `log(exposure)` is calculated separately for each new observation.
#'
#' **Offset supplied separately**
#'
#' If the offset was supplied through the `offset` argument during fitting,
#' supply the new offset values through `newoffset`. For example, if the
#' fitting offset was `log(exposure)`:
#'
#' ```r
#' predict(
#'   fit_offset,
#'   newdata = newdata["age"],
#'   newoffset = log(newdata$exposure),
#'   type = "response"
#' )
#' ```
#'
#' Supply the calculated offset, such as `log(exposure)`, rather than
#' `exposure` itself. An explicit `newoffset` overrides a formula offset;
#' the variables used only by that formula offset can then be omitted.
#'
#' For most models, supply one finite `newoffset` value per prediction row.
#' For example, `newoffset = c(0.1, 0.2)` supplies offsets for two observations.
#' `newoffset` requires `newdata` and a model fitted with an offset.
#'
#' **Multinomial offsets**
#'
#' For multinomial models, offsets must form a matrix with one row per
#' observation. Columns can represent all response categories or only
#' the non-reference categories. Values are offsets on the logit scale.
#'
#' A formula offset must produce such a matrix. For example, for response
#' categories A, B, and C, with A as the reference category:
#'
#' ```r
#' fit_multinomial <- mfp2(
#'   outcome ~ age + offset(cbind(A = off_A, B = off_B, C = off_C)),
#'   data = training,
#'   family = multinomial_family(reference = "A")
#' )
#'
#' newdata <- data.frame(
#'   age   = c(50, 60),
#'   off_A = c(0, 0),
#'   off_B = c(0.2, 0.4),
#'   off_C = c(-0.1, 0.3)
#' )
#'
#' predict(fit_multinomial, newdata = newdata, type = "response")
#' ```
#'
#' The offset matrix is calculated automatically from the columns in
#' `newdata`. The offset for reference category A is subtracted from
#' those for categories B and C.
#'
#' Alternatively, supply the matrix directly through `newoffset`.
#' This overrides the formula offset, so only `age` is needed in `newdata`:
#'
#' ```r
#' newoffset <- cbind(
#'   A = newdata$off_A,
#'   B = newdata$off_B,
#'   C = newdata$off_C
#' )
#'
#' predict(
#'   fit_multinomial,
#'   newdata = newdata["age"],
#'   newoffset = newoffset,
#'   type = "response"
#' )
#' ```
#'
#' Because the offsets for reference category A are zero here, supplying
#' only columns B and C gives the same predictions:
#'
#' ```r
#' predict(
#'   fit_multinomial,
#'   newdata = newdata["age"],
#'   newoffset = newoffset[, c("B", "C"), drop = FALSE],
#'   type = "response"
#' )
#' ```
#'
#' With columns for only the non-reference categories, supply offsets
#' relative to the reference category. Name the columns using the
#' corresponding response-category labels. The same matrix formats apply
#' when the offset was supplied separately during fitting.
#'
#' @section Proportional-hazards predictions:
#' Cox and Fine--Gray `type = "lp"` predictions are relative log-hazards;
#' `type = "risk"` returns their exponentials. For Fine--Gray models,
#' these refer to the subdistribution hazard. Neither type returns an
#' absolute hazard or a survival probability.
#'
#' **Reference for relative predictions**
#'
#' Use `cox_reference` with `"lp"` or `"risk"`:
#'
#' - `"zero"` uses zero on the fitted predictor scale and is the default.
#' - `"sample"` uses the average fitted predictor values in the fitting sample.
#' - `"strata"` uses separate averages within each fitted stratum.
#'
#' For example:
#'
#' ```r
#' predict(
#'   fit_cox,
#'   newdata = data.frame(age = c(50, 60)),
#'   type = "risk",
#'   cox_reference = "sample"
#' )
#' ```
#'
#' Changing the reference changes the reported log-hazards and risk scores.
#' Ratios of risk scores remain unchanged when a common reference is used.
#' With `"strata"`, make comparisons within the same stratum.
#' `cox_reference` is used only with `"lp"` and `"risk"`.
#'
#' **Cox survival probabilities and cumulative hazards**
#'
#' Use `type = "survival"` for survival probabilities and
#' `type = "expected"` for cumulative hazards.
#' When supplying `newdata`, also include the variables used in the fitted
#' `Surv()` response. The time values specify when predictions are evaluated.
#'
#' For example, suppose the model was fitted with
#' `Surv(time, event) ~ age`, with time measured in days and event coded 0/1.
#' To predict survival at 365 days:
#'
#' ```r
#' newdata <- data.frame(
#'   age = c(50, 60),
#'   time = c(365, 365),
#'   event = c(0, 0)
#' )
#'
#' predict(fit_cox, newdata = newdata, type = "survival")
#' ```
#'
#' The event-status column is required to construct the `Surv()` response,
#' but its values do not affect these predictions. Use valid status values.
#' With `type = "expected"`, the same data give cumulative hazards at 365 days.
#'
#' Alternatively, include exactly one `Surv` column:
#'
#' ```r
#' newdata <- data.frame(
#'   age = c(50, 60),
#'   followup = I(survival::Surv(c(365, 365), c(0, 0)))
#' )
#' ```
#'
#' For matrix-interface fits, supply follow-up information this way.
#' When `newdata = NULL`, the follow-up times used during fitting are reused.
#'
#' **Stratified models**
#'
#' Each observation must belong to a stratum available during fitting.
#' If the formula contains `strata(...)`, include its variables in `newdata`.
#' For example, for a model fitted with
#' `Surv(time, event) ~ age + strata(sex)`:
#'
#' ```r
#' newdata <- data.frame(
#'   age = c(50, 60),
#'   sex = c("F", "M")
#' )
#'
#' predict(fit_stratified, newdata = newdata, type = "risk")
#' ```
#'
#' Here, `"F"` and `"M"` must match fitted stratum labels.
#' Do not also supply the `strata` argument for formula-level strata.
#'
#' If strata were supplied separately during matrix-interface fitting,
#' supply them separately for the new observations:
#'
#' ```r
#' predict(
#'   fit_matrix,
#'   newdata = as.matrix(newdata["age"]),
#'   strata = newdata$sex,
#'   type = "risk"
#' )
#' ```
#'
#' Strata are matched by their labels, not by numeric positions.
#' For several stratification variables, supply the same variables in
#' the same column order used during fitting. New strata are not allowed.
#'
#' **Fine--Gray cumulative incidence**
#'
#' Use `type = "response"` with `times` to predict cumulative incidence.
#' Supply times in the same units used during fitting. For example:
#'
#' ```r
#' predict(
#'   fit_finegray,
#'   newdata = data.frame(age = c(50, 60)),
#'   type = "response",
#'   times = c(365, 730)
#' )
#' ```
#'
#' @param object A fitted object of class `"mfp2"`.
#'
#' @param newdata An optional data frame or matrix containing observations to
#'   predict. See **Supplying prediction data**.
#'
#' @param type A character value selecting the prediction type. See
#'   **Choosing the prediction type**.
#'
#' @param se.fit Logical; return standard errors for complete-model predictions.
#'   Not available for multinomial class labels. Ignored for term and contrast
#'   predictions, which always include standard errors and confidence limits.
#'   See **Choosing the prediction type**.
#'
#' @param terms An optional character vector naming terms to return with
#'   `type = "terms"` or `type = "contrasts"`. Use original variable or group
#'   names. If `NULL`, all terms included in the final model are used.
#'
#' @param terms_seq For term and contrast predictions, choose which values
#'   to predict at. `"equidistant"` (the default) uses `nseq` evenly spaced
#'   values between the variable's minimum and maximum. `"data"` uses the
#'   observed values. Both use `newdata` when supplied, otherwise the fitting
#'   data. Factors and grouped terms always use the supplied or fitted rows.
#'
#' @param alpha A single finite number strictly between 0 and 1 used to
#'   calculate confidence intervals for term and contrast results. The confidence level is `1 - alpha`. The
#'   default `alpha = 0.05` gives 95 percent confidence intervals.
#'
#' @param ref A named list of reference values used with
#'   `type = "contrasts"`. Numeric references must be supplied on the original
#'   variable scale. See **Contrasts and reference values**.
#'
#' @param strata Strata for new rows of a Cox or `survreg` model fitted with
#'   separately supplied strata. Supply a vector for one stratification
#'   variable, or a matrix or data frame for several. For formula strata,
#'   include the variables in `newdata` instead. See **Proportional-hazards
#'   predictions** and **Parametric survival predictions**.
#'
#' @param newoffset Offsets for new rows: one finite numeric value per row,
#'   or, for multinomial models, a matrix with one column per category or
#'   non-reference category. See **Offsets**.
#'
#' @param times Time points for Fine--Gray cumulative-incidence predictions,
#'   in the same units as the fitted follow-up time. Required for
#'   `type = "response"` in Fine--Gray models; unavailable otherwise.
#'   See **Proportional-hazards predictions**.
#'
#' @param nseq A single finite positive integer giving the number of equally
#'   spaced values used when `terms_seq = "equidistant"`. The default is 100.
#'
#' @param add_intercept Logical; defaults to `TRUE`. Includes the fitted
#'   intercept in term predictions where applicable and the reference-cutoff
#'   intercept in full-model ordinal link predictions. See **Choosing the
#'   prediction type**.
#' @param ordinal_cutoff Cutoff number for ordinal term predictions. Includes
#'   the chosen cutoff's intercept and accounts for its uncertainty.
#'   The default, `NULL`, omits cutoff intercepts. See **Choosing the
#'   prediction type**.
#'
#' @param cox_reference Reference for Cox or Fine--Gray `"lp"` and `"risk"`
#'   predictions: `"zero"` (default), `"sample"`, or `"strata"`.
#'   See **Proportional-hazards predictions**.
#'
#' @param ... Additional arguments for the underlying prediction method,
#'   such as `p = 0.5` for `survreg` quantiles. Internal prediction methods
#'   ignore these arguments; negative-binomial predictions reject unknown
#'   arguments. See **Choosing the prediction type**.
#'
#' @return
#' Complete-model predictions usually return a numeric vector, or a list
#' containing predictions and standard errors when `se.fit = TRUE`.
#' Where a native prediction method is used, its return format is retained.
#'
#' GEE predictions with `se.fit = TRUE` return a list containing `fit`,
#' `se.fit`, and `residual.scale`.
#'
#' Multinomial link and response predictions return matrices; class predictions
#' return a factor. Ordinal link and mean predictions return numeric vectors;
#' response predictions return a matrix with one column per ordered category.
#' For these numeric predictions, `se.fit = TRUE` returns a list containing
#' `fit` and `se.fit` with matching dimensions and observation names.
#' See **Choosing the prediction type**.
#'
#' Term and contrast predictions return a named list with one data frame per
#' requested term, containing:
#'
#' - `variable`: the evaluated value or factor level;
#' - `variable_pre`: for a single numeric term, the shifted value before the
#'   final transformation; for an indicator-only spike term, the zero indicator;
#' - `value`: the fitted term value or contrast;
#' - `se`: its standard error;
#' - `lower`, `upper`: its confidence limits.
#'
#' Grouped terms may also include their member columns. Multinomial term and
#' contrast results include a `logit` column identifying the non-reference
#' category.
#'
#' @examples
#' ## Gaussian models: complete predictions, terms, and contrasts
#' data("prostate")
#'
#' fit_gaussian <- mfp2(
#'   lpsa ~ fp(age) + fp(cavol) + fp(weight) + svi,
#'   data = prostate,
#'   family = "gaussian",
#'   keep = c("age", "cavol"),
#'   verbose = FALSE
#' )
#'
#' # Predictions for the fitting data.
#' predict(fit_gaussian)
#'
#' # Predictions and standard errors for five rows.
#' predict(
#'   fit_gaussian,
#'   newdata = prostate[1:5, ],
#'   type = "response",
#'   se.fit = TRUE
#' )
#'
#' # For this Gaussian model, the response is lpsa.
#' # These predictions are therefore on the log-PSA scale.
#'
#' # Term contributions, without the model intercept.
#' predict(
#'   fit_gaussian,
#'   type = "terms",
#'   terms = c("age", "cavol"),
#'   add_intercept = FALSE
#' )
#'
#' # Evaluate exactly these ages, rather than an equally spaced grid.
#' # Only age is needed for an age term prediction.
#' predict(
#'   fit_gaussian,
#'   newdata = data.frame(age = c(40, 50, 60)),
#'   type = "contrasts",
#'   terms = "age",
#'   terms_seq = "data",
#'   ref = list(age = 50)
#' )
#'
#' # References use the original variable scale, even for fitted FP functions.
#' predict(
#'   fit_gaussian,
#'   type = "contrasts",
#'   terms = c("age", "cavol"),
#'   ref = list(age = 50, cavol = 1)
#' )
#'
#' \dontrun{
#' ## Gaussian models: a predictor with a spike at zero
#' set.seed(1)
#' d_spike <- data.frame(
#'   exposure = c(rep(0, 35), runif(85, 0.1, 8))
#' )
#' d_spike$y <- with(
#'   d_spike,
#'   2 * (exposure > 0) + log1p(exposure) + rnorm(120)
#' )
#'
#' fit_spike <- mfp2(
#'   y ~ fp(exposure, spike = TRUE),
#'   data = d_spike,
#'   family = "gaussian",
#'   keep = "exposure",
#'   verbose = FALSE
#' )
#'
#' # Compare exactly 0, 1, and 4 with exposure = 0.
#' predict(
#'   fit_spike,
#'   newdata = data.frame(exposure = c(0, 1, 4)),
#'   type = "contrasts",
#'   terms = "exposure",
#'   terms_seq = "data",
#'   ref = list(exposure = 0)
#' )
#'
#' ## Poisson models: formula offsets and separately supplied offsets
#' set.seed(2)
#' d_count <- data.frame(
#'   age = runif(80, 30, 75),
#'   exposure = runif(80, 0.5, 3)
#' )
#' d_count$count <- rpois(
#'   80,
#'   lambda = with(d_count, exposure * exp(0.02 * (age - 50)))
#' )
#'
#' fit_poisson <- mfp2(
#'   count ~ fp(age) + offset(log(exposure)),
#'   data = d_count,
#'   family = "poisson",
#'   keep = "age",
#'   verbose = FALSE
#' )
#'
#' # The formula calculates log(exposure) from the new rows.
#' nd_count <- d_count[1:4, c("age", "exposure")]
#' predict(fit_poisson, newdata = nd_count, type = "response")
#'
#' # Response predictions are expected counts for the supplied exposure.
#' # Set exposure to 1 to obtain expected counts per unit of exposure.
#' nd_unit <- nd_count
#' nd_unit$exposure <- 1
#' predict(fit_poisson, newdata = nd_unit, type = "response")
#'
#' fit_external_offset <- mfp2(
#'   count ~ fp(age),
#'   data = d_count,
#'   family = "poisson",
#'   offset = log(d_count$exposure),
#'   keep = "age",
#'   verbose = FALSE
#' )
#'
#' # A separately supplied offset must be supplied again for new rows.
#' predict(
#'   fit_external_offset,
#'   newdata = nd_count["age"],
#'   newoffset = log(nd_count$exposure),
#'   type = "response"
#' )
#'
#' ## GEE models: predictions and standard errors
#' # Illustrate groups of four observations.
#' d_count$id <- rep(seq_len(20), each = 4)
#' fit_gee <- mfp2(
#'   count ~ fp(age) + offset(log(exposure)),
#'   data = d_count,
#'   family = gee_family(poisson(), corstr = "independence"),
#'   id = id,
#'   keep = "age",
#'   verbose = FALSE
#' )
#'
#' # New group identifiers are not needed for these mean predictions.
#' predict(
#'   fit_gee,
#'   newdata = nd_count,
#'   type = "response",
#'   se.fit = TRUE
#' )
#'
#' ## Multinomial models: category probabilities and matrix offsets
#' d_multi <- iris
#' d_multi$off_A <- 0.1
#' d_multi$off_B <- 0.05 * d_multi$Sepal.Width
#' d_multi$off_C <- -0.05 * d_multi$Sepal.Width
#'
#' fit_multinomial <- mfp2(
#'   Species ~ fp(Sepal.Length) +
#'     offset(cbind(
#'       setosa = off_A,
#'       versicolor = off_B,
#'       virginica = off_C
#'     )),
#'   data = d_multi,
#'   family = multinomial_family(reference = "setosa"),
#'   df = 1,
#'   select = 1,
#'   verbose = FALSE
#' )
#'
#' nd_multi <- d_multi[
#'   c(1, 51, 101),
#'   c("Sepal.Length", "off_A", "off_B", "off_C")
#' ]
#'
#' # The formula rebuilds the offset matrix from these columns.
#' predict(
#'   fit_multinomial,
#'   newdata = nd_multi,
#'   type = "response",
#'   se.fit = TRUE
#' )
#'
#' # Equivalent offsets for only the non-reference categories.
#' # Subtract the setosa offset; do not simply drop its column.
#' relative_offsets <- with(
#'   nd_multi,
#'   cbind(
#'     versicolor = off_B - off_A,
#'     virginica = off_C - off_A
#'   )
#' )
#'
#' # Explicit newoffset replaces the formula offset, so its source columns
#' # are no longer needed in newdata.
#' predict(
#'   fit_multinomial,
#'   newdata = nd_multi["Sepal.Length"],
#'   newoffset = relative_offsets,
#'   type = "response"
#' )
#'
#' # Class predictions return labels; do not request their standard errors.
#' predict(fit_multinomial, newdata = nd_multi, type = "class")
#'
#' ## Ordinal models: category probabilities and cutoff intercepts
#' set.seed(3)
#' d_ordinal <- data.frame(age = runif(150, 30, 75))
#' d_ordinal$severity <- cut(
#'   0.04 * (d_ordinal$age - 50) + rlogis(150),
#'   breaks = c(-Inf, -0.5, 0.5, Inf),
#'   labels = c("low", "medium", "high"),
#'   ordered_result = TRUE
#' )
#'
#' fit_ordinal <- mfp2(
#'   severity ~ fp(age),
#'   data = d_ordinal,
#'   family = ordinal_family(),
#'   df = 1,
#'   select = 1,
#'   verbose = FALSE
#' )
#'
#' nd_ordinal <- data.frame(age = c(40, 50, 60))
#' predict(
#'   fit_ordinal,
#'   newdata = nd_ordinal,
#'   type = "response",
#'   se.fit = TRUE
#' )
#'
#' # The common predictor without a cutoff intercept.
#' predict(
#'   fit_ordinal,
#'   newdata = nd_ordinal,
#'   type = "link",
#'   add_intercept = FALSE
#' )
#'
#' # Include cutoff 1: medium/high versus low.
#' predict(
#'   fit_ordinal,
#'   newdata = nd_ordinal,
#'   type = "terms",
#'   terms = "age",
#'   terms_seq = "data",
#'   ordinal_cutoff = 1
#' )
#'
#' ## Cox models: relative risk and survival at a chosen time
#' data("gbsg")
#' fit_cox <- mfp2(
#'   survival::Surv(rectime, censrec) ~
#'     fp(age) + fp(nodes) + strata(grade),
#'   data = gbsg,
#'   family = "cox",
#'   df = 1,
#'   select = 1,
#'   verbose = FALSE
#' )
#'
#' # Relative predictions need predictors and the strata variable.
#' nd_cox <- gbsg[1:3, c("age", "nodes", "grade")]
#' predict(fit_cox, newdata = nd_cox, type = "lp")
#' predict(
#'   fit_cox,
#'   newdata = nd_cox,
#'   type = "risk",
#'   cox_reference = "strata"
#' )
#'
#' # Predict survival at 365 days for each profile.
#' # Use the response-variable names from the fitted formula.
#' nd_survival <- nd_cox
#' nd_survival$rectime <- 365
#' nd_survival$censrec <- 0L
#'
#' # censrec is required to reconstruct Surv(); its value is a placeholder
#' # here, not knowledge of whether a future event will occur.
#' predict(fit_cox, newdata = nd_survival, type = "survival")
#'
#' # "expected" returns cumulative hazard, not an event probability.
#' predict(fit_cox, newdata = nd_survival, type = "expected")
#'
#' ## Parametric survival models: predicted survival-time quantiles
#' fit_survreg <- mfp2(
#'   survival::Surv(rectime, censrec) ~ fp(age) + fp(nodes),
#'   data = gbsg,
#'   family = survreg_family(dist = "weibull"),
#'   df = 1,
#'   select = 1,
#'   verbose = FALSE
#' )
#'
#' # Median survival time, in days. New follow-up values are not needed.
#' predict(
#'   fit_survreg,
#'   newdata = gbsg[1:3, c("age", "nodes")],
#'   type = "quantile",
#'   p = 0.5
#' )
#'
#' ## Fine--Gray models: cumulative incidence at chosen times
#' set.seed(4)
#' d_fg <- data.frame(
#'   time = rexp(180, rate = 1 / 1000),
#'   age = runif(180, 30, 75),
#'   event = factor(
#'     sample(c("censored", "relapse", "death"), 180, replace = TRUE),
#'     levels = c("censored", "relapse", "death")
#'   )
#' )
#'
#' # The first event level represents censoring.
#' fit_finegray <- mfp2(
#'   survival::Surv(time, event) ~ fp(age, select = 1),
#'   data = d_fg,
#'   family = finegray_family(etype = "relapse"),
#'   verbose = FALSE
#' )
#'
#' # Probability of relapse by 365 and 730 days, allowing for competing death.
#' # Unlike Cox survival prediction, supply times through the times argument.
#' predict(
#'   fit_finegray,
#'   newdata = data.frame(age = c(40, 50, 60)),
#'   type = "response",
#'   times = c(365, 730)
#' )
#' }
#' @seealso
#' [mfp2()], [get_selected_variable_names()], [stats::predict.glm()],
#' [survival::predict.coxph()], [geepack::geeglm()]
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
                         times = NULL,
                         nseq = 100,
                         add_intercept = TRUE,
                         cox_reference = NULL,
                         ordinal_cutoff = NULL,
                         ...) {

  # Record whether the caller supplied a non-NULL Cox covariate reference.
  # NULL follows predict.mfpi() semantics and means "use the method default".
  # The explicit-supply distinction is essential for expected/survival
  # predictions, where the underlying Cox method fixes its own sample reference.
  cox_reference_supplied <- !missing(cox_reference) && !is.null(cox_reference)

  # `newy` was considered during development but is intentionally not part of
  # the public API. Absolute Cox prediction data belong in `newdata`, matching
  # the data-oriented convention used by predict.coxph(). Catch legacy draft
  # calls here so the argument cannot disappear silently into `...`.
  dots_call <- match.call(expand.dots = FALSE)[["..."]]
  dots_names <- names(as.list(dots_call))
  if (!is.null(dots_names) && "reference" %in% dots_names) {
    stop(
      "'reference' has been renamed to 'cox_reference' in predict.mfp2().",
      call. = FALSE
    )
  }
  if (!is.null(dots_names) && "newy" %in% dots_names) {
    stop(
      "'newy' has been removed. Supply the Cox prediction response in ",
      "'newdata': include the original formula response variables or one ",
      "Surv column.",
      call. = FALSE
    )
  }

  terms_seq <- match.arg(terms_seq)

  # Validate scalar prediction controls before any fitted-object work. Reusing
  # shared validators gives NA-, length-, finiteness-, and integer-safe errors
  # without changing valid prediction behavior.
  validate_open_probability_scalar(alpha, "alpha")
  validate_positive_integer_scalar(nseq, "nseq")
  validate_logical_vector(se.fit, "se.fit", allowed_lengths = 1L)
  validate_logical_vector(add_intercept, "add_intercept", allowed_lengths = 1L)
  nseq <- as.integer(nseq)

  # Prediction metadata and fitted-model methods are defined only for mfp2 objects.
  if (!inherits(object, "mfp2")) {
    stop("The object is not an mfp2 object.", call. = FALSE)
  }

  # Resolve the prediction scale. Formula fits store the same conceptual term
  # names used by selection and printing, including source-variable names for
  # simple wrappers such as factor(x9).
  if (is.null(type)) {
    type <- ifelse(mfp2_family_is_ph(object$family_string), "lp", "link")
  }

  if (!is.character(type) || length(type) != 1L || is.na(type)) {
    stop("'type' must be a single non-missing character value.", call. = FALSE)
  }

  # `link` and `lp` are family-specific names for the same linear-predictor
  # quantity. Normalize that unambiguous alias before method-specific type
  # matching so users receive the intended prediction rather than a downstream
  # match.arg() error from the base model method.
  type <- normalize_prediction_type(
    type = type,
    family_string = object$family_string
  )

  if (!is.null(ordinal_cutoff)) {
    n_cutoffs <- length(object$mfp2_ordinal_intercepts)
    if (!identical(object$family_string, "ordinal") ||
        !identical(type, "terms") ||
        !is.numeric(ordinal_cutoff) || length(ordinal_cutoff) != 1L ||
        is.na(ordinal_cutoff) || !is.finite(ordinal_cutoff) ||
        ordinal_cutoff != floor(ordinal_cutoff) ||
        ordinal_cutoff < 1L || ordinal_cutoff > n_cutoffs) {
      stop("`ordinal_cutoff` must select an available cutoff for ordinal `type = 'terms'`.",
           call. = FALSE)
    }
  }

  # Cox prediction types are resolved here rather than left to
  # predict.coxph(). This lets mfp2 distinguish its own term/contrast paths from
  # the relative and absolute full-model paths and report argument conflicts
  # before any transformed data are constructed.
  resolved_cox_reference <- NULL
  if (mfp2_family_is_ph(object$family_string)) {
    type <- if (identical(object$family_string, "finegray")) {
      mfp2_match_finegray_prediction_type(type)
    } else {
      mfp2_match_cox_prediction_type(type)
    }
    resolved_cox_reference <- mfp2_validate_cox_reference(
      type = type,
      cox_reference = cox_reference,
      cox_reference_supplied = cox_reference_supplied
    )
  } else if (identical(object$family_string, "survreg")) {
    type <- mfp2_match_survreg_prediction_type(type)
    if (cox_reference_supplied) {
      stop("'cox_reference' is only available for proportional-hazards models.", call. = FALSE)
    }
  } else if (identical(object$family_string, "multinomial")) {
    type <- mfp2_match_multinomial_prediction_type(type)
    if (cox_reference_supplied) {
      stop("'cox_reference' is only available for Cox models.", call. = FALSE)
    }
  } else if (identical(object$family_string, "ordinal")) {
    type <- mfp2_match_ordinal_prediction_type(type)
    if (cox_reference_supplied) {
      stop("'cox_reference' is only available for Cox models.", call. = FALSE)
    }
  } else {
    type <- mfp2_match_glm_prediction_type(type)
    if (cox_reference_supplied) {
      stop(
        "'cox_reference' is only available for Cox models.",
        call. = FALSE
      )
    }
  }

  if (identical(object$family_string, "finegray") && identical(type, "response")) {
    if (!is.numeric(times) || length(times) < 1L || anyNA(times) ||
        any(!is.finite(times)) || any(times < 0)) {
      stop(
        "Fine--Gray `type = \"response\"` requires `times` to contain one or ",
        "more finite non-negative values.",
        call. = FALSE
      )
    }
    times <- as.numeric(times)
  } else if (!is.null(times)) {
    stop("`times` is used only for Fine--Gray `type = \"response\"` prediction.",
         call. = FALSE)
  }

  # `strata` and `newoffset` alter only full-model predictions on
  # supplied rows. Reject irrelevant combinations before formula evaluation or
  # transformation so that an accepted argument is never silently ignored.
  full_model_prediction <- !type %in% c("terms", "contrasts")

  if (!is.null(strata)) {
    if (!object$family_string %in% c("cox", "survreg", "finegray")) {
      stop("'strata' is available only for Cox, Fine--Gray, or stratified survreg predictions.", call. = FALSE)
    }
    if (!full_model_prediction) {
      stop(
        "'strata' is not used for mfp2 term or contrast predictions.",
        call. = FALSE
      )
    }
    if (is.null(newdata)) {
      stop(
        "'strata' can be supplied only together with 'newdata'.",
        call. = FALSE
      )
    }
    if (prediction_has_formula_strata(object)) {
      stop(
        "`strata` cannot be supplied when predicting from an `mfp2` model ",
        "fitted with formula-level strata. Include the original strata ",
        "variable(s) in `newdata` instead.",
        call. = FALSE
      )
    }
    if (!prediction_fit_has_strata(object)) {
      stop(
        "'strata' was supplied, but the retained survival model is not stratified.",
        call. = FALSE
      )
    }
  }

  if (!is.null(newoffset)) {
    if (!full_model_prediction) {
      stop(
        "'newoffset' is not used for mfp2 term or contrast predictions.",
        call. = FALSE
      )
    }
    if (is.null(newdata)) {
      stop(
        "'newoffset' can be supplied only together with 'newdata'.",
        call. = FALSE
      )
    }
    if (!isTRUE(object$has_offset)) {
      stop(
        "'newoffset' can be supplied only when the fitted model used an offset.",
        call. = FALSE
      )
    }
    if (!is.numeric(newoffset) || anyNA(newoffset) ||
        any(!is.finite(newoffset))) {
      stop(
        "'newoffset' must contain finite numeric values.",
        call. = FALSE
      )
    }
    offset_rows <- if (is.matrix(newoffset)) nrow(newoffset) else length(newoffset)
    if (offset_rows != NROW(newdata)) {
      stop(
        "The number of observation rows in 'newoffset' must equal the rows in 'newdata'.",
        call. = FALSE
      )
    }
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

  # Data-frame subclasses such as tibbles return data frames from `[, j]`,
  # which the term and contrast code expects to be vectors. Use a plain
  # data.frame throughout.
  if (is.data.frame(newdata) && !identical(class(newdata), "data.frame")) {
    newdata <- as.data.frame(newdata)
    class(newdata) <- "data.frame"
  }

  newdata_raw <- newdata

  # Check row count before constructing ranges or references, without inspecting
  # unused columns. Required predictors and formula inputs are validated below;
  # offsets, strata, and absolute Cox responses are checked when they are used.
  if (!is.null(newdata_raw)) {
    mfp2_validate_prediction_newdata(
      newdata = newdata_raw[, integer(0L), drop = FALSE]
    )
  }

  # Absolute Cox predictions with newdata need a Surv response because
  # predict.coxph() evaluates the cumulative baseline hazard at the supplied
  # time or interval. Reconstruct it before formula newdata are reduced to the
  # active predictor design columns.
  cox_prediction_response <- if (
    !is.null(newdata_raw) &&
    identical(object$family_string, "cox") &&
    type %in% c("expected", "survival")
  ) {
    reconstruct_cox_prediction_response(
      object = object,
      fit_obj = object,
      newdata = newdata_raw
    )
  } else {
    NULL
  }

  prediction_terms <- if (type %in% c("terms", "contrasts")) {
    intersect(terms, selected_internal_terms)
  } else {
    selected_internal_terms
  }

  # Zero-handling flags also exist for unrequested or eliminated terms. Those
  # columns must be ignored just like any other unused prediction input.
  zero_option_terms <- intersect(
    prediction_terms,
    active_zero_option_names(object$zero, object$catzero, object$spike)
  )

  if (!is.null(newdata_raw)) {
    validate_zero_option_covariates(newdata_raw, zero_option_terms)
  }

  # Explicit contrast references use the same exact-zero domain as newdata.
  # Validate them before adding fitted shifts or constructing transformed rows.
  if (!is.null(ref) && length(ref) > 0L) {
    for (term in intersect(names(ref), zero_option_terms)) {
      value <- ref[[term]]
      if (!is.null(value) && is.numeric(value)) {
        validate_zero_option_covariates(
          matrix(value, ncol = 1L, dimnames = list(NULL, term)),
          term
        )
      }
    }
  }

  if (!is.null(newdata)) {
    newdata <- reconstruct_formula_newdata(
      object,
      newdata,
      terms = prediction_terms
    )
  }

  if (!is.null(newdata)) {
    # Validate only the reconstructed predictors used by the requested
    # prediction. Irrelevant extra columns, including missing or infinite values,
    # are intentionally ignored for both fitting interfaces.
    validation_lookup <- prediction_term_to_columns(
      object,
      terms = prediction_terms
    )
    validation_columns <- unique(unlist(validation_lookup, use.names = FALSE))
    validation_columns <- intersect(validation_columns, colnames(newdata))
    mfp2_validate_prediction_newdata(
      newdata[, validation_columns, drop = FALSE]
    )
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
    if (identical(object$family_string, "multinomial")) {
      return(mfp2_predict_multinomial_terms(
        object,
        terms = terms,
        type = type,
        ref = ref,
        terms_seq = terms_seq,
        nseq = nseq,
        newdata = newdata,
        alpha = alpha,
        add_intercept = add_intercept,
        selected_internal_terms = selected_internal_terms
      ))
    }
    requested_terms <- terms
    internal_terms <- requested_terms
    selected <- internal_terms %in% selected_internal_terms

    # Name each omitted term, including unknown names and eliminated terms.
    omitted_terms <- unique(requested_terms[!selected])
    if (length(omitted_terms) > 0L) {
      action <- if (any(selected)) {
        "Omitting these terms."
      } else {
        "Returning an empty list."
      }
      warning(
        "predict(): requested term(s) not retained in the final model: ",
        paste(paste0("`", omitted_terms, "`"), collapse = ", "), ".\n",
        action,
        call. = FALSE
      )
    }
    if (!any(selected)) return(list())

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

    # Ordinal predictors share a slope but have one intercept per cutoff.
    # With no selected cutoff they return the shared, intercept-free effect.
    cf <- coef(object)
    intercept_name <- if (!is.null(ordinal_cutoff)) {
      names(object$mfp2_ordinal_intercepts)[ordinal_cutoff]
    } else {
      "(Intercept)"
    }
    intercept <- if (intercept_name %in% names(cf)) cf[[intercept_name]] else 0
    if (is.na(intercept) ||
        (is.null(ordinal_cutoff) &&
         (identical(object$family_string, "ordinal") || !add_intercept)) ||
        type == "contrasts") {
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
          terms = t,
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
            terms = t,
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
            terms = t,
            apply_pre = FALSE,
            check_binary = FALSE,
            allow_missing_predictors = TRUE
          ))
        } else {
          x_seq <- if (!is.null(newdata)) {
            # Use column indexing rather than newdata[[t]] so the documented
            # matrix and data-frame prediction interfaces follow the same path.
            matrix(
              newdata[, t] + object$transformations[t, "shift"],
              ncol = 1L,
              dimnames = list(NULL, t)
            )
          } else {
            object$x_original[, t, drop = FALSE]
          }

          # check_binary = FALSE: one or two distinct evaluation values must
          # still receive the fitted FP transformation.
          x_trafo <- as.matrix(prepare_newdata_for_predict(
            object,
            x_seq,
            terms = t,
            apply_pre = FALSE,
            check_binary = FALSE,
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

        # Legacy serialized fits may predate `spike_dec`. Treat absent or
        # incomplete decision metadata as not binary-only; the helper still
        # recognizes the explicit code stored by current SAZ fits.
        if (prediction_term_is_binary_only(object, t)) {
          variable <- as.integer(x_seq[, t] == 0)
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
            terms = t,
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
        include_intercept = type == "terms" &&
          (!is.null(ordinal_cutoff) ||
           (add_intercept && !identical(object$family_string, "ordinal"))),
        intercept_name = intercept_name
      )
      mult <- stats::qnorm(1 - alpha / 2)
      res$lower <- res$value - mult * res$se
      res$upper <- res$value + mult * res$se
      res_list[[display_term]] <- res
    }

    return(res_list)
  }

  # Full-model prediction reconstructs the required raw design columns
  # and fitted transformations, then uses the appropriate family-specific
  # prediction path.

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

      if (!is.numeric(newoffset) || anyNA(newoffset) ||
          any(!is.finite(newoffset))) {
        stop(
          "! newoffset must contain finite numeric values.",
          call. = FALSE
        )
      }

      offset_rows <- if (is.matrix(newoffset)) nrow(newoffset) else length(newoffset)
      if (offset_rows != nrow(newdata)) {
        stop(
          "The number of observation rows in newoffset must equal the rows of newdata",
          call. = FALSE
        )
      }

    } else {
      newoffset <- NULL
    }

    if (
      object$family_string %in% c("cox", "survreg", "finegray") &&
      !is.null(newdata_raw) &&
      prediction_has_formula_strata(object) &&
      !is.null(object$mfp2_internal_names$strata)
    ) {
      strata <- reconstruct_formula_strata_newdata(object, newdata_raw)
    }

    newdata <- prepare_newdata_for_predict(
      object,
      newdata,
      terms = selected_internal_terms,
      strata = strata,
      # Multinomial and ordinal predictors apply `newoffset` themselves and
      # have no formula offset helper column to attach it to.
      offset = if (object$family_string %in% c("multinomial", "ordinal")) {
        NULL
      } else {
        newoffset
      },
      check_binary = FALSE
    )

    if (!is.null(cox_prediction_response)) {
      newdata <- attach_cox_prediction_response(
        fit_obj = object,
        newdata = newdata,
        response = cox_prediction_response
      )
    }

    if (identical(object$family_string, "multinomial")) {
      return(mfp2_predict_multinomial(
        object = object,
        transformed = newdata,
        type = type,
        se.fit = se.fit,
        newoffset = newoffset
      ))
    }

    if (identical(object$family_string, "ordinal")) {
      return(mfp2_predict_ordinal(
        object = object,
        transformed = newdata,
        type = type,
        se.fit = se.fit,
        newoffset = newoffset,
        add_intercept = add_intercept
      ))
    }

    if (mfp2_family_is_ph(object$family_string)) {
      if (identical(object$family_string, "finegray") &&
          identical(type, "response")) {
        return(mfp2_predict_finegray_cif(
          object = object,
          newdata = newdata,
          times = times,
          se.fit = se.fit
        ))
      }
      return(
        mfp2_predict_cox_base(
          object = object,
          newdata = newdata,
          type = type,
          se.fit = se.fit,
          cox_reference = resolved_cox_reference,
          ...
        )
      )
    }

    if (inherits(object, "fastglm")) {
      prediction_data <- as.data.frame(newdata, check.names = FALSE)
      # fastglm final fits are matrix based, so their offset is passed directly
      # and no temporary formula column or internal-name metadata is needed.
      prediction_offset <- if (is.null(newoffset)) {
        rep(0, nrow(prediction_data))
      } else {
        as.numeric(newoffset)
      }
      prediction_matrix <- as.matrix(prediction_data)
      model_columns <- prediction_model_column_names(
        object,
        colnames(prediction_matrix)
      )
      colnames(prediction_matrix) <- model_columns

      return(
        mfp2_predict_fastglm_matrix(
          object = object,
          newx = prediction_matrix,
          offset = prediction_offset,
          type = type,
          se.fit = se.fit,
          ...
        )
      )
    }

    # GEE standard errors must use the covariance selected by `std.err`, which the
    # base predict.glm() path does not (and cannot) provide. Route se.fit for
    # GEE through the GEE matrix predictor; point predictions still use the
    # base delegation below, which is correct because it depends only on beta.
    if (identical(object$family_string, "gee") && isTRUE(se.fit)) {
      prediction_data <- as.data.frame(newdata, check.names = FALSE)
      prediction_offset <- if (is.null(newoffset)) {
        rep(0, nrow(prediction_data))
      } else {
        as.numeric(newoffset)
      }
      # prepare_newdata_for_predict() also attaches the offset as the internal
      # formula helper column used by predict.glm(). The matrix predictor
      # receives the same values through `offset`, so drop that helper column
      # rather than treating it as a coefficient column.
      offset_column <- object$mfp2_internal_names$offset
      if (!is.null(newoffset) && is.character(offset_column) &&
          length(offset_column) == 1L) {
        prediction_data[[offset_column]] <- NULL
      }
      prediction_matrix <- as.matrix(prediction_data)
      colnames(prediction_matrix) <- prediction_model_column_names(
        object, colnames(prediction_matrix)
      )
      prediction <- mfp2_predict_gee_matrix(
        object = object,
        newx = prediction_matrix,
        offset = prediction_offset,
        type = type,
        se.fit = TRUE
      )
      # Match the row names returned by the se.fit = FALSE (predict.glm) path.
      names(prediction$fit) <- names(prediction$se.fit) <-
        row.names(prediction_data)
      return(prediction)
    }

    # Strip "mfp2" so S3 dispatch reaches the underlying prediction method
    # rather than recursively entering predict.mfp2().
    obj_base <- object
    class(obj_base) <- setdiff(class(obj_base), "mfp2")
    obj_base <- mfp2_restore_gee_family(obj_base)

    return(
      stats::predict(
        obj_base,
        newdata = newdata,
        type = type,
        se.fit = se.fit,
        ...
      )
    )
  }

  # With no newdata, the fitted object contains the training response,
  # transformed design, strata, and offsets needed for prediction.
  if (identical(object$family_string, "multinomial")) {
    return(mfp2_predict_multinomial(
      object = object,
      transformed = NULL,
      type = type,
      se.fit = se.fit,
      newoffset = NULL
    ))
  }

  if (identical(object$family_string, "ordinal")) {
    return(mfp2_predict_ordinal(
      object = object,
      transformed = NULL,
      type = type,
      se.fit = se.fit,
      newoffset = NULL,
      add_intercept = add_intercept
    ))
  }

  if (mfp2_family_is_ph(object$family_string)) {
    if (identical(object$family_string, "finegray") &&
        identical(type, "response")) {
      return(mfp2_predict_finegray_cif(
        object = object,
        newdata = NULL,
        times = times,
        se.fit = se.fit
      ))
    }
    return(
      mfp2_predict_cox_base(
        object = object,
        type = type,
        se.fit = se.fit,
        cox_reference = resolved_cox_reference,
        ...
      )
    )
  }

  if (inherits(object, "fastglm")) {
    return(
      mfp2_predict_fastglm_matrix(
        object = object,
        type = type,
        se.fit = se.fit,
        ...
      )
    )
  }

  # GEE standard errors (no newdata) likewise use the selected covariance.
  if (identical(object$family_string, "gee") && isTRUE(se.fit)) {
    return(
      mfp2_predict_gee_matrix(
        object = object,
        newx = NULL,
        type = type,
        se.fit = TRUE
      )
    )
  }

  obj_base <- object
  class(obj_base) <- setdiff(class(obj_base), "mfp2")
  obj_base <- mfp2_restore_gee_family(obj_base)
  stats::predict(
    obj_base,
    type = type,
    se.fit = se.fit,
    ...
  )
}

#' Prepare a GEE Family for GLM Prediction
#'
#' Ensures that a GEE object stores the response family expected by
#' [stats::predict.glm()]. If its `family` slot contains an `mfp2_gee_family`
#' wrapper, replaces it with the inner GLM family. Otherwise, returns the
#' object unchanged.
#'
#' @param object A fitted object with the `mfp2` class removed, ready for
#'   prediction through its underlying method.
#'
#' @return The object with its inner GLM family restored when needed.
#'
#' @keywords internal
#' @noRd
mfp2_restore_gee_family <- function(object) {
  if (identical(object$family_string, "gee")) {
    inner <- object$family$response_family
    if (!inherits(inner, "family") &&
        inherits(object$mfp2_family, "mfp2_gee_family")) {
      inner <- object$mfp2_family$response_family
    }
    if (inherits(inner, "family")) {
      object$family <- inner
    }
  }
  object
}


#' Validate and Match a Cox Prediction Type
#'
#' Matches an exact or unambiguous abbreviated prediction type to its full
#' name before selecting the prediction method or transforming prediction
#' data. Invalid or ambiguous values produce an error.
#'
#' @param type A single character value supplied to `predict.mfp2()`.
#'
#' @return One of `"lp"`, `"risk"`, `"expected"`, `"survival"`,
#'   `"terms"`, or `"contrasts"`.
#'
#' @keywords internal
#' @noRd
mfp2_match_cox_prediction_type <- function(type) {
  choices <- c("lp", "risk", "expected", "survival", "terms", "contrasts")

  tryCatch(
    match.arg(type, choices),
    error = function(e) {
      stop(
        "For Cox models, 'type' must be one of: ",
        paste(shQuote(choices), collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  )
}

#' Validate and Match a GLM Prediction Type
#'
#' Matches an exact or unambiguous abbreviated prediction type to its full
#' name before transforming prediction data. Invalid or ambiguous values
#' produce an error. Aliases are handled by `normalize_prediction_type()`
#' before this function is called.
#'
#' @param type A single character value after alias normalization.
#'
#' @return One of `"link"`, `"response"`, `"terms"`, or `"contrasts"`.
#'
#' @keywords internal
#' @noRd
mfp2_match_glm_prediction_type <- function(type) {
  choices <- c("link", "response", "terms", "contrasts")

  if (identical(type, "risk")) {
    stop(
      "For GLM models, type = 'risk' is not supported. Use type = ",
      "'response' only when prediction on the model response scale is ",
      "intended.",
      call. = FALSE
    )
  }

  tryCatch(
    match.arg(type, choices),
    error = function(e) {
      stop(
        "For GLM models, 'type' must be one of: ",
        paste(shQuote(choices), collapse = ", "),
        ". The alias 'lp' is also accepted for 'link'.",
        call. = FALSE
      )
    }
  )
}

#' Validate and Match a Multinomial Prediction Type
#'
#' Matches an exact or unambiguous abbreviated prediction type to its full
#' name before prediction. Invalid or ambiguous values produce an error.
#'
#' @param type A single character value supplied to `predict.mfp2()`.
#'
#' @return One of `"link"`, `"response"`, `"class"`, `"terms"`,
#'   or `"contrasts"`.
#'
#' @keywords internal
#' @noRd
mfp2_match_multinomial_prediction_type <- function(type) {
  choices <- c("link", "response", "class", "terms", "contrasts")
  tryCatch(
    match.arg(type, choices),
    error = function(e) {
      stop(
        "For multinomial models, 'type' must be one of: ",
        paste(shQuote(choices), collapse = ", "), ".",
        call. = FALSE
      )
    }
  )
}

#' Validate and Match an Ordinal Prediction Type
#'
#' Matches an exact or unambiguous abbreviated prediction type to its full
#' name before prediction. Invalid or ambiguous values produce an error.
#' The alias `"lp"` is converted to `"link"` before this function is called.
#'
#' @param type A single character value after alias normalization.
#'
#' @return One of `"link"`, `"response"`, `"mean"`, `"terms"`,
#'   or `"contrasts"`.
#'
#' @keywords internal
#' @noRd
mfp2_match_ordinal_prediction_type <- function(type) {
  # `lp` is normalized to `link` upstream for non-PH families.
  choices <- c("link", "response", "mean", "terms", "contrasts")
  tryCatch(
    match.arg(type, choices),
    error = function(e) {
      stop(
        "For ordinal models, 'type' must be one of: ",
        paste(shQuote(choices), collapse = ", "),
        ". The alias 'lp' is also accepted for 'link'.",
        call. = FALSE
      )
    }
  )
}

#' Compute Full-Model Ordinal Predictions
#'
#' Returns link predictions, category probabilities, or expected category
#' values. Link predictions can include the reference-cutoff intercept;
#' probabilities and means always use all fitted cutoffs.
#' Standard errors account for uncertainty in the relevant coefficients
#' and their covariances.
#'
#' @param object A fitted ordinal `mfp2` model.
#' @param transformed Prepared prediction data with the fitted transformations
#'   applied, or `NULL` to use the fitting data.
#' @param type One of `"link"`, `"response"`, or `"mean"`.
#' @param se.fit Logical; whether to return prediction standard errors.
#' @param newoffset Offsets for new rows. When `transformed = NULL`, the
#'   fitting offsets are used instead.
#' @param add_intercept Logical; whether link predictions include the
#'   reference-cutoff intercept. Does not affect probabilities or means.
#'
#' @return A numeric vector for link or mean predictions, or a matrix with
#'   one column per response category for probabilities. With `se.fit = TRUE`,
#'   a list containing `fit` and `se.fit` with matching dimensions and names.
#'
#' @keywords internal
#' @noRd
mfp2_predict_ordinal <- function(object, transformed = NULL,
                                 type = c("link", "response", "mean"),
                                 se.fit = FALSE, newoffset = NULL,
                                 add_intercept = TRUE) {
  type <- match.arg(type)

  intercepts <- object$mfp2_ordinal_intercepts
  n_int <- length(intercepts)
  link <- object$mfp2_ordinal_link
  levels_labels <- object$mfp2_ordinal_levels
  all_coef <- object$coefficients
  betas <- if (length(all_coef) > n_int) {
    all_coef[seq.int(n_int + 1L, length(all_coef))]
  } else {
    numeric(0L)
  }

  if (is.null(transformed)) {
    xx <- object$x
    prediction_rownames <- rownames(xx)
    if (isTRUE(object$mfp2_family$prepared$has_offset)) {
      newoffset <- object$mfp2_family$prepared$offset
    }
  } else {
    transformed <- as.data.frame(transformed, check.names = FALSE)
    prediction_rownames <- row.names(transformed)
    model_columns <- prediction_model_column_names(object, colnames(transformed))
    names(model_columns) <- colnames(transformed)
    needed <- names(betas)
    selected <- match(needed, unname(model_columns))
    if (anyNA(selected)) {
      stop("Internal error: transformed ordinal prediction columns are incomplete.", call. = FALSE)
    }
    xx <- as.matrix(transformed[, names(model_columns)[selected], drop = FALSE])
    colnames(xx) <- needed
  }

  n <- if (is.null(transformed)) {
    if (!is.null(xx)) nrow(xx) else object$nobs
  } else {
    nrow(transformed)
  }
  if (length(betas) == 0L) {
    xx <- matrix(numeric(0L), nrow = n, ncol = 0L,
                 dimnames = list(prediction_rownames, character(0L)))
  } else {
    if (is.null(xx)) {
      stop("Internal error: ordinal prediction design is unavailable.",
           call. = FALSE)
    }
    xx <- xx[, names(betas), drop = FALSE]
  }
  eta <- if (length(betas) > 0L) {
    as.numeric(xx %*% betas)
  } else {
    rep(0, n)
  }
  if (!is.null(newoffset)) eta <- eta + as.numeric(newoffset)
  names(eta) <- prediction_rownames
  if (type == "link") {
    fit <- eta
    X <- xx
    if (isTRUE(add_intercept)) {
      # predictrms() uses interceptRef, falling back to the first cutoff for
      # older fits. Do not choose a cutoff from the new prediction data.
      cutoff <- object$interceptRef
      if (is.null(cutoff) || length(cutoff) == 0L) cutoff <- 1L
      if (!is.numeric(cutoff) || length(cutoff) != 1L ||
          is.na(cutoff) || !is.finite(cutoff) || cutoff != floor(cutoff) ||
          cutoff < 1L || cutoff > n_int) {
        stop("The fitted ordinal reference cutoff is invalid.", call. = FALSE)
      }
      fit <- fit + unname(intercepts[cutoff])
      if (isTRUE(se.fit)) {
        V <- mfp2_ordinal_vcov_all(object)
        columns <- c(names(intercepts)[cutoff], names(betas))
        V <- V[columns, columns, drop = FALSE]
        X <- cbind(1, X)
      }
    } else if (isTRUE(se.fit)) {
      V <- mfp2_ordinal_vcov_slopes(object)
    }
    if (!isTRUE(se.fit)) return(fit)
    se <- mfp2_prediction_standard_error(X, V)
    names(se) <- names(fit)
    return(list(fit = fit, se.fit = se))
  }

  ginv <- switch(link,
    logistic = stats::plogis,
    probit   = stats::pnorm,
    loglog   = function(z) exp(-exp(-z)),
    cloglog  = function(z) -expm1(-exp(z)),
    cauchit  = stats::pcauchy
  )
  # Derivatives of the inverse links. The exponential forms avoid Inf * 0
  # when a finite but extreme predictor makes exp(z) or exp(-z) overflow.
  density <- switch(link,
    logistic = stats::dlogis,
    probit   = stats::dnorm,
    loglog   = function(z) exp(-z - exp(-z)),
    cloglog  = function(z) exp(z - exp(z)),
    cauchit  = stats::dcauchy
  )
  k <- n_int + 1L
  cum <- matrix(1, n, k)                      # cum[, j] = P(Y >= j); P(Y>=1)=1
  for (j in seq_len(n_int)) cum[, j + 1L] <- ginv(intercepts[j] + eta)
  probs <- matrix(0, n, k)
  for (j in seq_len(k)) {
    probs[, j] <- cum[, j] - (if (j < k) cum[, j + 1L] else 0)
  }
  dimnames(probs) <- list(prediction_rownames, levels_labels)
  if (type == "response" && !isTRUE(se.fit)) return(probs)

  # type == "mean": expected response using numeric level values when available,
  # otherwise the 1..k category index.
  numeric_levels <- suppressWarnings(as.numeric(levels_labels))
  vals <- if (!anyNA(numeric_levels)) numeric_levels else seq_len(k)
  fit <- if (type == "response") probs else {
    stats::setNames(as.numeric(probs %*% vals), prediction_rownames)
  }
  if (!isTRUE(se.fit)) return(fit)

  V <- mfp2_ordinal_vcov_all(object)
  derivatives <- matrix(0, nrow = n, ncol = n_int)
  for (j in seq_len(n_int)) {
    derivatives[, j] <- density(intercepts[j] + eta)
  }
  if (type == "mean") {
    # E[Y] = value[1] + sum_j diff(value)[j] * P(Y >= j + 1).
    # Slopes affect every cutoff; include their covariance with all intercepts.
    threshold_gradient <- sweep(derivatives, 2L, diff(vals), "*")
    gradient <- cbind(threshold_gradient, xx * rowSums(threshold_gradient))
    se <- mfp2_prediction_standard_error(gradient, V)
    names(se) <- names(fit)
  } else {
    se <- matrix(0, nrow = n, ncol = k, dimnames = dimnames(probs))
    for (j in seq_len(k)) {
      # P(Y = j) is the difference of its two adjacent cumulative probabilities.
      # Only those cutoff coefficients enter this gradient. Restricting the
      # covariance block avoids multiplying a mostly zero, potentially large
      # matrix for ordinal responses with many categories.
      cutoffs <- seq.int(max(1L, j - 1L), min(n_int, j))
      threshold_gradient <- sweep(derivatives[, cutoffs, drop = FALSE], 2L,
                                  ifelse(cutoffs == j, -1, 1), "*")
      gradient <- cbind(threshold_gradient, xx * rowSums(threshold_gradient))
      columns <- c(cutoffs, n_int + seq_along(betas))
      se[, j] <- mfp2_prediction_standard_error(gradient, V[columns, columns, drop = FALSE])
    }
  }
  list(fit = fit, se.fit = se)
}

#' Compute Multinomial Term and Contrast Predictions
#'
#' Computes each requested term's contribution to the log-odds of each
#' non-reference category versus the reference category. Terms share the same
#' fitted transformations across categories but have separate coefficients.
#'
#' Term predictions include each category's intercept when
#' `add_intercept = TRUE`. Contrasts subtract the term's value at the
#' reference specified by `ref` and omit intercepts. Standard errors use
#' the covariance block for the corresponding category comparison.
#'
#' @param object A fitted multinomial `mfp2` model.
#' @param terms,type,ref,terms_seq,nseq,newdata,alpha,add_intercept Arguments
#'   passed from [predict.mfp2()], with the same meaning.
#' @param selected_internal_terms Internal names of terms retained in the
#'   final model, as determined by the caller.
#'
#' @return A named list with one data frame per requested term. Each contains
#'   `variable`, `variable_pre`, `value`, `se`, `lower`, and `upper`, plus a
#'   `logit` factor identifying the non-reference category. Factor levels
#'   follow the coefficient-matrix row order.
#'
#' @keywords internal
#' @noRd
mfp2_predict_multinomial_terms <- function(object,
                                           terms,
                                           type,
                                           ref,
                                           terms_seq,
                                           nseq,
                                           newdata,
                                           alpha,
                                           add_intercept,
                                           selected_internal_terms) {
  requested_terms <- terms
  internal_terms <- requested_terms
  selected <- internal_terms %in% selected_internal_terms

  # Name each omitted term, including unknown names and eliminated terms.
  omitted_terms <- unique(requested_terms[!selected])
  if (length(omitted_terms) > 0L) {
    action <- if (any(selected)) {
      "Omitting these terms."
    } else {
      "Returning an empty list."
    }
    warning(
      "predict(): requested term(s) not retained in the final model: ",
      paste(paste0("`", omitted_terms, "`"), collapse = ", "), ".\n",
      action,
      call. = FALSE
    )
  }
  if (!any(selected)) return(list())

  requested_terms <- requested_terms[selected]
  internal_terms <- internal_terms[selected]

  coefficient_matrix <- object$mfp2_coefficient_matrix
  if (!is.matrix(coefficient_matrix)) {
    stop("Internal error: multinomial coefficient matrix is unavailable.",
         call. = FALSE)
  }
  logits <- rownames(coefficient_matrix)
  covariance <- as.matrix(stats::vcov(object))
  mult <- stats::qnorm(1 - alpha / 2)

  # Standard error of a per-logit partial predictor: sqrt(diag(X V_q X')), where
  # V_q is the block of the multinomial covariance for logit `q`, addressed by
  # the "<logit>:<col>" naming that vcov() uses for the flattened coefficients.
  se_for_logit <- function(x_trafo, x_ref_trafo, q, include_intercept) {
    model_columns <- prediction_model_column_names(object, colnames(x_trafo))
    X <- x_trafo
    if (!is.null(x_ref_trafo)) {
      X <- sweep(X, 2, x_ref_trafo, "-")
    }
    cov_columns <- paste0(q, ":", model_columns)
    if (include_intercept && is.null(x_ref_trafo)) {
      X <- cbind("(Intercept)" = 1, X)
      cov_columns <- c(paste0(q, ":(Intercept)"), cov_columns)
    }
    ind <- match(cov_columns, colnames(covariance))
    if (anyNA(ind)) {
      stop(
        "Cannot compute multinomial SE; covariance matrix lacks column(s): ",
        paste(cov_columns[is.na(ind)], collapse = ", "),
        call. = FALSE
      )
    }
    block <- covariance[ind, ind, drop = FALSE]
    # Use the same finite-value and roundoff checks as scalar term predictions.
    mfp2_prediction_standard_error(X, block)
  }

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

    x_ref_trafo <- NULL

    if (block_term) {
      # Grouped/factor terms are predicted as complete design blocks, exactly as
      # in the scalar path; only the coefficient row and covariance block differ.
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
        object, source_block, terms = t,
        apply_pre = source_needs_preprocessing,
        allow_missing_predictors = TRUE, check_binary = FALSE
      ))
      labels <- decode_grouped_term_values(object, term = t, block = source_block)
      base_frame <- data.frame(
        variable = labels, variable_pre = labels, source_block,
        check.names = FALSE, stringsAsFactors = FALSE
      )
      if (type == "contrasts") {
        reference_block <- resolve_grouped_reference(
          object, term = t,
          reference = if (display_term %in% names(ref)) ref[[display_term]] else ref[[t]],
          block = source_block
        )
        x_ref_trafo <- as.matrix(prepare_newdata_for_predict(
          object, reference_block, terms = t, apply_pre = TRUE,
          allow_missing_predictors = TRUE, check_binary = FALSE
        ))
      }
    } else {
      # Singleton numeric terms: evaluate the FP transformation over the data or
      # an equidistant grid, on the shifted-but-unscaled fitting scale.
      if (terms_seq == "equidistant") {
        if (!is.null(newdata)) {
          x_range <- range(newdata[, t]) + object$transformations[t, "shift"]
        } else {
          x_range <- range(object$x_original[, t])
        }
        x_seq <- matrix(
          seq(x_range[1L], x_range[2L], length.out = nseq),
          ncol = 1L, dimnames = list(NULL, t)
        )
      } else {
        x_seq <- if (!is.null(newdata)) {
          matrix(
            newdata[, t] + object$transformations[t, "shift"],
            ncol = 1L, dimnames = list(NULL, t)
          )
        } else {
          object$x_original[, t, drop = FALSE]
        }
      }
      x_trafo <- as.matrix(prepare_newdata_for_predict(
        object, x_seq, terms = t, apply_pre = FALSE, check_binary = FALSE,
        allow_missing_predictors = TRUE
      ))
      variable <- as.numeric(x_seq) - object$transformations[t, "shift"]
      variable_pre <- as.numeric(x_seq)
      if (prediction_term_is_binary_only(object, t)) {
        variable <- as.integer(x_seq[, t] == 0)
        variable_pre <- variable
      }
      base_frame <- data.frame(variable = variable, variable_pre = variable_pre)

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
          object, x_ref, terms = t, apply_pre = FALSE, check_binary = FALSE,
          reset_zero = FALSE, allow_missing_predictors = TRUE
        ))
      }
    }

    # Map the transformed columns to coefficient-matrix columns once; the
    # mapping is shared across logits.
    model_columns <- prediction_model_column_names(object, colnames(x_trafo))
    include_intercept <- add_intercept && type == "terms"

    per_logit <- vector("list", length(logits))
    for (li in seq_along(logits)) {
      q <- logits[[li]]
      term_coef <- coefficient_matrix[q, model_columns]
      term_coef[is.na(term_coef)] <- 0
      value <- as.numeric(x_trafo %*% term_coef)
      if (type == "contrasts") {
        value <- value - as.numeric(x_ref_trafo %*% term_coef)
      } else if (include_intercept) {
        intercept_q <- coefficient_matrix[q, "(Intercept)"]
        if (!is.na(intercept_q)) value <- value + intercept_q
      }
      se <- se_for_logit(x_trafo, x_ref_trafo, q, include_intercept)
      df_q <- base_frame
      df_q$logit <- q
      df_q$value <- value
      df_q$se <- se
      df_q$lower <- value - mult * se
      df_q$upper <- value + mult * se
      per_logit[[li]] <- df_q
    }
    res <- do.call(rbind, per_logit)
    res$logit <- factor(res$logit, levels = logits)
    rownames(res) <- NULL
    res_list[[display_term]] <- res
  }

  res_list
}

#' Compute Full-Model Multinomial Predictions
#'
#' Returns log-odds relative to the reference category, category probabilities,
#' or predicted category labels. Predictions include any fitted or supplied
#' offsets. The reference category's logit is fixed at zero.
#'
#' @param object A fitted multinomial `mfp2` model.
#' @param transformed Prepared prediction data with the fitted transformations
#'   applied, or `NULL` to use the fitting data.
#' @param type One of `"link"`, `"response"`, or `"class"`.
#' @param se.fit Logical; whether to return standard errors. Not supported
#'   for `"class"` predictions.
#' @param newoffset Offset matrix for new rows. When `transformed = NULL`,
#'   the fitting offsets are used instead.
#'
#' @return A matrix of non-reference-category log-odds or category
#'   probabilities, or a factor of predicted category labels. For numeric
#'   predictions with `se.fit = TRUE`, a list containing `fit` and `se.fit`
#'   matrices with matching dimensions and names.
#'
#' @keywords internal
#' @noRd
mfp2_predict_multinomial <- function(object, transformed = NULL,
                                      type = c("link", "response", "class"),
                                      se.fit = FALSE, newoffset = NULL) {
  type <- match.arg(type)
  if (isTRUE(se.fit) && type == "class") {
    stop(
      "`se.fit = TRUE` is not available for multinomial class labels; use `type = 'response'` for probability standard errors.",
      call. = FALSE
    )
  }
  coefficient_matrix <- object$mfp2_coefficient_matrix
  if (!is.matrix(coefficient_matrix)) {
    stop("Internal error: multinomial coefficient matrix is unavailable.", call. = FALSE)
  }

  if (is.null(transformed)) {
    xx <- object$mfp2_design
    offset_matrix <- object$mfp2_offset_matrix
    prediction_rownames <- rownames(xx)
  } else {
    transformed <- as.data.frame(transformed, check.names = FALSE)
    # Match predict.multinom(): returned matrices use the prediction data's row
    # names, including automatic "1", "2", ... names for an unnamed matrix.
    prediction_rownames <- row.names(transformed)
    needed <- setdiff(colnames(coefficient_matrix), "(Intercept)")
    model_columns <- prediction_model_column_names(object, colnames(transformed))
    names(model_columns) <- colnames(transformed)
    selected <- match(needed, unname(model_columns))
    if (anyNA(selected)) {
      stop("Internal error: transformed multinomial prediction columns are incomplete.", call. = FALSE)
    }
    xx <- cbind(
      "(Intercept)" = 1,
      as.matrix(transformed[, names(model_columns)[selected], drop = FALSE])
    )
    colnames(xx) <- c("(Intercept)", needed)
    offset_matrix <- mfp2_multinomial_offset(
      newoffset,
      object$mfp2_family,
      nobs = nrow(xx),
      has_offset = isTRUE(object$has_offset)
    )
  }

  if (is.null(offset_matrix)) {
    offset_matrix <- mfp2_multinomial_offset(
      NULL,
      object$mfp2_family,
      nobs = nrow(xx),
      has_offset = FALSE
    )
  }

  xx <- xx[, colnames(coefficient_matrix), drop = FALSE]
  eta <- xx %*% t(coefficient_matrix)
  eta <- eta + offset_matrix[, -1L, drop = FALSE]
  dimnames(eta) <- list(prediction_rownames, rownames(coefficient_matrix))
  if (isTRUE(se.fit)) {
    # Flatten in logit-major order, matching vcov.multinom(). Match both axes
    # by name so neither coefficient nor covariance storage order is assumed.
    covariance_names <- unlist(lapply(rownames(coefficient_matrix), function(q) {
      paste0(q, ":", colnames(coefficient_matrix))
    }), use.names = FALSE)
    V <- as.matrix(stats::vcov(object))
    if (!all(covariance_names %in% rownames(V)) ||
        !all(covariance_names %in% colnames(V))) {
      stop("The multinomial covariance matrix lacks required coefficient names.",
           call. = FALSE)
    }
    V <- V[covariance_names, covariance_names, drop = FALSE]
    p <- ncol(xx)
    n_logits <- nrow(coefficient_matrix)
  }
  if (type == "link") {
    if (!isTRUE(se.fit)) return(eta)
    se <- matrix(0, nrow = nrow(eta), ncol = ncol(eta), dimnames = dimnames(eta))
    for (q in seq_len(n_logits)) {
      block <- (q - 1L) * p + seq_len(p)
      se[, q] <- mfp2_prediction_standard_error(xx, V[block, block, drop = FALSE])
    }
    return(list(fit = eta, se.fit = se))
  }

  eta_full <- cbind(0, eta)
  colnames(eta_full) <- object$class_levels
  row_max <- apply(eta_full, 1L, max)
  exp_eta <- exp(eta_full - row_max)
  probabilities <- exp_eta / rowSums(exp_eta)

  if (isTRUE(se.fit)) {
    se <- matrix(0, nrow = nrow(probabilities), ncol = ncol(probabilities),
                 dimnames = dimnames(probabilities))
    for (j in seq_len(ncol(probabilities))) {
      # Softmax derivative: d P_j / d eta_q = P_j * (I[j = q] - P_q).
      # Probability SEs need cross-logit covariance, including for the reference
      # category. Keep one n-by-parameter gradient at a time, not an n-by-n matrix.
      gradient <- do.call(cbind, lapply(seq_len(n_logits), function(q) {
        xx * (probabilities[, j] * ((j == q + 1L) - probabilities[, q + 1L]))
      }))
      se[, j] <- mfp2_prediction_standard_error(gradient, V)
    }
  }

  original_levels <- object$mfp2_family$prepared$original_levels
  if (!is.null(original_levels) && all(original_levels %in% colnames(probabilities))) {
    probabilities <- probabilities[, original_levels, drop = FALSE]
    if (isTRUE(se.fit)) se <- se[, original_levels, drop = FALSE]
  }
  if (type == "response") {
    if (isTRUE(se.fit)) return(list(fit = probabilities, se.fit = se))
    return(probabilities)
  }

  factor(
    colnames(probabilities)[max.col(probabilities, ties.method = "first")],
    levels = colnames(probabilities)
  )
}

#' Validate and Match a survreg Prediction Type
#'
#' Matches prediction types supported by [survival::predict.survreg()],
#' together with the `"terms"` and `"contrasts"` predictions handled by
#' `predict.mfp2()` on the linear-predictor scale.
#'
#' @param type A single character value supplied to `predict.mfp2()`.
#'
#' @return The matched prediction type used internally.
#'
#' @keywords internal
#' @noRd
mfp2_match_survreg_prediction_type <- function(type) {
  choices <- c("link", "response", "quantile", "uquantile", "terms", "contrasts")
  tryCatch(
    match.arg(type, choices),
    error = function(e) {
      stop(
        "For survreg models, 'type' must be one of: ",
        paste(shQuote(choices), collapse = ", "),
        ". The aliases 'lp' and 'linear' are also accepted for 'link'.",
        call. = FALSE
      )
    }
  )
}

#' Validate and Match a Fine--Gray Prediction Type
#'
#' Matches a supported prediction type to the name used internally before
#' selecting the prediction method. Unsupported values produce an error.
#'
#' @param type A single character value supplied to `predict.mfp2()`.
#'
#' @return One of `"lp"`, `"risk"`, `"response"`, `"terms"`,
#'   or `"contrasts"`.
#'
#' @keywords internal
#' @noRd
mfp2_match_finegray_prediction_type <- function(type) {
  choices <- c("lp", "risk", "response", "terms", "contrasts")
  tryCatch(
    match.arg(type, choices),
    error = function(e) {
      stop(
        "For Fine--Gray models, 'type' must be one of: ",
        paste(shQuote(choices), collapse = ", "), ".",
        call. = FALSE
      )
    }
  )
}

#' Compute Fine--Gray Cumulative-Incidence Predictions
#'
#' Predicts the probability of the selected event occurring by each requested
#' time, accounting for competing events. Uses [survival::survfit()] on the
#' retained weighted Cox fit and calculates cumulative incidence as one minus
#' the fitted subdistribution survival curve.
#'
#' @param object A fitted Fine--Gray `mfp2` model.
#' @param newdata Prediction data, or `NULL` to use the stored fitting data.
#' @param times Finite, non-negative prediction times, in the same units as
#'   the fitted follow-up time. Sorted and reduced to unique values.
#' @param se.fit Logical; whether to return prediction standard errors.
#'
#' @return A named numeric vector for one time, or a matrix with prediction
#'   rows and one column per time. With `se.fit = TRUE`, a list containing
#'   `fit`, `se.fit`, `times`, and the selected `event`.
#'
#' @details
#' Predictions are computed in batches. If a batch result cannot be aligned
#' reliably, its rows are evaluated individually. Output rows retain the
#' order of the prediction data.
#'
#' @keywords internal
#' @noRd
mfp2_predict_finegray_cif <- function(object, newdata = NULL, times,
                                      se.fit = FALSE) {
  if (!inherits(object, "coxph")) {
    stop("Internal error: Fine--Gray CIF prediction requires a retained coxph fit.",
         call. = FALSE)
  }
  if (!is.numeric(times) || length(times) < 1L || anyNA(times) ||
      any(!is.finite(times)) || any(times < 0)) {
    stop("`times` must contain one or more finite non-negative values.",
         call. = FALSE)
  }
  times <- sort(unique(as.numeric(times)))

  prediction_data <- if (is.null(newdata)) {
    object$mfp2_finegray_training_newdata
  } else {
    as.data.frame(newdata, check.names = FALSE)
  }
  if (is.null(prediction_data) || !is.data.frame(prediction_data) ||
      nrow(prediction_data) < 1L) {
    stop(
      "Fine--Gray prediction data are unavailable; refit with the current ",
      "package version or supply `newdata`.",
      call. = FALSE
    )
  }

  base_object <- object
  class(base_object) <- setdiff(class(base_object), "mfp2")
  n <- nrow(prediction_data)
  nt <- length(times)
  cif <- matrix(NA_real_, n, nt)
  cif_se <- if (isTRUE(se.fit)) matrix(NA_real_, n, nt) else NULL
  batch_size <- 256L
  for (first in seq.int(1L, n, by = batch_size)) {
    rows <- seq.int(first, min(n, first + batch_size - 1L))
    batch <- prediction_data[rows, , drop = FALSE]
    nb <- length(rows)
    curve <- tryCatch(
      survival::survfit(base_object, newdata = batch,
                        se.fit = isTRUE(se.fit)),
      error = function(e) NULL
    )
    point <- if (!is.null(curve)) {
      tryCatch(summary(curve, times = times, extend = TRUE),
               error = function(e) NULL)
    } else {
      NULL
    }

    # With baseline strata, survfit concatenates curves in newdata row order;
    # without strata, the curves occupy columns of a time-by-row matrix.
    # Stratum labels are implementation-dependent and need not be newdata row
    # names. Verify the curve counts and time grid instead of comparing labels.
    strata_aligned <- is.null(curve$strata) ||
      (length(curve$strata) == nb && length(point$strata) == nb &&
       !anyNA(point$strata) && all(point$strata == nt))
    expected_times <- if (is.matrix(point$surv)) times else rep(times, nb)
    has_prediction_shape <- function(values) {
      if (is.matrix(values)) {
        identical(dim(values), c(nt, nb))
      } else {
        is.null(dim(values)) && length(values) == nt * nb
      }
    }

    if (!isTRUE(strata_aligned) ||
        !identical(as.numeric(point$time), expected_times) ||
        !has_prediction_shape(point$surv) ||
        (isTRUE(se.fit) && !has_prediction_shape(point$std.err))) {
      # Some survival versions return curve shapes that cannot safely be
      # mapped to the requested rows. Reuse the original one-row procedure.
      for (i in rows) {
        single <- summary(survival::survfit(
          base_object, newdata = prediction_data[i, , drop = FALSE],
          se.fit = isTRUE(se.fit)
        ), times = times, extend = TRUE)
        values <- as.numeric(single$surv)
        if (length(values) != nt) {
          stop("Internal error: Fine--Gray CIF time grid could not be aligned.",
               call. = FALSE)
        }
        cif[i, ] <- 1 - values
        if (isTRUE(se.fit)) {
          errors <- as.numeric(single$std.err)
          if (length(errors) != nt) {
            stop("Internal error: Fine--Gray CIF standard errors are incomplete.",
                 call. = FALSE)
          }
          cif_se[i, ] <- errors
        }
      }
      next
    }

    cif[rows, ] <- 1 - t(matrix(as.numeric(point$surv), nt, nb))
    if (isTRUE(se.fit)) {
      cif_se[rows, ] <- t(matrix(as.numeric(point$std.err), nt, nb))
    }
  }

  time_names <- paste0("time=", format(times, trim = TRUE, scientific = FALSE))
  dimnames(cif) <- list(row.names(prediction_data), time_names)
  if (isTRUE(se.fit)) dimnames(cif_se) <- dimnames(cif)

  fit <- if (nt == 1L) stats::setNames(as.numeric(cif[, 1L]), row.names(prediction_data)) else cif
  if (!isTRUE(se.fit)) return(fit)
  se <- if (nt == 1L) stats::setNames(as.numeric(cif_se[, 1L]), row.names(prediction_data)) else cif_se
  list(fit = fit, se.fit = se, times = times,
       event = object$mfp2_finegray_event)
}

#' Validate Prediction Data
#'
#' Checks the selected prediction columns before transformations or model
#' matrix construction. Callers supply only columns needed for the requested
#' prediction. Errors identify missing values or columns containing infinities.
#'
#' @param newdata Prediction data to check.
#' @param check_missing Logical; whether to reject missing values.
#'   Numeric infinities are always rejected.
#'
#' @return Invisibly `TRUE` when validation succeeds.
#'
#' @keywords internal
#' @noRd
mfp2_validate_prediction_newdata <- function(newdata, check_missing = TRUE) {
  if (!is.logical(check_missing) || length(check_missing) != 1L ||
      is.na(check_missing)) {
    stop("Internal error: `check_missing` must be TRUE or FALSE.",
         call. = FALSE)
  }

  newdata <- as.data.frame(newdata, check.names = FALSE)

  # Every prediction mode requires at least one row. Reject an empty object at
  # the shared validation boundary, before term ranges, equidistant sequences,
  # formula reconstruction, or contrast reference means are calculated. This
  # prevents downstream Inf/-Inf ranges and NaN default reference values.
  if (nrow(newdata) == 0L) {
    stop(
      "! `newdata` must contain at least one row.",
      call. = FALSE
    )
  }

  if (isTRUE(check_missing) && anyNA(newdata)) {
    stop(
      "! newdata must not contain any NA (missing data).\n",
      "i Please remove any missing data before passing newdata to this function.",
      call. = FALSE
    )
  }

  non_finite_columns <- names(newdata)[vapply(
    newdata,
    function(column) {
      if (!is.numeric(column)) {
        return(FALSE)
      }
      bad <- if (isTRUE(check_missing)) {
        !is.finite(column)
      } else {
        is.infinite(column)
      }
      any(bad)
    },
    logical(1L)
  )]

  if (length(non_finite_columns) > 0L) {
    stop(
      "! newdata must contain only finite numeric values.\n",
      "i Infinite value(s) were found in numeric column(s): ",
      paste(non_finite_columns, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Validate the Cox Covariate Reference for the Requested Prediction Type
#'
#' `survival::predict.coxph()` uses `reference` only for relative
#' linear-predictor, risk, and native term predictions. `mfp2` computes its own
#' term and contrast output, while the Cox expected-event/survival path
#' unconditionally uses the sample reference internally. This helper enforces
#' those distinctions and returns a value only when the delegated Cox method can
#' use it.
#'
#' The explicit-supply flag is kept separate from the argument value so that
#' `NULL` means "use the method default", matching `predict.mfpi()`, while a
#' non-NULL value can be rejected on prediction paths where references do not
#' apply.
#'
#' @param type Resolved Cox prediction type.
#' @param cox_reference Candidate Cox reference argument.
#' @param cox_reference_supplied Logical scalar indicating whether the caller
#'   explicitly supplied `cox_reference` to `predict.mfp2()`.
#'
#' @return A single validated reference string for `type = "lp"` or
#'   `type = "risk"`; otherwise `NULL`.
#'
#' @keywords internal
#' @noRd
mfp2_validate_cox_reference <- function(type,
                                        cox_reference,
                                        cox_reference_supplied) {
  if (!is.logical(cox_reference_supplied) ||
      length(cox_reference_supplied) != 1L ||
      is.na(cox_reference_supplied)) {
    stop("Internal error: invalid cox_reference-supply flag.", call. = FALSE)
  }

  if (type %in% c("terms", "contrasts")) {
    if (cox_reference_supplied) {
      stop(
        "'cox_reference' is not used for mfp2 term or contrast predictions. ",
        "Use 'ref' to choose variable-specific contrast reference values.",
        call. = FALSE
      )
    }
    return(NULL)
  }

  if (type %in% c("expected", "survival", "response")) {
    if (cox_reference_supplied) {
      stop(
        "'cox_reference' does not apply to Cox predictions with type = '",
        type,
        "'. The baseline-hazard calculation uses the fitted sample ",
        "reference internally; remove the 'cox_reference' argument.",
        call. = FALSE
      )
    }
    return(NULL)
  }

  # The public NULL default resolves to "zero" for relative Cox predictions.
  # A non-NULL user value must be scalar and is validated below.
  if (!cox_reference_supplied) {
    return("zero")
  }

  match_cox_reference(
    value = cox_reference,
    default = "zero",
    argument = "cox_reference"
  )
}


#' Delegate a Full-Model Prediction to the Stored Cox Fit
#'
#' Removes only the `"mfp2"` class and then calls the native Cox prediction
#' method through `stats::predict()`. The underlying object remains a complete
#' `coxph` fit with its fitted response, transformed design matrix, residuals,
#' means, covariance matrix, strata metadata, and tie method intact.
#'
#' Relative predictions (`"lp"` and `"risk"`) receive the validated covariate
#' reference. Absolute predictions (`"expected"` and `"survival"`) deliberately
#' omit that argument because `predict.coxph()` resets it to `"sample"` before
#' calculating the baseline cumulative hazard. Omitting it also avoids implying
#' that an ignored user value controls absolute predictions.
#'
#' When the fitted model contains an offset, the stored training-offset origin
#' is validated before delegation. Offset centering itself remains the
#' responsibility of `predict.coxph()` and is independent of the covariate
#' reference.
#'
#' @param object Fitted `mfp2` Cox object.
#' @param newdata Optional fully transformed prediction data frame. For absolute
#'   predictions it must also contain the internal `Surv` response column.
#' @param type Resolved Cox prediction type.
#' @param se.fit Logical scalar passed to `predict.coxph()`.
#' @param cox_reference Validated reference string for relative predictions, or
#'   `NULL` for absolute predictions.
#' @param ... Further arguments passed to the native Cox prediction method.
#'
#' @return The unmodified result returned by `predict.coxph()`.
#'
#' @keywords internal
#' @noRd
mfp2_predict_cox_base <- function(object,
                                  newdata = NULL,
                                  type,
                                  se.fit,
                                  cox_reference = NULL,
                                  ...) {
  if (!inherits(object, "mfp2") ||
      !mfp2_family_is_ph(object$family_string)) {
    stop("Internal error: a fitted mfp2 proportional-hazards object is required.", call. = FALSE)
  }

  allowed_types <- if (identical(object$family_string, "finegray")) {
    c("lp", "risk")
  } else {
    c("lp", "risk", "expected", "survival")
  }
  if (!type %in% allowed_types) {
    stop("Internal error: unsupported delegated Cox prediction type.",
         call. = FALSE)
  }

  if (isTRUE(object$has_offset)) {
    mfp2_cox_offset_reference(object)
  }

  obj_base <- object
  class(obj_base) <- setdiff(class(obj_base), "mfp2")

  # Keep the native Cox method authoritative for both covariate-reference and
  # offset origins. In particular, do not add or subtract
  # object$cox_offset_reference here: predict.coxph() reconstructs the fitted
  # model frame when an offset is present, subtracts the mean training offset,
  # and applies the requested covariate reference in the same calculation.
  # Reimplementing either adjustment in mfp2 would double-center one component
  # and can diverge from the exact fitted coxph design.
  if (type %in% c("lp", "risk")) {
    if (is.null(cox_reference)) {
      stop("Internal error: relative Cox prediction lacks a reference.",
           call. = FALSE)
    }

    if (is.null(newdata)) {
      prediction <- stats::predict(
        obj_base,
        type = type,
        se.fit = se.fit,
        reference = cox_reference,
        ...
      )
      if (identical(object$family_string, "finegray")) {
        prediction <- collapse_finegray_training_prediction(
          prediction,
          object$mfp2_finegray_row_map,
          object$nobs
        )
      }
      return(prediction)
    }

    return(stats::predict(
      obj_base,
      newdata = newdata,
      type = type,
      se.fit = se.fit,
      reference = cox_reference,
      ...
    ))
  }

  if (is.null(newdata)) {
    return(stats::predict(
      obj_base,
      type = type,
      se.fit = se.fit,
      ...
    ))
  }

  stats::predict(
    obj_base,
    newdata = newdata,
    type = type,
    se.fit = se.fit,
    ...
  )
}


#' Resolve Conceptual Terms to Raw Prediction Columns
#'
#' Returns the term lookup stored by grouped-term fits. The names are conceptual
#' terms and each value contains the corresponding raw design-matrix columns in
#' fitted term order. The lookup can be restricted to active prediction terms.
#'
#' @param object A fitted \code{mfp2} object.
#' @param terms Optional character vector of conceptual terms to retain.
#'
#' @return Named list mapping conceptual terms to raw input columns.
#'
#' @keywords internal
#' @noRd
prediction_term_to_columns <- function(object, terms = NULL) {
  fitted_terms <- rownames(object$transformations)

  if (is.null(fitted_terms) || length(fitted_terms) == 0L) {
    stop("The fitted object does not contain transformation term names.",
         call. = FALSE)
  }

  lookup <- object$term_to_columns

  if (is.null(lookup)) {
    # Objects created before grouped categorical support did not store an
    # explicit conceptual-term lookup. Formula fits from newer releases may
    # still carry the equivalent formula-level lookup, so prefer it when
    # available. Otherwise a continuous-only legacy object has an identity
    # relationship between conceptual terms and raw predictor columns.
    lookup <- object$formula_term_to_columns

    if (is.null(lookup)) {
      selected_terms <- tryCatch(
        get_selected_variable_names(object),
        error = function(e) character(0L)
      )
      selected_raw_columns <- if (is.null(object$x_original)) {
        character(0L)
      } else {
        colnames(object$x_original)
      }

      # A retained grouped term cannot be reconstructed safely without its
      # stored member-column mapping. Fail explicitly rather than silently
      # treating the conceptual name as a raw design column.
      unresolved_selected <- setdiff(selected_terms, selected_raw_columns)
      if (length(unresolved_selected) > 0L) {
        stop(
          "The fitted object lacks grouped-term mapping metadata for selected term(s): ",
          paste(unresolved_selected, collapse = ", "),
          ". Refit the model with the current mfp2 version.",
          call. = FALSE
        )
      }

      lookup <- stats::setNames(as.list(fitted_terms), fitted_terms)
    }
  }

  if (!is.list(lookup) || is.null(names(lookup)) || anyDuplicated(names(lookup))) {
    stop("The fitted grouped-term lookup is malformed.", call. = FALSE)
  }

  missing_terms <- setdiff(fitted_terms, names(lookup))
  if (length(missing_terms) > 0L) {
    stop(
      "The fitted grouped-term lookup is missing term(s): ",
      paste(missing_terms, collapse = ", "),
      call. = FALSE
    )
  }

  lookup <- lookup[fitted_terms]

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

  if (!is.null(terms)) {
    terms <- unique(as.character(terms))
    unknown <- setdiff(terms, names(lookup))
    if (length(unknown) > 0L) {
      stop(
        "Unknown fitted prediction term(s): ",
        paste(unknown, collapse = ", "),
        call. = FALSE
      )
    }
    lookup <- lookup[terms]
  }

  lookup
}

#' Resolve Transformed Design Columns to Fitted Coefficient Names
#'
#' Maps a vector of transformed design-column names to the exact
#' coefficient names stored by the final formula-based fit. The
#' authoritative map is created once by `fit_mfp()` from the transformed
#' matrix and the fitted coefficient order, so prediction never has to
#' quote, unquote, sanitise, or otherwise guess model column names.
#'
#' @param object An `"mfp2"` model object.
#' @param transformed_columns Character vector of transformed
#'   design-column names.
#'
#' @return Character vector of fitted coefficient names, in the same order
#'   as `transformed_columns`.
#'
#' @keywords internal
#' @noRd
prediction_model_column_names <- function(object, transformed_columns) {
  column_map <- object$transformed_to_model_columns

  if (is.null(column_map)) {
    # Legacy continuous-only objects predate the explicit map. Recreate the
    # complete transformed source-column order from the stored training data,
    # then align it positionally with the fitted coefficient vector. Final
    # glm()/coxph() fits preserve formula-column order, including aliased
    # coefficients, so this recovers exact quoted or non-syntactic names.
    selected_terms <- get_selected_variable_names(object)

    if (length(selected_terms) == 0L) {
      if (length(transformed_columns) > 0L) {
        stop(
          "Cannot map transformed columns for an intercept-only fitted object.",
          call. = FALSE
        )
      }
      return(character(0L))
    }

    if (is.null(object$x_original)) {
      stop(
        "The legacy fitted object lacks training predictors needed to reconstruct coefficient-column metadata.",
        call. = FALSE
      )
    }

    reconstructed <- prepare_newdata_for_predict(
      object,
      object$x_original,
      terms = selected_terms,
      apply_pre = FALSE,
      allow_missing_predictors = FALSE
    )
    source_columns <- colnames(reconstructed)
    if (is.null(source_columns)) {
      source_columns <- character(0L)
    }

    fitted_columns <- names(stats::coef(object))
    if (is.null(fitted_columns)) {
      fitted_columns <- character(0L)
    }
    fitted_columns <- fitted_columns[fitted_columns != "(Intercept)"]

    if (length(source_columns) != length(fitted_columns)) {
      stop(
        "The legacy fitted object does not contain enough metadata to align transformed columns with fitted coefficients.",
        call. = FALSE
      )
    }

    column_map <- stats::setNames(fitted_columns, source_columns)
  }

  if (!is.character(column_map) || is.null(names(column_map)) ||
      anyDuplicated(names(column_map))) {
    stop(
      "The fitted transformed-to-model column mapping is malformed.",
      call. = FALSE
    )
  }

  model_columns <- unname(column_map[transformed_columns])
  if (anyNA(model_columns)) {
    stop(
      "The fitted object lacks coefficient mappings for transformed column(s): ",
      paste(transformed_columns[is.na(model_columns)], collapse = ", "),
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
  lookup <- prediction_term_to_columns(object)
  spike_source <- if (!is.null(object$spike)) {
    object$spike
  } else {
    object$fp_terms[, "spike"]
  }

  expanded <- expand_term_metadata_to_columns(
    term_to_columns = lookup,
    powers = object$fp_powers,
    raw_columns = raw_columns,
    acdx = object$acd,
    zero = object$zero,
    catzero = object$catzero,
    spike = spike_source,
    spike_decision = object$spike_dec,
    acd_parameter = object$acd_parameter
  )

  raw_terms <- unname(expanded$terms)
  shifts <- stats::setNames(
    vapply(
      raw_terms,
      function(term) as.numeric(object$transformations[term, "shift"]),
      numeric(1L)
    ),
    raw_columns
  )

  list(
    terms = expanded$terms,
    powers = expanded$powers,
    shift = shifts,
    acdx = expanded$acdx,
    zero = expanded$zero,
    catzero = expanded$catzero,
    spike = expanded$spike,
    spike_decision = expanded$spike_decision,
    acd_parameter = expanded$acd_parameter
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

#' Is a Prediction Term Explicitly Stored as Binary-Only SAZ?
#'
#' Current SAZ fits store one named decision code per term in
#' `object$spike_dec`. Older serialized fits may not contain this field, and
#' partially reconstructed legacy objects may contain an unnamed vector, omit a
#' term, or retain an `NA` value. None of those states proves that the fitted
#' term is binary-only, so this helper returns `FALSE` rather than allowing a
#' zero-length/indeterminate comparison to reach `if()`.
#'
#' @param object A fitted `mfp2` object.
#' @param v Character scalar giving the conceptual term name.
#'
#' @return A single logical value.
#'
#' @keywords internal
#' @noRd
prediction_term_is_binary_only <- function(object, v) {
  decisions <- object$spike_dec

  if (is.null(decisions) ||
      is.null(names(decisions)) ||
      !v %in% names(decisions)) {
    return(FALSE)
  }

  decision <- decisions[[v]]

  if (!is.numeric(decision) || length(decision) != 1L || is.na(decision)) {
    return(FALSE)
  }

  isTRUE(
    as.integer(decision) == saz_decision_codes[["binary_only"]]
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
#' Variables fitted with zero-handling are exempt because exact-zero values are
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

  # Reuse the same legacy-safe metadata interpretation as the term/contrast
  # display branch so missing `spike_dec` cannot produce inconsistent behavior.
  if (prediction_term_is_binary_only(object, v)) {
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

#' Resolve Formula Labels for Active Prediction Terms
#'
#' Formula fits retain original formula labels for evaluation and conceptual
#' term names for MFP selection. This helper links the two representations so
#' prediction can evaluate only the terms present in the final model.
#'
#' @keywords internal
#' @noRd
prediction_formula_source_name <- function(term_label) {
  factor_source <- formula_factor_source_name(term_label)
  if (!is.null(factor_source)) {
    return(factor_source)
  }

  formula_fp_source_name(term_label)
}


#' Map Formula-Term Labels to Prediction Names for a Fitted `mfp2` Model
#'
#' Returns the map from formula-term labels to the prediction names used
#' internally by [predict.mfp2()]. The map is stored on a formula-based
#' fit so that prediction can look up terms by the same identifiers users
#' see in the printed model, without re-parsing formula labels at
#' prediction time.
#'
#' @param object An `"mfp2"` model object fitted through the formula
#'   interface.
#'
#' @return A named character vector: names are formula-term labels,
#'   values are internal prediction names.
#'
#' @keywords internal
#' @noRd
prediction_formula_term_map <- function(object) {
  map <- object$formula_prediction_term_names
  labels <- attr(object$formula_terms, "term.labels")

  if (!is.null(map)) {
    if (!is.character(map) || is.null(names(map)) || anyDuplicated(names(map))) {
      stop("The fitted formula prediction-term map is malformed.",
           call. = FALSE)
    }
    missing_labels <- setdiff(labels, names(map))
    if (length(missing_labels) > 0L) {
      stop(
        "The fitted formula prediction-term map is missing formula term(s): ",
        paste(missing_labels, collapse = ", "),
        call. = FALSE
      )
    }
    return(map[labels])
  }

  # Backward-compatible fallback for objects fitted before the explicit map was
  # stored. Bare terms and simple factor wrappers can be resolved safely. More
  # complex legacy formula objects retain the previous complete-formula replay.
  map <- stats::setNames(labels, labels)
  lookup <- prediction_term_to_columns(object)
  for (label in labels) {
    source_name <- prediction_formula_source_name(label)
    if (!is.null(source_name) && source_name %in% names(lookup)) {
      map[[label]] <- source_name
    }
  }
  map
}


#' Restrict Stored Formula Terms to Active Conceptual Terms
#'
#' @keywords internal
#' @noRd
prediction_formula_terms <- function(object, terms) {
  terms_object <- object$formula_terms
  term_map <- prediction_formula_term_map(object)
  labels <- attr(terms_object, "term.labels")
  active_labels <- names(term_map)[unname(term_map) %in% terms]
  keep_labels <- active_labels

  # Preserve lower-order terms that use only variables from an active term.
  # model.matrix() uses the term hierarchy when deciding factor contrast coding;
  # retaining these support terms reproduces the fit-time columns for active
  # interactions without reintroducing unrelated predictor dependencies.
  factors <- attr(terms_object, "factors")
  if (!is.null(factors) && length(active_labels) > 0L) {
    for (active_label in active_labels) {
      active_variables <- rownames(factors)[factors[, active_label] != 0]
      support_labels <- colnames(factors)[
        vapply(
          seq_len(ncol(factors)),
          function(index) {
            term_variables <- rownames(factors)[factors[, index] != 0]
            length(term_variables) > 0L &&
              all(term_variables %in% active_variables)
          },
          logical(1L)
        )
      ]
      keep_labels <- union(keep_labels, support_labels)
    }
  }

  drop_indices <- which(!labels %in% keep_labels)
  if (length(drop_indices) > 0L) {
    terms_object <- terms_object[-drop_indices]
  }

  terms_object
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
reconstruct_formula_newdata <- function(object, newdata, terms = NULL) {
  # Matrix-interface fits already receive raw design columns and need no formula
  # evaluation. Formula fits are reconstructed with stored levels and contrasts.
  if (!isTRUE(object$formula_interface)) {
    return(newdata)
  }

  if (is.null(terms)) {
    terms <- get_selected_variable_names(object)
  }
  lookup <- prediction_term_to_columns(object, terms = terms)
  expected <- unlist(lookup, use.names = FALSE)

  # A final intercept-only model has no predictor dependencies. Preserve the
  # requested row count without evaluating any original formula expression.
  if (length(expected) == 0L) {
    return(data.frame(row.names = row.names(as.data.frame(newdata))))
  }

  # Accept callers that already supplied the exact active model-matrix columns.
  if (all(expected %in% colnames(newdata))) {
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
  active_formula_terms <- prediction_formula_terms(object, terms)
  active_factors <- attr(active_formula_terms, "factors")
  active_frame_variables <- if (is.null(active_factors)) {
    character(0L)
  } else {
    rownames(active_factors)[rowSums(active_factors != 0) > 0L]
  }
  active_xlevels <- object$formula_xlevels[
    intersect(names(object$formula_xlevels), active_frame_variables)
  ]

  # Validate only the inputs to the active formula. Checking the resulting
  # model matrix alone would miss inputs masked by an expression, e.g. 1 / Inf.
  source_columns <- intersect(all.vars(active_formula_terms), names(newdata_df))
  mfp2_validate_prediction_newdata(newdata_df[, source_columns, drop = FALSE])

  # model.frame() evaluates only active formula expressions and enforces their
  # fitted factor levels; model.matrix() recreates the fitted contrast coding.
  mf <- stats::model.frame(
    active_formula_terms,
    data = newdata_df,
    na.action = stats::na.pass,
    xlev = active_xlevels
  )
  active_contrasts <- object$formula_contrasts[
    intersect(names(object$formula_contrasts), names(mf))
  ]
  mm <- stats::model.matrix(
    active_formula_terms,
    data = mf,
    contrasts.arg = active_contrasts
  )

  keep_cols <- colnames(mm) != "(Intercept)"
  mm <- mm[, keep_cols, drop = FALSE]

  # Restore the exact column names used by the fitted object after formula
  # sanitisation or internal fp()/fp2() renaming.
  column_map <- object$formula_column_map
  if (!is.null(column_map)) {
    mapped <- unname(column_map[colnames(mm)])
    colnames(mm) <- ifelse(is.na(mapped), colnames(mm), mapped)
  }

  missing <- setdiff(expected, colnames(mm))
  if (length(missing) > 0L) {
    stop(
      "! Could not reconstruct required model-matrix column(s) from newdata: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  mm[, expected, drop = FALSE]
}

#' Reconstruct Formula-Level Survival Strata from Prediction Data
#'
#' Reconstructs the original formula-level survival stratification term from
#' `newdata` for `mfp2` or `mfpi` objects fitted through a formula interface.
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
#'   a survival strata object suitable for passing as the `strata` argument to
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

  # Formula strata are needed only for full-model survival predictions. Keep
  # their input checks local so term predictions can ignore these columns.
  source_columns <- intersect(all.vars(object$formula_strata_terms), names(newdata_df))
  mfp2_validate_prediction_newdata(newdata_df[, source_columns, drop = FALSE])

  mf <- tryCatch(
    stats::model.frame(
      object$formula_strata_terms,
      data = newdata_df,
      na.action = stats::na.pass,
      xlev = object$formula_strata_xlevels
    ),
    error = function(e) {
      stop(
        "! This fitted object uses formula-level survival strata, ",
        "but the strata term could not be reconstructed from `newdata`.\n",
        "i Include the original strata variable(s) in `newdata`; do not ",
        "supply a separate `strata` argument for a formula-stratified fit.\n",
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

#' Retrieve and Validate the Cox Offset Origin Stored at Fit Time
#'
#' `survival::predict.coxph()` handles offsets independently of the covariate
#' `reference`. For both training and new-data predictions it subtracts the mean
#' training offset, so a supplied prediction offset contributes as
#' `new_offset - mean(training_offset)`.
#'
#' The final `fit_mfp()` call stores that mean once as `cox_offset_reference`.
#' This helper validates the metadata before delegation. It intentionally does
#' not alter prediction offsets or combine them with covariate centering; the
#' native Cox method remains the single implementation of that arithmetic.
#'
#' A malformed value normally indicates a Cox object created by an older package
#' version or an object whose internals were modified after fitting. Failing here
#' gives a direct remediation message rather than allowing an unexplained shift
#' in the predicted linear predictor.
#'
#' @param object A fitted `mfp2` Cox object containing
#'   `cox_offset_reference`.
#'
#' @return Unnamed finite numeric scalar equal to the mean offset used in the
#'   final Cox fitting data.
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


#' Rebuild a Formula-Level Offset from Prediction Newdata
#'
#' Formula-interface fits store a response-free terms object containing only the
#' original `offset()` expression and the factor levels required to evaluate it.
#' This helper evaluates that expression in ordinary user-facing `newdata` and
#' returns the raw offset vector expected by the internal fitted Cox or GLM
#' formula.
#'
#' No centering is performed here. For Cox models the raw vector is attached as
#' the collision-safe internal offset column stored by the fitted model, after
#' which `predict.coxph()` subtracts the mean training offset.
#' Keeping reconstruction separate from centering prevents the offset origin
#' from being applied twice and keeps GLM and Cox delegation consistent with
#' their native prediction methods.
#'
#' The helper returns `NULL` for matrix-interface fits and for formula fits that
#' did not contain an offset. Evaluation errors identify the missing original
#' variable and direct callers to the explicit `newoffset` alternative.
#'
#' @param object Fitted `mfp2` object.
#' @param newdata User-facing prediction data before formula reconstruction.
#'
#' @return `NULL` when no formula offset is stored, otherwise a finite numeric
#'   vector with one value per prediction row.
#'
#' @keywords internal
#' @noRd
reconstruct_formula_offset_newdata <- function(object, newdata) {
  if (!isTRUE(object$formula_interface) ||
      is.null(object$formula_offset_terms)) {
    return(NULL)
  }

  newdata_df <- as.data.frame(newdata)

  # Validate the raw numeric variables used by the stored offset expression
  # before model.frame() evaluates transformations such as log(). This gives
  # missing and non-finite inputs a deterministic prediction error and prevents
  # transformation warnings (for example, "NaNs produced") from escaping.
  offset_variables_call <- attr(object$formula_offset_terms, "variables")
  offset_source_variables <- if (is.null(offset_variables_call)) {
    character()
  } else {
    unique(all.vars(offset_variables_call))
  }

  # Do not let model.frame() resolve a missing offset variable from the stored
  # formula environment. That environment can still contain the training-data
  # variable, whose length then differs from the prediction rows and produces a
  # misleading "wrong number of rows" error. Formula-level offsets must be
  # reconstructed from the supplied newdata; otherwise callers must provide
  # newoffset explicitly.
  missing_offset_variables <- setdiff(
    offset_source_variables,
    names(newdata_df)
  )
  if (length(missing_offset_variables) > 0L) {
    stop(
      "! This `mfp2` object was fitted with a formula-level offset, ",
      "but the offset could not be reconstructed from `newdata`.\n",
      "i Include the original offset variable(s) in `newdata`, or supply ",
      "`newoffset = ...` explicitly to `predict()`.\n",
      "i Missing offset variable(s): ",
      paste(missing_offset_variables, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  offset_source_columns <- intersect(offset_source_variables, names(newdata_df))
  if (length(offset_source_columns) > 0L) {
    mfp2_validate_prediction_newdata(
      newdata_df[, offset_source_columns, drop = FALSE]
    )
  }

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

  if (!identical(object$family_string, "multinomial")) out <- as.vector(out)

  out_rows <- if (is.matrix(out)) nrow(out) else length(out)
  if (out_rows != NROW(newdata_df)) {
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
#'   here because the stored Cox formula evaluates `strata()` on the collision-
#'   safe internal column name allocated during the final refit.
#' @param terms Character vector of conceptual terms required for this prediction.
#'   Defaults to the terms selected in the final model.
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
                                        terms = NULL,
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
  prediction_rownames <- row.names(newdata)

  if (is.null(colnames(newdata)) || any(colnames(newdata) == "")) {
    stop("! newdata must have non-empty column names.", call. = FALSE)
  }

  # Resolve the conceptual-term lookup once. Full prediction requires every raw
  # member column; term-specific prediction may include complete term subsets.
  if (is.null(terms)) {
    terms <- get_selected_variable_names(object)
  }
  lookup <- prediction_term_to_columns(object, terms = terms)
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
    newdata <- data.frame(row.names = seq_len(n_newdata))
  } else {
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

    # Keep this guard immediately before shifting and transformation. It also
    # protects internal callers that bypass the public method's earlier check.
    mfp2_validate_prediction_newdata(newdata)

    newdata <- as.matrix(newdata)
    storage.mode(newdata) <- "double"
    metadata <- expand_prediction_metadata(object, raw_columns)

    # Apply the fitted shift only to raw user data. Stored x_original data and
    # internally generated term grids are already on the shifted prediction scale.
    if (apply_pre) {
      newdata <- sweep(newdata, 2L, metadata$shift[raw_columns], "+")
    }

    # Domain validation is independent of who applied the shift. Full-model
    # prediction reaches this helper with raw data and apply_pre = TRUE, whereas
    # term/contrast prediction applies the shift before calling this helper and
    # therefore uses apply_pre = FALSE. In both cases `newdata` is now on the
    # shifted scale consumed by the FP transformation and must be checked before
    # log(), negative-power, fractional-power, or repeated-power evaluation.
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
          zero = x_trans$zero_expanded,
          zero_rows = x_trans$structural_zero_rows
        )
      }
    }

    if (NCOL(newdata) == 0L) {
      newdata <- data.frame(row.names = seq_len(n_newdata))
    } else {
      newdata <- data.frame(newdata, check.names = FALSE)
    }
  }

  # Transformation builders may drop matrix row names. Restore the names from
  # the supplied rows, including for intercept-only designs, before prediction
  # methods use them to label fitted values and standard errors.
  row.names(newdata) <- prediction_rownames

  # Attach high-level strata and offset variables expected by the stored model
  # formula after the predictor design has been reconstructed. The exact helper
  # names are recovered from fit-time metadata so user predictors named
  # `strata_`, `offset_`, or even `..mfp2_*` cannot be overwritten.
  if (object$family_string %in% c("cox", "survreg", "finegray") &&
      !is.null(strata) &&
      !is.null(object$mfp2_internal_names$strata)) {
    strata_name <- mfp2_internal_fit_name(
      object,
      component = "strata"
    )
    newdata[[strata_name]] <- normalize_cox_prediction_strata(
      strata = strata,
      fit_obj = object,
      nobs = nrow(newdata)
    )
  }

  if (!is.null(offset)) {
    if (!is.numeric(offset) || anyNA(offset) || any(!is.finite(offset))) {
      stop("! offset must be a finite numeric vector.", call. = FALSE)
    }
    if (length(offset) != nrow(newdata)) {
      stop("! offset must have one value per observation.", call. = FALSE)
    }
    # Formula-based fits store the exact offset helper name. Matrix-based
    # fastglm fits consume the offset separately and need no helper column.
    if (!inherits(object, "fastglm")) {
      offset_name <- mfp2_internal_fit_name(
        object,
        component = "offset"
      )
      newdata[[offset_name]] <- offset
    }
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
#' covariance should be included when no reference value is supplied. Shared
#' ordinal terms set this to `FALSE`; cutoff-specific terms set it to `TRUE`.
#' @param intercept_name Exact fitted coefficient name of the intercept when
#'   `include_intercept = TRUE`. Ordinal fits use a different name per cutoff.
#'
#' @details
#' See pages 91-92 and following in the book by Royston and Sauerbrei 2008
#' for the formulas and mathematical details.
#'
#' Only the covariance block required by `X` (and, when requested, the
#' intercept) is validated. Non-finite entries and materially negative
#' quadratic-form variances cause an error. A negative variance whose magnitude
#' is within a scale-aware floating-point tolerance is treated as numerical
#' roundoff and replaced by zero before the square root is taken.
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
                                     include_intercept = TRUE,
                                     intercept_name = "(Intercept)") {

  # vcov.orm() defaults to the middle cutoff only. Request every intercept
  # when a cutoff is included, or just the common slope block otherwise.
  vcovx <- if (identical(model$family_string, "ordinal")) {
    if (isTRUE(include_intercept) && is.null(xref)) {
      mfp2_ordinal_vcov_all(model)
    } else {
      mfp2_ordinal_vcov_slopes(model)
    }
  } else {
    vcov(object = model)
  }

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
  if (isTRUE(include_intercept) &&
      mfp2_family_has_intercept(model$family_string) && is.null(xref)) {
    intercept_ind <- match(intercept_name, colnames(vcovx))
    if (is.na(intercept_ind)) {
      stop(
        "Cannot compute SE with intercept; covariance matrix lacks `",
        intercept_name, "`.",
        call. = FALSE
      )
    }
    X <- cbind(1, X)
    ind <- c(intercept_ind, ind)
  }

  # Select before validation: an undefined covariance for an unrelated, unused
  # coefficient must not prevent an otherwise estimable prediction. All custom
  # prediction paths share the same finite-value and negative-variance checks.
  vcovx <- vcovx[ind, ind, drop = FALSE]
  mfp2_prediction_standard_error(X, vcovx)
}
