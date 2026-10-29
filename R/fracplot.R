#' Plot fitted term contributions and contrasts from an mfp2 model
#'
#' Shows how selected variables contribute to the fitted model, with confidence
#' intervals. Term plots can also include component-plus-residual points to
#' help assess the fitted relationship.
#'
#' @param x A fitted `mfp2` model.
#' @param terms Names of the variables to plot, for example,
#'   `c("age", "cavol")`. By default (`NULL`), requests all variables kept in
#'   the final model. See **Choosing terms**.
#' @param partial_only Whether to show only the fitted effect and confidence
#'   interval. Defaults to `FALSE`. See **Component-plus-residual plots**.
#' @param type Either `"terms"` (default) or `"contrasts"`.
#'   See **Choosing the plot type**.
#' @param ref Named list of reference values for contrast plots, for example,
#'   `list(age = 50)`. See **Contrasts and reference values**.
#' @param terms_seq Either `"data"` or `"equidistant"`.
#'   See **Evaluation values** for their meaning and defaults.
#' @param alpha Significance level for confidence intervals. The default,
#'   `0.05`, gives 95 percent intervals. See **Confidence intervals**.
#' @param show_titles Whether to show titles describing the fitted
#'   transformation or variable representation. Defaults to `TRUE`.
#' @param show_grid Whether to show the major and minor panel grid lines of
#'   [ggplot2::theme_bw()]. Defaults to `FALSE`, which draws plots without
#'   grid lines.
#' @param color_points Colour of component-plus-residual points, for example,
#'   `"grey60"` or `"#AAAAAA"`. All colour arguments accept ggplot2 colour
#'   names or hexadecimal codes. See [ggplot2::aes_colour_fill_alpha].
#' @param color_line Colour of fitted curves and fitted points for binary
#'   predictors and factors, for example, `"black"` or `"blue"`. Also used
#'   for the zero endpoint of a term fitted with `zero = TRUE` alone.
#' @param color_line_spike Colour of the separate zero endpoint for a
#'   `catzero` or spike-at-zero term with a continuous component,
#'   for example, `"red"`.
#' @param color_fill Fill colour of confidence bands, for example,
#'   `"lightblue"` or `"#56B4E9"`.
#' @param shape Shape of component-plus-residual points. Common ggplot2
#'   shape codes include `1` for open circles and `16` for solid circles.
#'   The default is `16`. See [ggplot2::scale_shape()].
#' @param size_points Numeric size of component-plus-residual points,
#'   in millimetres, for example, `1` or `2`. The default is `1`.
#'   See [ggplot2::geom_point()].
#' @param size_points_spike Numeric size of fitted points for binary
#'   predictors, factors, and zero endpoints.
#'   For example, use `3` for larger points. The default is `2`.
#'   See [ggplot2::geom_point()].
#' @param linetype Line style of fitted curves, for example, `"solid"`,
#'   `"dashed"`, or `"dotted"`. The default is `"solid"`.
#'   See [ggplot2::scale_linetype()].
#' @param linewidth Numeric width of fitted curves, for example, `0.5`
#'   or `1`. Larger values produce thicker lines. The default is `1`.
#'   See [ggplot2::geom_line()].
#' @param alpha_fill Confidence-band transparency, from `0` (transparent)
#'   to `1` (opaque).
#' @param ... Additional arguments are not used. Supplying `newdata` produces
#'   an error; other arguments are ignored with a warning.
#'
#' @section Choosing the plot type:
#' Use `type` to choose what the fitted curve or fitted points represent:
#'
#' - `"terms"` (default) shows a variable's fitted contribution, using the
#'   relationship selected during fitting.
#' - `"contrasts"` shows the difference from a reference value. See
#'   **Contrasts and reference values**.
#'
#' Fitted contributions use the following model scales:
#'
#' - Gaussian with an identity link: contributions to the mean response,
#'   in the same units as the response used for fitting. For example, fitting
#'   log-PSA gives contributions to mean log-PSA.
#' - Binomial with a logit link: log-odds.
#' - Poisson and negative-binomial with a log link: log mean.
#' - Negative-binomial with a square-root link: square root of the mean;
#'   with an identity link: the mean.
#' - Other GLM and GEE models: the scale defined by the fitted link function.
#' - Cox: log relative hazard.
#' - Fine--Gray: log relative subdistribution hazard.
#' - `survreg`: the fitted location scale, such as log time for Weibull models.
#'
#' Term plots include the intercept for GLM, negative-binomial, GEE, and
#' `survreg` models. Cox and Fine--Gray terms have no intercept. Offsets are
#' excluded from fitted term contributions. See **Multinomial and ordinal
#' models** for their scales, panels, and intercepts.
#'
#' @section Choosing terms:
#' Use `terms` to select variables kept in the final model. Use their original
#' names: `"age"` rather than `"fp(age)"`, and `"race"` rather than individual
#' contrast-column names such as `"raceB"` and `"raceC"`.
#'
#' For example, for a model containing age, cancer volume, and the factor race:
#'
#' ```r
#' # Plot age only.
#' age_plots <- plot(fit, terms = "age")
#' age_plots[["age"]]
#'
#' # Create separate plots and combine them using patchwork, if installed.
#' plots <- plot(fit, terms = c("age", "cavol"))
#' patchwork::wrap_plots(plots, ncol = 2)
#'
#' # Show all fitted race levels together.
#' race_plots <- plot(fit, terms = "race")
#' race_plots[["race"]]
#' ```
#'
#' @section Evaluation values:
#' Use `terms_seq` to choose where a numeric variable's fitted effect is
#' calculated:
#'
#' - `"data"` uses the observed predictor values in the data used to fit
#'   the model.
#' - `"equidistant"` uses equally spaced values between the smallest and
#'   largest observed predictor values.
#'
#' The default is `"data"` for term plots and `"equidistant"` for contrast
#' plots. Both settings use the same fitted model.
#'
#' ```r
#' # Evaluate the curve at observed ages.
#' plots <- plot(fit, terms = "age", terms_seq = "data")
#'
#' # Evaluate it at equally spaced ages over the same range.
#' curves <- plot(
#'   fit, terms = "age", terms_seq = "equidistant", partial_only = TRUE
#' )
#' curves[["age"]]
#' ```
#'
#' See **Binary and categorical predictors** for discrete variables.
#' To obtain estimates at specific new values, use [predict.mfp2()] with
#' `newdata` and `terms_seq = "data"`, then plot its output. An example is
#' included below.
#'
#' @section Contrasts and reference values:
#' A contrast is:
#'
#' `contribution at the plotted value - contribution at the reference`.
#'
#' For example, compare each age with age 50:
#'
#' ```r
#' plots <- plot(
#'   fit, terms = "age", type = "contrasts", ref = list(age = 50)
#' )
#' plots[["age"]]
#' ```
#'
#' The contrast is zero at age 50. Above zero, the fitted contribution is
#' larger than at age 50; below zero, it is smaller. Intercepts cancel from
#' the comparison.
#'
#' The interpretation follows the model scale: a Gaussian identity-link
#' contrast is a difference in mean response; a binomial logit contrast is
#' a log odds ratio; and a log-link count-model contrast is a log ratio of
#' mean counts. Cox and Fine--Gray contrasts are log hazard ratios and log
#' subdistribution hazard ratios, respectively. These comparisons hold the
#' other predictors and any offset fixed and, where applicable, use the same
#' stratum.
#'
#' Supply references on the original variable scale. The fitted shift,
#' transformation, and centering are applied automatically. For example, if
#' age was fitted using `log(age + 10)`, supply `ref = list(age = 50)`, not
#' `list(age = 60)` or `list(age = log(60))`.
#'
#' For a formula factor, use a category label available during fitting.
#' For example, if race has categories A, B, and C, compare each with race A:
#'
#' ```r
#' plots <- plot(
#'   fit, terms = "race", type = "contrasts", ref = list(race = "A")
#' )
#' plots[["race"]]
#' ```
#'
#' If a reference is omitted, the function uses:
#'
#' - the smaller value for a numeric predictor with two distinct values;
#' - the observed mean in the fitting data for other numeric predictors;
#' - the first fitted level for a formula factor.
#'
#' See **Contrasts and reference values** in [predict.mfp2()] for further
#' details.
#'
#' @section Component-plus-residual plots:
#' By default, term plots add one component-plus-residual point for each
#' observation used to fit the model. Each point has:
#'
#' - a horizontal position given by the observed predictor value;
#' - a height equal to the fitted contribution at that value plus a residual
#'   from the complete model.
#'
#' These points leave out the fitted contributions of the other predictors,
#' while keeping the plotted predictor's contribution and the residual from
#' the complete model. This helps assess the fitted relationship while
#' accounting for the other predictors. Residuals are calculated from the
#' complete fitted model, including any offsets.
#'
#' The type of residual added depends on the model:
#'
#' - GLM, negative-binomial, and `survreg`: deviance residuals.
#' - GEE: deviance residuals based on fitted marginal means and prior weights.
#' - Cox: martingale residuals.
#' - Fine--Gray: martingale residuals from the weighted Cox fit, summed over
#'   each subject's expanded rows.
#' - Multinomial: binary deviance residuals comparing each non-reference
#'   category with all other categories. This differs from the comparison
#'   with the reference category shown by the fitted curve.
#' - Ordinal: deviance residuals for the cutoff shown in each panel.
#'
#' For non-Gaussian models, these points provide a diagnostic display rather
#' than an exact adjustment of the observed response for other predictors.
#' Points always use observed predictor values, even when the fitted curve
#' uses `terms_seq = "equidistant"`.
#'
#' Set `partial_only = TRUE` to hide the component-plus-residual points:
#'
#' ```r
#' curves <- plot(fit, terms = "age", partial_only = TRUE)
#' curves[["age"]]
#' ```
#'
#' Contrast plots always omit component-plus-residual points.
#'
#' @section Confidence intervals:
#' Confidence bands and vertical intervals describe uncertainty in the fitted
#' contribution or contrast, rather than the spread of component-plus-residual
#' points. They use the estimated uncertainty in the final model's coefficients
#' and apply to each plotted value individually. They are not prediction
#' intervals for future outcomes.
#'
#' The intervals treat the selected model as fixed. They do not account for
#' uncertainty from selecting variables, FP powers, or spike-at-zero
#' components. For GEE models, they use the covariance selected for the fit.
#'
#' @section Multinomial and ordinal models:
#' Multinomial models return one plot per variable, with a panel for each
#' non-reference category. A panel labelled `"B vs A"` shows the fitted
#' contribution to the log-odds of category B versus reference category A,
#' including that comparison's intercept for term plots. The fitted
#' transformation is shared across panels, but the coefficients differ.
#' Contrasts are shown separately for each category comparison.
#'
#' Ordinal term plots have one panel per response cutoff. For ordered
#' categories low, medium, and high, the panels compare medium/high with low,
#' and high with low/medium. Each curve and confidence interval includes its
#' cutoff intercept. The default logistic link uses cumulative log-odds;
#' other ordinal links use their corresponding link scale.
#'
#' Ordinal contrast plots show the shared predictor effect without cutoff
#' intercepts and therefore use a single panel.
#'
#' @section Binary and categorical predictors:
#' Binary predictors and formula factors are displayed using fitted points
#' and vertical confidence intervals rather than a continuous curve.
#' Numeric binary predictors use their observed values on the horizontal
#' axis; formula factors use their fitted category labels.
#'
#' For numeric predictors with only two values in multinomial or ordinal
#' models, use `terms_seq = "data"` to avoid displaying intermediate values.
#'
#' @section Predictors containing zeros:
#' For a term fitted with `zero = TRUE`, the continuous function is defined
#' for positive values. Its curve and confidence band are drawn only above
#' zero, with the fitted value at exactly zero shown as a separate point.
#' Centering can make that contribution nonzero. This point comes from the
#' continuous term and does not represent a separately estimated zero effect.
#'
#' A `catzero` term also includes a separately estimated zero indicator.
#' A spike-at-zero term follows the representation selected in the final model:
#'
#' - a positive-value function with a zero indicator;
#' - the positive-value function alone;
#' - the zero indicator alone.
#'
#' When a continuous component is retained, the positive-value curve is not
#' connected to the zero endpoint. For `catzero` and spike-at-zero terms,
#' that endpoint uses `color_line_spike`; for `zero = TRUE` alone, it uses
#' `color_line`.
#'
#' A binary-only spike term uses labels such as `"exposure = 0"` and
#' `"exposure > 0"`, with fitted points in `color_line`.
#'
#' @return A named list of `ggplot2` objects, one per plotted variable.
#'   Individual plots can be printed by name or position and modified using
#'   `ggplot2`.
#'
#' @examples
#' ## Gaussian models: terms and contrasts
#' data("prostate")
#' fit_gaussian <- mfp2(
#'   lpsa ~ fp(age) + fp(cavol) + fp(weight) + svi,
#'   data = prostate,
#'   family = "gaussian",
#'   keep = c("age", "cavol"),
#'   verbose = FALSE
#' )
#'
#' # Term contributions with component-plus-residual points.
#' plots <- plot(fit_gaussian, terms = c("age", "cavol"))
#' plots[["age"]]
#'
#' # The response is lpsa, so the vertical scale remains log-PSA.
#' # Show only the fitted contribution and its confidence band.
#' curves <- plot(
#'   fit_gaussian,
#'   terms = "age",
#'   partial_only = TRUE,
#'   terms_seq = "equidistant"
#' )
#' curves[["age"]]
#'
#' # Compare each plotted age with age 50.
#' # Supply the reference before any shift or transformation.
#' contrasts <- plot(
#'   fit_gaussian,
#'   terms = "age",
#'   type = "contrasts",
#'   ref = list(age = 50)
#' )
#' contrasts[["age"]]
#'
#' # Modify an individual plot.
#' curves[["age"]] + ggplot2::labs(
#'   x = "Age (years)",
#'   y = "Age contribution to log-PSA"
#' )
#'
#' # To evaluate specific new ages, use predict() and plot its output.
#' age_predictions <- predict(
#'   fit_gaussian,
#'   newdata = data.frame(age = c(40, 50, 60)),
#'   type = "contrasts",
#'   terms = "age",
#'   terms_seq = "data",
#'   ref = list(age = 50)
#' )[["age"]]
#'
#' ggplot2::ggplot(
#'   age_predictions,
#'   ggplot2::aes(x = variable, y = value)
#' ) +
#'   ggplot2::geom_point() +
#'   ggplot2::geom_errorbar(
#'     ggplot2::aes(ymin = lower, ymax = upper),
#'     width = 1
#'   ) +
#'   ggplot2::labs(y = "Difference from age 50")
#'
#' \dontrun{
#' ## Gaussian models: categorical predictors
#' fit_factor <- mfp2(
#'   Sepal.Length ~ Species,
#'   data = iris,
#'   family = "gaussian",
#'   keep = "Species",
#'   verbose = FALSE
#' )
#'
#' # Use the factor name, not individual dummy-column names.
#' # Each species has a fitted point and vertical confidence interval.
#' factor_plots <- plot(
#'   fit_factor,
#'   terms = "Species",
#'   partial_only = TRUE
#' )
#' factor_plots[["Species"]]
#'
#' factor_contrasts <- plot(
#'   fit_factor,
#'   terms = "Species",
#'   type = "contrasts",
#'   ref = list(Species = "setosa")
#' )
#' factor_contrasts[["Species"]]
#'
#' ## Gaussian models: a predictor containing zeros
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
#' # The display follows the selected spike representation.
#' # If a continuous component is retained, zero is shown separately
#' # from the positive-value curve.
#' spike_plots <- plot(
#'   fit_spike,
#'   terms = "exposure",
#'   partial_only = TRUE,
#'   color_line_spike = "blue"
#' )
#' spike_plots[["exposure"]]
#'
#' spike_contrasts <- plot(
#'   fit_spike,
#'   terms = "exposure",
#'   type = "contrasts",
#'   ref = list(exposure = 0)
#' )
#' spike_contrasts[["exposure"]]
#'
#' ## Binomial models: log-odds contributions
#' set.seed(2)
#' d_binary <- data.frame(age = runif(150, 30, 75))
#' d_binary$event <- rbinom(
#'   150,
#'   size = 1,
#'   prob = plogis(0.06 * (d_binary$age - 50))
#' )
#'
#' fit_binomial <- mfp2(
#'   event ~ fp(age, df = 1, select = 1),
#'   data = d_binary,
#'   family = "binomial",
#'   verbose = FALSE
#' )
#'
#' # The curve is on the log-odds scale, not the probability scale.
#' binomial_plots <- plot(fit_binomial, partial_only = TRUE)
#' binomial_plots[["age"]]
#'
#' ## Poisson models: term plots with a formula offset
#' set.seed(3)
#' d_count <- data.frame(
#'   age = runif(100, 30, 75),
#'   exposure = runif(100, 0.5, 3)
#' )
#' d_count$count <- rpois(
#'   100,
#'   lambda = with(d_count, exposure * exp(0.03 * (age - 50)))
#' )
#'
#' fit_poisson <- mfp2(
#'   count ~ fp(age, select = 1) + offset(log(exposure)),
#'   data = d_count,
#'   family = "poisson",
#'   verbose = FALSE
#' )
#'
#' # The curve includes the intercept and age contribution on the log scale.
#' # It does not include the exposure offset or show expected counts.
#' poisson_plots <- plot(fit_poisson, partial_only = TRUE)
#' poisson_plots[["age"]]
#'
#' ## Negative-binomial models: overdispersed counts
#' if (requireNamespace("fastglm", quietly = TRUE)) {
#'   d_count$overdispersed <- rnbinom(
#'     nrow(d_count),
#'     mu = with(d_count, exposure * exp(0.03 * (age - 50))),
#'     size = 2
#'   )
#'
#'   fit_negbin <- mfp2(
#'     overdispersed ~ fp(age, df = 1, select = 1) +
#'       offset(log(exposure)),
#'     data = d_count,
#'     family = "negbin",
#'     verbose = FALSE
#'   )
#'
#'   # Contrasts are differences on the log-mean scale.
#'   negbin_plots <- plot(
#'     fit_negbin,
#'     type = "contrasts",
#'     ref = list(age = 50)
#'   )
#'   print(negbin_plots[["age"]])
#' }
#'
#' ## GEE models: component-plus-residual plots for grouped observations
#' # Illustrate 25 groups with four observations each.
#' d_count$id <- rep(seq_len(25), each = 4)
#' fit_gee <- mfp2(
#'   count ~ fp(age, df = 1, select = 1) + offset(log(exposure)),
#'   data = d_count,
#'   family = gee_family(poisson(), corstr = "independence"),
#'   id = id,
#'   verbose = FALSE
#' )
#'
#' # Confidence intervals use the covariance selected for the GEE fit.
#' # Points add marginal deviance residuals to the age contribution.
#' gee_plots <- plot(fit_gee, terms = "age")
#' gee_plots[["age"]]
#'
#' ## Multinomial models: category-specific panels
#' set.seed(5)
#' d_multi <- data.frame(age = runif(200, 30, 75), wide = rbinom(200, 1, 0.5))
#' eta_b <- with(d_multi, -0.5 + 0.04 * (age - 50) + 0.6 * wide)
#' eta_c <- with(d_multi, -0.3 + 0.07 * (age - 50) - 0.4 * wide)
#' p_multi <- cbind(1, exp(eta_b), exp(eta_c))
#' p_multi <- p_multi / rowSums(p_multi)
#' d_multi$type <- factor(
#'   apply(p_multi, 1, function(p) sample(c("A", "B", "C"), 1, prob = p))
#' )
#'
#' fit_multinomial <- mfp2(
#'   type ~ fp(age, df = 1, select = 1) + wide,
#'   data = d_multi,
#'   family = multinomial_family(reference = "A"),
#'   keep = "wide",
#'   verbose = FALSE
#' )
#'
#' # Each panel compares a non-reference category with category A.
#' multinomial_plots <- plot(
#'   fit_multinomial,
#'   terms = "age",
#'   partial_only = TRUE
#' )
#' multinomial_plots[["age"]]
#'
#' # For a numeric binary predictor, use observed values explicitly.
#' # Compare wide = 1 with wide = 0 within each category comparison.
#' binary_contrasts <- plot(
#'   fit_multinomial,
#'   terms = "wide",
#'   type = "contrasts",
#'   terms_seq = "data",
#'   ref = list(wide = 0)
#' )
#' binary_contrasts[["wide"]]
#'
#' ## Ordinal models: cutoff panels and shared contrasts
#' if (requireNamespace("rms", quietly = TRUE) &&
#'     utils::packageVersion("rms") >= "7.0-0") {
#'   set.seed(4)
#'   d_ordinal <- data.frame(age = runif(150, 30, 75))
#'   d_ordinal$severity <- cut(
#'     0.04 * (d_ordinal$age - 50) + rlogis(150),
#'     breaks = c(-Inf, -0.5, 0.5, Inf),
#'     labels = c("low", "medium", "high"),
#'     ordered_result = TRUE
#'   )
#'
#'   fit_ordinal <- mfp2(
#'     severity ~ fp(age, df = 1, select = 1),
#'     data = d_ordinal,
#'     family = ordinal_family(),
#'     verbose = FALSE
#'   )
#'
#'   # Panels compare medium/high with low, and high with low/medium.
#'   # Each panel includes its own cutoff intercept.
#'   ordinal_plots <- plot(fit_ordinal, partial_only = TRUE)
#'   print(ordinal_plots[["age"]])
#'
#'   # Contrasts omit cutoff intercepts and show one shared age effect.
#'   ordinal_contrasts <- plot(
#'     fit_ordinal,
#'     terms = "age",
#'     type = "contrasts",
#'     ref = list(age = 50)
#'   )
#'   print(ordinal_contrasts[["age"]])
#' }
#'
#' ## Cox models: contributions to log relative hazard
#' data("gbsg")
#' fit_cox <- mfp2(
#'   survival::Surv(rectime, censrec) ~
#'     fp(age, df = 1, select = 1) +
#'     fp(nodes, df = 1, select = 1) + strata(grade),
#'   data = gbsg,
#'   family = "cox",
#'   verbose = FALSE
#' )
#'
#' # Points add martingale residuals to the fitted age contribution.
#' # The strata term is not plotted as a predictor effect.
#' cox_plots <- plot(fit_cox, terms = "age")
#' cox_plots[["age"]]
#'
#' # Differences from age 50 are on the log relative-hazard scale.
#' # These are not survival probabilities.
#' cox_contrasts <- plot(
#'   fit_cox,
#'   terms = "age",
#'   type = "contrasts",
#'   ref = list(age = 50)
#' )
#' cox_contrasts[["age"]]
#'
#' ## Parametric survival models: fitted location contributions
#' fit_survreg <- mfp2(
#'   survival::Surv(rectime, censrec) ~
#'     fp(age, df = 1, select = 1) +
#'     fp(nodes, df = 1, select = 1),
#'   data = gbsg,
#'   family = survreg_family(dist = "weibull"),
#'   verbose = FALSE
#' )
#'
#' # For Weibull models, this contribution is on the log-time scale.
#' # It is not a survival probability or a predicted time quantile.
#' survreg_plots <- plot(
#'   fit_survreg,
#'   terms = "age",
#'   partial_only = TRUE
#' )
#' survreg_plots[["age"]]
#'
#' ## Fine--Gray models: log subdistribution-hazard contributions
#' set.seed(5)
#' d_fg <- data.frame(age = runif(180, 30, 75))
#' relapse_time <- rexp(
#'   180,
#'   rate = exp(0.04 * (d_fg$age - 50)) / 1000
#' )
#' death_time <- rexp(180, rate = 1 / 1800)
#' censor_time <- runif(180, 300, 1800)
#' d_fg$time <- pmin(relapse_time, death_time, censor_time)
#' d_fg$event <- factor(
#'   ifelse(
#'     d_fg$time == relapse_time, "relapse",
#'     ifelse(d_fg$time == death_time, "death", "censored")
#'   ),
#'   levels = c("censored", "relapse", "death")
#' )
#'
#' fit_finegray <- mfp2(
#'   survival::Surv(time, event) ~ fp(age, df = 1, select = 1),
#'   data = d_fg,
#'   family = finegray_family(etype = "relapse"),
#'   verbose = FALSE
#' )
#'
#' # This shows the age contribution to log subdistribution hazard.
#' # Cumulative-incidence probabilities require predict(type = "response",
#' # times = ...), rather than plot().
#' finegray_plots <- plot(
#'   fit_finegray,
#'   terms = "age",
#'   partial_only = TRUE
#' )
#' finegray_plots[["age"]]
#' }
#'
#' @seealso [predict.mfp2()], [mfp2()]
#' @importFrom ggplot2 .data
#' @method plot mfp2
#' @export
plot.mfp2 <- function(x,
                      terms = NULL,
                      partial_only = FALSE,
                      type = c("terms", "contrasts"),
                      ref = NULL,
                      terms_seq = c("data", "equidistant"),
                      alpha = 0.05,
                      show_titles = TRUE,
                      show_grid = FALSE,
                      color_points = "#AAAAAA",
                      color_line = "black",
                      color_line_spike = "black",
                      color_fill = "#000000",
                      shape = 16,
                      size_points = 1,
                      size_points_spike = 2,
                      linetype = "solid",
                      linewidth = 1,
                      alpha_fill = 0.1,
                      ...) {
  terms_seq_missing <- missing(terms_seq)

  # `...` exists only for the plot() generic. Do not let a supplied argument
  # disappear silently: plots are always drawn for the fitted data.
  dot_names <- names(match.call(expand.dots = FALSE)[["..."]])
  n_dots <- ...length()
  if (n_dots > 0L) {
    if (is.null(dot_names)) dot_names <- rep("", n_dots)
    if ("newdata" %in% dot_names) {
      stop(
        "`plot()` for mfp2 models draws the fitted data and has no ",
        "`newdata` argument. Use `predict(x, type = \"terms\", newdata = ...)` ",
        "or `type = \"contrasts\"` to evaluate terms at new values.",
        call. = FALSE
      )
    }
    warning(
      "Ignoring argument(s) not used by `plot()` for mfp2 models: ",
      paste(ifelse(nzchar(dot_names), dot_names, "<unnamed>"), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  plot_mfp2_impl(
    model = x,
    terms = terms,
    partial_only = partial_only,
    type = type,
    ref = ref,
    terms_seq = terms_seq,
    alpha = alpha,
    show_titles = show_titles,
    show_grid = show_grid,
    color_points = color_points,
    color_line = color_line,
    color_line_spike = color_line_spike,
    color_fill = color_fill,
    shape = shape,
    size_points = size_points,
    size_points_spike = size_points_spike,
    linetype = linetype,
    linewidth = linewidth,
    alpha_fill = alpha_fill,
    terms_seq_missing = terms_seq_missing
  )
}

#' Choose a Scale Label for a Fitted-Only Term Plot
#'
#' Uses the fitted link or survival distribution. Component-plus-residual and
#' contrast plots keep their existing labels and do not call this helper.
#'
#' @param model A fitted `mfp2` model.
#' @return A character value for the vertical axis.
#' @keywords internal
#' @noRd
mfp2_plot_partial_ylabel <- function(model) {
  family <- model$family_string
  if (identical(family, "cox")) return("Contribution to log hazard")
  if (identical(family, "finegray")) return("Contribution to log subhazard")
  if (identical(family, "multinomial")) return("Contribution to log odds")
  if (identical(family, "negbin") &&
      !(is.list(model$family) && is.character(model$family$link))) {
    # Negative-binomial fits store their link in the fastglm family; a fit
    # without that metadata can only have used the original log link.
    return("Contribution to log mean")
  }

  if (identical(family, "survreg")) {
    distribution <- model$dist
    if (is.null(distribution) && is.list(model$family)) {
      distribution <- model$family$dist
    }
    log_time <- c("weibull", "exponential", "lognormal", "loggaussian",
                  "loglogistic")
    if (is.character(distribution) && length(distribution) == 1L &&
        !is.na(distribution) && tolower(distribution) %in% log_time) {
      return("Contribution to log time")
    }
    # Custom or other distributions need not model log time.
    return("Contribution to fitted location")
  }

  response_family <- model$family
  if (identical(family, "gee") && !inherits(response_family, "family")) {
    response_family <- if (is.list(response_family)) {
      response_family$response_family
    } else {
      NULL
    }
    if (is.null(response_family) && is.list(model$mfp2_family)) {
      response_family <- model$mfp2_family$response_family
    }
  }
  link <- if (is.list(response_family)) response_family$link else NULL
  if (identical(family, "ordinal")) {
    if (!is.null(model$mfp2_ordinal_link)) link <- model$mfp2_ordinal_link
    if (is.null(link) && is.list(model$mfp2_family)) {
      link <- model$mfp2_family$link
    }
    if (identical(link, "logistic") || identical(link, "logit")) {
      return("Contribution to cumulative log odds")
    }
  }
  if (!is.character(link) || length(link) != 1L || is.na(link) ||
      !nzchar(link)) return("Fitted contribution")

  switch(link,
    identity = "Contribution to mean response",
    log = "Contribution to log mean",
    logit = "Contribution to log odds",
    inverse = "Contribution to inverse mean",
    "1/mu^2" = "Contribution to 1/mean\u00b2",
    sqrt = "Contribution to sqrt(mean)",
    probit = "Probit contribution",
    cloglog = "Cloglog contribution",
    loglog = "Log-log contribution",
    cauchit = "Cauchit contribution",
    paste0("Contribution (", link, " scale)")
  )
}


#' Is a Model Term Binary for Plotting Purposes?
#'
#' Identifies fitted terms that should be displayed as two-level effects
#' rather than as continuous curves. A formula factor with exactly two levels
#' is treated as binary; a singleton numeric term is treated as binary when
#' its fitted raw design values contain only two distinct non-missing values.
#'
#' @param model An `"mfp2"` model object.
#' @param term Character scalar naming a model term.
#'
#' @return `TRUE` when the term should be displayed as a two-level effect,
#'   `FALSE` otherwise.
#'
#' @keywords internal
#' @noRd
mfp2_plot_term_is_binary <- function(model, term) {
  factor_info <- if (!is.null(model$formula_factor_info)) {
    model$formula_factor_info[[term]]
  } else {
    NULL
  }

  if (!is.null(factor_info) && length(factor_info$levels) == 2L) {
    return(TRUE)
  }

  term_columns <- if (!is.null(model$term_to_columns) &&
                      term %in% names(model$term_to_columns)) {
    model$term_to_columns[[term]]
  } else {
    term
  }

  if (length(term_columns) != 1L ||
      is.null(model$x_original) ||
      is.null(colnames(model$x_original)) ||
      !term_columns %in% colnames(model$x_original)) {
    return(FALSE)
  }

  values <- model$x_original[, term_columns, drop = TRUE]
  values <- values[!is.na(values)]

  length(unique(values)) == 2L
}

#' Is a Model Term a Formula-Factor Effect for Plotting?
#'
#' Reports whether a term was created from a formula factor and should
#' therefore be plotted with discrete-level geoms regardless of the number of
#' fitted levels. Kept separate from binary detection so that factors with
#' three or more levels are never routed through the continuous line/ribbon
#' plotting path.
#'
#' @param model An `"mfp2"` model object.
#' @param term Character scalar naming a model term.
#'
#' @return `TRUE` when the term originated from a formula factor, `FALSE`
#'   otherwise.
#'
#' @keywords internal
#' @noRd
mfp2_plot_term_is_factor <- function(model, term) {
  !is.null(model$formula_factor_info) &&
    !is.null(model$formula_factor_info[[term]])
}

#' Predictions on the Observed Levels of a Discrete Term
#'
#' Re-evaluates one discrete model term on its fitted observations. This
#' preserves the general [predict.mfp2()] API while ensuring that
#' `plot(..., terms_seq = "equidistant")` does not display artificial
#' intermediate binary or factor-level values.
#'
#' @param model An `"mfp2"` model object.
#' @param term Character scalar naming the term to be plotted.
#' @param type Prediction type passed to [predict.mfp2()], e.g. `"terms"` or
#'   `"contrasts"`.
#' @param ref Named list of reference values used by
#'   `predict(type = "contrasts")`; may be `NULL`.
#' @param alpha Numeric confidence level for prediction intervals passed to
#'   [predict.mfp2()].
#'
#' @return A data frame of predictions on the fitted observations for
#'   `term`.
#'
#' @keywords internal
#' @noRd
mfp2_plot_binary_prediction <- function(model, term, type, ref, alpha) {
  term_ref <- NULL
  if (identical(type, "contrasts") &&
      !is.null(ref) &&
      term %in% names(ref) &&
      !is.null(ref[[term]])) {
    term_ref <- stats::setNames(list(ref[[term]]), term)
  }

  out <- predict(
    model,
    type = type,
    terms = term,
    ref = term_ref,
    terms_seq = "data",
    alpha = alpha
  )

  out[[term]]
}

#' Prepare Prediction Data for a Discrete-Term Plot
#'
#' Collapses repeated observed-value predictions to one fitted estimate per
#' level, preserving the fitted factor-level order for discrete axes,
#' including factors with more than two levels. For a spike-at-zero
#' binary-only term, translates the internal indicator `I(x == 0)` into
#' semantic labels on the original covariate scale (`"x = 0"` and
#' `"x > 0"`).
#'
#' @param model An `"mfp2"` model object.
#' @param term Character scalar naming the term to be plotted.
#' @param df Data frame of per-observation predictions for `term`.
#' @param residual_df Optional data frame of partial residuals matched to the
#'   original observations, used to overlay per-level residuals.
#'
#' @return A data frame with one row per fitted level of `term`, in fitted
#'   order, plus optional residual annotations.
#'
#' @keywords internal
#' @noRd
mfp2_plot_prepare_binary_data <- function(model, term, df, residual_df = NULL) {
  binary_df <- df[!duplicated(df$variable), , drop = FALSE]

  is_saz_binary_only <- FALSE
  if (!is.null(model$spike_dec) &&
      term %in% names(model$spike_dec) &&
      !is.na(model$spike_dec[[term]])) {
    is_saz_binary_only <- identical(
      as.integer(model$spike_dec[[term]]),
      as.integer(saz_decision_codes[["binary_only"]])
    )
  }

  if (is_saz_binary_only) {
    zero_label <- paste0(term, " = 0")
    positive_label <- paste0(term, " > 0")
    level_order <- c(zero_label, positive_label)

    saz_label <- function(x) {
      indicator <- suppressWarnings(as.numeric(as.character(x)))
      out <- rep(NA_character_, length(indicator))
      out[indicator == 1] <- zero_label
      out[indicator == 0] <- positive_label
      out
    }

    binary_df$variable <- factor(
      saz_label(binary_df$variable),
      levels = level_order
    )
    binary_df <- binary_df[order(binary_df$variable), , drop = FALSE]

    if (!is.null(residual_df)) {
      residual_df$variable <- factor(
        saz_label(residual_df$variable),
        levels = level_order
      )
    }

    return(list(
      fitted = binary_df,
      residual = residual_df,
      continuous_axis = FALSE,
      breaks = level_order
    ))
  }

  if (is.numeric(binary_df$variable)) {
    binary_df <- binary_df[order(binary_df$variable), , drop = FALSE]
    return(list(
      fitted = binary_df,
      residual = residual_df,
      continuous_axis = TRUE,
      breaks = as.numeric(binary_df$variable)
    ))
  }

  factor_info <- if (!is.null(model$formula_factor_info)) {
    model$formula_factor_info[[term]]
  } else {
    NULL
  }

  level_order <- if (!is.null(factor_info) &&
                     length(factor_info$levels) > 0L) {
    # Use the fitted factor metadata rather than observed ordering so the
    # reference and non-reference levels appear in their original order.
    as.character(factor_info$levels)
  } else {
    unique(as.character(binary_df$variable))
  }

  binary_df$variable <- factor(
    as.character(binary_df$variable),
    levels = level_order
  )
  binary_df <- binary_df[order(binary_df$variable), , drop = FALSE]

  if (!is.null(residual_df)) {
    residual_df$variable <- factor(
      as.character(residual_df$variable),
      levels = level_order
    )
  }

  list(
    fitted = binary_df,
    residual = residual_df,
    continuous_axis = FALSE,
    breaks = level_order
  )
}


#' Attach Model Residuals to an Observed-Value Prediction Frame
#'
#' Attaches fitted-model residuals to the observed-value prediction frame
#' produced by `predict(model, terms_seq = "data")` for each requested term.
#' Those predictions are constructed from `model$x_original`, so
#' fitted-observation order is the governing invariant. Row-count agreement
#' is validated explicitly to prevent R from recycling a shorter residual
#' vector when its length happens to divide the number of rows.
#'
#' @param pred_data Data frame of observed-value predictions for one term.
#' @param resid Numeric vector of fitted-model residuals aligned with the
#'   original fitting rows.
#' @param model An `"mfp2"` model object.
#'
#' @return `pred_data` with a `residual` column appended.
#'
#' @keywords internal
#' @noRd
mfp2_plot_attach_residuals <- function(pred_data, resid, model) {
  if (is.null(model$x_original) || is.null(nrow(model$x_original))) {
    stop(
      "Cannot construct component-plus-residual plots: the fitted object ",
      "does not contain its original predictor rows.",
      call. = FALSE
    )
  }

  if (!is.atomic(resid) || !is.null(dim(resid))) {
    stop(
      "Cannot construct component-plus-residual plots: model residuals must ",
      "be a vector.",
      call. = FALSE
    )
  }

  expected_n <- nrow(model$x_original)
  if (length(resid) != expected_n) {
    stop(
      "Cannot construct component-plus-residual plots: the model contains ",
      expected_n, " fitted observations but ", length(resid), " residuals.",
      call. = FALSE
    )
  }

  observation_ids <- rownames(model$x_original)
  residual_ids <- names(resid)

  if (!is.null(observation_ids) && !is.null(residual_ids)) {
    # When identifiers are available, use them rather than assuming the
    # residual method retained fitted-data order. Duplicate, missing, or
    # different identifiers make a one-to-one alignment impossible.
    invalid_ids <- anyNA(observation_ids) || anyNA(residual_ids) ||
      anyDuplicated(observation_ids) > 0L ||
      anyDuplicated(residual_ids) > 0L ||
      !setequal(observation_ids, residual_ids)

    if (invalid_ids) {
      stop(
        "Cannot align model residuals with the fitted observations: their ",
        "row identifiers are missing, duplicated, or different.",
        call. = FALSE
      )
    }

    resid <- resid[match(observation_ids, residual_ids)]
  }

  for (term in names(pred_data)) {
    term_data <- pred_data[[term]]
    term_n <- if (is.data.frame(term_data)) nrow(term_data) else NA_integer_

    if (is.na(term_n) || term_n != expected_n) {
      term_n_label <- if (is.na(term_n)) "a non-data-frame result" else term_n
      stop(
        "Cannot attach residuals for term `", term, "`: its prediction ",
        "frame contains ", term_n_label, " but ", expected_n,
        " fitted observations were expected.",
        call. = FALSE
      )
    }

    term_resid <- resid
    term_ids <- rownames(term_data)
    default_term_ids <- as.character(seq_len(expected_n))

    if (!is.null(observation_ids) && !is.null(term_ids) &&
        !identical(term_ids, observation_ids)) {
      if (!anyNA(term_ids) && !anyDuplicated(term_ids) &&
          setequal(term_ids, observation_ids)) {
        # Some grouped prediction frames preserve source row names. If such a
        # frame has been reordered, align residuals to that frame explicitly.
        term_resid <- resid[match(term_ids, observation_ids)]
      } else if (!identical(term_ids, default_term_ids)) {
        # Default sequential row names carry no observation identity and use
        # the documented x_original order. Non-default incompatible names,
        # however, are evidence that positional attachment is unsafe.
        stop(
          "Cannot attach residuals for term `", term, "`: its prediction ",
          "row identifiers do not match the fitted observations.",
          call. = FALSE
        )
      }
    }

    # Remove residual names after alignment. The prediction frame already
    # carries the row order, and retaining unrelated names on this column can
    # obscure that the assignment is deliberately positional at this point.
    term_data$resid <- unname(term_resid)
    pred_data[[term]] <- term_data
  }

  pred_data
}

#' Internal Implementation for `plot.mfp2()` and `fracplot()`
#'
#' Shared implementation of the diagnostic-plot output built by both
#' [plot.mfp2()] and the deprecated [fracplot()] wrapper. Keeping the
#' plotting logic here allows the preferred `plot()` method to run without
#' passing through the deprecated wrapper and its deprecation warning.
#'
#' @param model An `"mfp2"` model object.
#' @param terms Character vector of term names to plot, or `NULL` for all.
#' @param partial_only Logical. If `TRUE`, plot only partial residuals
#'   without the fitted curve.
#' @param ... Additional graphical parameters forwarded to the plotting
#'   layer.
#'
#' @return A list of `ggplot` objects, one per plotted term.
#'
#' @keywords internal
#' @noRd
plot_mfp2_impl <- function(model,
                           terms = NULL,
                           partial_only = FALSE,
                           type = c("terms", "contrasts"),
                           ref = NULL,
                           terms_seq = c("data", "equidistant"),
                           alpha = 0.05,
                           show_titles = TRUE,
                           show_grid = FALSE,
                           color_points = "#AAAAAA",
                           color_line = "red",
                           color_line_spike = "red",
                           color_fill = "#000000",
                           shape = 1,
                           size_points = 1,
                           size_points_spike = 2,
                           linetype = "solid",
                           linewidth = 1,
                           alpha_fill = 0.1,
                           terms_seq_missing = FALSE) {


  # Fail clearly if ggplot2 is unavailable at runtime.
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  # Validate the fitted object before prediction.
  if (!inherits(model, "mfp2"))
    stop("The model entered is not an mfp2 object.", call. = FALSE)

  # Contrast references must be supplied as a named list.
  if (!is.null(ref) && !is.list(ref))
    stop("`ref` must be a list or NULL.", call. = FALSE)

  if (!is.logical(show_titles) || length(show_titles) != 1L || is.na(show_titles)) {
    stop("show_titles must be a single non-missing logical value.", call. = FALSE)
  }
  if (!is.logical(show_grid) || length(show_grid) != 1L || is.na(show_grid)) {
    stop("show_grid must be a single non-missing logical value.", call. = FALSE)
  }

  # set defaults depending on type
  type <- match.arg(type)
  if (type == "contrasts") {
    partial_only <- TRUE
    if (terms_seq_missing)
      terms_seq <- "equidistant"
  }

  terms_seq <- match.arg(terms_seq)

  # Select once before repeated prediction calls (grids, observed values,
  # or ordinal cutoffs), and let this warning identify every omitted term.
  retained_terms <- get_selected_variable_names(model)
  terms <- if (is.null(terms)) retained_terms else terms
  terms <- mfp2_select_prediction_terms(terms, retained_terms, caller = "plot()")
  if (length(terms) == 0L) return(list())

  # Multinomial fits carry a separate coefficient row per non-reference logit,
  # so each term yields one partial predictor per logit. These are shown in a
  # single faceted plot per predictor (one facet per logit) rather than the
  # single-curve layout used for scalar-response families.
  if (identical(model$family_string, "multinomial") ||
      (identical(model$family_string, "ordinal") && type == "terms")) {
    return(mfp2_plot_multinomial_impl(
      model = model,
      terms = terms,
      partial_only = partial_only,
      type = type,
      ref = ref,
      terms_seq = terms_seq,
      alpha = alpha,
      show_titles = show_titles,
      show_grid = show_grid,
      color_points = color_points,
      color_line = color_line,
      color_line_spike = color_line_spike,
      color_fill = color_fill,
      shape = shape,
      size_points = size_points,
      size_points_spike = size_points_spike,
      linetype = linetype,
      linewidth = linewidth,
      alpha_fill = alpha_fill,
      ordinal = identical(model$family_string, "ordinal")
    ))
  }

  pred <- predict(model,
                  type = type,
                  terms = terms,
                  ref = ref,
                  terms_seq = terms_seq,
                  alpha = alpha)
  # for points also need the data point predictions
  pred_data <- pred
  if (!partial_only && terms_seq != "data") {
    pred_data <- predict(model,
                         type = "terms",
                         terms = terms,
                         terms_seq = "data")
  }

  if (length(pred) == 0L) return(list())

  ylab <- if (isTRUE(partial_only) && identical(type, "terms")) {
    mfp2_plot_partial_ylabel(model)
  } else {
    "Partial Predictor"
  }

  if (!partial_only) {
    # compute residuals to plot points
    # for glm, deviance residuals are required
    # while for cox martingale residuals
    resid <- if (mfp2_family_is_ph(model$family_string)) {
      model$residuals
    } else {
      # Deviance residuals for the component-plus-residual overlay, matching
      # Stata's `fracplot`. GEE uses the same convention: residuals.mfp2()
      # computes deviance residuals for GEE from the marginal fitted mean even
      # though geepack does not expose them directly.
      stats::residuals(model, type = "deviance")
    }
    if (identical(model$family_string, "finegray") &&
        !is.null(model$mfp2_finegray_row_map) &&
        length(resid) == length(model$mfp2_finegray_row_map)) {
      row_map <- model$mfp2_finegray_row_map
      aggregated <- rowsum(
        resid,
        row_map,
        reorder = FALSE
      )
      n_original <- if (is.null(model$nobs)) max(row_map) else model$nobs
      resid_original <- rep(NA_real_, n_original)
      resid_original[unique(row_map)] <- as.numeric(aggregated)
      resid <- resid_original
    }
    # Prediction frames and residuals must describe the same fitted rows.
    # The helper checks counts and uses observation names for alignment when
    # both the stored predictors and residual vector provide them.
    pred_data <- mfp2_plot_attach_residuals(pred_data, resid, model)

    # y label for the plot
    ylab <- "Partial Predictor + residuals"
  }

  # Preallocate the list for plots
  plots <- setNames(vector("list", length(names(pred))), names(pred))

  # in the calls to ggplot2::aes use the .data pronoun to avoid notes
  # generated by R CMD CHECK about missing bindings of global variables
  for (v in names(pred)) {
    df <- pred[[v]]
    # Check if variable was assessed by the spike-at-zero algorithm.
    # Use stored spike metadata rather than catzero_list: for positive-part-only
    # spike terms, the binary indicator is intentionally suppressed, but the
    # variable is still a spike-at-zero covariate and must be plotted as a
    # positive-value function.
    is_spike <- FALSE
    if (!is.null(model$spike) && v %in% names(model$spike)) {
      is_spike <- isTRUE(model$spike[[v]])
    } else if (!is.null(model$fp_terms) &&
               v %in% rownames(model$fp_terms) &&
               "spike" %in% colnames(model$fp_terms)) {
      is_spike <- isTRUE(model$fp_terms[v, "spike"])
    }

    spike_dec_v <- saz_decision_codes[["continuous_only"]]
    if (is_spike &&
        !is.null(model$spike_dec) &&
        v %in% names(model$spike_dec) &&
        !is.na(model$spike_dec[[v]])) {
      spike_dec_v <- as.integer(model$spike_dec[[v]])
    }

    is_acd <- !is.null(model$fp_terms) &&
      v %in% rownames(model$fp_terms) &&
      "acd" %in% colnames(model$fp_terms) &&
      isTRUE(model$fp_terms[v, "acd"])

    # zero = TRUE defines the uncentered continuous basis only on x > 0 and
    # assigns zero to that basis at x = 0. Ordinary centering may subsequently
    # move the zero-row design value away from zero. This is distinct from
    # catzero/spike modelling, which can include a separately estimated zero
    # indicator, but every zero-handled continuous curve must avoid drawing an
    # interpolating segment through the exact-zero point.
    has_zero_indicator <- !is.null(model$catzero) &&
      v %in% names(model$catzero) &&
      isTRUE(model$catzero[[v]])

    is_zero_handled <- (!is.null(model$zero) &&
      v %in% names(model$zero) &&
      isTRUE(model$zero[[v]])) || is_spike || has_zero_indicator

    is_binary <- (is_spike &&
                    spike_dec_v == saz_decision_codes[["binary_only"]]) ||
      mfp2_plot_term_is_binary(model, v)

    # A formula factor is discrete even when it has more than two levels.
    # Treat it like a binary effect for plotting purposes: predict only at the
    # fitted levels and display estimates with confidence intervals, rather
    # than drawing a continuous curve between category labels.
    is_factor <- mfp2_plot_term_is_factor(model, v)
    is_discrete <- is_binary || is_factor

    binary_plot_data <- NULL
    if (is_discrete) {
      binary_df <- mfp2_plot_binary_prediction(
        model = model,
        term = v,
        type = type,
        ref = ref,
        alpha = alpha
      )

      if (!is.null(binary_df)) {
        binary_plot_data <- mfp2_plot_prepare_binary_data(
          model = model,
          term = v,
          df = binary_df,
          residual_df = if (!partial_only) pred_data[[v]] else NULL
        )
        df <- binary_plot_data$fitted
        if (!partial_only) {
          pred_data[[v]] <- binary_plot_data$residual
        }
      }
    }

    power_label <- paste0(model$fp_powers[[v]], collapse = ", ")
    continuous_label <- if (is_acd) {
      sprintf("ACD FP(%s)", power_label)
    } else {
      sprintf("FP(%s)", power_label)
    }

    # Title with FP/ACD powers and explicit spike-at-zero representation.
    # The internal spike_dec code is not exposed in the plot title.
    # Factor variables get "Categorical" instead of a misleading FP label.
    plot_title <- if (is_factor) {
      "Categorical"
    } else if (is_spike) {
      saz_decision_label(
        spike_dec_v,
        style = "plot_title",
        continuous_label = continuous_label
      )
    } else if (is_zero_handled && has_zero_indicator) {
      paste0(continuous_label, " (x > 0) + zero indicator")
    } else if (is_zero_handled) {
      paste0(continuous_label, " (positive part; no zero indicator)")
    } else {
      continuous_label
    }

    p <- ggplot2::ggplot(data = df, ggplot2::aes(x = .data$variable, y = .data$value)) +
      ggplot2::xlab(v) + ggplot2::ylab(ylab) +
      ggplot2::theme_bw()
    if (!show_grid) {
      p <- p + ggplot2::theme(panel.grid = ggplot2::element_blank())
    }

    if (show_titles) {
      p <- p + ggplot2::ggtitle(plot_title)
    }

    # Residuals go first
    if (!partial_only) {
      p <- p + ggplot2::geom_point(data = pred_data[[v]],
                                   ggplot2::aes(y = .data$value + .data$resid),
                                   color = color_points,
                                   size = size_points,
                                   shape = shape)
    }

    # Then fitted line/ribbon on top.
    # A zero-handled positive-value FP/ACD curve must not extend through x = 0.
    # Binary-only SAZ terms were handled by the discrete branch above. All
    # remaining zero-handled terms have a continuous positive component; only
    # catzero or spike models use the special zero-component styling.
    if (is_discrete && !is.null(binary_plot_data)) {
      # Discrete terms (binary variables and formula factors of any size) use
      # one estimate and confidence interval per fitted level. Do not imply
      # values between unordered/ordered factor levels with a continuous line.
      bin_df <- binary_plot_data$fitted

      p <- p + ggplot2::geom_point(
        data = bin_df,
        ggplot2::aes(x = .data$variable, y = .data$value),
        color = color_line,
        size = size_points_spike
      ) +
        ggplot2::geom_errorbar(
          data = bin_df,
          ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
          width = 0.2,
          color = color_line
        )

      if (isTRUE(binary_plot_data$continuous_axis)) {
        p <- p + ggplot2::scale_x_continuous(
          breaks = binary_plot_data$breaks,
          labels = format(binary_plot_data$breaks, trim = TRUE),
          minor_breaks = NULL
        )
      } else {
        p <- p + ggplot2::scale_x_discrete(
          breaks = binary_plot_data$breaks,
          limits = binary_plot_data$breaks
        )
      }

    } else if (is_zero_handled) {
      # The continuous transformation is defined only for x > 0. At x = 0,
      # zero-only handling shows the endpoint from the same centered continuous
      # term with ordinary styling; catzero and spike models retain
      # zero-component styling.
      pos_df <- df[df$variable > 0, , drop = FALSE]
      pos_df <- pos_df[order(pos_df$variable), , drop = FALSE]
      zero_df <- df[df$variable == 0, , drop = FALSE]
      if (nrow(zero_df) > 1L) {
        zero_df <- zero_df[1L, , drop = FALSE]
      }

      p <- p + ggplot2::geom_line(data = pos_df,
                                  ggplot2::aes(x = .data$variable, y = .data$value),
                                  linewidth = linewidth,
                                  linetype = linetype,
                                  color = color_line) +
        ggplot2::geom_ribbon(data = pos_df,
                             ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                             alpha = alpha_fill, fill = color_fill)

      if (nrow(zero_df) > 0) {
        zero_color <- if (is_spike || has_zero_indicator) {
          color_line_spike
        } else {
          color_line
        }

        p <- p + ggplot2::geom_point(data = zero_df,
                                     ggplot2::aes(y = .data$value),
                                     color = zero_color,
                                     size = size_points_spike) +
          ggplot2::geom_errorbar(data = zero_df,
                                 ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                                 width = 0.2,
                                 color = zero_color)
      }

    } else {
      # Regular continuous plot (non-spike variable)
      p <- p + ggplot2::geom_line(linewidth = linewidth,
                                  linetype = linetype,
                                  color = color_line) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                             alpha = alpha_fill, fill = color_fill)
    }

    # Add residual points if needed
    # if (!partial_only) {
    #   p <- p + ggplot2::geom_point(data = pred_data[[v]],
    #                                ggplot2::aes(y = .data$value + .data$resid),
    #                                color = color_points,
    #                                size = size_points,
    #                                shape = shape)
    # }
    plots[[v]] <- p
  }

  plots
}

#' Faceted Partial-Predictor Plots for Categorical `mfp2` Models
#'
#' A baseline-category multinomial model shares the fractional-polynomial powers
#' across logits but estimates a separate coefficient for each non-reference
#' logit, so every predictor has one partial predictor per logit. Ordinal
#' models have a separate intercept for every cumulative-response cutoff.
#' This helper draws a single \pkg{ggplot2} object per predictor, faceted by
#' logit or cutoff.
#'
#' Component-plus-residual plots use a deviance residual for each logit or
#' ordinal cutoff. Each ordinal facet includes its fitted cutoff intercept,
#' the corresponding term coefficient covariance, and the matching residual.
#'
#' @param model,terms,partial_only,type,ref,terms_seq,alpha,show_titles,show_grid Plot
#'   inputs forwarded from [plot.mfp2()] with the same meaning.
#' @param ordinal Logical; facet by ordinal cutoff instead of multinomial logit.
#' @param color_points,color_line,color_line_spike,color_fill,shape,size_points,size_points_spike,linetype,linewidth,alpha_fill
#'   Styling controls forwarded from [plot.mfp2()].
#'
#' @return A named list of \pkg{ggplot2} objects, one per plotted term.
#'
#' @keywords internal
#' @noRd
mfp2_plot_multinomial_impl <- function(model,
                                       terms,
                                       partial_only,
                                       type,
                                       ref,
                                       terms_seq,
                                       alpha,
                                       show_titles,
                                       show_grid,
                                       color_points,
                                       color_line,
                                       color_line_spike,
                                       color_fill,
                                       shape,
                                       size_points,
                                       size_points_spike,
                                       linetype,
                                       linewidth,
                                       alpha_fill,
                                       ordinal = FALSE) {
  cutoff_ids <- if (ordinal) seq_along(model$mfp2_ordinal_intercepts) else NULL
  predict_facets <- function(sequence) {
    if (!ordinal) {
      return(predict(model, type = type, terms = terms, ref = ref,
                     terms_seq = sequence, alpha = alpha))
    }
    per_cutoff <- lapply(cutoff_ids, function(j) {
      predict(model, type = "terms", terms = terms, terms_seq = sequence,
              alpha = alpha, ordinal_cutoff = j)
    })
    if (length(per_cutoff[[1L]]) == 0L) return(list())
    lapply(names(per_cutoff[[1L]]), function(term) {
      do.call(rbind, lapply(seq_along(cutoff_ids), function(j) {
        df <- per_cutoff[[j]][[term]]
        df$logit <- as.character(j)
        df
      }))
    }) |> stats::setNames(names(per_cutoff[[1L]]))
  }
  pred <- predict_facets(terms_seq)
  if (length(pred) == 0L) return(list())

  ylab <- if (!partial_only) {
    "Partial Predictor + residuals"
  } else if (identical(type, "terms")) {
    mfp2_plot_partial_ylabel(model)
  } else {
    "Partial Predictor"
  }

  # Facet strips read "<non-reference class> vs <reference class>".
  if (ordinal) {
    logit_levels <- as.character(cutoff_ids)
    levels <- model$mfp2_ordinal_levels
    facet_labels <- stats::setNames(vapply(cutoff_ids, function(j) {
      paste0(paste(levels[seq.int(j + 1L, length(levels))], collapse = "/"),
             " vs ", paste(levels[seq_len(j)], collapse = "/"))
    }, character(1L)), logit_levels)
  } else {
    reference <- model$reference_class
    logit_levels <- rownames(model$mfp2_coefficient_matrix)
    facet_labels <- stats::setNames(
      paste0(logit_levels, " vs ", reference), logit_levels
    )
  }

  # Component-plus-residual overlays need data-scale per-logit predictions and
  # the per-logit deviance residuals. Both are ordered by fitted observation.
  resid_matrix <- NULL
  pred_data <- NULL
  if (!partial_only) {
    pred_data <- if (terms_seq == "data") {
      pred
    } else {
      predict_facets("data")
    }
    resid_matrix <- if (ordinal) {
      values <- lapply(cutoff_ids, function(j) {
        stats::residuals(model, type = "deviance", kint = j)
      })
      if (any(lengths(values) != model$nobs)) {
        stop("Ordinal cutoff residuals do not align with fitted observations.",
             call. = FALSE)
      }
      matrix(unlist(values, use.names = FALSE), ncol = length(cutoff_ids),
             dimnames = list(NULL, logit_levels))
    } else {
      stats::residuals(model, type = "deviance")
    }
  }

  plots <- stats::setNames(vector("list", length(names(pred))), names(pred))

  for (v in names(pred)) {
    df <- pred[[v]]
    df$logit <- factor(df$logit, levels = logit_levels)

    is_spike <- FALSE
    if (!is.null(model$spike) && v %in% names(model$spike)) {
      is_spike <- isTRUE(model$spike[[v]])
    } else if (!is.null(model$fp_terms) &&
               v %in% rownames(model$fp_terms) &&
               "spike" %in% colnames(model$fp_terms)) {
      is_spike <- isTRUE(model$fp_terms[v, "spike"])
    }
    spike_dec_v <- if (is_spike && !is.null(model$spike_dec) &&
                       v %in% names(model$spike_dec) &&
                       !is.na(model$spike_dec[[v]])) {
      as.integer(model$spike_dec[[v]])
    } else {
      saz_decision_codes[["continuous_only"]]
    }
    is_binary_only <- is_spike &&
      spike_dec_v == saz_decision_codes[["binary_only"]]
    is_factor <- mfp2_plot_term_is_factor(model, v)
    is_discrete <- is_binary_only || is_factor ||
      mfp2_plot_term_is_binary(model, v)
    has_zero_indicator <- !is.null(model$catzero) &&
      isTRUE(model$catzero[[v]])
    is_zero_handled <- (!is.null(model$zero) &&
      isTRUE(model$zero[[v]])) || is_spike || has_zero_indicator
    is_acd <- !is.null(model$fp_terms) &&
      v %in% rownames(model$fp_terms) &&
      "acd" %in% colnames(model$fp_terms) &&
      isTRUE(model$fp_terms[v, "acd"])

    power_label <- paste0(model$fp_powers[[v]], collapse = ", ")
    continuous_label <- if (is_acd) {
      sprintf("ACD FP(%s)", power_label)
    } else {
      sprintf("FP(%s)", power_label)
    }
    plot_title <- if (is_factor) {
      "Categorical"
    } else if (is_spike) {
      saz_decision_label(spike_dec_v, style = "plot_title",
                         continuous_label = continuous_label)
    } else if (is_zero_handled && has_zero_indicator) {
      paste0(continuous_label, " (x > 0) + zero indicator")
    } else if (is_zero_handled) {
      paste0(continuous_label, " (positive part; no zero indicator)")
    } else {
      continuous_label
    }

    # A binary-only spike is predicted on the internal indicator scale:
    # 1 means x = 0 and 0 means x > 0. Keep the original-scale labels in
    # every facet and in the observed residual layer.
    if (is_binary_only) {
      zero_label <- paste0(v, " = 0")
      positive_label <- paste0(v, " > 0")
      spike_levels <- c(zero_label, positive_label)
      label_spike <- function(indicator) {
        factor(ifelse(indicator == 1, zero_label, positive_label),
               levels = spike_levels)
      }
      df$variable <- label_spike(df$variable)
    } else if (is_factor) {
      factor_levels <- model$formula_factor_info[[v]]$levels
      df$variable <- factor(as.character(df$variable), levels = factor_levels)
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$variable, y = .data$value)) +
      ggplot2::xlab(v) + ggplot2::ylab(ylab) +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(
        ggplot2::vars(.data$logit), nrow = 1L,
        labeller = ggplot2::as_labeller(facet_labels)
      )
    if (!show_grid) {
      p <- p + ggplot2::theme(panel.grid = ggplot2::element_blank())
    }
    if (show_titles) {
      p <- p + ggplot2::ggtitle(plot_title)
    }

    # Residual points first, so the fitted curve is drawn on top.
    if (!partial_only && !is.null(resid_matrix)) {
      resid_df <- do.call(rbind, lapply(logit_levels, function(q) {
        pd <- pred_data[[v]]
        pd <- pd[as.character(pd$logit) == q, , drop = FALSE]
        if (nrow(pd) != nrow(resid_matrix)) {
          stop("Cutoff predictions and residuals have different row counts.",
               call. = FALSE)
        }
        if (is_binary_only) pd$variable <- label_spike(pd$variable)
        if (is_factor) {
          pd$variable <- factor(as.character(pd$variable),
                                levels = factor_levels)
        }
        pd$resid <- as.numeric(resid_matrix[, q])
        pd$logit <- factor(q, levels = logit_levels)
        pd
      }))
      p <- p + ggplot2::geom_point(
        data = resid_df,
        ggplot2::aes(y = .data$value + .data$resid),
        color = color_points, size = size_points, shape = shape
      )
    }

    if (is_discrete) {
      # One estimate and interval per fitted level; no interpolating line.
      level_df <- df[!duplicated(df[c("variable", "logit")]), , drop = FALSE]
      if (is_binary_only) {
        level_df$variable <- factor(as.character(level_df$variable),
                                    levels = spike_levels)
      }
      p <- p + ggplot2::geom_point(
        data = level_df,
        ggplot2::aes(x = .data$variable, y = .data$value),
        color = color_line, size = size_points_spike
      ) +
        ggplot2::geom_errorbar(
          data = level_df,
          ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
          width = 0.2, color = color_line
        )
      if (is_binary_only) {
        p <- p + ggplot2::scale_x_discrete(breaks = spike_levels,
                                           limits = spike_levels)
      } else if (is_factor) {
        p <- p + ggplot2::scale_x_discrete(breaks = factor_levels,
                                           limits = factor_levels)
      }
    } else if (is_zero_handled) {
      # The continuous basis is defined for x > 0. Its value at exactly zero
      # is a separate endpoint, potentially with an estimated zero indicator.
      pos_df <- df[df$variable > 0, , drop = FALSE]
      pos_df <- pos_df[order(pos_df$logit, pos_df$variable), , drop = FALSE]
      zero_df <- df[df$variable == 0, , drop = FALSE]
      zero_df <- zero_df[!duplicated(zero_df$logit), , drop = FALSE]
      p <- p + ggplot2::geom_line(
        data = pos_df, linewidth = linewidth, linetype = linetype,
        color = color_line
      ) + ggplot2::geom_ribbon(
        data = pos_df,
        ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
        alpha = alpha_fill, fill = color_fill
      )
      if (nrow(zero_df) > 0L) {
        zero_color <- if (is_spike || has_zero_indicator) {
          color_line_spike
        } else {
          color_line
        }
        p <- p + ggplot2::geom_point(
          data = zero_df, color = zero_color, size = size_points_spike
        ) + ggplot2::geom_errorbar(
          data = zero_df,
          ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
          width = 0.2, color = zero_color
        )
      }
    } else {
      ord <- order(df$logit, df$variable)
      df_line <- df[ord, , drop = FALSE]
      p <- p + ggplot2::geom_line(
        data = df_line, linewidth = linewidth, linetype = linetype,
        color = color_line
      ) +
        ggplot2::geom_ribbon(
          data = df_line,
          ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
          alpha = alpha_fill, fill = color_fill
        )
    }

    plots[[v]] <- p
  }

  plots
}

#' Deprecated Fractional-Polynomial Plotting Interface
#'
#' `fracplot()` is retained temporarily for compatibility with code written
#' for earlier versions of `mfp2`. Use [plot.mfp2()] or `plot()` instead.
#'
#' The function forwards its arguments to [plot.mfp2()] and returns the same
#' plot objects, except that its own defaults for `color_line` and
#' `color_line_spike` are `"red"`, whereas [plot.mfp2()] uses `"black"`.
#'
#' @inheritParams plot.mfp2
#' @param model A fitted object of class `"mfp2"`.
#'
#' @return
#' A named list of `ggplot2` objects, one for each plotted term, identical to
#' the result returned by [plot.mfp2()] called with the same arguments
#' (including the colours).
#'
#' @seealso [plot.mfp2()]
#'
#' @examples
#' data("prostate")
#'
#' fit <- mfp2(
#'   lpsa ~ fp(age) + fp(cavol) + fp(weight) + svi,
#'   data = prostate,
#'   verbose = FALSE
#' )
#'
#' # Recommended interface
#' plots <- plot(fit, terms = "cavol")
#' plots
#'
#' \donttest{
#' # Deprecated compatibility interface
#' old_plots <- fracplot(fit, terms = "cavol")
#' }
#'
#' @export
fracplot <- function(model,
                     terms = NULL,
                     partial_only = FALSE,
                     type = c("terms", "contrasts"),
                     ref = NULL,
                     terms_seq = c("data", "equidistant"),
                     alpha = 0.05,
                     show_titles = TRUE,
                     show_grid = FALSE,
                     color_points = "#AAAAAA",
                     color_line = "red",
                     color_line_spike = "red",
                     color_fill = "#000000",
                     shape = 1,
                     size_points = 1,
                     size_points_spike = 2,
                     linetype = "solid",
                     linewidth = 1,
                     alpha_fill = 0.1) {
  terms_seq_missing <- missing(terms_seq)

  .Deprecated(
    new = "plot",
    package = "mfp2",
    old = "fracplot",
    msg = paste0(
      "`fracplot()` is deprecated and will be removed in a future version ",
      "of mfp2. Use `plot()` with a fitted `mfp2` object instead, for example ",
      "`plot(fit)`."
    )
  )

  plot_args <- list(
    x = model,
    terms = terms,
    partial_only = partial_only,
    type = type,
    ref = ref,
    alpha = alpha,
    show_titles = show_titles,
    show_grid = show_grid,
    color_points = color_points,
    color_line = color_line,
    color_line_spike = color_line_spike,
    color_fill = color_fill,
    shape = shape,
    size_points = size_points,
    size_points_spike = size_points_spike,
    linetype = linetype,
    linewidth = linewidth,
    alpha_fill = alpha_fill
  )

  # Preserve plot.mfp2()'s type-dependent default when terms_seq was omitted.
  if (!terms_seq_missing) {
    plot_args$terms_seq <- terms_seq
  }

  do.call(plot.mfp2, plot_args)
}
