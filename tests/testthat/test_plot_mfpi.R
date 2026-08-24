# Public contract tests for the exported plot.mfpi() method.

.plot_mfpi_fit_cache_v <- new.env(parent = emptyenv())

get_plot_mfpi_fit_v <- function() {
  if (!exists("fit", envir = .plot_mfpi_fit_cache_v, inherits = FALSE)) {
    dat <- make_mfpi_factor_data(n = 180L)
    fit <- mfpi(
      y ~ trt + fp(x, df = 1),
      data = dat,
      group_var = "trt",
      cont_vars = "x",
      cont_var_forms = c(x = "linear"),
      flex = "flex1",
      cycles = 1,
      df = 1,
      select = 1,
      alpha = 1,
      p_interact = 1,
      shift = 0,
      scale = 1,
      center = FALSE,
      verbose = FALSE
    )
    assign("fit", fit, envir = .plot_mfpi_fit_cache_v)
  }

  get("fit", envir = .plot_mfpi_fit_cache_v, inherits = FALSE)
}


test_that("plot.mfpi returns documented fitted and difference plot structures", {
  skip_if_not_installed("ggplot2")
  fit <- get_plot_mfpi_fit_v()

  fitted <- plot(
    fit,
    terms = "x",
    plot_type = "fitted",
    auto_print = FALSE,
    show_ci_fitted = TRUE
  )
  difference <- plot(
    fit,
    terms = "x",
    plot_type = "difference",
    auto_print = FALSE,
    show_null_line = TRUE,
    show_maineffect_line = FALSE
  )

  expect_type(fitted, "list")
  expect_identical(names(fitted), "x")
  expect_length(fitted$x, 1L)
  expect_s3_class(fitted$x[[1L]], "ggplot")

  expect_type(difference, "list")
  expect_identical(names(difference), "x")
  expect_identical(names(difference$x), names(fitted$x))
  expect_s3_class(difference$x[[1L]], "ggplot")

  fitted_geoms <- vapply(
    fitted$x[[1L]]$layers,
    function(layer) class(layer$geom)[[1L]],
    character(1L)
  )
  difference_geoms <- vapply(
    difference$x[[1L]]$layers,
    function(layer) class(layer$geom)[[1L]],
    character(1L)
  )

  expect_true("GeomLine" %in% fitted_geoms)
  expect_true("GeomRibbon" %in% fitted_geoms)
  expect_true("GeomLine" %in% difference_geoms)
  expect_true("GeomHline" %in% difference_geoms)
  expect_match(fitted$x[[1L]]$labels$subtitle, "treated vs control", fixed = TRUE)
})


test_that("plot.mfpi returns patchwork objects for combined plots", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  fit <- get_plot_mfpi_fit_v()

  combined <- plot(
    fit,
    terms = "x",
    plot_type = "both",
    auto_print = FALSE
  )

  expect_type(combined, "list")
  expect_length(combined$x, 1L)
  expect_s3_class(combined$x[[1L]], "patchwork")
})


test_that("plot.mfpi validates public plotting arguments", {
  skip_if_not_installed("ggplot2")
  fit <- get_plot_mfpi_fit_v()

  expect_error(
    plot(fit, terms = character(), auto_print = FALSE),
    "`terms` must be a non-empty character vector"
  )
  expect_error(
    plot(fit, terms = "x", linewidth = 0, auto_print = FALSE),
    "`linewidth` must be a single positive numeric value"
  )
  expect_error(
    plot(fit, terms = "x", ribbon_alpha = 2, auto_print = FALSE),
    "`ribbon_alpha` must be a single numeric value in [0, 1]",
    fixed = TRUE
  )
  expect_error(
    plot(fit, terms = "x", legend_inside = c(-0.1, 0.5), auto_print = FALSE),
    "`legend_inside` must be a numeric vector of length 2",
    fixed = TRUE
  )
  expect_warning(
    plot(fit, terms = "x", auto_print = FALSE, unused_argument = TRUE),
    "Unused arguments in `...`: unused_argument"
  )
})
