# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# =============================================================================
# 21. plot()
# =============================================================================

# Test purpose: Checks that plot() can be called on a Gaussian mfp2 fit without
# error.
test_that("plot() runs without error for Gaussian model", {
  fit <- mfp2(
    x_prostate,
    y_prostate,
    verbose = FALSE
  )

  expect_error(
    plot(fit),
    NA
  )
})


# Test purpose: Checks that plot() returns a list and does not emit warnings for
# a Gaussian mfp2 model.
test_that("plot() runs without warning for Gaussian models", {
  fit <- mfp2(
    x_prostate,
    y_prostate,
    verbose = FALSE
  )

  expect_warning(
    plots <- plot(fit),
    NA
  )

  expect_type(plots, "list")
})


# Test purpose: Checks that the deprecated fracplot() wrapper still delegates to
# the plotting implementation and returns a list.
test_that("fracplot() is deprecated but still works for Gaussian models", {
  fit <- mfp2(
    x_prostate,
    y_prostate,
    verbose = FALSE
  )

  expect_warning(
    plots <- fracplot(fit),
    regexp = "fracplot.*deprecated.*plot",
    ignore.case = TRUE
  )

  expect_type(plots, "list")
})


# -----------------------------------------------------------------------------
# 21.1 Binary term plotting
# -----------------------------------------------------------------------------

# Test purpose: Numeric binary terms must be plotted as two fitted point
# estimates with vertical confidence intervals, even when an equidistant
# sequence is requested. No interpolating line or confidence ribbon is drawn.
test_that("plot() displays numeric binary terms as two point estimates", {
  skip_if_not_installed("ggplot2")

  set.seed(21011)
  n <- 120L
  binary <- rep(c(2, 5), each = n / 2L)
  y <- 1 + 1.8 * (binary == 5) + stats::rnorm(n, sd = 0.4)

  fit <- mfp2(
    x = cbind(binary = binary),
    y = y,
    df = 1,
    keep = "binary",
    center = FALSE,
    verbose = FALSE
  )

  p <- plot(
    fit,
    terms = "binary",
    partial_only = TRUE,
    terms_seq = "equidistant"
  )[["binary"]]

  geom_classes <- vapply(
    p$layers,
    function(layer) class(layer$geom)[[1L]],
    character(1L)
  )

  expect_true("GeomPoint" %in% geom_classes)
  expect_true("GeomErrorbar" %in% geom_classes)
  expect_false("GeomLine" %in% geom_classes)
  expect_false("GeomRibbon" %in% geom_classes)

  point_layer <- p$layers[[which(geom_classes == "GeomPoint")[[1L]]]]
  errorbar_layer <- p$layers[[which(geom_classes == "GeomErrorbar")[[1L]]]]

  expect_equal(nrow(point_layer$data), 2L)
  expect_equal(nrow(errorbar_layer$data), 2L)
  expect_equal(as.numeric(point_layer$data$variable), c(2, 5))

  x_scale <- p$scales$get_scales("x")
  expect_equal(as.numeric(x_scale$breaks), c(2, 5))
})


# Test purpose: Two-level formula factors retain their fitted level labels and
# use the same point-and-confidence-interval presentation as numeric binaries.
test_that("plot() displays two-level factors as binary effects", {
  skip_if_not_installed("ggplot2")

  set.seed(21012)
  n <- 120L
  group <- factor(
    rep(c("control", "treated"), each = n / 2L),
    levels = c("control", "treated")
  )
  y <- 0.5 + 1.4 * (group == "treated") + stats::rnorm(n, sd = 0.4)
  dat <- data.frame(y = y, group = group)

  fit <- mfp2(
    y ~ group,
    data = dat,
    keep = "group",
    verbose = FALSE
  )

  p <- plot(
    fit,
    terms = "group",
    partial_only = TRUE,
    terms_seq = "equidistant"
  )[["group"]]

  geom_classes <- vapply(
    p$layers,
    function(layer) class(layer$geom)[[1L]],
    character(1L)
  )

  expect_true("GeomPoint" %in% geom_classes)
  expect_true("GeomErrorbar" %in% geom_classes)
  expect_false("GeomLine" %in% geom_classes)
  expect_false("GeomRibbon" %in% geom_classes)

  point_layer <- p$layers[[which(geom_classes == "GeomPoint")[[1L]]]]
  expect_equal(nrow(point_layer$data), 2L)
  expect_equal(
    as.character(point_layer$data$variable),
    c("control", "treated")
  )

  x_scale <- p$scales$get_scales("x")
  expect_equal(x_scale$breaks, c("control", "treated"))
  expect_equal(x_scale$limits, c("control", "treated"))
})


# Test purpose: A SAZ binary-only plot must not expose the internal indicator
# coding, where 1 denotes the zero group. The x-axis instead displays the two
# original-scale groups in the intuitive order zero, then positive.
test_that("plot() labels SAZ binary-only groups on the original scale", {
  skip_if_not_installed("ggplot2")

  data("prostate", package = "mfp2")
  fit <- mfp2(
    lpsa ~ fp(
      pgg45,
      df = 4,
      select = 0.05,
      alpha = 0.05,
      spike = TRUE
    ),
    data = prostate,
    family = "gaussian",
    criterion = "pvalue",
    verbose = FALSE
  )

  expect_identical(
    as.integer(fit$spike_dec[["pgg45"]]),
    as.integer(saz_decision_codes[["binary_only"]])
  )

  p <- plot(
    fit,
    terms = "pgg45",
    partial_only = TRUE
  )[["pgg45"]]

  expected_labels <- c("pgg45 = 0", "pgg45 > 0")
  x_scale <- p$scales$get_scales("x")

  expect_equal(x_scale$breaks, expected_labels)
  expect_equal(x_scale$limits, expected_labels)

  geom_classes <- vapply(
    p$layers,
    function(layer) class(layer$geom)[[1L]],
    character(1L)
  )
  point_layer <- p$layers[[which(geom_classes == "GeomPoint")[[1L]]]]

  expect_equal(
    as.character(point_layer$data$variable),
    expected_labels
  )
  expect_false(any(as.character(point_layer$data$variable) %in% c("0", "1")))
})


# -----------------------------------------------------------------------------
# 21.2 Positive-part ACD plotting
# -----------------------------------------------------------------------------

# Test purpose: zero_vars without catzero or spike fits no separate zero
# coefficient. The ACD curve and ribbon must cover only x > 0, while x = 0 is
# shown once with the ordinary curve colour rather than spike styling.
test_that("plot() separates the zero endpoint from a positive-part ACD curve", {
  skip_if_not_installed("ggplot2")

  set.seed(21013)
  n <- 160L
  exposure <- stats::rexp(n)
  exposure[sample.int(n, size = 40L)] <- 0
  y <- 1 + 4 * (exposure > 0) * stats::pnorm(log(pmax(exposure, 0.05))) +
    stats::rnorm(n, sd = 0.1)

  fit <- mfp2(
    x = matrix(exposure, ncol = 1L, dimnames = list(NULL, "exposure")),
    y = y,
    acd_vars = "exposure",
    zero_vars = "exposure",
    keep = "exposure",
    force_max_fp_vars = "exposure",
    select = 1,
    alpha = 1,
    verbose = FALSE
  )

  expect_true(fit$acd[["exposure"]])
  expect_true(fit$zero[["exposure"]])
  expect_false(fit$catzero[["exposure"]])
  expect_false(fit$spike[["exposure"]])

  p <- plot(
    fit,
    terms = "exposure",
    partial_only = TRUE,
    terms_seq = "equidistant",
    color_line = "navy",
    color_line_spike = "orange"
  )[["exposure"]]

  geom_classes <- vapply(
    p$layers,
    function(layer) class(layer$geom)[[1L]],
    character(1L)
  )

  expect_equal(sum(geom_classes == "GeomLine"), 1L)
  expect_equal(sum(geom_classes == "GeomRibbon"), 1L)
  expect_equal(sum(geom_classes == "GeomPoint"), 1L)
  expect_equal(sum(geom_classes == "GeomErrorbar"), 1L)

  line_layer <- p$layers[[which(geom_classes == "GeomLine")[[1L]]]]
  ribbon_layer <- p$layers[[which(geom_classes == "GeomRibbon")[[1L]]]]
  point_layer <- p$layers[[which(geom_classes == "GeomPoint")[[1L]]]]
  errorbar_layer <- p$layers[[which(geom_classes == "GeomErrorbar")[[1L]]]]

  expect_true(all(line_layer$data$variable > 0))
  expect_true(all(ribbon_layer$data$variable > 0))
  expect_equal(point_layer$data$variable, 0)
  expect_equal(errorbar_layer$data$variable, 0)
  expect_identical(point_layer$aes_params$colour, "navy")
  expect_identical(errorbar_layer$aes_params$colour, "navy")
  expect_match(p$labels$title, "no zero indicator", fixed = TRUE)
  expect_false(grepl("spike", p$labels$title, ignore.case = TRUE))

  # Plot layers must use the same centered term predictions returned by the
  # public prediction method, including the exact-zero endpoint.
  predicted <- predict(
    fit,
    type = "terms",
    terms = "exposure",
    terms_seq = "equidistant"
  )[["exposure"]]
  predicted_zero <- predicted[predicted$variable == 0, , drop = FALSE]
  predicted_positive <- predicted[predicted$variable > 0, , drop = FALSE]

  expect_equal(point_layer$data$value, predicted_zero$value, tolerance = 1e-12)
  expect_equal(
    unname(line_layer$data$value),
    unname(predicted_positive$value[order(predicted_positive$variable)]),
    tolerance = 1e-12
  )
})
