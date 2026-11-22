test_that("rcs retains nested maximum and predicts with fitted knots", {
  x <- seq(-2, 2, length.out = 120)
  dat <- data.frame(x = x, y = x + 0.4 * x^3 + 0.05 * sin(13 * x))
  fit <- mfp2(y ~ rcs(x, df = 4), data = dat,
              select = 1, alpha = 1, verbose = FALSE)

  expect_equal(unname(fit$rcs_selected_df[["x"]]), 4L)
  expect_length(fit$rcs_specs$x$knots, 5L)
  expect_equal(length(fit$term_to_columns$x), 4L)
  expect_length(predict(fit, newdata = data.frame(x = 0.3)), 1L)
  expect_equal(as.numeric(predict(fit, newdata = data.frame(x = 0.3))),
               as.numeric(predict(fit, newdata = data.frame(x = c(0.3, 0.6))))[[1L]])
})

test_that("a forced, uncentered RCS agrees with rms at identical knots", {
  skip_if_not_installed("rms", minimum_version = "7.0-0")
  old_options <- options(rcspc = FALSE)
  on.exit(options(old_options), add = TRUE)

  x <- seq(0.2, 4, length.out = 120)
  dat <- data.frame(x = x, y = 0.3 + sin(2 * x) + 0.1 * cos(5 * x))
  knots <- c(0.3, 1.1, 2, 2.8, 3.8)
  fit <- mfp2(y ~ rcs(x, df = 4, knots = knots, center = FALSE),
              data = dat, select = 1, alpha = 1, verbose = FALSE)

  # rms::rcs() delegates to Hmisc::rcspline.eval() with the same normalized
  # truncated-cubic basis. A local formula environment prevents mfp2::rcs()
  # from masking rms::rcs() in the independent reference fit.
  rms_formula <- y ~ rcs(x, knots)
  rms_env <- new.env(parent = environment(rms_formula))
  rms_env$rcs <- rms::rcs
  environment(rms_formula) <- rms_env
  reference <- rms::ols(rms_formula, data = dat)

  rms_basis <- rms::rcs(x, knots)
  own_basis <- mfp2:::rcs_basis(x, knots)
  expect_equal(fit$rcs_selected_df[["x"]], 4L)
  expect_equal(fit$rcs_specs$x$knots, knots)
  expect_equal(dim(own_basis), dim(rms_basis))
  expect_equal(as.numeric(own_basis), as.numeric(rms_basis), tolerance = 1e-10)
  expect_equal(unname(stats::coef(fit)), unname(stats::coef(reference)),
               tolerance = 1e-7)

  query <- data.frame(x = c(0.2, 0.9, 2.3, 3.9))
  expect_equal(as.numeric(predict(fit, newdata = query)),
               as.numeric(predict(reference, newdata = query)),
               tolerance = 1e-7)
  centered <- mfp2(y ~ rcs(x, df = 4, knots = knots), data = dat,
                   select = 1, alpha = 1, verbose = FALSE)
  expect_equal(as.numeric(predict(centered, newdata = query)),
               as.numeric(predict(reference, newdata = query)),
               tolerance = 1e-7)
})

test_that("selecting a spline preserves the initial df of other FP terms", {
  x <- seq(0.1, 3, length.out = 120)
  z <- 2 + cos(5 * x)
  dat <- data.frame(x = x, z = z, y = sin(4 * x) + log(z))
  for (spline in c("rcs", "bspline")) {
    spline_call <- if (spline == "rcs") quote(rcs(x, df = 3)) else
      quote(bspline(x, df = 6))
    formula <- stats::reformulate(c(deparse(spline_call), "fp(z, df = 2)"), "y")
    fit <- mfp2(formula, data = dat, select = 1, alpha = 1,
                verbose = FALSE)
    expect_equal(as.numeric(fit$fp_terms["z", "df_setting"]), 2)
    expect_equal(as.numeric(fit$fp_terms["z", "df_initial"]), 2)
  }
})

test_that("formula spline columns inherit a TRUE default and honor individual centering", {
  x <- seq(0.1, 3, length.out = 120)
  z <- 2 + sin(5 * x)
  dat <- data.frame(x = x, z = z, y = sin(4 * x) + log(z))
  forms <- list(
    rcs = list(default = y ~ rcs(x, df = 3) + fp(z, df = 1),
               raw = y ~ rcs(x, df = 3, center = FALSE) + fp(z, df = 1)),
    bs = list(default = y ~ bspline(x, df = 6) + fp(z, df = 1),
              raw = y ~ bspline(x, df = 6, center = FALSE) + fp(z, df = 1))
  )
  for (kind in names(forms)) {
    centered <- mfp2(forms[[kind]]$default, data = dat, select = 1,
                     alpha = 1, verbose = FALSE)
    uncentered <- mfp2(forms[[kind]]$raw, data = dat, select = 1,
                       alpha = 1, verbose = FALSE)
    centered_info <- mfp2_design_column_info(centered)
    raw_info <- mfp2_design_column_info(uncentered)
    centered_info <- centered_info[centered_info$variable == "x", ]
    raw_info <- raw_info[raw_info$variable == "x", ]
    expect_true(all(centered_info$centered))
    expect_equal(centered_info$center[[1L]], mean(x), tolerance = 1e-7)
    expect_true(any(abs(centered_info$center) > 1e-7))
    expect_false(any(raw_info$centered))
    expect_true(all(raw_info$center == 0))
    query <- dat[c(5L, 42L, 103L), c("x", "z")]
    expect_equal(unname(predict(centered, newdata = query)),
                 unname(predict(uncentered, newdata = query)), tolerance = 1e-6)
  }
  expect_error(rcs(x, center = NA), "`center` must be TRUE or FALSE")
  expect_error(bspline(x, center = 1), "`center` must be TRUE or FALSE")
})

test_that("rcs compares the maximum spline with linear and removal", {
  x <- seq(-2, 2, length.out = 120)
  dat <- data.frame(x = x, y = x + 0.4 * x^3 + 0.05 * sin(13 * x))
  linear <- mfp2(y ~ rcs(x, df = 4), data = dat,
                 select = 1, alpha = 0, verbose = FALSE)
  expect_equal(unname(linear$rcs_selected_df[["x"]]), 1L)
  expect_length(linear$term_to_columns$x, 1L)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_equal(plot(linear, terms = "x", partial_only = TRUE)$x$labels$title,
                 "x: linear")
  }

  dat$y <- rep(c(-1, 1), 60)
  removed <- mfp2(y ~ rcs(x, df = 4), data = dat,
                  select = 0, verbose = FALSE)
  expect_true(is.na(removed$rcs_selected_df[["x"]]))
})

test_that("rcs resolves automatic knots on fitted subset", {
  dat <- data.frame(x = seq(-2, 2, length.out = 120), y = seq_len(120))
  fit <- mfp2(y ~ rcs(x, df = 3), data = dat,
              subset = x > 0, select = 1, alpha = 1, verbose = FALSE)
  fitted_x <- sort(dat$x[dat$x > 0])
  expected <- as.numeric(stats::quantile(
    fitted_x, probs = seq(0.05, 0.95, length.out = 4), names = FALSE
  ))
  expected[c(1, 4)] <- fitted_x[c(5, length(fitted_x) - 4)]
  expect_equal(fit$rcs_specs$x$knots, expected)
  expect_error(rcs(dat$x, df = 3, knots = 1:3), "df \\+ 1")
  expect_error(rcs(dat$x, df = 2^32), "integer `df`")
})

test_that("tied values reduce the effective maximum spline df", {
  x <- rep(0:4, each = 20)
  dat <- data.frame(x = x, y = x + 0.2 * x^2 + sin(seq_along(x)))
  expect_warning(
    fit <- mfp2(y ~ rcs(x, df = 4), data = dat,
                select = 1, alpha = 1, verbose = FALSE),
    "maximum df reduced from 4 to 2"
  )
  expect_equal(fit$rcs_specs$x$df, 2L)
  expect_equal(length(fit$rcs_specs$x$knots), 3L)
  expect_equal(fit$rcs_selected_df[["x"]], 2L)
})

test_that("rcs catzero keeps the exact-zero indicator with the selected spline", {
  x <- c(rep(0, 30), seq(0.1, 3, length.out = 90))
  dat <- data.frame(x = x, y = 2 * (x == 0) + x + 0.1 * sin(x))
  fit <- mfp2(y ~ rcs(x, df = 3, catzero = TRUE), data = dat,
              criterion = "aic", verbose = FALSE)
  expect_true(isTRUE(fit$rcs_specs$x$catzero))
  expect_true(tail(fit$term_to_columns$x, 1L) %in%
              fit$formula_term_to_columns$x)
  expect_equal(length(predict(fit, newdata = data.frame(x = c(0, 0.5)))), 2L)
})

test_that("rcs spike compares both component choices under AIC", {
  x <- c(rep(0, 30), seq(0.1, 3, length.out = 90))
  dat <- data.frame(x = x, y = 3 * (x == 0) + 0.01 * x)
  fit <- mfp2(y ~ rcs(x, df = 3, spike = TRUE), data = dat,
              criterion = "aic", verbose = FALSE)
  expect_true(isTRUE(fit$rcs_specs$x$spike))
  expect_true(fit$rcs_selected_df[["x"]] %in% 0:3)
  expect_length(predict(fit, newdata = data.frame(x = 0)), 1L)
})

test_that("rcs spike p-value selection can retain both components", {
  x <- c(rep(0, 30), seq(0.1, 3, length.out = 90))
  dat <- data.frame(x = x, y = 2 * (x == 0) + x + 0.1 * sin(x))
  fit <- mfp2(y ~ rcs(x, df = 3, spike = TRUE), data = dat,
              select = 1, alpha = 1, verbose = FALSE)
  expect_equal(fit$rcs_selected_df[["x"]], 3L)
  expect_length(fit$term_to_columns$x, 4L)
  expect_equal(length(predict(fit, newdata = data.frame(x = c(0, 0.5)))), 2L)
})

test_that("explicit catzero survives an ineligible spike request", {
  x <- c(0, seq(0.1, 3, length.out = 99))
  dat <- data.frame(x = x, y = x + as.integer(x == 0))
  expect_warning(
    fit <- mfp2(y ~ rcs(x, df = 3, spike = TRUE, catzero = TRUE),
                data = dat, criterion = "aic", verbose = FALSE),
    "spike was reset"
  )
  expect_false(fit$rcs_specs$x$spike)
  expect_true(fit$rcs_specs$x$catzero)
  expect_error(rcs(c(-1, 0, 1, 2, 3, 4), df = 2, spike = TRUE),
               "nonnegative")
})

test_that("RCS term predictions use the fitted basis on data and grid values", {
  x <- seq(-2, 2, length.out = 120)
  dat <- data.frame(x = x, y = x + 0.4 * x^3 + 0.05 * sin(13 * x))
  fit <- mfp2(y ~ rcs(x, df = 4), data = dat,
              select = 1, alpha = 1, verbose = FALSE)
  query <- data.frame(x = c(-0.4, 0.3, 1.1))
  term <- predict(fit, type = "terms", terms = "x", newdata = query,
                  terms_seq = "data", add_intercept = FALSE)$x
  full <- predict(fit, newdata = query, type = "link")
  expect_equal(term$variable, query$x)
  expect_equal(term$value, as.numeric(full) - unname(coef(fit)[["(Intercept)"]]),
               tolerance = 1e-6)
  grid <- predict(fit, type = "terms", terms = "x", nseq = 17)$x
  expect_equal(nrow(grid), 17L)
  expect_equal(range(grid$variable), range(dat$x))
  contrast <- predict(fit, type = "contrasts", terms = "x",
                      newdata = data.frame(x = 0.3), terms_seq = "data",
                      ref = list(x = 0.3))$x
  expect_equal(contrast$value, 0, tolerance = 1e-7)
})

test_that("RCS zero indicator remains separate in term prediction and plots", {
  x <- c(rep(0, 30), seq(0.1, 3, length.out = 90))
  dat <- data.frame(x = x, y = 2 * (x == 0) + x + 0.1 * sin(x))
  fit <- mfp2(y ~ rcs(x, df = 3, catzero = TRUE), data = dat,
              select = 1, alpha = 1, verbose = FALSE)
  values <- predict(fit, type = "terms", terms = "x", nseq = 11)$x
  expect_equal(values$variable[[1L]], 0)
  expect_equal(nrow(values), 12L)
  expect_true(all(values$variable[-1L] > 0))
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    graph <- plot(fit, terms = "x", partial_only = TRUE)$x
    expect_s3_class(graph, "ggplot")
    expect_equal(graph$labels$title, "x: rcs(df=3) + zero indicator")
  }
})

test_that("verbose RCS p-value progress prints comparisons within each cycle", {
  x <- seq(-2, 2, length.out = 120)
  dat <- data.frame(x = x, y = x + 0.4 * x^3)
  printed <- paste(capture.output(
    fit <- suppressMessages(mfp2(y ~ rcs(x, df = 4), data = dat,
                                 select = 1, alpha = 1, verbose = TRUE))
  ), collapse = "\n")
  expect_match(printed, "Variable: x\n", fixed = TRUE)
  expect_match(printed, "rcs\\(df=4\\)")
  expect_match(printed, "Deviance")
  expect_match(printed, "Deviance diff.", fixed = TRUE)
  expect_match(printed, "P-value")
  expect_match(printed, "Selected: rcs\\(df=4\\)")
})

test_that("RCS spike indicator-only terms use the original zero scale", {
  x <- c(rep(0, 30), seq(0.1, 3, length.out = 90))
  # A nonzero residual variance makes Gaussian AIC comparisons well defined.
  dat <- data.frame(x = x, y = 5 * (x == 0) +
                     0.1 * rep(c(-1, 1), length.out = length(x)))
  fit <- mfp2(y ~ rcs(x, df = 3, spike = TRUE), data = dat,
              criterion = "aic", verbose = FALSE)
  expect_equal(fit$rcs_selected_df[["x"]], 0L)
  values <- predict(fit, type = "terms", terms = "x", terms_seq = "data",
                    newdata = data.frame(x = c(0, 2)))$x
  expect_equal(values$variable, c(1, 0))
  contrast <- predict(fit, type = "contrasts", terms = "x",
                      terms_seq = "data", newdata = data.frame(x = c(0, 2)),
                      ref = list(x = 0))$x
  expect_equal(contrast$value[[1L]], 0, tolerance = 1e-7)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    graph <- plot(fit, terms = "x", partial_only = TRUE)$x
    expect_s3_class(graph, "ggplot")
    expect_true(any(vapply(graph$layers, function(layer) {
      inherits(layer$geom, "GeomErrorbar")
    }, logical(1L))))
  }
})

test_that("print and summary show RCS components and their fitted basis", {
  x <- c(rep(0, 30), seq(0.1, 3, length.out = 90))
  dat <- data.frame(x = x, y = 2 * (x == 0) + x^2 + 0.1 * sin(x))
  fit <- mfp2(y ~ rcs(x, df = 3, spike = TRUE), data = dat,
              select = 1, alpha = 1, verbose = FALSE)
  selected <- fit$rcs_selected_df[["x"]]
  expect_equal(selected, 3L)
  printed <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(printed, "RCS\\(df=3\\)")
  expect_match(printed, "Spike-at-Zero")

  result <- summary(fit, basis = TRUE, formulas = TRUE)
  expect_match(result$function_table$Function[result$function_table$Variable == "x"],
               "RCS\\(df=3\\)")
  expect_match(result$nonlinear_terms$form[result$nonlinear_terms$variable == "x"],
               "RCS\\(df=3\\)")
  expect_true(any(grepl("RCS_1\\(x; knots", result$basis$term)))
  expect_true(any(grepl("I\\(x = 0\\)", result$basis$term)))
  expect_true(any(grepl("RCS_1\\(x; knots", result$formulas)))
})

test_that("summary calls an indicator-only RCS spike binary", {
  x <- c(rep(0, 30), seq(0.1, 3, length.out = 90))
  fit <- mfp2(y ~ rcs(x, df = 3, spike = TRUE),
              data = data.frame(x = x, y = 5 * (x == 0) +
                                  0.1 * rep(c(-1, 1), length.out = length(x))),
              criterion = "aic", verbose = FALSE)
  expect_equal(fit$rcs_selected_df[["x"]], 0L)
  result <- summary(fit)
  expect_equal(result$function_table$Function[result$function_table$Variable == "x"],
               "zero indicator only")
  expect_true(any(result$linear_terms$variable == "I(x = 0)"))
  printed <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(printed, "zero indicator only")
})
