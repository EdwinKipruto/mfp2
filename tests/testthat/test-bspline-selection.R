test_that("cubic B-spline prefixes span nested fixed-knot spaces", {
  x <- seq(0, 1, length.out = 121)
  knots <- c(0.25, 0.5, 0.75)
  basis <- mfp2:::bspline_basis(x, knots, c(0, 1))
  one <- splines::bs(x, knots = 0.5, Boundary.knots = c(0, 1))
  three <- splines::bs(x, knots = knots, Boundary.knots = c(0, 1))

  expect_equal(dim(basis), c(length(x), 6L))
  expect_equal(qr(cbind(1, basis[, 1:4]))$rank, 5L)
  expect_equal(qr(cbind(1, basis))$rank, 7L)
  expect_lt(max(abs(qr.resid(qr(cbind(1, basis[, 1:4])), one))), 1e-8)
  expect_lt(max(abs(qr.resid(qr(cbind(1, basis)), three))), 1e-8)
})

test_that("forced cubic B-spline agrees with splines::bs at identical knots", {
  x <- seq(0, 1, length.out = 121)
  dat <- data.frame(x = x, y = sin(8 * x) + 0.05 * cos(27 * x))
  knots <- c(0.25, 0.5, 0.75)
  boundaries <- c(0, 1)
  fit <- mfp2(y ~ bs(x, df = 6, knots = knots,
                     Boundary.knots = boundaries, center = FALSE),
              data = dat, select = 1, alpha = 1, verbose = FALSE)
  reference <- stats::lm(
    y ~ splines::bs(x, knots = knots, Boundary.knots = boundaries,
                    degree = 3L, intercept = FALSE),
    data = dat
  )

  # mfp2 uses nested basis columns, so coefficients need not match the
  # standard B-spline coefficients even though both span the same model.
  own_basis <- mfp2:::bspline_basis(x, knots, boundaries)
  reference_basis <- splines::bs(x, knots = knots,
                                  Boundary.knots = boundaries,
                                  degree = 3L, intercept = FALSE)
  expect_equal(fit$bspline_selected_df[["x"]], 6L)
  expect_equal(fit$bspline_specs$x$knots, knots)
  expect_equal(fit$bspline_specs$x$boundary_knots, boundaries)
  expect_equal(qr(cbind(1, own_basis))$rank, 7L)
  expect_equal(qr(cbind(1, reference_basis))$rank, 7L)
  expect_lt(max(abs(qr.resid(qr(cbind(1, own_basis)), reference_basis))), 1e-8)
  expect_lt(max(abs(qr.resid(qr(cbind(1, reference_basis)), own_basis))), 1e-8)

  query <- data.frame(x = c(0, 0.08, 0.37, 0.68, 1))
  expect_equal(as.numeric(predict(fit, newdata = query)),
               as.numeric(predict(reference, newdata = query)),
               tolerance = 1e-7)
  expect_equal(as.numeric(stats::deviance(fit)),
               as.numeric(stats::deviance(reference)), tolerance = 1e-7)
  centered <- mfp2(y ~ bs(x, df = 6, knots = knots,
                          Boundary.knots = boundaries),
                   data = dat, select = 1, alpha = 1, verbose = FALSE)
  expect_equal(as.numeric(predict(centered, newdata = query)),
               as.numeric(predict(reference, newdata = query)),
               tolerance = 1e-7)
})

test_that("bs selects the same nested models as bspline in mixed FP formulas", {
  x <- seq(0, 1, length.out = 120)
  z <- 2 + cos(5 * x)
  dat <- data.frame(x = x, z = z, y = sin(8 * x) + 0.1 * z)
  fit <- mfp2(y ~ fp(z, df = 1) + bs(x), data = dat,
              select = 1, alpha = 1, verbose = FALSE)
  legacy <- mfp2(y ~ fp(z, df = 1) + bspline(x), data = dat,
                 select = 1, alpha = 1, verbose = FALSE)
  expect_equal(fit$bspline_selected_df[["x"]], 6L)
  expect_equal(fit$bspline_specs$x$knots, legacy$bspline_specs$x$knots)
  expect_equal(fit$bspline_specs$x$boundary_knots,
               legacy$bspline_specs$x$boundary_knots)
  expect_match(summary(fit)$function_table$Function[
    summary(fit)$function_table$Variable == "x"], "B-spline")
  expect_match(paste(capture.output(print(fit)), collapse = "\n"), "B-spline")
  query <- dat[c(5L, 48L, 105L), c("x", "z")]
  expect_equal(unname(predict(fit, newdata = query)),
               unname(predict(legacy, newdata = query)), tolerance = 1e-7)
  expect_equal(fit$fp_terms["z", "df_initial"],
               legacy$fp_terms["z", "df_initial"])
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_equal(plot(fit, terms = "x", partial_only = TRUE)$x$labels$title,
                 "x: bs(df=6)")
  }

  fixed <- mfp2(y ~ mfp2::bs(x, df = 6, knots = c(0.25, 0.5, 0.75),
                             Boundary.knots = c(0, 1), center = FALSE),
                data = dat, select = 1, alpha = 1, verbose = FALSE)
  expect_equal(fixed$bspline_specs$x$knots, c(0.25, 0.5, 0.75))
  expect_equal(fixed$bspline_specs$x$boundary_knots, c(0, 1))
  expect_false(any(mfp2_design_column_info(fixed)$centered))
  expect_length(predict(fixed, newdata = data.frame(x = 0.4)), 1L)
  expect_null(formula_transform_source_name("splines::bs(x, df = 6)"))

  masked <- y ~ bs(x)
  masked_env <- new.env(parent = environment(masked))
  masked_env$bs <- function(...) stop("masked bs was evaluated")
  environment(masked) <- masked_env
  masked_fit <- mfp2(masked, data = dat, select = 1, alpha = 1,
                     verbose = FALSE)
  expect_equal(masked_fit$bspline_selected_df[["x"]], 6L)
  expect_length(predict(masked_fit, newdata = data.frame(x = 0.4)), 1L)
})

test_that("bspline selects the nested six-column form and reuses fitted knots", {
  x <- seq(0, 1, length.out = 120)
  dat <- data.frame(x = x, y = sin(8 * x) + 0.05 * cos(27 * x))
  fit <- mfp2(y ~ bspline(x, df = 6), data = dat,
              select = 1, alpha = 1, verbose = FALSE)
  spec <- fit$bspline_specs$x

  expect_equal(fit$bspline_selected_df[["x"]], 6L)
  expect_equal(spec$knots,
               as.numeric(stats::quantile(x, c(0.25, 0.5, 0.75))))
  expect_equal(spec$boundary_knots, range(x))
  expect_equal(spec$candidate_df, c(4L, 6L))
  expect_length(fit$term_to_columns$x, 6L)
  expect_equal(length(predict(fit, newdata = data.frame(x = c(0.2, 0.8)))), 2L)
  expect_equal(as.numeric(predict(fit, newdata = data.frame(x = 0.2))),
               as.numeric(predict(fit, newdata = data.frame(x = c(0.2, 0.8))))[[1L]])
  expect_match(paste(capture.output(print(fit)), collapse = "\n"), "B-spline")
  expect_match(summary(fit)$function_table$Function[
    summary(fit)$function_table$Variable == "x"], "B-spline")
  expect_true(any(grepl("BS_1\\(x; knots", summary(fit, basis = TRUE)$basis$term)))
  term <- predict(fit, type = "terms", terms = "x", nseq = 15)$x
  expect_equal(nrow(term), 15L)
  contrast <- predict(fit, type = "contrasts", terms = "x",
                      newdata = data.frame(x = 0.4), terms_seq = "data",
                      ref = list(x = 0.4))$x
  expect_equal(contrast$value, 0, tolerance = 1e-7)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    graph <- plot(fit, terms = "x", partial_only = TRUE)$x
    expect_s3_class(graph, "ggplot")
    expect_equal(graph$labels$title, "x: bs(df=6)")
  }
})

test_that("bspline compares exclusion and linearity without a knotless cubic", {
  x <- seq(0, 1, length.out = 120)
  dat <- data.frame(x = x, y = x + 0.1 * sin(5 * x))
  linear <- mfp2(y ~ bspline(x, df = 6), data = dat,
                 select = 1, alpha = 0, verbose = FALSE)
  expect_equal(linear$bspline_selected_df[["x"]], 1L)
  expect_length(linear$term_to_columns$x, 1L)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_equal(plot(linear, terms = "x", partial_only = TRUE)$x$labels$title,
                 "x: linear")
  }

  dat$y <- rep(c(-1, 1), 60)
  removed <- mfp2(y ~ bspline(x, df = 6), data = dat,
                  select = 0, verbose = FALSE)
  expect_true(is.na(removed$bspline_selected_df[["x"]]))
})

test_that("automatic B-spline knots come from the fitted subset", {
  x <- seq(-2, 2, length.out = 120)
  dat <- data.frame(x = x, y = sin(3 * x))
  fit <- mfp2(y ~ bspline(x, df = 6), data = dat,
              subset = x > 0, select = 1, alpha = 1, verbose = FALSE)
  fitted_x <- dat$x[dat$x > 0]
  expect_equal(fit$bspline_specs$x$knots,
               as.numeric(stats::quantile(fitted_x, c(0.25, 0.5, 0.75))))
  expect_equal(fit$bspline_specs$x$boundary_knots, range(fitted_x))
})

test_that("bspline retains a separate zero indicator and validates knots", {
  x <- c(rep(0, 30), seq(0.1, 2, length.out = 90))
  dat <- data.frame(x = x, y = 2 * (x == 0) + sin(6 * x))
  fit <- mfp2(y ~ bspline(x, df = 6, catzero = TRUE), data = dat,
              select = 1, alpha = 1, verbose = FALSE)
  expect_true(isTRUE(fit$bspline_specs$x$indicator))
  expect_true(all(fit$bspline_specs$x$knots > 0))
  expect_equal(fit$bspline_specs$x$boundary_knots[[1L]], 0)
  expect_length(predict(fit, newdata = data.frame(x = c(0, 0.5))), 2L)

  expect_error(bspline(x, df = 6, knots = c(0.2, 0.2, 0.8)),
               "increasing interior knots")
  expect_error(bspline(c(-1, x), spike = TRUE), "nonnegative")
})

test_that("B-spline spike can retain only the zero indicator", {
  x <- c(rep(0, 30), seq(0.1, 2, length.out = 90))
  fit <- mfp2(y ~ bspline(x, df = 6, spike = TRUE),
              data = data.frame(x = x, y = 5 * (x == 0) +
                                  0.1 * rep(c(-1, 1), length.out = length(x))),
              criterion = "aic", verbose = FALSE)
  expect_equal(fit$bspline_selected_df[["x"]], 0L)
  expect_equal(length(fit$term_to_columns$x), 1L)
  expect_equal(length(predict(fit, newdata = data.frame(x = c(0, 1)))), 2L)
})

test_that("verbose B-spline AIC progress shows ordered candidates and total DF", {
  x <- seq(0, 1, length.out = 120)
  dat <- data.frame(x = x, y = sin(8 * x))
  printed <- paste(capture.output(
    fit <- suppressMessages(mfp2(y ~ bspline(x, df = 6), data = dat,
                                 criterion = "aic", verbose = TRUE))
  ), collapse = "\n")
  expect_match(printed, "Variable: x\n", fixed = TRUE)
  expect_match(printed, "DF")
  expect_match(printed, "Deviance")
  expect_match(printed, "AIC")
  expect_lt(regexpr("bs\\(df=6\\)", printed)[[1L]],
            regexpr("bs\\(df=4\\)", printed)[[1L]])
  expect_lt(regexpr("bs\\(df=4\\)", printed)[[1L]],
            regexpr("null", printed)[[1L]])
  expect_match(printed, "Selected: ")
})

test_that("spline progress uses the requested descending table layout", {
  fit <- list(
    candidate_keys = c("null", "6C", "1C", "4C"),
    metrics = cbind(df = c(4, 10, 5, 8),
                    deviance_rs = c(248.481, 200.036, 223.401, 202.094)),
    pvalue = c("maximum vs null" = 0.00001,
               "maximum vs linear" = 0.00002,
               "maximum vs df 4" = 0.393),
    stage1_key = "4C", model_best = "4C"
  )
  before <- fit
  output <- capture.output(mfp2:::print_spline_step(
    "cavol", "pvalue", fit,
    list(type = "bspline", df = 6L, indicator = FALSE), "gaussian"
  ))
  expect_true("Variable: cavol" %in% output)
  expect_match(output[grepl("DF", output)],
               "DF +Deviance +Versus +Deviance diff\\. +P-value")
  rows <- output[grepl("^(bs\\(df=|linear +|null +)", output)]
  fields <- strsplit(trimws(rows), "[[:space:]]+")
  expect_identical(vapply(fields, `[[`, character(1L), 1L),
                   c("bs(df=6)", "bs(df=4)", "linear", "null"))
  expect_identical(fields[[1L]], c("bs(df=6)", "10", "200.036", ".", ".", "."))
  expect_identical(fields[[2L]], c("bs(df=4)", "8", "202.094", "bs(df=6)",
                                  "2.058", "0.393"))
  expect_true("Selected: bs(df=4)" %in% output)
  expect_false(any(grepl("Model|excluded|\\[B-spline\\]|–", output)))
  expect_identical(fit, before)
})
