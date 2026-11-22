test_that("matrix spline specifications retain fitted knots and reuse them for prediction", {
  age <- seq(0, 1, length.out = 121)
  z <- 2 + cos(6 * pi * age)
  x <- cbind(age = age, z = z)
  y <- sin(8 * age) + 0.2 * z
  fit <- mfp2(x, y, splines = c(age = "bspline"),
              df = c(age = 6, z = 4), select = 1, alpha = 1,
              verbose = FALSE)
  expect_equal(fit$bspline_selected_df[["age"]], 6L)
  expect_equal(fit$bspline_specs$age$knots,
               as.numeric(stats::quantile(age, c(0.25, 0.5, 0.75))))
  centered_basis <- mfp2_design_column_info(fit)
  centered_basis <- centered_basis[centered_basis$variable == "age", ]
  expect_true(all(centered_basis$centered))
  expect_equal(centered_basis$center[[1L]], mean(age), tolerance = 1e-7)
  raw_fit <- mfp2(x, y, splines = c(age = "bspline"),
                  df = c(age = 6, z = 4), center = c(age = FALSE, z = TRUE),
                  select = 1, alpha = 1, verbose = FALSE)
  raw_basis <- mfp2_design_column_info(raw_fit)
  raw_basis <- raw_basis[raw_basis$variable == "age", ]
  expect_false(any(raw_basis$centered))
  expect_true(all(raw_basis$center == 0))
  expect_equal(unname(predict(fit, newdata = x[1:3, , drop = FALSE])),
               unname(predict(raw_fit, newdata = x[1:3, , drop = FALSE])),
               tolerance = 1e-6)
  query <- cbind(age = c(0.25, 0.8), z = c(2.1, 1.8))
  expect_length(predict(fit, newdata = query), 2L)
  expect_equal(unname(predict(fit, newdata = query[1L, , drop = FALSE])),
               unname(predict(fit, newdata = query))[[1L]], tolerance = 1e-7)
  expect_equal(nrow(predict(fit, type = "terms", terms = "age", nseq = 17)$age),
               17L)
  expect_match(paste(capture.output(print(fit)), collapse = "\n"), "B-spline")
  expect_match(summary(fit)$function_table$Function[
    summary(fit)$function_table$Variable == "age"], "B-spline")
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_equal(plot(fit, terms = "age", partial_only = TRUE)$age$labels$title,
                 "age: bs(df=6)")
  }
})

test_that("matrix spline constructors accept fixed knots and coexist with FP", {
  age <- seq(0, 1, length.out = 120)
  bmi <- 1 + 2 * age + 0.15 * sin(8 * pi * age)
  z <- 2 + 2 * age + 0.2 * cos(7 * pi * age)
  x <- cbind(age = age, bmi = bmi, z = z)
  y <- sin(8 * age) + bmi^2 + log(z)
  fit <- mfp2(x, y, splines = list(
    age = bs_spec(knots = c(0.25, 0.5, 0.75), boundary_knots = c(0, 1)),
    bmi = rcs_spec(knots = c(1, 1.5, 2.5, 3))
  ), df = c(z = 4), select = 1, alpha = 1, verbose = FALSE)
  expect_equal(fit$bspline_specs$age$df, 6L)
  expect_equal(fit$rcs_specs$bmi$df, 3L)
  expect_equal(fit$bspline_specs$age$knots, c(0.25, 0.5, 0.75))
  expect_equal(fit$rcs_specs$bmi$knots, c(1, 1.5, 2.5, 3))
  rcs_basis_info <- mfp2_design_column_info(fit)
  rcs_basis_info <- rcs_basis_info[rcs_basis_info$variable == "bmi", ]
  expect_true(all(rcs_basis_info$centered))
  expect_equal(length(predict(fit, newdata = x[1:3, , drop = FALSE])), 3L)
  expect_equal(nrow(predict(fit, type = "terms", terms = "bmi", nseq = 11)$bmi), 11L)
  expect_match(summary(fit)$function_table$Function[
    summary(fit)$function_table$Variable == "bmi"], "RCS")
  contrast <- predict(fit, type = "contrasts", terms = "age",
                      newdata = x[1:2, , drop = FALSE], terms_seq = "data",
                      ref = list(age = 0.5))$age
  expect_equal(nrow(contrast), 2L)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_equal(plot(fit, terms = "bmi", partial_only = TRUE)$bmi$labels$title,
                 "bmi: rcs(df=3)")
  }
  expect_error(mfp2(x, y, splines = list(age = bs_spec(df = 5,
    knots = c(0.25, 0.5, 0.75))), verbose = FALSE), "knot count")
})

test_that("matrix spline knots use fitted subset and existing zero options", {
  age <- c(rep(0, 30), seq(0.1, 3, length.out = 90), 30)
  x <- cbind(age = age)
  y <- sin(age) + 3 * (age == 0)
  fit <- mfp2(x, y, subset = seq_len(120),
              splines = list(age = rcs_spec(df = 3)),
              catzero_vars = "age", select = 1, alpha = 1, verbose = FALSE)
  expect_true(fit$rcs_specs$age$indicator)
  expect_equal(fit$rcs_specs$age$knots,
               mfp2:::rcs_default_knots(age[31:120], 3L))
  expect_length(predict(fit, newdata = cbind(age = c(0, 1))), 2L)
  expect_error(predict(fit, newdata = cbind(age = -1)), "nonnegative")

  spike <- mfp2(x[1:120, , drop = FALSE], y[1:120],
                splines = c(age = "rcs"), spike_vars = "age",
                criterion = "aic", verbose = FALSE)
  expect_true(spike$rcs_specs$age$spike)
  expect_true(spike$rcs_specs$age$indicator)
})

test_that("matrix splines can share a model with a grouped linear term", {
  age <- seq(0, 1, length.out = 120)
  group <- rep(1:3, length.out = length(age))
  x <- cbind(age = age, raceB = as.integer(group == 2L),
             raceC = as.integer(group == 3L))
  y <- sin(8 * age) + (group == 2L) - (group == 3L)
  fit <- mfp2(x, y, splines = c(age = "bspline"),
              term_groups = list(race = c("raceB", "raceC")),
              select = 1, alpha = 1, verbose = FALSE)
  expect_equal(fit$bspline_specs$age$df, 6L)
  expect_length(fit$term_to_columns$race, 2L)
  expect_length(predict(fit, newdata = x[1:2, , drop = FALSE]), 2L)
})
