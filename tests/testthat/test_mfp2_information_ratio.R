# Regression tests for the low-information diagnostic. These tests deliberately
# separate ratio arithmetic from full model fitting so failures identify whether
# the complexity calculation or the public-interface behavior regressed.

test_that("small matrix and formula fits warn instead of failing at five rows", {
  dat <- data.frame(
    y = c(1.0, 2.2, 2.7, 4.1),
    x = c(1.0, 2.0, 3.0, 4.0)
  )

  expect_warning(
    matrix_fit <- mfp2(
      x = as.matrix(dat["x"]),
      y = dat$y,
      family = "gaussian",
      df = 1,
      xorder = "original",
      verbose = FALSE
    ),
    "4.00 observations per initial model degree of freedom"
  )
  expect_s3_class(matrix_fit, "mfp2")

  # The formula method delegates the same four fitting rows to mfp2.default(),
  # so it must produce the same warning rather than its former absolute stop.
  expect_warning(
    formula_fit <- mfp2(
      y ~ fp(x, df = 1),
      data = dat,
      family = "gaussian",
      xorder = "original",
      verbose = FALSE
    ),
    "4.00 observations per initial model degree of freedom"
  )
  expect_s3_class(formula_fit, "mfp2")
})


test_that("public interfaces can disable only the low-information diagnostic", {
  # The diagnostic remains enabled by default for both public dispatch paths.
  expect_identical(formals(mfp2.default)$warn_low_information, TRUE)
  expect_identical(formals(mfp2.formula)$warn_low_information, TRUE)

  dat <- data.frame(
    y = c(1.0, 2.2, 2.7, 4.1),
    x = c(1.0, 2.0, 3.0, 4.0)
  )

  # Both interfaces would warn at 4 / 1 = 4 by default. The explicit FALSE
  # must reach fit_mfp() and suppress that diagnostic without changing the fit.
  expect_no_warning(
    matrix_fit <- mfp2(
      x = as.matrix(dat["x"]),
      y = dat$y,
      family = "gaussian",
      df = 1,
      xorder = "original",
      verbose = FALSE,
      warn_low_information = FALSE
    )
  )
  expect_s3_class(matrix_fit, "mfp2")

  expect_no_warning(
    formula_fit <- mfp2(
      y ~ fp(x, df = 1),
      data = dat,
      family = "gaussian",
      xorder = "original",
      verbose = FALSE,
      warn_low_information = FALSE
    )
  )
  expect_s3_class(formula_fit, "mfp2")
})


test_that("warn_low_information must be one non-missing logical value", {
  x <- matrix(seq_len(6), ncol = 1L, dimnames = list(NULL, "x"))
  y <- seq_len(6)

  expect_error(
    mfp2(x, y, df = 1, warn_low_information = NA, verbose = FALSE),
    "`warn_low_information`"
  )
  expect_error(
    mfp2(x, y, df = 1, warn_low_information = c(TRUE, FALSE), verbose = FALSE),
    "`warn_low_information`"
  )
  expect_error(
    mfp2(x, y, df = 1, warn_low_information = 1, verbose = FALSE),
    "`warn_low_information`"
  )
})


test_that("the warning uses effective term df rather than predictor count", {
  x <- matrix(seq_len(20), ncol = 1L, dimnames = list(NULL, "x"))

  # One FP2 term has four initial MFP df. Twenty observations therefore give
  # exactly five information units per df and should meet the warning boundary.
  expect_warning(
    info <- mfp2:::warn_mfp_information_ratio(
      x = x,
      y = seq_len(20),
      weights = rep(1, 20),
      family_string = "gaussian",
      df = c(x = 4),
      term_to_columns = list(x = "x"),
      catzero = c(x = FALSE),
      threshold = 5
    ),
    "5.00 observations per initial model degree of freedom"
  )

  expect_equal(info$information, 20)
  expect_equal(info$initial_df, 4)
  expect_equal(info$ratio, 5)
})


test_that("catzero indicators and grouped terms contribute their fitted df", {
  x <- cbind(
    biomarker = seq_len(30),
    group_b = rep(c(0, 1), 15),
    group_c = rep(c(0, 0, 1), 10)
  )

  # biomarker contributes FP2 df = 4 plus one catzero indicator (5 total);
  # the grouped factor contributes its two design columns despite its
  # fixed-linear search setting of df = 1. Total initial df is therefore 7.
  expect_warning(
    info <- mfp2:::warn_mfp_information_ratio(
      x = x,
      y = seq_len(30),
      weights = rep(1, 30),
      family_string = "gaussian",
      df = c(biomarker = 4, group = 1),
      term_to_columns = list(
        biomarker = "biomarker",
        group = c("group_b", "group_c")
      ),
      catzero = c(biomarker = TRUE, group = FALSE),
      threshold = 5
    ),
    "4.29 observations per initial model degree of freedom"
  )

  expect_equal(info$initial_df, 7)
  expect_equal(info$ratio, 30 / 7)
})


test_that("Cox and binomial fits use outcome-specific information units", {
  cox_x <- matrix(
    seq_len(600),
    ncol = 3L,
    dimnames = list(NULL, c("x1", "x2", "x3"))
  )
  cox_y <- cbind(time = seq_len(200), status = c(rep(1, 30), rep(0, 170)))

  # Two FP2 terms and one linear term contribute 9 df. Cox adequacy is based on
  # the 30 observed events, not on all 200 participant rows.
  expect_warning(
    cox_info <- mfp2:::warn_mfp_information_ratio(
      x = cox_x,
      y = cox_y,
      weights = rep(1, 200),
      family_string = "cox",
      df = c(x1 = 4, x2 = 4, x3 = 1),
      term_to_columns = list(x1 = "x1", x2 = "x2", x3 = "x3"),
      catzero = c(x1 = FALSE, x2 = FALSE, x3 = FALSE),
      threshold = 5
    ),
    "3.33 events per initial model degree of freedom"
  )
  expect_equal(cox_info$information, 30)
  expect_equal(cox_info$ratio, 30 / 9)

  # Binomial information is the weighted minority-outcome total. This catches
  # severe imbalance that an ordinary row-count ratio would conceal.
  expect_warning(
    binomial_info <- mfp2:::warn_mfp_information_ratio(
      x = matrix(seq_len(20), ncol = 1L, dimnames = list(NULL, "x")),
      y = c(rep(1, 3), rep(0, 17)),
      weights = rep(1, 20),
      family_string = "binomial",
      df = c(x = 1),
      term_to_columns = list(x = "x"),
      catzero = c(x = FALSE),
      threshold = 5
    ),
    "3.00 minority outcome units per initial model degree of freedom"
  )
  expect_equal(binomial_info$information, 3)
})


test_that("ratios above the threshold remain quiet", {
  expect_no_warning(
    info <- mfp2:::warn_mfp_information_ratio(
      x = matrix(seq_len(6), ncol = 1L, dimnames = list(NULL, "x")),
      y = seq_len(6),
      weights = rep(1, 6),
      family_string = "gaussian",
      df = c(x = 1),
      term_to_columns = list(x = "x"),
      catzero = c(x = FALSE),
      threshold = 5
    )
  )
  expect_equal(info$ratio, 6)
})
