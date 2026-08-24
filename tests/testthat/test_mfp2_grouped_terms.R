# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# -----------------------------------------------------------------------------
# 1.1 Grouped terms in the default matrix interface
# -----------------------------------------------------------------------------

# Test purpose: Checks that manually supplied dummy columns can be represented
# and fitted as one conceptual linear term.
test_that("mfp2.default() groups manually supplied dummy columns", {
  set.seed(2101)
  n <- 180

  group <- factor(rep(c("A", "B", "C"), length.out = n))
  group_mm <- stats::model.matrix(~ group)[, -1L, drop = FALSE]
  x <- cbind(
    x1 = runif(n, 1, 10),
    group_mm
  )
  y <- 0.4 * x[, "x1"] +
    0.8 * x[, "groupB"] -
    0.5 * x[, "groupC"] +
    rnorm(n, sd = 0.3)

  fit <- mfp2(
    x,
    y,
    df = 1,
    select = 1,
    term_groups = list(group = c("groupB", "groupC")),
    verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_identical(
    fit$term_to_columns[["group"]],
    c("groupB", "groupC")
  )
  expect_true(fit$fp_terms["group", "selected"])
  expect_true(all(c("groupB", "groupC") %in% colnames(fit$x_original)))
})


# Test purpose: Checks that selection operates on the complete grouped block;
# the final fitting matrix must contain either every member column or none.
test_that("grouped matrix columns are retained or removed together", {
  set.seed(2102)
  n <- 180

  group <- factor(rep(c("A", "B", "C"), length.out = n))
  group_mm <- stats::model.matrix(~ group)[, -1L, drop = FALSE]
  x <- cbind(
    x1 = runif(n, 1, 10),
    group_mm
  )
  y <- 0.5 * x[, "x1"] + rnorm(n)

  fit <- mfp2(
    x,
    y,
    term_groups = list(group = c("groupB", "groupC")),
    verbose = FALSE
  )

  present <- c("groupB", "groupC") %in% colnames(fit$x_original)
  expect_true(all(present) || !any(present))
  expect_identical(all(present), isTRUE(fit$fp_terms["group", "selected"]))
})


# Test purpose: A scalar df is a continuous-variable default and must not
# send a supplied ordered-factor contrast block into the FP search.
test_that("grouped ordered-factor columns are linear under scalar df default", {
  set.seed(21021)
  n <- 160
  stage <- ordered(rep(LETTERS[1:4], length.out = n))
  stage_mm <- stats::model.matrix(~ stage)[, -1L, drop = FALSE]
  x <- cbind(
    x1 = seq(0.5, 8, length.out = n),
    stage_mm
  )
  y <- 0.4 * x[, "x1"] + 0.8 * stage_mm[, 1L] -
    0.5 * stage_mm[, 2L] + stats::rnorm(n, sd = 0.3)

  fit <- mfp2(
    x,
    y,
    term_groups = list(stage = colnames(stage_mm)),
    keep = "stage",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_equal(as.numeric(fit$fp_terms["stage", "df_setting"]), 1)
  expect_equal(
    as.numeric(fit$fp_terms["stage", "df_initial"]),
    ncol(stage_mm)
  )
  expect_true(fit$fp_terms["stage", "selected"])
})


# Test purpose: Automatic scaling is a continuous-variable default and must not
# independently rescale the member columns of a supplied categorical block.
test_that("grouped matrix columns default to scale one", {
  set.seed(21022)
  n <- 160
  stage <- ordered(rep(LETTERS[1:4], length.out = n))
  stage_mm <- stats::model.matrix(~ stage)[, -1L, drop = FALSE]
  stage_mm[, 1L] <- 1000 * stage_mm[, 1L]
  x <- cbind(
    x1 = seq(0.5, 8, length.out = n),
    stage_mm
  )
  y <- 0.4 * x[, "x1"] + 0.001 * stage_mm[, 1L] -
    0.5 * stage_mm[, 2L] + stats::rnorm(n, sd = 0.3)

  fit <- mfp2(
    x,
    y,
    term_groups = list(stage = colnames(stage_mm)),
    keep = "stage",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_equal(as.numeric(fit$transformations["stage", "scale"]), 1)
  expect_true(fit$fp_terms["stage", "selected"])
})


# Test purpose: A retained manually grouped block reports its fitted rank
# contribution rather than the single linear power used by the selection engine.
test_that("grouped matrix terms report one final df per estimable coefficient", {
  set.seed(2103)
  n <- 180
  group <- factor(rep(c("A", "B", "C"), length.out = n))
  x <- stats::model.matrix(~ group)[, -1L, drop = FALSE]
  y <- 0.9 * x[, "groupB"] - 0.6 * x[, "groupC"] +
    stats::rnorm(n, sd = 0.2)

  fit <- mfp2(
    x,
    y,
    term_groups = list(group = c("groupB", "groupC")),
    keep = "group",
    select = 0,
    df = 1,
    verbose = FALSE
  )

  transformed_columns <- paste0(fit$term_to_columns[["group"]], ".1")
  coefficient_columns <- unname(
    fit$transformed_to_model_columns[transformed_columns]
  )
  expected_df <- sum(!is.na(stats::coef(fit)[coefficient_columns]))

  expect_equal(expected_df, 2)
  expect_equal(as.numeric(fit$fp_terms["group", "df_final"]), expected_df)
})


# Test purpose: The shared metadata expansion applies grouped-term invariants
# while preserving all fitted settings for ordinary singleton terms.
test_that("term metadata expansion is consistent for grouped and singleton terms", {
  acd_fit <- list(beta0 = 0.2, beta1 = 0.8, power = 1)
  expanded <- expand_term_metadata_to_columns(
    term_to_columns = list(
      x1 = "x1",
      group = c("groupB", "groupC"),
      omitted = c("omittedB", "omittedC")
    ),
    powers = list(x1 = c(1, 2), group = 1, omitted = NA_real_),
    raw_columns = c("groupC", "x1", "groupB", "omittedB"),
    center = c(x1 = TRUE, group = TRUE, omitted = FALSE),
    acdx = c(x1 = TRUE, group = TRUE, omitted = TRUE),
    zero = c(x1 = TRUE, group = TRUE, omitted = TRUE),
    catzero = c(x1 = TRUE, group = TRUE, omitted = TRUE),
    spike = c(x1 = TRUE, group = TRUE, omitted = TRUE),
    spike_decision = c(x1 = 1L, group = 1L, omitted = 1L),
    acd_parameter = list(x1 = acd_fit, group = acd_fit, omitted = acd_fit)
  )

  expect_identical(
    unname(expanded$terms),
    c("group", "x1", "group", "omitted")
  )
  expect_identical(expanded$powers[["groupC"]], 1)
  expect_identical(expanded$powers[["groupB"]], 1)
  expect_true(is.na(expanded$powers[["omittedB"]]))
  expect_identical(expanded$powers[["x1"]], c(1, 2))
  expect_true(expanded$center[["groupC"]])
  expect_true(expanded$center[["x1"]])

  grouped_columns <- c("groupC", "groupB", "omittedB")
  expect_false(any(expanded$acdx[grouped_columns]))
  expect_false(any(expanded$zero[grouped_columns]))
  expect_false(any(expanded$catzero[grouped_columns]))
  expect_false(any(expanded$spike[grouped_columns]))
  expect_true(all(
    expanded$spike_decision[grouped_columns] ==
      saz_decision_codes[["continuous_only"]]
  ))
  expect_true(all(vapply(
    expanded$acd_parameter[grouped_columns],
    is.null,
    logical(1L)
  )))

  expect_true(expanded$acdx[["x1"]])
  expect_true(expanded$zero[["x1"]])
  expect_true(expanded$catzero[["x1"]])
  expect_true(expanded$spike[["x1"]])
  expect_identical(expanded$spike_decision[["x1"]], 1L)
  expect_identical(expanded$acd_parameter[["x1"]], acd_fit)
})


# Test purpose: Grouped-term df uses fitted estimability, so aliased member
# coefficients do not count toward the final rank contribution.
test_that("grouped final df excludes non-estimable coefficients", {
  fp_terms <- create_fp_terms(
    fp_powers = list(group = 1),
    acdx = c(group = FALSE),
    df = c(group = 1),
    select = c(group = 1),
    alpha = c(group = 1),
    criterion = "pvalue",
    zero = c(group = FALSE),
    catzero = c(group = FALSE),
    spike = c(group = FALSE),
    spike_decision = c(group = 2),
    term_to_columns = list(group = c("groupB", "groupC")),
    transformed_to_model_columns = c(
      "groupB.1" = "groupB.1",
      "groupC.1" = "groupC.1"
    ),
    coefficients = c(
      "(Intercept)" = 0,
      "groupB.1" = 0.4,
      "groupC.1" = NA_real_
    )
  )

  expect_equal(as.numeric(fp_terms["group", "df_setting"]), 1)
  expect_equal(as.numeric(fp_terms["group", "df_initial"]), 2)
  expect_equal(as.numeric(fp_terms["group", "df_final"]), 1)
})


# Test purpose: Eligible SAZ terms enter selection with both the continuous
# component and the structural-zero indicator. Initial df must therefore count
# the binary indicator even when Stage 2 subsequently removes it.
test_that("eligible SAZ initial df includes the binary indicator", {
  fp_terms <- create_fp_terms(
    fp_powers = list(
      saz_full = c(1, 3),
      saz_continuous = c(1, 3),
      ordinary = c(1, 3)
    ),
    acdx = c(
      saz_full = FALSE,
      saz_continuous = FALSE,
      ordinary = FALSE
    ),
    df = c(saz_full = 4, saz_continuous = 4, ordinary = 4),
    select = c(saz_full = 1, saz_continuous = 1, ordinary = 1),
    alpha = c(saz_full = 1, saz_continuous = 1, ordinary = 1),
    criterion = "pvalue",
    zero = c(saz_full = TRUE, saz_continuous = TRUE, ordinary = FALSE),
    # catzero is final/effective metadata: Stage 2 removed the binary indicator
    # from saz_continuous, but its initial SAZ candidate still contained it.
    catzero = c(saz_full = TRUE, saz_continuous = FALSE, ordinary = FALSE),
    spike = c(saz_full = TRUE, saz_continuous = TRUE, ordinary = FALSE),
    spike_decision = c(
      saz_full = saz_decision_codes[["cont_binary"]],
      saz_continuous = saz_decision_codes[["continuous_only"]],
      ordinary = saz_decision_codes[["continuous_only"]]
    )
  )

  expect_equal(fp_terms$df_setting, c(4, 4, 4))
  expect_equal(fp_terms$df_initial, c(5, 5, 4))
  expect_equal(fp_terms$df_final, c(5, 4, 4))
})
