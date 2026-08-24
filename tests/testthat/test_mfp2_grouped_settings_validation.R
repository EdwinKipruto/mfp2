# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# =============================================================================
# 10. Edge cases and input validation
# =============================================================================

# Test purpose: Defines the shared scalar-or-named-override contract used by
# df, select, and alpha in the matrix interfaces.
test_that("named override settings fill defaults and match by name", {
  columns <- c("age", "bmi", "weight")

  global <- normalize_named_override_setting(
    value = 2,
    column_names = columns,
    default = 4,
    argument_name = "df"
  )
  expect_equal(global$value, c(age = 2, bmi = 2, weight = 2))
  expect_true(global$global_scalar)
  expect_true(all(global$supplied))

  partial <- normalize_named_override_setting(
    value = c(weight = 1, age = 2),
    column_names = columns,
    default = 4,
    argument_name = "df"
  )
  expect_equal(partial$value, c(age = 2, bmi = 4, weight = 1))
  expect_false(partial$global_scalar)
  expect_equal(partial$supplied, c(age = TRUE, bmi = FALSE, weight = TRUE))
})


# Test purpose: Prevents the previous positional assignment behavior and
# rejects malformed or unknown names.
test_that("named override settings reject positional and malformed vectors", {
  columns <- c("age", "bmi", "weight")

  expect_error(
    normalize_named_override_setting(
      value = c(1, 2),
      column_names = columns,
      default = 4,
      argument_name = "df"
    ),
    "single unnamed numeric value or a named numeric vector"
  )

  expect_error(
    normalize_named_override_setting(
      value = c(age = 1, unknown = 2),
      column_names = columns,
      default = 4,
      argument_name = "df"
    ),
    "unknown column name.*unknown"
  )

  expect_error(
    normalize_named_override_setting(
      value = stats::setNames(c(1, 2), c("age", "age")),
      column_names = columns,
      default = 4,
      argument_name = "df"
    ),
    "names must be unique"
  )

  expect_error(
    normalize_named_override_setting(
      value = stats::setNames(c(1, 2), c("age", "")),
      column_names = columns,
      default = 4,
      argument_name = "df"
    ),
    "single unnamed numeric value or a named numeric vector"
  )
})


# Test purpose: Verifies partial named df/select/alpha settings in
# mfp2.default(), including default filling and cardinality reduction for an
# omitted five-level predictor.
test_that("mfp2.default() accepts partial named df select and alpha overrides", {
  set.seed(10001)
  n <- 80L
  x <- cbind(
    age = stats::runif(n, 1, 8),
    bmi = rep(1:5, length.out = n),
    weight = stats::runif(n, 2, 10)
  )
  y <- 0.5 * x[, "age"] + stats::rnorm(n, sd = 0.2)

  fit <- mfp2(
    x,
    y,
    cycles = 5,
    df = c(age = 1),
    select = c(age = 1),
    alpha = c(weight = 1),
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )

  expect_equal(as.numeric(fit$fp_terms["age", "df_initial"]), 1)
  expect_equal(as.numeric(fit$fp_terms["bmi", "df_initial"]), 2)
  expect_equal(as.numeric(fit$fp_terms["weight", "df_initial"]), 4)
  expect_equal(as.numeric(fit$fp_terms["age", "select"]), 1)
  expect_equal(as.numeric(fit$fp_terms["bmi", "select"]), 0.05)
  expect_equal(as.numeric(fit$fp_terms["age", "alpha"]), 0.05)
  expect_equal(as.numeric(fit$fp_terms["weight", "alpha"]), 1)
})


# Test purpose: Checks the public error path rather than only the internal
# normalizer when an unnamed multi-value vector is supplied.
test_that("mfp2.default() rejects unnamed multi-value df select and alpha", {
  x <- cbind(age = 1:30, bmi = seq(2, 8, length.out = 30))
  # Response values are irrelevant because validation fails before fitting.
  y <- seq_len(30)

  expect_error(
    mfp2(x, y, df = c(1, 4), verbose = FALSE),
    "single unnamed numeric value or a named numeric vector"
  )
  expect_error(
    mfp2(x, y, select = c(1, 0.05), verbose = FALSE),
    "single unnamed numeric value or a named numeric vector"
  )
  expect_error(
    mfp2(x, y, alpha = c(1, 0.05), verbose = FALSE),
    "single unnamed numeric value or a named numeric vector"
  )
})


# Test purpose: Checks that every manually grouped member names an existing raw
# design-matrix column.
test_that("term_groups rejects unknown columns", {
  x <- cbind(x1 = 1:20, groupB = rep(0:1, 10))
  # Response values are irrelevant because validation fails before fitting.
  y <- seq_len(20)

  expect_error(
    mfp2(
      x,
      y,
      term_groups = list(group = c("groupB", "groupC"))
    ),
    "unknown column.*groupC|groupC.*unknown column",
    ignore.case = TRUE
  )
})


# Test purpose: Checks that one raw column cannot belong to two conceptual terms.
test_that("term_groups rejects columns in more than one group", {
  x <- cbind(
    groupB = rep(0:1, 10),
    groupC = rep(c(0, 0, 1, 0), 5),
    other = 1:20
  )
  # Response values are irrelevant because validation fails before fitting.
  y <- seq_len(20)

  expect_error(
    mfp2(
      x,
      y,
      term_groups = list(
        group1 = c("groupB", "groupC"),
        group2 = c("groupC", "other")
      )
    ),
    "only one.*Duplicated column|Duplicated column.*groupC",
    ignore.case = TRUE
  )
})


# Test purpose: Checks that a grouped term cannot enter the FP search with a
# nonlinear df setting; grouped blocks are fixed linear terms.
test_that("grouped terms require df = 1 for every member column", {
  x <- cbind(
    group1 = 1:20,
    group2 = (1:20)^2,
    x1 = seq(0.5, 10, length.out = 20)
  )
  # Response values are irrelevant because validation fails before fitting.
  y <- seq_len(20)

  expect_error(
    mfp2(
      x,
      y,
      df = c(group1 = 1, group2 = 4, x1 = 1),
      term_groups = list(group = c("group1", "group2"))
    ),
    paste0(
      "Grouped term 'group'.*`df`.*1.*fixed linear design blocks.*",
      "fractional-polynomial transformation search"
    ),
    ignore.case = TRUE
  )
})


# Test purpose: Checks that continuous-only extensions cannot be assigned to
# any member of a grouped categorical block.
test_that("grouped terms reject continuous-only processing options", {
  x <- cbind(
    groupB = rep(c(0, 1, 0), length.out = 24),
    groupC = rep(c(0, 0, 1), length.out = 24),
    x1 = seq(1, 12, length.out = 24)
  )
  # Response values are irrelevant because validation fails before fitting.
  y <- seq_len(24)
  grouped <- list(group = c("groupB", "groupC"))

  cases <- list(
    acdx = list(acdx = "groupB"),
    zero_vars = list(zero_vars = "groupB"),
    catzero_vars = list(catzero_vars = "groupB"),
    spike_vars = list(spike_vars = "groupB")
  )

  for (setting in names(cases)) {
    args <- c(
      list(x = x, y = y, term_groups = grouped),
      cases[[setting]]
    )
    expect_error(
      do.call(mfp2, args),
      paste0(
        "Grouped term 'group'.*`", setting, "`.*FALSE.*",
        "fixed categorical design blocks.*singleton continuous predictors"
      ),
      ignore.case = TRUE,
      info = paste("setting =", setting)
    )
  }
})


# Test purpose: Confirms that conceptual-term mappings are supplied explicitly
# throughout the fitting chain and are not read from an attribute on x.
test_that("grouped mappings use explicit arguments instead of matrix attributes", {
  x <- cbind(
    groupB = c(0, 1, 0, 1),
    groupC = c(0, 0, 1, 0),
    x1 = c(1, 2, 3, 4)
  )
  attr(x, "mfp2_term_to_columns") <- list(stale = "missing_column")

  term_to_columns <- list(
    group = c("groupB", "groupC"),
    x1 = "x1"
  )
  term_names <- names(term_to_columns)

  out <- build_adjustment_step(
    x = x,
    xi = "x1",
    powers_current = list(group = 1, x1 = 1),
    powers = list(group = 1, x1 = 1),
    acdx = stats::setNames(rep(FALSE, 2L), term_names),
    zero = stats::setNames(rep(FALSE, 2L), term_names),
    catzero = stats::setNames(vector("list", 2L), term_names),
    spike = stats::setNames(rep(FALSE, 2L), term_names),
    spike_decision = stats::setNames(rep(0L, 2L), term_names),
    acd_parameter = stats::setNames(vector("list", 2L), term_names),
    prev_adj_params = stats::setNames(vector("list", 2L), term_names),
    term_to_columns = term_to_columns
  )

  expect_identical(colnames(out$data_adj), c("groupB", "groupC"))
  expect_equal(
    unname(out$data_adj),
    unname(x[, c("groupB", "groupC"), drop = FALSE])
  )
  expect_false(any(grepl(
    "mfp2_term_to_columns",
    deparse(body(fit_mfp), width.cutoff = 500L),
    fixed = TRUE
  )))
})


# Test purpose: Covers the MFPI names that use the shared grouped-setting
# validator and confirms that both option classes receive an actionable reason.
test_that("grouped-setting errors explain why MFPI options are unsupported", {
  term_to_columns <- list(group = c("groupB", "groupC"))

  expect_error(
    validate_grouped_term_setting(
      term_to_columns = term_to_columns,
      values = c(groupB = TRUE, groupC = FALSE),
      setting = "acd_vars",
      predicate = function(v) !v,
      requirement = "FALSE"
    ),
    paste0(
      "fixed categorical design blocks.*",
      "singleton continuous predictors"
    ),
    ignore.case = TRUE
  )

  expect_error(
    validate_grouped_term_setting(
      term_to_columns = term_to_columns,
      values = c(groupB = TRUE, groupC = FALSE),
      setting = "force_max_fp_vars",
      predicate = function(v) !v,
      requirement = "FALSE"
    ),
    paste0(
      "fixed linear design blocks.*",
      "fractional-polynomial transformation search"
    ),
    ignore.case = TRUE
  )
})
