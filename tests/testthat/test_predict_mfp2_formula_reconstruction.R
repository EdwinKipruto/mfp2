# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# Test purpose: Checks that formula-interface prediction preserves the number
# of rows in newdata even when model selection drops all predictors.
test_that("predict.mfp2() returns newdata-length predictions for intercept-only formula models", {
  set.seed(109)
  n <- 120

  dat <- data.frame(
    y = rnorm(n),
    x = runif(n, 1, 10),
    group = factor(rep(c("A", "B", "C"), length.out = n))
  )

  fit <- mfp2(
    y ~ fp(x) + group,
    data = dat,
    verbose = FALSE
  )

  pred <- predict(fit, newdata = dat[1:10, ])

  expect_length(pred, 10)
  expect_true(all(is.finite(pred)))
})


# Test purpose: Checks that matrix-interface prediction preserves newdata row
# count when the final selected model is intercept-only.
test_that("predict.mfp2() returns newdata-length predictions for intercept-only matrix models", {
  set.seed(110)
  n <- 120

  x <- cbind(
    x1 = runif(n, 1, 10),
    x2 = runif(n, 1, 10)
  )
  y <- rnorm(n)

  fit <- mfp2(
    x,
    y,
    verbose = FALSE
  )

  pred <- predict(fit, newdata = x[1:10, , drop = FALSE])

  expect_length(pred, 10)
  expect_true(all(is.finite(pred)))
})


# Regression: full formula prediction must not require terms eliminated during
# selection. This also exercises the stored fp() formula-label mapping.
test_that("predict.mfp2() formula prediction requires only selected terms", {
  set.seed(1101)
  n <- 180
  dat <- data.frame(
    x = runif(n, 1, 10),
    group = factor(rep(c("A", "B", "C"), length.out = n))
  )
  dat$y <- 1 + 0.8 * dat$x + rnorm(n, sd = 0.15)

  fit <- mfp2(
    y ~ fp(x) + group,
    data = dat,
    keep = "x",
    select = 0,
    verbose = FALSE
  )

  expect_true(fit$fp_terms["x", "selected"])
  expect_false(fit$fp_terms["group", "selected"])
  expect_identical(fit$formula_prediction_term_names[["fp(x)"]], "x")

  nd_minimal <- dat[1:12, "x", drop = FALSE]
  pred_minimal <- predict(fit, newdata = nd_minimal)
  pred_complete <- predict(fit, newdata = dat[1:12, c("x", "group")])
  nd_with_unseen_dropped_level <- data.frame(
    x = nd_minimal$x,
    group = factor(rep("D", nrow(nd_minimal)), levels = c("A", "B", "C", "D"))
  )
  pred_unseen_dropped <- predict(fit, newdata = nd_with_unseen_dropped_level)
  term_minimal <- predict(
    fit,
    newdata = nd_minimal,
    type = "terms",
    terms = "x",
    terms_seq = "data"
  )

  expect_equal(pred_minimal, pred_complete, tolerance = 1e-12)
  expect_equal(pred_minimal, pred_unseen_dropped, tolerance = 1e-12)
  expect_named(term_minimal, "x")
  expect_equal(nrow(term_minimal$x), nrow(nd_minimal))
})


# Regression: retaining a categorical block must not require an eliminated
# continuous candidate in formula-style newdata.
test_that("predict.mfp2() reconstructs only a selected factor block", {
  set.seed(1102)
  n <- 180
  dat <- data.frame(
    x = runif(n, 1, 10),
    group = factor(rep(c("A", "B", "C"), length.out = n))
  )
  group_effect <- c(A = 0, B = 1, C = -0.7)[as.character(dat$group)]
  dat$y <- group_effect + rnorm(n, sd = 0.15)

  fit <- mfp2(
    y ~ x + group,
    data = dat,
    keep = "group",
    select = 0,
    verbose = FALSE
  )

  expect_false(fit$fp_terms["x", "selected"])
  expect_true(fit$fp_terms["group", "selected"])

  nd_minimal <- dat[1:12, "group", drop = FALSE]
  pred_minimal <- predict(fit, newdata = nd_minimal)
  pred_complete <- predict(fit, newdata = dat[1:12, c("x", "group")])

  expect_equal(pred_minimal, pred_complete, tolerance = 1e-12)
})


# Regression: selected factor interactions retain the fit-time contrast coding
# while unrelated eliminated formula terms are not required.
test_that("predict.mfp2() preserves active interaction design with minimal newdata", {
  set.seed(11021)
  n <- 210
  dat <- data.frame(
    x = runif(n, 1, 5),
    z = runif(n, -2, 2),
    group = factor(rep(c("A", "B", "C"), length.out = n))
  )
  group_slope <- c(A = 0.2, B = 0.8, C = -0.5)[as.character(dat$group)]
  dat$y <- group_slope * dat$x + rnorm(n, sd = 0.1)

  fit <- mfp2(
    y ~ x + group + x:group + z,
    data = dat,
    keep = "x:group",
    select = 0,
    verbose = FALSE
  )

  expect_true(fit$fp_terms["x:group", "selected"])
  expect_false(fit$fp_terms["z", "selected"])

  nd_minimal <- dat[1:15, c("x", "group"), drop = FALSE]
  pred_minimal <- predict(fit, newdata = nd_minimal)
  pred_complete <- predict(
    fit,
    newdata = dat[1:15, c("x", "group", "z"), drop = FALSE]
  )

  expect_equal(pred_minimal, pred_complete, tolerance = 1e-12)
})


# Regression: matrix prediction expands only selected conceptual terms, so a
# grouped block eliminated during selection is not a newdata dependency.
test_that("predict.mfp2() matrix prediction ignores eliminated grouped columns", {
  set.seed(1103)
  n <- 180
  x <- cbind(
    x1 = runif(n, 1, 10),
    groupB = rep(c(0, 1, 0), length.out = n),
    groupC = rep(c(0, 0, 1), length.out = n)
  )
  y <- 1 + 0.7 * x[, "x1"] + rnorm(n, sd = 0.15)

  fit <- mfp2(
    x,
    y,
    term_groups = list(group = c("groupB", "groupC")),
    keep = "x1",
    select = 0,
    df = 1,
    verbose = FALSE
  )

  expect_true(fit$fp_terms["x1", "selected"])
  expect_false(fit$fp_terms["group", "selected"])
  expect_equal(as.numeric(fit$fp_terms["group", "df_final"]), 0)

  nd_minimal <- x[1:12, "x1", drop = FALSE]
  pred_minimal <- predict(fit, newdata = nd_minimal)
  pred_complete <- predict(fit, newdata = x[1:12, , drop = FALSE])

  expect_equal(pred_minimal, pred_complete, tolerance = 1e-12)
})


# Regression: an intercept-only final model has no predictor dependencies and
# therefore accepts a zero-column data frame carrying only the requested rows.
test_that("predict.mfp2() intercept-only models accept zero-column newdata", {
  set.seed(1104)
  n <- 160
  dat <- data.frame(
    y = rnorm(n),
    x = runif(n, 1, 10),
    group = factor(rep(c("A", "B", "C"), length.out = n))
  )

  fit_formula <- mfp2(
    y ~ x + group,
    data = dat,
    select = 0,
    verbose = FALSE
  )
  fit_matrix <- mfp2(
    as.matrix(dat[, "x", drop = FALSE]),
    dat$y,
    select = 0,
    verbose = FALSE
  )

  expect_length(get_selected_variable_names(fit_formula), 0L)
  expect_length(get_selected_variable_names(fit_matrix), 0L)

  nd_empty <- data.frame(row.names = seq_len(9L))
  expect_length(predict(fit_formula, newdata = nd_empty), 9L)
  expect_length(predict(fit_matrix, newdata = nd_empty), 9L)
  omitted_terms <- NULL
  expect_warning(
    omitted_terms <- predict(
      fit_formula,
      newdata = nd_empty,
      type = "terms",
      terms = "group"
    ),
    "All the terms supplied are not in the final model"
  )
  expect_identical(omitted_terms, list())
})


# Test purpose: Checks that formula-fitted models can predict from ordinary
# newdata containing original factor variables rather than expanded dummy columns.
test_that("predict.mfp2() reconstructs formula-interface newdata with retained factors", {
  set.seed(111)
  n <- 120

  x <- runif(n, 1, 10)
  group <- factor(rep(c("A", "B", "C"), length.out = n))
  group_effect <- c(A = 0, B = 1, C = 2)[as.character(group)]

  dat <- data.frame(
    y = 0.2 * x + 0.5 * group_effect + rnorm(n, sd = 0.2),
    x = x,
    group = group
  )

  fit <- mfp2(
    y ~ fp(x) + group,
    data = dat,
    keep = c("x", "group"),
    verbose = FALSE
  )

  pred <- predict(fit, newdata = dat[1:10, ])

  expect_length(pred, 10)
  expect_true(all(is.finite(pred)))
})


# Test purpose: Checks conceptual-term prediction for an unordered factor,
# including level labels, raw dummy columns, fitted values, and uncertainty.
test_that("predict.mfp2() returns a complete unordered-factor term block", {
  set.seed(2201)
  n <- 150

  dat <- data.frame(
    group = factor(rep(c("A", "B", "C"), length.out = n)),
    x = runif(n, 1, 10)
  )
  effect <- c(A = 0, B = 1, C = -0.5)[as.character(dat$group)]
  dat$y <- 0.3 * dat$x + effect + rnorm(n, sd = 0.2)

  fit <- mfp2(
    y ~ x + group,
    data = dat,
    keep = "group",
    verbose = FALSE
  )
  nd <- dat[1:12, c("x", "group"), drop = FALSE]

  out <- predict(
    fit,
    newdata = nd,
    type = "terms",
    terms = "group"
  )

  expect_named(out, "group")
  expect_s3_class(out$group, "data.frame")
  expect_equal(nrow(out$group), nrow(nd))
  expect_true(
    all(c(
      "variable", "variable_pre", "groupB", "groupC",
      "value", "se", "lower", "upper"
    ) %in% names(out$group))
  )
  expect_equal(out$group$variable, as.character(nd$group))
  expect_true(all(is.finite(out$group$value)))
  expect_true(all(is.finite(out$group$se)))
})


# Test purpose: Checks grouped factor contrasts against a named fitted level.
test_that("predict.mfp2() computes factor contrasts from a level reference", {
  set.seed(2202)
  n <- 150

  dat <- data.frame(
    group = factor(rep(c("A", "B", "C"), length.out = n)),
    x = runif(n, 1, 10)
  )
  effect <- c(A = 0, B = 1, C = -0.5)[as.character(dat$group)]
  dat$y <- 0.3 * dat$x + effect + rnorm(n, sd = 0.2)

  fit <- mfp2(
    y ~ x + group,
    data = dat,
    keep = "group",
    verbose = FALSE
  )
  nd <- dat[1:12, c("x", "group"), drop = FALSE]

  out <- predict(
    fit,
    newdata = nd,
    type = "contrasts",
    terms = "group",
    ref = list(group = "A")
  )

  expect_named(out, "group")
  expect_equal(nrow(out$group), nrow(nd))
  expect_true(all(is.finite(out$group$value)))
  expect_true(all(is.finite(out$group$se)))
  expect_equal(
    out$group$value[nd$group == "A"],
    rep(0, sum(nd$group == "A")),
    tolerance = 1e-10
  )
})


# Test purpose: Checks that ordered-factor contrasts can be reconstructed from
# ordinary formula-style newdata for full and conceptual-term prediction.
test_that("predict.mfp2() reconstructs ordered-factor contrasts", {
  set.seed(2203)
  n <- 180

  dat <- data.frame(
    severity = ordered(
      rep(c("low", "medium", "high"), length.out = n),
      levels = c("low", "medium", "high")
    ),
    x = runif(n, 1, 10)
  )
  severity_effect <- c(low = 0, medium = 0.5, high = 1.2)[
    as.character(dat$severity)
  ]
  dat$y <- 0.2 * dat$x + severity_effect + rnorm(n, sd = 0.2)

  fit <- mfp2(
    y ~ x + severity,
    data = dat,
    keep = "severity",
    verbose = FALSE
  )
  nd <- dat[1:15, c("x", "severity"), drop = FALSE]

  ordinary <- predict(fit, newdata = nd)
  term_output <- predict(
    fit,
    newdata = nd,
    type = "terms",
    terms = "severity"
  )

  expect_length(ordinary, nrow(nd))
  expect_true(all(is.finite(ordinary)))
  expect_named(term_output, "severity")
  expect_equal(term_output$severity$variable, as.character(nd$severity))
  expect_true(all(is.finite(term_output$severity$value)))
})


# Test purpose: Ensures every fitted model stores a complete conceptual-term
# lookup, including ordinary identity mappings required by prediction.
test_that("mfp2 stores complete identity term mappings", {
  set.seed(22031)
  dat <- data.frame(
    y = rnorm(80),
    x1 = runif(80, 1, 4),
    x2 = runif(80, 2, 6)
  )

  fit <- mfp2(
    y ~ x1 + x2,
    data = dat,
    df = 1,
    select = 1,
    verbose = FALSE
  )

  expect_identical(fit$term_to_columns$x1, "x1")
  expect_identical(fit$term_to_columns$x2, "x2")
  expect_length(predict(fit, newdata = dat[1:5, ]), 5L)
})


# Test purpose: Checks that alternative simple factor wrappers follow the same
# source-variable conceptual naming contract as factor().
test_that("simple factor wrappers use source-variable conceptual names", {
  set.seed(22032)
  dat <- data.frame(
    y = rnorm(90),
    x = rep(1:3, length.out = 90)
  )

  fit <- mfp2(
    y ~ as.factor(x),
    data = dat,
    keep = "x",
    verbose = FALSE
  )

  expect_identical(
    fit$term_to_columns$x,
    c("as.factor(x)2", "as.factor(x)3")
  )
  expect_true("x" %in% rownames(fit$fp_terms))
  expect_false("as.factor(x)" %in% rownames(fit$fp_terms))
})


# Test purpose: Checks that a simple inline factor() wrapper uses the source
# variable name as the conceptual term while retaining model.matrix() column and
# coefficient names.
test_that("formula interface names inline factors by their source variable", {
  set.seed(2204)
  n <- 150

  dat <- data.frame(
    x1 = runif(n, 1, 10),
    x2 = rep(1:3, length.out = n)
  )
  x2_effect <- c(`1` = 0, `2` = 0.8, `3` = -0.4)[as.character(dat$x2)]
  dat$y <- 0.3 * dat$x1 + x2_effect + rnorm(n, sd = 0.2)

  fit <- mfp2(
    y ~ x1 + factor(x2),
    data = dat,
    keep = "x2",
    verbose = FALSE
  )

  expect_true("x2" %in% names(fit$term_to_columns))
  expect_false("factor(x2)" %in% names(fit$term_to_columns))
  expect_identical(
    fit$term_to_columns[["x2"]],
    c("factor(x2)2", "factor(x2)3")
  )
  expect_true("x2" %in% rownames(fit$fp_terms))
  expect_false("factor(x2)" %in% rownames(fit$fp_terms))

  coefficient_names <- gsub("`", "", names(stats::coef(fit)), fixed = TRUE)
  expect_true(all(c("factor(x2)2.1", "factor(x2)3.1") %in% coefficient_names))

  nd <- dat[1:10, c("x1", "x2"), drop = FALSE]
  pred <- predict(fit, newdata = nd)
  term_pred <- predict(
    fit,
    newdata = nd,
    type = "terms",
    terms = "x2"
  )

  expect_length(pred, 10L)
  expect_true(all(is.finite(pred)))
  expect_named(term_pred, "x2")
  expect_true(all(is.finite(term_pred[["x2"]]$value)))
  expect_true(all(is.finite(term_pred[["x2"]]$se)))
})


# Test purpose: Checks that a binary inline factor remains a mapped categorical
# term even though model.matrix() generates only one dummy column.
test_that("binary inline factors retain a source-name singleton mapping", {
  set.seed(22041)
  n <- 140

  dat <- data.frame(
    x1 = runif(n, 1, 10),
    x2 = rep(1:2, length.out = n)
  )
  dat$y <- 0.4 * dat$x1 + 0.9 * (dat$x2 == 2) + rnorm(n, sd = 0.2)

  fit <- mfp2(
    y ~ x1 + factor(x2),
    data = dat,
    keep = "x2",
    verbose = FALSE
  )

  expect_identical(fit$term_to_columns[["x2"]], "factor(x2)2")
  expect_true(fit$fp_terms["x2", "selected"])
  expect_equal(as.numeric(fit$fp_terms["x2", "df_final"]), 1)

  nd <- dat[1:12, c("x1", "x2"), drop = FALSE]
  out <- predict(fit, newdata = nd, type = "terms", terms = "x2")
  expect_named(out, "x2")
  expect_true(all(is.finite(out$x2$value)))
  expect_true(all(is.finite(out$x2$se)))
})


# Test purpose: Prevents ambiguous formulas that include the same source
# variable both directly and through a simple factor wrapper.
test_that("formula interface rejects duplicate conceptual source variables", {
  dat <- data.frame(
    # The duplicate-term validation occurs before fitting; no RNG is needed.
    y = seq_len(60),
    x2 = rep(1:3, length.out = 60)
  )

  expect_error(
    mfp2(y ~ x2 + factor(x2), data = dat, verbose = FALSE),
    "same conceptual variable"
  )
})


# Test purpose: Checks that formula reconstruction rejects factor levels that
# were not present when the model matrix and contrasts were fitted.
test_that("predict.mfp2() rejects unseen factor levels", {
  set.seed(2205)
  n <- 120

  dat <- data.frame(
    y = rnorm(n),
    group = factor(rep(c("A", "B", "C"), length.out = n))
  )
  fit <- mfp2(
    y ~ group,
    data = dat,
    keep = "group",
    verbose = FALSE
  )
  nd <- data.frame(
    group = factor("D", levels = c("A", "B", "C", "D"))
  )

  expect_error(
    predict(fit, newdata = nd),
    "new level|factor level|levels",
    ignore.case = TRUE
  )
})


# Test purpose: Checks that matrix-interface prediction rejects a partial grouped
# block rather than silently changing the fitted parameterisation.
test_that("matrix prediction requires every grouped-term column", {
  set.seed(2206)
  n <- 120

  x <- cbind(
    groupB = rep(c(0, 1, 0), length.out = n),
    groupC = rep(c(0, 0, 1), length.out = n),
    x1 = runif(n, 1, 10)
  )
  y <- x[, "groupB"] - x[, "groupC"] +
    0.2 * x[, "x1"] + rnorm(n)

  fit <- mfp2(
    x,
    y,
    df = 1,
    select = 1,
    term_groups = list(group = c("groupB", "groupC")),
    verbose = FALSE
  )
  incomplete <- x[1:5, c("groupB", "x1"), drop = FALSE]

  expect_error(
    predict(fit, newdata = incomplete),
    "groupC|grouped term|grouped-term",
    ignore.case = TRUE
  )
})
