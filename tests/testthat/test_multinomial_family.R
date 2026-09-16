# Multinomial family, common-power selection, and nnet equivalence -----------

make_multinomial_test_data <- function(n = 360L, seed = 1609L) {
  set.seed(seed)
  dat <- data.frame(
    x1 = stats::runif(n, 0.5, 4.5),
    x2 = stats::rnorm(n)
  )
  eta_b <- -0.45 + 0.35 * dat$x1 - 0.30 * dat$x2
  eta_c <- 0.25 - 0.20 * dat$x1 + 0.45 * dat$x2
  denom <- 1 + exp(eta_b) + exp(eta_c)
  probabilities <- cbind(A = 1 / denom, B = exp(eta_b) / denom,
                         C = exp(eta_c) / denom)
  draw <- apply(probabilities, 1L, function(p) sample.int(3L, 1L, prob = p))
  dat$y <- factor(colnames(probabilities)[draw], levels = colnames(probabilities))
  dat
}


test_that("multinomial family normalization and response validation are strict", {
  expect_identical(
    normalize_family_argument("multinomial")$family_string,
    "multinomial"
  )
  expect_s3_class(multinomial_family(), "mfp2_multinomial_family")
  expect_error(multinomial_family(reference = c("A", "B")), "reference")
  expect_error(
    validate_family_response(factor(c("A", "B", "A")), "multinomial", 3L),
    "at least three"
  )
  expect_error(
    validate_family_response(cbind(A = c(1, 0), B = c(0, 1)),
                             "multinomial", 2L),
    "at least three"
  )
  expect_error(
    validate_family_response(
      cbind(A = c(1, 0), B = c(0, 0), C = c(0, 0)),
      "multinomial", 2L
    ),
    "at least one trial"
  )
  expect_invisible(
    validate_family_response(c(1L, 2L, 3L), "multinomial", 3L)
  )
  expect_invisible(
    validate_family_response(c("A", "B", "C"), "multinomial", 3L)
  )
  expect_identical(mfp2_multinomial_n_classes(c(1L, 2L, 3L)), 3L)
  expect_error(
    validate_family_response(c(1, 2, Inf), "multinomial", 3L),
    "finite"
  )
})


test_that("integer-coded multinomial labels are converted to factor classes", {
  dat <- make_multinomial_test_data(n = 240L, seed = 1617L)
  dat$y <- as.integer(dat$y)

  fit <- mfp2(
    y ~ x1 + x2,
    data = dat,
    family = multinomial_family(),
    cycles = 1,
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )
  reference <- nnet::multinom(
    factor(y) ~ x1 + x2,
    data = dat,
    Hess = TRUE,
    trace = FALSE
  )

  expect_identical(fit$class_levels, c("1", "2", "3"))
  expect_identical(fit$reference_class, "1")
  expect_equal(
    unname(stats::coef(fit)),
    unname(stats::coef(reference)),
    tolerance = 1e-6
  )
  expect_equal(
    predict(fit, newdata = dat[1:20, ], type = "response"),
    predict(reference, newdata = dat[1:20, ], type = "probs"),
    tolerance = 1e-6
  )
})


test_that("forced-linear formula fit matches nnet::multinom", {
  dat <- make_multinomial_test_data()
  w <- stats::runif(nrow(dat), 0.7, 1.8)

  fit <- mfp2(
    y ~ x1 + x2,
    data = dat,
    weights = w,
    family = multinomial_family(reference = "A"),
    cycles = 1,
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )
  reference <- nnet::multinom(y ~ x1 + x2, data = dat, weights = w,
                              Hess = TRUE, trace = FALSE)

  expect_s3_class(fit, "mfp2")
  expect_s3_class(fit, "multinom")
  expect_identical(fit$family_string, "multinomial")
  expect_identical(fit$reference_class, "A")
  expect_identical(fit$n_logits, 2L)
  expect_equal(unname(stats::coef(fit)), unname(stats::coef(reference)),
               tolerance = 1e-6)
  expect_equal(unname(stats::vcov(fit)), unname(stats::vcov(reference)),
               tolerance = 1e-5)
  expect_equal(as.numeric(stats::logLik(fit)),
               as.numeric(stats::logLik(reference)), tolerance = 1e-6)

  nd <- dat[1:35, c("x1", "x2"), drop = FALSE]
  probabilities <- predict(fit, newdata = nd, type = "response")
  probabilities_reference <- predict(reference, newdata = nd, type = "probs")
  classes <- predict(fit, newdata = nd, type = "class")
  classes_reference <- predict(reference, newdata = nd, type = "class")
  logits <- predict(fit, newdata = nd, type = "link")
  logits_reference <- log(
    probabilities_reference[, c("B", "C"), drop = FALSE] /
      probabilities_reference[, "A"]
  )

  expect_identical(rownames(probabilities), row.names(nd))
  expect_identical(rownames(logits), row.names(nd))
  expect_equal(probabilities, probabilities_reference, tolerance = 1e-6)
  expect_equal(logits, logits_reference, tolerance = 1e-6)
  expect_identical(as.character(classes), as.character(classes_reference))

  manual_x <- cbind(`(Intercept)` = 1, x1 = nd$x1, x2 = nd$x2)
  manual_logits <- manual_x %*% t(stats::coef(fit))
  expect_equal(unname(logits), unname(manual_logits), tolerance = 1e-8)
  expect_equal(
    unname(rowSums(probabilities)),
    rep(1, nrow(nd)),
    tolerance = 1e-12
  )
})


test_that("matrix interface and non-first reference preserve public class order", {
  dat <- make_multinomial_test_data(seed = 1610L)
  x <- as.matrix(dat[, c("x1", "x2")])
  y_reference <- stats::relevel(dat$y, ref = "C")

  fit <- mfp2(
    x = x,
    y = dat$y,
    family = multinomial_family(reference = "C"),
    cycles = 1,
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )
  reference <- nnet::multinom(y_reference ~ x1 + x2, data = dat,
                              Hess = TRUE, trace = FALSE)

  expect_equal(unname(stats::coef(fit)), unname(stats::coef(reference)),
               tolerance = 1e-6)
  p <- predict(fit, newdata = x[1:20, , drop = FALSE], type = "response")
  p_reference <- predict(reference, newdata = dat[1:20, ], type = "probs")
  p_reference <- p_reference[, levels(dat$y), drop = FALSE]
  expect_identical(colnames(p), levels(dat$y))
  expect_equal(p, p_reference, tolerance = 1e-6)
})


test_that("grouped multinomial counts match nnet::multinom when linear", {
  dat <- make_multinomial_test_data(n = 180L, seed = 1611L)
  base_prob <- predict(
    nnet::multinom(y ~ x1 + x2, data = dat, trace = FALSE),
    type = "probs"
  )
  set.seed(1612L)
  counts <- t(vapply(seq_len(nrow(dat)), function(i) {
    as.vector(stats::rmultinom(1L, size = 12L, prob = base_prob[i, ]))
  }, numeric(3L)))
  colnames(counts) <- c("A", "B", "C")
  x <- as.matrix(dat[, c("x1", "x2")])

  fit <- mfp2(
    x, counts,
    family = "multinomial",
    cycles = 1,
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )
  reference_data <- data.frame(dat[, c("x1", "x2")], counts)
  reference <- nnet::multinom(
    cbind(A, B, C) ~ x1 + x2,
    data = reference_data,
    Hess = TRUE,
    trace = FALSE
  )

  expect_equal(unname(stats::coef(fit)), unname(stats::coef(reference)),
               tolerance = 1e-6)
  expect_equal(as.numeric(stats::logLik(fit)),
               as.numeric(stats::logLik(reference)), tolerance = 1e-6)
  expect_equal(
    predict(fit, newdata = x[1:20, , drop = FALSE], type = "response"),
    predict(reference, newdata = reference_data[1:20, ], type = "probs"),
    tolerance = 1e-6
  )
})


test_that("multinomial class offsets match nnet::multinom", {
  dat <- make_multinomial_test_data(n = 240L, seed = 1616L)
  offset_matrix <- cbind(
    A = 0,
    B = 0.10 * sin(dat$x1),
    C = -0.12 * dat$x2
  )

  fit <- mfp2(
    x = as.matrix(dat[, c("x1", "x2")]),
    y = dat$y,
    family = "multinomial",
    offset = offset_matrix,
    cycles = 1,
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )
  reference <- nnet::multinom(
    y ~ x1 + x2 + offset(offset_matrix),
    data = dat,
    Hess = FALSE,
    trace = FALSE
  )

  expect_equal(unname(stats::coef(fit)), unname(stats::coef(reference)),
               tolerance = 1e-6)
  expect_equal(as.numeric(stats::logLik(fit)),
               as.numeric(stats::logLik(reference)), tolerance = 1e-6)
  expect_equal(
    unname(predict(fit, type = "response")),
    unname(stats::fitted(reference)),
    tolerance = 1e-6
  )

  # nnet's own multinomHess() cannot align a matrix-valued fixed offset with
  # the free coefficients. Compare the package covariance with the exact
  # baseline-category logit information matrix instead.
  fitted_probabilities <- stats::fitted(reference)[, c("B", "C"), drop = FALSE]
  design <- stats::model.matrix(~ x1 + x2, data = dat)
  information <- matrix(0, nrow = 6L, ncol = 6L)
  for (j in seq_len(2L)) {
    rows <- (j - 1L) * 3L + seq_len(3L)
    for (k in seq_len(2L)) {
      columns <- (k - 1L) * 3L + seq_len(3L)
      working_weights <- fitted_probabilities[, j] *
        ((j == k) - fitted_probabilities[, k])
      information[rows, columns] <- crossprod(
        design,
        design * working_weights
      )
    }
  }
  expect_equal(
    unname(stats::vcov(fit)),
    unname(solve(information)),
    tolerance = 1e-5
  )
})


test_that("one nonlinear FP power is shared across all logits", {
  dat <- make_multinomial_test_data(n = 260L, seed = 1613L)
  fit <- mfp2(
    x = as.matrix(dat["x1"]),
    y = dat$y,
    family = "multinomial",
    cycles = 5,
    df = 2,
    powers = list(x1 = c(0.5, 2)),
    force_max_fp_vars = "x1",
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )

  expect_length(fit$fp_powers$x1, 1L)
  expect_true(fit$fp_powers$x1 %in% c(0.5, 2))
  expect_equal(nrow(stats::coef(fit)), 2L)
  expect_equal(fit$fp_terms["x1", "df_final"], 3)
  printed <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(printed, "FP powers: common across logits")
  expect_match(printed, "Outcome classes: A, B, C")
  expect_match(printed, "B vs A")
  expect_match(printed, "C vs A")

  summarized <- summary(fit)
  expect_s3_class(summarized, "summary.mfp2")
  expect_true(isTRUE(summarized$multinomial))
  expect_equal(unique(summarized$coefficients$outcome), c("B", "C"))
  expect_match(
    paste(capture.output(print(summarized)), collapse = "\n"),
    "Coefficient Tests"
  )
})


test_that("multinomial prediction rejects unsupported uncertainty and term paths", {
  dat <- make_multinomial_test_data(n = 180L, seed = 1614L)
  fit <- mfp2(
    y ~ x1 + x2, data = dat, family = "multinomial",
    cycles = 1, df = 1, select = 1, alpha = 1,
    shift = 0, scale = 1, center = FALSE, verbose = FALSE
  )
  expect_error(predict(fit, dat[1:3, ], se.fit = TRUE), "not available")
  expect_error(predict(fit, type = "terms"), "not yet available")
})
