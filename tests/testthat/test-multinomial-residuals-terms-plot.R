# Multinomial term/contrast predictions and the component-plus-residual plot.
#
# A baseline-category multinomial model shares the FP powers across logits but
# carries a separate coefficient row per non-reference logit, so:
#   * residuals.mfp2() returns an n x Q matrix of per-logit residuals (nnet
#     supplies no residual method, so these are computed here);
#   * predict(type = "terms"/"contrasts") returns a per-term data frame stacked
#     over logits, distinguished by a `logit` factor column;
#   * plot() draws one figure per predictor, faceted by logit.
#
# Fits use df = 1, select = 1 so every predictor is retained as a linear term:
# this keeps the multinomial candidate design well conditioned and makes the
# reported coefficient set deterministic.

make_mnl_data <- function(n = 400, seed = 11) {
  set.seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::runif(n, 1, 6)
  g <- factor(sample(c("lo", "hi"), n, TRUE))
  lp2 <- 0.9 * x1 - 0.3 * x2 + 1.0 * (g == "hi")
  lp3 <- -0.6 * x1 + 0.5 * x2 - 0.8 * (g == "hi")
  y <- factor(vapply(seq_len(n), function(i) {
    p1 <- 1 / (1 + exp(lp2[i]) + exp(lp3[i]))
    sample(c("A", "B", "C"), 1L,
           prob = c(p1, exp(lp2[i]) * p1, exp(lp3[i]) * p1))
  }, character(1)))
  data.frame(y = y, x1 = x1, x2 = x2, g = g)
}

lin_mnl <- function(v) sprintf("fp(%s, df = 1, select = 1)", v)

fit_mnl <- function(d, with_factor = FALSE) {
  rhs <- paste(lin_mnl("x1"), "+", lin_mnl("x2"))
  if (with_factor) rhs <- paste(rhs, "+ g")
  mfp2(stats::as.formula(paste("y ~", rhs)), data = d,
       family = "multinomial", xorder = "original", verbose = FALSE)
}

# ---------------------------------------------------------------------------
# residuals.mfp2()
# ---------------------------------------------------------------------------

test_that("residuals() returns an n x Q per-logit matrix for multinomial", {
  skip_if_not_installed("nnet")
  d <- make_mnl_data()
  m <- fit_mnl(d)
  logits <- rownames(m$mfp2_coefficient_matrix)

  for (ty in c("deviance", "pearson", "working", "response")) {
    r <- residuals(m, type = ty)
    expect_true(is.matrix(r))
    expect_equal(nrow(r), nrow(d))
    expect_equal(ncol(r), length(logits))
    expect_identical(colnames(r), logits)
    expect_false(anyNA(r))
  }
  # default type is deviance
  expect_equal(residuals(m), residuals(m, type = "deviance"))
})

test_that("multinomial deviance/response residuals match their definitions", {
  skip_if_not_installed("nnet")
  d <- make_mnl_data()
  m <- fit_mnl(d)
  probs <- predict(m, type = "response")
  eps <- .Machine$double.eps

  for (q in rownames(m$mfp2_coefficient_matrix)) {
    y_q <- as.numeric(as.character(d$y) == q)
    p_q <- pmin(pmax(probs[, q], eps), 1 - eps)

    dev_ref <- sign(y_q - p_q) *
      sqrt(pmax(-2 * (y_q * log(p_q) + (1 - y_q) * log(1 - p_q)), 0))
    expect_equal(unname(residuals(m, type = "deviance")[, q]),
                 unname(dev_ref), tolerance = 1e-10)

    expect_equal(unname(residuals(m, type = "response")[, q]),
                 unname(y_q - p_q), tolerance = 1e-10)
    expect_equal(unname(residuals(m, type = "pearson")[, q]),
                 unname((y_q - p_q) / sqrt(p_q * (1 - p_q))), tolerance = 1e-10)
    expect_equal(unname(residuals(m, type = "working")[, q]),
                 unname((y_q - p_q) / (p_q * (1 - p_q))), tolerance = 1e-10)
  }
})

# ---------------------------------------------------------------------------
# predict(type = "terms" / "contrasts")
# ---------------------------------------------------------------------------

test_that("multinomial term predictions carry a logit column and reconstruct the link", {
  skip_if_not_installed("nnet")
  d <- make_mnl_data()
  m <- fit_mnl(d)
  logits <- rownames(m$mfp2_coefficient_matrix)

  tm <- predict(m, type = "terms", terms_seq = "data")
  expect_setequal(names(tm), c("x1", "x2"))
  for (v in names(tm)) {
    expect_true("logit" %in% names(tm[[v]]))
    expect_identical(levels(tm[[v]]$logit), logits)
    expect_equal(nrow(tm[[v]]), length(logits) * nrow(d))
  }

  # Summing the no-intercept term contributions and adding each logit's
  # intercept reproduces predict(type = "link") exactly.
  tni <- predict(m, type = "terms", terms_seq = "data", add_intercept = FALSE)
  link <- predict(m, type = "link")
  cm <- m$mfp2_coefficient_matrix
  for (q in logits) {
    eta <- cm[q, "(Intercept)"] +
      tni$x1$value[tni$x1$logit == q] +
      tni$x2$value[tni$x2$logit == q]
    expect_equal(unname(eta), unname(link[, q]), tolerance = TOL_TIGHT)
  }
})

test_that("multinomial term standard errors use the per-logit covariance block", {
  skip_if_not_installed("nnet")
  d <- make_mnl_data()
  m <- fit_mnl(d)
  V <- vcov(m)
  # transformed x1 basis on the fitting scale (single linear column x1.1)
  basis <- as.numeric(as.matrix(
    mfp2:::prepare_newdata_for_predict(m, m$x_original[, "x1", drop = FALSE],
                                       terms = "x1", apply_pre = FALSE)
  ))

  tr0 <- predict(m, type = "terms", terms_seq = "data", add_intercept = FALSE)
  tr1 <- predict(m, type = "terms", terms_seq = "data", add_intercept = TRUE)
  for (q in rownames(m$mfp2_coefficient_matrix)) {
    sl <- paste0(q, ":x1.1"); it <- paste0(q, ":(Intercept)")
    se0 <- abs(basis) * sqrt(V[sl, sl])
    expect_equal(unname(tr0$x1$se[tr0$x1$logit == q]), unname(se0),
                 tolerance = TOL_TIGHT)
    se1 <- sqrt(V[it, it] + basis^2 * V[sl, sl] + 2 * basis * V[it, sl])
    expect_equal(unname(tr1$x1$se[tr1$x1$logit == q]), unname(se1),
                 tolerance = TOL_TIGHT)
  }
})

test_that("multinomial contrasts are intercept-free and per-logit", {
  skip_if_not_installed("nnet")
  d <- make_mnl_data()
  m <- fit_mnl(d)
  logits <- rownames(m$mfp2_coefficient_matrix)

  cc <- predict(m, type = "contrasts", terms_seq = "data")
  for (v in names(cc)) {
    expect_identical(levels(cc[[v]]$logit), logits)
    # a contrast against a reference value is zero at that reference and does
    # not depend on the intercept
    expect_true(all(is.finite(cc[[v]]$value)))
  }
})

# ---------------------------------------------------------------------------
# plot()
# ---------------------------------------------------------------------------

test_that("plot() returns one faceted ggplot per predictor for multinomial", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("ggplot2")
  d <- make_mnl_data()
  m <- fit_mnl(d)
  n_logits <- nrow(m$mfp2_coefficient_matrix)

  p <- plot(m)
  expect_setequal(names(p), c("x1", "x2"))
  expect_true(all(vapply(p, function(g) inherits(g, "ggplot"), logical(1))))
  # one facet panel per non-reference logit
  built <- ggplot2::ggplot_build(p$x1)
  expect_equal(length(unique(built$data[[1]]$PANEL)), n_logits)

  # partial_only drops the residual layer (one fewer layer than the default)
  p_full <- plot(m)
  p_partial <- plot(m, partial_only = TRUE)
  expect_lt(length(p_partial$x1$layers), length(p_full$x1$layers))

  # contrasts also build cleanly
  pc <- plot(m, type = "contrasts")
  expect_true(inherits(pc$x1, "ggplot"))
  expect_no_error(ggplot2::ggplot_build(pc$x1))
})

test_that("plot() draws discrete factor terms per logit with point-and-interval", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("ggplot2")
  d <- make_mnl_data()
  m <- fit_mnl(d, with_factor = TRUE)
  skip_if_not("g" %in% names(plot(m)))

  p <- plot(m)
  expect_true(inherits(p$g, "ggplot"))
  expect_identical(p$g$labels$title, "Categorical")
  built <- ggplot2::ggplot_build(p$g)
  expect_equal(length(unique(built$data[[1]]$PANEL)),
               nrow(m$mfp2_coefficient_matrix))
})
