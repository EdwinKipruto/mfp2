# Area 1a: multinomial_family() fitting correctness vs nnet::multinom()

test_that("mfp2 multinomial forced-linear fit matches nnet::multinom", {
  skip_on_cran()
  dat <- make_multinomial_data(n = 400, seed = 1)

  fit <- mfp2(y ~ fp(x1, df = 1, select = 1, scale = 1, shift = 0) +
                  fp(x2, df = 1, select = 1, scale = 1, shift = 0),
              data = dat, family = multinomial_family(),
              center = FALSE, control = list(maxit = 300), verbose = FALSE)
  ref <- nnet::multinom(y ~ x1 + x2, data = dat, trace = FALSE, maxit = 500)

  expect_equal(fit$family_string, "multinomial")
  expect_true(inherits(fit, "multinom"))

  # Log-likelihood and parameter count
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)), tolerance = TOL_NNET)
  expect_equal(attr(logLik(fit), "df"), ref$edf)

  # Slope coefficients: same non-reference logits, aligned by name. The
  # intercept is on mfp2's centred convention and is validated through the
  # matching fitted probabilities below, not compared directly.
  cf <- coef(fit); cr <- coef(ref)
  expect_equal(dim(cf), dim(cr))
  slopes <- setdiff(colnames(cf), "(Intercept)")
  ref_slopes <- setdiff(colnames(cr), "(Intercept)")
  expect_equal(max_abs_diff(cf[rownames(cf), slopes, drop = FALSE],
                            cr[rownames(cf), ref_slopes, drop = FALSE]),
               0, tolerance = 1e-2)

  # Fitted probabilities and predicted classes
  pm <- predict(fit, type = "response")
  pn <- predict(ref, type = "probs")
  expect_equal(max_abs_diff(pm, pn[, colnames(pm)]), 0, tolerance = TOL_NNET)

  cm <- as.character(predict(fit, type = "class"))
  cn <- as.character(predict(ref, type = "class"))
  expect_gt(mean(cm == cn), 0.99)
})

test_that("custom reference level makes coefficients relative to that class", {
  skip_on_cran()
  dat <- make_multinomial_data(n = 400, seed = 2)

  fit_ref_c <- mfp2(y ~ fp(x1, df = 1, select = 1, scale = 1, shift = 0),
                    data = dat, family = multinomial_family(reference = "C"),
                    center = FALSE, control = list(maxit = 300), verbose = FALSE)
  # nnet uses the first factor level as reference; relevel to match
  dat_c <- dat; dat_c$y <- relevel(dat$y, ref = "C")
  ref <- nnet::multinom(y ~ x1, data = dat_c, trace = FALSE, maxit = 500)

  # non-reference logits are A and B (C is the reference)
  expect_setequal(rownames(coef(fit_ref_c)), c("A", "B"))
  expect_equal(fit_ref_c$mfp2_reference_class, "C")
  expect_equal(as.numeric(logLik(fit_ref_c)), as.numeric(logLik(ref)),
               tolerance = TOL_NNET)
  # slopes relative to the chosen reference match nnet (intercept excluded)
  cf <- coef(fit_ref_c); cr <- coef(ref)
  slp <- setdiff(colnames(cf), "(Intercept)")
  crp <- setdiff(colnames(cr), "(Intercept)")
  expect_equal(max_abs_diff(cf[rownames(cf), slp, drop = FALSE],
                            cr[rownames(cf), crp, drop = FALSE]),
               0, tolerance = 1e-2)
})

test_that("count-matrix response produces the same fit as the factor response", {
  skip_on_cran()
  dat <- make_multinomial_data(n = 300, seed = 3)

  fit_factor <- mfp2(y ~ fp(x1, df = 1, select = 1, scale = 1, shift = 0),
                     data = dat, family = multinomial_family(),
                     center = FALSE, control = list(maxit = 300), verbose = FALSE)

  # one-hot count matrix with the same class ordering
  ymat <- nnet::class.ind(dat$y)
  colnames(ymat) <- levels(dat$y)
  dat2 <- data.frame(x1 = dat$x1)
  dat2$Y <- ymat
  fit_counts <- mfp2(Y ~ fp(x1, df = 1, select = 1, scale = 1, shift = 0),
                     data = dat2, family = multinomial_family(),
                     center = FALSE, control = list(maxit = 300), verbose = FALSE)

  expect_equal(as.numeric(logLik(fit_counts)), as.numeric(logLik(fit_factor)),
               tolerance = TOL_NNET)
  expect_equal(max_abs_diff(coef(fit_counts), coef(fit_factor)), 0,
               tolerance = 1e-2)
})

test_that("multinomial predict returns correct link/response/class shapes", {
  skip_on_cran()
  dat <- make_multinomial_data(n = 300, seed = 4)
  fit <- mfp2(y ~ fp(x1, df = 1, select = 1),
              data = dat, family = multinomial_family(),
              control = list(maxit = 300), verbose = FALSE)

  lnk <- predict(fit, type = "link")
  expect_true(is.matrix(lnk))
  expect_equal(ncol(lnk), nlevels(dat$y) - 1L)          # n x (C-1) logits
  expect_equal(nrow(lnk), nrow(dat))

  rsp <- predict(fit, type = "response")
  expect_true(is.matrix(rsp))
  expect_equal(ncol(rsp), nlevels(dat$y))               # n x C probabilities
  expect_true(all(abs(rowSums(rsp) - 1) < TOL_TIGHT))   # rows sum to 1

  cls <- predict(fit, type = "class")
  expect_length(cls, nrow(dat))
  expect_true(all(as.character(cls) %in% levels(dat$y)))
})
