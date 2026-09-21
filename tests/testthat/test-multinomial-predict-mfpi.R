# BUG 2 fix: predict.mfpi() now supports the multinomial family through a
# dedicated per-logit path. These tests verify not just that prediction runs,
# but that the reconstructed fitted functions and delta-method standard errors
# are numerically correct against the stored interaction model itself.

# Build one fitted multinomial mfpi model with a genuine group x covariate
# interaction, reused across the tests below.
mn_pred_fit <- function() {
  set.seed(7); n <- 900; a <- runif(n, 1, 10)
  g <- factor(sample(c("0", "1"), n, TRUE))
  lp2 <- -0.5 + 0.10 * a + 0.9 * (g == "1") * a   # slope of a differs by group
  lp3 <-  0.3 - 0.15 * a + 0.5 * (g == "1")
  y <- factor(vapply(seq_len(n), function(i) {
    p1 <- 1 / (1 + exp(lp2[i]) + exp(lp3[i]))
    sample(c("A", "B", "C"), 1L,
           prob = c(p1, exp(lp2[i]) * p1, exp(lp3[i]) * p1))
  }, character(1)))
  dm <- data.frame(y = y, a = a, g = g)
  list(
    fit = mfpi(y ~ g + fp(a), data = dm, group_var = "g",
               interaction_vars = "a", interaction_forms = c(a = "fp1"),
               family = multinomial_family(), control = list(maxit = 300),
               verbose = FALSE),
    data = dm
  )
}

test_that("fitted functions reproduce the interaction model's linear predictors", {
  skip_on_cran()
  obj <- mn_pred_fit(); fit <- obj$fit; dm <- obj$data

  im <- fit$var_winners[["a"]]$fit$test_results$interaction_model
  D  <- im$fit$mfp2_design               # centred design used by the inner fit
  cf <- im$coefficient_matrix            # Q x P per-logit coefficient matrix
  eta_gold <- D %*% t(cf)                # n x Q model linear predictors
  colnames(eta_gold) <- rownames(cf)

  pr <- predict(fit, newdata = dm, se.fit = TRUE)
  fu <- pr$functions
  disp <- mfp2:::mfpi_prediction_group_display_labels(fit, c("0", "1"))
  names(disp) <- c("0", "1")
  gi <- as.character(dm$g)

  maxerr <- 0
  for (cc in rownames(cf)) {
    for (grp in c("0", "1")) {
      rows <- which(gi == grp)
      sub <- fu[fu$class == cc & fu$group == disp[[grp]], ]
      # Each group's curve is evaluated at every newdata x; compare rows that
      # actually belong to that group against the model's own linear predictor.
      maxerr <- max(maxerr, max(abs(sub$fit[rows] - eta_gold[rows, cc])))
    }
  }
  expect_lt(maxerr, 1e-8)
})

test_that("fitted-function delta-method SEs match a direct computation", {
  skip_on_cran()
  obj <- mn_pred_fit(); fit <- obj$fit; dm <- obj$data
  fr <- fit$var_winners[["a"]]$fit
  im <- fr$test_results$interaction_model
  V  <- vcov(im$fit)

  # Reconstruct the dense fitted-function basis exactly as predict() does.
  pd <- mfp2:::mfpi_prepare_prediction_data(fit, "a", newdata = dm, grid = FALSE,
                                            n_grid = 200L, purpose = "function")
  basis <- mfp2:::mfpi_build_function_basis(fit, "a", fr,
                                            pd$cont_var_scaled, pd$x_display)
  pr <- predict(fit, newdata = dm, se.fit = TRUE)
  disp <- mfp2:::mfpi_prediction_group_display_labels(fit, c("0", "1"))
  names(disp) <- c("0", "1")

  i <- 5L
  # Baseline group ("0"): design row is [intercept, a01]; class B block of V.
  xrow <- c(1, basis$x[i, "a01"]); vn <- paste0("B:", c("(Intercept)", "a01"))
  se_hand <- sqrt(as.numeric(t(xrow) %*% V[vn, vn] %*% xrow))
  subf <- pr$functions[pr$functions$class == "B" &
                         pr$functions$group == disp[["0"]], ]
  expect_equal(subf$se.fit[i], se_hand, tolerance = 1e-9)

  # Non-baseline group ("1"): design row is [intercept, a11, g1 dummy].
  xrow1 <- c(1, basis$x[i, "a11"], 1)
  vn1 <- paste0("C:", c("(Intercept)", "a11", "g1"))
  se_hand1 <- sqrt(as.numeric(t(xrow1) %*% V[vn1, vn1] %*% xrow1))
  subf1 <- pr$functions[pr$functions$class == "C" &
                          pr$functions$group == disp[["1"]], ]
  expect_equal(subf1$se.fit[i], se_hand1, tolerance = 1e-9)
})

test_that("difference curves and their SEs are correct per logit", {
  skip_on_cran()
  obj <- mn_pred_fit(); fit <- obj$fit; dm <- obj$data
  fr <- fit$var_winners[["a"]]$fit
  im <- fr$test_results$interaction_model
  V  <- vcov(im$fit)
  pd <- mfp2:::mfpi_prepare_prediction_data(fit, "a", newdata = dm, grid = FALSE,
                                            n_grid = 200L, purpose = "function")
  basis <- mfp2:::mfpi_build_function_basis(fit, "a", fr,
                                            pd$cont_var_scaled, pd$x_display)
  pr <- predict(fit, newdata = dm, se.fit = TRUE)
  disp <- mfp2:::mfpi_prediction_group_display_labels(fit, c("0", "1"))
  names(disp) <- c("0", "1")

  # Point estimate: difference == curve(group1) - curve(group0), per logit.
  for (cc in c("B", "C")) {
    f1 <- pr$functions[pr$functions$class == cc & pr$functions$group == disp[["1"]], ]$fit
    f0 <- pr$functions[pr$functions$class == cc & pr$functions$group == disp[["0"]], ]$fit
    d  <- pr$differences[pr$differences$class == cc, ]$fit
    expect_equal(max(abs(d - (f1 - f0))), 0, tolerance = 1e-10)
  }

  # SE via delta method with the dense basis: gradient is -a01 (ref group0),
  # +a11 (group1), +g1 (dummy). Intercept cancels in the difference.
  i <- 5L
  dx <- c(-basis$x[i, "a01"], basis$x[i, "a11"], 1)
  vn <- paste0("C:", c("a01", "a11", "g1"))
  sed_hand <- sqrt(as.numeric(t(dx) %*% V[vn, vn] %*% dx))
  subd <- pr$differences[pr$differences$class == "C", ]
  expect_equal(subd$se.fit[i], sed_hand, tolerance = 1e-9)
})

test_that("type = 'function' and 'difference' return only the requested piece", {
  skip_on_cran()
  fit <- mn_pred_fit()$fit
  pf <- predict(fit, type = "function")
  expect_false(is.null(pf$functions))
  expect_null(pf$differences)
  pd <- predict(fit, type = "difference")
  expect_null(pd$functions)
  expect_false(is.null(pd$differences))
})

test_that("multinomial ordinary predictions return logits and probabilities", {
  skip_on_cran()
  obj <- mn_pred_fit(); fit <- obj$fit; dm <- obj$data
  im <- fit$var_winners[["a"]]$fit$test_results$interaction_model
  eta <- im$fit$mfp2_design %*% t(im$coefficient_matrix) +
    im$fit$mfp2_offset_matrix[, -1L, drop = FALSE]

  link <- predict(fit, type = "link")$predictions
  expect_equal(unname(as.matrix(link[paste0("fit.", colnames(eta))])),
               unname(eta), tolerance = 1e-8)
  se_columns <- grep("^se.fit\\.", names(link), value = TRUE)
  expect_length(se_columns, ncol(eta))
  expect_true(all(vapply(link[se_columns],
                         function(x) all(is.finite(x)), logical(1L))))

  eta_full <- cbind(0, eta)
  colnames(eta_full) <- im$fit$mfp2_family$prepared$levels
  exp_eta <- exp(eta_full - apply(eta_full, 1L, max))
  probability <- exp_eta / rowSums(exp_eta)
  response <- predict(fit, type = "response")$predictions
  expect_equal(unname(as.matrix(response[paste0("fit.", colnames(probability))])),
               unname(probability), tolerance = 1e-8)
  expect_equal(
    unname(rowSums(response[paste0("fit.", colnames(probability))])),
    rep(1, nrow(dm)),
    tolerance = 1e-10
  )

  response_new <- predict(fit, newdata = dm[1:7, ],
                          type = "response")$predictions
  expect_equal(response_new[paste0("fit.", colnames(probability))],
               response[1:7, paste0("fit.", colnames(probability))],
               tolerance = 1e-8)
})

test_that("se.fit = FALSE drops SE/CI columns but keeps per-logit fits", {
  skip_on_cran()
  fit <- mn_pred_fit()$fit
  pr <- predict(fit, se.fit = FALSE)
  expect_false("se.fit" %in% names(pr$functions))
  expect_true("class" %in% names(pr$functions))
  expect_true(all(is.finite(pr$functions$fit)))
})
