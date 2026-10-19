# Zero-handled predictors (`zero`, `catzero`, `spike`):
# * automatic scale is 1, so log terms are not shifted on the positive rows
#   only and FP selection compares the same models that are finally fitted;
# * in mfpi(), centring follows mfp2(): the basis is 0 at structural zeros and
#   the mean of the complete column is subtracted from every row, so centring
#   never changes the fitted model.

make_zero_data <- function(n = 400, seed = 11) {
  set.seed(seed)
  d <- data.frame(
    g = factor(sample(c("A", "B", "C"), n, TRUE)),
    w = rexp(n) * 30,            # automatic scale would be 100
    x = runif(n, 0.5, 80),       # ordinary predictor, automatic scale 10
    z = runif(n, 1, 5)
  )
  d$w[sample(n, 80)] <- 0
  d$y <- ifelse(d$w > 0, log(d$w), 0) * (1 + (d$g == "B")) / 3 +
    0.2 * d$z + rnorm(n)
  d
}

test_that("fix_zero_handled_scale() only resets automatic scales of zero-handled predictors with zeros", {
  x <- cbind(w = c(0, 1, 50), v = c(1, 2, 3), u = c(0, 2, 3), p = c(1, 5, 9))
  scale <- c(w = NA, v = NA, u = 10, p = NA)
  zero <- c(w = TRUE, v = TRUE, u = TRUE, p = FALSE)
  out <- mfp2:::fix_zero_handled_scale(scale, zero, x)
  expect_identical(out[["w"]], 1)         # zero-handled, has zeros, automatic
  expect_true(is.na(out[["v"]]))          # zero-handled but no exact zeros
  expect_identical(out[["u"]], 10)        # user-supplied scale is kept
  expect_true(is.na(out[["p"]]))          # ordinary predictor untouched
})

test_that("mfp2() does not scale zero-handled predictors automatically", {
  d <- make_zero_data()
  auto <- mfp2(y ~ fp(w, zero = TRUE) + fp(x) + fp(z), data = d, verbose = FALSE)
  expect_identical(auto$transformations["w", "scale"], 1)
  expect_true(auto$transformations["x", "scale"] != 1)   # ordinary: still scaled

  explicit <- mfp2(y ~ fp(w, zero = TRUE, scale = 1) + fp(x) + fp(z),
                   data = d, verbose = FALSE)
  expect_equal(auto$fp_powers, explicit$fp_powers)
  expect_equal(deviance(auto), deviance(explicit))

  user <- mfp2(y ~ fp(w, zero = TRUE, scale = 100) + fp(z), data = d,
               verbose = FALSE)
  expect_identical(user$transformations["w", "scale"], 100)
})

test_that("mfpi() does not scale zero-handled predictors automatically", {
  d <- make_zero_data()
  fit_auto <- suppressWarnings(mfpi(y ~ g + fp(w, zero = TRUE) + fp(z), data = d,
                                    group_var = "g", interaction_vars = "w",
                                    p_interact = 1, verbose = FALSE))
  fit_one <- suppressWarnings(mfpi(y ~ g + fp(w, zero = TRUE, scale = 1) + fp(z),
                                   data = d, group_var = "g",
                                   interaction_vars = "w", p_interact = 1,
                                   verbose = FALSE))
  expect_identical(fit_auto$scale[["w"]], 1)
  expect_equal(fit_auto$all_model_metrics$pvalue,
               fit_one$all_model_metrics$pvalue)
})

test_that("centring does not change zero-handled mfpi interaction results", {
  d <- make_zero_data()
  for (flex in c("flex1", "flex2", "flex3", "flex4")) {
    for (form in c("fp1", "fp2")) {
      fit_with <- function(center, center_type) {
        suppressWarnings(mfpi(
          stats::as.formula(paste0("y ~ g + fp(w, zero = TRUE, center = ",
                                   center, ") + fp(z)")),
          data = d, group_var = "g", interaction_vars = "w",
          interaction_forms = c(w = form), flex = flex,
          center_type = center_type, p_interact = 1, verbose = FALSE
        ))
      }
      ref <- fit_with(FALSE, "grand")
      for (variant in list(c("TRUE", "grand"), c("TRUE", "group"))) {
        fit <- fit_with(variant[1], variant[2])
        info <- paste(flex, form, variant[2])
        expect_equal(fit$all_model_metrics$deviance_int,
                     ref$all_model_metrics$deviance_int, tolerance = 1e-8,
                     info = info)
        expect_equal(fit$all_model_metrics$pvalue,
                     ref$all_model_metrics$pvalue, tolerance = 1e-8, info = info)
        expect_equal(fit$all_model_metrics$fp_powers_int,
                     ref$all_model_metrics$fp_powers_int, info = info)
        p_fit <- predict(fit, newdata = d, type = "link", model = "all",
                         se.fit = FALSE)$predictions$fit
        p_ref <- predict(ref, newdata = d, type = "link", model = "all",
                         se.fit = FALSE)$predictions$fit
        expect_equal(p_fit, p_ref, tolerance = 1e-8, info = info)
      }
    }
  }
})

test_that("structural zeros are centred like every other row in the interaction design", {
  d <- make_zero_data()
  for (ct in c("grand", "group")) {
    fit <- suppressWarnings(mfpi(y ~ g + fp(w, zero = TRUE) + fp(z), data = d,
                                 group_var = "g", interaction_vars = "w",
                                 flex = "flex2", center_type = ct,
                                 p_interact = 1, verbose = FALSE))
    res <- fit$var_winners$w$fit
    group_a <- res$coefficient_groups[[1L]]
    zero_a <- d$w == 0 & d$g == "A"
    expect_equal(unname(res$xinteraction[zero_a, group_a[1L]]),
                 rep(-unname(res$center_vals[group_a[1L]]), sum(zero_a)),
                 info = ct)
  }
})

test_that("fitted functions and differences do not depend on centring for zero variables", {
  d <- make_zero_data()
  fits <- lapply(c(TRUE, FALSE), function(cen) {
    suppressWarnings(mfpi(
      stats::as.formula(paste0("y ~ g + fp(w, zero = TRUE, center = ", cen,
                               ") + fp(z)")),
      data = d, group_var = "g", interaction_vars = "w", flex = "flex3",
      interaction_forms = c(w = "fp2"), p_interact = 1, verbose = FALSE
    ))
  })
  p <- lapply(fits, function(f) predict(f, model = "all", type = "both",
                                        newdata = data.frame(w = c(0, 1, 10, 100))))
  expect_equal(p[[1]]$functions$fit, p[[2]]$functions$fit, tolerance = 1e-8)
  expect_equal(p[[1]]$functions$se.fit, p[[2]]$functions$se.fit, tolerance = 1e-8)
  expect_equal(p[[1]]$differences$fit, p[[2]]$differences$fit, tolerance = 1e-8)
})
