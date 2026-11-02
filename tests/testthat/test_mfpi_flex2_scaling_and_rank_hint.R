# flex2 selects FP powers on the scaled covariate with uncentred candidates,
# like flex1, flex4 and mfp2(), and builds its final design once from the
# selected powers. Rank-deficient designs raise a classed error that explains
# the typical automatic-shift cause.

make_flex2_data <- function(n = 400, seed = 11) {
  set.seed(seed)
  d <- data.frame(
    g = factor(sample(c("A", "B", "C"), n, TRUE)),
    x = runif(n, 0.5, 80),   # automatic scale = 10
    z = runif(n, 1, 5)
  )
  d$y <- log(d$x) * (1 + (d$g == "B")) / 2 + 0.2 * d$z + rnorm(n)
  d$yb <- rbinom(n, 1, plogis(d$y - 2))
  d
}

flex2_fit <- function(d, lhs = "y", family = "gaussian", x_term = "fp(x)", ...) {
  suppressWarnings(mfpi(
    stats::as.formula(paste(lhs, "~ g +", x_term, "+ fp(z)")),
    data = d, family = family, group_var = "g", interaction_vars = "x",
    p_interact = 1, verbose = FALSE, ...
  ))
}

test_that("flex2/flex3 selection is invariant to scaling for ordinary variables", {
  d <- make_flex2_data()
  for (flex in c("flex2", "flex3")) {
    for (form in c("fp1", "fp2")) {
      for (fam in c("gaussian", "binomial")) {
        lhs <- if (fam == "gaussian") "y" else "yb"
        auto <- flex2_fit(d, lhs, fam, flex = flex, interaction_forms = c(x = form))
        unscaled <- flex2_fit(d, lhs, fam, x_term = "fp(x, scale = 1)",
                              flex = flex, interaction_forms = c(x = form))
        info <- paste(flex, form, fam)
        expect_true(auto$scale[["x"]] != 1, info = info)
        expect_equal(auto$all_model_metrics$fp_powers_int,
                     unscaled$all_model_metrics$fp_powers_int, info = info)
        expect_equal(auto$all_model_metrics$pvalue,
                     unscaled$all_model_metrics$pvalue, tolerance = 1e-8,
                     info = info)
      }
    }
  }
})

test_that("flex2 selection does not depend on the centring type for ordinary variables", {
  d <- make_flex2_data()
  for (form in c("fp1", "fp2")) {
    grand <- flex2_fit(d, flex = "flex2", interaction_forms = c(x = form),
                       center_type = "grand")
    group <- flex2_fit(d, flex = "flex2", interaction_forms = c(x = form),
                       center_type = "group")
    expect_equal(grand$all_model_metrics$fp_powers_int,
                 group$all_model_metrics$fp_powers_int, info = form)
    expect_equal(grand$all_model_metrics$pvalue,
                 group$all_model_metrics$pvalue, tolerance = 1e-8, info = form)
  }
})

test_that("flex2 final design matches the shared create_z_variables() builder", {
  d <- make_flex2_data()
  for (ct in c("grand", "group")) {
    fit <- flex2_fit(d, flex = "flex2", interaction_forms = c(x = "fp2"),
                     center_type = ct)
    res <- fit$var_winners$x$fit
    x_int <- fit$x_train_internal
    ref <- mfp2:::create_z_variables(
      cont_var = x_int[, "x", drop = FALSE],
      group_var = x_int[, "g", drop = FALSE],
      power = res$bestfp_main, shift = 0, scale = 1, center = TRUE,
      center_type = ct, scale_var = fit$scale[["x"]]
    )
    focal <- unlist(res$coefficient_groups, use.names = FALSE)
    expect_equal(unname(res$xinteraction[, focal]), unname(ref$z[, focal]),
                 tolerance = 1e-10, info = ct)
    expect_equal(res$center_vals, ref$center_vals, tolerance = 1e-10, info = ct)
  }
})

test_that("predictions from newdata match training predictions after the change", {
  d <- make_flex2_data()
  fit <- flex2_fit(d, flex = "flex2", interaction_forms = c(x = "fp2"))
  p_new <- predict(fit, newdata = d, type = "link", model = "all",
                   se.fit = FALSE)$predictions$fit
  p_train <- predict(fit, type = "link", model = "all",
                     se.fit = FALSE)$predictions$fit
  expect_equal(p_new, p_train, tolerance = 1e-10)
})

# Rank-deficiency errors -------------------------------------------------------

make_rank_data <- function() {
  set.seed(1)
  n <- 300
  d <- data.frame(g = rbinom(n, 1, 0.5), z = runif(n, 1, 5))
  d$x3 <- sample(c(-1, 0, 1), n, TRUE)
  d$xn <- rnorm(n)
  d$w <- rexp(n)
  d$w[sample(n, 60)] <- 0
  d$y <- d$xn * (1 + d$g) + d$x3 + log(d$z) + rnorm(n)
  d
}

test_that("automatic-shift rank failure gives a classed error with guidance", {
  d <- make_rank_data()
  err <- tryCatch(
    mfpi(y ~ g + fp(xn) + fp(z), data = d, group_var = "g",
         interaction_vars = "xn", interaction_forms = c(xn = "fp2"),
         flex = "flex1", verbose = FALSE),
    error = function(e) e
  )
  expect_s3_class(err, "mfpi_rank_deficient_error")
  msg <- conditionMessage(err)
  expect_match(msg, "rank deficient")
  expect_match(msg, "shifted automatically", fixed = TRUE)
  expect_match(msg, "fp(xn, shift = ...)", fixed = TRUE)
  expect_match(msg, "interaction_forms = c(xn = \"fp1\")", fixed = TRUE)
  expect_match(msg, "another `flex` level", fixed = TRUE)
})

test_that("a user-supplied shift never gets the automatic-shift guidance", {
  d <- make_rank_data()
  s_auto <- mfp2:::find_shift_factor(d$xn)
  err <- tryCatch(
    mfpi(y ~ g + fp(xn, shift = s_auto) + fp(z), data = d, group_var = "g",
         interaction_vars = "xn", interaction_forms = c(xn = "fp2"),
         flex = "flex1", verbose = FALSE),
    error = function(e) e
  )
  expect_s3_class(err, "mfpi_rank_deficient_error")
  expect_false(grepl("shifted automatically", conditionMessage(err), fixed = TRUE))
})

test_that("guidance is only produced for the automatic-shift signature", {
  spec <- list(kind = "continuous")
  x <- cbind(xn = c(1e-6, 1, 2, 3, 5))          # internal values (scale 1)
  args <- list(var_name = "xn", interaction_spec = spec,
               interaction_form = "fp2", x = x, scale = c(xn = 1),
               shift = c(xn = 2.5), shift_automatic = c(xn = TRUE))
  hint <- function(...) do.call(mfp2:::mfpi_rank_deficient_hint,
                                utils::modifyList(args, list(...)))
  expect_type(hint(powers = c(-2, 1)), "character")
  expect_null(hint(powers = c(1, 2)))                      # no negative/log power
  expect_null(hint(powers = c(-2, 1), shift_automatic = c(xn = FALSE)))
  expect_null(hint(powers = c(-2, 1), shift = c(xn = 0)))  # not shifted
  expect_null(hint(powers = c(-2, 1), interaction_spec = list(kind = "binary")))
  expect_null(hint(powers = c(-2, 1), x = cbind(xn = c(1, 2, 3, 5))))  # min not near 0
  expect_match(hint(powers = c(0)), "c(xn = \"fp1\")", fixed = TRUE)
  expect_match(hint(powers = c(0), interaction_form = "fp1"),
               "c(xn = \"linear\")", fixed = TRUE)
})

test_that("rank-deficiency helper signals the documented class", {
  err <- tryCatch(mfp2:::mfpi_stop_rank_deficient("a", "b"), error = function(e) e)
  expect_s3_class(err, "mfpi_rank_deficient_error")
  expect_identical(conditionMessage(err), "ab")
})
