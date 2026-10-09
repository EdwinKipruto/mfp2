test_that("pre-scaled GEE FP candidates preserve candidate fits and scores", {
  skip_if_not_installed("geepack")
  set.seed(9024)
  id <- rep(seq_len(30L), each = 3L)
  n <- length(id)
  x <- seq(0.5, 2.4, length.out = n) + stats::runif(n, 0, 0.15)
  z <- stats::rnorm(n)
  y <- 1 + 0.5 * log(x) + 0.3 * z +
    stats::rnorm(30L, sd = 0.25)[id] + stats::rnorm(n, sd = 0.3)
  family <- prepare_gee_family(
    gee_family(stats::gaussian(), corstr = "exchangeable"),
    y = y, id = id
  )

  basis <- cbind(x^(-2), log(x), x, x * log(x))
  template <- cbind("(Intercept)" = 1, V1 = 0, V2 = 0, z = z)
  candidate_map <- rbind(c(1L, 2L), c(3L, 4L))
  prepared <- prepare_gee_candidate_scaling(template, basis)

  for (criterion in c("aic", "pvalue")) {
    for (i in seq_len(nrow(candidate_map))) {
      source <- as.integer(candidate_map[i, ])
      raw <- template
      raw[, 2:3] <- basis[, source]
      old <- fit_gee(
        x = raw, y = y, family = family, x_has_intercept = TRUE,
        selection_criterion = criterion
      )

      scaled <- copy_fp_basis_candidate_cpp(
        target = prepared$design, basis = prepared$basis,
        source_cols = source, target_cols = 2:3
      )
      scales <- prepared$column_scales
      scales[2:3] <- prepared$basis_scales[source]
      expect_equal(scaled, sweep(raw, 2L, scales, "/"), tolerance = 1e-12)

      new <- fit_gee(
        x = scaled, y = y, family = family, x_has_intercept = TRUE,
        column_scales = scales, selection_criterion = criterion
      )
      expect_equal(new$logl, old$logl, tolerance = 1e-7)
      expect_equal(new$selection_deviance, old$selection_deviance,
                   tolerance = 1e-7)
      expect_equal(new$coefficients, old$coefficients, tolerance = 1e-7)
      expect_equal(new$robust_vcov, old$robust_vcov, tolerance = 1e-7)
      expect_equal(
        calculate_gee_metrics(new, n_obs = 30L, df_additional = 2),
        calculate_gee_metrics(old, n_obs = 30L, df_additional = 2),
        tolerance = 1e-7
      )
    }
  }
})

test_that("C++ GEE preparation matches R column scaling and leaves inputs intact", {
  design <- cbind("(Intercept)" = c(1, 1, 1),
                  focal = c(0, 0, 0), adjustment = c(-12, 3, 6))
  basis <- cbind(power = c(-20, 0, 10),
                 nonfinite = c(1, Inf, 2), missing = c(NA, 2, 3))
  original_design <- design
  original_basis <- basis

  expected_scales <- function(x) {
    apply(x, 2L, function(col) {
      s <- max(abs(col))
      if (!is.finite(s) || s <= 0) 1 else s
    })
  }
  design_scales <- expected_scales(design)
  design_scales["(Intercept)"] <- 1
  basis_scales <- expected_scales(basis)

  prepared <- prepare_gee_candidate_scaling(design, basis)
  expect_equal(prepared$column_scales, design_scales)
  expect_equal(prepared$basis_scales, basis_scales)
  expect_equal(prepared$design, sweep(design, 2L, design_scales, "/"))
  expect_equal(prepared$basis, sweep(basis, 2L, basis_scales, "/"))
  expect_identical(design, original_design)
  expect_identical(basis, original_basis)
})
