# summary.mfp2() keeps the components that summary() returned in mfp2 1.0.x
# (summary.glm() for GLM-type models), so that existing code reading them
# keeps working. They are not printed.

test_that("GLM summaries expose summary.glm() components", {
  data("prostate", package = "mfp2", envir = environment())
  fit <- mfp2(lpsa ~ fp(age) + svi, data = prostate, verbose = FALSE)
  s <- summary(fit)

  base <- fit
  class(base) <- setdiff(class(base), "mfp2")
  ref <- summary(base)

  expect_identical(s$coefficients, ref$coefficients)
  expect_identical(s$cov.scaled, ref$cov.scaled)
  expect_identical(s$cov.unscaled, ref$cov.unscaled)
  # Partial matching used by older code still works.
  expect_identical(s$coef, ref$coefficients)
  # The covariance agrees with vcov().
  expect_equal(s$cov.scaled, stats::vcov(fit))
})

test_that("compatibility components do not change the printed summary", {
  data("prostate", package = "mfp2", envir = environment())
  fit <- mfp2(lpsa ~ fp(age) + svi, data = prostate, verbose = FALSE)
  s <- summary(fit)
  trimmed <- s
  trimmed[c("coefficients", "cov.scaled", "cov.unscaled")] <- NULL
  class(trimmed) <- class(s)
  expect_identical(capture.output(print(s)), capture.output(print(trimmed)))
})

test_that("compatibility helper copies only components that exist", {
  expect_identical(mfp2_summary_compat_components(NULL), list())
  raw <- structure(list(coefficients = 1, other = 2), class = "summary.foo")
  expect_identical(mfp2_summary_compat_components(raw), list(coefficients = 1))
})
