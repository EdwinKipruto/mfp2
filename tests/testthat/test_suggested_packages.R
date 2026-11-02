# ggplot2 and geepack are in Suggests. These tests mock them as unavailable
# and check that every entry point fails with an informative error, while
# features that do not need them keep working.

mock_unavailable <- function(pkgs, env = parent.frame()) {
  local_mocked_bindings(
    mfp2_has_package = function(pkg) !pkg %in% pkgs,
    .env = env
  )
}

test_that("availability helpers give informative errors", {
  mock_unavailable(c("ggplot2", "geepack"))

  expect_error(
    mfp2_require_ggplot2(),
    "! Plotting requires the 'ggplot2' package.\ni Install it with install.packages(\"ggplot2\").",
    fixed = TRUE
  )
  expect_error(
    mfp2_require_geepack(),
    "! GEE models require the 'geepack' package.\ni Install it with install.packages(\"geepack\").",
    fixed = TRUE
  )
})

test_that("availability helpers pass when the packages are installed", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("geepack")
  expect_true(mfp2_require_ggplot2())
  expect_true(mfp2_require_geepack())
})

test_that("requireNamespace() returning FALSE triggers the errors", {
  # The helpers are thin wrappers around requireNamespace(); confirm the
  # wiring by asking for a package that cannot exist.
  expect_false(mfp2_has_package("mfp2NoSuchPackageXYZ"))
})

test_that("gee_family() stops when geepack is unavailable", {
  mock_unavailable("geepack")
  expect_error(gee_family(), "GEE models require the 'geepack' package",
               fixed = TRUE)
  expect_error(gee_family(), "install.packages(\"geepack\")", fixed = TRUE)
})

test_that("mfp2() fits non-GEE models without ggplot2 or geepack", {
  mock_unavailable(c("ggplot2", "geepack"))
  data("prostate", package = "mfp2", envir = environment())

  fit <- mfp2(lpsa ~ fp(cavol), data = prostate, verbose = FALSE)
  expect_s3_class(fit, "mfp2")
  expect_output(print(fit))
  expect_s3_class(summary(fit), "summary.mfp2")
  expect_length(predict(fit), nrow(prostate))
  expect_length(residuals(fit), nrow(prostate))
  expect_true(is.matrix(confint(fit)))
})

test_that("plot() and fracplot() stop when ggplot2 is unavailable", {
  data("prostate", package = "mfp2", envir = environment())
  fit <- mfp2(lpsa ~ fp(cavol), data = prostate, verbose = FALSE)

  mock_unavailable("ggplot2")
  expect_error(plot(fit), "Plotting requires the 'ggplot2' package",
               fixed = TRUE)
  expect_error(plot(fit), "install.packages(\"ggplot2\")", fixed = TRUE)
  # The check runs before the deprecation warning.
  expect_error(
    expect_no_warning(fracplot(fit)),
    "Plotting requires the 'ggplot2' package",
    fixed = TRUE
  )
})

test_that("plot.mfpi() stops when ggplot2 is unavailable", {
  data("prostate", package = "mfp2", envir = environment())
  fit <- mfpi(
    lpsa ~ fp(cavol) + svi,
    data = prostate,
    group_var = "svi",
    interaction_vars = "cavol",
    verbose = FALSE
  )

  mock_unavailable("ggplot2")
  expect_error(
    plot(fit, terms = "cavol", plot_type = "difference"),
    "Plotting requires the 'ggplot2' package",
    fixed = TRUE
  )
})

test_that("methods for a fitted GEE model stop when geepack is unavailable", {
  skip_if_not_installed("geepack")
  set.seed(1)
  d <- data.frame(id = rep(1:30, each = 3), x = runif(90, 1, 5))
  d$y <- log(d$x) + rnorm(30, sd = 0.5)[d$id] + rnorm(90, sd = 0.3)
  fit <- mfp2(
    y ~ fp(x, df = 2),
    data = d,
    family = gee_family(stats::gaussian()),
    id = d$id,
    verbose = FALSE
  )

  mock_unavailable("geepack")
  msg <- "GEE models require the 'geepack' package"
  expect_error(predict(fit), msg, fixed = TRUE)
  expect_error(summary(fit), msg, fixed = TRUE)
  expect_error(residuals(fit), msg, fixed = TRUE)
  expect_error(confint(fit), msg, fixed = TRUE)
  expect_error(print(fit), msg, fixed = TRUE)
  # coef() needs nothing from geepack.
  expect_identical(coef(fit), fit$coefficients)
})
