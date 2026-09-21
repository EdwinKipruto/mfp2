# BUG 1 fix: multinomial fits now use an nnet-native control object with a
# sensible default iteration cap (maxit = 300), instead of being routed through
# stats::glm.control() (maxit = 25). These tests pin the new control contract
# and the end-to-end behaviour it fixes.

# The normalizer is internal; reach it through the package namespace.
mn_control <- function(...) mfp2:::normalize_multinomial_control(...)

test_that("default multinomial control is nnet-native with maxit = 300", {
  skip_on_cran()
  ctl <- mn_control()
  expect_equal(ctl$maxit, 300L)
  # nnet-only field set: no glm.control() fields such as `epsilon`.
  expect_setequal(names(ctl), c("maxit", "reltol", "abstol", "trace", "MaxNWts"))
  expect_null(ctl$epsilon)
  # nnet's own defaults are preserved for the remaining knobs.
  expect_equal(ctl$reltol, 1.0e-8)
  expect_equal(ctl$abstol, 1.0e-4)
  expect_false(ctl$trace)
  expect_equal(ctl$MaxNWts, 1000L)
})

test_that("multinomial control does not depend on stats::glm.control()", {
  skip_on_cran()
  # The old default maxit came from glm.control() (25). The new default must be
  # decoupled from it: even if glm.control()'s maxit changed, mn_control() stays
  # at 300.
  expect_false(identical(mn_control()$maxit, stats::glm.control()$maxit))
})

test_that("user maxit is honoured and epsilon is accepted as a reltol alias", {
  skip_on_cran()
  expect_equal(mn_control(list(maxit = 50))$maxit, 50L)
  # glm-style `epsilon` maps onto nnet's `reltol` for backward-friendly control
  # lists, and does not survive as a stray field.
  ctl <- mn_control(list(epsilon = 1e-6))
  expect_equal(ctl$reltol, 1e-6)
  expect_null(ctl$epsilon)
})

test_that("unknown multinomial control fields are rejected", {
  skip_on_cran()
  expect_error(mn_control(list(foo = 1)), "Unknown multinomial control field")
  expect_error(mn_control(list(maxit = -1)), "positive integer")
  expect_error(mn_control(list(maxit = 2.5)), "positive integer")
  expect_error(mn_control(list(maxit = Inf)), "positive integer")
  expect_error(mn_control(list(MaxNWts = 1000.5)), "positive integer")
  expect_error(mn_control(list(MaxNWts = Inf)), "positive integer")
  expect_error(mn_control(list(trace = "yes")), "TRUE.*FALSE|trace")
})

test_that("multinomial FP selection succeeds out of the box (default controls)", {
  # This is the headline BUG 1 regression: with default controls an ordinary
  # multinomial FP model (df > 1) previously aborted as unconverged because the
  # candidate fast path inherited maxit = 25. It must now fit without any user
  # control tuning.
  skip_on_cran()
  dat <- make_multinomial_data(n = 500, seed = 33)
  expect_error(
    mfp2(y ~ fp(x1, df = 4) + fp(x2, df = 4), data = dat,
         family = multinomial_family(), criterion = "aic", verbose = FALSE),
    NA
  )
})

test_that("mfpi multinomial fits out of the box without a maxit override", {
  # Every mfpi() multinomial run used to require control = list(maxit = ...);
  # with the fixed default it works unaided.
  skip_on_cran()
  dm <- make_mfpi_multinomial()
  fit <- mfpi(y ~ g + fp(a), data = dm, group_var = "g",
              interaction_vars = "a", interaction_forms = c(a = "fp1"),
              family = multinomial_family(), p_interact = 0.1, verbose = FALSE)
  expect_s3_class(fit, "mfpi")
  expect_equal(fit$family_string, "multinomial")
})

test_that("raised default gives the same fit as an explicit maxit = 300", {
  skip_on_cran()
  dat <- make_multinomial_data(n = 500, seed = 41)
  f_default <- mfp2(y ~ fp(x1, df = 4) + fp(x2, df = 4), data = dat,
                    family = multinomial_family(), criterion = "aic",
                    verbose = FALSE)
  f_explicit <- mfp2(y ~ fp(x1, df = 4) + fp(x2, df = 4), data = dat,
                     family = multinomial_family(), criterion = "aic",
                     control = list(maxit = 300), verbose = FALSE)
  expect_equal(as.numeric(logLik(f_default)),
               as.numeric(logLik(f_explicit)), tolerance = TOL_NNET)
})
