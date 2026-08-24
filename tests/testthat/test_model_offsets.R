# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# Test purpose: Verifies that negative-binomial weights and offsets are passed
# consistently through fitting and prediction by comparison with glm.nb().
test_that("negative-binomial weights and offsets agree with MASS::glm.nb()", {
  skip_if_not_installed("fastglm")
  skip_if_not_installed("MASS")

  set.seed(1201)
  n <- 280L
  dat <- data.frame(
    x1 = stats::runif(n, 0.5, 2.5),
    x2 = stats::rnorm(n),
    exposure = stats::runif(n, 0.6, 2.2),
    w = sample(c(1, 2), n, replace = TRUE)
  )
  dat$log_exposure <- log(dat$exposure)
  dat$y <- stats::rnbinom(
    n,
    mu = exp(0.15 + 0.4 * dat$x1 - 0.2 * dat$x2 + dat$log_exposure),
    size = 3.2
  )

  fit_mfp2 <- mfp2(
    y ~ x1 + x2,
    data = dat,
    family = "negbin",
    fitter = "fastglm",
    weights = dat$w,
    offset = dat$log_exposure,
    cycles = 1,
    df = 1,
    select = 1,
    alpha = 1,
    xorder = "original",
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = FALSE
  )
  fit_mass <- MASS::glm.nb(
    y ~ x1 + x2 + offset(log_exposure),
    data = dat,
    weights = w,
    link = log,
    control = stats::glm.control(maxit = 100L)
  )

  expect_negbin_mfp2_mass_equal(
    fit_mfp2,
    fit_mass,
    coefficient_tolerance = 5e-4,
    likelihood_tolerance = 2e-3
  )
  expect_equal(unname(fit_mfp2$prior.weights), dat$w)

  nd <- dat[seq_len(20L), c("x1", "x2", "log_exposure"), drop = FALSE]
  pred_mfp2 <- predict(
    fit_mfp2,
    newdata = nd[, c("x1", "x2"), drop = FALSE],
    newoffset = nd$log_exposure,
    type = "link",
    se.fit = TRUE
  )
  pred_mass <- stats::predict(
    fit_mass,
    newdata = nd,
    type = "link",
    se.fit = TRUE
  )

  expect_equal(unname(pred_mfp2$fit), unname(pred_mass$fit), tolerance = 5e-4)
  expect_equal(
    unname(pred_mfp2$se.fit),
    unname(pred_mass$se.fit),
    tolerance = 5e-4
  )

  shifted <- predict(
    fit_mfp2,
    newdata = nd[, c("x1", "x2"), drop = FALSE],
    newoffset = nd$log_exposure + 0.35,
    type = "link"
  )
  expect_equal(
    unname(shifted - pred_mfp2$fit),
    rep(0.35, nrow(nd)),
    tolerance = 1e-10
  )
})
