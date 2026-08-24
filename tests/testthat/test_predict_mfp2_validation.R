test_that("predict.mfp2 validates scalar prediction controls", {
  dummy <- structure(list(), class = "mfp2")

  bad_alpha <- list(
    NA_real_, Inf, -Inf, 0, 1, -0.1, 1.1,
    c(0.05, 0.10), numeric(0), "0.05"
  )
  for (value in bad_alpha) {
    expect_error(
      stats::predict(dummy, alpha = value),
      "alpha"
    )
  }

  bad_nseq <- list(
    NA_real_, Inf, -Inf, 0, -1, 2.5,
    c(10, 20), numeric(0), 1e20, "100"
  )
  for (value in bad_nseq) {
    expect_error(
      stats::predict(dummy, nseq = value),
      "nseq"
    )
  }

  bad_logical <- list(NA, c(TRUE, FALSE), logical(0), 1L, "TRUE", NULL)
  for (value in bad_logical) {
    expect_error(
      stats::predict(dummy, se.fit = value),
      "se.fit"
    )
    expect_error(
      stats::predict(dummy, add_intercept = value),
      "add_intercept"
    )
  }
})


test_that("shared scalar prediction validators accept valid boundary-safe values", {
  expect_invisible(mfp2:::validate_open_probability_scalar(0.05, "alpha"))
  expect_invisible(mfp2:::validate_open_probability_scalar(.Machine$double.eps, "alpha"))
  expect_invisible(mfp2:::validate_open_probability_scalar(1 - .Machine$double.eps, "alpha"))

  expect_invisible(mfp2:::validate_positive_integer_scalar(1, "nseq"))
  expect_invisible(mfp2:::validate_positive_integer_scalar(100L, "nseq"))

  expect_invisible(
    mfp2:::validate_logical_vector(TRUE, "se.fit", allowed_lengths = 1L)
  )
  expect_invisible(
    mfp2:::validate_logical_vector(FALSE, "add_intercept", allowed_lengths = 1L)
  )
})


# Migrated coverage from the former test_mfp2.R

# =============================================================================
# End of tests
# =============================================================================

# =============================================================================
# Prediction API regression tests maintained for 1.1.0.9003
# =============================================================================

test_that("predict.mfp2 normalizes linear-predictor aliases", {
  data("prostate", package = "mfp2")
  fit_glm <- mfp2(
    lpsa ~ fp(age) + svi,
    data = prostate,
    verbose = FALSE
  )
  nd_glm <- prostate[1:12, c("age", "svi"), drop = FALSE]
  expect_equal(
    predict(fit_glm, newdata = nd_glm, type = "lp"),
    predict(fit_glm, newdata = nd_glm, type = "link")
  )
  expect_equal(
    predict(fit_glm, newdata = nd_glm, type = c(alias = "lp")),
    predict(fit_glm, newdata = nd_glm, type = "link")
  )

  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[complete.cases(dat[, c("time", "status", "age", "sex")]), ]
  fit_cox <- mfp2(
    survival::Surv(time, status) ~ fp(age) + sex,
    data = dat,
    family = "cox",
    verbose = FALSE
  )
  nd_cox <- dat[1:12, c("age", "sex"), drop = FALSE]
  expect_equal(
    predict(fit_cox, newdata = nd_cox, type = "link", cox_reference = "zero"),
    predict(fit_cox, newdata = nd_cox, type = "lp", cox_reference = "zero")
  )
  expect_equal(
    predict(
      fit_cox,
      newdata = nd_cox,
      type = c(alias = "link"),
      cox_reference = "zero"
    ),
    predict(fit_cox, newdata = nd_cox, type = "lp", cox_reference = "zero")
  )
})


test_that("predict.mfp2 rejects infinite numeric newdata columns", {
  data("prostate", package = "mfp2")
  fit <- mfp2(lpsa ~ fp(age) + svi, data = prostate, verbose = FALSE)

  nd_inf <- prostate[1:6, c("age", "svi"), drop = FALSE]
  nd_inf$age[2] <- Inf
  expect_error(
    predict(fit, newdata = nd_inf, type = "link"),
    "finite numeric values|Infinite"
  )

  nd_ninf <- prostate[1:6, c("age", "svi"), drop = FALSE]
  nd_ninf$age[2] <- -Inf
  expect_error(
    predict(fit, newdata = nd_ninf, type = "link"),
    "finite numeric values|Infinite"
  )
})


# =============================================================================
# Prediction argument and non-finite validation regressions
# =============================================================================

test_that("predict.mfp2 validates only required missing predictors", {
  data("prostate", package = "mfp2")
  fit <- mfp2(
    lpsa ~ fp(age) + svi,
    data = prostate,
    keep = "age",
    select = 1,
    alpha = 1,
    verbose = FALSE
  )
  nd <- prostate[1:8, c("age", "svi"), drop = FALSE]
  expected <- predict(fit, newdata = nd, type = "link")

  nd_extra <- nd
  nd_extra$unused_na <- NA_real_
  nd_extra$unused_nan <- NaN
  expect_equal(
    predict(fit, newdata = nd_extra, type = "link"),
    expected,
    tolerance = 1e-12
  )

  for (bad_value in list(NA_real_, NaN, Inf, -Inf)) {
    nd_bad <- nd
    nd_bad$age[2] <- bad_value
    expect_error(
      predict(fit, newdata = nd_bad, type = "link"),
      "missing|finite"
    )
  }

  fit_matrix <- mfp2(
    x = as.matrix(prostate[, c("age", "svi")]),
    y = prostate$lpsa,
    verbose = FALSE
  )
  nd_matrix <- data.frame(
    age = prostate$age[1:8],
    svi = prostate$svi[1:8],
    unused_na = NA_real_,
    unused_nan = NaN
  )
  expect_equal(
    predict(fit_matrix, newdata = nd_matrix, type = "link"),
    predict(
      fit_matrix,
      newdata = nd_matrix[, c("age", "svi"), drop = FALSE],
      type = "link"
    ),
    tolerance = 1e-12
  )
})


test_that("predict.mfp2 rejects irrelevant strata arguments", {
  data("prostate", package = "mfp2")
  fit_glm <- mfp2(lpsa ~ fp(age) + svi, data = prostate, verbose = FALSE)
  nd_glm <- prostate[1:6, c("age", "svi"), drop = FALSE]
  expect_error(
    predict(fit_glm, newdata = nd_glm, strata = rep(1, nrow(nd_glm))),
    "only for Cox"
  )

  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[complete.cases(dat[, c("time", "status", "age", "sex")]), ]
  fit_cox <- mfp2(
    survival::Surv(time, status) ~ age + sex,
    data = dat,
    family = "cox",
    df = 1,
    select = 1,
    alpha = 1,
    verbose = FALSE
  )
  nd_cox <- dat[1:8, c("age", "sex"), drop = FALSE]

  expect_error(
    predict(
      fit_cox,
      newdata = nd_cox,
      type = "lp",
      strata = rep(1, nrow(nd_cox))
    ),
    "not stratified"
  )
  expect_error(
    predict(
      fit_cox,
      newdata = nd_cox,
      type = "terms",
      strata = rep(1, nrow(nd_cox))
    ),
    "not used.*term or contrast"
  )
  expect_error(
    predict(fit_cox, type = "lp", strata = rep(1, nrow(dat))),
    "together with.*newdata"
  )

  dat$stratum <- factor(rep(c("A", "B"), length.out = nrow(dat)))
  fit_stratified <- mfp2(
    survival::Surv(time, status) ~ age + sex + strata(stratum),
    data = dat,
    family = "cox",
    df = 1,
    select = 1,
    alpha = 1,
    verbose = FALSE
  )
  nd_stratified <- dat[1:8, c("age", "sex", "stratum"), drop = FALSE]
  bad_strata <- as.numeric(nd_stratified$stratum)
  bad_strata[2] <- Inf
  expect_error(
    predict(
      fit_stratified,
      newdata = nd_stratified,
      type = "lp",
      strata = bad_strata
    ),
    "strata.*finite"
  )
})


test_that("predict.mfp2 rejects irrelevant and non-finite replacement offsets", {
  data("prostate", package = "mfp2")
  fit_plain <- mfp2(lpsa ~ fp(age) + svi, data = prostate, verbose = FALSE)
  nd_plain <- prostate[1:7, c("age", "svi"), drop = FALSE]
  expect_error(
    predict(
      fit_plain,
      newdata = nd_plain,
      type = "link",
      newoffset = rep(0, nrow(nd_plain))
    ),
    "used an offset"
  )

  set.seed(21001)
  dat <- prostate
  dat$exposure <- stats::runif(nrow(dat), 0.7, 2.2)
  fit_offset <- mfp2(
    lpsa ~ fp(age) + svi + offset(log(exposure)),
    data = dat,
    verbose = FALSE
  )
  nd <- dat[1:7, c("age", "svi", "exposure"), drop = FALSE]
  good_offset <- log(nd$exposure)

  expect_error(
    predict(
      fit_offset,
      newdata = nd,
      type = "terms",
      newoffset = good_offset
    ),
    "not used.*term or contrast"
  )
  expect_error(
    predict(fit_offset, type = "link", newoffset = log(dat$exposure)),
    "together with.*newdata"
  )

  for (bad_value in list(NA_real_, NaN, Inf, -Inf)) {
    bad_offset <- good_offset
    bad_offset[2] <- bad_value
    expect_error(
      predict(
        fit_offset,
        newdata = nd,
        type = "link",
        newoffset = bad_offset
      ),
      "finite numeric"
    )
  }
})
