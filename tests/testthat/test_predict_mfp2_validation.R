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


test_that("predict.mfp2 rejects infinite required predictor columns", {
  data("prostate", package = "mfp2")
  fit <- mfp2(lpsa ~ fp(age) + svi, data = prostate, keep = "age",
              select = 1, alpha = 1, verbose = FALSE)

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

test_that("predict.mfp2 ignores missing and infinite unused columns", {
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
  nd_extra$unused_inf <- Inf
  nd_extra$unused_ninf <- -Inf
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
    unused_nan = NaN,
    unused_inf = Inf,
    unused_ninf = -Inf
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
  expect_equal(
    predict(fit_matrix, newdata = as.matrix(nd_matrix), type = "link"),
    predict(fit_matrix, newdata = as.matrix(nd_matrix[, c("age", "svi")]),
            type = "link"),
    tolerance = 1e-12
  )
})


test_that("term predictions ignore unrequested predictors and formula inputs", {
  dat <- data.frame(x = seq(1, 5, length.out = 90L),
                    z = rep(seq(1, 3, length.out = 30L), 3L))
  dat$y <- 0.4 * dat$x + 0.8 / dat$z + 0.1 * sin(seq_len(nrow(dat)))
  fit_formula <- mfp2(y ~ x + I(1 / z), data = dat, df = 1,
                      select = 1, alpha = 1, shift = 0, scale = 1,
                      center = FALSE, verbose = FALSE)
  fit_matrix <- mfp2(as.matrix(dat[, c("x", "z")]), dat$y, df = 1,
                     select = 1, alpha = 1, shift = 0, scale = 1,
                     center = FALSE, verbose = FALSE)
  for (fit in list(fit_formula, fit_matrix)) {
    for (type in c("terms", "contrasts")) {
      prediction <- function(nd) {
        predict(fit, newdata = nd, type = type, terms = "x",
                terms_seq = "data", ref = list(x = 1))
      }
      nd <- dat[1:4, "x", drop = FALSE]
      expected <- prediction(nd)
      for (bad in c(NA_real_, NaN, Inf, -Inf)) {
        extra <- nd
        extra$z <- extra$y <- bad
        expect_equal(prediction(extra), expected)
      }
    }
  }

  # Force formula reconstruction (the design column I(1 / z) is not supplied).
  # Checking only its finite result would incorrectly accept z = Inf.
  nd <- dat[1:4, c("x", "z")]
  nd$z[1L] <- Inf
  expect_error(predict(fit_formula, newdata = nd), "finite numeric values")
})


test_that("unrequested zero-handled terms do not validate unused values", {
  dat <- data.frame(x = seq(1, 5, length.out = 90L),
                    z = rep(c(0, 1, 2), 30L))
  dat$y <- 0.4 * dat$x + 0.2 * dat$z + 0.1 * sin(seq_len(nrow(dat)))
  fit <- mfp2(as.matrix(dat[, c("x", "z")]), dat$y, zero_vars = "z",
               df = 1, select = 1, alpha = 1, shift = 0, scale = 1,
               verbose = FALSE)
  for (type in c("terms", "contrasts")) {
    prediction <- function(nd) {
      predict(fit, newdata = nd, type = type, terms = "x",
              terms_seq = "data", ref = list(x = 1))
    }
    nd <- dat[1:4, "x", drop = FALSE]
    expected <- prediction(nd)
    # Negative values are also invalid for a zero-handled term, but irrelevant
    # when that term is not requested. The same values must fail full prediction.
    for (bad in c(-1, NA_real_, NaN, Inf, -Inf)) {
      extra <- nd
      extra$z <- bad
      expect_equal(prediction(extra), expected)
      expect_error(predict(fit, newdata = extra), "nonnegative|missing|finite")
    }
  }
})


test_that("eliminated terms and intercept-only fits ignore unused values", {
  dat <- data.frame(x = seq(1, 5, length.out = 90L),
                    z = rep(c(0, 1, 2), 30L))
  dat$y <- 0.8 * dat$x + 0.1 * sin(seq_len(nrow(dat)))
  # Formula-level select is scalar; set each variable's threshold in fp().
  fit <- mfp2(y ~ fp(x, df = 1, select = 1) + fp(z, df = 1, select = 0),
               data = dat, df = 1, select = 1, alpha = 1,
               shift = 0, scale = 1, verbose = FALSE)
  expect_identical(get_selected_variable_names(fit), "x")
  nd <- dat[1:4, "x", drop = FALSE]
  expected <- predict(fit, newdata = nd)
  for (bad in c(NA_real_, NaN, Inf, -Inf)) {
    extra <- nd
    extra$z <- bad
    expect_equal(predict(fit, newdata = extra), expected)
  }

  # A retained factor forces model-matrix reconstruction. Its eliminated
  # companion must still be excluded from the raw formula-input checks.
  dat$group <- factor(rep(c("a", "b", "c"), each = 30L))
  dat$y <- as.numeric(dat$group) + 0.1 * sin(seq_len(nrow(dat)))
  fit_factor <- mfp2(y ~ group + fp(z, df = 1, select = 0),
                      data = dat, df = 1, select = 1, alpha = 1,
                      shift = 0, scale = 1, verbose = FALSE)
  nd_factor <- dat[1:4, "group", drop = FALSE]
  expected_factor <- predict(fit_factor, newdata = nd_factor)
  for (bad in c(NA_real_, NaN, Inf, -Inf)) {
    extra <- nd_factor
    extra$z <- bad
    expect_equal(predict(fit_factor, newdata = extra), expected_factor)
  }

  fit_null <- mfp2(y ~ x + z, data = dat, df = 1, select = 0,
                    alpha = 1, verbose = FALSE)
  expect_length(get_selected_variable_names(fit_null), 0L)
  empty <- data.frame(row.names = seq_len(4L))
  expect_equal(
    predict(fit_null, newdata = data.frame(x = rep(Inf, 4L), z = NA_real_)),
    predict(fit_null, newdata = empty)
  )
  expect_error(predict(fit_null, newdata = empty[FALSE, , drop = FALSE]),
               "must contain at least one row")
})


test_that("offset inputs are checked only when the formula offset is used", {
  dat <- data.frame(x = seq(1, 5, length.out = 90L),
                    exposure = rep(c(1, 2, 3), 30L))
  dat$y <- 0.4 * dat$x + 1 / dat$exposure + 0.1 * sin(seq_len(nrow(dat)))
  fit <- mfp2(y ~ x + offset(1 / exposure), data = dat, df = 1,
               select = 1, alpha = 1, shift = 0, scale = 1, verbose = FALSE)
  nd <- dat[1:4, c("x", "exposure")]
  expected <- predict(fit, newdata = nd, newoffset = rep(0, nrow(nd)))
  for (bad in c(NA_real_, NaN, Inf, -Inf)) {
    extra <- nd
    extra$exposure <- bad
    expect_error(predict(fit, newdata = extra), "missing|finite")
    # Explicit newoffset takes precedence, so exposure is now unused.
    expect_equal(predict(fit, newdata = extra, newoffset = rep(0, nrow(nd))),
                 expected)
    for (type in c("terms", "contrasts")) {
      expect_equal(
        predict(fit, newdata = extra, type = type, terms_seq = "data"),
        predict(fit, newdata = nd["x"], type = type, terms_seq = "data")
      )
    }
  }
})


test_that("Cox validation uses only the response and strata inputs requested", {
  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[stats::complete.cases(dat[, c("time", "status", "age", "sex")]), ]
  fit <- mfp2(survival::Surv(1 / time, status) ~ age + strata(1 / sex),
               data = dat, family = "cox", df = 1, select = 1, alpha = 1,
               shift = 0, scale = 1, verbose = FALSE)
  nd <- dat[1:4, c("time", "status", "age", "sex")]
  nd$unused <- Inf
  expect_equal(predict(fit, newdata = nd, type = "survival"),
               predict(fit, newdata = nd[, names(nd) != "unused"], type = "survival"))
  expected_lp <- predict(fit, newdata = nd, type = "lp")
  expected_term <- predict(fit, newdata = nd["age"], type = "terms",
                            terms = "age", terms_seq = "data")
  for (bad in c(NA_real_, NaN, Inf, -Inf)) {
    extra <- nd
    extra$time <- extra$status <- bad
    expect_equal(predict(fit, newdata = extra, type = "lp"), expected_lp)
    expect_error(predict(fit, newdata = extra, type = "survival"), "missing|finite")
    extra$sex <- bad
    expect_error(predict(fit, newdata = extra, type = "lp"), "missing|finite")
    expect_equal(predict(fit, newdata = extra, type = "terms", terms = "age",
                         terms_seq = "data"), expected_term)
  }
})


test_that("full-model predictions preserve supplied row names through transformation", {
  dat <- data.frame(x = seq(1, 5, length.out = 90L),
                    z = rep(c(0, 1, 2), 30L))
  dat$y <- 0.8 * dat$x + 0.2 * dat$z + 0.1 * sin(seq_len(nrow(dat)))
  formula_fit <- mfp2(y ~ x + z, data = dat, df = 1, select = 1,
                       alpha = 1, verbose = FALSE)
  matrix_fit <- mfp2(as.matrix(dat[, c("x", "z")]), dat$y, df = 1,
                      select = 1, alpha = 1, verbose = FALSE)
  for (fit in list(formula_fit, matrix_fit)) {
    nd <- dat[c(40L, 3L, 70L), c("x", "z")]
    row.names(nd) <- c("subject-C", "subject-A", "subject-F")
    for (input in list(nd, as.matrix(nd))) {
      for (type in c("link", "response")) {
        point <- predict(fit, newdata = input, type = type)
        result <- predict(fit, newdata = input, type = type, se.fit = TRUE)
        expect_identical(names(point), row.names(nd))
        expect_identical(names(result$fit), row.names(nd))
        expect_identical(names(result$se.fit), row.names(nd))
        expect_equal(result$fit, point)
      }
    }
  }

  # Neither formula reconstruction nor an empty transformed design may reset
  # row names when selection removes every predictor.
  null_fit <- mfp2(y ~ x + z, data = dat, df = 1, select = 0, verbose = FALSE)
  nd <- data.frame(row.names = c("subject-C", "subject-A", "subject-F"))
  reconstructed <- mfp2:::reconstruct_formula_newdata(null_fit, nd)
  prepared <- mfp2:::prepare_newdata_for_predict(null_fit, reconstructed)
  expect_identical(row.names(reconstructed), row.names(nd))
  expect_identical(row.names(prepared), row.names(nd))
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
    x = data.matrix(dat[, c("age", "sex"), drop = FALSE]),
    y = survival::Surv(dat$time, dat$status),
    family = "cox",
    strata = dat$stratum,
    df = 1,
    select = 1,
    alpha = 1,
    verbose = FALSE
  )
  nd_stratified <- data.matrix(dat[1:8, c("age", "sex"), drop = FALSE])
  bad_strata <- as.numeric(dat$stratum[1:8])
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
