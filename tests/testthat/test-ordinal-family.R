# Ordinal (proportional-odds) family via rms::orm. Fitting equivalence is
# checked against a native rms::orm() fit; the candidate path uses orm.fit()
# (integer-coded response) for speed, and the retained model is a native orm
# object so coef/vcov/logLik/AIC/predict work.

skip_if_no_rms <- function() skip_if_not_installed("rms")

make_ordinal_data <- function(n = 500, seed = 21) {
  set.seed(seed)
  age <- runif(n, 20, 80)
  sex <- rbinom(n, 1, 0.5)
  eta <- 0.05 * age + 0.6 * sex
  y <- as.integer(cut(eta + rlogis(n), 4))
  data.frame(y = ordered(y), yint = y, age = age, sex = sex)
}

test_that("ordinal_family constructor validates the link and stores it", {
  skip_on_cran(); skip_if_no_rms()
  fam <- ordinal_family()
  expect_s3_class(fam, "mfp2_ordinal_family")
  expect_equal(fam$link, "logistic")
  expect_equal(ordinal_family(link = "probit")$link, "probit")
  expect_error(ordinal_family(link = "foo"))
})

test_that("forced-linear mfp2 ordinal reproduces native rms::orm", {
  skip_on_cran(); skip_if_no_rms()
  dm <- make_ordinal_data()
  fit <- mfp2(y ~ fp(age, df = 1, select = 1, scale = 1, shift = 0) +
                fp(sex, df = 1, select = 1),
              data = dm, family = ordinal_family(), center = FALSE, verbose = FALSE)
  ref <- rms::orm(yint ~ age + sex, data = dm, family = "logistic")

  expect_equal(fit$family_string, "ordinal")
  expect_true(inherits(fit, "orm"))
  expect_s3_class(fit$mfp2_family, "mfp2_ordinal_family")
  expect_equal(fit$mfp2_family$prepared$y, as.integer(dm$y))
  cf <- coef(fit); names(cf) <- sub("\\.1$", "", names(cf))
  expect_equal(unname(cf[["age"]]), unname(coef(ref)[["age"]]), tolerance = 1e-5)
  expect_equal(unname(cf[["sex"]]), unname(coef(ref)[["sex"]]), tolerance = 1e-5)
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)), tolerance = 1e-4)
  expect_equal(
    unname(fit$mfp2_ordinal_intercepts),
    unname(fit$coefficients[seq_len(fit$non.slopes)]),
    tolerance = 0
  )
  expect_true(all(c("maxit", "eps", "tol", "trace") %in%
                    names(as.list(fit$call))))
  expect_true(is.matrix(vcov(fit)))
  expect_true(is.finite(AIC(fit)))
})

test_that("FP selection works out of the box with default controls", {
  skip_on_cran(); skip_if_no_rms()
  dm <- make_ordinal_data()
  expect_error(
    mfp2(y ~ fp(age, df = 4) + fp(sex, df = 1), data = dm,
         family = ordinal_family(), criterion = "aic", verbose = FALSE),
    NA
  )
})

test_that("all five links fit", {
  skip_on_cran(); skip_if_no_rms()
  dm <- make_ordinal_data()
  for (lk in c("logistic", "probit", "loglog", "cloglog", "cauchit")) {
    expect_error(
      mfp2(y ~ fp(age, df = 1, select = 1) + fp(sex, df = 1, select = 1),
           data = dm, family = ordinal_family(link = lk), verbose = FALSE),
      NA, info = lk
    )
  }
})

test_that("ordinal control uses explicit orm candidate settings", {
  skip_on_cran(); skip_if_no_rms()
  ctl <- mfp2:::normalize_ordinal_control()
  expect_equal(ctl$maxit, 30L)
  expect_setequal(names(ctl), c("maxit", "eps", "tol", "trace"))
  expect_null(ctl$epsilon)
  # glm-style epsilon maps onto eps
  expect_equal(mfp2:::normalize_ordinal_control(list(epsilon = 1e-6))$eps, 1e-6)
  expect_error(mfp2:::normalize_ordinal_control(list(foo = 1)), "Unknown ordinal control")
  expect_error(mfp2:::normalize_ordinal_control(list(maxit = 2.5)),
               "positive integer")
  expect_error(mfp2:::normalize_ordinal_control(list(maxit = Inf)),
               "positive integer")
})

test_that("response contract matches rms with a guard message for unordered input", {
  skip_on_cran(); skip_if_no_rms()
  set.seed(3); n <- 400; age <- runif(n, 20, 80)
  # ordered factor / numeric: no message
  yo <- ordered(as.integer(cut(0.05 * age + rlogis(n), 4)))
  expect_silent(
    mfp2(y ~ fp(age, df = 1, select = 1), data = data.frame(y = yo, age = age),
         family = ordinal_family(), verbose = FALSE)
  )
  # character response: informational order message
  yc <- sample(c("low", "med", "high"), n, TRUE)
  msg <- character(0)
  withCallingHandlers(
    mfp2(y ~ fp(age, df = 1, select = 1), data = data.frame(y = yc, age = age),
         family = ordinal_family(), verbose = FALSE),
    message = function(m) { msg <<- c(msg, conditionMessage(m)); invokeRestart("muffleMessage") }
  )
  expect_true(any(grepl("order inferred", msg)))
})

test_that("two-level response and case weights are rejected", {
  skip_on_cran(); skip_if_no_rms()
  set.seed(5); n <- 300; age <- runif(n, 20, 80)
  expect_error(
    mfp2(y ~ fp(age, df = 1), data = data.frame(y = ordered(rbinom(n, 1, .5)), age = age),
         family = ordinal_family(), verbose = FALSE)
  )
  dm <- make_ordinal_data()
  expect_error(
    mfp2(y ~ fp(age, df = 1, select = 1), data = dm, family = ordinal_family(),
         weights = runif(nrow(dm), 1, 2), verbose = FALSE),
    "does not support case weights"
  )
  expect_error(
    mfp2(y ~ fp(age, df = 1, select = 1), data = dm, family = ordinal_family(),
         weights = rep(2, nrow(dm)), verbose = FALSE),
    "does not support case weights"
  )
})

test_that("ordinal preparation drops unused response levels", {
  skip_on_cran(); skip_if_no_rms()
  y <- ordered(c("low", "middle", "high", "middle"),
               levels = c("unused", "low", "middle", "high"))
  prepared <- mfp2:::prepare_ordinal_family(
    ordinal_family(), y, weights = rep(1, length(y))
  )$prepared
  expect_identical(prepared$levels, c("low", "middle", "high"))
  expect_identical(sort(unique(prepared$y)), 1:3)
  expect_equal(prepared$n_intercepts, 2L)
  expect_false(prepared$has_offset)
  expect_null(prepared$offset)

  expect_error(
    mfp2:::prepare_ordinal_family(
      ordinal_family(), y, weights = rep(2, length(y))
    ),
    "does not support case weights"
  )
})

test_that("ordinal preparation caches a supplied offset once", {
  skip_on_cran(); skip_if_no_rms()
  dm <- make_ordinal_data(n = 120, seed = 24)
  off <- seq(-0.2, 0.2, length.out = nrow(dm))
  prepared <- mfp2:::prepare_ordinal_family(
    ordinal_family(), dm$y, weights = rep.int(1, nrow(dm)),
    offset = off, has_offset = TRUE
  )$prepared

  expect_true(prepared$has_offset)
  expect_identical(prepared$offset, off)
})

test_that("intercept-only ordinal fitting retains a native orm model", {
  skip_on_cran(); skip_if_no_rms()
  dm <- make_ordinal_data(n = 500, seed = 25)
  ycodes <- as.integer(dm$y)
  family <- mfp2:::prepare_ordinal_family(
    ordinal_family(), dm$y, weights = rep.int(1, nrow(dm)),
    has_offset = FALSE
  )

  fitted <- mfp2:::fit_ordinal(
    x = matrix(numeric(0L), nrow = nrow(dm), ncol = 0L),
    family = family,
    control = mfp2:::normalize_ordinal_control(),
    fast = FALSE,
    keep_fit = TRUE,
    x_has_intercept = FALSE
  )
  reference <- rms::orm(ycodes ~ 1, data = data.frame(ycodes),
                        x = TRUE, y = TRUE)

  expect_s3_class(fitted$fit, "orm")
  expect_length(fitted$coefficients, 0L)
  expect_equal(fitted$fit$coefficients, reference$coefficients,
               tolerance = 1e-8)
  expect_true(all(is.finite(diag(stats::vcov(fitted$fit)))))
  expect_equal(as.numeric(stats::logLik(fitted$fit)),
               as.numeric(stats::logLik(reference)), tolerance = 1e-8)
  expect_equal(fitted$logl, as.numeric(stats::logLik(reference)),
               tolerance = 1e-8)
  expect_true(is.finite(stats::AIC(fitted$fit)))
})

test_that("intercept-only ordinal fitting honors a varying offset", {
  skip_on_cran(); skip_if_no_rms()
  dm <- make_ordinal_data(n = 500, seed = 26)
  ycodes <- as.integer(dm$y)
  off <- seq(-0.4, 0.4, length.out = nrow(dm))
  family <- mfp2:::prepare_ordinal_family(
    ordinal_family(), dm$y, weights = rep.int(1, nrow(dm)),
    offset = off, has_offset = TRUE
  )

  fitted <- mfp2:::fit_ordinal(
    x = matrix(numeric(0L), nrow = nrow(dm), ncol = 0L),
    family = family,
    control = mfp2:::normalize_ordinal_control(),
    fast = FALSE,
    keep_fit = TRUE,
    x_has_intercept = FALSE
  )
  reference_data <- data.frame(ycodes = ycodes, off = off)
  reference <- rms::orm(
    ycodes ~ offset(off), data = reference_data,
    family = "logistic", x = TRUE, y = TRUE
  )

  expect_equal(fitted$fit$coefficients, reference$coefficients,
               tolerance = 1e-8)
  expect_equal(as.numeric(stats::logLik(fitted$fit)),
               as.numeric(stats::logLik(reference)), tolerance = 1e-8)
  expect_equal(fitted$logl, as.numeric(stats::logLik(reference)),
               tolerance = 1e-8)
  expect_true(all(is.finite(diag(stats::vcov(fitted$fit)))))
})

test_that("intercept-only ordinal prediction returns every requested row", {
  skip_on_cran(); skip_if_no_rms()
  intercepts <- c("y>=2" = 0.7, "y>=3" = -0.4)
  object <- list(
    coefficients = intercepts,
    mfp2_ordinal_intercepts = intercepts,
    mfp2_ordinal_link = "logistic",
    mfp2_ordinal_levels = c("low", "middle", "high"),
    x = NULL,
    nobs = 6L
  )
  probabilities <- mfp2:::mfp2_predict_ordinal(
    object, type = "response", se.fit = FALSE
  )
  expect_equal(dim(probabilities), c(6L, 3L))
  expect_equal(rowSums(probabilities), rep(1, 6L), tolerance = 1e-12)
})

test_that("predict returns correct link, probabilities, and mean", {
  skip_on_cran(); skip_if_no_rms()
  dm <- make_ordinal_data()
  fit <- mfp2(y ~ fp(age, df = 1, select = 1, scale = 1, shift = 0) +
                fp(sex, df = 1, select = 1),
              data = dm, family = ordinal_family(), center = FALSE, verbose = FALSE)
  ref <- rms::orm(yint ~ age + sex, data = dm, family = "logistic", x = TRUE, y = TRUE)

  probs <- predict(fit, type = "response")
  expect_equal(dim(probs), c(nrow(dm), 4L))
  expect_equal(max(abs(rowSums(probs) - 1)), 0, tolerance = 1e-10)
  expect_equal(max(abs(probs - predict(ref, type = "fitted.ind"))), 0, tolerance = 1e-6)

  # covariate linear predictor differs from orm's centred lp only by a constant
  lp <- predict(fit, type = "link")
  expect_equal(sd(lp - predict(ref, type = "lp")), 0, tolerance = 1e-8)

  # newdata prediction agrees with the training-data subset
  expect_equal(
    max(abs(predict(fit, newdata = dm[1:10, ], type = "response") - probs[1:10, ])),
    0, tolerance = 1e-10
  )
  expect_true(all(is.finite(predict(fit, type = "mean"))))
  # lp alias
  expect_equal(unname(predict(fit, type = "lp")), unname(lp), tolerance = 1e-12)
})

test_that("mfpi() supports the ordinal family end to end", {
  skip_on_cran(); skip_if_no_rms()
  # mfpi() now supports proportional-odds ordinal models via rms::orm; the
  # previous guard ("does not yet support") has been removed. Use a group
  # with a genuine effect and modification so the interaction model is
  # well-conditioned (a random group produced a singular Hessian).
  set.seed(202)
  n <- 400L
  g <- factor(rep(c("c", "t"), each = n / 2L))
  age <- runif(n, 20, 80)
  eta <- 0.04 * age + 0.5 * (g == "t") + 0.03 * age * (g == "t")
  y <- ordered(as.integer(cut(eta + rlogis(n), 4)))
  dm <- data.frame(y = y, age = age, g = g)

  fit <- mfpi(
    y ~ g + fp(age, df = 1) + 1,
    data = dm,
    group_var = "g",
    interaction_vars = "age",
    interaction_forms = c(age = "linear"),
    flex = "flex1",
    family = ordinal_family(),
    cycles = 1, df = 1, select = 1, alpha = 1, p_interact = 1,
    shift = 0, scale = 1, center = FALSE, verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_equal(fit$family_string, "ordinal")
  expect_true(is.character(fit$ordinal_levels))
  expect_equal(fit$n_intercepts, length(fit$ordinal_levels) - 1L)
  sm <- summary(fit)
  expect_identical(sm$family_string, "ordinal")
  expect_identical(sm$ordinal_levels, fit$ordinal_levels)
  expect_identical(sm$ordinal_link, fit$ordinal_link)
  expect_identical(sm$n_intercepts, fit$n_intercepts)

  # Ordinary response prediction yields a valid class-probability matrix.
  nd <- dm[1:10, c("g", "age"), drop = FALSE]
  pr <- predict(fit, newdata = nd, terms = "age", model = "all",
                type = "response", se.fit = FALSE)
  fit_cols <- grep("^fit(\\.|$)", names(pr$predictions), value = TRUE)
  probs <- as.matrix(pr$predictions[, fit_cols, drop = FALSE])
  expect_equal(ncol(probs), length(fit$ordinal_levels))
  expect_equal(rowSums(probs), rep(1, nrow(nd)), tolerance = 1e-8)

  # Training-data prediction must retain all rows rather than returning the
  # single placeholder row formerly used when no design was reconstructed.
  training <- predict(fit, terms = "age", model = "all",
                      type = "response", se.fit = FALSE)$predictions
  training_cols <- grep("^fit(\\.|$)", names(training), value = TRUE)
  expect_equal(nrow(training), nrow(dm))
  expect_equal(rowSums(training[, training_cols, drop = FALSE]),
               rep(1, nrow(dm)), tolerance = 1e-8)
})
