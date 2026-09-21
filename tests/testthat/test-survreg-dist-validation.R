# Item 3: survreg_family() now validates a character `dist` against
# survival::survreg.distributions at construction time (fail fast), instead of
# deferring the error to fit time. A list-valued `dist` (a custom distribution
# object) is still passed through unchecked.

test_that("an unrecognized distribution name errors at construction", {
  skip_on_cran()
  expect_error(survreg_family(dist = "invalid_dist"),
               "not a recognized")
  expect_error(survreg_family(dist = "gauss"),      # near-miss typo
               "not a recognized")
})

test_that("the fail-fast error names the valid built-in choices", {
  skip_on_cran()
  err <- tryCatch(survreg_family(dist = "nope"), error = function(e) conditionMessage(e))
  # A couple of the canonical built-ins should be listed to guide the user.
  expect_true(grepl("weibull", err))
  expect_true(grepl("lognormal", err))
})

test_that("all built-in survreg distributions still construct", {
  skip_on_cran()
  for (d in names(survival::survreg.distributions)) {
    expect_s3_class(survreg_family(dist = d), "mfp2_survreg_family")
  }
})

test_that("a custom list-valued distribution is passed through unchecked", {
  skip_on_cran()
  custom <- survival::survreg.distributions$weibull
  expect_s3_class(survreg_family(dist = custom), "mfp2_survreg_family")
})

test_that("construction-time validation does not disturb a valid end-to-end fit", {
  skip_on_cran()
  data(gbsg, package = "mfp2")
  fit <- mfp2(Surv(rectime, censrec) ~ fp(age, df = 1, select = 1) +
                fp(nodes, df = 1, select = 1),
              data = gbsg, family = survreg_family(dist = "loglogistic"),
              verbose = FALSE)
  expect_equal(fit$family_string, "survreg")
  expect_true(inherits(fit, "survreg"))
})
