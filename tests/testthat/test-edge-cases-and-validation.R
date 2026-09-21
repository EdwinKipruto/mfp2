# Area 5: edge cases and input validation.

test_that("multinomial_family with a 2-class response errors", {
  skip_on_cran()
  d2 <- data.frame(y = factor(sample(c("A", "B"), 200, TRUE)), x = runif(200, 1, 5))
  expect_error(
    mfp2(y ~ fp(x, df = 1, select = 1), data = d2,
         family = multinomial_family(), verbose = FALSE),
    regexp = "three classes|at least three", ignore.case = TRUE
  )
})

test_that("finegray_family with a non-multi-state Surv errors", {
  skip_on_cran()
  data(gbsg, package = "mfp2")
  expect_error(
    mfp2(Surv(rectime, censrec) ~ fp(age, df = 1, select = 1), data = gbsg,
         family = finegray_family(etype = "relapse"), verbose = FALSE),
    regexp = "multi-state|multistate", ignore.case = TRUE
  )
})

test_that("survreg with an invalid distribution errors informatively", {
  skip_on_cran()
  data(gbsg, package = "mfp2")
  # NOTE: survreg_family(dist = "invalid_dist") does not error at construction;
  # the informative error is raised when the model is fitted.
  expect_error(
    mfp2(Surv(rectime, censrec) ~ fp(age, df = 1, select = 1), data = gbsg,
         family = survreg_family(dist = "invalid_dist"), verbose = FALSE),
    regexp = "distribution", ignore.case = TRUE
  )
})

test_that("survreg_family(scale = -1) errors at construction", {
  expect_error(survreg_family(scale = -1),
               regexp = "non-negative|finite", ignore.case = TRUE)
})

test_that("finegray_family(strata_action = 'invalid') errors", {
  expect_error(finegray_family(strata_action = "invalid"),
               regexp = "both|censoring|baseline", ignore.case = TRUE)
})

test_that("mfpi multinomial with a single-level group errors or warns", {
  skip_on_cran()
  dm <- make_mfpi_multinomial(n = 400, seed = 91)
  dm$g <- factor(rep("0", nrow(dm)))          # only one group level
  # Acceptable outcomes: an error, or a warning about the group having < 2 levels.
  res <- tryCatch(
    mfpi(y ~ g + fp(a), data = dm, group_var = "g",
         interaction_vars = "a", interaction_forms = c(a = "fp1"),
         family = multinomial_family(), control = list(maxit = 300),
         verbose = FALSE),
    error   = function(e) "error",
    warning = function(w) "warning"
  )
  expect_true(identical(res, "error") || identical(res, "warning"))
})
