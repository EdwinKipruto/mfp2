# Scoping of the `select` argument in the formula interface.
#
# Design intent (now documented explicitly in ?mfp2 and ?fp):
#   * The top-level `select` given to mfp2() sets the selection level for plain
#     (non-fp()) terms only.
#   * Terms wrapped in fp() carry their OWN `select`, which defaults to 0.05 and
#     is NOT inherited from the top-level value.
#
# So a top-level `select = 1` does not force an fp() term in; only a per-term
# fp(x, select = 1) (or listing the variable in `keep`) does. These tests pin
# that documented contract for the new survival families.

test_that("top-level select = 1 does NOT force an fp() term (documented scoping)", {
  skip_on_cran()
  data(gbsg, package = "mfp2")
  # age is a weak predictor of recurrence time here; at the fp() default
  # select = 0.05 it is dropped, and top-level select = 1 must not change that.
  fit_global <- mfp2(Surv(rectime, censrec) ~ fp(age, df = 1) +
                       fp(nodes, df = 1) + hormon,
                     data = gbsg, family = survreg_family(dist = "weibull"),
                     select = 1, verbose = FALSE)
  expect_false("age.1" %in% names(coef(fit_global)))
})

test_that("per-term fp(select = 1) does force the predictor in", {
  skip_on_cran()
  data(gbsg, package = "mfp2")
  fit_perterm <- mfp2(Surv(rectime, censrec) ~ fp(age, df = 1, select = 1) +
                        fp(nodes, df = 1, select = 1) + hormon,
                      data = gbsg, family = survreg_family(dist = "weibull"),
                      verbose = FALSE)
  expect_true("age.1" %in% names(coef(fit_perterm)))
})

test_that("keep also forces an fp() term in, matching the documented equivalence", {
  skip_on_cran()
  data(gbsg, package = "mfp2")
  fit_keep <- mfp2(Surv(rectime, censrec) ~ fp(age, df = 1) +
                     fp(nodes, df = 1) + hormon,
                   data = gbsg, family = survreg_family(dist = "weibull"),
                   keep = "age", verbose = FALSE)
  expect_true("age.1" %in% names(coef(fit_keep)))
})
