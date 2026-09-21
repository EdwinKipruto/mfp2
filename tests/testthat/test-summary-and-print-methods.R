# Area 6: summary() and print() methods for each family.

capture_ok <- function(expr) {
  out <- tryCatch(utils::capture.output(expr), error = function(e) e)
  expect_false(inherits(out, "error"))
  paste(out, collapse = "\n")
}

test_that("summary/print work and identify the family (survreg)", {
  skip_on_cran()
  data(gbsg, package = "mfp2")
  fit <- mfp2(Surv(rectime, censrec) ~ fp(age, df = 1, select = 1),
              data = gbsg, family = survreg_family(dist = "weibull"), verbose = FALSE)
  txt_p <- capture_ok(print(fit))
  txt_s <- capture_ok(print(summary(fit)))
  expect_match(paste(txt_p, txt_s), "survreg|Weibull|weibull", ignore.case = TRUE)
})

test_that("summary/print work and identify the family (finegray)", {
  skip_on_cran()
  dfg <- make_finegray_data(n = 500, seed = 61)
  fit <- mfp2(Surv(obs, ev) ~ fp(z1, df = 1, select = 1),
              data = dfg, family = finegray_family(etype = "relapse"), verbose = FALSE)
  txt_p <- capture_ok(print(fit))
  txt_s <- capture_ok(print(summary(fit)))
  expect_match(paste(txt_p, txt_s), "fine|gray|grey|relapse|subdist",
               ignore.case = TRUE)
})

test_that("summary/print work and identify the family (multinomial)", {
  skip_on_cran()
  dm <- make_multinomial_data(n = 300, seed = 62)
  fit <- mfp2(y ~ fp(x1, df = 1, select = 1), data = dm,
              family = multinomial_family(), control = list(maxit = 300),
              verbose = FALSE)
  txt_p <- capture_ok(print(fit))
  txt_s <- capture_ok(print(summary(fit)))
  expect_match(paste(txt_p, txt_s), "multinom", ignore.case = TRUE)
})

test_that("mfpi summary/print include interaction test statistics", {
  skip_on_cran()
  dfg <- make_mfpi_finegray()
  fit <- mfpi(Surv(obs, ev) ~ grp + fp(z1), data = dfg, group_var = "grp",
              interaction_vars = "z1", interaction_forms = c(z1 = "fp1"),
              family = finegray_family(etype = "relapse"),
              p_interact = 0.1, verbose = FALSE)
  txt_p <- capture_ok(print(fit))
  txt_s <- capture_ok(print(summary(fit)))
  # interaction output should mention a p-value / interaction test
  expect_match(paste(txt_p, txt_s), "p.?value|interact", ignore.case = TRUE)
})
