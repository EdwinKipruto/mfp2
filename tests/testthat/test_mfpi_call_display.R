test_that("MFPI calls use the public name for both interfaces", {
  dat <- make_mfpi_factor_data()

  calls <- list(
    default = quote(mfpi(
      x = dat[, c("trt", "x")], y = dat$y,
      group_var = "trt", interaction_vars = "x",
      interaction_forms = c(x = "linear"), flex = "flex1",
      df = 1, select = 1, verbose = FALSE
    )),
    formula = quote(mfpi(
      formula = y ~ trt + x, data = dat,
      group_var = "trt", interaction_vars = "x",
      interaction_forms = c(x = "linear"), flex = "flex1",
      df = 1, select = 1, verbose = FALSE
    ))
  )

  for (interface in names(calls)) {
    expected_call <- calls[[interface]]
    for (entry in c("mfpi", paste0("mfpi.", interface))) {
      input_call <- expected_call
      input_call[[1L]] <- as.name(entry)
      if (interface == "formula" && entry == "mfpi") {
        # The generic dispatches on x; the formula method records it as formula.
        names(input_call)[2L] <- ""
      }
      fit <- eval(input_call)

      expect_identical(fit$call, expected_call)
      expect_identical(summary(fit)$call, expected_call)
      printed <- paste(capture.output(print(fit)), collapse = "\n")
      expect_match(printed, "mfpi(", fixed = TRUE)
      expect_false(grepl("mfpi.default(", printed, fixed = TRUE))
      expect_false(grepl("mfpi.formula(", printed, fixed = TRUE))
    }
  }
})
