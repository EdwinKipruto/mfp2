test_that("Fine-Gray summary refits preserve the fitted baseline risk sets", {
  skip_on_cran()
  dat <- make_finegray_data(n = 600, seed = 918, strong_strata = TRUE)
  dat$off <- 0.05 * (dat$z1 - 5)

  for (mode in c("both", "baseline", "censoring")) {
    fit <- mfp2(
      Surv(obs, ev) ~ fp(z1, df = 4, select = 1, alpha = 1) +
        strata(sx) + offset(off),
      data = dat,
      family = finegray_family(etype = "relapse", strata_action = mode),
      verbose = FALSE
    )
    context <- mfp2_summary_refit_context(fit)
    expect_true(context$valid, info = mode)
    expect_identical(is.null(context$strata), mode == "censoring", info = mode)
    if (!is.null(context$strata)) {
      expect_equal(length(context$strata), NROW(context$y), info = mode)
      expect_identical(as.integer(context$strata), as.integer(fit$strata),
                       info = mode)
    }

    classified <- mfp2_summary_classify_terms(fit)
    expect_true(classified$is_nonlinear[match("z1", classified$variable_names)],
                info = mode)
    keep <- setdiff(colnames(fit$x), classified$cols_by_var[["z1"]])
    reduced_x <- fit$x[, keep, drop = FALSE]
    reduced <- survival::agreg.fit(
      x = reduced_x, y = context$y,
      strata = if (is.null(context$strata)) NULL else as.integer(context$strata),
      offset = context$offset, init = NULL, control = context$control,
      weights = context$weights, method = context$method,
      rownames = NULL, resid = FALSE, nocenter = context$nocenter
    )
    expected_lr <- max(0, 2 * (context$full_logl - tail(reduced$loglik, 1L)))
    term_test <- mfp2_summary_lrt_drop_variable(fit, classified, "z1", context)
    expect_equal(term_test$lr, expected_lr, tolerance = 1e-6, info = mode)
    summary_row <- summary(fit)$nonlinear_terms
    expect_equal(summary_row$lr_chisq[summary_row$variable == "z1"],
                 expected_lr, tolerance = 1e-6, info = mode)

    # An independent, row-wise native survival calculation is the reference
    # for both a mixed-stratum batch and a pooled-baseline batch.
    by_stratum <- split(seq_len(nrow(dat)), dat$sx)
    rows <- c(by_stratum[[1L]][1:2], by_stratum[[2L]][1:2])[c(1, 3, 4, 2)]
    newdata <- fit$mfp2_finegray_training_newdata[rows, , drop = FALSE]
    row.names(newdata) <- c("first_f", "first_m", "second_m", "second_f")
    tmid <- as.numeric(stats::quantile(dat$obs, c(0.3, 0.7)))
    times <- c(0, tmid, tmid[1L], max(dat$obs) + 1)
    time_grid <- sort(unique(times))
    base_fit <- fit
    class(base_fit) <- setdiff(class(base_fit), "mfp2")
    reference <- lapply(seq_len(nrow(newdata)), function(i) {
      curve <- survival::survfit(base_fit,
                                 newdata = newdata[i, , drop = FALSE],
                                 se.fit = TRUE)
      summary(curve, times = time_grid, extend = TRUE)
    })
    expected_fit <- t(vapply(reference, function(x) 1 - as.numeric(x$surv),
                             numeric(length(time_grid))))
    expected_se <- t(vapply(reference, function(x) as.numeric(x$std.err),
                            numeric(length(time_grid))))

    batch <- mfp2_predict_finegray_cif(fit, newdata, times, se.fit = TRUE)
    expect_identical(dimnames(batch$fit),
                     list(row.names(newdata),
                          paste0("time=", format(time_grid, trim = TRUE,
                                                 scientific = FALSE))),
                     info = mode)
    expect_equal(unname(batch$fit), unname(expected_fit),
                 tolerance = 1e-8, info = mode)
    expect_equal(unname(batch$se.fit), unname(expected_se),
                 tolerance = 1e-8, info = mode)
    expect_equal(mfp2_predict_finegray_cif(fit, newdata, times), batch$fit,
                 info = mode)

    one_time <- mfp2_predict_finegray_cif(fit, newdata[2L, , drop = FALSE],
                                          tmid[1L], se.fit = TRUE)
    expect_identical(names(one_time$fit), row.names(newdata)[2L], info = mode)
    expect_equal(unname(one_time$fit), expected_fit[2L, 2L],
                 tolerance = 1e-8, info = mode)
    expect_equal(unname(one_time$se.fit), expected_se[2L, 2L],
                 tolerance = 1e-8, info = mode)

    if (mode == "both") {
      many_rows <- fit$mfp2_finegray_training_newdata[1:257, , drop = FALSE]
      chunked <- mfp2_predict_finegray_cif(fit, many_rows, tmid[1L])
      expect_length(chunked, 257L)
      for (i in c(1L, 256L, 257L)) {
        native <- summary(survival::survfit(
          base_fit, newdata = many_rows[i, , drop = FALSE], se.fit = FALSE
        ), times = tmid[1L], extend = TRUE)
        expect_equal(unname(chunked[i]), 1 - as.numeric(native$surv),
                     tolerance = 1e-8, info = paste(mode, i))
      }
    }
  }
})
