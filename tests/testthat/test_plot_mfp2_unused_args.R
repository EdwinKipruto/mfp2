# plot.mfp2() has no newdata argument. Arguments caught by `...` used to be
# ignored silently, so plot(fit, newdata = nd) looked as if it plotted nd.

test_that("plot() rejects newdata and warns about other unused arguments", {
  skip_if_not_installed("ggplot2")
  data("prostate", package = "mfp2")
  fit <- mfp2(lpsa ~ fp(cavol, df = 2, select = 1) + age, data = prostate,
              select = 1, verbose = FALSE)

  expect_error(plot(fit, terms = "cavol", newdata = prostate[1:5, ]),
               "has no `newdata` argument")
  expect_warning(plot(fit, terms = "cavol", main = "x"),
                 "Ignoring argument\\(s\\) not used by `plot\\(\\)`.*main")
  expect_no_warning(p <- plot(fit, terms = "cavol"))
  expect_s3_class(p$cavol, "ggplot")
})

test_that("plotted term values equal predict(type = 'terms')", {
  skip_if_not_installed("ggplot2")
  data("prostate", package = "mfp2")
  fit <- mfp2(lpsa ~ fp(cavol, df = 4, select = 1, force_max_fp = TRUE) + age,
              data = prostate, select = 1, verbose = FALSE)
  p <- plot(fit, terms = "cavol", partial_only = TRUE)$cavol
  built <- ggplot2::ggplot_build(p)
  geoms <- vapply(p$layers, function(l) class(l$geom)[1L], character(1L))
  line <- built$data[[which(geoms == "GeomLine")]]
  pred <- predict(fit, type = "terms", terms = "cavol", terms_seq = "data")$cavol
  pred <- pred[order(pred$variable), ]
  expect_equal(line$y, pred$value)
  # Intercept is included by default (add_intercept = TRUE).
  no_int <- predict(fit, type = "terms", terms = "cavol", terms_seq = "data",
                    add_intercept = FALSE)$cavol
  no_int <- no_int[order(no_int$variable), ]
  expect_equal(line$y - no_int$value,
               rep(unname(coef(fit)[["(Intercept)"]]), nrow(no_int)))
})
