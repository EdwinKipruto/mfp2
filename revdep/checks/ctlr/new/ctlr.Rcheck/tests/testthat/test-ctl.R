test_that("ctl works with basic inputs", {
  # Test that ctl runs without error
  expect_no_error(
    ctl(
      ctl_dataset1,
      idvar = "id",
      ynew = "y1",
      yref = "y2",
      intercept = 5,
      slope = 0
    )
  )
})

test_that("ctl returns correct structure when results = TRUE", {
  result <- ctl(
    ctl_dataset1,
    idvar = "id",
    ynew = "y1",
    yref = "y2",
    intercept = 5,
    slope = 0
  )

  expect_type(result, "list")
  expect_named(result, c("data", "overall_agreement"))
  expect_s3_class(result$data, "data.frame")
  expect_type(result$overall_agreement, "list")
})

test_that("ctl works with cpaplot = TRUE", {
  testthat::skip_if_not_installed("withr")
  withr::with_tempdir({
    expect_no_error(
      ctl(
        ctl_dataset1,
        idvar = "id",
        ynew = "y1",
        yref = "y2",
        intercept = 0,
        slope = 0.15,
        cpaplot = TRUE,
        nbsimul = 100,
        plots_to_file = TRUE,
        outputs_path = "."
      )
    )
    expect_true(dir.exists("ctl_outputs"))
    expect_true(file.exists("ctl_outputs/cpa_plot.pdf"))
  })
})

test_that("ctl works with tlplot = TRUE", {
  testthat::skip_if_not_installed("withr")
  withr::with_tempdir({
    expect_no_error(
      ctl(
        ctl_dataset1,
        idvar = "id",
        ynew = "y1",
        yref = "y2",
        intercept = 5,
        slope = 0,
        tlplot = TRUE,
        plots_to_file = TRUE,
        outputs_path = "."
      )
    )

    expect_true(dir.exists("ctl_outputs"))
    expect_true(file.exists("ctl_outputs/tl_plot.pdf"))
  })
})

test_that("ctl validates overall agreement calculation", {
  testthat::skip_if_not_installed("withr")
  withr::with_tempdir({
    result <- ctl(
      ctl_dataset1,
      idvar = "id",
      ynew = "y1",
      yref = "y2",
      intercept = 5,
      slope = 0,
      results_to_file = TRUE,
      outputs_path = "."
    )

    expect_true(dir.exists("ctl_outputs"))
    expect_true(file.exists("ctl_outputs/ctl_results.rda"))

    oa <- result$overall_agreement
    expect_true(oa$overall_agreement >= 0 && oa$overall_agreement <= 1)
    expect_true(oa$overall_agreement_lo >= 0 && oa$overall_agreement_lo <= 1)
    expect_true(oa$overall_agreement_up >= 0 && oa$overall_agreement_up <= 1)
    expect_true(oa$overall_agreement_lo <= oa$overall_agreement)
    expect_true(oa$overall_agreement <= oa$overall_agreement_up)
  })
})

test_that("ctl handles dataset2", {
  expect_no_error(
    ctl(
      ctl_dataset2,
      idvar = "id",
      ynew = "y1",
      yref = "y2",
      intercept = 1,
      slope = 0.2
    )
  )
})

test_that("ctl throws error if outputs_path is missing when writing files", {
  expect_error(
    ctl(
      ctl_dataset1,
      idvar = "id",
      ynew = "y1",
      yref = "y2",
      intercept = 5,
      slope = 0.15,
      results_to_file = TRUE
    ),
    "Please provide the `outputs_path` argument"
  )

  expect_error(
    ctl(
      ctl_dataset1,
      idvar = "id",
      ynew = "y1",
      yref = "y2",
      intercept = 5,
      slope = 0.15,
      plots_to_file = TRUE
    ),
    "Please provide the `outputs_path` argument"
  )
})
