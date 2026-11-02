test_that("convert_powers_list_to_matrix uses safe power-slot sequences", {
  expect_error(
    mfp2:::convert_powers_list_to_matrix(list()),
    "must contain at least one variable"
  )
  expect_error(
    mfp2:::convert_powers_list_to_matrix(list(x = numeric())),
    "must contain at least one power"
  )

  out <- mfp2:::convert_powers_list_to_matrix(
    list(x1 = 1, x2 = c(-2, 0.5), x3 = NA_real_)
  )

  expect_identical(dim(out), c(3L, 2L))
  expect_identical(colnames(out), c("power1", "power2"))
  expect_equal(unname(out["x1", ]), c(1, NA_real_))
  expect_equal(unname(out["x2", ]), c(-2, 0.5))
  expect_equal(unname(out["x3", ]), c(NA_real_, NA_real_))
})
