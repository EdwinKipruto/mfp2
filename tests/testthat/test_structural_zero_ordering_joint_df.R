# Tests for blockwise likelihood-ratio ordering of catzero/retained-spike terms.

# The binary indicator and positive-part continuous column form one conceptual
# term. Their ordering test must use the fitted rank difference for the whole
# block, not one df and not raw block width unconditionally.
test_that("structural-zero ordering uses joint fitted rank difference", {
  x <- cbind(
    exposure = c(0, 0, 1, 2),
    exposure_bin = c(1, 1, 0, 0),
    z = c(1, 2, 3, 4)
  )
  mapping <- list(
    exposure = c("exposure", "exposure_bin"),
    z = "z"
  )

  testthat::local_mocked_bindings(
    fit_model = function(x, ...) {
      remaining <- colnames(x)
      if (identical(remaining, "z")) {
        return(list(logl = 95, df = 1))
      }
      if (identical(remaining, c("exposure", "exposure_bin"))) {
        return(list(logl = 96, df = 2))
      }
      stop("Unexpected reduced design in structural-zero ordering test.")
    },
    .package = "mfp2"
  )

  ordered <- order_variables_by_significance(
    xorder = "ascending",
    x = x,
    y = rep(0, nrow(x)),
    family = stats::gaussian(),
    family_string = "gaussian",
    weights = NULL,
    offset = NULL,
    strata = NULL,
    method = NULL,
    control = NULL,
    nocenter = NULL,
    full_reference = list(logl = 100, df = 3),
    term_to_columns = mapping
  )

  # exposure: LR = 10 on 2 df; z: LR = 8 on 1 df. The 1-df z term is more
  # significant, so it is visited first. Treating exposure as one df reverses
  # this order and would fail the test.
  expect_identical(ordered, c("z", "exposure"))
})

# If an augmented term contributes no fitted rank after aliasing, it must remain
# in the visiting order but receive NA significance and sort after valid terms.
test_that("non-identifiable structural-zero block is ordered after valid terms", {
  x <- cbind(
    exposure = c(0, 0, 1, 2),
    exposure_bin = c(1, 1, 0, 0),
    z = c(1, 2, 3, 4)
  )
  mapping <- list(
    exposure = c("exposure", "exposure_bin"),
    z = "z"
  )

  testthat::local_mocked_bindings(
    fit_model = function(x, ...) {
      remaining <- colnames(x)
      if (identical(remaining, "z")) {
        # Dropping exposure changes no fitted rank: invalid LRT for exposure.
        return(list(logl = 100, df = 3))
      }
      if (identical(remaining, c("exposure", "exposure_bin"))) {
        return(list(logl = 96, df = 2))
      }
      stop("Unexpected reduced design in structural-zero alias test.")
    },
    .package = "mfp2"
  )

  ordered <- order_variables_by_significance(
    xorder = "ascending",
    x = x,
    y = rep(0, nrow(x)),
    family = stats::gaussian(),
    family_string = "gaussian",
    weights = NULL,
    offset = NULL,
    strata = NULL,
    method = NULL,
    control = NULL,
    nocenter = NULL,
    full_reference = list(logl = 100, df = 3),
    term_to_columns = mapping
  )

  expect_identical(ordered, c("z", "exposure"))
})
