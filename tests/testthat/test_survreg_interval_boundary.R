# Regression tests for log-time survreg censoring-boundary preparation.


prepare_log_time_response <- function(y, dist) {
  prepare_survreg_family(
    family = survreg_family(dist = dist),
    y = y,
    weights = rep(1, NROW(y))
  )$prepared$y
}

prepare_loglogistic_response <- function(y) {
  prepare_log_time_response(y, dist = "loglogistic")
}


test_that("zero-lower log-time intervals equal left censoring at upper", {
  upper <- c(2, 4, 8)
  manual_left <- survival::Surv(
    time = upper,
    event = rep(0, length(upper)),
    type = "left"
  )
  boundary_inputs <- list(
    interval = survival::Surv(
      time = rep(0, length(upper)),
      time2 = upper,
      event = rep(3, length(upper)),
      type = "interval"
    ),
    interval2 = survival::Surv(
      time = rep(0, length(upper)),
      time2 = upper,
      type = "interval2"
    )
  )

  for (dist in c("weibull", "exponential", "rayleigh", "lognormal", "loglogistic")) {
    expected <- prepare_log_time_response(manual_left, dist = dist)

    for (input_name in names(boundary_inputs)) {
      prepared <- prepare_log_time_response(boundary_inputs[[input_name]], dist = dist)
      context <- paste(dist, input_name)
      expect_true(all(is.finite(prepared)), info = context)
      expect_equal(
        unname(prepared),
        unname(expected),
        tolerance = 0,
        info = context
      )
    }
  }
})


test_that("log-logistic left and positive interval responses retain correct coding", {
  left <- survival::Surv(
    time = c(2, 5),
    event = c(0, 1),
    type = "left"
  )
  expect_equal(
    unname(prepare_loglogistic_response(left)),
    unname(cbind(log(c(2, 5)), c(2, 1))),
    tolerance = 0
  )

  lower <- c(1, 2)
  upper <- c(3, 5)
  expected_interval <- cbind(log(lower), log(upper), c(3, 3))
  interval_inputs <- list(
    interval = survival::Surv(
      time = lower,
      time2 = upper,
      event = c(3, 3),
      type = "interval"
    ),
    interval2 = survival::Surv(
      time = lower,
      time2 = upper,
      type = "interval2"
    )
  )

  for (input_name in names(interval_inputs)) {
    expect_equal(
      unname(prepare_loglogistic_response(interval_inputs[[input_name]])),
      unname(expected_interval),
      tolerance = 0,
      info = input_name
    )
  }
})
