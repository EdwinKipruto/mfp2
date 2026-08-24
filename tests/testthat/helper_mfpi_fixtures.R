# Shared MFPI fitting fixtures.

make_mfpi_factor_data <- function(n = 240L, ordered_stage = FALSE) {
  set.seed(9101)

  trt <- factor(rep(c("control", "treated"), length.out = n))
  stage_values <- rep(c("I", "II", "III"), length.out = n)
  stage <- if (ordered_stage) {
    ordered(stage_values, levels = c("I", "II", "III"))
  } else {
    factor(stage_values, levels = c("I", "II", "III"))
  }

  x <- runif(n, 1, 8)
  z <- rnorm(n)
  stage_effect <- c(I = 0, II = 0.7, III = -0.5)[as.character(stage)]

  data.frame(
    y = 1 + 0.4 * x + 0.9 * x * (trt == "treated") +
      stage_effect + 0.2 * z + rnorm(n, sd = 0.25),
    trt = trt,
    x = x,
    stage = stage,
    z = z
  )
}


# Test purpose: Verifies that mfpi.default() uses the same partial named
