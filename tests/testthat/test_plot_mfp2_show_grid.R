grid_is_blank <- function(p) {
  inherits(ggplot2::calc_element("panel.grid.major.x", ggplot2::theme_bw() + p$theme),
           "element_blank")
}

test_that("plot.mfp2() hides grid lines by default and shows them with show_grid = TRUE", {
  set.seed(11)
  d <- data.frame(x1 = runif(150, 1, 5), x2 = runif(150, 1, 5))
  d$y <- log(d$x1) + d$x2 + rnorm(150)
  fit <- mfp2(y ~ fp(x1, df = 2) + fp(x2, df = 2, select = 1), data = d, verbose = FALSE)

  default_plots <- plot(fit)
  grid_plots <- plot(fit, show_grid = TRUE)
  expect_true(all(vapply(default_plots, grid_is_blank, logical(1))))
  expect_false(any(vapply(grid_plots, grid_is_blank, logical(1))))
})

test_that("multinomial plots honour show_grid", {
  set.seed(12)
  d <- data.frame(x1 = runif(300, 1, 5))
  eta2 <- 0.8 * log(d$x1); eta3 <- -0.5 + 0.3 * d$x1
  pr <- cbind(1, exp(eta2), exp(eta3)); pr <- pr / rowSums(pr)
  d$y <- factor(apply(pr, 1, function(p) sample(1:3, 1, prob = p)))
  fit <- mfp2(y ~ fp(x1, df = 2, select = 1), data = d,
              family = multinomial_family(), verbose = FALSE)

  expect_true(grid_is_blank(plot(fit)[[1]]))
  expect_false(grid_is_blank(plot(fit, show_grid = TRUE)[[1]]))
})

test_that("show_grid must be a single non-missing logical value", {
  set.seed(13)
  d <- data.frame(x1 = runif(80, 1, 5)); d$y <- log(d$x1) + rnorm(80)
  fit <- mfp2(y ~ fp(x1, df = 2, select = 1), data = d, verbose = FALSE)
  expect_error(plot(fit, show_grid = NA), "show_grid")
  expect_error(plot(fit, show_grid = c(TRUE, FALSE)), "show_grid")
  expect_error(plot(fit, show_grid = "yes"), "show_grid")
})

test_that("deprecated fracplot() forwards show_grid", {
  set.seed(14)
  d <- data.frame(x1 = runif(80, 1, 5)); d$y <- log(d$x1) + rnorm(80)
  fit <- mfp2(y ~ fp(x1, df = 2, select = 1), data = d, verbose = FALSE)
  p <- suppressWarnings(fracplot(fit, show_grid = TRUE))
  expect_false(grid_is_blank(p[[1]]))
})
