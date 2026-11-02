make_discrete_mfpi_data <- function(n = 240L, seed = 867L) {
  set.seed(seed)
  group <- factor(rep(c("control", "treated"), each = n / 2L))
  binary <- rep(c(0, 1, 1, 0), length.out = n)
  stage <- factor(rep(c("I", "II", "III", "II"), length.out = n),
                  levels = c("I", "II", "III"))
  x <- stats::runif(n, 0.5, 4)
  y <- 1 + 0.4 * (group == "treated") - 0.6 * binary +
    1.2 * (group == "treated") * binary +
    0.5 * (stage == "II") - 0.3 * (stage == "III") +
    0.8 * (group == "treated") * (stage == "II") -
    0.5 * (group == "treated") * (stage == "III") +
    0.25 * x + stats::rnorm(n, sd = 0.35)
  data.frame(y, group, binary, stage, x)
}

test_that("binary linear MFPI matches an explicit Gaussian interaction", {
  dat <- make_discrete_mfpi_data()

  fit <- mfpi(
    y ~ group + binary + x,
    data = dat,
    group_var = "group",
    interaction_vars = "binary",
    select = 1,
    p_interact = 1,
    center = TRUE,
    verbose = FALSE
  )
  reference <- stats::lm(y ~ group * binary + x, data = dat)
  stored <- fit$var_winners$binary$fit

  expect_identical(fit$interaction_specs$binary$kind, "binary")
  expect_equal(unname(fit$scale["binary"]), 1)
  expect_equal(unname(fit$shift["binary"]), 0)
  expect_null(stored$center_vals)
  expect_equal(stored$test_results$df$interaction_only, 1)
  expect_equal(stored$test_results$interaction_model$logl,
               as.numeric(stats::logLik(reference)), tolerance = 1e-7)

  nd <- dat[1:30, c("group", "binary", "x")]
  observed <- predict(
    fit, terms = "binary", model = "all", type = "response", newdata = nd
  )$predictions$fit
  expected <- stats::predict(reference, newdata = nd)
  expect_equal(observed, unname(expected), tolerance = 1e-7)

  beta <- stored$test_results$interaction_model$coefficients
  groups <- stored$coefficient_groups
  expect_equal(unname(beta[groups[[1L]]]),
               unname(stats::coef(reference)["binary"]), tolerance = 1e-7)
  expect_equal(
    unname(beta[groups[[2L]]]),
    unname(stats::coef(reference)["binary"] +
             stats::coef(reference)["grouptreated:binary"]),
    tolerance = 1e-7
  )
})

test_that("categorical linear MFPI uses the complete contrast block", {
  dat <- make_discrete_mfpi_data()

  fit <- mfpi(
    y ~ group + stage + x,
    data = dat,
    group_var = "group",
    interaction_vars = "stage",
    interaction_forms = c(stage = "linear"),
    select = 1,
    p_interact = 1,
    verbose = FALSE
  )
  reference <- stats::lm(y ~ group * stage + x, data = dat)
  stored <- fit$var_winners$stage$fit

  expect_identical(fit$interaction_specs$stage$kind, "categorical")
  expect_equal(fit$interaction_specs$stage$width, 2)
  expect_equal(fit$interaction_specs$stage$level_labels,
               c("I", "II", "III"))
  expect_null(stored$center_vals)
  expect_equal(stored$test_results$df$interaction_only, 2)
  expect_equal(stored$test_results$interaction_model$logl,
               as.numeric(stats::logLik(reference)), tolerance = 1e-7)

  nd <- dat[1:36, c("group", "stage", "x")]
  observed <- predict(
    fit, terms = "stage", model = "all", type = "response", newdata = nd
  )$predictions$fit
  expected <- stats::predict(reference, newdata = nd)
  expect_equal(observed, unname(expected), tolerance = 1e-7)

  functions <- predict(
    fit, terms = "stage", model = "all", type = "function", grid = TRUE
  )$functions
  expect_setequal(unique(functions$x), levels(dat$stage))
})

test_that("matrix categorical MFPI matches an explicit Gaussian interaction", {
  dat <- make_discrete_mfpi_data()
  stage_matrix <- stats::model.matrix(~ stage, dat)[, -1L, drop = FALSE]
  x_matrix <- cbind(
    group = as.integer(dat$group) - 1L,
    stage_matrix,
    x = dat$x
  )

  fit <- mfpi(
    x_matrix,
    dat$y,
    group_var = "group",
    interaction_vars = "stage",
    term_groups = list(stage = colnames(stage_matrix)),
    select = 1,
    p_interact = 1,
    verbose = FALSE
  )
  reference <- stats::lm(y ~ group * stage + x, data = dat)
  stored <- fit$var_winners$stage$fit

  expect_identical(fit$interaction_specs$stage$kind, "categorical")
  expect_equal(fit$interaction_specs$stage$columns, colnames(stage_matrix))
  expect_equal(stored$test_results$df$interaction_only, 2)
  expect_equal(
    stored$test_results$interaction_model$logl,
    as.numeric(stats::logLik(reference)),
    tolerance = 1e-7
  )

  nd_rows <- 1:36
  nd <- x_matrix[nd_rows, , drop = FALSE]
  observed <- predict(
    fit, terms = "stage", model = "all", type = "response", newdata = nd
  )$predictions$fit
  expected <- stats::predict(reference, newdata = dat[nd_rows, ])
  expect_equal(observed, unname(expected), tolerance = 1e-7)
})

test_that("categorical FP requests are reset to linear with a warning", {
  dat <- make_discrete_mfpi_data()

  expect_warning(
    fit <- mfpi(
      y ~ group + stage + x,
      data = dat,
      group_var = "group",
      interaction_vars = "stage",
      interaction_forms = c(stage = "fp1"),
      verbose = FALSE
    ),
    "reset to 'linear': stage"
  )
  expect_identical(fit$interaction_forms, c(stage = "linear"))
  expect_identical(fit$interaction_specs$stage$form, "linear")
  expect_equal(fit$all_model_metrics$df_int, 2)
})

test_that("omitted and partial forms match explicit mixed-variable models", {
  dat <- make_discrete_mfpi_data()
  args <- list(
    y ~ group + binary + stage + x,
    data = dat,
    group_var = "group",
    interaction_vars = c("binary", "stage", "x"),
    df = 1,
    select = 1,
    powers = list(x = c(0, 1)),
    flex = "flex1",
    verbose = FALSE
  )
  expected_forms <- c(binary = "linear", stage = "linear", x = "fp1")
  explicit <- do.call(mfpi, c(args, list(interaction_forms = expected_forms)))

  for (forms in list(NULL, c(x = "fp1"))) {
    expect_no_warning(
      fit <- do.call(mfpi, c(args, list(interaction_forms = forms)))
    )
    expect_identical(fit$interaction_forms, expected_forms)
    expect_identical(
      vapply(fit$interaction_specs, `[[`, character(1L), "form"),
      expected_forms
    )
    expect_equal(fit$all_model_metrics, explicit$all_model_metrics)
    expect_equal(unname(fit$shift[c("binary", "stageII", "stageIII")]),
                 c(0, 0, 0))
    expect_equal(unname(fit$scale[c("binary", "stageII", "stageIII")]),
                 c(1, 1, 1))
    # Linear adjustment does not force a continuous interaction to be linear.
    expect_equal(fit$adjustment_model$fp_terms["x", "df_final"], 1)
    expect_identical(
      fit$all_model_metrics$type[fit$all_model_metrics$variable == "x"],
      "fp1"
    )
  }
})

test_that("explicit FP1 and FP2 requests reset all discrete terms together", {
  dat <- make_discrete_mfpi_data()

  for (form in c("fp1", "fp2")) {
    expect_warning(
      fit <- mfpi(
        y ~ group + binary + stage + x,
        data = dat,
        group_var = "group",
        interaction_vars = c("binary", "stage", "x"),
        interaction_forms = c(stage = form, binary = form, x = "linear"),
        df = 1,
        select = 1,
        verbose = FALSE
      ),
      "reset to 'linear': stage, binary"
    )
    expect_identical(
      fit$interaction_forms,
      c(binary = "linear", stage = "linear", x = "linear")
    )
    expect_equal(
      fit$all_model_metrics$df_int[
        match(c("binary", "stage", "x"), fit$all_model_metrics$variable)
      ],
      c(1, 2, 1)
    )
  }
})

test_that("form validation rejects malformed requests before resetting", {
  specs <- list(binary = list(discrete = TRUE), x = list(discrete = FALSE))
  resolve <- function(forms) resolve_mfpi_interaction_forms(forms, specs)

  expect_error(resolve(list(binary = "fp1")), "named character vector")
  expect_error(resolve("fp1"), "Every entry.*must be named")
  expect_error(resolve(setNames("fp1", NA_character_)), "must be named")
  expect_error(resolve(c(binary = "fp1", binary = "fp2")), "must be unique")
  expect_error(resolve(c(other = "fp1")), "not in `interaction_vars`")
  expect_error(resolve(c(binary = NA_character_)), "Invalid value")
  expect_error(resolve(c(binary = "spline")), "Invalid value")
})

test_that("categorical interaction df scales across multinomial logits", {
  expect_equal(
    mfp2:::interaction_model_df(
      n_groups = 2,
      degree = 0,
      n_logits = 3,
      flex = "flex0",
      linear_width = 2
    )$dfint,
    6
  )
})

test_that("binary multinomial MFPI matches nnet::multinom", {
  set.seed(868L)
  n <- 360L
  dat <- data.frame(
    group = factor(rep(c("control", "active"), each = n / 2L),
                   levels = c("control", "active")),
    binary = rep(c(0, 1, 1, 0), length.out = n),
    x = stats::rnorm(n)
  )
  eta_b <- with(
    dat,
    -0.35 + 0.55 * (group == "active") - 0.45 * binary +
      0.80 * (group == "active") * binary + 0.25 * x
  )
  eta_c <- with(
    dat,
    0.20 - 0.30 * (group == "active") + 0.35 * binary -
      0.50 * (group == "active") * binary - 0.20 * x
  )
  denominator <- 1 + exp(eta_b) + exp(eta_c)
  probabilities <- cbind(
    A = 1 / denominator,
    B = exp(eta_b) / denominator,
    C = exp(eta_c) / denominator
  )
  draw <- apply(
    probabilities,
    1L,
    function(probability) sample.int(3L, 1L, prob = probability)
  )
  dat$y <- factor(colnames(probabilities)[draw], levels = colnames(probabilities))

  fit <- mfpi(
    y ~ group + binary + x,
    data = dat,
    family = multinomial_family(),
    group_var = "group",
    interaction_vars = "binary",
    interaction_forms = c(binary = "linear"),
    df = 1,
    select = 1,
    p_interact = 1,
    verbose = FALSE
  )
  reference <- nnet::multinom(
    y ~ group * binary + x,
    data = dat,
    trace = FALSE
  )
  stored <- fit$var_winners$binary$fit$test_results$interaction_model
  stored_coef <- stored$coefficient_matrix
  reference_coef <- stats::coef(reference)

  expect_equal(
    fit$var_winners$binary$fit$test_results$df$interaction_only,
    2
  )
  expect_equal(stored$logl, as.numeric(stats::logLik(reference)), tolerance = 1e-6)
  expect_equal(
    unname(stored$fit$fitted.values),
    unname(stats::predict(reference, type = "probs")),
    tolerance = 1e-4
  )
  expect_equal(
    unname(stored_coef[, "binary01"]),
    unname(reference_coef[, "binary"]),
    tolerance = 1e-4
  )
  expect_equal(
    unname(stored_coef[, "binary11"]),
    unname(reference_coef[, "binary"] +
             reference_coef[, "groupactive:binary"]),
    tolerance = 1e-4
  )
})

test_that("discrete MFPI plots use points instead of continuous curves", {
  skip_if_not_installed("ggplot2")
  dat <- make_discrete_mfpi_data()
  fit <- mfpi(
    y ~ group + stage + x,
    data = dat,
    group_var = "group",
    interaction_vars = "stage",
    interaction_forms = c(stage = "linear"),
    select = 1,
    p_interact = 1,
    verbose = FALSE
  )

  plots <- plot(
    fit,
    terms = "stage",
    plot_type = "fitted",
    show_ci_fitted = FALSE,
    auto_print = FALSE
  )
  p <- plots$stage[[1L]]
  geoms <- vapply(p$layers, function(layer) class(layer$geom)[1L], character(1L))
  expect_true("GeomPoint" %in% geoms)
  expect_false("GeomLine" %in% geoms)

  difference <- plot(
    fit, terms = "stage", plot_type = "difference",
    show_ci_maineffect = TRUE, auto_print = FALSE
  )$stage[[1L]]
  difference_geoms <- vapply(difference$layers, function(layer) {
    class(layer$geom)[1L]
  }, character(1L))
  expect_true("GeomRect" %in% difference_geoms)
  expect_true("GeomHline" %in% difference_geoms)
  expect_no_error(ggplot2::ggplot_build(difference))
})
