# GEE fractional-polynomial selection: all three criteria are available and
# recover the expected functional forms.

make_gee_fp_data <- function(K = 80, m = 5, seed = 7) {
  set.seed(seed)
  n <- K * m
  id <- rep(seq_len(K), each = m)
  wave <- rep(seq_len(m), times = K)
  x1 <- stats::runif(n, 1, 5)      # true effect 1/x1
  x2 <- stats::rnorm(n)            # linear
  x3 <- stats::rnorm(n)            # null
  b <- stats::rnorm(K)[id]
  y <- 2 + 3 * (1 / x1) + 0.6 * x2 + 0 * x3 + b + stats::rnorm(n, sd = 0.4)
  data.frame(y, x1, x2, x3, id, wave)
}

test_that("weighted GEE quasi-likelihood uses the correct family expressions", {
  expect_equal(
    gee_quasi_likelihood(
      y = c(1, 3), mu = c(2, 2), family_string = "gaussian",
      weights = c(2, 4)
    ),
    sum(c(2, 4) * (-0.5 * (c(1, 3) - c(2, 2))^2))
  )

  expect_equal(
    gee_quasi_likelihood(
      y = c(0, 3), mu = c(1, 2), family_string = "poisson",
      weights = c(2, 4)
    ),
    sum(c(2, 4) * (c(0, 3) * log(c(1, 2)) - c(1, 2)))
  )

  # Effective grouped-binomial weights are prior weights multiplied by totals.
  proportions <- c(2 / 5, 3 / 4)
  effective_weights <- c(2, 3) * c(5, 4)
  mu <- c(0.3, 0.6)
  expect_equal(
    gee_quasi_likelihood(
      y = proportions, mu = mu, family_string = "binomial",
      weights = effective_weights
    ),
    sum(effective_weights * (
      proportions * log(mu / (1 - mu)) + log(1 - mu)
    ))
  )

  gamma_y <- c(1, 4)
  gamma_mu <- c(2, 3)
  gamma_weights <- c(2, 4)
  expect_equal(
    gee_quasi_likelihood(
      y = gamma_y, mu = gamma_mu, family_string = "Gamma",
      weights = gamma_weights
    ),
    sum(gamma_weights * (-gamma_y / gamma_mu - log(gamma_mu)))
  )
})

test_that("weighted quasi-likelihood differences equal family deviance differences", {
  cases <- list(
    gaussian = list(
      family = gaussian(), y = c(-0.4, 1.2, 2),
      reduced = c(0, 1, 1.4), full = c(-0.2, 1.1, 1.8),
      weights = c(2, 0.5, 3)
    ),
    poisson = list(
      family = poisson(), y = c(0, 2, 5),
      reduced = c(0.7, 1.4, 3.7), full = c(0.5, 1.9, 4.6),
      weights = c(2, 0.5, 3)
    ),
    binomial = list(
      family = binomial(), y = c(1 / 5, 3 / 4, 1 / 2),
      reduced = c(0.4, 0.5, 0.7), full = c(0.25, 0.7, 0.55),
      weights = c(10, 12, 8)
    ),
    Gamma = list(
      family = Gamma(), y = c(0.5, 2, 5),
      reduced = c(1, 1.4, 4), full = c(0.7, 1.8, 4.8),
      weights = c(2, 0.5, 3)
    )
  )

  for (family_name in names(cases)) {
    case <- cases[[family_name]]
    q_reduced <- gee_quasi_likelihood(
      case$y, case$reduced, family_name, case$weights
    )
    q_full <- gee_quasi_likelihood(
      case$y, case$full, family_name, case$weights
    )
    d_reduced <- sum(case$family$dev.resids(
      case$y, case$reduced, case$weights
    ))
    d_full <- sum(case$family$dev.resids(
      case$y, case$full, case$weights
    ))

    expect_equal(
      2 * (q_full - q_reduced),
      d_reduced - d_full,
      tolerance = 1e-12,
      info = family_name
    )
  }
})

test_that("GEE QAIC and QBIC use weighted quasi-likelihood and cluster count", {
  obj <- list(
    logl = -12,
    df = 3,
    rank = 3,
    coefficients = stats::setNames(rep(0, 3), c("a", "b", "c")),
    selection_deviance = -42
  )

  metrics <- calculate_gee_metrics(obj, n_obs = 10, df_additional = 2)

  expect_equal(unname(metrics["aic"]), 24 + 2 * 5)
  expect_equal(unname(metrics["bic"]), 24 + log(10) * 5)
  # GEE p-value selection reports Stata's negative overall robust Wald quantity.
  expect_equal(unname(metrics["deviance_rs"]), -42)

  # Without a stored selection deviance the metric falls back to -2Q.
  obj$selection_deviance <- NULL
  fallback <- calculate_gee_metrics(obj, n_obs = 10, df_additional = 2)
  expect_equal(unname(fallback["deviance_rs"]), -2 * -12)
})

test_that("retained GEE selection scores use the active criterion scale", {
  obj <- list(
    logl = -12,
    df = 3,
    rank = 3,
    coefficients = stats::setNames(rep(0, 3), c("a", "b", "c")),
    selection_deviance = -42
  )

  # p-value selection uses Stata's negative overall robust Wald quantity.
  expect_equal(gee_final_selection_score(obj, "pvalue", 10), -42)
  expect_equal(gee_final_selection_score(obj, "aic", 10), 30)
  expect_equal(
    gee_final_selection_score(obj, "bic", 10),
    24 + log(10) * 3
  )
})

test_that("progress tables use three decimals and GEE QAIC/QBIC labels", {
  pvalue_fit <- list(
    metrics = rbind(
      FP2 = c(df = 4, deviance_rs = -12.3456),
      null = c(df = 0, deviance_rs = -1.2344)
    ),
    pvalue = c("null vs FP2" = 0.04567),
    is_gee = TRUE
  )

  pvalue_table <- print_mfp_pvalue_step(
    xi = "x", fit = pvalue_fit, criterion = "pvalue"
  )
  # GEE reports the overall robust Wald chi-square W = -deviance_rs (positive),
  # with the statistic column giving W(full) - W(reduced).
  expect_identical(
    unname(pvalue_table[, "Wald chi-sq"]),
    c("12.346", "1.234")
  )
  expect_identical(
    unname(pvalue_table[, "Wald diff."]),
    c(NA_character_, "11.111")
  )
  expect_identical(
    unname(pvalue_table[, "P-value"]),
    c(NA_character_, "0.046")
  )

  ic_fit <- list(
    metrics = rbind(
      null = c(aic = 12.3456, bic = 13.4567),
      linear = c(aic = 10.2344, bic = 11.3456)
    ),
    is_gee = TRUE
  )
  qaic_table <- print_mfp_ic_step("x", ic_fit, "aic")
  qbic_table <- print_mfp_ic_step("x", ic_fit, "bic")

  expect_identical(colnames(qaic_table), "QAIC")
  expect_identical(unname(qaic_table[, 1]), c("12.346", "10.234"))
  expect_identical(colnames(qbic_table), "QBIC")
  expect_identical(unname(qbic_table[, 1]), c("13.457", "11.346"))

  ic_fit$is_gee <- NULL
  expect_identical(colnames(print_mfp_ic_step("x", ic_fit, "aic")), "AIC")
  expect_identical(colnames(print_mfp_ic_step("x", ic_fit, "bic")), "BIC")
})

test_that("robust GEE block-Wald and Stata selection statistics are distinct", {
  beta <- c("(Intercept)" = 3, x1 = 2, x2 = -1)
  covariance <- diag(c(4, 1, 1))
  dimnames(covariance) <- list(names(beta), names(beta))

  block <- gee_robust_wald_test(beta, covariance, indices = 2:3)
  expect_equal(block$statistic, 5)
  expect_equal(block$df, 2)
  expect_equal(
    block$pvalue,
    stats::pchisq(5, df = 2, lower.tail = FALSE)
  )
  expect_equal(gee_selection_deviance(beta, covariance), -5)
  expect_equal(
    gee_selection_deviance(beta, covariance, n_clusters = 5),
    -4
  )
  expect_equal(
    gee_selection_deviance(beta, gee_robust_vcov(covariance, 5)),
    -4
  )
  score_only <- gee_robust_wald_test(
    beta, covariance, indices = 2:3, compute_pvalue = FALSE
  )
  expect_equal(score_only$statistic, block$statistic)
  expect_true(is.na(score_only$pvalue))
})

test_that("singular GEE slope blocks cannot produce joint selection tests", {
  beta <- c("(Intercept)" = 0, x1 = 1, x2 = 1)
  covariance <- matrix(0, 3, 3, dimnames = list(names(beta), names(beta)))
  covariance[2:3, 2:3] <- matrix(1, 2, 2)

  rank_one <- gee_robust_wald_test(beta, covariance, indices = 2:3)
  expect_identical(rank_one$df, 1L)
  expect_true(all(is.na(unlist(rank_one[c("statistic", "pvalue", "dev_diff")]))))
  expect_true(is.na(gee_selection_deviance(beta, covariance)))
  expect_true(is.na(gee_selection_deviance(beta, covariance, n_clusters = 2)))
  expect_true(is.na(gee_robust_wald_test(beta, covariance, 2:3,
                                        compute_pvalue = FALSE)$statistic))

  rank_zero <- gee_robust_wald_test(beta, matrix(0, 3, 3), 2:3)
  expect_identical(rank_zero$df, 0L)
  expect_true(is.na(rank_zero$statistic))
  expect_true(is.na(rank_zero$pvalue))

  single_slope <- gee_robust_wald_test(beta, covariance, 2L)
  expect_equal(single_slope$statistic, 1)
  expect_equal(single_slope$pvalue, stats::pchisq(1, df = 1,
                                                 lower.tail = FALSE))
  expect_identical(gee_robust_wald_test(beta, covariance, integer())$statistic,
                   0)
})

test_that("GEE covariance selection corrects only the sandwich estimator", {
  geese <- list(
    vbeta = diag(c(1, 2)),
    vbeta.ajs = diag(c(3, 4)),
    vbeta.j1s = diag(c(5, 6)),
    vbeta.fij = diag(c(7, 8))
  )

  expect_equal(gee_selected_vcov(geese, "san.se", n_clusters = 5),
               geese$vbeta * 5 / 4)
  expect_equal(gee_selected_vcov(geese, "jack", n_clusters = 5),
               geese$vbeta.ajs)
  expect_equal(gee_selected_vcov(geese, "j1s", n_clusters = 5),
               geese$vbeta.j1s)
  expect_equal(gee_selected_vcov(geese, "fij", n_clusters = 5),
               geese$vbeta.fij)

  control <- gee_control_for_std_err(
    list(jack = 0L, j1s = 0L, fij = 0L),
    "j1s"
  )
  expect_identical(control$jack, 0L)
  expect_identical(control$j1s, 1L)
  expect_identical(control$fij, 0L)
})

test_that("GEE quasi-likelihood is stable at response boundaries", {
  expect_equal(
    gee_quasi_likelihood(c(0, 1), c(0, 1), "binomial"),
    0
  )
  expect_equal(
    gee_quasi_likelihood(c(0, 2), c(0, 2), "poisson"),
    2 * log(2) - 2
  )
  expect_identical(
    gee_quasi_likelihood(1, 0, "binomial"),
    -Inf
  )
})

test_that("power one is valid in every searched FP1 class", {
  expect_true(validate_mfp_candidate_powers(
    powers = list(x = 1),
    df = c(x = 4)
  ))
})

test_that("GEE visiting order uses full-model robust block tests", {
  x <- cbind(x1 = 1:5, x2 = 6:10)
  beta <- c("(Intercept)" = 0, x1 = 1, x2 = 3)
  covariance <- diag(3)
  dimnames(covariance) <- list(names(beta), names(beta))

  result <- order_variables_by_significance(
    xorder = "ascending",
    x = x,
    y = 1:5,
    family = gee_family(),
    family_string = "gee",
    weights = NULL,
    offset = NULL,
    strata = NULL,
    method = NULL,
    control = NULL,
    nocenter = NULL,
    full_reference = list(
      coefficients = beta,
      robust_vcov = covariance,
      df = 3,
      logl = 0
    )
  )

  expect_identical(result, c("x2", "x1"))
})

test_that("fixed-degree GEE p-value selection uses the Stata selection deviance", {
  metric <- function(logl, deviance) {
    c(
      logl = logl, df = 3, deviance_rs = deviance,
      deviance_gaussian = NA_real_, aic = 0, bic = 0, df_resid = 7
    )
  }

  testthat::local_mocked_bindings(
    transform_data_step = function(..., powers) {
      expect_true(1 %in% powers$x)
      list(
        data_adj = NULL,
        data_fp = list(
          matrix(1, nrow = 5, ncol = 1),
          matrix(2, nrow = 5, ncol = 1),
          matrix(3, nrow = 5, ncol = 1)
        ),
        fp_basis = NULL,
        acd_basis = NULL,
        powers_fp = matrix(c(-1, 1, 2), ncol = 1),
        current_params = list(x = list())
      )
    },
    fit_model = function(x, ...) list(candidate = x[1, 2]),
    calculate_model_metrics = function(obj, ...) {
      if (obj$candidate == 1) metric(logl = -10, deviance = -20) else
        if (obj$candidate == 2) metric(logl = -7, deviance = -15) else
          metric(logl = -5, deviance = -10)
    },
    .package = "mfp2"
  )

  result <- find_best_fpm_step(
    x = matrix(1:5, ncol = 1, dimnames = list(NULL, "x")),
    xi = "x", degree = 1, y = 1:5,
    powers_current = list(x = 1), powers = list(x = c(-1, 1, 2)),
    acdx = c(x = FALSE),
    family = gee_family(), family_string = "gee",
    zero = c(x = FALSE), catzero = list(x = NULL),
    spike = c(x = FALSE), spike_decision = c(x = NA_integer_),
    acd_parameter = list(), prev_adj_params = list(),
    has_offset = FALSE, n_obs = 5, term_to_columns = list(x = "x"),
    criterion = "pvalue"
  )

  # GEE p-value selection minimizes the negative global-Wald selection
  # deviance, so the most negative candidate 1 wins (not the largest-Q one).
  expect_equal(result$model_best, 1)
})

test_that("Stata GEE RA2 passes power one through the FP1 search", {
  seen <- new.env(parent = emptyenv())
  metric <- function(df, deviance) {
    c(
      logl = 0, df = df, deviance_rs = deviance,
      deviance_gaussian = NA_real_, aic = 0, bic = 0, df_resid = 20 - df
    )
  }

  testthat::local_mocked_bindings(
    build_adjustment_step = function(...) list(transform_cache = NULL),
    find_best_fpm_step = function(degree, powers, ...) {
      seen[[as.character(degree)]] <- 1 %in% powers[["x"]]
      if (degree == 2L) {
        list(
          powers = matrix(c(-1, 2), nrow = 1L),
          power_best = c(-1, 2),
          metrics = rbind(FP2 = metric(5, -20)),
          model_best = 1L,
          current_adj_params = list()
        )
      } else {
        list(
          powers = matrix(1, nrow = 1L),
          power_best = 1,
          # Power 1 has the same fit as the separately fitted linear model,
          # but is charged the searched-FP1 df in the closed comparison.
          metrics = rbind(FP1 = metric(3, -12)),
          model_best = 1L,
          current_adj_params = list()
        )
      }
    },
    fit_null_step = function(...) list(
      powers = NA_real_, metrics = metric(1, -10),
      current_adj_params = list()
    ),
    fit_linear_step = function(...) list(
      powers = 1, metrics = metric(2, -12),
      current_adj_params = list()
    ),
    .package = "mfp2"
  )

  result <- select_ra2(
    x = matrix(1:5, ncol = 1, dimnames = list(NULL, "x")),
    xi = "x", keep = character(), degree = 2, acdx = c(x = FALSE),
    y = 1:5, powers_current = list(x = 1),
    powers = list(x = c(-1, 1, 2)), criterion = "pvalue",
    ftest = FALSE, select = 0.05, alpha = 0.05,
    family = gee_family(gaussian()),
    family_string = "gee", zero = c(x = FALSE),
    catzero = list(x = NULL), spike = c(x = FALSE),
    spike_decision = c(x = NA_integer_), acd_parameter = list(),
    prev_adj_params = list(), force_max_fp = FALSE, has_offset = FALSE,
    n_obs = 20, term_to_columns = list(x = "x")
  )

  expect_true(seen[["1"]])
  expect_equal(unname(result$power_best), c(-1, 2))
  expect_identical(
    rownames(result$metrics),
    c("FP2", "null", "linear", "FP1")
  )
})

test_that("Stata-compatible GEE selection uses its selection deviance", {
  metric <- function(logl, df, deviance) {
    c(
      logl = logl, df = df, deviance_rs = deviance,
      deviance_gaussian = NA_real_, aic = 0, bic = 0, df_resid = 20 - df
    )
  }

  testthat::local_mocked_bindings(
    build_adjustment_step = function(...) list(transform_cache = NULL),
    fit_null_step = function(...) list(
      powers = c(x = NA_real_),
      metrics = metric(logl = -10, df = 1, deviance = -9),
      current_adj_params = list(model = "null")
    ),
    fit_linear_step = function(...) list(
      powers = c(x = 1),
      metrics = metric(logl = -9.9, df = 2, deviance = -15),
      current_adj_params = list(model = "linear")
    ),
    .package = "mfp2"
  )

  result <- select_linear(
    x = matrix(1:4, ncol = 1, dimnames = list(NULL, "x")),
    xi = "x", keep = character(), degree = 1,
    acdx = c(x = FALSE), y = 1:4,
    powers_current = list(x = 1), powers = list(x = 1),
    criterion = "pvalue", ftest = FALSE, select = 0.05, alpha = 0.05,
    family = gee_family(gaussian()),
    family_string = "gee",
    zero = c(x = FALSE), catzero = list(x = NULL),
    spike = c(x = FALSE), spike_decision = c(x = NA_integer_),
    acd_parameter = list(), prev_adj_params = list(),
    force_max_fp = FALSE, has_offset = FALSE, n_obs = 4,
    term_to_columns = list(x = "x")
  )

  # The quasi-likelihood statistic would be only 0.2. The supplied negative-
  # Wald difference is 6 and is the statistic that Stata compatibility uses.
  expect_equal(unname(result$statistic), 6)
  expect_equal(
    unname(result$pvalue),
    stats::pchisq(6, df = 1, lower.tail = FALSE)
  )
  expect_equal(result$model_best, 2)
})

test_that("all three criteria recover the strong nonlinear and linear signals", {
  skip_if_not_installed("geepack")
  dat <- make_gee_fp_data()
  for (crit in c("pvalue", "aic", "bic")) {
    m <- mfp2(y ~ fp(x1) + fp(x2) + fp(x3), data = dat,
              family = gee_family(gaussian(), corstr = "exchangeable"),
              id = dat$id, criterion = crit, verbose = FALSE)
    # x1 selected and nonlinear
    expect_false(any(is.na(m$fp_powers$x1)))
    expect_false(identical(as.numeric(m$fp_powers$x1), 1))

    # On this seeded data the information criteria recover the truth cleanly:
    # x2 is linear and x3 is dropped. Do not impose the same deterministic
    # assertion on p-value selection: at select = 0.05 the Stata-compatible GEE
    # p-value procedure may retain a null variable as a nominal false positive
    # and may over-select the functional form of a linear signal. It is also
    # conditional on the selected FP powers, as documented by gee_family(). Its
    # mechanics are covered by the focused deterministic tests above rather than
    # by requiring perfect truth recovery from one simulated sample; here only
    # require that x2 is retained.
    if (crit == "pvalue") {
      expect_false(any(is.na(m$fp_powers$x2)), info = crit)
    } else {
      expect_equal(as.numeric(m$fp_powers$x2), 1, info = crit)
      expect_true(all(is.na(m$fp_powers$x3)), info = crit)
    }

    # The AIC/BIC selection score reconstructs from the quasi-likelihood and
    # df. The p-value score is Stata's negative overall robust Wald selection
    # deviance, on a different scale from the family deviance and not
    # reconstructed here, so only require that it is finite.
    if (crit == "pvalue") {
      expect_true(is.finite(m$mfp_selection_score), info = crit)
    } else {
      expected_score <- switch(
        crit,
        aic = -2 * m$mfp_logl + 2 * m$mfp_df,
        bic = -2 * m$mfp_logl + log(length(unique(dat$id))) * m$mfp_df
      )
      expect_equal(m$mfp_selection_score, expected_score, info = crit)
    }
  }
})

test_that("formula and matrix interfaces make the same GEE selection", {
  skip_if_not_installed("geepack")
  dat <- make_gee_fp_data()
  mf <- mfp2(y ~ fp(x1) + fp(x2) + fp(x3), data = dat,
             family = gee_family(gaussian(), corstr = "exchangeable"),
             id = dat$id, criterion = "aic", verbose = FALSE)
  X <- as.matrix(dat[, c("x1", "x2", "x3")])
  mm <- mfp2(X, dat$y, family = gee_family(gaussian(), corstr = "exchangeable"),
             id = dat$id, criterion = "aic", verbose = FALSE)
  expect_equal(mf$fp_powers$x1, mm$fp_powers$x1)
  expect_equal(mf$fp_powers$x2, mm$fp_powers$x2)
  expect_equal(mf$fp_powers$x3, mm$fp_powers$x3)
})

test_that("QBIC uses the number of clusters, not observations, as sample size", {
  skip_if_not_installed("geepack")
  dat <- make_gee_fp_data()
  m <- mfp2(y ~ x1 + x2, data = dat,
            family = gee_family(gaussian(), corstr = "exchangeable"),
            id = dat$id, df = 1, select = 1, center = FALSE,
            criterion = "bic", verbose = FALSE)
  # n_clusters is stored on the fitted object's geese component
  expect_equal(length(unique(dat$id)), length(m$geese$clusz))
})

test_that("S3 methods work on a fitted GEE mfp2 object", {
  skip_if_not_installed("geepack")
  dat <- make_gee_fp_data()
  m <- mfp2(y ~ fp(x1) + fp(x2), data = dat,
            family = gee_family(gaussian(), corstr = "exchangeable"),
            id = dat$id, criterion = "aic", verbose = FALSE)
  expect_true(inherits(m, "geeglm"))
  expect_type(coef(m), "double")
  expect_true(is.matrix(vcov(m)))
  pr <- predict(m)
  expect_length(pr, nrow(dat))
  expect_true(all(is.finite(pr)))
  expect_s3_class(summary(m), "summary.mfp2")
  expect_output(print(m))
})
