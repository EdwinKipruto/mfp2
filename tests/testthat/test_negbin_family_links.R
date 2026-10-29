# Regression tests for negbin_family() and its log, square-root, and identity
# links. Expected values are computed directly from fastglm::fastglm_nb() fits
# on independently constructed design matrices, or from closed-form inverse
# links and their derivatives, never from the mfp2 helpers under test.

skip_if_no_fastglm_nb <- function() {
  skip_if_not_installed("fastglm", minimum_version = "0.1.1")
  skip_if_not(
    "fastglm_nb" %in% getNamespaceExports("fastglm"),
    "fastglm does not export fastglm_nb()"
  )
  skip_if_not(
    "link" %in% names(formals(fastglm::fastglm_nb)),
    "fastglm_nb() does not accept a link"
  )
}

negbin_links <- c("log", "sqrt", "identity")

# Closed-form inverse links and their derivatives d mu / d eta.
nb_linkinv <- list(
  log = function(eta) exp(eta),
  sqrt = function(eta) eta^2,
  identity = function(eta) eta
)
nb_mu_eta <- list(
  log = function(eta) exp(eta),
  sqrt = function(eta) 2 * eta,
  identity = function(eta) rep(1, length(eta))
)

# Linear predictors are kept well inside each link's feasible domain: positive
# covariates, positive coefficients on the dominant term, and small positive
# offsets give a positive linear predictor (and mean) for every observation.
make_nb_link_data <- function(link, n = 320L, seed = 7301L) {
  set.seed(seed)
  dat <- data.frame(
    x1 = stats::runif(n, 1, 4),
    x2 = stats::runif(n, 0, 1),
    w = sample(c(1, 2), n, replace = TRUE),
    o = stats::runif(n, 0, 0.3)
  )
  eta <- switch(
    link,
    log = 0.30 + 0.35 * dat$x1 - 0.40 * dat$x2,
    sqrt = 0.80 + 0.45 * dat$x1 - 0.30 * dat$x2,
    identity = 1.00 + 1.20 * dat$x1 - 0.80 * dat$x2
  )
  dat$y <- stats::rnbinom(n, mu = nb_linkinv[[link]](eta + dat$o), size = 4)
  dat
}

# Fit a forced-linear model so that the comparison with fastglm_nb() isolates
# fitting from FP selection.
fit_nb_fixed_linear <- function(dat, link, center = FALSE) {
  mfp2(
    as.matrix(dat[, c("x1", "x2")]),
    dat$y,
    family = negbin_family(link = link),
    weights = dat$w,
    offset = dat$o,
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = center,
    cycles = 1,
    xorder = "original",
    verbose = FALSE
  )
}

# Direct fastglm_nb() fit with the controls mfp2 passes to the backend.
fit_nb_direct <- function(X, dat, link) {
  control <- stats::glm.control()
  fastglm::fastglm_nb(
    x = cbind("(Intercept)" = 1, X),
    y = dat$y,
    weights = dat$w,
    offset = dat$o,
    link = link,
    method = 0L,
    tol = control$epsilon,
    maxit = control$maxit
  )
}

first_ggplot <- function(x) {
  if (inherits(x, "ggplot")) return(x)
  if (is.list(x)) {
    for (element in x) {
      found <- first_ggplot(element)
      if (!is.null(found)) return(found)
    }
  }
  NULL
}


# ---------------------------------------------------------------------------
# Constructor and family normalisation (no fastglm needed)
# ---------------------------------------------------------------------------

test_that("negbin_family() defaults to the log link and accepts all three links", {
  default <- negbin_family()
  expect_s3_class(default, c("mfp2_negbin_family", "mfp2_family"), exact = TRUE)
  expect_identical(default$family, "negbin")
  expect_identical(default$link, "log")

  for (link in negbin_links) {
    expect_identical(negbin_family(link = link)$link, link)
    expect_identical(negbin_family(link)$link, link)
  }
})

test_that("negbin_family() rejects unsupported link specifications", {
  expect_error(negbin_family("logit"), "Unsupported link.*'logit'")
  expect_error(negbin_family("inverse"), "Supported links are")
  expect_error(negbin_family("Log"), "Unsupported link")
  expect_error(negbin_family(c("log", "sqrt")), "single character string")
  expect_error(negbin_family(NA_character_), "single character string")
  expect_error(negbin_family(""), "single character string")
  expect_error(negbin_family(1), "single character string")
  expect_error(negbin_family(log), "by name, not as a function")
  expect_error(negbin_family(stats::make.link("log")), "by name")
})

test_that("family = 'negbin' normalises to negbin_family(link = 'log')", {
  shorthand <- mfp2:::normalize_family_argument("negbin")
  explicit <- mfp2:::normalize_family_argument(negbin_family(link = "log"))
  expect_identical(shorthand$family_string, "negbin")
  expect_identical(explicit$family_string, "negbin")
  expect_identical(shorthand$family, negbin_family(link = "log"))
  expect_identical(explicit$family, negbin_family(link = "log"))

  sqrt_info <- mfp2:::normalize_family_argument(negbin_family("sqrt"))
  expect_identical(sqrt_info$family$link, "sqrt")
})

test_that("a modified negbin_family() object with an invalid link is rejected", {
  dat <- make_nb_link_data("log", n = 60L)
  bad <- negbin_family()
  bad$link <- "logit"
  expect_error(
    mfp2(y ~ x1, data = dat, family = bad, verbose = FALSE),
    "Unsupported link.*'logit'"
  )
  expect_error(
    mfpi(
      y ~ x1 + x2,
      data = transform(dat, g = factor(rep(c("a", "b"), length.out = 60L))),
      family = bad,
      group_var = "g",
      interaction_vars = "x1",
      verbose = FALSE
    ),
    "Unsupported link"
  )
})

test_that("GEE does not accept negative-binomial families", {
  expect_error(gee_family(negbin_family()), "must be one of")
  expect_error(gee_family(negbin_family("sqrt")), "must be one of")
})


# ---------------------------------------------------------------------------
# Backward compatibility and interface equivalence
# ---------------------------------------------------------------------------

test_that("'negbin' and negbin_family('log') select and fit the same model", {
  skip_if_no_fastglm_nb()

  dat <- make_negbin_test_data(n = 240L, seed = 7302L)
  # Keep x2 linear: plain numeric formula terms are also FP-tested by default.
  fit_shorthand <- mfp2(
    y ~ fp(x1) + fp(x2, df = 1) + x3, data = dat, family = "negbin",
    verbose = FALSE
  )
  fit_explicit <- mfp2(
    y ~ fp(x1) + fp(x2, df = 1) + x3, data = dat,
    family = negbin_family(link = "log"), verbose = FALSE
  )

  expect_identical(fit_explicit$family$link, "log")
  expect_identical(fit_shorthand$family$link, "log")
  expect_identical(fit_explicit$fp_terms, fit_shorthand$fp_terms)
  expect_equal(coef(fit_explicit), coef(fit_shorthand), tolerance = 1e-12)
  expect_equal(fit_explicit$theta, fit_shorthand$theta, tolerance = 1e-12)
  nd <- dat[1:15, ]
  for (type in c("link", "response")) {
    expect_equal(
      predict(fit_explicit, newdata = nd, type = type),
      predict(fit_shorthand, newdata = nd, type = type),
      tolerance = 1e-12
    )
  }
})

test_that("formula and matrix interfaces agree for every negative-binomial link", {
  skip_if_no_fastglm_nb()

  for (link in negbin_links) {
    dat <- make_nb_link_data(link, seed = 7303L)
    X <- as.matrix(dat[, c("x1", "x2")])

    fit_formula <- mfp2(
      y ~ fp(x1, df = 2) + fp(x2, df = 2) + offset(o),
      data = dat,
      weights = w,
      family = negbin_family(link = link),
      verbose = FALSE
    )
    fit_matrix <- mfp2(
      X, dat$y,
      weights = dat$w,
      offset = dat$o,
      df = 2,
      family = negbin_family(link = link),
      verbose = FALSE
    )

    expect_identical(fit_formula$family$link, link, info = link)
    expect_identical(fit_matrix$family$link, link, info = link)
    # expect_equal(): the interfaces store df settings as double vs integer.
    expect_equal(fit_formula$fp_terms, fit_matrix$fp_terms, info = link)
    expect_equal(unname(coef(fit_formula)), unname(coef(fit_matrix)),
                 tolerance = 1e-10, info = link)
    expect_equal(fit_formula$theta, fit_matrix$theta, tolerance = 1e-10,
                 info = link)
    expect_equal(
      unname(predict(fit_formula, type = "response")),
      unname(predict(fit_matrix, type = "response")),
      tolerance = 1e-10,
      info = link
    )
  }
})


# ---------------------------------------------------------------------------
# Agreement with fastglm::fastglm_nb()
# ---------------------------------------------------------------------------

test_that("fixed linear fits match fastglm_nb() for every link", {
  skip_if_no_fastglm_nb()

  for (link in negbin_links) {
    dat <- make_nb_link_data(link, seed = 7304L)
    X <- as.matrix(dat[, c("x1", "x2")])
    fit <- fit_nb_fixed_linear(dat, link, center = FALSE)
    direct <- fit_nb_direct(X, dat, link)

    expect_identical(fit$family$link, link, info = link)
    expect_identical(direct$family$link, link, info = link)
    expect_equal(unname(coef(fit)), unname(direct$coefficients),
                 tolerance = 1e-7, info = link)
    expect_equal(fit$theta, direct$theta, tolerance = 1e-6, info = link)
    expect_equal(fit$mfp_logl, direct$twologlik / 2, tolerance = 1e-8,
                 info = link)
    expect_equal(as.numeric(fit$fitted.values),
                 as.numeric(direct$fitted.values),
                 tolerance = 1e-7, info = link)
    expect_equal(unname(as.matrix(stats::vcov(fit))),
                 unname(as.matrix(direct$cov.unscaled)),
                 tolerance = 1e-6, info = link)
  }
})

test_that("centred fits match fastglm_nb() on a design centred at the column means", {
  skip_if_no_fastglm_nb()

  for (link in negbin_links) {
    dat <- make_nb_link_data(link, seed = 7305L)
    X <- as.matrix(dat[, c("x1", "x2")])
    fit <- fit_nb_fixed_linear(dat, link, center = TRUE)

    # Linear terms are centred at their unweighted means.
    expected_centers <- colMeans(X)
    expect_equal(unname(fit$centers), unname(expected_centers),
                 tolerance = 1e-12, info = link)

    direct <- fit_nb_direct(sweep(X, 2L, expected_centers), dat, link)
    expect_equal(unname(coef(fit)), unname(direct$coefficients),
                 tolerance = 1e-7, info = link)
    expect_equal(as.numeric(fit$fitted.values),
                 as.numeric(direct$fitted.values),
                 tolerance = 1e-7, info = link)
  }
})


# ---------------------------------------------------------------------------
# Prediction, standard errors, and confidence intervals
# ---------------------------------------------------------------------------

test_that("link and response predictions use the selected link with and without newdata", {
  skip_if_no_fastglm_nb()

  for (link in negbin_links) {
    dat <- make_nb_link_data(link, seed = 7306L)
    fit <- mfp2(
      y ~ fp(x1, df = 1, select = 1, shift = 0, scale = 1, center = FALSE) +
        fp(x2, df = 1, select = 1, shift = 0, scale = 1, center = FALSE) +
        offset(o),
      data = dat,
      weights = w,
      family = negbin_family(link = link),
      verbose = FALSE
    )
    direct <- fit_nb_direct(as.matrix(dat[, c("x1", "x2")]), dat, link)
    beta <- unname(direct$coefficients)

    # Training data.
    eta_train <- drop(cbind(1, dat$x1, dat$x2) %*% beta) + dat$o
    expect_equal(unname(predict(fit, type = "link")), eta_train,
                 tolerance = 1e-7, info = link)
    expect_equal(unname(predict(fit, type = "response")),
                 nb_linkinv[[link]](eta_train),
                 tolerance = 1e-7, info = link)

    # New data, including a new offset.
    nd <- data.frame(x1 = c(1.2, 2.5, 3.8), x2 = c(0.1, 0.5, 0.9),
                     o = c(0, 0.1, 0.25))
    eta_new <- drop(cbind(1, nd$x1, nd$x2) %*% beta) + nd$o
    expect_equal(unname(predict(fit, newdata = nd, type = "link")), eta_new,
                 tolerance = 1e-7, info = link)
    expect_equal(unname(predict(fit, newdata = nd, type = "response")),
                 nb_linkinv[[link]](eta_new),
                 tolerance = 1e-7, info = link)
  }
})

test_that("response-scale standard errors follow the delta method for each link", {
  skip_if_no_fastglm_nb()

  for (link in negbin_links) {
    dat <- make_nb_link_data(link, seed = 7307L)
    fit <- fit_nb_fixed_linear(dat, link, center = FALSE)
    direct <- fit_nb_direct(as.matrix(dat[, c("x1", "x2")]), dat, link)

    newx <- cbind(x1 = c(1.1, 2.0, 3.9), x2 = c(0.9, 0.4, 0.05))
    newoffset <- c(0.05, 0, 0.2)
    Xn <- cbind(1, newx)
    V <- as.matrix(direct$cov.unscaled)
    eta <- drop(Xn %*% direct$coefficients) + newoffset
    se_link <- sqrt(rowSums((Xn %*% V) * Xn))
    se_response <- se_link * abs(nb_mu_eta[[link]](eta))

    pred_link <- predict(fit, newdata = newx, newoffset = newoffset,
                         type = "link", se.fit = TRUE)
    pred_response <- predict(fit, newdata = newx, newoffset = newoffset,
                             type = "response", se.fit = TRUE)

    expect_equal(unname(pred_link$fit), eta, tolerance = 1e-7, info = link)
    expect_equal(unname(pred_link$se.fit), se_link, tolerance = 1e-6,
                 info = link)
    expect_equal(unname(pred_response$fit), nb_linkinv[[link]](eta),
                 tolerance = 1e-7, info = link)
    expect_equal(unname(pred_response$se.fit), se_response,
                 tolerance = 1e-6, info = link)
  }
})

test_that("confint() gives the documented Wald intervals for each link", {
  skip_if_no_fastglm_nb()

  for (link in negbin_links) {
    dat <- make_nb_link_data(link, seed = 7308L)
    fit <- fit_nb_fixed_linear(dat, link, center = FALSE)
    direct <- fit_nb_direct(as.matrix(dat[, c("x1", "x2")]), dat, link)
    se <- sqrt(diag(as.matrix(direct$cov.unscaled)))

    for (level in c(0.95, 0.90)) {
      z <- stats::qnorm(1 - (1 - level) / 2)
      expected <- cbind(direct$coefficients - z * se,
                        direct$coefficients + z * se)
      ci <- confint(fit, level = level)
      expect_equal(unname(ci), unname(expected), tolerance = 1e-6,
                   info = paste(link, level))
      expect_identical(rownames(ci), names(coef(fit)))
    }
  }
})


# ---------------------------------------------------------------------------
# FP selection, nonlinear fits, and plotting
# ---------------------------------------------------------------------------

make_nb_nonlinear_data <- function(n = 400L, seed = 7309L) {
  set.seed(seed)
  dat <- data.frame(x = stats::runif(n, 0.5, 4), z = stats::runif(n))
  # sqrt(mu) is a strongly nonlinear, positive function of x.
  dat$y <- stats::rnbinom(n, mu = (0.6 + 2.4 / dat$x + 0.2 * dat$z)^2,
                          size = 5)
  dat
}

test_that("the chosen link is used for every backend fit during FP selection", {
  skip_if_no_fastglm_nb()

  dat <- make_nb_nonlinear_data()
  real_fastglm_nb <- fastglm::fastglm_nb
  links_used <- character()
  local_mocked_bindings(
    fastglm_nb = function(..., link) {
      links_used <<- c(links_used, link)
      real_fastglm_nb(..., link = link)
    },
    .package = "fastglm"
  )

  fit <- mfp2(y ~ fp(x) + z, data = dat, family = negbin_family("sqrt"),
              verbose = FALSE)

  # Null, linear, all FP1 and FP2 candidates, and the final fit.
  expect_gt(length(links_used), 40L)
  expect_true(all(links_used == "sqrt"))
  expect_identical(fit$family$link, "sqrt")
  expect_identical(fit$family_string, "negbin")
})

test_that("a nonlinear square-root-link fit predicts and plots on the right scales", {
  skip_if_no_fastglm_nb()

  dat <- make_nb_nonlinear_data()
  fit <- mfp2(y ~ fp(x) + z, data = dat, family = negbin_family("sqrt"),
              verbose = FALSE)

  power_columns <- grep("^power", names(fit$fp_terms), value = TRUE)
  powers <- unlist(fit$fp_terms["x", power_columns])
  powers <- powers[!is.na(powers)]
  expect_true(length(powers) >= 1L && !identical(unname(powers), 1))

  nd <- data.frame(x = c(0.6, 1, 2, 3.5), z = 0.5)
  eta <- predict(fit, newdata = nd, type = "link", se.fit = TRUE)
  mu <- predict(fit, newdata = nd, type = "response", se.fit = TRUE)
  expect_true(all(is.finite(eta$fit)) && all(eta$fit > 0))
  expect_equal(unname(mu$fit), unname(eta$fit)^2, tolerance = 1e-10)
  expect_equal(unname(mu$se.fit), unname(eta$se.fit) * 2 * unname(eta$fit),
               tolerance = 1e-10)
  # The fitted mean is decreasing in x, as in the generating model.
  expect_true(all(diff(mu$fit[order(nd$x)]) < 0))

  partial <- first_ggplot(plot(fit, terms = "x", partial_only = TRUE))
  expect_s3_class(partial, "ggplot")
  expect_identical(partial$labels$y, "Contribution to sqrt(mean)")

  for (link in c("log", "identity")) {
    dat_link <- make_nb_link_data(link, seed = 7310L)
    fit_link <- fit_nb_fixed_linear(dat_link, link)
    plotted <- first_ggplot(plot(fit_link, terms = "x1", partial_only = TRUE))
    expect_identical(
      plotted$labels$y,
      if (identical(link, "log")) "Contribution to log mean" else
        "Contribution to mean response",
      info = link
    )
  }

  summary_text <- paste(capture.output(print(summary(fit))), collapse = "\n")
  expect_match(summary_text, "sqrt link", fixed = TRUE)
})


# ---------------------------------------------------------------------------
# MFPI propagation
# ---------------------------------------------------------------------------

test_that("mfpi() propagates the link to functions and between-group contrasts", {
  skip_if_no_fastglm_nb()

  set.seed(7311L)
  n <- 300L
  dat <- data.frame(
    group = factor(rep(c("control", "treated"), each = n / 2L)),
    x = stats::runif(n, 1, 4),
    z = stats::runif(n)
  )
  treated <- as.numeric(dat$group == "treated")
  eta_true <- 0.9 + 0.3 * dat$x + 0.2 * treated + 0.25 * dat$x * treated -
    0.2 * dat$z
  dat$y <- stats::rnbinom(n, mu = eta_true^2, size = 4)

  fit <- mfpi(
    dat[, c("group", "x", "z")], dat$y,
    family = negbin_family("sqrt"),
    group_var = "group",
    interaction_vars = "x",
    interaction_forms = c(x = "linear"),
    flex = "flex1",
    cycles = 1,
    df = 1,
    select = 1,
    alpha = 1,
    p_interact = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = FALSE
  )
  expect_identical(fit$family_string, "negbin")
  stored <- fit$var_winners[["x"]]$fit$test_results$interaction_model$fit
  expect_identical(stored$family$link, "sqrt")

  # Independent interaction model: separate x slopes in each group.
  X <- cbind(
    treated = treated,
    x_control = dat$x * (1 - treated),
    x_treated = dat$x * treated,
    z = dat$z
  )
  direct <- fastglm::fastglm_nb(
    cbind("(Intercept)" = 1, X), dat$y, link = "sqrt", method = 0L,
    tol = stats::glm.control()$epsilon, maxit = stats::glm.control()$maxit
  )
  b <- unname(direct$coefficients)
  V <- as.matrix(direct$cov.unscaled)

  x_new <- c(1.5, 2.5, 3.5)
  pred <- predict(
    fit,
    newdata = data.frame(group = "control", x = x_new, z = 0.5),
    terms = "x",
    type = "both",
    se.fit = TRUE
  )

  # Between-group contrast on the square-root-mean scale and its Wald SE.
  A <- cbind(0, 1, -x_new, x_new, 0)
  expected_diff <- drop(A %*% b)
  expected_diff_se <- sqrt(rowSums((A %*% V) * A))
  differences <- pred$differences[order(pred$differences$x), ]
  expect_equal(differences$fit, expected_diff, tolerance = 1e-6)
  expect_equal(differences$se.fit, expected_diff_se, tolerance = 1e-5)

  # Group functions differ by exactly the reported contrast.
  functions <- pred$functions
  f_control <- functions$fit[functions$group == "control"][order(x_new)]
  f_treated <- functions$fit[functions$group == "treated"][order(x_new)]
  expect_equal(f_treated - f_control, expected_diff, tolerance = 1e-6)

  # Individual response predictions square the linear predictor.
  nd <- dat[1:10, c("group", "x", "z")]
  link_pred <- predict(fit, newdata = nd, terms = "x", type = "link",
                       se.fit = FALSE)
  response_pred <- predict(fit, newdata = nd, terms = "x", type = "response",
                           se.fit = FALSE)
  eta_direct <- drop(cbind(1, X[1:10, , drop = FALSE]) %*% b)
  expect_equal(link_pred$predictions$fit, eta_direct, tolerance = 1e-6)
  expect_equal(response_pred$predictions$fit, eta_direct^2, tolerance = 1e-6)

  difference_plot <- first_ggplot(
    plot(fit, terms = "x", plot_type = "difference")
  )
  expect_s3_class(difference_plot, "ggplot")
  expect_match(difference_plot$labels$y, "^Link-scale difference")
})


# ---------------------------------------------------------------------------
# Unsuccessful fits
# ---------------------------------------------------------------------------

test_that("infeasible identity-link FP candidates are skipped without changing the link", {
  skip_if_no_fastglm_nb()

  dat <- make_nb_link_data("identity", n = 260L, seed = 7312L)
  real_fastglm_nb <- fastglm::fastglm_nb
  links_used <- character()
  corrupted <- 0L
  local_mocked_bindings(
    fastglm_nb = function(x, ..., link) {
      links_used <<- c(links_used, link)
      fit <- real_fastglm_nb(x, ..., link = link)
      # FP2 candidates for x1 have an intercept and two FP columns. Make those
      # whose first FP column decreases in x1 (negative first power) report an
      # infeasible mean, as an identity-link fit outside its domain would.
      if (ncol(x) == 3L && stats::cor(x[, 2L], dat$x1) < 0) {
        corrupted <<- corrupted + 1L
        fit$fitted.values[1L] <- -0.5
      }
      fit
    },
    .package = "fastglm"
  )

  expect_warning(
    fit <- mfp2(
      dat[, "x1", drop = FALSE] |> as.matrix(), dat$y,
      family = negbin_family("identity"),
      shift = 0, scale = 1,
      verbose = FALSE
    ),
    "could not be fitted and were skipped.*valid domain.*identity link"
  )
  expect_gt(corrupted, 0L)
  expect_true(all(links_used == "identity"))
  expect_identical(fit$family$link, "identity")
  expect_true(all(fit$fitted.values > 0))
})

test_that("an infeasible required fit stops with an informative error", {
  skip_if_no_fastglm_nb()

  dat <- make_nb_link_data("sqrt", n = 120L, seed = 7313L)
  real_fastglm_nb <- fastglm::fastglm_nb
  local_mocked_bindings(
    fastglm_nb = function(...) {
      fit <- real_fastglm_nb(...)
      fit$linear.predictors[1L] <- -0.1
      fit$fitted.values[1L] <- 0.01
      fit
    },
    .package = "fastglm"
  )

  expect_error(
    mfp2(y ~ fp(x1) + x2, data = dat, family = negbin_family("sqrt"),
         verbose = FALSE),
    "sqrt link requires a positive linear predictor"
  )
})

test_that("a backend that ignores the requested link is rejected", {
  skip_if_no_fastglm_nb()

  dat <- make_nb_link_data("sqrt", n = 120L, seed = 7314L)
  real_fastglm_nb <- fastglm::fastglm_nb
  local_mocked_bindings(
    fastglm_nb = function(..., link) real_fastglm_nb(..., link = "log"),
    .package = "fastglm"
  )

  expect_error(
    mfp2(y ~ x1 + x2, data = dat, family = negbin_family("sqrt"),
         verbose = FALSE),
    "link 'log' instead of the requested 'sqrt' link"
  )
})

test_that("backend errors for restricted links are reported with the link", {
  skip_if_no_fastglm_nb()

  dat <- make_nb_link_data("identity", n = 120L, seed = 7315L)
  local_mocked_bindings(
    fastglm_nb = function(...) stop("no valid set of coefficients"),
    .package = "fastglm"
  )

  expect_error(
    mfp2(y ~ x1 + x2, data = dat, family = negbin_family("identity"),
         verbose = FALSE),
    "\\(identity link\\): no valid set of coefficients"
  )
})
