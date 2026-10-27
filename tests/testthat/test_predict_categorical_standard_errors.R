# Fixed covariance fixtures make the numerical checks independent of optimizer
# tolerances. Their nonzero off-diagonal entries detect omitted covariance terms.
categorical_ordinal_fixture <- function(link = "logistic") {
  X <- matrix(c(-0.5, 0.2, 0, 0.5, 0.8, -0.3), ncol = 2L, byrow = TRUE,
              dimnames = list(c("third", "first", "seventh"), c("x1", "x2")))
  intercepts <- c("y>=5" = 0.7, "y>=9" = -0.8)
  cf <- c(intercepts, x1 = 0.4, x2 = -0.2)
  A <- matrix(seq_len(16L) / 20, 4L)
  V <- crossprod(A) + diag(0.2, 4L)
  dimnames(V) <- list(names(cf), names(cf))
  list(
    family_string = "ordinal", coefficients = cf,
    mfp2_ordinal_intercepts = intercepts, mfp2_ordinal_link = link,
    mfp2_ordinal_levels = c("2", "5", "9"), interceptRef = 2L,
    x = X, nobs = nrow(X),
    mfp2_family = list(prepared = list(has_offset = TRUE, offset = c(0.1, -0.2, 0.3))),
    transformed_to_model_columns = c(x1 = "x1", x2 = "x2"),
    prediction_test_covariance = V
  )
}

categorical_multinomial_fixture <- function() {
  X <- cbind("(Intercept)" = 1, x = c(-0.5, 0.2, 1))
  rownames(X) <- c("third", "first", "seventh")
  B <- matrix(c(0.2, 0.4, -0.3, 0.1), 2L, byrow = TRUE,
              dimnames = list(c("A", "C"), colnames(X)))
  cf_names <- unlist(lapply(rownames(B), function(q) paste0(q, ":", colnames(B))))
  A <- matrix(seq_len(16L) / 25, 4L)
  V <- crossprod(A) + diag(0.15, 4L)
  dimnames(V) <- list(cf_names, cf_names)
  offsets <- cbind(B = 0, A = c(0.1, -0.2, 0.3), C = c(-0.1, 0.2, 0.1))
  list(
    mfp2_coefficient_matrix = B, mfp2_design = X,
    mfp2_offset_matrix = offsets, class_levels = c("B", "A", "C"),
    has_offset = TRUE,
    mfp2_family = list(prepared = list(
      levels = c("B", "A", "C"), original_levels = c("A", "B", "C"),
      nonreference = c("A", "C"), n_classes = 3L, n_logits = 2L
    )),
    transformed_to_model_columns = c(x = "x"),
    prediction_test_covariance = V
  )
}

test_that("ordinal link SEs include exactly the requested cutoff and slopes", {
  object <- categorical_ordinal_fixture()
  # These internal helpers deliberately use stats::vcov(), as native fits do.
  testthat::local_mocked_bindings(vcov = function(object, intercepts = "all", ...) {
    V <- object$prediction_test_covariance
    if (intercepts == "none") V <- V[c("x1", "x2"), c("x1", "x2")]
    V
  }, .package = "stats")
  eta <- drop(object$x %*% object$coefficients[c("x1", "x2")]) +
    object$mfp2_family$prepared$offset
  for (add in c(TRUE, FALSE)) {
    result <- mfp2:::mfp2_predict_ordinal(object, se.fit = TRUE, add_intercept = add)
    X <- if (add) cbind(1, object$x) else object$x
    columns <- if (add) c("y>=9", "x1", "x2") else c("x1", "x2")
    V <- object$prediction_test_covariance[columns, columns]
    expected <- eta + if (add) -0.8 else 0
    expect_equal(result$fit, expected)
    expect_equal(unname(result$se.fit), unname(sqrt(diag(X %*% V %*% t(X)))))
    expect_identical(names(result$fit), names(result$se.fit))
    # A single new observation preserves vector shape, names, and the offset.
    one <- mfp2:::mfp2_predict_ordinal(
      object, transformed = object$x[2L, , drop = FALSE],
      newoffset = -0.2, se.fit = TRUE, add_intercept = add
    )
    expect_equal(one$fit, result$fit[2L])
    expect_equal(one$se.fit, result$se.fit[2L])
  }
})

test_that("ordinal probability and mean SEs agree with numerical derivatives for all links", {
  testthat::local_mocked_bindings(vcov = function(object, ...) {
    # Intentionally permute both axes to check coefficient-name alignment.
    V <- object$prediction_test_covariance
    V[rev(rownames(V)), rev(colnames(V)), drop = FALSE]
  }, .package = "stats")
  cdfs <- list(logistic = stats::plogis, probit = stats::pnorm,
               loglog = function(z) exp(-exp(-z)),
               cloglog = function(z) 1 - exp(-exp(z)), cauchit = stats::pcauchy)
  for (link in names(cdfs)) {
    object <- categorical_ordinal_fixture(link)
    V <- object$prediction_test_covariance
    result <- mfp2:::mfp2_predict_ordinal(object, type = "response", se.fit = TRUE)
    mean <- mfp2:::mfp2_predict_ordinal(object, type = "mean", se.fit = TRUE)
    for (i in seq_len(nrow(object$x))) {
      probabilities <- function(theta) {
        eta <- sum(object$x[i, ] * theta[3:4]) +
          object$mfp2_family$prepared$offset[i]
        -diff(c(1, cdfs[[link]](theta[1:2] + eta), 0))
      }
      G <- numDeriv::jacobian(probabilities, object$coefficients)
      expected_se <- sqrt(diag(G %*% V %*% t(G)))
      expect_equal(unname(result$fit[i, ]), unname(probabilities(object$coefficients)))
      expect_equal(unname(result$se.fit[i, ]), unname(expected_se), tolerance = 1e-6)
      mean_gradient <- numDeriv::grad(function(theta) {
        sum(c(2, 5, 9) * probabilities(theta))
      }, object$coefficients)
      expect_equal(unname(mean$se.fit[i]),
                   sqrt(drop(t(mean_gradient) %*% V %*% mean_gradient)),
                   tolerance = 1e-6)
    }
    expect_identical(dimnames(result$fit), dimnames(result$se.fit))
    expect_equal(unname(rowSums(result$fit)), rep(1, 3L))
    expect_equal(unname(mean$fit), unname(drop(result$fit %*% c(2, 5, 9))))
    for (type in c("response", "mean")) {
      expect_equal(
        mfp2:::mfp2_predict_ordinal(object, type = type, se.fit = TRUE,
                                   add_intercept = FALSE),
        mfp2:::mfp2_predict_ordinal(object, type = type, se.fit = TRUE)
      )
    }
    # A reference category has a genuine probability SE, even though its
    # corresponding cumulative endpoint is fixed at one.
    expect_true(all(result$se.fit[, 1L] > 0))
    one <- mfp2:::mfp2_predict_ordinal(
      object, transformed = object$x[1L, , drop = FALSE], newoffset = 0.1,
      type = "response", se.fit = TRUE
    )
    expect_equal(one$fit, result$fit[1L, , drop = FALSE])
    expect_equal(one$se.fit, result$se.fit[1L, , drop = FALSE])
  }
})

test_that("intercept-only ordinal models and category scores have valid SEs", {
  object <- categorical_ordinal_fixture()
  object$coefficients <- object$mfp2_ordinal_intercepts
  object$x <- NULL
  object$prediction_test_covariance <- object$prediction_test_covariance[1:2, 1:2]
  object$mfp2_ordinal_levels <- c("low", "middle", "high")
  testthat::local_mocked_bindings(vcov = function(object, ...) {
    object$prediction_test_covariance
  }, .package = "stats")
  bare <- mfp2:::mfp2_predict_ordinal(object, se.fit = TRUE, add_intercept = FALSE)
  expect_equal(bare$fit, object$mfp2_family$prepared$offset)
  expect_equal(bare$se.fit, rep(0, object$nobs))
  probabilities <- mfp2:::mfp2_predict_ordinal(object, type = "response", se.fit = TRUE)
  mean <- mfp2:::mfp2_predict_ordinal(object, type = "mean", se.fit = TRUE)
  expect_equal(dim(probabilities$se.fit), c(3L, 3L))
  expect_true(all(is.finite(probabilities$se.fit)))
  expect_equal(unname(mean$fit), unname(drop(probabilities$fit %*% 1:3)))
  expect_true(all(is.finite(mean$se.fit) & mean$se.fit > 0))
})

test_that("multinomial link and probability SEs use the full named covariance", {
  object <- categorical_multinomial_fixture()
  testthat::local_mocked_bindings(vcov = function(object, ...) {
    V <- object$prediction_test_covariance
    V[c(3, 1, 4, 2), c(2, 4, 1, 3), drop = FALSE]
  }, .package = "stats")
  V <- object$prediction_test_covariance
  theta <- as.vector(t(object$mfp2_coefficient_matrix))
  link <- mfp2:::mfp2_predict_multinomial(object, type = "link", se.fit = TRUE)
  response <- mfp2:::mfp2_predict_multinomial(object, type = "response", se.fit = TRUE)
  for (i in seq_len(nrow(object$mfp2_design))) {
    logits <- function(theta) {
      B <- matrix(theta, 2L, byrow = TRUE)
      drop(B %*% object$mfp2_design[i, ]) + object$mfp2_offset_matrix[i, -1L]
    }
    probabilities <- function(theta) {
      eta <- c(0, logits(theta))
      p <- exp(eta - max(eta))
      # Internal order B,A,C is restored to original response order A,B,C.
      (p / sum(p))[c(2, 1, 3)]
    }
    for (type in c("link", "response")) {
      f <- if (type == "link") logits else probabilities
      result <- if (type == "link") link else response
      G <- numDeriv::jacobian(f, theta)
      expect_equal(unname(result$fit[i, ]), unname(f(theta)))
      expect_equal(unname(result$se.fit[i, ]), unname(sqrt(diag(G %*% V %*% t(G)))),
                   tolerance = 1e-6)
    }
  }
  expect_identical(colnames(response$fit), c("A", "B", "C"))
  expect_identical(dimnames(response$fit), dimnames(response$se.fit))
  expect_true(all(response$se.fit[, "B"] > 0))
  # Newdata offsets can be supplied on the logit scale; names preserve order.
  for (type in c("link", "response")) {
    many <- if (type == "link") link else response
    one <- mfp2:::mfp2_predict_multinomial(
      object, transformed = object$mfp2_design[2L, "x", drop = FALSE],
      newoffset = object$mfp2_offset_matrix[2L, -1L, drop = FALSE],
      type = type, se.fit = TRUE
    )
    expect_equal(one$fit, many$fit[2L, , drop = FALSE])
    expect_equal(one$se.fit, many$se.fit[2L, , drop = FALSE])
  }
  expect_error(mfp2:::mfp2_predict_multinomial(object, type = "class", se.fit = TRUE),
               "class labels")
})

test_that("binary ordinal and intercept-only multinomial predictions preserve SE shapes", {
  testthat::local_mocked_bindings(vcov = function(object, ...) {
    object$prediction_test_covariance
  }, .package = "stats")
  ordinal <- categorical_ordinal_fixture()
  ordinal$mfp2_ordinal_intercepts <- ordinal$mfp2_ordinal_intercepts[1L]
  ordinal$coefficients <- ordinal$coefficients[-2L]
  ordinal$prediction_test_covariance <- ordinal$prediction_test_covariance[-2L, -2L]
  ordinal$mfp2_ordinal_levels <- c("0", "1")
  ordinal$interceptRef <- 1L
  response <- mfp2:::mfp2_predict_ordinal(ordinal, type = "response", se.fit = TRUE)
  mean <- mfp2:::mfp2_predict_ordinal(ordinal, type = "mean", se.fit = TRUE)
  expect_identical(dim(response$se.fit), c(3L, 2L))
  expect_equal(mean$fit, response$fit[, 2L])
  expect_equal(mean$se.fit, response$se.fit[, 2L])

  multinomial <- categorical_multinomial_fixture()
  multinomial$mfp2_design <- multinomial$mfp2_design[, 1L, drop = FALSE]
  multinomial$mfp2_coefficient_matrix <- multinomial$mfp2_coefficient_matrix[, 1L, drop = FALSE]
  multinomial$prediction_test_covariance <- multinomial$prediction_test_covariance[c(1, 3), c(1, 3)]
  for (type in c("link", "response")) {
    result <- mfp2:::mfp2_predict_multinomial(multinomial, type = type, se.fit = TRUE)
    expect_identical(dim(result$fit), dim(result$se.fit))
    expect_true(all(is.finite(result$se.fit) & result$se.fit > 0))
    one <- mfp2:::mfp2_predict_multinomial(
      multinomial, transformed = data.frame(row.names = "first"),
      newoffset = multinomial$mfp2_offset_matrix[2L, -1L, drop = FALSE],
      type = type, se.fit = TRUE
    )
    expect_equal(one$fit, result$fit[2L, , drop = FALSE])
    expect_equal(one$se.fit, result$se.fit[2L, , drop = FALSE])
  }
})

test_that("categorical SEs validate required covariance and negative variances", {
  ordinal <- categorical_ordinal_fixture()
  multinomial <- categorical_multinomial_fixture()
  testthat::local_mocked_bindings(vcov = function(object, intercepts = "all", ...) {
    V <- object$prediction_test_covariance
    if (intercepts == "none") V <- V[c("x1", "x2"), c("x1", "x2")]
    V
  }, .package = "stats")
  ordinal$prediction_test_covariance[1L, 1L] <- NA_real_
  expect_no_error(mfp2:::mfp2_predict_ordinal(ordinal, se.fit = TRUE))
  expect_no_error(mfp2:::mfp2_predict_ordinal(ordinal, se.fit = TRUE, add_intercept = FALSE))
  for (bad in c(NA_real_, NaN, Inf, -Inf)) {
    ordinal$prediction_test_covariance[1L, 1L] <- bad
    for (type in c("response", "mean")) {
      expect_error(mfp2:::mfp2_predict_ordinal(ordinal, type = type, se.fit = TRUE),
                   "non-finite entries")
    }
    multinomial$prediction_test_covariance[1L, 1L] <- bad
    for (type in c("link", "response")) {
      expect_error(mfp2:::mfp2_predict_multinomial(multinomial, type = type, se.fit = TRUE),
                   "non-finite entries")
    }
  }
  for (type in c("link", "response", "mean")) {
    ordinal$prediction_test_covariance <- -diag(4L)
    dimnames(ordinal$prediction_test_covariance) <- list(names(ordinal$coefficients),
                                                        names(ordinal$coefficients))
    expect_error(mfp2:::mfp2_predict_ordinal(ordinal, type = type, se.fit = TRUE),
                 "materially negative")
  }
  V <- multinomial$prediction_test_covariance
  V[,] <- 0
  diag(V) <- -1
  multinomial$prediction_test_covariance <- V
  for (type in c("link", "response")) {
    expect_error(mfp2:::mfp2_predict_multinomial(multinomial, type = type, se.fit = TRUE),
                 "materially negative")
  }
})

test_that("linear ordinal full-model predictions and SEs match rms", {
  skip_if_not_installed("rms", minimum_version = "7.0-0")
  set.seed(824)
  dat <- data.frame(x1 = runif(300L, 0.5, 3), x2 = rnorm(300L))
  latent <- 0.4 * dat$x1 - 0.3 * dat$x2 + rlogis(nrow(dat))
  dat$y <- ordered(cut(latent, c(-Inf, -0.5, 0.7, Inf), labels = 1:3))
  fit <- mfp2(y ~ fp(x1, df = 1, select = 1, center = FALSE, shift = 0, scale = 1) +
                fp(x2, df = 1, select = 1, center = FALSE, shift = 0, scale = 1),
              data = dat, family = ordinal_family(), verbose = FALSE)
  native <- rms::orm(y ~ x1 + x2, data = dat, x = TRUE, y = TRUE)
  for (rows in list(seq_len(nrow(dat)), c(2L, 9L, 20L), 9L)) {
    nd <- dat[rows, c("x1", "x2"), drop = FALSE]
    result <- predict(fit, newdata = nd, type = "link", se.fit = TRUE)
    expected <- predict(native, newdata = nd, type = "lp", se.fit = TRUE)
    expect_equal(unname(result$fit), unname(expected$linear.predictors), tolerance = 1e-6)
    expect_equal(unname(result$se.fit), unname(expected$se.fit), tolerance = 1e-6)
    response <- predict(fit, newdata = nd, type = "response", se.fit = TRUE)
    expect_equal(as.numeric(response$fit),
                 as.numeric(predict(native, newdata = nd, type = "fitted.ind")),
                 tolerance = 1e-6)
  }
  training <- predict(fit, type = "link", se.fit = TRUE)
  reference <- predict(native, type = "lp", se.fit = TRUE)
  expect_equal(unname(training$fit), unname(reference$linear.predictors), tolerance = 1e-6)
  expect_equal(unname(training$se.fit), unname(reference$se.fit), tolerance = 1e-6)
})
