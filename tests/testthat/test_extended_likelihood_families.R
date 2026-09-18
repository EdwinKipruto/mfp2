# Constructor, validation, and MFPI integration coverage for the likelihood
# families added to the public family interface. Forced-linear equivalence with
# the native model functions lives in the model-specific test files.


test_that("survreg distribution parameters follow the native interface", {
  y <- survival::Surv(c(1, 2, 3, 4), c(1, 1, 0, 1))
  weights <- rep(1, NROW(y))

  positional <- prepare_survreg_family(
    survreg_family(dist = "t", parms = 5),
    y,
    weights
  )
  expect_equal(positional$prepared$parms[["df"]], 5)

  named <- prepare_survreg_family(
    survreg_family(dist = "t", parms = c(df = 5)),
    y,
    weights
  )
  expect_equal(named$prepared$parms[["df"]], 5)
})


test_that("Fine--Gray id is required only for start-stop responses", {
  n_subject <- 18L
  id <- rep(seq_len(n_subject), each = 2L)
  start <- rep(c(0, 1), n_subject)
  stop <- start + 1
  terminal <- rep(c("relapse", "death", "censor"), length.out = n_subject)
  event <- factor(
    c(rbind(rep("censor", n_subject), terminal)),
    levels = c("censor", "relapse", "death")
  )
  y <- survival::Surv(start, stop, event)
  weights <- rep(1, length(id))
  family <- finegray_family(etype = "relapse")

  expect_error(
    prepare_finegray_family(family, y, weights),
    "require the `id` argument"
  )
  expect_s3_class(
    prepare_finegray_family(family, y, weights, id = id),
    "mfp2_finegray_family"
  )

  ordinary <- survival::Surv(seq_len(n_subject), factor(
    terminal,
    levels = c("censor", "relapse", "death")
  ))
  expect_s3_class(
    prepare_finegray_family(
      family,
      ordinary,
      rep(1, n_subject)
    ),
    "mfp2_finegray_family"
  )
})


test_that("Fine--Gray family validates endpoint specifications", {
  expect_error(finegray_family(etype = list("relapse")), "event-state")
  expect_error(finegray_family(etype = matrix("relapse")), "event-state")
  expect_warning(
    family <- finegray_family(etype = c("relapse", "death")),
    "first value"
  )
  expect_identical(family$etype, "relapse")
})


test_that("MFPI accepts the expanded likelihood-family objects", {
  set.seed(4104)
  n <- 150L
  dat <- data.frame(
    group = factor(rep(c("control", "treated"), each = n / 2L)),
    x = stats::runif(n, 0.5, 3),
    z = stats::runif(n, 0.5, 2)
  )
  mu <- exp(0.2 + 0.2 * dat$x + 0.15 * dat$z +
    0.1 * (dat$group == "treated") * dat$x)
  dat$y <- stats::rgamma(n, shape = 6, scale = mu / 6)

  fit <- mfpi(
    y ~ group + fp(x, df = 1) + z,
    data = dat,
    family = stats::Gamma(link = "log"),
    group_var = "group",
    interaction_vars = "x",
    interaction_forms = c(x = "linear"),
    flex = "flex1",
    cycles = 1,
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_identical(fit$family_string, "Gamma")
  expect_s3_class(fit$adjustment_model, "glm")
})
