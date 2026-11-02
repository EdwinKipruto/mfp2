# ---------------------------------------------------------------------------
# Independent test suite for the three new mfp2 model families
#   - multinomial_family()
#   - finegray_family()
#   - survreg_family()
#
# These helpers build well-conditioned synthetic data and provide small
# utilities used across the test files. They are intentionally independent of
# the package's own tests/testthat fixtures: the goal is to re-derive the
# expected behaviour from first principles and native survival / nnet fits.
#
# Design notes that the tests rely on (verified against mfp2 1.1.0):
#   * fp() applies a per-variable shift/scale and mfp2() centres transformed
#     terms by default. Slope coefficients are invariant to centring but the
#     reported intercept is not. To compare *coefficients* directly against a
#     native fit we therefore turn centring off (center = FALSE) and pin the
#     FP transform to the identity via fp(x, scale = 1, shift = 0). To compare
#     *predictions* we can leave the defaults on -- predictions are invariant.
#   * The formula interface does not force a variable into the model via the
#     top-level `select = 1`; forcing must be requested per term with
#     fp(x, select = 1). All equivalence fits below force their terms in.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(mfp2)
  library(survival)
  library(nnet)
})

# Tolerances -----------------------------------------------------------------
TOL_TIGHT  <- 1e-6   # linear-algebra-exact quantities (predictions, exact refits)
TOL_MED    <- 1e-4   # quantities that pass through a second optimiser
# nnet's BFGS converges to a looser tolerance than survreg/coxph Newton steps;
# multinomial log-likelihood / probability comparisons use TOL_NNET.
TOL_NNET   <- 1e-3

# Multinomial: 3-class factor response with two positive continuous predictors
make_multinomial_data <- function(n = 400, seed = 1) {
  set.seed(seed)
  x1 <- runif(n, 1, 8)
  x2 <- runif(n, 1, 8)
  lp2 <- -0.4 + 0.20 * x1 - 0.05 * x2
  lp3 <-  0.3 - 0.15 * x2 + 0.05 * x1
  cls <- vapply(seq_len(n), function(i) {
    p1 <- 1 / (1 + exp(lp2[i]) + exp(lp3[i]))
    sample(c("A", "B", "C"), 1L,
           prob = c(p1, exp(lp2[i]) * p1, exp(lp3[i]) * p1))
  }, character(1))
  data.frame(y = factor(cls), x1 = x1, x2 = x2)
}

# Competing-risks data with a multi-state Surv() response and a binary stratum.
# `strong_strata` inflates the stratum effect on the event/censoring processes
# so that the three strata_action modes yield visibly different coefficients.
make_finegray_data <- function(n = 600, seed = 9, strong_strata = FALSE) {
  set.seed(seed)
  z1 <- runif(n, 1, 10)
  grp <- factor(sample(c("t", "c"), n, TRUE))
  sx  <- factor(sample(c("m", "f"), n, TRUE))
  sx_eff <- if (strong_strata) 0.9 else 0.3
  eta <- 0.12 * (z1 - mean(z1)) + 0.4 * (grp == "t") + sx_eff * (sx == "f")
  tm  <- rexp(n, rate = 0.05 * exp(eta))
  # censoring time depends on stratum so stratified IPCW weights differ by mode
  cn  <- runif(n, 0, if (strong_strata) 20 else 30) *
    ifelse(sx == "f" & strong_strata, 0.6, 1)
  status <- ifelse(tm <= cn, sample(1:2, n, TRUE, prob = c(0.6, 0.4)), 0L)
  ev <- factor(status, levels = 0:2, labels = c("censor", "relapse", "death"))
  data.frame(obs = pmin(tm, cn), ev = ev, z1 = z1, grp = grp, sx = sx)
}

# Right-censored survival data (single event) for survreg tests. gbsg is the
# natural choice for real data; this generator is used where a controlled
# stratum with a per-stratum scale is needed.
make_survreg_data <- function(n = 500, seed = 5) {
  set.seed(seed)
  x1 <- runif(n, 1, 10)
  sx <- factor(sample(c("m", "f"), n, TRUE))
  # different scale (log-time SD) per stratum
  scale_s <- ifelse(sx == "f", 1.3, 0.7)
  lp <- 3 + 0.10 * x1
  logt <- lp + scale_s * (log(-log(runif(n))))     # Weibull-ish, per-stratum scale
  tt <- exp(logt)
  cn <- quantile(tt, 0.85)
  status <- as.integer(tt <= cn)
  data.frame(time = pmin(tt, cn), status = status, x1 = x1, sx = sx)
}

# Return the value of a testthat skip-friendly numeric max abs difference
max_abs_diff <- function(a, b) max(abs(as.numeric(a) - as.numeric(b)))

# ---- mfpi data generators (real group x covariate interaction) -------------

make_mfpi_multinomial <- function(n = 800, seed = 88) {
  set.seed(seed)
  a <- runif(n, 1, 10)
  g <- factor(sample(c("0", "1"), n, TRUE))
  bz <- ifelse(g == "1", 0.35, -0.10)          # z-effect differs by group
  lp2 <- -0.3 + bz * (a - mean(a))
  lp3 <-  0.2 - 0.10 * (a - mean(a))
  cls <- vapply(seq_len(n), function(i) {
    p1 <- 1 / (1 + exp(lp2[i]) + exp(lp3[i]))
    sample(c("A", "B", "C"), 1L, prob = c(p1, exp(lp2[i]) * p1, exp(lp3[i]) * p1))
  }, character(1))
  data.frame(y = factor(cls), a = a, g = g)
}

make_mfpi_survreg <- function(n = 700, seed = 5) {
  set.seed(seed)
  x1 <- runif(n, 1, 10)
  grp <- factor(sample(c("t", "c"), n, TRUE))
  sx  <- factor(sample(c("m", "f"), n, TRUE))
  bz <- ifelse(grp == "t", 0.15, -0.02)
  sc <- ifelse(sx == "f", 1.3, 0.7)
  logt <- 3 + bz * (x1 - mean(x1)) + sc * log(-log(runif(n)))
  tt <- exp(logt); cn <- quantile(tt, 0.85)
  data.frame(time = pmin(tt, cn), status = as.integer(tt <= cn),
             x1 = x1, grp = grp, sx = sx)
}

# recursively locate the first fitted object of a given class inside a list
find_fit <- function(x, cls) {
  if (inherits(x, cls)) return(x)
  if (is.list(x)) for (el in x) { r <- find_fit(el, cls); if (!is.null(r)) return(r) }
  NULL
}
