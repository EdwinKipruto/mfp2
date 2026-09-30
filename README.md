# `mfp2`

<!-- badges: start -->
[![R-CMD-check](https://github.com/EdwinKipruto/mfp2/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EdwinKipruto/mfp2/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

`mfp2` implements multivariable fractional polynomial (MFP) models and related
extensions. It performs variable selection and functional-form selection for
continuous covariates. Supported responses include Gaussian, binomial,
Poisson, Gamma, inverse Gaussian, and negative-binomial outcomes; unordered
multinomial and ordinal outcomes; Cox, parametric survival, and Fine--Gray
competing-risks models; and clustered marginal models fitted by generalized
estimating equations (GEE). Multinomial fitting uses `nnet`; one FP
transformation is selected per predictor and shared across logits, with
logit-specific coefficients. Negative-binomial models use the optional
`fastglm` backend, selected automatically for `family = "negbin"`.
Ordinal models use the optional `rms` package. GEE models use `geepack` and
require a cluster identifier supplied through `id` (with optional `waves`).

In addition to standard MFP modelling, `mfp2` provides:

- approximate cumulative distribution (ACD) transformations for selected
  sigmoid-shaped covariate effects;
- spike-at-zero modelling for semi-continuous covariates with a distinct zero
  component and a positive continuous component; and
- `mfpi()` for investigating interactions between a categorical grouping
  variable, such as treatment group, and continuous covariates using
  fractional polynomial submodels.

Both formula and matrix interfaces are available.

## Compatibility with existing software packages

`mfp2` closely emulates the functionality of the `mfp` and `mfpa` packages in
Stata. It extends the existing `mfp` package in R by providing:

- both matrix and formula interfaces for input data;
- sigmoid transformations through the ACD transformation;
- support for covariates with a spike at zero;
- fractional-polynomial interaction analyses between categorical grouping
  variables and continuous covariates;
- estimation and plotting of contrasts and partial linear predictors for
  investigating nonlinear effects; and
- computational optimizations to improve speed and usability.

## Installation

```r
# Install the development version from GitHub with pak
# install.packages("pak")
pak::pak("EdwinKipruto/mfp2")

# Alternatively, install it with remotes
# install.packages("remotes")
remotes::install_github("EdwinKipruto/mfp2")
```

## Quick start

### Fit an MFP model

For most users, the formula interface is the simplest way to fit an MFP model.
The `fp()` terms identify continuous covariates for functional-form selection.

```r
library(mfp2)

data("prostate")

fit <- mfp2(
  lpsa ~ fp(age) + fp(cavol) + fp(weight) + svi,
  data = prostate,
  family = "gaussian",
  verbose = FALSE
)
```

Use the standard model methods to inspect and use the fitted model.

```r
fit
summary(fit)
coef(fit)
get_selected_variable_names(fit)

predict(
  fit,
  newdata = prostate[1:5, ]
)

plots <- plot(fit)
plots[[1]]
```

The matrix interface is available when the predictors have already been
prepared as a numeric matrix.

```r
x <- as.matrix(prostate[, c("age", "cavol", "weight", "svi")])
y <- prostate$lpsa

fit_matrix <- mfp2(
  x = x,
  y = y,
  family = "gaussian",
  verbose = FALSE
)

fit_matrix
```

### Supported model families

Choose `family` according to the response. These specifications work with
`mfp2()`; the formula examples below show their required inputs.

| Response or model | `family` argument | Requirement |
| --- | --- | --- |
| Continuous Gaussian | `"gaussian"` | Numeric response |
| Binary or grouped binomial | `"binomial"` | Binary response or successes/failures |
| Counts | `"poisson"` | Nonnegative integer response |
| Positive continuous | `stats::Gamma(link = "log")` | Positive response |
| Positive continuous | `stats::inverse.gaussian(link = "log")` | Positive response |
| Overdispersed counts | `"negbin"` | Optional `fastglm`; backend chosen automatically |
| Unordered classes | `multinomial_family()` | At least three classes; `nnet` |
| Ordered categories | `ordinal_family()` | At least three ordered categories; optional `rms` |
| Right-censored survival | `"cox"` | `survival::Surv(time, status)` |
| Parametric survival | `survreg_family(dist = "weibull")` | `survival::Surv()` response |
| Competing risks | `finegray_family(etype = "cause1")` | Multi-state `survival::Surv()` response |
| Clustered marginal mean | `gee_family(gaussian())` | Contiguous clusters through `id` |

GEE also accepts binomial, Poisson, and Gamma marginal response families.
Quasi families are not supported by the likelihood-based GLM selection.

### Fit binomial, Poisson, Gamma, and inverse-Gaussian models

This reproducible data set illustrates each GLM response. The two
positive-response fits use different variance assumptions for the same
illustrative outcome.

```r
set.seed(2026)
n_family <- 240L
family_data <- data.frame(
  x = runif(n_family, 0.5, 3),
  z = rnorm(n_family)
)
eta <- with(family_data, 0.2 + 0.5 * log(x) + 0.25 * z)
mu <- exp(eta)
family_data$binary <- rbinom(n_family, 1, plogis(eta - 0.7))
family_data$count <- rpois(n_family, mu)
family_data$positive <- rgamma(n_family, shape = 5, scale = mu / 5)

fit_binomial <- mfp2(
  binary ~ fp(x, df = 2) + z, data = family_data,
  family = "binomial", verbose = FALSE
)
fit_poisson <- mfp2(
  count ~ fp(x, df = 2) + z, data = family_data,
  family = "poisson", verbose = FALSE
)
fit_gamma <- mfp2(
  positive ~ fp(x, df = 2) + z, data = family_data,
  family = stats::Gamma(link = "log"), verbose = FALSE
)
fit_inverse_gaussian <- mfp2(
  positive ~ fp(x, df = 2) + z, data = family_data,
  family = stats::inverse.gaussian(link = "log"), verbose = FALSE
)
```

### Fit a negative-binomial model

`fastglm` is optional; when installed, `family = "negbin"` selects that
backend automatically.

```r
if (requireNamespace("fastglm", quietly = TRUE)) {
  family_data$overdispersed_count <- rnbinom(n_family, mu = mu, size = 2)
  fit_negbin <- mfp2(
    overdispersed_count ~ fp(x, df = 2) + z,
    data = family_data, family = "negbin", verbose = FALSE
  )
}
```

### Fit a multinomial MFP model

Use a factor response with at least three levels, or a matrix of class counts.
The first level is the reference by default; choose another with
`multinomial_family()`.

```r
fit_multinomial <- mfp2(
  Species ~ fp(Sepal.Length) + fp(Petal.Length) + Sepal.Width,
  data = iris,
  family = multinomial_family(reference = "setosa"),
  verbose = FALSE
)

# One row per non-reference logit; the selected powers are shared.
coef(fit_multinomial)

# Class probabilities and predicted classes.
predict(fit_multinomial, iris[1:6, ], type = "response")
predict(fit_multinomial, iris[1:6, ], type = "class")
```

### Fit an ordinal MFP model

For an ordered response, `ordinal_family()` fits a proportional-odds
model. Install the optional `rms` package first. Supplying an ordered
factor makes the category order explicit.

```r
if (requireNamespace("rms", quietly = TRUE)) {
  ordinal_data <- family_data[, c("x", "z")]
  ordinal_data$rating <- ordered(cut(
    eta + rlogis(n_family),
    breaks = c(-Inf, 0.3, 1.3, Inf),
    labels = c("low", "middle", "high")
  ))
  fit_ordinal <- mfp2(
    rating ~ fp(x, df = 2) + z,
    data = ordinal_data, family = ordinal_family(),
    verbose = FALSE
  )
  predict(fit_ordinal, ordinal_data[1:3, ], type = "response")
}
```

### Fit Cox and parametric survival models

Cox models use `family = "cox"`. For an accelerated failure-time model,
use `survreg_family()`; Weibull is its default distribution. Both accept a
`survival::Surv()` response.

```r
event_time <- rexp(n_family, rate = 0.12 * exp(0.3 * eta))
censor_time <- rexp(n_family, rate = 0.08)
survival_data <- family_data[, c("x", "z")]
survival_data$time <- pmin(event_time, censor_time)
survival_data$status <- as.integer(event_time <= censor_time)

fit_cox <- mfp2(
  survival::Surv(time, status) ~ fp(x, df = 2) + z,
  data = survival_data, family = "cox", verbose = FALSE
)
fit_weibull <- mfp2(
  survival::Surv(time, status) ~ fp(x, df = 2) + z,
  data = survival_data,
  family = survreg_family(dist = "weibull"),
  verbose = FALSE
)
```

### Fit a Fine--Gray competing-risks model

The event response must be a multi-state `Surv()` object: the first
factor level represents censoring, and `etype` identifies the event of
interest. With `strata_action = "both"` (the default), the formula stratum
affects both censoring weights and the baseline subdistribution hazard.

```r
competing_data <- family_data[, c("x", "z")]
competing_data$time <- rexp(n_family, rate = 0.12 * exp(0.2 * eta))
competing_data$event <- factor(
  sample(c("censor", "cause1", "cause2"), n_family, replace = TRUE,
         prob = c(0.35, 0.40, 0.25)),
  levels = c("censor", "cause1", "cause2")
)
competing_data$site <- factor(rep(c("A", "B"), length.out = n_family))

fit_finegray <- mfp2(
  survival::Surv(time, event) ~ fp(x, df = 2) +
    survival::strata(site),
  data = competing_data,
  family = finegray_family(etype = "cause1", strata_action = "both"),
  verbose = FALSE
)
```

### Fit a GEE MFP model for clustered data

Supply a cluster identifier through `id` (clusters must occupy contiguous
rows). `gee_family()` defaults to an exchangeable working correlation and the
robust sandwich covariance; `std.err` also accepts `"jack"`, `"j1s"`, and
`"fij"`. The retained model is a native `geepack::geeglm` object.

```r
set.seed(2027)
clustered_data <- data.frame(
  id = rep(seq_len(50), each = 4),
  x1 = runif(200, 1, 4),
  x2 = rnorm(200)
)
cluster_effect <- rnorm(50, sd = 0.7)
clustered_data$y <- with(clustered_data, 1 + log(x1) + 0.3 * x2) +
  cluster_effect[clustered_data$id] + rnorm(nrow(clustered_data), sd = 0.4)

fit_gee <- mfp2(
  y ~ fp(x1, df = 2) + fp(x2, df = 1),
  data = clustered_data,
  family = gee_family(gaussian(), corstr = "exchangeable"),
  id = clustered_data$id,
  verbose = FALSE
)

# Robust (sandwich) covariance and coefficients.
coef(fit_gee)
vcov(fit_gee)
```

### Model a spike-at-zero covariate

A spike-at-zero covariate has a distinct group of zero observations and a
positive continuous component. Request spike-at-zero modelling with
`spike = TRUE` inside `fp()`.

```r
set.seed(1)

n <- 200
exposure <- numeric(n)
positive_rows <- sample(seq_len(n), size = 150)
exposure[positive_rows] <- rgamma(length(positive_rows), shape = 2, rate = 0.5)
age <- runif(n, 20, 80)
outcome <-
  1.5 * (exposure == 0) +
  2 * log(ifelse(exposure > 0, exposure, 1)) +
  0.02 * age +
  rnorm(n)

spike_data <- data.frame(outcome, exposure, age)

fit_spike <- mfp2(
  outcome ~ fp(exposure, spike = TRUE) + fp(age),
  data = spike_data,
  family = "gaussian",
  verbose = FALSE
)

fit_spike
plot(fit_spike, terms = "exposure")
```

With the matrix interface, identify spike-at-zero variables by name through
`spike_vars`.

```r
x_spike <- as.matrix(spike_data[, c("exposure", "age")])

fit_spike_matrix <- mfp2(
  x = x_spike,
  y = spike_data$outcome,
  spike_vars = "exposure",
  family = "gaussian",
  verbose = FALSE
)
```

### Investigate interactions with MFPI

`mfpi()` investigates whether the effects of selected continuous covariates
differ across the levels of a categorical grouping variable. In this example,
the effects of `cavol` and `age` are evaluated across the two `svi` groups.

```r
fit_interaction <- mfpi(
  lpsa ~ fp(age) + svi + fp(pgg45) + fp(cavol) + fp(weight) +
    fp(bph) + fp(cp),
  data = prostate,
  group_var = "svi",
  cont_vars = c("cavol", "age"),
  flex = "flex1",
  center = FALSE,
  include_group_var = TRUE,
  show_models = FALSE,
  verbose = FALSE
)

fit_interaction
summary(fit_interaction)
```

Plot group-specific fitted functions and their differences. Difference curves
are shown on the model's linear-predictor scale.

```r
if (requireNamespace("patchwork", quietly = TRUE)) {
  plot(
    fit_interaction,
    terms = "cavol",
    plot_type = "both"
  )
}
```

## Documentation

The package vignettes provide detailed guidance:

```r
# Practical introduction to fitting, prediction, and plotting
vignette("mfp2_introduction", package = "mfp2")

# MFP methodology and practical considerations
vignette("MFP_Introduction", package = "mfp2")

# Approximate cumulative distribution transformations
vignette("mfp2_ACD", package = "mfp2")

# Spike-at-zero modelling
vignette("mfp2_spike", package = "mfp2")

# MFPI interaction analysis
vignette("mfpi", package = "mfp2")
```

Function-level help is available through:

```r
?mfp2
?fp
?mfpi
?multinomial_family
?ordinal_family
?survreg_family
?finegray_family
?gee_family
?predict.mfp2
?plot.mfp2
```

## References

To learn more about the MFP algorithm, see Royston and Sauerbrei (2008),
*Multivariable Model-Building: A Pragmatic Approach to Regression Analysis
based on Fractional Polynomials for Modelling Continuous Variables*.
John Wiley & Sons.

For details on the ACD transformation, see Royston (2014),
*A smooth covariate rank transformation for use in regression models with a
sigmoid dose-response function*. The Stata Journal.

For the spike-at-zero algorithm, see Becher et al. (2012),
*Analysing covariates with spike at zero: a modified FP procedure and conceptual
issues*. Biometrical Journal, 54(5), 686-700.

For fractional-polynomial interaction analyses between treatment and continuous
covariates, see Royston and Sauerbrei (2004), *A new approach to modelling
interactions between treatment and continuous covariates in clinical trials by
using fractional polynomial submodels*. Statistics in Medicine, 23(16),
2509-2525. doi: [10.1002/sim.1815](https://doi.org/10.1002/sim.1815)

See also Royston and Sauerbrei (2013), *Interaction of treatment with a
continuous variable: simulation study of significance level for several methods
of analysis*. Statistics in Medicine, 32(22), 3788-3803. doi:
[10.1002/sim.5813](https://doi.org/10.1002/sim.5813)

For simulation results on power, see Royston and Sauerbrei (2014), *Interaction
of treatment with a continuous variable: simulation study of power for several
methods of analysis*. Statistics in Medicine, 33(27), 4695-4708. doi:
[10.1002/sim.6308](https://doi.org/10.1002/sim.6308)
