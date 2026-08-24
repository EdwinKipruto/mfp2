# `mfp2`

<!-- badges: start -->
[![R-CMD-check](https://github.com/EdwinKipruto/mfp2/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EdwinKipruto/mfp2/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

`mfp2` implements multivariable fractional polynomial (MFP) models and related
extensions. It performs variable selection and functional-form selection for
continuous covariates. The package supports generalized linear models with
families `"gaussian"`, `"binomial"`, `"poisson"`, and `"negbin"`, as well as
Cox proportional hazards models for right-censored survival outcomes.
Negative-binomial models are specified with `family = "negbin"` and require
`fitter = "fastglm"`.

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
