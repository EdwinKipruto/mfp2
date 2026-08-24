# Minimal MFPI objects used by printing tests.

make_mfpi_print_metrics <- function(n_groups = 2L) {
  group_names <- LETTERS[seq_len(n_groups)]
  int_powers <- setNames(
    lapply(seq_len(n_groups), function(i) c(i)),
    group_names
  )

  out <- data.frame(
    variable = "x",
    type = "fp1",
    deviance_int = 10.1234,
    deviance_diff = 2.3456,
    df_int = 2,
    pvalue = 0.01234,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out$fp_powers_main <- I(list(list(x = c(1, 2))))
  out$fp_powers_int <- I(list(int_powers))
  out
}

make_minimal_mfpi_print_object <- function() {
  structure(
    list(
      group_var = "group",
      nobs = 10L,
      flex = "flex3",
      criterion = "pvalue",
      p_adjust_method = "holm",
      p_interact = 0.05,
      min_improvement = 2,
      digits = 3L,
      adjust_terms = NULL,
      all_model_metrics = NULL,
      best_model_metrics = NULL,
      var_winners = NULL
    ),
    class = "mfpi"
  )
}


make_mfpi_adjustment_object <- function(criterion = "pvalue") {
  x <- make_minimal_mfpi_print_object()
  x$criterion <- criterion
  x$adjust_terms <- data.frame(
    df_setting = c(1, 4),
    df_initial = c(1, 4),
    select = if (criterion == "pvalue") c(0.05, 0.10) else toupper(criterion),
    alpha = if (criterion == "pvalue") c(0.05, 0.05) else toupper(criterion),
    df_final = c(1, 4),
    power1 = c(1, 3),
    power2 = c(NA_real_, 3),
    selected = c(TRUE, TRUE),
    row.names = c("hx", "age"),
    check.names = FALSE
  )
  x
}

make_mfpi_interaction_summary_object <- function(criterion = "pvalue") {
  x <- make_minimal_mfpi_print_object()
  x$criterion <- criterion
  x$p_adjust_method <- "holm"
  x$p_interact <- 0.05
  x$min_improvement <- 2

  if (criterion == "pvalue") {
    age_metric <- data.frame(pvalue = 0.01, p_adjusted = 0.02)
    wt_metric <- data.frame(pvalue = 0.20, p_adjusted = 0.20)
  } else if (criterion == "aic") {
    age_metric <- data.frame(AIC_main_minus_int = 3.5)
    wt_metric <- data.frame(AIC_main_minus_int = 1.5)
  } else {
    age_metric <- data.frame(BIC_main_minus_int = 4.0)
    wt_metric <- data.frame(BIC_main_minus_int = 0.5)
  }

  x$var_winners <- list(
    age = list(fit = list(ok = TRUE), metric = age_metric, type = "fp2"),
    wt = list(fit = list(ok = TRUE), metric = wt_metric, type = "fp1")
  )
  x$best_model_metrics <- data.frame(
    variable = "age",
    stringsAsFactors = FALSE
  )
  x
}
