# Degrees-of-freedom calculation for MFPI interaction tests
#
# interaction_model_df() is called by test_interaction() to obtain the three df
# quantities needed for the LRT, AIC, and BIC:
#
#   dfmain    -- parameters in the main-effects model (excl. intercept & adj.)
#   total_df  -- parameters in the interaction model  (excl. intercept & adj.)
#   dfint     -- total_df - dfmain  (LRT degrees of freedom)
#
# Reference: Royston & Sauerbrei (2009), The Stata Journal 9(2), 230-251,
#            Table 1 and Section 7.1.


#' Degrees of Freedom for the Main-Effects and Interaction Models
#'
#' Computes the MFPI degrees-of-freedom accounting for the main-effects and
#' interaction models. These counts combine estimable regression coefficients
#' with the method's FP power-search allowances, while excluding the intercept
#' and adjustment-variable coefficients (which are common to both models and
#' therefore cancel in likelihood-ratio and AIC/BIC differences).
#'
#' @section Parameter counts:
#' Let \eqn{K} be the number of groups (`n_groups`), \eqn{m} the FP degree,
#' \eqn{R} the width of a prespecified linear design block (`linear_width`),
#' and \eqn{Q} the number of outcome logits (`n_logits`). For
#' scalar-response models, \eqn{Q = 1}. For multinomial models, regression
#' coefficients are logit specific but FP powers are shared across logits.
#'
#' The **main-effects model** contains:
#' \itemize{
#'   \item \eqn{K - 1} group dummy coefficients \eqn{\gamma_1, \ldots,
#'     \gamma_{K-1}},
#'   \item \eqn{J} parameters for the interaction variable, where \eqn{J}
#'     depends on the term type:
#'     \itemize{
#'       \item \strong{FP} (\eqn{m \geq 1}): \eqn{J = 2m}, counting \eqn{m}
#'         regression coefficients and \eqn{m} estimated FP powers (the
#'         **mfp2** convention of counting both).
#'       \item \strong{Linear} (\eqn{m = 0}): \eqn{J = R}, counting the
#'         coefficients in the complete linear design block. Thus \eqn{R=1}
#'         for continuous or binary variables, while a categorical factor with
#'         \eqn{L} levels normally has \eqn{R=L-1}. There are no FP powers to
#'         estimate, so the general formula
#'         \eqn{J = 2m} does not apply; the linear term is a special case.
#'     }
#' }
#' Hence
#' \deqn{
#'   p_{\text{main}} = (K - 1) + J,
#' }
#' giving \eqn{p_{\text{main}} = (K-1) + 2m} for FP terms and
#' \eqn{p_{\text{main}} = (K-1)+R} for a linear design block.
#'
#' The **interaction model** adds group-specific FP slopes, and depending on
#' the `flex` level, may also add group-specific FP powers:
#'
#' \describe{
#'   \item{`flex0` (linear, \eqn{m = 0})}{Each group gets its own
#'     \eqn{R}-coefficient block.
#'     \deqn{df_{\text{int}} = (K - 1)R, \quad
#'           p_{\text{int}} = (K-1)+KR.}}
#'
#'   \item{`flex1` / `flex2` / `flex3`}{Powers are constrained to be equal
#'     across groups (either taken from the pooled model or estimated jointly).
#'     Only regression coefficients differ: \eqn{K \times m} regression
#'     coefficients plus \eqn{m} shared powers \eqn{= Km + m = (K+1)m},
#'     but the power count is already absorbed into the main-effects df, so
#'     the extra parameters relative to the main model are the
#'     \eqn{(K-1)m} additional regression coefficients.
#'     \deqn{df_{\text{int}} = (K-1)\,m, \quad
#'           p_{\text{int}} = (K-1) + (K+1)m.}}
#'
#'   \item{`flex4`}{Each group has its own power set. The interaction model
#'     has \eqn{K \times m} regression coefficients and \eqn{K \times m}
#'     power parameters (versus \eqn{m} shared in the main model), giving
#'     \eqn{(K-1)m} extra regression parameters and \eqn{(K-1)m} extra power
#'     parameters.
#'     \deqn{df_{\text{int}} = 2\,(K-1)\,m, \quad
#'           p_{\text{int}} = (K-1) + 2m + 2\,(K-1)\,m.}}
#' }
#'
#' For \eqn{Q > 1}, the implemented common-power multinomial formulas are
#' \deqn{p_{main}=Q(K-1+R)} for a linear design block and
#' \deqn{p_{main}=Q(K-1+m)+m} for an FP term. The interaction increment is
#' \eqn{Q(K-1)R} for a linear design block, \eqn{Q(K-1)m} for
#' `flex1`--`flex3`, and
#' \eqn{(K-1)m(Q+1)} for `flex4`. Setting \eqn{Q=1} reduces these expressions
#' to the scalar-response formulas above.
#'
#' @section Verification table (K = 2):
#' \tabular{llrrr}{
#'   `flex`  \tab `degree` \tab \eqn{p_{\text{main}}} \tab
#'   \eqn{p_{\text{int}}} \tab \eqn{df_{\text{int}}} \cr
#'   flex0   \tab 0        \tab 2                     \tab 3               \tab 1 \cr
#'   flex1/2/3 \tab 1      \tab 3                     \tab 4               \tab 1 \cr
#'   flex1/2/3 \tab 2      \tab 5                     \tab 7               \tab 2 \cr
#'   flex4   \tab 1        \tab 3                     \tab 5               \tab 2 \cr
#'   flex4   \tab 2        \tab 5                     \tab 9               \tab 4 \cr
#' }
#'
#' @section Verification table (K = 3):
#' \tabular{llrrr}{
#'   `flex`  \tab `degree` \tab \eqn{p_{\text{main}}} \tab
#'   \eqn{p_{\text{int}}} \tab \eqn{df_{\text{int}}} \cr
#'   flex0   \tab 0        \tab 3                     \tab 5               \tab 2 \cr
#'   flex1/2/3 \tab 1      \tab 4                     \tab 6               \tab 2 \cr
#'   flex1/2/3 \tab 2      \tab 6                     \tab 10              \tab 4 \cr
#'   flex4   \tab 1        \tab 4                     \tab 8               \tab 4 \cr
#'   flex4   \tab 2        \tab 6                     \tab 14              \tab 8 \cr
#' }
#'
#' @section Note on non-nested models:
#' For `flex3` and `flex4`, the main-effects and interaction models may use
#' different FP families. The models are therefore not strictly nested and the
#' LRT p-value should be treated as indicative rather than exact (Royston &
#' Sauerbrei 2014). The degrees of freedom returned here follow the same
#' formulas as `flex1`/`flex2` for `flex3`, and the `flex4` formula for
#' `flex4`, as tabulated above.
#'
#' @param n_groups Positive integer \eqn{\geq 2}. Number of distinct levels of
#'   the grouping variable (`group_var`).
#' @param degree Non-negative integer. FP degree: `0` = linear (no FP
#'   transformation), `1` = FP1, `2` = FP2.
#' @param n_logits Positive integer. Number of independently parameterized
#'   outcome logits. Use `1` for scalar-response models and `C - 1` for a
#'   `C`-class multinomial model. FP powers remain common across these logits.
#' @param linear_width Positive integer. Number of columns in the prespecified
#'   linear focal-variable block. It is `1` for continuous and binary terms and
#'   typically one fewer than the number of levels for a treatment-coded factor.
#' @param flex Character string; one of `"flex0"`, `"flex1"`, `"flex2"`,
#'   `"flex3"`, or `"flex4"`. Controls the flexibility of the interaction
#'   model. See [mfp2::mfpi()] for a full description.
#'
#' @return A named list with four elements:
#' \describe{
#'   \item{`dfmain`}{\eqn{p_{\text{main}}}: main-model df contribution,
#'     including any FP power-search allowance and excluding intercept and
#'     adjustment terms.}
#'   \item{`total_df`}{\eqn{p_{\text{int}}}: interaction-model df
#'     contribution, including any FP power-search allowance and excluding
#'     intercept and adjustment terms.}
#'   \item{`dfint`}{\eqn{df_{\text{int}} = p_{\text{int}} - p_{\text{main}}}:
#'     degrees of freedom for the likelihood-ratio test.}
#'   \item{`n_groups`}{The value of `n_groups` passed in, returned for
#'     downstream use.}
#' }
#'
#' @examples
#' \dontrun{
#' # FP1, two groups (K=2): dfint = (K-1)*m = 1
#' interaction_model_df(n_groups = 2, degree = 1, flex = "flex1")
#'
#' # FP2, three groups (K=3): dfint = (K-1)*m = 4
#' interaction_model_df(n_groups = 3, degree = 2, flex = "flex1")
#'
#' # flex4, FP1, two groups: dfint = 2*(K-1)*m = 2
#' interaction_model_df(n_groups = 2, degree = 1, flex = "flex4")
#'}
#' @references
#' Royston, P. and Sauerbrei, W. (2009). Two techniques for investigating
#' interactions between treatment and continuous covariates in clinical trials.
#' *The Stata Journal*, 9(2), 230–251. (Table 1 and Section 7.1.)
#'
#' Royston, P. and Sauerbrei, W. (2014). Interaction of treatment with a
#' continuous variable: simulation study of power for several methods of
#' analysis. *Statistics in Medicine*, 33, 4695–4708.
#'
#' @keywords internal
#' @noRd
interaction_model_df <- function(n_groups,
                                 degree,
                                 n_logits = 1L,
                                 flex = c("flex0", "flex1", "flex2",
                                          "flex3", "flex4"),
                                 linear_width = 1L) {
  
  flex <- match.arg(flex)
  
  if (!is.numeric(n_groups) || length(n_groups) != 1L || n_groups < 2L) {
    stop("`n_groups` must be a single integer >= 2.", call. = FALSE)
  }
  if (!is.numeric(degree) || length(degree) != 1L || degree < 0L) {
    stop("`degree` must be a single non-negative integer.", call. = FALSE)
  }
  
  k <- as.integer(n_groups)
  m <- as.integer(degree)
  if (!is.numeric(n_logits) || length(n_logits) != 1L || anyNA(n_logits) ||
      n_logits < 1L || n_logits != as.integer(n_logits)) {
    stop("`n_logits` must be a single integer >= 1.", call. = FALSE)
  }
  q <- as.integer(n_logits)
  if (!is.numeric(linear_width) || length(linear_width) != 1L ||
      anyNA(linear_width) || linear_width < 1L ||
      linear_width != as.integer(linear_width)) {
    stop("`linear_width` must be a single integer >= 1.", call. = FALSE)
  }
  r <- as.integer(linear_width)
  
  # ---------------------------------------------------------------------------
  # p_main: parameters in the main-effects model
  #   (K-1) group dummies + J FP terms, where J = 2m (or 1 for linear)
  # ---------------------------------------------------------------------------
  dfmain <- if (m == 0L) {
    q * (k - 1L + r)
  } else {
    q * (k - 1L + m) + m
  }
  
  # ---------------------------------------------------------------------------
  # df_int: extra parameters in the interaction model vs the main model
  #
  #   flex0           : (K-1) extra slopes,                 df_int = K-1
  #   flex1/flex2/flex3: (K-1)*m extra regression coefs,    df_int = (K-1)*m
  #   flex4           : (K-1)*m regression + (K-1)*m powers, df_int = 2*(K-1)*m
  #
  # (flex0 forces degree = 0 in flex_fit(), so the degree == 0 branch below
  #  also covers flex0 correctly.)
  # ---------------------------------------------------------------------------
  if (m == 0L) {
    # Linear / flex0: one extra slope per non-reference group; no powers
    dfint <- q * (k - 1L) * r
  } else if (flex == "flex4") {
    # Group-specific powers remain common across outcome logits.
    dfint <- (k - 1L) * m * (q + 1L)
  } else {
    # Shared powers cancel; only group-specific regression coefficients differ.
    dfint <- q * (k - 1L) * m
  }
  
  total_df <- dfmain + dfint
  
  list(
    dfmain   = dfmain,
    dfint    = dfint,
    total_df = total_df,
    n_groups = k
  )
}
