#' Artificial dataset with continuous response
#'
#' The ART data set mimics the GBSG breast cancer study in terms of the
#' distribution of predictors and correlation structure.
#'
#' @name art
#' @docType data
#' @usage data(art)
#' @keywords data
#' @format The dataset has 250 observations and 11 variables (1 response and 10 covariates)
#' \describe{
#'   \item{y}{Continuous response variable.}
#'   \item{x1, x3, x5-x7, x10}{Continuous covariates.}
#'   \item{x2}{Binary variable.}
#'   \item{x4}{Ordinal variable with 3 levels.}
#'   \item{x8}{Binary variable.}
#'   \item{x9}{Nominal variable with 3 levels.}
#' }
"art"

#' Breast cancer dataset used in Royston and Sauerbrei (2008) book.
#'
#' @name gbsg
#' @docType data
#' @usage data(gbsg)
#' @keywords data
#' @format A dataset with 686 observations and 12 variables.
#' \describe{
#'   \item{id}{Patient identifier.}
#'   \item{age}{Age in years.}
#'   \item{meno}{Menopausal status (0 = premeno, 1 = postmeno).}
#'   \item{size}{Tumor size (mm).}
#'   \item{grade}{Tumor grade.}
#'   \item{nodes}{Number of positive lymph nodes.}
#'   \item{enodes}{exp(-0.12*nodes).}
#'   \item{pgr}{Progesterone receptor status.}
#'   \item{er}{Estrogen receptor status.}
#'   \item{hormon}{Tamoxifen treatment.}
#'   \item{rectime}{Time (days) to death or cancer recurrence.}
#'   \item{censrec}{Censoring (0 = censored, 1 = event).}
#' }
"gbsg"


#' Pima Indians dataset used in Royston and Sauerbrei (2008) book.
#'
#' The dataset arises from an investigation of potential predictors of
#' the onset of diabetes in a cohort of 768 female Pima Indians of whom 268
#' developed diabetes. Missing values were imputed using the ice procedure
#' for Stata.
#'
#' @name pima
#' @docType data
#' @usage data(pima)
#' @keywords data
#' @format A dataset with 768 observations and 10 variables.
#' \describe{
#'   \item{id}{Patient identifier.}
#'   \item{pregnant}{Number of times pregnant.}
#'   \item{glucose}{Plasma glucose concentration at 2h in an oral glucose tolerance test.}
#'   \item{diastolic}{Diastolic blood pressure in mmHg.}
#'   \item{triceps}{Triceps skin fold thickness in mm.}
#'   \item{insulin}{2-h serum insulin.}
#'   \item{bmi}{Body mass index.}
#'   \item{diabetes}{Diabetes pedigree function.}
#'   \item{age}{Age in years.}
#'   \item{y}{Binary outcome variable (diabetes, yes/no).}
#' }
"pima"


#' Prostate cancer dataset used in Royston and Sauerbrei (2008) book.
#'
#' @name prostate
#' @docType data
#' @usage data(prostate)
#' @keywords data
#' @format A dataset with 97 observations and 9 variables.
#' \describe{
#'   \item{obsno}{Observation number.}
#'   \item{age}{Age in years.}
#'   \item{svi}{Seminal vessel invasion (yes/no).}
#'   \item{pgg45}{Percentage Gleason score 4 or 5.}
#'   \item{cavol}{Cancer volume (mm).}
#'   \item{weight}{Prostate weight (g).}
#'   \item{bph}{Amount of benign prostatic hyperplasia (g).}
#'   \item{cp}{Amount of capsular penetration (g).}
#'   \item{lpsa}{Log PSA concentration (outcome variable).}
#' }
"prostate"

#' Advanced prostate cancer dataset used in Royston and Sauerbrei (2008) book.
#'
#' @name advanced_prostate_cancer
#' @docType data
#' @usage data(advanced_prostate_cancer)
#' @keywords data
#' @format A dataset with 475 observations (338 deaths) and 17 variables.
#' \describe{
#'   \item{patnr}{Patient number.}
#'   \item{age}{Age at diagnosis in years (continuous).}
#'   \item{wt}{Standardized weight (continuous).}
#'   \item{sbp}{Systolic blood pressure (continuous).}
#'   \item{dbp}{Diastolic blood pressure (continuous).}
#'   \item{sz}{Size of primary tumour (continuous).}
#'   \item{ap}{Serum acid phosphatase (continuous).}
#'   \item{hg}{Haemoglobin (g/100 ml) (continuous).}
#'   \item{sg}{Gleason stage-grade category (continuous).}
#'   \item{pf}{Performance status (binary).}
#'   \item{hx}{History of cardiovascular disease (binary).}
#'   \item{bm}{Presence of bone metastases (binary).}
#'   \item{stage}{Stage 4 vs stage 3 (binary).}
#'   \item{ekg}{Abnormal electrocardiogram (binary).}
#'   \item{rx}{Treatment group (binary).}
#'   \item{survtime}{Time to death (overall survival).}
#'   \item{cens}{Censoring (0 = censored, 1 = event).}
#'   }
"advanced_prostate_cancer"