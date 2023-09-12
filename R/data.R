#' Artificial dataset with continuous response
#'
#' The ART data set mimic the GBSG breast cancer study in terms of the 
#' distribution of predictors and correlation structure
#'
#' @name art
#' @docType data
#' @usage data(art)
#' @keywords data
#' @format The dataset has 250 observations and 10 covariates
#' \describe{
#'   \item{y}{continuous response variable.}
#'   \item{x1, x3, x5-x7, x10}{continous covariates.}
#'   \item{x2}{binary variable.}
#'   \item{x4}{ordinal variable with 3 levels.}
#'   \item{x8}{binary variable.}
#'   \item{x9}{nominal variable with 3 levels}
#' }
"art"

#' Breast cancer data set used in Royston and Sauerbrei (2008) book
#'
#' @name gbsg
#' @docType data
#' @usage data(gbsg)
#' @keywords data
#' @format A dataset with 686 observations and 11 variables
#' \describe{
#'   \item{id}{patient identifier.}
#'   \item{age}{age in years.}
#'   \item{meno}{menopausal status (0 = premeno, 1 = postmeno).}
#'   \item{size}{tumor size(mm)}
#'   \item{grade}{tummor grade}
#'   \item{nodes}{number of positive lymph nodes}
#'   \item{enodes}{exp(-0.12*nodes)}
#'   \item{pgr}{progestrone receptor status}
#'   \item{er}{oestrogen receptor status}
#'   \item{hormon}{tamoxifen treatment}
#'   \item{rectime}{time(days) to death or cancer recurrence}
#'   \item{censrec}{censoring(0 = censored, 1 = event)}
#' }
"gbsg"

#' Pima Indians dataset used in Royston and Sauerbrei (2008) book
#'
#' The dataset arises from an investigation of potential predictors of
#' the onset of diabetes in a cohort of 768 female Pima Indians of whom 268
#' developed diabetes. Missing values were imputed using the ice procedure
#' for stata
#'
#' @name pima
#' @docType data
#' @usage data(pima)
#' @keywords data
#' @format A data set with 768 observations and 9 variables.
#' \describe{
#'   \item{id}{patient identifier.}
#'   \item{pregnant}{number of times pregnant.}
#'   \item{glucose}{plasma glucose concentration at 2h in an oral glucose tolerance test.}
#'   \item{diastolic}{diastolic blood pressure in mmHg}
#'   \item{triceps}{triceps skin fold thickness in mm}
#'   \item{insulin}{2-h serum insulin}
#'   \item{bmi}{body mass index}
#'   \item{diabetes}{diabetes pedigree function}
#'   \item{age}{age in years}
#'   \item{y}{binary outcome variable (diabetes, yes/no)}
#' }
"pima"

#' Prostate cancer dataset used in Royston and Sauerbrei (2008)
#'
#' @name prostate
#' @docType data
#' @usage data(prostate)
#' @keywords data
#' @format A data set with 97 observations and 8 variables.
#' \describe{
#'   \item{obsno}{observation number.}
#'   \item{age}{age in years.}
#'   \item{svi}{seminal vessel invasion (yes/no)}
#'   \item{pgg45}{percentage gleason score 4 or 5}
#'   \item{cavol}{cancer volume (mm)}
#'   \item{weight}{prostate weight (g)}
#'   \item{bph}{amount of benign prostatic hyperplasia (g)}
#'   \item{cp}{amount of capsular penetration (g)}
#'   \item{lpsa}{log PSA concentration-outcome variable}
#' }
"prostate"