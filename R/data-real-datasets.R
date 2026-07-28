#' German Breast Cancer Study Group 2 data
#'
#' A numeric recoding of recurrence-free survival data for 686 women with node-positive breast cancer in the German Breast Cancer Study Group 2 trial.
#'
#' @format A data frame with 686 rows and 12 variables:
#' \describe{
#'   \item{\code{X}}{Row index in the supplied data file.}
#'   \item{\code{id}}{Patient identifier.}
#'   \item{\code{time}}{Recurrence-free survival time, in days.}
#'   \item{\code{status}}{Event indicator: 0 for censored and 1 for recurrence or death.}
#'   \item{\code{treat}}{Tamoxifen indicator: 0 for no tamoxifen and 1 for tamoxifen.}
#'   \item{\code{age}}{Age in years.}
#'   \item{\code{men}}{Menopausal status: 1 for premenopausal and 2 for postmenopausal.}
#'   \item{\code{size}}{Tumor size, in millimeters.}
#'   \item{\code{grade}}{Tumor grade, coded 1 to 3.}
#'   \item{\code{nodes}}{Number of positive lymph nodes.}
#'   \item{\code{prog}}{Progesterone receptor level, in fmol.}
#'   \item{\code{oest}}{Estrogen receptor level, in fmol.}
#' }
#' @source Adapted from the \code{GBSG2} data distributed in the \pkg{TH.data} package.
#' @references Schumacher, M., et al. (1994). Randomized 2 x 2 trial evaluating hormonal treatment and the duration of chemotherapy in node-positive breast cancer patients. \emph{Journal of Clinical Oncology}, 12, 2086-2093. \doi{10.1200/JCO.1994.12.10.2086}.
#' @docType data
#' @keywords datasets
#' @usage data(BreastCancer)
"BreastCancer"

#' Diabetic Retinopathy Study data
#'
#' Paired-eye survival data from 197 high-risk diabetic patients. One randomly selected eye received laser treatment and the other eye was untreated; blindness was the event of interest.
#'
#' @format A data frame with 394 rows and 8 variables, with two rows for each patient:
#' \describe{
#'   \item{\code{subject_id}}{Unique patient identifier.}
#'   \item{\code{eye}}{Eye indicator: 1 for right and 2 for left.}
#'   \item{\code{time}}{Observed follow-up time.}
#'   \item{\code{status}}{Outcome indicator: 0 for censored and 1 for blindness.}
#'   \item{\code{treated}}{Treatment indicator: 0 for untreated and 1 for treated.}
#'   \item{\code{age_at_onset}}{Age at onset of diabetes, in years.}
#'   \item{\code{laser_type}}{Laser type: 1 for xenon and 2 for argon.}
#'   \item{\code{diabetes_type}}{Diabetes type: 1 for juvenile and 2 for adult.}
#' }
#' @source Adapted from the \code{drs} data distributed in the \pkg{frailtySurv} package.
#' @references Huster, W. J., Brookmeyer, R., and Self, S. G. (1989). Modelling paired survival data with covariates. \emph{Biometrics}, 45, 145-156. \doi{10.2307/2532041}.
#' @docType data
#' @keywords datasets
#' @usage data(drs)
"drs"

#' Kidney catheter infection data
#'
#' Recurrence times to infection at the catheter insertion point for 38 kidney patients using portable dialysis equipment. Each patient has exactly two observations.
#'
#' @format A data frame with 76 rows and 7 variables:
#' \describe{
#'   \item{\code{id}}{Patient identifier.}
#'   \item{\code{time}}{Observed recurrence or censoring time.}
#'   \item{\code{status}}{Event indicator: 0 for censored and 1 for infection.}
#'   \item{\code{age}}{Age in years.}
#'   \item{\code{sex}}{Sex: 1 for male and 2 for female.}
#'   \item{\code{disease}}{Disease type, a factor with levels \code{Other}, \code{GN}, \code{AN}, and \code{PKD}.}
#'   \item{\code{frail}}{Frailty estimate reported in the original study.}
#' }
#' @source Adapted from the \code{kidney} data distributed in the \pkg{survival} package.
#' @references McGilchrist, C. A. and Aisbett, C. W. (1991). Regression with frailty in survival analysis. \emph{Biometrics}, 47, 461-466. \doi{10.2307/2532138}.
#' @docType data
#' @keywords datasets
#' @usage data(kidney)
"kidney"

#' Acute myeloid leukemia survival data
#'
#' Survival data for 1,043 adult acute myeloid leukemia patients in northwest England, including residence coordinates, administrative district, and subject-specific prognostic factors.
#'
#' @format A data frame with 1,043 rows and 9 variables:
#' \describe{
#'   \item{\code{time}}{Survival time, in days.}
#'   \item{\code{cens}}{Right-censoring status: 0 for censored and 1 for dead.}
#'   \item{\code{xcoord}}{Residence coordinate on the x-axis.}
#'   \item{\code{ycoord}}{Residence coordinate on the y-axis.}
#'   \item{\code{age}}{Age in years.}
#'   \item{\code{sex}}{Sex: 0 for female and 1 for male.}
#'   \item{\code{wbc}}{White blood cell count at diagnosis, truncated at 500.}
#'   \item{\code{tpi}}{Townsend deprivation score; higher values indicate less affluent areas.}
#'   \item{\code{district}}{Administrative district of residence, among 24 districts.}
#' }
#' @source Adapted from the \code{LeukSurv} data distributed in the \pkg{spBayesSurv} package.
#' @references Henderson, R., Shimakura, S., and Gorst, D. (2002). Modeling spatial variation in leukemia survival data. \emph{Journal of the American Statistical Association}, 97, 965-972. \doi{10.1198/016214502388618753}.
#' @docType data
#' @keywords datasets
#' @usage data(LeukSurv)
"LeukSurv"
