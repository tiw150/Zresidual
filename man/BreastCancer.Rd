\name{BreastCancer}
\alias{BreastCancer}
\docType{data}
\title{German Breast Cancer Study Group 2 data}
\description{
A numeric recoding of recurrence-free survival data for 686 women with node-positive breast cancer in the German Breast Cancer Study Group 2 trial.
}
\usage{data(BreastCancer)}
\format{
A data frame with 686 rows and 12 variables:
\describe{
  \item{\code{X}}{Row index in the supplied data file.}
  \item{\code{id}}{Patient identifier.}
  \item{\code{time}}{Recurrence-free survival time, in days.}
  \item{\code{status}}{Event indicator: 0 for censored and 1 for recurrence or death.}
  \item{\code{treat}}{Tamoxifen indicator: 0 for no tamoxifen and 1 for tamoxifen.}
  \item{\code{age}}{Age in years.}
  \item{\code{men}}{Menopausal status: 1 for premenopausal and 2 for postmenopausal.}
  \item{\code{size}}{Tumor size, in millimeters.}
  \item{\code{grade}}{Tumor grade, coded 1 to 3.}
  \item{\code{nodes}}{Number of positive lymph nodes.}
  \item{\code{prog}}{Progesterone receptor level, in fmol.}
  \item{\code{oest}}{Estrogen receptor level, in fmol.}
}
}
\source{
Adapted from the \code{GBSG2} data distributed in the \pkg{TH.data} package.
}
\references{
Schumacher, M., et al. (1994). Randomized 2 x 2 trial evaluating hormonal treatment and the duration of chemotherapy in node-positive breast cancer patients. \emph{Journal of Clinical Oncology}, 12, 2086-2093. \doi{10.1200/JCO.1994.12.10.2086}.
}
\keyword{datasets}
