\name{LeukSurv}
\alias{LeukSurv}
\docType{data}
\title{Acute myeloid leukemia survival data}
\description{
Survival data for 1,043 adult acute myeloid leukemia patients in northwest England, including residence coordinates, administrative district, and subject-specific prognostic factors.
}
\usage{data(LeukSurv)}
\format{
A data frame with 1,043 rows and 9 variables:
\describe{
  \item{\code{time}}{Survival time, in days.}
  \item{\code{cens}}{Right-censoring status: 0 for censored and 1 for dead.}
  \item{\code{xcoord}}{Residence coordinate on the x-axis.}
  \item{\code{ycoord}}{Residence coordinate on the y-axis.}
  \item{\code{age}}{Age in years.}
  \item{\code{sex}}{Sex: 0 for female and 1 for male.}
  \item{\code{wbc}}{White blood cell count at diagnosis, truncated at 500.}
  \item{\code{tpi}}{Townsend deprivation score; higher values indicate less affluent areas.}
  \item{\code{district}}{Administrative district of residence, among 24 districts.}
}
}
\source{
Adapted from the \code{LeukSurv} data distributed in the \pkg{spBayesSurv} package.
}
\references{
Henderson, R., Shimakura, S., and Gorst, D. (2002). Modeling spatial variation in leukemia survival data. \emph{Journal of the American Statistical Association}, 97, 965-972. \doi{10.1198/016214502388618753}.
}
\keyword{datasets}
