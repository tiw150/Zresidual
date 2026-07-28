\name{drs}
\alias{drs}
\docType{data}
\title{Diabetic Retinopathy Study data}
\description{
Paired-eye survival data from 197 high-risk diabetic patients. One randomly selected eye received laser treatment and the other eye was untreated; blindness was the event of interest.
}
\usage{data(drs)}
\format{
A data frame with 394 rows and 8 variables, with two rows for each patient:
\describe{
  \item{\code{subject_id}}{Unique patient identifier.}
  \item{\code{eye}}{Eye indicator: 1 for right and 2 for left.}
  \item{\code{time}}{Observed follow-up time.}
  \item{\code{status}}{Outcome indicator: 0 for censored and 1 for blindness.}
  \item{\code{treated}}{Treatment indicator: 0 for untreated and 1 for treated.}
  \item{\code{age_at_onset}}{Age at onset of diabetes, in years.}
  \item{\code{laser_type}}{Laser type: 1 for xenon and 2 for argon.}
  \item{\code{diabetes_type}}{Diabetes type: 1 for juvenile and 2 for adult.}
}
}
\source{
Adapted from the \code{drs} data distributed in the \pkg{frailtySurv} package.
}
\references{
Huster, W. J., Brookmeyer, R., and Self, S. G. (1989). Modelling paired survival data with covariates. \emph{Biometrics}, 45, 145-156. \doi{10.2307/2532041}.
}
\keyword{datasets}
