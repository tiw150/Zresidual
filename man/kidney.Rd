\name{kidney}
\alias{kidney}
\docType{data}
\title{Kidney catheter infection data}
\description{
Recurrence times to infection at the catheter insertion point for 38 kidney patients using portable dialysis equipment. Each patient has exactly two observations.
}
\usage{data(kidney)}
\format{
A data frame with 76 rows and 7 variables:
\describe{
  \item{\code{id}}{Patient identifier.}
  \item{\code{time}}{Observed recurrence or censoring time.}
  \item{\code{status}}{Event indicator: 0 for censored and 1 for infection.}
  \item{\code{age}}{Age in years.}
  \item{\code{sex}}{Sex: 1 for male and 2 for female.}
  \item{\code{disease}}{Disease type, a factor with levels \code{Other}, \code{GN}, \code{AN}, and \code{PKD}.}
  \item{\code{frail}}{Frailty estimate reported in the original study.}
}
}
\source{
Adapted from the \code{kidney} data distributed in the \pkg{survival} package.
}
\references{
McGilchrist, C. A. and Aisbett, C. W. (1991). Regression with frailty in survival analysis. \emph{Biometrics}, 47, 461-466. \doi{10.2307/2532138}.
}
\keyword{datasets}
