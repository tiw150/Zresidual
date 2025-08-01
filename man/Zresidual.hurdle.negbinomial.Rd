% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Zresidual.hurdle.negbinomial.R
\name{Zresidual.hurdle.negbinomial}
\alias{Zresidual.hurdle.negbinomial}
\title{A function to calculate z-residuals of a 'brm' fit.
This function is to be used when the user needs to calculate the Z-residuals of TNB/HNB}
\usage{
Zresidual.hurdle.negbinomial(fit, type, method = "iscv", nrep = 1)
}
\arguments{
\item{fit}{is the fit object are from 'brms' package.}
}
\value{
\itemize{
 \item{Zresid}{Z-residual}
}
}
\description{
A function to calculate z-residuals of a 'brm' fit.
This function is to be used when the user needs to calculate the Z-residuals of TNB/HNB
}
