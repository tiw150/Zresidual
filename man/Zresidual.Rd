\name{Zresidual}
\alias{Zresidual}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Z-residual
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Zresidual(fit.object, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit.object}{The fit object are ’coxph’, ’survreg’.
%%     ~~Describe \code{x} here~~
}
\item{data}{A data that used for fitting the survival model.
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{Zresid}{Z-residual}
  \item{censored.Zresid}{censored Z-residual}
  \item{SP}{Survival Probabilities}
  \item{ucs}{unmodified CS residuals}
  \item{mcs}{modified CS residuals}
  \item{martg}{Martingale residuals}
  \item{dev}{Deviance residuals}

%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
