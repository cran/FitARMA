\name{residuals.FitARMA}
\alias{residuals.FitARMA}
\title{ residuals method for class FitARMA }
\description{
The innovation residuals are obtained.
}

\usage{
\method{residuals}{FitARMA}(object, ...)
}

\arguments{
  \item{object}{ class FitARMA object }
  \item{\dots}{ auxiliary parameters}
}

\value{vector or ts object}

\author{ A.I. McLeod }
\seealso{ \code{\link{FitARMA}} }
\examples{
data(SeriesA)
out<-FitARMA(SeriesA, c(1,0,1))
resid(out)
}
\keyword{ ts }


