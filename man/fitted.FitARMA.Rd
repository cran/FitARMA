\name{fitted.FitARMA}
\alias{fitted.FitARMA}
\title{ fitted method for class FitARMA }
\description{
  The fitted values are the observed minus residuals.
  If there is differencing, the observed values are those corresponding to the
  differenced time series.
}

\usage{
\method{fitted}{FitARMA}(object, ...)
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
fitted(out)
}
\keyword{ ts }


