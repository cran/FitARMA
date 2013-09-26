\name{summary.FitARMA}
\alias{summary.FitARMA}
\title{ print method for class FitARMA }
\description{
a summary is printed out of the fitted model
}

\usage{
\method{summary}{FitARMA}(object, ...)
}

\arguments{
  \item{object}{ object, class FitARMA }
  \item{\dots}{ optional arguments }
}

\value{the result is displayed}

\author{ A.I. McLeod }
\seealso{ \code{\link{FitARMA}} }
\examples{
data(SeriesA)
out<-FitARMA(SeriesA, c(1,0,1))
summary(out)
}
\keyword{ ts }

