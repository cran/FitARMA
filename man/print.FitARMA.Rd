\name{print.FitARMA}
\alias{print.FitARMA}
\title{ print method for class FitARMA }
\description{
a brief summary is printed out of the fitted model
}

\usage{
\method{print}{FitARMA}(x, ...)
}


\arguments{
  \item{x}{ object, class FitARMA }
  \item{\dots}{ optional arguments }
}

\value{the result is displayed}

\author{ A.I. McLeod }
\seealso{ \code{\link{FitARMA}} }
\examples{
data(SeriesA)
FitARMA(SeriesA, c(1,0,1))
}
\keyword{ ts }


