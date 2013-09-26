\name{coef.FitARMA}
\alias{coef.FitARMA}
\title{ coef method for class FitARMA }
\description{
produces table showing parameters, standard errors and Z-ratios
}

\usage{
\method{coef}{FitARMA}(object, ...)
}
\arguments{
  \item{object}{ class FitARMA object }
  \item{\dots}{ auxiliary parameters}
}

\value{matrix with 3 columns}

\author{ A.I. McLeod }
\seealso{ \code{\link{FitARMA}} }
\examples{
data(SeriesA)
out<-FitARMA(SeriesA, c(1,0,1))
coef(out)
}
\keyword{ ts }
