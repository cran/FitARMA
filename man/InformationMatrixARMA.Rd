\name{InformationMatrixARMA}
\alias{InformationMatrixARMA}
\title{ Expected large-sample information matrix for ARMA  }
\description{
The expected large-sample information matrix per observation for ARMA(p,q) models
is computed.
}
\usage{
InformationMatrixARMA(phi = numeric(0), theta = numeric(0))
}
\arguments{
  \item{phi}{ AR coefficients }
  \item{theta}{ MA coefficients }
}
\details{
The information matrix is derived by Box and Jenkins (1970).

}
\value{
a matrix of order (p+q)
}

\references{Box and Jenkins (1970). Time Series Analysis: Forecasting and Control.} 

\author{ A.I. McLeod }

\seealso{ \code{\link{FitARMA}} }

\examples{
#The covariance matrix estimates of the parameters phi and theta in an ARMA(1,1)
#with phi=0.9 and theta=0.5 and n=200 is
v<-solve(InformationMatrixARMA(0.9,0.5))/200
v
#and the standard errors are
sqrt(diag(v))
  }

\keyword{ ts }
