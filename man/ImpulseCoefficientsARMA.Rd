\name{ImpulseCoefficientsARMA}
\alias{ImpulseCoefficientsARMA}
\title{ Impulse coefficients of ARMA}
\description{
The coefficients in the infinite MA expansion of the ARMA model are determined.
}
\usage{
ImpulseCoefficientsARMA(phi, theta, lag.max)
}
\arguments{
  \item{phi}{ AR coefficients }
  \item{theta}{ MA coefficients }
  \item{lag.max}{ lags 0,...,lag.max }
}
\value{
vector length lag.max+1
}

\author{ A.I. McLeod }

\examples{
ImpulseCoefficientsARMA(0.9,0.5,20)
}
\keyword{ ts }
