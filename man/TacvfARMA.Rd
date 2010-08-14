\name{TacvfARMA}
\alias{TacvfARMA}

\title{ Theoretical Autocovariance Function of ARMA }
\description{
The theoretical autocovariance function of an ARMA(p,q) with unit
variance is computed. This algorithm has many applications.
In this package it is used for the computation of the information
matrix, in simulating p initial starting values for AR
simulations and in the computation of the exact mle for the mean.
}

\usage{
TacvfARMA(phi = numeric(0), theta = numeric(0), lag.max = 20)
}

\arguments{
  \item{phi}{ AR coefficients }
  \item{theta}{ MA coefficients }
  \item{lag.max}{ computes autocovariances lags 0,1,...,lag.max }
}

\details{
The algorithm given by McLeod (1975) is used.

The built-in R function ARMAacf could also be used but it is quite
complicated and apart from the source code, the precise algorithm
used is not described.  The only reference given for ARMAacf is the
Brockwell and Davis (1991) but this text does not give any detailed
exact algorithm for the general case.

Another advantage of TacvfARMA over ARMAacf is that it will be easier
for to translate and implement this algorithm in other computing
environments such as MatLab etc.
}

\value{
vector of length lag.max+1 containing the autocovariances is returned
}

\references{
McLeod, A.I. (1975), 
Derivation of the theoretical autocorrelation function of autoregressive 
moving-average time series, 
Applied Statistics 24, 255-256.
} 

\author{ A.I. McLeod }

\seealso{ 
\code{\link{ARMAacf}}, 
\code{\link{InformationMatrixARMA}}
 }

\examples{
#calculate and plot the autocorrelations from an ARMA(1,1) model
# with parameters phi=0.9 and theta=0.5
 g<-TacvfARMA(0.9,0.5,20)
 AcfPlot(g/g[1], LagZeroQ=FALSE)
}

\keyword{ ts }


