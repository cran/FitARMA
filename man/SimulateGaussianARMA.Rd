\name{SimulateGaussianARMA}
\alias{SimulateGaussianARMA}
\title{ Simulate Gaussian ARMA model }
\description{
  An exact simulation method is used to simulate Gaussian ARMA models.
}
\usage{
SimulateGaussianARMA(phi, theta, n, InnovationVariance = 1, UseC = TRUE)
}
\arguments{
  \item{phi}{ AR coefficients }
  \item{theta}{ MA coefficients }
  \item{n}{ length of series }
  \item{InnovationVariance}{ innovation variable, default is 1 }
  \item{UseC}{ if \code{UseC}=TRUE, use C code. Otherwise, use slower R code. }
}
\details{
The detailed description is given in Hipel and McLeod (1994, 2006).
}
\value{
a vector containing the time series
}
\references{
Hipel, K.W. and McLeod, A.I. (2006).
Time Series Modelling of Water Resources and Environmental Systems.
}
\author{ A.I. McLeod }


\seealso{ \code{\link{arima.sim}} }

\examples{
z<-SimulateGaussianARMA(0.9, 0.5, 200)
FitARMA(z, c(1,0,1))

}
\keyword{ ts }
