\name{tccfAR}
\alias{tccfAR}
\title{ Theoretical cross-covariances of auxilary AR process in ARMA(p,q) }
\description{
The auxilary AR processes in the ARMA(p,q) model phi(B)z(t)=theta(B)a(t) are defined by
phi(B)u(t)=-a(t) and theta(B)v(t)=a(t).
The upper off-diagonal p-by-q block of the ARMA information matrix is obtained
from the cross-covariances of u(t) and v(t).
This function obtains these covariances.
}
\usage{
tccfAR(phi, theta)
}
\arguments{
  \item{phi}{ AR coefficients in ARMA }
  \item{theta}{ MA coefficients in ARMA }
}
\details{
A set of linear equations which determine the covariances is solved.
The algorithm is similar in spirit to that for the autocovariances (McLeod, 1975).
}
\value{ vector of cross-covariances
}

\references{
McLeod, A.I. (1975), 
Derivation of the theoretical autocorrelation function of autoregressive 
moving-average time series, 
Applied Statistics 24, 255-256.
} 

\author{ A.I. McLeod }

\seealso{ \code{\link{InformationMatrixARMA}} }

\examples{
tccfAR(0.9,0.5)
}

\keyword{ ts }
