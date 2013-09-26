\name{GetFitARMA}
\alias{GetFitARMA}
\title{ Fit ARMA(p,q) model to mean zero time series. }
\description{
The algorithm of McLeod and Zhang (2007) is used.
}
\usage{
GetFitARMA(y, p, q, pApprox = 30, init = 0)
}
\arguments{
  \item{y}{ time series }
  \item{p}{ AR order }
  \item{q}{ MA order }
  \item{pApprox}{ AR approximation  }
  \item{init}{ initial parameter estimates }
}
\details{
  See McLeod and Zhang (2006).
}
\value{
  \item{loglikelihood}{value of maximized loglikelihood}
  \item{phiHat}{estimated phi parameters}
  \item{thetaHat}{estimated theta parameters}
  \item{convergence}{result from optim}
  \item{algorithm}{ indicates "L-BFGS-B" or "Nelder-Mead" according as
  which algorithm was used in optim}
}

\references{
A.I. McLeod andY. Zhang (2008), Faster ARMA maximum likelihood estimation, 
Computational Statistics & Data Analysis, 52-4, 2166-2176.   
DOI link: http://dx.doi.org/10.1016/j.csda.2007.07.020}

\author{ A.I. McLeod, aimcleod@uwo.ca }

\seealso{ 
\code{\link{arima}}, 
\code{\link{FitARMA}} 
}
\examples{
data(SeriesA)
z<-SeriesA-mean(SeriesA)
GetFitARMA(z, 1, 1)
w<-diff(z, differences=1)
GetFitARMA(w, 0, 1)
}
\keyword{ ts }
