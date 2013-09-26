\name{FitARMA}
\alias{FitARMA}
\title{ Fit ARMA/ARIMA using fast MLE algorithm }
\description{
Fits an ARIMA(p,d,q) model using the algorithm given in McLeod and Zhang (2007). 
}
\usage{
FitARMA(z, order = c(0, 0, 0), demean = TRUE, MeanMLEQ = FALSE, pApprox = 30, MaxLag = 30)
}
\arguments{
  \item{z}{ time series }
  \item{order}{ model order, c(p,d,q) }
  \item{demean}{ if TRUE, mean parameter included otherwise assumed zero }
  \item{MeanMLEQ}{ exact MLE for mean, ignored unless demean=TRUE }
  \item{pApprox}{ order of approximation to be used }
  \item{MaxLag}{ maximum number of lags for portmanteau test }
}
\details{
  See McLeod and Ying (2007).
}

\value{
  A list with class name "FitARMA" and components: 
  \item{loglikelihood }{value of the loglikelihood}
  \item{phiHat }{AR coefficients }
  \item{thetaHat }{MA coefficients }
  \item{sigsqHat }{innovation variance estimate}
  \item{muHat }{estimate of the mean}
  \item{covHat }{covariance matrix of the coefficient estimates}
  \item{racf }{residual autocorrelations }
  \item{LjungBox}{table of Ljung-Box portmanteau test statistics}
  \item{res}{innovation residuals, same length as z}
  \item{fits}{fitted values, same length as z}
  \item{demean }{TRUE if mean estimated otherwise assumed zero}
  \item{IterationCount }{number of iterations in mean mle estimation}
  \item{convergence }{value returned by optim -- should be 0}
  \item{MLEMeanQ }{TRUE if mle for mean algorithm used}
  \item{tsp}{tsp(z)}
  \item{call}{result from match.call() showing how the function was called}
  \item{ModelTitle}{description of model}
  \item{DataTitle}{returns attr(z,"title") }
}

\references{
A.I. McLeod andY. Zhang (2008), Faster ARMA maximum likelihood estimation, 
Computational Statistics & Data Analysis, 52-4, 2166-2176.   
DOI link: http://dx.doi.org/10.1016/j.csda.2007.07.020
}

\author{ A.I. McLeod, aimcleod@uwo.ca }

\note{  
When d>0 and demean=TRUE, the mean of the differenced series is estimated.
This corresponds to including a polynomial of degree d.

When d>0, the AIC/BIC are computed for the differenced series and so
they are not comparable to the values obtained for models with d=0.

}

\seealso{ 
\code{\link{GetFitARMA}}, 
\code{\link{print.FitARMA}}, 
\code{\link{coef.FitARMA}}, 
\code{\link{residuals.FitARMA}}, 
\code{\link{fitted.FitARMA}}, 
\code{\link{arima}}
 }

\examples{
data(SeriesA) #in datasets()
out1<-FitARMA(SeriesA, c(1,0,1))
out1
coef(out1)
out2<-FitARMA(SeriesA, c(0,1,1))
out2
coef(out2)
}
\keyword{ ts }
