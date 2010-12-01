\name{FitARMA-package}
\alias{FitARMA-package}
\docType{package}
\title{ FitARMA: Fit ARMA or ARIMA using fast MLE algorithm}

\description{
Fit ARMA/ARIMA time series model using fast algorithm.
All MLE computations in R.
Two estimation functions: 'FitARMA' and 'GetFitARMA' are provided.
'FitARMA' provides more options including
an option for exact MLE estimation of the intercept term.
'GetFitARMA' is provided for bootstrapping and simulation experiments.
S3 Methods 'print', 'summary', 'coef', 'residuals', 'fitted' provided. 
Fast exact Gaussian ARMA simultation using C.
}

\details{
\tabular{ll}{
Package: \tab FitARMA\cr
Type: \tab Package\cr
Version: \tab 1.4\cr
Date: \tab 2010-12-01\cr
License: \tab GLP 2.0 or greater\cr
LazyLoad: \tab yes\cr
}

FitARMA is the main function. 
}

\author{A.I. McLeod}

\references{
A. I. McLeod, Ying Zhang (2007). Faster ARMA maximum likelihood estimation,
Computational Statistics & Data Analysis 52(4),
URL http://dx.doi.org/10.1016/j.csda.2007.07.020
}


\keyword{ package }
\keyword{ ts }

\seealso{
\code{\link{arima}},
\code{\link[FitAR]{AcfPlot}}

}

\examples{
data(SeriesA)
#ARIMA(0,1,1) with exact estimation of mean of differenced series
ans<-FitARMA(SeriesA, order=c(0,1,1), MeanMLEQ=TRUE)
ans
coef(ans)
#ARIMA(0,1,1) with sample-mean estimation of mean of differenced series
ans<-FitARMA(SeriesA, order=c(0,1,1))
ans
coef(ans)
#ARIMA(0,1,1) with mean of differenced series set to zero
#as in 'arima'
ans<-FitARMA(SeriesA, order=c(0,1,1), demean=FALSE)
ans
coef(ans)
# illustrating methods
summary(ans)
resid(ans)
fitted(ans)
ans$racf
#Simulate and fit Gaussian ARMA
z<-SimulateGaussianARMA(0.9, 0.5, 200)
#GetFitARMA is faster than FitARMA.
#Use GetFitARMA for parametric bootstrap and simulation experiments
GetFitARMA(z, p=1, q=1)

}
