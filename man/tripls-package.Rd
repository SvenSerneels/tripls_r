\name{tripls-package}
\alias{tripls-package}
\alias{tripls}
\docType{package}
\title{
Trilinear partial Least squares Regression
}
\description{
This packge provides an R implementation of trilinear partial least squares regression. it includes the Tri-PLS estimator, a predict() function, a summary() function, a function to compute the Jacobian matrix, a function to compute confidence and prediction intervals as well as a versatile plot() function. 
}
\details{
\tabular{ll}{
Package: \tab tripls\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-04-18\cr
License: \tab GNU \cr
}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Sven Serneels 

Maintainer: sven.serneels@basf.com
}
\references{
See individual functions
}

\keyword{TriPLS}
\keyword{Regression}
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
data(data_3)
res.tripls <- tripls(data_3$X3,data_3$y3,2)
plot(res.tripls, type = "coefficients")
plot(res.tripls, type = "yyp")
}
