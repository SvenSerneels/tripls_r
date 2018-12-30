\name{intervals.tripls}
\alias{intervals.tripls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Inference for Tri-PLS estimates
%%  ~~function to do ... ~~
}
\description{Depending if a data set of new cases is provided, this function provides either confidence intervals for the vector of Tri-PLS regression coefficients, or prediction intervals for the set of new cases. The intervals are computed using error propagation based on a local linearization of order 1. } 
\usage{
intervals.tripls(res.tripls, alpha = 0.05, Xnew)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res.tripls}{A "tripls" class regression object containing Tri-PLS regression estimates
%%     ~~Describe \code{res.tripls} here~~
}
  \item{alpha}{The significance level at which to compute the intervals (positive numeric)
%%     ~~Describe \code{alpha} here~~
}
  \item{Xnew}{Optional. If missing, intervals for the regression coefficients will be calculated. If provided, it should be a three-dimensional numeric array of dimensions n' x p x q, where the dimenions p and q should match those in res.tripls. In that case, predictions and prediction intervals for the n' new cases will be returned.   
%%     ~~Describe \code{Xnew} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Returns a list object containing 
  \item{intervals}{A data frame containing the prediction interval's lower limit (llim), the predicted values (either called coefficients or y_predicted, depending on the inputs) and the prediction interval's upper limit (ulim)}
  \item{df}{The estimated number of degrees of freedom}
%% ...
}
\references{
S. Serneels, K. Faber, T. Verdonck, P.J. Van Espen, Case specific prediction intervals for tri-PLS1: The full local linearization. Chemometrics and Intelligent Laboratory Systems, 108 (2011), 93-99.
}
\author{Sven Serneels, BASF Corp. 
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
 data(data_3)
 intervals <- intervals.tripls(tripls(data_3$X3,data_3$y3,1),alpha=.05)
 intervals <- intervals.tripls(tripls(data_3$X3,data_3$y3,1),alpha=.05,data_3$X3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
