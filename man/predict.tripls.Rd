\name{predict.tripls}
\alias{predict.tripls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Tri-PLS Predictions

}
\description{
Make predictions based on a Trilinear Partial Least Squares regression model.  
}
\usage{
predict.tripls(res.tripls, data, casenames = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res.tripls}{A "tripls" class regression object containing Tri-PLS regression estimates
%%     ~~Describe \code{res.tripls} here~~
}
  \item{data}{An optional (n' x p x q) three-way data frame containing n' cases for which to calculate tri-PLS predictions
%%     ~~Describe \code{data} here~~
}
  \item{casenames}{An optional n' x 1 character array containing the corresponding case names.  
%%     ~~Describe \code{casenames} here~~
}
}
\details{This function implicitly adapts to the type of scaling used to construct the tri-PLS model, and rescales the predictions accordingly. 
%%  ~~ If necessary, more details than the description above ~~
}
\value{The predicted responses as a length n' numeric, named if names have been provided. 
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
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
res.tripls <- tripls(data_3$X3,data_3$y3,2)
predict(res.tripls)
predict(res.tripls,data=data_3$X3)
predict(res.tripls,data=data_3$X3,casename=paste("C",18:34,sep=""))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
