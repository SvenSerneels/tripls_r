\name{summary.tripls}
\alias{summary.tripls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summary of Tri-PLS regression results 
%%  ~~function to do ... ~~
}
\description{
Some summary statistics of a Tri-PLS regression model 
}
\usage{
summary.tripls(res.tripls)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res.tripls}{A "tripls" class regression object containing Tri-PLS regression estimates
%%     ~~Describe \code{res.tripls} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Does not return output. 
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

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(data_3)
summary(tripls(data_3$X3,data_3$y3,2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.