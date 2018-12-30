\name{dtrids}
\alias{dtrids}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Jacobian Matrix of the Tri-PLS regression coefficients
}
\description{
This function computes the Jacobian matrix of the Tri-PLS vector of regression coefficients with respect to the vector of cross products. 
}
\usage{
dtrids(X, y, p, q, h)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
The unfolded X data as a matrix (n x pq). 
}
  \item{y}{
The corresponding y data as an n x 1 vector or a length n numeric
}
  \item{p}{
The number of variables in the p mode (integer)
}
  \item{q}{
The number of variables in the q mode (integer)
}
  \item{h}{
The number of latent variables to compute (integer)
}
}
\details{
}
\value{
The Jacobian matrix db/ds (pq x pq)
}
\references{S. Serneels, K. Faber, T. Verdonck, P.J. Van Espen, Case specific prediction intervals for tri-PLS1: The full local linearization. Chemometrics and Intelligent Laboratory Systems, 108 (2011), 93-99.

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
J <- dtrids(array(data_3$X3,dim=c(17,33)),data_3$y3,11,3,2)}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

