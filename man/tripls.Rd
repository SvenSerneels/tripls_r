\name{tripls}
\alias{tripls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Trilinear Partial Least Squares Regression
%%  ~~function to do ... ~~
}
\description{This function computes trilinear partial least squares (Tri-PLS) regression estimates 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
tripls(X, y, A, scaling = TRUE, nnames = NULL, pnames = NULL, qnames = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{The X data as a 3 way (n x p x q) array.
%%     ~~Describe \code{X} here~~
}
  \item{y}{The corresponding y data as an n x 1 vector or a length n numeric
%%     ~~Describe \code{y} here~~
}
  \item{A}{The number of latent components to estimate (integer)
%%     ~~Describe \code{A} here~~
}
  \item{scaling}{Logical flag. If TRUE, the data are internally centered and scaled to unit variance. If FALSE, the data are only centered. 
%%     ~~Describe \code{scaling} here~~
}
  \item{nnames}{An optional n x 1 character matrix containing the case names. 
%%     ~~Describe \code{nnames} here~~
}
  \item{pnames}{An optional p x 1 character matrix containing p mode variable names.
%%     ~~Describe \code{pnames} here~~
}
  \item{qnames}{An optional q x 1 character matrix containing q mode variable names.
%%     ~~Describe \code{qnames} here~~
}
}
\details{
The actual algorithm is described in Reference [4], which is a minor modification of the implementations described in references [1-3]. 
}
\value{
Returns a class "tripls" regression object containing the indivdual Tri-PLS results, i.e.: 
  \item{coefficients}{The vector of regression coeficients (pq x 1)}
  \item{intercept}{The intercept (n x 1)}
  \item{scores}{The latent variables (or scores, n x A)}
  \item{fitted.values}{The fitted responses from an A component model }
  \item{W}{The combined weighting vectors (pq x A)}
  \item{WJ}{The p mode weighting vectors (p x A)}
  \item{WK}{The q mode weighting vectors (q x A)}
  \item{YMeans}{The y mean (length 1 numeric) )}
  \item{YScales}{The y scale (length 1 numeric, 1 if scaling=FALSE)}
  \item{XMeans}{The X columnwise means (length pq numeric)}
  \item{XScales}{The X columnwise scales (length pq numeric, all ones if scaling=FALSE)}
  \item{X.scaled}{The scaled, unfolded predictor matrix (n x pq)}
  \item{y.scaled}{The scaled response (length n numeric)}
  \item{sev}{The percentage of explained covariance}
  \item{rmsec}{The root mean squared error of calibration}
  \item{inputs}{A list object containing the input data}
  
  
%% ...
}
\references{
[1] L. Staahle, Aspects of the analysis of three-way data. Chemometrics and Intelligent Laboratory Systems, 7 (1989), 95-100.

[2] R. Bro, Multiway calibration. Multilinear PLS. Journal of Chemometrics, 10 (1996), 47-61.

[3] S. de Jong, Regression coefficients in multilinear PLS. Journal of Chemometrics, 12 (1998) 77-81.

[4] S. Serneels, K. Faber, T. Verdonck, P.J. Van Espen, Case specific prediction intervals for tri-PLS1: The full local linearization. Chemometrics and Intelligent Laboratory Systems, 108 (2011), 93-99.
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
res.tripls <- tripls(data_3$X3,data_3$y3,2,scaling=FALSE)
res.tripls <- tripls(data_3$X3,data_3$y3,1,qnames=paste("B",1:3,sep=""))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.