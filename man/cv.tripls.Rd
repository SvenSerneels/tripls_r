\name{cv.tripls}
\alias{cv.tripls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cross-validation for trilinear PLS regression
}
\description{Provides random subset k-fold cross validation for the tripls function. If plot=TRUE, a root mean squared error of prediction plot is drawn.
}
\usage{
cv.tripls(X, y, Amax, scaling = TRUE, fold = 10, repeats = 100, plot = TRUE, colors = list(line = "#0000AA", background = "#CCCCFF"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{The X data as a 3 way (n x p x q) array. 
}
  \item{y}{
The corresponding y data as an n x 1 vector or a length n numeric
}
  \item{Amax}{
The maximal number of latent components to estimate (integer)
}
  \item{scaling}{
Logical flag. If TRUE, the data are internally centered and scaled to unit variance. If FALSE, the data are only centered. 
}
  \item{fold}{The number of cases to leave out in each random cross-validation repeat (integer)
}
  \item{repeats}{
The number of times to repeat the cross validation process (integer)
}
  \item{plot}{
Logical flag. If TRUE, an RMSEP plot vs. the number of components is plotted. 
}
  \item{colors}{
A List object specifying the colors for the RMSEP plot. It contains the elements "lines" and "background". 
}
}
\details{
}
\value{
The function returns a list object containing
  \item{cvres}{A data frame containing a vector of the number of components (1:Amax), and the corresponding RMSEC and RMSEP values}
  \item{mincv}{The number of components that corresponds to the minimal RMSEP}
}
\references{
}
\author{
Sven Serneels, BASF Corp.
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
cv.tripls(data_3$X3, data_3$y3, 3, scaling = TRUE, fold = 4, repeats = 100, 
    plot = TRUE, colors = list(line = "#0000AA", background = "#CCCCFF")) }
