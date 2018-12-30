\name{data_3}
\alias{data_3}
\docType{data}
\title{Three-way data set
%%   ~~ data name/kind ... ~~
}
\description{A generic three-way data set to test the algorithms on. It follows a 1 component tri-PLS model if scaled to unit variance. 
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(data_3)}
\format{The data are presented as a list object containing
X3 (Predictor data of dimensionality 17 x 11 x 3), y3 (Predictand of dimensionality 17 x 1),
pnames (Character arry of length 11 containin p mode names)
qnames (Character arry of length 3 containin q mode names).}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(data_3)
X3 <- data_3$X3
y3 <- data_3$y3
res.tripls <- tripls(X3,y3,2)
}