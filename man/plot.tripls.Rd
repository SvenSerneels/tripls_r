\name{plot.tripls}
\alias{plot.tripls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Graphical visualization for Tri-PLS
}
\description{This function provides three types of plots from Trilinear Partial Least Squares regression estimates. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.tripls(res.tripls, type = "yyp", alpha = 0.05, comp = NULL, colors = list(bars = "#0000AA", errorbars = "red", background = "#BBBBEE", abline = "#21A0D2"), textsize = 6, errorbar_width = 0.3, intervals, Xnew, ynew, yscale, nam, range = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res.tripls}{A "tripls" class regression object containing Tri-PLS regression estimates
%%     ~~Describe \code{res.tripls} here~~
}
  \item{type}{A flag (character) that determines the plot type. If set to "weights", plot.tripls yields a level plot of the matrix of weighting vectors, showing how individual variables contribute to latent variables. In that case, the input variable comp controls which latent variable's weighting matrix is being plotted. If set to "coefficients", the individual regression coefficients, as well as their confidence intervals, will be plotted. If set to "yyp", the predicted responses will be plotted versus the original responses (if known), together with their prediction intervals. In the "yyp" case, if no new X cases are supplied through the variable "Xnew", the y vs. y predicted plot will be shown for the cases used to construct the Tri-PLS model. Otherwise, the predicted responses and their prediction intervals will be plotted for the new cases (in Xnew), versus either their index (if the responses are unknown) or versus the corresponding true responses (to be provided in ynew). In both cases, alpha controls the significance level. 
%%     ~~Describe \code{type} here~~
}
  \item{alpha}{Optional positive numeric. Sinificance level at which to compute confidence intervals and prediction intervals. If unspecified, set to 0.05. 
%%     ~~Describe \code{alpha} here~~
}
  \item{comp}{Optional integer. The number of tri-PLS components at which to plot the weighting matrix plot. If unspecified, it is set to the maximal number of components at which res.tripls was computed. Only has an effect if type="weights". 
%%     ~~Describe \code{comp} here~~
}
  \item{colors}{An optional list object containing the colors in which to plot bars (bars), errorbars (errorbars), background (background) and the diagonal line in a y vs. y predicted plot (abline). 
%%     ~~Describe \code{colors} here~~
}
  \item{textsize}{An optional positive numeric value that specifies the font size for individual case names in the plot. Only relevant in case type="coefficients" or if type="yyp". 
%%     ~~Describe \code{textsize} here~~
}
  \item{errorbar_width}{An optional positive numeric that specifies the errorbar width. Only relevant in case type="coefficients" or if type="yyp". 
%%     ~~Describe \code{errorbar_width} here~~
}
  \item{intervals}{An optional data frame that contains confidence or prediction intrevals, returned from intervals.tripls (if data dimensions are big, it may be useful to store these results). Only relevant in case type="coefficients" or if type="yyp". 
%%     ~~Describe \code{intervals} here~~
}
  \item{Xnew}{An optional three-way array (n' x p x q) that contains a set of new cases for which to plot the predicted vs. the actual responses. Only relevant in case type="yyp". 
%%     ~~Describe \code{Xnew} here~~
}
  \item{ynew}{An optional length n' numeric or n' x 1 matrix or data frame containing the responses that correspond to Xnew. 
%%     ~~Describe \code{ynew} here~~
}
  \item{yscale}{An optional length 2 numeric specifying the y scale in the plots (useful if several models have to be plotted on a comparable scale). 
%%     ~~Describe \code{yscale} here~~
}
  \item{nam}{An optional lenth n or n' character array containing the names for the individual cases. 
%%     ~~Describe \code{nam} here~~
}
  \item{range}{An optional input that specifies the fraction of the output to be plotted. If either the number of cases or the number of regression coefficients is big, the range can be specified to e.g. 1:100, in which case only the first 100 elelments will be plotted. 
%%     ~~Describe \code{range} here~~
}
}
\details{Requires ggplot2. 
%%  ~~ If necessary, more details than the description above ~~
}
\value{Does not return output. 
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
plot(res.tripls, type = "weights")
plot(res.tripls, type = "weights",comp=1)

plot(res.tripls, type = "coefficients")
plot(res.tripls, type = "coefficients",alpha=.025)
plot(res.tripls, type = "coefficients",range=1:10)
intervals <- intervals.tripls(tripls(data_3$X3,data_3$y3,1),alpha=.05)
plot(res.tripls, type = "coefficients", alpha = 0.05, colors = list(bars = "#0000AA", errorbars = "red", background = "#BBBBEE", abline = "#21A0D2"), errorbar_width = 0.3, intervals=intervals$intervals, range = 1:10)

plot(res.tripls, type = "yyp")
plot(res.tripls, type = "yyp",alpha=.025)
plot(res.tripls, type = "yyp",range=1:10)
plot(res.tripls, type = "yyp",Xnew=data_3$X3)
plot(res.tripls, type = "yyp",Xnew=data_3$X3, ynew=data_3$y3, nam=paste("C",18:34,sep=""))
plot(res.tripls, type = "yyp",Xnew=data_3$X3, ynew=data_3$y3, range=1:10)
intervals <- intervals.tripls(tripls(data_3$X3,data_3$y3,1),alpha=.05,data_3$X3)
plot(res.tripls, type = "yyp", alpha = 0.05, colors = list(bars = "#0000AA", errorbars = "red", background = "#BBBBEE", abline = "#21A0D2"), textsize = 6, errorbar_width = 0.3, intervals=intervals, Xnew=data_3$X3, ynew=data_3$y3, yscale=c(50,90), nam=paste("C",18:34,sep=""), range = 1:7)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
