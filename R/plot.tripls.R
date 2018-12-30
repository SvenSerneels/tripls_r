plot.tripls <-
function(res.tripls,type="yyp",alpha=.05,comp=NULL,colors=list(bars="#0000AA",errorbars="red",background="#BBBBEE",abline="#21A0D2"),textsize=6,errorbar_width=.3,intervals,Xnew,ynew,yscale,nam,range=NULL)
  
  # plot.tripls automatically generates Trilinear Partial Least Squares regression Weights and Regression Coefficients Plots
  # inputs : res.tripls, a "tripls" class regression object 
  #          type, either "weights", "coefficients" or "yyp" 
  #          colors: a list containing the desired bar, errorbars, diagnal line and background colors as strings
  #          textsize: the text size in which to print the scores and loading names
  #          errorbar_width: a numeric containing the width of the errorbars on the regcoefficients  
  #          data (optional) a data frame containing new cases to predict and plot
  #          yscale, an optional scale vector for the yscale in the y vs y predicted plot (e.g. if two different regression plots have to be on the same scale)
  
  # written by Sven serneels, BASF Corp., April 2014. 

{
  require(ggplot2)
  require(grid)
  require(reshape)
  
if(!(class(res.tripls)=="tripls")){stop("The tripls plot function can only be applied to tripls class objects")}
     optcomp <- attr(res.tripls,"ncomp") 
     if(is.null(comp)){comp=optcomp}
     w <- res.tripls$W[,comp]
     wJ <- res.tripls$WJ[,comp]
     p <- length(wJ)
     wK <- res.tripls$WK[,comp]
     q <- length(wK)
  yp <- predict(res.tripls)
  
if(type=="weights"){
  Wmat <- array(w,dim=c(p,q))
  rownames(Wmat) <- names(wJ)
  colnames(Wmat) <- names(wK)
  plotty <- ggplot(melt(Wmat),aes(X1,X2,fill=value))
  plotty <- plotty + geom_tile()
  plotty <- plotty + labs(title=paste("Tri-PLS component",comp,"weighting matrix plot"))
  plotty <- plotty + theme(panel.background=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_text(angle=-90,face="bold"),
                           axis.text.y=element_text(face="bold"),plot.title=element_text(face="bold"),panel.grid=element_blank())
  print(plotty)
} else { 
  if(type=="coefficients"){
    if(is.null(range)){range <- 1:(p*q)}
    if(missing(intervals)){plotcoefs <- intervals.tripls(res.tripls,alpha)$intervals}
    else {
      DI <- dim(intervals)
      if(!DI[1]==(p*q)){stop("The set of supplied intervals' dimernsionality does not match")}
      plotcoefs <- intervals
    }
    plotcoefs$nam <- rownames(plotcoefs)
    plotty <- ggplot(plotcoefs[range,],aes(nam,coefficients)) + geom_bar(stat="identity",size=3,fill=colors$bars) + labs(title=paste(names(res.tripls$YMeans)," Trilinear PLS Regression Coefficients Plot"))
    plotty <- plotty + coord_flip() + geom_errorbar(aes(ymin=llim,ymax=ulim),width=errorbar_width,color=colors$errorbars)
    plotty <- plotty + theme(panel.background=element_rect(fill=colors$background),plot.title=element_text(size=rel(1),face="bold"),
                axis.text.x=element_text(angle=-90),axis.title.x=element_blank(),axis.title.y=element_blank(),panel.grid=element_blank())
    print(plotty)
    } else {
      if(type=="yyp"){
        ynames <- names(res.tripls$YMeans)
      if(missing(Xnew)){
        Xnew <- res.tripls$inputs$X0
        ynew <- res.tripls$inputs$y0
        if(missing(nam)){nam <- rownames(res.tripls$y.scaled)}
      }
      yp <- predict(res.tripls,Xnew)
      n <- length(yp) 
        if(is.null(range)){range <- 1:n}  
      ifelse(missing(ynew),y0 <- array(1:n,c(n,1)),y0 <- ynew )
      if(missing(nam)){nam <- as.character(1:n)}
          Xn <- attr(yp,"Xn.scaled")[,]
        if(missing(intervals)){plotyyp <- intervals.tripls(res.tripls,alpha,Xnew)}
        else {
          plotyyp <- intervals
        }
      df <- plotyyp$df 
      plotyyp <- plotyyp$intervals
      DI <- dim(plotyyp)
      if(!DI[1]==n){stop("The set of supplied intervals' dimernsionality does not match")}  
      plotyyp$nam <- nam
      plotyyp$y_original <- as.numeric(y0[,])
      rsq <- 1 - sum((yp-plotyyp$y_original)^2)/((n-1)*var(plotyyp$y_original))
      plotty <- ggplot(plotyyp[range,])
      plotty <- plotty + geom_text(aes(x=y_original,y=y_predicted,label=nam),size=textsize,color=colors$bars) 
      plotty <- plotty + geom_point(aes(x=y_original,y=y_predicted,label=nam),size=textsize,color=NA) 
      plotty <- plotty + geom_errorbar(aes(x=y_original,y=y_predicted,ymin=llim,ymax=ulim),width=errorbar_width,color=colors$errorbars)
      plotty <- plotty + labs(title=paste("Trilinear PLS Regression ", ynames," vs. ", ynames, " Predicted -", optcomp," LV,",ifelse(y0==1:n," ",paste("R? =",round(rsq,2),",")), "df = ",round(df,1))) 
      plotty <- plotty + geom_abline(intercept=0,slope=1,color=colors$abline) 
      if(!missing(yscale)){plotty <- plotty + ylim(yscale)}
      plotty <- plotty + theme(panel.background=element_rect(fill=colors$background),plot.title=element_text(size=rel(1)),
                           axis.text.x=element_text(),axis.title.x=element_blank(),axis.title.y=element_blank())
      print(plotty)
      } else {stop("Plot type can be 'weights', 'coefficients', or 'yyp'")}
    }
}}
