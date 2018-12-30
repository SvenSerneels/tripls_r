cv.tripls <-
function(X,y,Amax,scaling=TRUE,fold=10,repeats=100,plot=TRUE,colors=list(line="#0000AA",background="#CCCCFF")){
  
  DX <- dim(X) 
  n <- length(y)
  if(!(n==DX[1])){stop("Number of cases in X and y must match")}
  set <- 1:n
  rc <- array(dim=c(repeats,Amax+1))
  rv <- array(dim=c(repeats,Amax+1))
  for(i in 1:repeats){
    newset <- sample(set)
    XCV <- X[newset,,]
    yCV <- y[newset]
    calsamps <- ((fold+1):n)
    valsamps <- 1:fold
    XC <- XCV[calsamps,,]
    XV <- XCV[valsamps,,]
    yc <- yCV[calsamps] 
    yv <- yCV[valsamps]
    for(j in 0:Amax){
      if(j==0){
        ycp <- mean(yc)
        rc[i,j+1] <- sum((yc-ycp)^2)/(n-fold)
        yvp <- mean(yv)
        rv[i,j+1] <- sum((yv-yvp)^2)/fold
      } else {
        cvres <- tripls(XC,yc,j,scaling)
        rc[i,j+1] <- cvres$rmsec^2 
        yvp <- predict(cvres,XV)
        rv[i,j+1] <- sum((yv-yvp)^2)/fold
      }
    }
  }  
rc <- sqrt(apply(rc,2,mean))
rv <- sqrt(apply(rv,2,mean))
optcomp <- which(rv==min(rv))-1  
plotcv <- data.frame(RMSEC=rc,RMSEP=rv,n=0:Amax)  
if(plot==TRUE){
  require(ggplot2)
  print(ggplot(plotcv,aes(n,RMSEP))+geom_line(size=1.5,color=colors$line)+labs(x="Number of components",y="RMSEP",title=paste("TriPLS",names(y)," leave ",fold," out random cross validation",sep="")) +
          geom_text(aes(x=4,y=1),label=paste("Minimum = ",optcomp),size=7,color=colors$line) +
          theme(panel.background=element_rect(fill=colors$background),title=element_text(size=rel(1.5),face="bold")))
}
return(list(cvres=plotcv,mincv=optcomp))  
}
