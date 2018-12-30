intervals.tripls <-
function(res.tripls,alpha=.05,Xnew){

tr <- function(X){return(sum(diag(X)))}  
optcomp <- attr(res.tripls,"ncomp") 
wJ <- res.tripls$WJ[,optcomp]
p <- length(wJ)
wK <- res.tripls$WK[,optcomp]
q <- length(wK)
X <- as.matrix(res.tripls$X.scaled)
y <- res.tripls$y.scaled
n <- length(y)
zb <- qnorm(1-alpha)
J <- dtrids(X,y,p,q,optcomp) 
bw <- coef(res.tripls)
r <- y - X%*%bw[,optcomp]
dYh <-  X%*%J%*%t(X) + matrix(1,n,n)/n
dr <- (diag(n) - dYh)
df <- (tr(t(dr)%*%dr))
sigma2 <- as.numeric(t(r)%*%r)/df
covb <- J%*%t(X)%*%X%*%t(J)
covb2 <- diag(covb*sigma2)
df <- (n+1-tr(dYh))

if(missing(Xnew)){
  type <- "Coefficients"
  intervals <- data.frame(llim=bw[,optcomp] -
                            zb*sqrt(covb2),coefficients=bw[,optcomp],ulim=bw[,optcomp]+zb*sqrt(covb2))  
} else {
  type <- "yyp"
  yp <- predict(res.tripls,Xnew)
  Xn <- attr(yp,"Xn.scaled")[,]
  sigma2 <- sigma2
  covb <- covb*sigma2
  var.yp <- (sigma2*(1+1/n) + Xn%*%covb%*%t(Xn))
  sd.yp <- sqrt(diag(var.yp))
  ulim <- yp + qt(1-alpha/2,df)*sd.yp*res.tripls$YScales 
  llim <- yp - qt(1-alpha/2,df)*sd.yp*res.tripls$YScales  
  intervals <- data.frame(llim=llim,y_predicted=yp,ulim=ulim)  
}
inputs <- list(alpha=alpha,optcomp=optcomp,type=type)
output <- list(intervals=intervals,df=df)
attr(output,"inputs") <- inputs
return(output)
}
