tripls <-
function(X,y,A,scaling=TRUE,nnames=NULL,pnames=NULL,qnames=NULL){
# Adapted from S. de Jong, Journal of Chemometrics, 12 (1998) 77-81.
# Alternative deflation by S. Serneels, November 2010. 
# Converted into R by S. Serneels, 1.2012

# Initialization   
DX <- dim(X);
I <- DX[1]
J <- DX[2]
K <- DX[3]
if(is.na(K)){stop("The X data have to be presented as a 3 way array.")}
T <- matrix(0,DX[1],A) ; 
W <- matrix(0,J*K,A); 
WJ <- matrix(0,J,A); 
WK <- matrix(0,K,A); 
bNPLS <- matrix(0,J*K,A);
if(is.null(nnames)){nnames <- 1:I}
if(is.null(pnames)){pnames <- 1:J}
if(is.null(qnames)){qnames <- 1:K}

# Centering and scaling
X0 <- X
X <- as.matrix(as.data.frame(X))
centering <- TRUE
X <- scale(X,center=centering,scale=scaling)
mX <- attr(X,"scaled:center")
ifelse(scaling,sX <- attr(X,"scaled:scale"),sX <- rep(1,J*K))
DY <- dim(y)
if(is.null(DY[2])){
  n <- length(y)
  if(!n==I){stop("Number of cases in y must equal number of cases in X")}
  y <- as.matrix(y,nrow=n,ncol=1)
} else {
  if(DY[2]>1){stop("This routine only works for univariate y")}
  n <- DY[1]
  if(!n==I){stop("Number of cases in y must equal number of cases in X")}
}
y0 <- y
y <- scale(y,center=centering,scale=scaling)
my <- attr(y,"scaled:center")
ifelse(scaling,sy <- attr(y,"scaled:scale"),sy <- 1)

s0 <- t(X)%*%y;
s <- s0;                        # initialize
for(a in 1:A){                   # for each factor
   Z <- array(t(s),dim=c(J,K));    # vector of covariances -> matrix
   wJK <- svd(Z);        # find weights
   w <- kronecker(wJK$v[,1],wJK$u[,1]); # maximizing covariance
   W[,a] <- w;                 # save weights
   WJ[,a] <- wJK$u[,1];         # J-mode
   WK[,a] <- wJK$v[,1];         # K-mode
   T[,a] <- X%*%w;             # save scores I-mode
   qrT <- qr(T[,1:a]);
   b <- qr.coef(qrT,y);                 # y loadings wrt T
   bNPLS[,a] <- W[,1:a] %*% b;        # regression coefficients
   s <- s0 - t(X)%*%T[,1:a]%*%b;             # residual s
}

b0 <- mean(y-X%*%bNPLS[,A])
yfit <- (X %*% bNPLS[,A] + b0)*sy + my

sev <- 1 - sum((s0 - t(X)%*%T[,1:a]%*%b)^2)/sum((s0)^2)
rmsec <-  sqrt(mean((y0-yfit)^2))

cnames <- paste("Comp",1:A)
colnames(T) <- cnames
colnames(bNPLS) <- cnames
colnames(WJ) <- cnames
colnames(WK) <- cnames
colnames(W) <- cnames
rownames(yfit) <- nnames
rownames(T) <- nnames 
rownames(y) <- nnames  
rownames(WJ) <- pnames
rownames(WK) <- qnames
Cnames <- NULL
for (i in 1:K){
  Cnames <- c(Cnames,paste(qnames[i],"_",pnames,sep=""))
}
rownames(W) <- Cnames
rownames(bNPLS) <- Cnames 
names(mX) <- Cnames
names(sX) <- Cnames
X.scaled <- as.data.frame(X) 
colnames(X.scaled) <- Cnames
rownames(X.scaled) <- nnames

inputs <- list(X0=X0,y0=y0)
output <- list(coefficients=bNPLS,intercept=b0,scores=T,fitted.values=yfit,W=W,WJ=WJ,WK=WK,YMeans=my,YScales=sy,XMeans=mX,XScales=sX,X.scaled=X.scaled,y.scaled=y,sev=sev,rmsec=rmsec,inputs=inputs)
class(output) <- "tripls"
attr(output,"ncomp") <- A
return(output)
}
