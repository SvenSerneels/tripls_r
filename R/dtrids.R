dtrids <-
function(X,y,p,q,h){

# Sven Serneels, 3. 2012.  


library(MASS);
DX <- dim(X);
n <- DX[1];
T <- matrix(0,DX[1],h) ; 
W <- matrix(0,p*q,h); 
WJ <- matrix(0,p,h); 
WK <- matrix(0,q,h); 
bNPLS <- matrix(0,p*q,h);
s0 <- t(X)%*%y;
s <- s0;                        # initialize
dbds <- matrix(0,h,p*q);
dTds <- matrix(0,h*n,p*q);
dWds <- matrix(0,h*p*q,p*q);
dbNPLSds <- matrix(0,p*q,p*q);
e <- y;

for(i in 1:h){
    Z <- array(t(s),dim=c(p,q));
    if(i==1){
        dZds <- diag(1,p*q);}
    else{
        dZds <- diag(1,p*q) - kronecker(t(b[1:(i-1)]),t(X))%*%dTds[1:((i-1)*n),] - t(X)%*%T[,1:(i-1)]%*%(dbds);
    }
    wJK <- svd(Z);   
    sigma <- wJK$d[1]; 
    wJ <- wJK$u[,1]; 
    wK <- wJK$v[,1];
    dwJds <- (kronecker(t(wJ)%*%Z,ginv(sigma^2*diag(1,p)-Z%*%t(Z),tol=1e-10))
             +kronecker(ginv(sigma^2*diag(1,p)-Z%*%t(Z),tol=1e-10)%*%Z,t(wJ)))%*%dZds;
    WJ[,i] <- wJ;
    dwKds <- (kronecker(ginv(sigma^2*diag(1,q)-t(Z)%*%Z,tol=1e-10),t(wK)%*%t(Z))+
              kronecker(t(wK),ginv(sigma^2*diag(1,q)-t(Z)%*%Z,tol=1e-10)%*%t(Z)))%*%dZds;
    WK[,i] <- wK;
    w <- kronecker(wK,wJ);
    W[,i] <- w;
    dwds <- kronecker(diag(1,q),wJ)%*%dwKds+kronecker(wK,diag(1,p))%*%dwJds;
    dWds[((i-1)*p*q+1):(i*p*q),] <- dwds;
    t <- X%*%w;
    T[,i] <- t;
    dTds[((i-1)*n+1):(i*n),] <- X%*%dwds;
    qrT <- qr(T[,1:i]);
    b <- qr.coef(qrT,y);   
    TTT <- ginv(t(T[,1:i])%*%T[,1:i],tol=1e-10);
    dbds1 <- TTT%*%t(W[,1:i]); 
    dbds2 <- kronecker(TTT,t(s0))%*%dWds[1:(i*p*q),];
    dbds3a <- kronecker(t(s0)%*%W[,1:i]%*%TTT,TTT%*%t(T[,1:i]))
    dbds3b <- kronecker(TTT,t(s0)%*%W[,1:i]%*%TTT%*%t(T[,1:i]))
    dbds <- dbds1 + dbds2 - (dbds3a + dbds3b)%*%dTds[1:(i*n),];
    bNPLS[,i] <- W[,1:i] %*% b;
    dbNPLSds <- dbNPLSds + matrix(b[i],(p*q),(p*q))*dwds + w%*%t(dbds[i,]);
    s <- s0 - t(X) %*% T[,1:i] %*% b;
    e <- y - T[,1:i]%*%b;
} #for
return(dbNPLSds)
}
