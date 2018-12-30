predict.tripls <-
function(res.tripls,data,casenames=NULL){
  
# PREDICT.TRIPLS predicts response based on a trilinear partial least squares model 
#     Inputs: res.tripls, a "tripls" class regression object 
#             data, an optional data matrix or frame containing a new set of cases 
#
# Written by Sven Serneels, BASF Corp, December 2013. 
  
  b <- coef(res.tripls)
  db <- dim(b)
  b <- b[,db[2]]
  if(missing(data)){fitted.values <- res.tripls$fitted.values}
  else{
    DD <- dim(data)
    if(DD[2]*DD[3]==length(b)){Xn <- as.matrix(as.data.frame(data))} 
    else{stop("New caswes must be of the same dimensionality as original cases")}  
    Xn <- scale(Xn,center=res.tripls$XMeans,scale=res.tripls$XScales)
    fitted.values <- (Xn%*%b + res.tripls$intercept)*res.tripls$YScales + res.tripls$YMeans
    Cnames <- names(res.tripls$XMeans)
    if(!is.null(casenames)){
      rownames(fitted.values) <- casenames
      rownames(Xn) <- casenames}
    colnames(Xn) <- Cnames
    attr(fitted.values,"Xn.scaled")<- Xn
  }
  return(fitted.values)
}
