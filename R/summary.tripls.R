summary.tripls <-
function(res.tripls){
  
  if(!(class(res.tripls)=="tripls")){stop("The tripls summary function can only be applied to tripls class objects")}
  
  
  cat(" ------------------------------------------","\n",
      "Trilinear Partial Least Squares Regression","\n",
      "------------------------------------------","\n",
      c("# Components:",attr(res.tripls,"ncomp")),"\n",
      c("Explained covariance:", round(res.tripls$sev*100,2),"%"),"\n",
      c("RMSEC:",round(res.tripls$rmsec,2)),"\n",
      "------------------------------------------")
  
}
