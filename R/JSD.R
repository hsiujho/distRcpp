


JSD <- function(inMatrix, taxa_are_rows = TRUE, pseudocount=0.000001,numThreads=3, ...) {
  inMatrix=apply(inMatrix,1:2,function(x) ifelse(x==0,pseudocount,x))
  library(RcppParallel)

  setThreadOptions(numThreads = numThreads)

  if(taxa_are_rows){
    resultsMatrix <- rcpp_parallel_js_distance(t(inMatrix))
  } else {
    resultsMatrix <- rcpp_parallel_js_distance(inMatrix)
  }
  colnames(resultsMatrix)=rownames(resultsMatrix)=colnames(inMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix)
}
