
distRcpp_DEMO=function(){
  require(BacteriaIden)
  require(distRcpp)

  k0=subset_samples(MS16049_gg_phylo,group=="ST") %>>%
    (prune_taxa(taxa_sums(.)>0,.))

  #weighted unifrac distance, pow=1為weigthed unifrac, numThreads為執行緒數目
  aa=GUniFrac(k0,pow=1,numThreads=3)
  aa1=as.matrix(aa)
  dimnames(aa1)=NULL
  #與其他套件比較結果是否一致
  bb=MiSPU::GUniFrac(otu_table(k0)@.Data %>>% t(),phy_tree(k0))
  bb1=bb$dw
  dimnames(bb1)=NULL
  all.equal(aa1,bb1)

  #Jensen-Shannon distance

  otu=otu_table(k0)@.Data %>>% t() %>>% wisconsin() %>>% t()
  system.time(a3<-JSD(otu))

  dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
    KLD <- function(x,y) sum(x *log(x/y))
    JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
    matrixColSize <- length(colnames(inMatrix))
    matrixRowSize <- length(rownames(inMatrix))
    colnames <- colnames(inMatrix)
    resultsMatrix <- matrix(0, matrixColSize, matrixColSize)

    inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))

    for(i in 1:matrixColSize) {
      for(j in 1:matrixColSize) {
        resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                               as.vector(inMatrix[,j]))
      }
    }
    colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
    as.dist(resultsMatrix)->resultsMatrix
    attr(resultsMatrix, "method") <- "dist"
    return(resultsMatrix)
  }

  a4=dist.JSD(otu)
  all.equal(a3,a4)

}
