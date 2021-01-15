rowVars <- function(x, na.rm=FALSE, dims=1, unbiased=TRUE,
                    SumSquares=FALSE, twopass=FALSE) {
  if (SumSquares) return(rowSums(x^2, na.rm, dims))
  N <- rowSums(!is.na(x), FALSE, dims)
  Nm1 <- if (unbiased) N-1 else N
  if (twopass) {x <- if (dims==0) x - mean(x, na.rm=na.rm) else
    sweep(x, 1:dims, rowMeans(x,na.rm,dims))}
  (rowSums(x^2, na.rm, dims) - rowSums(x, na.rm, dims)^2/N) / Nm1
}

expit2 = function(x) 2^x/(1+2^x)

newgapfinder<-function(betamatrix=NULL,cutoffnum=5){

  NAwhenGap<-apply(betamatrix,1,function(x){

    originalbeta<-x
    perc25<-quantile(originalbeta,0.25,na.rm=TRUE)
    perc75<-quantile(originalbeta,0.75,na.rm=TRUE)
    gapsize<-abs(perc75-perc25)*3
    tempbeta<-sort(originalbeta)
    diffbetas<-diff(tempbeta)
    gapindexes <- which(diffbetas>=gapsize)
    numbergaps<-length(gapindexes)

    if(numbergaps>0){

      for(i in 1:numbergaps){
        tempind<-gapindexes[i]

        ## If a low outlier and number of outliers is less than the allowed cutoff number
        if(tempind<ncol(betamatrix)/2 & tempind <= cutoffnum){
          OutlierIDs<-names(tempbeta)[1:tempind]
          originalbeta[OutlierIDs]<-NA
        }

        ## If a high outlier and number of outliers is less than the allowed cutoff number
        if (tempind>ncol(betamatrix)/2 & tempind >= (ncol(betamatrix) - cutoffnum)){
          OutlierIDs<-names(tempbeta)[(tempind+1):ncol(betamatrix)]
          originalbeta[OutlierIDs]<-NA
        }
      }
    }
    return(originalbeta)
  })

  NAwhenGap<-t(NAwhenGap)
  colnames(NAwhenGap)<-colnames(betamatrix)
  rownames(NAwhenGap)<-rownames(betamatrix)

  return(NAwhenGap)

}
