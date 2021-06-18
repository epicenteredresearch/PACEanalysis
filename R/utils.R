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

outliermethod <- function(betamatrix=NULL,quantilemethod=NULL,trimming=FALSE,pct=0.005,NumberOutliers=5) {

  ## If trimming, using the IQR; percentile is 0.25
  if(trimming) pct=0.25
  betamatrixoriginal<-betamatrix

  if(quantilemethod=="EmpiricalBeta"){

    ## Using Method of Moments approach; faster than MLE
    sample.means <- rowMeans(betamatrix, na.rm=T)
    sample.vars <- matrixStats::rowVars(betamatrix, na.rm=T)

    v <- sample.means * (1 - sample.means)

    # initialize
    shape1<-rep(NA,length(v))
    shape2<-rep(NA,length(v))

    ## less than v
    ind1 <- which(sample.means < v , arr.ind=T)
    shape2[ind1] <- sample.means[ind1] * (v[ind1]/sample.vars[ind1] - 1)
    shape1[ind1]  <- (1 - sample.means[ind1]) * (v[ind1]/sample.vars[ind1] - 1)

    ## greater than v; switch shape parameters
    ind2 <- which(sample.means >= v, arr.ind=T)
    shape1[ind2] <- sample.means[ind2] * (v[ind2]/sample.vars[ind2] - 1)
    shape2[ind2] <- (1 - sample.means[ind2]) * (v[ind2]/sample.vars[ind2] - 1)

    low<-qbeta(pct,shape1,shape2,lower.tail=TRUE)
    upper<-qbeta((1-pct),shape1,shape2,lower.tail=TRUE)
    quantiles<-cbind(low,upper)

  }

  if(quantilemethod=="Quantile"){
    quantiles <- matrixStats::rowQuantiles(betamatrix, probs=c(pct,(1-pct)), na.rm=T)
  }

  ## Trimming is based on 3*IQR
  if(trimming){

    IQRtemp<-quantiles[,2]-quantiles[,1]
    low <- quantiles[,1]-3*IQRtemp
    upper <- quantiles[,2]+3*IQRtemp

    outliers.lower <- rowSums(betamatrix < low, na.rm=T)
    outliers.upper <- rowSums(betamatrix > upper, na.rm=T)

    betamatrix[which(betamatrix <= low, arr.ind=T)] <- NA
    betamatrix[which(betamatrix >= upper, arr.ind=T)] <- NA

  ## If we don't trim, we winsorize
  } else {

    low <- quantiles[,1]
    upper <- quantiles[,2]

    outliers.lower <- rowSums(betamatrix < low, na.rm=T)
    outliers.upper <- rowSums(betamatrix > upper, na.rm=T)

    idx <- which(betamatrix < low, arr.ind=T)
    betamatrix[idx] <- low[idx[,1]]

    idx <- which(betamatrix > upper, arr.ind=T)
    betamatrix[idx] <- upper[idx[,1]]

  }

  ## If more outliers than a certain specified number (NumberOutliers)
  ## assuming they are not outliers, but reflective of the data distribution
  idx <- which(outliers.lower >= NumberOutliers)
  betamatrix[idx,]<-betamatrixoriginal[idx,]

  idx <- which(outliers.upper >= NumberOutliers)
  betamatrix[idx,]<-betamatrixoriginal[idx,]

  n <- rowSums(!is.na(betamatrix))
  log <- data.frame(outliers.lower, outliers.upper, n)

  return(list(methylation=betamatrix, log=log))

}



