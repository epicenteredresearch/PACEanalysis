#'Correcting Outliers
#'
#'@description Correcting outliers from the beta-value matrix based on specified
#'  approach
#'@param processedBetas The preprocessed beta-value matrix
#'@param quantilemethod A character string indicating the approach to estimate
#'  quantiles. Two options: 'EmpiricalBeta' or 'Quantile'. To estimate the
#'  quantiles based on the empirical beta distribution, specify 'EmpiricalBeta'.
#'  To estimate the quantiles using the quantiles function (type7), specify
#'  'Quantile'.
#'@param trimming Logic statement (TRUE/FALSE; default is FALSE). If FALSE,
#'  winsorize the outliers based on the specified percentile. If TRUE, trims the
#'  outliers based on gaps larger than 3*IQR; the cutoff for the number of
#'  outliers in a group is either a maximum of 5 or 0.0025 of the total number
#'  of samples (whichever is larger). Detected outliers are recoded as NA.
#'@param pct Numeric value indicating the percentile winsorized on each side of
#'  the distribution; default is 0.005, indicating winsorizing 1 percent.
#'@param destinationfolder A character string indicating the location where
#'  files should be saved, e.g. "C:\\Home\\PACE\\BirthSize"
#'@param cohort A character string for the cohort's acronym (e.g. "HEBC")
#'@param analysisdate A character string indicating the date the analysis was
#'  run. Please specify in the form: YEARMONTHDAY, e.g. "20200205" for February
#'  5th 2020
#'@details  If estimating the quantiles based on the empirical beta
#'  distribution, the shape parameters are estimated using the method of moments
#'  approach. If trimming is TRUE, remove extreme outliers based on gaps that
#'  must be at least 3*IQR; the cutoff for the number of outliers in a group is
#'  either a maximum of 5 or 0.0025 of the total number of samples (whichever is
#'  larger). Detected outliers are recoded as NA.
#'@return Returns and saves an RData file of the preprocessed beta-value matrix
#'  with the outliers removed or corrected. Also saves a csv of the CpG
#'  distributions after removing or correcting the outliers, a csv of the number
#'  of outliers per CpG locus, and a csv of the proportion of outliers per
#'  sample.
#'@examples
#'\dontrun{
#'Betasnooutliers<-outlierprocess(processedBetas=betasabovedetection,
#'                                  quantilemethod="EmpiricalBeta",
#'                                  trimming=FALSE,
#'                                  pct=0.005,
#'                                  destinationfolder="H:/UCLA/PACE/Birthweight-placenta",
#'                                  cohort="HEBC",analysisdate="20200710")
#'}

outlierprocess<-function(processedBetas=NULL,quantilemethod=NULL,trimming=FALSE,pct=0.005,
                         cohort=NULL,analysisdate=NULL,destinationfolder=NULL){

  if(!is.matrix(processedBetas)) stop("Please make sure processedBetas is a matrix")
  if(!(quantilemethod %in% c("EmpiricalBeta","Quantile"))) stop("Quantile method must be  'EmpiricalBeta' or 'Quantile'")
  if(!is.numeric(pct)) stop("pct must be numeric")
  if(is.null(destinationfolder)) stop("Please specify a destination folder")

  Outputname<-paste(cohort,"_",analysisdate,"_Output",sep="")
  dir.create(file.path(destinationfolder, Outputname), showWarnings = FALSE)
  destinationfolder<-paste(destinationfolder,Outputname,sep="/")
  setwd(destinationfolder)

  #################################################################################
  # R-code for detecting outliers

  cat("Detecting beta-value outliers...","\n")
  # Required input & parameters
  N <- ncol(processedBetas) ## No. of samples in your data.
  cutoff <- ifelse(0.0025*N<5,5,ceiling(N*0.0025))
  ## This cutoff is chosen for detecting probes with outliers. We have chosen this
  ## cutoff such that a probe can have a maximum of 5 or 0.0025 of the total number
  ## of samples (whichever is larger) as outliers. Can change it if required.

  betafinal.nooutlier<-outliermethod(betamatrix=processedBetas,quantilemethod=quantilemethod,trimming=trimming,pct=pct,NumberOutliers=cutoff)
  write.csv(betafinal.nooutlier$log,paste(cohort,"_",analysisdate,"_OutliersSummary_CpGs.csv",sep=""))

  betafinal.nooutlier<-betafinal.nooutlier$methylation
  betafinal.nooutlier.sampleperc<-apply(betafinal.nooutlier,2,function(x) (sum(is.na(x))/nrow(betafinal.nooutlier))*100)
  write.csv(betafinal.nooutlier.sampleperc,paste(cohort,"_",analysisdate,"_OutliersSummary_Individuals.csv",sep=""))

  summarybetas<-apply(betafinal.nooutlier,1,function(x){
    tempx<-na.omit(x)
    Ntemp<-length(tempx)
    Nmisstemp<-length(x)-Ntemp
    Mintemp<-min(tempx)
    Maxtemp<-max(tempx)
    Mediantemp<-median(tempx)
    Meantemp<-mean(tempx)
    sdtemp<-sd(tempx)
    c(n=Ntemp,Nmissing=Nmisstemp,min=Mintemp,max=Maxtemp,median=Mediantemp,mean=Meantemp,sd=sdtemp)
  })

  summarybetas<-t(summarybetas)
  write.csv(summarybetas,paste(cohort,"_",analysisdate,"_Summarize_Beta_Values_nooutliers.csv",sep=""))

  save(file=paste(cohort,"_",analysisdate,"_PreprocessedBetas_nooutliers.RData",sep=""),compress=TRUE, list=c("betafinal.nooutlier"))
  return(betafinal.nooutlier)

}
