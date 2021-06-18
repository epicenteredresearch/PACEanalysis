#'Detection Masking
#'
#'@description Masks beta-values based on poor intensity values
#'@param processedBetas The preprocessed beta-value matrix
#'@param DetectionPvals An optional matrix of detection p-values
#'@param DetectionPvalCutoff a numeric value indicating the detection p-value
#'  cut-off (default=0.05)
#'@param IndicatorGoodIntensity An optional matrix used to mask methylation
#'  values with a poor intensity value based on the number of beads and/or
#'  intensity values of zero; NA if poor intensity value and 1 otherwise. Output
#'  from ExploratoryDataAnalysis function
#'@param destinationfolder A character string indicating the location where
#'  files should be saved, e.g. "C:\\Home\\PACE\\BirthSize"
#'@param cohort A character string for the cohort's acronym (e.g. "HEBC")
#'@param analysisdate A character string indicating the date the analysis was
#'  run. Please specify in the form: YEARMONTHDAY, e.g. "20200205" for February
#'  5th 2020
#'@details  Masks beta-values based on the detection p-value; the default
#'  p-value cut-off is 0.05.
#'@return Returns and saves an RData file of the preprocessed beta-value matrix
#'  with beta-values masked based on poor intensity values. Also saves a csv of
#'  the CpG distributions after this masking
#'@examples
#'\dontrun{
#'betasabovedetection<-detectionMask(processedBetas=processedOut$processedBetas,
#'                                  DetectionPvals=EDAresults$DetectionPval,
#'                                  DetectionPvalCutoff=0.05,
#'                                  IndicatorGoodIntensity=EDAresults$IndicatorGoodIntensity,
#'                                  destinationfolder="H:/UCLA/PACE/Birthweight-placenta",
#'                                  cohort="HEBC",analysisdate="20200710")
#'}

detectionMask<-function(processedBetas=NULL,DetectionPvals=NULL,DetectionPvalCutoff=0.05,
                        IndicatorGoodIntensity=NULL,cohort=NULL,analysisdate=NULL,destinationfolder=NULL){

  if(!is.matrix(processedBetas)) stop("Please make sure processedBetas is a matrix")
  if(is.null(destinationfolder)) stop("Please specify a destination folder")

  ## If masking based on detection p-value
  if(!is.null(DetectionPvals)){

    if(!is.matrix(DetectionPvals)) stop("Please make sure DetectionPvals is a matrix")
    if(!is.numeric(DetectionPvalCutoff)) stop("DetectionPvalCutoff must be numeric")

    DetectionPvals[DetectionPvals>DetectionPvalCutoff]<-NA
    DetectionPvals[!is.na(DetectionPvals)]<-1

    ## first restricting to the same probes and samples
    DetectionPvals<-DetectionPvals[rownames(processedBetas),]
    DetectionPvals<-DetectionPvals[,colnames(processedBetas)]

    processedBetas<-processedBetas*DetectionPvals

  }

  ## If masking based on N beads and zero intensity value
  if(!is.null(IndicatorGoodIntensity)){

    if(!is.matrix(IndicatorGoodIntensity)) stop("Please make sure IndicatorGoodIntensity is a matrix")

    ## first restricting to the same probes and samples
    IndicatorGoodIntensity<-IndicatorGoodIntensity[rownames(processedBetas),]
    IndicatorGoodIntensity<-IndicatorGoodIntensity[,colnames(processedBetas)]

    processedBetas<-processedBetas*IndicatorGoodIntensity

  }

  summarybetas<-apply(processedBetas,1,function(x){
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
  write.csv(summarybetas,paste(cohort,"_",analysisdate,"_Summarize_Beta_Values_detectionmasked.csv",sep=""))

  save(file=paste(cohort,"_",analysisdate,"_PreprocessedBetas_detectionmasked.RData",sep=""),compress=TRUE, list=c("processedBetas"))
  return(processedBetas)

}
