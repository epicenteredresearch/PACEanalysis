#'Load Data and Rename Phenotype Data to Standardize
#'
#'@description Reorganize data to standardize for analysis
#'@param SamplePlacement An data.frame that includes 'Basename' and sample ID.
#'  Not required if 'Basename'is in PhenoData; might also include an indicator
#'  for batch.
#'@param PhenoData A required data.frame that includes the participant
#'  characteristics and a sample ID
#'@param IDlink A required character string indicating the column name for the
#'  sample ID in the data.frame PhenoData
#'@param BWTvar An optional character string indicating the column name for the
#'  birthweight in the data.frame PhenoData; if there is no birthweight variable
#'  in data.frame PhenoData, can specify argument value is NULL. Birthweight
#'  values should be in grams
#'@param BATCHvar An optional character string indicating the column name for
#'  the batch variable in the data.frame PhenoData; if there is no batch
#'  variable in the data.frame PhenoData, can specify argument value is NULL.
#'  Batch variable is expected to be categorical
#'@param SEXvar A required character string indicating the column name for the
#'  infant in the data.frame PhenoData; assumed only two sexes are in the
#'  dataset
#'@param FemaleInd A character string indicating how female sex is coded in the
#'  dataset; e.g. "1"
#'@param MaleInd A character string indicating how male sex is coded in the
#'  dataset; e.g. "2"
#'@param ETHNICvar An optional character string indicating the column name for
#'  the maternal race/ethnicity variable in the data.frame PhenoData; if there
#'  is no maternal race/ethnicity variable in the data.frame PhenoData, can
#'  specify argument value is NULL. Maternal race/ethnicity variable is expected
#'  to be categorical
#'@param GESTvar An optional character string indicating the column name for the
#'  gestational age variable in the data.frame PhenoData; if there is no
#'  gestational age variable in the data.frame PhenoData, can specify argument
#'  value is NULL. Gestational age variable is expected to be in days
#'@param BIRTHLENGTHvar An optional character string indicating the column name
#'  for the birth length variable in the data.frame PhenoData; if there is no
#'  birth length variable in the data.frame PhenoData, can specify argument
#'  value is NULL. Birth length variable is expected to be in centimeters
#'@param HEADCIRCUMvar An optional character string indicating the column name
#'  for the head circumference variable in the data.frame PhenoData; if there is
#'  no head circumference variable in the data.frame PhenoData, can specify
#'  argument value is NULL. Head circumference variable is expected to be in
#'  centimeters
#'@param IDATdir A character string indicating the location of the IDAT files,
#'  e.g. "C:/Home/PACE/BirthSize/IDATs"
#'@param destinationfolder A character string indicating the location where
#'  files should be saved, e.g. "C:/Home/PACE/BirthSize"
#'@param savelog logic; TRUE indicates to save a log of the functions run
#'  (default)
#'@param cohort A character string for the cohort's acronym (e.g. "HEBC")
#'@param analysisdate A character string indicating the date the analysis was
#'  run. Please specify in the form: YEARMONTHDAY, e.g. "20200205" for February
#'  5th 2020
#'@details This function Renames the variables, puts the variables in the right
#'  mode for the analysis, and loads the raw array data.
#'@return A RGset (of class RGChannelSetExtended) with raw microarray intensity
#'  values. In addition to the columns included in the input PhenoData
#'  data.frame, the pData for this RGset will include the following columns:
#'  \item{ID}{Sample ID} \item{BWT}{Birthweight variable if BWTvar argument was
#'  not NULL; numeric} \item{Age}{Maternal age variable if AGEvar argument was
#'  not NULL; numeric} \item{BMI}{Maternal pre-pregnancy variable if BMIvar
#'  argument was not NULL; numeric} \item{Sex}{Infant sex variable if SEXvar
#'  argument was not NULL; factor, 'Female' or 'Male'} \item{Ethnic}{Maternal
#'  race/ethnicity variable if ETHNICvar argument was not NULL; factor}
#'  \item{Gestage}{Gestational age variable if GESTvar argument was not NULL;
#'  numeric}\item{Batch}{Batch variable if BATCHvar argument was not NULL;
#'  factor} \item{BirthLength}{Birth length variable if BIRTHLENGTHvar argument
#'  was not NULL; numeric} \item{HeadCircum}{Head circumference variable if
#'  HEADCIRCUMvar argument was not NULL; numeric} \item{BWTkg}{Birthweight in kg
#'  if BWTvar argument was not NULL} \item{wlr}{Birthweight birth length ratio
#'  if BWTvar and BirthLength arguments were not NULL}
#'  \item{BWT_Zscore}{Birthweight gestational age and sex-specific Z-scores
#'  based on Newborn Cross-Sectional Study of the INTERGROWTH-21st Project}
#'  \item{BirthLength_Zscore}{Birth length gestational age and sex-specific
#'  Z-scores based on Newborn Cross-Sectional Study of the INTERGROWTH-21st
#'  Project} \item{HeadCircum_Zscore}{Head circumference gestational age and
#'  sex-specific Z-scores based on Newborn Cross-Sectional Study of the
#'  INTERGROWTH-21st Project} \item{WLR_Zscore}{Birthweight birth length ratio
#'  gestational age and sex-specific Z-scores based on Newborn Cross-Sectional
#'  Study of the INTERGROWTH-21st Project}
#'@examples
#'\dontrun{
#'exampledat<-loadingSamples(SamplePlacement=NULL,PhenoData=allphenodata,IDlink="ID",
#'                           BWTvar="BWT",BATCHvar="Study",
#'                           SEXvar="GENDER_A",FemaleInd="1",MaleInd="2",
#'                           ETHNICvar="ETHNIC",GESTvar="gestAge",
#'                           BIRTHLENGTHvar="BirthLength",HEADCIRCUMvar=NULL,
#'                           IDATdir="H:/UCLA/PACE/Birthweight-placenta/IDATfiles",
#'                           destinationfolder="H:/UCLA/PACE/Birthweight-placenta",
#'                           savelog=TRUE,
#'                           cohort="HEBC",analysisdate="20200710")
#'}

loadingSamples<-function(SamplePlacement=NULL,PhenoData=PhenoData,IDlink="ID",
                         BWTvar=NULL,BATCHvar=NULL,
                         SEXvar=NULL,FemaleInd="1",MaleInd="2",
                         ETHNICvar=NULL,GESTvar=NULL,
                         BIRTHLENGTHvar=NULL,HEADCIRCUMvar=NULL,
                         IDATdir=NULL,destinationfolder=NULL,savelog=TRUE,
                         cohort=NULL,analysisdate=NULL){

  if(is.null(IDATdir)) stop("Please specify IDAT directory")

  Outputname<-paste(cohort,"_",analysisdate,"_Output",sep="")
  dir.create(file.path(destinationfolder, Outputname), showWarnings = FALSE)
  destinationfolder<-paste(destinationfolder,Outputname,sep="/")
  setwd(destinationfolder)

  if(savelog){

    sink()
    sink(type="message")

    con <- file("log_step1.log")
    sink(con, append=TRUE)
    sink(con, append=TRUE, type="message")
  }

  if(!is.null(SamplePlacement)){
    ## checking if there is an ID linking
    if (!(IDlink %in% colnames(SamplePlacement))) stop("Please make sure IDlink is in the SamplePlacement file")
    if (!(IDlink %in% colnames(PhenoData))) stop("Please make sure IDlink is in the PhenoData file")

    cat("PhenoData Number of Samples:",nrow(PhenoData),"\n")
    cat("SamplePlacement Number of Samples:",nrow(SamplePlacement),"\n")
    PhenoData$ID<-as.character(PhenoData[,IDlink])
    SamplePlacement$ID<-as.character(SamplePlacement[,IDlink])
    sampledat<-merge(PhenoData,SamplePlacement,by="ID")
  }

  if(is.null(SamplePlacement)){
    ## checking if there is an ID
    if (!(IDlink %in% colnames(PhenoData))) stop("Please make sure IDlink is in the PhenoData file")
    PhenoData$ID<-as.character(PhenoData[,IDlink])
    sampledat<-PhenoData
  }

  ## checking if Basename is in sampledat
  if(!("Basename" %in% colnames(sampledat))) {

    if(("Sentrix_ID" %in% colnames(sampledat)) & ("Sentrix_Position" %in% colnames(sampledat))){

      sampledat$Sentrix_ID<-as.character(sampledat$Sentrix_ID)
      sampledat$Sentrix_Position<-as.character(sampledat$Sentrix_Position)
      sampledat$Basename<-paste(sampledat$Sentrix_ID,sampledat$Sentrix_Position,sep="_")
      sampledat$Basename<-paste(sampledat$Sentrix_ID,sampledat$Sentrix_Position,sep="_")

    } else if("Basenames" %in% colnames(sampledat)){
      sampledat$Basename<-as.character(sampledat$Basenames)
    } else {
      stop("Please make sure 'Basename' is a column in the SamplePlacement csv")
    }
  }
  sampledat$Basename<-as.character(sampledat$Basename)
  setwd(IDATdir)
  IDATdircontents<-dir()
  IDATdircontents<-gsub("_Red.idat","",IDATdircontents)
  IDATdircontents<-gsub("_Grn.idat","",IDATdircontents)
  IDATdircontents<-unique(IDATdircontents)
  setwd(destinationfolder)

  IDATsMissing<-sampledat$Basename[!(sampledat$Basename %in% IDATdircontents)]
  if(length(IDATsMissing)>0){
    cat("IDATs missing from IDATdir folder:",paste(IDATsMissing,collapse = ", "),"\n")
  } else {
    cat("No IDATs missing from IDATdir folder","\n")
  }

  ## Checking which variables specified
  if(is.null(BWTvar)) cat("No birth weight variable indicated")
  if(is.null(BATCHvar)) cat("No batch indicated","\n")
  if(is.null(SEXvar)) stop("No infant sex variable indicated","\n")
  if(is.null(ETHNICvar)) cat("No ethnicity variable indicated","\n")
  if(is.null(GESTvar)) cat("No gestational age variable indicated","\n")
  if(is.null(BIRTHLENGTHvar)) cat("No birth length variable indicated","\n")
  if(is.null(HEADCIRCUMvar)) cat("No head circumference variable indicated","\n")

  ## If variable is not NULL, checking if it is in the sample data
  if(!is.null(BATCHvar)) if(!(BATCHvar %in% colnames(sampledat))) stop("Batch variable not in column names")
  if(!is.null(BWTvar)) if(!(BWTvar %in% colnames(sampledat)))  stop("Birth weight variable not in column names")
  if(!is.null(SEXvar)) if(!(SEXvar %in% colnames(sampledat))) stop("Infant sex variable not in column names")
  if(!is.null(ETHNICvar)) if(!(ETHNICvar %in% colnames(sampledat))) stop("Ethnicity variable not in column names")
  if(!is.null(GESTvar)) if(!(GESTvar %in% colnames(sampledat))) stop("Gestational age variable not in column names")
  if(!is.null(BIRTHLENGTHvar)) if(!(BIRTHLENGTHvar %in% colnames(sampledat))) stop("Birth length variable not in column names")
  if(!is.null(HEADCIRCUMvar)) if(!(HEADCIRCUMvar %in% colnames(sampledat))) stop("Head circumference variable not in column names")

  ## Changing column names
  sampledat$ID<-as.character(sampledat$ID)
  if(!is.null(BWTvar)) sampledat$BWT<-as.numeric(as.character(sampledat[,BWTvar]))
  if(!is.null(SEXvar)) sampledat$Sex<-sampledat[,SEXvar]
  if(!is.null(ETHNICvar)) sampledat$Ethnic<-as.factor(sampledat[,ETHNICvar])
  if(!is.null(GESTvar)) sampledat$Gestage<-as.numeric(as.character(sampledat[,GESTvar]))
  if(!is.null(BATCHvar)) sampledat$Batch<-as.factor(sampledat[,BATCHvar])
  if(!is.null(BIRTHLENGTHvar)) sampledat$BirthLength<-as.numeric(as.character(sampledat[,BIRTHLENGTHvar]))
  if(!is.null(HEADCIRCUMvar)) sampledat$HeadCircum<-as.numeric(as.character(sampledat[,HEADCIRCUMvar]))

  if("Sex" %in% colnames(sampledat)) {

    if(length(unique(sampledat$Sex))>2) stop("More than two sexes")

    sampledat$Sex<-as.character(sampledat$Sex)
    sampledat$Sex[which(sampledat$Sex==FemaleInd)]<-"Female"
    sampledat$Sex[which(sampledat$Sex==MaleInd)]<-"Male"
    sampledat$Sex<-as.factor(sampledat$Sex)

  }


  ## Calculating Z-scores based on Newborn Cross-Sectional Study of the INTERGROWTH-21st Project

  if("Sex" %in% colnames(sampledat) & "Gestage" %in% colnames(sampledat)) {
    if("BWT" %in% colnames(sampledat)){
      sampledat$BWTkg<-sampledat$BWT/1000
      sampledat$BWT_Zscore<-growthstandards::igb_wtkg2zscore(gagebrth=sampledat$Gestage*7, wtkg=sampledat$BWTkg, sex = as.character(sampledat$Sex))
    }
    if("BirthLength" %in% colnames(sampledat)){
      sampledat$BirthLength_Zscore<-growthstandards::igb_lencm2zscore(gagebrth=sampledat$Gestage*7, lencm=sampledat$BirthLength, sex = as.character(sampledat$Sex))
    }
    if("HeadCircum" %in% colnames(sampledat)){
      sampledat$HeadCircum_Zscore<-growthstandards::igb_hcircm2zscore(gagebrth=sampledat$Gestage*7, hcircm=sampledat$HeadCircum, sex = as.character(sampledat$Sex))
    }
    if("BWT" %in% colnames(sampledat) & "BirthLength" %in% colnames(sampledat)){
      sampledat$wlr<-sampledat$BWTkg/(sampledat$BirthLength/100) ## must be kg/m
      sampledat$wlr_Zscore<-growthstandards::igb_wlr2zscore(gagebrth=sampledat$Gestage*7, wlr=sampledat$wlr, sex = as.character(sampledat$Sex))
    }
  }

  cat("Loading IDAT files...","\n")
  msetraw<-minfi::read.metharray.exp(base=IDATdir,targets=sampledat,recursive=TRUE,extended = TRUE)
  msetraw

  if(savelog){
    sink()
    sink(type="message")
  }

  return(msetraw)

}
