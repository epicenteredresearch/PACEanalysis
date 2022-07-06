## ----eval=FALSE---------------------------------------------------------------
#  
#  ## First need to install required packages if you don't have them already
#  install.packages(c("ggplot2","gplots","reshape","RPMM","pvclust","doParallel",
#                     "GGally","Hmisc","MASS","sandwich", "lmtest","plyr","remotes","devtools","parallel","dplyr"))
#  
#  remotes::install_version("RefFreeEWAS", "2.2")
#  remotes::install_version("heatmap.plus", "1.3")
#  
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#  BiocManager::install(c("methylumi","minfi","sva","sesame","wateRmelon","EpiDISH",
#                         "IlluminaHumanMethylationEPICmanifest",
#                         "IlluminaHumanMethylation450kmanifest",
#                         "IlluminaHumanMethylation450kanno.ilmn12.hg19",
#                         "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
#                         "TxDb.Hsapiens.UCSC.hg19.knownGene",
#                         "org.Hs.eg.db","FDb.InfiniumMethylation.hg19",
#                         "FlowSorted.CordBloodCombined.450k",
#                         "FlowSorted.Blood.EPIC",
#                         "FlowSorted.CordBlood.450k",
#                         "FlowSorted.Blood.450k",
#                         "illuminaio"))
#  
#  remotes::install_github("bokeh/rbokeh")
#  remotes::install_github("ki-tools/growthstandards")
#  remotes::install_github("hhhh5/ewastools")
#  
#  ## If ExperimentHub (>1.17.2), need to update caching location
#  moveFiles<-function(package){
#         olddir <- path.expand(rappdirs::user_cache_dir(appname=package))
#         newdir <- tools::R_user_dir(package, which="cache")
#         dir.create(path=newdir, recursive=TRUE)
#         files <- list.files(olddir, full.names =TRUE)
#         moveres <- vapply(files,
#             FUN=function(fl){
#             filename = basename(fl)
#             newname = file.path(newdir, filename)
#             file.rename(fl, newname)
#             },
#             FUN.VALUE = logical(1))
#         if(all(moveres)) unlink(olddir, recursive=TRUE)
#         }
#  
#  package="ExperimentHub"
#  moveFiles(package)
#  
#  ## If recently installed sesame, need to cache the associated annotation data
#  ## This only needs to be done once per new installation of sesame
#  sesameData::sesameDataCacheAll()
#  eh<-ExperimentHub()
#  eh[["EH6019"]] ## one that isn't automatically downloaded
#  
#  ## Need to then install package, specifying path to the source package
#  install.packages("G:\\PACE\\PACEanalysis_0.1.8.tar.gz",
#                   repos = NULL, type="source")
#  

## ----eval=FALSE---------------------------------------------------------------
#  ## Attach package
#  library(PACEanalysis)
#  
#  setwd("G:\\PACE\\Birthweight-placenta")
#  allphenodata<-read.csv("All_Pheno_and_Basenames.csv",header=TRUE)
#  dim(allphenodata)
#  
#  ## If the data.frame you are going to specify as PhenoData does not include
#  ## the column "Basename", you will need to also load a data.frame you are
#  ## going to specify for the SamplePlacement argument. PhenoData and SamplePlacement
#  ## are merged on the column name specified by the IDlink argument
#  
#  exampledat<-loadingSamples(SamplePlacement=NULL,PhenoData=allphenodata,IDlink="ID",
#                    BWTvar="BWT",BATCHvar="Study",
#                    SEXvar="GENDER_A",FemaleInd="1",MaleInd="2",
#                    ETHNICvar="ETHNIC",GESTvar="gestAge",
#                    BIRTHLENGTHvar="BirthLength",HEADCIRCUMvar=NULL,
#                    IDATdir="H:\\UCLA\\PACE\\Birthweight-placenta\\IDATfiles",
#                    destinationfolder="H:\\UCLA\\PACE\\Birthweight-placenta",
#                    savelog=TRUE,
#                    cohort="HEBC",analysisdate="20201229")
#  
#  EDAresults<-ExploratoryDataAnalysis(RGset=exampledat,
#                    globalvarexplore=c("BWT","Sex"),
#                    DetectionPvalMethod="SeSAMe",
#                    DetectionPvalCutoff=0.05,
#                    minNbeads=3,
#                    FilterZeroIntensities=TRUE,
#                    destinationfolder="H:\\UCLA\\PACE\\Birthweight-placenta",
#                    savelog=TRUE,
#                    cohort="HEBC",analysisdate="20201229")

## ----eval=FALSE---------------------------------------------------------------
#  
#  processedOut<-preprocessingofData(RGset=exampledat,
#                    SamplestoRemove=EDAresults$SamplestoRemove,
#                    ProbestoRemove=EDAresults$ProbestoRemove,
#                    destinationfolder="H:\\UCLA\\PACE\\Birthweight-placenta",
#                    compositeCellType="Placenta",
#                    KchooseManual=NULL,
#                    savelog=TRUE,
#                    cohort="HEBC",analysisdate="20201229")
#  
#  betasabovedetection<-detectionMask(processedBetas=processedOut$processedBetas,
#                                    DetectionPvals=EDAresults$DetectionPval,
#                                    DetectionPvalCutoff=0.05,
#                                    IndicatorGoodIntensity=EDAresults$IndicatorGoodIntensity,
#                                    destinationfolder="H:\\UCLA\\PACE\\Birthweight-placenta",
#                                    cohort="HEBC",analysisdate="20201229")
#  
#  Betasnooutliers<-outlierprocess(processedBetas=betasabovedetection,
#                                    quantilemethod="Quantile",
#                                    trimming=TRUE,
#                                    pct=0.25,
#                                    destinationfolder="H:\\UCLA\\PACE\\Birthweight-placenta",
#                                    cohort="HEBC",analysisdate="20201229")
#  

## ----eval=FALSE---------------------------------------------------------------
#  
#  ## If you closed prior R session, you can load list of pre-processed objects
#  ## that is automatically saved by the preprocessingofData function
#  
#  setwd("H:\\UCLA\\PACE\\Birthweight-placenta\\HEBC_20201229_Output")
#  load("HEBC_20201229_Preprocessed.RData")
#  
#  phenodataframe<-as.data.frame(pData(processedOut$mset))
#  phenodataframe$LBWbin<-ifelse(phenodataframe$BWT<2500,1,0)
#  phenodataframe$HBWbin<-ifelse(phenodataframe$BWT>4000,1,0)
#  
#  modelstorun<-data.frame(varofinterest=c("BWT","LBWbin","HBWbin","BirthLength",
#                                          "WeightLengthRatio","HeadCircum"))
#  modelstorun$varofinterest<-as.character(modelstorun$varofinterest)
#  modelstorun$vartype<-"OutcomeCont"
#  modelstorun$vartype[modelstorun$varofinterest %in% c("LBWbin","HBWbin")]<-"OutcomeBin"
#  
#  ## You can reduce this dataframe to whatever variables you have.
#  ## For example, if you only have birthweight, you would specify:
#  ## modelstorun<-data.frame(varofinterest=c("BWT","LBWbin","HBWbin"))
#  
#  ## Make sure that you have sufficient numbers of samples between categories
#  ## of exposure/outcome of interest as well as adjustment variables;
#  ## you are assumed to have at least 10
#  
#  table(phenodataframe$LBWbin)
#  table(phenodataframe$HBWbin)
#  table(phenodataframe$Sex)
#  table(phenodataframe$Parity)
#  table(phenodataframe$MaternalEd)
#  table(phenodataframe$Smoke)
#  table(phenodataframe$Ethnic)
#  
#  ## Make sure all categorical adjustment variables are coded as factors or characters
#  ## 'Sex' is already coded as a character
#  
#  phenodataframe$Parity<-as.factor(phenodataframe$Parity)
#  phenodataframe$MaternalEd<-as.factor(phenodataframe$MaternalEd)
#  phenodataframe$Smoke<-as.factor(phenodataframe$Smoke)
#  phenodataframe$Ethnic<-as.factor(phenodataframe$Ethnic)
#  

## ----eval=FALSE---------------------------------------------------------------
#  
#  for (i in 1:nrow(modelstorun)){
#  
#    cat("OutcomG:",modelstorun$varofinterest[i],"\n")
#    tempresults<-dataAnalysis(phenofinal=phenodataframe,
#                    betafinal=Betasnooutliers[1:100,], ## restricting to first 100 loci
#                    array="450K",
#                    maxit=100,
#                    Omega=processedOut$Omega,
#                    vartype=modelstorun$vartype[i],
#                    robust=TRUE,
#                    varofinterest=modelstorun$varofinterest[i],
#                    Table1vars=c("Gestage","Sex","Age","Parity","MaternalEd",
#                                     "Smoke","preBMI","Ethnic"),
#                    StratifyTable1=FALSE,
#                    StratifyTable1var=NULL,
#                    adjustmentvariables=c("Gestage","Sex","Age","Parity","MaternalEd",
#                                     "Smoke","preBMI","Ethnic"),
#                    RunUnadjusted=TRUE,
#                    RunAdjusted=TRUE,
#                    RunCellTypeAdjusted=TRUE,
#                    RunSexSpecific=TRUE,
#                    RunCellTypeInteract=TRUE,
#                    RestrictToSubset=FALSE,
#                    RestrictionVar=NULL,
#                    RestrictToIndicator=NULL,
#                    destinationfolder="H:\\UCLA\\PACE\\Birthweight-placenta",
#                    savelog=TRUE,
#                    cohort="HEBC",analysisdate="20210103",analysisname="main")
#  
#  }
#  

## ----eval=FALSE---------------------------------------------------------------
#  for (i in 1:nrow(modelstorun)){
#  
#    cat("OutcomG:",modelstorun$varofinterest[i],"\n")
#    tempresults<-dataAnalysis(phenofinal=phenodataframe,
#                    betafinal=Betasnooutliers,
#                    array="450K",
#                    maxit=100,
#                    Omega=processedOut$Omega,
#                    vartype=modelstorun$vartype[i],
#                    robust=TRUE,
#                    varofinterest=modelstorun$varofinterest[i],
#                    Table1vars=c("Gestage","Sex","Age","Parity","MaternalEd",
#                                     "Smoke","preBMI","Ethnic"),
#                    StratifyTable1=FALSE,
#                    StratifyTable1var=NULL,
#                    adjustmentvariables=c("Gestage","Sex","Age","Parity","MaternalEd",
#                                     "Smoke","preBMI","Ethnic"),
#                    RunUnadjusted=TRUE,
#                    RunAdjusted=TRUE,
#                    RunCellTypeAdjusted=TRUE,
#                    RunSexSpecific=TRUE,
#                    RunCellTypeInteract=TRUE,
#                    RestrictToSubset=FALSE,
#                    RestrictionVar=NULL,
#                    RestrictToIndicator=NULL,
#                    destinationfolder="H:\\UCLA\\PACE\\Birthweight-placenta",
#                    savelog=TRUE,
#                    cohort="HEBC",analysisdate="20210103")
#  
#      tempresultsNonHispanicWhite<-dataAnalysis(phenofinal=phenodataframe,
#                    betafinal=Betasnooutliers,
#                    array="450K",
#                    maxit=100,
#                    Omega=processedOut$Omega,
#                    vartype=modelstorun$vartype[i],
#                    robust=TRUE,
#                    varofinterest=modelstorun$varofinterest[i],
#                    Table1vars=c("Gestage","Sex","Age","Parity","MaternalEd",
#                                     "Smoke","preBMI"),
#                    StratifyTable1=FALSE,
#                    StratifyTable1var=NULL,
#                    adjustmentvariables=c("Gestage","Sex","Age","Parity","MaternalEd",
#                                     "Smoke","preBMI"),
#                    RunUnadjusted=TRUE,
#                    RunAdjusted=TRUE,
#                    RunCellTypeAdjusted=TRUE,
#                    RunSexSpecific=TRUE,
#                    RunCellTypeInteract=TRUE,
#                    RestrictToSubset=TRUE,
#                    RestrictionVar="Ethnic",
#                    RestrictToIndicator="1",
#                    destinationfolder="H:\\UCLA\\PACE\\Birthweight-placenta",
#                    savelog=TRUE,
#                    cohort="HEBC",analysisdate="20210103",analysisname="main")
#  
#  }
#  

## ----eval=FALSE---------------------------------------------------------------
#  
#  for (i in 1:nrow(modelstorun)){
#  
#    cat("OutcomG:",modelstorun$varofinterest[i],"\n")
#    tempresults<-dataAnalysis(phenofinal=phenodataframe,
#                    betafinal=Betasnooutliers,
#                    array="450K",
#                    maxit=100,
#                    Omega=processedOut$Omega,
#                    vartype=modelstorun$vartype[i],
#                    robust=TRUE,
#                    varofinterest=modelstorun$varofinterest[i],
#                    Table1vars=c("Gestage","Sex","Age","Parity","MaternalEd",
#                                     "Smoke","preBMI"),
#                    StratifyTable1=FALSE,
#                    StratifyTable1var=NULL,
#                    adjustmentvariables=c("Gestage","Sex","Age","Parity","MaternalEd",
#                                     "Smoke","preBMI"),
#                    RunUnadjusted=TRUE,
#                    RunAdjusted=TRUE,
#                    RunCellTypeAdjusted=TRUE,
#                    RunSexSpecific=TRUE,
#                    RunCellTypeInteract=TRUE,
#                    RestrictToSubset=FALSE,
#                    RestrictionVar=NULL,
#                    RestrictToIndicator=NULL,
#                    destinationfolder="H:\\UCLA\\PACE\\Birthweight-placenta",
#                    savelog=TRUE,
#                    cohort="HEBC",analysisdate="20210103",analysisname="main")
#  
#  }
#  

