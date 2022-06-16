#' "fix" an RGChannelSet (for which IDATs may be unavailable) with Sesame
#' The input is an RGSet and the output is a sesamized minfi::GenomicRatioSet
#' 
#' \code{HDF5Array} package required.
#' 
#' @param rgSet an RGChannelSet, perhaps with colData of various flavors
#' @param naFrac maximum NA fraction for a probe before it gets dropped (1)
#' @param BPPARAM get parallel with MulticoreParam(n)
#' @param HDF5 is the rgSet HDF5-backed? if so, avoid eating RAM (perhaps)
#' @param HDF5SEdestination character(1) path to where the
#' HDF5-backed GenomicRatioSet will be stored
#' @param replace logical(1) passed to saveHDF5SummarizedExperiment
#' @note We employ BPREDO for a second chance if bplapply hits an error.
#' @return a sesamized GenomicRatioSet
#' @import BiocParallel
#' @importFrom S4Vectors metadata
#' @importFrom S4Vectors metadata<-
#' @importFrom SummarizedExperiment assays
#' @examples
#' 
#' \dontrun{
#'
#' library(FlowSorted.CordBloodNorway.450k)
#' library(BiocParallel)
#' sesamize(FlowSorted.CordBloodNorway.450k[,1:2],
#'     BPPARAM=MulticoreParam(2))
#'
#' }
#' @export 
sesamize <- function(
    rgSet, naFrac=1, BPPARAM=SerialParam(), HDF5=NULL,
    HDF5SEdestination=paste0(tempdir(check=TRUE), "/sesamize_HDF5_scratch"),
    replace=FALSE) {
  
  stopifnot(is(rgSet, "RGChannelSet"))
  
  #pkgTest('minfi')
  #pkgTest('SummarizedExperiment')
  indices = seq_len(ncol(rgSet))
  
  if (is.null(HDF5)) {
    ## are we working on an HDF5-backed RGChannelSet?
    HDF5 <- is(assays(rgSet)[[1]], "DelayedMatrix")
  }
  t1 = BiocParallel::bptry(BiocParallel::bplapply(indices, function(i) {
    message("Sesamizing ", i, "...")
    sdf <- RGChannelSet1ToSigDF(rgSet[,i])
    sdf <- sesame::dyeBiasNL(sesame::noob(sdf))
    SigDFToRatioSet(sdf)}, BPPARAM=BPPARAM))
  lk = vapply(t1, inherits, logical(1), "bperror")  # second try?
  if (any(lk)) {
    t1 = BiocParallel::bptry(BiocParallel::bplapply(indices, function(i) {
      message("Sesamizing ", i, "...")
      sdf <- RGChannelSet1ToSigDF(rgSet[,i])
      sdf <- sesame::dyeBiasNL(sesame::noob(sdf))
      SigDFToRatioSet(sdf)}, BPREDO=t1, BPPARAM=BPPARAM))
  }
  
  ratioSet <- do.call(SummarizedExperiment::cbind, t1)
  colnames(ratioSet) = colnames(rgSet)
  if (HDF5) {
    # pkgTest('HDF5Array')
    ##td <- paste(tempdir(check=TRUE), "sesamize_HDF5_scratch", sep="/")
    ratioSet <- HDF5Array::saveHDF5SummarizedExperiment(
      ratioSet, dir=HDF5SEdestination, replace=replace) #td, replace=TRUE)
  }
  ## mapping occurs first, SNPs get separated here
  ratioSet <- minfi::mapToGenome(ratioSet)
  
  ## keep only probes surviving naFrac
  kept <- seq_len(nrow(ratioSet))
  if (naFrac < 1) { 
    kept <- which((
      rowSums(is.na(minfi::getBeta(ratioSet)))/ncol(ratioSet)) <= naFrac)
    if (length(kept) < 1) 
      stop("No probes survived with naFrac <= ",naFrac,".")
  } 
  
  ## put back colData(), @processMethod, $SNPs
  mfst <- packageVersion(paste(minfi::annotation(ratioSet), collapse="anno."))
  ratioSet@preprocessMethod <- c(
    rg.norm="SeSAMe (type I)",
    p.value="SeSAMe (pOOBAH)",
    sesame=as.character(packageVersion("sesame")),
    minfi=as.character(packageVersion("minfi")),
    manifest=as.character(mfst))
  
  ## SNP not adjusted in minfi, so keep them that way
  metadata(ratioSet)$SNPs <- minfi::getSnpBeta(rgSet)
  SummarizedExperiment::assays(ratioSet)[["M"]] <- NULL 
  SummarizedExperiment::colData(ratioSet) <- colData(rgSet)
  
  return(ratioSet[kept, ])
}

platformMinfiToSm <- function(platform) {
  plf <- sub("HMEPIC", "EPIC", 
             sub("IlluminaHumanMethylation", "HM", 
                 sub("k$", "", platform)))
  stopifnot(plf %in% c('EPIC','HM450','HM27'))
  plf
}

platformSmToMinfi <- function(platform) {
  plf <- sub(
    "HM", "IlluminaHumanMethylation",
    sub("EPIC", "HMEPIC",
        sub("(450|27)$", "\\1k", platform)))
  stopifnot(plf %in% c(
    'IlluminaHumanMethylationEPIC',
    'IlluminaHumanMethylation450k',
    'IlluminaHumanMethylation27k'
  ))
  plf
}

## reverse of chipAddressToSignal
SigDFToRGChannel <- function(sdf, manifest = NULL, controls = NULL) {
  
  if (is.null(manifest)) {
    tempplatform<-paste0(sesame::sdfPlatform(sdf),'.address')
    dfAddress <- sesameData::sesameDataGet(tempplatform)
    manifest <- dfAddress$ordering
    controls <- dfAddress$controls
  }
  
  d1 = InfI(sdf)
  SSRed <- c(
    setNames(d1$MR, manifest$M[match(d1$Probe_ID, manifest$Probe_ID)]),
    setNames(sdf$UR, manifest$U[match(sdf$Probe_ID, manifest$Probe_ID)]))
  SSGrn <- c(
    setNames(d1$MG, manifest$M[match(d1$Probe_ID, manifest$Probe_ID)]),
    setNames(sdf$UG, manifest$U[match(sdf$Probe_ID, manifest$Probe_ID)]))
  
  ## controls
  if (!is.null(controls)) {
    ctl = controls(sdf)
    control.names <- make.names(controls$Name, unique = TRUE)
    SSGrn <- c(SSGrn, setNames(ctl[match(
      control.names, rownames(ctl)),'G'],
      as.character(controls$Address)))
    SSRed <- c(SSRed, setNames(ctl[match(
      control.names, rownames(ctl)),'R'],
      as.character(controls$Address)))
  } ## else TODO controls obtained from manifest
  
  list(grn=SSGrn, red=SSRed)
}

## annotation, if not given is guessed
guessMinfiAnnotation <- function(ptf, annotation = NA) {
  if (ptf == "MM285") {
    ## 20211109: still waiting for the BioC official code
    stop("SigDFsToRGChannelSet does not support mouse array.")
  }
  if (is.na(annotation)) {
    if (ptf %in% c("HM450", "HM27")) {
      'ilmn12.hg19'
    } else { # EPIC
      'ilm10b4.hg19'
    }
  } else {
    annotation
  }
}

#' Convert sesame::SigDF to minfi::RGChannelSet
#'
#' This function does not support the mouse array.
#' 
#' @param sdfs a list of sesame::SigDF
#' @param BPPARAM get parallel with MulticoreParam(n)
#' @param annotation the minfi annotation string, guessed if not given
#' @return a minfi::RGChannelSet
#' @import BiocParallel
#' @examples
#'
#' \donttest{
#'
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' rgSet <- SigDFsToRGChannelSet(sdf)
#'
#' }
#' 
#' @export 
SigDFsToRGChannelSet <- function(sdfs, BPPARAM=SerialParam(), annotation=NA) {
  
  if (is(sdfs, 'SigDF')) {
    sdfs <- list(sample=sdfs)
  }
  
  pt <- sesame::sdfPlatform(sdfs[[1]])
  annotation <- guessMinfiAnnotation(pt, annotation)
  
  ss_all <- BiocParallel::bplapply(sdfs, SigDFToRGChannel, BPPARAM=BPPARAM)
  rgset <- minfi::RGChannelSet(
    Green=do.call(cbind, lapply(ss_all, function(ss) ss$grn)), 
    Red=do.call(cbind, lapply(ss_all, function(ss) ss$red)), 
    annotation=c(
      array=unname(platformSmToMinfi(pt)),
      annotation=annotation))
}


#' Lookup address in one sample
#'
#' Lookup address and transform address to probe
#'
#' Translate data in chip address to probe address.
#' Type I probes can be separated into Red and Grn channels. The
#' methylated allele and unmethylated allele are at different
#' addresses. For type II probes methylation allele and unmethylated allele are
#' at the same address. Grn channel is for methylated allele and Red channel is
#' for unmethylated allele. The out-of-band signals are type I probes measured
#' using the other channel.
#'
#' @param dm data frame in chip address, 2 columns: cy3/Grn and cy5/Red
#' @param mft a data frame with columns Probe_ID, M, U and col
#' @return a SigDF, indexed by probe ID address
chipAddressToSignal <- function(dm, mft) {
  
  ## Infinium-I
  mft1 <- mft[!is.na(mft$col),]
  tmpM <- dm[match(mft1$M, rownames(dm)),]
  tmpU <- dm[match(mft1$U, rownames(dm)),]
  sdf <- data.frame(
    Probe_ID=mft1$Probe_ID,
    MG=tmpM[,"G"], MR=tmpM[,"R"], UG=tmpU[,"G"], UR=tmpU[,"R"],
    col=mft1$col, mask=FALSE)
  ## if ("mask" %in% colnames(mft1)) { sdf$mask <- mft1$mask; }
  
  ## Infinium-II
  mft2 <- mft[is.na(mft$col),]
  if (nrow(mft2) > 0) {
    tmp <- dm[match(mft2$U, rownames(dm)),]
    s2 <- data.frame(
      Probe_ID=mft2$Probe_ID,
      MG=NA, MR=NA, UG=tmp[,"G"], UR=tmp[,"R"], col="2", mask=FALSE)
    ## if ("mask" %in% colnames(mft2)) { s2$mask <- mft2$mask; }
    sdf <- rbind(sdf, s2)
  }
  sdf$col <- factor(sdf$col, levels=c("G","R","2"))
  sdf <- sdf[match(mft$Probe_ID, sdf$Probe_ID),] # always the mft order
  sdf <- structure(sdf, class=c("SigDF", "data.frame"))
  rownames(sdf) <- NULL
  sdf
}

readControls <- function(dm, controls) {
  if ("Color_Channel" %in% colnames(controls)) { # legacy control data
    ctl <- as.data.frame(dm[match(controls$Address, rownames(dm)),])
    rownames(ctl) <- make.names(controls$Name, unique=TRUE)
    ctl <- cbind(ctl, controls[, c("Color_Channel","Type")])
    colnames(ctl) <- c('G','R','col','type')
    ctl <- ctl[!(is.na(ctl$G)|is.na(ctl$R)),] # no NA in controls
  } else {
    ctl <- as.data.frame(chipAddressToSignal(dm, controls))
  }
  ctl
}

## helper: convert RGChannelSet of one sample
RGChannelSet1ToSigDF <- function(rgSet1, manifest = NULL, controls = NULL) {
  
  stopifnot(ncol(rgSet1) == 1)
  
  ## chipaddress/rownames are automatically the same
  dm <- cbind(
    G=as.matrix(minfi::getGreen(rgSet1)),
    R=as.matrix(minfi::getRed(rgSet1)))
  
  colnames(dm) <- c('G','R') # just in case..
  if (is.null(manifest)) {
    attr(dm, 'platform') <- platformMinfiToSm(
      minfi::annotation(rgSet1)['array'])
    df_address <- sesameData::sesameDataGet(paste0(
      attr(dm, 'platform'), '.address'))
    manifest <- df_address$ordering
    controls <- df_address$controls
  }
  
  tempsdf<-chipAddressToSignal(dm, manifest)
  sdf = sesame::pOOBAH(sesame::qualityMask(tempsdf))
  if (!is.null(controls)) {
    attr(sdf, "controls") = readControls(dm, controls)
  }
  sdf
}

#' Convert RGChannelSet (minfi) to a list of SigDF (SeSAMe)
#'
#' Notice the colData() and rowData() is lost.
#'
#' @param rgSet a minfi::RGChannelSet
#' @param BPPARAM get parallel with MulticoreParam(n)
#' @param manifest manifest file
#' @return a list of sesame::SigDF
#' @import BiocParallel
#' @examples
#' 
#' \donttest{
#'
#' library(FlowSorted.Blood.450k)
#' rgSet <- FlowSorted.Blood.450k[,1:2]
#' sdfs <- RGChannelSetToSigDFs(rgSet)
#' ## set names(sdfs) yourself. each data has a different place
#' ## to store sample name
#'
#' }
#' 
#' @export
RGChannelSetToSigDFs <- function(
    rgSet, manifest=NULL, BPPARAM=SerialParam()) {
  
  samples <- colnames(rgSet)
  BiocParallel::bplapply(
    seq_len(ncol(rgSet)), function(i) {
      RGChannelSet1ToSigDF(rgSet[,i], manifest=manifest)
    }, BPPARAM=BPPARAM)
}

#' Convert one sesame::SigDF to minfi::RatioSet
#'
#' @param sdf a sesame::SigDF
#' @param annotation minfi annotation string
#' @return a minfi::RatioSet
#' @examples
#' 
#' \donttest{
#' 
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' ratioSet <- SigDFToRatioSet(sdf)
#'
#' }
#'
#' @export
SigDFToRatioSet <- function(sdf, annotation = NA) {
  Beta <- as.matrix(sesame::getBetas(sdf))
  CN <- as.matrix(log2(sesame::totalIntensities(sdf))[rownames(Beta)])
  annotation <- guessMinfiAnnotation(sesame::sdfPlatform(sdf), annotation)
  platform <- platformSmToMinfi(sesame::sdfPlatform(sdf))
  minfi::RatioSet(Beta = Beta, CN = CN, annotation = c(
    array = unname(platform), annotation = annotation))
}
