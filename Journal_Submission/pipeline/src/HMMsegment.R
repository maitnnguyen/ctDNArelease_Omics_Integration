suppressMessages({
  #library(HMMcopy)
  library(GenomicRanges)
  library(optparse)
  library(ichorCNA)
  library(PureCN)
  library(dplyr)
})

options(stringsAsFactors=FALSE, scipen=0)
options(bitmapType='cairo')

option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "normalize coverage file by patient"),
  make_option(c("-o", "--outfile"), type = "character", help = "Output file."),
  make_option(c("-O", "--outDir"), type = "character", help = "Output folder for saving visualisations"),
  make_option(c("-g", "--gcmap"), type="character", help = "File containing GC and mappability info for correction"),
  make_option(c("-c", "--centromere"), type="character", help = "File containing Centromere locations"),
  make_option(c("-p", "--PoN"), type =  "character", help = "Panel of normal file"),
  make_option(c("--chrs"), type="character", default="c(1:22,\"X\")", help = "Specify chromosomes to analyze."),
  make_option(c("--method"), type = "character", default="median", help="Median or Mean."),
  make_option(c("--chrNormalize"), type="character", default="c(1:22)", help = "Specify chromosomes to normalize GC/mappability biases"),
  make_option(c("--chrTrain"), type="character", default="c(1:22)", help = "Specify chromosomes to estimate params. Default: [%default]"),
  make_option(c("--maleChrXLogRThres"), type="numeric", default=-0.80, help = "ChrX Log ratio threshold to confirm as male gender."),
  make_option(c("--normal"), type="character", default="seq(0.5:0.9,by=.1)", help = "Initial normal contamination; can be more than one value if additional normal initializations are desired. Default: [%default]"),
  make_option(c("--coverage"), type="numeric", default=.1, help = "PICARD sequencing coverage. Default: [%default]"),
  make_option(c("--scStates"), type="character", default="c(1,3)", help = "Subclonal states to consider. Default: [%default]"),
  make_option(c("--lambda"), type="character", default="NULL", help="Initial Student's t precision; must contain 4 values (e.g. c(1500,1500,1500,1500)); if not provided then will automatically use based on variance of data. Default: [%default]"),
  make_option(c("--lambdaScaleHyperParam"), type="numeric", default=3, help="Hyperparameter (scale) for Gamma prior on Student's-t precision. Default: [%default]"),
  make_option(c("--ploidy"), type="character", default="c(2:4)", help = "Initial tumour ploidy; can be more than one value if additional ploidy initializations are desired. Default: [%default]"),
  make_option(c("--maxCN"), type="numeric", default=10, help = "Total clonal CN states. Default: [%default]"),
  make_option(c("--estimateNormal"), type="logical", default=TRUE, help = "Estimate normal. Default: [%default]"),
  make_option(c("--estimateScPrevalence"), type="logical", default=TRUE, help = "Estimate subclonal prevalence. Default: [%default]"),
  make_option(c("--estimatePloidy"), type="logical", default=TRUE, help = "Estimate tumour ploidy. Default: [%default]"),
  make_option(c("--maxFracCNASubclone"), type="numeric", default=0.7, help="Exclude solutions with fraction of subclonal events greater than this value. Default: [%default]"),
  make_option(c("--maxFracGenomeSubclone"), type="numeric", default=0.5, help="Exclude solutions with subclonal genome fraction greater than this value. Default: [%default]"),
  make_option(c("--minSegmentBins"), type="numeric", default=50, help="Minimum number of bins for largest segment threshold required to estimate tumor fraction; if below this threshold, then will be assigned zero tumor fraction."),
  make_option(c("--altFracThreshold"), type="numeric", default=0.01, help="Minimum proportion of bins altered required to estimate tumor fraction; if below this threshold, then will be assigned zero tumor fraction. Default: [%default]"),
  make_option(c("--minTumFracToCorrect"), type="numeric", default=0.05, help = "Tumor-fraction correction of bin and segment-level CNA if sample has minimum estimated tumor fraction. [Default: %default]")
  )
parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
print(opt)

# source for supporting import coverage file
source("~/src/R/utils.R")

# input import and parameters setup
{
  # infor of input
  inputFile <- opt$input
  outfile <- opt$outfile 
  outDir <- opt$outDir # save visualisations
  #patientID <- opt$patientID
  # chr infor
  chrs <- as.character(eval(parse(text = opt$chrs)))
  chrNormalize <- as.character(eval(parse(text=opt$chrNormalize)))
  chrTrain <- as.character(eval(parse(text=opt$chrTrain)))
  gc.map.file <- opt$gcmap
  normDB <- opt$PoN
  centromereFile <- opt$centromere
  method <- opt$method
  maleChrXLogRThres <- opt$maleChrXLogRThres

  # parameters
  normal <- eval(parse(text = opt$normal))
  scStates <- eval(parse(text = opt$scStates))
  lambda <- eval(parse(text = opt$lambda))
  lambdaScaleHyperParam <- opt$lambdaScaleHyperParam
  estimateNormal <- opt$estimateNormal
  estimatePloidy <- opt$estimatePloidy
  estimateScPrevalence <- opt$estimateScPrevalence
  maxFracCNASubclone <- opt$maxFracCNASubclone
  maxFracGenomeSubclone <- opt$maxFracGenomeSubclone
  minSegmentBins <- opt$minSegmentBins
  altFracThreshold <- opt$altFracThreshold
  minTumFracToCorrect <- opt$minTumFracToCorrect
  ploidy <- eval(parse(text = opt$ploidy))
  maxCN <- opt$maxCN
  coverage <- opt$coverage
}

# proceed file

{
  patientID = 'patient'
  files = read.delim(inputFile)
  # normalize to PoN
  numSamples <- nrow(files) #length(inputFile)
  tumour_copy <- list()
  for (i in 1:numSamples) {
  id <- files[i,'SampleID']
  # load raw coverage and correct bias
  rawCounts <- loadReadCount(files[i,'File'], map.gc.AnotFile = gc.map.file,
                            centromere = centromereFile)
  normCounts <- correctReadCounts(rawCounts, chrNormalize=chrNormalize)
  # get gender
  gender <- extractGender(rawReads = rawCounts, normReads = normCounts)

  tumour_copy[[id]] <- ichorCNA::normalizeByPanelOrMatchedNormal(normCounts, chrs = chrs, 
                                    normal_panel = normDB, gender = gender$gender)

  outMat <- as.data.frame(tumour_copy[[id]])
  outMat <- outMat[,c("seqnames","start","end","copy")]
  colnames(outMat) <- c("chr","start","end","log2_TNratio_corrected")
  outFile <- paste0(outDir,"/",id,".correctedDepth.txt")
 # message(paste("Outputting to:", outFile))
  write.table(outMat, file=outFile, row.names=F, col.names=T, quote=F, sep="\t")
  }
}
  chrInd <- as.character(seqnames(tumour_copy[[1]])) %in% chrTrain
  ## get positions that are valid
  valid <- tumour_copy[[1]]$valid
  if (length(tumour_copy) >= 2) {
  for (i in 2:length(tumour_copy)){ 
  valid <- valid & tumour_copy[[i]]$valid 
    }
  }


  # set parameter for segmentation model
{
ptmTotalSolutions <- proc.time() # start total timer
results <- list()
loglik <- as.data.frame(matrix(NA, nrow = length(normal) * length(ploidy), ncol = 7, 
                 dimnames = list(c(), c("init", "n_est", "phi_est", "BIC", 
                                        "Frac_genome_subclonal", "Frac_CNA_subclonal", "loglik"))))
counter <- 1
compNames <- rep(NA, nrow(loglik))
mainName <- rep(NA, length(normal) * length(ploidy))
}
  
for (n in normal){
  for (p in ploidy){
    if (n == 0.95 & p != 2) {
        next
    }
    logR <- as.data.frame(lapply(tumour_copy, function(x) { x$copy })) # NEED TO EXCLUDE CHR X #
    param <- ichorCNA::getDefaultParameters(logR[valid & chrInd, , drop=F], maxCN = maxCN,  
                ct.sc=scStates, ploidy = floor(p), e.same = 50)
    param$phi_0 <- rep(p, numSamples)
    param$n_0 <- rep(n, numSamples)
    
    ############################################
    ######## CUSTOM PARAMETER SETTINGS #########
    ############################################
    # 0.1x cfDNA #
    logR.var <- 1 / ((apply(logR, 2, sd, na.rm = TRUE) / sqrt(length(param$ct))) ^ 2)
    param$lambda <- rep(logR.var, length(param$ct))
    param$lambda[param$ct %in% c(2)] <- logR.var 
    param$lambda[param$ct %in% c(1,3)] <- logR.var 
    param$lambda[param$ct >= 4] <- logR.var / 5
    param$lambda[param$ct == max(param$ct)] <- logR.var / 15
    param$lambda[param$ct.sc.status] <- logR.var / 10
    param$alphaLambda <- rep(lambdaScaleHyperParam, length(param$ct)) 

    # HMM segmentation
    hmmResults.cor <- ichorCNA::HMMsegment(tumour_copy, valid, dataType = "copy", 
                                 param = param, chrTrain = chrTrain, maxiter = 20,
                                 estimateNormal = estimateNormal, estimatePloidy = estimatePloidy,
                                 estimateSubclone = estimateScPrevalence, verbose = TRUE)

    for (s in 1:numSamples){
      iter <- hmmResults.cor$results$iter
      id <- names(hmmResults.cor$cna)[s]

      ## convert full diploid solution (of chrs to train) to have 1.0 normal or 0.0 purity
      ## check if there is an altered segment that has at least a minimum # of bins
      segsS <- hmmResults.cor$results$segs[[s]]
      segsS <- segsS[segsS$chr %in% chrTrain, ]
      segAltInd <- which(segsS$event != "NEUT")
      maxBinLength = -Inf
      if (sum(segAltInd) > 0){
        maxInd <- which.max(segsS$end[segAltInd] - segsS$start[segAltInd] + 1)
        maxSegRD <- GRanges(seqnames=segsS$chr[segAltInd[maxInd]], 
                  ranges=IRanges(start=segsS$start[segAltInd[maxInd]], end=segsS$end[segAltInd[maxInd]]))
        hits <- findOverlaps(query=maxSegRD, subject=tumour_copy[[s]][valid, ])
        maxBinLength <- length(subjectHits(hits))
      }
      ## check if there are proportion of total bins altered 
      # if segment size smaller than minSegmentBins, but altFrac > altFracThreshold, then still estimate TF
      cnaS <- hmmResults.cor$cna[[s]]
      altInd <- cnaS[cnaS$chr %in% chrTrain, "event"] == "NEUT"
      altFrac <- sum(!altInd, na.rm=TRUE) / length(altInd)
      if ((maxBinLength <= minSegmentBins) & (altFrac <= altFracThreshold)){
        hmmResults.cor$results$n[s, iter] <- 1.0
      }

      # correct integer copy number based on estimated purity and ploidy
      correctedResults <- ichorCNA::correctIntegerCN(cn = hmmResults.cor$cna[[s]],
            segs = hmmResults.cor$results$segs[[s]], 
            purity = 1 - hmmResults.cor$results$n[s, iter], ploidy = hmmResults.cor$results$phi[s, iter],
            cellPrev = 1 - hmmResults.cor$results$sp[s, iter], 
            maxCNtoCorrect.autosomes = maxCN, maxCNtoCorrect.X = maxCN, minPurityToCorrect = minTumFracToCorrect, 
            gender = gender$gender, chrs = chrs)
      hmmResults.cor$results$segs[[s]] <- correctedResults$segs
      hmmResults.cor$cna[[s]] <- correctedResults$cn

        ## plot solution ##
      outPlotFile <- paste0(outDir, "/", id, "/", id, "_genomeWide_", "n", n, "-p", p)
      #mainName[counter] <- paste0(id, ", n: ", n, ", p: ", p, ", log likelihood: ", 
      #  signif(hmmResults.cor$results$loglik[hmmResults.cor$results$iter], digits = 4))
      #plotGWSolution(hmmResults.cor, s=s, outPlotFile=outPlotFile, plotFileType='pdf', 
      #      logR.column = "logR", call.column = "Corrected_Call",
      #       estimateScPrevalence=estimateScPrevalence, main=mainName[counter])
    }
    iter <- hmmResults.cor$results$iter
    results[[counter]] <- hmmResults.cor
    loglik[counter, "loglik"] <- signif(hmmResults.cor$results$loglik[iter], digits = 4)
    subClonalBinCount <- unlist(lapply(hmmResults.cor$cna, function(x){ sum(x$subclone.status) }))
    fracGenomeSub <- subClonalBinCount / unlist(lapply(hmmResults.cor$cna, function(x){ nrow(x) }))
    fracAltSub <- subClonalBinCount / unlist(lapply(hmmResults.cor$cna, function(x){ sum(x$copy.number != 2) }))
    fracAltSub <- lapply(fracAltSub, function(x){if (is.na(x)){0}else{x}})
    loglik[counter, "Frac_genome_subclonal"] <- paste0(signif(fracGenomeSub, digits=2), collapse=",")
    loglik[counter, "Frac_CNA_subclonal"] <- paste0(signif(as.numeric(fracAltSub), digits=2), collapse=",")
    loglik[counter, "init"] <- paste0("n", n, "-p", p)
    loglik[counter, "n_est"] <- paste(signif(hmmResults.cor$results$n[, iter], digits = 2), collapse = ",")
    loglik[counter, "phi_est"] <- paste(signif(hmmResults.cor$results$phi[, iter], digits = 4), collapse = ",")

    counter <- counter + 1
  }
}
## get total time for all solutions ##
elapsedTimeSolutions <- proc.time() - ptmTotalSolutions
message("Total ULP-WGS HMM Runtime: ", format(elapsedTimeSolutions[3] / 60, digits = 2), " min.")

### SAVE R IMAGE ###
#save.image(outImage)
#save(tumour_copy, results, loglik, file=paste0(outDir,"/",id,".RData"))

### SELECT SOLUTION WITH LARGEST LIKELIHOOD ###
loglik <- loglik[!is.na(loglik$init), ]
if (estimateScPrevalence){ ## sort but excluding solutions with too large % subclonal 
  fracInd <- which(loglik[, "Frac_CNA_subclonal"] <= maxFracCNASubclone & 
                   loglik[, "Frac_genome_subclonal"] <= maxFracGenomeSubclone)
  if (length(fracInd) > 0){ ## if there is a solution satisfying % subclonal
    ind <- fracInd[order(loglik[fracInd, "loglik"], decreasing=TRUE)]
  }else{ # otherwise just take largest likelihood
    ind <- order(as.numeric(loglik[, "loglik"]), decreasing=TRUE) 
  }
}else{#sort by likelihood only
  ind <- order(as.numeric(loglik[, "loglik"]), decreasing=TRUE) 
}
write.table(loglik, paste0(outDir,'/loglik.csv'), col.names=T, row.names=F, sep='\t', quote=F)
#new loop by order of solutions (ind)

hmmResults.cor <- results[[ind[1]]]
hmmResults.cor$results$loglik <- as.data.frame(loglik)
hmmResults.cor$results$gender <- gender$gender
hmmResults.cor$results$chrYCov <- gender$chrYCovRatio
hmmResults.cor$results$chrXMedian <- gender$chrXMedian
hmmResults.cor$results$coverage <- coverage

outputHMM(cna = hmmResults.cor$cna, segs = hmmResults.cor$results$segs, 
                      results = hmmResults.cor$results, outDir=outDir)
outFile <- paste0(outDir, "/params.txt")
outputParametersToFile(hmmResults.cor, file = outFile)

saveRDS(hmmResults.cor, outfile)

## plot solutions for all samples 
plotSolutions(hmmResults.cor, tumour_copy, chrs, outDir, numSamples=numSamples,
              logR.column = "logR", call.column = "Corrected_Call",
              plotFileType='pdf', #plotYLim=plotYLim, seqinfo = seqinfo,
              estimateScPrevalence=estimateScPrevalence, maxCN=maxCN)








