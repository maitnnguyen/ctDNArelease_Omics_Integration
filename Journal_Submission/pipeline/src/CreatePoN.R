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
  make_option(c("-f", "--filelist"), type = "character", help = "List of of corrected bias read count file."),
  make_option(c("-o", "--outfile"), type = "character", help = "Output file."),
  make_option(c("-g", "--gcmap"), type="character", help = "File containing GC and mappability info for correction"),
  make_option(c("-c", "--centromere"), type="character", help = "File containing Centromere locations"),
  make_option(c("--chrs"), type="character", default="c(1:22,\"X\")", help = "Specify chromosomes to analyze."),
  make_option(c("--chrNormalize"), type="character", default="c(1:22)", help = "Specify chromosomes to normalize GC/mappability biases"),
  make_option(c("--maleChrXLogRThres"), type="numeric", default=-0.80, help = "ChrX Log ratio threshold to confirm as male gender."),
  make_option(c("-r", "--hitologyFile"), type = "character", default=NULL, help = "File from eeduni to annotate histology of sample"),
  make_option(c("--tp53"), type = "character", default=NULL, help = "File for mutation calling from mutation pipeline"),
  make_option(c("--method"), type = "character", default="median", help="Median or Mean."),
  make_option(c("--pon"), type = "character", default='BDNA', help="Using BDNA or 0 TP53 plasma as PoN")
)
parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
print(opt)

# source for supporting import coverage file
source("~/src/R/utils.R")

# infor of input
filelist <- opt$filelist
histologyFile <- opt$hitologyFile
tp53File <- opt$tp53
gc.map.file <- opt$gcmap
centromereFile <- opt$centromere
method <- opt$method
pon_method <- opt$pon
outfile <- opt$outfile
maleChrXLogRThres <- opt$maleChrXLogRThres
chrs <- as.character(eval(parse(text = opt$chrs)))
chrNormalize <- as.character(eval(parse(text=opt$chrNormalize))); 

# gc and mappability info and centromere
gc_map = read.delim(gc.map.file) %>%
    dplyr::rename(Target = Targets)

centromere = read.delim(centromereFile)

# proceeding
if (is.null(histologyFile)){
  if (pon_method == 'BDNA'){
        files = read.delim(filelist) %>%
        dplyr::rename(SampleID = Key) %>%
        filter(grepl('BDNA', SampleID))
      }
  else {
      tp53 = read.delim(tp53File) %>%
        dplyr::rename(SampleID = sample)

      files = read.delim(filelist) %>%
        dplyr::rename(SampleID = Key) %>%
        left_join(tp53, by = c('SampleID')) %>%
        filter(vaf <= 0.001, grepl('TPL', SampleID))
    }
  } else {
   files = read.delim(filelist) %>%
   dplyr::rename(SampleID = Key) %>%
   inner_join(read.delim(histologyFile) %>%
    filter(Histology=='benign'), by = 'SampleID') 
  }


normalGR <- NULL

for (i in 1:dim(files)[1]){
  sid <- files[i,'SampleID']
  rawCounts <- loadReadCount(files[i,'File'], map.gc.AnotFile = gc.map.file,
                            centromere = centromereFile) #readRDS(files[i])
  normCounts <- correctReadCounts(rawCounts, chrNormalize=chrNormalize)
  # get gender
  gender <- extractGender(rawReads = rawCounts, normReads = normCounts)
  message("Correcting ", sid)
  if (is.null(normalGR)){
    normalGR <- normCounts
    values(normalGR) <- values(normalGR)$copy
    colnames(values(normalGR)) <- sid
  }else{
    values(normalGR)[[sid]] <- normCounts$copy
  }
  
  chrXMedian <- gender$chrXMedian
  chrXStr <- grep("X", chrs, value = TRUE)
  chrXInd <- as.character(seqnames(normalGR)) == chrXStr
  ## Normalize chrX ##
  values(normalGR)[[sid]][chrXInd] <- values(normalGR)[[sid]][chrXInd] - chrXMedian
}

mat <- values(normalGR)
if (method == "median"){
   medianVal <- apply(mat, 1, median, na.rm = TRUE)
}else if (method == "mean"){
   medianVal <- apply(mat, 1, mean, na.rm = TRUE)
}else{
  stop("method is not specified as median or mean.")
}
values(normalGR)[["Median"]] <- medianVal

## output to text file ##
write.table(as.data.frame(normalGR[,"Median"]), file=paste0(outfile, "_", method, ".txt"), col.names=TRUE, row.names=F, quote=F, sep="\t")

## save GR object ##
saveRDS(normalGR, file = outfile)
