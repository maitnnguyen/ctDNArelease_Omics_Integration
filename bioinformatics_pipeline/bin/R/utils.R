suppressMessages({
  library(dplyr)
  library(PureCN)
  library(HMMcopy)
  library(GenomicRanges)
})

# centromere info for excluding
excludeCentromere <- function(x, centromere, flankLength = 0){
  require(GenomeInfoDb)
  colnames(centromere)[1:3] <- c("seqnames","start","end")
  centromere$start <- centromere$start - flankLength
  centromere$end <- centromere$end + flankLength
  centromere <- as(centromere, "GRanges")
  #seqlevelsStyle(centromere) <- genomeStyle
  centromere <- sort(centromere)	
  hits <- findOverlaps(query = x, subject = centromere)
  ind <- queryHits(hits)
  message("Removed ", length(ind), " bins near centromeres.")
  if (length(ind) > 0){
    x <- x[-ind, ]
  }
  return(x)
}

# function to load read counts from PureCN
# exclude centromere region
# and then correct for gc and mappability

loadReadCount <- function(rawReadCountFile, map.gc.AnotFile, map_thresh = 0.85, centromereFile, count_thresh = 0){
  message(paste0("Removed regions with mappability lower than ", map_thresh))
  # import mappability and gc content annotation
  map_gc <- read.delim(map.gc.AnotFile) %>%
    dplyr::rename(Target = Targets)
  # centromore
  centromere = read.delim(centromereFile) %>%
    mutate(Chr=sub('chr','',Chr))
  # import readcound and add mappability and gc
  rawReadCount = read.delim(rawReadCountFile, sep = ' ') %>%
    mutate(chr = sub('chr','',sub(':.*','',Target)),
           ranges = sub('.*:','', Target)) %>%
    left_join(map_gc) %>%
    # filter chr in 1:22,X,Y only
    filter(chr %in% c(seq(1:22),'X','Y'),
           mappability>=map_thresh,
           counts >= count_thresh)
  
  readCount <- GenomicRanges::GRanges(seqnames = rawReadCount$chr,
                                      ranges= rawReadCount$ranges,
                                      reads = rawReadCount$counts,
                                      map = rawReadCount$mappability,
                                      gc = rawReadCount$gc_bias)
  out <- excludeCentromere(readCount, centromere)
  return(out)
}

# normalize readcount, function copy from ichorCNA
# https://github.com/broadinstitute/ichorCNA/blob/master/R/utils.R
correctReadCounts <- function(x, chrNormalize = c(1:22), mappability = 0.85, samplesize = 500000,
                              verbose = TRUE) {
  if (length(x$reads) == 0 | length(x$gc) == 0) {
    stop("Missing one of required columns: reads, gc")
  }
  chrInd <- sub('chr','',as.character(seqnames(x))) %in% chrNormalize
  if(verbose) { message("Applying filter on data...") }
  x$valid <- TRUE
  x$valid[x$reads <= 0 | x$gc < 0] <- FALSE
  x$ideal <- TRUE
  routlier <- 0.01
  range <- quantile(x$reads[x$valid & chrInd], prob = c(0, 1 - routlier), na.rm = TRUE)
  doutlier <- 0.001
  domain <- quantile(x$gc[x$valid & chrInd], prob = c(doutlier, 1 - doutlier), na.rm = TRUE)
  if (length(x$map) != 0) {
    x$ideal[!x$valid | x$map < mappability | x$reads <= range[1] |
              x$reads > range[2] | x$gc < domain[1] | x$gc > domain[2]] <- FALSE
  } else {
    x$ideal[!x$valid | x$reads <= range[1] |
              x$reads > range[2] | x$gc < domain[1] | x$gc > domain[2]] <- FALSE
  }
  
  if (verbose) { message("Correcting for GC bias...") }
  set <- which(x$ideal & chrInd)
  select <- sample(set, min(length(set), samplesize))
  rough = loess(x$reads[select] ~ x$gc[select], span = 0.03)
  i <- seq(0, 1, by = 0.001)
  final = loess(predict(rough, i) ~ i, span = 0.3)
  x$cor.gc <- x$reads / predict(final, x$gc)
  
  if (length(x$map) != 0) {
    if (verbose) { message("Correcting for mappability bias...") }
    coutlier <- 0.01
    range <- quantile(x$cor.gc[which(x$valid & chrInd)], prob = c(0, 1 - coutlier), na.rm = TRUE)
    set <- which(x$cor.gc < range[2] & chrInd)
    select <- sample(set, min(length(set), samplesize))
    final = approxfun(lowess(x$map[select], x$cor.gc[select]))
    x$cor.map <- x$cor.gc / final(x$map)
  } else {
    x$cor.map <- x$cor.gc
  }
  x$copy <- x$cor.map
  x$copy[x$copy <= 0] = NA
  x$copy <- log(x$copy, 2)
  return(x)
}

# get Gender from sample, copy from ichorCNA
extractGender <- function(rawReads, normReads, fracReadsInChrYForMale = 0.002, chrXMedianForMale = -0.5, useChrY = TRUE,
                      centromere=NULL, flankLength=1e5, targetedSequences=NULL, genomeStyle="NCBI"){
  chrXStr <- grep("X", runValue(seqnames(normReads)), value = TRUE)
  chrYStr <- grep("Y", runValue(seqnames(rawReads)), value = TRUE)
  chrXInd <- as.character(seqnames(normReads)) == chrXStr
  chrYInd <- as.character(seqnames(rawReads)) == chrYStr
  if (sum(chrXInd) > 1){ ## if no X 
    chrXMedian <- median(normReads[chrXInd, ]$copy, na.rm = TRUE)
    # proportion of reads in chrY #
    tumY <- rawReads[chrYInd]
    chrYCov <- sum(tumY$reads) / sum(rawReads$reads)
    if (chrXMedian < chrXMedianForMale){
      if (useChrY && (chrYCov < fracReadsInChrYForMale)){ #trumps chrX if using chrY
        gender <- "female"  
      }else{
        gender <- "male" # satisfies decreased chrX log ratio and/or increased chrY coverage
      }
    }else{
      gender <- "female" # chrX is provided but does not satisfies male critera
    }
  }else{
    gender <- "unknown" # chrX is not provided
    chrYCov <- NA
    chrXMedian <- NULL
  }
  return(list(gender=gender, chrYCovRatio=chrYCov, chrXMedian=chrXMedian))
}


# create PoN
param_export <- function(res_path, server = F){
  res_files = list.files(path = res_path, pattern = 'segment_')
  samples = sub('segment_','',res_files)
  res_links = paste(res_path, res_files,'folder1/params.txt', sep = '/')
  out = data.frame(matrix(nrow=0, ncol = 3))
  colnames(out) = c('sample', 'TF', 'N') # sample, TF: purity, N: ploidy
  
  # run loop to collect parameters
  for (i in 1:length(samples)){
    df <- read.delim(res_links[i])#[1,2:3]
    out[i,] = c(samples[i], df[1,2:3])
  }
  out <- out %>%
    mutate(TF = as.numeric(TF), N = as.numeric(N))
  # return result
  return(out)
}
