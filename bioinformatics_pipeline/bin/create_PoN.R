#!/usr/bin/env Rscript

suppressMessages({
  library(GenomicRanges)
  library(optparse)
  library(ichorCNA)
  library(PureCN)
  library(dplyr)
})

options(stringsAsFactors = FALSE, scipen = 0)
options(bitmapType = "cairo")

option_list <- list(
  make_option(c("--bin_dir"), type="character", help="Path to bin folder"),
  make_option(c("-f", "--filelist"), type="character",
              help="File containing list of corrected bias count files OR comma-separated list."),
  make_option(c("-o", "--outfile"), type="character", help="Output prefix."),
  make_option(c("-g", "--gcmap"), type="character", help="GC/mappability annotation file."),
  make_option(c("-c", "--centromere"), type="character", help="Centromere locations file."),
  make_option(c("--chrs"), type="character", default="c(1:22,'X')",
              help="Chromosomes to analyze."),
  make_option(c("--chrNormalize"), type="character", default="c(1:22)",
              help="Chromosomes used for GC/mappability normalization."),
  make_option(c("--maleChrXLogRThres"), type="numeric", default=-0.80,
              help="ChrX LogR threshold for male detection."),
  make_option(c("-r", "--hitologyFile"), type="character", default=NULL,
              help="Histology annotation file."),
  make_option(c("--tp53"), type="character", default=NULL,
              help="TP53 mutation annotation file."),
  make_option(c("--method"), type="character", default="median",
              help="Aggregation method: median or mean."),
  make_option(c("--pon"), type="character", default="BDNA",
              help="BDNA or 0-TP53 plasma PoN")
)

parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
print(opt)

# -------------------------------------------------------------------------
# Load utility functions
# -------------------------------------------------------------------------
bin_dir <- opt$bin_dir
source(file.path(bin_dir, "R/utils.R"))

# -------------------------------------------------------------------------
# Parse parameters
# -------------------------------------------------------------------------
gc.map.file    <- opt$gcmap
centromereFile <- opt$centromere
histologyFile  <- opt$hitologyFile
tp53File       <- opt$tp53
method         <- opt$method
pon_method     <- opt$pon
outfile        <- opt$outfile
maleChrXLogRThres <- opt$maleChrXLogRThres

chrs <- as.character(eval(parse(text = opt$chrs)))
chrNormalize <- as.character(eval(parse(text = opt$chrNormalize)))


# -------------------------------------------------------------------------
# Parse filelist (file OR comma-separated)
# -------------------------------------------------------------------------
if (file.exists(opt$filelist)) {
    files_df <- read.delim(opt$filelist)
} else {
    files_df <- data.frame(File = strsplit(opt$filelist, ",")[[1]])
    files_df$Key <- basename(files_df$File)
}

files_df <- files_df %>%
    rename(SampleID = Key)


# -------------------------------------------------------------------------
# Select samples for PoN
# -------------------------------------------------------------------------
if (is.null(histologyFile)) {

    if (pon_method == "BDNA") {
        files <- files_df %>% filter(grepl("BDNA", SampleID))

    } else {
        tp53 <- read.delim(tp53File) %>% rename(SampleID = sample)

        files <- files_df %>%
            left_join(tp53, by = "SampleID") %>%
            filter(vaf <= 0.001 & grepl("TPL", SampleID))
    }

} else {

    histology <- read.delim(histologyFile) %>% filter(Histology %in% c("benign", "control", "healthy"))

    files <- files_df %>%
        inner_join(histology, by = "SampleID")
}

# -------------------------------------------------------------------------
# Load reference data
# -------------------------------------------------------------------------
gc_map <- read.delim(gc.map.file) %>% rename(Target = Targets)
centromere <- read.delim(centromereFile)

normalGR <- NULL

# -------------------------------------------------------------------------
# Process each sample
# -------------------------------------------------------------------------
for (i in seq_len(nrow(files))) {

    sid <- files$SampleID[i]

    rawCounts <- loadReadCount(
        files$File[i],
        map.gc.AnotFile = gc.map.file,
        centromere      = centromereFile
    )

    normCounts <- correctReadCounts(rawCounts, chrNormalize = chrNormalize)
    gender <- extractGender(rawReads = rawCounts, normReads = normCounts)

    message("Correcting ", sid)

    if (is.null(normalGR)) {
        normalGR <- normCounts
        values(normalGR) <- normCounts$copy
        colnames(values(normalGR)) <- sid
    } else {
        values(normalGR)[[sid]] <- normCounts$copy
    }

    chrXMedian <- gender$chrXMedian
    chrXInd <- seqnames(normalGR) == grep("X", chrs, value = TRUE)
    values(normalGR)[[sid]][chrXInd] <- values(normalGR)[[sid]][chrXInd] - chrXMedian
}

# -------------------------------------------------------------------------
# Compute PoN (median or mean)
# -------------------------------------------------------------------------
mat <- values(normalGR)

if (method == "median") {
    consensus <- apply(mat, 1, median, na.rm = TRUE)
} else if (method == "mean") {
    consensus <- apply(mat, 1, mean, na.rm = TRUE)
} else {
    stop("Method must be 'median' or 'mean'")
}

values(normalGR)[["Consensus"]] <- consensus

# -------------------------------------------------------------------------
# Write outputs
# -------------------------------------------------------------------------
write.table(
    as.data.frame(normalGR[, "Consensus"]),
    file = paste0(outfile, "_", method, ".txt"),
    col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t"
)

saveRDS(normalGR, file = outfile)

