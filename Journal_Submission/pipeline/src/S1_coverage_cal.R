#!/usr/bin/env Rscript

# 1 bam.file, 2 interval.file, 3 output
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 2) {
  stop("Check input, need at least 2 inputs: bam.file, interval.file, and output (optional)", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  args[3] = "out.txt"
}

suppressWarnings(suppressMessages({
	library(PureCN) 
	library(dplyr)
	library(stringr)
}))

## running command
calculateBamCoverageByInterval(bam.file = args[1],
                                 interval.file = args[2], 
                                 output.file = args[3]#,
                                 #index.file = args[1],#gsub('bam','bai',args[1]), 
                                 #keep.duplicates = FALSE
                                 )