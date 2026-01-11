#!/usr/bin/env nextflow

// CNA and tumor fraction estimation from plasma sWGS data
// Author: Mai T.N. Nguyen

// Input (mandatory):
// --sample_info: tsv file with tumor and normal labels and paths to bam files
// --pubDir: output directory

// TODO: update documentation also inside the workflow

nextflow.enable.dsl=2
 
// list of processes to be included
include { coverage_calc } from '../../modules/coverage_calc.nf'
include { segment_sample } from '../../modules/segment_sample.nf'
include { fragment_distribution } from '../../modules/fragment_dist.nf'

workflow {
    // Load BAM files from the specified directory
    bam_files = Channel.fromPath("${params.bam_dir}/*.bam")
                     .map { 
                        bam -> 
                        tuple(bam.baseName, bam, 
                        file("${bam.parent}/${bam.baseName}.bai")) }

    // option: check the DNA fragment size distribution
    fragment_distribution(bam_files)

    coverage = coverage_calc(bam_files)

    segment_results = segment_sample(coverage)

    // Aggregate parameters into separate TSVs per sample
    param_files = aggregate_params(segment_results)

    //param_files.view()
    // Merge all sample TSVs into a single summary
    merge_summary(param_files.collect())
}

process aggregate_params {
    publishDir "${params.out_dir}/summary", mode: 'copy'

    input:
    tuple val(sample), path(segment_dir)

    output:
    path "param_summary_${sample}.tsv"

    script:
    """
    Rscript -e '
      sample <- "${sample}"
      dir <- "${segment_dir}"

      param_file <- paste0(dir, "/params.txt")
      p <- read.delim(param_file)
      p[2,1:3]
      res <- data.frame(sample=sample, purity=p[1,2], ploidy=p[1,3])

      write.table(res, "param_summary_${sample}.tsv", 
                col.names = T, row.names = F, quote = F, sep = "\\t")
    '
    """
}

process merge_summary {
    publishDir "${params.out_dir}/summary", mode: 'copy'
    input:
    path summary_files

    output:
    path "params_summary.tsv"

    script:
    """
    awk 'FNR==1 && NR!=1{next} {print}' ${summary_files} > params_summary.tsv
    """
}
