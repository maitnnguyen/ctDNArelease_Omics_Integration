#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// --------------------------------------------------
// LOAD SAMPLE METADATA
// --------------------------------------------------
Channel
    .fromPath(params.sample_info)
    .splitCsv(header: true, sep: '\t')   // header = true, so we can use column names
    .filter { row -> row.normal == 'TRUE' && row.usable == 'TRUE' && row.normalSample != 'NA' }
    .map { row -> tuple(row.normalSample, file(row.normalBamFile)) }
    .distinct()
    .set { bam_channel }

// --------------------------------------------------
// MODULES
// --------------------------------------------------
include { coverage_calc } from '../modules/coverage_calc'
include { create_PoN    } from '../modules/create_PoN'

// --------------------------------------------------
// STEP1: DOWNSAMPLING FROM BDNA GENOME DATA THEN CALCULATE COVERAGE
// --------------------------------------------------
process DOWNSAMPLE_AND_COVERAGE {

    publishDir "${params.out_dir}/coverage", mode: 'copy'
    tag { sample }

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${sample}_coverage.txt")

    script:
    """
    # Downsample BAM (local to work dir)
    samtools view -s ${params.downsampling_rate} -b $bam > ${sample}.downsampled.bam
    samtools index ${sample}.downsampled.bam

    # Calculate coverage
    Rscript -e "
    suppressWarnings(suppressMessages({
        library(PureCN)
        library(dplyr)
    }))

    calculateBamCoverageByInterval(
        bam.file = '${sample}.downsampled.bam',
        interval.file = '${params.interval_file}',
        output.file = '${sample}_coverage.txt'
    )
    "

    # Optional cleanup (not strictly needed)
    rm -f ${sample}.downsampled.bam ${sample}.downsampled.bam.bai
    """
}   

// --------------------------------------------------
// MAIN WORKFLOW
// --------------------------------------------------
workflow create_PoN_workflow {

    take:
      bam_channel
      gc_map
      centromere

    main:
      // 1. Coverage calculation (your existing module) 
      // downsampled bam files are removed to save space
      cov = DOWNSAMPLE_AND_COVERAGE(bam_channel)
        .map { sample, covfile -> covfile }
        .collect()

      // 2. Create PoN
      pon = create_PoN(cov, gc_map, centromere)

    emit:
      pon
}

workflow {

    create_PoN_workflow(
        bam_channel,
        file(params.gc_map_file),
        file(params.centromere_file)
    )
}
