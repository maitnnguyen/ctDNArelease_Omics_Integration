
process coverage_calc {
    publishDir "${params.out_dir}/coverage", mode: 'copy'
    tag { sample }

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    tuple val(sample), path("${sample}_coverage.txt")  // emit sample name + coverage file

    script:
    """
    #mkdir -p ${params.out_dir}/coverage/

    Rscript -e '
    suppressWarnings(suppressMessages({
        library(PureCN)
        library(dplyr)
    }))

    calculateBamCoverageByInterval(
        bam.file = "${bam}", 
        interval.file = "${params.interval_file}", 
        output.file = "${sample}_coverage.txt"
    )
    '

    """
}