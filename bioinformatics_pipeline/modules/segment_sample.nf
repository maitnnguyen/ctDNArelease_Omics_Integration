
process segment_sample {
    tag { sample }
    publishDir "${params.out_dir}/segmentation", mode: 'copy'

    input:
    tuple val(sample), path(covfile)

    output:
    tuple val(sample), path("${sample}_segment") , emit: segment_dir

    script:
    """
    mkdir -p ${sample}_segment
    Rscript ${params.bin_dir}/segment_sample.R \
        --bin_dir ${params.bin_dir} \
        --outDir ${sample}_segment \
        --input ${covfile} \
        --PoN ${params.pon_file} \
        --gcmap ${params.gc_map_file} \
        --centromere ${params.centromere_file} \
        --mapscore ${params.mapscore}
    """
}