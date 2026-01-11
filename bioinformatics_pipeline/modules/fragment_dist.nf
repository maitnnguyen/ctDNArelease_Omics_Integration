process fragment_distribution {
    publishDir "${params.out_dir}/fragment_distribution", mode: 'copy'
    tag { sample }

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    	tuple val(sample), 
	path("${sample}_fragment_distribution.txt"), 
	path("${sample}_distribution_histogram.pdf") // emit sample name + fragment distribution file

    // note: bam file seq style is chrM, chr1, chr2, ..., chrX, chrY
    // please adjust if needed (style 1,2,...,X,Y)

    script:
    """
    samtools view -b ${bam} chr{1..22} chr{X,Y} > temp.bam

    java -jar /opt/share/picard-2.20.2/picard.jar CollectInsertSizeMetrics \
        I=temp.bam \
        O=${sample}_fragment_distribution.txt \
        H=${sample}_distribution_histogram.pdf

    rm temp.bam
    """
}
