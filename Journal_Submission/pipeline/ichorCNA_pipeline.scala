#!/usr/bin/env anduril
//$OPT --threads 20
//$OPT --pipe "tee $ANDURIL_EXECUTION_DIR/_log"
//$OPT --pipe "anduril-pager --ls"
//$OPT --wrapper slurm-prefix
//$PRE export ANDURIL_SELECTNODE="mintaken"
//$PRE export ANDURIL_NODELIST="evm01 evm02 evm03 evm04 evm05 evm06 evm07 evm08 evm09 evm10 evmfull01 evmfull02"
//$PRE export PATH="/mnt/csc-gc8/work/joikkone/miniconda2/bin:/mnt/csc-gc8/work/joikkone/miniconda2/envs/snakemake/bin:$PATH"

import anduril.builtin._
import anduril.tools._
import org.anduril.runtime._
import anduril.microarray._
import anduril.sequencing._

object ichorCNA {
	val bam_files = INPUT(path="~/result_alignment/outBamCSV/out.csv")
    val intervals = INPUT(path="~/int/binSizeInterval/interval_1Mb.intervals")
    val interval_GCannot = INPUT(path = "~/int/binSizeInterval/gcmap_1Mb.intervals")
    val centromere = INPUT(path ="~/int/GRCh38_centromere_acen.txt")

    // WGS normal bams
    val WGSbams = INPUT(path = "~/[WGS normal bam files]")
    val DownBDNASample = NamedMap[BashEvaluate]("DownBDNASample")
    val DownSampleArray = NamedMap[BinaryFile]("DownSampleArray")

    // step1: calculate raw coverage from bamfiles
    val Coverage = NamedMap[BashEvaluate]("Coverage")
    val CoverageArray = NamedMap[BinaryFile]("CoverageArray")
    val InsertMatric = NamedMap[BashEvaluate]("InsertMatric")
    
    // step2: correct GC and mappability and normalizePON
    val normalizePON = NamedMap[BashEvaluate]("normalizePON")
    // step3: making PoN
    val CNVcall = NamedMap[BashEvaluate]("CNVcall")
    val LogRArray = NamedMap[BinaryFile]("LogRArray")
    val SegmArray = NamedMap[BinaryFile]("SegmArray")
    val segment = NamedMap[BashEvaluate]("segment")
    
    // PoN from previous run set1
    val PoN = INPUT(path = "/makePoNfromBDNA/normDB.rds")


    // downsampling WGS normal samples
    for ( row <- iterCSV(WGSbams)) {
            //val patient = row("patient")
            val bamfile = INPUT(row("File"))
            //val vcffile = row("vcf")
            val sample = row("sample")

            DownBDNASample(sample) = BashEvaluate(var1 = bamfile,
                                         script = """
                                         samtools view -bs 50.01 @var1@ > @out1@
                                         samtools index -b @out1@

                                         """
            )
                                         
            DownBDNASample(sample)._filename("out1","out1.bam")
            //DownSampling(sample)._filename("out2","out2.bam")
            DownBDNASample(sample)._custom("memory") = "5G"
            DownBDNASample(sample)._custom("cpu") = "5"

            DownSampleArray(sample) = DownBDNASample(sample).out1
    } 
    val DownSampleList = Array2CSV(in = DownSampleArray)

    // combine sWGS and downsampling WGS normal samples
    val bamfiles = REvaluate(table1 = bam_files,
                            table2 = DownSampleList.out,
                            script = """
                                suppressMessages({
                                    library(dplyr)
                                })
                                table1 = table1 %>%
                                    dplyr::rename(sample=Key) %>% 
                                    mutate(patient=sub('_.*', '', sample))
 
                                table2 <- table2 %>%
                                    dplyr::rename(sample=Key) %>% 
                                    mutate(patient=sub('_.*', '', sample))

                                table.out = data.frame(rbind.data.frame(table1, table2))
                                """)

    // Step 1: calculate coverage by bin size       
    for ( row <- iterCSV(bam_files)) {
            val bamfile = INPUT(row("File"))
            val sample = row("sample")

            Coverage(sample) = BashEvaluate(var1 = bamfile,
                                         var2 = intervals,
                                         script = """
                                         Rscript /src/S1_coverage_cal.R @var1@ @var2@ @out1@
                                         """
            )
            Coverage(sample)._filename("out1","out1.txt")
            Coverage(sample)._custom("memory") = "4G"
            CoverageArray(sample) = Coverage(sample).out1
            
    }
    val CoverageList = Array2CSV(in = CoverageArray)

    val covFileByPatient = REvaluate(table1 = CoverageList, //Array2CSV(in = Coverage(sample).out1),
                            script = """
                                suppressMessages({
                                    library(dplyr)
                                })

                                table1 <- table1 %>% 
                                        mutate(patient = sub('_.*','',Key)) %>%
                                        dplyr::rename(sample = Key) %>%
                                        filter(grepl('TPL', sample))

                                table.out = table1
                                array.out = split(table.out[,c("sample","File")],table.out$patient)
                                """)
    val covFileByPatientList = Array2CSV(covFileByPatient.outArray)

    // Step 2: making Panel of Normal 
    val makePoNfromBDNA = BashEvaluate(var1 = CoverageList.out,
                                var2 = interval_GCannot,
                                var3 = centromere,
                                script = """
                                Rscript ~/src/CreatePoN.R \
                                -f @var1@ -g @var2@ -c @var3@ -o @out1@ 

                                """)
    makePoNfromBDNA._filename("out1","normDB.rds")
    makePoNfromBDNA._custom("memory") = "3G"

    // Step 3: segmentation with ichorCNA
    for ( row <- iterCSV(covFileByPatient.table)){
        val sample = row("sample")
        val covfile = INPUT(row("File"))
        
        segment(sample) = BashEvaluate(var1 = covfile,
                                        var2 = makePoNfromBDNA.out1,
                                        var3 = interval_GCannot,
                                        var4 = centromere,
                                        script = """
                                            Rscript /src/segmentforSample.R \
                                            -i @var1@ -g @var3@ -c @var4@ -p @var2@ -o @out1@ -O @folder1@
                                            """)
        segment(sample)._filename("out1","out1.rds")
        segment(sample)._custom("memory") = "3G"
        SegmArray(sample) = segment(sample).out1

    }

    val Segmlist = Array2CSV(in = SegmArray)
    
}
