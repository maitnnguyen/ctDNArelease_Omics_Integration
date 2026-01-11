# Copy number calling pipeline for ctDNA sWGS

This project directory contains pipelines for CNA calling from ctDNA sWGS samples.

### Main Workflow  
Pipeline info:
- Name: main.nf
- Config file: nextflow.config

Workflow includes:
- Calculate coverage per bin size (1Mb). Of note, the bin size depends on the coverage of the data, the lower coverage, the longer bin size. Bin size needs to be annotated for GC content and mappability, which can be done with GATK tools. 
- Coverage then was corrected for GC and mappability values with HMM package.
- Panel of normals (PoNs) was then used for normalizing as well as minimizing noises.
    - In this specific pipeline, PoNs was built from genomic WGS while-blood control samples.
    - Public has other PoNs built mostly from simulated healthy cfDNA donor samples.
- Segmentation was done with *ichorCNA* and the results returned including copy-number profiles (segmentation file) and sample parameters estimation (sample's tumor fraction and ploidy).

### Subworkflow
There is the pipeline for generating PoN from WGS BDNA samples by downsampling to similar depth of coverage of sWGS data. This pipeline costs heavy computational power due to running `samtools` for downsampling deep coverage.
- Name: `pon.nf`
- Input file is configed in the file name `pon.config`.
PoN has been built from 314 BDNA samples, no need to run this pipeline if not needed. PoN file is deposited into `asset` folder.

### Key Features
- Input files: bam file produced from Pre-processing pipeline (DNAseq-preprocessing)
- Helper files: all are located in `assets/` director
	- Interval files (regions of interest for calculating coverage) with GC and mappability annotation
	- Centromere
	- Previous constructed PoN: `assets/bdna_hg38_wgs_downsample_1M_median_PoN.rds`
	- PoN from ichorCNA (simulated from healthy cfDNA samples) => this is for reference if needed. If this is used for the pipelinem update config file specifying PoN: `HD_ULP_PoN_hg38_1Mb_median_normAutosome_median.rds`
- Centralized configuration: Most of parameters are default except mappability threshold for filtering (recommended for ctDNA sWGS data). This parameter is defined in `params.json`. If the user wants to apply the pipeine for deeper coverage sequencing (≥ 2x, try to tune parameter - not implementeed in this pipeline).
- R script for modules/(sub)workflows: located in `bin/` directory

*Additional module (**recommended**) *:
- Module `fragment_dist.nf` is recommended to run and check how DNA fragment size of the ctDNA sample, this will help evaluate quality of the sequencing. The peak of DNA fragment size should be around 140-170bp, if the peak of the DNA fragment size is bigger than 200bp (look at the distribution histogram).
- This module is added into the main pipeline, can be commented out if not needed.

### Configuration
In the configuration file for each pipeline, adjust information if needed such as:
- Path to bam files to run the main pipeline. Usually, after the NGS data is processed and the output files (bam files) are moved to `/mnt/storageBig8/resources/processed_data/ctDNA/[folder contain bam files]`. This path should be updated when running the pipeline.
- PoN file: the reference file for the pipeline.
	- Not recommended to rerun the pipeline unless needed. 
	- Other PoN can be used is PoN provided by ichorCNA that simulated from healthy cfDNA samples. This file is deposited in `assets` folder, file name: ``
- Other configurations: `output folder` and some parameters for the pipeline, for example *mappability score* for filtering genomic regions that are not reliable for identifying copy-number segmentation.
	- The pipeline for now is using default parameters. If the user is interested in investigating other parameters, please check script `bin/segment_sample.R` or run `Rscript bin/segment_sample.R --help` to see parameters used for the algorithm.

### Usage
#### Basic Usage
Access to folder `workflow/swgs/` then run command below:

```/opt/bin/share/nextflow run main.nf```

### Scripts
All scripts are deposited in folder `bin` in which, folder `R` contains backbone script for the algorithm. Scripts in `bin/R` are modified from ichorCNA (https://github.com/broadinstitute/ichorCNA/tree/master/R).

### Output Structure
The main aim of the pipeline is to estimate the tumor fraction of ctDNA sample. The summary file in `summary` folder collect purity and ploidy parameters from all ctDNA samples - `params_summary.tsv` while the file name with sample name at the end is for each sample. Output folder structure is below:

```
results/
├── coverage
│   ├── {sample}_coverage.txt
├── fragment_distribution
│   ├── {sample}_fragment_distribution.txt
│   ├── {sample}_distrbution_histogram.pdf
├── segmentation
│   ├── {sample}_segment
│   │   ├── normalizeCount.rds
│   │   ├── params.txt
│   │   ├── sample.cna.seg
│   │   ├── sample.seg
│   │   └── sample.seg.txt
└── summary
    ├── params_summary.tsv
    └── param_summary_{sample}.tsv
```

