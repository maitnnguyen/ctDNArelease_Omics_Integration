
# Multi-Omics Characterization of ctDNA Release Mechanisms in Ovarian Cancer

## Back ground:
- High grade serous carcinoma (HGSC) is usually diagnosed at late stages but from our data, we observed that there is considerable proportion of HGSC patients showed very low ctDNA fraction. Thus, we integrated genomics and transcriptomics information from tissue biopsies to understand biological features behind ctDNA level phenotypes.

## Material and method
- Liquid biopsy: plasma samples
  - shallow WGS
  - tumor fraction estimation with ichorCNA, the bioinformatic pipeline is in the subfolder `bioinformatics_pipeline`. The pipeline was written with `nextflow`.
- Genomics:
  - mutational burden
  - mutational signatures: COSMIC database
  - number of breakpoints from segmentation profiles
- Transcriptomics:
  - cancer cell component decomposed from bulkRNA expression
  - DESeq2 for differential expression gene
  - Gene set enrichment analysis with gsea
- Single-cell transcriptomics for validation:
  - General cell types were annotated based on this paper: https://www.biorxiv.org/content/10.1101/2025.06.13.659489v1
  - SingleR, celldex are then used for annotating subtypes of T cells
  - UCell is used to calculate gene signature scores
 
## Source code
- Pipeline for data from plasma samples
  - `~/bioinformatics_pipeline`
  - `nextflow` pipeline
- Analysis and figures
  - `~/Journal_Submission`
 
## Graphical abstracts
[TBA]
