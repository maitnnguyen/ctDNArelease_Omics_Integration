
# Multi-Modal Characterization of ctDNA Release Mechanisms in Ovarian Cancer Associate ctDNA Levels to Treatment Response

## Back ground:
- High grade serous carcinoma (HGSC) is usually diagnosed at late stages but from our data, we observed that there is considerable proportion of HGSC patients showed very low ctDNA fraction. Thus, we integrated genomics and transcriptomics information from tissue biopsies to understand biological features behind ctDNA level phenotypes.

## Material and method
- Liquid biopsy: plasma samples
  - shallow WGS
  - tumor fraction estimation with ichorCNA
- Genomics:
  - mutational burden
  - mutational signatures
  - number of breakpoints from segmentation profiles
- Transcriptomics:
  - cancer cell component decomposed from bulkRNA expression
  - DESeq2 for differential expression gene
  - Gene set enrichment analysis with gsea
  - over representation analysis with pathfindR
- Single-cell RNA:
  - Manual annotation of 3 main cell types: epithelial, fibroblast, and immune
  - After that, immune cells were extracted and annotated with celldex package, using 2 reference sets: MonacoImmuneData and DatabaseImmuneCellExpressionData
 
## Source code
- Pipeline for data from plasma samples
  - ~/Journal_Submission/pipeline
- Analysis and figures
  - ~/Journal_Submission/analysis
 
## Graphical abstracts
[TBA]
