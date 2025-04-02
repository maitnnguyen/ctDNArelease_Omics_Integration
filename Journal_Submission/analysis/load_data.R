source('/int/00_functions.R')

# color setting:
list_color <- c('#cc3333', '#339933', "#003399")

##################### load data
# 1. ctdna data
# Clinical data
ctdna_df_prim <- readxl::read_excel('data/Supplementary_Tables.xlsx', 
                                      sheet = 'S1_ctDNA_samples') %>%
  as.data.frame() 

# note: details of clinical data will be provided on request
# this information is not included in this source code
clinical_data <- read.delim('~/clinical_data.tsv')

ctdna_df_prim <- ctdna_df_prim |>
  left_join(clinical_data)

# threshold:
thresh1 <- max(ctdna_df_prim[ctdna_df_prim$histology=='benign',]$ctDNA_fraction)
thresh2 <- 0.065

ctdna_cohort <- ctdna_df_prim |>
  filter(histology!='benign') |>
  dplyr::rename(concentration = `DNA_concentration (ng/ml)`)

# 2. RNA data
# note: details of processed RNA data is not included in this repository
# data including meta data file + deconvoluted expression matrix of cancer cell
load('~/OmicsIntegration/data/processed/RNA/pretrt_RNA.RData')


# 3. WGS data
# note: details of processed RNA data is not included in this repository
genome_size <- 2900077904/1e6 # (Megabase)
# load WGS tissue mutation data
load('~/OmicsIntegration/data/WGS/mutation_call_4_6.RData')

# mutational signature on sample level
load('~/OmicsIntegration/data/WGS/mutational_signature_sample_lev_231005.RData')

# purity estimation based on WGS CNV profile with ASCAT
load('~/OmicsIntegration/data/WGS/ascat_results.RData')
