################################################
# Figure 1: cohort info
# author: Mai T.N. Nguyen
# email: mai.tn.nguyen@helsinki.fi/ ntnmai303@gmail.com
# git: https://github.com/maitnnguyen/ctDNArelease_Omics_Integration/
################################################
source('src/00_functions.R')
list_color <- c('#cc3333', '#339933', "#003399")

# load data
source('src/final/01_load_data.R')

#------------------------------- Figure 1a -------------------------------#
#### visualisation - ctdna TF between stages
comp1 <- list(c('Benign\n(n=10)','Early\n(n=14)'),
              c('Benign\n(n=10)','III\n(n=85)'),
              c('Benign\n(n=10)','IV\n(n=19)'),
              c('Early\n(n=14)','IV\n(n=19)'),
              c('III\n(n=85)','IV\n(n=19)')
)

fig1a <- cfdna_samples |>
  mutate(
    Stage = ifelse(is.na(FIGO_stage), 'Benign', 
                   sub('[A-C].*', '', FIGO_stage)),
    type = ifelse(is.na(FIGO_stage), 'Benign', 
                  gsub("[^ABC]", "", FIGO_stage)),
    figo = ifelse(Stage %in% c('I', 'II'), 'Early', Stage),
    ctdna_lev = case_when(
      TF <= thresh1 ~ "low",
      TF > thresh1 & TF <= thresh2 ~ "med",
      TF > thresh2 ~ "high"
    )
  ) |>
  ggplot(aes(x=factor(figo, 
                      labels = c('Benign\n(n=10)',
                                 'Early\n(n=14)',
                                 'III\n(n=85)',
                                 'IV\n(n=19)')),
             y = TF)) +
  geom_boxplot(width=.3,outlier.shape = NA,
               alpha=.4, show.legend = F) + 
  geom_jitter(aes(color=ctdna_lev), 
              size=3, position=position_jitter(0.2)) + 
  stat_compare_means(comparisons = comp1, size = 5,
                     label.y = c(.12,.16,.20,.24, .28) ) +
  geom_hline(yintercept = thresh1, color=list_color[1]) + 
  geom_hline(yintercept = thresh2, color=list_color[3]) + 
  scale_colour_manual(values = c("high" = list_color[3], 
                                 "med" = list_color[2],
                                 "low"=list_color[1])) +
  labs(x='Stages at diagnosis', y ='ctDNA fraction',
       color = 'ctDNA level') + 
  theme_classic() +
  theme(axis.text.x = element_text(size=14,  family = 'Arial', colour = 'black'),
        axis.title = element_text(size=16,  family = 'Arial', colour = 'black'),
        axis.text = element_text(size=14, family = 'Arial', colour = 'black'),
        legend.position = 'none', 
        legend.text = element_text(size=11, family = 'Arial', colour = 'black'),
        legend.title = element_text(size=13, family = 'Arial', colour = 'black'))

svg('results/figures/main/01_ctDNA_fraction_stratification.svg', width = 4, height = 4)
fig1a
dev.off()

#------------------------------- Figure 1b -------------------------------#
# pie plot between HRD or treatment administered
fig1b1 <- ctdna_cohort |>
  mutate(HRD = ifelse(is.na(HRD), 'NA', HRD)) %>%
  dplyr::group_by(HRD) %>%
  dplyr::summarise(n = n()) %>%
  ungroup() %>%
  mutate(prop = n/sum(n),
         label = paste0(n, " (", round(prop*100,0), "%)")) %>%
  as.data.frame() |>
  ggplot(aes(x='',y = prop, fill=HRD)) +
  geom_col(width=1, color='white') +
  geom_label(aes(label = label), size=5,
             position = position_stack(vjust = 0.5),
             show.legend = FALSE,
             color = "black") +
  coord_polar(theta = "y") +
  labs(x = "", y='') +
  theme_void() +
  scale_fill_manual(values = c("HRD" = "#ffa500", "HRP" = "#c39797", 'NA'='grey')) +
  theme(axis.text.x = element_blank(),#element_text(size=14,  family = 'Arial', colour = 'black'),
        axis.title = element_text(size=16,  family = 'Arial', colour = 'black'),
        axis.text = element_text(size=14, family = 'Arial', colour = 'black'),
        legend.position = 'none', 
        legend.text = element_text(size=11, family = 'Arial', colour = 'black'),
        legend.title = element_text(size=13, family = 'Arial', colour = 'black'))

fig1b2 <- ctdna_cohort |>
  filter(!is.na(HRD)) |>
  ggplot(aes(x = HRD, y = TF)) +
  geom_violin(aes(fill = factor(HRD, levels = c('HRD','HRP'))), 
              alpha = .8) +
  geom_boxplot(width=.1) +
  scale_fill_manual(values = c("HRD" = "#ffa500", "HRP" = "#c39797"))+
  stat_compare_means(label.x.npc = .35, size = 6, label.y.npc = .8,
                     label = '..p..') +
  theme_pubr()+
  labs(x = 'Homologous recombination', y = 'ctDNA fraction') +
  theme(axis.text.x = element_text(size=14,  family = 'Arial', colour = 'black'),
        axis.title = element_text(size=16,  family = 'Arial', colour = 'black'),
        axis.text = element_text(size=14, family = 'Arial', colour = 'black'),
        legend.position = 'none', 
        legend.text = element_text(size=11, family = 'Arial', colour = 'black'),
        legend.title = element_text(size=13, family = 'Arial', colour = 'black'))

fig1b3 <- ctdna_cohort |>
  #mutate(HRD = ifelse(is.na(HRD), 'NA', HRD)) %>%
  dplyr::group_by(treatment) %>%
  dplyr::summarise(n = n()) %>%
  ungroup() %>%
  mutate(prop = n/sum(n),
         label = paste0(n, " (", round(prop*100,0), "%)")) %>%
  as.data.frame() |>
  ggplot(aes(x='',y = prop, fill=factor(treatment, levels = c('PDS','NACT','Other')))) +
  geom_col(width=1, color='white') +
  geom_label(aes(label = label), size=5,
             position = position_stack(vjust = 0.5),
             show.legend = FALSE,
             color = "black") +
  coord_polar(theta = "y") +
  theme_void() +
  labs(x = '', y='') +
  scale_fill_manual(values = c("NACT" = "#744700", "PDS" = "#76a5af", 'Other'='gray')) +
  theme(axis.text.x = element_blank(),
        axis.title = element_text(size=16,  family = 'Arial', colour = 'black'),
        axis.text = element_text(size=14, family = 'Arial', colour = 'black'),
        legend.position = 'none', 
        legend.text = element_text(size=11, family = 'Arial', colour = 'black'),
        legend.title = element_text(size=13, family = 'Arial', colour = 'black'))

fig1b4 <- ctdna_cohort |>
  filter(treatment!='Other') |>
  ggplot(aes(x = treatment, y = TF)) +
  geom_violin(aes(fill = factor(treatment, levels = c('PDS','NACT','Other'))), 
              alpha = .8) +
  geom_boxplot(width=.1) +
  scale_fill_manual(values = c("NACT" = "#744700", "PDS" = "#76a5af", 'Other'='gray'))+
  stat_compare_means(label.x.npc = .35, size = 6, label.y.npc = .8,
                     label = '..p..') +
  theme_pubr()+
  labs(x = 'Treatment modality', y = '') +
  theme(axis.text.x = element_text(size=14,  family = 'Arial', colour = 'black'),
        axis.title = element_text(size=16,  family = 'Arial', colour = 'black'),
        axis.text = element_text(size=14, family = 'Arial', colour = 'black'),
        legend.position = 'none', 
        legend.text = element_text(size=11, family = 'Arial', colour = 'black'),
        legend.title = element_text(size=13, family = 'Arial', colour = 'black'))

svg('results/figures/main/1b_HRD_trtm_ctDNA.svg', height = 5, width = 11)
ggarrange(
  fig1a, ggarrange(fig1b1, fig1b2, fig1b3, fig1b4, ncol = 2, nrow = 2),
  ncol = 2, widths = c(.36,.64)
)
dev.off()

#------------------------------- Figure 1d -------------------------------#
meta_mat <- ctdna_cohort |>
  left_join(
    clinical_data |>
      dplyr::select(patient, eoc_code)
  ) |>
  mutate(tissue_DNA_avail = ifelse(patient %in% dna_cohort$patient, 'yes', 'no'),
         tissue_RNA_avail = ifelse(patient %in% rna_cohort$patient, 'yes', 'no'),
         stage = gsub("^(IV|III|II|I).*", "\\1", FIGO_stage),
         chronic_ill_atDg = tolower(chronic_ill_atDg),
         prev_cancer = tolower(prev_cancer)) |>
  dplyr::select(eoc_code, patient, TF, ctdna_lev, stage, bmi_gr, 
                age_gr, prev_cancer, treatment, 
                tissue_RNA_avail, tissue_DNA_avail, HRD) |>
  left_join(
    purple_est |>
      filter(sample %in% highest_purity_dna$sample) |>
      dplyr::select(patient, WGD) |>
      mutate(WGD = ifelse(WGD=='TRUE', 'yes', 'no'))
  ) |>
  tibble::column_to_rownames(var='eoc_code')

# CNV info: MYC, MECOM, KRAS, CCNE1, and manual curated RB1 
rb1_info <- readxl::read_excel('data/input/clinicalInfo/RB1_curated.xlsx',
                               sheet='calls') |>
  as.data.frame() |>
  dplyr::select(patient, RB1=RB1_class) |>
  mutate(RB1 = ifelse(RB1 == 'Loss', 'deletion', 
                      ifelse(RB1 == 'Normal', 'normal', 'loss')))

# adding CN gene status to the heatmap of the cohort
cnv_df_wide <- cna_df |>
  dplyr::rename(CN = CNA) |>
  mutate(status = ifelse(CN >= 8, 'amplification',
                         ifelse(CN > 4, 'gain',
                                ifelse(CN < .9, 'deletion',
                                       ifelse(CN < 1.3, 'loss', 'normal'))))) |>
  dplyr::select(patient, ctdna_lev, gene=gene_name, status) |>
  tidyr::spread(key=gene, value = status) |>
  inner_join(rb1_info)

meta <- meta_mat |>
  tibble::rownames_to_column(var='eoc') |>
  mutate(prev_cancer = ifelse(is.na(prev_cancer), 'no', prev_cancer)) |>
  left_join(cnv_df_wide |>
              dplyr::select(patient, MECOM, MYC, KRAS, CCNE1, MAGEB10, RB1)) |>
  dplyr::select(-patient) |>
  tibble::column_to_rownames(var='eoc')
# Convert all character columns to factors
meta <- as.data.frame(meta)
meta[] <- lapply(meta, function(x) {
  if (is.character(x)) factor(x) else x
})

# Create annotation heatmap using dummy matrix
{
  dummy_mat <- matrix(1, nrow = 1, ncol = nrow(meta))
  colnames(dummy_mat) <- rownames(meta)
  meta <- meta[colnames(dummy_mat), , drop = FALSE]
  colnames(meta) <- c(
    "TF" = "ctDNA fraction",
    "ctdna_lev" = "ctDNA level",
    "stage" = "Diagnosis stage",
    "bmi_gr" = "BMI group",
    "age_gr" = "Age group",
    "prev_cancer" = "Previous cancer", 
    "treatment" = "Treatment",
    "tissue_RNA_avail" = "Avai. tissue RNA",
    "tissue_DNA_avail" = "Avai. tissue DNA",
    "WGD" = "WGD",
    "HRD" = "HR status",
    "MECOM" = "MECOM",
    "MYC" = "MYC",
    "KRAS" = "KRAS",
    "CCNE1" = "CCNE1",
    "MAGEB10" = "MAGEB10",
    "RB1" = "RB1"
  )[colnames(meta)]
  
  meta$`CNA status` = rep('amp', 118)
  meta$`matched tissue` = rep('yes', 118)
  
  col_list <- list(
    "ctDNA level" = c("low" = "#cc3333", "med" = "#339933", "high" = "#003399"),
    "BMI group" = c("normal" = "#9fc5e8", "overweight" = "#6987a2", "obese" = "#2d4861"),
    "Age group" = c("> 70" = "#8e7cc3", "<= 70" = "#6ec44a"),
    "Treatment" = c("NACT" = "#744700", "PDS" = "#76a5af", "Other" = "gray"),
    "Avai. tissue RNA" = c("yes" = "#8da0cb", "no" = "gray"),
    "Avai. tissue DNA" = c("yes" = "#8da0cb", "no" = "gray"),
    "WGD" = c("yes" = "#8da0cb", "no" = "gray"),
    "Previous cancer" = c("yes" = "#8da0cb", "no" = "#e78ac3"), 
    #"Chronic disease" = c("yes" = "#8da0cb", "no" = "#e78ac3"),
    "Diagnosis stage" = c("I" = "#efe599", "II" = "#f1c232", "III" = "#ce7e00", "IV" = "#783f04"),
    "ctDNA fraction" = circlize::colorRamp2(c(0, 0.315), c("lightblue", "red")),  # adjust range as needed
    "HR status" = c("HRD" = "#ffa500", "HRP" = "#c39797"),
    "MECOM" = c("amplification" = "#ff0800", "gain" = "#e4717a", "normal" = "#e0eeee"),
    "MYC" = c("amplification" = "#ff0800", "gain" = "#e4717a", "normal" = "#e0eeee"),
    "KRAS" = c("amplification" = "#ff0800", "gain" = "#e4717a", "normal" = "#e0eeee", "loss"="#318ce7"),
    "CCNE1" = c("amplification" = "#ff0800", "gain" = "#e4717a", "normal" = "#e0eeee", "loss"="#318ce7"),
    "MAGEB10" = c("deletion"="#00008b", "loss"="#318ce7","normal"= "#e0eeee","gain" = "#e4717a"),
    "RB1" = c("deletion"="#00008b", "loss"="#318ce7","normal"= "#e0eeee")
  )
  
  ha <- ComplexHeatmap::HeatmapAnnotation(
    df = meta,
    which = "column",
    col = col_list,
    annotation_name_side = "left",
    annotation_legend_param = list(nrow = 2,
                                   title_gp = gpar(fontsize = 15, fontfamily = "Arial"),  # legend title
                                   labels_gp = gpar(fontsize = 13, fontfamily = "Arial")) ,
    annotation_name_gp = gpar(fontsize = 15, fontfamily = "Arial")
  )
}
# Draw the heatmap (cohort annotation only)
ht <- ComplexHeatmap::Heatmap(
  dummy_mat, 
  name = "Patient Metadata",
  top_annotation = ha,
  column_names_gp = gpar(fontsize = 14, fontfamily = "Arial"),  # ✅ defined once
  rect_gp = gpar(col = "white", lwd = 1.2),
  #width = unit(ncol(dummy_mat) * 5, "mm"),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_column_names = TRUE,
  column_names_side = "top",
  show_heatmap_legend = FALSE
)

svg('results/figures/main/01_cohort_info.svg', width = 22, height = 7)
ComplexHeatmap::draw(ht,
                     heatmap_legend_side = "bottom",
                     annotation_legend_side = "bottom")
dev.off()
