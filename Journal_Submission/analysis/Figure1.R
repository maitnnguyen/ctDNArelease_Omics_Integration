source('src/00_functions.R')

# color setting:
list_color <- c('#cc3333', '#339933', "#003399")
# load data
# 1. ctdna data
ca125 <- read.delim('data/input/clinicalInfo/ca125.tsv') |>
  dplyr::select(patient, ca125=result)
# Clinical export
{  
  clinical_export <- readxl::read_excel('data/input/clinicalInfo/230929_Clinical export.xlsx', 
                                        sheet = 'DATA') %>%
    as.data.frame() %>%
    #filter(Histology == 'high grade serous') %>%
    dplyr::rename(patient = `Patient card::Patient cohort code_Patient Card`,
                  diag_date = `Date of diagnosis`,
                  relapse_date = `Date of 1st progression`,
                  FIGO_stage = Stage_FIGO2014,
                  age = `Age at Diagnosis`,
                  BMI = `BMI at Dg`,
                  chronic_ill_atDg = `Chronic illnesses at Dg`,
                  BRCA_Decider = `BRCA any in Decider`,
                  treatment = `Treatment strategy`,
                  prim_outcome = `Primary therapy outcome`,
                  prev_cancer = `Previous cancer yes no`,
                  maintain_therapy = `Maintenance therary after 1st line`,
                  oper_date = `Oper 1 date`,
                  ids_date = `Oper 2 IDS date`,
                  dissem_idx1 = `Oper1_Dissemination_Index`,
                  dissem_idx2 = `Oper2_Dissemination_Index`,
                  death_date = `Date of death`,
                  relapse2death = `Time from 1st prog to Death_Days Post progression survival`,
                  progress = `Progression Yes_No`,
                  CRS = `CRS Omental`,
                  outcome_update_date = `Date of Outcome update`,
                  last_chemo_prim_date = `Date of last chemoherapy of primary ther`,
                  eoc_code = `Patient card::Publication code`,
                  OS = OS_KaplanM_allHGSC,
                  PFS = PFS_KaplanM_allHGSC,
                  PFI = PFI_KaplanM_allHGSC,
                  relapse2death = `Time from 1st prog to Death_Days Post progression survival`,
                  DeathCause = `Cause of death`) %>%
    
    dplyr::select(patient, eoc_code, Histology, FIGO_stage, age, BMI, Attention, CRS, 
                  HRD = `HR signature SBS3 per patient`,BRCA_Decider,
                  treatment, prim_outcome, prev_cancer, dissem_idx1, dissem_idx2, 
                  PFI, PFS, OS, progress, alive = Survival, relapse2death, DeathCause) 
  dim(clinical_export)
}  
ctdna_df_prim <- read.delim('data/input/ctdnaData/ctdna_data.tsv') 

dna_concentration <- readxl::read_excel('data/input/ctdnaData/ctDNA_sWGS_QC_combined.xlsx') |>
  as.data.frame() |>
  mutate(sample = sub('-.*', '', `Sample Name`),
         patient = sub('_.*', '', sample)) |>
  dplyr::select(patient, sample, concentration = `Concentration(ng/ul)`, 
                volume = `Volume(ul)`, DNA_amount = `Total amount(ug)`) |>
  dplyr::slice_max(concentration, by = sample)

ctdna_df_prim <- ctdna_df_prim |>
  dplyr::select(-concentration, -amount, -volumn) |>
  dplyr::slice_max(cov, by = sample) |>
  left_join(dna_concentration) |>
  as.data.frame()

length(unique(ctdna_df_prim$patient)) == nrow(ctdna_df_prim)

# threshold:
thresh1 <- max(ctdna_df_prim[ctdna_df_prim$histology=='benign',]$TF)
thresh2 <- 0.065

ctdna_df_prim <- ctdna_df_prim |>
  left_join(
    # ca125 data
    ca125 |>
      dplyr::select(patient, ca125)
  ) |>
  left_join(
    clinical_export) |>
  mutate(age_gr = ifelse(age > 70, '>70', '≤70'),
         ctdna_lev = ifelse(TF < thresh1, 'low',
                        ifelse(TF < thresh2, 'med', 'high')),
         PFI_gr = ifelse(PFI>=365, '≥12M',
                         ifelse(PFI<=180, '≤6M', '6-12M')),
         Stage = ifelse(is.na(FIGO_stage), 'Benign', 
                        sub('[A-C].*', '', FIGO_stage)),
         type = ifelse(is.na(FIGO_stage), 'Benign', 
                       gsub("[^ABC]", "", FIGO_stage)),
         bmi_gr = ifelse(BMI < 25, 'norm', 'obese'),
         outcome=ifelse(prim_outcome %in% c('Complete Response', 'Partial Response'), 
                        prim_outcome, 'No Response')) 

dim(ctdna_df_prim)  
table(ctdna_df_prim$Stage)

# ctdna cohort
hrd_status <- read.delim('data/input/WGSdata/HRD_status.tsv')

ctdna_cohort <- ctdna_df_prim |>
  filter(histology!='benign') |>
  left_join(hrd_status |> 
              dplyr::select(patient, hrd_sig = hrd))

############# 1a - concentration ~ diagnosis stages
{
  fig1a <- ctdna_df_prim %>%
    mutate(figo = ifelse(Stage %in% c('I', 'II'), 'Early', Stage)) |>
    ggplot(aes(x=factor(figo, 
                        labels = c('Benign\n(n=10)',
                                   'Early\n(n=14)',
                                   'III\n(n=85)',
                                   'IV\n(n=20)')),
               y = concentration)) +
    geom_boxplot(aes(fill=figo),width=.5,#outlier.shape = NA,
                 alpha=.7, show.legend = F) +
    #ylim(c(-2.5,2.5))+
    scale_fill_brewer(palette = 'Blues') +
    labs(y='DNA concentration (ng/ml)', x='') +
    stat_compare_means(label.y.npc = .9, size = 5) +
    theme_classic() +
    theme(axis.text.x = element_text(size=13,  family = 'Arial'),
          axis.title = element_text(size=15,  family = 'Arial'),
          axis.text = element_text(size=13, family = 'Arial'))
  fig1a
  svg('results/figures/final/Supp/s2b.svg', width=6, height = 5)
  fig1a
  dev.off()
}

############# 1b - concentration ~ CA125
{
  fig1b <- ctdna_cohort |> 
    mutate(y= TF,
           x = concentration) |> 
    ggscatter(x = 'x', y= 'y', 
              add = 'reg.line', cor.coef = T, 
              cor.coef.size = 5) + 
    labs(y='ctDNA fraction', x = 'cfDNA concentration (ng/ml)') +
    theme(axis.title = element_text(size=15,  family = 'Arial'),
          axis.text = element_text(size=13, family = 'Arial'))
  fig1b
  svg('results/figures/final/Supp/s2d.svg', width=5, height = 5)
  fig1b
  dev.off()
}

############# 1c - concentration ~ treatment
{
  fig1c <- ctdna_cohort |>
    filter(treatment!='Other',
           #Stage %in% c('III', 'IV'))
    ) |>
    ggplot(aes(x=factor(treatment,
                        labels=c('NACT (n=45)', 'PDS (n=73)')),
               y = concentration)) +
    geom_violin(aes(fill=treatment), 
                show.legend = F, alpha = .9)+
    geom_boxplot(width=.1)+
    geom_jitter(size=.3) +
    scale_y_break(c(2,5))+
    labs(y='DNA concentration (ng/ml)', x='', fill='Treatment') +
    stat_compare_means(label.y = 5.3, label.x = 1.3,
                       size = 5)+
    scale_fill_brewer(palette = 'Set2')+
    theme_classic() +
    theme(axis.title = element_text(size=15,  family = 'Arial'),
          axis.text.x = element_blank(),
          axis.text = element_text(size=13, family = 'Arial'),
          legend.position = 'top')
  fig1c
  svg('results/figures/final/Supp/s2a.svg', width=4, height = 3.5)
  fig1c
  dev.off()
  
  fig1c1 <- ctdna_cohort |>
    filter(treatment!='Other',
           #Stage %in% c('III', 'IV'))
    ) |>
    ggplot(aes(x=factor(treatment,
                        labels=c('NACT (n=45)', 'PDS (n=73)')),
               y = TF)) +
    geom_violin(aes(fill=treatment), 
                show.legend = F, alpha = .9)+
    geom_boxplot(width=.1)+
    geom_jitter(size=.3) +
    scale_y_break(c(.1,.25))+
    labs(y='ctDNA fraction', x='', fill='Treatment') +
    stat_compare_means(label.y = .28, label.x = 1.2,
                       size = 5)+
    scale_fill_brewer(palette = 'Set2')+
    theme_classic() +
    theme(axis.title = element_text(size=15,  family = 'Arial'),
          axis.text = element_text(size=13, family = 'Arial'),
          legend.position = 'top')
  fig1c1
  svg('results/figures/final/Supp/s2a1.svg', width=4, height = 3.5)
  fig1c1
  dev.off()
}

############# 1d - TF ~ Stage
{
  comps <- list(c('Benign\n(n=10)','Early\n(n=14)'),
                c('Benign\n(n=10)','III\n(n=85)'),
                c('Benign\n(n=10)','IV\n(n=20)'),
                c('Early\n(n=14)','III\n(n=85)'),
                c('III\n(n=85)','IV\n(n=20)'))
  fig1d <- ctdna_df_prim %>%
    mutate(figo = ifelse(Stage %in% c('I', 'II'), 'Early', Stage)) |>
    ggplot(aes(x=factor(figo, 
                        labels = c('Benign\n(n=10)',
                                   'Early\n(n=14)',
                                   'III\n(n=85)',
                                   'IV\n(n=20)')),
               y = TF)) +
    geom_boxplot(aes(fill=figo),width=.5,#outlier.shape = NA,
                 alpha=.7, show.legend = F) + 
    scale_y_break(c(.12, .28)) +
    stat_compare_means(comparisons = comps, 
                       tip.length = 0.01, 
                       size=5,
                       label = 'p.signif',
                       label.y = c(.055, .075, .1, .09, .28)) +
    stat_compare_means(label.y = .29, size = 5, label.x = 1.1) +
    scale_fill_brewer(palette = 'Blues') +
    labs(y='ctDNA Fraction', x='') +
    theme_classic() +
    theme(axis.title = element_text(size=15,  family = 'Arial'),
          axis.text.x = element_text(size=13, family = 'Arial'),
          axis.text.y = element_text(size=13, 
                                     family = 'Arial'))
  fig1d
  svg('results/figures/final/Supp/s2c.svg', width=6, height = 5)
  fig1d
  dev.off()
}

############# 1ef - concentration/TF ~ prim response
ctdna_cohort <- ctdna_cohort |>
  mutate(prim_outcome = ifelse(prim_outcome %in% c('died during chemotherapy','ND', 'No Chemotherapy'),
                               'Other', 
                          ifelse(prim_outcome %in% c('Stable Disease', 'Progressive Disease'),
                                 'Stable/Progressive', prim_outcome)),
         prev_cancer = ifelse(prev_cancer=='no', 'No', prev_cancer))

{
  fig1e <- ctdna_cohort |>
    filter(prim_outcome != 'Other') |>
    ggplot(aes(x=factor(prim_outcome,
                        levels=c('Complete Response', 'Partial Response',
                                 'Stable/Progressive'),
                        labels=c('Complete\nResponse\n(n=71)', 'Partial\nResponse\n(n=31)',
                                 'Progressive/\nStable\n(n=11)')),
               y = concentration)) +
    geom_boxplot(width=.5, 
                 alpha=.7, show.legend = F)+
    geom_jitter(size=.7) +
    scale_y_break(c(2,5))+
    labs(y='DNA concentration (ng/ml))', x='') +
    theme_classic() +
    theme(axis.text.x = element_text(size=13,  family = 'Arial'), 
          axis.title = element_text(size=15,  family = 'Arial'),
          axis.text = element_text(size=13, family = 'Arial'))
  fig1e
  
  fig1f <- ctdna_cohort |>
    filter(prim_outcome != 'Other') |>
    ggplot(aes(x=factor(prim_outcome,
                        levels=c('Complete Response', 'Partial Response',
                                 'Stable/Progressive'),
                        labels=c('Complete\nResponse\n(n=71)', 'Partial\nResponse\n(n=31)',
                                 'Progressive/\nStable\n(n=11)')),
               y = TF)) +
    scale_y_break(c(.12,.28))+
    geom_boxplot(width=.5, outlier.shape = NA,
                 alpha=.7, show.legend = F)+
    geom_jitter(size=.7) +
    labs(y='ctDNA Fraction', x='') +
    theme_classic() +
    theme(axis.title = element_text(size=15,  family = 'Arial'),
          axis.text = element_text(size=13, family = 'Arial'))
  fig1f
}

{
  svg('results/figures/final/Supp/fig1a.svg', width = 5, height = 5)
  fig1a
  dev.off()
  svg('results/figures/final/Supp/fig1b.svg', width = 5, height = 5)
  fig1b
  dev.off()
  
  svg('results/figures/final/Supp/fig1c.svg', width = 3.5, height = 5)
  fig1c
  dev.off()
  
  svg('results/figures/final/Supp/fig1d.svg', width = 5, height = 5)
  fig1d
  dev.off()
  
  svg('results/figures/final/Main/fig1e.svg', width = 5, height = 5)
  fig1e
  dev.off()
  
  svg('results/figures/final/Main/fig1f.svg', width = 5, height = 5)
  fig1f
  dev.off()
}
  
################ done ##################