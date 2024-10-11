source('/int/00_functions.R')

# color setting:
list_color <- c('#cc3333', '#339933', "#003399")
# load data
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

# 2. Figures
##### 1ab - Kaplan Meier (KM) curves
surv_df <- ctdna_cohort %>%
  mutate(event = ifelse(progress=='Yes', 1, 0),
         died = ifelse(alive == 'Alive', 0, 1),
         PFI = ifelse(PFI < 0, 0, PFI),
         # 0.0142: max of benign samples
         ctdna_lev = factor(ctdna_lev, levels = c('low' , 'med', 'high')),
         # BMI group
         bmi_gr = ifelse(BMI < 25, 'normal', 
                         ifelse(BMI < 30, 'overweight', 'obese')),
         stage = ifelse(Stage %in% c('I', 'II'), 'early', Stage),
         age_gr = ifelse(age > 70, '>70', '≤70') |>
  filter(!is.na(PFI))

### perform on subset of NACT & PDS patients
surv_df_nact <- surv_df |>
  filter(treatment=='NACT') |>
  mutate(operable = ifelse(residual_ids == 'no_surgery','No','Yes'))

surv_df_pds <- surv_df |>
  filter(treatment=='PDS') |>
  mutate(residual = ifelse(resi_pds == 'res_0', 'good', 'bad'))

{
  # PDS
  pfi_object1 <- Surv(time = surv_df_pds$PFI/30, event = surv_df_pds$event)
  fit1 <- survfit(pfi_object1 ~ ctdna_lev, data = surv_df_pds)
  
  pfi_object2 <- Surv(time = surv_df_nact$PFI/30, event = surv_df_nact$event)
  fit2 <- survfit(pfi_object2 ~ ctdna_lev, data = surv_df_nact)
  
  #check pairwise in KM with BH correction for p values
  #pairwise_survdiff(Surv(time = PFI, event = event) ~ ctdna_lev, data = surv_df_nact)
  
  # Figure 1a: KM for PDS patients
  fig1a <- ggsurvplot(fit1, data = surv_df_pds, pval = TRUE, 
                    title = 'ctDNA levels in PDS Patients',
                    # Change legends: title & labels
                    legend.title = "",
                    conf.int = F, 
                    legend.labs = c("low ctDNA", "med ctDNA","high ctDNA"),
                    risk.table = T,
                    legend = 'none',#c(.8,.8),
                    ylab='PFI Probability',
                    xlab='Time (Months)',
                    risk.table.col = "strata",
                    palette = list_color,
                    #risk.table.y.text=F,
                    show.legend=F)
  
  # Figure 1b: KM for NACT patients
  fig1b <- ggsurvplot(fit2, data = surv_df_nact, pval = TRUE, 
                    title = 'ctDNA levels in NACT Patients',
                    # Change legends: title & labels
                    legend.title = "",
                    conf.int = F, 
                    legend.labs = c("low ctDNA", "med ctDNA","high ctDNA"),
                    risk.table = T,
                    legend = 'none',#c(.8,.8),
                    ylab='PFI Probability',
                    xlab='Time (Months)',
                    risk.table.col = "strata",
                    palette = list_color,
                    #risk.table.y.text=F,
                    show.legend=F)
}
##### 1cd - Cox regression
{
  # PDS
  surv_df_pds$ctdna_level = relevel(surv_df_pds$pheno, ref = "med")
  coxph1 <- coxph(Surv(PFI, event) ~ ctdna_level + 
                  residual + 
                  stage +
                  HRD +
                  bmi_gr +
                  age_gr , 
                  data = surv_df_pds )
  
  fig1c <- ggforest(coxph1, data = surv_df_pds , fontsize = 1)

  
  # NACT
  surv_df_nact$ctdna_level = relevel(surv_df$ctdna_level, ref = "med")
  coxph2 <- coxph(Surv(PFI, event) ~ ctdna_level + 
                  operable +
                  stage +
                  HRDconcensus +
                  bmi_gr +
                  age_gr,
                  data = surv_df_nact )
  
  fig1d <- ggforest(coxph2, data = surv_df_nact , fontsize = 1)
}

################ done ##################
