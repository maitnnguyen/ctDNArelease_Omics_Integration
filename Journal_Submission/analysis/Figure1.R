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

##### 1a - concentration/ctDNA fraction ~ treatment
{
  fig1a <- ggarrange(
    ctdna_cohort |>
      ggplot(aes(x=factor(treatment,
                          labels=c('NACT', 'PDS')),
                 y = concentration)) +
      geom_violin(aes(fill=treatment), 
                  show.legend = F, alpha = .9)+
      geom_boxplot(width=.1)+
      geom_jitter(size=.3) +
      scale_y_break(c(2,5))+
      labs(y='DNA concentration (ng/ml)', 
           x='', fill='Treatment') +
      stat_compare_means(label.y = 5, label.x = 1.2,
                         size = 5)+
      scale_fill_brewer(palette = 'Set2')+
      theme_classic() +
      theme(axis.title = element_text(size=15,  family = 'Arial'),
            #axis.text.x = element_blank(),
            axis.text = element_text(size=13, family = 'Arial'),
            legend.position = 'top'), 
  
    ctdna_cohort |>
      ggplot(aes(x=factor(treatment,
                          labels=c('NACT', 'PDS')),
                 y = ctdna_fraction)) +
      geom_violin(aes(fill=treatment), 
                  show.legend = F, alpha = .9)+
      geom_boxplot(width=.1)+
      geom_jitter(size=.3) +
      scale_y_break(c(2,5))+
      labs(y='ctDNA fraction', 
           x='', fill='Treatment') +
      stat_compare_means(label.y = 5, label.x = 1.2,
                         size = 5)+
      scale_fill_brewer(palette = 'Set2')+
      theme_classic() +
      theme(axis.title = element_text(size=15,  family = 'Arial'),
            #axis.text.x = element_blank(),
            axis.text = element_text(size=13, family = 'Arial'),
            legend.position = 'top'),
    nrow = 2
  )
}

##### 1bc - primary treatment response ~ DNA concentration/ctDNA fraction
dat1 <- ctdna_cohort |>
  mutate(prim_outcome = ifelse(prim_outcome %in% c('died during chemotherapy','ND', 'No Chemotherapy'),
                               'Other', 
                               ifelse(prim_outcome %in% c('Stable Disease', 'Progressive Disease'),
                                      'Stable/Progressive', prim_outcome)))
{
  fig1b <- dat1 |>
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
  
  fig1c <- dat1 |>
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
}

##### 1de - survival analysis
surv_df <- ctdna_cohort %>%
  mutate(event = ifelse(progress=='Yes', 1, 0),
         died = ifelse(alive == 'Alive', 0, 1),
         PFI = ifelse(PFI < 0, 0, PFI),
         # 0.01418: max of benign samples
         pheno = factor(ctdna_lev, levels = c('low' , 'med', 'high')),
         # BMI group
         bmi_gr = ifelse(BMI < 25, 'normal', 
                         ifelse(BMI < 30, 'overweight', 'obese'))) |>
  filter(!is.na(PFI))

{
  pfi_object <- Surv(time = surv_df$PFI/30, event = surv_df$event)
  fit1 <- survfit(pfi_object ~ pheno, data = surv_df)
  #pairwise_survdiff(Surv(time = PFI, event = event) ~ pheno, data = surv_df)
  
  # Figure 1d: Kaplan Meier for ctdna_group
  fig1d <- ggsurvplot(fit1, data = surv_df, pval = TRUE, 
                      # Change legends: title & labels
                      legend.title = "",
                      conf.int = F, 
                      legend.labs = c("low ctDNA", "med ctDNA","high ctDNA"),
                      risk.table = T,
                      legend = 'none',
                      ylab='PFI Probability',
                      xlab='Time (Months)',
                      risk.table.col = "strata",
                      palette = list_color,
                      show.legend=F)
  
  # Fit a Cox proportional hazards model
  surv_df$ctdna_level = relevel(surv_df$pheno, ref = "med")
  coxph1 <- coxph(Surv(PFI, event) ~ ctdna_level + 
                    treatment + 
                    prev_cancer + 
                    `BRCA1/2` +
                    bmi_gr +
                    age_gr , 
                  data = surv_df )
  
  fig1e <- ggforest(coxph1, data = surv_df , fontsize = 1)
  }

################ done ##################
