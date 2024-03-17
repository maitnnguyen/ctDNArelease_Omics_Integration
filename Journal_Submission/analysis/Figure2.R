source('src/00_functions.R')

# survival analysis
surv_df <- ctdna_cohort %>%
  mutate(event = ifelse(progress=='Yes', 1, 0),
         died = ifelse(alive == 'Alive', 0, 1),
         PFI = ifelse(PFI < 0, 0, PFI),
         HRD = sub('NA', 'ND', ifelse(is.na(HRD), 'NA', HRD)),
         # 0.01418: max of benign samples
         pheno = factor(ctdna_lev, levels = c('low' , 'med', 'high')),
         treatment = ifelse(treatment == 'Other', 'NACT', treatment)) |>
  filter(!is.na(PFI)) 
dim(surv_df)
pfi_object <- Surv(time = surv_df$PFI/30, event = surv_df$event)
pfs_object <- Surv(time = surv_df$PFS/30, event = surv_df$died)
os_object <- Surv(time = surv_df$OS/30, event = surv_df$died)

fit1 <- survfit(pfi_object ~ pheno, data = surv_df)
#pairwise_survdiff(Surv(time = PFI, event = event) ~ pheno, data = surv_df)

# Figure 2a: Kaplan Meier for ctdna_group
fig2a <- ggsurvplot(fit1, data = surv_df, pval = TRUE, 
                    #surv.median.line = "hv",
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

svg('results/figures/final/Main/fig2a.svg', width = 6, height = 7)
fig2a
dev.off()

# Fit a Cox proportional hazards model
dat <- surv_df |>
  mutate(prev_cancer = tolower(prev_cancer),
         `BRCA1/2` = tolower(ifelse(BRCA_Decider=='NA', 'no',BRCA_Decider)),
         figo = ifelse(Stage %in% c('I', 'II'), 'early', Stage),
         bmi_gr = ifelse(BMI < 25, 'normal', 
                      ifelse(BMI < 30, 'overweight', 'obese')))
dim(dat)
dat$ctdna_level = relevel(dat$pheno, ref = "med")
coxph1 <- coxph(Surv(PFI, event) ~ ctdna_level + 
                  treatment + 
                  prev_cancer + 
                  `BRCA1/2` +
                  bmi_gr +
                  #figo +
                  age_gr , 
                data = dat )
fig2b <- ggforest(coxph1, data = dat , fontsize = 1)
fig2b

ggsurvplot(survfit(Surv(OS, died) ~ pheno, data = dat), data = dat, pval = TRUE, 
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

svg('results/figures/final/Main/fig2b.svg', width = 6, height = 7)
fig2b
dev.off()

########## done ##########
coxph2 <- coxph(Surv(PFI, event) ~ ctdna_level + 
                  #hrd_sig + 
                  treatment + 
                  #prev_cancer + 
                  `BRCA1/2` , 
                data = dat)
ggforest(coxph2, data = dat , fontsize = 1)

############## NACT Patients ###############
dat1 <- dat %>% filter(treatment != 'PDS')
dim(dat1)
coxph3 <- coxph(Surv(PFI, event) ~ ctdna_level + 
                  hrd_sig + 
                  #treatment + 
                  bmi_gr +
                  Stage +
                  age_gr , 
                data = dat1 )
ggforest(coxph3, data = dat1 , fontsize = 1)

ggsurvplot(survfit(Surv(PFI, event) ~ pheno, data = dat1), data = dat1, pval = TRUE, 
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

################ done ##################