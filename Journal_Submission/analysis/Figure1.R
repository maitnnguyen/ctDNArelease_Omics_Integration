########################################## Figures ##########################################
##### 4ab - Kaplan Meier (KM) curves
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
##### 4cd - Univariate + Multivariate Cox regression
{
  # PDS
  surv_df_pds$ctdna_level = relevel(surv_df_pds$pheno, ref = "med")
  
  dependent_pfi <- "Surv(PFI, event)"
  expl_var <- c("ctdna_level", "residual", "stage", "HRD", "bmi_gr", "age_gr")

  # regression results below were used for generating figure 4c
  cox_res_pds <- surv_df_pds %>% 
    finalfit::finalfit(dependent_pfi, 
                       expl_var, 
                       add_dependent_label = FALSE)

  # multivariate plot
  coxph2 <- coxph(Surv(PFI, event) ~ ctdna_level + 
                  operable +
                  stage +
                  HRD +
                  bmi_gr +
                  age_gr,
                  data = surv_df_nact )
  
  # combining univariate + multivariate plot
  fig4c <- ggplot(cox_res_pds, aes(x = HR, xmax = high_conf, xmin = low_conf, 
                   y = label, color = model, fill = model)) +
          geom_point(size = 3, shape = 18, position = position_dodge(width = 0.5)) +
          geom_linerange(position = position_dodge(width = 0.5), size = 1) +
          geom_vline(xintercept = 1, size = 1) +
          #geom_hline(yintercept = 0, size = 1) +
          facet_grid(variable ~ ., scales = "free_y", space = "free_y") +
          scale_alpha_identity() +
          scale_fill_manual(values = barCOLS) +
          scale_color_manual(values = dotCOLS) +
          scale_y_discrete(name = "") +
          scale_x_continuous(name = "Hazard ratio") +
          geom_text(aes(label = variable), x = -Inf, y = Inf, hjust = 1, 
                    vjust = .5,size = 6,
                    check_overlap = TRUE, color = "darkgreen") +
          coord_trans(x = 'log10', clip = 'off') +
          theme(
            panel.background = element_blank(),
            panel.spacing = unit(0, "pt"),
            axis.line.x.bottom = element_line(size = 1),
            axis.text.y.left =  element_text(margin = margin(l = 15, unit = "pt"),
                                             family = 'Arial', size=13, color = 'black',
                                             vjust = 0),
            strip.background = element_blank(), 
            strip.text = element_blank(),
            legend.position = 'none',
            axis.title.x = element_text(size = 15),
            axis.text.x = element_text(family = 'Arial', color = 'black',size = 13),
            axis.text.y = element_text(family = 'Arial', color = 'black',size = 13)
          )
  
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

  # settings
  dependent_pfi <- "Surv(PFI, event)"
  expl_var <- c("ctdna_level", "operable", "stage", "HRD", "bmi_gr", "age_gr")

  # regression results below were used for generating figure 4d
  cox_res_nact <- dat %>% 
    finalfit::finalfit(dependent_pfi, 
                       expl_var, 
                       add_dependent_label = FALSE)

  # multivariate plot
  coxph2 <- coxph(Surv(PFI, event) ~ ctdna_level + 
                  operable +
                  stage +
                  HRD +
                  bmi_gr +
                  age_gr,
                  data = surv_df_nact )

  # combining univariate + multivariate plot
  fig4d <- ggplot(cox_res_nact, aes(x = HR, xmax = high_conf, xmin = low_conf, 
                   y = label, color = model, fill = model)) +
          geom_point(size = 3, shape = 18, position = position_dodge(width = 0.5)) +
          geom_linerange(position = position_dodge(width = 0.5), size = 1) +
          geom_vline(xintercept = 1, size = 1) +
          #geom_hline(yintercept = 0, size = 1) +
          facet_grid(variable ~ ., scales = "free_y", space = "free_y") +
          scale_alpha_identity() +
          scale_fill_manual(values = barCOLS) +
          scale_color_manual(values = dotCOLS) +
          scale_y_discrete(name = "") +
          scale_x_continuous(name = "Hazard ratio") +
          geom_text(aes(label = variable), x = -Inf, y = Inf, hjust = 1, 
                    vjust = .5,size = 6,
                    check_overlap = TRUE, color = "darkgreen") +
          coord_trans(x = 'log10', clip = 'off') +
          theme(
            panel.background = element_blank(),
            panel.spacing = unit(0, "pt"),
            axis.line.x.bottom = element_line(size = 1),
            axis.text.y.left =  element_text(margin = margin(l = 15, unit = "pt"),
                                             family = 'Arial', size=13, color = 'black',
                                             vjust = 0),
            strip.background = element_blank(), 
            strip.text = element_blank(),
            legend.position = 'none',
            axis.title.x = element_text(size = 15),
            axis.text.x = element_text(family = 'Arial', color = 'black',size = 13),
            axis.text.y = element_text(family = 'Arial', color = 'black',size = 13)
          )
}

################ done ##################
