source('src/00_functions.R')

genome_size <- 2900077904/1e6 # (Megabase)
# load data
load('~/mnt/storageBig8/work/nguyenma/projects/RNA/OmicsIntegration/data/processed/mut_call/mutation_call_4_6.RData')

# mutational signature on sample level
load('data/input/WGSdata/mutational_signature_sample_lev_231005.RData')

load('data/input/WGSdata/ascat_results.RData')

dna_cohort <- ascat_est |>
  filter(patient %in% ctdna_cohort$patient,
         !(tissue %in% c('Asc', 'Plf'))#, purity >= 0.1
         ) |>
  filter(!is.na(purity), aberrant == TRUE) |>
  mutate(tissue_gr = ifelse(tissue %in% c('Ome', 'Per'), tissue, 
                            ifelse(tissue %in% c('Ova', 'Tub', 'Adn', 'And'),
                                   'Ova_Tub','Other'))) |>
  dplyr::slice_max(purity, by = c(patient, tissue_gr)) |>
  dplyr::slice_max(goodnessOfFit, by = c(patient, tissue_gr)) |>
  left_join(ctdna_cohort |>
              dplyr::select(patient, TF, Stage,
                            ctdna_lev, treatment)) 
  
length(unique(dna_cohort$patient))
summary(dna_cohort$purity)
dim(dna_cohort)
# - mutation call
{
  mut_call <- snv_res %>% 
    mutate(tissue = stringr::str_match(sample, "^[^_]+_[piro]\\d?(LN|[A-Z][a-z]+)")[, 2],
           tissue_gr = ifelse(tissue %in% c('Ome', 'Per'), tissue, 
                              ifelse(tissue %in% c('Ova', 'Tub', 'Adn', 'And'),
                                     'Ova_Tub','Other')),
           patient = sub('_.*', '', sample)) |>
    filter(sample %in% dna_cohort$sample) |>
    inner_join(ctdna_cohort |>
                 dplyr::select(-sample)) 
  #length(unique(mut_call$patient))==length(unique(dna_cohort$patient))
  
  mut_genes <- mut_call |>
    dplyr::select(chr, pos, ref, alt, Func, gene,VAF, 
                  TF, patient, ctdna_lev, hrd_sig) |>
    dplyr::slice_max(VAF, by=c(chr, pos, ref, alt, Func, gene, 
                               TF, patient, ctdna_lev, hrd_sig))
  #length(unique(mut_genes$patient))==length(unique(dna_cohort$patient))
  mut_pat_lev <- mut_call |>
    dplyr::select(chr, pos, ref, alt, Func, gene, TF, patient, ctdna_lev, hrd_sig) |>
    distinct()
  dim(mut_pat_lev)
  #length(unique(mut_pat_lev$patient))==length(unique(dna_cohort$patient))
}

# - mutational signature data - sample level
{
  mut_sbs_sample_prop <- single_sample_sbs_prop %>%
    dplyr::rename(sample = sample_id) |>
    filter(sample %in% dna_cohort$sample) |>
    left_join(ctdna_cohort %>%
                dplyr::select(patient, TF, eoc_code, 
                              ctdna_lev)) |>
    mutate(tissue = stringr::str_match(sample, "^[^_]+_[piro]\\d?(LN|[A-Z][a-z]+)")[, 2],
           tissue_gr = ifelse(tissue %in% c('Ome', 'Per'), tissue, 
                              ifelse(tissue %in% c('Ova', 'Tub', 'Adn', 'And'),
                                     'Ova_Tub','Other'))) 
  
  mut_sbs_sample_count <- single_sample_sbs_count %>%
    dplyr::rename(sample = sample_id) |>
    filter(sample %in% dna_cohort$sample) |>
    left_join(ctdna_cohort %>%
                dplyr::select(patient, TF, eoc_code, 
                              ctdna_lev)) |>
    mutate(tissue = stringr::str_match(sample, "^[^_]+_[piro]\\d?(LN|[A-Z][a-z]+)")[, 2],
           tissue_gr = ifelse(tissue %in% c('Ome', 'Per'), tissue, 
                              ifelse(tissue %in% c('Ova', 'Tub', 'Adn', 'And'),
                                     'Ova_Tub','Other'))) 
  
  mut_id_sample_prop <- single_sample_id_prop %>%
    dplyr::rename(sample = sample_id) |>
    filter(sample %in% dna_cohort$sample) |>
    left_join(ctdna_cohort %>%
                dplyr::select(patient, TF, eoc_code, 
                              ctdna_lev)) |>
    mutate(tissue = stringr::str_match(sample, "^[^_]+_[piro]\\d?(LN|[A-Z][a-z]+)")[, 2],
           tissue_gr = ifelse(tissue %in% c('Ome', 'Per'), tissue, 
                              ifelse(tissue %in% c('Ova', 'Tub', 'Adn', 'And'),
                                     'Ova_Tub','Other'))) 
  
  mut_id_sample_count <- single_sample_id_count %>%
    dplyr::rename(sample = sample_id) |>
    filter(sample %in% dna_cohort$sample) |>
    left_join(ctdna_cohort %>%
                dplyr::select(patient, TF, eoc_code, 
                              ctdna_lev)) |>
    mutate(tissue = stringr::str_match(sample, "^[^_]+_[piro]\\d?(LN|[A-Z][a-z]+)")[, 2],
           tissue_gr = ifelse(tissue %in% c('Ome', 'Per'), tissue, 
                              ifelse(tissue %in% c('Ova', 'Tub', 'Adn', 'And'),
                                     'Ova_Tub','Other'))) 
  
  mut_dbs_sample_prop <- single_sample_dbs_prop %>%
    dplyr::rename(sample = sample_id) |>
    filter(sample %in% dna_cohort$sample) |>
    left_join(ctdna_cohort %>%
                dplyr::select(patient, TF, eoc_code, 
                              ctdna_lev)) |>
    mutate(tissue = stringr::str_match(sample, "^[^_]+_[piro]\\d?(LN|[A-Z][a-z]+)")[, 2],
           tissue_gr = ifelse(tissue %in% c('Ome', 'Per'), tissue, 
                              ifelse(tissue %in% c('Ova', 'Tub', 'Adn', 'And'),
                                     'Ova_Tub','Other'))) 
  
  mut_dbs_sample_count <- single_sample_dbs_count %>%
    dplyr::rename(sample = sample_id) |>
    filter(sample %in% dna_cohort$sample) |>
    left_join(ctdna_cohort %>%
                dplyr::select(patient, TF, eoc_code, 
                              ctdna_lev)) |>
    mutate(tissue = stringr::str_match(sample, "^[^_]+_[piro]\\d?(LN|[A-Z][a-z]+)")[, 2],
           tissue_gr = ifelse(tissue %in% c('Ome', 'Per'), tissue, 
                              ifelse(tissue %in% c('Ova', 'Tub', 'Adn', 'And'),
                                     'Ova_Tub','Other'))) 
}
genome_size <- 2900
########### Figure supplementary - tissue purity and ctDNA TF

{
  s3 <- dna_cohort |>
    mutate(tissue_gr = factor(tissue_gr, levels=c('Ova_Tub', 'Ome',
                                                  'Per', 'Other'))) |>
    ggscatter(y = 'TF', x = 'purity', add = 'reg.line', 
              cor.coef = T, facet.by = 'tissue_gr') +
    labs(x='Tissue Purity', y = 'ctDNA Fraction') +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size=15), 
          axis.text.y = element_text(size=13),
          #axis.text.x = element_blank(),
          legend.position = c(.85,.75),
          legend.text = element_text(size=10),
          legend.title = element_text(size=15),
          strip.text.x = element_text(size = 15))
  s3
}
svg('results/figures/final/Supp/s3a.svg', width = 5, height = 5)
s3
dev.off()

# fig 3a: copy-number breakpoint
{
  ascat_break <- ascat_break |>
    mutate(patient = sub('_.*', '', sample),
           tissue = stringr::str_match(sample, '^[^_]+_[piro]\\d?(LN|[A-Z][a-z]+)')[, 2],
           tissue_gr = ifelse(tissue %in% c('Ome', 'Per'), tissue, 
                              ifelse(tissue %in% c('Ova', 'Tub', 'Adn', 'And'),
                                     'Ova_Tub','Other'))) |>
    filter(sample %in% dna_cohort$sample) |>
    left_join(ctdna_cohort |>
                dplyr::select(patient, ctdna_lev))
    
  fig3a <- ascat_break |> 
    left_join(dna_cohort |>
                dplyr::select(sample, purity, TP53.VAF)) |>
    filter(sample %in% highest_purity_dna$sample) |> 
    #dplyr::slice_max(purity, by = patient) |> 
    #dplyr::slice_max(TP53.VAF, by = patient) |>
    as.data.frame() |> 
    ggplot(aes(x=factor(ctdna_lev, 
                        levels = c('low', 'med', 'high')), 
               y = breaks/23)) +
    geom_boxplot(aes(color=factor(ctdna_lev, 
                              levels = c('low', 'med', 'high'))),
                 show.legend = F) +
    geom_jitter(aes(color=factor(ctdna_lev, 
                                  levels = c('low', 'med', 'high'))),
                 show.legend = F) +
    scale_color_manual(values = list_color) +
    stat_compare_means(size=5, label.y.npc = .85) +
    stat_compare_means(comparisons = list(c('low', 'med'),
                                          c('high', 'med'),
                                          c('high', 'low')),
                       label = 'p.signif', size=5,
                       label.y = c(32,34,36), 
                       tip.length = .01) +
    labs(fill='ctDNA level', x='', y = 'Avg # Break Points per Chromosome') +
    theme_pubr() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size=12), 
          axis.text.y = element_text(size=10),
          axis.text.x = element_blank(),
          legend.position = 'none',
          legend.text = element_text(size=10),
          legend.title = element_text(size=15),
          strip.text.x = element_text(size = 15))

  fig3a
}
svg('results/figures/final/Main/fig3a.svg', width = 5, height = 5)
fig3a
dev.off()
## Figure 3b: mutational burden (N/Mb) by ctdna groups
{
  fig3b <- mut_pat_lev |>
    dplyr::group_by(patient, ctdna_lev) %>%
    dplyr::summarise(n_mut = n()/genome_size) %>%
    as.data.frame() %>%
    ggplot(aes(x = factor(ctdna_lev, 
                          levels = c('low', 'med', 'high')), 
               y = n_mut)) +
    geom_violin(aes(color=factor(ctdna_lev, 
                                 levels = c('low', 'med', 'high')))) +
    geom_jitter(
      aes(color = factor(ctdna_lev, 
                         levels = c('low', 'med', 'high'))), 
      position = position_jitterdodge(jitter.width = 1.4, dodge.width = .4),
      size = 2, show.legend = F
    ) +
    stat_summary(color = 'black', fun.data="mean_sdl",  
                 fun.args = list(mult=1), 
                 geom = "pointrange",  size = 1,
                 position = position_dodge(1))+
    scale_color_manual(values = list_color) +
    stat_compare_means(size=5) +
    stat_compare_means(comparisons = list(c('high', 'med'),
                                          c('low', 'med'),
                                          c('high', 'low')),
                       label = 'p.signif', size=5,
                       label.y = c(13,11,15), 
                       tip.length = .01) +
    labs(x = '', y='N somatic mutations per Megabase', fill='ctDNA Level') +
    #ylim(c(0,.4)) +
    theme_pubr() +
    theme(axis.title = element_text(size=15), 
          axis.text.y = element_text(size=13),
          axis.text.x = element_blank(),
          legend.position = 'none',
          legend.text = element_text(size=10),
          legend.title = element_text(size=12),
          strip.text.x = element_text(size = 15))
  fig3b
}

svg('results/figures/final/Main/fig3b.svg', width = 5, height = 5)
fig3b
dev.off()

# mutation information integration
### part c,d,e
# get highest purity samples
highest_purity_dna <- dna_cohort |>
  dplyr::slice_max(purity, by = c(patient)) |>
  dplyr::slice_max(TP53.VAF, by = c(patient))
dim(highest_purity_dna)

# figure 3d: Mutation type ~ ctdna_group
{
  mut_sigs <- mut_sbs_sample_count |>
    filter(sample %in% highest_purity_dna$sample) |>
    dplyr::select(patient, mut_count, TF, ctdna_lev) |>
    #dplyr::group_by(patient, TF, ctdna_lev) |>
    #dplyr::summarise(N = median(mut_count)) |>
    mutate(type='SBS') |>
    rbind(
      mut_dbs_sample_count |>
        filter(sample %in% highest_purity_dna$sample) |>
        dplyr::select(patient, mut_count, TF, ctdna_lev) |>
        #dplyr::group_by(patient, TF, ctdna_lev) |>
        #dplyr::summarise(N = median(mut_count)) |>
        mutate(type='DBS')
    ) |>
    rbind(
      mut_id_sample_count |>
        filter(sample %in% highest_purity_dna$sample) |>
        dplyr::select(patient, mut_count, TF, ctdna_lev) |>
        #dplyr::group_by(patient, TF, ctdna_lev) |>
        #dplyr::summarise(N = median(mut_count)) |>
        mutate(type='ID')
    ) 
  
  fig3c <- mut_sigs |>
    ggplot(aes(x = factor(ctdna_lev, levels = c('low', 'med', 'high')),
               y = mut_count/genome_size)) + #geom_boxplot()
    geom_jitter(
      aes(color = factor(ctdna_lev, 
                         levels = c('low', 'med', 'high'))), 
      position = position_jitterdodge(jitter.width = 1, dodge.width = .4),
      size = 2, show.legend = F
    ) +
    stat_summary(color = 'black', fun.data="mean_sdl",  
                 fun.args = list(mult=1), 
                 geom = "pointrange",  size = 1,
                 position = position_dodge(1))+
    scale_color_manual(values = list_color) +
    facet_wrap(~factor(type, levels=c('SBS', 'DBS', 'ID')), scales = 'free_y') +
    stat_compare_means(label.sep = '\n',
                       size=5, label.y.npc = .9) +
    scale_color_manual(values = list_color) +
    labs(y='Mutation Count per Megabase', x ='') +
    theme_pubr() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size=15), 
          axis.text.x = element_blank(),
          #axis.text.y = element_blank(),
          legend.position = c(.85,.75),
          legend.text = element_text(size=10),
          legend.title = element_text(size=15),
          strip.text.x = element_text(size = 15))
  fig3c
}
svg('results/figures/final/Main/fig3c.svg', width = 10, height = 4)
fig3c
dev.off()


sig_df1 <- mut_sbs_sample_count |>
  filter(sample %in% highest_purity_dna$sample) |>
  tidyr::gather(key = 'signature', value = 'count', SBS1:SBS40) |>
  as.data.frame() |>
  rbind(
    mut_dbs_sample_count |>
      filter(sample %in% highest_purity_dna$sample) |>
      tidyr::gather(key = 'signature', value = 'count', DBS1:DBS11) |>
      as.data.frame()
  )


sig_df2 <- mut_id_sample_count |>
  filter(sample %in% highest_purity_dna$sample) |>
  tidyr::gather(key = 'signature', value = 'count', ID1:ID10) |>
  as.data.frame()

# figure 2d: significant signatures
{
  fig3d1 <- sig_df1 |>
      filter(signature %in% c('SBS40', 'DBS9')) |>
      ggplot(aes(x = factor(ctdna_lev, 
                            levels = c('low', 'med', 'high')),
                 y = count/genome_size)) +
    geom_jitter(
      aes(color = factor(ctdna_lev, 
                         levels = c('low', 'med', 'high'))), 
      position = position_jitterdodge(jitter.width = 1, dodge.width = .4),
      size = 1, show.legend = F
    ) +
    stat_summary(color = '#696969', fun.data="mean_sdl",  
                 fun.args = list(mult=1), 
                 geom = "pointrange",  size = .5,
                 position = position_dodge(1))+
    scale_color_manual(values = list_color) +
    stat_compare_means(size = 4, 
                         label.sep = '\n',
                         label.y.npc = .9) +
      facet_wrap(~factor(signature, levels=c('SBS40', 'DBS9')), 
                 scale = 'free', ncol = 3) +
      labs(x='', y ='Mutation Count per Mb') +
      theme_pubr() +
      theme(panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_text(size=12), 
            axis.text.y = element_text(size=10),
            axis.text.x = element_blank(),
            legend.position = c(.85,.75),
            legend.text = element_text(size=10),
            legend.title = element_text(size=15),
            strip.text.x = element_text(size = 15)) 
  fig3d1
  
  fig3d2 <- sig_df2 |>
      filter(signature %in% c('ID6', 'ID8')) |>
      ggplot(aes(x = factor(ctdna_lev, 
                            levels = c('low', 'med', 'high')),
                 y = count/genome_size)) +
    geom_jitter(
      aes(color = factor(ctdna_lev, 
                         levels = c('low', 'med', 'high'))), 
      position = position_jitterdodge(jitter.width = 1, dodge.width = .4),
      size = 1, show.legend = F
    ) +
    stat_summary(color = '#696969', fun.data="mean_sdl",  
                 fun.args = list(mult=1), 
                 geom = "pointrange",  size = .5,
                 position = position_dodge(1))+
    stat_compare_means(size = 4, 
                       label.sep = '\n',
                       label.y.npc = .9) +
     scale_color_manual(values = list_color) +
      facet_wrap(~signature, scale = 'free', ncol = 3) +
      labs(x='', y ='Mutation Count per Mb') +
      theme_pubr() +
      theme(panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_text(size=12), 
            axis.text.y = element_text(size=10),
            axis.text.x = element_blank(),
            legend.position = c(.85,.75),
            legend.text = element_text(size=10),
            legend.title = element_text(size=15),
            strip.text.x = element_text(size = 15)) 
  fig3d2
}

svg('results/figures/final/Main/fig3d1.svg',width=6, height=3.5)
fig3d1
dev.off()

svg('results/figures/final/Main/fig3d2.svg',width=6, height=3.5)
fig3d2
dev.off()

### figure 3e: sbs1 
sig_df3 <- mut_sbs_sample_prop |>
  filter(sample %in% highest_purity_dna$sample) |>
  as.data.frame() |>
  tidyr::gather(key = 'signature', value = 'prop', SBS1:SBS40) |>
  as.data.frame() |>
  rbind(
    mut_id_sample_prop |>
      filter(sample %in% highest_purity_dna$sample) |>
      as.data.frame() |>
      tidyr::gather(key = 'signature', value = 'prop', ID1:ID10) |>
      as.data.frame()
  )

fig3e <- sig_df3 |>
  filter(signature %in% c('SBS1')) |>
  ggplot(aes(x = factor(ctdna_lev, 
                        levels = c('low', 'med', 'high')),
             y = prop)) +
  geom_violin(aes(fill = factor(ctdna_lev, 
                                levels = c('low', 'med', 'high'))), 
              show.legend = F, alpha=.9) +
  geom_boxplot(width=.2, outlier.shape = NA,
               show.legend = F) +
  #geom_jitter(size=.2, color="#999999") +
  stat_compare_means(size = 3.5, label.y.npc = .9,
                     label.sep = '\n', label.x = 1.5) +
  scale_fill_manual(values = list_color) +
  labs(x='', y ='SBS1 Contribution') +
  theme_pubr() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size=12), 
        axis.text.y = element_text(size=10),
        axis.text.x = element_blank(),
        legend.position = c(.85,.75),
        legend.text = element_text(size=10),
        legend.title = element_text(size=15),
        strip.text.x = element_text(size = 15))
fig3e

svg('results/figures/final/Main/fig3e.svg',width=3.1, height=3.5)
fig3e
dev.off()

# 3f: ID5 with prop-test
sig_df4 <- mut_id_sample_prop |>
  filter(sample %in% highest_purity_dna$sample) |>
  mutate(id5 = ifelse(ID5>0, 'non-zero', 'zero'))
stat.test <- xchisq.test(ctdna_lev ~id5, data=sig_df4)
fig3f <- sig_df4 |>
  ggplot(aes(x = factor(ctdna_lev, 
                        levels = c('low', 'med', 'high')))) +
  geom_bar(aes(fill = id5,
               color = factor(ctdna_lev, 
                              levels = c('low', 'med', 'high'))), 
           position = 'fill', linewidth=1.2) +
  annotate("text", x=2, y =1.1, size = 4, family='Arial',
           label = paste0('X-squared: ',
                          round((chisq.test(sig_df4$ctdna_lev, sig_df4$id5))$stat,4), 
                         ',\np-value: ', 
                         round((chisq.test(sig_df4$ctdna_lev, 
                                            sig_df4$id5))$p.value,4))) +
  scale_fill_brewer(palette = 'PuBu', direction = -1) +
  scale_color_manual(values=list_color) +
  labs(x='', y ='ID5', color='ctDNA level',
       fill = 'ID5') +
  ylim(0,1.15)+
  theme_pubr() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size=12), 
        axis.text.y = element_text(size=10),
        axis.text.x = element_blank(),
        legend.position = 'right',
        legend.text = element_text(size=10),
        legend.title = element_text(size=15),
        strip.text.x = element_text(size = 15))
fig3f
svg('results/figures/final/Main/fig3f.svg',width=4, height=3.5)
fig3f
dev.off()

# SBS3
dat <- single_sample_sbs_prop |> 
  filter(patient %in% dna_cohort$patient) |>
  dplyr::group_by(patient) |> 
  dplyr::summarise(sbs3 = sum(SBS3 > 0), N = n()) |> 
  as.data.frame() |> 
  mutate(hrd = ifelse(sbs3 == N, 'HRD', 'HRP')) |>
  left_join(ctdna_cohort |>
              dplyr::select(patient, ctdna_lev)) 
length(unique(dat$patient))==nrow(dat)
ggplot(dat %>% dplyr::group_by(hrd) %>% dplyr::summarise(n=n()), 
       aes(x = hrd, y = n)) + 
  stat_pvalue_manual(data = data.frame(
    group1 = "HRD", group2 = "HRP",
    p = prop.test(table(dat$hrd))$p.value,
    y.position = max(table(dat$hrd)) * 1.2),
    label = "prop.test, p = {round(p, 3)}", bracket.size = 1,
    size = 6, tip.length = 0.1) +
  geom_col(width = 0.5, aes(fill = hrd)) + 
  scale_y_continuous(limits = c(0, 80)) +
  theme_light(base_size = 20) +
  scale_fill_manual(values = c("deepskyblue4", "orange"), guide = "none")


####### done

#-f DBS
{
  fig2g <- ggarrange(
    mut_dbs_count |>
      ggplot(aes(x = factor(ctdna_lev, levels = c('neg', 'low', 'high')),
                 y = DBS6/genome_size)) +
      geom_violin(aes(fill=factor(ctdna_lev, levels = c('neg', 'low', 'high'))), 
                  show.legend = F, alpha=.9) +
      geom_boxplot(width=.2, 
                   show.legend = F, alpha=.7) +
      stat_compare_means(size = 4, label.y = .055) +
      stat_compare_means(comparisons = list(c('neg', 'low'),
                                            c('low', 'high'),
                                            c('neg', 'high')),
                         label = 'p.signif',
                         size=4, tip.length = .01,
                         label.y.npc = .7,
                         label.x.npc = 'center') +
      scale_fill_manual(values = list_color) +
      labs(x='DBS6 count per Megabase', y ='') +
      theme_pubr() +
      theme(panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_text(size=12), 
            axis.text.y = element_text(size=10),
            axis.text.x = element_blank(),
            legend.position = c(.85,.75),
            legend.text = element_text(size=10),
            legend.title = element_text(size=15),
            strip.text.x = element_text(size = 15)) ,
    
    mut_dbs_count |>
      ggplot(aes(x = factor(ctdna_lev, levels = c('neg', 'low', 'high')),
                 y = DBS4/genome_size)) +
      geom_violin(aes(fill=factor(ctdna_lev, levels = c('neg', 'low', 'high'))), 
                  show.legend = F, alpha=.9) +
      geom_boxplot(width=.2, 
                   show.legend = F, alpha=.7) +
      #ylim(c(0,.1)) +
      stat_compare_means(size = 4, label.y = .05) +
      stat_compare_means(comparisons = list(c('neg', 'low'),
                                            c('low', 'high'),
                                            c('neg', 'high')),
                         label = 'p.signif',
                         size=4, tip.length = .01,
                         label.y.npc = .7,
                         label.x.npc = 'center') +
      scale_fill_manual(values = list_color) +
      labs(x='DBS9 count per Megabase', y ='') +
      theme_pubr() +
      theme(panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_text(size=12), 
            axis.text.y = element_text(size=10),
            axis.text.x = element_blank(),
            legend.position = c(.85,.75),
            legend.text = element_text(size=10),
            legend.title = element_text(size=15),
            strip.text.x = element_text(size = 15)) ,
    ncol=2)
fig2g
}

svg('results/figures/submission/Main/fig2e.svg', width = 6, height = 5)
fig2e
dev.off()

svg('results/figures/submission/Main/fig2f.svg', width = 7, height = 3)
fig2f
dev.off()

svg('results/figures/submission/Main/fig2g.svg', width = 7, height = 3)
fig2g
dev.off()



######### ------------- CNV integration
{
ascat_break <- ascat_break |>
  mutate(patient = sub('_.*', '', sample),
         tissue = stringr::str_match(sample, '^[^_]+_[piro]\\d?(LN|[A-Z][a-z]+)')[, 2],
         tissue_gr = ifelse(tissue %in% c('Ome', 'Per'), tissue, 
                            ifelse(tissue %in% c('Ova', 'Tub', 'Adn', 'And'),
                                   'Ova_Tub','Other'))) |>
  filter(sample %in% dna_cohort$sample) |>
  left_join(ctdna_cohort |>
              dplyr::select(patient, ctdna_lev))

ascat_break |> 
  left_join(dna_cohort |>
              dplyr::select(sample, purity, TP53.VAF)) |>
  left_join(ctdna_cohort |>
              dplyr::select(patient, Stage)) |>
  #filter(tissue_gr!='Other') |>
  dplyr::slice_max(purity, by = patient) |> 
  dplyr::slice_max(TP53.VAF, by = patient) |> #dim()
  #filter(sample %in%
  #         (dna_cohort |>
  #            dplyr::slice_max(purity, by = c(patient, tissue_gr)) |>
  #            dplyr::slice_max(goodnessOfFit, by = c(patient, tissue_gr)))$sample) |>
  dplyr::group_by(patient, tissue_gr, ctdna_lev, Stage) |>
  dplyr::summarise(breaks = median(ascatBreaks)) |>
  #ggscatter(y='TF', x = 'purity', cor.coef = T, add='reg.line') |>
  ggplot(aes(x=factor(ctdna_lev), 
             y = breaks)) +
  geom_boxplot(aes(fill=factor(ctdna_lev))) +
  geom_jitter() +
  stat_compare_means(#label.x.npc = .5, 
                     #label.y.npc = .8,
                     size=5, label.sep = '\n') +
  scale_fill_manual(values = list_color) +
  #facet_wrap(~tissue_gr, scales = 'free') +
  labs(fill='ctDNA level', x='', y = 'Number of Break Points') +
  theme_pubr() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size=12), 
        axis.text.y = element_text(size=10),
        axis.text.x = element_blank(),
        #legend.position = c(.85,.75),
        legend.text = element_text(size=10),
        legend.title = element_text(size=15),
        strip.text.x = element_text(size = 15))
}
# mutation information integration
{
  library("rstatix")  # https://github.com/kassambara/rstatix
    mut_sig_mat <- mut_sbs_prop |>
      dplyr::select(patient, SBS1, SBS5, SBS40) |>
      left_join(mut_dbs_prop |>
                  dplyr::select(patient, DBS1:DBS11)) |>
      left_join(mut_id_prop |>
                  dplyr::select(patient, ID1:ID10)) |>
      tibble::column_to_rownames(var='patient') |>
      dplyr::select(SBS1,SBS5,SBS40,
                    DBS2,DBS4,DBS6,DBS9,
                    ID5:ID9)
    
    annot <- mut_sig_mat |>
      tibble::rownames_to_column(var='patient') |>
      dplyr::select(patient) |>
      left_join(ctdna_cohort |>
                  dplyr::select(patient, TF)) |>
      tibble::column_to_rownames(var='patient')
      
    pheatmap::pheatmap(t(mut_sig_mat), scale = 'row')
    
}

{
  # heatmap for SBS and DBS signatures
  heat_annot <- mut_sbs_prop |>
    left_join(ctdna_cohort |>
                dplyr::select(patient, age, Stage)) |>
    left_join(df_remain_surv |>
                dplyr::select(patient, hrd)) |>
    dplyr::arrange(by=TF) |>
    tibble::column_to_rownames(var='patient') |>
    dplyr::select(TF, age, hrd, Stage)
  
  sbs_heat <- mut_sbs_prop |>
    dplyr::arrange(by=TF) |>
    tibble::column_to_rownames(var='patient') |>
    dplyr::select(SBS1:SBS13, SBS40)
  
  pat_lev <- factor(mut_sbs_prop$patient)
  sbs_dat1 <- mut_sbs_prop |>
    dplyr::arrange(by=TF) |>
    dplyr::select(patient, TF, ctdna_lev, SBS1, SBS3, SBS5, SBS40) |>
    tidyr::gather(key = 'signature', value = 'proportion', SBS1:SBS40) 
  
  comp2 <- list(c('< 3%', '≥ 3%'))
  
  fig2b <- ggarrange(
    # signature
    ggarrange(
      mut_sbs_prop |>
        ggplot(aes(x = factor(ctdna_lev, levels = c('< 3%', '≥ 3%')), 
                   y = SBS5)) +
        geom_boxplot(aes(color=factor(ctdna_lev, 
                                      levels = c('< 3%', '≥ 3%'))),
                     width = .4,
                     alpha=.4) +
        geom_jitter(aes(color=factor(ctdna_lev, levels = c('< 3%', '≥ 3%')))) +
        stat_compare_means(aes(label=..p.adj..), 
                           comparisons=comp2,
                           label.y = .5, 
                           tip.length = 0) +
        scale_fill_brewer(palette = 'Set1', direction = -1) +
        scale_color_brewer(palette = 'Set1', direction = -1) +
        labs(x='', y ='SBS5 Prop', color='') +
        #ylim(c(0.,.55)) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_text(size=15), 
              axis.text.y = element_text(size=10),
              axis.text.x = element_blank(),
              legend.position = 'none',
              legend.text = element_text(size=10),
              legend.title = element_text(size=10),
              strip.text.x = element_text(size=12)) ,
      
      mut_sbs_count |>
        ggplot(aes(x = factor(ctdna_lev, levels = c('< 3%', '≥ 3%')), 
                   y = SBS5/1e3)) +
        geom_boxplot(aes(color=factor(ctdna_lev, 
                                      levels = c('< 3%', '≥ 3%'))),
                     width = .5,
                     alpha=.4) +
        geom_jitter(aes(color=factor(ctdna_lev, levels = c('< 3%', '≥ 3%')))) +
        stat_compare_means(aes(label=..p.adj..), 
                           comparisons=comp2,
                           label.y = 5) +
        scale_fill_brewer(palette = 'Set1', direction = -1) +
        scale_color_brewer(palette = 'Set1', direction = -1) +
        labs(x='', y ='N SBS5/1kb', color='') +
        #ylim(c(0.025,.25)) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_text(size=15), 
              axis.text.y = element_text(size=10),
              axis.text.x = element_blank(),
              legend.position = 'none',
              legend.text = element_text(size=10),
              legend.title = element_text(size=10),
              strip.text.x = element_text(size=12)) ,
      
      mut_sbs_prop |>
        ggplot(aes(x = factor(ctdna_lev, levels = c('< 3%', '≥ 3%')), 
                   y = SBS1)) +
        geom_boxplot(aes(color=factor(ctdna_lev, 
                                      levels = c('< 3%', '≥ 3%'))),
                     width = .5,
                     alpha=.4) +
        geom_jitter(aes(color=factor(ctdna_lev, levels = c('< 3%', '≥ 3%')))) +
        stat_compare_means(aes(label=..p.adj..), 
                           comparisons=comp2,
                           label.y = .16) +
        scale_fill_brewer(palette = 'Set1', direction = -1) +
        scale_color_brewer(palette = 'Set1', direction = -1) +
        labs(x='', y ='SBS1 Prop', color='') +
        #ylim(c(0,.5)) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_text(size=15), 
              axis.text.y = element_text(size=10),
              axis.text.x = element_blank(),
              legend.position = 'none',
              legend.text = element_text(size=10),
              legend.title = element_text(size=10),
              strip.text.x = element_text(size=12)) ,
      
      mut_sbs_count |>
        ggplot(aes(x = factor(ctdna_lev, levels = c('< 3%', '≥ 3%')), 
                   y = SBS1/1e3)) +
        geom_boxplot(aes(color=factor(ctdna_lev, 
                                      levels = c('< 3%', '≥ 3%'))),
                     width = .5,
                     alpha=.4) +
        geom_jitter(aes(color=factor(ctdna_lev, levels = c('< 3%', '≥ 3%')))) +
        stat_compare_means(aes(label=..p.adj..), 
                           comparisons=comp2,
                           label.y = 1.6) +
        scale_fill_brewer(palette = 'Set1', direction = -1) +
        scale_color_brewer(palette = 'Set1', direction = -1) +
        labs(x='', y ='N SBS1/1kb', color='') +
        #ylim(c(0,.5)) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_text(size=15), 
              axis.text.y = element_text(size=10),
              axis.text.x = element_blank(),
              legend.position = 'none',
              legend.text = element_text(size=10),
              legend.title = element_text(size=10),
              strip.text.x = element_text(size=12)) ,
      
      mut_sbs_count |>
        ggplot(aes(x = factor(ctdna_lev, levels = c('< 3%', '≥ 3%')), 
                   y = SBS40/1e3)) +
        geom_boxplot(aes(color=factor(ctdna_lev, 
                                      levels = c('< 3%', '≥ 3%'))),
                     width = .5,
                     alpha=2) +
        geom_jitter(aes(color=factor(ctdna_lev, levels = c('< 3%', '≥ 3%')))) +
        stat_compare_means(aes(label=..p.adj..), 
                           comparisons=comp2,
                           label.y = 12) +
        scale_fill_brewer(palette = 'Set1', direction = -1) +
        scale_color_brewer(palette = 'Set1', direction = -1) +
        labs(x='', y ='N SBS40/1kb', color='') +
        #ylim(c(0,.5)) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_text(size=15), 
              axis.text.y = element_text(size=10),
              axis.text.x = element_blank(),
              legend.position = 'none',
              legend.text = element_text(size=10),
              legend.title = element_text(size=10),
              strip.text.x = element_text(size=12)) ,
      
      mut_dbs_count |>
        ggplot(aes(x = factor(ctdna_lev, levels = c('< 3%', '≥ 3%')), 
                   y = DBS2/1e3)) +
        geom_boxplot(aes(color=factor(ctdna_lev, 
                                      levels = c('< 3%', '≥ 3%'))),
                     width = .5,
                     alpha=2) +
        geom_jitter(aes(color=factor(ctdna_lev, levels = c('< 3%', '≥ 3%')))) +
        stat_compare_means(aes(label=..p.adj..), 
                           comparisons=comp2,
                           label.y = .07) +
        scale_fill_brewer(palette = 'Set1', direction = -1) +
        scale_color_brewer(palette = 'Set1', direction = -1) +
        labs(x='', y ='N DBS2/1kb', color='') +
        #ylim(c(0,.5)) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_text(size=15), 
              axis.text.y = element_text(size=10),
              axis.text.x = element_blank(),
              legend.position = 'none',
              legend.text = element_text(size=10),
              legend.title = element_text(size=10),
              strip.text.x = element_text(size=12)) ,
      
      mut_dbs_count |>
        ggplot(aes(x = factor(ctdna_lev, levels = c('< 3%', '≥ 3%')), 
                   y = DBS4/1e3)) +
        geom_boxplot(aes(color=factor(ctdna_lev, 
                                      levels = c('< 3%', '≥ 3%'))),
                     width = .5,
                     alpha=2) +
        geom_jitter(aes(color=factor(ctdna_lev, levels = c('< 3%', '≥ 3%')))) +
        stat_compare_means(aes(label=..p.adj..), 
                           comparisons=comp2,
                           label.y = .07) +
        scale_fill_brewer(palette = 'Set1', direction = -1) +
        scale_color_brewer(palette = 'Set1', direction = -1) +
        labs(x='', y ='N DBS4/1kb', color='') +
        #ylim(c(0,.5)) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_text(size=15), 
              axis.text.y = element_text(size=10),
              axis.text.x = element_blank(),
              legend.position = 'none',
              legend.text = element_text(size=10),
              legend.title = element_text(size=10),
              strip.text.x = element_text(size=12)) ,
      
      mut_dbs_count |>
        ggplot(aes(x = factor(ctdna_lev, levels = c('< 3%', '≥ 3%')), 
                   y = DBS6/1e3)) +
        geom_boxplot(aes(color=factor(ctdna_lev, 
                                      levels = c('< 3%', '≥ 3%'))),
                     width = .5,
                     alpha=2) +
        geom_jitter(aes(color=factor(ctdna_lev, levels = c('< 3%', '≥ 3%')))) +
        stat_compare_means(aes(label=..p.adj..), 
                           comparisons=comp2,
                           label.y = .07) +
        scale_fill_brewer(palette = 'Set1', direction = -1) +
        scale_color_brewer(palette = 'Set1', direction = -1) +
        labs(x='', y ='N DBS6/1kb', color='') +
        #ylim(c(0,.5)) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_text(size=15), 
              axis.text.y = element_text(size=10),
              axis.text.x = element_blank(),
              legend.position = 'none',
              legend.text = element_text(size=10),
              legend.title = element_text(size=10),
              strip.text.x = element_text(size=12)) ,
      nrow=2, ncol=4
    ) ,
    # mutation count/vaf
    ggarrange(
      mut_gene_pat_lev |>
        filter(gene != 'TP53') |>
        dplyr::group_by(patient, Func, TF, ctdna_lev) |>
        dplyr::summarise(N=n()) |>
        as.data.frame() |>
        filter(Func %in% c('synonymous', 'stopgain', 'frameshift_indel', 'splicing', 'inframe_indel')) |>
        mutate(ctdna_lev = factor(ctdna_lev, levels = c('< 3%', '≥ 3%'))) |>
        ggscatter(x = 'TF', y = 'N', add='reg.line', 
                  cor.coef = T, cor.method = 'spearman',
                  conf.int = T,
                  cor.coef.size = 5) + 
        geom_point(aes(color=ctdna_lev),
                   show.legend = F,
                   alpha=2) +
        facet_grid(~factor(Func,levels = c('synonymous', 'frameshift_indel', 'stopgain', 'splicing', 'inframe_indel'),
                           labels = c('Synonymous', 'Frameshift', 'Stopgain', 'Splicing', 'Inframe')), 
                   scales = 'free') +
        geom_vline(xintercept = .03, linetype='dotted', color='blue', linewidth=1)+
        #geom_vline(xintercept = .07, linetype='dotted', color='blue', linewidth=1) +
        labs( x = 'ctDNA Tumor Fraction') +
        scale_color_brewer(palette = 'Pastel1', direction = -1) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_text(size=15), 
              axis.text.y = element_text(size=10),
              axis.text.x = element_blank(),
              legend.position = 'none',
              legend.text = element_text(size=10),
              legend.title = element_text(size=15),
              strip.text.x = element_text(size = 15)) ,
      
      mut_call |> 
        filter(gene!='TP53',
               Func %in% c('synonymous', 'stopgain', 'frameshift_indel', 'splicing', 'inframe_indel')
               ) |>
        ggplot(aes(x=VAF)) + 
        geom_density(aes(fill=factor(ctdna_lev, levels = c('< 3%', '≥ 3%'))), alpha=.4) + 
        labs(fill='', x='VAF calling from Tissue Biopsy') + 
        scale_fill_brewer(palette = 'Set1', direction = -1) + 
        facet_grid(~factor(Func,levels = c('synonymous', 'frameshift_indel', 'stopgain', 'splicing', 'inframe_indel'),
                           labels = c('Synonymous', 'Frameshift', 'Stopgain', 'Splicing', 'Inframe')), 
                   scales = 'free') +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_text(size=15), 
              axis.text.y = element_text(size=10),
              axis.text.x = element_blank(),
              legend.position = 'none',
              legend.text = element_text(size=10),
              legend.title = element_text(size=15),
              strip.text = element_blank()#element_text(size = 20)
        ) ,
      ncol=1, nrow=2
    ),
    ncol = 1, nrow = 2, heights = c(.45,.55) , 
    labels = c('d', 'e'), font.label = list(size = 20)
  )
}
svg('results/figures/submission/Supp_3/fig2b.svg', width = 12, height = 10)
fig2b
dev.off()


{
  fig2b <- ggarrange(
    mut_gene_pat_lev |>
      filter(gene != 'TP53') |>
      dplyr::group_by(patient, TF, ctdna_lev) |>
      dplyr::summarise(N=n()/1e3) |>
      as.data.frame() |>
      mutate(ctdna_lev = factor(ctdna_lev, levels = c('< 3%', '≥ 3%'))) |>
      ggscatter(x = 'TF', y = 'N', add='reg.line', 
                cor.coef = T, cor.method = 'spearman',
                conf.int = T,
                cor.coef.size = 8) + 
      geom_point(aes(color=ctdna_lev),
                 show.legend = F,
                 alpha=2) +
      geom_vline(xintercept = .03, linetype='dotted', color='blue', linewidth=1)+
      #geom_vline(xintercept = .07, linetype='dotted', color='blue', linewidth=1) +
      labs( x = 'ctDNA Tumor Fraction') +
      scale_color_brewer(palette = 'Set1', direction = -1) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_text(size=15), 
            axis.text.y = element_text(size=10),
            axis.text.x = element_blank(),
            legend.position = 'none',
            legend.text = element_text(size=10),
            legend.title = element_text(size=15),
            strip.text.x = element_text(size = 15)) ,
    
    # proportion of low vaf mutations
    mut_call |>
      dplyr::select(chr, pos, ref, alt, gene, VAF, patient, TF) |>
      distinct() |>
      mutate(low_vaf = ifelse(VAF < .05, 'yes', 'no')) |>
      dplyr::group_by(patient, TF, low_vaf) |>
      dplyr::summarise(N = n()) |>
      as.data.frame() |>
      dplyr::group_by(patient, TF) |>
      mutate(p = N/sum(N)) |>
      as.data.frame() |>
      filter(low_vaf=='yes') |>
      ggplot(aes(x=ctdna_detection, y = p)) +
      geom_violin() +
      geom_jitter() +
      stat_compare_means()
      
    
  )
}


# figure 2c: 
{
  eoc.score <- read.delim('~/mnt/storageBig8/work/nguyenma/projects/RNA/OmicsIntegration/data/processed/Progeny_scores_ctDNA_x_RNA_EOC.tsv') |>
    dplyr::rename(sample=X) |>
    mutate(patient = sub('_.*', '', sample)) |>
    left_join(rna_cohort |>
                dplyr::select(sample, tissue, EOC, Fibroblast, Immune)) |>
    inner_join(ctdna_cohort |>
                 dplyr::select(patient, TF)) |>
    mutate(ctdna_lev = ifelse(TF>=.03, 'high', 'low')) |>
    dplyr::arrange(TF)
  
  heatmap_df <- eoc.score %>%
    tibble::column_to_rownames(var='sample') |>
    dplyr::select(Androgen:p53) 
  annot_heat <- eoc.score %>%
    tibble::column_to_rownames(var='sample') |>
    dplyr::select(TF, ctdna_lev)
  
  ComplexHeatmap::Heatmap(heatmap_df, cluster_rows = F)
  
  eoc.score |>
    ggplot(aes(x=ctdna_lev, y = PI3K)) +
    geom_boxplot() +
    geom_jitter() +
    #geom_errorbar() +
    stat_compare_means()
  
  eoc.score.df <- eoc.score |>
    mutate(tissue = stringr::str_match(sample, '^[^_]+_[piro]\\d?(LN|[A-Z][a-z]+)')[, 2]) |>
    filter(!(tissue %in% c('Asc', 'Plf', 'LN', 'Lnr'))) |>
    dplyr::select(patient, Hypoxia, PI3K, TNFa, WNT, JAK.STAT) |>
    dplyr::group_by(patient) |>
    dplyr::summarise_all(median) |>
    as.data.frame()
  
  intg_df <- mut_sbs_prop |> 
    dplyr::select(patient, TF, SBS1.p=SBS1, 
                  SBS3.p=SBS3, SBS5.p=SBS5, 
                  SBS40.p=SBS40) |> 
    left_join(
      mut_sbs_count %>% 
        mutate_at(colnames(.)[grepl('SBS', colnames(.))], function(x) x/genome_size) %>%
        dplyr::select(patient, TF, SBS1.c=SBS1, 
                      SBS3.c=SBS3, SBS5.c=SBS5, 
                      SBS40.c=SBS40)
    ) |>
    left_join(
      mut_dbs_prop |> 
        dplyr::select(patient, DBS2.p=DBS2, DBS9.p=DBS9,
                      DBS4.p=DBS4, DBS6.p=DBS6)) |>
    left_join(
      mut_dbs_count %>%
        mutate_at(colnames(.)[grepl('DBS', colnames(.))], function(x) x/genome_size) %>%
        dplyr::select(patient, DBS2.c=DBS2, DBS9.c=DBS9,
                      DBS4.c=DBS4, DBS6.c=DBS6)) |>
    #inner_join(eoc.score.df) |>
    tibble::column_to_rownames(var='patient') #%>%
    filter(rownames(.)!='H084')
  
  ann.intg <- intg_df |>
    tibble::rownames_to_column(var='patient') |>
    dplyr::select(patient) |>
    left_join(ctdna_cohort |>
                mutate(ctdna_lev = ifelse(TF>=.03, 'pos', 'neg')) %>%
                dplyr::select(patient, ctdna_lev, progress, HRD, Stage, age, PFI, OS)) |>
    tibble::column_to_rownames(var='patient')
  
  pheatmap::pheatmap(intg_df[,2:ncol(intg_df)], annotation_row = ann.intg)  
  #### kmeans
  k=3
  kmean_res <- kmeans(intg_df[,-1], centers = k, iter.max = 1000)
  res <- intg_df |>
    mutate(cluster = paste0('Cluster', kmean_res$cluster))
  res |> ggboxplot(x='cluster', y = 'TF', add='jitter')
  
  {
    # correlation matrix
    corrplot::corrplot(cor(intg_df), type = 'upper', insig = 'label_sig')
    
    ComplexHeatmap::Heatmap(as.matrix(scale(intg_df)), row_km = 3)

  pheatmap::pheatmap(intg_df[,c('TF','DBS6.c','SBS1.c','SBS40.c')], scale = 'column', annotation_row = ann.intg, 
                     colorRampPalette(c("#000080", "white", "#cc0000"))(80))
  corrplot::corrplot(cor(intg_df))
  }
  
  intg_df[1,]
  library(psych)
  pairs.panels(intg_df[, 1:ncol(intg_df)],smooth = TRUE,      # If TRUE, draws loess smooths
               scale = FALSE,      # If TRUE, scales the correlation text font
               density = TRUE,     # If TRUE, adds density plots and histograms
               ellipses = TRUE,    # If TRUE, draws ellipses
               method = "spearman", # Correlation method (also "spearman" or "kendall")
               pch = 21,           # pch symbol
               lm = TRUE,         # If TRUE, plots linear fit rather than the LOESS (smoothed) fit
               cor = TRUE,         # If TRUE, reports correlations
               jiggle = FALSE,     # If TRUE, data points are jittered
               factor = 2,         # Jittering factor
               hist.col = 4,       # Histograms color
               stars = TRUE,       # If TRUE, adds significance level with stars
               ci = TRUE)
  
  # linear mixed model
  {
    library(nlme)
    library(lme4)
    
    lm1 <- lm(TF ~ DBS6.c, 
              data=intg_df)
    
    lm2 <- lm(TF ~ SBS5 + clonality + SBS5 * clonality + 
                DBS2 + DBS4 + DBS6 + sig3 + SBS1, data=intg_df)
    
    lm3 <- lm(TF ~ SBS5 + DBS6 + sig3 + 
                clonality +
                SBS5 * clonality + 
                DBS2 + DBS4 , 
              data=intg_df)
    
    lm4 <- lm(TF ~ SBS5.p + SBS1.p + DBS2.c + DBS9.c + SBS1.p +
                DBS6.c + DBS4.c, data=intg_df)
    summary(lm4)
    anova(lm1, lm4)
    tab_model(lm4, 
              show.se = T,
              #show.p = T,
              show.aic = T,
              show.loglik = T,
              show.fstat = T,
              show.stat = T)
    lm5 <- lm(TF ~ PI3K + Hypoxia + TNFa +
                DBS2 + DBS4 , 
              data=intg_df)
    summary(lm5)
    glm1 <- gls(TF ~ SBS5 + clonality + SBS5 * clonality, data = intg_df,
                method = "ML")
    
    glm2 <- gls(TF ~ SBS5 + clonality + SBS5 * clonality + 
                  DBS2 + DBS4 + DBS6 + sig3 + SBS1, data = intg_df,
                method = "ML")
    
    glm3 <- gls(TF ~ SBS5 + DBS6 + sig3 +
                  clonality + DBS2 + DBS4,
                data = intg_df,
                method = "ML")
    tab_model(glm3, 
              show.std = T,
              show.p = T)
    summary(glm3)
    lme1 <- lme(TF ~ SBS5 + clonality + SBS5 * clonality, 
                data = intg_df,
                random = ~1|sig3, method = "ML")
    
    lme2 <- lme(TF ~ SBS5 + clonality + SBS5 * clonality + 
                  DBS2 + DBS4 + DBS6  + SBS1, 
                data = intg_df,
                random = ~1|sig3, method = "ML")
    
    lme3 <- lme(TF ~ SBS5 + 
                  clonality + DBS2 + DBS4 +
                  SBS5 * clonality + 
                  DBS6, 
                data = intg_df,
                random = ~1|sig3, method = "ML")
    summary(lme3)
    anova(glm3, glm1)
    anova(glm1, glm2, glm3, lme1, lme2, lme3)
    plot(ranef(lme1, level = 1))
    summary(glm3)
    devtools::install_github("strengejacke/sjPlot")
    library(sjPlot)
    tab_model(glm3)
  }
  
  svg('results/figures/submission/Supp_3/fig2c.svg', width = 12, height = 12)
  pairs.panels(intg_df[, 2:ncol(intg_df)],
               smooth = TRUE,      # If TRUE, draws loess smooths
               scale = FALSE,      # If TRUE, scales the correlation text font
               density = TRUE,     # If TRUE, adds density plots and histograms
               ellipses = F,    # If TRUE, draws ellipses
               method = "spearman", # Correlation method (also "spearman" or "kendall")
               pch = 21,           # pch symbol
               lm = TRUE,         # If TRUE, plots linear fit rather than the LOESS (smoothed) fit
               cor = TRUE,         # If TRUE, reports correlations
               jiggle = FALSE,     # If TRUE, data points are jittered
               factor = 1,         # Jittering factor
               hist.col = 4,       # Histograms color
               stars = TRUE,       # If TRUE, adds significance level with stars
               ci = TRUE)
  dev.off()
}


# figure 2b: evolution integration -- skip for now
{
  evol_df <- evol_gr_annot %>%
    inner_join(ctdna_cohort %>%
                 dplyr::select(patient, TF)) %>%
    mutate(clonality = as.numeric(clonality),
           divergence = as.numeric(divergence)) %>%
    mutate(ctdna_lev = ifelse(TF < .01, '< 1%',
                              ifelse(TF < 0.07, '1-7%', '> 7%')),
           ctdna_cat = ifelse(TF <= .03, '≤ 3%', '< 3%')) |>
    filter(patient %in% dna_cohort$patient)
  
  nrow(evol_df)
  comp1 <- list(#c('adaptive', 'evolving'),
                c('adaptive', 'maintaining'),
                c('maintaining', 'evolving'))
  comp2 <- list(c('< 1%', '> 7%'),
                c('1-7%', '> 7%'))
  fig2a <- ggarrange(
    # - Association between evolution groups
    evol_df |>
      #ggplot(aes(x=state, y = TF)) +
      #geom_boxplot(aes(color=state), show.legend = F) +
      ggplot(aes(x = factor(state, 
                            levels = c('adaptive', 'evolving', 'maintaining')))) +
      geom_bar(aes(fill = factor(ctdna_lev, levels=c('< 3%', '≥ 3%'))), 
               alpha=2, position = 'fill', show.legend = F) +
      scale_fill_brewer(palette = 'Pastel1', direction = -1) +
      labs(x='Evolution State', y = 'Prop') +
      theme_classic() +
      theme(axis.title = element_text(size=20), 
            axis.text = element_text(size=15),
            legend.position = c(.85,.8),
            legend.text = element_text(size=10),
            legend.title = element_text(size=15)),
    
    evol_df |>
      ggplot(aes(x=factor(ctdna_lev, levels=c('< 3%', '≥ 3%')), 
                 y=clonality)) +
      geom_violin(aes(color=factor(ctdna_lev, levels=c('< 3%', '≥ 3%'))), 
                  show.legend = F,
                  trim = FALSE)+
      geom_dotplot(aes(fill=factor(ctdna_lev,
                                   levels=c('< 3%', '≥ 3%'))),
                   show.legend = F, alpha=2,
                   binaxis='y', stackdir='center',
                   stackratio=1.2, dotsize=.7) +
      #stat_summary(fun=mean, geom="point", shape=15,
      #             size=5, color="red") +
      stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                   geom="pointrange", color="#ff0000") +
      stat_compare_means(comparisons = comp2,
                         label = 'p.signif',
                         label.y = c(6.6, 6.1)) +
      scale_fill_brewer(palette = 'Pastel1', direction = -1) +
      scale_color_brewer(palette = 'Pastel1', direction = -1) +
      labs(x='ctDNA Level', y='Clonality Score') +
      lims(y=c(0,7))+
      theme_classic() +
      theme(axis.title = element_text(size=20), 
            axis.text = element_text(size=15),
            legend.position = c(.85,.8),
            legend.text = element_text(size=10),
            legend.title = element_text(size=15)),
  
    # - Association of ctDNA TF and evolution features
    evol_df |>
      ggscatter(x = 'TF', y='clonality', add='reg.line', cor.coef = T,cor.method = 'spearman',
                conf.int = T, cor.coef.size = 6, cor.coef.coord = c(.05, 6.1),
                color='#999999') +
      labs(x = 'ctDNA Fraction', y = '') +
      lims(y=c(0,7))+
      theme(axis.title = element_text(size=20), 
            axis.text = element_text(size=15)) +
      #geom_vline(xintercept = .01, linetype='dotted', color='blue', linewidth=1)+
      geom_vline(xintercept = .03, linetype='dotted', color='blue', linewidth=1),
    
    nrow = 1, ncol = 3, widths = c(.35,.3,.35),
    labels = c('c', 'd', 'e'), 
    font.label = list(size = 20, color = "black")
  )
}
svg('results/figures/submission/Supp_3/fig2b.svg', width = 12, height = 5)
fig2a
dev.off()

#### radar chart for signature proportions - skip
{
  nb.cols <- 3
  mycolors <- colorRampPalette(RColorBrewer::brewer.pal(nb.cols, "Pastel1"), alpha=.6)(nb.cols)
  mycolors <- c("#4DAF4A", "#377EB8", "#E41A1C")
  area_col <- c("#CCEBC5FF", "#B3CDE3FF", "#FBB4AEFF")
  radar_df <- mut_sbs_prop |>
    dplyr::select(patient, ctdna_lev, SBS1, SBS5, SBS40) |>
    left_join(mut_dbs_prop |>
                dplyr::select(patient, ctdna_lev, DBS2, DBS4, DBS6)) |>
    dplyr::group_by(ctdna_lev) |>
    summarise_at(vars(SBS1:DBS6), mean) |>
    as.data.frame() |>
    tibble::column_to_rownames(var = 'ctdna_lev')
  max_min <- radar_df |>
    mutate_all(max) |>
    head(1) |>
    rbind(radar_df |>
            mutate_all(min) |>
            head(1))
  x=min(max_min)
  y=max(max_min)
  rownames(max_min) <- c('max', 'min')
  radarchart(max_min |> rbind(ttt[c(2,1,3),]),maxmin = T, 
             pcol = mycolors, plwd = 3, plty = 1,
             pfcol = c(alpha("#CCEBC5FF", 0.3), alpha("#B3CDE3FF", 0.3), alpha("#FBB4AEFF",.3)),
             axislabcol = 1, caxislabels = c(8,23,38,54))
}


{
  tt <- mut_call |> 
    dplyr::select(chr, start=pos, end=pos, vaf=VAF, patient) |> 
    distinct() |>
    #filter(vaf < .05) |>
    mutate(chrom = sub('chr','',chr)) |>
    dplyr::arrange(chrom, start)
  
  thresh <- median(ctdna_cohort$TF)
  #tt1 <- tt |> filter(patient %in% ctdna_cohort[ctdna_cohort$TF>thresh,]$patient)
  #tt2 <- tt |> filter(patient %in% ctdna_cohort[ctdna_cohort$TF<=thresh,]$patient)
  #circos.initializeWithIdeogram(species='hg38', chromosome.index = paste0("chr", 1:22))
  #circos.genomicDensity(tt1, col = c("#ffdab9"), track.height = 0.1)
  #circos.genomicDensity(tt2, col = c("#c6e2ff"), track.height = 0.1)
  
  tt1 <- tt |> filter(patient %in% ctdna_cohort[ctdna_cohort$TF>=.03,]$patient)
  tt2 <- tt |> filter(patient %in% ctdna_cohort[ctdna_cohort$TF<.03,]$patient)
  
  circos.initializeWithIdeogram(species='hg38', chromosome.index = paste0("chr", 1:22))
  circos.genomicDensity(tt1, col = '#FF000080')
  circos.genomicDensity(tt2, col = c("#0000FF80"), track.height = 0.1)
  #circos.genomicDensity(tt1, col = c("#66CDAA"), track.height = 0.1)
}

# figure 2d -- skip for now:
{
  # - CNA burden
  cna_wgs <- wgs_data_prim |>
    filter(sample %in% dna_cohort$sample) |>
    left_join(ctdna_cohort |>
                dplyr::select(-sample)) |>
    mutate(ctdna_lev = ifelse(TF < .03, '< 3%', '≥ 3%'))
  #length(unique(cna_df$sample))==length(unique(dna_cohort$sample))
  #length(unique(cna_df$patient))==100
  
  # remove outlier
  remove_outliers <- function(x, na.rm = TRUE, ...) {
    qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
    H <- 1.5 * IQR(x, na.rm = na.rm)
    y <- x
    y[x < (qnt[1] - H)] <- NA
    y[x > (qnt[2] + H)] <- NA
    y
  }
  
  cna_dat1 <- cna_wgs %>%
    dplyr::slice_max(purity, by = c(patient, tissue)) |>
    filter(!(grepl('seq', sample))) |> 
    as.data.frame() |>
    mutate(tissue_gr = ifelse(tissue %in% c('Tub', 'Ova'), 'Ova_Tub',
                          ifelse(tissue %in% c('Ome', 'Per'), tissue, 'Others')))
  
  cna_dat2 <- cna_wgs %>%
    dplyr::slice_max(purity, by = patient) |>
    dplyr::slice_max(ploidy, by = patient) |>
    filter(!(grepl('seq', sample))) |> 
    as.data.frame()
  
  my_func1 <- function(df) {
    new <- remove_outliers(df$logR)
    return(cbind(new,df))
  }
  
  logR_rm_outlier <- cna_dat2 %>% 
    group_by(chromosome, ctdna_lev) %>% 
    nest() %>%
    mutate(data = purrr::map(data, my_func1)) %>%
    unnest(cols = c(data)) |>
    as.data.frame()
  
  my_func2 <- function(df) {
    new <- remove_outliers(df$CN)
    return(cbind(new,df))
  }
  
  CN_rm_outlier <- cna_dat1 %>% 
    mutate(CN=minorAlleleCopyNumber+majorAlleleCopyNumber) |>
    group_by(chromosome, ctdna_lev) %>% 
    nest() %>%
    mutate(data = purrr::map(data, my_func2)) %>%
    unnest(cols = c(data)) |>
    as.data.frame()
  
  CN_rm_outlier |>
    #filter(chromosome %in% c('chr1', 'chr3', 'chr6', 'chr8',
    #                        'chr12', 'chr17', 'chr19')) |>
    ggplot(aes(x = factor(ctdna_lev,
                          levels=c('< 3%', '≥ 3%')),
               y = new)) + 
    geom_boxplot(aes(color=factor(ctdna_lev,
                                  levels=c('< 3%', '≥ 3%'))), #outlier.shape = NA,
                 show.legend = F) +
    #ylim(c(0,20))+
    facet_wrap(~chromosome, scales = 'free') +
    labs(x='', y='Copy Number Value')+
    stat_compare_means(label.y.npc = .3, size=3, 
                       comparisons = comp2,
                       label = '..p..')
  
  logR_rm_outlier |>
    #filter(chromosome %in% c('chr1', 'chr3', 'chr6', 'chr8',
    #                        'chr12', 'chr17', 'chr19')) |>
    ggplot(aes(x = factor(ctdna_lev,
                          levels=c('< 3%', '≥ 3%')),
               y = new)) + 
    geom_boxplot(aes(color=factor(ctdna_lev,
                                  levels=c('< 3%', '≥ 3%'))), #outlier.shape = NA,
                 show.legend = F) +
    #ylim(c(0,20))+
    facet_wrap(~chromosome, scales = 'free') +
    labs(x='', y='Copy Number Value')+
    stat_compare_means(label.y.npc = .3, size=3, 
                       comparisons = comp2,
                       label = '..p..')
  
  CN_rm_outlier |>
    mutate(CN=minorAlleleCopyNumber+majorAlleleCopyNumber) |> 
    filter(tissue_gr %in% c('Per', 'Ome', 'Ova_Tub'), 
           chromosome %in% c('chr2', 'chr3', 'chr5', 'chr8', 'chr12', 'chr17', 'chr20')) |>
    ggplot(aes(x = CN)) + 
    geom_density(aes(color=factor(ctdna_lev,
                                  levels=c('< 3%', '≥ 3%'))), #outlier.shape = NA,
                 show.legend = T) + 
    facet_grid(tissue_gr~chromosome, scales = 'free') + 
    scale_color_brewer(palette = 'Pastel1', direction = -1) + 
    labs(fill='ctdna_lev') +
    theme(legend.position = 'bottom')
  
  ###
  {
    species = "hg38"
    circos.initializeWithIdeogram(species = species)
  }
  fig2c <- ggarrange(
    
  )
  # - gistic landscape comparison

}

#### Top driver amplifcations

driv_cna <- data.frame(
  geneSymbol = c('MECOM', 'MYC', 'KRAS', 'CCNE1'),
  chr = c('chr3', 'chr8', 'chr12', 'chr19'),
  start = c(169083507, 27735434, 25209178, 19811991),
  end = c(169663712, 127740477, 25250936, 29824312)
) 
driv_cna_gr <- top_cna |>
  GenomicRanges::GRanges()

# infer CN values of those genes 
wgs_segment <- wgs_data_prim |>
  filter(sample %in% dna_cohort$sample) |>
  dplyr::select(chr=chromosome, start, end, copyNumber, 
                minor=minorAlleleCopyNumber,
                major=majorAlleleCopyNumber,
                logR, Loh,
                sample, patient)
segment_gr <- with(wgs_segment, GRanges(chr, IRanges(start, end)))
overlaps <- findOverlaps(driv_cna_gr,segment_gr)
df <- data.frame(
  driv_cna[from(overlaps), ],
  wgs_segment[to(overlaps), ] %>%
    dplyr::select(-chr, -start, -end)
) |>
  dplyr::group_by(patient, sample, geneSymbol, chr, start, end) |>
  dplyr::summarise(CN = min(copyNumber), 
                   minor = min(minor),
                   major = min(major),
                   logR=min(logR),
                   Loh = min(Loh)) |>
  as.data.frame() 
df |> write.table('data/input/WGSdata/driver_CNA.tsv', sep = '\t', 
                  col.names = T, row.names = F, quote = F)




## CCNE1 status - skip
{
  ccne_pat <- readxl::read_excel('data/input/clinicalInfo/CCNE1amp_20230510.xlsx',
                                 sheet = 'CCNE1pat', skip=1) |>
    dplyr::rename(patient = Patient) |>
    as.data.frame() |> inner_join(ctdna_cohort) |>
    filter(patient %in% dna_cohort$patient)
  
  ccne_pat |>
    ggplot(aes(x = CCNE1amp_hetStatus, y = TF)) +
    geom_boxplot() +
    stat_compare_means()
}
}
