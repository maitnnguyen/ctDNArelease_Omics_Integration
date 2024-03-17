source('00_functions.R')

genome_size <- 2900077904/1e6 # (Megabase)
# load WGS tissue mutation data
load('~/OmicsIntegration/data/WGS/mutation_call_4_6.RData')

# mutational signature on sample level
load('~/OmicsIntegration/data/WGS/mutational_signature_sample_lev_231005.RData')

load('~/OmicsIntegration/data/WGS/ascat_results.RData')

dna_cohort <- ascat_est |>
  filter(patient %in% ctdna_cohort$patient) |>
  mutate(tissue_gr = ifelse(tissue %in% c('Ome', 'Per'), tissue, 
                            ifelse(tissue %in% c('Ova', 'Tub', 'Adn', 'And'),
                                   'Ova_Tub','Other'))) |>
  left_join(ctdna_cohort |>
              dplyr::select(patient, TF=ctDNA_fraction, Stage,
                            ctdna_lev, treatment)) 
  
##### fig2a: number of breakpoints
{
  ascat_break <- ascat_break_point |>
    inner_join(ctdna_cohort |>
                dplyr::select(patient, ctdna_lev))
  
  fig2a <- ascat_break |> 
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
    labs(fill='ctDNA level', x='', 
         y = 'Avg # Break Points per Chromosome') +
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
  
  fig2a
}

##### fig2b: mutational burden
{
  fig2b <- mut_pat_lev |>
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
    labs(fill='ctDNA Level', x = '', 
         y='N somatic mutations per Megabase') +
    theme_pubr() +
    theme(axis.title = element_text(size=15), 
          axis.text.y = element_text(size=13),
          axis.text.x = element_blank(),
          legend.position = 'none',
          legend.text = element_text(size=10),
          legend.title = element_text(size=12),
          strip.text.x = element_text(size = 15))
}

##### fig2c: mutation types
{
  fig2c <- mut_sigs |>
    ggplot(aes(x = factor(ctdna_lev, levels = c('low', 'med', 'high')),
               y = mut_count/genome_size)) + 
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
          legend.position = c(.85,.75),
          legend.text = element_text(size=10),
          legend.title = element_text(size=15),
          strip.text.x = element_text(size = 15))
}

##### fig2d: mutational signatures
{
  fig2d1 <- signatures_df |>
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
  
  fig2d2 <- signatures_df |>
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
}

##### fig2e: mutational signatures
{
  fig2e <- signatures_df |>
    filter(signature %in% c('SBS1')) |>
    ggplot(aes(x = factor(ctdna_lev, 
                          levels = c('low', 'med', 'high')),
               y = prop)) +
    geom_violin(aes(fill = factor(ctdna_lev, 
                                  levels = c('low', 'med', 'high'))), 
                show.legend = F, alpha=.9) +
    geom_boxplot(width=.2, outlier.shape = NA,
                 show.legend = F) +
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
}

##### fig2f: mutational signatures
{
  stat.test <- xchisq.test(ctdna_lev ~id5, data=signatures_df)
  fig2f <- signatures_df |>
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
}

################ done ##################
