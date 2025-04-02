source('load_data.R')

##### fig3a: number of breakpoints
{
  ascat_break <- ascat_break_point |>
    inner_join(ctdna_cohort |>
                dplyr::select(patient, ctdna_lev))
  
  fig3a <- ascat_break |> 
    inner_join(dna_cohort |>
                dplyr::select(sample, purity, TP53.VAF)) |>
    filter(sample %in% highest_purity_dna$sample) |> 
    as.data.frame() |> 
    ggplot(aes(x=factor(ctdna_lev, 
                        levels = c('low', 'med', 'high')), 
               y = breaks/genome_size)) +
    geom_boxplot(aes(color=factor(ctdna_lev, 
                              levels = c('low', 'med', 'high'))),
                 outlier.shape = NA,
                 show.legend = F, width=.5) +
    geom_dotplot(aes(fill=factor(ctdna_lev, 
                                  levels = c('low', 'med', 'high'))),
                 show.legend = F,
                 binaxis='y', stackdir='center', dotsize=.7) +
    scale_color_manual(values = list_color) +
    scale_fill_manual(values = list_color) +
    stat_compare_means(comparisons = list(c('low', 'med'),
                                          c('high', 'med'),
                                          c('high', 'low')),
                       label = '.p.', size=5,
                       label.y = c(.27,.30,.33), 
                       tip.length = .01) +
    labs(fill='ctDNA level', x='', y = 'Avg N Breakpoints per Megabase') +
    theme_pubr() +
    theme(axis.title = element_text(size=15), 
          axis.text.y = element_text(size=13),
          axis.text.x = element_blank(),
          legend.position = 'none',
          legend.text = element_text(size=10),
          legend.title = element_text(size=12),
          strip.text.x = element_text(size = 15))
  
  fig3a
}

##### fig3b: mutational burden
{
  fig3b <- mut_pat_lev |>
    dplyr::group_by(patient, ctdna_lev) %>%
    dplyr::summarise(n_mut = n()/genome_size) %>%
    as.data.frame() %>%
    ggplot(aes(x = factor(ctdna_lev, 
                          levels = c('low', 'med', 'high')), 
               y = n_mut)) +
    geom_boxplot(aes(color=factor(ctdna_lev, 
                                 levels = c('low', 'med', 'high'))),
                 outlier.shape = NA, width=.5) +
    geom_dotplot(aes(fill=factor(ctdna_lev, 
                                 levels = c('low', 'med', 'high'))),
                 show.legend = F,
                 binaxis='y', stackdir='center', dotsize=.7) +
    scale_color_manual(values = list_color) +
    scale_fill_manual(values = list_color) +
    stat_compare_means(comparisons = list(c('high', 'med'),
                                          c('low', 'med'),
                                          c('high', 'low')),
                       label = '..p..', size=5,
                       label.y = c(13,11,15), 
                       tip.length = .01) +
    labs(x = '', y='Mutation Count per Megabase', fill='ctDNA Level') +
    theme_pubr() +
    theme(axis.title = element_text(size=15), 
          axis.text.y = element_text(size=13),
          axis.text.x = element_blank(),
          legend.position = 'none',
          legend.text = element_text(size=10),
          legend.title = element_text(size=12),
          strip.text.x = element_text(size = 15))
}

##### fig3c: mutation types
{
  fig3c <- mut_sigs |>
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
          axis.text.y = element_text(size=13),
          axis.text.x = element_blank(),
          legend.position = c(.85,.75),
          legend.text = element_text(size=10),
          legend.title = element_text(size=12),
          strip.text.x = element_text(size = 15))
}

##### fig3d: mutational signatures
{  
  fig3d <- signatures_df |>
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
    labs(x='', y ='Mutation Count per Megabase') +
    theme_pubr() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size=12), 
          axis.text.y = element_text(size=10),
          axis.text.x = element_blank(),
          legend.position = c(.85,.75),
          legend.text = element_text(size=10),
          legend.title = element_text(size=12),
          strip.text.x = element_text(size = 15)) 
}

##### fig3e: mutational signatures
{
  stat.test <- xchisq.test(ctdna_lev ~id5, data=signatures_df)
  fig3e <- signatures_df |>
    ggplot(aes(x = factor(ctdna_lev, 
                          levels = c('low', 'med', 'high')))) +
    geom_bar(aes(fill = factor(id5, level=c('zero', 'non-zero')),
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
          axis.text.y = element_text(size=13),
          axis.text.x = element_blank(),
          legend.position = 'right',
          legend.text = element_text(size=10),
          legend.title = element_text(size=12),
          strip.text.x = element_text(size = 15))
}

################ done ##################
