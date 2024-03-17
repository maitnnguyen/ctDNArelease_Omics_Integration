############## GSVA scores for hallmark geneset
# I. Sample level
# 1.1. EOC Component
{
  samples <- rna_cohort |>
    dplyr::slice_max(EOC, by = c(patient, tissue_gr)) |>
    dplyr::arrange(desc(TF))
  
  exp_mat <- rna_z_prim#[rowSums(rna_z_prim)>0,] 
  exp_mat <- exp_mat |>
    as.data.frame() |>
    tibble::rownames_to_column(var='geneID') |>
    left_join(geneID2name |>
                dplyr::select(geneID, geneName)) |>
    dplyr::select(-geneID) |>
    filter(!(is.na(geneName)),
           !duplicated(geneName)) |>
    as.data.frame() |>
    tibble::column_to_rownames(var='geneName') %>%
    dplyr::select(colnames(.)[grepl('EOC', colnames(.))]) %>%
    dplyr::rename_all(.funs = function(x) sub('.EOC','',x))
  
  dim(exp_mat)
  exp_mat[1:2,1:4]
  EOC.Z <- exp_mat[, samples$sample]
  ncol(EOC.Z) == nrow(samples)  
  
  G.EOC <- as.matrix(rna_g[, samples$sample])
  W <- as.matrix(samples %>% 
                   tibble::column_to_rownames(var='sample') |>
                   dplyr::select(EOC, Fibroblast, Immune, Unknown))
  W.EOC <- W[ which(W[,1] > 0.25) , grep('EOC', colnames(W)) ]
  
  # normalisation count data with size factor
  # this is used for heatmap
  Z.n <- diag.mul.r( as.matrix(EOC.Z), 1./ (G.EOC * W.EOC )) 
  max(Z.n)
  Z.n[1:3,1:7]
  count_mtx <- Z.n[rowSums(Z.n)>0,] #countdata#[rowSum > 10,] |>
  
  data_mtx <- as.matrix(count_mtx)
  
  hsi_NES <- GSVA::gsva(as.matrix(data_mtx),
                        pathways.hallmark,
                        method="gsva",
                        min.sz=20,
                        max.sz=1000, 
                        kcdf="Poisson")
}

# visualisation config
{
  ann.df <- colnames(count_mtx) |>
    as.data.frame() |>
    dplyr::rename_with(~'sample') |>
    left_join(samples) |>
    left_join(ctdna_cohort %>%
                dplyr::select(patient, treatment, ctdna_lev)) |>
    dplyr::arrange(TF) |>
    tibble::column_to_rownames(var='sample') |>
    dplyr::select(ctDNA=TF, ctdna_lev, EOC, tissue_gr, treatment, Stage)
  
  nb.cols <- 5
  stage <- colorRampPalette(RColorBrewer::brewer.pal(4, "Blues"))(4)#[2:5]
  tis <- colorRampPalette(RColorBrewer::brewer.pal(4, "Dark2"))(4)#[2:5]
  ann.col=list(Stage=c('I'=stage[1],'II'=stage[2], 'III'=stage[3],'IV'=stage[4]),
               ctdna_lev = c('neg'='maroon', 'low'='darkgreen', 'high'='darkblue'),
               treatment=c('NACT'='#794044','PDS'='#088da5'),
               tissue_gr=c('Ome'=tis[1], 'Per'=tis[2], 'Ova_Tub'=tis[3], 'Other'=tis[4]))
  
}

{
  hsi_sample <- t(hsi_NES) %>%
    as.data.frame() |>
    tibble::rownames_to_column(var='sample') %>%
    mutate(patient = sub('_.*', '', sample)) |>
    left_join(samples) 
  
  mean_score_pathway <- hsi_sample %>% 
    tidyr::gather(key = 'pathway', value = 'score', 2:51) |>
    mutate(pathway = sub('HALLMARK_', '', pathway)) |>
    dplyr::group_by(patient, pheno, pathway) |>
    dplyr::summarise(Score = mean(score))
  
  # 4 hallmarks proliferation-related gene sets
  proliferate_pathways <- c( 'E2F_TARGETS', 'G2M_CHECKPOINT',
                            'MYC_TARGETS_V1','MYC_TARGETS_V2')
  proliferation_comp <- mean_score_pathway |>
    filter(pathway %in% proliferate_pathways) |>
    ggplot(aes(x=factor(pheno, 
                        levels=c('low', 'med', 'high')), 
               y = Score)) + 
    geom_violin(aes(fill = pheno),
                alpha=.4,
                show.legend = F)+ 
    scale_fill_manual(values = c("low"=list_color[1],
                                 "med" = list_color[2],
                                 "high" = list_color[3])) +
    geom_boxplot(aes(fill = pheno),
                 width=.2, show.legend = F)+
    labs(x = '', y = 'pathway score') +
    stat_compare_means(comparisons = list(c('high','low'),
                                          c('med', 'low'),
                                          c('high', 'med')),
                       label = 'p.signif', size=2.5) +
    facet_wrap(~factor(pathway, levels = proliferate_pathways), ncol = 4) +
    theme_pubr() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size=15), 
          axis.text.y = element_text(size=13),
          axis.text.x = element_blank(),
          legend.position = c(.85,.75),
          legend.text = element_text(size=10),
          legend.title = element_text(size=15),
          strip.text.x = element_text(size = 7))
  proliferation_comp
  svg('results/figures/final/Supp/fig5a.svg', width = 5, height = 3)
  proliferation_comp
  dev.off()
  
  # metabolism-related pathway
  metabolic_pathways <- c('OXIDATIVE_PHOSPHORYLATION', 'BILE_ACID_METABOLISM',
                          #'CHOLESTEROL_HOMEOSTASIS',
                          'HEME_METABOLISM',
                          #'XENOBIOTIC_METABOLISM',
                          'FATTY_ACID_METABOLISM')
  metabolism_comp <- mean_score_pathway |>
    filter(pathway %in% metabolic_pathways) |>
    ggplot(aes(x=factor(pheno, 
                        levels=c('low', 'med', 'high')), 
               y = Score)) + 
    geom_violin(aes(fill = pheno),
                alpha=.4,
                show.legend = F)+ 
    scale_fill_manual(values = c("low"=list_color[1],
                                 "med" = list_color[2],
                                 "high" = list_color[3])) +
    geom_boxplot(aes(fill = pheno),
                 width=.2, show.legend = F)+
    labs(x = '', y = 'pathway score') +
    stat_compare_means(comparisons = list(c('high','low'),
                                          c('med', 'low'),
                                          c('high', 'med')),
                       label = 'p.signif', size=2.5) +
    facet_wrap(~factor(pathway, levels = metabolic_pathways), ncol = 4) +
    theme_pubr() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size=15), 
          axis.text.y = element_text(size=13),
          axis.text.x = element_blank(),
          legend.position = c(.85,.75),
          legend.text = element_text(size=10),
          legend.title = element_text(size=15),
          strip.text.x = element_text(size = 7))
  metabolism_comp
  svg('results/figures/final/Supp/fig5b.svg', width = 5, height = 3)
  metabolism_comp
  dev.off()
  
  # immune-related pathway
  imm_pathways <- c('IL6_JAK_STAT3_SIGNALING',
                    'ALLOGRAFT_REJECTION', 
                    'COMPLEMENT',#'COAGULATION',
                    #'INTERFERON_ALPHA_RESPONSE','INTERFERON_GAMMA_RESPONSE', 
                    'INFLAMMATORY_RESPONSE')
  imm_comp <- mean_score_pathway |>
    filter(pathway %in% imm_pathways) |>
    ggplot(aes(x=factor(pheno, 
                        levels=c('low', 'med', 'high')), 
               y = Score)) + 
    geom_violin(aes(fill = pheno),
                alpha=.4,
                show.legend = F)+ 
    scale_fill_manual(values = c("low"=list_color[1],
                                 "med" = list_color[2],
                                 "high" = list_color[3])) +
    geom_boxplot(aes(fill = pheno),
                 width=.2, show.legend = F)+
    labs(x = '', y = 'pathway score') +
    stat_compare_means(comparisons = list(c('high','low'),
                                          c('med', 'low'),
                                          c('high', 'med')),
                       label = 'p.signif', size=2.5) +
    facet_wrap(~factor(pathway, levels = imm_pathways), ncol = 4) +
    theme_pubr() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size=15), 
          axis.text.y = element_text(size=13),
          axis.text.x = element_blank(),
          legend.position = c(.85,.75),
          legend.text = element_text(size=10),
          legend.title = element_text(size=15),
          strip.text.x = element_text(size = 7))
  imm_comp
  svg('results/figures/final/Supp/fig5c.svg', width = 5, height = 3)
  imm_comp
  dev.off()
  # signaling-related pathway
  signaling_pathways <- c('MTORC1_SIGNALING',
                          'HEDGEHOG_SIGNALING',#'PI3K_AKT_MTOR_SIGNALING',
                          'TNFA_SIGNALING_VIA_NFKB','TGF_BETA_SIGNALING')
  signal_comp <- mean_score_pathway |>
    filter(pathway %in% signaling_pathways) |>
    ggplot(aes(x=factor(pheno, 
                        levels=c('low', 'med', 'high')), 
               y = Score)) + 
    geom_violin(aes(fill = pheno),
                alpha=.4,
                show.legend = F)+ 
    scale_fill_manual(values = c("low"=list_color[1],
                                 "med" = list_color[2],
                                 "high" = list_color[3])) +
    geom_boxplot(aes(fill = pheno),
                 width=.2, show.legend = F)+
    labs(x = '', y = 'pathway score') +
    stat_compare_means(comparisons = list(c('high','low'),
                                          c('med', 'low'),
                                          c('high', 'med')),
                       label = 'p.signif', size=2.5) +
    facet_wrap(~factor(pathway, levels = signaling_pathways), ncol = 4) +
    theme_pubr() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size=15), 
          axis.text.y = element_text(size=13),
          axis.text.x = element_blank(),
          legend.position = c(.85,.75),
          legend.text = element_text(size=10),
          legend.title = element_text(size=15),
          strip.text.x = element_text(size = 7))
  signal_comp
  svg('results/figures/final/Supp/fig5d.svg', width = 5, height = 3)
  signal_comp
  dev.off()
  
  
  pathways <- unique(mean_score_pathway$pathway) |>
    as.data.frame() |>
    dplyr::rename_with(~'pathway') |>
    mutate(Kruskal_Wallis_stat = NA,
           p_val = NA)
  
  for (i in pathways$pathway){ 
    dat <- mean_score_pathway |>
      filter(pathway == i)
    
    comp_res <- kruskal.test(Score ~ pheno, data = dat)
    
    pathways[pathways$pathway==i,2:3] = c(comp_res$statistic[[1]], comp_res$p.value)
  }
  
  sig.pathways <- pathways[pathways$p_val<.05,]$pathway
  # heatmap for pathway scores
  #tt <- hsi_sample[,c('sample', 
  #                    hallmark.resTidy[hallmark.resTidy$padj<.01,]$pathway)]
  interested_pathways <- hsi_sample[,c('patient',
                      colnames(hsi_sample)[grepl('HALLMARK',colnames(hsi_sample))])]
  #tt <- hsi_sample[,c('patient',
  #                    colnames(hsi_sample)[sub('HALLMARK_', '',colnames(hsi_sample)) %in%
  #                                           sig.pathways])]
  names(interested_pathways) <- sub('HALLMARK_', '', names(interested_pathways))
  interested_pathways <- interested_pathways |>
    dplyr::group_by(patient) |>
    dplyr::summarise_all(.funs = mean) |>
    left_join(ctdna_cohort |>
                dplyr::select(patient, TF)) |>
    dplyr::arrange(TF) |>
    tibble::column_to_rownames(var='patient') |>
    dplyr::select(-TF)
  interested_pathways <- interested_pathways[,c(proliferate_pathways, metabolic_pathways, imm_pathways, signaling_pathways)]
  ann.row <- interested_pathways |>
    tibble::rownames_to_column(var='patient') |>
    dplyr::select(patient) |>
    left_join(ctdna_cohort |>
                dplyr::select(patient, ctdna_tf=TF, ctdna_lev)) |>
    dplyr::arrange(ctdna_tf) |>
    tibble::column_to_rownames(var='patient') 
  ann.col=list(ctdna_lev = c('low'=list_color[1],
                             'med'=list_color[2], 
                             'high'=list_color[3]),
               treatment=c('NACT'='#794044','PDS'='#088da5'))
  breaksList = seq(-3, 3, by = .4)
  
  svg('results/figures/final/Main/fig5a.svg', width = 10, height = 6)
  pheatmap::pheatmap(interested_pathways, annotation_row = ann.row, 
                     annotation_colors = ann.col, 
                     cluster_rows = F, show_rownames = F, 
                     scale = 'column',
                     color = colorRampPalette(c("#3399ff", "white", "#ff6666"))(length(breaksList)),
                     breaks = breaksList, gaps_row = c(cumsum(c(48,31,24)),cumsum(c(48,31,24))), border_color = 'black',
                     fontsize_col = 8, cutree_rows = 4, cutree_cols = 2,
                    height = 3, width = 6)
  dev.off()
  
  order_names <- c('HEDGEHOG_SIGNALING','TGF_BETA_SIGNALING',
                   'IL6_JAK_STAT3_SIGNALING', 'INFLAMMATORY_RESPONSE',
                   'TNFA_SIGNALING_VIA_NFKB', 'ALLOGRAFT_REJECTION',
                   'HEME_METABOLISM', 'COMPLEMENT',
                   'BILE_ACID_METABOLISM', 'MTORC1_SIGNALING', 
                   'OXIDATIVE_PHOSPHORYLATION','FATTY_ACID_METABOLISM',
                   'E2F_TARGETS', 'G2M_CHECKPOINT', 'MYC_TARGETS_V1',
                   'MYC_TARGETS_V2')
  svg('results/figures/final/Main/fig5b.svg', width = 10, height = 10)
  corrplot::corrplot(cor(interested_pathways), type = 'upper', 
                   tl.cex = .6, tl.col = c(rep('black',4),
                                           rep('#660066',4),
                                           rep('#000080',4),
                                           rep('#794044',4)))
  #corrplot::corrplot(cor(tt), type = 'upper', 
  #                   tl.cex = .8, tl.col = 'black')
  dev.off()
  # pca
  {
    data <- t(hsi_NES) %>%
      as.data.frame()
    colnames(data) = sub('HALLMARK_','',colnames(data))
    pathway2check <- sub('HALLMARK_','',
                         hallmark.resTidy$pathway[hallmark.resTidy$padj<.05])
    data <- data %>%  dplyr::select(colnames(.)[colnames(.) %in%  #pathway2check]) |>
                                                  c(proliferate_pathways, 
                                                    metabolic_pathways, 
                                                    imm_pathways, 
                                                    signaling_pathways)]) |>
      tibble::rownames_to_column(var='sample') |>
      mutate(patient = sub('_.*', '', sample)) |>
      dplyr::select(-sample) |>
      dplyr::group_by(patient) |>
      dplyr::summarise_all(.funs = median) |>
      as.data.frame() |>
      left_join(ctdna_cohort %>%
                  dplyr::select(patient, TF, ctdna_lev)) |>
      tibble::column_to_rownames(var = 'patient') 
    dim(data)
    
    res.pca <- prcomp(data[, 1:(ncol(data)-2)], scale=F)
    # Visualize eigenvalues/variances
    factoextra::fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 70))
    # Extract the results for variables
    var <- factoextra::get_pca_var(res.pca)
    # Contributions of variables to PC1
    factoextra::fviz_contrib(res.pca, choice = "var", axes = 1, top = 10, )
    # Contributions of variables to PC2
    factoextra::fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
    #fviz_contrib(res.pca, choice = "var", axes = 3, top = 10)
    #fviz_contrib(res.pca, choice = "var", axes = 4, top = 10)
    
    # Control variable colors using their contributions to the principle axis
    svg('results/figures/final/Main/fig5c.svg')
    factoextra::fviz_pca_var(res.pca, col.var="contrib", axes = c(1, 2),
                             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                             repel = TRUE, # Avoid text overlapping
                             labelsize = 4
    ) + theme_pubr() +
      labs(color = 'Contribution', x = 'PCA1 (43.9%)', y = 'PCA2 (17.6%)') +
      theme(title = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_text(size=12, family = 'Arial'), 
            #axis.text.y = element_blank(),
            axis.text.x = element_text(size=10, family = 'Arial'),
            legend.position = c(.2,.85),
            legend.text = element_text(size=10, family = 'Arial'),
            legend.title = element_text(size=15, family = 'Arial'),
            strip.text.x = element_text(size = 15, family = 'Arial'))
    dev.off()
    
    factoextra::fviz_pca_biplot(res.pca, label = "contrib", #axes = c(1, 3),
                                habillage=factor(data$ctdna_lev, 
                                                 levels = c('low', 'med', 'high')),
                                addEllipses=TRUE, ellipse.level=0.3,
                                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                xlim=c(-1.25,1.25), ylim=c(-1,1),
                                palette = list_color, legend=NULL
                                #ggtheme = theme_minimal()
    ) 
    
    
    svg('results/figures/final/Main/fig5d.svg')
    factoextra::fviz_pca_ind(res.pca, label="contrib",#axes = c(1, 3), 
                             habillage=factor(data$ctdna_lev, 
                                              levels = c('low', 'med', 'high')),
                             #addEllipses=F, 
                             alpha=.9,
                             #addEllipses=TRUE, ellipse.level=0.7,
                             pointsize = 2, 
                             invisible="quali",
                             palette = list_color) +
      labs(x = 'PCA1 (43.9%)', y = 'PCA2 (17.6%)') +
      #geom_point(size=3)
      theme(title = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_text(size=12, family = 'Arial'), 
            axis.text.y = element_text(size=10, family = 'Arial'),
            axis.text.x = element_text(size=10, family = 'Arial'),
            legend.position = 'none',
            strip.text.x = element_text(size = 15, family = 'Arial'))
    dev.off()
  }
  }
  
# export cluster of 4 and run surv
{
  
}


  hsi_sample %>% 
    ggplot(aes(x=ctdna_lev, y = HALLMARK_HYPOXIA)) +
    geom_violin()+ 
    geom_boxplot(width=.2)+
    stat_compare_means(comparisons = list(c('high','low'),
                                          c('neg', 'low'),
                                          c('high', 'neg')),
                       label = '..p.adj..')
  
  hsi_sample %>% ggplot(aes(x=ctdna_lev, y = HALLMARK_HYPOXIA)) + 
    geom_violin()+ 
    geom_boxplot(width=.2)+
    stat_compare_means(comparisons = list(c('high','low'),
                                          c('neg', 'low'),
                                          c('high', 'neg')),
                       label = '..p.adj..')
  
  hsi_sample %>% ggplot(aes(x=ctdna_lev, y = HALLMARK_SPERMATOGENESIS)) + 
    geom_violin()+ 
    geom_boxplot(width=.2)+
    stat_compare_means(comparisons = list(c('high','low'),
                                          c('neg', 'low'),
                                          c('high', 'neg')),
                       label = '..p.adj..') 
  
  hsi_sample %>% ggplot(aes(x=ctdna_lev, y = HALLMARK_MYC_TARGETS_V1)) + 
    geom_violin(aes(fill = ctdna_lev), alpha=.7, show.legend = F)+ 
    geom_boxplot(fill='white', width=.2)+
    stat_compare_means(comparisons = list(c('high','low'),
                                          c('neg', 'low'),
                                          c('high', 'neg')),
                       label = '..p.adj..') +
    scale_fill_manual(values = c('darkblue', 'darkgreen', 'maroon')) +
    labs(x = '') +
    theme()
  
  hsi_sample %>% ggplot(aes(x=ctdna_lev, y = HALLMARK_G2M_CHECKPOINT)) + 
    geom_violin(aes(fill = ctdna_lev), alpha=.7, show.legend = F)+ 
    geom_boxplot(fill='white', width=.2)+
    stat_compare_means(comparisons = list(c('high','low'),
                                          c('neg', 'low'),
                                          c('high', 'neg')),
                       label = '..p.adj..') +
    scale_fill_manual(values = c('darkblue', 'darkgreen', 'maroon')) +
    labs(x = '') +
    theme()
  
  hsi_sample %>% ggplot(aes(x=ctdna_lev, y = HALLMARK_TNFA_SIGNALING_VIA_NFKB)) + 
    geom_violin(aes(fill = ctdna_lev), alpha=.7, show.legend = F)+ 
    geom_boxplot(fill='white', width=.2)+
    stat_compare_means(comparisons = list(c('high','low'),
                                          c('neg', 'low'),
                                          c('high', 'neg')),
                       label = '..p.adj..') +
    scale_fill_manual(values = c('darkblue', 'darkgreen', 'maroon')) +
    labs(x = '') +
    theme()
  
  hsi_sample %>% ggplot(aes(x=ctdna_lev, y = HALLMARK_TGF_BETA_SIGNALING)) + 
    geom_violin(aes(fill = ctdna_lev), alpha=.7, show.legend = F)+ 
    geom_boxplot(fill='white', width=.2)+
    stat_compare_means(comparisons = list(c('high','low'),
                                          c('neg', 'low'),
                                          c('high', 'neg')),
                       label = '..p.adj..') +
    scale_fill_manual(values = c('darkblue', 'darkgreen', 'maroon')) +
    labs(x = '') +
    theme()
  
  hsi_sample %>% ggplot(aes(x=ctdna_lev, y = HALLMARK_IL6_JAK_STAT3_SIGNALING)) + 
    geom_violin(aes(fill = ctdna_lev), alpha=.7, show.legend = F)+ 
    geom_boxplot(fill='white', width=.2)+
    stat_compare_means(comparisons = list(c('high','low'),
                                          c('neg', 'low'),
                                          c('high', 'neg')),
                       label = '..p.adj..') +
    scale_fill_manual(values = c('darkblue', 'darkgreen', 'maroon')) +
    labs(x = '') +
    theme()
  
  hsi_sample %>% ggplot(aes(x=ctdna_lev, y = HALLMARK_TNFA_SIGNALING_VIA_NFKB)) + 
    geom_violin(aes(fill = ctdna_lev), alpha=.7, show.legend = F)+ 
    geom_boxplot(fill='white', width=.2)+
    stat_compare_means(comparisons = list(c('high','low'),
                                          c('neg', 'low'),
                                          c('high', 'neg')),
                       label = '..p..') +
    scale_fill_manual(values = c('darkblue', 'darkgreen', 'maroon')) +
    labs(x = '') +
    theme()
  
  hsi_sample %>% ggplot(aes(x=ctdna_lev, y = HALLMARK_IL6_JAK_STAT3_SIGNALING)) + 
    geom_violin(aes(fill = ctdna_lev), alpha=.7, show.legend = F)+ 
    geom_boxplot(fill='white', width=.2)+
    stat_compare_means(comparisons = list(c('high','low'),
                                          c('neg', 'low'),
                                          c('high', 'neg')),
                       label = '..p.adj..') +
    scale_fill_manual(values = c('darkblue', 'darkgreen', 'maroon')) +
    labs(x = '') +
    theme()
  
  hsi_sample %>% ggplot(aes(x=ctdna_lev, y = HALLMARK_TGF_BETA_SIGNALING)) + 
    geom_violin(aes(fill = ctdna_lev), alpha=.7, show.legend = F)+ 
    geom_boxplot(fill='white', width=.2)+
    stat_compare_means(comparisons = list(c('high','low'),
                                          c('neg', 'low'),
                                          c('high', 'neg')),
                       label = '..p.adj..') +
    scale_fill_manual(values = c('darkblue', 'darkgreen', 'maroon')) +
    labs(x = '') +
    theme()
  
  hsi_sample %>% ggplot(aes(x=ctdna_lev, y = HALLMARK_E2F_TARGETS)) + 
    geom_violin(aes(fill = ctdna_lev), alpha=.7, show.legend = F)+ 
    geom_boxplot(fill='white', width=.2)+
    stat_compare_means(comparisons = list(c('high','low'),
                                          c('neg', 'low'),
                                          c('high', 'neg')),
                       label = '..p.adj..') +
    scale_fill_manual(values = c('darkblue', 'darkgreen', 'maroon')) +
    labs(x = '') +
    theme()
  
  hsi_sample %>% ggplot(aes(x=ctdna_lev, y = HALLMARK_KRAS_SIGNALING_DN)) + 
    geom_violin(aes(fill = ctdna_lev), alpha=.7, show.legend = F)+ 
    geom_boxplot(fill='white', width=.2)+
    stat_compare_means(comparisons = list(c('high','low'),
                                          c('neg', 'low'),
                                          c('high', 'neg')),
                       label = '..p.adj..') +
    scale_fill_manual(values = c('darkblue', 'darkgreen', 'maroon')) +
    labs(x = '') +
    theme()
  
  hsi_sample %>% ggplot(aes(x=ctdna_lev, y = HALLMARK_INTERFERON_GAMMA_RESPONSE)) + 
    geom_violin()+
    geom_boxplot(width=.2)+ 
    stat_compare_means(comparisons = list(c('high','low'),
                                          c('neg', 'low'),
                                          c('high', 'neg')),
                       label = '..p.adj..') +
    facet_wrap(~tissue_gr)
  
  hsi_sample %>% ggplot(aes(x=ctdna_lev, y = HALLMARK_FATTY_ACID_METABOLISM)) + 
    geom_violin()+
    geom_boxplot(width=.2)+ 
    stat_compare_means(comparisons = list(c('high','low'),
                                          c('neg', 'low'),
                                          c('high', 'neg')),
                       label = '..p.adj..') +
    facet_wrap(~tissue_gr)
  
  hsi_sample %>% ggplot(aes(x=ctdna_lev, y = HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY)) + 
    geom_violin()+ 
    geom_boxplot(width=.1)+
    stat_compare_means(comparisons = list(c('high','low'),
                                          c('neg', 'low'),
                                          c('high', 'neg')),
                       label = '..p.adj..')  +
    facet_wrap(~tissue_gr)
  
  hsi_sample %>% 
    dplyr::group_by(patient, ctdna_lev) |>
    dplyr::summarise(score = median(HALLMARK_SPERMATOGENESIS)) |>
    ggplot(aes(x=ctdna_lev, y = score)) + 
    #ggplot(aes(x=ctdna_lev, y = HALLMARK_PEROXISOME)) + 
    geom_violin()+ 
    geom_boxplot(width=.1)+
    stat_compare_means(comparisons = list(c('high','low'),
                                          c('neg', 'low'),
                                          c('high', 'neg')),
                       label = '..p.adj..') 
  
  
  hsi_sample %>% ggplot(aes(x=ctdna_lev, y = HALLMARK_E2F_TARGETS)) + 
    geom_violin()+ 
    geom_boxplot(width=.1)+
    stat_compare_means(comparisons = list(c('high','low'),
                                          c('neg', 'low'),
                                          c('high', 'neg')),
                       label = '..p.adj..')  +
    facet_wrap(~tissue_gr)
  
  
  hsi_sample %>% #ggplot(aes(x=TF, y = HALLMARK_OXIDATIVE_PHOSPHORYLATION)) +
    #geom_point() +
    ggscatter(x='TF', y='HALLMARK_SPERMATOGENESIS', cor.coef = T, add='reg.line') +
    facet_grid(HRD ~.)
    
  hsi_sample %>% 
    filter(TF>thresh1) |>
    ggscatter(x='TF', y = 'HALLMARK_SPERMATOGENESIS', cor.coef = T,
                   add = 'reg.line', facet.by = 'tissue_gr')
  
  lm1 <- glm(TF ~ (HALLMARK_DNA_REPAIR+HALLMARK_E2F_TARGETS+
              HALLMARK_OXIDATIVE_PHOSPHORYLATION)*HRD, data=tt)
  summary(lm1)
  tt <- tt %>%
    mutate(y = ifelse(pheno=='neg', 1, 0),
           HRD=ifelse(is.na(HRD), 'ND', HRD))
  model <- glm( y ~ (HALLMARK_DNA_REPAIR + HALLMARK_E2F_TARGETS +
                  HALLMARK_OXIDATIVE_PHOSPHORYLATION +
                  HALLMARK_FATTY_ACID_METABOLISM ) * HRD , 
                data = tt, family = binomial)
  summary(model)
  
  hsi_pat <- hsi_sample |>
    dplyr::group_by(patient, TF, ctdna_lev, HRD) |>
    dplyr::summarise(score = mean(HALLMARK_GLYCOLYSIS)) |>
    as.data.frame() |>
    dplyr::arrange(TF)
  pat <- factor(hsi_pat$patient, levels = hsi_pat$patient)
  hsi_pat |> 
    ggplot(aes(x = ctdna_lev, 
               y = score)) + 
    geom_violin() +
    geom_boxplot(width = .1) +
    stat_compare_means(comparisons = list(c('high', 'low'),
                                          c('high', 'neg'),
                                          c('low', 'neg'))) +
    facet_wrap(~HRD)
    
    ggplot(aes(x = factor(patient, levels = pat), 
               y = score)) + 
    geom_col(aes(fill=ctdna_lev))+ 
    
  #  scale_fill_continuous() +
    theme(axis.text.x = element_text(angle = 45))
  # PCA
  hsi_NES
  #pat <- samples[samples$pheno!='neg',]$patient
  


# 1.1.1 EOC Component - each patient with 1 highestt EOC sample
{
  samples <- rna_cohort |>
    mutate(pheno = ifelse(TF < thresh, 'low', 'high'),
           tissue_gr = ifelse(tissue %in% c('Ome', 'Per'), tissue, 
                              ifelse(tissue %in% c('Ova', 'Tub', 'Adn', 'And'),
                                     'Ova_Tub','Other'))) |>
    dplyr::slice_max(EOC, by = c(patient)) |>
    dplyr::arrange(desc(TF))

  dim(samples)
  
  EOC.Z <- rna_z_prim[,grep('EOC', colnames(rna_z_prim))]
  colnames(EOC.Z) = sub('.EOC', '', colnames(EOC.Z))
  EOC.Z <- EOC.Z[, samples$sample]
  ncol(EOC.Z) == nrow(samples)  
  
  G.EOC <- as.matrix(rna_g_prim[, samples$sample])
  W <- as.matrix(samples %>% 
                   tibble::column_to_rownames(var='sample') |>
                   dplyr::select(EOC, Fibroblast, Immune, Unknown))
  W.EOC <- W[ which(W[,1] >0.25) , grep('EOC', colnames(W)) ]
  
  Z.n <- diag.mul.r( as.matrix(EOC.Z), 1./ (G.EOC * W.EOC )) 
  max(Z.n)
  
  # count matrix for EOC component
  countdata <- EOC.Z
  #rowSum <- apply(countdata, 1, sum)
  
  count_mtx <- countdata#[rowSum > 10,] |>
  
  ann.df <- colnames(count_mtx) |>
    as.data.frame() |>
    dplyr::rename_with(~'sample') |>
    left_join(samples) |>
    dplyr::arrange(TF) |>
    tibble::column_to_rownames(var='sample') |>
    dplyr::select(ctDNA=TF, EOC, tissue_gr, treatment, Stage)
  
  stage <- colorRampPalette(RColorBrewer::brewer.pal(4, "Blues"))(4)#[2:5]
  tis <- colorRampPalette(RColorBrewer::brewer.pal(4, "Dark2"))(4)#[2:5]
  ann.col=list(Stage=c('I'=stage[1],'II'=stage[2], 'III'=stage[3],'IV'=stage[4]),
               treatment=c('NACT'='#794044','PDS'='#088da5'),
               tissue_gr=c('Ome'=tis[1], 'Per'=tis[2], 'Ova_Tub'=tis[3], 'Other'=tis[4]))
  data_mtx <- as.matrix(count_mtx[,rownames(ann.df)])
  
  hsi_NES <- GSVA::gsva(as.matrix(data_mtx),
                        pathways.hallmark,
                        method="ssgsea",
                        min.sz=10,
                        max.sz=1000)
  svg('results/figures/supplementary/3_2hallmark_pathway_score_max_eoc_pat_lev.svg', width = 15, height = 15)
  pheatmap::pheatmap(hsi_NES, scale = 'row',
                     annotation_colors = ann.col, 
                     annotation_col = ann.df, cluster_cols = F, 
                     color = colorRampPalette(c("navy", "white", "red"))(70),
                     show_colnames = F, fontsize_row = 10, cutree_rows = 10)
  dev.off()
}

# 1.1.2 EOC Component - each patient with 1 highest EOC Ova_Tub sample
{
  samples <- rna_cohort |>
    dplyr::slice_max(EOC, by = c(patient)) |>
    dplyr::arrange(desc(TF))
  
  dim(samples)
  
  EOC.Z <- rna_z_prim[,grep('EOC', colnames(rna_z_prim))]
  colnames(EOC.Z) = sub('.EOC', '', colnames(EOC.Z))
  EOC.Z <- EOC.Z[, samples$sample]
  ncol(EOC.Z) == nrow(samples)  
  
  G.EOC <- as.matrix(rna_g_prim[, samples$sample])
  W <- as.matrix(samples %>% 
                   tibble::column_to_rownames(var='sample') |>
                   dplyr::select(EOC, Fibroblast, Immune, Unknown))
  W.EOC <- W[ which(W[,1] >0.25) , grep('EOC', colnames(W)) ]
  
  Z.n <- diag.mul.r( as.matrix(EOC.Z), 1./ (G.EOC * W.EOC )) 
  max(Z.n)
  
  # count matrix for EOC component
  countdata <- Z.n |>
    as.data.frame()
  countdata$symbol = geneID2name$uniqName[match(sub( '[.].*$', '', rownames(countdata)), 
                                                sub( '[.].*$', '', geneID2name$geneID ))]
  rownames(countdata) <- NULL
  rownames(countdata) <- countdata$symbol
  
  rowSum <- apply(countdata, 1, sum)
  
  count_mtx <- as.matrix(countdata[rowSum>0,] |> dplyr::select(-symbol))
    
  #count_mtx <- counts(dds)
  ann.df <- colnames(count_mtx) |>
    as.data.frame() |>
    dplyr::rename_with(~'sample') |>
    left_join(samples) |>
    dplyr::arrange(TF) |>
    tibble::column_to_rownames(var='sample') |>
    dplyr::select(ctDNA=TF, EOC, tis_gr)
  
  stage <- colorRampPalette(RColorBrewer::brewer.pal(4, "Blues"))(4)#[2:5]
  tis <- colorRampPalette(RColorBrewer::brewer.pal(4, "Dark2"))(4)#[2:5]
  ann.col=list(Stage=c('I'=stage[1],'II'=stage[2], 'III'=stage[3],'IV'=stage[4]),
               treatment=c('NACT'='#794044','PDS'='#088da5'),
               tissue_gr=c('Ome'=tis[1], 'Per'=tis[2], 'Ova_Tub'=tis[3], 'Other'=tis[4]))
  data_mtx <- as.matrix(count_mtx[,rownames(ann.df)])
  
  hsi_NES <- GSVA::gsva(as.matrix(count_mtx),
                        pathways.hallmark,
                        method="ssgsea",
                        min.sz=10,
                        max.sz=1000)
  hsi_NES_scale <- hsi_NES/apply(hsi_NES, 1, median)
  
  hall_mark_score <- hsi_NES_scale |>
    as.data.frame() |>
    tibble::rownames_to_column(var='hallmark') |>
    tidyr::gather(key='sample', value='score', 2:98) |>
    left_join(samples)
  
  hall_mark_score |>
    filter(tis_gr=='Ome',
           hallmark=='HALLMARK_E2F_TARGETS') |>
    ggplot(aes(x=ctdna_detection, y=score)) +
    geom_violin() +
    geom_boxplot(width=.1) +
    stat_compare_means()
  
  
  
  svg('results/figures/supplementary/3_2hallmark_pathway_score_max_eoc_Ova_lev.svg', width = 15, height = 15)
  pheatmap::pheatmap(hsi_NES, scale = 'row',
                     annotation_colors = ann.col, 
                     annotation_col = ann.df, cluster_cols = F, 
                     color = colorRampPalette(c("navy", "white", "red"))(70),
                     show_colnames = F, fontsize_row = 10, cutree_rows = 10) +
    scale_fill_continuous(limits = c(0,50), breaks = c(0, 5, 10, 15, 20))
                          
  dev.off()
}

# 1.1.3 EOC Component - each patient with 1 highest EOC Ome sample
{
  samples <- rna_cohort |>
    mutate(pheno = ifelse(TF < thresh, 'low', 'high'),
           tissue_gr = ifelse(tissue %in% c('Ome', 'Per'), tissue, 
                              ifelse(tissue %in% c('Ova', 'Tub', 'Adn', 'And'),
                                     'Ova_Tub','Other'))) |>
    filter(tissue_gr=='Ome') |>
    dplyr::slice_max(EOC, by = c(patient)) |>
    dplyr::arrange(desc(TF))
  
  dim(samples)
  
  EOC.Z <- rna_z_prim[,grep('EOC', colnames(rna_z_prim))]
  colnames(EOC.Z) = sub('.EOC', '', colnames(EOC.Z))
  EOC.Z <- EOC.Z[, samples$sample]
  ncol(EOC.Z) == nrow(samples)  
  
  G.EOC <- as.matrix(rna_g_prim[, samples$sample])
  W <- as.matrix(samples %>% 
                   tibble::column_to_rownames(var='sample') |>
                   dplyr::select(EOC, Fibroblast, Immune, Unknown))
  W.EOC <- W[ which(W[,1] >0.25) , grep('EOC', colnames(W)) ]
  
  Z.n <- diag.mul.r( as.matrix(EOC.Z), 1./ (G.EOC * W.EOC )) 
  max(Z.n)
  
  # count matrix for EOC component
  countdata <- EOC.Z
  #rowSum <- apply(countdata, 1, sum)
  
  count_mtx <- countdata#[rowSum > 10,] |>
  
  ann.df <- colnames(count_mtx) |>
    as.data.frame() |>
    dplyr::rename_with(~'sample') |>
    left_join(samples) |>
    dplyr::arrange(TF) |>
    tibble::column_to_rownames(var='sample') |>
    dplyr::select(ctDNA=TF, EOC, tissue_gr, treatment, Stage)
  
  stage <- colorRampPalette(RColorBrewer::brewer.pal(4, "Blues"))(4)#[2:5]
  tis <- colorRampPalette(RColorBrewer::brewer.pal(4, "Dark2"))(4)#[2:5]
  ann.col=list(Stage=c('I'=stage[1],'II'=stage[2], 'III'=stage[3],'IV'=stage[4]),
               treatment=c('NACT'='#794044','PDS'='#088da5'),
               tissue_gr=c('Ome'=tis[1], 'Per'=tis[2], 'Ova_Tub'=tis[3], 'Other'=tis[4]))
  data_mtx <- as.matrix(count_mtx[,rownames(ann.df)])
  
  hsi_NES <- GSVA::gsva(as.matrix(data_mtx),
                        pathways.hallmark,
                        method="ssgsea",
                        min.sz=10,
                        max.sz=1000)
  hsi_NES_scale <- hsi_NES/apply(hsi_NES, 1, median)
  svg('results/figures/supplementary/3_2hallmark_pathway_score_max_eoc_Ome_lev.svg', width = 15, height = 15)
  pheatmap::pheatmap(hsi_NES, scale = 'row',
                     annotation_colors = ann.col, 
                     annotation_col = ann.df, cluster_cols = F, 
                     color = colorRampPalette(c("navy", "white", "red"))(70),
                     show_colnames = F, fontsize_row = 10, cutree_rows = 10)
  dev.off()
}

# 1.2. Fibroblast Component
{
  samples <- rna_cohort |>
    mutate(pheno = ifelse(TF < thresh, 'low', 'high'),
           tissue_gr = ifelse(tissue %in% c('Ome', 'Per'), tissue, 
                              ifelse(tissue %in% c('Ova', 'Tub', 'Adn', 'And'),
                                     'Ova_Tub','Other'))) |>
    #dplyr::slice_max(EOC, by = c(patient, tissue_gr)) |>
    dplyr::arrange(desc(TF))
   
  EOC.Z <- rna_z_prim[,grep('Fibroblast', colnames(rna_z_prim))]
  colnames(EOC.Z) = sub('.Fibroblast', '', colnames(EOC.Z))
  EOC.Z <- EOC.Z[, samples$sample]
  ncol(EOC.Z) == nrow(samples)  
  
  G.EOC <- as.matrix(rna_g_prim[, samples$sample])
  W <- as.matrix(samples %>% 
                   tibble::column_to_rownames(var='sample') |>
                   dplyr::select(EOC, Fibroblast, Immune, Unknown))
  W.EOC <- W[ which(W[,1] >0.25) , grep('Fibroblast', colnames(W)) ]
  
  Z.n <- diag.mul.r( as.matrix(EOC.Z), 1./ (G.EOC * W.EOC )) 
  max(Z.n)
  
  # count matrix for EOC component
  countdata <- EOC.Z
  #rowSum <- apply(countdata, 1, sum)
  
  count_mtx <- countdata#[rowSum > 10,] |>
  
  ann.df <- colnames(count_mtx) |>
    as.data.frame() |>
    dplyr::rename_with(~'sample') |>
    left_join(samples) |>
    dplyr::arrange(TF) |>
    tibble::column_to_rownames(var='sample') |>
    dplyr::select(ctDNA=TF, EOC, tissue_gr, treatment, Stage)
  
  stage <- colorRampPalette(RColorBrewer::brewer.pal(4, "Blues"))(4)#[2:5]
  tis <- colorRampPalette(RColorBrewer::brewer.pal(4, "Dark2"))(4)#[2:5]
  ann.col=list(Stage=c('I'=stage[1],'II'=stage[2], 'III'=stage[3],'IV'=stage[4]),
               treatment=c('NACT'='#794044','PDS'='#088da5'),
               tissue_gr=c('Ome'=tis[1], 'Per'=tis[2], 'Ova_Tub'=tis[3], 'Other'=tis[4]))
  data_mtx <- as.matrix(count_mtx[,rownames(ann.df)])
  
  hsi_NES <- GSVA::gsva(as.matrix(data_mtx),
                        pathways.hallmark,
                        method="ssgsea",
                        min.sz=10,
                        max.sz=1000)
  svg('results/figures/supplementary/3_3hallmark_pathway_score_fib_sample_lev.svg', width = 15, height = 15)
  pheatmap::pheatmap(t(scale(t(hsi_NES))), #scale = 'row',
                     annotation_colors = ann.col, 
                     annotation_col = ann.df, cluster_cols = F, 
                     color = colorRampPalette(c("navy", "white", "red"))(70),
                     show_colnames = F, fontsize_row = 10, cutree_rows = 10)
  dev.off()
}

# 1.2. count of mut on gene level
{
  #data_mtx <- as.matrix(countdata)
  mut_count <- mut_call |>
    dplyr::arrange(TF) |>
    mutate(mut=paste0(gene,'_',Func)) |>
    dplyr::group_by(gene, sample) |>
    dplyr::summarise(N=n()) |>
    tidyr::spread(key=sample, value=N) |>
    as.data.frame() |>
    tibble::column_to_rownames(var='gene')
  mut_count[is.na(mut_count)]=0
  
  ann.df <- colnames(mut_count) |>
    as.data.frame() |>
    dplyr::rename_with(~'sample') |>
    left_join(dna_cohort) |>
    dplyr::arrange(TF) |>
    tibble::column_to_rownames(var='sample') |>
    dplyr::select(ctDNA=TF, purity, tissue_gr, treatment, Stage)
  
  stage <- colorRampPalette(RColorBrewer::brewer.pal(4, "Blues"))(4)#[2:5]
  tis <- colorRampPalette(RColorBrewer::brewer.pal(4, "Dark2"))(4)#[2:5]
  ann.col=list(Stage=c('I'=stage[1],'II'=stage[2], 'III'=stage[3],'IV'=stage[4]),
               treatment=c('NACT'='#794044','PDS'='#088da5'),
               tissue_gr=c('Ome'=tis[1], 'Per'=tis[2], 'Ova_Tub'=tis[3], 'Other'=tis[4]))
  mut_count <- as.matrix(mut_count[,rownames(ann.df)])
  
  hsi_NES <- GSVA::gsva(as.matrix(mut_count),
                        pathways.hallmark,
                        #method='zscore',
                        min.sz=0,
                        max.sz=1000)
  svg('results/figures/supplementary/hallmark_pathway_score_mutation_count.svg', width = 15, height = 15)
  pheatmap::pheatmap(hsi_NES, scale = 'row',annotation_colors = ann.col, 
                     annotation_col = ann.df, cluster_cols = T, 
                     color = colorRampPalette(c("navy", "white", "red"))(70),
                     show_colnames = F, fontsize_row = 10, cutree_rows = 10)
  dev.off()
}

# II. Annotation of driver gene CNV status
# MECOM, CCNE1, MYC, KRAS
driv_cna <- read.delim('data/input/WGSdata/driver_CNA.tsv') |>
  dplyr::group_by(patient, geneSymbol) |>
  dplyr::summarise(CN = mean(CN),
                   minor=mean(minor),
                   major=mean(major),
                   logR=mean(logR)) |>
  as.data.frame() |>
  dplyr::select(patient, geneSymbol, CN) |>
  tidyr::spread(key=geneSymbol, value=CN) |>
  as.data.frame()
  
cna <- read.delim('data/input/WGSdata/driver_CNA.tsv') |>
  left_join(dna_cohort |>
              dplyr::select(sample, TF, ctdna_lev)) |>
  mutate(status=ifelse(CN>=8,'amp', 'norm'))
  

cna |>
  ggplot(aes(x = status, y = TF)) +
  geom_boxplot(aes(fill=status), alpha=.7, show.legend = F) +
  scale_fill_brewer(palette = 'Set1', direction = -1) +
  facet_wrap(~geneSymbol, scales = 'free')+
  stat_compare_means() +
  theme_minimal()

cna |> 

length(unique(driv_cna$patient))

ann.df <- colnames(Z.n) |>
  as.data.frame() |>
  dplyr::rename_with(~'sample') |>
  left_join(samples) |>
  dplyr::arrange(TF) |>
  mutate(tmp=sub('_RNA.*', '', sample)) |>
  inner_join(driv_cna) |>
  tibble::column_to_rownames(var='sample') |>
  dplyr::select(ctDNA=TF, EOC, tissue_gr, treatment, Stage, CCNE1:MYC)

pheatmap::pheatmap(Z.n, scale = 'row', annotation_col = ann.df)



{
  #________ Heatmap Plot_____________#
  # Tidy the results:
  fgseaRes <- hallmark.res
  fgseaResTidy <- hallmark.res %>%
    as_tibble() %>%
    arrange(desc(NES)) # order by normalized enrichment score (NES)
  
  # To see what genes are in each of these pathways:
  gene.in.pathway <- pathways.hallmark %>% 
    tibble::enframe("pathway", "SYMBOL") %>% 
    unnest(cols = c(SYMBOL)) %>% 
    as.data.frame() %>%
    inner_join(res %>% as.data.frame(), by="SYMBOL")
  
  # pathways with significant enrichment score
  fgseaResTidy$adjPvalue <- ifelse(fgseaResTidy$padj <= 0.01, "significant", "non-significant")
  cols <- c("non-significant" = "grey", "significant" = "red")
  
  sig.path <- fgseaResTidy$pathway[fgseaResTidy$adjPvalue == "significant"]
  sig.gen <- unique(na.omit(gene.in.pathway$SYMBOL[gene.in.pathway$pathway %in% sig.path]))
  
  ### create a new data-frame that has '1' for when a gene is part of a term, and '0' when not
  h.dat <- dcast(gene.in.pathway[, c(1,2)], SYMBOL~pathway)
  rownames(h.dat) <- h.dat$SYMBOL
  h.dat <- h.dat[, -1]
  
  h.dat <- h.dat[rownames(h.dat) %in% sig.gen, ]
  h.dat <- h.dat[, colnames(h.dat) %in% sig.path]
  
  tt <- h.dat
  tt[!is.na(tt)]=1
  tt[is.na(tt)]=0
  h.dat <- tt
  h.dat <- h.dat %>% mutate_all(as.numeric)
  # keep those genes with 3  or more occurnes
  table(data.frame(rowSums(h.dat)))
  
  h.dat <- h.dat[data.frame(rowSums(h.dat)) >= 4, ]
  
  #
  topTable <- res[res$SYMBOL %in% rownames(h.dat), ]
  rownames(topTable) <- topTable$SYMBOL
  # match the order of rownames in toptable with that of h.dat
  topTableAligned <- topTable[which(rownames(topTable) %in% rownames(h.dat)),]
  topTableAligned <- topTableAligned[match(rownames(h.dat), rownames(topTableAligned)),]
  all(rownames(topTableAligned) == rownames(h.dat))
  
  # colour bar for -log10(adjusted p-value) for sig.genes
  dfMinusLog10FDRGenes <- data.frame(-log10(
    topTableAligned[which(rownames(topTableAligned) %in% rownames(h.dat)), 'padj']))
  dfMinusLog10FDRGenes[dfMinusLog10FDRGenes == 'Inf'] <- 0
  
  # colour bar for fold changes for sigGenes
  dfFoldChangeGenes <- data.frame(
    topTableAligned[which(rownames(topTableAligned) %in% rownames(h.dat)), 'log2FoldChange'])
  
  # merge both
  dfGeneAnno <- data.frame(dfMinusLog10FDRGenes, dfFoldChangeGenes)
  colnames(dfGeneAnno) <- c('Gene score', 'Log2FC')
  dfGeneAnno$direction <- ifelse(dfGeneAnno$Log2FC > 0, 'Up-regulated',
                           ifelse(dfGeneAnno$Log2FC < 0, 'Down-regulated', 'Unchanged'))
  colours <- list('Log2FC' = colorRamp2(c(min(dfGeneAnno$Log2FC)*1.2, 0, 
                                          max(dfGeneAnno$Log2FC)*1.2), 
                                        c("darkblue", "white", "maroon")),
                  'direction' = c('Up-regulated' = 'maroon', 
                                  'Down-regulated' = 'darkblue', 
                                  'Unchanged'='grey'))
  
  haGenes <- rowAnnotation(
    df = dfGeneAnno,
    col = colours,
    width = unit(1,'cm'),
    annotation_name_side = 'top')
  
  # Now a separate color bar for the GSEA enrichment padj. This will 
  # also contain the enriched term names via annot_text()
  
  # colour bar for enrichment score from fgsea results
  dfEnrichment <- fgseaRes[, c("pathway", "NES")]
  dfEnrichment <- dfEnrichment[dfEnrichment$pathway %in% colnames(h.dat)]
  dd <- dfEnrichment$pathway
  dfEnrichment <- dfEnrichment[, -1]
  rownames(dfEnrichment) <- dd
  #colnames(dfEnrichment) <- 'Normalized\n Enrichment score'
  col_fun = colorRamp2(c(-4, 0, 4), c("darkblue", "white", "maroon"))
  haTerms <- HeatmapAnnotation(
    df=dfEnrichment,
    col = list(NES = col_fun),
    Term = anno_text(
      colnames(h.dat),
      rot = 45, #location = unit(1, 'npc'),
      just = 'right',
      gp = gpar(fontsize = 12)),
    annotation_height = unit.c(unit(.5, 'cm'), unit(1, 'cm')),
    annotation_name_side = 'left')
  # now generate the heatmap
  hmapGSEA <- Heatmap(h.dat, km = 4, #name = '', split
                      name = 'GSEA hallmark pathways enrichment',
                      split = dfGeneAnno$direction,
                      col = c('0' = 'white', '1' = 'forestgreen'),
                      rect_gp = gpar(col = 'grey85'),
                      cluster_rows = TRUE,
                      show_row_dend = TRUE,
                      row_title = 'Top Genes',
                      row_title_side = 'left',
                      row_title_gp = gpar(fontsize = 11, fontface = 'bold'),
                      row_title_rot = 90,
                      show_row_names = TRUE,
                      row_names_gp = gpar(fontsize = 11, fontface = 'bold'),
                      row_names_side = 'left',
                      row_dend_width = unit(35, 'mm'),
                      cluster_columns = TRUE,
                      show_column_dend = TRUE,
                      column_title = 'Enriched terms',
                      column_title_side = 'top', 
                      column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
                      column_title_rot = 0,
                      show_column_names = FALSE,
                      show_heatmap_legend = FALSE,
                      clustering_distance_columns = 'euclidean',
                      clustering_method_columns = 'ward.D2',
                      clustering_distance_rows = 'euclidean',
                      clustering_method_rows = 'ward.D2',
                      bottom_annotation = haTerms)
  
  svg("results/figures/others/GSEA_enrichment_2.svg", width=20, height=40)
  draw(hmapGSEA + haGenes,
       heatmap_legend_side = 'right',
       annotation_legend_side = 'right')
  dev.off()
}




# string DB
library(STRINGdb)
string_db <- STRINGdb$new( version="11.5", species=9606,score_threshold=200, network_type="full", input_directory="")

diff_exp <- as.data.frame(res)
diff_exp <- diff_exp %>%
  dplyr::arrange(padj)

example1_mapped <- string_db$map( diff_exp, "SYMBOL", removeUnmappedRows = TRUE )

hits <- example1_mapped$STRING_id[1:500]
p1 <- string_db$plot_network( hits )

example1_mapped_pval05 <- string_db$add_diff_exp_color( subset(example1_mapped, padj<0.05),logFcColStr="log2FoldChange")
hits05 <- example1_mapped_pval05$STRING_id[1:200]
payload_id <- string_db$post_payload( example1_mapped_pval05$STRING_id,colors=example1_mapped_pval05$color )
string_db$plot_network( hits05, payload_id=payload_id )

clustersList <- string_db$get_clusters(example1_mapped_pval05$STRING_id[1:200])
par(mfrow=c(2,3))
for(i in seq(1:6)){
  string_db$plot_network(clustersList[[i]])}

hits <- example1_mapped$STRING_id[1:500]
enrichment <- string_db$get_enrichment(hits, category = "KEGG")
