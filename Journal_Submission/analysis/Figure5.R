############## GSVA scores for hallmark geneset

pathways.hallmark <- readRDS(paste0(~/int/,'hallmark_pathway.rds'))

# I. Sample level - highest cancer fraction sample from 1 tissue group 
# Calculate pathway score by sample
{
  samples <- rna_cohort |>
    dplyr::slice_max(EOC, by = c(patient, tissue_gr)) |>
    dplyr::arrange(desc(TF))
  
  exp_mat <- rna_z_prim
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
                      
  count_mtx <- Z.n[rowSums(Z.n)>0,] 
  
  data_mtx <- as.matrix(count_mtx)
  
  hsi_NES <- GSVA::gsva(as.matrix(data_mtx),
                        pathways.hallmark,
                        method="gsva",
                        min.sz=20,
                        max.sz=1000, 
                        kcdf="Poisson")
}

# II. Analysis from pathway scores
{
  hsi_sample <- t(hsi_NES) %>%
    as.data.frame() |>
    tibble::rownames_to_column(var='sample') %>%
    mutate(patient = sub('_.*', '', sample)) |>
    left_join(samples) 
  
  # 4 hallmarks proliferation-related gene sets
  proliferate_pathways <- c( 'E2F_TARGETS', 'G2M_CHECKPOINT',
                            'MYC_TARGETS_V1','MYC_TARGETS_V2')
  
  # metabolism-related pathway
  metabolic_pathways <- c('OXIDATIVE_PHOSPHORYLATION', 'BILE_ACID_METABOLISM',
                          'HEME_METABOLISM', 'FATTY_ACID_METABOLISM')
  
  # immune-related pathway
  imm_pathways <- c('IL6_JAK_STAT3_SIGNALING',
                    'ALLOGRAFT_REJECTION', 
                    'COMPLEMENT', 'INFLAMMATORY_RESPONSE')
  
  # signaling-related pathway
  signaling_pathways <- c('MTORC1_SIGNALING',
                          'HEDGEHOG_SIGNALING',
                          'TNFA_SIGNALING_VIA_NFKB','TGF_BETA_SIGNALING')
  
  # filter pathways of interested listed above
  interested_pathways <- hsi_sample[,c('patient',
                      colnames(hsi_sample)[grepl('HALLMARK',colnames(hsi_sample))])]
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
  dev.off()
  
  # pca
  {
    data <- t(hsi_NES) %>%
      as.data.frame()
    colnames(data) = sub('HALLMARK_','',colnames(data))
    data <- data %>%  dplyr::select(colnames(.)[colnames(.) %in%  
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
  
