############## load data
source('./load_data.R')

### fig 1a
fig1a <- ctdna_df_prim %>%
    mutate(figo = ifelse(Stage %in% c('I', 'II'), 'Early', Stage)) |>
    ggplot(aes(x=factor(figo, 
                        labels = c('Benign\n(n=10)',
                                   'Early\n(n=14)',
                                   'III\n(n=85)',
                                   'IV\n(n=19)')),
               y = TF)) +
    geom_boxplot(width=.3,outlier.shape = NA,
                 alpha=.4, show.legend = F) + 
  geom_jitter(aes(color=ctdna_lev), 
              size=3, position=position_jitter(0.2)) + 
  stat_compare_means(comparisons = comp1,
                     label.y = c(.12,.16,.20,.24) ) +
  geom_hline(yintercept = thresh1, color=list_color[1]) + 
  geom_hline(yintercept = thresh2, color=list_color[3]) + 
  scale_colour_manual(values = c("high" = list_color[3], 
                                 "med" = list_color[2],
                                 "low"=list_color[1])) +
  labs(x='Stages at diagnosis', y ='ctDNA fraction',
       color = 'ctDNA level') + 
  theme_classic() +
  theme(axis.text.x = element_text(size=13,  family = 'Arial', colour = 'black'),
        axis.title = element_text(size=15,  family = 'Arial', colour = 'black'),
        axis.text = element_text(size=13, family = 'Arial', colour = 'black'),
        legend.position = 'none', 
        legend.text = element_text(size=11, family = 'Arial', colour = 'black'),
        legend.title = element_text(size=13, family = 'Arial', colour = 'black'))

############## GSVA scores for hallmark geneset
## MDBSig hallmark
pathways.hallmark <- readRDS('int/hallmark_pathway.rds')

# EOC Component
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
  
  EOC.Z <- exp_mat[, samples$sample]
  
  G.EOC <- as.matrix(rna_g[, samples$sample])
  W <- as.matrix(samples %>% 
                   tibble::column_to_rownames(var='sample') |>
                   dplyr::select(EOC, Fibroblast, Immune, Unknown))
  W.EOC <- W[ which(W[,1] > 0.25) , grep('EOC', colnames(W)) ]
  
  # normalisation count data with size factor
  Z.n <- diag.mul.r( as.matrix(EOC.Z), 1./ (G.EOC * W.EOC )) 
  count_mtx <- Z.n[rowSums(Z.n)>0,] 
  
  data_mtx <- as.matrix(count_mtx)
  
  ## build GSVA parameter object
  gsvapar <- GSVA::gsvaParam(data_mtx, pathways.hallmark, maxDiff=TRUE,
                       minSize = 20, maxSize = 1000, kcdf = 'Poisson')
  hsi_NES <- GSVA::gsva(gsvapar)
}

## analysis for constructing figure 1b
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
  
  ## get significant pathways
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

  # significant pathways below are also provided in supplementary table S5
  sig.pathways <- pathways[pathways$p_val<.05,]$pathway

  interested_pathways <- hsi_sample[,c('patient',
                                       colnames(hsi_sample)[grepl('HALLMARK',colnames(hsi_sample))])]
  
  names(interested_pathways) <- sub('HALLMARK_', '', names(interested_pathways))
  interested_pathways <- interested_pathways |>
    dplyr::group_by(patient) |>
    dplyr::summarise_all(.funs = mean) |>
    left_join(ctdna_cohort |>
                dplyr::select(patient, 
                              TF=ctDNA_fraction)) |>
    dplyr::arrange(TF) |>
    tibble::column_to_rownames(var='patient') |>
    dplyr::select(-TF)
  
  interested_pathways <- interested_pathways[,sig.pathways]
  
  ann.row <- interested_pathways |>
    tibble::rownames_to_column(var='patient') |>
    dplyr::select(patient) |>
    left_join(ctdna_cohort |>
                dplyr::select(patient, ctdna_tf=TF, ctdna_lev)) |>
    dplyr::arrange(ctdna_tf) |>
    tibble::column_to_rownames(var='patient') 
  
  ann.col=list(ctdna_lev = c('low'=list_color[1],
                             'med'=list_color[2], 
                             'high'=list_color[3]))
  breaksList = seq(-3, 3, by = .4)
  
  ##### fig1b: heatmap
  fig1b <- pheatmap::pheatmap(interested_pathways, annotation_row = ann.row, 
                     annotation_colors = ann.col, 
                     cluster_rows = F, show_rownames = F, 
                     scale = 'column',
                     color = colorRampPalette(c("#3399ff", "white", "#ff6666"))(length(breaksList)),
                     breaks = breaksList, gaps_row = c(cumsum(c(49,27,26)),cumsum(c(49,27,26))), border_color = 'black',
                     fontsize_col = 8, cutree_rows = 4, cutree_cols = 2,
                     height = 3, width = 6)
}

## running PCA
{
  data <- t(hsi_NES) %>%
    as.data.frame()
  colnames(data) = sub('HALLMARK_','',colnames(data))
  pathway2check <- sub('HALLMARK_','',
                       hallmark.resTidy$pathway[hallmark.resTidy$padj<.05])
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

  res.pca <- prcomp(data[, 1:(ncol(data)-2)], scale=F)
  # Visualize eigenvalues/variances
  factoextra::fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 70))
  # Extract the results for variables
  var <- factoextra::get_pca_var(res.pca)
  # Contributions of variables to PC1
  factoextra::fviz_contrib(res.pca, choice = "var", axes = 1, top = 10, )
  # Contributions of variables to PC2
  factoextra::fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)

  ### figure 1c
  fig1c <- factoextra::fviz_pca_biplot(res.pca, geom.ind = "point",axes = c(1, 2),
                                fill.ind=factor(data$ctdna_lev, 
                                                 levels = c('low', 'med', 'high')),
                                pointshape=21, pointsize=5, addEllipses = F,
                                col.var="contrib",repel = T, 
                                geom.var = c("arrow"),arrowsize = 1.3,
                                labelsize = 5,alpha.ind=.7,
                                gradient.cols = "BuPu",# c("#00AFBB", "#E7B800", "#FC4E07"),
                                palette = list_color) + 
      labs(x = 'PCA1 (41.4%)', y = 'PCA2 (17.8%)') +
      theme_classic()+
      theme(#title = element_blank(),
            axis.title = element_text(size=15, family = 'Arial'), 
            axis.text.y = element_text(size=12, family = 'Arial'),
            axis.text.x = element_text(size=12, family = 'Arial'),
            legend.position = "right", 
            legend.text = element_text(size=20, family = 'Arial'),
            legend.title = element_text(size=20, family = 'Arial'),
            strip.text.x = element_text(size = 15, family = 'Arial'))
}

################ done ##################
