source('./load_data.R')

library(fmsb)
library(ggradar)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(igraph)

# load data
## MDBSig hallmark
pathways.hallmark <- readRDS('int/hallmark_pathway.rds')

### processing - all samples
{
  rna_cohort <- rna_samples |>
    mutate(patient = sub('_.*', '', sample),
           tissue = stringr::str_match(sample, "^[^_]+_[piro]\\d?(LN|[A-Z][a-z]+)")[, 2]) |>
    inner_join(ctdna_cohort |>
                dplyr::select(patient, TF=ctdna_fraction, Stage)) |>
    mutate(pheno=ifelse(TF<thresh1, 'low', ifelse(TF<thresh2, 'med', 'high')))
  
  # EOC expression analysis
  # solid tumor only
  samples <- rna_cohort |> 
    # get one sample per patient only
    dplyr::slice_max(EOC, by = c(patient)) |>
    dplyr::arrange(desc(TF))
  nrow(samples) == length(unique(samples$patient))
  
  {
    EOC.Z <- rna_z_prim[,grep('EOC', colnames(rna_z_prim))]
    colnames(EOC.Z) = sub('.EOC', '', colnames(EOC.Z))
    EOC.Z <- EOC.Z[, samples$sample]
    ncol(EOC.Z) == nrow(samples)  
    EOC.Z[1:2,1:4]
    rna_g_prim <- rna_g
    G.EOC <- as.matrix(rna_g_prim[, samples$sample])
    W <- as.matrix(samples %>% 
                     tibble::column_to_rownames(var='sample') |>
                     dplyr::select(EOC, Fibroblast, Immune, Unknown))
    W.EOC <- W[ which(W[,1] > 0.25) , grep('EOC', colnames(W)) ]
    
    # normalisation count data with size factor
    Z.n <- diag.mul.r( as.matrix(EOC.Z), 1./ (G.EOC * W.EOC )) 
    max(Z.n)
    
    # assign count matrix
    countdata <- EOC.Z
    
    # Run DESeq2
    nrow(samples) == ncol(countdata)
    col_data <- samples |>
      dplyr::select(sample, tissue_gr, pheno, TF) %>%
      tibble::column_to_rownames(var='sample') |>
      mutate(detect=ifelse(pheno!='low','pos',pheno))
    dim(col_data)
    dds <- DESeq2::DESeqDataSetFromMatrix(countData=round(countdata[,rownames(col_data)]), 
                                                 colData=col_data, 
                                                 design= ~ 1 + pheno)
    # filtering
    #table(col_data$pheno)
    smallestGroupSize <- 26
    keep <- rowSums(DESeq2::counts(dds) > 0) >= smallestGroupSize
    dds <- dds[keep,]
    dim(dds)
    # normalisation factor
    normFactors <- (G.EOC * W.EOC )
    dds$sizeFactor = unlist(normFactors)
    dds$pheno = relevel(dds$pheno, "low")
    # run DESeq
    dds <- DESeq2::DESeq(dds)
    resultsNames(dds)
    
    ########## high versus low groups
    res_high_low <- DESeq2::lfcShrink(dds, coef='pheno_high_vs_low', type="apeglm")
    
    res_high_low$SYMBOL = geneInfo$geneName[match(sub( '[.].*$', '', rownames(res_high_low)), 
                                                  sub( '[.].*$', '', geneInfo$geneID ))]
    res_high_low_Sig <- res_high_low[ which(res_high_low$padj < 0.05 ), ]
  }
}

# list of DEGs: Supplementary table S4
res_high_low_Sig |>
  as.data.frame() |>
  write.table('results/outputs/DGE_highvslow_shrink.tsv', col.names = T, row.names = F, 
              quote = F, sep = '\t')
##### fig2a: Volcano plot
{
  vol.df.high.low <- res_high_low %>% 
    as.data.frame() %>%
    filter(!is.na(padj)) |>
    arrange(padj) %>%
    dplyr::rename(geneName = SYMBOL) %>%
    mutate(thresh=ifelse(log2FoldChange > 1 & padj < .05, 'high ctDNA', 
                         ifelse(log2FoldChange < -1 & padj < .05, 'low ctDNA', '')),
           gene = ifelse(abs(log2FoldChange) > 1 & padj < .05, 
                         geneName, NA))
  
  fig2a <- vol.df.high.low %>% 
    mutate(y = -log10(padj)) %>% 
    ggplot(., aes(x = log2FoldChange, y = y, 
                  color=thresh)) + 
    geom_point(size=.4)+
    geom_text_repel(aes(label=gene), 
                    size=4, family = 'Arial') +
    geom_vline(xintercept = -1, linetype="dotdash", color=list_color[3]) +
    geom_vline(xintercept = 1, linetype="dotdash", color=list_color[3])+
    geom_hline(yintercept = -log10(.05), linetype="dotdash", color=list_color[3])+
    scale_colour_manual(values = c("high ctDNA" = list_color[3], "low ctDNA"=list_color[1])) +
    labs(y='Log10 padj', color='Expressed in') + 
    theme_classic()+
    theme(title = element_text(size = 13),
          legend.position = 'none',
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size=12), 
          axis.text.x = element_text(size=10),
          legend.text = element_text(size=10),
          legend.title = element_text(size=15),
          strip.text.x = element_text(size = 15))
}

##### fig3b: heatmap of 3 genes CIITA, MUC4, and MUC1
{
  gene_markers <- c('CIITA', 'MUC4', 'MUC1')
  gene_markers_df <- Z.n %>% 
    as.data.frame() %>%
    tibble::rownames_to_column(var = 'gene') |>
    filter(gene %in% gene_markers) %>%
    left_join(ctdna_cohort |>
                dplyr::select(patient, ctdna_lev, Stage, treatment)) |>
    dplyr::select(-sample) |>
    dplyr::group_by(patient, ctdna_lev, Stage, treatment) |>
    dplyr::summarise_all(.funs = function(x) {mean(x)}) |>
    as.data.frame()

  gene_mx <- gene_markers_df[, c('CIITA','MUC4','MUC1')]
  annot_mx <- gene_markers_df[,c(1:4,ncol(gene_markers_df))] |> 
    tibble::column_to_rownames(var='patient')

  row_ha = ComplexHeatmap::rowAnnotation('ctdna level' = annot_mx$ctdna_lev, 
                                       'ctdna fraction' = annot_mx$TF, 
                                       treatment = annot_mx$treatment,
                                       #FIGO = annot_mx$Stage,
                                       col = list('ctdna level' = c("low" = list_color[1], "med" = list_color[2], "high" = list_color[3]),
                                                  'ctdna fraction' = colorRamp2::colorRamp2(c(0, .06, .4), c("white", "#C5D4B9", "#7ABC6D")),
                                                  treatment=c('NACT'='#794044','PDS'='#088da5')))

  fig2b <- ComplexHeatmap::Heatmap(scale(as.matrix(log2(gene_mx+1e-3))), 
                        left_annotation = row_ha,
                        name = "z-score") 
}

#### figure 2c: gene set enrichment analysis
# hallmark pathway analysis
{
  pathways.hallmark <- readRDS('int/hallmark_pathway.rds')
  
  res <- results(dds, contrast = c('pheno', 'high', 'low'))
  res$SYMBOL = geneInfo$geneName[match(sub( '[.].*$', '', rownames(res)), 
                                       sub( '[.].*$', '', geneInfo$geneID ))]
  
  res2 <- as.data.frame(res) %>% 
    filter(!is.na(padj),
           !is.na(SYMBOL),
           !duplicated(SYMBOL)) %>% 
    dplyr::arrange(padj) %>% 
    dplyr::select(symbol=SYMBOL, stat) 
  ranks <- tibble::deframe(res2)
  
  # hallmark database
  hallmark.res <- fgsea::fgsea(pathways=pathways.hallmark, stats=ranks, 
                               minSize = 20, maxSize = 500,
                               nperm=1000)
  
  hallmark.resTidy <- hallmark.res %>%
    as_tibble() %>%
    arrange(desc(NES)) |>
    as.data.frame()
  
  fig2c <- hallmark.resTidy %>%
    filter(padj < 0.05) %>% 
    mutate(pathway=sub('HALLMARK_', '', pathway),
           gr = ifelse(NES > 1, 'high', 'low'),
           logp = -log10(padj)) |>
    ggplot( aes(x=NES,
                y=reorder(pathway, NES))) +
    geom_col(aes(fill= factor(gr, levels=c('low','high'),
                              labels = c('low ctDNA',
                                         'high ctDNA'))),
             width = .5, show.legend = F) +
    scale_fill_manual(values = c(list_color[1], list_color[3])) +
    scale_color_manual(values = c('grey', 'black')) +
    labs(x="Normalised Enrichment Score", 
         y="MSigDB Gene Sets",
         title="") + 
    theme_bw()+
    theme(axis.title = element_text(size=12, family = 'Arial'), 
          axis.text.y = element_text(size=9, family = 'Arial'),
          axis.text.x = element_text(size=9, family = 'Arial'),
          #legend.position = 'none',# c(.2,.9),
          legend.text = element_text(size=8, family = 'Arial'),
          legend.title = element_text(size=8, family = 'Arial'),
          strip.text.x = element_text(size = 15, family = 'Arial'))
}
################ done ##################
