source('src/00_functions.R')

library(fmsb)
library(ggradar)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(igraph)

# load data
load('~/OmicsIntegration/data/RNA/pretrt_RNA.RData')
# contains related data tables: 
# - bulk_rna: Bulk RNA data from primary samples
# - rna_samples: proportion of cell populations decomposed from primary bulk
# - rna_z: weight for normalization
# - rna_g: decomposed expression

rna_z_prim <- rna_z
rownames(rna_z_prim) = rna_z_prim$geneID

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
    smallestGroupSize <- 24
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
##### fig3a: Volcano plot
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
  
  fig3a <- vol.df.high.low %>% 
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

##### fig3b,c: ORA analysis
{
  hallmark_desc <- names(pathways.hallmark)
  names(hallmark_desc) <- names(pathways.hallmark)
  hallmark_sets <- pathfindR::fetch_gene_set(gene_sets = 'Custom',
                                             custom_genes = pathways.hallmark,
                                             custom_descriptions = hallmark_desc)
  ############# pathfindR
  ########### high versus low ctDNA #############
  highvslow_pathfindRinput <- res_high_low |>
    as.data.frame() |>
    dplyr::select(SYMBOL, log2FoldChange, padj) |>
    filter(!duplicated(SYMBOL),
           !is.na(padj)) |>
    as.data.frame()
  dim(highvslow_pathfindRinput)
  
  ##
  output_highvslow <- pathfindR::run_pathfindR(highvslow_pathfindRinput, 
                                               gene_sets = 'Custom',
                                               custom_genes = pathways.hallmark,
                                               custom_descriptions = hallmark_desc, 
                                               min_gset_size = 20,
                                               p_val_threshold = .05,
                                               iterations = 25) |>
    mutate(ID = sub('HALLMARK_', '', ID),
           Term_Description = sub('HALLMARK_', '', Term_Description))
  
  highvslow_processed <- pathfindR::input_processing(highvslow_pathfindRinput)
  
  fig3b <- pathfindR::term_gene_heatmap(output_highvslow, 
                               num_terms = 9, low = list_color[1],high = list_color[3]) +
    theme(axis.text.x = element_text(size=7))
  
  fig3c <- pathfindR::term_gene_graph(
    output_highvslow,
    num_terms = 9, 
    layout = "stress",
    use_description = T,
    node_size = "num_genes"
  ) +
    scale_color_manual(values = c('grey', list_color[3], list_color[1]),
                       labels = c("enriched term", "up-regulated genes in high ctDNA", 
                                  "up-regulated genes in med ctDNA")) +
    labs(color='', `# genes` = '')+
    guides(fill=guide_legend(ncol=2)) +
    theme(title = element_text(family = 'Arial', size = 15),
          legend.text = element_text(size = 10, family = 'Arial')) 
}

#### figure 3d: gene set enrichment analysis
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
  
  fig3d <- hallmark.resTidy %>%
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
