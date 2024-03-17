source('src/00_functions.R')

library(fmsb)
library(ggradar)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(igraph)
#devtools::install_github("ricardo-bion/ggradar")
#devtools::install_github('KChen-lab/METAFlux')

# load data
load('~/mnt/storageBig8/work/nguyenma/projects/RNA/OmicsIntegration/data/processed/RNA/pretrt_RNA.RData')
# contains related data tables: 
# - bulk_rna: Bulk RNA data from primary samples
# - rna_samples: proportion of cell populations decomposed from primary bulk
# - rna_z: weight for normalization
# - rna_g: decomposed expression

rna_z_prim <- rna_z
rna_z_prim[1:2,1:5]
dim(rna_z_prim)
rownames(rna_z_prim) = rna_z_prim$geneID

# check row name of count matrix
rna_z_prim <- rna_z_prim[,-1]
rna_z_prim[1:2,1:3]
### processing - all samples
{
  rna_cohort <- rna_samples |>
    mutate(patient = sub('_.*', '', sample),
           tissue = stringr::str_match(sample, "^[^_]+_[piro]\\d?(LN|[A-Z][a-z]+)")[, 2]) |>
    filter(grepl('_p', sample),
           # filter patients in ctDNA cohort only
           patient %in% ctdna_cohort$patient,
           !(tissue %in% c('Asc', 'Plf'))) |>
    inner_join(ctdna_cohort |>
                dplyr::select(patient, TF, Stage)) |>
    mutate(pheno=ifelse(TF<thresh1, 'low', ifelse(TF<thresh2, 'med', 'high')),
           tissue_gr = ifelse(tissue %in% c('Ome', 'Per'), tissue, 
                              ifelse(tissue %in% c('Ova', 'Tub', 'Adn', 'And'),
                                     'Ova_Tub','Other'))) |>
    # exclude samples Daria mentioned
    filter(sample!='H284_pOvaR1_RNA1')
  
  # EOC expression analysis
  # solid tumor only
  samples <- rna_cohort |> 
    # filter advanced patients only
    #filter(Stage %in% c('III','IV')) |>
    #mutate(pheno=ifelse(TF>=thresh, 'detected', 'undetected')) |>
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
    # this is used for heatmap
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
    svg('results/figures/final/Supp/maplot_high_low.svg', width = 6, height = 5)
    plotMA(res_high_low)
    dev.off()
    
    res_high_low$SYMBOL = geneInfo$geneName[match(sub( '[.].*$', '', rownames(res_high_low)), 
                                                  sub( '[.].*$', '', geneInfo$geneID ))]
    res_high_low_Sig <- res_high_low[ which(res_high_low$padj < 0.05 ), ]
    dim(res_high_low_Sig)
    res_high_low_Sig |>
      as.data.frame() |>
      write.table('results/outputs/DGE_highvslow_shrink.tsv', col.names = T, row.names = F, 
                  quote = F, sep = '\t')
    # Volcano plot
    vol.df.high.low <- res_high_low %>% 
      as.data.frame() %>%
      filter(!is.na(padj)) |>
      arrange(padj) %>%
      dplyr::rename(geneName = SYMBOL) %>%
      mutate(thresh=ifelse(log2FoldChange > 1 & padj < .05, 'high ctDNA', 
                           ifelse(log2FoldChange < -1 & padj < .05, 'low ctDNA', '')),
             gene = ifelse(abs(log2FoldChange) > 1 & padj < .05, 
                           geneName, NA))

    fig_high_low <- vol.df.high.low %>% 
      #filter(log2FoldChange < 5)|>
      mutate(y = -log10(padj)) %>% 
      ggplot(., aes(x = log2FoldChange, y = y, 
                    color=thresh)) + 
      geom_point(size=.4)+
      geom_text_repel(aes(label=gene), 
                      size=4, family = 'Arial') +
      #geom_text_repel() + 
      geom_vline(xintercept = -1, linetype="dotdash", color=list_color[3]) +
      geom_vline(xintercept = 1, linetype="dotdash", color=list_color[3])+
      geom_hline(yintercept = -log10(.05), linetype="dotdash", color=list_color[3])+
      scale_colour_manual(values = c("high ctDNA" = list_color[3], "low ctDNA"=list_color[1])) +
      labs(y='Log10 padj', color='Expressed in') + 
      #lims( x = c(-4,6)) +
      theme_classic()+
      theme(title = element_text(size = 13),
            legend.position = 'none',
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_text(size=12), 
            axis.text.x = element_text(size=10),
            #legend.position = c(.85,.75),
            legend.text = element_text(size=10),
            legend.title = element_text(size=15),
            strip.text.x = element_text(size = 15))
    svg('results/figures/final/Main/fig4a.svg', width = 7, height = 7)
    fig_high_low
    dev.off()

    # figure 3c: pathway enrichment analysis
    # pathway analysis
    {
      dat_path <- '~/mnt/storageBig8/work/nguyenma/projects/RNA/OmicsIntegration/data/temp/'
      #hallmark.db <- read.gmt('data/input/pathways/h.all.v2023.2.Hs.symbols.gmt')
      pathways.hallmark <- readRDS(paste0(dat_path,'hallmark_pathway.rds'))
      
      #pathways.kegg <- readRDS(paste0(dat_path,'kegg_pathway.rds'))
      res <- results(dds, contrast = c('pheno', 'high', 'low'))
      res$SYMBOL = geneInfo$geneName[match(sub( '[.].*$', '', rownames(res)), 
                                           sub( '[.].*$', '', geneInfo$geneID ))]
      res2 <- as.data.frame(res) %>% #head(2)
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
      
      fig4e <- hallmark.resTidy %>%
        filter(padj < 0.05) %>% #head(2)
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
      
      svg('results/figures/final/Main/fig4e.svg', width = 5, height = 6.5)
      fig4e
      dev.off()
      
      
      
      fig3d <- ggarrange(
        plotEnrichment(pathways.hallmark[["HALLMARK_E2F_TARGETS"]], 
                       ranks) + labs(title="HALLMARK_E2F_TARGETS"),
        plotEnrichment(pathways.hallmark[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]],
                       ranks) + labs(title="HALLMARK_OXIDATIVE_PHOSPHORYLATION"),
        plotEnrichment(pathways.hallmark[["HALLMARK_INFLAMMATORY_RESPONSE"]],
                       ranks) + labs(title="HALLMARK_INFLAMMATORY_RESPONSE"),
        plotEnrichment(pathways.hallmark[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]],
                       ranks) + labs(title="HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
      nrow=4)
      svg('results/figures/submission/Main/fig34.svg', width = 9, height = 7)
      fig3d
      dev.off()
      
      # gene signatures for pathway heatmap
      {
        #________ Heatmap Plot_____________#
        # Tidy the results:
        fgseaRes <- hallmark.res
        fgseaResTidy <- hallmark.resTidy
        
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
        h.dat <- reshape2::dcast(gene.in.pathway[, c(1,2)], SYMBOL~pathway)
        heat_df <- gene.in.pathway |>
          #mutate(pathway = sub('HALLMARK_', '', pathway)) |>
          #tibble::column_to_rownames(var='pathway') |>
          dplyr::select(SYMBOL, pathway, log2FoldChange) |>
          filter(SYMBOL %in% sig.gen,
                 pathway %in% sig.path) |>
          tidyr::spread(key=pathway, value=log2FoldChange) |>
          tibble::column_to_rownames(var='SYMBOL') |>
          mutate_all(function(x) ifelse(is.na(x), 0, x)) 
        
        dim(heat_df)
        rownames(h.dat) <- h.dat$SYMBOL
        h.dat <- h.dat[, -1]
        
        h.dat <- h.dat[rownames(h.dat) %in% sig.gen, ]
        h.dat <- h.dat[, colnames(h.dat) %in% sig.path]
        
        tt <- h.dat
        tt[!is.na(tt)]=1
        tt[is.na(tt)]=0
        h.dat <- tt
        h.dat <- h.dat %>% mutate_all(as.numeric)
        # keep those genes with 4 or more occurnes
        idx <- apply(heat_df, 1, function(x) sum(x!=0))
        #table(data.frame(rowSums(h.dat)))
        
        h.dat <- h.dat[data.frame(rowSums(h.dat!=0)) >= 4, ]
        colnames(h.dat) = sub('HALLMARK_','',colnames(h.dat))
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
        colnames(dfGeneAnno) <- c('Gene_score', 'Log2FC')
        dfGeneAnno$direction <- ifelse(dfGeneAnno$Log2FC > 0, 'Up-regulated',
                                       ifelse(dfGeneAnno$Log2FC < 0, 'Down-regulated', 'Unchanged'))
        colours <- list('Log2FC' = colorRamp2(c(min(dfGeneAnno$Log2FC)*1.2, 0, 
                                                max(dfGeneAnno$Log2FC)*1.2), 
                                              c("darkblue", "white", "maroon")),
                        'direction' = c('Up-regulated' = '#ff6666', 
                                        'Down-regulated' = '#6897bb', 
                                        'Unchanged'='grey'))
        
        #'Gene_score' = colorRamp2(c(min(dfGeneAnno$Gene_score), 1, 
        #                            max(dfGeneAnno$Gene_score)*1.1), 
        #                          c("white", "#b4eeb4", "#008000"))
        haGenes <- rowAnnotation(
          df = dfGeneAnno,
          col = colours,
          width = unit(1,'cm'),
          annotation_name_side = 'top'
          )
        
        # Now a separate color bar for the GSEA enrichment padj. This will 
        # also contain the enriched term names via annot_text()
        
        # colour bar for enrichment score from fgsea results
        fgseaRes <- hallmark.res
        dfEnrichment <- fgseaRes[, c("pathway", "NES")]
        dfEnrichment$pathway = sub('HALLMARK_', '', dfEnrichment$pathway)
        dfEnrichment <- dfEnrichment[dfEnrichment$pathway %in% colnames(h.dat),]
        dfEnrichment <- dfEnrichment |> 
          mutate(pathway = sub('HALLMARK_', '', pathway),
                 pathway = sub('.* .','',pathway)) |>
          filter(pathway %in% colnames(h.dat)) |>
          tibble::column_to_rownames(var='pathway')
        dd <- dfEnrichment$pathway
        dfEnrichment <- dfEnrichment[, -1]
        rownames(dfEnrichment) <- dd
        #colnames(dfEnrichment) <- 'Normalized\n Enrichment score'
        col_fun = list('NES'=colorRamp2(c(-2, 0, 3.1), 
                             c("darkblue", "white", "maroon")))
        haTerms <- HeatmapAnnotation(
          df=dfEnrichment,
          col = col_fun,
          Term = anno_text(
            colnames(h.dat),
            rot = 45, #location = unit(1, 'npc'),
            just = 'right',
            gp = gpar(fontsize = 12)),
          annotation_height = unit.c(unit(.5, 'cm'), unit(1, 'cm')),
          annotation_name_side = 'left')
        # generate the heatmap
        hmapGSEA <- Heatmap(h.dat, km = 4,  #name = '', split
                            name = 'GSEA hallmark pathways enrichment',
                            split = dfGeneAnno$direction, 
                            col = c('0' = 'white', '1' = '#b0e0e6'),
                            rect_gp = gpar(col = 'grey85'), 
                            cluster_rows = TRUE,
                            show_row_dend = TRUE,
                            row_title = 'Top Genes',
                            row_title_side = 'left',
                            row_title_gp = gpar(fontsize = 15, fontface = 'bold', fontfamily='Arial'),
                            row_title_rot = 90,
                            show_row_names = TRUE,
                            row_names_gp = gpar(fontsize = 15, fontface = 'bold', fontfamily='Arial'),
                            row_names_side = 'left',
                            row_dend_width = unit(35, 'mm'),
                            cluster_columns = TRUE,
                            show_column_dend = TRUE,
                            column_title = 'Enriched terms',
                            column_title_side = 'top', 
                            column_title_gp = gpar(fontsize = 15, fontface = 'bold', fontfamily='Arial'),
                            column_title_rot = 0,
                            show_column_names = FALSE,
                            show_heatmap_legend = FALSE,
                            clustering_distance_columns = 'euclidean',
                            clustering_method_columns = 'ward.D2',
                            clustering_distance_rows = 'euclidean',
                            clustering_method_rows = 'ward.D2', 
                            bottom_annotation = haTerms, 
                            heatmap_legend_param = list(fontsize = 20, fontfamily='Arial'))
        
        svg("results/figures/submission/Main/fig35.svg", width=15, height=20)
        svg("results/figures/submission/Supp/hallmark_7.svg", width=15, height=20)
        draw(hmapGSEA + haGenes,
             heatmap_legend_side = 'right',
             annotation_legend_side = 'right')
        dev.off()
      }
      
      # KEGG
      {
        fgseaKEGGRes <- fgsea(pathways=pathways.kegg, 
                              nperm = 10000,
                              stats=ranks, minSize = 10)
        
        fgseaResKEGGTidy <- fgseaKEGGRes %>%
          as_tibble() %>%
          arrange(desc(NES))
        
        s3b <- fgseaResKEGGTidy %>%
          filter(padj <= 0.01) %>%
          ggplot( aes(reorder(pathway, NES), NES)) +
          geom_col(aes(fill= factor(NES>1, 
                                    levels = c(TRUE, FALSE),
                                    labels=c('ctDNA ≥ 3%', 'ctDNA < 3%')))) +
          scale_fill_manual(values = c('Maroon', 'Dark blue'))+
          coord_flip() +
          labs(x="Pathway", y="Normalized Enrichment Score",
               fill='ctDNA Level',
               title="KEGG pathways with GSEA - padj ≤ 0.01") + 
          theme_minimal()+
          theme(legend.position = c(.2,.9))
        
        #________ Heatmap Plot_____________#
        # Tidy the results:
        fgseaResTidy <- fgseaResKEGGTidy
        
        # To see what genes are in each of these pathways:
        gene.in.pathway <- pathways.kegg %>% 
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
        # keep those genes with 4 or more occures
        table(colSums(h.dat))
        h.dat <- h.dat[data.frame(rowSums(h.dat!=0)) >= 4, colSums(h.dat)>=30]
        h.dat <- h.dat[, colSums(h.dat)!=0]
        dim(h.dat)
        colnames(h.dat) = sub('KEGG_','',colnames(h.dat))
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
        colnames(dfGeneAnno) <- c('Gene_score', 'Log2FC')
        dfGeneAnno$direction <- ifelse(dfGeneAnno$Log2FC > 0, 'Up-regulated',
                                       ifelse(dfGeneAnno$Log2FC < 0, 'Down-regulated', 'Unchanged'))
        colours <- list('Log2FC' = colorRamp2(c(min(dfGeneAnno$Log2FC)*1.2, 0, 
                                                max(dfGeneAnno$Log2FC)*1.2), 
                                              c("darkblue", "white", "maroon")),
                        'direction' = c('Up-regulated' = '#ff6666', 
                                        'Down-regulated' = '#6897bb', 
                                        'Unchanged'='grey'))
        
        haGenes <- rowAnnotation(
          df = dfGeneAnno,
          col = colours,
          width = unit(1,'cm'),
          annotation_name_side = 'top')
        
        # Now a separate color bar for the GSEA enrichment padj. This will 
        # also contain the enriched term names via annot_text()
        
        # colour bar for enrichment score from fgsea results
        dfEnrichment <- fgseaKEGGRes[, c("pathway", "NES")]
        dfEnrichment$pathway = sub('KEGG_', '', dfEnrichment$pathway)
        dfEnrichment <- dfEnrichment[dfEnrichment$pathway %in% colnames(h.dat)]
        dd <- dfEnrichment$pathway
        dfEnrichment <- dfEnrichment[, -1]
        rownames(dfEnrichment) <- dd
        #colnames(dfEnrichment) <- 'Normalized\n Enrichment score'
        col_fun = list('NES'=colorRamp2(c(min(dfEnrichment$NES)*1.1, 0, max(dfEnrichment$NES)*1.1), 
                                        c("darkblue", "white", "maroon")))
        haTerms <- HeatmapAnnotation(
          df=dfEnrichment,
          col = col_fun,
          Term = anno_text(
            colnames(h.dat),
            rot = 45, #location = unit(1, 'npc'),
            just = 'right',
            gp = gpar(fontsize = 12)),
          annotation_height = unit.c(unit(.5, 'cm'), unit(1, 'cm')),
          annotation_name_side = 'left')
        # now generate the heatmap
        hmapGSEA_kegg <- Heatmap(h.dat, km = 4, #name = '', split
                            name = 'GSEA hallmark pathways enrichment',
                            split = dfGeneAnno$direction, 
                            col = c('0' = 'white', '1' = '#b0e0e6'),
                            rect_gp = gpar(col = 'grey85'), 
                            cluster_rows = TRUE,
                            show_row_dend = TRUE,
                            row_title = 'Top Genes',
                            row_title_side = 'left',
                            row_title_gp = gpar(fontsize = 15, fontface = 'bold', fontfamily='Arial'),
                            row_title_rot = 90,
                            show_row_names = TRUE,
                            row_names_gp = gpar(fontsize = 15, fontface = 'bold', fontfamily='Arial'),
                            row_names_side = 'left',
                            row_dend_width = unit(35, 'mm'),
                            cluster_columns = TRUE,
                            show_column_dend = TRUE,
                            column_title = 'Enriched terms',
                            column_title_side = 'top', 
                            column_title_gp = gpar(fontsize = 15, fontface = 'bold', fontfamily='Arial'),
                            column_title_rot = 0,
                            show_column_names = FALSE,
                            show_heatmap_legend = FALSE,
                            clustering_distance_columns = 'euclidean',
                            clustering_method_columns = 'ward.D2',
                            clustering_distance_rows = 'euclidean',
                            clustering_method_rows = 'ward.D2',
                            bottom_annotation = haTerms, 
                            heatmap_legend_param = list(fontsize = 20, fontfamily='Arial'))
        
        svg("results/figures/submission/Supp/kegg_cluster.svg", width=10, height=15)
        draw(hmapGSEA_kegg + haGenes,
             heatmap_legend_side = 'right',
             annotation_legend_side = 'right')
        dev.off()
      }
      
      
      # fig3c
      
      
      # we want the log2 fold change 
      original_gene_list <- res$log2FoldChange
      
      # name the vector
      names(original_gene_list) <- res$SYMBOL
      
      # omit any NA values 
      gene_list<-na.omit(original_gene_list)
      
      # sort the list in decreasing order (required for clusterProfiler)
      gene_list = sort(gene_list, decreasing = TRUE)
      
      ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
      # remove duplicate IDS (here I use "SYMBOL", but it should be whatever was selected as keyType)
      dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
      # Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
      df2 = res[res$SYMBOL %in% dedup_ids$SYMBOL,]
      
      # Create a new column in df2 with the corresponding ENTREZ IDs
      df2$Y = dedup_ids$ENTREZID
      
      # Create a vector of the gene unuiverse
      kegg_gene_list <- df2$log2FoldChange
      
      # Name vector with ENTREZ ids
      names(kegg_gene_list) <- df2$Y
      
      # omit any NA values 
      kegg_gene_list<-na.omit(kegg_gene_list)
      
      # sort the list in decreasing order (required for clusterProfiler)
      kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
      set.seed(123)
      kk2 <- gseKEGG(geneList     = kegg_gene_list,
                     organism     = 'hsa',
                     nPerm        = 10000, 
                     minGSSize = 20,
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH",
                     keyType       = "ncbi-geneid")
      
      dotplot(kk2, title = "Enriched Pathways" , split=".sign") + 
        facet_grid(.~.sign)
      emapplot(kk2)
      cnetplot(kk2, categorySize="pvalue", foldChange=gene_list)
      ridgeplot(kk2) + labs(x = "enrichment distribution")
      
      x <- enrichKEGG(gene=names(kegg_gene_list), 
                      organism = 'hsa',
                      pAdjustMethod = 'BH',
                      pvalueCutoff = 0.05, 
                      minGSSize = 20,
                      qvalueCutoff = 0.05)
      dotplot(x)
      
      x <- enrichPathway(gene=names(kegg_gene_list), 
                         pvalueCutoff = 0.05, readable=TRUE
                      pAdjustMethod = 'BH',
                      pvalueCutoff = 0.05, 
                      minGSSize = 20,
                      qvalueCutoff = 0.05)
      dotplot(x)
      
      #https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
      
      
      
      # kegg database
      fgseaKEGGRes <- fgsea(pathways=pathways.kegg, stats=ranks)
      
      fgseaResKEGGTidy <- fgseaKEGGRes %>%
        as_tibble() %>%
        arrange(desc(NES))
      
      kegg <- ggplot(fgseaResKEGGTidy %>% 
                       as.data.frame()%>%
                       filter(abs(NES)>1.5)%>%
                       mutate(pathway = tomeder(sub('KEGG_','', pathway))), 
                     aes(reorder(pathway, NES), NES)) + 
        geom_col(aes(fill=padj<0.01)) +
        coord_flip() +labs(x="Pathway", y="Normalized Enrichment Score",
                           title="KEGG pathways NES from GSEA\n on Cancer Cells Population") + 
        theme_minimal(base_size = 10) + theme(title = element_text(size=12))
    }
    
    
    # additional investigation
    {
      # (unevaluated code chunk)
      library("IHW")
      resIHW <- results(dds, filterFun=ihw)
      summary(resIHW)
      sum(resIHW$padj < 0.1, na.rm=TRUE)
      metadata(resIHW)$ihwResult
    }
    
    # evaluatio shrinkage
    resLFC <- lfcShrink(dds, coef=2)
    resNorm <- lfcShrink(dds, coef=2, type="normal")
    resAsh <- lfcShrink(dds, coef=2, type="ashr")
    par(mfrow=c(1,3), mar=c(4,4,2,1))
    xlim <- c(1,1e5); ylim <- c(-3,3)
    plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
    plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
    plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")
    
    resOrdered <- res[order(res$pvalue),]
    resSig <- subset(resOrdered, padj < 0.05)
    
    }    
    # transformation
  {}
    {
      ntd <- normTransform(dds)
      library("vsn")
      vsn::meanSdPlot(assay(ntd))
      vsd <- vst(dds, blind=FALSE)
      vsn::meanSdPlot(assay(vsd))
    }
    
    conv_tbl <- read.delim('data/input/RNAdata/entrezid_conv.tsv')
    
    resLFC <- lfcShrink(dds, type = 'ashr', contrast = c('pheno', 'TRUE', 'FALSE'), res = res)# type = 'ashr', contrast = c('pheno', 'high', 'med'))
    resLFC$SYMBOL = geneID2name$uniqName[match(sub( '[.].*$', '', rownames(resLFC)), 
                                         sub( '[.].*$', '', geneID2name$geneID ))]
    resLFC$geneType = geneID2name$geneType[match(sub( '[.].*$', '', rownames(resLFC)), 
                                           sub( '[.].*$', '', geneID2name$geneID ))]
    resLFC$ENTREZID = conv_tbl$ENTREZID[match(sub( '[.].*$', '', resLFC$SYMBOL), 
                                              sub( '[.].*$', '', conv_tbl$SYMBOL ))]
    
    sigresLFC <- resLFC |>
      as.data.frame() |>
      filter(padj < .05, abs(log2FoldChange) > 2)
    dim(sigresLFC)
    
    heat_mat <- as.data.frame(Z.n[,rownames(col_data)]) %>% 
                  filter(rownames(.)%in% rownames(sigresLFC)) %>%
                  tibble::rownames_to_column(var='geneID') |>
                  left_join(geneInfo |>
                              dplyr::select(geneID, symbol=uniqName)) |>
                  dplyr::select(-geneID) |>
                  tibble::column_to_rownames(var='symbol')
    pheatmap::pheatmap(heat_mat,
                       cluster_cols = F, annotation_col = col_data, scale = 'row')
    scale_mat <- t(scale(t(heat_mat)))
    sigGenes <- resLFC$ENTREZID[ resLFC$padj < 0.05 & 
                                    !is.na(resLFC$padj) &
                                    abs(resLFC$log2FoldChange) > 1.5 ]
    length(sigGenes)
    sigGenes <- na.exclude(sigGenes)
    geneUniverse <- resLFC$ENTREZID
    kk <- clusterProfiler::enrichKEGG(gene = sigGenes, universe = geneUniverse, 
                                      organism = 'hsa',
                                    pvalueCutoff = .05)
    #dotplot(kk)
  }
  


### processing - slice_max(EOC)
{
  # EOC expression analysis
  # solid tumor only
  samples_max_eoc <- rna_cohort |>
    mutate(pheno = ifelse(TF < .07, 'med', 'high')) |>
    dplyr::slice_max(EOC, by=patient) |>
    dplyr::arrange(TF)
  
  {
    EOC.Z <- rna_z_prim %>%
      #tibble::column_to_rownames(var='geneID') |>
      as.data.frame() %>%
      dplyr::select(colnames(.)[grepl('EOC', colnames(.))])
    colnames(EOC.Z) <- sub('.EOC', '', colnames(EOC.Z))
    EOC.Z <- EOC.Z[, samples_max_eoc$sample]
    ncol(EOC.Z) == nrow(samples_max_eoc)  
    sum(colnames(EOC.Z) %in% samples_max_eoc$sample)
    
    G.EOC <- as.matrix(rna_g_prim[, samples_max_eoc$sample])  
    W <- as.matrix(samples_max_eoc %>% 
                     tibble::column_to_rownames(var='sample') |>
                     dplyr::select(EOC, Fibroblast, Immune, Unknown))
    W.EOC <- W[ which(W[,1] >0.25) , grep('EOC', colnames(W)) ]
    
    Z.n <- diag.mul.r( EOC.Z, 1./ (G.EOC *W.EOC ) )
    
    # running DEG
    countdata = EOC.Z 
    dim(countdata)
    sum(is.na(countdata))
    
    # DEG
    nrow(samples_max_eoc) == ncol(countdata)
    
    dds <- DESeq2::DESeqDataSetFromMatrix(countData=round(countdata), 
                                          colData=samples_max_eoc, 
                                          design= ~ pheno)
    keep <- rowSums( counts(dds) ) > 100
    dds <- dds[keep,]
    
    sizeFactor <- ( G.EOC *W.EOC ) #normFactor(rna.g, samples_max_eoc$EOC)
    dds$sizeFactor = unlist(sizeFactor)
    
    dds <- DESeq2::DESeq(dds)
    resultsNames(dds)
    
    res <- DESeq2::results( dds )
    summary(res)
  }
}

# figure 3a: volcano plot
{
  geneInfo = geneID2name
  res$symbol = geneInfo$uniqName[match(sub( '[.].*$', '', rownames(res)), 
                                       sub( '[.].*$', '', geneInfo$geneID ))]
  res$geneType = geneInfo$geneType[match(sub( '[.].*$', '', rownames(res)), 
                                         sub( '[.].*$', '', geneInfo$geneID ))]
  resSig <- res[ which(res$padj < 0.01 ), ]
  dim(resSig)
  # Volcano plot
  vol.df <- res %>% 
    as.data.frame() %>%
    filter(!is.na(padj)) |>
    arrange(padj) %>%
    dplyr::rename(geneName = symbol) %>%
    mutate(thresh=ifelse(log2FoldChange > 1.5 & padj < .05, 'ctDNA ≥ 3%', 
                      ifelse(log2FoldChange < -1.5 & padj < .05, 'ctDNA < 3%', '')),
           gene = ifelse(abs(log2FoldChange) >= 1.5 & padj < .05, 
                         geneName, ''))
  
  sig.gene <- vol.df |>
    filter(padj < .05, abs(log2FoldChange) > 2,
           geneType=='protein_coding')
  
  fig3c <- vol.df %>% 
    #filter(log2FoldChange < 5)|>
    mutate(y = -log10(padj)) %>% 
    ggplot(., aes(x = log2FoldChange, y = y, label = geneName)) + 
    geom_point(aes(color=thresh))+
    geom_vline(xintercept = -1.5, linetype="dotdash", color='navy') +
    geom_vline(xintercept = 1.5, linetype="dotdash", color='navy')+
    geom_hline(yintercept = -log10(.05), linetype="dotdash", color='navy')+
    scale_colour_manual(values = c("ctDNA < 3%" = "dark blue", "ctDNA ≥ 3%"="maroon")) +
    labs(y='Log10 padj', color='Expressed in') + 
    theme_classic()+
    theme(title = element_text(size = 13, color = 'navy'),
          legend.position = c(.8,.8))
  
  column_ha = ComplexHeatmap::HeatmapAnnotation(TF = samples$TF,
                                                ctdna_lev = samples$pheno)
  heat_mat <- Z.n[,samples$sample] |> as.data.frame() |>
    tibble::rownames_to_column(var='geneID') |>
    filter(geneID %in% rownames(sig.gene)) |>
    left_join(geneInfo |>
                dplyr::select(geneID, uniqName)) |>
    dplyr::select(-geneID) |>
    tibble::column_to_rownames(var='uniqName') |>
    as.matrix()
  dim(heat_mat)
  heat_mat <- t(scale(t(heat_mat)))
  fig3b <- ComplexHeatmap::Heatmap(heat_mat, show_row_names = T, show_column_names = F,
                                   #column_split = col_data$pheno,
                                   #row_km = 4, 
                                  cluster_columns = F, 
                                  top_annotation = column_ha) 
          
  svg('results/figures/submission/Supp_3/fig3b.svg', width = 8, height = 5)
  fig3b
  dev.off()
  sig.gene <- vol.df %>% 
    filter(padj < .05,
           geneType=='protein_coding',
           abs(log2FoldChange) > 2)
  
  dim(sig.gene)
  heatmap_df <- data.frame(Z.n)[,samples$sample] %>% 
    filter(rownames(.)%in% rownames(sig.gene)) |>
    tibble::rownames_to_column(var='geneID') |>
    left_join(geneInfo %>%
                dplyr::select(geneID, symbol=uniqName)) |>
    dplyr::select(-geneID) |>
    tibble::column_to_rownames(var='symbol')
  
  
  #pheatmap::pheatmap(heatmap_df, scale = 'row', cluster_cols = F)
  scaled_mat = t(scale(t(heatmap_df)))
  col_fun = circlize::colorRamp2(c(min(scaled_mat), 0, max(scaled_mat)), c("#3399ff", "white", "#ff6666"))
  
  column_ha = ComplexHeatmap::HeatmapAnnotation(TF = samples$TF,
                                                ctdna_lev = samples$pheno)
  
  ComplexHeatmap::Heatmap(scaled_mat, 
                          cluster_columns = F, 
                          col = col_fun,
                          row_km =3,
                          top_annotation = column_ha,
                          show_column_names = F) 
}


# figure 3b: Top significant differently expressed genes



# heatmap
{
  sample_order <- samples |>
    dplyr::arrange(by = TF) |>
    dplyr::select(sample, ctdna_lev, TF)
  
  
  # Create a DGEList object
  dge <- DGEList(counts = countdata)
  
  # Perform TMM normalization
  dge <- calcNormFactors(dge)
  
  # Get normalized counts
  DEG_normalized <- cpm(dge, log = TRUE)
  
  DEG_ordered <- DEG_normalized[, sample_order$sample] |>
    as.data.frame() %>% 
    filter(rownames(.) %in% rownames(resSig))
  gplots::heatmap.2(as.matrix(t(DEG_ordered)),
                    trace = "none",  # Remove row/column names on the sides
                    margins = c(10, 10),  # Add margins for annotations
                    Colv = NA,  # Do not cluster columns
                    dendrogram = "none",  # Do not show dendrograms
                    #RowSideColors = sample_order$ctdna_lev,  # Color rows by group
                    key = TRUE,  # Add color key
                    keysize = 1,  # Set size of color key
                    key.title = "Groups",  # Set color key title
                    key.xlab = "Group",  # Set color key x-axis label
                    main = "Gene Expression Heatmap",  # Set heatmap title
                    labRow = FALSE  # Do not show row labels
  )
}

# heterogeneity scores
{
  # https://github.com/WangX-Lab/DEPTH
  # devtools::install_github("WangX-Lab/DEPTH")
  library(DEPTH)
  rna_exp <- countdata 
  sample_cat <- samples |>
    mutate(Identification = 'Tumor') |>
    dplyr::select(State=sample, Identification) 
  
  het_score_res <- DEPTH(rna_exp, sample_cat)
  
  # load configurations
  source("bin/GSECA/Scripts/config.R")
  
}

# by tissue
{
  ### processing
  {
    # EOC expression analysis
    # solid tumor only
    samples_ome <- rna_cohort |>
      mutate(pheno = ifelse(TF <= .07, 'med', 'high')) |>
      filter(tissue == 'Ome')
    {
      EOC.Z.ome <- rna_z_prim %>%
        as.data.frame() %>%
        dplyr::select(colnames(.)[sub('.EOC', '', colnames(.)) %in% samples_ome$sample]) 
      colnames(EOC.Z.ome) = sub('.EOC', '', colnames(EOC.Z.ome))
      ncol(EOC.Z.ome) == nrow(samples_ome)  
      sum(colnames(EOC.Z.ome) %in% samples_ome$sample)
      
      G.EOC.ome <- rna_g_prim %>%
        dplyr::select(colnames(.)[colnames(.) %in% samples_ome$sample]) 
      W.EOC.ome <- samples_ome$EOC
      Z.n <- diag.mul.r( EOC.Z.ome, 1./G.EOC.ome *W.EOC.ome )
      
      # filter gene count matrix for DEG
      countdata = EOC.Z.ome
      dim(countdata)
      sum(is.na(countdata))
      # filter gene count matrix for DEG
      rowSum = apply(countdata, 1, sum)
      rowVar = apply(countdata, 1, var)
      countdata <- countdata[rowSum > 100,]
      Z.r <- t(unit.ranks( t(countdata) )) # Z is the expression matrix (gene by cell)
      dim(Z.r)
      sum(is.na(Z.r))
      
      # DEG
      nrow(samples_ome) == ncol(countdata)
      rna.g <- rna_g_prim %>%
        dplyr::select(colnames(.)[colnames(.) %in% samples_ome$sample])
      
      deseq2Data <- DESeq2::DESeqDataSetFromMatrix(countData=round(countdata), 
                                                   colData=samples_ome, 
                                                   design= ~ factor(pheno))
      dds <- DESeq2::DESeq(deseq2Data)
      resultsNames(dds)
      
      sizeFactor <- normFactor(rna.g, samples_ome$EOC)
      dds$sizeFactor = sizeFactor
      #resLFC <- lfcShrink(dds, coef = 2)
      
      res <- DESeq2::results( dds )
      summary(res)
      }
    
  }
  
  # figure 3a: volcano plot
  {
    geneInfo = geneID2name
    res$symbol = geneInfo$uniqName[match(sub( '[.].*$', '', rownames(res)), 
                                         sub( '[.].*$', '', geneInfo$geneID ))]
    res$geneType = geneInfo$geneType[match(sub( '[.].*$', '', rownames(res)), 
                                           sub( '[.].*$', '', geneInfo$geneID ))]
    resSig <- res[ which(res$padj < 0.01 ), ]
    # Volcano plot
    vol.df <- res %>% 
      as.data.frame() %>%
      arrange(padj) %>%
      dplyr::rename(geneName = symbol) %>%
      mutate(thresh=ifelse(log2FoldChange > 2 & padj < .05 & geneType=='protein_coding', 'Up', 
                           ifelse(log2FoldChange < -2 & padj < .05 & geneType=='protein_coding', 'Down', 'Norm')),
             gene = ifelse(log2FoldChange > 2 & padj < .05 & geneType=='protein_coding', 
                           geneName, ''))
    
    vol.df %>% 
      filter(log2FoldChange<20)|>
      mutate(y = -log(padj)) %>% 
      ggplot(., aes(x = log2FoldChange, y = y, color=thresh, label = geneName)) + 
      geom_point(aes(colort=thresh))+
      #geom_label(aes(label=gene)) + ggrepel::geom_text_repel(size = 2) + 
      scale_colour_manual(values = c("Up"= "red", "Down"="blue",  "Norm"= "gray")) +
      labs(y='log2 padj', title = 'Vocalno Plot for DEG in all samples') + 
      theme_classic()+
      theme(title = element_text(size = 13, color = 'navy'))
    
    vol.df %>% 
      filter(padj < .01) %>%
      dplyr::arrange(desc(log2FoldChange)) #%>% dim()
  }
  
  
  # figure 3b: Top significant differently expressed genes
  
  
  # figure 3c: pathway enrichment analysis
  # pathway analysis
  
  {
    dat_path <- '~/mnt/storageBig8/work/nguyenma/projects/RNA/OmicsIntegration/data/temp/'
    pathways.hallmark <- readRDS(paste0(dat_path,'hallmark_pathway.rds'))
    
    pathways.kegg <- readRDS(paste0(dat_path,'kegg_pathway.rds'))
    
    res2 <- as.data.frame(res) %>% 
      dplyr::arrange(padj) %>%
      dplyr::select(symbol, stat) %>% 
      na.omit() %>% 
      distinct() %>% 
      group_by(symbol) %>% 
      summarize(stat=mean(stat))
    
    ranks <- tibble::deframe(res2)
    
    # hallmark database
    hallmark.res <- fgsea(pathways=pathways.hallmark, stats=ranks)
    
    hallmark.resTidy <- hallmark.res %>%
      as_tibble() %>%
      arrange(desc(NES))
    hallmark.resTidy %>%
      filter(padj < 0.05) %>%
      ggplot( aes(reorder(pathway, NES), NES)) +
      geom_col(aes(fill=NES<0), alpha=.7) +
      scale_fill_brewer(palette='Set1', direction = -1)+
      coord_flip() +
      labs(x="Pathway", y="Normalized Enrichment Score",
           fill='ctDNA TF > 7%',
           title="Hallmark pathways NES>1.5 from GSEA") + 
      theme_minimal()+
      theme(legend.position = 'bottom')
    
    #https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
    plotEnrichment(pathways.hallmark[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]],
                   ranks) + labs(title="HALLMARK_TNFA_SIGNALING_VIA_NFKB")
    plotEnrichment(pathways.hallmark[["HALLMARK_INFLAMMATORY_RESPONSE"]],
                   ranks) + labs(title="HALLMARK_INFLAMMATORY_RESPONSE")
    plotEnrichment(pathways.hallmark[["HALLMARK_E2F_TARGETS"]],
                   ranks) + labs(title="HALLMARK_E2F_TARGETS")
    
    
    # kegg database
    fgseaKEGGRes <- fgsea(pathways=pathways.kegg, stats=ranks)
    
    fgseaResKEGGTidy <- fgseaKEGGRes %>%
      as_tibble() %>%
      arrange(desc(NES))
    
    kegg <- ggplot(fgseaResKEGGTidy %>% 
                     as.data.frame()%>%
                     filter(padj<0.05)%>%
                     mutate(pathway = tomeder(sub('KEGG_','', pathway))), 
                   aes(reorder(pathway, NES), NES)) + 
      geom_col(aes(fill=padj<0.01)) +
      coord_flip() +labs(x="Pathway", y="Normalized Enrichment Score",
                         title="KEGG pathways NES from GSEA\n on Cancer Cells Population") + 
      theme_minimal(base_size = 10) + theme(title = element_text(size=12))
  }
  
  # progeny score from Daria
  {
    bulk.score <- read.delim('~/mnt/storageBig8/work/nguyenma/projects/RNA/OmicsIntegration/data/processed/Progeny_scores_ctDNA_x_RNA_bulk.tsv') |>
      dplyr::rename(sample=X) |>
      mutate(patient = sub('_.*', '', sample)) |>
      left_join(rna_cohort |>
                  dplyr::select(sample, tissue)) |>
      inner_join(ctdna_cohort |>
                  dplyr::select(patient, TF)) |>
      left_join(evol_df) |>
      mutate(ctdna_lev = ifelse(TF>.07, 'high', ifelse(TF > .01, 'med', 'med')),
             ctdna_cat = ifelse(TF > .07, 'high', 'med'))
    bulk.score |> 
      ggplot(aes(x = ctdna_cat, y = Androgen)) +
      geom_boxplot() +
      stat_compare_means()
    
    bulk.score |> 
      #filter(TF >= .03) |>
      ggscatter(x = 'clonality', y = 'JAK.STAT', add = 'reg.line', cor.coef = T, facet.by = 'tissue')
    
      ggplot(aes(x = TF, y = TNFa)) +
      geom_point() +
      geom_smooth(method = 'lm')
      
      
      eoc.score <- read.delim('~/mnt/storageBig8/work/nguyenma/projects/RNA/OmicsIntegration/data/processed/Progeny_scores_ctDNA_x_RNA_EOC.tsv') |>
        dplyr::rename(sample=X) |>
        mutate(patient = sub('_.*', '', sample)) |>
        left_join(rna_cohort |>
                    dplyr::select(sample, tissue, EOC, Fibroblast, Immune)) |>
        inner_join(ctdna_cohort |>
                     dplyr::select(patient, TF)) |>
        mutate(ctdna_lev = ifelse(TF>=.03, '≥ 3%', '< 3%'))
      
      eoc.score |> 
        dplyr::group_by(patient, ctdna_lev) |>
        dplyr::summarise(eoc=median(PI3K)) |>
        dplyr::ungroup() |>
        as.data.frame() |> #head(3)
        ggplot(aes(x = ctdna_lev, y = eoc)) +
        geom_boxplot(aes(color = ctdna_lev), show.legend = F) +
        geom_jitter(aes(color = ctdna_lev), show.legend = F) +
        scale_color_brewer(palette = 'Set1') +
        stat_compare_means(comparisons = comp2) +
        theme_pubr()
      
      eoc.score %>% 
        dplyr::slice_max(EOC, by = patient) |>
        #filter(TF >= .07) |>
        ggscatter(x = 'TF', y = 'PI3K', add = 'reg.line', cor.coef = T)
      
      ggplot(aes(x = TF, y = PI3K)) +
        geom_point() +
        geom_smooth(method = 'lm')
    
  }
}
### heatmap by gene signatures
Z.n[1:3,1:5]
head(resSig)
exp.mat <- Z.n %>%
  filter(rownames(.) %in% rownames(resSig))



