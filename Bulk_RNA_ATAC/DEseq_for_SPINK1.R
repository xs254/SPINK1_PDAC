#Load libraries
library(tximport)
library(DESeq2)
library(DOSE)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(ggridges)
library(patchwork)
library(fs)
library(RColorBrewer)
library(tidyverse)
library(scales)
library(decoupleR)
library(aracne.networks)
library(viper)
#Specify the folder you want all results to be exportd to.
quants_dir <- 'C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\DEseq_SPINK1_KO_v2\\salmon_quants\\'
output_dir <- 'C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\DEseq_SPINK1_KO_v2\\Results\\'
meta_data <- data.frame(row.names = c('Scramble_1','Scramble_2','KO_1','KO_2'),
                               condition = factor(c('Scramble_Control','Scramble_Control','Knock_Out','Knock_Out')))
dir.create(output_dir)
# Read in quant.sf files and convert to gene level abundance
setwd(quants_dir)
tx2gene <- read.csv('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\tx2gene\\tx2gene_gencode_v45_basic.csv',header = TRUE)
txi <- tximport(c('SCR1_257_320_S11_L002_R1_001_quant.sf','SCR2_245_332_S12_L002_R1_001_quant.sf',
                      '2G-6-1_221_356_S14_L002_R1_001_quant.sf','2G-6-2_209_368_S15_L002_R1_001_quant.sf'), 
                    type = 'salmon', tx2gene = tx2gene)
setwd(output_dir)
txi_cts_matrix <- txi$counts
write.csv(txi_cts_matrix,file = 'estimated_gene_counts.csv')
#Step 1: Pre-processing and QC of data
dds <- DESeqDataSetFromTximport(txi,colData = meta_data,design = ~ condition)
dds$condition <- relevel(dds$condition,ref = 'Scramble_Control')
rownames(dds) <- gsub('\\..*','',rownames(dds))
#PCA for QC
dds_norm <- normTransform(dds)
bmp(file = 'QC by PCA.bmp',width = 512,height = 512,units = "px")
plotPCA(dds_norm, intgroup=c("condition"))
dev.off()
#Export normalized count matrix from dds for reference
dds <- estimateSizeFactors(dds)
cts_matrix_normalized <- counts(dds,normalized = TRUE)
cts_matrix_normalized  <- as.data.frame(cts_matrix_normalized)
gene_symbols <- mapIds(org.Hs.eg.db,keys = rownames(dds),column = 'SYMBOL',
                       keytype = 'ENSEMBL',multiVals = 'first')
cts_matrix_normalized <- cbind(gene_symbols,cts_matrix_normalized)
write.csv(cts_matrix_normalized,file = 'normalized_counts_matrix.csv')
#Setp 2: Run DEseq, shrink the LFC and export the results
dds <- DESeq(dds)
res <- results(dds)
resultsNames(dds) #Check the names of the estimated effects
res <- lfcShrink(dds,coef = 2,type='apeglm') #coef = 2 means the target coef to shrink is the 2nd in resultNames.
#Export the DEseq results
gene_symbols <- mapIds(org.Hs.eg.db,keys = rownames(res),column = 'SYMBOL',
                       keytype = 'ENSEMBL',multiVals = 'first')
res_df <- as.data.frame(res)
res_df <- cbind(gene_symbols,res_df)
write.csv(res_df,file = 'DEseq_results.csv')
#Visualize DEseq results
bmp(file = 'EnhancedVolcano plot.bmp',width = 2048,height = 1024,units = "px")
EnhancedVolcano(res,lab=gene_symbols,x='log2FoldChange',y='padj',
                pCutoff = 0.01, FCcutoff = 1)
dev.off()
#Visualize DEseq results focusing on EMT markers
EMT_markers <- read.csv('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\DEseq_SPINK1_KO_v2\\Total EMT gene list.csv',header = FALSE)
EMT_markers <- EMT_markers[1:50,c(-1,-8)]
colnames(EMT_markers) <- c('Gene_Symbols','Gene_Length','Type','Chrosome','Full_Name','Description')
#Plot 1 volcano plot with EMT markers labelled
bmp(file = 'EnhancedVolcano plot EMT.bmp',width = 2048,height = 1024,units = "px")
EnhancedVolcano(res,selectLab = EMT_markers$Gene_Symbols,lab=gene_symbols,x='log2FoldChange',y='padj',
                pCutoff = 0.01, FCcutoff = 1)
dev.off()
#Plot 2 heatmap of valid EMT markers
res_EMT_markers <- res_df[res_df$gene_symbols %in% EMT_markers$Gene_Symbols,]
res_EMT_markers_valid <- na.omit(res_EMT_markers)
res_EMT_markers_valid <- res_EMT_markers_valid[order(res_EMT_markers_valid$log2FoldChange,decreasing = TRUE),]
EMT_DEGs <- res_EMT_markers_valid$gene_symbols[res_EMT_markers_valid$padj<=0.05]
cts_matrix_EMT_markers_valid <- cts_matrix_normalized[rownames(res_EMT_markers_valid),]
rownames(cts_matrix_EMT_markers_valid) <- cts_matrix_EMT_markers_valid$gene_symbols
cts_matrix_EMT_markers_valid <- cts_matrix_EMT_markers_valid[,-1]
cts_matrix_EMT_markers_valid_scaled <- as.data.frame(t(scale(t(cts_matrix_EMT_markers_valid))))
heatmap_input <- cbind(rownames(cts_matrix_EMT_markers_valid_scaled),cts_matrix_EMT_markers_valid_scaled)
colnames(heatmap_input) <- c('gene_symbols','Scramble Control 1','Scramble Control 2','SPINK1 KO 1','SPINK1 KO 2')
heatmap_input[EMT_DEGs,'gene_symbols'] <- paste0('*',heatmap_input[EMT_DEGs,'gene_symbols'])
heatmap_input <- tidyr::pivot_longer(heatmap_input,cols = !gene_symbols,names_to = 'condition',values_to = 'Z_score')
heatmap_input$gene_symbols <- factor(heatmap_input$gene_symbols,levels = unique(heatmap_input$gene_symbols))
heatmap_input$condition <- factor(heatmap_input$condition, levels = unique(heatmap_input$condition))
bmp(file = 'Heatmap for EMT markers.bmp',width = 800,height = 1000,units = "px",res = 150)
ggplot(heatmap_input,aes(x = condition, y = gene_symbols, fill = Z_score))+geom_tile()+
  scale_fill_distiller(palette = 'Spectral',limits = c(-1.5,1.5), oob = squish)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  xlab('')+ylab('')
dev.off()
# Plot 3 Overall statistic EMT markers
nEMTs_invalid <- sum(is.na(res_EMT_markers$padj))
nEMTs_notSig <- sum(res_EMT_markers_valid$padj>0.05)
nEMTs_upSig <- sum(res_EMT_markers_valid$padj<=0.05 & res_EMT_markers_valid$log2FoldChange>0)
nEMTs_downSig <- sum(res_EMT_markers_valid$padj<=0.05 & res_EMT_markers_valid$log2FoldChange<0)
EMT_markers_statistics <- data.frame(outcome = c('Invalid','Not significant','Significantlly upregulated','Significantly downregulated'),
                                     frequency = c(nEMTs_invalid,nEMTs_notSig,nEMTs_upSig,nEMTs_downSig))
EMT_markers_statistics['percentage'] <- EMT_markers_statistics$frequency/sum(EMT_markers_statistics$frequency)
bmp(file = 'Piechart for overall stastics of EMT markers.bmp',width = 750,height = 500, units = 'px',res = 150)
ggplot(EMT_markers_statistics, aes(x="", y=percentage, fill=outcome))+
  geom_bar(stat='identity', width=1, color="white")+xlab('')+ylab('')+
  coord_polar("y", start=0)+theme_minimal()
dev.off()
#Step 3: Prepare ranked genelist for GSEA
padj_nozero <- res$padj
padj_nozero[padj_nozero==0] <- min(res$padj[res$padj>0])/10 #To avoid InF due to log0.
ranked_genelist <- res$log2FoldChange
names(ranked_genelist) <- rownames(res)
ranked_genelist <- sort(ranked_genelist,decreasing=TRUE)
#Step 4a: Perform GSEA with GO
gse_from_GO <- gseGO(geneLis = ranked_genelist, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             minGSSize = 10, 
             maxGSSize = 200, 
             pvalueCutoff = 0.01, 
             verbose = TRUE, 
             OrgDb = 'org.Hs.eg.db', 
             pAdjustMethod = 'BH')
# Visualize GO results
gse_from_GO <- pairwise_termsim(gse_from_GO,showCategory = dim(gse_from_GO)[1])
set.seed(1) # For reproducibility of emapplot
bmp(file = 'enrichment map for GO.bmp', width = 2560, height = 1440, units = "px")
emapplot(gse_from_GO,showCategory = dim(gse_from_GO)[1],
         cex.params = list(category_node = 0.5,category_label = 0.75),
         group_category = T, group_legend = T, node_label = 'category')
dev.off()
bmp(file = 'dotplot for GO.bmp', width = 1024, height = 1440, units = "px")
dotplot(gse_from_GO,showCategory = 100,font.size = 12,label_format = 200)
dev.off()
#Step 4b: Perform GSEA with KEGG
ranked_genelist_uniprot <- ranked_genelist
uniprot_symbols <- mapIds(org.Hs.eg.db,keys = names(ranked_genelist),column = 'UNIPROT',
                          keytype = 'ENSEMBL',multivals = 'first')
names(ranked_genelist_uniprot) <- uniprot_symbols
gse_from_KEGG <- gseKEGG(geneLis = ranked_genelist_uniprot,
                         organism = 'hsa',
                         keyType = 'uniprot',
                         pvalueCutoff = 0.05, 
                         verbose = TRUE,
                         pAdjustMethod = 'BH')
# Visualize KEGG results
gse_from_KEGG <- pairwise_termsim(gse_from_KEGG,showCategory = dim(gse_from_KEGG)[1])
gse_from_KEGG_readable <- setReadable(gse_from_KEGG,org.Hs.eg.db,keyType = 'UNIPROT')
set.seed(1) # For reproducibility of emapplot
bmp(file = 'enrichment map for KEGG.bmp', width = 2560, height = 1440, units = "px")
emapplot(gse_from_KEGG,showCategory = dim(gse_from_KEGG)[1],
         cex.params = list(category_node = 0.5,category_label = 0.75),
         group_category = T, group_legend = T, node_label = 'category')
dev.off()
bmp(file = 'dotplot for KEGG.bmp', width = 1024, height = 1440, units = "px")
dotplot(gse_from_KEGG,showCategory = 100,font.size = 12,label_format = 200)
dev.off()
# Leading edge analysis
leading_edge_dir = paste(output_dir,'leading edge analysis\\',sep = '')
dir.create(leading_edge_dir)
for(to_plot in 1:dim(gse_from_KEGG)[1]){
  core_enrichment_genes <- gse_from_KEGG@result$core_enrichment[to_plot]
  core_enrichment_genes <- unlist(strsplit(core_enrichment_genes,'/'))
  core_enrichment_genes_readbale <- gse_from_KEGG_readable@result$core_enrichment[to_plot]
  core_enrichment_genes_readbale <- unlist(strsplit(core_enrichment_genes_readbale,'/'))
  
  core_enrichment_lfc <- data.frame(genes = core_enrichment_genes_readbale,
    lfc = ranked_genelist_uniprot[core_enrichment_genes])
  core_enrichment_lfc$genes <- make.unique(core_enrichment_lfc$genes,sep = '_')
  core_enrichment_lfc$genes <- factor(core_enrichment_lfc$genes,
                                      levels = core_enrichment_lfc$genes[order(core_enrichment_lfc$lfc,decreasing = TRUE)])
  
  p1 <- gseaplot(gse_from_KEGG,geneSetID = to_plot,by = 'runningScore',
           title = sprintf('%s, q = %.2g',gse_from_KEGG@result$Description[to_plot],gse_from_KEGG@result$qvalue[to_plot]))
  p1 <- p1+labs(subtitle = gse_from_KEGG@result$leading_edge[to_plot])+
    theme(plot.title = element_text(hjust = 0.5,size = 15),plot.subtitle = element_text(hjust = 0.5,size = 12),
                 axis.title = element_text(size = 15))
  p2 <- ggplot(core_enrichment_lfc,aes(x = genes, y = lfc))+
    geom_col()+xlab('Genes')+ylab('LFC')+
    theme(axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1,size = 12),
          axis.title = element_text(size = 15),
          plot.title = element_text(size = 15))
  
  p2 <- p2+labs(title = 'log2fold change of core enriched genes')+theme(plot.title = element_text(hjust = 0.5,size = 15))
  
  filename_bmp <- path_sanitize(gse_from_KEGG@result$Description[to_plot],replacement = ' ')
  filename_bmp <- paste(leading_edge_dir,filename_bmp,'.bmp',sep='')
  bmp(filename = filename_bmp,
      width = 1440,height = 800, unit = 'px')
  print(p1/p2)
  dev.off()
}
# Transcription factor activity inference
cts_for_inference <- read.csv('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\DEseq_SPINK1_KO_v2\\Results\\normalized_counts_matrix.csv')
cts_for_inference <- na.omit(cts_for_inference)
cts_for_inference <- cts_for_inference[!duplicated(cts_for_inference$gene_symbols),]
row.names(cts_for_inference) <- cts_for_inference$gene_symbols
cts_for_inference <- cts_for_inference[,-c(1,2)]
cts_for_inference <- log1p(cts_for_inference)
cts_for_inference <- as.matrix(cts_for_inference)
TF_net <- get_collectri(organism='human', split_complexes=FALSE)
TF_inferred <- run_ulm(mat=cts_for_inference, net=TF_net, .source='source', .target='target',.mor='mor', minsize = 5)
TF_inferred_wide <- pivot_wider(TF_inferred,id_cols = source,names_from = condition,values_from = score)
# Visualize results from TF inference
#TF_statistics <- TF_inferred %>% group_by(source) %>% summarise(std = sd(score))
#TF_statistics <- TF_variability[order(TF_variability$std,decreasing = TRUE),]
n_TFs <- length(TF_inferred_wide$source)
TF_statistics <- data.frame(p = rep(NA,n_TFs),p_adjusted = rep(NA,n_TFs),change = rep(NA,n_TFs),row.names = TF_inferred_wide$source)
for(n in 1:n_TFs)
{
  to_analyze <- as.matrix(TF_inferred_wide[n,2:5])
  test_res <- t.test(to_analyze[1:2],to_analyze[3:4]) 
  TF_statistics$p[n] <- test_res$p.value
  TF_statistics$change[n] <- mean(to_analyze[1:2])-mean(to_analyze[3:4])
}
TF_statistics$p_adjusted <- p.adjust(TF_statistics$p,'BH')
TF_statistics <- TF_statistics[order(TF_statistics$p_adjusted,decreasing = FALSE),]
write.csv(TF_statistics,'Differential_TF_statstics_inferred_from_collectri.csv')
TF_to_plot <- c(TF_variability$source[1:50])
heatmap_input <- TF_inferred[TF_inferred$source %in% TF_to_plot,]
heatmap_input <- as.data.frame(pivot_wider(heatmap_input,id_cols = 'condition',names_from = 'source',values_from = 'score'))
heatmap_input[,-1] <- as.data.frame(scale(heatmap_input[,-1]))
heatmap_input <- pivot_longer(heatmap_input,cols = !condition,names_to = 'Source',values_to = 'Z_score')
heatmap_input$condition <- factor(heatmap_input$condition,levels = c('Scramble_1','Scramble_2','KO_1','KO_2'))
heatmap_input$Source <- factor(heatmap_input$Source,levels = TF_to_plot)
ggplot(heatmap_input,aes(x = condition, y = Source, fill = Z_score))+geom_tile()+
  scale_fill_distiller(palette = 'Spectral',limits = c(-1.5,1.5), oob = squish)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  xlab('')+ylab('')
#
CDH1_regulator <- TF_net[TF_net$target=='CDH1',]
CDH2_regulator <- TF_net[TF_net$target=='CDH2',]
#
omni_path <- OmnipathR::import_transcriptional_interactions(resources = c('DoRothEA'),dorothea_levels = c('A','B'),references_by_resource = FALSE)
CDH1_omni_path <- omni_path[omni_path$target_genesymbol=='CDH1',]
#
cts_for_inference <- read.csv('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\DEseq_SPINK1_KO_v2\\Results\\normalized_counts_matrix.csv')
cts_for_inference <- na.omit(cts_for_inference)
cts_for_inference <- cts_for_inference[!duplicated(cts_for_inference$gene_symbols),]
row.names(cts_for_inference) <- mapIds(org.Hs.eg.db,keys = cts_for_inference$gene_symbols,column = 'ENTREZID',
                                       keytype = 'SYMBOL',multiVals = 'first')
cts_for_inference <- cts_for_inference[,-c(1,2)]
cts_for_inference <- as.matrix(cts_for_inference)
signature <- rowTtest(cts_for_inference[,3:4],cts_for_inference[,1:2])
signature <- (qnorm(signature$p.value/2, lower.tail = FALSE)*sign(signature$statistic))[, 1]
null_model <- ttestNull(cts_for_inference[,3:4],cts_for_inference[,1:2],per = 1000,repos = TRUE)
mrs <- msviper(signature,regulonpaad,null_model,verbose = TRUE)
mrs_es <- data.frame(row.names = names(mrs$es$nes),
                        symbol = mapIds(org.Hs.eg.db,keys = names(mrs$es$nes),column = 'SYMBOL',keytype = 'ENTREZID',multiVals = 'first'),
                        size = mrs$es$size,
                        nes = mrs$es$nes,
                        p_value = mrs$es$p.value)
mrs_es <- mrs_es[order(mrs_es$nes,decreasing = TRUE),]
mrs_es['p_adjusted'] <- p.adjust(mrs_es$p_value,'BH')
write.regulon(regulonpaad,file = 'paad_regulon.csv',header = FALSE,sep = ',')
# 