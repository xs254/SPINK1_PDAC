library(tximport)
library(DESeq2)
library(EnhancedVolcano)
library(ggplot2)
library(org.Hs.eg.db)
# Create the metadata for the KO and OE datasets
KO_metadata <- data.frame(row.names = c(sprintf('Scramble_%d',1:5),sprintf('KO_%d',1:5)),
                          file_name = c(sprintf('TM59Z7_%d_SCR_%d_quant.sf',1:5,1:5),
                                        sprintf('TM59Z7_%d_KO_%d_quant.sf',6:10,1:5)),
                          condition = factor(rep(c('Scramble_Control','Knockout'),each = 5),
                                             levels = c('Scramble_Control','Knockout'))) 
OE_metadata <- data.frame(row.names = c(sprintf('OE_Control_%d',1:5),sprintf('OE_%d',1:5)),
                          file_name = c(sprintf('TM59Z7_%d_CONT_%d_quant.sf',11:15,1:5),
                                        sprintf('TM59Z7_%d_OE_%d_quant.sf',16:20,1:5)),
                          condition = factor(rep(c('OE_Control','OE'),each = 5),
                                             levels = c('OE_Control','OE')))
merged_metadata <- rbind(KO_metadata,OE_metadata)
# Convert transcript counts to gene counts
setwd('C:/Users/xshan/OneDrive/Desktop/R code/DEseq_SPINK1_KO_v2/data/TM59Z7_quant_files')
tx2gene <- read.csv('C:/Users/xshan/OneDrive/Desktop/R code/tx2gene/tx2gene_gencode_v49_basic.csv',header = TRUE)
txi_KO <- tximport(KO_metadata$file_name,type = 'salmon', tx2gene = tx2gene)
txi_OE <- tximport(OE_metadata$file_name,type = 'salmon', tx2gene = tx2gene)
txi_merged <- tximport(merged_metadata$file_name,type = 'salmon', tx2gene = tx2gene)
# Pre-processing and QC of data
setwd('C:/Users/xshan/OneDrive/Desktop/R code/DEseq_SPINK1_KO_v2/results')
dds_KO <- DESeqDataSetFromTximport(txi_KO,colData = KO_metadata,design = ~ condition)
dds_KO$condition <- relevel(dds_KO$condition,ref = 'Scramble_Control')
rownames(dds_KO) <- gsub('\\..*','',rownames(dds_KO))
dds_OE <- DESeqDataSetFromTximport(txi_OE,colData = OE_metadata,design = ~ condition)
dds_OE$condition <- relevel(dds_OE$condition,ref = 'OE_Control')
rownames(dds_OE) <- gsub('\\..*','',rownames(dds_OE))
# Merged dds is only for plotting PCR plots. DEseq will be run separately on dds_KO and dds_OE
dds_merged <- DESeqDataSetFromTximport(txi_merged,colData = merged_metadata,design = ~ condition) 
#PCA for QC
dds_KO_norm <- normTransform(dds_KO)
dds_OE_norm <- normTransform(dds_OE)
dds_merged_norm <- normTransform(dds_merged)
png('KO_QC_by_PCA.png',width = 2000,height = 1500,units = 'px',res = 300)
plotPCA(dds_KO_norm, intgroup=c('condition'))
dev.off()
png('OE_QC_by_PCA.png',width = 2000,height = 1500,units = 'px',res = 300)
plotPCA(dds_OE_norm, intgroup=c('condition'))
dev.off()
png('OE_KO_merged_QC_by_PCA.png',width = 2000,height = 1500,units = 'px',res = 300)
plotPCA(dds_merged_norm, intgroup=c('condition'))
dev.off()
# Run DESeq
dds_KO <- DESeq(dds_KO)
res_KO <- results(dds_KO)
resultsNames(dds_KO) #Check the names of the estimated effects
res_KO <- lfcShrink(dds_KO,coef = 2,type='apeglm') #coef = 2 means the target coef to shrink is the 2nd in resultNames.
dds_OE <- DESeq(dds_OE)
res_OE <- results(dds_OE)
resultsNames(dds_OE) #Check the names of the estimated effects
res_OE <- lfcShrink(dds_OE,coef = 2,type='apeglm') #coef = 2 means the target coef to shrink is the 2nd in resultNames.
#Export the DEseq results
KO_gene_symbols <- mapIds(org.Hs.eg.db,keys = rownames(res_KO),column = 'SYMBOL',
                       keytype = 'ENSEMBL',multiVals = 'first')
res_KO_df <- as.data.frame(res_KO)
res_KO_df <- cbind(KO_gene_symbols,res_KO_df)
colnames(res_KO_df)[1] <- 'gene_symbols' 
OE_gene_symbols <- mapIds(org.Hs.eg.db,keys = rownames(res_OE),column = 'SYMBOL',
                       keytype = 'ENSEMBL',multiVals = 'first')
res_OE_df <- as.data.frame(res_OE)
res_OE_df <- cbind(OE_gene_symbols,res_OE_df)
colnames(res_OE_df)[1] <- 'gene_symbols' 
write.csv(res_KO_df,'KO_DEseq_results.csv')
write.csv(res_OE_df,'OE_DEseq_results.csv')
# Visualize DESeq results
png('KO_enhancedVolcano_plot.png',width = 2500,height = 4000,res = 300,units = 'px')
EnhancedVolcano(res_KO,lab=KO_gene_symbols,x='log2FoldChange',y='padj',
                pCutoff = 0.05, FCcutoff = 1)
dev.off()
png('OE_enhancedVolcano_plot.png',width = 2500,height = 4000,res = 300,units = 'px')
EnhancedVolcano(res_OE,lab=OE_gene_symbols,x='log2FoldChange',y='padj',
                pCutoff = 0.05, FCcutoff = 1)
dev.off()
#Export normalized count matrix from dds for reference
KO_cts_matrix_normalized <- counts(dds_KO,normalized = TRUE)
KO_cts_matrix_normalized  <- as.data.frame(KO_cts_matrix_normalized)
KO_cts_matrix_normalized <- cbind(KO_gene_symbols,KO_cts_matrix_normalized)
colnames(KO_cts_matrix_normalized)[1] <- 'gene_symbols'
write.csv(KO_cts_matrix_normalized,file = 'KO_normalized_counts_matrix.csv')
OE_cts_matrix_normalized <- counts(dds_OE,normalized = TRUE)
OE_cts_matrix_normalized <- as.data.frame(OE_cts_matrix_normalized)
OE_cts_matrix_normalized <- cbind(OE_gene_symbols,OE_cts_matrix_normalized)
colnames(OE_cts_matrix_normalized)[1] <- 'gene_symbols'
write.csv(OE_cts_matrix_normalized,file = 'OE_normalized_counts_matrix.csv')