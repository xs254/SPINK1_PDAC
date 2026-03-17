library(ggplot2)
library(scales)
setwd('C:/Users/xshan/OneDrive/Desktop/R code/DEseq_SPINK1_KO_v2')
# Read unfiltered results from DESeq
res_KO <- read.csv('results/KO_DEseq_results.csv',row.names = 1)
res_OE <- read.csv('results/OE_DEseq_results.csv',row.names = 1)
# Correct the gene symbol for ENSG00000269226 from TMSB15B to TMSB15C
res_OE['ENSG00000269226','gene_symbols'] <- 'TMSB15C'
res_KO['ENSG00000269226','gene_symbols'] <- 'TMSB15C'
# Filter out significant genes by p_adj<=0.05
res_KO_filtered <- na.omit(res_KO)
res_OE_filtered <- na.omit(res_OE)
res_KO_filtered <- res_KO_filtered[res_KO_filtered$padj<=0.05,]
res_OE_filtered <- res_OE_filtered[res_OE_filtered$padj<=0.05,]

# Identify consensually regulated genes in both OE and KO
OE_KO_DEGs <- data.frame(gene_symbols = union(res_OE_filtered$gene_symbols,res_KO_filtered$gene_symbols),
                         consensus = NA)
for(i_gene in 1:nrow(OE_KO_DEGs))
{
  target_gene <- OE_KO_DEGs$gene_symbols[i_gene]
  is_in_KO <- target_gene %in% res_KO_filtered$gene_symbols
  is_in_OE <- target_gene %in% res_OE_filtered$gene_symbols
  if(is_in_KO & is_in_OE)
  {
    log2FC_in_KO <- res_KO_filtered$log2FoldChange[res_KO_filtered$gene_symbols==target_gene]
    log2FC_in_OE <- res_OE_filtered$log2FoldChange[res_OE_filtered$gene_symbols==target_gene]
    if(log2FC_in_KO<0 & log2FC_in_OE>0)
    {
      OE_KO_DEGs$consensus[i_gene] <- 'Consistent promotion'
    }else if(log2FC_in_KO>0 & log2FC_in_OE<0)
    {
      OE_KO_DEGs$consensus[i_gene] <- 'Consistent suppression'
    }else
    {
      OE_KO_DEGs$consensus[i_gene] <- 'Contradictory'
    }
    
  }else if(is_in_KO)
  {
    OE_KO_DEGs$consensus[i_gene] <- 'Only in KO'
  }else if(is_in_OE)
  {
    OE_KO_DEGs$consensus[i_gene] <- 'Only in OE'
  }
}
OE_KO_DEGs$consensus <- factor(OE_KO_DEGs$consensus,levels = c('Consistent promotion','Consistent suppression',
                                                               'Only in OE','Only in KO','Contradictory'))
OE_KO_DEGs_consensus_summary <- data.frame(types = names(summary(OE_KO_DEGs$consensus)),
                                counts = summary(OE_KO_DEGs$consensus))
OE_KO_DEGs_consensus_summary$types <- factor(OE_KO_DEGs_consensus_summary$types,
                                             levels = c('Consistent promotion','Consistent suppression',
                                                                     'Only in OE','Only in KO','Contradictory'))
png('results/Summary of consensus between KO and OE RNA seq.png',width = 3000,height = 2000,res = 300)
ggplot(OE_KO_DEGs_consensus_summary, aes(x='',y = counts,fill = types))+
  geom_bar(stat = 'identity',width=1) +coord_polar('y',start=0)+
  scale_fill_manual(values = c('salmon','dodgerblue2','gold','palegreen','grey'),name = NULL)+
  theme_void()+
  theme(legend.text = element_text(size = 20))
dev.off()
write.csv(OE_KO_DEGs,'results/consensus_of_DEGs_from_KO_and_OE.csv',row.names = FALSE)
# Integrate with results from ATACseq
ATAC_diffbind_integrated <- read.csv('C:/Users/xshan/OneDrive/Desktop/R code/ATACseq_SPINK1_KO/Results/diffbind_integrated.csv')
OE_KO_DEGs <- read.csv('results/consensus_of_DEGs_from_KO_and_OE.csv')
OE_KO_DEGs_consistent <- OE_KO_DEGs[OE_KO_DEGs$consensus %in% c('Consistent promotion','Consistent suppression'),]
ATAC_OE_KO_DEGs <- data.frame(gene_symbols = intersect(OE_KO_DEGs_consistent$gene_symbols,ATAC_diffbind_integrated$gene_symbols),
                              consensus = NA)
for(i_gene in 1:nrow(ATAC_OE_KO_DEGs))
{
  target_gene <- ATAC_OE_KO_DEGs$gene_symbols[i_gene]
  ATAC_consensus <- ATAC_diffbind_integrated$integrated_accessibility[ATAC_diffbind_integrated$gene_symbols==target_gene]
  OE_KO_consensus <- OE_KO_DEGs_consistent$consensus[OE_KO_DEGs_consistent$gene_symbols==target_gene]
  if(ATAC_consensus==0)
  {
    ATAC_OE_KO_DEGs$consensus[i_gene] <- 'ATAC inconclusive'
  }else
  {
    if(ATAC_consensus==1 & OE_KO_consensus=='Consistent suppression')
    {
      ATAC_OE_KO_DEGs$consensus[i_gene] <- 'Consistent epigenetic suppression'
    }else if(ATAC_consensus==-1 & OE_KO_consensus=='Consistent promotion')
    {
      ATAC_OE_KO_DEGs$consensus[i_gene] <- 'Consistent epigenetic promotion'
    }else
    {
      ATAC_OE_KO_DEGs$consensus[i_gene] <- 'ATAC RNA inconsistent'
    }   
  }
}
ATAC_OE_KO_DEGs$consensus <- factor(ATAC_OE_KO_DEGs$consensus)
ATAC_OE_KO_DEGs_cosensus_summary <- data.frame(types = names(summary(ATAC_OE_KO_DEGs$consensus)),
                                               counts = summary(ATAC_OE_KO_DEGs$consensus))
ATAC_OE_KO_DEGs_cosensus_summary$types <- factor(ATAC_OE_KO_DEGs_cosensus_summary$types,
                                                 levels = c('Consistent epigenetic promotion','Consistent epigenetic suppression',
                                                            'ATAC RNA inconsistent','ATAC inconclusive'))
png('results/Summary of consensus between ATAC and RNA seq.png',width = 3500,height = 2000,res = 300)
ggplot(ATAC_OE_KO_DEGs_cosensus_summary, aes(x='',y = counts,fill = types))+
  geom_bar(stat = 'identity',width=1) +coord_polar('y',start=0)+
  scale_fill_manual(values = c('salmon','dodgerblue2','palegreen','grey'),name = NULL)+
  theme_void()+
  theme(legend.text = element_text(size = 20))
dev.off()
writeLines(ATAC_OE_KO_DEGs$gene_symbols[ATAC_OE_KO_DEGs$consensus=='Consistent epigenetic suppression'],
           'results/Consistent_epigenetic_supression_genes_by_SPINK1.txt')
writeLines(ATAC_OE_KO_DEGs$gene_symbols[ATAC_OE_KO_DEGs$consensus=='Consistent epigenetic promotion'],
           'results/Consistent_epigenetic_promotion_genes_by_SPINK1.txt')
# Changes in EMT related genes
EMT_markers <- read.csv('Total EMT gene list.csv',header = FALSE)
EMT_markers <- EMT_markers[1:50,c(-1,-8)]
colnames(EMT_markers) <- c('gene_symbols','gene_length','type','chrosome','full_name','description')
epigenetic_consistent_EMTs <- intersect(EMT_markers$gene_symbols,
                                        ATAC_OE_KO_DEGs$gene_symbols[ATAC_OE_KO_DEGs$consensus %in% c('Consistent epigenetic suppression','Consistent epigenetic promotion')])
res_KO_EMT_markers <- res_KO[res_KO$gene_symbols %in% EMT_markers$gene_symbols,]
res_OE_EMT_markers <- res_OE[res_OE$gene_symbols %in% EMT_markers$gene_symbols,]
na_genes <- union(rownames(res_KO_EMT_markers)[is.na(res_KO_EMT_markers$padj)],
                  rownames(res_OE_EMT_markers)[is.na(res_OE_EMT_markers$padj)])
res_KO_EMT_markers <- res_KO_EMT_markers[!(rownames(res_KO_EMT_markers) %in% na_genes),]
res_OE_EMT_markers <- res_OE_EMT_markers[!(rownames(res_OE_EMT_markers) %in% na_genes),]

res_KO_EMT_markers <- res_KO_EMT_markers[order(res_KO_EMT_markers$log2FoldChange,decreasing = TRUE),]
res_OE_EMT_markers <- res_OE_EMT_markers[rownames(res_KO_EMT_markers),]

KO_EMT_DEGs <- res_KO_EMT_markers$gene_symbols[res_KO_EMT_markers$padj<=0.05]
OE_EMT_DEGs <- res_OE_EMT_markers$gene_symbols[res_OE_EMT_markers$padj<=0.05]

KO_EMT_cts <- read.csv('results/KO_normalized_counts_matrix.csv',row.names = 1)
KO_EMT_cts <- KO_EMT_cts[rownames(res_KO_EMT_markers),]
OE_EMT_cts <- read.csv('results/OE_normalized_counts_matrix.csv',row.names = 1)
OE_EMT_cts <- OE_EMT_cts[rownames(res_OE_EMT_markers),]

rownames(KO_EMT_cts) <- KO_EMT_cts$gene_symbols
KO_EMT_cts <- KO_EMT_cts[,-1]
KO_EMT_scaled <- as.data.frame(t(scale(t(KO_EMT_cts))))
rownames(OE_EMT_cts) <- OE_EMT_cts$gene_symbols
OE_EMT_cts <- OE_EMT_cts[,-1]
OE_EMT_scaled <- as.data.frame(t(scale(t(OE_EMT_cts))))

KO_heatmap_input <- cbind(rownames(KO_EMT_scaled),KO_EMT_scaled)
colnames(KO_heatmap_input)[1] <- 'gene_symbols'
KO_heatmap_input[KO_EMT_DEGs,'gene_symbols'] <- paste0('*',KO_heatmap_input[KO_EMT_DEGs,'gene_symbols'])
KO_heatmap_input <- tidyr::pivot_longer(KO_heatmap_input,cols = !gene_symbols,
                                        names_to = 'condition',values_to = 'Z_score')
KO_heatmap_input$gene_symbols <- factor(KO_heatmap_input$gene_symbols,levels = unique(KO_heatmap_input$gene_symbols))
KO_heatmap_input$condition <- factor(KO_heatmap_input$condition,levels = unique(KO_heatmap_input$condition))

OE_heatmap_input <- cbind(rownames(OE_EMT_scaled),OE_EMT_scaled)
colnames(OE_heatmap_input)[1] <- 'gene_symbols'
OE_heatmap_input[OE_EMT_DEGs,'gene_symbols'] <- paste0('*',OE_heatmap_input[OE_EMT_DEGs,'gene_symbols'])
OE_heatmap_input <- tidyr::pivot_longer(OE_heatmap_input,cols = !gene_symbols,
                                        names_to = 'condition',values_to = 'Z_score')
OE_heatmap_input$gene_symbols <- factor(OE_heatmap_input$gene_symbols,levels = unique(OE_heatmap_input$gene_symbols))
OE_heatmap_input$condition <- factor(OE_heatmap_input$condition,levels = unique(OE_heatmap_input$condition))
# Use red color to highlight genes that are in the ATAC-RNA consistent group.
y_tick_colors <- rep('black',nrow(KO_EMT_scaled))
names(y_tick_colors) <- rownames(KO_EMT_scaled)
y_tick_colors[epigenetic_consistent_EMTs] <- 'red'

png(file = 'results/KO_heatmap_for_EMT_markers.png',width = 1600,height = 2000,units = "px",res = 300)
ggplot(KO_heatmap_input,aes(x = condition, y = gene_symbols, fill = Z_score))+geom_tile()+
  scale_fill_distiller(palette = 'Spectral',limits = c(-1.5,1.5), oob = squish)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.text.y=ggtext::element_markdown(color=y_tick_colors))+
  xlab('')+ylab('')
dev.off()
png(file = 'results/OE_heatmap_for_EMT_markers.png',width = 1600,height = 2000,units = "px",res = 300)
ggplot(OE_heatmap_input,aes(x = condition, y = gene_symbols, fill = Z_score))+geom_tile()+
  scale_fill_distiller(palette = 'Spectral',limits = c(-1.5,1.5), oob = squish)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.text.y=ggtext::element_markdown(color=y_tick_colors))+
  xlab('')+ylab('')
dev.off()