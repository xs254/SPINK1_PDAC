library(ggplot2)
library(eulerr)
library(Gviz)
# Define helper function
get_inconclusive_genes <- function(annotated_resgions)
{
  inconclusive_genes <- c()
  for(gene in unique(annotated_resgions$geneId))
  {
    if(sum(gene==annotated_resgions$geneId)>1)
    {
      gene_directions <- annotated_resgions$direction[gene==annotated_resgions$geneId]
      if(length(unique(gene_directions))>1)
      {
        inconclusive_genes <- c(inconclusive_genes,gene)
      }
    }
  }
  return(inconclusive_genes)
}
setwd('C:/Users/xshan/OneDrive/Desktop/R code/CutNTag_SPINK1_KO')
# Import and filter H3K4Me3 and H3K27Me3 data
H3K4Me3 <- read.csv('results/H3K4Me3_csaw_edgeR_annotated_results.csv',header = TRUE)
H3K4Me3_filtered <- H3K4Me3[H3K4Me3$FDR<=0.05 & H3K4Me3$direction!='mixed',]
H3K27Me3 <- read.csv('results/H3K27Me3_csaw_edgeR_annotated_results.csv',header = TRUE)
H3K27Me3_filtered <- H3K27Me3[H3K27Me3$FDR<=0.05 & H3K27Me3$direction!='mixed',]
# For H3K4Me3, focus on regions close to TSS
H3K4Me3_filtered_TSS <- H3K4Me3_filtered[abs(H3K4Me3_filtered$distanceToTSS)<=1000,]
# Remove inconclusive genes in H3K4Me3
H3K4Me3_inconclusive <- get_inconclusive_genes(H3K4Me3_filtered_TSS)
H3K4Me3_filtered_TSS_conclusive <- H3K4Me3_filtered_TSS[!(H3K4Me3_filtered_TSS$geneId %in% H3K4Me3_inconclusive),]
# Remove inconclusive genes in H3K27Me3
H3K27Me3_inconclusive <- get_inconclusive_genes(H3K27Me3_filtered)
H3K27Me3_filtered_conclusive <- H3K27Me3_filtered[!(H3K27Me3_filtered$geneId %in% H3K27Me3_inconclusive),]
# Get up/downregulated genes
logFC_threshold <- 1
H3K4Me3_up_genes <- H3K4Me3_filtered_TSS_conclusive$SYMBOL[H3K4Me3_filtered_TSS_conclusive$rep.logFC>=logFC_threshold]
H3K4Me3_up_genes <- na.omit(unique(H3K4Me3_up_genes))
H3K4Me3_down_genes <- H3K4Me3_filtered_TSS_conclusive$SYMBOL[H3K4Me3_filtered_TSS_conclusive$rep.logFC<=-logFC_threshold]
H3K4Me3_down_genes <- na.omit(unique(H3K4Me3_down_genes))

H3K27Me3_up_genes <- H3K27Me3_filtered_conclusive$SYMBOL[H3K27Me3_filtered_conclusive$rep.logFC>=logFC_threshold]
H3K27Me3_up_genes <- na.omit(unique(H3K27Me3_up_genes))
H3K27Me3_down_genes <- H3K27Me3_filtered_conclusive$SYMBOL[H3K27Me3_filtered_conclusive$rep.logFC<=-logFC_threshold]
H3K27Me3_down_genes <- na.omit(unique(H3K27Me3_down_genes))
# Check overlap between up and down genes
for_venn_diagram <- list(H3K4Me3_up_genes,H3K4Me3_down_genes,H3K27Me3_up_genes,H3K27Me3_down_genes)
names(for_venn_diagram) <- c('H3K4Me3 up','H3K4Me3 down','H3K27Me3 up','H3K27Me3 down')
for_venn_diagram <- euler(for_venn_diagram)
png('results/Overview_of_H3K4Me3_H3K27Me3_CutNTag.png',width = 1200,height = 1500,units = 'px',res = 300)
plot(for_venn_diagram,edges = TRUE,quantities = TRUE,
     fills = c('salmon','gold2','palegreen3','dodgerblue2'))
dev.off()
# Check overlap with SPINK1 promotion gene identified from RNA and ATACseq
SPINK1_promo_genes <- readLines('C:/Users/xshan/OneDrive/Desktop/R code/DEseq_SPINK1_KO_v2/results/Epigenetic_promotion_genes_by_SPINK1.txt')
SPINK1_suppression_genes <- readLines('C:/Users/xshan/OneDrive/Desktop/R code/DEseq_SPINK1_KO_v2/results/Epigenetic_supression_genes_by_SPINK1.txt')

for_venn_diagram <- list(SPINK1_promo_genes,
                         H3K4Me3_up_genes,H3K4Me3_down_genes,H3K27Me3_up_genes,H3K27Me3_down_genes)
names(for_venn_diagram) <- c('Promotion',
                             'H3K4Me3 up','H3K4Me3 down','H3K27Me3 up','H3K27Me3 down')
for_venn_diagram <- euler(for_venn_diagram)
png('results/SPINK1_promotion_genes_overlap.png',width = 1200,height = 1500,units = 'px',res = 300)
plot(for_venn_diagram,edges = TRUE,quantities = FALSE,labels = FALSE,
     fills = c('salmon',rep('grey',4)))
dev.off()
png('results/SPINK1_promotion_genes_overlap_labelled.png',width = 3000,height = 2000,units = 'px',res = 300)
plot(for_venn_diagram,edges = TRUE,quantities = TRUE,labels = TRUE,
     fills = c('salmon',rep('grey',4)))
dev.off()

for_venn_diagram <- list(SPINK1_suppression_genes,
                         H3K4Me3_up_genes,H3K4Me3_down_genes,H3K27Me3_up_genes,H3K27Me3_down_genes)
names(for_venn_diagram) <- c('Suppression',
                             'H3K4Me3 up','H3K4Me3 down','H3K27Me3 up','H3K27Me3 down')
for_venn_diagram <- euler(for_venn_diagram)
png('results/SPINK1_suppression_genes_overlap.png',width = 1200,height = 1500,units = 'px',res = 300)
plot(for_venn_diagram,edges = TRUE,quantities = FALSE,labels = FALSE,
     fills = c('dodgerblue2',rep('grey',4)))
dev.off()
png('results/SPINK1_suppression_genes_overlap_labelled.png',width = 3000,height = 2000,units = 'px',res = 300)
plot(for_venn_diagram,edges = TRUE,quantities = TRUE,labels = TRUE,
     fills = c('dodgerblue2',rep('grey',4)))
dev.off()
# Check overlap with stemness related genes
all_genes <- read.delim('C:/Users/xshan/OneDrive/Desktop/R code/DEseq_SPINK1_KO_v2/results/Epigenetic_promotion_CellMarker_Augmented_2021_table.txt',sep = '\t')
CellMarker_genes <- c(strsplit(all_genes$Genes[1],';')[[1]],
                      strsplit(all_genes$Genes[4],';')[[1]],
                      strsplit(all_genes$Genes[5],';')[[1]],
                      strsplit(all_genes$Genes[6],';')[[1]])
CellMarker_genes <- unique(CellMarker_genes)
all_genes <- read.delim('C:/Users/xshan/OneDrive/Desktop/R code/DEseq_SPINK1_KO_v2/results/Epigenetic_promotion_PanglaoDB_Augmented_2021_table.txt',sep = '\t')
PanglaoDB_genes <- c(strsplit(all_genes$Genes[1],';')[[1]],
                     strsplit(all_genes$Genes[11],';')[[1]],
                     strsplit(all_genes$Genes[12],';')[[1]])
PanglaoDB_genes <- unique(PanglaoDB_genes)
stem_related <- union(CellMarker_genes,PanglaoDB_genes)
writeLines(stem_related,'results/stem_related_genes_from_enrichr.txt')
for_venn_diagram <- list(stem_related,
                         H3K4Me3_up_genes,H3K4Me3_down_genes,H3K27Me3_up_genes,H3K27Me3_down_genes)
names(for_venn_diagram) <- c('Stem or progenitor',
                             'H3K4Me3 up','H3K4Me3 down','H3K27Me3 up','H3K27Me3 down')
for_venn_diagram <- euler(for_venn_diagram)
png('results/Stemness_genes_overlap.png',width = 1200,height = 1500,units = 'px',res = 300)
plot(for_venn_diagram,edges = TRUE,quantities = FALSE,labels = FALSE,
     fills = c('skyblue2',rep('grey',4)))
dev.off()
png('results/Stemness_genes_overlap_labelled.png',width = 3000,height = 2000,units = 'px',res = 300)
plot(for_venn_diagram,edges = TRUE,quantities = TRUE,labels = TRUE,
     fills = c('skyblue2',rep('grey',4)))
dev.off()
# Get stem or progenitor genes down in H3K4Me3 after KO
H3K4Me3_down_stem <- stem_related[stem_related %in% H3K4Me3_down_genes]
H3K27Me3_up_stem <- stem_related[stem_related %in% H3K27Me3_up_genes]

# Integrated visualization of RNA seq, ATACseq and CutNTag seq for examplery genes
# 
# H3K4Me3_win_counts <- readRDS('results/H3K4Me3_csaw_normalized_counts.RDS')
# H3K4Me3_cpm <- csaw::asDGEList(H3K4Me3_win_counts)
# H3K4Me3_cpm <- edgeR::cpm(H3K4Me3_win_counts,log = FALSE)
# H3K4Me3_cpm_aggregated <- data.frame(KO = rowMeans(H3K4Me3_cpm[,1:4]),
#                                      SCR = rowMeans(H3K4Me3_cpm[,5:8]))
# H3K4Me3_gr <- rowRanges(H3K4Me3_win_counts)
# mcols(H3K4Me3_gr)$KO_mean_CPM <- H3K4Me3_cpm_aggregated$KO
# mcols(H3K4Me3_gr)$SCR_mean_CPM <- H3K4Me3_cpm_aggregated$SCR
# # For H3K27Me3
# H3K27Me3_win_counts <- readRDS('results/H3K27Me3_csaw_normalized_counts.RDS')
# H3K27Me3_cpm <- csaw::asDGEList(H3K27Me3_win_counts)
# H3K27Me3_cpm <- edgeR::cpm(H3K27Me3_win_counts,log = FALSE)
# H3K27Me3_cpm_aggregated <- data.frame(KO = rowMeans(H3K27Me3_cpm[,1:4]),
#                                      SCR = rowMeans(H3K27Me3_cpm[,5:7]))
# H3K27Me3_gr <- rowRanges(H3K27Me3_win_counts)
# mcols(H3K27Me3_gr)$KO_mean_CPM <- H3K27Me3_cpm_aggregated$KO
# mcols(H3K27Me3_gr)$SCR_mean_CPM <- H3K27Me3_cpm_aggregated$SCR
# # Plot tracks
# gtrack <- GenomeAxisTrack()
# # For NES H3K4Me3
# GOI_chr <- na.omit(unique(H3K4Me3$geneChr[H3K4Me3$SYMBOL=='NES']))
# GOI_start <- na.omit(unique(H3K4Me3$geneStart[H3K4Me3$SYMBOL=='NES']))
# GOI_end <- na.omit(unique(H3K4Me3$geneEnd[H3K4Me3$SYMBOL=='NES']))
# 
# KO_track <- DataTrack(range = H3K4Me3_gr,data = 'KO_mean_CPM',
#                       type = 'histogram',name = 'KO',
#                       cex.title = 2,cex.axis = 1,ylim = c(0,1))
# SCR_track <- DataTrack(range = H3K4Me3_gr,data = 'SCR_mean_CPM',
#                       type = 'histogram',name = 'SCR',
#                       cex.title = 2,cex.axis = 1,ylim = c(0,1))
# png('results/NES_H3K4Me3_track.png',width = 3000,height = 1500,units = 'px',res = 300)
# plotTracks(list(gtrack,SCR_track, KO_track),
#   chromosome = GOI_chr,from = GOI_start,to = GOI_end)
# dev.off()
# # For NES H3K27Me3
# GOI_chr <- na.omit(unique(H3K27Me3$geneChr[H3K27Me3$SYMBOL=='NES']))
# GOI_start <- na.omit(unique(H3K27Me3$geneStart[H3K27Me3$SYMBOL=='NES']))
# GOI_end <- na.omit(unique(H3K27Me3$geneEnd[H3K27Me3$SYMBOL=='NES']))
# 
# KO_track <- DataTrack(range = H3K27Me3_gr,data = 'KO_mean_CPM',
#                       type = 'histogram',name = 'KO',
#                       cex.title = 2,cex.axis = 1,ylim = c(0,1))
# SCR_track <- DataTrack(range = H3K27Me3_gr,data = 'SCR_mean_CPM',
#                        type = 'histogram',name = 'SCR',
#                        cex.title = 2,cex.axis = 1,ylim = c(0,1))
# png('results/NES_H3K27Me3_track.png',width = 3000,height = 1500,units = 'px',res = 300)
# plotTracks(list(gtrack,SCR_track, KO_track),
#            chromosome = GOI_chr,from = GOI_start,to = GOI_end)
# dev.off()
# # For YES1
# GOI_chr <- na.omit(unique(H3K4Me3$geneChr[H3K4Me3$SYMBOL=='YES1']))
# GOI_start <- na.omit(unique(H3K4Me3$geneStart[H3K4Me3$SYMBOL=='YES1']))
# GOI_end <- na.omit(unique(H3K4Me3$geneEnd[H3K4Me3$SYMBOL=='YES1']))
# 
# KO_track <- DataTrack(range = H3K4Me3_gr,data = 'KO_mean_CPM',
#                       type = 'histogram',name = 'KO',
#                       cex.title = 2,cex.axis = 1,ylim = c(0,5))
# SCR_track <- DataTrack(range = H3K4Me3_gr,data = 'SCR_mean_CPM',
#                        type = 'histogram',name = 'SCR',
#                        cex.title = 2,cex.axis = 1,ylim = c(0,5))
# png('results/YES1_H3K4Me3_track.png',width = 3000,height = 1500,units = 'px',res = 300)
# plotTracks(list(gtrack,SCR_track, KO_track),
#            chromosome = GOI_chr,from = GOI_start,to = GOI_end[2])
# dev.off()