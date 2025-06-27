library(DiffBind)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(tidyverse)
library(clusterProfiler)
library(patchwork)
# Load the peak files from MACS2 to DiffBind
setwd('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\ATACseq_SPINK1_KO')
ATAC_peakset<-dba.peakset(NULL,peaks = '01_0FJU_025ZYale_SCR-P-6_ATAC_hs_i20-48_MACS2_peaks.xls',
                          peak.caller = 'macs',sampID = "scramble_1",treatment = 'scramble',replicate = 1,
                          bamReads = '01_0FJU_025ZYale_SCR-P-6_ATAC_hs_i20-48_peak_calling_ready.bam')

ATAC_peakset<-dba.peakset(ATAC_peakset,peaks = '05_0FJY_025ZYale_SCR_ATAC_hs_i24-4_MACS2_peaks.xls',
                          peak.caller = 'macs',sampID = "scramble_2",treatment = 'scramble',replicate = 2,
                          bamReads = '05_0FJY_025ZYale_SCR_ATAC_hs_i24-4_peak_calling_ready.bam')

ATAC_peakset<-dba.peakset(ATAC_peakset,peaks = '02_0FJV_025ZYale_2g-6-SP-P-4_ATAC_hs_i21-1_MACS2_peaks.xls',
                          peak.caller = 'macs',sampID = "knockout_1",treatment = 'knockout',replicate = 1,
                          bamReads = '02_0FJV_025ZYale_2g-6-SP-P-4_ATAC_hs_i21-1_peak_calling_ready.bam')

ATAC_peakset<-dba.peakset(ATAC_peakset,peaks = '06_0FJZ_025ZYale_2G-6-SP_ATAC_hs_i25-521_MACS2_peaks.xls',
                          peak.caller = 'macs',sampID = "knockout_2",treatment = 'knockout',replicate = 2,
                          bamReads = '06_0FJZ_025ZYale_2G-6-SP_ATAC_hs_i25-521_peak_calling_ready.bam')
# Plot Venn diagrams of replicates
setwd('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\ATACseq_SPINK1_KO\\Results')
png(file = 'Venn_diagram_for_replicates_of_scrmble_control.png',width = 1200,height = 800,units = "px",res = 150)
dba.plotVenn(ATAC_peakset,ATAC_peakset$masks$scramble, main = "Open chromatic region overlaps in scramble control replicates")
dev.off()
png(file = 'Venn_diagram_for_replicates_of_knockout.png',width = 1200,height = 800,units = "px",res = 150)
dba.plotVenn(ATAC_peakset,ATAC_peakset$masks$knockout, main = "Open chromatic region overlaps in knockout replicates")
dev.off()
# Calculate peak counts, this step can take tens of minutes to finish
setwd('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\ATACseq_SPINK1_KO')
ATAC_peakset<-dba.count(ATAC_peakset,bParallel = TRUE)
# Plot heatmap for correlation
setwd('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\ATACseq_SPINK1_KO\\Results')
png(file = 'Heatmap_of_sample_correlation.png',width = 1500,height = 1500,units = "px",res = 150)
plot(ATAC_peakset, main="Correlation plot of studied samples")
dev.off()
# Normalize 
ATAC_peakset <- dba.normalize(ATAC_peakset)
# Set up contrast for diffbind analysis
ATAC_peakset <- dba.contrast(ATAC_peakset, minMembers = 2, reorderMeta = list(Treatment = 'scramble'))
# Run diffbind analysis
# No need to apply grey/black list as the bam has already been filtered
ATAC_peakset <- dba.analyze(ATAC_peakset,bGreylist = FALSE,bBlacklist = FALSE)
# Retrieve differentially up/down sites and export it
diffbind_sites <- dba.report(ATAC_peakset)
write.csv(as.data.frame(diffbind_sites),'diffbind_results.csv',quote = FALSE,row.names = TRUE)
# Prepare bed file for HOMER motif analysis. HOMER and MEME needs slightly different bed files
# MEME only requires chrom, chromStart and chromEnd.
# HOMER requires chrom, chromStart, chromEnd, ID,Fold,Strand
# HOMER does not actually use Fold, but there just needs to be a column there to ensure Strand is the 6th column
# For ATACseq, strand does not matter, just put + to all to let HOMER work.
diffbind_sites_bed <- as.data.frame(diffbind_sites)
diffbind_sites_bed <- diffbind_sites_bed[,c(1:3)]
colnames(diffbind_sites_bed)[1:3] <-  c('chrom','chromStart','chromEnd')
diffbind_sites_bed['ID'] <- row.names(diffbind_sites_bed)
diffbind_sites_bed['Fold'] <- diffbind_sites$Fold
diffbind_sites_bed['Strand'] <- rep('+',dim(diffbind_sites_bed)[1])
diffbind_sites_bed <- diffbind_sites_bed[diffbind_sites_bed$chrom %in% c(paste('chr',1:22,sep = ''),'chrX'),]
diffbind_sites_bed_up_homer <- diffbind_sites_bed[diffbind_sites_bed$Fold>0,]
diffbind_sites_bed_up_meme <- diffbind_sites_bed_up_homer[,1:3]
diffbind_sites_bed_down_homer <- diffbind_sites_bed[diffbind_sites_bed$Fold<0,]
diffbind_sites_bed_down_meme <- diffbind_sites_bed_down_homer[,1:3]
dir.create('bed files for HOMER motif analysis')
write.table(diffbind_sites_bed_up_homer, file = 'bed files for HOMER motif analysis\\SPINK1_KO_diffbind_sites_bed_up_HOMER.bed',sep = '\t',row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(diffbind_sites_bed_down_homer, file = 'bed files for HOMER motif analysis\\SPINK1_KO_diffbind_sites_bed_down_HOMER.bed',sep = '\t',row.names = FALSE,col.names = FALSE,quote = FALSE)
dir.create('bed files for MEME motif analysis')
write.table(diffbind_sites_bed_up_meme, file = 'bed files for MEME motif analysis\\SPINK1_KO_diffbind_sites_bed_up_MEME.bed',sep = '\t',row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(diffbind_sites_bed_down_meme, file = 'bed files for MEME motif analysis\\SPINK1_KO_diffbind_sites_bed_down_MEME.bed',sep = '\t',row.names = FALSE,col.names = FALSE,quote = FALSE)
# Visualize diffbind results
setwd('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\ATACseq_SPINK1_KO')
profiles <- dba.plotProfile(ATAC_peakset)
setwd('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\ATACseq_SPINK1_KO\\Results')
png(file = 'Gain_and_loss_profiles.png',width = 1000,height = 1500,units = "px",res = 150)
dba.plotProfile(profiles)
dev.off()
# Visualize coverage of differentially up/down sites
chrosomes_to_plot <- c()
for(n in 1:22)
{
  chrosomes_to_plot <- c(chrosomes_to_plot,paste0('chr',n))
}
chrosomes_to_plot <- c(chrosomes_to_plot,'chrX')
png(file = 'ATAC_seq_coverage_over_chrosomes.png',width = 750,height = 2000,units = "px",res = 150)
covplot(diffbind_sites,weightCol = 'Conc',title = 'ATACseq coverage over chromosomes',chrs = chrosomes_to_plot)
dev.off()
# Annotate the differentially up/down sites to the closest genes
# If a peak lies within tssRegion, it will be annotated as promoter
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
diffbind_sites_annotated <- annotatePeak(diffbind_sites,tssRegion = c(-1000, 1000),TxDb = txdb, annoDb = "org.Hs.eg.db")
diffbind_sites_annotated_df <- as.data.frame(diffbind_sites_annotated)
write.csv(diffbind_sites_annotated_df,'diffbind_annotated_results.csv',
          quote = TRUE,row.names = FALSE)
# Visualize distribution of annotation
png(file = 'ATAC_seq_diffbind_annotated_peaks_location_pie_chart.png',width = 1200,height = 600,units = "px",res = 150)
plotAnnoPie(diffbind_sites_annotated)
dev.off()
# One gene can have multiple diffbind sites, integrate them into an overall up/down/contradictory status
setwd('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\ATACseq_SPINK1_KO\\Results')
diffbind_sites_annotated_df <- read.csv('diffbind_annotated_results.csv')
diffbind_sites_annotated_df <- diffbind_sites_annotated_df[!is.na(diffbind_sites_annotated_df$SYMBOL),]
diffbind_genes <- unique(diffbind_sites_annotated_df$SYMBOL)
diffbind_integrated <- data.frame(gene_symbols = diffbind_genes,
                                  integrated_accessibility = rep(NA,length(diffbind_genes)),
                                  n_diffbind_sites = rep(NA,length(diffbind_genes)),
                                  shortest_TSS_distance = rep(NA,length(diffbind_genes)))
for(n in 1:length(diffbind_genes))
{
  gene_profile <- diffbind_sites_annotated_df[diffbind_sites_annotated_df$SYMBOL==diffbind_integrated$gene_symbols[n],]
  diffbind_integrated$n_diffbind_sites[n] <- dim(gene_profile)[1]
  if(dim(gene_profile)[1]==1)
  {
    diffbind_integrated$shortest_TSS_distance[n] <- gene_profile$distanceToTSS
    if(gene_profile$Fold>0){diffbind_integrated$integrated_accessibility[n] = 1}
    else if(gene_profile$Fold<0){diffbind_integrated$integrated_accessibility[n] = -1}
  }
  else
  {
    diffbind_integrated$shortest_TSS_distance[n] <- min(gene_profile$distanceToTSS)
    all_up <- unique(gene_profile$Fold > 0)
    if(length(all_up)==1)
    {
      if(all_up == TRUE){diffbind_integrated$integrated_accessibility[n] = 1}
      else if(all_up == FALSE){diffbind_integrated$integrated_accessibility[n] = -1}
    }
    else if(length(all_up)!=1){diffbind_integrated$integrated_accessibility[n] = 0}
  }
}
sum(diffbind_integrated$integrated_accessibility==1)
# Integrate with RNAseq data
KO_deseq_results <- read.csv('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\DEseq_SPINK1_KO_v2\\Results\\DEseq_results.csv')
KO_deseq_DEGs <- KO_deseq_results[KO_deseq_results$padj<=0.05,]
KO_deseq_DEGs_diffbind <- KO_deseq_DEGs[KO_deseq_DEGs$gene_symbols %in% diffbind_integrated$gene_symbols,]
diffbind_integrated_DEGs <- diffbind_integrated[diffbind_integrated$gene_symbols %in% KO_deseq_DEGs$gene_symbols,]
ATAC_RNA_DEGs <- merge(diffbind_integrated_DEGs,KO_deseq_DEGs_diffbind,by = 'gene_symbols')
ATAC_RNA_DEGs[,'ATAC_RNA_seq_consensus'] <- ATAC_RNA_DEGs$integrated_accessibility*sign(ATAC_RNA_DEGs$log2FoldChange)
write.csv(ATAC_RNA_DEGs,'ATAC_RNA_DEGs.csv',row.names = FALSE)
# Visualize the consensus between RNAseq and ATACseq
ATAC_RNA_DEGs <- read.csv('ATAC_RNA_DEGs.csv')
ATAC_RNA_seq_consensus_stat <- data.frame(results = factor(c('Consensually upregulated','Consensually downregulated','Disagreed','ATACseq inconclusive'),
                                                           levels = c('Consensually upregulated','Consensually downregulated','Disagreed','ATACseq inconclusive')),
                                          frequency = c(sum(ATAC_RNA_DEGs$ATAC_RNA_seq_consensus==1 & ATAC_RNA_DEGs$integrated_accessibility==1),
                                                        sum(ATAC_RNA_DEGs$ATAC_RNA_seq_consensus==1 & ATAC_RNA_DEGs$integrated_accessibility==-1),
                                                        sum(ATAC_RNA_DEGs$ATAC_RNA_seq_consensus==-1),
                                                        sum(ATAC_RNA_DEGs$ATAC_RNA_seq_consensus==0)))
ATAC_RNA_seq_consensus_stat['probability'] = ATAC_RNA_seq_consensus_stat$frequency/sum(ATAC_RNA_seq_consensus_stat$frequency)
png(file = 'Piechart for RNAseq ATACseq integration consensus statistics.png',width = 1500,height = 1000,units = 'px',res = 300)
ggplot(ATAC_RNA_seq_consensus_stat,aes(x = '',y = probability,fill = results))+
  geom_bar(stat='identity', width=1, color='black',linewidth=0.1)+xlab('')+ylab('')+
  coord_polar("y", start=0)+theme_minimal()+ scale_fill_manual(name = '',values = c('Consensually upregulated'='coral2','Consensually downregulated'='dodgerblue2','Disagreed'='palegreen','ATACseq inconclusive'='gray'))+
  theme(axis.ticks = element_blank(),axis.text.y = element_blank(),axis.text.x = element_blank(),panel.grid = element_blank()) 
dev.off()
# Visualize binding site distribution of RNAseq ATACseq agreed genes
bmp(file = 'Histogram for binding site distribution of RNAseq ATACseq agreed genes.bmp', width = 1500,height = 1000,units = 'px',res = 300)
ggplot(ATAC_RNA_DEGs[ATAC_RNA_DEGs$ATAC_RNA_seq_consensus==1,],aes(x = shortest_TSS_distance))+
  geom_histogram()+xlim(c(-5000,5000))+xlab('shortest distance to TSS')+theme_minimal()
dev.off()
# enrichR analysis of RNAseq ATACseq agreed genes
genes_to_enrich_up <- ATAC_RNA_DEGs$gene_symbols[ATAC_RNA_DEGs$ATAC_RNA_seq_consensus==1 & ATAC_RNA_DEGs$log2FoldChange > 0]
writeLines(genes_to_enrich_up, 'ATAC_RNA_both_up_gene_lists_threshold_0.txt') # Export to txt for enrichR
genes_to_enrich_down <- ATAC_RNA_DEGs$gene_symbols[ATAC_RNA_DEGs$ATAC_RNA_seq_consensus==1 & ATAC_RNA_DEGs$log2FoldChange < 0]
writeLines(genes_to_enrich_down, 'ATAC_RNA_both_down_gene_lists_threshold_0.txt') # Export to txt for enrichR
# After enrichR analysis, H3K27Ac and H3K4Me3 and embryonic stem cell related gene set were found to be enriched
# Visualize the degree of overlapping of these three gene sets
setwd('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\ATACseq_SPINK1_KO\\Results')
epiRoadMap_enrichR <- read.table('ATAC_RNA_down_Epigenomics_Roadmap_HM_ChIP-seq_table.txt',sep = '\t', header = TRUE)
PDB_enrichR <- read.table('ATAC_RNA_down_ProteomicsDB_2020_table.txt', sep = '\t', header = TRUE)
TabulaMuris_enrichR <- read.table('ATAC_RNA_down_Tabula_Muris_table.txt', sep = '\t', header = TRUE)
CellMarker2021_enrichR <- read.table('ATAC_RNA_down_CellMarker_Augmented_2021_table.txt', sep = '\t', header = TRUE)

H3K4Me3_genes <- strsplit(epiRoadMap_enrichR[3,'Genes'], split = ';')
H3K4Me3_genes <- H3K4Me3_genes[[1]]
H3K27Ac_genes <- strsplit(epiRoadMap_enrichR[1,'Genes'], split = ';')
H3K27Ac_genes <- H3K27Ac_genes[[1]]
PDB_stem_genes <- vector(mode = 'list',length = 3)
for(n in 1:3)
{
  temp <- strsplit(PDB_enrichR[n,'Genes'], split = ';')
  PDB_stem_genes[[n]] <- temp[[1]]
}
PDB_stem_genes <- Reduce(union,PDB_stem_genes)
TabulaMuris_stem_genes <- strsplit(TabulaMuris_enrichR[1,'Genes'], split = ';')
TabulaMuris_stem_genes <- TabulaMuris_stem_genes[[1]]
CellMarker2021_stem_genes <- strsplit(CellMarker2021_enrichR[1,'Genes'], split = ';')
CellMarker2021_stem_genes <- CellMarker2021_stem_genes[[1]]
all_stem_genes <- Reduce(union,list(PDB_stem_genes,TabulaMuris_stem_genes,CellMarker2021_stem_genes))
c('ACACA','BCL7C','Yes1','ASPM','HMMR','NES') %in% union(H3K4Me3_genes,H3K27Ac_genes)
VennDiagram::venn.diagram(x = list(H3K4Me3_genes, H3K27Ac_genes, all_stem_genes),
  category.names = c('H3K4Me3', 'H3K27Ac', 'Stem cell related genes'), filename = 'enrichR_venn_diagramm.png',output=TRUE,
  height = 2500, width = 2500, resolution = 300,lwd = 0, cat.cex = 2, cex = 2,
  fill = c('springgreen','lightskyblue','coral'))
stem_genes_epi_regulated <- intersect(stem_genes, union(H3K27Ac_genes, H3K4Me3_genes))
print(c('SRSF8','EIF4EBP2','SH3BP4','PTK2','DNAJB6','DNAJC10') %in% stem_genes_epi_regulated)
# Integrated RNA-seq ATAC-seq EMT panel analysis
EMT_markers <- read.csv('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\DEseq_SPINK1_KO_v2\\Total EMT gene list.csv',header = FALSE)
EMT_markers <- EMT_markers[1:50,c(-1,-8)]
colnames(EMT_markers) <- c('Gene_Symbols','Gene_Length','Type','Chrosome','Full_Name','Description')
EMT_markers_RNAseq <- KO_deseq_results[KO_deseq_results$gene_symbols %in% EMT_markers$Gene_Symbols,]
EMT_markers_RNAseq_valid <- na.omit(EMT_markers_RNAseq)
EMT_markers_consensus_DEGs <- ATAC_RNA_DEGs[ATAC_RNA_DEGs$gene_symbols %in% EMT_markers$Gene_Symbols &
                                              ATAC_RNA_DEGs$ATAC_RNA_seq_consensus == 1,]
nEMTs_invalid <- sum(is.na(EMT_markers_RNAseq))
nEMTs_notSig <- sum(EMT_markers_RNAseq_valid$padj>0.05)
nEMTs_upSigConsensus <- sum(EMT_markers_RNAseq_valid$padj<=0.05 & EMT_markers_RNAseq_valid$log2FoldChange>0 & EMT_markers_RNAseq_valid$gene_symbols %in% EMT_markers_consensus_DEGs$gene_symbols)
nEMTs_upSigRest <- sum(EMT_markers_RNAseq_valid$padj<=0.05 & EMT_markers_RNAseq_valid$log2FoldChange>0 & !(EMT_markers_RNAseq_valid$gene_symbols %in% EMT_markers_consensus_DEGs$gene_symbols))
nEMTs_downSigConsensus <- sum(EMT_markers_RNAseq_valid$padj<=0.05 & EMT_markers_RNAseq_valid$log2FoldChange<0 & EMT_markers_RNAseq_valid$gene_symbols %in% EMT_markers_consensus_DEGs$gene_symbols)
nEMTs_downSigRest <- sum(EMT_markers_RNAseq_valid$padj<=0.05 & EMT_markers_RNAseq_valid$log2FoldChange<0 & !(EMT_markers_RNAseq_valid$gene_symbols %in% EMT_markers_consensus_DEGs$gene_symbols))
EMT_markers_statistics <- data.frame(results = factor(c('Invalid','Not significant','Consensually upregulated','RNA-seq upregulated','Consensually downregulated','RNA-seq downregulated'),
                                                      levels = c('Consensually upregulated','RNA-seq upregulated','Consensually downregulated','RNA-seq downregulated','Not significant','Invalid')),
                                     frequency = c(nEMTs_invalid,nEMTs_notSig,nEMTs_upSigConsensus,nEMTs_upSigRest,nEMTs_downSigConsensus,nEMTs_downSigRest))
EMT_markers_statistics['probability'] <- EMT_markers_statistics['frequency']/sum(EMT_markers_statistics['frequency'])
# Visualize RNAseq ATACseq consensus in EMT panel
bmp(file = 'Piechart for EMT panel RNAseq ATACseq integration consensus statistics.bmp',width = 750,height = 500,units = 'px',res = 150)
ggplot(EMT_markers_statistics,aes(x = '',y = probability,fill = results))+
  geom_bar(stat='identity', width=1, color='black',linewidth=0.1)+xlab('')+ylab('')+
  coord_polar("y", start=0)+theme_minimal()+ 
  scale_fill_manual(name = '',values = c('Invalid' = 'gray','Not significant' = 'palegreen','Consensually upregulated'='red3','RNA-seq upregulated'='lightcoral','Consensually downregulated'='blue3','RNA-seq downregulated' = 'lightblue'))+
  theme(axis.ticks = element_blank(),axis.text.y = element_blank(),axis.text.x = element_blank(),panel.grid = element_blank()) 
dev.off()
# Heatmap of EMT panel in RNAseq with ATACseq results supplemented
RNAseq_cts_matrix_normalized <- read.csv('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\DEseq_SPINK1_KO_v2\\Results\\normalized_counts_matrix.csv')
EMT_markers_RNAseq_valid <- EMT_markers_RNAseq_valid[order(EMT_markers_RNAseq_valid$log2FoldChange,decreasing = TRUE),]
EMT_DEGs <- EMT_markers_RNAseq_valid$gene_symbols[EMT_markers_RNAseq_valid$padj<=0.05]
RNAseq_cts_matrix_EMT_markers_valid <- RNAseq_cts_matrix_normalized[RNAseq_cts_matrix_normalized$gene_symbols %in% EMT_markers_RNAseq_valid$gene_symbols,]
rownames(RNAseq_cts_matrix_EMT_markers_valid) <- RNAseq_cts_matrix_EMT_markers_valid$gene_symbols
RNAseq_cts_matrix_EMT_markers_valid <- RNAseq_cts_matrix_EMT_markers_valid[,c(-1,-2)]
RNAseq_cts_matrix_EMT_markers_valid_scaled <- as.data.frame(t(scale(t(RNAseq_cts_matrix_EMT_markers_valid))))
RNAseq_cts_matrix_EMT_markers_valid_scaled <- RNAseq_cts_matrix_EMT_markers_valid_scaled[EMT_markers_RNAseq_valid$gene_symbols,]
heatmap_input <- cbind(rownames(RNAseq_cts_matrix_EMT_markers_valid_scaled),RNAseq_cts_matrix_EMT_markers_valid_scaled)
colnames(heatmap_input) <- c('gene_symbols','Scramble Control 1','Scramble Control 2','SPINK1 KO 1','SPINK1 KO 2')
y_tick_colors <- heatmap_input$gene_symbols %in% EMT_markers_consensus_DEGs$gene_symbols
y_tick_colors[y_tick_colors==TRUE] = 'darkgoldenrod3' #Highlight consensually up/down genes in RNA/ATAC-seq
y_tick_colors[y_tick_colors==FALSE] = 'black'
heatmap_input[EMT_DEGs,'gene_symbols'] <- paste0('*',heatmap_input[EMT_DEGs,'gene_symbols']) #Highlight degs
heatmap_input <- tidyr::pivot_longer(heatmap_input,cols = !gene_symbols,names_to = 'condition',values_to = 'Z_score')
heatmap_input$gene_symbols <- factor(heatmap_input$gene_symbols,levels = unique(heatmap_input$gene_symbols))
heatmap_input$condition <- factor(heatmap_input$condition, levels = unique(heatmap_input$condition))
bmp(file = 'Heatmap for EMT markers with RNAseq and ATACseq integrated.bmp',width = 800,height = 1000,units = "px",res = 150)
ggplot(heatmap_input,aes(x = condition, y = gene_symbols, fill = Z_score))+geom_tile()+
  scale_fill_distiller(palette = 'Spectral',limits = c(-1.5,1.5), oob = scales::squish)+
  theme(axis.text.y=ggtext::element_markdown(color=y_tick_colors),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  xlab('')+ylab('')
dev.off()
# Visualize the distribution of peaks of selected genes/regions by Lollipop plot
genes_to_plot <- c('CDH1','CDH2','VIM','ZEB1')
sample_names <- c('Scramble control 1','Scramble control 2','Knockout 1','Knockout 2')
genes_to_plot_ID <- mapIds(org.Hs.eg.db,keys = genes_to_plot,column = 'ENTREZID',
                           keytype = 'SYMBOL',multivals = 'first')
genes_info <- as.data.frame(genes(txdb,filter = list(gene_id = genes_to_plot_ID))) #get the genomic location from txdb
genes_info['Symbol'] <- mapIds(org.Hs.eg.db,keys = genes_info$gene_id,column = 'SYMBOL',
                               keytype = 'ENTREZID',multivals = 'first')
for(i_gene in 1:dim(genes_info)[1])
{
  max_peak_height <- 0 # Store the largest peak height (RPKM) to set a global ylim for all samples
  for(i_sample in 1:length(sample_names))
  {
    peaks <- ATAC_peakset$peaks[[i_sample]]
    peaks <- peaks[peaks$Chr == genes_info$seqnames[i_gene],]
    # If the center of a peak is within 3000 bp upstream to the end of the gene, then it is regarded
    # as affecting that gene, and will be ploted out. 
    peaks['peaks_center'] <- (peaks$Start+peaks$End)/2-genes_info$start[i_gene]
    peaks <- peaks[peaks$peaks_center>=-3000 & 
                     peaks$peaks_center<=genes_info$end[i_gene]-genes_info$start[i_gene],]
    max_peak_height <- max(max(peaks$RPKM),max_peak_height)
    if(i_sample==1)
    {
      gene_plot <- ggplot(peaks, aes(x = peaks_center, y = RPKM))+ 
        geom_segment(aes(x = peaks_center, xend = peaks_center, y = 0, yend = RPKM))+
        xlim(-3000,genes_info$end[i_gene]-genes_info$start[i_gene])+xlab(sample_names[i_sample])+
        ggtitle(genes_info$Symbol[i_gene])+theme_minimal()+theme(plot.title = element_text(hjust=0.5))
    }
    else
    {
      new_plot <- ggplot(peaks, aes(x = peaks_center, y = RPKM))+ 
        geom_segment(aes(x = peaks_center, xend = peaks_center, y = 0, yend = RPKM))+
        xlim(-3000,genes_info$end[i_gene]-genes_info$start[i_gene])+xlab(sample_names[i_sample])+theme_minimal()
      gene_plot <- gene_plot/new_plot
    }
  }
  gene_plot <- gene_plot & ylim(0,max_peak_height) # For easy visual comparison, set all subplots to have the same ylim
  bmp(file = sprintf('Peak visualization of %s.bmp',genes_info$Symbol[i_gene]),width = 1000,height = 1000,units = "px",res = 150)
  print(gene_plot)
  dev.off()
}
################################################################################################################################
# Visualize HOMER and MEME enrichment results and picks motif of interest by elbow plot
setwd('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\ATACseq_SPINK1_KO\\Results')
HOMER_up_report <- homerkit::read_homerResults_html('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\ATACseq_SPINK1_KO\\Results\\bed files for HOMER motif analysis\\SPINK1_KO_up_motifs_HOMER\\homerResults.html')
HOMER_down_report <- homerkit::read_homerResults_html('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\ATACseq_SPINK1_KO\\Results\\bed files for HOMER motif analysis\\SPINK1_KO_down_motifs_HOMER\\homerResults.html')
MEME_up_report <- read_tsv('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\ATACseq_SPINK1_KO\\Results\\bed files for MEME motif analysis\\centrimo_MEME_up.tsv')
MEME_down_report <- read_tsv('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\ATACseq_SPINK1_KO\\Results\\bed files for MEME motif analysis\\centrimo_MEME_down.tsv')
elbow_plot_data_HOMER <- data.frame(group = factor(c(rep('upregulated motif',dim(HOMER_up_report)[1]),rep('downregulated motif',dim(HOMER_down_report)[1])),levels = c('upregulated motif','downregulated motif')),
                                    rank = c(1:dim(HOMER_up_report)[1],1:dim(HOMER_down_report)[1]),
                                    minus_log_p = -c(HOMER_up_report$log_pvalue,HOMER_down_report$log_pvalue))
elbow_plot_data_MEME <- data.frame(group = factor(c(rep('upregulated motif',dim(MEME_up_report)[1]),rep('downregulated motif',dim(MEME_down_report)[1])),levels = c('upregulated motif','downregulated motif')),
                                   rank = c(1:dim(MEME_up_report)[1],1:dim(MEME_down_report)[1]),
                                   minus_log_p = -c(MEME_up_report$`log_adj_p-value`,MEME_down_report$`log_adj_p-value`))
bmp(file = 'Elbow_plot_for_HOMER_motif_enrichment.bmp',width = 500,height = 500,units = "px",res = 150)
ggplot(elbow_plot_data_HOMER, aes(x = rank,y = minus_log_p,color = group))+geom_line()+geom_point()+theme_classic()+
  theme(legend.position = c(0.75, 0.8),legend.title = element_blank())+ylab('-log of p value')
dev.off()
bmp(file = 'Elbow_plot_for_MEME_motif_enrichment.bmp',width = 500,height = 500,units = "px",res = 150)
ggplot(elbow_plot_data_MEME, aes(x = rank,y = minus_log_p,color = group))+geom_line()+geom_point()+theme_classic()+
  theme(legend.position = c(0.75, 0.8),legend.title = element_blank())+ylab('-log of p value')
dev.off()
# Analyze peaks assigned to motif of interest by integrating with RNAseq data
KO_deseq_results <- read.csv('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\DEseq_SPINK1_KO_v2\\Results\\DEseq_results.csv')
peaks_motif_1 <- read.table('output_motif1_up_to_peaks.txt',sep = '\t',header = TRUE)
peaks_motif_2 <- read.table('output_motif2_up_to_peaks.txt',sep = '\t',header = TRUE)
peaks_motif_1_info <-  diffbind_sites[as.character(unique(peaks_motif_1$PositionID)),]
peaks_motif_2_info <-  diffbind_sites[as.character(unique(peaks_motif_2$PositionID)),]
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peaks_motif_info_1_annotated <- annotatePeak(peaks_motif_1_info,tssRegion = c(-1000, 1000),TxDb = txdb, annoDb = "org.Hs.eg.db")
peaks_motif_info_2_annotated <- annotatePeak(peaks_motif_2_info,tssRegion = c(-1000, 1000),TxDb = txdb, annoDb = "org.Hs.eg.db")
plotAnnoPie(peaks_motif_info_1_annotated)
plotAnnoPie(peaks_motif_info_2_annotated)
peaks_motif_info_1_annotated_df <- as.data.frame(peaks_motif_info_1_annotated)
peaks_motif_info_2_annotated_df <- as.data.frame(peaks_motif_info_2_annotated)
peaks_motif_info_1_annotated_subset_df <- peaks_motif_info_1_annotated_df[peaks_motif_info_1_annotated_df$annotation %in% c('Promoter','Distal Intergenic'),]
peaks_motif_info_2_annotated_subset_df <- peaks_motif_info_2_annotated_df[peaks_motif_info_2_annotated_df$annotation %in% c('Promoter','Distal Intergenic'),]
motif_1_DEGs <- KO_deseq_DEGs[KO_deseq_DEGs$gene_symbols %in% unique(peaks_motif_info_1_annotated_subset_df$SYMBOL),]
motif_2_DEGs <- KO_deseq_DEGs[KO_deseq_DEGs$gene_symbols %in% unique(peaks_motif_info_2_annotated_subset_df$SYMBOL),]
motif_1_DEGs <- motif_1_DEGs[abs(motif_1_DEGs$log2FoldChange)>=1,]
motif_2_DEGs <- motif_2_DEGs[abs(motif_2_DEGs$log2FoldChange)>=1,]
VennDiagram::venn.diagram(list(unique(peaks_motif_info_1_annotated_subset_df$SYMBOL),
                               unique(peaks_motif_info_2_annotated_subset_df$SYMBOL)),
                          category.names = c('',''),fill = c('lightgreen','lightcoral'),lty = 'blank',
                          cex = 2,fontface = 'bold',
                          filename = 'genes_overlap_between_AP1_green_FOX_red_motif.png',output = TRUE)
VennDiagram::venn.diagram(list(unique(motif_1_DEGs$gene_symbols),
                               unique(motif_2_DEGs$gene_symbols)),
                          category.names = c('',''),fill = c('lightgreen','lightcoral'),lty = 'blank',
                          cex = 2,fontface = 'bold',
                          filename = 'DEGs_overlap_between_AP1_green_FOX_red_motif.png',output = TRUE)
motif_1_2_overlapped <- intersect(peaks_motif_info_1_annotated_subset_df$SYMBOL,peaks_motif_info_2_annotated_subset_df$SYMBOL)
motif_1_2_overlapped_DEGs <- KO_deseq_DEGs[KO_deseq_DEGs$gene_symbols %in% unique(motif_1_2_overlapped),]
motif_1_2_overlapped_DEGs <- motif_1_2_overlapped_DEGs[order(motif_1_2_overlapped_DEGs$log2FoldChange,decreasing = TRUE),]
# Narrows down target from family level to individual protein level by analyzing RNAseq data
family_cts_plot <- function(family_to_plot,RNAseq_cts)
{
  family_cts <- RNAseq_cts[RNAseq_cts$gene_symbols %in% family_to_plot,2:6]
  family_cts['Scramble_mean'] <- rowMeans(family_cts[c('Scramble_1','Scramble_2')])
  family_cts['KO_mean'] <- rowMeans(family_cts[c('KO_1','KO_2')])
  family_cts_long <- pivot_longer(family_cts[c('gene_symbols','Scramble_mean','KO_mean')],!gene_symbols,
                                  names_to = 'condition',values_to = 'mean_counts')
  family_cts_long$condition <- factor(family_cts_long$condition,levels = c('KO_mean','Scramble_mean'))
  family_cts_long$gene_symbols <- factor(family_cts_long$gene_symbols,levels = rev(family_to_plot))
  ggplot(family_cts_long,aes(x = gene_symbols,y = mean_counts, fill = condition))+
    geom_bar(position = 'dodge',stat = 'identity')+coord_flip()+
    scale_fill_manual(values = c('Scramble_mean' = "lightgreen", 'KO_mean' = 'lightcoral'))+
    theme_minimal()+theme(legend.position = 'none')+xlab('')+ylab('')
}
normalized_cts <- read.csv('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\DEseq_SPINK1_KO_v2\\Results\\normalized_counts_matrix.csv')
family_list <- vector('list', length = 5)
names(family_list) <- c('JUN','FOS','MAF','ATF','FOX')
family_list$JUN <- c('JUN','JUNB','JUND')
family_list$FOS <- c('FOS','FOSB','FOSL1','FOSL2')
family_list$MAF <- c('MAF','MAFA','MAFB','MAFF','MAFG','MAFK','NRL')
family_list$ATF <- c(paste0('ATF',1:5),'ATF6','ATF6B','ATF7','BATF','BATF2','BATF3','JDP1','JDP2')
family_list$FOX <- sort(normalized_cts$gene_symbols[grep('^FOX',normalized_cts$gene_symbols)])
family_plots <- family_list
for (n in 1:length(family_list))
{
  family_plots[[n]] <- family_cts_plot(family_list[[n]],RNAseq_cts = normalized_cts)
}
combined_plot <- plot_spacer()+family_plots[[1]]+plot_spacer()+family_plots[[2]]+plot_spacer()+family_plots[[3]]+plot_spacer()+family_plots[[4]]+plot_layout(ncol = 1, heights = c(2,3,2,4,2,7,2,13))
combined_plot <- combined_plot|family_plots[[5]]
bmp(file = 'AP1_FOX_family_overview.bmp',width = 1200,height = 1500,units = "px",res = 150)
print(combined_plot)
dev.off()
CCLE <- read.csv('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\CCLE_RNA_seq\\OmicsExpressionProteinCodingGenesTPMLogp1.csv')
CCLE_subset <- CCLE[CCLE$X %in% c('ACH-000164','ACH-000601','ACH-000535','ACH-000222','ACH-000138','ACH-000107','ACH-000094','ACH-000354'),]
row.names(CCLE_subset) <- CCLE_subset$X
CCLE_subset <- CCLE_subset[c('ACH-000164','ACH-000601','ACH-000535','ACH-000222','ACH-000138','ACH-000107','ACH-000094','ACH-000354'),]
CCLE_subset$X <- c('PANC1','MiaPaCa2','BxPc3','AsPC1','CFPAC1','CAPAN2','HPAF2','CAPAN1')
colnames(CCLE_subset) <- gsub('\\.\\..*','', colnames(CCLE_subset))
write.csv(CCLE_subset,'C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\CCLE_RNA_seq\\OmicsExpressionProteinCodingGenesTPMLogp1_PDAC_subset.csv')
genes_to_plot <- c('KRT19','SOX9','MMP7','TSPAN8','VTCN1','REG1A')#c('CDH1','CDH2','VIM','ZEB1','ELF3','KLF5','SPINK1','COL18A1')
CCLE_subset_to_plot <- CCLE_subset[,c('X',genes_to_plot)]
CCLE_subset_to_plot <- pivot_longer(CCLE_subset_to_plot,!X,names_to = 'Gene',values_to = 'log1p_RSEM')
CCLE_subset_to_plot$X <- factor(CCLE_subset_to_plot$X,levels = c('PANC1','MiaPaCa2','BxPc3','AsPC1','CFPAC1','CAPAN2','HPAF2','CAPAN1'))
CCLE_subset_to_plot$Gene <- factor(CCLE_subset_to_plot$Gene,levels = genes_to_plot)
bmp('Distribution of genes of interest set 2 in different PDAC cell lines.bmp',width = 1000,height=1700,units='px',res=150)
ggplot(CCLE_subset_to_plot,aes(x = Gene, y = log1p_RSEM, fill = X))+xlab('')+ylab('log1p RSEM')+
  geom_bar(position = 'dodge',stat = 'identity')+coord_flip()+
  scale_fill_manual('',values = c('coral','lightcoral','gold3','gold1','palegreen','palegreen2','palegreen3','palegreen4'))
dev.off()  
