library(csaw)
library(Gviz)
library(edgeR)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
setwd('C:/Users/xshan/OneDrive/Desktop/R code/CutNTag_SPINK1_KO')
# Set bin size for each sequencing
RNA_bin_size <- 150
ATAC_bin_size <- 150
H3K4Me3_bin_size <-150
H3K27Me3_bin_size <- 150
# Count BAM based on bin size
## For RNA seq
param <- readParam(minq = 30)
bam_files <- c(sprintf('C:/Users/xshan/OneDrive/Desktop/R code/DEseq_SPINK1_KO_v2/data/TM59Z7_bam/TM59Z7_%d_KO_%d_dedup-mapped-reads.bam',6:10,1:5),
               sprintf('C:/Users/xshan/OneDrive/Desktop/R code/DEseq_SPINK1_KO_v2/data/TM59Z7_bam/TM59Z7_%d_SCR_%d_dedup-mapped-reads.bam',1:5,1:5))
RNA_KO_win_counts <- windowCounts(bam_files,width = RNA_bin_size,spacing = RNA_bin_size,param = param)
bam_files <- c(sprintf('C:/Users/xshan/OneDrive/Desktop/R code/DEseq_SPINK1_KO_v2/data/TM59Z7_bam/TM59Z7_%d_OE_%d_dedup-mapped-reads.bam',16:20,1:5),
               sprintf('C:/Users/xshan/OneDrive/Desktop/R code/DEseq_SPINK1_KO_v2/data/TM59Z7_bam/TM59Z7_%d_CONT_%d_dedup-mapped-reads.bam',11:15,1:5))
RNA_OE_win_counts <- windowCounts(bam_files,width = RNA_bin_size,spacing = RNA_bin_size,param = param)
## For ATAC
param <- readParam(minq = 30,pe = 'both')
bam_files <- c('C:/Users/xshan/OneDrive/Desktop/R code/ATACseq_SPINK1_KO/02_0FJV_025ZYale_2g-6-SP-P-4_ATAC_hs_i21-1_peak_calling_ready.bam',
               'C:/Users/xshan/OneDrive/Desktop/R code/ATACseq_SPINK1_KO/06_0FJZ_025ZYale_2G-6-SP_ATAC_hs_i25-521_peak_calling_ready.bam',
               'C:/Users/xshan/OneDrive/Desktop/R code/ATACseq_SPINK1_KO/01_0FJU_025ZYale_SCR-P-6_ATAC_hs_i20-48_peak_calling_ready.bam',
               'C:/Users/xshan/OneDrive/Desktop/R code/ATACseq_SPINK1_KO/05_0FJY_025ZYale_SCR_ATAC_hs_i24-4_peak_calling_ready.bam')
ATAC_win_counts <- windowCounts(bam_files,width = ATAC_bin_size,spacing = ATAC_bin_size,param = param)
## For H3K4Me3
param <- readParam(minq = 30,pe = 'both')
bam_files <- c(sprintf('C:/Users/xshan/OneDrive/Desktop/R code/CutNTag_SPINK1_KO/data/ko_%d_k4_bowtie2_sorted.bam',1:4),
               sprintf('C:/Users/xshan/OneDrive/Desktop/R code/CutNTag_SPINK1_KO/data/scr_%d_k4_bowtie2_sorted.bam',1:4))
H3K4Me3_win_counts <- windowCounts(bam_files,width = H3K27Me3_bin_size,spacing = H3K27Me3_bin_size,param = param)
## For H3K27Me3
param <- readParam(minq = 30,pe = 'both')
bam_files <- c(sprintf('C:/Users/xshan/OneDrive/Desktop/R code/CutNTag_SPINK1_KO/data/ko_%d_27_bowtie2_sorted.bam',1:4),
               sprintf('C:/Users/xshan/OneDrive/Desktop/R code/CutNTag_SPINK1_KO/data/scr_%d_27_bowtie2_sorted.bam',2:4))
H3K27Me3_win_counts <- windowCounts(bam_files,width = H3K4Me3_bin_size,spacing = H3K4Me3_bin_size,param = param)
# Create grange with aggregated mean of KO and SCR
## For RNA
### For RNA KO
RNA_KO_cpm <- asDGEList(RNA_KO_win_counts)
RNA_KO_cpm <- cpm(RNA_KO_win_counts,log = FALSE)
RNA_KO_gr <- rowRanges(RNA_KO_win_counts)
mcols(RNA_KO_gr)$KO_mean_CPM <- rowMeans(RNA_KO_cpm[,1:5])
mcols(RNA_KO_gr)$SCR_mean_CPM <- rowMeans(RNA_KO_cpm[,6:10])
saveRDS(RNA_KO_gr,'results/integrated_visualization/RNA_KO_gr.rds')
### For RNA OE
RNA_OE_cpm <- asDGEList(RNA_OE_win_counts)
RNA_OE_cpm <- cpm(RNA_OE_win_counts,log = FALSE)
RNA_OE_gr <- rowRanges(RNA_OE_win_counts)
mcols(RNA_OE_gr)$OE_mean_CPM <- rowMeans(RNA_OE_cpm[,1:5])
mcols(RNA_OE_gr)$CONT_mean_CPM <- rowMeans(RNA_OE_cpm[,6:10])
saveRDS(RNA_OE_gr,'results/integrated_visualization/RNA_OE_gr.rds')
## For ATAC
ATAC_cpm <- asDGEList(ATAC_win_counts)
ATAC_cpm <- cpm(ATAC_win_counts,log = FALSE)
ATAC_gr <- rowRanges(ATAC_win_counts)
mcols(ATAC_gr)$KO_mean_CPM <- rowMeans(ATAC_cpm[,1:2])
mcols(ATAC_gr)$SCR_mean_CPM <- rowMeans(ATAC_cpm[,3:4])
saveRDS(ATAC_gr,'results/integrated_visualization/ATAC_gr.rds')
## For H3K4Me3
H3K4Me3_cpm <- asDGEList(H3K4Me3_win_counts)
H3K4Me3_cpm <- cpm(H3K4Me3_win_counts,log = FALSE)
H3K4Me3_gr <- rowRanges(H3K4Me3_win_counts)
mcols(H3K4Me3_gr)$KO_mean_CPM <- rowMeans(H3K4Me3_cpm[,1:4])
mcols(H3K4Me3_gr)$SCR_mean_CPM <- rowMeans(H3K4Me3_cpm[,5:8])
saveRDS(H3K4Me3_gr,'results/integrated_visualization/H3K4Me3_gr.rds')
## For H3K27Me3
H3K27Me3_cpm <- asDGEList(H3K27Me3_win_counts)
H3K27Me3_cpm <- cpm(H3K27Me3_win_counts,log = FALSE)
H3K27Me3_gr <- rowRanges(H3K27Me3_win_counts)
mcols(H3K27Me3_gr)$KO_mean_CPM <- rowMeans(H3K27Me3_cpm[,1:4])
mcols(H3K27Me3_gr)$SCR_mean_CPM <- rowMeans(H3K27Me3_cpm[,5:7])
saveRDS(H3K27Me3_gr,'results/integrated_visualization/H3K27Me3_gr.rds')
########################### End of preparing gene range files##########################
setwd('C:/Users/xshan/OneDrive/Desktop/R code/CutNTag_SPINK1_KO')
RNA_KO_gr <- readRDS('results/integrated_visualization/RNA_KO_gr.rds')
RNA_OE_gr <- readRDS('results/integrated_visualization/RNA_OE_gr.rds')
ATAC_gr <- readRDS('results/integrated_visualization/ATAC_gr.rds')
H3K4Me3_gr <- readRDS('results/integrated_visualization/H3K4Me3_gr.rds')
H3K27Me3_gr <- readRDS('results/integrated_visualization/H3K27Me3_gr.rds')
# Add chr to RNA related gr seqnames to align with other gr
seqlevelsStyle(RNA_KO_gr) <- seqlevelsStyle(ATAC_gr)
seqlevelsStyle(RNA_OE_gr) <- seqlevelsStyle(ATAC_gr)
# Remove levels not started with chr
chr_levels <- grep('^chr', seqlevels(RNA_KO_gr), value = TRUE)
RNA_KO_gr <- keepSeqlevels(RNA_KO_gr,chr_levels,pruning.mode = 'coarse')
chr_levels <- grep('^chr', seqlevels(RNA_OE_gr), value = TRUE)
RNA_OE_gr <- keepSeqlevels(RNA_OE_gr,chr_levels,pruning.mode = 'coarse')
# Set reference transcriptome
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
ex_by_tx <- exonsBy(txdb, by='tx', use.names=TRUE)
genes_gr <- genes(txdb)
# Get gene(s) of interest to plot
GOI_collection_symbol <- readLines('results/stem_related_genes_from_enrichr.txt')
GOI_collection_symbol <- c(GOI_collection_symbol,'CDH1')
GOI_collection_id <- mapIds(org.Hs.eg.db,keys = GOI_collection_symbol,
                            column = 'ENTREZID',keytype = 'SYMBOL',multiVals = 'first')
# Get genome coordinate of a specific gene
for(i_GOI in 1:length(GOI_collection_id))
{
  GOI <- GOI_collection_id[i_GOI]
  GOI_tx_ref <- 'ENST00000336868.8' #Specify the transcript reference to use
  GOI_coord <- genes_gr[names(genes_gr) == GOI]
  GOI_chr <- as.character(seqnames(GOI_coord))
  # Define TSS
  if(length(GOI_coord)==0){next}
  if(as.character(strand(GOI_coord))=='+')
  {
    GOI_TSS <- start(GOI_coord)
  }else{
    GOI_TSS <- end(GOI_coord)  
  }
  # Calculate limit of y axis
  GOI_query <- GRanges(seqnames = seqnames(GOI_coord),
                       ranges = IRanges(start = start(GOI_coord)-1000,end = end(GOI_coord)))
  GOI_gr <- subsetByOverlaps(H3K4Me3_gr, GOI_query)
  H3K4Me3_ylim <- ceiling(max(mcols(GOI_gr)$KO_mean_CPM,mcols(GOI_gr)$SCR_mean_CPM))
  H3K4Me3_ylim <- max(H3K4Me3_ylim,1)
  GOI_gr <- subsetByOverlaps(H3K27Me3_gr, GOI_query)
  H3K27Me3_ylim <- ceiling(max(mcols(GOI_gr)$KO_mean_CPM,mcols(GOI_gr)$SCR_mean_CPM))
  H3K27Me3_ylim <- max(H3K27Me3_ylim,1)
  GOI_gr <- subsetByOverlaps(ATAC_gr, GOI_query)
  ATAC_ylim <- ceiling(max(mcols(GOI_gr)$KO_mean_CPM,mcols(GOI_gr)$SCR_mean_CPM))
  ATAC_ylim <- max(ATAC_ylim,1)
  GOI_gr <- subsetByOverlaps(RNA_KO_gr, GOI_query)
  RNA_KO_ylim <- ceiling(max(mcols(GOI_gr)$KO_mean_CPM,mcols(GOI_gr)$SCR_mean_CPM))
  RNA_KO_ylim <- max(RNA_KO_ylim,1)
  GOI_gr <- subsetByOverlaps(RNA_OE_gr, GOI_query)
  RNA_OE_ylim <- ceiling(max(mcols(GOI_gr)$OE_mean_CPM,mcols(GOI_gr)$CONT_mean_CPM))
  RNA_OE_ylim <- max(RNA_OE_ylim,1)
  # Visualization by gene track
  track_colors <- c('orange2','gold2','palegreen2','skyblue2','orchid2')
  #gtrack <- GeneRegionTrack(ex_by_tx[GOI_tx_ref],fill = 'black')
  gtrack <- GenomeAxisTrack()
  H3K4Me3_KO_track <- DataTrack(range = H3K4Me3_gr,data = 'KO_mean_CPM',
                                type = 'histogram',name = 'KO',
                                ylim = c(0,H3K4Me3_ylim),yTicksAt = c(0,H3K4Me3_ylim/2,H3K4Me3_ylim),
                                col.histogram = track_colors[1],fill = track_colors[1],col.axis = track_colors[1])
  H3K4Me3_SCR_track <- DataTrack(range = H3K4Me3_gr,data = 'SCR_mean_CPM',
                                 type = 'histogram',name = 'SCR',
                                 ylim = c(0,H3K4Me3_ylim),yTicksAt = c(0,H3K4Me3_ylim/2,H3K4Me3_ylim),
                                 col.histogram = track_colors[1],fill = track_colors[1],col.axis = track_colors[1])
  H3K27Me3_KO_track <- DataTrack(range = H3K27Me3_gr,data = 'KO_mean_CPM',
                                 type = 'histogram',name = 'KO',
                                 ylim = c(0,H3K27Me3_ylim),yTicksAt = c(0,H3K27Me3_ylim/2,H3K27Me3_ylim),
                                 col.histogram = track_colors[2],fill = track_colors[2],col.axis = track_colors[2])
  H3K27Me3_SCR_track <- DataTrack(range = H3K27Me3_gr,data = 'SCR_mean_CPM',
                                  type = 'histogram',name = 'SCR',
                                  ylim = c(0,H3K27Me3_ylim),yTicksAt = c(0,H3K27Me3_ylim/2,H3K27Me3_ylim),
                                  col.histogram = track_colors[2],fill = track_colors[2],col.axis = track_colors[2])
  ATAC_KO_track <- DataTrack(range = ATAC_gr,data = 'KO_mean_CPM',
                             type = 'histogram',name = 'KO',
                             ylim = c(0,ATAC_ylim),yTicksAt = c(0,ATAC_ylim/2,ATAC_ylim),
                             col.histogram = track_colors[3],fill = track_colors[3],col.axis = track_colors[3])
  ATAC_SCR_track <- DataTrack(range = ATAC_gr,data = 'SCR_mean_CPM',
                              type = 'histogram',name = 'SCR',
                              ylim = c(0,ATAC_ylim),yTicksAt = c(0,ATAC_ylim/2,ATAC_ylim),
                              col.histogram = track_colors[3],fill = track_colors[3],col.axis = track_colors[3])
  RNA_KO_track <- DataTrack(range = RNA_KO_gr,data = 'KO_mean_CPM',
                            type = 'histogram',name = 'KO',
                            ylim = c(0,RNA_KO_ylim),yTicksAt = c(0,RNA_KO_ylim/2,RNA_KO_ylim),
                            col.histogram = track_colors[4],fill = track_colors[4],col.axis = track_colors[4])
  RNA_SCR_track <- DataTrack(range = RNA_KO_gr,data = 'SCR_mean_CPM',
                             type = 'histogram',name = 'SCR',
                             ylim = c(0,RNA_KO_ylim),yTicksAt = c(0,RNA_KO_ylim/2,RNA_KO_ylim),
                             col.histogram = track_colors[4],fill = track_colors[4],col.axis = track_colors[4])
  RNA_OE_track <- DataTrack(range = RNA_OE_gr,data = 'OE_mean_CPM',
                            type = 'histogram',name = 'OE',
                            ylim = c(0,RNA_OE_ylim),yTicksAt = c(0,RNA_OE_ylim/2,RNA_OE_ylim),
                            col.histogram = track_colors[5],fill = track_colors[5],col.axis = track_colors[5])
  RNA_CONT_track <- DataTrack(range = RNA_OE_gr,data = 'CONT_mean_CPM',
                              type = 'histogram',name = 'CONT',
                              ylim = c(0,RNA_OE_ylim),yTicksAt = c(0,RNA_OE_ylim/2,RNA_OE_ylim),
                              col.histogram = track_colors[5],fill = track_colors[5],col.axis = track_colors[5])
  TSS_highlighted <- HighlightTrack(list(gtrack,H3K4Me3_SCR_track,H3K4Me3_KO_track,
                                         H3K27Me3_SCR_track,H3K27Me3_KO_track,
                                         ATAC_SCR_track,ATAC_KO_track,
                                         RNA_KO_track,RNA_SCR_track,
                                         RNA_OE_track,RNA_CONT_track),
                                    start = GOI_TSS-1000,end = GOI_TSS+1000,chromosome = GOI_chr)
  png(sprintf('results/integrated_visualization/plots/%s_integrated_visualization.png',names(GOI)),
      width = 2500,height = 7500,units = 'px',res = 300)
  plotTracks(TSS_highlighted,
             chromosome = GOI_chr,from = start(GOI_coord)-1000,to = end(GOI_coord)+1000,
             showTitle = FALSE,background.title = 'white',cex.axis = 1.8)
  dev.off()
}







GOI <- GOI_collection_id['CDH1']
GOI_tx_ref <- 'ENST00000261769.10' #Specify the transcript reference to use
GOI_coord <- genes_gr[names(genes_gr) == GOI]
GOI_chr <- as.character(seqnames(GOI_coord))
# Define TSS
if(length(GOI_coord)==0){next}
if(as.character(strand(GOI_coord))=='+')
{
  GOI_TSS <- start(GOI_coord)
}else{
  GOI_TSS <- end(GOI_coord)  
}
# Calculate limit of y axis
GOI_query <- GRanges(seqnames = seqnames(GOI_coord),
                     ranges = IRanges(start = start(GOI_coord)-1000,end = end(GOI_coord)))
GOI_gr <- subsetByOverlaps(H3K4Me3_gr, GOI_query)
H3K4Me3_ylim <- ceiling(max(mcols(GOI_gr)$KO_mean_CPM,mcols(GOI_gr)$SCR_mean_CPM))
H3K4Me3_ylim <- max(H3K4Me3_ylim,1)
GOI_gr <- subsetByOverlaps(H3K27Me3_gr, GOI_query)
H3K27Me3_ylim <- ceiling(max(mcols(GOI_gr)$KO_mean_CPM,mcols(GOI_gr)$SCR_mean_CPM))
H3K27Me3_ylim <- max(H3K27Me3_ylim,1)
GOI_gr <- subsetByOverlaps(ATAC_gr, GOI_query)
ATAC_ylim <- ceiling(max(mcols(GOI_gr)$KO_mean_CPM,mcols(GOI_gr)$SCR_mean_CPM))
ATAC_ylim <- max(ATAC_ylim,1)
GOI_gr <- subsetByOverlaps(RNA_KO_gr, GOI_query)
RNA_KO_ylim <- ceiling(max(mcols(GOI_gr)$KO_mean_CPM,mcols(GOI_gr)$SCR_mean_CPM))
RNA_KO_ylim <- max(RNA_KO_ylim,1)
GOI_gr <- subsetByOverlaps(RNA_OE_gr, GOI_query)
RNA_OE_ylim <- ceiling(max(mcols(GOI_gr)$OE_mean_CPM,mcols(GOI_gr)$CONT_mean_CPM))
RNA_OE_ylim <- max(RNA_OE_ylim,1)
# Visualization by gene track
track_colors <- c('orange2','gold2','palegreen2','skyblue2','orchid2')
gtrack <- GeneRegionTrack(ex_by_tx[GOI_tx_ref],fill = 'black')
H3K4Me3_KO_track <- DataTrack(range = H3K4Me3_gr,data = 'KO_mean_CPM',
                              type = 'histogram',name = 'KO',
                              ylim = c(0,H3K4Me3_ylim),yTicksAt = c(0,H3K4Me3_ylim/2,H3K4Me3_ylim),
                              col.histogram = track_colors[1],fill = track_colors[1],col.axis = track_colors[1])
H3K4Me3_SCR_track <- DataTrack(range = H3K4Me3_gr,data = 'SCR_mean_CPM',
                               type = 'histogram',name = 'SCR',
                               ylim = c(0,H3K4Me3_ylim),yTicksAt = c(0,H3K4Me3_ylim/2,H3K4Me3_ylim),
                               col.histogram = track_colors[1],fill = track_colors[1],col.axis = track_colors[1])
H3K27Me3_KO_track <- DataTrack(range = H3K27Me3_gr,data = 'KO_mean_CPM',
                               type = 'histogram',name = 'KO',
                               ylim = c(0,H3K27Me3_ylim),yTicksAt = c(0,H3K27Me3_ylim/2,H3K27Me3_ylim),
                               col.histogram = track_colors[2],fill = track_colors[2],col.axis = track_colors[2])
H3K27Me3_SCR_track <- DataTrack(range = H3K27Me3_gr,data = 'SCR_mean_CPM',
                                type = 'histogram',name = 'SCR',
                                ylim = c(0,H3K27Me3_ylim),yTicksAt = c(0,H3K27Me3_ylim/2,H3K27Me3_ylim),
                                col.histogram = track_colors[2],fill = track_colors[2],col.axis = track_colors[2])
ATAC_KO_track <- DataTrack(range = ATAC_gr,data = 'KO_mean_CPM',
                           type = 'histogram',name = 'KO',
                           ylim = c(0,ATAC_ylim),yTicksAt = c(0,ATAC_ylim/2,ATAC_ylim),
                           col.histogram = track_colors[3],fill = track_colors[3],col.axis = track_colors[3])
ATAC_SCR_track <- DataTrack(range = ATAC_gr,data = 'SCR_mean_CPM',
                            type = 'histogram',name = 'SCR',
                            ylim = c(0,ATAC_ylim),yTicksAt = c(0,ATAC_ylim/2,ATAC_ylim),
                            col.histogram = track_colors[3],fill = track_colors[3],col.axis = track_colors[3])
RNA_KO_track <- DataTrack(range = RNA_KO_gr,data = 'KO_mean_CPM',
                          type = 'histogram',name = 'KO',
                          ylim = c(0,RNA_KO_ylim),yTicksAt = c(0,RNA_KO_ylim/2,RNA_KO_ylim),
                          col.histogram = track_colors[4],fill = track_colors[4],col.axis = track_colors[4])
RNA_SCR_track <- DataTrack(range = RNA_KO_gr,data = 'SCR_mean_CPM',
                           type = 'histogram',name = 'SCR',
                           ylim = c(0,RNA_KO_ylim),yTicksAt = c(0,RNA_KO_ylim/2,RNA_KO_ylim),
                           col.histogram = track_colors[4],fill = track_colors[4],col.axis = track_colors[4])
RNA_OE_track <- DataTrack(range = RNA_OE_gr,data = 'OE_mean_CPM',
                          type = 'histogram',name = 'OE',
                          ylim = c(0,RNA_OE_ylim),yTicksAt = c(0,RNA_OE_ylim/2,RNA_OE_ylim),
                          col.histogram = track_colors[5],fill = track_colors[5],col.axis = track_colors[5])
RNA_CONT_track <- DataTrack(range = RNA_OE_gr,data = 'CONT_mean_CPM',
                            type = 'histogram',name = 'CONT',
                            ylim = c(0,RNA_OE_ylim),yTicksAt = c(0,RNA_OE_ylim/2,RNA_OE_ylim),
                            col.histogram = track_colors[5],fill = track_colors[5],col.axis = track_colors[5])
TSS_highlighted <- HighlightTrack(list(gtrack,
                                       H3K4Me3_KO_track,H3K4Me3_SCR_track,
                                       H3K27Me3_KO_track,H3K27Me3_SCR_track,
                                       ATAC_KO_track,ATAC_SCR_track,
                                       RNA_KO_track,RNA_SCR_track,
                                       RNA_OE_track,RNA_CONT_track),
                                  start = GOI_TSS-1000,end = GOI_TSS+1000,chromosome = GOI_chr,
                                  fill = '#fff992',col = 'white')
png(sprintf('results/integrated_visualization/plots/%s_integrated_visualization.png',names(GOI)),
    width = 3500,height = 7500,units = 'px',res = 300)
plotTracks(TSS_highlighted,
           chromosome = GOI_chr,from = start(GOI_coord)-1000,to = end(GOI_coord)+1000,
           showTitle = FALSE,background.title = 'white',cex.axis = 1.8)
dev.off()
