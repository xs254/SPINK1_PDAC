library(csaw)
library(edgeR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
setwd('C:/Users/xshan/OneDrive/Desktop/R code/CutNTag_SPINK1_KO')
# For H3K4Me3 data
target_name <- 'H3K4Me3'
bam_files <- c(sprintf('C:/Users/xshan/OneDrive/Desktop/R code/CutNTag_SPINK1_KO/data/ko_%d_k4_bowtie2_sorted.bam',1:4),
               sprintf('C:/Users/xshan/OneDrive/Desktop/R code/CutNTag_SPINK1_KO/data/scr_%d_k4_bowtie2_sorted.bam',1:4))
condition <- factor(c(rep('KO',4),rep('SCR',4)))
condition <- relevel(condition, ref = 'SCR')
design <- model.matrix(~condition)
param <- readParam(minq = 30,pe = 'both')
# Specify key parameters
win_width <- 150
bin_width <- 2000
merge_width <- 1000
# Counting for windows
win_counts <- windowCounts(bam_files,width = win_width,param = param)
# Counting for binned regions for normalization and filtering
bin_counts <- windowCounts(bam_files,width = bin_width, param=param)
# Filter out low abundance windows that are <5X enriched over background
win_to_keep <- filterWindowsGlobal(win_counts,bin_counts)$filter>log2(5)
win_counts <- win_counts[win_to_keep,]
# Normalization
win_counts <- normFactors(bin_counts,se.out=win_counts)
saveRDS(win_counts,sprintf('results/%s_csaw_normalized_counts.RDS',target_name))
# PCA plot for QC
win_counts <- readRDS(sprintf('results/%s_csaw_normalized_counts.RDS',target_name))
y <- asDGEList(win_counts)
y_logCPM <- cpm(y, log = TRUE, prior.count = 2)
y_pca <- prcomp(t(y_logCPM), scale. = TRUE)
pca_df <- data.frame(PC1 = y_pca$x[,1],PC2 = y_pca$x[,2],
                     condition = factor(condition,levels = c('KO','SCR')))
png(sprintf('results/PCA_plot_%s.png',target_name),
    width = 2000,height = 1500,units = 'px',res = 300)
ggplot(pca_df, aes(PC1, PC2,color = condition))+geom_point(size = 4)+
  labs(title = target_name,color = NULL)+
  xlab(sprintf('PC1 %d%%',round(100 * summary(y_pca)$importance[2,1])))+
  ylab(sprintf('PC2 %d%%',round(100 * summary(y_pca)$importance[2,2])))+
  theme_classic()+
  theme(axis.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.text = element_text(size = 15),
        plot.title = element_text(size = 20,hjust = 0.5))
dev.off()
# Identifying DB regions by edgeR
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design,robust = TRUE)
results <- glmQLFTest(fit,coef = 2)
# Merge significant windows to regions
merged <- mergeResults(win_counts,results$table,tol = merge_width)
# Annotate regions
regions <- merged$regions
mcols(regions) <- cbind(mcols(regions),merged$combined)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
regions_annotated <- annotatePeak(regions,TxDb = txdb,annoDb = 'org.Hs.eg.db',tssRegion = c(-1000, 1000))
saveRDS(regions_annotated,sprintf('results/%s_annotated_edger_results.rds',target_name))
regions_annotated <- readRDS(sprintf('results/%s_annotated_edger_results.rds',target_name))
regions_annotated_df <- as.data.frame(regions_annotated)
write.csv(regions_annotated_df,
          sprintf('results/%s_csaw_edgeR_annotated_results.csv',target_name),row.names = FALSE)
# Visualization of statistically significant regions
regions_filtered <- regions[regions$FDR<=0.05]
regions_filtered_annotated <- annotatePeak(regions,TxDb = txdb,annoDb = 'org.Hs.eg.db',tssRegion = c(-1000, 1000))
# Plot distance to TSS
png(sprintf('results/%s_distance_of_diff_peaks_to_TSS.png',target_name),
    width = 2000,height = 600,units = 'px',res = 300)
plotDistToTSS(regions_filtered_annotated,title = target_name)
dev.off()
# Pie chart of diff peaks location
png(sprintf('results/%s_distribution_of_diff_peak_resgions.png',target_name),
    width = 2000,height = 1000,units = 'px',res = 300)
plotAnnoPie(regions_filtered_annotated)
dev.off()

# For H3K27Me3 data
target_name <- 'H3K27Me3'
bam_files <- c(sprintf('C:/Users/xshan/OneDrive/Desktop/R code/CutNTag_SPINK1_KO/data/ko_%d_27_bowtie2_sorted.bam',1:4),
               sprintf('C:/Users/xshan/OneDrive/Desktop/R code/CutNTag_SPINK1_KO/data/scr_%d_27_bowtie2_sorted.bam',2:4))
condition <- factor(c(rep('KO',4),rep('SCR',3)))
condition <- relevel(condition, ref = 'SCR')
design <- model.matrix(~condition)
param <- readParam(minq = 30,pe = 'both')
# Specify key parameters
win_width <- 2000
bin_width <- 20000
merge_width <- 1000
# Counting for windows
win_counts <- windowCounts(bam_files,width = win_width,param = param)
# Counting for binned regions for normalization and filtering
bin_counts <- windowCounts(bam_files,width = bin_width, param=param)
# Filter out low abundance windows that are <5X enriched over background
win_to_keep <- filterWindowsGlobal(win_counts,bin_counts)$filter>log2(5)
win_counts <- win_counts[win_to_keep,]
# Normalization
win_counts <- normFactors(bin_counts,se.out=win_counts)
saveRDS(win_counts,sprintf('results/%s_csaw_normalized_counts.RDS',target_name))
# PCA plot for QC
win_counts <- readRDS(sprintf('results/%s_csaw_normalized_counts.RDS',target_name))
y <- asDGEList(win_counts)
y_logCPM <- cpm(y, log = TRUE, prior.count = 2)
y_pca <- prcomp(t(y_logCPM), scale. = TRUE)
pca_df <- data.frame(PC1 = y_pca$x[,1],PC2 = y_pca$x[,2],
                     condition = factor(condition,levels = c('KO','SCR')))
png(sprintf('results/PCA_plot_%s.png',target_name),
    width = 2000,height = 1500,units = 'px',res = 300)
ggplot(pca_df, aes(PC1, PC2,color = condition))+geom_point(size = 4)+
  labs(title = target_name,color = NULL)+
  xlab(sprintf('PC1 %d%%',round(100 * summary(y_pca)$importance[2,1])))+
  ylab(sprintf('PC2 %d%%',round(100 * summary(y_pca)$importance[2,2])))+
  theme_classic()+
  theme(axis.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.text = element_text(size = 15),
        plot.title = element_text(size = 20,hjust = 0.5))
dev.off()
# Identifying DB regions by edgeR
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design,robust = TRUE)
results <- glmQLFTest(fit,coef = 2)
# Merge significant windows to regions
merged <- mergeResults(win_counts,results$table,tol = merge_width)
# Annotate the regions
regions <- merged$regions
mcols(regions) <- cbind(mcols(regions),merged$combined)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
regions_annotated <- annotatePeak(regions,TxDb = txdb,annoDb = 'org.Hs.eg.db',tssRegion = c(-1000, 1000))
saveRDS(regions_annotated,sprintf('results/%s_annotated_edger_results.rds',target_name))
regions_annotated <- readRDS(sprintf('results/%s_annotated_edger_results.rds',target_name))
regions_annotated_df <- as.data.frame(regions_annotated)
write.csv(regions_annotated_df,
          sprintf('results/%s_csaw_edgeR_annotated_results.csv',target_name),row.names = FALSE)
# Visualization of statistically significant regions
regions_filtered <- regions[regions$FDR<=0.05]
regions_filtered_annotated <- annotatePeak(regions,TxDb = txdb,annoDb = 'org.Hs.eg.db',tssRegion = c(-1000, 1000))
# Plot distance to TSS
png(sprintf('results/%s_distance_of_diff_peaks_to_TSS.png',target_name),
    width = 2000,height = 600,units = 'px',res = 300)
plotDistToTSS(regions_filtered_annotated,title = target_name)
dev.off()
# Pie chart of diff peaks location
png(sprintf('results/%s_distribution_of_diff_peak_resgions.png',target_name),
    width = 2000,height = 1000,units = 'px',res = 300)
plotAnnoPie(regions_filtered_annotated)
dev.off()
