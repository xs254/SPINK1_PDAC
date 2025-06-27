library(Seurat)
library(ggplot2)
library(patchwork)
library(monocle3)
# Define helper functions
plot_gene_trajectory <- function(data_set,genes_to_plot,bin_size = 10,method_to_use = 'median',pseudotime_cap = NULL,x_lim = c(-5,85),output_width = 2000,output_height = 1000,violin_width = 0.5)
{
  to_plot <- FetchData(data_set,vars = c(genes_to_plot,'pseudotime'),layer = 'data')
  # Assign a small number to cells with pseudotime 0 such that during binning by ceiling, 
  # these cells can bin into 1 rather than form a separate bin 0
  to_plot$pseudotime[to_plot$pseudotime==0] <- 0.001
  to_plot_high <- to_plot[high_SPINK1_trajectory,]
  to_plot_low <- to_plot[low_SPINK1_trajectory,]
  if(!is.null(pseudotime_cap))
  {
    to_plot_high$pseudotime[to_plot_high$pseudotime>=pseudotime_cap[1]] <- pseudotime_cap[1]
    to_plot_low$pseudotime[to_plot_low$pseudotime>=pseudotime_cap[2]] <- pseudotime_cap[2]
  }
  to_plot_high$pseudotime_bin <- ceiling(to_plot_high$pseudotime/bin_size)
  to_plot_low$pseudotime_bin <- ceiling(to_plot_low$pseudotime/bin_size)
  for(gene in genes_to_plot)
  {
    high_stats <- aggregate(get(gene) ~ pseudotime_bin,to_plot_high,get(method_to_use))
    colnames(high_stats)[2] <- gene
    low_stats <- aggregate(get(gene) ~ pseudotime_bin,to_plot_low,get(method_to_use))
    colnames(low_stats)[2] <- gene
    # Calculate a consensus y axis upper limit for plotting
    y_upper <- max(max(to_plot_high[,gene],max(to_plot_low[,gene])))+0.5
    p1 <- ggplot(to_plot_high, aes(bin_size*pseudotime_bin-bin_size/2,get(gene),group = pseudotime_bin))+geom_violin(scale = 'width',width = violin_width)+
      geom_point(data = high_stats, aes(bin_size*pseudotime_bin-bin_size/2,get(gene)))+
      geom_line(data = high_stats, aes(bin_size*pseudotime_bin-bin_size/2,get(gene),group = 1))+xlim(x_lim[1],x_lim[2])+ylim(0,y_upper)+xlab('')+ylab(gene)+
      theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20))
    p2 <- ggplot(to_plot_low, aes(bin_size*pseudotime_bin-bin_size/2,get(gene),group = pseudotime_bin))+geom_violin(scale = 'width',width = violin_width)+
      geom_point(data = low_stats, aes(bin_size*pseudotime_bin-bin_size/2,get(gene)))+
      geom_line(data = low_stats, aes(bin_size*pseudotime_bin-bin_size/2,get(gene),group = 1))+xlim(x_lim[1],x_lim[2])+ylim(0,y_upper)+xlab('Pseudotime')+ylab(gene)+
      theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20))
    png(paste0(gene,' pseudotime trend two trajectories.png'),width = output_width,height = output_height, units = 'px', res = 300)
    print(p1/p2)
    dev.off()
  }
}
mean_for_log1p <- function(input)
{
  return(log1p(mean(exp(input)-1)))
}
setwd('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\reference_dataset_core')
# 1st way of integration: RPCA by project
PDAC_all_cells_RPCA <- readRDS('PDAC_all_cells_15430_RPCA_by_project.rds')
PDAC_all_cells_RPCA <- RunUMAP(PDAC_all_cells_RPCA, reduction = 'integrated_rpca',dims = 1:30)#30
setwd('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\final_results_v4')
png('PDAC_all_cells_RPCA_dim30_umap.png',width = 2000,height = 2000,res = 300)
DimPlot(PDAC_all_cells_RPCA, group.by = 'Cell_type',label = TRUE)+NoLegend()+NoAxes()
dev.off()
# One subcluster of cells originally identified as ductal cell type 2 and macrophages are likely neutrophil
png('neutrophil_markers.png',width = 4000,height = 4000,res = 300)
FeaturePlot(PDAC_all_cells_RPCA, features = c('FCGR3B','CSF3R','G0S2','CXCR2'))&NoAxes()&NoLegend()
dev.off()
# Clustering
PDAC_all_cells_RPCA <- FindNeighbors(PDAC_all_cells_RPCA,reduction = 'integrated_rpca',dims = 1:30)
PDAC_all_cells_RPCA <- FindClusters(PDAC_all_cells_RPCA, resolution = 0.5)
PDAC_all_cells_RPCA$seurat_clusters_v2 <- Idents(PDAC_all_cells_RPCA)
png('PDAC_all_cells_louvian_clustering.png',width = 2000,height = 2000,res = 300)
DimPlot(PDAC_all_cells_RPCA,group.by = 'seurat_clusters_v2',label = TRUE)+NoLegend()+NoAxes()
dev.off()
# Find markers for cluster 19 and 21.
# Although originally labelled as T cell, 
# cluster 19 is most likely antibody-secreting cells (ASC) due to high expression of XBP1,JCHAIN,MZB1,PRDM1
Idents(PDAC_all_cells_RPCA) <- PDAC_all_cells_RPCA$seurat_clusters_v2
cluster_19_markers <- FindMarkers(PDAC_all_cells_RPCA,ident.1 = 19)
cluster_19_markers <- cluster_19_markers[order(cluster_19_markers$avg_log2FC,decreasing = TRUE),]
cluster_19_markers_filtered <- cluster_19_markers[cluster_19_markers$pct.1>=0.5,]
# cluster 21 is most likely mast cells due to high expression of KIT, CPA3, TSAB1
cluster_21_markers <- FindMarkers(PDAC_all_cells_RPCA,ident.1 = 21)
cluster_21_markers <- cluster_21_markers[order(cluster_21_markers$avg_log2FC,decreasing = TRUE),]
cluster_21_markers_filtered <- cluster_21_markers[cluster_21_markers$pct.1>=0.5,]
# Plot markers of all types of cells
cell_type_markers <-  list(c('PRSS1','CTRB1','CTRB2','REG1B'),
                           c('XBP1','JCHAIN','MZB1','PRDM1'),
                           c('MS4A1','CD79A','CD79B','CD52'),
                           c('AMBP','CFTR','MMP7'),
                           c('KRT19','KRT7','TSPAN8','SLPI'),
                           c('CHGB','INS','IAPP','SCGN','ERO1B','RBP4','SLC30A8','TTR','SST','NEUROD1','PCSK1','PAX6'),
                           c('CDH5','PLVAP','VWF','CLDN5'),
                           c('LUM','DCN','COL1A1'),
                           c('AIF1','FCGR1A','CD14','CD68'),
                           c('KIT','TPSAB1','CPA3','HDC'),
                           c('FCGR3B','CSF3R','G0S2','CXCR2'),
                           c('RGS5','ACTA2','PDGFRB','ADIRF'),
                           c('CD3D','CD3E','CD4','CD8A'))
names(cell_type_markers) <- c('Acinar cell',
                              'Antibody-secreting cell',
                              'B cell',
                              'Ductal cell type 1',
                              'Ductal cell type 2',
                              'Endocrine cell',
                              'Endothelial cell',
                              'Fibroblast',
                              'Macrophage',
                              'Mast cell',
                              'Neutrophil',
                              'Stellate cell',
                              'T cell')
cell_type_colors <- c('steelblue',
                      'gold',
                      'gold4',
                      'lightgreen',
                      'coral',
                      'aquamarine1',
                      'cyan2',
                      'plum',
                      'yellowgreen',
                      'magenta2',
                      'burlywood',
                      'darkseagreen',
                      'pink',
                      'grey')
names(cell_type_colors) <- c(names(cell_type_markers),'Unknown')
setwd('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\final_results_v4\\cell_markers')
for(n in 1:length(cell_type_markers))
{
  genes_to_plot <- cell_type_markers[[n]]
  if(length(genes_to_plot)<=3)
  {
    png(paste0(names(cell_type_markers)[n],'.png'),width = 2000*length(genes_to_plot),height = 2000,res = 300)
    print(FeaturePlot(PDAC_all_cells_RPCA,features = genes_to_plot,ncol = length(genes_to_plot))&
          NoLegend()&NoAxes()&theme(plot.title = element_text(size=40)))
    dev.off()
  }
  else if(length(genes_to_plot)==4)
  {
    png(paste0(names(cell_type_markers)[n],'.png'),width = 4000,height = 4000,res = 300)
    print(FeaturePlot(PDAC_all_cells_RPCA,features = genes_to_plot,ncol = 2)&
          NoLegend()&NoAxes()&theme(plot.title = element_text(size=40)))
    dev.off()
  }
  else
  {
    png(paste0(names(cell_type_markers)[n],'.png'),width = 8000,height = 2000*ceiling(length(genes_to_plot)/4),res = 300)
    print(FeaturePlot(PDAC_all_cells_RPCA,features = genes_to_plot,ncol = 4)&
            NoLegend()&NoAxes()&theme(plot.title = element_text(size=40)))
    dev.off()
  }
}
# Reassign labels
PDAC_all_cells_RPCA$cell_type_v2 <- NaN
PDAC_all_cells_RPCA$cell_type_v2[PDAC_all_cells_RPCA$seurat_clusters_v2 %in% c(15)] <- 'Acinar cell'
PDAC_all_cells_RPCA$cell_type_v2[PDAC_all_cells_RPCA$seurat_clusters_v2 %in% c(19)] <- 'Antibody-secreting cell'
PDAC_all_cells_RPCA$cell_type_v2[PDAC_all_cells_RPCA$seurat_clusters_v2 %in% c(12,28)] <- 'B cell'
PDAC_all_cells_RPCA$cell_type_v2[PDAC_all_cells_RPCA$seurat_clusters_v2 %in% c(4,14)] <- 'Ductal cell type 1'
PDAC_all_cells_RPCA$cell_type_v2[PDAC_all_cells_RPCA$seurat_clusters_v2 %in% c(3,5,6,16,17,22,29,32)] <- 'Ductal cell type 2'
PDAC_all_cells_RPCA$cell_type_v2[PDAC_all_cells_RPCA$seurat_clusters_v2 %in% c(23)] <- 'Endocrine cell'
PDAC_all_cells_RPCA$cell_type_v2[PDAC_all_cells_RPCA$seurat_clusters_v2 %in% c(2)] <- 'Endothelial cell'
PDAC_all_cells_RPCA$cell_type_v2[PDAC_all_cells_RPCA$seurat_clusters_v2 %in% c(10,13)] <- 'Fibroblast'
PDAC_all_cells_RPCA$cell_type_v2[PDAC_all_cells_RPCA$seurat_clusters_v2 %in% c(7,9,18,20,24,30)] <- 'Macrophage'
PDAC_all_cells_RPCA$cell_type_v2[PDAC_all_cells_RPCA$seurat_clusters_v2 %in% c(21)] <- 'Mast cell'
PDAC_all_cells_RPCA$cell_type_v2[PDAC_all_cells_RPCA$seurat_clusters_v2 %in% c(11)] <- 'Neutrophil'
PDAC_all_cells_RPCA$cell_type_v2[PDAC_all_cells_RPCA$seurat_clusters_v2 %in% c(8)] <- 'Stellate cell'
PDAC_all_cells_RPCA$cell_type_v2[PDAC_all_cells_RPCA$seurat_clusters_v2 %in% c(0,1,25,26,34)] <- 'T cell'
# Cluster 27,31,33 are positive in markers for two types of cells
# They might be doublet or cells in hydrid state, we label them as unkonwn
PDAC_all_cells_RPCA$cell_type_v2[PDAC_all_cells_RPCA$seurat_clusters_v2 %in% c(27,31,33)] <- 'Unknown'
PDAC_all_cells_RPCA$cell_type_v2 <- factor(PDAC_all_cells_RPCA$cell_type_v2,levels = c(names(cell_type_markers),'Unknown'))
# A small cluster near fibroblast, were somehow clustered as Ductal cell type 2, 
# but based on results from markers, it is more likely fibroblast, so we would manually correct it to fibroblast
to_choose_from <- DimPlot(PDAC_all_cells_RPCA,group.by = 'cell_type_v2')+NoLegend()+xlim(-12,-6)+ylim(2.5,12)
cells_to_relabel <- CellSelector(to_choose_from)
cells_to_relabel <- intersect(colnames(PDAC_all_cells_RPCA)[PDAC_all_cells_RPCA$cell_type_v2 == 'Ductal cell type 2'],cells_to_relabel)
PDAC_all_cells_RPCA$cell_type_v2[cells_to_relabel] <- 'Fibroblast'
# A small cluster near steallate, were somehow clustered as Ductal cell type 2, 
# but based on results from markers, it is more likely stellate cell, so we would manually correct it to fibroblast
to_choose_from <- DimPlot(PDAC_all_cells_RPCA,group.by = 'cell_type_v2')+NoLegend()+xlim(-10,-5)+ylim(-5,5)
cells_to_relabel <- CellSelector(to_choose_from)
cells_to_relabel <- intersect(colnames(PDAC_all_cells_RPCA)[PDAC_all_cells_RPCA$cell_type_v2 == 'Ductal cell type 2'],cells_to_relabel)
PDAC_all_cells_RPCA$cell_type_v2[cells_to_relabel] <- 'Stellate cell'
# A small cluster near steallate, were somehow clustered as Ductal cell type 2, 
# but based on results from markers, it is not clear, so we would label it as unknown
to_choose_from <- DimPlot(PDAC_all_cells_RPCA,group.by = 'cell_type_v2')+NoLegend()+xlim(5,10)+ylim(-10,0)
cells_to_relabel <- CellSelector(to_choose_from)
cells_to_relabel <- intersect(colnames(PDAC_all_cells_RPCA)[PDAC_all_cells_RPCA$cell_type_v2 == 'Ductal cell type 2'],cells_to_relabel)
PDAC_all_cells_RPCA$cell_type_v2[cells_to_relabel] <- 'Unknown'
setwd('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\final_results_v4')
png('PDAC_all_cells_v2_RPCA_dim30_nn30_umap.png',width = 2500,height = 2000,res = 300)
DimPlot(PDAC_all_cells_RPCA,group.by = 'cell_type_v2')+
  NoAxes()+scale_color_manual(values = cell_type_colors)
dev.off()
saveRDS(PDAC_all_cells_RPCA,'C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\reference_dataset_core\\PDAC_all_cells_15430_RPCA_by_project.rds')
# Aggregate based on cell type and study the distribution of candidate genes from PDAC patient microarray data
PDAC_all_cells_RPCA <- readRDS('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\reference_dataset_core\\PDAC_all_cells_15430_RPCA_by_project.rds')
PDAC_all_cells_RPCA_exclude_unknown <- subset(PDAC_all_cells_RPCA, cell_type_v2!='Unknown')
png('PDAC_all_cells_v2_exclude_unknownRPCA_dim30_nn30_umap.png',width = 2500,height = 2000,res = 300)
DimPlot(PDAC_all_cells_RPCA_exclude_unknown,group.by = 'cell_type_v2',label = TRUE,label.size = 5)+
  NoLegend()+NoAxes()+ggtitle(NULL)+scale_color_manual(values = cell_type_colors)
dev.off()
# For comparison, also plot the original cell type label
PDAC_all_cells_RPCA_exclude_unknown$Cell_type[PDAC_all_cells_RPCA_exclude_unknown$Cell_type == 'Fibroblast cell'] <- 'Fibroblast'
PDAC_all_cells_RPCA_exclude_unknown$Cell_type[PDAC_all_cells_RPCA_exclude_unknown$Cell_type == 'Macrophage cell'] <- 'Macrophage'
cell_types_v1 <- c('Acinar cell','B cell','Ductal cell type 1','Ductal cell type 2',
                   'Endocrine cell','Endothelial cell','Fibroblast','Macrophage','Stellate cell','T cell')
PDAC_all_cells_RPCA_exclude_unknown$Cell_type <- factor(PDAC_all_cells_RPCA_exclude_unknown$Cell_type,
                                                        levels = cell_types_v1)
png('PDAC_all_cells_v1_exclude_unknownRPCA_dim30_nn30_umap.png',width = 2500,height = 2000,res = 300)
DimPlot(PDAC_all_cells_RPCA_exclude_unknown,group.by = 'Cell_type',label = TRUE,label.size = 5)+
  NoLegend()+NoAxes()+ggtitle(NULL)+scale_color_manual(values = cell_type_colors[cell_types_v1])
dev.off()

PDAC_all_cells_RPCA_tumor_exclude_unknown <- subset(PDAC_all_cells_RPCA_exclude_unknown,Type=='Tumor')
cell_type_v2_pseudobulk <- AggregateExpression(PDAC_all_cells_RPCA_tumor_exclude_unknown,
                                               assays = 'RNA',
                                                group.by = 'cell_type_v2')
cell_type_v2_pseudobulk <- cell_type_v2_pseudobulk$RNA
STS_LTS_genes <- read.csv('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_clinical\\results\\limma_results_filtered_protein_coding.csv',row.names = 1)
#STS_LTS_genes <- STS_LTS_genes[!(rownames(STS_LTS_genes) %in% c('SLCO1B3','C4orf51','GPR182','SCG3','SHOX')),]
STS_LTS_genes_cell_type_v2_distribution <- data.frame(matrix(data = NA,ncol = length(rownames(STS_LTS_genes)),nrow = length(colnames(cell_type_v2_pseudobulk))),
                                                      row.names = colnames(cell_type_v2_pseudobulk))
colnames(STS_LTS_genes_cell_type_v2_distribution) <- rownames(STS_LTS_genes)

STS_LTS_genes_total_counts <- rep(NaN,length(rownames(STS_LTS_genes)))
names(STS_LTS_genes_total_counts) <- rownames(STS_LTS_genes)
for(n_gene in 1:length(colnames(STS_LTS_genes_cell_type_v2_distribution)))
{
  gene_to_calculate <- colnames(STS_LTS_genes_cell_type_v2_distribution)[n_gene]
  if(gene_to_calculate %in% rownames(cell_type_v2_pseudobulk))
  {
    gene_to_calcuate_pseudobulk <- cell_type_v2_pseudobulk[gene_to_calculate,]
    STS_LTS_genes_cell_type_v2_distribution[,gene_to_calculate] <- gene_to_calcuate_pseudobulk/sum(gene_to_calcuate_pseudobulk)
    STS_LTS_genes_total_counts[gene_to_calculate] <- sum(gene_to_calcuate_pseudobulk)
  }
  else
  {
    print(sprintf('%s is not in the scRNAseq dataset.',gene_to_calculate))  
  }
}
STS_LTS_genes_total_counts_filtered <- na.omit(STS_LTS_genes_total_counts)
STS_LTS_genes_total_counts_filtered <- STS_LTS_genes_total_counts_filtered[STS_LTS_genes_total_counts_filtered>=100]
STS_LTS_genes_filtered <- STS_LTS_genes[names(STS_LTS_genes_total_counts_filtered),]
STS_LTS_genes_cell_type_v2_distribution_filtered <- STS_LTS_genes_cell_type_v2_distribution[names(STS_LTS_genes_total_counts_filtered)]
# Visualize STS gene distribution across different cell types
STS_genes <- rownames(STS_LTS_genes_filtered)[STS_LTS_genes_filtered$logFC>0]
STS_genes_plots <- vector(mode = 'list',length = length(STS_genes))
for(n_gene in 1:length(STS_genes))
{
  gene_to_plot <- STS_genes[n_gene]
  data_to_plot <- data.frame('percentage' = STS_LTS_genes_cell_type_v2_distribution_filtered[,gene_to_plot],
                             'cell_type_v2' = rownames(STS_LTS_genes_cell_type_v2_distribution_filtered))
  STS_genes_plots[[n_gene]] <- ggplot(data_to_plot, aes(x = '',y = percentage,fill = cell_type_v2))+
    geom_bar(stat="identity", width=1)+coord_polar('y', start=0)+xlab(NULL)+ylab(NULL)+theme_minimal()+NoLegend()+
    ggtitle(gene_to_plot)+theme(plot.title = element_text(hjust = 0.5,size = 20,face = 'bold'))+
    scale_fill_manual(values = cell_type_colors)+NoLegend()
  
}
png('Pseudobulk_distribution_of_STS_genes.png',width = 1000*5,height = 1000*3,res = 300)
print(wrap_plots(STS_genes_plots,nrow = 3,ncol = 5))
dev.off()
png('Feature_plot_of_STS_genes.png',width = 1000*5,height = 1000*3,res = 300)
FeaturePlot(PDAC_all_cells_RPCA_exclude_unknown,features = STS_genes,ncol = 5)&
  NoAxes()&NoLegend()&theme(plot.title = element_text(size=20))
dev.off()
# Plot LTS genes
LTS_genes <- rownames(STS_LTS_genes_filtered)[STS_LTS_genes_filtered$logFC<0]
LTS_genes_plots <- vector(mode = 'list',length = length(LTS_genes))
for(n_gene in 1:length(LTS_genes))
{
  gene_to_plot <- LTS_genes[n_gene]
  data_to_plot <- data.frame('percentage' = STS_LTS_genes_cell_type_v2_distribution_filtered[,gene_to_plot],
                             'cell_type_v2' = rownames(STS_LTS_genes_cell_type_v2_distribution_filtered))
  LTS_genes_plots[[n_gene]] <- ggplot(data_to_plot, aes(x = '',y = percentage,fill = cell_type_v2))+
    geom_bar(stat="identity", width=1)+coord_polar('y', start=0)+xlab(NULL)+ylab(NULL)+theme_minimal()+NoLegend()+
    ggtitle(gene_to_plot)+theme(plot.title = element_text(hjust = 0.5,size = 20,face = 'bold'))+
    scale_fill_manual(values = cell_type_colors)+NoLegend()
}
png('Pseudobulk_distribution_of_LTS_genes.png',width = 1000*5,height = 1000,res = 300)
print(wrap_plots(LTS_genes_plots,nrow = 1,ncol = 5))
dev.off()
png('Feature_plot_of_LTS_genes.png',width = 1000*5,height = 1000,res = 300)
FeaturePlot(PDAC_all_cells_RPCA_exclude_unknown,features = LTS_genes,ncol = 5,order = TRUE,max.cutoff = 1)&
  NoAxes()&NoLegend()&theme(plot.title = element_text(size=20))
dev.off()
######################################Integration of acinar and ductal cells#####################################
setwd('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\reference_dataset_core')
PDAC_AD_cells <- readRDS('PDAC_all_cells_15430_RPCA_by_project.rds')
PDAC_AD_cells <- subset(PDAC_AD_cells, cell_type_v2 %in% c('Acinar cell','Ductal cell type 1','Ductal cell type 2'))
sample_list <- unique(PDAC_AD_cells$Patient2) # 72 samples in total
normal_sample_list <- unique(PDAC_AD_cells$Patient2[PDAC_AD_cells$Type == 'Normal'])# 14 normal samples in total
tumor_sample_list <- unique(PDAC_AD_cells$Patient2[PDAC_AD_cells$Type == 'Tumor'])# 58 tumor samples in total
hvg_info_normal <- vector('list', length = length(normal_sample_list))
hvg_info_tumor <- vector('list',length = length(tumor_sample_list))
# Collect HVGs for each sample
for(i_sample in 1:length(sample_list))
{
  sample_to_analyze <- subset(PDAC_AD_cells,Patient2 == sample_list[i_sample])
  sample_to_analyze <- NormalizeData(sample_to_analyze)
  sample_to_analyze <- FindVariableFeatures(sample_to_analyze)
  
  if(sample_list[i_sample] %in% normal_sample_list)
  {
    sample_index <- which(sample_list[i_sample] == normal_sample_list)
    hvg_info_normal[[sample_index]] <- VariableFeatures(sample_to_analyze)
  }
  else
  {
    sample_index <- which(sample_list[i_sample] == tumor_sample_list)
    hvg_info_tumor[[sample_index]] <- VariableFeatures(sample_to_analyze)
  }
}
HVG_occurance <- data.frame(row.names = rownames(PDAC_AD_cells),
                            total_occurance = rep(0,length(rownames(PDAC_AD_cells))),
                            tumor_occurance = rep(0,length(rownames(PDAC_AD_cells))),
                            normal_occurance = rep(0,length(rownames(PDAC_AD_cells))))
# Calculate occurance as HVG for each gene
for(gene in rownames(PDAC_AD_cells))
{
  tumor_occurance <- 0
  for(n in 1:length(hvg_info_tumor))
  {
    tumor_occurance <- tumor_occurance+(gene %in% hvg_info_tumor[[n]])
  }
  normal_occurance <- 0 
  for(n in 1:length(hvg_info_normal))
  {
    normal_occurance <- normal_occurance+(gene %in% hvg_info_normal[[n]])
  }
  HVG_occurance[gene,'total_occurance'] <- tumor_occurance+normal_occurance
  HVG_occurance[gene,'tumor_occurance'] <- tumor_occurance
  HVG_occurance[gene,'normal_occurance'] <- normal_occurance
}
HVG_occurance <- HVG_occurance[order(HVG_occurance$total_occurance,decreasing = TRUE),]
saveRDS(HVG_occurance,'C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\final_results_v4\\HVG_occurance_acinar_ductal.rds')
# Second, HVGs for Harmony integration were identified by picking genes pop up as HVGs in at least N patients
setwd('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\reference_dataset_core')
PDAC_AD_cells <- readRDS('PDAC_all_cells_15430_RPCA_by_project.rds')
PDAC_AD_cells <- DietSeurat(PDAC_AD_cells)
PDAC_AD_cells <- subset(PDAC_AD_cells, cell_type_v2 %in% c('Acinar cell','Ductal cell type 1','Ductal cell type 2'))
PDAC_AD_cells <- NormalizeData(PDAC_AD_cells)
HVG_occurance <- readRDS('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\final_results_v4\\HVG_occurance_acinar_ductal.rds')
VariableFeatures(PDAC_AD_cells) <- rownames(HVG_occurance)[HVG_occurance$total_occurance>=35]#35
PDAC_AD_cells <- ScaleData(PDAC_AD_cells)
PDAC_AD_cells <- RunPCA(PDAC_AD_cells)
# Harmony has some randomness in the algorithm which can affect the output.
# This effect can somhow significantly affect trajectory influence by monocle3.
# Try a few random seeds to pick the best one.
setwd('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\final_results_v4\\harmony random seed optimization HVG35 ndim 25')
for(n_seed in 1:10)
{
  set.seed(n_seed)
  PDAC_AD_cells <- harmony::RunHarmony(PDAC_AD_cells,c('Patient2'),dims.use = 1:25)
  PDAC_AD_cells <- RunUMAP(PDAC_AD_cells,dims = 1:6,n.neighbors = 5,reduction = 'harmony',reduction.name = 'umap_harmony',seed.use = 42)
  PDAC_AD_cells_cds <- SeuratWrappers::as.cell_data_set(PDAC_AD_cells)
  reducedDimNames(PDAC_AD_cells_cds)[3] <- 'UMAP' # We need to use the UMAP_HARMONY for further analysis, so we need to change its name to UMAP for monocle3 to recognize it as umap
  PDAC_AD_cells_cds <- cluster_cells(PDAC_AD_cells_cds,reduction_method = 'UMAP', k = 75)
  PDAC_AD_cells_cds <- learn_graph(PDAC_AD_cells_cds,close_loop = FALSE,
                                   learn_graph_control = list(minimal_branch_len = 5))
  png(sprintf('seed %d SPINK1.png',n_seed),width = 600,height = 500,units = 'px',res = 150)
  print(FeaturePlot(PDAC_AD_cells,'SPINK1')+NoAxes())
  dev.off()
  png(sprintf('seed %d trajectory.png',n_seed),width = 600,height = 500,units = 'px',res = 150)
  print(plot_cells(PDAC_AD_cells_cds))
  dev.off()
}
set.seed(5) # Seed 5 ndim 25 and uamp ndim 6 gives best results
PDAC_AD_cells <- harmony::RunHarmony(PDAC_AD_cells,c('Patient2'),dims.use = 1:25) #25
PDAC_AD_cells <- RunUMAP(PDAC_AD_cells,dims = 1:6,n.neighbors = 5,reduction = 'harmony',reduction.name = 'umap_harmony',seed.use = 42)
png('SPINK1 overview harmony seed 5 ndim 25 umap ndim 6.png',width = 1250,height = 1000,units = 'px',res = 300)
FeaturePlot(PDAC_AD_cells,'SPINK1')+NoAxes()+ggtitle(NULL)
dev.off()
png('Cell type markers.png',width = 3000, height = 1000,res = 300)
FeaturePlot(PDAC_AD_cells, features = c('PRSS1','SOX9','KRT19'),ncol = 3)&NoAxes()
dev.off()

# Relabel cells based on the new integration results
PDAC_AD_cells <- FindNeighbors(PDAC_AD_cells,reduction = 'harmony',dims = 1:25)
PDAC_AD_cells <- FindClusters(PDAC_AD_cells,resolution = 0.3)
PDAC_AD_cells$seurat_clusters_AD_harmony <- Idents(PDAC_AD_cells)
png('Seurat clusters after harmony by patients on acinar ductal cells.png',width = 1000,height = 1000,res = 300)
DimPlot(PDAC_AD_cells,group.by = 'seurat_clusters_AD_harmony',label = TRUE)+NoAxes()+NoLegend()
dev.off()
PDAC_AD_cells$cell_type_v2_AD_harmony <- NaN
PDAC_AD_cells$cell_type_v2_AD_harmony[PDAC_AD_cells$seurat_clusters_AD_harmony %in% c(4,7)] <- 'Acinar cell'
PDAC_AD_cells$cell_type_v2_AD_harmony[PDAC_AD_cells$seurat_clusters_AD_harmony %in% c(1)] <- 'Ductal cell type 1'
PDAC_AD_cells$cell_type_v2_AD_harmony[PDAC_AD_cells$seurat_clusters_AD_harmony %in% c(0,2,3,5,6,8,9)] <- 'Ductal cell type 2'
PDAC_AD_cells$cell_type_v2_AD_harmony <- factor(PDAC_AD_cells$cell_type_v2_AD_harmony,
                                                levels = c('Acinar cell','Ductal cell type 1','Ductal cell type 2'))
png('acinar ductal cells after harmony by patients.png',width = 1500,height = 1000,res = 300)
DimPlot(PDAC_AD_cells,group.by = 'cell_type_v2_AD_harmony')+
  scale_color_manual(values = c('steelblue','lightgreen','coral'))+
  ggtitle(NULL)+NoAxes()
dev.off()
png('Tumor vs normal for acinar ductal cells.png',width = 2000,height = 1000,res = 300)
DimPlot(PDAC_AD_cells,split.by = 'Type',group.by = 'cell_type_v2_AD_harmony')+
  scale_color_manual(values = c('steelblue','lightgreen','coral'))+
  ggtitle(NULL)+NoAxes()+NoLegend()
dev.off()

# Trajectory inference by monocle3
PDAC_AD_cells_cds <- SeuratWrappers::as.cell_data_set(PDAC_AD_cells)
reducedDimNames(PDAC_AD_cells_cds)[3] <- 'UMAP' # We need to use the UMAP_HARMONY for further analysis, so we need to change its name to UMAP for monocle3 to recognize it as umap
PDAC_AD_cells_cds <- cluster_cells(PDAC_AD_cells_cds,reduction_method = 'UMAP', k = 75)
PDAC_AD_cells_cds <- learn_graph(PDAC_AD_cells_cds,close_loop = FALSE,
                                 learn_graph_control = list(minimal_branch_len = 5))
PDAC_AD_cells_cds <- order_cells(PDAC_AD_cells_cds)
PDAC_AD_cells$pseudotime <- pseudotime(PDAC_AD_cells_cds, reduction_method = 'UMAP')
setwd('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\final_results_v4')
png('Pesudotime trajectory of PDAC k 75 minimal branch 5.png',width = 2500,height = 2000, units = 'px', res = 300)
plot_cells(PDAC_AD_cells_cds, color_cells_by = 'pseudotime', label_cell_groups=FALSE,
           label_leaves=FALSE, label_branch_points=FALSE,
           graph_label_size=3, trajectory_graph_color = 'green', 
           trajectory_graph_segment_size = 1.5)+NoAxes()+
           theme(legend.title= element_blank(),legend.text = element_text(size = 20),legend.key.size = unit(1,'cm'))
dev.off()
saveRDS(PDAC_AD_cells,'C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\reference_dataset_core\\PDAC_acinar_ductal_cells_15430_genes_Harmony_by_patients_v4.rds')
saveRDS(PDAC_AD_cells_cds,'C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\reference_dataset_core\\PDAC_acinar_ductal_cells_15430_genes_Harmony_by_patients_v4_cds.rds')

# Select high spink1 and low spink1 trajectory
setwd('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\final_results_v4')
high_SPINK1_trajectory <- colnames(choose_graph_segments(PDAC_AD_cells_cds))
saveRDS(high_SPINK1_trajectory,'high_SPINK1_trajectory.rds')
low_SPINK1_trajectory <- colnames(choose_graph_segments(PDAC_AD_cells_cds))
saveRDS(low_SPINK1_trajectory,'low_SPINK1_trajectory.rds')
# Plot trajectories for SPINK1 along these two trajectories
PDAC_AD_cells <- readRDS('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\reference_dataset_core\\PDAC_acinar_ductal_cells_15430_genes_Harmony_by_patients_v4.rds')
PDAC_AD_cells_cds <- readRDS('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\reference_dataset_core\\PDAC_acinar_ductal_cells_15430_genes_Harmony_by_patients_v4_cds.rds')
high_SPINK1_trajectory <- readRDS('high_SPINK1_trajectory.rds')
low_SPINK1_trajectory <- readRDS('low_SPINK1_trajectory.rds')
p1 <- DimPlot(PDAC_AD_cells,cells.highlight = high_SPINK1_trajectory,reduction = 'umap_harmony',cols.highlight = 'gold3')+NoLegend()+NoAxes()
p2 <- DimPlot(PDAC_AD_cells,cells.highlight = low_SPINK1_trajectory, reduction = 'umap_harmony',cols.highlight = 'deepskyblue')+NoLegend()+NoAxes()
png('Cells on high and low SPINK1 trajectories.png',width = 1000,height = 2000,res = 300)
print(p1/p2)
dev.off()
plot_gene_trajectory(PDAC_AD_cells,c('SPINK1'),method_to_use = 'median',pseudotime_cap = c(70,75),x_lim = c(0,80),
                     output_width = 2000,output_height = 1200, bin_size = 10,violin_width = 5)

# Partition cells into start, high SPINK1 end and low SPINK1 end
PDAC_DT2_cells <- readRDS('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\reference_dataset_core\\PDAC_acinar_ductal_cells_15430_genes_Harmony_by_patients_v4.rds')
PDAC_DT2_cells <- subset(PDAC_DT2_cells,cell_type_v2_AD_harmony == 'Ductal cell type 2')
PDAC_DT2_cells$trajectory_position <- 'Other'
transition_cells <- intersect(colnames(PDAC_DT2_cells)[PDAC_DT2_cells$pseudotime>=28 & PDAC_DT2_cells$pseudotime<=42],
                         low_SPINK1_trajectory)
PDAC_DT2_cells$trajectory_position[transition_cells] <- 'Transition'
bifurcation_cells <- intersect(colnames(PDAC_DT2_cells)[PDAC_DT2_cells$pseudotime>=60 & PDAC_DT2_cells$pseudotime<=66.5],
                               low_SPINK1_trajectory)
PDAC_DT2_cells$trajectory_position[bifurcation_cells] <- 'Bifurcation'
high_SPINK1_end_cells <- colnames(PDAC_DT2_cells)[PDAC_DT2_cells$seurat_clusters_AD_harmony == 5]
PDAC_DT2_cells$trajectory_position[high_SPINK1_end_cells] <- 'High SPINK1 end'
low_SPINK1_end_cells <- intersect(colnames(PDAC_DT2_cells)[PDAC_DT2_cells$pseudotime >=72],low_SPINK1_trajectory)
PDAC_DT2_cells$trajectory_position[low_SPINK1_end_cells] <- 'Low SPINK1 end'
PDAC_DT2_cells$trajectory_position <- factor(PDAC_DT2_cells$trajectory_position,
                                             levels = c('Transition','Other','Bifurcation','High SPINK1 end','Low SPINK1 end'))
png('transition_28_42 bifurcation_60_66.5 low_72 high_cluster_5 for DEG testing with legend.png',width = 1500,height = 1000,res = 300)
DimPlot(PDAC_DT2_cells,group.by = 'trajectory_position')+
  ggtitle(NULL)+NoAxes()+scale_color_manual(values = c('chartreuse2','grey','chartreuse4','gold3','deepskyblue'))
dev.off()
png('transition_28_42 bifurcation_60_66.5 low_72 high_cluster_5 for DEG testing.png',width = 750,height = 1000,res = 300)
DimPlot(PDAC_DT2_cells,group.by = 'trajectory_position')+
  ggtitle(NULL)+NoAxes()+NoLegend()+scale_color_manual(values = c('chartreuse2','grey','chartreuse4','gold3','deepskyblue'))
dev.off()
saveRDS(PDAC_DT2_cells,'C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\reference_dataset_core\\PDAC_ductal_type_2_cells_15430_genes_Harmony_by_patients_v4.rds')

# Identify genes differentially expressed at the end point of the high SPINK1 and low SPINK1 trajectory
PDAC_DT2_cells <- readRDS('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\reference_dataset_core\\PDAC_ductal_type_2_cells_15430_genes_Harmony_by_patients_v4.rds')
Idents(PDAC_DT2_cells) <- 'trajectory_position'
bifurcation_vs_transition_markers <- FindMarkers(PDAC_DT2_cells,ident.1 = 'Bifurcation',ident.2 = 'Transition')
bifurcation_vs_transition_markers <- bifurcation_vs_transition_markers[order(bifurcation_vs_transition_markers$avg_log2FC,decreasing = TRUE),]
saveRDS(bifurcation_vs_transition_markers,'bifurcation_vs_transition_markers.rds')
high_vs_transition_markers <- FindMarkers(PDAC_DT2_cells,ident.1 = 'High SPINK1 end',ident.2 = 'Transition')
high_vs_transition_markers <- high_vs_transition_markers[order(high_vs_transition_markers$avg_log2FC,decreasing = TRUE),]
saveRDS(high_vs_transition_markers,'high_vs_transition_markers.rds')
low_vs_transition_markers <- FindMarkers(PDAC_DT2_cells,ident.1 = 'Low SPINK1 end',ident.2 = 'Transition')
low_vs_transition_markers <- low_vs_transition_markers[order(low_vs_transition_markers$avg_log2FC,decreasing = TRUE),]
saveRDS(low_vs_transition_markers,'low_vs_transition_markers.rds')
high_vs_low_markers <- FindMarkers(PDAC_DT2_cells,ident.1 = 'High SPINK1 end',ident.2 = 'Low SPINK1 end',logfc.threshold = 0)
high_vs_low_markers <- high_vs_low_markers[order(high_vs_low_markers$avg_log2FC,decreasing = TRUE),]
saveRDS(high_vs_low_markers,'high_vs_low_markers.rds')

# High SPINK1 end cells have significantly higher proliferative markers
setwd('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\final_results_v4')
high_vs_transition_markers <- readRDS('high_vs_transition_markers.rds')
low_vs_transition_markers <- readRDS('low_vs_transition_markers.rds')
high_vs_low_markers <- readRDS('high_vs_low_markers.rds')
# SPINK1
print(low_vs_transition_markers[c('SPINK1'),])
print(high_vs_transition_markers[c('SPINK1'),])
# Proliferation related genes
print(high_vs_transition_markers[c('MKI67','TOP2A','CCNB1','CCNB2'),])
print(high_vs_low_markers[c('MKI67','TOP2A','CCNB1','CCNB2'),])
# DNA damage related genes
print(low_vs_transition_markers[c('TPT1','NABP1','GADD45B','DDIT3'),])
print(high_vs_low_markers[c('TPT1','NABP1','GADD45B','DDIT3'),])
# EMT related genes
print(low_vs_transition_markers[c('ELF3','KLF5','CDH1','VIM','ZEB1','CDH2'),])
print(high_vs_transition_markers[c('ELF3','KLF5','CDH1','VIM','ZEB1','CDH2'),])

# Visualize proliferative markers upregulated in high SPINK1 end
png('Proliferative markers.png',width = 2000,height = 2000,res = 300)
FeaturePlot(PDAC_DT2_cells,features = c('MKI67','TOP2A','CCNB1','CCNB2'),ncol = 2)&NoAxes()
dev.off()
# Visualize DNA damage related genes upregulated in low SPINK1 end
p1 <- FeaturePlot(PDAC_DT2_cells,features = c('TPT1'),min.cutoff = 3,max.cutoff = 6)&NoAxes()
p2 <- FeaturePlot(PDAC_DT2_cells,features = c('NABP1'),order = TRUE)&NoAxes()
p3 <- FeaturePlot(PDAC_DT2_cells,features = c('GADD45B'),order = TRUE)&NoAxes()
p4 <- FeaturePlot(PDAC_DT2_cells,features = c('DDIT3'),order = TRUE)&NoAxes()
png('DNA damage related genes.png',width = 2000, height = 2000, res = 300)
print((p1|p2)/(p3|p4))
dev.off()
# Visualize EMT markers
png('EMT markers 1.png',width = 3000,height = 1000,res = 300)
FeaturePlot(PDAC_DT2_cells,features = c('ELF3','KLF5','CDH1'),ncol = 3)&NoAxes()
dev.off()
png('EMT markers 2.png',width = 3000,height = 1000,res = 300)
p1 <- FeaturePlot(PDAC_DT2_cells,features = c('VIM'),order = TRUE)&NoAxes()
p2 <- FeaturePlot(PDAC_DT2_cells,features = c('ZEB1'),max.cutoff = 2,order = TRUE)&NoAxes()
p3 <- FeaturePlot(PDAC_DT2_cells,features = c('CDH2'),max.cutoff = 2,order = TRUE)&NoAxes()
print(p1|p2|p3)
dev.off()

# Visualize patient composition of high spink1 end and low spink1 end cells
high_SPINK1_end_cells_patient <- summary(as.factor(PDAC_DT2_cells$Patient2[high_SPINK1_end_cells]))
high_SPINK1_end_cells_patient['GSE155698_T11'] <- high_SPINK1_end_cells_patient['GSE155698_T11A']+high_SPINK1_end_cells_patient['GSE155698_T11B']
high_SPINK1_end_cells_patient <- high_SPINK1_end_cells_patient[!(names(high_SPINK1_end_cells_patient) %in% c('GSE155698_T11A','GSE155698_T11B'))]
high_SPINK1_end_cells_patient <- data.frame(patient = names(high_SPINK1_end_cells_patient),
                                            n_cells = high_SPINK1_end_cells_patient,
                                            p_cells = high_SPINK1_end_cells_patient/length(high_SPINK1_end_cells))
high_SPINK1_end_cells_patient <- high_SPINK1_end_cells_patient[order(high_SPINK1_end_cells_patient$p_cells,decreasing = FALSE),]
high_SPINK1_end_cells_patient$patient <- factor(high_SPINK1_end_cells_patient$patient,levels = high_SPINK1_end_cells_patient$patient)

low_SPINK1_end_cells_patient <- summary(as.factor(PDAC_DT2_cells$Patient2[low_SPINK1_end_cells]))
low_SPINK1_end_cells_patient['GSE155698_T11'] <- low_SPINK1_end_cells_patient['GSE155698_T11A']+low_SPINK1_end_cells_patient['GSE155698_T11B']
low_SPINK1_end_cells_patient <- low_SPINK1_end_cells_patient[!(names(low_SPINK1_end_cells_patient) %in% c('GSE155698_T11A','GSE155698_T11B'))]
low_SPINK1_end_cells_patient <- data.frame(patient = names(low_SPINK1_end_cells_patient),
                                           n_cells = low_SPINK1_end_cells_patient,
                                           p_cells = low_SPINK1_end_cells_patient/length(low_SPINK1_end_cells))
low_SPINK1_end_cells_patient <- low_SPINK1_end_cells_patient[order(low_SPINK1_end_cells_patient$p_cells,decreasing = FALSE),]
low_SPINK1_end_cells_patient$patient <- factor(low_SPINK1_end_cells_patient$patient,levels = low_SPINK1_end_cells_patient$patient)

p1 <- ggplot(high_SPINK1_end_cells_patient[(nrow(high_SPINK1_end_cells_patient)-20):nrow(high_SPINK1_end_cells_patient),],
       aes(x = patient,y = p_cells))+geom_bar(stat = "identity")+coord_flip()+
       ggtitle(sprintf('High SPINK1 end cells (N= %d cells)',length(high_SPINK1_end_cells)))+xlab(NULL)+ylab('Percentage')+ylim(0,0.25)
p2 <- ggplot(low_SPINK1_end_cells_patient[(nrow(low_SPINK1_end_cells_patient)-20):nrow(low_SPINK1_end_cells_patient),],
       aes(x = patient,y = p_cells))+geom_bar(stat = "identity")+coord_flip()+
       ggtitle(sprintf('Low SPINK1 end cells (N = %d cells)',length(low_SPINK1_end_cells)))+xlab(NULL)+ylab('Percentage')+ylim(0,0.25)
png('patient distribution of high and low SPINK1 end 72 cells.png',width = 3000,height = 2000,res = 300)
print(p1|p2)
dev.off()

#Visualize percentage of high and low SPINK1 end cells in each patient's ductal type 2 cells
normal_sample_list <- unique(PDAC_AD_cells$Patient2[PDAC_AD_cells$Type == 'Normal'])# 14 normal samples in total
high_low_end_SPINK1_patient_summary <- data.frame(row.names = unique(PDAC_AD_cells$Patient2),
                                                  'N_ductal_type_2' = rep(NaN,length(unique(PDAC_AD_cells$Patient2))),
                                                  'N_high_SPINK1_end' = rep(NaN,length(unique(PDAC_AD_cells$Patient2))),
                                                  'N_low_SPINK1_end' = rep(NaN,length(unique(PDAC_AD_cells$Patient2))))
for(n in unique(PDAC_AD_cells$Patient2))
{
  ductal_type_2_of_patient <- colnames(PDAC_AD_cells)[(PDAC_AD_cells$Patient2 == n)&(PDAC_AD_cells$cell_type_v2_AD_harmony == 'Ductal cell type 2')]
  high_low_end_SPINK1_patient_summary[n,'N_ductal_type_2'] <- length(ductal_type_2_of_patient)
  high_low_end_SPINK1_patient_summary[n,'N_high_SPINK1_end'] <- length(intersect(ductal_type_2_of_patient,high_SPINK1_end_cells))
  high_low_end_SPINK1_patient_summary[n,'N_low_SPINK1_end'] <- length(intersect(ductal_type_2_of_patient,low_SPINK1_end_cells))
}
high_low_end_SPINK1_patient_summary['GSE155698_T11',] <- high_low_end_SPINK1_patient_summary['GSE155698_T11A',]+high_low_end_SPINK1_patient_summary['GSE155698_T11B',]
high_low_end_SPINK1_patient_summary <- high_low_end_SPINK1_patient_summary[!(rownames(high_low_end_SPINK1_patient_summary) %in% c('GSE155698_T11A','GSE155698_T11B')),]
high_low_end_SPINK1_patient_summary['p_high_SPINK1_end'] <- high_low_end_SPINK1_patient_summary['N_high_SPINK1_end']/high_low_end_SPINK1_patient_summary['N_ductal_type_2']
high_low_end_SPINK1_patient_summary['p_low_SPINK1_end'] <- high_low_end_SPINK1_patient_summary['N_low_SPINK1_end']/high_low_end_SPINK1_patient_summary['N_ductal_type_2']
high_low_end_SPINK1_patient_summary['p_high_SPINK1_end_normalized'] <- high_low_end_SPINK1_patient_summary['p_high_SPINK1_end']/(length(high_SPINK1_end_cells)/ncol(PDAC_DT2_cells))
high_low_end_SPINK1_patient_summary['p_low_SPINK1_end_normalized'] <- high_low_end_SPINK1_patient_summary['p_low_SPINK1_end']/(length(low_SPINK1_end_cells)/ncol(PDAC_DT2_cells))
high_low_end_SPINK1_patient_summary <- high_low_end_SPINK1_patient_summary[!(rownames(high_low_end_SPINK1_patient_summary) %in% normal_sample_list),]
# Assign each patient to a group based on their normalized percentage
high_low_end_SPINK1_patient_summary['trajectory_type'] <- 'Not determined'
for(n in rownames(high_low_end_SPINK1_patient_summary))
{
  p_high <- high_low_end_SPINK1_patient_summary[n,'p_high_SPINK1_end_normalized']
  p_low <- high_low_end_SPINK1_patient_summary[n,'p_low_SPINK1_end_normalized']
  if(high_low_end_SPINK1_patient_summary[n,'N_ductal_type_2']<200){next}
  if(p_high<0.5)
  {
    if(p_low<0.5){high_low_end_SPINK1_patient_summary[n,'trajectory_type'] <- 'Low proliferation low DNA damage'}
    else if(p_low>=0.5&p_low<=1.5){high_low_end_SPINK1_patient_summary[n,'trajectory_type'] <- 'Low proliferation normal DNA damage'}
    else if(p_low>1.5){high_low_end_SPINK1_patient_summary[n,'trajectory_type'] <- 'Low proliferation high DNA damage'}
  }
  else if(p_high>=0.5 & p_high<=1.5)
  {
    if(p_low<0.5){high_low_end_SPINK1_patient_summary[n,'trajectory_type'] <- 'Normal proliferation low DNA damage'}
    else if(p_low>=0.5&p_low<=1.5){high_low_end_SPINK1_patient_summary[n,'trajectory_type'] <- 'Normal proliferation normal DNA damage'}
    else if(p_low>1.5){high_low_end_SPINK1_patient_summary[n,'trajectory_type'] <- 'Normal proliferation high DNA damage'}
  }
  else if(p_high>1.5)
  {
    if(p_low<0.5){high_low_end_SPINK1_patient_summary[n,'trajectory_type'] <- 'High proliferation low DNA damage'}
    else if(p_low>=0.5&p_low<=1.5){high_low_end_SPINK1_patient_summary[n,'trajectory_type'] <- 'High proliferation normal DNA damage'}
    else if(p_low>1.5){high_low_end_SPINK1_patient_summary[n,'trajectory_type'] <- 'High proliferation high DNA damage'}
  }
}
saveRDS(high_low_end_SPINK1_patient_summary,'high_low_end_SPINK1_patient_summary.rds')
high_low_end_SPINK1_patient_summary_filtered <- high_low_end_SPINK1_patient_summary[high_low_end_SPINK1_patient_summary$N_ductal_type_2>=200,]

for_violin_plot <- data.frame('group' = c(rep('p_high_SPINK1_end',nrow(high_low_end_SPINK1_patient_summary_filtered)),rep('p_low_SPINK1_end',nrow(high_low_end_SPINK1_patient_summary_filtered))),
                              'percentage' = c(high_low_end_SPINK1_patient_summary_filtered$p_high_SPINK1_end_normalized,high_low_end_SPINK1_patient_summary_filtered$p_low_SPINK1_end_normalized))
setwd('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\final_results_v4')
png('Normalized percentgae of high and low end SPINK1 cells in each patient.png',width = 2000,height = 2000,res = 300)
set.seed(1)
ggplot(for_violin_plot,aes(x = group,y = percentage))+
  geom_violin()+
  geom_hline(yintercept = 1,linetype = 'dashed',color = 'red')+
  geom_jitter(width = 0.3,height = 0)+
  xlab(NULL)+ylab('Ratio of percentage in each patient to \nexpected percentage')+
  scale_x_discrete(label = c('High SPINK1 end','Low SPINK1 end'))+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20))
dev.off()
png('Scatter plot of normalized percentgae of high and low end SPINK1 cells in each patient.png',width = 2000,height = 2000,res = 300)
ggplot(high_low_end_SPINK1_patient_summary_filtered,aes(x = p_high_SPINK1_end_normalized,y = p_low_SPINK1_end_normalized))+
  geom_hline(yintercept = 0.5,linetype = 'dashed',color = 'red')+geom_hline(yintercept = 1.5,linetype = 'dashed',color = 'red')+
  geom_vline(xintercept = 0.5,linetype = 'dashed',color = 'red')+geom_vline(xintercept = 1.5,linetype = 'dashed',color = 'red')+
  geom_point()+xlab('High SPINK1 end')+ylab('Low SPINK1 end')+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20))
dev.off()

# Visualize the distribution of each patient's type 2 ductal cells along the trajectory
setwd('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\final_results_v4\\patient trajectory selection')
for(n_dir in unique(high_low_end_SPINK1_patient_summary_filtered$trajectory_type))
{
  dir.create(n_dir)
}
for(n in rownames(high_low_end_SPINK1_patient_summary_filtered))
{
  setwd(paste0('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\final_results_v4\\patient trajectory selection\\',
               high_low_end_SPINK1_patient_summary_filtered[n,'trajectory_type']))
  png(paste0(n,'.png'),width = 1000,height = 1000,res = 300)
  if(n=='GSE155698_T11')
  {
    print(DimPlot(PDAC_DT2_cells,cells.highlight = colnames(PDAC_DT2_cells)[PDAC_DT2_cells$Patient2 %in% c('GSE155698_T11A','GSE155698_T11B')])+
      ggtitle(sprintf('%s (N=%d)',n,high_low_end_SPINK1_patient_summary_filtered[n,'N_ductal_type_2']))+
      NoLegend()+NoAxes()+theme(plot.title = element_text(hjust = 0.5)))
  }
  else
  {
    print(DimPlot(PDAC_DT2_cells,cells.highlight = colnames(PDAC_DT2_cells)[PDAC_DT2_cells$Patient2 %in% (n)])+
      ggtitle(sprintf('%s (N=%d)',n,high_low_end_SPINK1_patient_summary_filtered[n,'N_ductal_type_2']))+
      NoLegend()+NoAxes()+theme(plot.title = element_text(hjust = 0.5)))
  }
  dev.off()
}

# Find genes correlates with SPINK1 level in ductal cell type 2
setwd('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\final_results_v4')
PDAC_DT2_cells <- readRDS('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\reference_dataset_core\\PDAC_acinar_ductal_cells_15430_genes_Harmony_by_patients_v4.rds')
PDAC_DT2_cells <- subset(PDAC_DT2_cells,cell_type_v2_AD_harmony == 'Ductal cell type 2')
lognormalized_matrix <- GetAssayData(object = PDAC_DT2_cells, layer = 'data')
pairs_to_test <- cbind(rep('SPINK1',nrow(lognormalized_matrix)),rownames(lognormalized_matrix))
SPINK1_corr_pairs <- scran::correlatePairs(lognormalized_matrix,pairings = pairs_to_test)
write.csv(SPINK1_corr_pairs,'SPINK1_corr_pairs_ductal_type_2.csv')
#
SPINK1_corr_pairs_DEG_up <- SPINK1_corr_pairs[SPINK1_corr_pairs$FDR<=0.05 & SPINK1_corr_pairs$rho>0.1,]
writeLines(SPINK1_corr_pairs_DEG_up$gene2,'SPINK1_corr_pairs_DEG_up.txt')
SPINK1_corr_pairs_DEG_down <- SPINK1_corr_pairs[SPINK1_corr_pairs$FDR<=0.05 & SPINK1_corr_pairs$rho<=-0.1,]
writeLines(SPINK1_corr_pairs_DEG_down$gene2,'SPINK1_corr_pairs_DEG_down.txt')

# Find the intersection between stem cell enriched ATAC-RNA down DEGs and genes correlated with SPINK1 in ductal cell type 2
setwd('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\final_results_v4')
SPINK1_corr_pairs <- read.csv('SPINK1_corr_pairs_ductal_type_2.csv',header = TRUE,row.names = 1)
PDB_enrichR <- read.table('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\ATACseq_SPINK1_KO\\Results\\ATAC_RNA_down_ProteomicsDB_2020_table.txt', sep = '\t', header = TRUE)
TabulaMuris_enrichR <- read.table('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\ATACseq_SPINK1_KO\\Results\\ATAC_RNA_down_Tabula_Muris_table.txt', sep = '\t', header = TRUE)
CellMarker2021_enrichR <- read.table('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\ATACseq_SPINK1_KO\\Results\\ATAC_RNA_down_CellMarker_Augmented_2021_table.txt', sep = '\t', header = TRUE)

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
# Plot the distribution of each gene of all_stem_genes in PDAC_DT2_cells
setwd('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\final_results_v4\\all_stem_genes')
for(n in 1:length(all_stem_genes))
{
  gene_to_plot <- all_stem_genes[n]
  if(gene_to_plot %in% rownames(PDAC_DT2_cells))
  {
    png(paste0(gene_to_plot,'.png'),width = 500,height = 500,res = 150)
    print(FeaturePlot(PDAC_DT2_cells,features = gene_to_plot,order = TRUE)+NoAxes())
    dev.off()
  }
}
high_vs_low_markers <- readRDS('high_vs_low_markers.rds')
# 5 genes in all_stem_genes are not present in high_vs_low_markers
sum(!(all_stem_genes %in% rownames(high_vs_low_markers)))
# 3 out of these 5 genes are not present in the atlas
# Therefore, another 2 are not tested due to low expression level 
sum(!(all_stem_genes %in% rownames(PDAC_DT2_cells)))
# Check how many genes in all_setm_genes are up/down in high SPINK1 end vs low SPINK1 end
high_vs_low_markers_all_stem_genes <- na.omit(high_vs_low_markers[all_stem_genes,])
sum(high_vs_low_markers_all_stem_genes$avg_log2FC>=1 & high_vs_low_markers_all_stem_genes$p_val_adj <=0.05)
sum(high_vs_low_markers_all_stem_genes$avg_log2FC<=-1 & high_vs_low_markers_all_stem_genes$p_val_adj <=0.05)
high_vs_low_markers_all_stem_genes[high_vs_low_markers_all_stem_genes$avg_log2FC>=1 
                                   & high_vs_low_markers_all_stem_genes$p_val_adj <=0.05,]
# Plot 6 genes with well-established functions in stem cell
p1 <- FeaturePlot(PDAC_DT2_cells,features = c('ACACA'),max.cutoff = 2,order = TRUE)&NoAxes()
p2 <- FeaturePlot(PDAC_DT2_cells,features = c('BCL7C'),order = TRUE)&NoAxes()
p3 <- FeaturePlot(PDAC_DT2_cells,features = c('YES1'),max.cutoff = 2,order = TRUE)&NoAxes()
png('Stem cell related genes set 1.png',width = 3000,height = 1000,unit = 'px',res = 300)
print(p1|p2|p3)
dev.off()
p1 <- FeaturePlot(PDAC_DT2_cells,features = c('ASPM'),order = TRUE)&NoAxes()
p2 <- FeaturePlot(PDAC_DT2_cells,features = c('HMMR'),order = TRUE)&NoAxes()
p3 <- FeaturePlot(PDAC_DT2_cells,features = c('NES'),max.cutoff = 2,order = TRUE)&NoAxes()
png('Stem cell related genes set 2.png',width = 3000,height = 1000,unit = 'px',res = 300)
print(p1|p2|p3)
dev.off()


