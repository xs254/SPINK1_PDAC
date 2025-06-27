library(Seurat)
library(monocle3)
library(ggplot2)
library(SeuratWrappers)
setwd('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\individual datasets collection')
# Load core genes and synonymous genes to be used
core_genes <- readRDS('Genes_shared_across_datasets.rds')
extra_core_genes <- readRDS('extra_core_genes.rds')
synonymous_symbol_group <- readRDS('synonymous_symbol_group.rds')
# Load raw counts data from PRJCA001063 and subset by core genes
PRJCA001063 <- readRDS('PRJCA001063_raw_counts.rds')
PRJCA001063 <- PRJCA001063[['RNA']]$counts
for(n in 1:length(synonymous_symbol_group))
{
  need_update <- synonymous_symbol_group[[n]]
  rownames(PRJCA001063)[rownames(PRJCA001063) %in% need_update] <- need_update[1] #The most updated gene symbol is the first one
}
PRJCA001063 <- PRJCA001063[c(core_genes,extra_core_genes),]
PRJCA001063 <- CreateSeuratObject(PRJCA001063,project = 'PRJCA001063')
saveRDS(PRJCA001063,'C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\reference_dataset_core\\PRJCA001063_15430_genes.rds')
# Load raw counts data from GSE111672 and subset by core genes
GSE111672 <- readRDS('GSE111672_raw_counts.rds')
GSE111672 <- JoinLayers(GSE111672)
GSE111672 <- GSE111672[['RNA']]$counts
for(n in 1:length(synonymous_symbol_group))
{
  need_update <- synonymous_symbol_group[[n]]
  rownames(GSE111672)[rownames(GSE111672) %in% need_update] <-  need_update[1]
}
GSE111672 <- GSE111672[c(core_genes,extra_core_genes),]
GSE111672 <- CreateSeuratObject(GSE111672,project = 'GSE111672')
saveRDS(GSE111672,'C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\reference_dataset_core\\GSE111672_15430_genes.rds')
# Load the original reference data set pk_all.rds and subset by core genes
PDAC_all_cells <- readRDS('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\reference_dataset\\pk_all.rds')
saveRDS(PDAC_all_cells@meta.data,'meta_data_of_all_cells_in_pk_all.rds')
PDAC_all_cells <- subset(PDAC_all_cells,Project!='CA001063' & Project!='GSE111672')
PDAC_all_cells_meta_data <- PDAC_all_cells@meta.data
PDAC_all_cells <- PDAC_all_cells[['RNA']]$counts
symbols_need_update <- unlist(synonymous_symbol_group)
symbols_absent <- symbols_need_update[!(symbols_need_update %in% rownames(PDAC_all_cells))]
symbols_absent_mtx <- Matrix::Matrix(0,nrow = length(symbols_absent),ncol = ncol(PDAC_all_cells),sparse = TRUE)
rownames(symbols_absent_mtx) <- symbols_absent
colnames(symbols_absent_mtx) <- colnames(PDAC_all_cells)
PDAC_all_cells <- rbind(PDAC_all_cells,symbols_absent_mtx)
extra_core_genes_mtx <- Matrix::Matrix(0,nrow = length(extra_core_genes),ncol = ncol(PDAC_all_cells),sparse = TRUE)
rownames(extra_core_genes_mtx) <- extra_core_genes
colnames(extra_core_genes_mtx) <- colnames(PDAC_all_cells)
for(n in 1:length(synonymous_symbol_group))
{
  need_update <- synonymous_symbol_group[[n]]
  # Sum up all gene counts to the row with the correct name
  # No need to delete the rows with outdated names because we only extract rows with most updated names later
  extra_core_genes_mtx[need_update[1],] <- colSums(PDAC_all_cells[need_update,])
  if(n%%10==0){print(n)}
}
PDAC_all_cells <- rbind(PDAC_all_cells[core_genes,],extra_core_genes_mtx)
PDAC_all_cells <- CreateSeuratObject(PDAC_all_cells,project = 'PDAC_all_cells_core',meta.data = PDAC_all_cells_meta_data)
saveRDS(PDAC_all_cells,'C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\reference_dataset_core\\pk_all_15430_genes.rds')
# The pk_all_core has two problems:
# 1) All raw counts from PRJCA001063 are from log normalized counts
# 2) All 'date' genes from GSE111672 have 0 counts.
# These two problems can be fixed by updating the raw counts
setwd('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\reference_dataset_core')
PRJCA001063_core <- readRDS('PRJCA001063_15430_genes.rds')
GSE111672_core <- readRDS('GSE111672_15430_genes.rds')
PDAC_all_cells_core <- readRDS('pk_all_15430_genes.rds')
pk_all_meta_data <- readRDS('meta_data_of_all_cells_in_pk_all.rds')
# Double check cell IDs across datasets are matched
PRJCA001063_cells <- rownames(pk_all_meta_data)[pk_all_meta_data$Project == 'CA001063']
GSE111672_cells <- rownames(pk_all_meta_data)[pk_all_meta_data$Project == 'GSE111672']
sum(!(PRJCA001063_cells %in% colnames(PRJCA001063_core)))
sum(!(GSE111672_cells %in% colnames(GSE111672_core)))
# Update the counts matrix of PDAC_all_cells_core
PRJCA001063_core$Cell_IDs <- colnames(PRJCA001063_core)
PRJCA001063_core <- subset(PRJCA001063_core,Cell_IDs %in% PRJCA001063_cells)
PRJCA001063_core@meta.data <- pk_all_meta_data[PRJCA001063_core$Cell_IDs,]
GSE111672_core$Cell_IDs <- colnames(GSE111672_core)
GSE111672_core <- subset(GSE111672_core,Cell_IDs %in% GSE111672_cells)
GSE111672_core@meta.data <- pk_all_meta_data[GSE111672_core$Cell_IDs,]
print(unique(PDAC_all_cells_core$Project))
PDAC_all_cells_core_final <- merge(PDAC_all_cells_core,y = c(PRJCA001063_core,GSE111672_core),project = 'PDAC_all_cells_core')
PDAC_all_cells_core_final <- JoinLayers(PDAC_all_cells_core_final)
print(unique(PDAC_all_cells_core_final$Project))
saveRDS(PDAC_all_cells_core_final,'PDAC_all_cells_15430_genes_final.rds')
# RPCA integration of PDAC_all_cells_core_final
setwd('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\reference_dataset_core')
PDAC_all_cells <- readRDS('PDAC_all_cells_15430_genes_final.rds')
PDAC_all_cells[['RNA']] <- split(PDAC_all_cells[["RNA"]], f = PDAC_all_cells$Project)
PDAC_all_cells <- NormalizeData(PDAC_all_cells, normalization.method = 'LogNormalize')
PDAC_all_cells <- FindVariableFeatures(PDAC_all_cells , selection.method = 'vst', nfeatures = 2000)
PDAC_all_cells <- ScaleData(PDAC_all_cells)
PDAC_all_cells <- RunPCA(PDAC_all_cells)
options(future.globals.maxSize = 5*1024^3)
PDAC_all_cells <- IntegrateLayers(object = PDAC_all_cells, method = RPCAIntegration,
                                      orig.reduction = 'pca', new.reduction = 'integrated_rpca',verbose = TRUE)
PDAC_all_cells <- JoinLayers(PDAC_all_cells)
saveRDS(PDAC_all_cells, 'PDAC_all_cells_15430_RPCA_by_project.rds')
