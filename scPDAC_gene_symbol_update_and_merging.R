library(Seurat)
library(scCustomize)
setwd('C:\\Users\\xshan\\OneDrive\\Desktop\\R code\\PDAC_reference\\individual datasets collection')
# Load raw counts for PRJCA001063
PRJCA001063 <- read.table('PRJCA001063_count_matrix.txt')
PRJCA001063 <- CreateSeuratObject(PRJCA001063,project = 'PRJCA001063')
saveRDS(PRJCA001063,file = 'PRJCA001063_raw_counts.rds')
# Load raw counts for GSE111672. 
# This dataset needs some extra correction of gene names because some gene symbols are converted to date by excel
GSE111672_A <- read.delim('GSE111672_PDAC-A-indrop-filtered-expMat.txt',sep = '\t',header = TRUE)
GSE111672_B <- read.delim('GSE111672_PDAC-B-indrop-filtered-expMat.txt',sep = '\t',header = TRUE)
GSM3036909 <- read.delim('GSM3036909.tsv',sep = '\t',header = TRUE) # Dataset with unconverted gene names for reference from same study
date_gene_index_A <- which(!(GSE111672_A$Genes %in% GSM3036909$Genes))
date_gene_index_B <- which(!(GSE111672_B$Genes %in% GSM3036909$Genes))
print(GSE111672_A$Genes[date_gene_index_A])
print(GSE111672_B$Genes[date_gene_index_B])
print(sum(date_gene_index_A!=date_gene_index_B)) # Not surprisingly, the location are the same 
# The correct gene names in in GSM3036909, we just need to delete some mitochondrial and ribosomal genes
correct_gene_list <- (GSM3036909$Genes[!(GSM3036909$Genes %in% GSE111672_A$Genes)])
correct_gene_list <- correct_gene_list[-c(2,16:106,120)]
# Make a converting table for final check
converting_table <- data.frame(wrong_name = GSE111672_A$Genes[date_gene_index_A],
                               correct_name = correct_gene_list)
GSE111672_A$Genes[date_gene_index_A] <- correct_gene_list
GSE111672_B$Genes[date_gene_index_B] <- correct_gene_list
row.names(GSE111672_A) <- GSE111672_A$Genes
row.names(GSE111672_B) <- GSE111672_B$Genes
GSE111672_A <- GSE111672_A[,-1]
GSE111672_B <- GSE111672_B[,-1]
GSE111672_A <- CreateSeuratObject(GSE111672_A,project = 'GSE111672_A')
GSE111672_A$Patient <- 'GSE111672_A'
GSE111672_B <- CreateSeuratObject(GSE111672_B,project = 'GSE111672_B')
GSE111672_A$Patient <- 'GSE111672_B'
GSE111672 <- merge(GSE111672_A,GSE111672_B, add.cell.ids = c('GSE111672_PDAC_A','GSE111672_PDAC_B'))
GSE111672$Project <- 'GSE111672'
GSE111672$Type <- 'Tumor'
saveRDS(GSE111672,'GSE111672_raw_counts.rds')
# Load GSE154778
GSE154778 <- read.csv('GSE154778_dgeMtx.csv',row.names = 1,header = TRUE)
GSE154778 <- CreateSeuratObject(GSE154778,project = 'GSE154778')
saveRDS(GSE154778,'GSE154778_raw_counts.rds')
# Load GSM4293555
GSM4293555 <- read.delim('GSM4293555_Human.csv', sep = '\t',row.names = 1,header = TRUE)
GSM4293555 <- CreateSeuratObject(GSM4293555,project = 'GSM4293555')
saveRDS(GSM4293555,'GSM4293555_raw_counts.rds')
# Establish a consensus set of genes measured in all datasets
# Load genes measured in each project
PRJCA001063_genes <- rownames(readRDS('PRJCA001063_raw_counts.rds'))
GSE111672_genes <- rownames(readRDS('GSE111672_raw_counts.rds'))
GSE155698_genes <- rownames(Read10X('GSE155698_P1_filtered_feature_bc_matrix'))
GSE154778_genes <- rownames(readRDS('GSE154778_raw_counts.rds'))
GSM4293555_genes <- rownames(readRDS('GSM4293555_raw_counts.rds'))
gene_list <- list(PRJCA001063_genes,GSE111672_genes,GSE155698_genes,GSE154778_genes,GSM4293555_genes)
shared_genes <- Reduce(intersect,gene_list)
shared_genes <- shared_genes[order(shared_genes)]
saveRDS(shared_genes,'core_genes_shared_across_datasets.rds')
not_shared_genes <- vector('list',length(gene_list))
# Extract genes that does not present in all datasets
for(n in 1:length(gene_list))
{
  not_shared_genes[[n]] <- gene_list[[n]][!(gene_list[[n]] %in% shared_genes)]
}
# Update the names of all not shared genes
not_shared_genes_updated <- vector('list',length(gene_list))
for(n in 1:length(gene_list))
{
  updated <- Updated_HGNC_Symbols(not_shared_genes[[n]],case_check_as_warn = TRUE,verbose = FALSE)
  updated <- updated[,c(1,5)]
  uncertain_genes <- updated$input_features[duplicated(updated$input_features)]
  updated <- updated[!(updated$input_features %in% uncertain_genes),]
  updated <- rbind(updated,data.frame('input_features' = uncertain_genes,'Output_Features' = uncertain_genes))
  not_shared_genes[[n]] <- updated$input_features
  not_shared_genes_updated[[n]] <- updated$Output_Features
}
# After updating the name, another 468 genes are identified to be shared
shared_genes_after_update <- Reduce(intersect,not_shared_genes_updated)
# Pair up the old name(s) and the updated names for these 468 genes
synonymous_symbol_group <- vector('list',length(shared_genes_after_update))
for(n_gene in 1:length(shared_genes_after_update))
{
  # The first symbol in the group is the most updated one
  new_name <- shared_genes_after_update[n_gene]
  synonymous_symbol_group[[n_gene]] <- new_name
  for(n_list in 1:length(not_shared_genes_updated))
  {
    old_name <- not_shared_genes[[n_list]][not_shared_genes_updated[[n_list]] == new_name]
    # If multiple genes in the old name map on to one single new name, this is not a very clear case and 
    # this gene will be discarded
    if(length(old_name)>1)
    {
      synonymous_symbol_group[[n_gene]] <- NA
      break
    }
    else
    {
      if(!(old_name %in% synonymous_symbol_group[[n_gene]]))
      {
        synonymous_symbol_group[[n_gene]] <- c(synonymous_symbol_group[[n_gene]],old_name)
      }
    }
  }
}
synonymous_symbol_group <- synonymous_symbol_group[!is.na(synonymous_symbol_group)]
saveRDS(synonymous_symbol_group,'synonymous_symbol_group.rds')
# Save shared genes after update as extra core genes
extra_core_genes <- rep(NaN,length(synonymous_symbol_group))
for(n in 1:length(synonymous_symbol_group))
{
  extra_core_genes[n] <- synonymous_symbol_group[[n]][1]
}
saveRDS(extra_core_genes,'extra_core_genes.rds')
