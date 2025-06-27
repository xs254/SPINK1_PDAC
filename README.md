# SPINK1_PDAC
This repository is a collection of all codes used for analysis of scRNAseq,RNAseq,ATACseq in the project of [SPINK1 PDAC](https://www.biorxiv.org/content/10.1101/2025.03.03.641251v1)
## A modified scRNAseq atlas of PDAC patients
The atlas is originally from [Ryota et al.](https://pubmed.ncbi.nlm.nih.gov/35847558/), but during the analysis of this ATLAS, we noticed multiple mistakes:
- log normalized counts is mistakenly used as raw counts for ~50% of cells. 
- Some genes have different names likely due to different datasets are aligned to different versions of reference (e.g. CXCL8 and IL8 both refer to the same gene). During merging, Seruat treats IL8 and CXCL8 as two different genes and a dataset with CXCL8 will be treated as having 0 IL8 for all cells and vice versa.
- The existence of "date" genes likely due to some datasets being opend and saved in Excel which triggers auto conversion for genes whose name look like dates (e.g. SEPT1 to 1-SEPT).
- Several clusters of cells are mislabelled by the label transfer algorithm . 
Therefore, we first corrected these mistakes (see [scPDAC_all_cell_correction_and_integration](scRNASeq/scPDAC_all_cell_correction_and_integration.R)
and [scPDAC_gene_symbol_update_and_merging](scRNASeq/scPDAC_gene_symbol_update_and_merging.R)) and then performed analysis focusing on SPINK1 (see [scPDAC_analysis_v4](scRNASeq/scPDAC_analysis_v4.R)) .
## Integrated analysis of RNA and ATAC seq for SPINK1 KO PANC-1 cell line
We also analyzed change in transcription and chromatin accessibility in a SPINK1 KO PDAC cell line (PANC-1) and identified a group of stemness related genes with coherent change in both transcription level and chromatin accessibility.
- See [DEseq_for_SPINK1](Bulk_RNA_ATAC/DEseq_for_SPINK1.R) for the analysis of RNAseq data.
- See [ATACseq_SPINK1_KO_v2](Bulk_RNA_ATAC/ATACseq_SPINK1_KO_v2.R) for the analysis of ATAC and the intergated analysis of RNA and ATAC. 
