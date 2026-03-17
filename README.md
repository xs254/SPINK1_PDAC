# SPINK1_PDAC
This repository is a collection of all codes used for analysis of scRNA-seq, RNA-seq, ATAC-seq, Cut&Tag-seq and microscopy images in the manuscript [SPINK1-COL18A1 Crosstalk Shapes Epigenome and Drives Cancer Stemness in Pancreatic Ductal Adenocarcinoma](https://www.biorxiv.org/content/10.1101/2025.03.03.641251v1)
## A modified scRNAseq atlas of PDAC patients
The atlas is originally from [Ryota et al.](https://pubmed.ncbi.nlm.nih.gov/35847558/), but during the analysis of this ATLAS, we noticed multiple mistakes:
- log normalized counts is mistakenly used as raw counts for ~50% of cells. 
- Some genes have different names likely due to different datasets are aligned to different versions of reference (e.g. CXCL8 and IL8 both refer to the same gene). During merging, Seruat treats IL8 and CXCL8 as two different genes and a dataset with CXCL8 will be treated as having 0 IL8 for all cells and vice versa.
- The existence of "date" genes likely due to some datasets being opend and saved in Excel which triggers auto conversion for genes whose name look like dates (e.g. SEPT1 to 1-SEPT).
- Several clusters of cells are mislabelled by the label transfer algorithm .

Therefore, we first corrected these mistakes (see [scPDAC_all_cell_correction_and_integration](scRNASeq/scPDAC_all_cell_correction_and_integration.R)
and [scPDAC_gene_symbol_update_and_merging](scRNASeq/scPDAC_gene_symbol_update_and_merging.R)) and then performed analysis focusing on SPINK1 (see [scPDAC_analysis_v4](scRNASeq/scPDAC_analysis_v4.R)) . 

## Integrated analysis of RNA-seq, ATAC-seq, Cut&Tag-Seq for SPINK1 KO PANC-1 cell line
We also analyzed change in transcription (RNA-seq), chromatin accessibility (ATAC-seq), and histone modifications (H3K4Me3 and H3K27Me3 Cut&Tag-seq) in SPINK1-KO cells (SPINK1-OE and SPINK1-KO for RNA-seq). The integrated analysis identified a group of stemness related genes with coherent change in histone modifications, chromatin accessibility and transcription. 
- See [DEseq_for_SPINK1_TM59Z7](Bulk_RNA_ATAC/DEseq_for_SPINK1_TM59Z7.R) for DEseq of RNAseq data.
- See [Analysis_for_SPINK1_TM59Z7](Bulk_RNA_ATAC/Analysis_for_SPINK1_TM59Z7.R) for DEseq results analysis and integrated analysis with ATAC-seq.
- See [ATACseq_SPINK1_KO_v2](Bulk_RNA_ATAC/ATACseq_SPINK1_KO_v2.R) for the analysis of ATAC-seq data.
- See [CutNTag_csaw](CutNTag/CutNTag_csaw.R) for the differential sites calling by csaw.
- See [CutNtag_analysis](CutNTag/CutNTag_analysis.R) for csaw results analysis and intergated analysis with RNA-seq and ATAC-seq.
- See [RNA_ATAC_CutNTag_integrated_visualization_v2](CutNTag/RNA_ATAC_CutNTag_integrated_visualization_v2.R) for integrated visualization of RNA, ATAC, Cut&Tag results for exemplary genes.

## Pearson correlation analysis of pixel intensities from individual cell
We analyzed the colocalization of SPINK1, COL18A1 and GOLGA2 (golgi apparatus marker) by quantifying Pearson correlation of pixel intensities of individual segmented cells. See [correlation_for_colocalization](Image_Analysis/correlation_for_colocalization.m) for details.
