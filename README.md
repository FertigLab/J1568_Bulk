# J1568_Bulk

### Author:
### Joe Tandurella

Repository created to house code pertaining to the J1568 MultiOmics Paper

# PDAC_Atlas_Immune_Domino

### Authors:
### Jacob Mitchell, Melissa Lyman

## Order of Operations

scripts used for the analysis are stored in the PDAC_Atlas_Immune_Domino/R/ folder

### 1) 20220622_Steele_TNFRSF9_classifier_domino_preprocessing.rmd

Subsets the atlas data to Steele et al and sets the cutoff value for the TNFRSF9 classifier at 1.282 after log-normalization of the counts data. Scales the expression data to recompute PCA and UMAP embeddings with the subset data. Includes plots that were remade in a later script and a loom file for domino that was not used in the final version.

### 2) 20220624_Steele_Immune_DE_TNFRSF9_hi-lo.rmd

Carries out differential expression analysis between TNFRSF9_hi and TNFRSF9_lo cells from CD8+, Effector CD8+, and Tregs clusters. Carries out gene set enrichment analysis of the DE results ranked by log2 fold change using HALLMARK, KEGG, REACTOME, IMMUNESIG, GOBP, and GOMF gene sets from MsigDB v. 7.5.1 as well as a custome immune gene set. Data frames of all test results and only significant results saved in figures/20220624/

### 3) 20220624_B_Plot_GSEA_Results.rmd

Plots waterfall plots of significantly enriched pathways from the GSEA results calculated in the previous script

### 4) 20220627_Steele_Immune_Domino_hg38_reference_rockfish.rmd

Subsets the Steele et al data set to remove cells annotated as NA or in the small clusters 32 and 33 that could not be confidently annotated as a singulat cell type. Counts matrix is also saved as a .loom file to be passed to domino for intercellular interaction analysis.

### 5) Scripts for pyscenic and domino run on the RockFish computing cluster
#### 20220627_pyscenic_Steele_immune.sh
#### 20220627_domino_Steele_immune.sh
#### 20220627_domino_Steele_immune.R

Runs pyscenic on the counts matrix loom file and scaled expression values from the Seurat object using the pyscyenc singularity instance from aertslab (version 0.11.0) and human transcription factor and motif references from the hg38 reference genome. Domino is then run with the resulting regulons file and CellPhoneDB annotations of ligand-receptor interactions (version 2.0.0, Teichlab). The directory with the scripts and outputs from RockFish are in processed_data/20220627_Steele_immune_NA_excluded_hg38.

### 6) 20220627_Steele_Immune_Results_Plotting.rmd

Calculates and plots the proportions of cells belonging to the TNFSF9 groups. Plots the jitter plot of cells in each group and the UMAP embedding of the cells used for the domino analysis.

### 7) 20220629_Domino_Results_hg38ref_noNA.rmd

Plots results from the domino analysis as heat maps of incoming signals, gene networks in grid orientation, the heatmap of TF features expressed by each cluster, and the correlation of receptor and TF expression heatmap.

### 8) 20220721_Supplement_Tables_DE_GSEA_Results.rmd

Save tables of significant gene set enrichment analysis and differential expression results from the MAST test between TNFRSF9 hi and lo cells by cell type.
