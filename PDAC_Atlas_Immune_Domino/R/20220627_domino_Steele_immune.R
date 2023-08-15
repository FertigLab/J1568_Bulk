library(Seurat)
library(domino)

sessionInfo()

# download cellphoneDB reference data
cell_db_path <- "data/cell_db"
if(!dir.exists(cell_db_path)){
  dir.create(cell_db_path, recursive = TRUE)
}

# check if the cellphone db files have been downloaded
complexes_check <- file.exists(paste0(cell_db_path, "/complexes.csv"))
interactions_check <- file.exists(paste0(cell_db_path, "/interactions.csv"))
proteins_check <- file.exists(paste0(cell_db_path, "/proteins.csv"))
genes_check <- file.exists(paste0(cell_db_path, "/genes.csv"))

if(!complexes_check){
  system(
    'curl -O https://raw.githubusercontent.com/Teichlab/cellphonedb-data/master/data/sources/complex_curated.csv \
    mv complex_curated.csv data/cell_db/complexes.csv')}
if(!interactions_check){
  system(
  'curl -O https://raw.githubusercontent.com/Teichlab/cellphonedb-data/master/data/sources/interaction_curated.csv \
  mv interaction_curated.csv data/cell_db/interactions.csv')}
if(!proteins_check){
  system(
  'curl -O https://raw.githubusercontent.com/Teichlab/cellphonedb-data/master/data/sources/protein_curated.csv \
  mv protein_curated.csv data/cell_db/proteins.csv')}
if(!genes_check){
  system(
  'curl -O https://raw.githubusercontent.com/Teichlab/cellphonedb-data/master/data/gene_input.csv \
  mv gene_input.csv data/cell_db/genes.csv')}

# cellphoneDB files downloaded 1:00 pm 6/16/22
# cellphoneDB V2.0.0

# create results directory
result_dir <- "results/20220627_Steele_immune_NA_excluded_hg38"
if(!dir.exists(result_dir)){
	dir.create(result_dir)
}

## Steele et al immune cell Seurat object
seurat <- readRDS("PDAC_immune_Steele/Steele_immune_seurat_NA_excluded.rds")
if(seurat@active.assay == "SCT"){
  seurat <- RenameAssays(object = seurat, SCT = "RNA")
}

## Pysenic results data
auc <- read.table(paste0(result_dir, "/auc.csv"),
                    header = TRUE, row.names = 1,
                    stringsAsFactors = FALSE, sep = ',')
regulons <- paste0(result_dir, "/regulons.csv")

# seurat object information
counts <- seurat@assays$RNA@counts
z_scores <- as.matrix(seurat@assays$RNA@scale.data)
clusters <- as.factor(seurat$cell_type_TNFRSF9_class)

# Create domino object
domino <- create_domino(signaling_db = "data/cell_db",
              features = t(auc),
              counts = counts,
              z_scores = z_scores,
              clusters = clusters,
              df = regulons,
	      remove_rec_dropout = FALSE)
# save the unbuilt domino object for changing build parameters
saveRDS(domino, file = paste0(result_dir, "/Steele_immune_domino_unbuilt.rds"))

# building the object will be done locally to allow parameter editting

sessionInfo()
