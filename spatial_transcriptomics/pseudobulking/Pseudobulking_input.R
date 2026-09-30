library(Matrix)
library(SeuratObject)
library(Seurat)
library(tidyr)
library(reshape2)
library(plyr)
library(dplyr)
library(pbapply)
library(textTinyR)

folder = "../pseudobulking/input_files"

## read in RDS file
spatial.sct <- readRDS("../spatial_rds/sct_processed.rds")

## pull out raw data (not normalized)
DefaultAssay(spatial.sct) <- "Spatial"

## split the whole object by cluster assignment
spatial_clusters <- SplitObject(spatial.sct, split.by = "seurat_clusters")

## loop through clusters and save as R file
for (i in 1:length(spatial_clusters)){
  name <- names(spatial_clusters)[i]
  saveRDS(spatial_clusters[[i]], file = paste0(folder, "/cluster_", name, ".rds"))
}

clusters <- c("0","1","2","3","4","5","6","7","8","9")

## calculate pseudobulk sums per sample for the clusters
for (i in 1:length(clusters)){
  name <- clusters[i]
  
  ## read in data (raw UMI counts stored in Spatial slot) for each cluster
  data <- readRDS(paste0(folder, "/cluster_", name, ".rds"))
  
  ## get raw UMI counts
  raw_data_sparse <- GetAssayData(data, assay = "Spatial", slot = "counts")
  meta_data <- data@meta.data
  sample_colname <- "orig.ident"
  
  ## sample names
  IDs <- as.data.frame(meta_data)[, sample_colname]
  unique_ID_list <- as.list(unique(IDs))
  
  ## calculate pseudobulk by summing across all clusters
  pseudobulk <- as.data.frame(pblapply(unique_ID_list, FUN = function(x){sparse_Sums(raw_data_sparse[,IDs == x, drop = FALSE], rowSums = TRUE)}))
  cellcount <- pblapply(unique_ID_list, FUN = function(x){ncol(raw_data_sparse[,IDs == x, drop = FALSE])})
  colnames(pseudobulk) <- names(cellcount) <- unique_ID_list
  rownames(pseudobulk) <- rownames(x = data)
  
  saveRDS(pseudobulk, file = paste0(folder, "/cluster_",name,"_pseudobulk.rds"))
}
