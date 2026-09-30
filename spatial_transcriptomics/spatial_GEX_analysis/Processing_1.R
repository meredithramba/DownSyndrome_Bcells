library(Seurat)
library(dplyr)
library(tidyverse)

#create list of samples in my folder - can also just make a list manually
batchNames_list <- list.dirs(path = '../data/', recursive = F, full.names = F)

#Create Seurat objects for each sample & filter
for (x in batchNames_list){
  seuratObject <- Load10X_Spatial(data.dir = paste0('../data/',x), filename = "filtered_feature_bc_matrix.h5", 
                                  assay = "Spatial")
  seuratObject[["percent.mt"]] <- PercentageFeatureSet(seuratObject, pattern = "^MT-") #add column with % of mito gene usage
  
  pdf(file = paste0("../figures/spatial_qc/", x, "_unfilt_vln.pdf"))
  print(VlnPlot(seuratObject, 
                features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3))
  dev.off()
  
  pdf(file = paste0("../figures/spatial_qc/", x, "_unfilt_scatter.pdf"))
  print(FeatureScatter(seuratObject, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial", cols = "black") + 
          ggtitle(x)) 
  dev.off()
  
  #prefiltered data
  print(seuratObject)
  
  #filtering (less than 200 unique genes, less than 10% mito genes)
  seuratObject <- subset(seuratObject, 
                         subset = nFeature_Spatial > 200 & percent.mt < 10)
  
  pdf(file = paste0("../figures/spatial_qc/", x, "_filt_vln.pdf"))
  print(VlnPlot(seuratObject, 
                features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3))
  dev.off()
  
  pdf(file = paste0("../figures/spatial_qc/", x, "_filt_scatter.pdf"))
  print(FeatureScatter(seuratObject, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial", cols = "black") + 
          ggtitle(x)) 
  dev.off()
  
  print(seuratObject)
  
  assign(paste0(x,"_filt"), seuratObject)
  
}

seuratObjectlist_filt = c(DUBO23_46127_0_filt, DUBO23_46131_0_filt, 
                          DUBO23_46142_0_filt, DUBO23_46152_0_filt, 
                          DUBO23_46156_0_filt, DUBO23_46167_0_filt)

#adding in label for each sample so they can be identified after merging
DUBO23_46127_0_filt$orig.ident <- "46127"
DUBO23_46131_0_filt$orig.ident <- "46131"
DUBO23_46142_0_filt$orig.ident <- "46142"
DUBO23_46152_0_filt$orig.ident <- "46152"
DUBO23_46156_0_filt$orig.ident <- "46156"
DUBO23_46167_0_filt$orig.ident <- "46167"

#merge samples 
mergespatial <- merge(DUBO23_46127_0_filt, y = c(DUBO23_46131_0_filt, 
                                                 DUBO23_46142_0_filt, DUBO23_46152_0_filt, 
                                                 DUBO23_46156_0_filt, DUBO23_46167_0_filt),
                      add.cell.ids = c("DUBO23_46127_filt", "DUBO23_46131_filt", "DUBO23_46142_filt",
                                       "DUBO23_46152_filt", "DUBO23_46156_filt", "DUBO23_46167_filt"), 
                      project = "Spatial")


#creating new column with group status
spatial_meta = mergespatial@meta.data

for (i in 1:nrow(spatial_meta)){
  if (spatial_meta[i, "orig.ident"] == 46127 | spatial_meta[i, "orig.ident"] == 46131 | spatial_meta[i, "orig.ident"] == 46142){
    spatial_meta[i, "Group"] <- "T21"}
  else{
    spatial_meta[i, "Group"] <- "D21"}
}

mergespatial <- AddMetaData(mergespatial, spatial_meta)

Idents(mergespatial) <- mergespatial$orig.ident

library(harmony)

#SCtransform and Harmony
spatial.list <- SplitObject(mergespatial, split.by="orig.ident") #split by each sample

spatial.list <- lapply(X = spatial.list, FUN = SCTransform, assay = "Spatial", new.assay.name = "SCT", return.only.var.genes = FALSE) #SCTransform

var.features <- SelectIntegrationFeatures(object.list = spatial.list, nfeatures = 3000) #integration features

spatial.sct <- merge(x = spatial.list[[1]], y = spatial.list[2:length(spatial.list)], merge.data=TRUE) #re-merge all samples

VariableFeatures(spatial.sct) <- var.features

spatial.sct <- RunPCA(spatial.sct, verbose = FALSE)
spatial.sct <- RunHarmony(spatial.sct, assay.use="SCT", group.by.vars = "Group") #integration
spatial.sct <- FindNeighbors(spatial.sct, reduction = "harmony", dims = 1:30)
spatial.sct <- FindClusters(spatial.sct, resolution = 0.5)
spatial.sct <- FindClusters(spatial.sct, resolution = 0.3)
spatial.sct <- RunUMAP(spatial.sct, reduction = "harmony", dims = 1:30)

saveRDS(spatial.sct, file = "../spatial_rds/sct_processed.rds")
