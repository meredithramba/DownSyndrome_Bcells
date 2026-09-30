library(patchwork)
library(tidyverse)
library(RColorBrewer)
library(Seurat)

batchNames_list <- list.dirs(path = "data/", recursive = F, full.names = F)

#Create Seurat objects for each sample
for (x in batchNames_list){
  
  print(x)
  
  cts <- ReadMtx(mtx = paste0('data/',x,'/sample_filtered_feature_bc_matrix/matrix.mtx.gz'),
                 features = paste0('data/',x,'/sample_filtered_feature_bc_matrix/features.tsv.gz'),
                 cells = paste0('data/',x,'/sample_filtered_feature_bc_matrix/barcodes.tsv.gz'))
  
  seuratObject <- CreateSeuratObject(counts = cts, min.cells = 5)
  seuratObject[["percent.mt"]] <- PercentageFeatureSet(seuratObject, pattern = "^MT-")  
  
  #unfiltered plots
  pdf(file = paste0('figures/filtering/', x,'_unfiltered_violin.pdf'))
  print(VlnPlot(seuratObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  dev.off() 
  
  pdf(paste0('figures/filtering/', x, '_unfiltered_correlation.pdf'))
  print(FeatureScatter(seuratObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm"))
  dev.off()
  
  #pre-processed data
  print(seuratObject)
  
  seuratObject <- subset(seuratObject, 
                         subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 20 & nCount_RNA < 25000 & nCount_RNA > 2000)
  
  #processed data
  print(seuratObject)
  
  #rename seuratObject with sample name
  assign(paste0(x,"_filt"), seuratObject)

  #filtered plots
  pdf(paste0('figures/filtering/', x, '_filtered_violin.pdf'))
  print(VlnPlot(seuratObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  dev.off()
  
  pdf(paste0('figures/filtering/', x, '_filtered_correlation.pdf'))
  print(FeatureScatter(seuratObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm"))
  dev.off()
  
}

#merge samples, adds in sample name in cell ID
mergeBcells <- merge(DM001_filt, y = c(DM002_filt, DM003_filt, DM004_filt, 
		     DM005_filt, DM006_filt, DM007_filt, DM008_filt, DM009_filt,
		     DM0011_filt),
                     add.cell.ids = c("DM001", "DM002", "DM003", "DM004", 
		     "DM005", "DM006","DM007", "DM008", "DM009", "DM0011"), 
                     project = "B_cells")

#add in metadata column with sample name
mergeBcells$sample <- rownames(mergeBcells@meta.data)
mergeBcells@meta.data <- separate(mergeBcells@meta.data, col = 'sample', into = c("Sample", "Barcode"), 
                                  sep = "_")
mergeBcells@meta.data <- select(mergeBcells@meta.data, -c("Barcode")) #remove unnecessary column

#add in group, age, sex, and batch information
sample_data <- read.csv("sample_data.csv")
mergeBcells@meta.data <- left_join(rownames_to_column(mergeBcells@meta.data), sample_data, by = "Sample")
rownames(mergeBcells@meta.data) <- mergeBcells@meta.data[,1]
mergeBcells@meta.data[,1] <- NULL

#normalization & integration on sample
mergeBcells <- NormalizeData(mergeBcells)
mergeBcells <- FindVariableFeatures(mergeBcells)
mergeBcells <- ScaleData(mergeBcells, vars.to.regress = "percent.mt")
mergeBcells <- RunPCA(mergeBcells)

options(future.globals.maxSize = 20000 * 1024^2)

#harmony integration
mergeBcells <- IntegrateLayers(object = mergeBcells, method = HarmonyIntegration, orig.reduction = "pca", 
                              new.reduction = "harmony", verbose = FALSE)

saveRDS(mergeBcells, file = "bcells_integrated.rds")
