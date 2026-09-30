library(Seurat)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)
library(data.table)
library(pbapply)
library(tidyr)

folder = "../pseudobulking"

## read in RDS file
spatial.sct <- readRDS("../spatial_rds/sct_processed.rds")

## adding metadata column with group and cluster
spatial.sct$group.sample <- paste(spatial.sct$Group, spatial.sct$seurat_clusters, sep = "_")

Idents(spatial.sct) <- "group.sample"
Idents(spatial.sct) <- factor(Idents(spatial.sct), levels=c("D21_0", "T21_0", "D21_1", "T21_1", "D21_2", "T21_2", "D21_3",
                                                            "T21_3", "D21_4", "T21_4", "D21_5", "T21_5", "D21_6", "T21_6",
                                                            "D21_7", "T21_7", "D21_8", "T21_8", "D21_9", "T21_9"))

## creating a dataframe with information of # of spots per cluster/sample
spatial.sct$sample.cluster <- paste(spatial.sct$orig.ident, spatial.sct$seurat_clusters, spatial.sct$Group, sep = "_")
spatial.counts <- as.data.frame(table(spatial.sct@meta.data$sample.cluster))

## order in desired order
list_samples <- c("46127_0_T21","46127_1_T21", "46127_2_T21", "46127_3_T21", "46127_4_T21", "46127_5_T21", 
                  "46127_6_T21", "46127_7_T21", "46127_8_T21", "46127_9_T21", 
                  "46131_0_T21","46131_1_T21", "46131_2_T21", "46131_3_T21", "46131_4_T21", "46131_5_T21", 
                  "46131_6_T21", "46131_7_T21", "46131_8_T21", "46131_9_T21",
                  "46142_0_T21","46142_1_T21", "46142_2_T21", "46142_3_T21", "46142_4_T21", "46142_5_T21", 
                  "46142_6_T21", "46142_7_T21", "46142_8_T21", "46142_9_T21", 
                  "46152_0_D21","46152_1_D21", "46152_2_D21", "46152_3_D21", "46152_4_D21", "46152_5_D21", 
                  "46152_6_D21", "46152_7_D21", "46152_8_D21", "46152_9_D21", 
                  "46156_0_D21","46156_1_D21", "46156_2_D21", "46156_3_D21", "46156_4_D21", "46156_5_D21", 
                  "46156_6_D21", "46156_7_D21", "46156_8_D21", "46156_9_D21",
                  "46167_0_D21","46167_1_D21", "46167_2_D21", "46167_3_D21", "46167_4_D21", "46167_5_D21", 
                  "46167_6_D21", "46167_7_D21", "46167_8_D21", "46167_9_D21")

## transform model column into factor and use level to organize
spatial.counts$Var1 <- factor(spatial.counts$Var1, levels = list_samples) 

spatial.counts <- spatial.counts[order(spatial.counts$Var1), ]
rownames(spatial.counts) <- spatial.counts[,1]
spatial.counts[,1] <- NULL

## summing nCount_Spatial by subcluster
metadata <- spatial.sct@meta.data
ncounts_Spatial <- metadata %>% group_by(sample.cluster) %>% summarise(nCount_Spatial = sum(nCount_Spatial)) %>% as.data.frame()

## transform model column into factor and use level to organize
ncounts_Spatial$sample.cluster <- factor(ncounts_Spatial$sample.cluster, levels = list_samples) 

ncounts_Spatial <- ncounts_Spatial[order(ncounts_Spatial$sample.cluster), ]
rownames(ncounts_Spatial) <- ncounts_Spatial[,1]

## split into 3 columns
ncounts_Spatial <- ncounts_Spatial %>% separate("sample.cluster", c('orig.ident', 'Cluster', 'Group'))

## bind together to make metadata file
bulk_meta <- cbind(spatial.counts, ncounts_Spatial)

write.csv(bulk_meta, paste0(folder,"/metadata.csv"))
