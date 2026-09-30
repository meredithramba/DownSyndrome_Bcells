library(dplyr)
library(patchwork)
library(tidyverse)
library(RColorBrewer)
library(Seurat)

#read in rds file from B_GEX_analysis 
Bcells_integrated <- readRDS("bcells_integrated.rds")

safe_colorblind_palette <- c("#88CCEE", "#af97c8", "#5570b6", "#A020F0","#AA4499",
                             "#44AA99", "#661100","#24b34b",
                             "#6699CC", "#025a48" , "#888888","#ffc81c")

# clustering
Bcells_integrated <- FindNeighbors(Bcells_integrated, dims = 1:30, reduction = "harmony")
Bcells_integrated <- FindClusters(Bcells_integrated, resolution = 0.3, cluster.name = "harmony_clusters")
Bcells_integrated <- RunUMAP(Bcells_integrated, dims = 1:30, reduction = "harmony", reduction.name = "umap.harmony")

#join layers for marker testing
Bcells_integrated <- JoinLayers(Bcells_integrated)

pdf("figures/UMAP_prefiltered.pdf")
DimPlot(Bcells_integrated, label = TRUE, cols = safe_colorblind_palette)
dev.off()

pdf("figures/UMAP_group_prefiltered.pdf")
DimPlot(Bcells_integrated, group.by = "Group", cols = safe_colorblind_palette)
dev.off()

pdf("figures/UMAP_sample_prefiltered.pdf")
DimPlot(Bcells_integrated, group.by = "Sample", cols = safe_colorblind_palette)
dev.off()

pdf("figures/UMAP_batch_prefiltered.pdf")
DimPlot(Bcells_integrated, group.by = "Batch", cols = safe_colorblind_palette)
dev.off()

# remove extreme points - they don't cluster with rest of cells
Bcells_integrated <- subset(x = Bcells_integrated, subset = umapharmony_1 < 10)

pdf("figures/UMAP_filtered.pdf")
DimPlot(Bcells_integrated, label = TRUE, cols = safe_colorblind_palette)
dev.off()

pdf("figures/UMAP_group_filtered.pdf")
DimPlot(Bcells_integrated, group.by = "Group", cols = safe_colorblind_palette)
dev.off()

pdf("figures/UMAP_sample_filtered.pdf")
DimPlot(Bcells_integrated, group.by = "Sample", cols = safe_colorblind_palette)
dev.off()

pdf("figures/UMAP_batch_filtered.pdf")
DimPlot(Bcells_integrated, group.by = "Batch", cols = safe_colorblind_palette)
dev.off()

#find markers of each cluster
Bcells_integrated@misc$markers <- FindAllMarkers(Bcells_integrated, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)

Bcells_integrated@misc$markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

pdf("figures/heatmap_top10.pdf", height = 20, width = 15)
DoHeatmap(Bcells_integrated, features = top10$gene) + NoLegend()
dev.off()

saveRDS(Bcells_integrated, file = "bcells_integrated_processed.rds")
