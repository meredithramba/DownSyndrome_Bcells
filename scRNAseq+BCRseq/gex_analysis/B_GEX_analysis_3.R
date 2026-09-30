library(dplyr)
library(patchwork)
library(tidyverse)
library(RColorBrewer)
library(Seurat)

# read in rds file from B_GEX_analysis_2
Bcells_integrated <- readRDS("bcells_integrated_processed.rds")

# Define color scheme
safe_colorblind_palette <- c("#88CCEE", "#5570b6","#A020F0", "#af97c8", 
                             "#AA4499", "#44AA99", "#24b34b",
                             "#6699CC", "#025a48" , "#888888","#ffc81c")

# dotplot
Idents(Bcells_integrated) <- factor(Idents(Bcells_integrated), levels=c("0","4","3","7","5","8","2","1","6","10"))

pdf("figures/dotplot_unlabeled.pdf", width = 14)
DotPlot(Bcells_integrated, features = c("BACH2", "TCL1A", "IL4R",
                                        "IGHD" , "IGHM", "STAT1",
                                        "MX1", "ISG15", "IFI44L",
                                        "CD69","MYC", "FOS", "CXCR3",
                                        "ITGAX", "FCRL5", "ZEB2", "MIR155HG", "NME1", "TNFRSF18",
                                        "CD27", "AIM2", "CD86", "IGHG2", "IGHA1","IGHG1",
				        "PRDM1", "XBP1", "JCHAIN", "MZB1"),  dot.scale = 10) +
  RotatedAxis() +
  theme(axis.text.x = element_text(face = 'italic')) +
  guides(color = guide_colorbar(title = 'Average Expression')) +
  scale_colour_gradient(low =c("white"), high =c("#332288"))
dev.off()

# label all clusters
new.cluster.ids <- c("Naive #1", "Naive #2", "ISG+ Naive", "Activated Naive",
                     "CXCR3+ Atypical B cells", "CXCR3+ Atypical B cells", "miR-155+ B cells", 
                     "Memory", "CD86+ Switched Memory", "Antibody Secreting Cells")

names(new.cluster.ids) <- levels(Bcells_integrated)
Bcells_integrated <- RenameIdents(Bcells_integrated, new.cluster.ids)

# add into metadata
Bcells_integrated[["cluster_label"]] <- Idents(Bcells_integrated)

## Supplemental S4 CD11c score
cd11c_features <- list(c("ITGAX", "ZEB2", "FCRL5"))

Bcells_integrated <- AddModuleScore(
  object = Bcells_integrated,
  features = cd11c_features,
  name = 'cd11c_score'
)

pdf("figures/cd11c.pdf")
FeaturePlot(Bcells_integrated, features = "cd11c_score1", min.cutoff = 0, max.cutoff = 1.5) +
  scale_colour_gradient(low =c("lightgrey"), high =c("#332288")) + ggtitle("CD11c score") +
  theme(plot.title = element_text(size = 20, face = "bold"))
dev.off()

# umap
pdf("figures/umap_precd11c_labeled.pdf")
DimPlot(Bcells_integrated, label = TRUE, cols = safe_colorblind_palette, label.size = 4, label.box = TRUE) +
        NoLegend()
dev.off()

# rename predicted cell type if CD11c score is high enough
for (i in 1:length(Cells(Bcells_integrated))){
  sample <- colnames(Bcells_integrated)[i]
  if (Bcells_integrated@meta.data[sample, ]["cd11c_score1"] > 0.75){
    Bcells_integrated@meta.data[sample, ]["cluster_label"] <- "CXCR3+ Atypical B cells"}
  
}

Idents(Bcells_integrated) <- "cluster_label"

### Figure 4C - labeled umap
pdf("figures/umap_labeled.pdf")
DimPlot(Bcells_integrated, label = TRUE, cols = safe_colorblind_palette, 
        label.size = 5, label.box = TRUE, repel = TRUE) +
  NoLegend()
dev.off()

## Supplemental S4 - umap by group
pdf("figures/umap_labeled_group.pdf", width = 10)
DimPlot(Bcells_integrated, label = TRUE, split.by = "Group", cols = safe_colorblind_palette, 
        label.size = 4, label.box = TRUE, repel = TRUE) +
        NoLegend()
dev.off()

### Figure 4B - dotplot labeled
pdf("figures/dotplot_labeled.pdf", width = 14)
DotPlot(Bcells_integrated, features = c("BACH2", "TCL1A", "IL4R",
                                        "IGHD" , "IGHM", "STAT1",
                                        "MX1", "ISG15", "IFI44L",
                                        "CD69","MYC", "FOS", "CXCR3",
                                        "ITGAX", "FCRL5", "ZEB2", "MIR155HG", "NME1", "TNFRSF18",
                                        "CD27", "AIM2", "CD86", "IGHG2", "IGHA1","IGHG1",
                                        "PRDM1", "XBP1", "JCHAIN", "MZB1"),  dot.scale = 10) +
  RotatedAxis() +
  theme(axis.text.x = element_text(face = 'italic')) +
  guides(color = guide_colorbar(title = 'Average Expression')) +
  scale_colour_gradient(low =c("white"), high =c("#332288"))
dev.off()

#cell %
Bcells_integrated_split <- SplitObject(Bcells_integrated, split.by = "Sample")

print("DM001")
prop.table(table(Idents(Bcells_integrated_split[["DM001"]])))
print("DM002")
prop.table(table(Idents(Bcells_integrated_split[["DM002"]])))
print("DM003")
prop.table(table(Idents(Bcells_integrated_split[["DM003"]])))
print("DM004")
prop.table(table(Idents(Bcells_integrated_split[["DM004"]])))
print("DM005")
prop.table(table(Idents(Bcells_integrated_split[["DM005"]])))
print("DM006")
prop.table(table(Idents(Bcells_integrated_split[["DM006"]])))
print("DM007")
prop.table(table(Idents(Bcells_integrated_split[["DM007"]])))
print("DM008")
prop.table(table(Idents(Bcells_integrated_split[["DM008"]])))
print("DM009")
prop.table(table(Idents(Bcells_integrated_split[["DM009"]])))
print("DM0011")
prop.table(table(Idents(Bcells_integrated_split[["DM0011"]])))

saveRDS(Bcells_integrated, file = "bcells_integrated_processed_labeled.rds")

