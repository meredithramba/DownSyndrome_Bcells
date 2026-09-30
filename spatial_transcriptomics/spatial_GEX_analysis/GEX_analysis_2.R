library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)
library(forcats)
library(scales)

library(pheatmap)
library(ggdendro)
library(cowplot)

spatial.sct <- readRDS("../spatial_rds/sct_processed.rds")

# Define color scheme
safe_colorblind_palette <- c("#88CCEE", "#af97c8", "#5570b6", "#ffc81c", "#44AA99","#24b34b", "#AA4499" ,
                             "#025a48" , "#A020F0", "#661100", "#6699CC", "#888888")

DefaultAssay(spatial.sct) <- "SCT"
Idents(spatial.sct) <- spatial.sct$seurat_clusters

# prep SCT values, as SCT was run on each individual sample
spatial.sct <- PrepSCTFindMarkers(spatial.sct)

# find markers
all_markers <- FindAllMarkers(spatial.sct, assay = "SCT")

# filter top 10 genes for each cluster
all_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

## Supplemental S2 - heatmap of clusters
pdf("../figures/heatmap.pdf", height = 10, width = 8)
DoHeatmap(spatial.sct, features = top10$gene, group.colors = safe_colorblind_palette) + 
  scale_fill_gradient(low = "white", high = "#332288") 
dev.off() 

# Combine clusters 3 and 8, label all clusters
Idents(spatial.sct) <- factor(Idents(spatial.sct), levels=c("0", "1", "2", "3", "4", "5", "6",
                                                            "7", "8", "9"))

new.cluster.ids <- c("Medullary sinus (Medulla #1)", "Capsule & Muscle", "Red blood cells/vascularity",
                     "Paracortex/T cell zones", "Subcapsular sinus", "Medullary cords (Medulla #2)",
                     "B cell follicles", "Medulla with neutrophil infiltration (Medulla #3)", 
                     "Paracortex/T cell zones", "Medulla #4")

names(new.cluster.ids) <- levels(spatial.sct)
spatial.sct <- RenameIdents(spatial.sct, new.cluster.ids)

# add labels into metadata
spatial.sct[["tissue_type"]] <- Idents(object = spatial.sct)

# adjust coordinates in images to integers so they graph properly
image_names <- c('slice1', 'slice1_DUBO23_46131_filt.1', 'slice1_DUBO23_46142_filt.2', 
                 'slice1_DUBO23_46152_filt.3', 'slice1_DUBO23_46156_filt.4', 
                 'slice1_DUBO23_46167_filt.5')

# change from character to int type
for (i in image_names) {
  
  # dynamic indexing
  coordinates <- spatial.sct@images[[i]]@coordinates
  
  # change type
  coordinates$imagerow <- as.numeric(coordinates$imagerow)
  coordinates$imagecol <- as.numeric(coordinates$imagecol)
  coordinates$row <- as.numeric(coordinates$row)
  coordinates$col <- as.numeric(coordinates$col)
  
  # assignment
  spatial.sct@images[[i]]@coordinates <- coordinates
}

### Fig 2A - example spatial data

pdf("../figures/46152_spatialplot.pdf", height = 6, width = 6)
SpatialDimPlot(spatial.sct, images = "slice1_DUBO23_46152_filt.3",  pt.size.factor = 3.5) + 
  scale_fill_manual(values=safe_colorblind_palette) +
  theme(legend.text = element_blank())
dev.off()

## Supplemental S3 cluster plots
pdf("../figures/46127_spatialplot.pdf", height = 6, width = 6)
SpatialDimPlot(spatial.sct, images = "slice1", pt.size.factor = 3.5) + 
  scale_fill_manual(values=safe_colorblind_palette) 
dev.off()

pdf("../figures/46131_spatialplot.pdf", height = 6, width = 6)
SpatialDimPlot(spatial.sct, images = "slice1_DUBO23_46131_filt.1",pt.size.factor = 3.5) + 
  scale_fill_manual(values=safe_colorblind_palette) 
dev.off()

pdf("../figures/46142_spatialplot.pdf", height = 6, width = 6)
SpatialDimPlot(spatial.sct, images = "slice1_DUBO23_46142_filt.2", pt.size.factor = 3.5) + 
  scale_fill_manual(values=safe_colorblind_palette) 
dev.off()

pdf("../figures/46152_spatialplot_2.pdf", height = 6, width = 6)
SpatialDimPlot(spatial.sct, images = "slice1_DUBO23_46152_filt.3",  pt.size.factor = 3.5) + 
  scale_fill_manual(values=safe_colorblind_palette)
dev.off()

pdf("../figures/46156_spatialplot.pdf", height = 6, width = 6)
SpatialDimPlot(spatial.sct, images = "slice1_DUBO23_46156_filt.4",  pt.size.factor = 4) + 
  scale_fill_manual(values=safe_colorblind_palette) 
dev.off()

pdf("../figures/46167_spatialplot.pdf", height = 6, width = 6)
SpatialDimPlot(spatial.sct, images = "slice1_DUBO23_46167_filt.5", pt.size.factor = 2) + 
  scale_fill_manual(values=safe_colorblind_palette) 
dev.off()

### Fig 2B - cluster identity

#reorder
Idents(spatial.sct) <- factor(Idents(spatial.sct), levels=c("Paracortex/T cell zones", "B cell follicles", 
                                                            "Medulla #4", "Medulla with neutrophil infiltration (Medulla #3)",
                                                            "Medullary sinus (Medulla #1)", "Medullary cords (Medulla #2)", 
                                                            "Red blood cells/vascularity", "Subcapsular sinus", "Capsule & Muscle"))

pdf("../figures/dotplot.pdf", height = 6, width = 14)
DotPlot(spatial.sct, features = c("MYH11", "MCAM", "ACTA2", "COL1A1", "COL1A2","COL3A1",
                                  "NTS", "TNFRSF9", "CXCL5", "ACKR1", 
                                  "HBD", "HBA2", "HBB", 
                                  "MARCO", "LYVE1", "CD209", "CLEC4M",
                                  "IGKC", "IGHG1", "IGHA1", 
                                  "S100A8", "S100A9", "MPO", "AZU1", "CTSG", 
                                  "CXCL13", "CD79A", "CD79B", "FCER2", "MS4A1", "CD22", "CD37",
                                  "IL7R","CCL21","CCL19", "CCR7", "TRBC1","TRBC2"), dot.scale = 10) + 
  RotatedAxis() + 
  theme(axis.text.x = element_text(face = 'italic')) +
  scale_y_discrete(labels = c("Paracortex/\n T cell zones", "B cell follicles", 
                              "Medulla #4", "Medulla with \nneutrophil infiltration\n(Medulla #3)",
                              "Medullary cords \n (Medulla #2)", "Medullary sinus \n (Medulla #1)",  
                              "Red blood cells/\nvascularity", "Subcapsular sinus", "Capsule & Muscle")) + 
  guides(color = guide_colorbar(title = 'Average Expression')) + 
  scale_colour_gradient(low =c("white"), high =c("#332288"))
dev.off()

### Figure 2C - UMAP

#reorder
Idents(spatial.sct) <- factor(Idents(spatial.sct), levels=c("Capsule & Muscle", "Subcapsular sinus", "Red blood cells/vascularity",
                                                            "Medullary sinus (Medulla #1)", "Medullary cords (Medulla #2)", 
                                                            "Medulla with neutrophil infiltration (Medulla #3)", "Medulla #4",
                                                            "B cell follicles","Paracortex/T cell zones" ))

safe_colorblind_palette_reordered <- c("#af97c8", "#44AA99", "#5570b6", "#88CCEE", "#24b34b", 
                                       "#025a48", "#A020F0", "#AA4499", "#ffc81c")

pdf("../figures/UMAP_group.pdf", width = 12)
DimPlot(spatial.sct, reduction = "umap", split.by = "Group", cols = safe_colorblind_palette_reordered)
dev.off()

### Figure 2D - BCL6 boxplot

# subset data only on B cell follicles
Idents(spatial.sct) <- "tissue_type"
bcell_follicles <- subset(x = spatial.sct, idents = "B cell follicles")
Idents(bcell_follicles) <- "orig.ident"

bcl6 <- DotPlot(bcell_follicles, features = c("BCL6")) 

for (i in 1:nrow(bcl6$data)){
  if (bcl6$data[i,"id"] == 46127 | bcl6$data[i,"id"] == 46131 | bcl6$data[i,"id"] == 46142){
    bcl6$data[i, "Group"] <- "T21"}
  else{
    bcl6$data[i, "Group"] <- "D21"}
}

bcl6$data$Group <- factor(bcl6$data$Group, levels = c("D21", "T21"))

print(bcl6$data)

pdf("../figures/BCL6_boxplot.pdf", width = 5, height = 4.5)
ggplot(bcl6$data, aes(x=features.plot,y=pct.exp, fill = Group)) + scale_fill_manual(values = c("#a3a6d2", "#332288")) +
  geom_boxplot() + geom_dotplot(binaxis='y', dotsize=1, stackdir='center', position = position_jitterdodge(jitter.width  = 0.01)) +
  ylab("% B Cell Follicle Spots Expressing BCL6") + ylim(0,50) + theme_bw() + xlab(NULL) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(size = 14, color = "black", face = "italic"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        panel.background = element_blank()) + theme(text = element_text(size = 15))
dev.off()

### Figure 2E - B cell migration/cytokine boxplot

bcell_migration <- DotPlot(bcell_follicles, features = c("CXCR5", "CXCL13")) 

for (i in 1:nrow(bcell_migration$data)){
  if (bcell_migration$data[i,"id"] == 46127 | bcell_migration$data[i,"id"] == 46131 | bcell_migration$data[i,"id"] == 46142){
    bcell_migration$data[i, "Group"] <- "T21"}
  else{
    bcell_migration$data[i, "Group"] <- "D21"}
}

bcell_migration$data$Group <- factor(bcell_migration$data$Group, levels = c("D21", "T21"))

print(bcell_migration$data)

pdf("../figures/B cell migration boxplot.pdf", width = 8, height = 4.5)
ggplot(bcell_migration$data,aes(x=features.plot,y=pct.exp, fill = Group)) + scale_fill_manual(values = c("#a3a6d2", "#332288")) +
  geom_boxplot() + geom_dotplot(binaxis='y', dotsize=1, stackdir='center', position = position_jitterdodge(jitter.width  = 0.01)) +
  ylab("% B Cell Follicle Spots Expressing Genes") + ylim(0,100) + theme_bw() + xlab(NULL) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(size = 14, color = "black", face = "italic"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        panel.background = element_blank()) + theme(text = element_text(size = 15))
dev.off()

### Figure 2F - T cell migration/cytokine boxplot

# subset data only on T cell zones
Idents(spatial.sct) <- "tissue_type"
tcell_zones <- subset(x = spatial.sct, idents = "Paracortex/T cell zones")
Idents(tcell_zones) <- "orig.ident"

tcell_migration <- DotPlot(tcell_zones, features = c("CCR7", "CCL19", "CCL21")) 

for (i in 1:nrow(tcell_migration$data)){
  if (tcell_migration$data[i,"id"] == 46127 | tcell_migration$data[i,"id"] == 46131 | tcell_migration$data[i,"id"] == 46142){
    tcell_migration$data[i, "Group"] <- "T21"}
  else{
    tcell_migration$data[i, "Group"] <- "D21"}
}

tcell_migration$data$Group <- factor(tcell_migration$data$Group, levels = c("D21", "T21"))

print(tcell_migration$data)

pdf("../figures/T cell migration boxplot.pdf", width = 8, height = 4.5)
ggplot(tcell_migration$data,aes(x=features.plot,y=pct.exp, fill = Group)) + scale_fill_manual(values = c("#a3a6d2", "#332288")) +
  geom_boxplot() + geom_dotplot(binaxis='y', dotsize=1, stackdir='center', position = position_jitterdodge(jitter.width  = 0.01)) +
  ylab("% T Cell Zone Spots Expressing Genes") + ylim(0,100) + theme_bw() + xlab(NULL) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(size = 14, color = "black", face = "italic"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        panel.background = element_blank()) + theme(text = element_text(size = 15))
dev.off()

## Supplemental S3 - frequency of each cluster, per sample

#order clusters 
cluster_order <- c(
  "Capsule & Muscle", "Subcapsular sinus", "Red blood cells/vascularity",
  "Medullary sinus (Medulla #1)", "Medullary cords (Medulla #2)", 
  "Medulla with neutrophil infiltration (Medulla #3)", "Medulla #4",
  "B cell follicles","Paracortex/T cell zones")

sample_order <- c(
  "46152", "46156", "46167",
  "46127", "46131", "46142")

#order data
spatial.freq <- as.data.frame(table(spatial.sct@meta.data$tissue_type))

Regions.df <- table(spatial.sct@meta.data$tissue_type, spatial.sct@meta.data$orig.ident) %>% 
  as.data.frame() %>% 
  dplyr::rename(Spatial_Region=Var1, Sample=Var2) %>%
  mutate(Spatial_Region = factor(Spatial_Region, levels = cluster_order)) %>%
  arrange(Spatial_Region) %>%
  mutate(Sample = factor(Sample, levels = sample_order)) %>%
  arrange(Sample)

brks <- c(0, 0.25, 0.5, 0.75, 1)

pdf("../figures/barplot.pdf", width = 7, height = 8)
ggplot(Regions.df,aes(x=Sample,y=Freq, fill = Spatial_Region)) + 
  geom_bar(stat="identity", position="fill") + scale_fill_manual(values = safe_colorblind_palette_reordered) + 
  ggplot2::theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank()) + 
  scale_y_continuous(breaks = brks, labels = scales::percent(brks)) + 
  ylab("% Composition of Sample by Spatial Region / Clusters") 
dev.off()

### Figure 3A

# adding metadata column with tissue type and cluster
spatial.sct$group.sample <- paste(spatial.sct$Group, spatial.sct$tissue_type, sep = " ")

Idents(spatial.sct) <- "group.sample"

Idents(spatial.sct) <- factor(Idents(spatial.sct), levels=c("D21 Medullary sinus (Medulla #1)", "T21 Medullary sinus (Medulla #1)", 
                                                            "D21 Capsule & Muscle", "T21 Capsule & Muscle", 
                                                            "D21 Red blood cells/vascularity", "T21 Red blood cells/vascularity", 
                                                            "D21 Subcapsular sinus", "T21 Subcapsular sinus", 
                                                            "D21 Medullary cords (Medulla #2)", "T21 Medullary cords (Medulla #2)", 
                                                            "D21 Medulla with neutrophil infiltration (Medulla #3)", "T21 Medulla with neutrophil infiltration (Medulla #3)", 
                                                            "D21 Medulla #4", "T21 Medulla #4",
                                                            "D21 Paracortex/T cell zones", "T21 Paracortex/T cell zones",
                                                            "D21 B cell follicles", "T21 B cell follicles"))

pdf("../figures/dotplot_stromal.pdf", width = 12, height = 6)
DotPlot(spatial.sct, features = c("ICAM1", "VCAM1", 
                                  "MARCO", "LYVE1", "CD209", "CLEC4M",
                                  "CR1", "CR2",
                                  "PDPN", "ATXN2", "OTUB1", "PAXBP1", "BACH1"), 
        idents = c("D21 Medullary sinus (Medulla #1)", "T21 Medullary sinus (Medulla #1)",
                   "D21 Medullary cords (Medulla #2)", "T21 Medullary cords (Medulla #2)", 
                   "D21 B cell follicles", "T21 B cell follicles",
                   "D21 Medulla with neutrophil infiltration (Medulla #3)", 
                   "T21 Medulla with neutrophil infiltration (Medulla #3)", 
                   "D21 Medulla #4", "T21 Medulla #4"), dot.scale = 10) + 
  RotatedAxis() + 
  scale_y_discrete(labels = label_wrap(40)) +
  theme(axis.text.x = element_text(face = 'italic')) +
  guides(color = guide_colorbar(title = 'Average Expression')) + 
  scale_colour_gradient(low =c("white"), high =c("#332288"))
dev.off()

### Figure 3B - Boxplot of adhesion markers
Idents(spatial.sct) <- "orig.ident"

stromal <- DotPlot(spatial.sct, features = c("ICAM1", "VCAM1")) 

for (i in 1:nrow(stromal$data)){
  if (stromal$data[i,"id"] == 46127 | stromal$data[i,"id"] == 46131 | stromal$data[i,"id"] == 46142){
    stromal$data[i, "Group"] <- "T21"}
  else{
    stromal$data[i, "Group"] <- "D21"}
}

stromal$data$Group <- factor(stromal$data$Group, levels = c("D21", "T21"))

print(stromal$data)

pdf("../figures/boxplot_adhesion.pdf", width = 4.5, height = 6)
ggplot(stromal$data,aes(x=features.plot,y=pct.exp, fill = Group)) + scale_fill_manual(values = c("#a3a6d2", "#332288")) +
  geom_boxplot() + geom_dotplot(binaxis='y', dotsize=1, stackdir='center', position = position_jitterdodge(jitter.width  = 0.01)) +
  ylab("% Spots Expressing Genes") + ylim(0,100) + theme_bw() + xlab(NULL) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(size = 14, color = "black", face = "italic"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        panel.background = element_blank()) + theme(text = element_text(size = 15))
dev.off()

### Figure 3C - Boxplot of FDC markers

gene1 <- "CR1"
gene2 <- "CR2"

# Create a data frame indicating if each gene is expressed in each spot
# This code has the same results as DotPlot() in terms of pct.exp, except this code only
# calculates spots that express both genes
expr_matrix <- data.frame(
  gene1_expr = FetchData(spatial.sct, vars = gene1)[[1]] > 0.5,
  gene2_expr = FetchData(spatial.sct, vars = gene2)[[1]] > 0.5
)

## add cell IDs to the data frame
expr_matrix$cell_IDs <- rownames(spatial.sct@meta.data)
expr_matrix$sample <- spatial.sct@meta.data$orig.ident

samples = c("46127", "46131", "46142", "46152", "46156", "46167")

results = NULL

for (i in 1:length(samples)){
  sample_num = samples[i]
  
  # Calculate number of spots expressing all genes
  spots <- filter(expr_matrix, sample == sample_num)
  
  count = 0
  
  for (i in 1:nrow(spots)){
    if (spots$gene1_expr[i] == TRUE & spots$gene2_expr[i] == TRUE) {
      count <- count + 1
    } 
  }
  percentage <- (count / nrow(spots)) * 100
  
  results = append(results, percentage)
}

Idents(spatial.sct) <- "orig.ident"

fdc_df = data.frame(pct.exp=results,
                    features.plot=c("CR1, CR2","CR1, CR2","CR1, CR2","CR1, CR2","CR1, CR2","CR1, CR2"),
                    id=c(46127, 46131, 46142, 46152, 46156, 46167))

for (i in 1:nrow(fdc_df)){
  if (fdc_df[i,"id"] == 46127 | fdc_df[i,"id"] == 46131 | fdc_df[i,"id"] == 46142){
    fdc_df[i, "Group"] <- "T21"}
  else{
    fdc_df[i, "Group"] <- "D21"}
}

fdc_df$Group <- factor(fdc_df$Group, levels = c("D21", "T21"))

print(fdc_df)

pdf("../figures/boxplot_fdc.pdf", width = 4.5, height = 6)
ggplot(fdc_df,aes(x=features.plot,y=pct.exp, fill = Group)) + scale_fill_manual(values = c("#a3a6d2", "#332288")) +
  geom_boxplot() + geom_dotplot(binaxis='y', dotsize=1, stackdir='center', position = position_jitterdodge(jitter.width  = 0.01)) +
  ylab("% Spots Expressing Genes") + ylim(0,15) + theme_bw() + xlab(NULL) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(size = 14, color = "black", face = "italic"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        panel.background = element_blank()) + theme(text = element_text(size = 15))
dev.off()

### Figure 3D - boxplot of MARCO LECs

gene1 <- "MARCO"
gene2 <- "LYVE1"
gene3 <- "CD209"
gene4 <- "CLEC4M"

# Create a data frame indicating if each gene is expressed in each spot
expr_matrix <- data.frame(
  gene1_expr = FetchData(spatial.sct, vars = gene1)[[1]] > 0.5,
  gene2_expr = FetchData(spatial.sct, vars = gene2)[[1]] > 0.5,
  gene3_expr = FetchData(spatial.sct, vars = gene3)[[1]] > 0.5,
  gene4_expr = FetchData(spatial.sct, vars = gene4)[[1]] > 0.5
)

## add cell IDs to the data frame
expr_matrix$cell_IDs <- rownames(spatial.sct@meta.data)
expr_matrix$sample <- spatial.sct@meta.data$orig.ident

samples = c("46127", "46131", "46142", "46152", "46156", "46167")

results = NULL

for (i in 1:length(samples)){
  sample_num = samples[i]
  
  # Calculate number of spots expressing all genes
  spots <- filter(expr_matrix, sample == sample_num)
  
  count = 0
  
  for (i in 1:nrow(spots)){
    if (spots$gene1_expr[i] == TRUE & spots$gene2_expr[i] == TRUE & 
        spots$gene3_expr[i] == TRUE & spots$gene4_expr[i] == TRUE) {
      count <- count + 1
    } 
  }
  percentage <- (count / nrow(spots)) * 100
  
  results = append(results, percentage)
}

Idents(spatial.sct) <- "orig.ident"

marco = data.frame(pct.exp=results,
                   features.plot=c("MARCO LECs (MARCO+ LYVE1+ CD209+ CLEC4M+)","MARCO LECs (MARCO+ LYVE1+ CD209+ CLEC4M+)",
                                   "MARCO LECs (MARCO+ LYVE1+ CD209+ CLEC4M+)","MARCO LECs (MARCO+ LYVE1+ CD209+ CLEC4M+)",
                                   "MARCO LECs (MARCO+ LYVE1+ CD209+ CLEC4M+)","MARCO LECs (MARCO+ LYVE1+ CD209+ CLEC4M+)"),
                   id=c(46127, 46131, 46142, 46152, 46156, 46167))

for (i in 1:nrow(marco)){
  if (marco[i,"id"] == 46127 | marco[i,"id"] == 46131 | marco[i,"id"] == 46142){
    marco[i, "Group"] <- "T21"}
  else{
    marco[i, "Group"] <- "D21"}
}

marco$Group <- factor(marco$Group, levels = c("D21", "T21"))

print(marco)

pdf("../figures/boxplot_marco.pdf", width = 4.5, height = 6)
ggplot(marco,aes(x=features.plot,y=pct.exp, fill = Group)) + scale_fill_manual(values = c("#a3a6d2", "#332288")) +
  geom_boxplot() + geom_dotplot(binaxis='y', dotsize=1, stackdir='center', position = position_jitterdodge(jitter.width  = 0.01)) +
  ylab("% Spots Expressing All Genes") + ylim(0,100) + theme_bw() + xlab(NULL) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(size = 14, color = "black", face = "italic"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        panel.background = element_blank()) + theme(text = element_text(size = 15))
dev.off()

### Figure 3E - MARCO LECs Spatial Plot
marco <- WhichCells(spatial.sct, expression= LYVE1 > 0.5 & MARCO > 0.5 & CD209 > 0.5 & CLEC4M > 0.5) 

pdf("../figures/46127_spatialplot_MARCO.pdf", height = 6, width = 6)
SpatialPlot(spatial.sct, images = "slice1", cells.highlight = marco,  pt.size.factor = 3.5,
            cols.highlight = c("lemonchiffon", "#332288"), label =FALSE)
dev.off()

pdf("../figures/46131_spatialplot_MARCO.pdf", height = 6, width = 6)
SpatialPlot(spatial.sct, images = "slice1_DUBO23_46131_filt.1", cells.highlight = marco,  pt.size.factor = 3.5,
            cols.highlight = c("lemonchiffon", "#332288"), label =FALSE) 
dev.off()

pdf("../figures/46142_spatialplot_MARCO.pdf", height = 6, width = 6)
SpatialPlot(spatial.sct, images = "slice1_DUBO23_46142_filt.2", cells.highlight = marco,  pt.size.factor = 3.5,
            cols.highlight = c("lemonchiffon", "#332288"), label =FALSE)
dev.off()

pdf("../figures/46152_spatialplot_MARCO.pdf", height = 6, width = 6)
SpatialPlot(spatial.sct, images = "slice1_DUBO23_46152_filt.3", cells.highlight = marco,  pt.size.factor = 3.5,
            cols.highlight = c("lemonchiffon", "#332288"), label =FALSE)
dev.off()

pdf("../figures/46156_spatialplot_MARCO.pdf", height = 6, width = 6)
SpatialPlot(spatial.sct, images = "slice1_DUBO23_46156_filt.4", cells.highlight = marco,  pt.size.factor = 4,
            cols.highlight = c("lemonchiffon", "#332288"), label =FALSE) 
dev.off()

pdf("../figures/46167_spatialplot_MARCO.pdf", height = 6, width = 6)
SpatialPlot(spatial.sct, images = "slice1_DUBO23_46167_filt.5", cells.highlight = marco,  pt.size.factor = 2,
            cols.highlight = c("lemonchiffon", "#332288"), label =FALSE) 
dev.off()

## Supplemental S3 - PDPN
frc <- DotPlot(spatial.sct, features = c("PDPN"))

for (i in 1:nrow(frc$data)){
  if (frc$data[i,"id"] == 46127 | frc$data[i,"id"] == 46131 | frc$data[i,"id"] == 46142){
    frc$data[i, "Group"] <- "T21"}
  else{
    frc$data[i, "Group"] <- "D21"}
}

frc$data$Group <- factor(frc$data$Group, levels = c("D21", "T21"))

print(frc$data)

pdf("../figures/boxplot_pdpn.pdf", width = 4.5, height = 6)
ggplot(frc$data,aes(x=features.plot,y=pct.exp, fill = Group)) + scale_fill_manual(values = c("#a3a6d2", "#332288")) +
  geom_boxplot() + geom_dotplot(binaxis='y', dotsize=1, stackdir='center', position = position_jitterdodge(jitter.width  = 0.01)) +
  ylab("% Spots Expressing Genes") + ylim(0,25) + theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(face = 'italic'),
        panel.background = element_blank()) + theme(text = element_text(size = 15))
dev.off()

## Supplemental S3 - stromal cell groups
Idents(spatial.sct) <- "group.sample"

Idents(spatial.sct) <- factor(Idents(spatial.sct), levels=c("D21 Medullary sinus (Medulla #1)", "T21 Medullary sinus (Medulla #1)", 
                                                            "D21 Capsule & Muscle", "T21 Capsule & Muscle", 
                                                            "D21 Red blood cells/vascularity", "T21 Red blood cells/vascularity", 
                                                            "D21 Paracortex/T cell zones", "T21 Paracortex/T cell zones", 
                                                            "D21 Subcapsular sinus", "T21 Subcapsular sinus", 
                                                            "D21 Medullary cords (Medulla #2)", "T21 Medullary cords (Medulla #2)", 
                                                            "D21 B cell follicles", "T21 B cell follicles",
                                                            "D21 Medulla with neutrophil infiltration (Medulla #3)", "T21 Medulla with neutrophil infiltration (Medulla #3)", 
                                                            "D21 Medulla #4", "T21 Medulla #4"))
# vLECs
pdf("../figures/dotplot_vLEC.pdf", width = 6, height = 4.5)
DotPlot(spatial.sct, features = c("FOXC2", "CLDN11"), 
        idents = c("D21 Medullary sinus (Medulla #1)", "T21 Medullary sinus (Medulla #1)",
                   "D21 Medullary cords (Medulla #2)", "T21 Medullary cords (Medulla #2)", 
                   "D21 Subcapsular sinus", "T21 Subcapsular sinus", 
                   "D21 Medulla with neutrophil infiltration (Medulla #3)", "T21 Medulla with neutrophil infiltration (Medulla #3)", 
                   "D21 Medulla #4", "T21 Medulla #4"), dot.scale = 10)  + 
  RotatedAxis() + 
  scale_y_discrete(labels = label_wrap(40)) +
  theme(axis.text.x = element_text(face = 'italic')) +
  guides(color = guide_colorbar(title = 'Average Expression')) + 
  scale_colour_gradient(low =c("white"), high =c("#332288"))
dev.off()

# cLECs
pdf("../figures/dotplot_cLEC.pdf", width = 8.5, height = 6)
DotPlot(spatial.sct, features = c("ACKR4", "MMRN1", "FOXC2", "PDGFB","EDN1", "CAV1", "TNFRSF11A", "NT5E"), 
        idents = c("D21 Medullary sinus (Medulla #1)", "T21 Medullary sinus (Medulla #1)",
                   "D21 Medullary cords (Medulla #2)", "T21 Medullary cords (Medulla #2)", 
                   "D21 Subcapsular sinus", "T21 Subcapsular sinus", 
                   "D21 Medulla with neutrophil infiltration (Medulla #3)", "T21 Medulla with neutrophil infiltration (Medulla #3)", 
                   "D21 Medulla #4", "T21 Medulla #4"), dot.scale = 10)  + 
  RotatedAxis() + 
  scale_y_discrete(labels = label_wrap(40)) +
  theme(axis.text.x = element_text(face = 'italic')) +
  guides(color = guide_colorbar(title = 'Average Expression')) + 
  scale_colour_gradient(low =c("white"), high =c("#332288"))
dev.off()

# fLECs
pdf("../figures/dotplot_fLEC.pdf", width = 8, height = 6)
DotPlot(spatial.sct, features = c("CCL20", "CD74", "CD274", "CSF1", "ACKR1"), 
        idents = c("D21 Medullary sinus (Medulla #1)", "T21 Medullary sinus (Medulla #1)",
                   "D21 Medullary cords (Medulla #2)", "T21 Medullary cords (Medulla #2)", 
                   "D21 Subcapsular sinus", "T21 Subcapsular sinus", 
                   "D21 Medulla with neutrophil infiltration (Medulla #3)", "T21 Medulla with neutrophil infiltration (Medulla #3)", 
                   "D21 Medulla #4", "T21 Medulla #4"), dot.scale = 10)  + 
  RotatedAxis() + 
  scale_y_discrete(labels = label_wrap(40)) +
  theme(axis.text.x = element_text(face = 'italic')) +
  guides(color = guide_colorbar(title = 'Average Expression')) + 
  scale_colour_gradient(low =c("white"), high =c("#332288"))
dev.off()

# PTX3-LECS
pdf("../figures/dotplot_PTX3_LEC.pdf", width = 9, height = 6)
DotPlot(spatial.sct, features = c("PTX3", "ITIH3", "LYVE1", "FLT4","NRP2", "PDPN", "CD36"), 
        idents = c("D21 Medullary sinus (Medulla #1)", "T21 Medullary sinus (Medulla #1)",
                   "D21 Medullary cords (Medulla #2)", "T21 Medullary cords (Medulla #2)", 
                   "D21 Subcapsular sinus", "T21 Subcapsular sinus", 
                   "D21 Medulla with neutrophil infiltration (Medulla #3)", "T21 Medulla with neutrophil infiltration (Medulla #3)", 
                   "D21 Medulla #4", "T21 Medulla #4"), dot.scale = 10)  + 
  RotatedAxis() + 
  scale_y_discrete(labels = label_wrap(40)) +
  theme(axis.text.x = element_text(face = 'italic')) +
  guides(color = guide_colorbar(title = 'Average Expression')) + 
  scale_colour_gradient(low =c("white"), high =c("#332288"))
dev.off()
