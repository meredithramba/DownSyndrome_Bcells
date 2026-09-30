library(Seurat)
library(msigdbr)
library(escape)
library(msigdb)
library(GSVA)
library(UCell)
library(ggplot2)
library(dplyr)

## create gene list from hallmark gene sets
GS.hallmark <- msigdbr(species  = "Homo sapiens",  
                       category = "H")

msigdbr_list <- split(x = GS.hallmark$gene_symbol, f = GS.hallmark$gs_name)

names(msigdbr_list) <- gsub("_", "-",
  names(msigdbr_list), fixed = TRUE)

## read in seurat data
Bcells_integrated <- readRDS("bcells_integrated_processed_labeled.rds")

## create matrix 
enrichment.scores <- escape.matrix(Bcells_integrated, 
                                   method = "UCell",
                                   gene.sets = msigdbr_list, 
                                   groups = 1000, 
                                   min.size = 5)

metadata <- Bcells_integrated@meta.data

# make row names into first column
metadata$cell_id <- rownames(metadata)
rownames(metadata) <- NULL
metadata <- metadata[, c("cell_id", setdiff(names(metadata), "cell_id"))]

# make scores into df and make row names into first column
enrichment.scores.df <- as.data.frame(enrichment.scores)
enrichment.scores.df$cell_id <- rownames(enrichment.scores.df)
rownames(enrichment.scores.df) <- NULL
enrichment.scores.df <- enrichment.scores.df[, c("cell_id", setdiff(names(enrichment.scores.df), "cell_id"))]

# merge scores and metadata
enrichment.scores.df <- merge(enrichment.scores.df, metadata, by.y = "cell_id")

write.table(enrichment.scores.df, "ucell_scores.csv", sep = ",", row.names = FALSE)

## subset on memory cells
enrichment.scores.memory <- subset(enrichment.scores.df, enrichment.scores.df$cluster_label == "Memory")

### Figure 4F and supplemental S4 - ROS pathway score in memory
pdf("figures/ucell_ros_memory.pdf")
ggplot(enrichment.scores.memory, aes(y = `HALLMARK.REACTIVE.OXYGEN.SPECIES.PATHWAY`, x = Sample, fill = Sample)) +
  geom_boxplot() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  labs(x = "", y = "HALLMARK-REACTIVE-OXYGEN-SPECIES-PATHWAY Enrichment Score", fill = "Sample", title = "Memory Cluster") +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 14, color = "black", face = "bold"),
        title = element_text(size = 14, color = "black", face = "bold")) +   
  scale_fill_manual(values = safe_colorblind_palette)
dev.off()

mean_memory <- enrichment.scores.memory %>%
  group_by(Sample) %>%
  summarise(mean_val = mean(`HALLMARK.REACTIVE.OXYGEN.SPECIES.PATHWAY`)) 

print('HALLMARK-REACTIVE-OXYGEN-SPECIES-PATHWAY, Memory')
print(mean_memory)

### Figure 4G and supplemental S4 - PI3K pathway score in memory
pdf("figures/ucell_pi3k_memory.pdf")
ggplot(enrichment.scores.memory, aes(y = `HALLMARK.PI3K.AKT.MTOR.SIGNALING`, x = Sample, fill = Sample)) +
  geom_boxplot() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  labs(x = "", y = "HALLMARK-PI3K-AKT-MTOR-SIGNALING Enrichment Score", fill = "Sample", title = "Memory Cluster") +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 14, color = "black", face = "bold"),
        title = element_text(size = 14, color = "black", face = "bold")) +   
  scale_fill_manual(values = safe_colorblind_palette)
dev.off()

mean_memory <- enrichment.scores.memory %>%
  group_by(Sample) %>%
  summarise(mean_val = mean(`HALLMARK.PI3K.AKT.MTOR.SIGNALING`))

print('HALLMARK-PI3K-AKT-MTOR-SIGNALING, Memory')
print(mean_memory)

## subset on CD86+ switched memory cells
enrichment.scores.memory86 <- subset(enrichment.scores.df, enrichment.scores.df$cluster_label == "CD86+ Switched Memory")

### Figure 4F and supplemental S4 - ROS pathway score in memory
pdf("figures/ucell_ros_memory86.pdf")
ggplot(enrichment.scores.memory86, aes(y = `HALLMARK.REACTIVE.OXYGEN.SPECIES.PATHWAY`, x = Sample, fill = Sample)) +
  geom_boxplot() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  labs(x = "", y = "HALLMARK-REACTIVE-OXYGEN-SPECIES-PATHWAY Enrichment Score", fill = "Sample", title = "CD86+ Switched Memory Cluster") +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 14, color = "black", face = "bold"),
        title = element_text(size = 14, color = "black", face = "bold")) +   
  scale_fill_manual(values = safe_colorblind_palette)
dev.off()

mean_memory86 <- enrichment.scores.memory86 %>%
  group_by(Sample) %>%
  summarise(mean_val = mean(`HALLMARK.REACTIVE.OXYGEN.SPECIES.PATHWAY`)) 

print('HALLMARK-REACTIVE-OXYGEN-SPECIES-PATHWAY, Switched Memory')
print(mean_memory86)

### Figure 4G and supplemental S4 - PI3K pathway score in memory
pdf("figures/ucell_pi3k_memory86.pdf")
ggplot(enrichment.scores.memory86, aes(y = `HALLMARK.PI3K.AKT.MTOR.SIGNALING`, x = Sample, fill = Sample)) +
  geom_boxplot() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  labs(x = "", y = "HALLMARK-PI3K-AKT-MTOR-SIGNALING Enrichment Score", fill = "Sample", title = "CD86+ Switched Memory Cluster") +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 14, color = "black", face = "bold"),
        title = element_text(size = 14, color = "black", face = "bold")) +   
  scale_fill_manual(values = safe_colorblind_palette)
dev.off()

mean_memory86 <- enrichment.scores.memory86 %>%
  group_by(Sample) %>%
  summarise(mean_val = mean(`HALLMARK.PI3K.AKT.MTOR.SIGNALING`)) 

print('HALLMARK-PI3K-AKT-MTOR-SIGNALING, Switched Memory')
print(mean_memory86)

