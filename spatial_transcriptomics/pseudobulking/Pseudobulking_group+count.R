library(Seurat)
library(limma)
library(edgeR)
library(tidyr)
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)
library(data.table)
library(pbapply)
library(Matrix)
library(textTinyR)
library(RColorBrewer) 
library(ggrepel) 
library(EnhancedVolcano)
library(qvalue)

folder = "../pseudobulking"

clusters <- c("0","1","2","3","4","5","6","7","8","9")

## read in meta data csv
md <- read.table(paste0("../pseudobulking/metadata.csv"), sep = ",", header = TRUE, row.names = 1)
md$orig.ident <- as.character(md$orig.ident)
md$Group <- factor(md$Group, levels = c("D21","T21"))
md$Cluster <- factor(md$Cluster, levels = c("0","1","2","3","4","5","6","7","8","9"))
md$nCount_Spatial <- scale(md$nCount_Spatial)

### ANALYSIS PART 1
for (i in 1:length(clusters)){
  
  cluster_i <- clusters[i]
  meta_data_i <- md
  
  ## read in pseudobulk
  reads <- readRDS(paste0("../pseudobulking/input_files/cluster_",cluster_i,"_pseudobulk.rds"))
  
  ## subset meta_data
  meta_data_i <- subset(meta_data_i, Cluster == cluster_i)
  row.names(meta_data_i) <- meta_data_i$orig.ident
  
  ## remove lowly-expressed genes
  dge <- DGEList(counts = reads)
  dge <- calcNormFactors(dge)
  design = model.matrix(~ Group + nCount_Spatial, data = meta_data_i)
  
  ## remove columns that are all 0s
  design <- design[, colSums(design != 0) > 0]
  v <- voom(dge, design, plot = FALSE)
  
  ## median logCPM filtering
  tab = data.frame(genes = rownames(reads), medians=apply(v$E, 1, median), order = 1:nrow(reads))
  tab = tab[order(-tab$medians), ]
  tab$order_by_median = 1:nrow(tab)
  tab = tab[order(tab$order), ]
  reads <- reads[which(tab$medians > 2), ]
  
  ## voom after removal of lowly expressed genes
  dge <- DGEList(counts = reads)
  dge <- calcNormFactors(dge)
  design = model.matrix(~ Group + nCount_Spatial, data = meta_data_i)
  v <- voom(dge, design, plot = FALSE)
  
  expression <- v$E
  weights <- v$weights
  colnames(weights) <- colnames(expression)
  rownames(weights) <- rownames(expression)
  
  ## write expression 
  write.table(expression, paste0(folder,"/pseudobulk_expression_cluster_",cluster_i,"_group+count.txt"), quote = FALSE, sep = ",")
  
  ## write weights
  write.table(weights, paste0(folder,"/pseudobulk_weights_cluster_",cluster_i,"_group+count.txt"), quote = FALSE, sep = ",")
  print(paste0(nrow(dge), " genes"))
}

### ANALYSIS PART 2

#removed cluster 9 from further analysis - too few spots per cluster/not represented in all samples
clusters <- c("0","1","2","3","4","5","6","7","8")

for (i in 1:length(clusters)){
  
  cluster_i <- clusters[i]
  meta_data_i <- md
  
  ## subset meta_data
  meta_data_i <- subset(meta_data_i, Cluster == cluster_i)
  row.names(meta_data_i) <- meta_data_i$orig.ident
  
  ## read in corrected expression
  residuals <- read.table(paste0(folder,"/pseudobulk_expression_cluster_",cluster_i,"_group+count.txt"), header = TRUE, sep = ",", check.names = FALSE)
  residuals <- residuals[colnames(residuals) %in% (rownames(meta_data_i))]
  
  ## read in weights
  weights <- read.table(paste0(folder,"/pseudobulk_weights_cluster_",cluster_i,"_group+count.txt"), header = TRUE, sep = ",", check.names = FALSE)
  weights <- weights[colnames(weights) %in% (rownames(meta_data_i))]
  
  print(length(which(rownames(meta_data_i)!=colnames(residuals))))
  print(length(which(colnames(weights)!=colnames(residuals))))
  
  ## model differential expression
  design = model.matrix(~ Group + nCount_Spatial, data = meta_data_i)
  design <- design[, colSums(design != 0) > 0]
  
  vfit <- lmFit(residuals, weights = weights, design)
  vfit <- eBayes(vfit)
  
  results <- topTable(vfit, coef = "GroupT21", num = Inf) 
  results$genes <- rownames(results)

  write.table(results, file = paste0(folder,"/cluster_",cluster_i,"_results_group+count.txt"), quote = FALSE, sep = ",")
  
  print(cluster_i)
  
}

cluster0_results <- read.table(file = paste0(folder,"/cluster_0_results_group+count.txt"), sep = ",")

pdf("../pseudobulking/cluster0_medullarysinus_pseudobulk.pdf", width = 8, height = 5)
EnhancedVolcano(cluster0_results, lab = cluster0_results$genes,
                x = "logFC", y = "P.Value", title = paste0("Medullary sinus (Medulla #1)"),
                subtitle = paste0("FDR 10%"),
                selectLab = c("MARCO"),
                pCutoff = 0.01912803,
                FCcutoff = 1.5,
                pointSize = 3.0,
                labSize = 6.0, 
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                col = c("#332288", "#332288","#332288", "#88CCEE"),
                boxedLabels = TRUE,
                colAlpha  = 0.75,
                legendPosition = "none")
dev.off()

cluster5_results <- read.table(file = paste0(folder,"/cluster_5_results_group+count.txt"), sep = ",")

pdf("../pseudobulking/cluster5_medullarycords_pseudobulk.pdf", width = 8, height = 5)
EnhancedVolcano(cluster5_results, lab = cluster5_results$genes,
                      x = "logFC", y = "P.Value", title = paste0("Medullary cords (Medulla #2)"),
                      subtitle = paste0("FDR 10%"),
                      selectLab = c("MARCO", "CLEC4M", "CD209", "SIGLEC1"),
                      pCutoff = 4.486682e-04,
                      FCcutoff = 1.5,
                      pointSize = 3.0,
                      labSize = 6.0, 
                      drawConnectors = TRUE,
                      widthConnectors = 0.75,
                      col = c("#332288", "#332288","#332288", "#88CCEE"),
                      boxedLabels = TRUE,
                      colAlpha  = 0.75,
                      legendPosition = "none")
dev.off()

cluster7_results <- read.table(file = paste0(folder,"/cluster_7_results_group+count.txt"), sep = ",")

pdf("../pseudobulking/cluster7_medulla+neutrophils_pseudobulk.pdf", width = 8, height = 5)
EnhancedVolcano(cluster7_results, lab = cluster7_results$genes,
                x = "logFC", y = "P.Value", title = paste0("Medulla with neutrophil infiltration (Medulla #3)"),
                subtitle = paste0("FDR 10%"),
                selectLab = c("MARCO", "SIGLEC1"),
                pCutoff = 9.38E-04,
                FCcutoff = 1.5,
                pointSize = 3.0,
                labSize = 6.0, 
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                col = c("#332288", "#332288","#332288", "#88CCEE"),
                boxedLabels = TRUE,
                colAlpha  = 0.75,
                legendPosition = "none")
dev.off()
