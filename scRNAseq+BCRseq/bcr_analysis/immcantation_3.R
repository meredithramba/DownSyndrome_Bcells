suppressPackageStartupMessages(library(airr))
suppressPackageStartupMessages(library(alakazam))
suppressPackageStartupMessages(library(dowser))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(scoper))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(shazam))
library(tidyverse)
library(rstatix)

data_mut <- read.csv("data_mut.csv")

#merge sample info
sample_code <- read.csv("sample_data.csv")
data_mut <- merge(sample_code, data_mut, by.y = "sample_id")

#order samples into groups
sample_order <- c(
  "DM001",
  "DM003",
  "DM006",
  "DM009",
  "DM002",
  "DM004",
  "DM005",
  "DM007",
  "DM008",
  "DM0011"
)

data_mut <- data_mut %>%
  mutate(sample_id = factor(sample_id, levels = sample_order)) %>%
  arrange(sample_id)

## frequency of each Ig class, per sample
safe_colorblind_palette <- c("#24b34b","#af97c8", "#5570b6", "#A020F0","#AA4499",
                             "#661100", "#44AA99", "#88CCEE", 
                             "#6699CC", "#025a48" )

#create table with frequencies
ccall.df.sample <- table(data_mut$c_call, data_mut$sample_id) %>% 
  as.data.frame() %>% 
  dplyr::rename(IgH=Var1, Sample=Var2) %>%
  mutate(Sample = factor(Sample, levels = sample_order)) %>%
  arrange(Sample) 

write.table(ccall.df.sample, "igclass_sample.csv", sep = ",", row.names = FALSE)

brks <- c(0, 0.25, 0.5, 0.75, 1)

pdf("figures/immcantation/ig_class_sample.pdf", width = 6, height = 5.5)
ggplot(ccall.df.sample, aes(x=Sample, y=Freq, fill = IgH)) + 
  geom_bar(stat="identity", position="fill") + scale_fill_manual(values = safe_colorblind_palette) + 
  ggplot2::theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank()) + 
  scale_y_continuous(breaks = brks, labels = scales::percent(brks)) + 
  ylab("% Composition of Sample by IgH") 
dev.off()

## Supplemental S4 - frequency of each Ig class, per cluster

#order clusters 
cluster_order <- c(
  "Naive #1",
  "Naive #2",
  "Activated Naive",
  "ISG+ Naive",
  "CXCR3+ Atypical B cells",
  "miR-155+ B cells",
  "Memory",
  "CD86+ Switched Memory",
  "Antibody Secreting Cells"
)

#order data
data_mut <- data_mut %>%
  mutate(gex_annotation = factor(gex_annotation, levels = cluster_order)) %>%
  arrange(gex_annotation)

#create table with frequencies
ccall.df.gex <- table(data_mut$c_call, data_mut$gex_annotation) %>% 
  as.data.frame() %>% 
  dplyr::rename(IgH=Var1, GEX=Var2) 

pdf("figures/immcantation/ig_class_cluster.pdf", width = 4.5, height = 5.5)
ggplot(ccall.df.gex, aes(x=GEX, y=Freq, fill = IgH)) + 
  geom_bar(stat="identity", position="fill") + scale_fill_manual(values = safe_colorblind_palette) + 
  ggplot2::theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black", face = "bold"),
                 axis.text.y = element_text(size = 14, color = "black"),
                 title = element_text(size = 12, color = "black", face = "bold")) + 
  scale_y_continuous(breaks = brks, labels = scales::percent(brks)) + 
  ylab("% Composition of Cluster by IgH") + xlab(NULL)
dev.off()

### Figure 4E - frequency of each IgH, per cluster, split between D21 and T21

data_mut_D21 <- subset(data_mut, data_mut$Group == "D21")
data_mut_T21 <- subset(data_mut, data_mut$Group == "T21")

ccall.df.d21 <- table(data_mut_D21$c_call, data_mut_D21$gex_annotation, data_mut_D21$sample_id) %>% 
  as.data.frame() %>% 
  dplyr::rename(IgH=Var1, GEX=Var2, Sample=Var3) 

write.table(ccall.df.d21, "igclass_sample_cluster_D21.csv", sep = ",", row.names = FALSE)

pdf("figures/immcantation/ig_class_cluster_D21.pdf", width = 4.5, height = 6.5)
ggplot(ccall.df.d21, aes(x=GEX, y=Freq, fill = IgH)) + 
  geom_bar(stat="identity", position="fill") + scale_fill_manual(values = safe_colorblind_palette) + 
  ggplot2::theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black", face = "bold"),
                 axis.text.y = element_text(size = 14, color = "black"),
                 title = element_text(size = 12, color = "black", face = "bold")) + 
  scale_y_continuous(breaks = brks, labels = scales::percent(brks)) + 
  labs(x = "", y = "% Composition of Cluster by IgH", title = "D21 Clusters by IgH Composition")
dev.off()

ccall.df.t21 <- table(data_mut_T21$c_call, data_mut_T21$gex_annotation, data_mut_T21$sample_id) %>% 
  as.data.frame() %>% 
  dplyr::rename(IgH=Var1, GEX=Var2, Sample=Var3) 

write.table(ccall.df.t21, "igclass_sample_cluster_T21.csv", sep = ",", row.names = FALSE)

pdf("figures/immcantation/ig_class_cluster_T21.pdf", width = 4.5, height = 6.5)
ggplot(ccall.df.t21, aes(x=GEX, y=Freq, fill = IgH)) + 
  geom_bar(stat="identity", position="fill") + scale_fill_manual(values = safe_colorblind_palette) + 
  ggplot2::theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black", face = "bold"),
                 axis.text.y = element_text(size = 14, color = "black"),
                 title = element_text(size = 12, color = "black", face = "bold")) + 
  scale_y_continuous(breaks = brks, labels = scales::percent(brks)) + 
  labs(x = "", y = "% Composition of Cluster by IgH", title = "T21 Clusters by IgH Composition")
dev.off()

### Figure 5A - SHM in memory cluster
memory <- data_mut %>%
  filter(gex_annotation == "Memory")

pdf("figures/immcantation/memory_SHM.pdf")
ggplot(memory, aes(y = mu_freq, x = sample_id, fill = sample_id)) +
  geom_boxplot() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  labs(x = "", y = "Mutation frequency", fill = "Sample", title = "Memory") +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 14, color = "black", face = "bold"),
        title = element_text(size = 14, color = "black", face = "bold")) + 
  scale_fill_manual(values = safe_colorblind_palette)
dev.off()

mean_sample_memory <- data_mut %>%
  group_by(sample_id) %>%
  filter(gex_annotation == "Memory") %>%
  summarise(mean_val = mean(mu_freq)) 

print(mean_sample_memory)
print("Mean mu_freq per sample for Memory")

### Figure 5B - SHM in combined memory and switched memory clusters
mean_sample_memory_combo <- data_mut %>%
  group_by(sample_id) %>%
  filter(gex_annotation == "CD86+ Switched Memory" | gex_annotation == "Memory") %>%
  summarise(mean_val = mean(mu_freq)) 

print(mean_sample_memory_combo)
print("Mean mu_freq per sample for Memory + CD86 Memory")

## Supplemental S5 - SHM by cluster

pdf("figures/immcantation/cluster_SHM.pdf", width = 4.5, height = 6)
ggplot(data_mut, aes(y = mu_freq, x = gex_annotation, fill = gex_annotation)) +
  geom_boxplot() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  labs(x = "", y = "Mutation frequency", fill = "Cluster") +
  theme(axis.text.x = element_text(size = 12, color = "black", 
                                   face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.text = element_blank(),
        legend.title = element_blank(),
        axis.title = element_text(size = 12, color = "black", face = "bold")) + 
  scale_fill_manual(values = safe_colorblind_palette)
dev.off()

### Figure 5C - V gene usage
vh_gene <- countGenes(data_mut, gene="v_call", groups="sample_id", mode="gene", clone = "clone_id")
vh_gene <- merge(sample_code, vh_gene, by.y = "sample_id")

# filter genes that do not appear in both D21 and T21 groups 
vh_gene_filtered <- vh_gene %>%
  filter(gene != "IGHV3-35" & gene != "IGHV3-69-1" & gene != "IGHV4/OR15-8")

# list of genes
unique_genes <- unique(vh_gene_filtered$gene)

# add in 0s for clones that don't appear in samples
sample_D21 <- c(
  "DM001",
  "DM003",
  "DM006",
  "DM009"
)

for (x in sample_D21){
  for (y in unique_genes){
    sample <- subset(vh_gene_filtered, vh_gene_filtered$sample_id==x)
    
    if (any(sample$gene == y) == FALSE){
      vh_gene_filtered <- rbind(vh_gene_filtered, list(sample_id = x, Group = "D21",
                                     Age = "", Sex = "", Batch = "", 
                                     gene = y, clone_count = 0, clone_freq = 0)) 
    }
  }
}

sample_T21 <- c(
  "DM002",
  "DM004",
  "DM005",
  "DM007",
  "DM008",
  "DM0011"
)

for (x in sample_T21){
  for (y in unique_genes){
    sample <- subset(vh_gene_filtered, vh_gene_filtered$sample_id==x)
    
    if (any(sample$gene == y) == FALSE){
      vh_gene_filtered <- rbind(vh_gene_filtered, list(sample_id = x, Group = "T21",
                                     Age = "", Sex = "", Batch = "", 
                                     gene = y, clone_count = 0, clone_freq = 0)) 
    }
  }
}

# run t tests for each gene, comparing clone freq between D21 and T21, and sort by statistic
stat.test <- vh_gene_filtered %>%
  group_by(gene) %>%
  t_test(clone_freq ~ Group, error.as.na = TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance() %>%
  arrange(statistic)

# extract list of genes, sorted by most elevated in T21 to most elevated in D21
gene_list <- pull(stat.test, gene)

# arrange v genes by statistics & sample
vh_gene_sorted_filtered <- vh_gene_filtered %>%
  mutate(sample_id = factor(sample_id, levels = sample_order)) %>%
  mutate(gene = factor(gene, levels = gene_list)) %>%
  arrange(gene)

pdf("figures/immcantation/IGV_heatmap.pdf", width = 12, height = 4)
ggplot(vh_gene_sorted_filtered, aes(x = gene, y = sample_id, fill = clone_freq)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="blue")  + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(fill = "Clone Freq")
dev.off()

### Figure 5E - IGVH434 usage per sample, in memory cluster
vh_usage_mem <- data.frame()

for (x in sample_order){
  sample <- subset(data_mut, data_mut$sample_id==x & data_mut$gex_annotation=="Memory")
  
  sample$VH4 <- ((sum(sample$v_call=="IGHV4-34*01", na.rm = TRUE)/nrow(sample))*100)
  
  #pull first row to get value
  vh_usage_sample <- sample[1,]
  
  #remove unnecessary columns
  vh_usage_sample <- subset(vh_usage_sample, select = c(sample_id, VH4))
  
  vh_usage_mem <- rbind(vh_usage_mem, vh_usage_sample)
  
}

print(vh_usage_mem)
print("IGHV4-34 usage in memory")

