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

safe_colorblind_palette <- c("#88CCEE", "#af97c8", "#5570b6", "#A020F0","#AA4499",
                             "#44AA99", "#661100","#24b34b",
                             "#6699CC", "#025a48" , "#888888","#ffc81c")

# read in merged tsv, see R script "tsv_merge.R"
bcr_data <- airr::read_rearrangement(file = 'BCR_DM001-DM011_data.tsv',
                                     aux_types = c("v_germline_length" = "i",
                                                   "d_germline_length" = "i",
                                                   "j_germline_length" = "i",
                                                   "day" = "i"))

cat(paste("There are", nrow(bcr_data), 
          "rows in `bcr_data`.\n"))

## There are 183184 rows in `bcr_data`.

# Filter out inconsistent sequences - edited this to c_call (from c_gene)
bcr_data <- bcr_data %>%
  dplyr::filter(
    (grepl("^IGHV", v_call) & grepl("^IGHJ", j_call) & grepl("^IGH[MGADE]", c_call)) |
      (grepl("^IGKV", v_call) & grepl("^IGKJ", j_call) & grepl("^IGKC", c_call)) |
      (grepl("^IGLV", v_call) & grepl("^IGLJ", j_call) & grepl("^IGLC", c_call))
  )
cat(paste("There are", nrow(bcr_data), 
          "rows in the data after filtering V/J/C calls inconsistent with the respective locus.\n"))

## There are 182613 rows in the data after filtering V/J/C calls inconsistent with the respective locus.

# read in the data
bcr_data <- bcr_data %>% dplyr::filter(productive)

cat(paste("There are", nrow(bcr_data), "rows in the data.\n"))

## There are 182611 rows in the data.

bcr_data %>% slice_sample(n = 5) # random examples

# remove cells with multiple heavy chain
multi_heavy <- table(dplyr::filter(bcr_data, locus == "IGH")$cell_id)
multi_heavy_cells <- names(multi_heavy)[multi_heavy > 1]

bcr_data <- dplyr::filter(bcr_data, !cell_id %in% multi_heavy_cells)
cat(paste("There are", nrow(bcr_data),
          "rows in the data after filtering out cells with multiple heavy chains.\n"))

## There are 172495 rows in the data after filtering out cells with multiple heavy chains.

# split cells by heavy and light chains
heavy_cells <- dplyr::filter(bcr_data, locus == "IGH")$cell_id
light_cells <- dplyr::filter(bcr_data, locus == "IGK" | locus == "IGL")$cell_id
no_heavy_cells <- light_cells[which(!light_cells %in% heavy_cells)]

bcr_data <- dplyr::filter(bcr_data, !cell_id %in% no_heavy_cells)
cat(paste("There are", nrow(bcr_data), "rows in the data after filtering out
          cells without heavy chains."))

## There are 163406 rows in the data after filtering out
## cells without heavy chains.

### read in GEX
gex_db <- readRDS("bcells_integrated_processed_labeled.rds")

# Object summary
# `print` can be used to obtain a general overview of the Seurat object
# (number of features, number of samples, etc.)
print(gex_db)

# cell type annotations
head(Idents(gex_db), 1)

## DM001_AAACCAAAGCCACAAC-1 
## Naive #1 
## 10 Levels: Naive #1 ISG+ Naive Activated Naive Naive #2 ... Antibody Secreting Cells

Cells(gex_db)[10000]

# make cell IDs in BCR match those in Seurat object
bcr_data$cell_id_unique <- paste0(bcr_data$sample_id, "_", bcr_data$cell_id)

# first id in the BCR data
bcr_data$cell_id_unique[1]

## "DM001_AAACCAAAGCCACAAC-1"

# first id in the GEX data
Cells(gex_db)[1]

## "DM001_AAACCAAAGCCACAAC-1"

# match indices to find the position of the BCR cells in the GEX data
# different from finding the position of the GEX cells in the BCR data!
match.index <- match(bcr_data$cell_id_unique, Cells(gex_db))

mean(is.na(match.index))

## [1] 0.05463692

# add cell annotations to BCR data
Idents(gex_db) <- "cluster_label"
cell.annotation <- as.character(Idents(gex_db))

bcr_data$gex_annotation <-
  unlist(lapply(match.index, function(x) {
    ifelse(is.na(x), NA, cell.annotation[x])
  }))

bcr_data$gex_annotation[1:5]

# remove cells that didn't match
bcr_data <- dplyr::filter(bcr_data, !is.na(gex_annotation))

cat(paste("There are", nrow(bcr_data), "that match with GEX data"))

## There are 154478 that match with GEX data

### CLONAL ANALYSIS

dist_heavy <- distToNearest(bcr_data, 
                            cellIdColumn = "cell_id",
                            first = FALSE,
                            onlyHeavy = TRUE,
                            fields = "sample_id",
                            nproc = 2, 
                            progress = TRUE)

pdf("figures/immcantation/distance_bcr.pdf", height = 12)
ggplot(subset(dist_heavy, !is.na(dist_nearest)), 
       aes(x=dist_nearest)) +
  geom_histogram(color="white", binwidth=0.02) +
  geom_vline(xintercept=0.12, color="firebrick", linetype=2) +
  labs(x = "Grouped Hamming distance", y = "Count") + 
  facet_grid(sample_id ~ ., scales="free_y") + 
  theme_bw()
dev.off()

threshold_output <- shazam::findThreshold(dist_heavy$dist_nearest,
                                          method = "gmm", model = "gamma-norm",
                                          cutoff = "user", spc = 0.995)
threshold <- threshold_output@threshold
threshold

### DEFINE CLONAL GROUPS

# call clones using hierarchicalClones
results <- hierarchicalClones(bcr_data,
                              cell_id = "cell_id_unique",
                              threshold = threshold,
                              only_heavy = TRUE, split_light = FALSE,
                              summarize_clones = FALSE,
                              fields = "sample_id")

# calculate and plot the rank-abundance curve
abund <- estimateAbundance(dplyr::filter(results, locus == "IGH"),
                           group = "sample_id", nboot = 100)

# plot by sample_id
pdf("figures/immcantation/rankabundance_sample_bcr.pdf")
plot(abund, silent=T, col = safe_colorblind_palette) + facet_wrap("sample_id", ncol = 3)
dev.off()

# get clone sizes using dplyr functions
clone_sizes <- countClones(dplyr::filter(results, locus == "IGH"),
                           groups = "sample_id")

# plot cells per clone
pdf("figures/immcantation/cellsperclone_bcr.pdf")
ggplot(clone_sizes, aes(x = seq_count)) +
  geom_bar() +
  facet_wrap("sample_id", ncol = 3) +
  labs(x = "Sequences per clone") +
  theme_bw()
dev.off()

# calculate and plot the rank-abundance curve
div <- alphaDiversity(dplyr::filter(results, locus == "IGH"),
                      group = "sample_id", nboot = 100)

pdf("figures/immcantation/alphadiversity_bcr.pdf")
plot(div, silent = TRUE, col = safe_colorblind_palette) + facet_wrap("sample_id", ncol = 3)
dev.off()

# read in IMGT files
references <- readIMGT(dir = "/groups/db3700_gp/mr4527/germlines/imgt/human/vdj")

## [1] "Read in 1200 from 17 fasta files"

results <- results[!is.na(results$clone_id), ]

# reconstruct germlines
results <- createGermlines(results, references, fields = "sample_id",
                           nproc = 4)

# this is typically only done on heavy chains, but can also be done on light chains
results_heavy <- dplyr::filter(results, locus == "IGH")

# calculate SHM frequency in the V gene
data_mut <-
  shazam::observedMutations(results_heavy,
                            sequenceColumn = "sequence_alignment",
                            germlineColumn = "germline_alignment_d_mask",
                            regionDefinition = IMGT_V,
                            frequency = TRUE,
                            combine = TRUE,
                            nproc = 4)

safe_colorblind_palette <- c("#88CCEE", "#af97c8", "#5570b6", "#A020F0","#AA4499",
                             "#44AA99", "#661100","#24b34b",
                             "#6699CC", "#025a48" )

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

# Sort by that custom order
data_mut <- data_mut %>%
  mutate(sample_id = factor(sample_id, levels = sample_order)) %>%
  arrange(sample_id)

pdf("figures/immcantation/mu_freq_sample.pdf")
ggplot(data_mut, aes(y = mu_freq, x = sample_id, fill = sample_id)) +
  geom_boxplot(col = safe_colorblind_palette) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(), axis.title.y = element_text(size = 20), 
        axis.text.y = element_text(size = 16)) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  labs(x = "", y = "Mutation frequency", fill = "sample_id") + 
  scale_fill_manual(values = safe_colorblind_palette)
dev.off()

# plotting mu_freq by isotype
pdf("figures/immcantation/mu_freq_isotype.pdf", width = 16)
ggplot(data_mut, aes(y = mu_freq, x = sample_id, fill = sample_id)) +
  geom_boxplot() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  labs(x = "", y = "Mutation frequency", fill = "Sample") +
  theme(axis.text.x = element_blank()) + facet_grid( . ~ c_call) + 
  scale_fill_manual(values = safe_colorblind_palette)
dev.off()

## Supplemental S5
# plotting mu_freq by cell annotation
pdf("figures/immcantation/mu_freq_cluster.pdf", width = 16)
ggplot(data_mut, aes(y = mu_freq, x = sample_id, fill = sample_id)) +
  geom_boxplot() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  labs(x = "", y = "Mutation frequency", fill = "sample_id") +
  theme(axis.text.x = element_blank()) + facet_grid( . ~ gex_annotation) + 
  scale_fill_manual(values = safe_colorblind_palette)
dev.off()

write_csv(data_mut, "data_mut.csv")
