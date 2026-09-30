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

## merged tsv, see R script "tsv_merge.R"
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

## filter for just heavy chains
bcr_data <- dplyr::filter(bcr_data, locus == "IGH")

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

# match indices
match.index <- match(Cells(gex_db), bcr_data$cell_id_unique)

mean(is.na(match.index))

## [1] 0.1982484

# add BCR heavy chain V gene to GEX metadata
v_call <- as.character(bcr_data$v_call)

gex_db@meta.data$v_call <-
  unlist(lapply(match.index, function(x) {
    ifelse(is.na(x), NA, v_call[x])
  }))

### Figure 5D - VH434 UMAP
meta_data = gex_db@meta.data

for (i in 1:nrow(meta_data)){
  if (is.na(meta_data[i, "v_call"])) {
    meta_data[i, "IGHV4-34"] <- 0
  }
  else if (meta_data[i, "v_call"] == "IGHV4-34*01"){
    meta_data[i, "IGHV4-34"] <- 1}
  else {
    meta_data[i, "IGHV4-34"] <- 0}
}

gex_db <- AddMetaData(gex_db, meta_data)

pdf("figures/UMAP_VH434.pdf", width = 10, height = 5)
FeaturePlot(gex_db, features = "IGHV4-34", split.by = "Group", order = TRUE) 
dev.off()
