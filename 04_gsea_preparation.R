############################################
# 04_gsea_preparation.R
# Prepare ranked gene list for GSEA
############################################

# Load libraries
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)

# Load differential expression results
res_table <- read.csv("res_annot_post_metformin_vs_placebo.csv")

# Map Ensembl IDs to gene symbols
res_table$symbol <- mapIds(
  org.Hs.eg.db,
  keys = res_table$ensembl_id,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Keep one entry per gene symbol (max |stat|)
res_table_unique <- res_table %>%
  filter(!is.na(symbol) & !is.na(stat)) %>%
  group_by(symbol) %>%
  slice_max(abs(stat), n = 1) %>%
  ungroup()

# Create ranked vector for GSEA
gene_ranks <- res_table_unique$stat
names(gene_ranks) <- res_table_unique$symbol

# Add tiny noise to break ties
set.seed(123) # for reproducibility
gene_ranks <- gene_ranks + rnorm(length(gene_ranks), mean = 0, sd = 1e-6)

# Sort decreasingly
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

# Quick checks
head(gene_ranks)
length(gene_ranks)
anyDuplicated(names(gene_ranks))  # should be 0

# Save gene ranks for later scripts
saveRDS(gene_ranks, file = "gene_ranks_metformin_vs_placebo.rds")
