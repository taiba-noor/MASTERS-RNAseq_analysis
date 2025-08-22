############################################
# 06_custom_pathways.R
# Custom GSEA (example: mTORC1 signaling)
############################################

# Load libraries
library(msigdbr)
library(fgsea)
library(ggplot2)

# Load ranked genes
gene_ranks <- readRDS("gene_ranks_metformin_vs_placebo.rds")

# Get mTORC1 hallmark gene set
mtorc1_geneset <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::filter(gs_name == "HALLMARK_MTORC1_SIGNALING") %>%
  dplyr::pull(gene_symbol)

# Wrap into a list for fgsea
mtorc1_pathway <- list(HALLMARK_MTORC1_SIGNALING = mtorc1_geneset)

# Run fgsea
fgsea_mtorc1 <- fgsea(
  pathways = mtorc1_pathway,
  stats = gene_ranks,
  minSize = 15,
  maxSize = 500,
  nperm = 10000
)

# Inspect results
fgsea_mtorc1

# Plot enrichment curve
plotEnrichment(
  mtorc1_pathway$HALLMARK_MTORC1_SIGNALING,
  gene_ranks
) + ggtitle("Hallmark: mTORC1 Signaling")
