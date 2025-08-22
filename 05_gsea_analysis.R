############################################
# 05_gsea_analysis.R
# Run GSEA on Hallmark pathways
############################################

# Load libraries
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# Load ranked genes
gene_ranks <- readRDS("gene_ranks_metformin_vs_placebo.rds")

# Get Hallmark gene sets
hallmark_df <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_terms <- hallmark_df[, c("gs_name", "gene_symbol")]

# Run GSEA
gsea_hallmark <- GSEA(
  geneList = gene_ranks,
  TERM2GENE = hallmark_terms,
  pvalueCutoff = 0.1,
  verbose = FALSE
)

# Inspect results
head(gsea_hallmark@result, 10)

# Plots
dotplot(gsea_hallmark, showCategory = 10, split = ".sign") +
  facet_grid(. ~ .sign) +
  ggtitle("Top 10 Hallmark Pathways (GSEA)")

ridgeplot(gsea_hallmark) + 
  ggtitle("Hallmark Pathways - Ridgeplot")

# Enrichment plot for Oxidative Phosphorylation
gseaplot2(
  gsea_hallmark,
  geneSetID = "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  pvalue_table = TRUE,
  title = "Hallmark: Oxidative Phosphorylation"
)

# Heatmap of top 5 pathways
gseaplot2(gsea_hallmark, geneSetID = 1:5, pvalue_table = TRUE)
