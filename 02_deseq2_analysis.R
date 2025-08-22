# 02_deseq2_analysis.R
# Differential expression analysis

library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Load preprocessed DESeq2 object
dds <- readRDS("dds_preprocessed.rds")

# Clean IDs
clean_ids <- sub("\\..*", "", rownames(dds))
clean_ids <- sub("_.*", "", clean_ids)
gene_symbols <- sub(".*_", "", rownames(dds))
rownames(dds) <- clean_ids

# Annotation helper
annotate_res <- function(res, gene_symbols) {
  res_df <- as.data.frame(res)
  res_df$EnsemblID <- rownames(res_df)
  res_df$Symbol <- gene_symbols[match(rownames(res_df), clean_ids)]
  res_df$EntrezID <- mapIds(
    org.Hs.eg.db,
    keys = res_df$EnsemblID,
    column = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  return(res_df)
}

# Define contrasts
contrasts <- list(
  placebo = c("condition", "PRT_Placebo", "Baseline_Placebo"),
  metformin = c("condition", "PRT_Metformin", "Baseline_Metformin"),
  post = c("condition", "PRT_Metformin", "PRT_Placebo")
)

# Run contrasts
res_placebo <- annotate_res(results(dds, contrast = contrasts$placebo), gene_symbols)
res_metformin <- annotate_res(results(dds, contrast = contrasts$metformin), gene_symbols)
res_post <- annotate_res(results(dds, contrast = contrasts$post), gene_symbols)

# Save results
saveRDS(res_placebo, "res_annot_placebo_pre_vs_post.rds")
write.csv(res_placebo, "res_annot_placebo_pre_vs_post.csv", row.names = FALSE)

saveRDS(res_metformin, "res_annot_metformin_pre_vs_post.rds")
write.csv(res_metformin, "res_annot_metformin_pre_vs_post.csv", row.names = FALSE)

saveRDS(res_post, "res_annot_post_metformin_vs_placebo.rds")
write.csv(res_post, "res_annot_post_metformin_vs_placebo.csv", row.names = FALSE)
