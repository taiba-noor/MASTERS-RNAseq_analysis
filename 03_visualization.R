# 03_visualization.R
# PCA plots, MA plots, Volcano plots, UpSet analysis

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(UpSetR)

# Load preprocessed object and DESeq2 results
dds <- readRDS("dds_preprocessed.rds")
res_placebo <- readRDS("res_annot_placebo_pre_vs_post.rds")
res_metformin <- readRDS("res_annot_metformin_pre_vs_post.rds")
res_post <- readRDS("res_annot_post_metformin_vs_placebo.rds")

# Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)
sample_info <- as.data.frame(colData(dds))

### --- PCA PLOTS --- ###

plot_pca_subset <- function(vsd, sample_info, groups, title, label_ids = FALSE) {
  keep <- rownames(sample_info) %in% colnames(vsd)[sample_info$condition %in% groups]
  vsd_sub <- vsd[, keep]
  info_sub <- sample_info[colnames(vsd_sub), ]
  
  pca <- prcomp(t(assay(vsd_sub)))
  pca_df <- as.data.frame(pca$x[, 1:2])
  pca_df$condition <- factor(info_sub$condition, levels = groups)
  pca_df$participant <- info_sub$participant
  
  var_expl <- round(summary(pca)$importance[2, 1:2] * 100, 1)
  
  custom_cols <- c(
    "Baseline_Placebo"   = "#1f78b4",
    "PRT_Placebo"        = "#a6cee3",
    "Baseline_Metformin" = "#33a02c",
    "PRT_Metformin"      = "#b2df8a"
  )
  
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 4, alpha = 0.9) +
    scale_color_manual(values = custom_cols, drop = TRUE) +
    labs(
      title = title,
      x = paste0("PC1 (", var_expl[1], "%)"),
      y = paste0("PC2 (", var_expl[2], "%)")
    ) +
    theme_minimal(base_size = 15) +
    theme(
      panel.grid.major = element_line(color = "gray85"),
      panel.grid.minor = element_line(color = "gray95"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
  
  if (label_ids) {
    p <- p + geom_text_repel(aes(label = participant), size = 3, max.overlaps = 20)
  }
  
  return(p)
}

# Example usage
plot_pca_subset(vsd, sample_info, c("Baseline_Metformin", "PRT_Metformin"),
                "PCA: Pre vs Post Metformin", label_ids = TRUE)


### --- MA PLOTS --- ###

plotMA(results(dds, contrast = c("condition", "PRT_Metformin", "Baseline_Metformin")),
       ylim = c(-3, 3), main = "MA Plot: Pre vs Post Metformin")

plotMA(results(dds, contrast = c("condition", "PRT_Placebo", "Baseline_Placebo")),
       ylim = c(-3, 3), main = "MA Plot: Pre vs Post Placebo")

plotMA(results(dds, contrast = c("condition", "PRT_Metformin", "PRT_Placebo")),
       ylim = c(-3, 3), main = "MA Plot: Post Metformin vs Post Placebo")


### --- VOLCANO PLOTS --- ###

plot_volcano <- function(res_df, title, top_labels = FALSE, top_n = 20) {
  res_df$negLog10Padj <- -log10(res_df$padj)
  res_df$category <- "NS"
  res_df$category[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Up"
  res_df$category[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down"
  
  p <- ggplot(res_df, aes(x = log2FoldChange, y = negLog10Padj, color = category)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = c("NS" = "grey70", "Up" = "darkgreen", "Down" = "darkorange")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
    coord_cartesian(xlim = c(-4, 4)) +
    theme_minimal() +
    labs(title = title, x = "Log2 Fold Change", y = "-Log10 Adjusted p-value")
  
  if (top_labels) {
    top_genes <- res_df[res_df$category != "NS", ]
    top_genes <- top_genes[order(top_genes$padj), ]
    top_genes <- head(top_genes, top_n)
    p <- p + geom_text_repel(data = top_genes, aes(label = rownames(top_genes)),
                             size = 3, max.overlaps = 15)
  }
  
  return(p)
}

# Example volcano
plot_volcano(as.data.frame(res_metformin), "Volcano: Pre vs Post Metformin", top_labels = TRUE)


### --- UPSET PLOT --- ###

deg_met <- rownames(res_metformin)[which(res_metformin$padj < 0.05)]
deg_placebo <- rownames(res_placebo)[which(res_placebo$padj < 0.05)]

deg_list <- list("Metformin PRT" = deg_met, "Placebo PRT" = deg_placebo)

upset(
  fromList(deg_list),
  order.by = "freq",
  keep.order = TRUE,
  main.bar.color = c("#f5a7a6", "#83c0df", "#a788b5"),
  sets.bar.color = c("#f5a7a6", "#83c0df"),
  matrix.color = "grey40",
  mainbar.y.label = "Number of DEGs",
  sets.x.label = "DEGs per condition",
  text.scale = c(1.5, 1.5, 1, 1, 1.5, 1.2)
)
