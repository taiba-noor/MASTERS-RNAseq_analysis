
# 01_preprocessing.R
# Load raw counts and create DESeq2 dataset

library(DESeq2)

# Load counts (download separately from GEO: GSE157585)
counts <- read.delim("data/GSE157585_Kulkarni_Peck_et_al_MASTERS_raw_counts.txt",
                     row.names = 1, check.names = FALSE)

# Build sample metadata
samples <- colnames(counts)
coldata <- data.frame(
  sample = samples,
  group = ifelse(grepl("^P", samples), "Placebo", "Metformin"),
  participant = sub("^[PM]_0*([0-9]+)_.*", "\\1", samples),
  time = ifelse(grepl("_1$", samples), "Baseline", "PRT"),
  batch = "batch1", # modify if multiple batches are known
  stringsAsFactors = FALSE
)

rownames(coldata) <- coldata$sample
coldata$participant <- factor(coldata$participant)
coldata$condition <- factor(
  paste0(coldata$time, "_", coldata$group),
  levels = c("Baseline_Placebo","PRT_Placebo",
             "Baseline_Metformin","PRT_Metformin")
)

# Construct DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ participant + condition
)

# Filter low counts
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]

# Run DESeq
dds <- DESeq(dds)

# Save preprocessed object
saveRDS(dds, "dds_preprocessed.rds")