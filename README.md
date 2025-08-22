# RNA-seq Analysis (MASTERS Trial)

This repository contains R scripts for RNA-seq analysis of skeletal muscle samples from the MASTERS trial (GEO: GSE157585). The analysis includes preprocessing, differential expression, visualization, and gene set enrichment analysis (GSEA).

## Folder Structure

MASTERS-RNAseq_analysis/
├── data/
│   └── GSE157585_Kulkarni_Peck_et_al_MASTERS_raw_counts.txt
├── 01_preprocessing.R
├── 02_deseq2_analysis.R
├── 03_visualization.R
├── 04_gsea_preparation.R
├── 05_gsea_analysis.R
├── 06_custom_pathways.R
└── .DS_Store

## How to Run
1. Open this repository in **RStudio**.  
2. Ensure the `data/` folder contains the counts file(s) needed for analysis.  
3. Run the scripts in order:  
   - `01_preprocessing.R` → `02_deseq2_analysis.R` → `03_visualization.R` → `04_gsea_preparation.R` → `05_gsea_analysis.R` → `06_custom_pathways.R`.  
4. Scripts use **relative paths**; no `setwd()` changes are required.  

## Data
- Counts file(s) included in the `data/` folder (can be a subset for testing).  
- Full dataset available from GEO: [GSE157585](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157585).  

## Notes
- Make sure required R packages (e.g., `DESeq2`, `ggplot2`, `fgsea`) are installed before running scripts.  
- For reproducibility, run scripts in order and do not modify the folder structure.
