#!/usr/bin/env Rscript

library(Seurat)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
data_dir <- file.path(args[1])
output_file <- args[2]
sample_name <- args[3]

seurat_obj_path <- args[1]
seurat_obj <- readRDS(seurat_obj_path)

n_cells_before <- ncol(seurat_obj)

thresholds <- list(
  nCount_ATAC_max = 10000,
  nCount_ATAC_min = 5000,
  nCount_RNA_max = 40000,
  nCount_RNA_min = 100,
  percent.mt_max = 20
)

# Prefilter plot generation
png("prefilter_vlnplot.png")
VlnPlot(seurat_obj, features = c("nCount_ATAC", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE, pt.size = 0) + NoLegend()
dev.off()

# Apply filter
seurat_obj <- subset(
  x = seurat_obj,
  subset =
    nCount_ATAC < thresholds$nCount_ATAC_max &
    nCount_ATAC > thresholds$nCount_ATAC_min &
    nCount_RNA < thresholds$nCount_RNA_max &
    nCount_RNA > thresholds$nCount_RNA_min &
    percent.mt < thresholds$percent.mt_max
)

n_cells_after <- ncol(seurat_obj)

# Postfilter
png("postfilter_vlnplot.png")
VlnPlot(seurat_obj, features = c("nCount_ATAC", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE, pt.size = 0) + NoLegend()
dev.off()

qc_summary <- data.frame(
  sample = sample_name,
  cells_preQC = n_cells_before,
  cells_postQC = n_cells_after,
  nCount_ATAC_min = thresholds$nCount_ATAC_min,
  nCount_ATAC_max = thresholds$nCount_ATAC_max,
  nCount_RNA_min = thresholds$nCount_RNA_min,
  nCount_RNA_max = thresholds$nCount_RNA_max,
  percent_mt_max = thresholds$percent.mt_max
)

write.csv(qc_summary, file = paste0(sample_name, "_qc_summary.csv"), row.names = FALSE)

saveRDS(seurat_obj, file = output_file)