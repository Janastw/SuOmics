#!/usr/bin/env Rscript

library(Seurat)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
data_dir <- file.path(args[1])
output_file <- args[2]
sample_name <- args[3]

seurat_obj_path <- args[1]
seurat_obj <- readRDS(seurat_obj_path)

# Prefilter plot generation
png("prefilter_vlnplot.png")
VlnPlot(seurat_obj, features = c("nCount_ATAC", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE, pt.size = 0) + NoLegend()
dev.off()

# Apply filter
seurat_obj <- subset(
  x = seurat_obj,
  subset =
    nCount_ATAC < 1e5 &
    nCount_ATAC > 5000 &
    nCount_RNA < 40000 &
    nCount_RNA > 100 &
    percent.mt < 20
)

# Postfilter
png("postfilter_vlnplot.png")
VlnPlot(seurat_obj, features = c("nCount_ATAC", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE, pt.size = 0) + NoLegend()
dev.off()

saveRDS(seurat_obj, file = output_file)