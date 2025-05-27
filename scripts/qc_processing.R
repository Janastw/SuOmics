#!/usr/bin/env Rscript

library(Seurat)


args <- commandArgs(trailingOnly = TRUE)
data_dir <- file.path(args[1])
output_file <- args[2]
sample_name <- args[3]

seurat_obj_path <- args[1]
seurat_obj <- readRDS(seurat_obj_path)

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")

png("prefilter_vlnplot.png", width = 1200, height = 800)
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
dev.off()

seurat_obj <- subset(seurat_obj,
  subset = nCount_RNA > 100 & nCount_RNA < 40000 & percent.mt < 20
)

png("postfilter_vlnplot.png", width = 1200, height = 800)
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
dev.off()

saveRDS(seurat_obj, file = output_file)