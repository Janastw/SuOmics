#!/usr/bin/env Rscript

library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
data_dir <- file.path(args[1], "outs", "filtered_feature_bc_matrix")
output_file <- args[2]
sample_name <- args[3]

data <- Read10X(data.dir = data_dir)
seurat_obj <- CreateSeuratObject(counts = data$`Gene Expression`, project = sample_name)
saveRDS(seurat_obj, file = output_file)