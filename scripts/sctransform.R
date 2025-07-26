#!/bin/usr/env Rscript

# Inputs needed: sample name, sample filepath
list.files()
args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]
seurat_obj_filename <- args[2]

# seurat_obj_filename <- base::file.path(paste0(sample_name, ".rds"))

library(Seurat)
seurat_obj <- base::readRDS(seurat_obj_filename)

seurat_obj <- Seurat::SCTransform(seurat_obj, verbose = FALSE)
seurat_obj <- Seurat::RunPCA(seurat_obj, verbose = FALSE)

png("ElbowPlot.png")
Seurat::ElbowPlot(seurat_obj)
dev.off()

seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = 0.5)

# Creates table to check if clusters exist
# table(seurat_obj$seurat_clusters)

seurat_obj <- Seurat::RunUMAP(
    seurat_obj,
    dims = 1:30,
    reduction.name = "umap.sct"
)

png("UMAP.png")
Seurat::DimPlot(seurat_obj, reduction = "umap.sct", label = TRUE, repel = TRUE)
dev.off()

library(dplyr)

markers_celltypes <- FindAllMarkers(
    seurat_obj,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
)
top5 <- markers_celltypes %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC)

png("DGEA.png")
DoHeatmap(
    seurat_obj,
    features = top5$gene
) + NoLegend()
dev.off()

base::saveRDS(seurat_obj, file = paste0(sample_name, "_sct.rds"))