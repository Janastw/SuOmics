#!/bin/usr/env Rscript
# 
# Description:
#   This script takes a Seurat object and performs SCTransform normalization.
#   It outputs a Seurat object normalized gene expression.
#
# Usage:
#   Rscript normalize_data.R sample_name seurat_obj_rds
#
# Arguments:
#   sample_name  - Name of sample
#   seurat_obj_rds - Path to the seurat object .rds of the sample
#
############### LOAD DEPENDENCIES ################

library(Seurat)
library(dplyr)

############### READ ARGUMENTS ##################

args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]
seurat_obj_filename <- args[2]

############### CORE SCRIPT ##############

seurat_obj <- base::readRDS(seurat_obj_filename)
seurat_obj <- Seurat::SCTransform(seurat_obj, verbose = FALSE)
seurat_obj <- Seurat::RunPCA(seurat_obj, verbose = FALSE)
seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- Seurat::RunUMAP(
    seurat_obj,
    dims = 1:30,
    reduction.name = "umap.sct"
)
base::saveRDS(seurat_obj, file = paste0(sample_name, "_sct.rds"))

################ PLOT GENERATION #####################

# Creates table to check if clusters exist
# table(seurat_obj$seurat_clusters)

png("ElbowPlot.png")
Seurat::ElbowPlot(seurat_obj)
dev.off()

png("UMAP.png")
Seurat::DimPlot(seurat_obj, reduction = "umap.sct", label = TRUE, repel = TRUE)
dev.off()

markers_celltypes <- Seurat::FindAllMarkers(
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


#################### EOF #####################