#!/usr/bin/env Rscript
# 
# Description:
#   This script takes a Seurat object and performs tf-idf normalization.
#   It outputs a Seurat object with normalized peaks.
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
library(Signac)
library(EnhancedVolcano)

############### READ ARGUMENTS ##################

args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]
seurat_obj_filename <- args[2]

############### CORE SCRIPT ##############

seurat_obj <- base::readRDS(seurat_obj_filename)
Seurat::DefaultAssay(seurat_obj) <- "ATAC"
seurat_obj <- Signac::RunTFIDF(seurat_obj)
seurat_obj <- Signac::FindTopFeatures(seurat_obj, min.cutoff = 'q0')
seurat_obj <- Signac::RunSVD(seurat_obj)
seurat_obj <- Seurat::FindNeighbors(seurat_obj, reduction = "lsi", dims = 2:30)
seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- Seurat::RunUMAP(seurat_obj, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
seurat_obj$celltype <- Seurat::Idents(seurat_obj)
base::saveRDS(seurat_obj, file = paste0(sample_name, "_tfidf.rds"))

################ PLOT GENERATION #####################

grDevices::png("depthCorr.png")
Signac::DepthCor(seurat_obj)
grDevices::dev.off()

png("UMAP.png")
Seurat::DimPlot(
  seurat_obj,
  reduction = "umap.atac",
  group.by = "celltype",
  label = TRUE,
  label.size = 2.5,
  repel = TRUE
) + ggplot2::ggtitle("ATAC")
dev.off()

atac_markers <- Seurat::FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25
)

png("volcano.png")
EnhancedVolcano::EnhancedVolcano(
  atac_markers,
  lab = rownames(atac_markers),
  x = 'avg_log2FC',
  y = 'p_val_adj',
  pCutoff = 0.05,
  FCcutoff = 0.5,
  selectLab = head(rownames(atac_markers), 15), # Label only top 15
  pointSize = 1.5,         # Small points
  labSize = 3.0,           # Small label
  max.overlaps = 20,       # Limiting labels
  title = "Differentially Accessible Peaks",
  subtitle = "Significant Peaks (ATAC-seq)",
  caption = paste0("total = ", nrow(atac_markers), " peaks")
)
dev.off()

# library(dplyr)

# top10_atac <- atac_markers %>%
#   group_by(cluster) %>%
#   top_n(n = 10, wt = avg_log2FC)

# png("diffAcc.png")
# DoHeatmap(
#   seurat_obj,
#   features = top10_atac$gene
# )
# dev.off()

#################### EOF #####################