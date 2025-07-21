#!/usr/bin/env Rscript

list.files()
args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]

output_file <- sample_name
seurat_obj_filename <- base::file.path(base::paste0(sample_name, ".rds"))

library(Seurat)
library(Signac)
seurat_obj <- base::readRDS(seurat_obj_filename)

Seurat::DefaultAssay(seurat_obj) <- "ATAC"
seurat_obj <- Signac::RunTFIDF(seurat_obj)
seurat_obj <- Signac::FindTopFeatures(seurat_obj, min.cutoff = 'q0')
seurat_obj <- Signac::RunSVD(seurat_obj)

grDevices::png("depthCorr.png")
Signac::DepthCor(seurat_obj)
grDevices::dev.off()

seurat_obj <- Seurat::FindNeighbors(seurat_obj, reduction = "lsi", dims = 2:30)
seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- Seurat::RunUMAP(seurat_obj, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

seurat_obj$celltype <- Seurat::Idents(seurat_obj)

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

library(EnhancedVolcano)

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

base::saveRDS(seurat_obj, file = output_file)