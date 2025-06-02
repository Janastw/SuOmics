#!/usr/bin/env Rscript

library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
data_dir <- file.path(args[1])
output_file <- args[2]
sample_name <- args[3]

seurat_obj_path <- args[1]
print(seurat_obj_path)
seurat_obj <- readRDS(seurat_obj_path)

library(Signac)
# library(EnsDb.Hsapiens.v86) For Human
library(EnsDb.Mmusculus.v79)
library(dplyr)
library(ggplot2)

# RNA Analysis
DefaultAssay(seurat_obj) <- "RNA"
# Normalization, dimensionality reduction, and vizualtion of RNA assay
seurat_obj <- SCTransform(seurat_obj, verbose = FALSE) %>%
    RunPCA() %>%
    RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')


DefaultAssay(seurat_obj) <- "ATAC"
# TF-IDF - Normalize Peak Accessibility for Downstream Analysis (SVD/UMAP)
seurat_obj <- RunTFIDF(seurat_obj)

# q0 keeps all features - Adjust to q25 for retaining top 25% or any other integer
seurat_obj <- FindTopFeatures(seurat_obj, min.cutoff = 'q0')
# SVD -> PCA for sparse data (treated as Latent Semantic Indexing (LSI))
seurat_obj <- RunSVD(seurat_obj)
# UMAP on LSI dims 2 to 50 ( 1st dim associated with sequencing depth/artifacts)
seurat_obj <- RunUMAP(seurat_obj, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# Start of WNN for RNA and ATAC modes
# Creates weighted scheme from neighbors in each space of PCA (RNA) and LSI (ATAC)
seurat_obj <- FindMultiModalNeighbors(seurat_obj, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
# UMAP - Create joint embedding
seurat_obj <- RunUMAP(seurat_obj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
# Use Leiden clustering on WNN for clustering
seurat_obj <- FindClusters(seurat_obj, graph.name = "wsnn", algorithm = 3, verbose = FALSE)


# Find subclustering done on cluster 6 here due to biologically distinct subtypes of CD8 T cells
seurat_obj <- FindSubCluster(seurat_obj, cluster = 6, graph.name = "wsnn", algorithm = 3)
# Assign and rename identities for immune cells
Idents(seurat_obj) <- "sub.cluster"
seurat_obj <- RenameIdents(seurat_obj,
  "0" = "CD14 Mono",
  "1" = "CD4 Naive",
  "2" = "CD8 Naive",
  "3" = "CD4 TCM",
  "4" = "CD4 TEM",
  "5" = "CD16 Mono",
  "6_0" = "CD8 TEM_1",
  "6_1" = "CD8 TEM_2",
  "7" = "NK",
  "8" = "Naive B",
  "9" = "CD14 Mono",
  "10" = "Intermediate B",
  "11" = "Memory B",
  "12" = "MAIT",
  "13" = "Treg"
)
seurat_obj$celltype <- Idents(seurat_obj)

p1 <- DimPlot(
  seurat_obj,
  reduction = "umap.rna",
  group.by = "celltype",
  label = TRUE,
  label.size = 2.5,
  repel = TRUE
) + ggtitle("RNA")

p2 <- DimPlot(
  seurat_obj,
  reduction = "umap.atac",
  group.by = "celltype",
  label = TRUE,
  label.size = 2.5,
  repel = TRUE
) + ggtitle("ATAC")

p3 <- DimPlot(
  seurat_obj,
  reduction = "wnn.umap",
  group.by = "celltype",
  label = TRUE,
  label.size = 2.5,
  repel = TRUE
) + ggtitle("WNN")

combined_plot <- (p1 + p2 + p3) &
  NoLegend() &
  theme(plot.title = element_text(hjust = 0.5))

ggsave("combined_umap_plot.png", combined_plot, width = 18, height = 6, dpi = 300)

saveRDS(seurat_obj, file = output_file)