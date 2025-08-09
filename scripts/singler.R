#!/usr/bin/env Rscript
# 
# Description:
#   This script takes a Seurat object and annotates clusters by
#   cell types using SingleR's mouse reference cell data and
#   log-normalized RNA.
#   It outputs a Seurat object with annotations.
#
# Usage:
#   Rscript singler.R sample_name seurat_obj_rds
#
# Arguments:
#   sample_name  - Name of sample
#   seurat_obj_rds - Path to the seurat object .rds of the sample
#
############### LOAD DEPENDENCIES ################

library(SingleR)
library(celldex)
library(ggplot2)
library(Seurat)

############### READ ARGUMENTS ##################

args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]
seurat_obj_path <- args[2]

############### CORE SCRIPT ##############

seurat_obj <- readRDS(seurat_obj_path)
# markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# SingleR doesn't use sctransformed data. It uses RNA log-normalized.
ref_general <- celldex::MouseRNAseqData()
# ref_general <- celldex::TabulaMurisData()
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- NormalizeData(seurat_obj)
# DO NOT SCALE DATA BEFORE ANNOTATING - SINGLER doesn't use that
# seurat_obj <- ScaleData(seurat_obj)
# Check overlap
# overlap <- length(intersect(rownames(seurat_obj), rownames(ref_general)))
# total_seurat <- length(rownames(seurat_obj))
# total_ref <- length(rownames(ref_general))
# Per cell annotation
singleR_annotations <- SingleR(test = GetAssayData(seurat_obj, slot = "data"), ref = ref_general, labels = ref_general$label.main, de.method = 'wilcox')
seurat_obj$SingleR_full <- singleR_annotations$labels
seurat_obj$celltype <- singleR_annotations$pruned.labels
# unique(singleR_annotations$pruned.labels)
# Per cluster annotation
# singleR_annotations <- SingleR(test = GetAssayData(seurat_obj, slot = "data"), ref = ref_general, labels = ref_general$label.main, clusters = seurat_obj$seurat_clusters)
# seurat_obj$celltype <- pred$labels[seurat_obj$seurat_clusters]
# seurat_obj <- AddMetaData(seurat_obj, singleR_annotations$pruned.labels, col.name = "SingleR_Annotations")
# seurat_obj <- SetIdent(seurat_obj, value = "SingleR_Annotations")
seurat_obj <- Seurat::FindVariableFeatures(seurat_obj)
seurat_obj <- Seurat::ScaleData(seurat_obj)
seurat_obj <- Seurat::RunPCA(seurat_obj, npcs = 30)
seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:30, reduction.name = "umap.rna", reduction.key = "rnaUMAP_")
base::saveRDS(seurat_obj, file = base::paste0(sample_name, "_singler.rds"))

################ PLOT GENERATION #####################

png("Annotation Score Heatmap.png")
plotScoreHeatmap(singleR_annotations)
dev.off()

png("Annotation Delta Distribution.png")
plotDeltaDistribution(singleR_annotations, ncol = 4, dots.on.top = FALSE)
dev.off()

png("UMAP.png")
DimPlot(seurat_obj,
        reduction = "umap.rna",
        group.by = "celltype",
        label = FALSE,
        repel = TRUE) + ggtitle("RNA") + labs( x = "UMAP1", y = "UMAP2")
dev.off()

#################### EOF #####################



