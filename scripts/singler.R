#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
seurat_obj_path <- args[1]
output_file <- args[2]
sample_name <- args[3]

seurat_obj <- readRDS(output_file)

library(SingleR)
library(celldex)
library(ggplot2)
library(Seurat)

# markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

seurat_obj <- SCTransform(seurat_obj, verbose = FALSE) |> RunPCA() |> RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
ref_general <- celldex::MouseRNAseqData()
# seurat_obj <- NormalizeData(seurat_obj)

# Check overlap
# overlap <- length(intersect(rownames(seurat_obj), rownames(ref_general)))
# total_seurat <- length(rownames(seurat_obj))
# total_ref <- length(rownames(ref_general))

DefaultAssay(seurat_obj) <- "SCT"

# Per cell annotation
singleR_annotations <- SingleR(test = GetAssayData(seurat_obj, slot = "data"), ref = ref_general, labels = ref_general$label.main, de.method = 'wilcox')
seurat_obj$celltype <- singleR_annotations$labels

# Per cluster annotation
# singleR_annotations <- SingleR(test = GetAssayData(seurat_obj, slot = "data"), ref = ref_general, labels = ref_general$label.main, clusters = seurat_obj$seurat_clusters)
# seurat_obj$celltype <- pred$labels[seurat_obj$seurat_clusters]

unique(singleR_annotations$pruned.labels)
png("Annotation Score Heatmap.png")
plotScoreHeatmap(singleR_annotations)
dev.off()

png("Annotation Delta Distribution.png")
plotDeltaDistribution(singleR_annotations, ncol = 4, dots.on.top = FALSE)
dev.off()


# seurat_obj <- AddMetaData(seurat_obj, singleR_annotations$pruned.labels, col.name = "SingleR_Annotations")
# seurat_obj <- SetIdent(seurat_obj, value = "SingleR_Annotations")

png("UMAP.png")
DimPlot(seurat_obj,
        reduction = "umap.rna",
        group.by = "celltype",
        label = TRUE,
        label.size = 2.5,
        repel = TRUE) + ggtitle("RNA")
dev.off()

base::saveRDS(seurat_obj, file = paste0(sample_name, "_singler.rds"))



