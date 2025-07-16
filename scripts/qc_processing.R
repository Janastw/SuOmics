#!/usr/bin/env Rscript

library(Seurat)
library(Signac)
library(ggplot2)
library(Matrix)
library(GenomeInfoDb)

args <- commandArgs(trailingOnly = TRUE)
data_dir <- base::file.path(args[1])
output_file <- args[2]
sample_name <- args[3]
utils_dir <- args[4]

seurat_obj_path <- args[1]
seurat_obj <- base::readRDS(seurat_obj_path)
n_cells_before <- base::ncol(seurat_obj)


thresholds <- list(
  nCount_ATAC_max = 10000,
  nCount_ATAC_min = 5000,
  nCount_RNA_max = 40000,
  nCount_RNA_min = 100,
  percent.mt_max = 20,
  pct_reads_in_peaks = 8,
  blacklist_fraction = 0.15,
  nucleosome_signal = 4,
  TSS.enrichment = 1.5
)
Seurat::DefaultAssay(seurat_obj) <- "ATAC"
seurat_obj$pct_reads_in_peaks <- seurat_obj$peak_region_fragments / seurat_obj$passed_filters * 100

blacklist_file <- list.files(path = utils_dir,
                            pattern = base::paste0("^", seurat_obj@misc$reference_genome, "-blacklist.*\\.bed\\.gz$"),
                            full.names = TRUE
                            )[1]

if (!base::is.na(blacklist_file)) {
  message("Using blacklist file: ", blacklist_file)
  library(rtracklayer)
  blacklist <- rtracklayer::import(blacklist_file)
  GenomeInfoDb::seqlevelsStyle(blacklist) <- "UCSC"

  blacklist_counts <- FeatureMatrix(
    fragments = Fragments(seurat_obj)[[1]],
    features = blacklist,
    cells = colnames(seurat_obj)
  )
  seurat_obj$blacklist_fraction <- Matrix::colSums(blacklist_counts) / seurat_obj$nCount_ATAC
} else {
  warning("No blacklist file found for reference genome:", seurat_obj@misc$reference_genome)
  thresholds$blacklist_fraction = 0
  seurat_obj$blacklist_fraction <- 1
}

# Generate needed qc metrics
seurat_obj <- Signac::TSSEnrichment(seurat_obj)
seurat_obj <- Signac::NucleosomeSignal(seurat_obj)

# Prefilter plot generation
png("prefilter_vlnplot.png")
VlnPlot(seurat_obj, features = c("nCount_ATAC", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE, pt.size = 0) + NoLegend()
dev.off()

# Apply filter
seurat_obj <- subset(
  x = seurat_obj,
  subset =
    nCount_ATAC < thresholds$nCount_ATAC_max &
    nCount_ATAC > thresholds$nCount_ATAC_min &
    nCount_RNA < thresholds$nCount_RNA_max &
    nCount_RNA > thresholds$nCount_RNA_min &
    percent.mt < thresholds$percent.mt_max &
    pct_reads_in_peaks > thresholds$pct_reads_in_peaks &
    blacklist_fraction < thresholds$blacklist_fraction #&
    # nucleosome_signal < thresholds$nucleosome_signal #&
    # TSS.enrichment < thresholds$TSS.enrichment
)

n_cells_after <- ncol(seurat_obj)

# Postfilter
png("postfilter_vlnplot.png")
VlnPlot(seurat_obj, features = c("nCount_ATAC", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE, pt.size = 0) + NoLegend()
dev.off()

qc_summary <- base::data.frame(
  sample = sample_name,
  cells_preQC = n_cells_before,
  cells_postQC = n_cells_after,
  nCount_ATAC_min = thresholds$nCount_ATAC_min,
  nCount_ATAC_max = thresholds$nCount_ATAC_max,
  nCount_RNA_min = thresholds$nCount_RNA_min,
  nCount_RNA_max = thresholds$nCount_RNA_max,
  percent_mt_max = thresholds$percent.mt_max
)

utils::write.csv(qc_summary, file = base::paste0(sample_name, "_qc_summary.csv"), row.names = FALSE)

base::saveRDS(seurat_obj, file = output_file)