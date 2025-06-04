#!/usr/bin/env Rscript

library(Seurat)
library(Signac)
# library(EnsDb.Hsapiens.v86) For Human
library(EnsDb.Mmusculus.v79)
library(GenomeInfoDb)


args <- commandArgs(trailingOnly = TRUE)
data_dir <- file.path(args[1], "outs", "filtered_feature_bc_matrix")
output_file <- args[2]
sample_name <- args[3]

data <- Read10X(data.dir = data_dir)
rna_counts <- data$`Gene Expression`
atac_counts <- data$Peaks
seurat_obj <- CreateSeuratObject(counts = data$`Gene Expression`, project = sample_name)

# Labeling refence genome in seurat object - Defaults to mm10 if unable
reference_genome <- "mm10"

summary_filepath <- file.path(args[1], "outs", "summary.csv")
tryCatch({
  genome_name <- read.csv(summary_filepath, stringsAsFactors = FALSE)$Genome[1]
}, error = function(e) {
  message("summary.csv or reference genome not found — defaulting to 'mm10'")
})

seurat_obj@misc$reference_genome <- reference_genome

# Calculate percentage of RNA reads mapped to mitochondrial genes - indicator for stress/damaged/lysed cells for exclusion
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")

# Convert rownames of ATAC count matrix to Genomic Ranges (coordinates) object
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
# %in% creates logical vector of peaks from standard chromosomes (filtering non-standard chromosomes like chrM) to reduce noise
grange.use <- seqnames(grange.counts) %in% GenomeInfoDb::standardChromosomes(grange.counts)
# Filters nonstandard peaks not in grange.use
atac_counts <- atac_counts[as.vector(grange.use), ]
# Annotations
# annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86) # for human
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79) # for mouse
# Convert naming style from "1" to "chr1" format
seqlevelsStyle(annotations) <- 'UCSC'
# Set genome tag 
# genome(annotations) <- "hg38" # for human
genome(annotations) <- "mm10" # for mouse

# Find fragment file
frag.file <- file.path(args[1], "outs", "atac_fragments.tsv.gz")

# Generate ATAC before filtering
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = 'mm10',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )

# Store ATAC information in object
seurat_obj[["ATAC"]] <- chrom_assay

saveRDS(seurat_obj, file = output_file)