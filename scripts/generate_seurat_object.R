#!/usr/bin/env Rscript
# 
# Description:
#   This script takes multimodal scRNAseq/scATACseq
#   (post-cellranger arc processed) and creates a Seurat object
#    Seurat object and performs SCTransform normalization.
#   It outputs a Seurat object normalized gene expression.
#   
#   Flexible for human or mouse data
#
# Usage:
#   Rscript generate_seurat_object.R sample_name
#
# Arguments:
#   sample_name  - Name of sample
#
############### LOAD DEPENDENCIES ################
library(Seurat)
library(Signac)
library(GenomeInfoDb)

############### READ ARGUMENTS ##################

args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]

############### CORE SCRIPT ##############

# Can create a fallback for reading in the regular data instead of relying on a .h5 file being present
sample_data_file <- base::file.path(sample_name, "outs", "filtered_feature_bc_matrix.h5")
data <- Seurat::Read10X_h5(filename = sample_data_file)
rna_counts <- data$`Gene Expression`
atac_counts <- data$Peaks
seurat_obj <- Seurat::CreateSeuratObject(counts = data$`Gene Expression`, project = sample_name)
# Labeling refence genome in seurat object - Defaults to mm10 if unable
reference_genome <- "mm10"
summary_filepath <- base::file.path(args[1], "outs", "summary.csv")
if (base::file.exists(summary_filepath)) {
  genome_val <- base::tolower(utils::read.csv(summary_filepath, stringsAsFactors = FALSE)$Genome[1])
  if (!is.na(genome_val) && genome_val != "") {
    if (genome_val == "grch38") {
      reference_genome <- "hg38"
    } else {
      reference_genome <- genome_val
    }
  }
} else {
  base::message("summary.csv or reference genome not found — defaulting to 'mm10'")
}
seurat_obj@misc$reference_genome <- reference_genome
annotations <- NULL
# Calculate percentage of RNA reads mapped to mitochondrial genes - indicator for stress/damaged/lysed cells for exclusion
if (reference_genome == "mm10") {
  library(EnsDb.Mmusculus.v79) # For Mouse
  annotations <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79) # for mouse
  GenomeInfoDb::seqlevelsStyle(annotations) <- 'UCSC'
  GenomeInfoDb::genome(annotations) <- "mm10" # for mouse
  seurat_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern = "^mt-")
} else if (reference_genome == "hg38") {
  library(EnsDb.Hsapiens.v86) # For Human
  annotations <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86) # for human
  GenomeInfoDb::seqlevelsStyle(annotations) <- 'UCSC'
  GenomeInfoDb::genome(annotations) <- "hg38" # for human
  seurat_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern = "^MT-")
} else {
  stop("No reference genome found for: ", reference_genome)
}
# TODO Address whether or not GRCH38 should stay as a reference genome later. Right now we are maintaining hg38
# } else if (reference_genome == "GRCh38") {
#   library(AnnotationHub)
#   library(ensembldb)
#   ah <- AnnotationHub()
#   annotations <- GetGRangesFromEnsDb(ensdb = ah[["AH104864"]]) # for human
#   cat("Genome string in annotation object:", unique(genome(annotations)), "-\n")
#   cat("Genome string being passed to CreateChromatinAssay:", reference_genome, "-\n")
#   # genome(annotations) <- "GRCh38" # for human
#   seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Convert rownames of ATAC count matrix to Genomic Ranges (coordinates) object
grange.counts <- Signac::StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
# %in% creates logical vector of peaks from standard chromosomes (filtering non-standard chromosomes like chrM) to reduce noise
grange.use <- GenomicRanges::seqnames(grange.counts) %in% GenomeInfoDb::standardChromosomes(grange.counts)
# Filters nonstandard peaks not in grange.use
atac_counts <- atac_counts[base::as.vector(grange.use), ]

# Find fragment file - USING RELATIVE FILE PATH
# frag.file <- base::file.path(args[1], "outs", "atac_fragments.tsv.gz")

library(fs)
frag.file.abs <- base::file.path(args[1], "outs", "atac_fragments.tsv.gz")
frag.file.rel <- fs::path_rel(frag.file.abs, start = getwd())

# Generate ATAC before filtering
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = base::unique(GenomeInfoDb::genome(annotations)),
   fragments = frag.file.rel,
   min.cells = 10,
   annotation = annotations
 )

# Store ATAC information in object
seurat_obj[["ATAC"]] <- chrom_assay
DefaultAssay(seurat_obj) <- "ATAC"

frag_counts <- Signac::CountFragments(fragments = frag.file.rel)
seurat_obj$passed_filters <- frag_counts$frequency_count[base::match(base::colnames(seurat_obj), frag_counts$CB)]
peak_counts <- Signac::FeatureMatrix(
  fragments = Signac::Fragments(seurat_obj),
  features = GenomicRanges::granges(seurat_obj[["ATAC"]], use.names = TRUE),
  cells = base::colnames(seurat_obj)
)
seurat_obj$peak_region_fragments <- Matrix::colSums(peak_counts)

base::saveRDS(seurat_obj, file = paste0(sample_name, ".rds"))

#################### EOF #####################
