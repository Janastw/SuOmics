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

library(ArchR)

############### READ ARGUMENTS ##################

args <- commandArgs(trailingOnly = TRUE)
data_file <- base::file.path(args[1], "outs", "atac_fragments.tsv.gz")
sample_name <- args[2]

############### CORE SCRIPT ##############

addArchRLocking(FALSE)
# addArchRThreads(threads = 1)
addArchRGenome("mm10")

ArrowFiles <- createArrowFiles(
    inputFiles = c(data_file),
    sampleNames = c(sample_name),
    minTSS = 4,
    minFrags = 1000, 
    addTileMat = TRUE,
    addGeneScoreMat = TRUE
)