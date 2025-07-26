#!/usr/bin/env Rscript


args <- commandArgs(trailingOnly = TRUE)
data_file <- base::file.path(args[1], "outs", "atac_fragments.tsv.gz")
sample_name <- args[2]

library(ArchR)

addArchRLocking(FALSE)
# addArchRThreads(threads = 1)
addArchRGenome("mm10")

ArrowFiles <- createArrowFiles(
    inputFiles = c(data_file),
    sampleNames = c(sample_name),
    minTSS = 4,        # Don't set this too high because you can always increase later
    minFrags = 1000, 
    addTileMat = TRUE,
    addGeneScoreMat = TRUE
)