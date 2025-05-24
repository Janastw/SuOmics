#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { generate_seurat_object } from "./modules/seurat_qc/generate_seurat_object.nf"

workflow {    
    def samples_ch = Channel.fromPath("data/*", type: 'dir')

    samples_ch | generate_seurat_object
}
