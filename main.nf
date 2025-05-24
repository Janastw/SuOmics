#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// include { generate_seurat_object } from "./modules/seurat_qc/generate_seurat_object.nf"

process generate_seurat_object {
    publishDir "results/${sample_name}/seurat_object", mode: 'copy'

    input:
    path sample_name

    output:
    path "${sample_name}.txt"

    script:
    """
    echo \$PWD > ${sample_name}.txt
    """
}

workflow {
    
    def samples_ch = Channel.fromPath("data/*", type: 'dir')

    samples_ch | generate_seurat_object
}
