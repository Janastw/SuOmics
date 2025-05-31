#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { generate_seurat_object } from "./modules/seurat_object/generate_seurat_object.nf"
include { qc_processing } from "./modules/seurat_qc/qc_processing.nf"
include { wnn_generation } from "./modules/wnn/wnn_generation.nf"

params.data_dir = "$baseDir/data"
params.results_dir = "$baseDir/results"

workflow {    
    def samples_ch = Channel.fromPath("data/*", type: 'dir')
    def script_ch = Channel.value(file("scripts/generate_seurat_object.R"))

    def seurat_objects = generate_seurat_object(samples_ch, script_ch, params.data_dir)    
    
    def qc_script_ch = Channel.value(file("scripts/qc_processing.R"))

    def qc_outputs = qc_processing(seurat_objects, qc_script_ch, params.results_dir)

    def wnn_script_ch = Channel.value(file("scripts/wnn_generation.R"))
    wnn_generation(qc_outputs, wnn_script_ch, params.results_dir)
}
