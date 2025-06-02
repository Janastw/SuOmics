#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { generate_seurat_object } from "./modules/seurat_object/generate_seurat_object.nf"
include { qc_processing } from "./modules/seurat_qc/qc_processing.nf"
include { wnn_generation } from "./modules/wnn/wnn_generation.nf"
include { aggregate_summaries } from "./modules/aggregate_summaries/aggregate_summaries.nf"

params.data_dir = "$baseDir/data"
params.results_dir = "$baseDir/results"

workflow {    
    // Create seurat object
    def samples_ch = Channel.fromPath("data/*", type: 'dir')
    def script_ch = Channel.value(file("scripts/generate_seurat_object.R"))
    def seurat_objects = generate_seurat_object(samples_ch, script_ch, params.data_dir)    
    
    // Perform QC on each sample
    def qc_script_ch = Channel.value(file("scripts/qc_processing.R"))
    def qc_outputs = qc_processing(seurat_objects, qc_script_ch, params.results_dir)

    // Perform UMAPS for ATAC, RNA, and WNN
    def wnn_script_ch = Channel.value(file("scripts/wnn_generation.R"))
    def wnn_outputs = wnn_generation(qc_outputs.outputs, wnn_script_ch, params.results_dir)

    // Create final summaries
    def aggregate_summaries_ch = Channel.value(file("scripts/aggregate_summaries.py"))
    aggregate_summaries(qc_outputs.sample_names.collect(), aggregate_summaries_ch, params.results_dir)

}
// Can add email notifications
workflow.onComplete {
    log.info(workflow.success
        ? "\nDone: View graphs in --> ${params.results_dir}/summary"
        : "\nPipeline failed. Check logs.")
}

