#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { generate_seurat_object } from "./modules/seurat_object/generate_seurat_object.nf"
include { generate_archr_object } from "./modules/archr_object/generate_archr_object.nf"
include { qc_processing } from "./modules/seurat_qc/qc_processing.nf"
include { wnn_generation } from "./modules/wnn/wnn_generation.nf"
include { aggregate_summaries } from "./modules/reports/aggregate_summaries.nf"
include { dl_blacklists } from "./modules/utils/dl_blacklists.nf"
include { singler } from "./modules/annotations/singler.nf"
include { sctransform } from './modules/normalization/sctransform.nf'
include { tfidf_lsi } from './modules/normalization/tfidf_lsi.nf'
// include { rpca } from './modules/rpca/rpca.nf'

params.data_dir = "$baseDir/data"
params.results_dir = "$baseDir/results"
params.utils_dir = "$baseDir/utils"

workflow {    
    // Download utilities
    def blacklist_script = Channel.value(file("scripts/dl_blacklists.sh"))
    dl_blacklists(blacklist_script)

    // Create sample channels
    def samples_ch = Channel.fromPath("data/*", type: 'dir')

    // Create seurat object
    def gen_seu_obj_ch = Channel.value(file("scripts/generate_seurat_object.R"))
    def seurat_objects = generate_seurat_object(samples_ch, gen_seu_obj_ch)    

    // // Create ArchR object
    def gen_archr_obj_ch = Channel.value(file("scripts/generate_archr_object.R"))
    def archr_objects = generate_archr_object(samples_ch, gen_archr_obj_ch, params.data_dir)    
    
    // Perform QC on each sample
    def qc_script_ch = Channel.value(file("scripts/qc_processing.R"))
    def qc_outputs = qc_processing(seurat_objects, qc_script_ch, params.utils_dir)

    // scRNA-Seq Normalization and Clustering
    def sctransform_script_ch = Channel.value(file("scripts/sctransform.R"))
    def sctransform_outputs = sctransform(qc_outputs.outputs, sctransform_script_ch)

    // scATAC-Seq Normalization and Clustering
    def tfidf_lsi_script_ch = Channel.value(file("scripts/tfidf_lsi.R"))
    def tfidf_lsi_outputs = tfidf_lsi(sctransform_outputs.outputs, tfidf_lsi_script_ch)
    
    // Perform Cell Type Annotations with SingleR
    def singler_script_ch = Channel.value(file("scripts/singler.R"))
    def singler_outputs = singler(tfidf_lsi_outputs.outputs, singler_script_ch)

    // // SCARLink
    // // def scarlink_script_ch = Channel.value(file("scripts/scarlink.py"))
    // // def scarlink_outputs = scarlink(sctransform_outputs.outputs, archr_objects.outputs, scarlink_script_ch)

    // Perform UMAPS for ATAC, RNA, and WNN
    def wnn_script_ch = Channel.value(file("scripts/wnn_generation.R"))
    def wnn_outputs = wnn_generation(tfidf_lsi_outputs.outputs, wnn_script_ch, params.results_dir)

    // Create RPCA
    // def rpca_ch = Channel.value(file("scripts/rpca.R"))
    // rpca(qc_outputs.sample_names.collect()[0], qc_outputs.sample_names.collect()[1], rpca_ch)
    // def samples_by_tumor_type = qc_outputs.sample_names.collect()

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

