#!/usr/bin/env nextflow

process qc_processing {
    publishDir "results/${sample_name}/qc_seurat", mode: 'copy'
    container 'seurat_qc'
    cache 'lenient'

    input:
    tuple val(sample_name), file(seurat_object)
    path script_file
    val utils_dir


    output:
    tuple val("${sample_name}"), file("${sample_name}_qc.rds"), emit: outputs
    val("${sample_name}"),                                      emit: sample_names
    file("prefilter_vlnplot.png")
    file("postfilter_vlnplot.png")
    file("${sample_name}_qc_summary.csv")

    script:
    """
    Rscript $script_file $sample_name $seurat_object $utils_dir
    """
}