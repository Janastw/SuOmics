#!/usr/bin/env nextflow



process generate_seurat_object {
    publishDir "results/${sample_name}/seurat_object", mode: 'copy'
    container 'seurat_qc'
    cache 'lenient'
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'terminate'}

    input:
    path sample_name
    path script_file

    output:
    tuple val("${sample_name}"), file("${sample_name}.rds"), emit: outputs

    script:
    """
    Rscript $script_file $sample_name
    """
}