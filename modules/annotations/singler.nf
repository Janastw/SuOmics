#!/usr/bin/env nextflow



process singler {
    publishDir "results/${sample_name}/annotations", mode: 'copy'
    container 'seurat_qc'
    cache 'lenient'

    input:
    tuple val(sample_name), file(seurat_object)
    path script_file

    output:
    tuple val("${sample_name}"), file("${sample_name}_singler.rds"), emit: outputs
    val("${sample_name}"),                                   emit: sample_names
    file("Annotation Score Heatmap.png")
    file("Annotation Delta Distribution.png")
    file("UMAP.png")

    script:
    """
    Rscript $script_file $sample_name $seurat_object
    """
}