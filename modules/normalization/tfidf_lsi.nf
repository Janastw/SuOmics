#!/usr/bin/env nextflow

process tfidf_lsi {
    publishDir "results/${sample_name}/tfidf_lsi", mode: 'copy'
    container 'seurat_qc'
    cache 'lenient'

    input:
    tuple val(sample_name), file(seurat_obj_rds)
    path script_file

    output:
    tuple val("${sample_name}"), file("${sample_name}_tfidf.rds"), emit: outputs
    file("depthCorr.png")
    file("volcano.png")
    file("UMAP.png")

    script:
    """
    Rscript $script_file $sample_name $seurat_obj_rds
    """
}