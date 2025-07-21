#!/usr/bin/env nextflow

process sctransform {
    publishDir "results/${sample_name}/sctransform", mode: 'copy'
    container 'seurat_qc'
    cache 'lenient'

    input:
    tuple val(sample_name), file(seurat_obj_rds)
    path script_file

    output:
    tuple val("${sample_name}"), file("${sample_name}.rds"), emit: outputs
    file("ElbowPlot.png")
    file("UMAP.png")
    file("DGEA.png")


    script:
    """
    Rscript $script_file $sample_name
    """
}