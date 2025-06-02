#!/usr/bin/env nextflow



process wnn_generation {
    publishDir "results/${sample_name}/wnn", mode: 'copy'
    container 'seurat_qc'
    cache 'lenient'

    input:
    tuple val(sample_name), file(seurat_object)
    path script_file
    val results_dir

    output:
    tuple val("${sample_name}"), file("${sample_name}.rds"), emit: outputs
    file("combined_umap_plot.png")

    script:
    /* 
     * PARAMETERS FOR SCRIPT.R -> results_dir, outfile_name, sample_name <- NOT NECESSARY
     */
    """
    Rscript $script_file $results_dir/$sample_name/qc_seurat/$seurat_object $seurat_object $sample_name
    """
}