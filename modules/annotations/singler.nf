#!/usr/bin/env nextflow



process singler {
    publishDir "results/${sample_name}/annotations", mode: 'copy'
    container 'seurat_qc'
    cache 'lenient'

    input:
    tuple val(sample_name), file(seurat_object)
    path script_file
    val results_dir

    output:
    tuple val("${sample_name}"), file("${sample_name}_singler.rds"), emit: outputs
    val("${sample_name}"),                                   emit: sample_names
    file("Annotation Score Heatmap.png")
    file("Annotation Delta Distribution.png")
    file("UMAP.png")

    script:
    /* 
     * PARAMETERS FOR SCRIPT.R -> results_dir, outfile_name, sample_name <- NOT NECESSARY
     */
    // name=\$(basename $sample_name .rds)
    """
    Rscript $script_file $results_dir/$sample_name/qc_seurat/$seurat_object $seurat_object $sample_name
    """
}