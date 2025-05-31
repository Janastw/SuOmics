#!/usr/bin/env nextflow



process qc_processing {
    publishDir "results/${sample_name}/qc_seurat", mode: 'copy'
    container 'seurat_qc'
    cache 'lenient'
    errorStrategy 'retry'

    input:
    tuple val(sample_name), file(seurat_object)
    path script_file
    val results_dir

    output:
    tuple val("${sample_name}"), file("${sample_name}.rds")
    file("prefilter_vlnplot.png")
    file("postfilter_vlnplot.png")
    file("${sample_name}_qc_summary.csv")

    script:
    /* 
     * PARAMETERS FOR SCRIPT.R -> results_dir, outfile_name, sample_name <- NOT NECESSARY
     */
    // name=\$(basename $sample_name .rds)
    """
    Rscript $script_file $results_dir/$sample_name/seurat_object/$seurat_object $seurat_object $sample_name
    """
}