#!/usr/bin/env nextflow



process generate_seurat_object {
    publishDir "results/${sample_name}/seurat_object", mode: 'copy'
    container 'seurat_qc'
    cache 'lenient'

    input:
    path sample_name
    path script_file
    val data_dir

    output:
    tuple val("${sample_name}"), file("${sample_name}.rds")

    script:
    /* 
     * PARAMETERS FOR SCRIPT.R -> data_dir, outfile_name, sample_name
     */

    """
    Rscript $script_file $data_dir/$sample_name ${sample_name}.rds ${sample_name}
    """
}