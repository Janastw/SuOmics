#!/usr/bin/env nextflow



process generate_seurat_object {
    publishDir "results/${sample_name}/seurat_object", mode: 'copy'

    input:
    path sample_name

    output:
    path "${sample_name}.rds"

    script:
    /* 
     * PARAMETERS FOR SCRIPT.R -> data_dir, outfile_name, sample_name
     */
    """
    Rscript ${projectDir}/scripts/generate_seurat_object.R ${projectDir}/data/${sample_name}/ ${sample_name}.rds ${sample_name}
    """
}