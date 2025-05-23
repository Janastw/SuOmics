#!/usr/bin/env nextflow

process test_process {

    publishDir 'results', mode: 'copy'

    input:
    val variable

    output:

    script:
    """
    
    """
}

/*
 * Pipeline Parameters
 */

params.samples = 'test.csv'

workflow {
    samples_array = ['sample1', 'sample2', 'sample3']
    samples_ch = Channel.fromPath(params.samples)
                        .splitCsv()
    test_process(samples_ch)

}
