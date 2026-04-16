process RESISTANCE_MODELING {
    tag "resistance"
    label 'process_low'
    publishDir "${params.outdir}/09_resistance_modeling", mode: 'copy'

    input:
    path moa_results
    path families

    output:
    path "resistance_predictions.tsv", emit: resistance

    script:
    template 'resistance_modeling.py'
}
