process MOA_PREDICTION {
    tag "moa"
    label 'process_low'
    publishDir "${params.outdir}/07_moa_prediction", mode: 'copy'

    input:
    path features

    output:
    path "moa_predictions.tsv", emit: moa_results
    path "moa_summary.json",    emit: summary

    script:
    template 'moa_prediction.py'
}
