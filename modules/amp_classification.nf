process AMP_CLASSIFICATION {
    tag "classification"
    label 'process_low'
    publishDir "${params.outdir}/04_classification", mode: 'copy'

    input:
    path candidates

    output:
    path "amp_families.tsv",    emit: families
    path "family_summary.json", emit: summary

    script:
    template 'amp_classification.py'
}
