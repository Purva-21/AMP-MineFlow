process ML_FEATURES {
    tag "ml"
    label 'process_medium'
    publishDir "${params.outdir}/13_ml_features", mode: 'copy'

    input:
    path features

    output:
    path "ml_feature_matrix.tsv", emit: ml_matrix
    path "ml_summary.json",       emit: summary

    script:
    template 'ml_features.py'
}
