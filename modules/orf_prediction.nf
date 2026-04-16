process ORF_PREDICTION {
    tag "${assembly.simpleName}"
    label 'process_medium'
    publishDir "${params.outdir}/02_orf_prediction", mode: 'copy'

    input:
    path assembly

    output:
    path "predicted_orfs.fasta", emit: orfs
    path "orf_stats.json",       emit: stats

    script:
    template 'orf_prediction.py'
}
