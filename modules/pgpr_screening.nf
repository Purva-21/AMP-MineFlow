process PGPR_SCREENING {
    tag "pgpr"
    label 'process_low'
    publishDir "${params.outdir}/12_pgpr_screening", mode: 'copy'

    input:
    path orfs

    output:
    path "pgpr_traits.tsv",   emit: pgpr_results
    path "pgpr_summary.json", emit: summary

    script:
    template 'pgpr_screening.py'
}
