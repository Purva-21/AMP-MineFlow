process EPS_PATHWAY {
    tag "eps"
    label 'process_low'
    publishDir "${params.outdir}/11_eps_pathway", mode: 'copy'

    input:
    path orfs
    path cazyme_results

    output:
    path "eps_operon_reconstruction.tsv", emit: eps
    path "eps_summary.json",              emit: summary

    script:
    template 'eps_pathway.py'
}
