process AMP_SCREENING {
    tag "screening"
    label 'process_medium'
    publishDir "${params.outdir}/03_amp_screening", mode: 'copy'

    input:
    path orfs

    output:
    path "amp_candidates.tsv",   emit: candidates
    path "amp_candidates.fasta", emit: candidates_fasta
    path "nrps_pks_genes.tsv",   emit: nrps_pks

    script:
    template 'amp_screening.py'
}
