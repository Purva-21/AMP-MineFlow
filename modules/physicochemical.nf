process PHYSICOCHEMICAL {
    tag "physicochemical"
    label 'process_medium'
    publishDir "${params.outdir}/05_physicochemical", mode: 'copy'

    input:
    path candidates

    output:
    path "physicochemical_features.tsv", emit: features
    path "top_amp_details.json",         emit: top_details

    script:
    template 'physicochemical.py'
}
