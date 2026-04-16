process VISUALIZE {
    tag "plots"
    label 'process_low'
    publishDir "${params.outdir}/14_report/plots", mode: 'copy'

    input:
    path summary
    path clusters
    path fam_summary
    path cazyme_summary
    path moa_summary
    path eps
    path pgpr

    output:
    path "*.png", emit: plots

    script:
    """
    python3 ${projectDir}/bin/generate_plots.py \
        --results ${params.outdir} \
        --outdir .
    """
}
