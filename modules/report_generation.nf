process REPORT_GENERATION {
    tag "report"
    label 'process_low'
    publishDir "${params.outdir}/14_report", mode: 'copy'

    input:
    path stats
    path candidates
    path families
    path features
    path clusters
    path moa
    path spectrum
    path resistance

    output:
    path "pipeline_summary.json", emit: report

    script:
    template 'report_generation.py'
}
