process CAZYME_ANNOTATION {
    tag "${assembly.simpleName}"
    label 'process_high'
    publishDir "${params.outdir}/10_cazyme_annotation", mode: 'copy'

    input:
    path assembly

    output:
    path "cazyme_annotations.tsv", emit: cazyme_results
    path "cazyme_summary.json",    emit: summary

    script:
    template 'cazyme_annotation.py'
}
