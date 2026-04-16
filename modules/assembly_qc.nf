process ASSEMBLY_QC {
    tag "${assembly.simpleName}"
    label 'process_low'
    publishDir "${params.outdir}/01_assembly_qc", mode: 'copy'

    input:
    path assembly

    output:
    path "assembly_stats.json",    emit: stats
    path "gc_content.tsv",         emit: gc
    path "coverage_analysis.tsv",  emit: coverage

    script:
    template 'assembly_qc.py'
}
