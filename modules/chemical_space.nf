process CHEMICAL_SPACE {
    tag "chemical_space"
    label 'process_medium'
    publishDir "${params.outdir}/06_chemical_space", mode: 'copy'

    input:
    path features

    output:
    path "chemical_space_clusters.tsv", emit: clusters
    path "pca_summary.json",            emit: pca

    script:
    template 'chemical_space.py'
}
