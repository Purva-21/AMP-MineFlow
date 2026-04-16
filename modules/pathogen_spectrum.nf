process PATHOGEN_SPECTRUM {
    tag "spectrum"
    label 'process_medium'
    publishDir "${params.outdir}/08_pathogen_spectrum", mode: 'copy'

    input:
    path features

    output:
    path "pathogen_susceptibility.tsv", emit: spectrum

    script:
    template 'pathogen_spectrum.py'
}
