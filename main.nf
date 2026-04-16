#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
========================================================================================
    AMP-MineFlow v1.0.0
    Antimicrobial Peptide Mining & Multi-Functional Genome Analysis Pipeline
========================================================================================
    GitHub  : https://github.com/Purva-21/AMP-MineFlow
    Author  : Purva Gohil
    License : MIT
========================================================================================
*/

params.input           = null
params.outdir          = "${launchDir}/results"
params.min_orf_aa      = 30
params.amp_min_score   = 4
params.amp_max_len     = 200
params.amp_min_len     = 10
params.tsne_perplexity = 30
params.kmeans_k        = 6
params.top_n           = 50
params.skip_cazyme     = false
params.skip_pgpr       = false
params.skip_ml         = false
params.help            = false

if (params.help) {
    log.info """
    =========================================
     AMP-MineFlow v1.0.0
    =========================================
    Usage:
      nextflow run main.nf --input assembly.fasta [options]

    Required:
      --input           Path to genome assembly FASTA file

    Optional:
      --outdir          Output directory [default: ./results]
      --min_orf_aa      Minimum ORF length in amino acids [default: 30]
      --amp_min_score   Minimum AMP candidate score [default: 4]
      --amp_max_len     Maximum AMP length [default: 200]
      --amp_min_len     Minimum AMP length [default: 10]
      --top_n           Top N AMPs for detailed analysis [default: 50]
      --tsne_perplexity t-SNE perplexity [default: 30]
      --kmeans_k        K-means clusters [default: 6]
      --skip_cazyme     Skip CAZyme annotation [default: false]
      --skip_pgpr       Skip PGPR screening [default: false]
      --skip_ml         Skip ML features [default: false]
    =========================================
    """.stripIndent()
    exit 0
}

if (!params.input) {
    error "ERROR: Please provide an input FASTA file with --input"
}

log.info """
=========================================
 AMP-MineFlow v1.0.0
=========================================
 input       : ${params.input}
 outdir      : ${params.outdir}
 min_orf_aa  : ${params.min_orf_aa}
 amp_score   : ${params.amp_min_score}
 amp_range   : ${params.amp_min_len}-${params.amp_max_len} aa
 top_n       : ${params.top_n}
 skip_cazyme : ${params.skip_cazyme}
 skip_pgpr   : ${params.skip_pgpr}
 skip_ml     : ${params.skip_ml}
=========================================
"""

include { ASSEMBLY_QC }        from './modules/assembly_qc'
include { ORF_PREDICTION }     from './modules/orf_prediction'
include { AMP_SCREENING }      from './modules/amp_screening'
include { AMP_CLASSIFICATION } from './modules/amp_classification'
include { PHYSICOCHEMICAL }    from './modules/physicochemical'
include { CHEMICAL_SPACE }     from './modules/chemical_space'
include { MOA_PREDICTION }     from './modules/moa_prediction'
include { PATHOGEN_SPECTRUM }  from './modules/pathogen_spectrum'
include { RESISTANCE_MODELING }from './modules/resistance_modeling'
include { CAZYME_ANNOTATION }  from './modules/cazyme_annotation'
include { EPS_PATHWAY }        from './modules/eps_pathway'
include { PGPR_SCREENING }    from './modules/pgpr_screening'
include { ML_FEATURES }        from './modules/ml_features'
include { REPORT_GENERATION }  from './modules/report_generation'
include { VISUALIZE }          from './modules/visualize'

workflow {
    ch_assembly = Channel.fromPath(params.input, checkIfExists: true)

    ASSEMBLY_QC(ch_assembly)
    ORF_PREDICTION(ch_assembly)
    AMP_SCREENING(ORF_PREDICTION.out.orfs)
    AMP_CLASSIFICATION(AMP_SCREENING.out.candidates)
    PHYSICOCHEMICAL(AMP_SCREENING.out.candidates)
    CHEMICAL_SPACE(PHYSICOCHEMICAL.out.features)
    MOA_PREDICTION(PHYSICOCHEMICAL.out.features)
    PATHOGEN_SPECTRUM(PHYSICOCHEMICAL.out.features)
    RESISTANCE_MODELING(MOA_PREDICTION.out.moa_results, AMP_CLASSIFICATION.out.families)

    if (!params.skip_cazyme) {
        CAZYME_ANNOTATION(ch_assembly)
        EPS_PATHWAY(ORF_PREDICTION.out.orfs, CAZYME_ANNOTATION.out.cazyme_results)
    }

    if (!params.skip_pgpr) {
        PGPR_SCREENING(ORF_PREDICTION.out.orfs)
    }

    if (!params.skip_ml) {
        ML_FEATURES(PHYSICOCHEMICAL.out.features)
    }

    REPORT_GENERATION(
        ASSEMBLY_QC.out.stats,
        AMP_SCREENING.out.candidates,
        AMP_CLASSIFICATION.out.families,
        PHYSICOCHEMICAL.out.features,
        CHEMICAL_SPACE.out.clusters,
        MOA_PREDICTION.out.moa_results,
        PATHOGEN_SPECTRUM.out.spectrum,
        RESISTANCE_MODELING.out.resistance
    )

    VISUALIZE(
        REPORT_GENERATION.out.report,
        CHEMICAL_SPACE.out.clusters,
        AMP_CLASSIFICATION.out.summary,
        CAZYME_ANNOTATION.out.summary,
        MOA_PREDICTION.out.summary,
        EPS_PATHWAY.out.eps,
        PGPR_SCREENING.out.pgpr_results
    )
}

workflow.onComplete {
    log.info """
    =========================================
     AMP-MineFlow Complete!
     Status   : ${workflow.success ? 'SUCCESS' : 'FAILED'}
     Duration : ${workflow.duration}
     Output   : ${params.outdir}
    =========================================
    """
}
