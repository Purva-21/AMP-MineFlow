# Output File Descriptions

## Phase I — Assembly QC (`01_assembly_qc/`)

| File | Format | Description |
|------|--------|-------------|
| `assembly_stats.json` | JSON | N50, N90, L50, GC%, total length, quality verdict |
| `gc_content.tsv` | TSV | Per-contig GC content |
| `coverage_analysis.tsv` | TSV | Coverage statistics and anomaly flags |

## Phase II — ORF Prediction (`02_orf_prediction/`)

| File | Format | Description |
|------|--------|-------------|
| `predicted_orfs.fasta` | FASTA | All predicted ORFs (amino acid sequences) |
| `orf_stats.json` | JSON | Total ORFs, mean length, contigs processed |

## Phase III — AMP Screening (`03_amp_screening/`)

| File | Format | Description |
|------|--------|-------------|
| `amp_candidates.tsv` | TSV | All candidates with 8-point scores, charge, hydrophobicity |
| `amp_candidates.fasta` | FASTA | AMP candidate sequences |
| `nrps_pks_genes.tsv` | TSV | Identified NRPS/PKS biosynthetic genes |

## Phase IV — Classification (`04_classification/`)

| File | Format | Description |
|------|--------|-------------|
| `amp_families.tsv` | TSV | Family assignment per AMP (9 families) |
| `family_summary.json` | JSON | Distribution and confidence breakdown |

## Phase V — Physicochemical (`05_physicochemical/`)

| File | Format | Description |
|------|--------|-------------|
| `physicochemical_features.tsv` | TSV | 18-D descriptor matrix (MW, charge, μH, Boman, GRAVY, etc.) |
| `top_amp_details.json` | JSON | Detailed profiles for top N AMPs |

## Phase VI — Chemical Space (`06_chemical_space/`)

| File | Format | Description |
|------|--------|-------------|
| `chemical_space_clusters.tsv` | TSV | PCA, t-SNE coordinates and K-means cluster assignments |
| `pca_summary.json` | JSON | Explained variance, cluster sizes |

## Phase VII — Mechanism of Action (`07_moa_prediction/`)

| File | Format | Description |
|------|--------|-------------|
| `moa_predictions.tsv` | TSV | Primary and secondary MOA predictions per AMP |
| `moa_summary.json` | JSON | MOA distribution summary |

## Phase VIII — Pathogen Spectrum (`08_pathogen_spectrum/`)

| File | Format | Description |
|------|--------|-------------|
| `pathogen_susceptibility.tsv` | TSV | Predicted activity against 12 ESKAPE+ pathogens |

## Phase IX — Resistance Modeling (`09_resistance_modeling/`)

| File | Format | Description |
|------|--------|-------------|
| `resistance_predictions.tsv` | TSV | Resistance frequency, risk level, cross-resistance |

## Phase X — CAZyme Annotation (`10_cazyme_annotation/`)

| File | Format | Description |
|------|--------|-------------|
| `cazyme_annotations.tsv` | TSV | CAZyme hits (class, family, substrate, e-value) |
| `cazyme_summary.json` | JSON | Class distribution, unique families |

## Phase XI — EPS Pathway (`11_eps_pathway/`)

| File | Format | Description |
|------|--------|-------------|
| `eps_operon_reconstruction.tsv` | TSV | epsA–O gene detection, positions, identity |
| `eps_summary.json` | JSON | Operon completeness verdict |

## Phase XII — PGPR Screening (`12_pgpr_screening/`)

| File | Format | Description |
|------|--------|-------------|
| `pgpr_traits.tsv` | TSV | 6 PGPR traits with marker gene detection |
| `pgpr_summary.json` | JSON | Overall PGPR potential rating |

## Phase XIII — ML Features (`13_ml_features/`)

| File | Format | Description |
|------|--------|-------------|
| `ml_feature_matrix.tsv` | TSV | 48-D feature matrix (AA comp + grouped + physicochemical + dipeptide) |
| `ml_summary.json` | JSON | Feature group breakdown |

## Phase XIV — Report (`14_report/`)

| File | Format | Description |
|------|--------|-------------|
| `pipeline_summary.json` | JSON | Consolidated results from all phases |
| `*.png` | PNG | Publication-quality visualizations |
