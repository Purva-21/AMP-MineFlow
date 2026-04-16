# AMP-MineFlow 🧬⚔️

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-%E2%89%A53.9-blue.svg)](https://python.org)

**A comprehensive Nextflow DSL2 pipeline for genome-scale antimicrobial peptide discovery, CAZyme annotation, EPS pathway reconstruction, PGPR trait identification, and AI/ML-ready feature engineering from bacterial genome assemblies.**

Designed for *Bacillus* and related genera isolated from diverse ecological niches.

---

## Pipeline Overview

AMP-MineFlow performs **14 analytical phases** in a reproducible, parallelized Nextflow workflow:

```
┌─────────────────────────────────────────────────────────────────┐
│  INPUT: Genome Assembly (FASTA)                                 │
├─────────────────────────────────────────────────────────────────┤
│  Phase I     │ Assembly QC (N50, GC, coverage anomaly detection)│
│  Phase II    │ Six-frame ORF Prediction                         │
│  Phase III   │ Multi-criteria AMP Screening (8-point scoring)   │
│  Phase IV    │ AMP Family Classification (9 families)           │
│  Phase V     │ Deep Physicochemical Characterization (18-D)     │
│  Phase VI    │ Chemical Space Analysis (t-SNE, PCA, K-means)    │
│  Phase VII   │ Mechanism of Action Prediction                   │
│  Phase VIII  │ Pathogen Susceptibility Spectrum (12 pathogens)  │
│  Phase IX    │ Resistance Frequency Modeling                    │
│  Phase X     │ CAZyme Annotation (6 classes)                    │
│  Phase XI    │ EPS Pathway Reconstruction (epsA-O operon)       │
│  Phase XII   │ PGPR Trait Identification (6 functional traits)  │
│  Phase XIII  │ ML Feature Engineering (48-D vectors)            │
│  Phase XIV   │ Summary Report Generation                        │
├─────────────────────────────────────────────────────────────────┤
│  OUTPUT: TSV, JSON, FASTA, PNG reports                          │
└─────────────────────────────────────────────────────────────────┘
```

## Quick Start

### Prerequisites

- [Nextflow](https://www.nextflow.io/) >= 23.04.0
- Python >= 3.9 with BioPython, NumPy, SciPy, scikit-learn, matplotlib
- OR: Docker / Singularity / Conda (recommended)

### Installation

```bash
git clone https://github.com/Purva-21/AMP-MineFlow.git
cd AMP-MineFlow
```

### Run with your genome

```bash
# Conda
nextflow run main.nf --input your_assembly.fasta -profile conda

# Docker
nextflow run main.nf --input your_assembly.fasta -profile docker

# Singularity (HPC)
nextflow run main.nf --input your_assembly.fasta -profile singularity

# SLURM + Singularity
nextflow run main.nf --input your_assembly.fasta -profile singularity,slurm
```

### Test run

```bash
nextflow run main.nf -profile test,conda
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | *required* | Path to genome assembly FASTA |
| `--outdir` | `./results` | Output directory |
| `--min_orf_aa` | `30` | Minimum ORF length (amino acids) |
| `--amp_min_score` | `4` | Minimum AMP candidate score (1–8) |
| `--amp_max_len` | `200` | Maximum AMP candidate length |
| `--amp_min_len` | `10` | Minimum AMP candidate length |
| `--top_n` | `50` | Top N AMPs for detailed analysis |
| `--tsne_perplexity` | `30` | t-SNE perplexity parameter |
| `--kmeans_k` | `6` | Number of K-means clusters |
| `--skip_cazyme` | `false` | Skip CAZyme annotation |
| `--skip_pgpr` | `false` | Skip PGPR trait screening |
| `--skip_ml` | `false` | Skip ML feature engineering |

## Output Structure

```
results/
├── 01_assembly_qc/          # Genome metrics, GC, coverage
├── 02_orf_prediction/       # Predicted ORFs (FASTA + stats)
├── 03_amp_screening/        # AMP candidates, NRPS/PKS genes
├── 04_classification/       # Family assignments (9 families)
├── 05_physicochemical/      # 18-D descriptor matrix
├── 06_chemical_space/       # PCA, t-SNE, K-means clusters
├── 07_moa_prediction/       # Mechanism of action predictions
├── 08_pathogen_spectrum/    # Activity vs 12 ESKAPE+ pathogens
├── 09_resistance_modeling/  # Resistance frequency estimates
├── 10_cazyme_annotation/    # CAZyme annotations (6 classes)
├── 11_eps_pathway/          # EPS operon reconstruction
├── 12_pgpr_screening/       # PGPR trait identification
├── 13_ml_features/          # 48-D ML-ready feature matrix
└── 14_report/               # Summary JSON + visualizations
```

## Example Output

The `example_output/` directory contains results from a demo run on a simulated *B. amyloliquefaciens* genome fragment (~200 kb, 6 contigs). Key results:

- **320 AMP candidates** identified (score ≥ 4/8)
- **6 AMP families**: surfactin (150), iturin (88), fengycin (58), subtilin (11), plantazolicin (8), subtilosin (5)
- **55 CAZymes** across 6 classes
- **EPS operon**: 14/15 genes detected (complete)
- **PGPR potential**: HIGH (6/6 traits detected)
- **ML matrix**: 320 × 48 features ready for downstream modeling

### Example Visualizations

See `example_output/14_report/` for:
- AMP family distribution
- Chemical space t-SNE clustering
- Mechanism of action breakdown
- CAZyme class distribution
- EPS operon map
- PGPR trait radar chart
- AMP score distribution

## AMP Scoring System

Each candidate is scored on 8 criteria (max score = 8):

1. **Length**: 10–200 amino acids
2. **Net charge**: Cationic (≥ +2)
3. **Hydrophobic ratio**: 30–70%
4. **Amphipathic moment**: μH > 0.25
5. **Disulfide potential**: Even cysteine pairs
6. **Sequence complexity**: Not low-complexity
7. **Signal peptide**: Hydrophobic N-terminus
8. **Motif match**: Known AMP family motifs

## Citation

If you use AMP-MineFlow, please cite:

```
Gohil P. (2026). AMP-MineFlow: Antimicrobial Peptide Mining & Multi-Functional
Genome Analysis Pipeline. GitHub: https://github.com/Purva-21/AMP-MineFlow
```

## License

MIT License — see [LICENSE](LICENSE).

## Contributing

Contributions welcome! Please open an issue or submit a pull request.
