# Installation Guide

## Option 1: Conda (Recommended)

```bash
git clone https://github.com/Purva-21/AMP-MineFlow.git
cd AMP-MineFlow

# Create environment
conda env create -f environment.yml
conda activate amp-mineflow

# Verify
bash bin/check_dependencies.sh
```

## Option 2: Docker

```bash
# Pull pre-built image
docker pull purva21/amp-mineflow:1.0.0

# Run pipeline
nextflow run main.nf --input assembly.fasta -profile docker
```

## Option 3: Singularity (HPC)

```bash
# Singularity will auto-pull from Docker Hub
nextflow run main.nf --input assembly.fasta -profile singularity
```

## Nextflow Installation

```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash
chmod +x nextflow
mv nextflow ~/bin/   # or /usr/local/bin/

# Verify
nextflow -version
```

## System Requirements

| Component | Minimum | Recommended |
|-----------|---------|-------------|
| CPU       | 4 cores | 8+ cores    |
| RAM       | 8 GB    | 16+ GB      |
| Storage   | 5 GB    | 20+ GB      |
| OS        | Linux / macOS | Linux  |
