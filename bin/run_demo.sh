#!/usr/bin/env bash
# Quick demo run using test data
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"

echo "Running AMP-MineFlow demo..."
echo "Input: ${PIPELINE_DIR}/test_data/test_assembly.fasta"
echo ""

nextflow run "${PIPELINE_DIR}/main.nf" \
    --input "${PIPELINE_DIR}/test_data/test_assembly.fasta" \
    --outdir "${PIPELINE_DIR}/test_results" \
    -profile conda \
    -resume

echo ""
echo "Demo complete! Results in: ${PIPELINE_DIR}/test_results/"
