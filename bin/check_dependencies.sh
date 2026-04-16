#!/usr/bin/env bash
# AMP-MineFlow dependency checker
set -e

echo "========================================="
echo " AMP-MineFlow Dependency Check"
echo "========================================="

check_cmd() {
    if command -v "$1" &>/dev/null; then
        echo "  [OK]  $1 $(command -v $1)"
    else
        echo "  [MISSING] $1 — required for $2"
        MISSING=1
    fi
}

MISSING=0

echo ""
echo "Core tools:"
check_cmd nextflow "pipeline execution"
check_cmd python3 "analysis scripts"

echo ""
echo "Bioinformatics tools:"
check_cmd prodigal "ORF prediction (Phase II)"
check_cmd hmmsearch "CAZyme annotation (Phase X)"
check_cmd diamond "sequence alignment"

echo ""
echo "Python packages:"
python3 -c "import Bio; print('  [OK]  biopython', Bio.__version__)" 2>/dev/null || echo "  [MISSING] biopython"
python3 -c "import numpy; print('  [OK]  numpy', numpy.__version__)" 2>/dev/null || echo "  [MISSING] numpy"
python3 -c "import scipy; print('  [OK]  scipy', scipy.__version__)" 2>/dev/null || echo "  [MISSING] scipy"
python3 -c "import sklearn; print('  [OK]  scikit-learn', sklearn.__version__)" 2>/dev/null || echo "  [MISSING] scikit-learn"
python3 -c "import matplotlib; print('  [OK]  matplotlib', matplotlib.__version__)" 2>/dev/null || echo "  [MISSING] matplotlib"
python3 -c "import pandas; print('  [OK]  pandas', pandas.__version__)" 2>/dev/null || echo "  [MISSING] pandas"

echo ""
if [ "$MISSING" -eq 0 ]; then
    echo "All dependencies found!"
else
    echo "Some dependencies are missing. Install via:"
    echo "  conda env create -f environment.yml"
    echo "  conda activate amp-mineflow"
fi
echo "========================================="
