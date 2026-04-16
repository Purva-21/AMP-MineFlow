#!/usr/bin/env python3
"""Validate AMP-MineFlow output directory structure and file integrity."""
import os, sys, json, csv

EXPECTED = {
    "01_assembly_qc": ["assembly_stats.json", "gc_content.tsv", "coverage_analysis.tsv"],
    "02_orf_prediction": ["predicted_orfs.fasta", "orf_stats.json"],
    "03_amp_screening": ["amp_candidates.tsv", "amp_candidates.fasta", "nrps_pks_genes.tsv"],
    "04_classification": ["amp_families.tsv", "family_summary.json"],
    "05_physicochemical": ["physicochemical_features.tsv", "top_amp_details.json"],
    "06_chemical_space": ["chemical_space_clusters.tsv", "pca_summary.json"],
    "07_moa_prediction": ["moa_predictions.tsv", "moa_summary.json"],
    "08_pathogen_spectrum": ["pathogen_susceptibility.tsv"],
    "09_resistance_modeling": ["resistance_predictions.tsv"],
    "10_cazyme_annotation": ["cazyme_annotations.tsv", "cazyme_summary.json"],
    "11_eps_pathway": ["eps_operon_reconstruction.tsv", "eps_summary.json"],
    "12_pgpr_screening": ["pgpr_traits.tsv", "pgpr_summary.json"],
    "13_ml_features": ["ml_feature_matrix.tsv", "ml_summary.json"],
    "14_report": ["pipeline_summary.json"],
}

def validate(outdir):
    errors = 0
    for phase, files in EXPECTED.items():
        phase_dir = os.path.join(outdir, phase)
        if not os.path.isdir(phase_dir):
            print(f"  [FAIL] Missing directory: {phase}")
            errors += 1
            continue
        for f in files:
            fp = os.path.join(phase_dir, f)
            if not os.path.isfile(fp):
                print(f"  [FAIL] Missing file: {phase}/{f}")
                errors += 1
            elif os.path.getsize(fp) == 0:
                print(f"  [WARN] Empty file: {phase}/{f}")
            else:
                # Validate JSON files
                if f.endswith('.json'):
                    try:
                        with open(fp) as jf:
                            json.load(jf)
                        print(f"  [OK]   {phase}/{f}")
                    except json.JSONDecodeError:
                        print(f"  [FAIL] Invalid JSON: {phase}/{f}")
                        errors += 1
                # Validate TSV files
                elif f.endswith('.tsv'):
                    with open(fp) as tf:
                        reader = csv.reader(tf, delimiter='\t')
                        rows = sum(1 for _ in reader)
                    print(f"  [OK]   {phase}/{f} ({rows} rows)")
                else:
                    print(f"  [OK]   {phase}/{f}")

    if errors == 0:
        print(f"\nAll {sum(len(v) for v in EXPECTED.values())} expected files present and valid!")
    else:
        print(f"\n{errors} error(s) found.")
    return errors

if __name__ == "__main__":
    outdir = sys.argv[1] if len(sys.argv) > 1 else "./results"
    print(f"Validating AMP-MineFlow output: {outdir}\n")
    sys.exit(validate(outdir))
