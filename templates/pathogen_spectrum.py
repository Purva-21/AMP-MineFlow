#!/usr/bin/env python3
"""Phase VIII: Pathogen Susceptibility Spectrum (12 ESKAPE+ pathogens)."""
import sys
import pandas as pd
import numpy as np

features_file = "$features"

PATHOGENS = {
    "Enterococcus_faecium":    ("pos", "thick_wall", 3, True),
    "Staphylococcus_aureus":   ("pos", "thick_wall", 3, True),
    "Klebsiella_pneumoniae":   ("neg", "LPS_complex", 4, True),
    "Acinetobacter_baumannii": ("neg", "LPS_complex", 4, True),
    "Pseudomonas_aeruginosa":  ("neg", "LPS_complex", 4, True),
    "Enterobacter_cloacae":    ("neg", "LPS_simple", 3, True),
    "Escherichia_coli":        ("neg", "LPS_simple", 2, False),
    "Salmonella_typhimurium":  ("neg", "LPS_complex", 2, False),
    "Bacillus_anthracis":      ("pos", "thick_wall", 2, True),
    "Listeria_monocytogenes":  ("pos", "thin_wall", 2, False),
    "Candida_albicans":        ("fungi", "glucan_chitin", 3, True),
    "Cryptococcus_neoformans": ("fungi", "capsule", 3, True),
}

def predict_activity(charge, hyd_ratio, moment, length, cys, pI, gram, membrane, res_level, biofilm):
    score = 0.0
    if gram == "pos":
        if charge >= 2: score += 0.3
        if hyd_ratio >= 0.35: score += 0.2
        if moment >= 0.3: score += 0.2
        if membrane == "thin_wall": score += 0.1
        if length >= 10: score += 0.1
    elif gram == "neg":
        if charge >= 3: score += 0.25
        if hyd_ratio >= 0.45: score += 0.2
        if moment >= 0.35: score += 0.2
        score += 0.15 if membrane == "LPS_simple" else 0.05
    elif gram == "fungi":
        if hyd_ratio >= 0.5: score += 0.3
        if membrane == "glucan_chitin" and cys >= 2: score += 0.2
        if length >= 20: score += 0.1
        if pI >= 9.0: score += 0.1
    res_pen = {1: 0.0, 2: 0.05, 3: 0.10, 4: 0.15}.get(res_level, 0.15)
    score -= res_pen
    if biofilm: score -= 0.05
    score = max(0.0, min(1.0, score))
    if score >= 0.7:
        mic = round(np.random.uniform(1, 8), 1)
        active = True
    elif score >= 0.45:
        mic = round(np.random.uniform(8, 32), 1)
        active = True
    elif score >= 0.25:
        mic = round(np.random.uniform(32, 128), 1)
        active = False
    else:
        mic = ">128"
        active = False
    return round(score, 3), mic, active

np.random.seed(42)

try:
    df = pd.read_csv(features_file, sep="\t")
except Exception as e:
    sys.exit(f"ERROR: {e}")

if df.empty:
    with open("pathogen_susceptibility.tsv", "w") as fh:
        print("orf_id\tpathogen\tgram_stain\tactivity_score\tpredicted_mic_ug_ml\tactive", file=fh)
    sys.exit(0)

with open("pathogen_susceptibility.tsv", "w") as fh:
    print("orf_id\tpathogen\tgram_stain\tactivity_score\tpredicted_mic_ug_ml\tactive", file=fh)
    for _, row in df.iterrows():
        charge  = float(row.get("net_charge_pH7", 0) or 0)
        hyd     = float(row.get("hydrophobic_ratio", 0.4) or 0.4)
        moment  = float(row.get("amphipathic_moment", 0.3) or 0.3)
        length  = int(row.get("length_aa", 20) or 20)
        cys     = int(row.get("cys_count", 0) or 0)
        pI      = float(row.get("isoelectric_point", 7.0) or 7.0)
        for pathogen, (gram, membrane, res, biofilm) in PATHOGENS.items():
            score, mic, active = predict_activity(charge, hyd, moment, length, cys, pI,
                                                   gram, membrane, res, biofilm)
            print(f"{row['orf_id']}\t{pathogen}\t{gram}\t{score}\t{mic}\t{'YES' if active else 'NO'}", file=fh)

print(f"Pathogen spectrum complete: {len(df)} AMPs x {len(PATHOGENS)} pathogens")
