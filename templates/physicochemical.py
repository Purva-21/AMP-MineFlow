#!/usr/bin/env python3
"""Phase V: Physicochemical Characterization (18-D descriptors)."""
import json, sys
import pandas as pd
import numpy as np

candidates_file = "$candidates"
top_n = int("$params.top_n")

HYDROPHOBICITY = {
    "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8, "G": -0.4, "H": -3.2,
    "I": 4.5, "K": -3.9, "L": 3.8, "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5,
    "R": -4.5, "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3
}
MOLWT = {
    "A": 89.09, "C": 121.16, "D": 133.10, "E": 147.13, "F": 165.19, "G": 75.03,
    "H": 155.16, "I": 131.17, "K": 146.19, "L": 131.17, "M": 149.20, "N": 132.12,
    "P": 115.13, "Q": 146.15, "R": 174.20, "S": 105.09, "T": 119.12, "V": 117.15,
    "W": 204.23, "Y": 181.19
}
PKA = {"D": 3.9, "E": 4.07, "H": 6.04, "C": 8.18, "Y": 10.46, "K": 10.54, "R": 12.48,
       "N_term": 8.0, "C_term": 3.1}
STD_AAS = set("ACDEFGHIKLMNPQRSTVWY")
HYD_SET = set("VILMFYWCA")

def clean_seq(seq):
    return "".join(c for c in seq.upper() if c in STD_AAS)

def mw(seq):
    return sum(MOLWT.get(aa, 110.0) for aa in seq) - (len(seq)-1)*18.02

def charge_pH7(seq, pH=7.4):
    c = 1.0/(1+10**(pH-PKA["N_term"])) - 1.0/(1+10**(PKA["C_term"]-pH))
    for aa in seq:
        if aa == "D": c -= 1.0/(1+10**(PKA["D"]-pH))
        elif aa == "E": c -= 1.0/(1+10**(PKA["E"]-pH))
        elif aa == "H": c += 1.0/(1+10**(pH-PKA["H"]))
        elif aa == "C": c -= 1.0/(1+10**(PKA["C"]-pH))
        elif aa == "Y": c -= 1.0/(1+10**(PKA["Y"]-pH))
        elif aa == "K": c += 1.0/(1+10**(pH-PKA["K"]))
        elif aa == "R": c += 1.0/(1+10**(pH-PKA["R"]))
    return round(c, 3)

def isoelectric_point(seq):
    for pH_x10 in range(0, 141):
        if charge_pH7(seq, pH_x10/10.0) <= 0:
            return round(pH_x10/10.0, 2)
    return 14.0

def amphipathic_moment(seq, angle=100):
    hyd = [HYDROPHOBICITY.get(aa, 0.0) for aa in seq]
    angles = [i * angle * np.pi / 180.0 for i in range(len(hyd))]
    hx = sum(h * np.cos(a) for h, a in zip(hyd, angles))
    hy = sum(h * np.sin(a) for h, a in zip(hyd, angles))
    return round(float(np.sqrt(hx**2 + hy**2) / max(len(seq), 1)), 3)

def boman(seq):
    BOMAN = {"L":-4.92,"I":-4.92,"V":-4.04,"F":-2.98,"M":-2.35,"W":-2.33,"A":-1.81,
             "C":-1.28,"G":-0.94,"Y":0.14,"T":2.57,"S":3.40,"H":4.66,"Q":5.54,
             "K":5.81,"N":6.64,"E":6.81,"D":8.72,"R":14.92,"P":0.0}
    return round(float(np.mean([BOMAN.get(aa, 0.0) for aa in seq])), 3) if seq else 0.0

try:
    df = pd.read_csv(candidates_file, sep="\t")
except Exception as e:
    sys.exit(f"ERROR: {e}")

if df.empty:
    with open("physicochemical_features.tsv", "w") as fh:
        print("orf_id", file=fh)
    with open("top_amp_details.json", "w") as fh:
        json.dump({}, fh)
    sys.exit(0)

results = []
for _, row in df.iterrows():
    seq = clean_seq(str(row.get("sequence", "")))
    if len(seq) < 5:
        continue
    total = max(len(seq), 1)
    results.append({
        "orf_id": row["orf_id"],
        "length_aa": len(seq),
        "amp_score": row.get("amp_score", 0),
        "molecular_weight_da": round(mw(seq), 2),
        "net_charge_pH7": charge_pH7(seq),
        "isoelectric_point": isoelectric_point(seq),
        "mean_hydrophobicity": round(float(np.mean([HYDROPHOBICITY.get(aa, 0.0) for aa in seq])), 3),
        "hydrophobic_ratio": round(sum(1 for aa in seq if aa in HYD_SET) / total, 3),
        "amphipathic_moment": amphipathic_moment(seq),
        "instability_index": round(10.0 / total * sum(1.0 for i in range(len(seq)-1)), 2),
        "boman_index": boman(seq),
        "gravy_score": round(float(np.mean([HYDROPHOBICITY.get(aa, 0.0) for aa in seq])), 3),
        "aromaticity": round(sum(1 for aa in seq if aa in "FYW") / total, 4),
        "cys_count": seq.count("C"),
        "disulfide_bonds": seq.count("C") // 2,
        "charged_pos_ratio": round(sum(1 for aa in seq if aa in "KR") / total, 4),
        "charged_neg_ratio": round(sum(1 for aa in seq if aa in "DE") / total, 4),
        "polar_ratio": round(sum(1 for aa in seq if aa in "STNQ") / total, 4),
        "aromatic_ratio": round(sum(1 for aa in seq if aa in "FYW") / total, 4),
        "aliphatic_ratio": round(sum(1 for aa in seq if aa in "AVILM") / total, 4),
        "pro_gly_ratio": round(sum(1 for aa in seq if aa in "PG") / total, 4),
        "motif_match": row.get("motif_match", "none"),
        "sequence": seq
    })

cols = ["orf_id","length_aa","amp_score","molecular_weight_da","net_charge_pH7","isoelectric_point",
        "mean_hydrophobicity","hydrophobic_ratio","amphipathic_moment","instability_index","boman_index",
        "gravy_score","aromaticity","cys_count","disulfide_bonds","charged_pos_ratio","charged_neg_ratio",
        "polar_ratio","aromatic_ratio","aliphatic_ratio","pro_gly_ratio","motif_match","sequence"]

with open("physicochemical_features.tsv", "w") as fh:
    print("\t".join(cols), file=fh)
    for r in results:
        print("\t".join(str(r[c]) for c in cols), file=fh)

top_amps = sorted(results, key=lambda x: float(x["amp_score"]), reverse=True)[:top_n]
with open("top_amp_details.json", "w") as fh:
    json.dump({
        "total_characterized": len(results),
        "top_n": top_n,
        "top_amps": [{k: v for k, v in a.items() if k != "sequence"} for a in top_amps]
    }, fh, indent=2)

print(f"Physicochemical characterization complete: {len(results)} AMPs with 18-D descriptors")
