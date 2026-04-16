#!/usr/bin/env python3
"""Phase XIII: ML Feature Engineering (48-D vectors)."""
import json, sys
import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler

features_file = "$features"

BASE_FEATURES = [
    "molecular_weight_da","net_charge_pH7","isoelectric_point","mean_hydrophobicity",
    "hydrophobic_ratio","amphipathic_moment","instability_index","boman_index",
    "gravy_score","aromaticity","cys_count","charged_pos_ratio","charged_neg_ratio",
    "polar_ratio","aromatic_ratio","aliphatic_ratio","pro_gly_ratio","length_aa"
]

TOP_DIPEPTIDES = [
    "KK","KR","RK","RR","LL","LI","IL","II","LK","KL","RL","LR",
    "GG","GA","AG","WF","FW","YF","SS","ST","TS","DE","ED","EE",
    "DD","CP","PC","CC","LV","VL"
]

def dipeptide_freqs(seq, dipeptides):
    n = max(len(seq) - 1, 1)
    return {dp: sum(1 for i in range(len(seq)-1) if seq[i:i+2] == dp) / n
            for dp in dipeptides}

try:
    df = pd.read_csv(features_file, sep="\t")
except Exception as e:
    sys.exit(f"ERROR: {e}")

if df.empty:
    with open("ml_feature_matrix.tsv", "w") as fh:
        print("orf_id", file=fh)
    with open("ml_summary.json", "w") as fh:
        json.dump({"error": "no data"}, fh)
    sys.exit(0)

feature_names = BASE_FEATURES + [f"dp_{dp}" for dp in TOP_DIPEPTIDES]
X_rows = []
orf_ids = []
amp_scores = []

for _, row in df.iterrows():
    seq = "".join(c for c in str(row.get("sequence","")).upper() if c in "ACDEFGHIKLMNPQRSTVWY")
    base_vals = []
    for col in BASE_FEATURES:
        try:
            base_vals.append(float(row.get(col, 0) or 0))
        except (ValueError, TypeError):
            base_vals.append(0.0)
    dp = dipeptide_freqs(seq, TOP_DIPEPTIDES)
    dp_vals = [round(dp[d], 5) for d in TOP_DIPEPTIDES]
    X_rows.append(base_vals + dp_vals)
    orf_ids.append(row["orf_id"])
    amp_scores.append(row.get("amp_score", 0))

X = np.array(X_rows, dtype=float)
scaler = MinMaxScaler()
X_norm = scaler.fit_transform(X) if X.shape[0] > 1 else X

with open("ml_feature_matrix.tsv", "w") as fh:
    print("\t".join(["orf_id","amp_score"] + feature_names), file=fh)
    for i, oid in enumerate(orf_ids):
        vals = "\t".join(f"{X_norm[i,j]:.6f}" for j in range(X_norm.shape[1]))
        print(f"{oid}\t{amp_scores[i]}\t{vals}", file=fh)

feat_var = np.var(X_norm, axis=0)
top_idx = np.argsort(feat_var)[::-1][:10]
top_feats = [(feature_names[i], round(float(feat_var[i]), 4)) for i in top_idx]

with open("ml_summary.json", "w") as fh:
    json.dump({
        "n_samples": len(orf_ids),
        "n_features": len(feature_names),
        "feature_dimensions": {
            "physicochemical": len(BASE_FEATURES),
            "dipeptide_composition": len(TOP_DIPEPTIDES),
            "total": len(feature_names)
        },
        "normalization": "MinMax [0,1]",
        "top_variance_features": top_feats
    }, fh, indent=2)

print(f"ML feature engineering complete: {len(orf_ids)} AMPs x {len(feature_names)} features")
