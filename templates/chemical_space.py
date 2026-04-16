#!/usr/bin/env python3
"""Phase VI: Chemical Space Analysis (PCA, t-SNE, K-means)."""
import json, sys
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE

features_file = "$features"
perplexity = int("$params.tsne_perplexity")
k = int("$params.kmeans_k")

NUMERIC_COLS = [
    "molecular_weight_da","net_charge_pH7","isoelectric_point","mean_hydrophobicity",
    "hydrophobic_ratio","amphipathic_moment","instability_index","boman_index",
    "gravy_score","aromaticity","cys_count","charged_pos_ratio","charged_neg_ratio",
    "polar_ratio","aromatic_ratio","aliphatic_ratio","pro_gly_ratio"
]

try:
    df = pd.read_csv(features_file, sep="\t")
except Exception as e:
    sys.exit(f"ERROR: {e}")

if df.empty or len(df) < 3:
    with open("chemical_space_clusters.tsv", "w") as fh:
        print("orf_id\tcluster\tpc1\tpc2\ttsne1\ttsne2", file=fh)
    with open("pca_summary.json", "w") as fh:
        json.dump({"error": "insufficient data"}, fh)
    sys.exit(0)

avail = [c for c in NUMERIC_COLS if c in df.columns]
X_raw = df[avail].fillna(0).values.astype(float)
X_scaled = StandardScaler().fit_transform(X_raw)

n_comp = min(X_scaled.shape[1], X_scaled.shape[0], 10)
pca = PCA(n_components=n_comp)
X_pca = pca.fit_transform(X_scaled)

eff_perp = min(perplexity, max(5, len(df) // 4))
tsne = TSNE(n_components=2, perplexity=eff_perp, random_state=42,
            max_iter=500, learning_rate="auto", init="pca" if len(df) >= 4 else "random")
X_tsne = tsne.fit_transform(X_pca[:, :min(5, X_pca.shape[1])])

eff_k = min(k, len(df))
kmeans = KMeans(n_clusters=eff_k, random_state=42, n_init=10)
clusters = kmeans.fit_predict(X_pca[:, :min(5, X_pca.shape[1])])

with open("chemical_space_clusters.tsv", "w") as fh:
    print("orf_id\tamp_score\tlength_aa\tcluster\tpc1\tpc2\ttsne1\ttsne2", file=fh)
    for i, (_, row) in enumerate(df.iterrows()):
        print(f"{row['orf_id']}\t{row.get('amp_score',0)}\t{row.get('length_aa',0)}\t"
              f"{clusters[i]}\t{X_pca[i,0]:.4f}\t{X_pca[i,1]:.4f}\t"
              f"{X_tsne[i,0]:.4f}\t{X_tsne[i,1]:.4f}", file=fh)

explained = pca.explained_variance_ratio_
cluster_stats = {}
for ci in range(eff_k):
    mask = clusters == ci
    cdf = df[mask]
    cluster_stats[f"cluster_{ci}"] = {
        "size": int(mask.sum()),
        "mean_amp_score": round(float(cdf["amp_score"].mean()), 2) if "amp_score" in cdf else 0,
    }

with open("pca_summary.json", "w") as fh:
    json.dump({
        "n_samples": len(df),
        "n_features": len(avail),
        "explained_variance_ratio": [round(float(v), 4) for v in explained],
        "tsne_perplexity_used": eff_perp,
        "kmeans_k_used": eff_k,
        "cluster_statistics": cluster_stats
    }, fh, indent=2)

print(f"Chemical space: {len(df)} AMPs, {eff_k} clusters, PC1+PC2={round(sum(explained[:2])*100,1)}%")
