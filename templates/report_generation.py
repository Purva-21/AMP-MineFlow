#!/usr/bin/env python3
"""Phase XIV: Summary Report Generation."""
import json, sys
import pandas as pd

stats_file      = "$stats"
candidates_file = "$candidates"
families_file   = "$families"
features_file   = "$features"
clusters_file   = "$clusters"
moa_file        = "$moa"
spectrum_file   = "$spectrum"
resistance_file = "$resistance"

def safe_json(path):
    try:
        with open(path) as f: return json.load(f)
    except Exception: return {}

def safe_tsv(path):
    try: return pd.read_csv(path, sep="\t")
    except Exception: return pd.DataFrame()

asm_stats   = safe_json(stats_file)
cand_df     = safe_tsv(candidates_file)
fam_df      = safe_tsv(families_file)
feat_df     = safe_tsv(features_file)
clust_df    = safe_tsv(clusters_file)
moa_df      = safe_tsv(moa_file)
spec_df     = safe_tsv(spectrum_file)
resist_df   = safe_tsv(resistance_file)

n_cand = len(cand_df)

score_dist = {}
if not cand_df.empty and "amp_score" in cand_df:
    score_dist = {str(int(k)): int(v) for k, v in cand_df["amp_score"].value_counts().items()}

fam_dist = {}
if not fam_df.empty and "family" in fam_df:
    fam_dist = {str(k): int(v) for k, v in fam_df["family"].value_counts().items()}

physchem = {}
if not feat_df.empty:
    for col in ["net_charge_pH7","hydrophobic_ratio","amphipathic_moment","molecular_weight_da"]:
        if col in feat_df:
            vals = pd.to_numeric(feat_df[col], errors="coerce").dropna()
            if len(vals):
                physchem[col] = {"mean": round(float(vals.mean()),3), "std": round(float(vals.std()),3)}

moa_dist = {}
dominant_moa = "unknown"
if not moa_df.empty and "primary_moa" in moa_df:
    moa_dist = {str(k): int(v) for k, v in moa_df["primary_moa"].value_counts().items()}
    if moa_dist:
        dominant_moa = max(moa_dist, key=moa_dist.get)

spec_summary = {}
if not spec_df.empty and "pathogen" in spec_df and "active" in spec_df:
    for path in spec_df["pathogen"].unique():
        sub = spec_df[spec_df["pathogen"] == path]
        ac = int((sub["active"] == "YES").sum())
        spec_summary[str(path)] = {"tested": len(sub), "active": ac,
                                    "rate": round(ac / max(len(sub),1), 3)}

resist_dist = {}
if not resist_df.empty and "overall_resistance_risk" in resist_df:
    resist_dist = {str(k): int(v) for k, v in resist_df["overall_resistance_risk"].value_counts().items()}

report = {
    "pipeline": "AMP-MineFlow v1.0.0",
    "phases_completed": 14,
    "assembly": {
        "num_contigs": asm_stats.get("num_contigs","NA"),
        "total_length_bp": asm_stats.get("total_length_bp","NA"),
        "n50_bp": asm_stats.get("n50_bp","NA"),
        "gc_percent": asm_stats.get("gc_percent","NA"),
        "mean_coverage": asm_stats.get("mean_coverage","NA"),
    },
    "amp_discovery": {
        "total_candidates": n_cand,
        "score_distribution": score_dist,
        "family_distribution": fam_dist,
        "amp_potential": "HIGH" if n_cand >= 100 else ("MODERATE" if n_cand >= 20 else "LOW")
    },
    "physicochemical_summary": physchem,
    "chemical_space": {"n_clusters": len(clust_df["cluster"].unique()) if not clust_df.empty and "cluster" in clust_df else 0},
    "mechanism_of_action": {"distribution": moa_dist, "dominant": dominant_moa},
    "pathogen_spectrum": spec_summary,
    "resistance_risk": resist_dist,
    "flags": {
        "high_family_diversity": len(fam_dist) >= 4,
        "broad_spectrum": len([p for p,s in spec_summary.items() if s.get("rate",0)>=0.3]) >= 6,
        "low_resistance_risk": (resist_dist.get("very_low",0)+resist_dist.get("low",0)) >
                                (resist_dist.get("moderate",0)+resist_dist.get("high",0))
    }
}

with open("pipeline_summary.json", "w") as fh:
    json.dump(report, fh, indent=2)

print("=" * 50)
print("AMP-MineFlow Complete!")
print("=" * 50)
print(f"Contigs:        {report['assembly']['num_contigs']}")
print(f"AMP candidates: {n_cand}")
print(f"AMP families:   {len(fam_dist)}")
print(f"Dominant MOA:   {dominant_moa}")
print("=" * 50)
