#!/usr/bin/env python3
"""Phase VII: Mechanism of Action Prediction."""
import json, sys
from collections import Counter
import pandas as pd

features_file = "$features"

MOA_MECHANISMS = [
    "membrane_disruption_carpet",
    "membrane_disruption_barrel_stave",
    "membrane_disruption_toroidal_pore",
    "cell_wall_synthesis_inhibition",
    "dna_rna_synthesis_inhibition",
    "protein_synthesis_inhibition",
    "enzymatic_activity",
    "immunomodulatory"
]

def predict_moa(charge, hyd_ratio, moment, length, mw, instab, boman, pI):
    s = {m: 0.0 for m in MOA_MECHANISMS}
    if charge >= 4 and 0.35 <= hyd_ratio <= 0.55: s["membrane_disruption_carpet"] += 3
    if charge >= 2: s["membrane_disruption_carpet"] += 1
    if hyd_ratio >= 0.55 and moment >= 0.35 and length >= 20: s["membrane_disruption_barrel_stave"] += 3
    if length >= 25 and hyd_ratio >= 0.50: s["membrane_disruption_barrel_stave"] += 2
    if 0.3 <= hyd_ratio <= 0.6 and moment >= 0.3 and 2 <= charge <= 6: s["membrane_disruption_toroidal_pore"] += 3
    if 15 <= length <= 35: s["membrane_disruption_toroidal_pore"] += 1
    if mw >= 2000 and 0 <= charge <= 4: s["cell_wall_synthesis_inhibition"] += 2
    if length >= 30 and hyd_ratio <= 0.45: s["cell_wall_synthesis_inhibition"] += 1
    if charge >= 6: s["dna_rna_synthesis_inhibition"] += 3
    if pI >= 10.0 and length <= 25: s["dna_rna_synthesis_inhibition"] += 2
    if 2 <= charge <= 5 and instab < 40: s["protein_synthesis_inhibition"] += 2
    if boman > 5.0: s["protein_synthesis_inhibition"] += 1
    if length >= 50 and abs(charge) <= 2: s["enzymatic_activity"] += 2
    if mw >= 5000: s["enzymatic_activity"] += 1
    if instab < 30 and 1 <= charge <= 4: s["immunomodulatory"] += 2
    if 10 <= length <= 40 and hyd_ratio <= 0.5: s["immunomodulatory"] += 1
    total = sum(s.values()) or 1
    probs = {m: round(v/total, 4) for m, v in s.items()}
    sorted_p = sorted(probs.items(), key=lambda x: x[1], reverse=True)
    return sorted_p[0][0], sorted_p[0][1], sorted_p[1][0], sorted_p[1][1], probs

try:
    df = pd.read_csv(features_file, sep="\t")
except Exception as e:
    sys.exit(f"ERROR: {e}")

if df.empty:
    with open("moa_predictions.tsv", "w") as fh:
        print("orf_id\tprimary_moa\tprimary_confidence\tsecondary_moa\tsecondary_confidence", file=fh)
    with open("moa_summary.json", "w") as fh:
        json.dump({"total": 0, "mechanisms": {}}, fh)
    sys.exit(0)

results = []
for _, row in df.iterrows():
    charge  = float(row.get("net_charge_pH7", 0) or 0)
    hyd     = float(row.get("hydrophobic_ratio", 0.4) or 0.4)
    moment  = float(row.get("amphipathic_moment", 0) or 0)
    length  = int(row.get("length_aa", 20) or 20)
    mw_val  = float(row.get("molecular_weight_da", 2200) or 2200)
    instab  = float(row.get("instability_index", 40) or 40)
    boman   = float(row.get("boman_index", 0) or 0)
    pI      = float(row.get("isoelectric_point", 7.0) or 7.0)
    pm, pc, sm, sc, probs = predict_moa(charge, hyd, moment, length, mw_val, instab, boman, pI)
    r = {"orf_id": row["orf_id"], "length_aa": length,
         "primary_moa": pm, "primary_confidence": pc,
         "secondary_moa": sm, "secondary_confidence": sc}
    r.update({f"prob_{m}": p for m, p in probs.items()})
    results.append(r)

prob_cols = [f"prob_{m}" for m in MOA_MECHANISMS]
cols = ["orf_id","length_aa","primary_moa","primary_confidence","secondary_moa","secondary_confidence"] + prob_cols

with open("moa_predictions.tsv", "w") as fh:
    print("\t".join(cols), file=fh)
    for r in results:
        print("\t".join(str(r.get(c, 0)) for c in cols), file=fh)

moa_dist = Counter(r["primary_moa"] for r in results)
with open("moa_summary.json", "w") as fh:
    json.dump({
        "total_amps": len(results),
        "mechanisms": dict(moa_dist),
        "dominant_mechanism": moa_dist.most_common(1)[0][0] if moa_dist else "unknown"
    }, fh, indent=2)

print(f"MOA prediction complete: {len(results)} AMPs")
