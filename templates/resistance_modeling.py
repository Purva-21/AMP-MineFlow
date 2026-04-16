#!/usr/bin/env python3
"""Phase IX: Resistance Frequency Modeling."""
import sys
from collections import Counter
import pandas as pd

moa_file      = "$moa_results"
families_file = "$families"

MOA_RESISTANCE = {
    "membrane_disruption_carpet":          ("1e-10", "lipid_remodeling", "very_low"),
    "membrane_disruption_barrel_stave":    ("5e-10", "membrane_composition_change", "very_low"),
    "membrane_disruption_toroidal_pore":   ("1e-9",  "surface_charge_modification", "very_low"),
    "cell_wall_synthesis_inhibition":      ("1e-7",  "target_modification", "low"),
    "dna_rna_synthesis_inhibition":        ("1e-6",  "efflux_pump_upregulation", "moderate"),
    "protein_synthesis_inhibition":        ("5e-7",  "ribosome_modification", "low"),
    "enzymatic_activity":                  ("1e-8",  "protease_secretion", "low"),
    "immunomodulatory":                    ("1e-11", "host_adaptation", "extremely_low"),
}

FAMILY_RESISTANCE = {
    "surfactin": ("rare", "none"), "iturin": ("rare", "low"),
    "fengycin": ("rare", "low"), "subtilin": ("moderate", "nisin_class"),
    "mersacidin": ("low", "vancomycin_partial"), "plantazolicin": ("low", "none"),
    "subtilosin": ("rare", "none"), "bacillibactin": ("moderate", "iron_chelation"),
    "bacilysin": ("moderate", "none"), "unclassified": ("unknown", "unknown"),
}

def classify_risk(freq_str):
    try:
        f = float(freq_str)
        if f <= 1e-10: return "very_low"
        if f <= 1e-8:  return "low"
        if f <= 1e-6:  return "moderate"
        return "high"
    except Exception:
        return "unknown"

try:
    moa_df = pd.read_csv(moa_file, sep="\t")
    fam_df = pd.read_csv(families_file, sep="\t")
except Exception as e:
    sys.exit(f"ERROR: {e}")

if moa_df.empty:
    with open("resistance_predictions.tsv", "w") as fh:
        print("orf_id\tprimary_moa\tfamily\tresistance_frequency\toverall_resistance_risk", file=fh)
    sys.exit(0)

merged = moa_df.merge(fam_df[["orf_id", "family"]], on="orf_id", how="left")
merged["family"] = merged["family"].fillna("unclassified")

cols = ["orf_id","primary_moa","family","length_aa","resistance_frequency",
        "resistance_rate_category","resistance_mechanism","acquired_resistance_risk",
        "cross_resistance","overall_resistance_risk"]

with open("resistance_predictions.tsv", "w") as fh:
    print("\t".join(cols), file=fh)
    for _, row in merged.iterrows():
        moa = str(row.get("primary_moa", "membrane_disruption_carpet"))
        fam = str(row.get("family", "unclassified"))
        length = int(row.get("length_aa", 20) or 20)
        freq, mech, rate = MOA_RESISTANCE.get(moa, ("1e-10", "unknown", "very_low"))
        acq_res, cross_res = FAMILY_RESISTANCE.get(fam, ("unknown", "unknown"))
        risk = classify_risk(freq)
        print(f"{row['orf_id']}\t{moa}\t{fam}\t{length}\t{freq}\t{rate}\t{mech}\t{acq_res}\t{cross_res}\t{risk}", file=fh)

risk_counts = Counter(merged.apply(
    lambda r: classify_risk(MOA_RESISTANCE.get(str(r.get("primary_moa",""))or"", ("1e-10","",""))[0]),
    axis=1
))
print(f"Resistance modeling complete. Risk: {dict(risk_counts)}")
