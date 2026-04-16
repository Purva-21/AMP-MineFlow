#!/usr/bin/env python3
"""Phase IV: AMP Family Classification (9 families)."""
import json, re, sys
from collections import Counter
import pandas as pd

candidates_file = "$candidates"

FAMILY_RULES = {
    "surfactin":     (re.compile(r'[LVI]{2,}[KR][LVI]', re.I), 0, 4, 7, 15, "Lipopeptide; membrane-active"),
    "iturin":        (re.compile(r'[FY][KR][FY]|[KR]{2}[FY]', re.I), 2, 6, 7, 12, "Antifungal lipopeptide"),
    "fengycin":      (re.compile(r'[VILM]{3,}[DE]|[FY]{2}[KR]', re.I), -2, 2, 10, 18, "Antifungal lipopeptide"),
    "subtilin":      (re.compile(r'[ACST]{2}[KR]|[KR]C[A-Z]{2}C', re.I), 3, 8, 30, 50, "Lantibiotic; cell wall inhibitor"),
    "mersacidin":    (re.compile(r'C[A-Z]{2,5}C[A-Z]{2,5}C', re.I), 0, 3, 20, 40, "Type B lantibiotic"),
    "plantazolicin": (re.compile(r'[AG]{2,}[ILV]|[GS]{3,}', re.I), 4, 12, 10, 25, "Bottromycin-type antibacterial"),
    "subtilosin":    (re.compile(r'C[A-Z]{5,12}C[A-Z]{3,8}C', re.I), -1, 2, 32, 42, "Sactipeptide"),
    "bacillibactin": (re.compile(r'[DE]{2}[GS]|[GS][DE]{2}', re.I), -4, 0, 8, 30, "Siderophore-derived"),
    "bacilysin":     (re.compile(r'[VILM]P[VILM]|PP[VILM]', re.I), 0, 3, 2, 10, "Dipeptide antibiotic"),
}

def classify_amp(seq, charge, length):
    best_fam, best_score = "unclassified", 0
    for fam, (pat, cmin, cmax, lmin, lmax, desc) in FAMILY_RULES.items():
        s = 0
        if pat.search(seq.upper()): s += 3
        if cmin <= charge <= cmax: s += 2
        if lmin <= length <= lmax: s += 2
        if s > best_score:
            best_score, best_fam = s, fam
    desc = FAMILY_RULES[best_fam][5] if best_fam != "unclassified" else "Unknown"
    return best_fam, best_score, desc

try:
    df = pd.read_csv(candidates_file, sep="\t")
except Exception as e:
    sys.exit(f"ERROR: {e}")

results = []
for _, row in df.iterrows():
    seq = str(row.get("sequence", ""))
    try:
        charge = float(row.get("net_charge", 0))
    except (ValueError, TypeError):
        charge = 0.0
    length = int(row.get("length_aa", 0))
    fam, conf, desc = classify_amp(seq, charge, length)
    results.append({
        "orf_id": row["orf_id"],
        "length_aa": length,
        "amp_score": row.get("amp_score", 0),
        "family": fam,
        "classification_confidence": conf,
        "description": desc
    })

cols = ["orf_id", "length_aa", "amp_score", "family", "classification_confidence", "description"]
with open("amp_families.tsv", "w") as fh:
    print("\t".join(cols), file=fh)
    for r in results:
        print("\t".join(str(r[c]) for c in cols), file=fh)

family_counts = Counter(r["family"] for r in results)
with open("family_summary.json", "w") as fh:
    json.dump({
        "total_candidates": len(results),
        "families": dict(family_counts),
        "top_family": family_counts.most_common(1)[0][0] if family_counts else "none"
    }, fh, indent=2)

print(f"AMP classification complete: {len(results)} AMPs in {len(set(r['family'] for r in results))} families")
