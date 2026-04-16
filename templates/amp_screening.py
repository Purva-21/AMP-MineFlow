#!/usr/bin/env python3
"""Phase III: Multi-criteria AMP Screening (8-point scoring)."""
import json, re, sys
from Bio import SeqIO
import numpy as np

orfs_file = "$orfs"
min_score = int("$params.amp_min_score")
max_len = int("$params.amp_max_len")
min_len = int("$params.amp_min_len")

AMP_MOTIFS = {
    "surfactin":     re.compile(r'[LVI]{2,}[KR]', re.IGNORECASE),
    "iturin":        re.compile(r'[FY][KR][FY]', re.IGNORECASE),
    "fengycin":      re.compile(r'[VILM]{3,}[DE]', re.IGNORECASE),
    "subtilin":      re.compile(r'[ACST]{2}[KR]', re.IGNORECASE),
    "mersacidin":    re.compile(r'C[A-Z]{2,4}C', re.IGNORECASE),
    "plantazolicin": re.compile(r'[AG]{2,}[ILV]', re.IGNORECASE),
    "subtilosin":    re.compile(r'C[A-Z]{6,10}C', re.IGNORECASE),
    "bacillibactin": re.compile(r'[DE]{2}[GS]', re.IGNORECASE),
    "bacilysin":     re.compile(r'[VILM]P[VILM]', re.IGNORECASE),
}

NRPS_PKS_MOTIFS = {
    "NRPS_A_domain": re.compile(r'[LIVMF]{2}[DE][LIVMF][ST]', re.IGNORECASE),
    "NRPS_T_domain": re.compile(r'[LIVMF]{2}[KR][LIVMF]{2}[DE]', re.IGNORECASE),
    "PKS_KS_domain": re.compile(r'[LIVMF]C[DN][LIVMF]{3}', re.IGNORECASE),
    "PKS_AT_domain": re.compile(r'[GS][LIVMF]{2}[HN][LIVMF]', re.IGNORECASE),
}

HYD_SET = set("VILMFYWCA")
POLAR_SET = set("STNQHKRDE")
STD_AAS = set("ACDEFGHIKLMNPQRSTVWY")

def score_amp(seq, seq_len):
    score = 0
    details = {}
    if min_len <= seq_len <= max_len:
        score += 1
        details["length_ok"] = True
    clean_seq = "".join(c for c in seq.upper() if c in STD_AAS)
    if len(clean_seq) < 5:
        return score, details
    charge = clean_seq.count("K") + clean_seq.count("R") - clean_seq.count("D") - clean_seq.count("E")
    details["net_charge"] = charge
    if charge >= 2:
        score += 1
    hyd_ratio = sum(1 for c in clean_seq if c in HYD_SET) / len(clean_seq)
    details["hydrophobic_ratio"] = round(hyd_ratio, 3)
    if 0.30 <= hyd_ratio <= 0.70:
        score += 1
    window = min(11, len(clean_seq))
    max_moment = 0.0
    for i in range(len(clean_seq) - window + 1):
        w = clean_seq[i:i+window]
        hyd_sin = sum(np.sin(j * 100 * np.pi / 180) * (1 if w[j] in HYD_SET else -0.5) for j in range(len(w)))
        hyd_cos = sum(np.cos(j * 100 * np.pi / 180) * (1 if w[j] in HYD_SET else -0.5) for j in range(len(w)))
        moment = np.sqrt(hyd_sin**2 + hyd_cos**2) / len(w)
        if moment > max_moment:
            max_moment = moment
    details["amphipathic_moment"] = round(float(max_moment), 3)
    if max_moment > 0.25:
        score += 1
    cys_count = clean_seq.count("C")
    details["cys_count"] = cys_count
    if cys_count > 0 and cys_count % 2 == 0:
        score += 1
    unique_aa = len(set(clean_seq))
    complexity = unique_aa / 20.0
    details["sequence_complexity"] = round(complexity, 3)
    if complexity > 0.35:
        score += 1
    n_term = clean_seq[:12]
    n_hyd = sum(1 for c in n_term if c in HYD_SET) / max(len(n_term), 1)
    details["n_term_hydrophobicity"] = round(n_hyd, 3)
    if n_hyd >= 0.5:
        score += 1
    matched_motif = None
    for family, pattern in AMP_MOTIFS.items():
        if pattern.search(seq):
            matched_motif = family
            score += 1
            break
    details["motif_match"] = matched_motif
    return score, details

records = list(SeqIO.parse(orfs_file, "fasta"))
candidates = []
nrps_pks_genes = []

for rec in records:
    seq = str(rec.seq).upper()
    seq_len = len(seq)
    score, details = score_amp(seq, seq_len)
    for sig, pattern in NRPS_PKS_MOTIFS.items():
        if pattern.search(seq):
            nrps_pks_genes.append({"orf_id": rec.id, "signature": sig, "length_aa": seq_len})
            break
    if score >= min_score and min_len <= seq_len <= max_len:
        candidates.append({
            "orf_id": rec.id,
            "length_aa": seq_len,
            "amp_score": score,
            "net_charge": details.get("net_charge", "NA"),
            "hydrophobic_ratio": details.get("hydrophobic_ratio", "NA"),
            "amphipathic_moment": details.get("amphipathic_moment", "NA"),
            "cys_count": details.get("cys_count", 0),
            "sequence_complexity": details.get("sequence_complexity", "NA"),
            "motif_match": details.get("motif_match", "none"),
            "sequence": seq,
        })

cols = ["orf_id", "length_aa", "amp_score", "net_charge", "hydrophobic_ratio",
        "amphipathic_moment", "cys_count", "sequence_complexity", "motif_match", "sequence"]

with open("amp_candidates.tsv", "w") as fh:
    print("\t".join(cols), file=fh)
    for c in candidates:
        print("\t".join(str(c[col]) for col in cols), file=fh)

with open("amp_candidates.fasta", "w") as fh:
    for c in candidates:
        print(f">{c['orf_id']} score={c['amp_score']} len={c['length_aa']}aa", file=fh)
        seq = c["sequence"]
        for chunk in [seq[i:i+60] for i in range(0, len(seq), 60)]:
            print(chunk, file=fh)

with open("nrps_pks_genes.tsv", "w") as fh:
    print("orf_id\tsignature\tlength_aa", file=fh)
    for g in nrps_pks_genes:
        print(f"{g['orf_id']}\t{g['signature']}\t{g['length_aa']}", file=fh)

print(f"AMP screening complete: {len(candidates)} candidates (score>={min_score}), {len(nrps_pks_genes)} NRPS/PKS")
