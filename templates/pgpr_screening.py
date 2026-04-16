#!/usr/bin/env python3
"""Phase XII: PGPR Trait Identification (6 traits)."""
import json, re, sys
from Bio import SeqIO

orfs_file = "$orfs"

PGPR_TRAITS = {
    "nitrogen_fixation": ("Biological nitrogen fixation", 2, {
        "nifH": (re.compile(r'[DE].{2}[LIVMF]{2}[KR].{3}[GS]', re.I), (290, 330)),
        "nifD": (re.compile(r'[HY].{2}[LIVMF]{3}[DE].{3}[GS]', re.I), (480, 560)),
        "nifK": (re.compile(r'[KR].{4}[LIVMF]{2}[DE]', re.I), (480, 530)),
        "nifA": (re.compile(r'[DE]{2}[LIVMF].{2}[KR]{2}', re.I), (400, 550)),
    }),
    "phosphate_solubilization": ("Inorganic phosphate solubilization", 2, {
        "pqqA": (re.compile(r'[DE].[ST][LIVMF]{2}', re.I), (20, 50)),
        "pqqB": (re.compile(r'[LIVMF]{3}[DE].[GS]', re.I), (300, 400)),
        "gcd":  (re.compile(r'[WF][LIVMF][DE][LIVMF]{2}[HKR]', re.I), (700, 900)),
    }),
    "iaa_production": ("IAA biosynthesis", 1, {
        "ipdC": (re.compile(r'[LIVMF]{2}[ST].{2}[DE][LIVMF]', re.I), (350, 430)),
        "trpE": (re.compile(r'[KR].{3}[LIVMF][DE].{2}[GS]', re.I), (480, 540)),
    }),
    "siderophore_production": ("Siderophore biosynthesis", 2, {
        "dhbA": (re.compile(r'[DE]{2}.[GS][LIVMF]{2}', re.I), (290, 350)),
        "dhbB": (re.compile(r'[LIVMF]{2}[HY][DE]', re.I), (420, 490)),
        "dhbC": (re.compile(r'[DE].{3}[WF].{2}[GS]', re.I), (430, 530)),
        "sfp":  (re.compile(r'[LIVMF][DE].{2}[GS].{3}[LIVMF]', re.I), (220, 280)),
    }),
    "acc_deaminase": ("ACC deaminase (ethylene reduction)", 1, {
        "acdS": (re.compile(r'[ST].{2}[LIVMF]{2}[KR].{3}[DE]', re.I), (330, 380)),
    }),
    "biofilm_formation": ("Biofilm formation / root colonization", 2, {
        "sinR": (re.compile(r'[HY].[LIVMF]{2}[DE][KR]', re.I), (100, 130)),
        "eps":  (re.compile(r'[LIVMF]{2}[KR].{3}[LIVMF]{2}', re.I), (400, 600)),
        "tasA": (re.compile(r'[LIVMF]{3}[ST][LIVMF]{2}', re.I), (260, 310)),
        "tapA": (re.compile(r'[DE].[LIVMF]{3}[KR]', re.I), (200, 280)),
    }),
}

records = list(SeqIO.parse(orfs_file, "fasta"))
results = []
positive_traits = []

for trait, (desc, min_genes, genes) in PGPR_TRAITS.items():
    found = {}
    for gene_name, (motif, (lmin, lmax)) in genes.items():
        for rec in records:
            seq = str(rec.seq).upper()
            if lmin * 0.5 <= len(seq) <= lmax * 1.5 and motif.search(seq):
                if gene_name not in found:
                    found[gene_name] = {"orf_id": rec.id, "length_aa": len(seq)}
                break
    positive = len(found) >= min_genes
    if positive:
        positive_traits.append(trait)
    conf = "HIGH" if len(found) >= len(genes) * 0.7 else ("MODERATE" if positive else "LOW")
    results.append({
        "trait": trait, "description": desc, "positive": "YES" if positive else "NO",
        "confidence": conf, "genes_found": len(found), "genes_total": len(genes),
        "min_required": min_genes,
        "detected_genes": ",".join(found.keys()) if found else "none"
    })

cols = ["trait","description","positive","confidence","genes_found","genes_total","min_required","detected_genes"]
with open("pgpr_traits.tsv", "w") as fh:
    print("\t".join(cols), file=fh)
    for r in results:
        print("\t".join(str(r[c]) for c in cols), file=fh)

n_pos = len(positive_traits)
pgpr_potential = "HIGH" if n_pos >= 4 else ("MODERATE" if n_pos >= 2 else "LOW")
with open("pgpr_summary.json", "w") as fh:
    json.dump({
        "total_traits_screened": len(PGPR_TRAITS),
        "positive_traits": n_pos,
        "positive_trait_list": positive_traits,
        "pgpr_potential": pgpr_potential
    }, fh, indent=2)

print(f"PGPR screening: {n_pos}/{len(PGPR_TRAITS)} traits positive, potential={pgpr_potential}")
