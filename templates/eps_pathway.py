#!/usr/bin/env python3
"""Phase XI: EPS Pathway Reconstruction (epsA-O operon)."""
import json, re, sys
from Bio import SeqIO
import pandas as pd

orfs_file    = "$orfs"
cazyme_file  = "$cazyme_results"

EPS_GENES = {
    "epsA": (re.compile(r'[LIVMF]{2}[DE].{2}[KR]|[KR].{3}[DE][LIVMF]', re.I), (350, 550), True,  "signal_transduction_kinase"),
    "epsB": (re.compile(r'[DE].{2}[ST].{3}[LIVMF]{2}', re.I),                  (250, 400), True,  "phosphatase_regulator"),
    "epsC": (re.compile(r'[GS].{2}[ST][GS][KR]', re.I),                        (200, 350), True,  "tyrosine_kinase"),
    "epsD": (re.compile(r'[LIVMF]{2}[DE].[GS][LIVMF]', re.I),                  (300, 500), True,  "glucosyltransferase"),
    "epsE": (re.compile(r'[DE].{3}[LIVMF].{2}[GS].{3}[KR]', re.I),             (350, 550), True,  "glucosyltransferase_GT2"),
    "epsF": (re.compile(r'[LIVMF]{2}[HY][DE].{2,3}[GS]', re.I),                (250, 450), True,  "glycosyltransferase_GT4"),
    "epsG": (re.compile(r'[KR].{2}[LIVMF]{3}[DE]', re.I),                      (150, 300), False, "polysaccharide_polymerase"),
    "epsH": (re.compile(r'[LIVMF][DE].[LIVMF]{2}[HY]', re.I),                  (150, 350), False, "acetyltransferase"),
    "epsI": (re.compile(r'[GS].{1,3}[ST].{2,4}[DE]', re.I),                    (200, 400), False, "glycosyltransferase"),
    "epsJ": (re.compile(r'[LIVMF]{3}[KR][DE].{2}[LIVMF]', re.I),               (350, 550), True,  "polysaccharide_export"),
    "epsK": (re.compile(r'[LIVMF]{2}[KR].{3}[LIVMF]{2}', re.I),                (400, 600), True,  "flippase_Wzx"),
    "epsL": (re.compile(r'[DE]{2}.{2}[LIVMF][KR]', re.I),                      (300, 500), True,  "polymerase_Wzy"),
    "epsM": (re.compile(r'[KR].{2}[DE][LIVMF]{2}[KR]', re.I),                  (250, 400), True,  "chain_length_regulator"),
    "epsN": (re.compile(r'[LIVMF][DE].{2}[GS].{3}[LIVMF]', re.I),              (200, 380), False, "methyltransferase"),
    "epsO": (re.compile(r'[KR].{3}[DE][GS][LIVMF]', re.I),                     (200, 350), False, "dehydratase"),
}

orf_records = list(SeqIO.parse(orfs_file, "fasta"))

try:
    cazyme_df = pd.read_csv(cazyme_file, sep="\t")
    gt_contigs = set(cazyme_df[cazyme_df["cazyme_class"] == "GT"]["contig_id"].tolist()) if not cazyme_df.empty else set()
except Exception:
    gt_contigs = set()

results = []
detected_genes = set()

for gene_name, (motif, (lmin, lmax), essential, func) in EPS_GENES.items():
    best_orf, best_len = None, None
    for rec in orf_records:
        seq = str(rec.seq).upper()
        slen = len(seq)
        if lmin * 0.5 <= slen <= lmax * 1.5 and motif.search(seq):
            if best_orf is None or (lmin <= slen <= lmax):
                best_orf = rec.id
                best_len = slen
    detected = best_orf is not None
    if detected:
        detected_genes.add(gene_name)
    conf = "HIGH" if (detected and best_len and lmin <= best_len <= lmax) else \
           "MODERATE" if detected else "ABSENT"
    results.append({
        "gene": gene_name, "function": func, "essential": "YES" if essential else "NO",
        "detected": "YES" if detected else "NO", "confidence": conf,
        "best_orf_id": best_orf or "NA", "orf_length_aa": best_len or "NA",
        "expected_length_range": f"{lmin}-{lmax}"
    })

cols = ["gene","function","essential","detected","confidence","best_orf_id","orf_length_aa","expected_length_range"]
with open("eps_operon_reconstruction.tsv", "w") as fh:
    print("\t".join(cols), file=fh)
    for r in results:
        print("\t".join(str(r[c]) for c in cols), file=fh)

total = len(EPS_GENES)
completeness = len(detected_genes) / total
with open("eps_summary.json", "w") as fh:
    json.dump({
        "total_eps_genes": total,
        "detected_genes": len(detected_genes),
        "detected_gene_list": sorted(list(detected_genes)),
        "missing_gene_list": sorted([g for g in EPS_GENES if g not in detected_genes]),
        "operon_completeness_pct": round(completeness * 100, 1),
        "eps_potential": "HIGH" if completeness >= 0.8 else ("MODERATE" if completeness >= 0.5 else "LOW")
    }, fh, indent=2)

print(f"EPS pathway: {len(detected_genes)}/{total} genes ({round(completeness*100,1)}% complete)")
