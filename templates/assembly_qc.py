#!/usr/bin/env python3
"""Phase I: Assembly QC."""
import json, re, sys
from Bio import SeqIO
import numpy as np

assembly_file = "$assembly"

records = list(SeqIO.parse(assembly_file, "fasta"))
if not records:
    sys.exit("ERROR: No sequences found in assembly FASTA")

lengths = [len(r) for r in records]
total_len = sum(lengths)
gc_counts = [r.seq.upper().count("G") + r.seq.upper().count("C") for r in records]
total_gc = sum(gc_counts)
gc_pct = round(100.0 * total_gc / total_len, 2) if total_len > 0 else 0.0

sorted_lens = sorted(lengths, reverse=True)
cumsum, n50, n90 = 0, 0, 0
for l in sorted_lens:
    cumsum += l
    if n50 == 0 and cumsum >= total_len / 2:
        n50 = l
    if n90 == 0 and cumsum >= total_len * 0.9:
        n90 = l

coverages = []
for r in records:
    hdr = r.description
    m = re.search(r'coverage[_ ]+([0-9.]+)', hdr, re.IGNORECASE) or \
        re.search(r'cov[_ ]+([0-9.]+)', hdr, re.IGNORECASE)
    coverages.append(float(m.group(1)) if m else None)

has_cov = [c for c in coverages if c is not None]
mean_cov = round(float(np.mean(has_cov)), 2) if has_cov else None

anomalies = []
if has_cov and mean_cov:
    for r, cov in zip(records, coverages):
        if cov is not None and (cov > 3 * mean_cov or cov < 0.3 * mean_cov):
            anomalies.append(r.id)

stats = {
    "num_contigs": len(records),
    "total_length_bp": total_len,
    "n50_bp": n50,
    "n90_bp": n90,
    "max_contig_bp": max(lengths),
    "min_contig_bp": min(lengths),
    "mean_contig_bp": round(float(np.mean(lengths)), 1),
    "gc_percent": gc_pct,
    "mean_coverage": mean_cov,
    "coverage_anomaly_contigs": len(anomalies),
    "anomalous_contig_ids": anomalies[:20]
}

with open("assembly_stats.json", "w") as fh:
    json.dump(stats, fh, indent=2)

with open("gc_content.tsv", "w") as fh:
    print("contig_id\tlength_bp\tgc_count\tgc_percent", file=fh)
    for r, gc in zip(records, gc_counts):
        pct = round(100.0 * gc / len(r), 2) if len(r) > 0 else 0.0
        print(f"{r.id}\t{len(r)}\t{gc}\t{pct}", file=fh)

with open("coverage_analysis.tsv", "w") as fh:
    print("contig_id\tlength_bp\tcoverage\tnormalized_coverage\tanomaly_flag", file=fh)
    for r, cov in zip(records, coverages):
        norm = round(cov / mean_cov, 3) if (cov is not None and mean_cov) else "NA"
        flag = "YES" if r.id in anomalies else "NO"
        print(f"{r.id}\t{len(r)}\t{cov if cov is not None else 'NA'}\t{norm}\t{flag}", file=fh)

print(f"Assembly QC complete: {len(records)} contigs, {total_len} bp, N50={n50}, GC={gc_pct}%")
