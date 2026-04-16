#!/usr/bin/env python3
"""Phase II: Six-frame ORF Prediction."""
import json, sys
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np

assembly_file = "$assembly"
min_aa = int("$params.min_orf_aa")

def find_orfs(seq, min_len, frame_offset, strand, seq_id):
    orfs = []
    start_codons = {"ATG", "GTG", "TTG"}
    stop_codons = {"TAA", "TAG", "TGA"}
    seq_str = str(seq).upper()
    i = 0
    while i + 2 < len(seq_str):
        codon = seq_str[i:i+3]
        if codon in start_codons:
            for j in range(i, len(seq_str) - 2, 3):
                stop = seq_str[j:j+3]
                if stop in stop_codons:
                    orf_nt = seq_str[i:j+3]
                    orf_aa = str(Seq(orf_nt).translate(table=11))[:-1]
                    if len(orf_aa) >= min_len:
                        if strand == "+":
                            start = frame_offset + i
                            end = frame_offset + j + 3
                        else:
                            total = frame_offset + len(seq_str)
                            start = total - (j + 3)
                            end = total - i
                        orfs.append({
                            "id": f"{seq_id}_{strand}{start+1}_{end}",
                            "seq_id": seq_id,
                            "start": start + 1,
                            "end": end,
                            "strand": strand,
                            "length_aa": len(orf_aa),
                            "protein": orf_aa
                        })
                    break
        i += 3
    return orfs

all_orfs = []
records = list(SeqIO.parse(assembly_file, "fasta"))
for rec in records:
    seq = rec.seq.upper()
    for frame in range(3):
        all_orfs.extend(find_orfs(seq[frame:], min_aa, frame, "+", rec.id))
    rev_seq = seq.reverse_complement()
    for frame in range(3):
        all_orfs.extend(find_orfs(rev_seq[frame:], min_aa, frame, "-", rec.id))

for idx, orf in enumerate(all_orfs):
    orf["orf_id"] = f"ORF_{idx+1:06d}"

with open("predicted_orfs.fasta", "w") as fh:
    for orf in all_orfs:
        print(f">{orf['orf_id']} {orf['seq_id']} {orf['start']}-{orf['end']}({orf['strand']}) len={orf['length_aa']}aa", file=fh)
        prot = orf["protein"]
        for chunk in [prot[i:i+60] for i in range(0, len(prot), 60)]:
            print(chunk, file=fh)

lens = [o["length_aa"] for o in all_orfs]
stats = {
    "total_orfs": len(all_orfs),
    "min_orf_aa_cutoff": min_aa,
    "mean_length_aa": round(float(np.mean(lens)), 1) if lens else 0,
    "median_length_aa": float(np.median(lens)) if lens else 0,
    "max_length_aa": max(lens) if lens else 0,
    "forward_strand_orfs": sum(1 for o in all_orfs if o["strand"] == "+"),
    "reverse_strand_orfs": sum(1 for o in all_orfs if o["strand"] == "-"),
    "contigs_with_orfs": len(set(o["seq_id"] for o in all_orfs))
}

with open("orf_stats.json", "w") as fh:
    json.dump(stats, fh, indent=2)

print(f"ORF prediction complete: {len(all_orfs)} ORFs predicted")
