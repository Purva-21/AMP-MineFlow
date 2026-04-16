#!/usr/bin/env python3
"""Phase X: CAZyme Annotation (6 classes)."""
import json, re, sys
from collections import Counter
from Bio import SeqIO

assembly_file = "$assembly"

CAZYME_SIGS = {
    "GH": {
        "GH13": (re.compile(r'[DIVN]VV[HY][NG][HQ]', re.I), "Alpha-amylase"),
        "GH18": (re.compile(r'[FY]D.{0,3}[GS]..[LIVMF]', re.I), "Chitinase"),
        "GH3":  (re.compile(r'[LIV]T[DN][TS][GS][LIVMF]', re.I), "Beta-glucosidase"),
        "GH5":  (re.compile(r'[LIVMF]{2}[NQ][EP][STA]', re.I), "Cellulase"),
        "GH73": (re.compile(r'[YW][LIVMF]H[TS]', re.I), "Peptidoglycan hydrolase"),
    },
    "GT": {
        "GT2":  (re.compile(r'[DE].{1,3}[LIVMF].{1,2}[GS].{1,3}[KR]', re.I), "Cellulose synthase"),
        "GT4":  (re.compile(r'[LIVMF]{2}[HY][DE].{2,3}[GS]', re.I), "Sucrose synthase"),
        "GT51": (re.compile(r'[ST].{2}[KR].{3}[DE]', re.I), "Murein polymerase"),
        "GT83": (re.compile(r'[KR]{2}[LIVMF][DE]', re.I), "EPS synthase"),
    },
    "PL": {
        "PL1":  (re.compile(r'[KR].{2}[YW][LIVMF]{2}', re.I), "Pectate lyase"),
        "PL6":  (re.compile(r'[DE].{3}[WF].{2}[GS]', re.I), "Alginate lyase"),
    },
    "CE": {
        "CE1":  (re.compile(r'[GS].{1,3}[ST].{2,4}[DE]', re.I), "Feruloyl esterase"),
        "CE4":  (re.compile(r'[DE].{2}H[LIVMF]', re.I), "Deacetylase"),
    },
    "AA": {
        "AA3":  (re.compile(r'[WF][LIVMF][DE][LIVMF]{2}[HKR]', re.I), "GMC oxidoreductase"),
        "AA10": (re.compile(r'H[LIVMF]{2}[DE].{2}H', re.I), "LPMO"),
    },
    "CBM": {
        "CBM5":  (re.compile(r'[WY].{2,4}[WY].{3,6}[GS]', re.I), "Chitin-binding"),
        "CBM50": (re.compile(r'[LIVMF]{2}[WY][LIVMF][KR]', re.I), "LysM domain"),
    },
}

CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

def translate(nt):
    aa = []
    for i in range(0, len(nt)-2, 3):
        c = CODON_TABLE.get(nt[i:i+3].upper(), 'X')
        if c == '*': break
        aa.append(c)
    return ''.join(aa)

def revcomp(seq):
    comp = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    return ''.join(comp.get(b,'N') for b in reversed(seq.upper()))

records = list(SeqIO.parse(assembly_file, "fasta"))
seen = set()
unique_cazymes = []

for rec in records:
    seq_nt = str(rec.seq).upper()
    for strand, nt in [("+", seq_nt), ("-", revcomp(seq_nt))]:
        for frame in range(3):
            pos = frame
            while pos + 90 < len(nt):
                prot = translate(nt[pos:pos+300])
                if len(prot) >= 30:
                    for cclass, families in CAZYME_SIGS.items():
                        for fam_id, (pat, desc) in families.items():
                            if pat.search(prot):
                                key = (rec.id, cclass, fam_id, pos // 300)
                                if key not in seen:
                                    seen.add(key)
                                    unique_cazymes.append({
                                        "contig_id": rec.id, "start": pos+1,
                                        "end": pos+len(prot)*3, "strand": strand,
                                        "cazyme_class": cclass, "family": fam_id,
                                        "description": desc, "length_aa": len(prot)
                                    })
                                break
                pos += 150

cols = ["contig_id","start","end","strand","cazyme_class","family","description","length_aa"]
with open("cazyme_annotations.tsv", "w") as fh:
    print("\t".join(cols), file=fh)
    for c in unique_cazymes:
        print("\t".join(str(c[col]) for col in cols), file=fh)

class_counts = Counter(c["cazyme_class"] for c in unique_cazymes)
with open("cazyme_summary.json", "w") as fh:
    json.dump({
        "total_cazymes": len(unique_cazymes),
        "classes": dict(class_counts),
        "contigs_with_cazymes": len(set(c["contig_id"] for c in unique_cazymes))
    }, fh, indent=2)

print(f"CAZyme annotation complete: {len(unique_cazymes)} CAZymes in {len(class_counts)} classes")
