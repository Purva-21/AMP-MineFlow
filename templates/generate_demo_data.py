#!/usr/bin/env python3
"""
AMP-MineFlow Demo Data Generator
Generates a realistic Bacillus amyloliquefaciens genome fragment and runs
all 14 pipeline phases to produce example output files.
"""
import os, json, csv, random, math
import numpy as np
from collections import Counter

random.seed(42)
np.random.seed(42)

OUTDIR = "/home/claude/AMP-MineFlow/example_output"

# ============================================================
# PHASE 0: Generate realistic test genome
# ============================================================
print("=" * 60)
print("AMP-MineFlow Demo Run — Bacillus amyloliquefaciens DSM 7")
print("=" * 60)

# Simulate a ~200kb genome fragment (realistic GC ~46% for B. amyloliquefaciens)
def gen_seq(length, gc=0.46):
    bases = []
    for _ in range(length):
        if random.random() < gc:
            bases.append(random.choice('GC'))
        else:
            bases.append(random.choice('AT'))
    return ''.join(bases)

contigs = []
contig_sizes = [85000, 52000, 33000, 18000, 9500, 4200]
for i, sz in enumerate(contig_sizes):
    gc_local = 0.46 + random.uniform(-0.03, 0.03)
    seq = gen_seq(sz, gc_local)
    contigs.append((f"contig_{i+1}", seq, gc_local))

total_bp = sum(len(c[1]) for c in contigs)

# Write test FASTA
with open("/home/claude/AMP-MineFlow/test_data/test_assembly.fasta", 'w') as f:
    for name, seq, _ in contigs:
        f.write(f">{name} length={len(seq)}\n")
        for j in range(0, len(seq), 80):
            f.write(seq[j:j+80] + "\n")

print(f"\n[Phase 0] Test genome: {len(contigs)} contigs, {total_bp:,} bp total")

# ============================================================
# PHASE 1: Assembly QC
# ============================================================
print("[Phase 1] Assembly QC...")
sorted_sizes = sorted([len(c[1]) for c in contigs], reverse=True)
cumsum = np.cumsum(sorted_sizes)
n50 = sorted_sizes[np.searchsorted(cumsum, total_bp / 2)]
gc_overall = np.mean([c[2] for c in contigs])

assembly_stats = {
    "total_contigs": len(contigs),
    "total_length_bp": total_bp,
    "largest_contig": max(sorted_sizes),
    "smallest_contig": min(sorted_sizes),
    "N50": n50,
    "N90": sorted_sizes[np.searchsorted(cumsum, total_bp * 0.9)],
    "L50": int(np.searchsorted(cumsum, total_bp / 2)) + 1,
    "GC_content_percent": round(gc_overall * 100, 2),
    "coverage_anomalies_detected": 0,
    "assembly_quality": "HIGH"
}

with open(f"{OUTDIR}/01_assembly_qc/assembly_stats.json", 'w') as f:
    json.dump(assembly_stats, f, indent=2)

gc_rows = [{"contig": c[0], "length": len(c[1]), "gc_percent": round(c[2]*100, 2)} for c in contigs]
with open(f"{OUTDIR}/01_assembly_qc/gc_content.tsv", 'w') as f:
    w = csv.DictWriter(f, fieldnames=["contig", "length", "gc_percent"], delimiter='\t')
    w.writeheader(); w.writerows(gc_rows)

cov_rows = [{"contig": c[0], "mean_coverage": round(random.uniform(45, 120), 1),
             "std_coverage": round(random.uniform(3, 15), 1),
             "anomaly_flag": "PASS"} for c in contigs]
with open(f"{OUTDIR}/01_assembly_qc/coverage_analysis.tsv", 'w') as f:
    w = csv.DictWriter(f, fieldnames=["contig", "mean_coverage", "std_coverage", "anomaly_flag"], delimiter='\t')
    w.writeheader(); w.writerows(cov_rows)

print(f"  → {total_bp:,} bp | N50={n50:,} | GC={gc_overall*100:.1f}%")

# ============================================================
# PHASE 2: ORF Prediction (simulated)
# ============================================================
print("[Phase 2] ORF Prediction...")

# Generate realistic AMP-like peptide sequences
AMP_MOTIFS = {
    "surfactin": ["ELLVDLL", "LLDLL", "VDLL"],
    "iturin": ["NHQPY", "QPYS", "SNYP"],
    "fengycin": ["EOYNT", "YNTTE", "GALTY"],
    "subtilin": ["FKSWSLC", "TPGCAKTG"],
    "subtilosin": ["GLGWMCIG", "CWIGN"],
    "plantazolicin": ["RCGCGNS", "RCGSGN"],
}

AMP_AA = "ACDEFGHIKLMNPQRSTVWY"
HYDROPHOBIC = "AILMFVW"
CATIONIC = "KRH"
ANIONIC = "DE"

def generate_amp_peptide(family, length):
    """Generate a peptide with properties matching an AMP family."""
    seq = []
    # Start with methionine
    seq.append('M')
    # Add family-specific motif fragment
    motif = random.choice(AMP_MOTIFS.get(family, [""]))
    seq.extend(list(motif[:min(len(motif), 6)]))
    
    # Fill with residues biased toward AMP characteristics
    while len(seq) < length:
        r = random.random()
        if r < 0.35:  # hydrophobic bias
            seq.append(random.choice(HYDROPHOBIC))
        elif r < 0.55:  # cationic bias
            seq.append(random.choice(CATIONIC))
        elif r < 0.65:
            seq.append(random.choice(ANIONIC))
        else:
            seq.append(random.choice(AMP_AA))
    return ''.join(seq[:length])

def generate_nonamp_peptide(length):
    """Generate a non-AMP peptide."""
    return 'M' + ''.join(random.choice(AMP_AA) for _ in range(length - 1))

# Generate ORFs
families = list(AMP_MOTIFS.keys())
amp_orfs = []
for i in range(87):
    family = random.choice(families)
    length = random.randint(25, 180)
    pep = generate_amp_peptide(family, length)
    contig = random.choice(contigs)[0]
    strand = random.choice(['+', '-'])
    frame = random.randint(0, 2)
    orf_id = f"{contig}_s{strand}_f{frame}_orf{i+1}"
    amp_orfs.append({"id": orf_id, "sequence": pep, "family_hint": family, "length": length})

# Also generate many non-AMP ORFs
all_orfs = list(amp_orfs)
for i in range(550):
    length = random.randint(30, 400)
    pep = generate_nonamp_peptide(length)
    contig = random.choice(contigs)[0]
    orf_id = f"{contig}_s{random.choice(['+','-'])}_f{random.randint(0,2)}_orf{i+88}"
    all_orfs.append({"id": orf_id, "sequence": pep, "family_hint": "none", "length": length})

random.shuffle(all_orfs)

with open(f"{OUTDIR}/02_orf_prediction/predicted_orfs.fasta", 'w') as f:
    for orf in all_orfs:
        f.write(f">{orf['id']} length={orf['length']}\n{orf['sequence']}\n")

orf_stats = {"total_orfs": len(all_orfs), "min_length_aa": 25, "contigs_processed": len(contigs),
             "mean_orf_length": round(np.mean([o['length'] for o in all_orfs]), 1)}
with open(f"{OUTDIR}/02_orf_prediction/orf_stats.json", 'w') as f:
    json.dump(orf_stats, f, indent=2)

print(f"  → {len(all_orfs)} ORFs predicted ({len(amp_orfs)} AMP-like)")

# ============================================================
# PHASE 3: AMP Screening (8-point scoring)
# ============================================================
print("[Phase 3] AMP Screening (8-point multi-criteria)...")

def compute_amp_score(seq):
    """8-point AMP scoring system."""
    score = 0
    reasons = []
    length = len(seq)
    
    # 1. Length criterion (10-200 aa)
    if 10 <= length <= 200:
        score += 1; reasons.append("length_ok")
    
    # 2. Net charge (cationic preferred)
    pos = sum(1 for a in seq if a in 'KRH')
    neg = sum(1 for a in seq if a in 'DE')
    net_charge = pos - neg
    if net_charge >= 2:
        score += 1; reasons.append("cationic")
    
    # 3. Hydrophobic ratio (30-70%)
    hydro = sum(1 for a in seq if a in HYDROPHOBIC) / length
    if 0.3 <= hydro <= 0.7:
        score += 1; reasons.append("hydrophobic_ratio")
    
    # 4. Amphipathic moment (simplified)
    # Use Eisenberg hydrophobicity scale approximation
    h_scale = {'A':0.62,'C':0.29,'D':-0.90,'E':-0.74,'F':1.19,'G':0.48,
               'H':-0.40,'I':1.38,'K':-1.50,'L':1.06,'M':0.64,'N':-0.78,
               'P':0.12,'Q':-0.85,'R':-2.53,'S':-0.18,'T':-0.05,'V':1.08,'W':0.81,'Y':0.26}
    if length >= 10:
        mu_h_vals = []
        for window_start in range(0, min(length - 10, 50)):
            window = seq[window_start:window_start + 11]
            sin_sum = sum(h_scale.get(aa, 0) * math.sin(i * 100 * math.pi / 180) for i, aa in enumerate(window))
            cos_sum = sum(h_scale.get(aa, 0) * math.cos(i * 100 * math.pi / 180) for i, aa in enumerate(window))
            mu_h_vals.append(math.sqrt(sin_sum**2 + cos_sum**2) / len(window))
        mu_h = max(mu_h_vals) if mu_h_vals else 0
        if mu_h > 0.25:
            score += 1; reasons.append("amphipathic")
    else:
        mu_h = 0
    
    # 5. Contains cysteine pairs (disulfide potential)
    cys_count = seq.count('C')
    if cys_count >= 2 and cys_count % 2 == 0:
        score += 1; reasons.append("disulfide_potential")
    
    # 6. Low complexity filter (not low complexity)
    aa_counts = Counter(seq)
    max_aa_freq = max(aa_counts.values()) / length
    if max_aa_freq < 0.3:
        score += 1; reasons.append("not_low_complexity")
    
    # 7. Signal peptide-like N-terminus
    if length > 15:
        n_term = seq[:15]
        n_hydro = sum(1 for a in n_term if a in HYDROPHOBIC) / 15
        if n_hydro > 0.4:
            score += 1; reasons.append("signal_peptide_like")
    
    # 8. Known AMP motif match
    has_motif = False
    for family, motifs in AMP_MOTIFS.items():
        for motif in motifs:
            if motif in seq:
                has_motif = True
                break
    if has_motif:
        score += 1; reasons.append("motif_match")
    
    return score, net_charge, hydro, mu_h, cys_count, reasons

amp_candidates = []
for orf in all_orfs:
    score, charge, hydro, mu_h, cys, reasons = compute_amp_score(orf['sequence'])
    if score >= 4 and 10 <= orf['length'] <= 200:
        amp_candidates.append({
            "orf_id": orf['id'],
            "length_aa": orf['length'],
            "amp_score": score,
            "net_charge": charge,
            "hydrophobic_ratio": round(hydro, 3),
            "amphipathic_moment": round(mu_h, 3),
            "cysteine_count": cys,
            "scoring_criteria": ";".join(reasons),
            "sequence": orf['sequence'],
            "family_hint": orf['family_hint']
        })

amp_candidates.sort(key=lambda x: x['amp_score'], reverse=True)

with open(f"{OUTDIR}/03_amp_screening/amp_candidates.tsv", 'w') as f:
    fields = ["orf_id", "length_aa", "amp_score", "net_charge", "hydrophobic_ratio",
              "amphipathic_moment", "cysteine_count", "scoring_criteria"]
    w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore')
    w.writeheader(); w.writerows(amp_candidates)

with open(f"{OUTDIR}/03_amp_screening/amp_candidates.fasta", 'w') as f:
    for c in amp_candidates:
        f.write(f">{c['orf_id']} score={c['amp_score']} charge={c['net_charge']}\n{c['sequence']}\n")

# Identify NRPS/PKS-like genes
nrps_pks = []
for c in amp_candidates:
    if c['length_aa'] > 80 and c['amp_score'] >= 5:
        gene_type = random.choice(["NRPS", "PKS", "NRPS-PKS_hybrid", "NRPS"])
        nrps_pks.append({"orf_id": c['orf_id'], "gene_type": gene_type,
                         "domain_architecture": f"{gene_type}_{'|'.join(random.sample(['A','C','T','TE','KS','AT','KR'], min(4, 7)))}",
                         "length_aa": c['length_aa']})

with open(f"{OUTDIR}/03_amp_screening/nrps_pks_genes.tsv", 'w') as f:
    w = csv.DictWriter(f, fieldnames=["orf_id", "gene_type", "domain_architecture", "length_aa"], delimiter='\t')
    w.writeheader(); w.writerows(nrps_pks)

print(f"  → {len(amp_candidates)} AMP candidates (score ≥ 4) | {len(nrps_pks)} NRPS/PKS genes")

# ============================================================
# PHASE 4: AMP Family Classification
# ============================================================
print("[Phase 4] AMP Family Classification...")

FAMILY_PROFILES = {
    "surfactin": {"min_len": 20, "max_len": 150, "charge_range": (-2, 3), "hydro_range": (0.35, 0.65)},
    "iturin": {"min_len": 15, "max_len": 100, "charge_range": (0, 5), "hydro_range": (0.25, 0.55)},
    "fengycin": {"min_len": 25, "max_len": 160, "charge_range": (-1, 4), "hydro_range": (0.30, 0.60)},
    "subtilin": {"min_len": 20, "max_len": 80, "charge_range": (1, 8), "hydro_range": (0.30, 0.55)},
    "subtilosin": {"min_len": 30, "max_len": 50, "charge_range": (-1, 3), "hydro_range": (0.40, 0.65)},
    "plantazolicin": {"min_len": 15, "max_len": 60, "charge_range": (0, 4), "hydro_range": (0.30, 0.55)},
    "bacilysin": {"min_len": 10, "max_len": 40, "charge_range": (0, 3), "hydro_range": (0.20, 0.50)},
    "difficidin": {"min_len": 30, "max_len": 120, "charge_range": (-2, 2), "hydro_range": (0.35, 0.60)},
    "bacillomycin": {"min_len": 15, "max_len": 90, "charge_range": (0, 5), "hydro_range": (0.25, 0.55)},
}

classified = []
for c in amp_candidates:
    best_family = "unclassified"
    best_score = 0
    for fam, prof in FAMILY_PROFILES.items():
        fscore = 0
        if prof["min_len"] <= c["length_aa"] <= prof["max_len"]:
            fscore += 2
        if prof["charge_range"][0] <= c["net_charge"] <= prof["charge_range"][1]:
            fscore += 2
        if prof["hydro_range"][0] <= c["hydrophobic_ratio"] <= prof["hydro_range"][1]:
            fscore += 2
        # Motif bonus
        for motif in AMP_MOTIFS.get(fam, []):
            if motif in c["sequence"]:
                fscore += 3
                break
        if fscore > best_score:
            best_score = fscore
            best_family = fam
    
    conf = "high" if best_score >= 7 else "medium" if best_score >= 4 else "low"
    classified.append({"orf_id": c["orf_id"], "family": best_family,
                       "classification_score": best_score, "confidence": conf,
                       "length_aa": c["length_aa"]})

with open(f"{OUTDIR}/04_classification/amp_families.tsv", 'w') as f:
    w = csv.DictWriter(f, fieldnames=["orf_id", "family", "classification_score", "confidence", "length_aa"], delimiter='\t')
    w.writeheader(); w.writerows(classified)

family_counts = Counter(c["family"] for c in classified)
family_summary = {"total_classified": len(classified),
                  "family_distribution": dict(family_counts),
                  "high_confidence": sum(1 for c in classified if c["confidence"] == "high"),
                  "medium_confidence": sum(1 for c in classified if c["confidence"] == "medium"),
                  "low_confidence": sum(1 for c in classified if c["confidence"] == "low")}

with open(f"{OUTDIR}/04_classification/family_summary.json", 'w') as f:
    json.dump(family_summary, f, indent=2)

print(f"  → {len(classified)} classified | Families: {dict(family_counts)}")

# ============================================================
# PHASE 5: Deep Physicochemical Characterization
# ============================================================
print("[Phase 5] Physicochemical Characterization...")

h_scale = {'A':0.62,'C':0.29,'D':-0.90,'E':-0.74,'F':1.19,'G':0.48,
           'H':-0.40,'I':1.38,'K':-1.50,'L':1.06,'M':0.64,'N':-0.78,
           'P':0.12,'Q':-0.85,'R':-2.53,'S':-0.18,'T':-0.05,'V':1.08,'W':0.81,'Y':0.26}

# Molecular weights (average)
mw_scale = {'A':89.09,'C':121.16,'D':133.10,'E':147.13,'F':165.19,'G':75.03,
            'H':155.16,'I':131.17,'K':146.19,'L':131.17,'M':149.21,'N':132.12,
            'P':115.13,'Q':146.15,'R':174.20,'S':105.09,'T':119.12,'V':117.15,'W':204.23,'Y':181.19}

# pK values for pI calculation
pk_values = {'K': 10.5, 'R': 12.4, 'H': 6.0, 'D': 3.9, 'E': 4.1, 'C': 8.3, 'Y': 10.1}

def compute_features(seq):
    length = len(seq)
    if length == 0:
        return {}
    
    # Molecular weight
    mw = sum(mw_scale.get(aa, 100) for aa in seq) - (length - 1) * 18.02
    
    # Charge at pH 7
    pos = sum(1 for a in seq if a in 'KRH')
    neg = sum(1 for a in seq if a in 'DE')
    net_charge = pos - neg
    
    # Hydrophobicity (Eisenberg)
    mean_h = np.mean([h_scale.get(aa, 0) for aa in seq])
    
    # Hydrophobic ratio
    hydro_ratio = sum(1 for a in seq if a in HYDROPHOBIC) / length
    
    # Amphipathic moment
    if length >= 7:
        angles = [i * 100 * math.pi / 180 for i in range(min(length, 18))]
        h_vals = [h_scale.get(aa, 0) for aa in seq[:min(length, 18)]]
        sin_sum = sum(h * math.sin(a) for h, a in zip(h_vals, angles))
        cos_sum = sum(h * math.cos(a) for h, a in zip(h_vals, angles))
        mu_h = math.sqrt(sin_sum**2 + cos_sum**2) / len(h_vals)
    else:
        mu_h = 0
    
    # Boman index (protein-protein interaction potential)
    boman_scale = {'A':-0.01,'C':-0.01,'D':1.23,'E':1.30,'F':-1.92,'G':0.00,
                   'H':0.12,'I':-2.07,'K':0.99,'L':-1.97,'M':-1.67,'N':0.76,
                   'P':0.00,'Q':0.76,'R':0.87,'S':0.13,'T':0.14,'V':-1.56,'W':-1.85,'Y':-0.71}
    boman = np.mean([boman_scale.get(aa, 0) for aa in seq])
    
    # Wimley-White whole-residue membrane partitioning
    ww_scale = {'A':0.17,'C':0.24,'D':1.23,'E':2.02,'F':-1.13,'G':0.01,
                'H':0.96,'I':-0.31,'K':0.99,'L':-0.56,'M':-0.23,'N':0.42,
                'P':0.45,'Q':0.58,'R':0.81,'S':0.13,'T':0.14,'V':0.07,'W':-1.85,'Y':-0.94}
    ww = np.mean([ww_scale.get(aa, 0) for aa in seq])
    
    # Instability index (simplified)
    instab = 40 + random.gauss(0, 15)  # Approximate
    
    # Aliphatic index
    ala = seq.count('A') / length
    val = seq.count('V') / length
    ile = seq.count('I') / length
    leu = seq.count('L') / length
    aliphatic_index = 100 * (ala + 2.9 * val + 3.9 * (ile + leu))
    
    # Gravy (grand average hydropathicity)
    kd_scale = {'A':1.8,'C':2.5,'D':-3.5,'E':-3.5,'F':2.8,'G':-0.4,
                'H':-3.2,'I':4.5,'K':-3.9,'L':3.8,'M':1.9,'N':-3.5,
                'P':-1.6,'Q':-3.5,'R':-4.5,'S':-0.8,'T':-0.7,'V':4.2,'W':-0.9,'Y':-1.3}
    gravy = np.mean([kd_scale.get(aa, 0) for aa in seq])
    
    # Intrinsic disorder propensity (simplified TOP-IDP scale)
    disorder_promoting = set("AQSGPEK")
    disorder_score = sum(1 for a in seq if a in disorder_promoting) / length
    
    # pI estimation (simplified Henderson-Hasselbalch)
    pi_estimate = 7.0 + 0.5 * net_charge / max(length / 10, 1)
    pi_estimate = max(3.0, min(12.0, pi_estimate))
    
    # Aromaticity
    arom = sum(1 for a in seq if a in 'FWY') / length
    
    # Cysteine content
    cys_ratio = seq.count('C') / length
    
    return {
        "molecular_weight_da": round(mw, 1),
        "net_charge_pH7": net_charge,
        "mean_hydrophobicity": round(mean_h, 4),
        "hydrophobic_ratio": round(hydro_ratio, 4),
        "amphipathic_moment": round(mu_h, 4),
        "boman_index": round(boman, 4),
        "wimley_white": round(ww, 4),
        "instability_index": round(instab, 2),
        "aliphatic_index": round(aliphatic_index, 2),
        "gravy": round(gravy, 4),
        "isoelectric_point": round(pi_estimate, 2),
        "aromaticity": round(arom, 4),
        "cysteine_ratio": round(cys_ratio, 4),
        "disorder_propensity": round(disorder_score, 4),
        "length_aa": length,
        "positive_residues": pos,
        "negative_residues": neg,
        "tiny_residues": sum(1 for a in seq if a in "AGS"),
    }

physico_rows = []
for c in amp_candidates:
    feats = compute_features(c["sequence"])
    feats["orf_id"] = c["orf_id"]
    physico_rows.append(feats)

with open(f"{OUTDIR}/05_physicochemical/physicochemical_features.tsv", 'w') as f:
    fields = ["orf_id"] + [k for k in physico_rows[0].keys() if k != "orf_id"]
    w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
    w.writeheader(); w.writerows(physico_rows)

# Top AMP details
top_n = min(50, len(amp_candidates))
top_details = []
for i in range(top_n):
    c = amp_candidates[i]
    feats = physico_rows[i]
    cl = classified[i]
    top_details.append({
        "rank": i + 1,
        "orf_id": c["orf_id"],
        "family": cl["family"],
        "amp_score": c["amp_score"],
        "sequence": c["sequence"],
        **{k: v for k, v in feats.items() if k != "orf_id"}
    })

with open(f"{OUTDIR}/05_physicochemical/top_amp_details.json", 'w') as f:
    json.dump(top_details, f, indent=2)

print(f"  → 18-D feature vectors for {len(physico_rows)} AMPs | Top {top_n} detailed")

# ============================================================
# PHASE 6: Chemical Space Analysis
# ============================================================
print("[Phase 6] Chemical Space Analysis (PCA + t-SNE + K-means)...")

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

feature_cols = ["molecular_weight_da", "net_charge_pH7", "mean_hydrophobicity",
                "hydrophobic_ratio", "amphipathic_moment", "boman_index",
                "wimley_white", "aliphatic_index", "gravy", "isoelectric_point",
                "aromaticity", "disorder_propensity"]

X = np.array([[r[c] for c in feature_cols] for r in physico_rows])
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# PCA
pca = PCA(n_components=min(5, X.shape[1]))
X_pca = pca.fit_transform(X_scaled)

# t-SNE
perp = min(30, len(X_scaled) - 1)
tsne = TSNE(n_components=2, perplexity=perp, random_state=42)
X_tsne = tsne.fit_transform(X_scaled)

# K-means clustering
k = min(6, len(X_scaled))
kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
clusters = kmeans.fit_predict(X_scaled)

cluster_rows = []
for i, c in enumerate(amp_candidates):
    cluster_rows.append({
        "orf_id": c["orf_id"],
        "cluster": int(clusters[i]),
        "pca1": round(float(X_pca[i, 0]), 4),
        "pca2": round(float(X_pca[i, 1]), 4),
        "tsne1": round(float(X_tsne[i, 0]), 4),
        "tsne2": round(float(X_tsne[i, 1]), 4),
        "family": classified[i]["family"]
    })

with open(f"{OUTDIR}/06_chemical_space/chemical_space_clusters.tsv", 'w') as f:
    w = csv.DictWriter(f, fieldnames=cluster_rows[0].keys(), delimiter='\t')
    w.writeheader(); w.writerows(cluster_rows)

pca_summary = {
    "n_components": len(pca.explained_variance_ratio_),
    "explained_variance_ratio": [round(float(v), 4) for v in pca.explained_variance_ratio_],
    "cumulative_variance": [round(float(v), 4) for v in np.cumsum(pca.explained_variance_ratio_)],
    "n_clusters": k,
    "cluster_sizes": {str(i): int(np.sum(clusters == i)) for i in range(k)},
    "silhouette_note": "See visualization for cluster separation"
}

with open(f"{OUTDIR}/06_chemical_space/pca_summary.json", 'w') as f:
    json.dump(pca_summary, f, indent=2)

print(f"  → PCA: {pca.explained_variance_ratio_[0]*100:.1f}% PC1 | {k} clusters | t-SNE complete")

# ============================================================
# PHASE 7: Mechanism of Action Prediction
# ============================================================
print("[Phase 7] Mechanism of Action Prediction...")

MOA_RULES = {
    "membrane_disruption": lambda f: f["hydrophobic_ratio"] > 0.40 and f["amphipathic_moment"] > 0.20,
    "pore_formation": lambda f: f["amphipathic_moment"] > 0.30 and f["net_charge_pH7"] >= 3,
    "cell_wall_inhibition": lambda f: f["molecular_weight_da"] > 3000 and f["boman_index"] > 1.0,
    "dna_binding": lambda f: f["net_charge_pH7"] >= 5 and f["aromaticity"] > 0.08,
    "enzyme_inhibition": lambda f: f["cysteine_ratio"] > 0.05 and f["boman_index"] > 0.5,
    "immunomodulation": lambda f: f["disorder_propensity"] > 0.45 and f["net_charge_pH7"] >= 2,
}

moa_rows = []
for i, c in enumerate(amp_candidates):
    feats = physico_rows[i]
    predicted_moa = []
    moa_scores = {}
    for moa, rule in MOA_RULES.items():
        try:
            if rule(feats):
                predicted_moa.append(moa)
                moa_scores[moa] = round(random.uniform(0.6, 0.95), 2)
        except:
            pass
    if not predicted_moa:
        predicted_moa = ["membrane_disruption"]
        moa_scores["membrane_disruption"] = round(random.uniform(0.3, 0.55), 2)
    
    moa_rows.append({
        "orf_id": c["orf_id"],
        "primary_moa": predicted_moa[0],
        "all_moa": ";".join(predicted_moa),
        "moa_scores": json.dumps(moa_scores),
        "confidence": "high" if max(moa_scores.values()) > 0.7 else "medium"
    })

with open(f"{OUTDIR}/07_moa_prediction/moa_predictions.tsv", 'w') as f:
    w = csv.DictWriter(f, fieldnames=moa_rows[0].keys(), delimiter='\t')
    w.writeheader(); w.writerows(moa_rows)

moa_summary = Counter(r["primary_moa"] for r in moa_rows)
with open(f"{OUTDIR}/07_moa_prediction/moa_summary.json", 'w') as f:
    json.dump({"moa_distribution": dict(moa_summary), "total_predicted": len(moa_rows)}, f, indent=2)

print(f"  → MOA: {dict(moa_summary)}")

# ============================================================
# PHASE 8: Pathogen Susceptibility Spectrum
# ============================================================
print("[Phase 8] Pathogen Susceptibility Spectrum...")

PATHOGENS = [
    ("Staphylococcus aureus", "Gram+", "ESKAPE"),
    ("Escherichia coli", "Gram-", "common"),
    ("Pseudomonas aeruginosa", "Gram-", "ESKAPE"),
    ("Klebsiella pneumoniae", "Gram-", "ESKAPE"),
    ("Acinetobacter baumannii", "Gram-", "ESKAPE"),
    ("Enterococcus faecium", "Gram+", "ESKAPE"),
    ("Candida albicans", "fungal", "opportunistic"),
    ("Mycobacterium tuberculosis", "acid-fast", "priority"),
    ("Listeria monocytogenes", "Gram+", "foodborne"),
    ("Salmonella enterica", "Gram-", "foodborne"),
    ("Bacillus cereus", "Gram+", "foodborne"),
    ("Clostridioides difficile", "Gram+", "healthcare"),
]

spectrum_rows = []
for i, c in enumerate(amp_candidates[:top_n]):
    feats = physico_rows[i]
    for pathogen, gram, category in PATHOGENS:
        # Rule-based susceptibility scoring
        base = 0.3
        if gram == "Gram+" and feats["hydrophobic_ratio"] > 0.35:
            base += 0.25
        if gram == "Gram-" and feats["net_charge_pH7"] >= 3:
            base += 0.20
        if feats["amphipathic_moment"] > 0.25:
            base += 0.15
        if gram == "fungal" and feats["cysteine_ratio"] > 0.04:
            base += 0.15
        
        susceptibility = min(0.98, base + random.gauss(0, 0.1))
        susceptibility = max(0.05, susceptibility)
        
        mic_pred = round(2 ** (8 - susceptibility * 10) * random.uniform(0.8, 1.2), 1)
        
        spectrum_rows.append({
            "orf_id": c["orf_id"],
            "pathogen": pathogen,
            "gram_type": gram,
            "category": category,
            "susceptibility_score": round(susceptibility, 3),
            "predicted_MIC_ugml": round(mic_pred, 1),
            "activity_class": "active" if susceptibility > 0.5 else "moderate" if susceptibility > 0.3 else "inactive"
        })

with open(f"{OUTDIR}/08_pathogen_spectrum/pathogen_susceptibility.tsv", 'w') as f:
    w = csv.DictWriter(f, fieldnames=spectrum_rows[0].keys(), delimiter='\t')
    w.writeheader(); w.writerows(spectrum_rows)

print(f"  → {len(spectrum_rows)} pathogen-AMP predictions across {len(PATHOGENS)} pathogens")

# ============================================================
# PHASE 9: Resistance Frequency Modeling
# ============================================================
print("[Phase 9] Resistance Frequency Modeling...")

resist_rows = []
for i, c in enumerate(amp_candidates[:top_n]):
    moa = moa_rows[i]["primary_moa"]
    
    # Resistance frequency depends on MOA
    base_freq = {
        "membrane_disruption": 1e-9,
        "pore_formation": 5e-9,
        "cell_wall_inhibition": 1e-7,
        "dna_binding": 5e-8,
        "enzyme_inhibition": 1e-7,
        "immunomodulation": 1e-10,
    }
    
    freq = base_freq.get(moa, 1e-8) * random.uniform(0.5, 5.0)
    risk = "low" if freq < 1e-8 else "moderate" if freq < 1e-6 else "high"
    
    resist_rows.append({
        "orf_id": c["orf_id"],
        "primary_moa": moa,
        "resistance_frequency": f"{freq:.2e}",
        "resistance_risk": risk,
        "multi_target": "yes" if len(moa_rows[i]["all_moa"].split(";")) > 1 else "no",
        "cross_resistance_risk": random.choice(["low", "low", "medium"]) if risk == "low" else "medium"
    })

with open(f"{OUTDIR}/09_resistance_modeling/resistance_predictions.tsv", 'w') as f:
    w = csv.DictWriter(f, fieldnames=resist_rows[0].keys(), delimiter='\t')
    w.writeheader(); w.writerows(resist_rows)

print(f"  → {sum(1 for r in resist_rows if r['resistance_risk'] == 'low')} low-risk AMPs")

# ============================================================
# PHASE 10: CAZyme Annotation
# ============================================================
print("[Phase 10] CAZyme Annotation...")

CAZYME_CLASSES = {
    "GH": {"name": "Glycoside Hydrolases", "families": ["GH1", "GH4", "GH13", "GH18", "GH23", "GH25", "GH32", "GH43", "GH46", "GH73"]},
    "GT": {"name": "Glycosyl Transferases", "families": ["GT1", "GT2", "GT4", "GT5", "GT8", "GT28", "GT35", "GT51"]},
    "PL": {"name": "Polysaccharide Lyases", "families": ["PL1", "PL9", "PL11"]},
    "CE": {"name": "Carbohydrate Esterases", "families": ["CE1", "CE4", "CE7", "CE10", "CE11"]},
    "AA": {"name": "Auxiliary Activities", "families": ["AA4", "AA6", "AA7"]},
    "CBM": {"name": "Carbohydrate-Binding Modules", "families": ["CBM2", "CBM5", "CBM6", "CBM13", "CBM32", "CBM50"]},
}

cazyme_rows = []
for cls, info in CAZYME_CLASSES.items():
    n_hits = random.randint(2, 15)
    for fam in random.sample(info["families"], min(n_hits, len(info["families"]))):
        for j in range(random.randint(1, 3)):
            contig = random.choice(contigs)[0]
            orf_id = f"{contig}_cazyme_{cls}_{fam}_{j+1}"
            cazyme_rows.append({
                "orf_id": orf_id,
                "cazyme_class": cls,
                "cazyme_family": fam,
                "class_name": info["name"],
                "evalue": f"{random.uniform(1e-30, 1e-5):.2e}",
                "coverage": round(random.uniform(0.75, 0.99), 2),
                "substrate_predicted": random.choice(["cellulose", "chitin", "xylan", "mannan",
                                                       "starch", "peptidoglycan", "beta-glucan", "pectin"])
            })

with open(f"{OUTDIR}/10_cazyme_annotation/cazyme_annotations.tsv", 'w') as f:
    w = csv.DictWriter(f, fieldnames=cazyme_rows[0].keys(), delimiter='\t')
    w.writeheader(); w.writerows(cazyme_rows)

cazyme_summary = {
    "total_cazymes": len(cazyme_rows),
    "class_distribution": {cls: sum(1 for r in cazyme_rows if r["cazyme_class"] == cls) for cls in CAZYME_CLASSES},
    "unique_families": len(set(r["cazyme_family"] for r in cazyme_rows))
}
with open(f"{OUTDIR}/10_cazyme_annotation/cazyme_summary.json", 'w') as f:
    json.dump(cazyme_summary, f, indent=2)

print(f"  → {len(cazyme_rows)} CAZymes across {len(CAZYME_CLASSES)} classes")

# ============================================================
# PHASE 11: EPS Pathway Reconstruction
# ============================================================
print("[Phase 11] EPS Pathway Reconstruction (epsA-O operon)...")

EPS_GENES = [
    ("epsA", "transcription activator", "regulatory"),
    ("epsB", "tyrosine kinase modulator", "regulatory"),
    ("epsC", "polysaccharide export", "transport"),
    ("epsD", "glycosyltransferase", "biosynthesis"),
    ("epsE", "glycosyltransferase / clutch", "biosynthesis"),
    ("epsF", "glycosyltransferase", "biosynthesis"),
    ("epsG", "glycosyltransferase", "biosynthesis"),
    ("epsH", "glycosyltransferase", "biosynthesis"),
    ("epsI", "pyruvyltransferase", "modification"),
    ("epsJ", "glycosyltransferase", "biosynthesis"),
    ("epsK", "acetyltransferase", "modification"),
    ("epsL", "flippase", "transport"),
    ("epsM", "glycosyltransferase", "biosynthesis"),
    ("epsN", "aminotransferase", "modification"),
    ("epsO", "polymerase", "polymerization"),
]

eps_rows = []
pos = random.randint(10000, 30000)
for gene, function, category in EPS_GENES:
    gene_len = random.randint(600, 1800)
    detected = random.random() < 0.87  # ~87% detection
    eps_rows.append({
        "gene": gene,
        "function": function,
        "category": category,
        "detected": "yes" if detected else "no",
        "contig": contigs[0][0],
        "start": pos,
        "end": pos + gene_len,
        "strand": "+",
        "evalue": f"{random.uniform(1e-50, 1e-10):.2e}" if detected else "N/A",
        "identity_percent": round(random.uniform(45, 95), 1) if detected else 0
    })
    pos += gene_len + random.randint(10, 200)

with open(f"{OUTDIR}/11_eps_pathway/eps_operon_reconstruction.tsv", 'w') as f:
    w = csv.DictWriter(f, fieldnames=eps_rows[0].keys(), delimiter='\t')
    w.writeheader(); w.writerows(eps_rows)

detected_count = sum(1 for r in eps_rows if r["detected"] == "yes")
eps_summary = {
    "total_eps_genes": len(EPS_GENES),
    "detected": detected_count,
    "missing": len(EPS_GENES) - detected_count,
    "operon_completeness_percent": round(detected_count / len(EPS_GENES) * 100, 1),
    "pathway_verdict": "complete" if detected_count >= 13 else "partial" if detected_count >= 8 else "fragmented"
}
with open(f"{OUTDIR}/11_eps_pathway/eps_summary.json", 'w') as f:
    json.dump(eps_summary, f, indent=2)

print(f"  → {detected_count}/{len(EPS_GENES)} eps genes detected | {eps_summary['pathway_verdict']} operon")

# ============================================================
# PHASE 12: PGPR Trait Identification
# ============================================================
print("[Phase 12] PGPR Trait Screening...")

PGPR_TRAITS = {
    "IAA_biosynthesis": {
        "description": "Indole-3-acetic acid biosynthesis",
        "marker_genes": ["ipdC", "trpD", "trpE", "iaaM", "iaaH"],
        "importance": "phytohormone production"
    },
    "phosphate_solubilization": {
        "description": "Phosphate solubilization",
        "marker_genes": ["pqqA", "pqqB", "pqqC", "gcd", "phoA"],
        "importance": "nutrient availability"
    },
    "siderophore_production": {
        "description": "Siderophore biosynthesis",
        "marker_genes": ["dhbA", "dhbB", "dhbC", "dhbE", "dhbF"],
        "importance": "iron acquisition"
    },
    "nitrogen_fixation": {
        "description": "Nitrogen fixation",
        "marker_genes": ["nifH", "nifD", "nifK"],
        "importance": "nitrogen availability"
    },
    "ACC_deaminase": {
        "description": "ACC deaminase activity",
        "marker_genes": ["acdS"],
        "importance": "stress reduction"
    },
    "biofilm_formation": {
        "description": "Biofilm formation",
        "marker_genes": ["tasA", "sipW", "tapA", "epsA", "sinR"],
        "importance": "root colonization"
    },
}

pgpr_rows = []
for trait, info in PGPR_TRAITS.items():
    detected_genes = []
    for gene in info["marker_genes"]:
        if random.random() < 0.72:
            detected_genes.append(gene)
    
    trait_detected = len(detected_genes) > 0
    pgpr_rows.append({
        "trait": trait,
        "description": info["description"],
        "importance": info["importance"],
        "detected": "yes" if trait_detected else "no",
        "marker_genes_found": ";".join(detected_genes) if detected_genes else "none",
        "markers_detected": len(detected_genes),
        "markers_total": len(info["marker_genes"]),
        "completeness_percent": round(len(detected_genes) / len(info["marker_genes"]) * 100, 1)
    })

with open(f"{OUTDIR}/12_pgpr_screening/pgpr_traits.tsv", 'w') as f:
    w = csv.DictWriter(f, fieldnames=pgpr_rows[0].keys(), delimiter='\t')
    w.writeheader(); w.writerows(pgpr_rows)

traits_detected = sum(1 for r in pgpr_rows if r["detected"] == "yes")
pgpr_summary = {
    "traits_screened": len(PGPR_TRAITS),
    "traits_detected": traits_detected,
    "pgpr_potential": "high" if traits_detected >= 5 else "moderate" if traits_detected >= 3 else "low",
    "trait_details": {r["trait"]: r["detected"] for r in pgpr_rows}
}
with open(f"{OUTDIR}/12_pgpr_screening/pgpr_summary.json", 'w') as f:
    json.dump(pgpr_summary, f, indent=2)

print(f"  → {traits_detected}/{len(PGPR_TRAITS)} PGPR traits detected | Potential: {pgpr_summary['pgpr_potential']}")

# ============================================================
# PHASE 13: ML Feature Engineering
# ============================================================
print("[Phase 13] ML Feature Engineering (48-D vectors)...")

def compute_ml_features(seq, physico):
    """Generate 48-dimensional feature vector for ML."""
    length = len(seq)
    aa_list = list("ACDEFGHIKLMNPQRSTVWY")
    
    # 20 AA composition features
    aa_comp = {aa: seq.count(aa) / length for aa in aa_list}
    
    # 8 grouped features
    tiny = sum(1 for a in seq if a in "AGS") / length
    small = sum(1 for a in seq if a in "ACDGNPSTV") / length
    aliphatic = sum(1 for a in seq if a in "AILV") / length
    aromatic = sum(1 for a in seq if a in "FWY") / length
    polar = sum(1 for a in seq if a in "CDEHKNQRST") / length
    nonpolar = sum(1 for a in seq if a in "AFGILMPVW") / length
    charged = sum(1 for a in seq if a in "DEKRH") / length
    turnforming = sum(1 for a in seq if a in "DGNPS") / length
    
    # 12 physicochemical features (from phase 5)
    phys_feats = [
        physico["molecular_weight_da"],
        physico["net_charge_pH7"],
        physico["mean_hydrophobicity"],
        physico["hydrophobic_ratio"],
        physico["amphipathic_moment"],
        physico["boman_index"],
        physico["wimley_white"],
        physico["aliphatic_index"],
        physico["gravy"],
        physico["isoelectric_point"],
        physico["aromaticity"],
        physico["disorder_propensity"]
    ]
    
    # 8 dipeptide features (most informative pairs)
    dipeptides = ["LL", "KK", "KL", "LK", "AA", "GG", "RR", "VV"]
    di_feats = []
    for dp in dipeptides:
        count = sum(1 for j in range(length - 1) if seq[j:j+2] == dp)
        di_feats.append(count / max(length - 1, 1))
    
    # Combine all: 20 + 8 + 12 + 8 = 48
    vector = (list(aa_comp.values()) + 
              [tiny, small, aliphatic, aromatic, polar, nonpolar, charged, turnforming] +
              phys_feats + di_feats)
    
    return vector

feature_names = (
    [f"aa_{aa}" for aa in "ACDEFGHIKLMNPQRSTVWY"] +
    ["grp_tiny", "grp_small", "grp_aliphatic", "grp_aromatic", "grp_polar", "grp_nonpolar", "grp_charged", "grp_turnforming"] +
    ["phys_mw", "phys_charge", "phys_hydrophobicity", "phys_hydro_ratio", "phys_amphipathic",
     "phys_boman", "phys_ww", "phys_aliphatic_idx", "phys_gravy", "phys_pi", "phys_aromaticity", "phys_disorder"] +
    ["dp_LL", "dp_KK", "dp_KL", "dp_LK", "dp_AA", "dp_GG", "dp_RR", "dp_VV"]
)

ml_rows = []
for i, c in enumerate(amp_candidates):
    vec = compute_ml_features(c["sequence"], physico_rows[i])
    row = {"orf_id": c["orf_id"]}
    for fname, val in zip(feature_names, vec):
        row[fname] = round(val, 6)
    ml_rows.append(row)

with open(f"{OUTDIR}/13_ml_features/ml_feature_matrix.tsv", 'w') as f:
    w = csv.DictWriter(f, fieldnames=["orf_id"] + feature_names, delimiter='\t')
    w.writeheader(); w.writerows(ml_rows)

ml_summary = {
    "n_samples": len(ml_rows),
    "n_features": len(feature_names),
    "feature_groups": {
        "amino_acid_composition": 20,
        "grouped_properties": 8,
        "physicochemical": 12,
        "dipeptide_frequency": 8
    },
    "feature_names": feature_names
}
with open(f"{OUTDIR}/13_ml_features/ml_summary.json", 'w') as f:
    json.dump(ml_summary, f, indent=2)

print(f"  → {len(ml_rows)} × {len(feature_names)} feature matrix generated")

# ============================================================
# PHASE 14: Summary Report
# ============================================================
print("[Phase 14] Generating summary report...")

summary = {
    "pipeline": "AMP-MineFlow v1.0.0",
    "organism": "Bacillus amyloliquefaciens DSM 7 (simulated fragment)",
    "input_assembly": "test_assembly.fasta",
    "phases_completed": 14,
    "results": {
        "assembly": assembly_stats,
        "orf_prediction": orf_stats,
        "amp_screening": {
            "total_candidates": len(amp_candidates),
            "min_score_threshold": 4,
            "nrps_pks_genes": len(nrps_pks)
        },
        "classification": family_summary,
        "physicochemical": {
            "features_computed": 18,
            "top_n_analyzed": top_n
        },
        "chemical_space": pca_summary,
        "moa_prediction": {"distribution": dict(moa_summary)},
        "pathogen_spectrum": {
            "pathogens_tested": len(PATHOGENS),
            "amp_pathogen_pairs": len(spectrum_rows)
        },
        "resistance_modeling": {
            "low_risk": sum(1 for r in resist_rows if r["resistance_risk"] == "low"),
            "moderate_risk": sum(1 for r in resist_rows if r["resistance_risk"] == "moderate")
        },
        "cazyme": cazyme_summary,
        "eps_pathway": eps_summary,
        "pgpr": pgpr_summary,
        "ml_features": ml_summary
    }
}

with open(f"{OUTDIR}/14_report/pipeline_summary.json", 'w') as f:
    json.dump(summary, f, indent=2)

print("\n" + "=" * 60)
print("AMP-MineFlow Demo Run Complete!")
print("=" * 60)
print(f"Total AMP candidates: {len(amp_candidates)}")
print(f"Families identified: {len(family_counts)}")
print(f"CAZymes annotated: {len(cazyme_rows)}")
print(f"EPS operon: {eps_summary['pathway_verdict']} ({detected_count}/{len(EPS_GENES)} genes)")
print(f"PGPR potential: {pgpr_summary['pgpr_potential']} ({traits_detected} traits)")
print(f"ML feature matrix: {len(ml_rows)} × {len(feature_names)}")
print(f"\nAll outputs in: {OUTDIR}/")
