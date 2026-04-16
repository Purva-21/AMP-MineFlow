#!/usr/bin/env python3
"""Generate publication-quality plots for AMP-MineFlow example output."""
import json, csv, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import Counter

OUTDIR = "/home/claude/AMP-MineFlow/example_output"
PLOTDIR = f"{OUTDIR}/14_report"

plt.rcParams.update({'font.size': 10, 'figure.dpi': 150, 'savefig.bbox': 'tight'})

# 1. AMP Family Distribution
with open(f"{OUTDIR}/04_classification/family_summary.json") as f:
    fam = json.load(f)
families = fam["family_distribution"]
fig, ax = plt.subplots(figsize=(8, 5))
colors = plt.cm.Set2(np.linspace(0, 1, len(families)))
bars = ax.barh(list(families.keys()), list(families.values()), color=colors, edgecolor='white')
ax.set_xlabel("Number of AMPs")
ax.set_title("AMP Family Classification Distribution")
for bar, val in zip(bars, families.values()):
    ax.text(bar.get_width() + 1, bar.get_y() + bar.get_height()/2, str(val), va='center', fontsize=9)
plt.tight_layout()
plt.savefig(f"{PLOTDIR}/amp_family_distribution.png")
plt.close()

# 2. Chemical Space t-SNE
clusters = []
with open(f"{OUTDIR}/06_chemical_space/chemical_space_clusters.tsv") as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        clusters.append(row)

fig, ax = plt.subplots(figsize=(8, 6))
unique_clusters = sorted(set(int(c["cluster"]) for c in clusters))
cmap = plt.cm.tab10
for cl in unique_clusters:
    pts = [(float(c["tsne1"]), float(c["tsne2"])) for c in clusters if int(c["cluster"]) == cl]
    x, y = zip(*pts)
    ax.scatter(x, y, c=[cmap(cl)], label=f"Cluster {cl}", alpha=0.7, s=40, edgecolors='white', linewidth=0.3)
ax.set_xlabel("t-SNE 1")
ax.set_ylabel("t-SNE 2")
ax.set_title("AMP Chemical Space (t-SNE + K-means Clustering)")
ax.legend(loc='best', fontsize=8)
plt.tight_layout()
plt.savefig(f"{PLOTDIR}/chemical_space_tsne.png")
plt.close()

# 3. MOA Distribution
with open(f"{OUTDIR}/07_moa_prediction/moa_summary.json") as f:
    moa = json.load(f)
moa_dist = moa["moa_distribution"]
fig, ax = plt.subplots(figsize=(7, 5))
labels = [k.replace("_", " ").title() for k in moa_dist.keys()]
sizes = list(moa_dist.values())
colors = plt.cm.Pastel1(np.linspace(0, 0.8, len(sizes)))
wedges, texts, autotexts = ax.pie(sizes, labels=labels, autopct='%1.1f%%', colors=colors, startangle=90)
ax.set_title("Predicted Mechanism of Action Distribution")
plt.tight_layout()
plt.savefig(f"{PLOTDIR}/moa_distribution.png")
plt.close()

# 4. CAZyme Class Distribution
with open(f"{OUTDIR}/10_cazyme_annotation/cazyme_summary.json") as f:
    caz = json.load(f)
cazyme_dist = caz["class_distribution"]
fig, ax = plt.subplots(figsize=(7, 5))
cls_names = {"GH": "Glycoside\nHydrolases", "GT": "Glycosyl\nTransferases", "PL": "Polysaccharide\nLyases",
             "CE": "Carbohydrate\nEsterases", "AA": "Auxiliary\nActivities", "CBM": "Carbohydrate\nBinding Modules"}
labels = [cls_names.get(k, k) for k in cazyme_dist.keys()]
vals = list(cazyme_dist.values())
colors = ['#E74C3C', '#3498DB', '#2ECC71', '#F39C12', '#9B59B6', '#1ABC9C']
ax.bar(labels, vals, color=colors[:len(vals)], edgecolor='white', width=0.6)
ax.set_ylabel("Count")
ax.set_title("CAZyme Class Distribution")
for i, v in enumerate(vals):
    ax.text(i, v + 0.3, str(v), ha='center', fontsize=9, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{PLOTDIR}/cazyme_distribution.png")
plt.close()

# 5. EPS Operon Map
with open(f"{OUTDIR}/11_eps_pathway/eps_operon_reconstruction.tsv") as f:
    reader = csv.DictReader(f, delimiter='\t')
    eps_genes = list(reader)

fig, ax = plt.subplots(figsize=(14, 3))
cat_colors = {"regulatory": "#E74C3C", "transport": "#3498DB", "biosynthesis": "#2ECC71",
              "modification": "#F39C12", "polymerization": "#9B59B6"}
for i, g in enumerate(eps_genes):
    color = cat_colors.get(g["category"], "#95A5A6")
    alpha = 1.0 if g["detected"] == "yes" else 0.25
    ax.barh(0, 1, left=i, color=color, alpha=alpha, edgecolor='white', height=0.5)
    ax.text(i + 0.5, 0, g["gene"], ha='center', va='center', fontsize=7, fontweight='bold',
            color='white' if alpha == 1 else 'gray')
ax.set_xlim(-0.2, len(eps_genes) + 0.2)
ax.set_ylim(-0.5, 0.8)
ax.set_yticks([])
ax.set_xticks([])
ax.set_title("EPS Operon Reconstruction (epsA–O)")
patches = [mpatches.Patch(color=c, label=l.replace("_", " ").title()) for l, c in cat_colors.items()]
patches.append(mpatches.Patch(facecolor='gray', alpha=0.25, label='Not Detected'))
ax.legend(handles=patches, loc='upper right', fontsize=7, ncol=3)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
plt.tight_layout()
plt.savefig(f"{PLOTDIR}/eps_operon_map.png")
plt.close()

# 6. PGPR Trait Radar Chart
with open(f"{OUTDIR}/12_pgpr_screening/pgpr_traits.tsv") as f:
    reader = csv.DictReader(f, delimiter='\t')
    pgpr = list(reader)

traits = [r["trait"].replace("_", " ").title() for r in pgpr]
completeness = [float(r["completeness_percent"]) for r in pgpr]

angles = np.linspace(0, 2 * np.pi, len(traits), endpoint=False).tolist()
completeness_plot = completeness + [completeness[0]]
angles += [angles[0]]

fig, ax = plt.subplots(figsize=(7, 7), subplot_kw=dict(polar=True))
ax.fill(angles, completeness_plot, alpha=0.25, color='#2ECC71')
ax.plot(angles, completeness_plot, 'o-', color='#2ECC71', linewidth=2)
ax.set_xticks(angles[:-1])
ax.set_xticklabels(traits, fontsize=8)
ax.set_ylim(0, 100)
ax.set_title("PGPR Trait Completeness (%)", pad=20)
plt.tight_layout()
plt.savefig(f"{PLOTDIR}/pgpr_radar.png")
plt.close()

# 7. AMP Score Distribution
scores = []
with open(f"{OUTDIR}/03_amp_screening/amp_candidates.tsv") as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        scores.append(int(row["amp_score"]))

fig, ax = plt.subplots(figsize=(7, 4))
score_counts = Counter(scores)
x = sorted(score_counts.keys())
y = [score_counts[s] for s in x]
colors_score = ['#FFC107' if s < 6 else '#4CAF50' if s < 7 else '#2196F3' for s in x]
ax.bar(x, y, color=colors_score, edgecolor='white')
ax.set_xlabel("AMP Score (out of 8)")
ax.set_ylabel("Number of Candidates")
ax.set_title("AMP Candidate Score Distribution")
ax.set_xticks(x)
plt.tight_layout()
plt.savefig(f"{PLOTDIR}/amp_score_distribution.png")
plt.close()

print(f"Generated 7 publication-quality plots in {PLOTDIR}/")
