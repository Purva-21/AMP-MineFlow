#!/usr/bin/env python3
"""
AMP-MineFlow Publication-Quality Plot Generator
Generates 9 visualizations from pipeline output.
Usage: python generate_plots.py --results <results_dir> --outdir <plot_dir>
"""
import json, sys, argparse, os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap
import warnings
warnings.filterwarnings('ignore')

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 11,
    'axes.titlesize': 13,
    'axes.titleweight': 'bold',
    'axes.spines.top': False,
    'axes.spines.right': False,
    'figure.dpi': 150,
    'savefig.bbox': 'tight',
    'savefig.dpi': 150,
})

PALETTE = ['#2E86AB', '#A23B72', '#F18F01', '#C73E1D', '#3B1F2B',
           '#44BBA4', '#E94F37', '#393E41', '#F5A623', '#7B2D8B']

def load_json(path):
    try:
        with open(path) as f:
            return json.load(f)
    except Exception:
        return {}

def load_tsv(path):
    try:
        return pd.read_csv(path, sep='\t')
    except Exception:
        return pd.DataFrame()

# ─────────────────────────────────────────────
# Plot 1: AMP Family Distribution (pie + bar)
# ─────────────────────────────────────────────
def plot_amp_families(fam_summary, outdir):
    fam = fam_summary.get('families', {})
    if not fam:
        return
    labels = [k.capitalize() for k in fam.keys()]
    values = list(fam.values())
    colors = PALETTE[:len(labels)]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle('AMP Family Distribution', fontsize=15, fontweight='bold', y=1.01)

    # Pie
    wedges, texts, autotexts = ax1.pie(
        values, labels=None, colors=colors,
        autopct=lambda p: f'{p:.1f}%' if p > 3 else '',
        startangle=140, pctdistance=0.75,
        wedgeprops=dict(linewidth=1.5, edgecolor='white')
    )
    for at in autotexts:
        at.set_fontsize(9)
    ax1.legend(wedges, [f'{l} ({v:,})' for l, v in zip(labels, values)],
               loc='center left', bbox_to_anchor=(-0.35, 0.5), fontsize=9)
    ax1.set_title('Proportion by Family')

    # Bar
    bars = ax2.barh(labels[::-1], values[::-1], color=colors[::-1], edgecolor='white', linewidth=0.5)
    for bar, val in zip(bars, values[::-1]):
        ax2.text(bar.get_width() + max(values)*0.01, bar.get_y() + bar.get_height()/2,
                 f'{val:,}', va='center', fontsize=9)
    ax2.set_xlabel('Number of AMP Candidates')
    ax2.set_title('Count per Family')
    ax2.set_xlim(0, max(values) * 1.15)

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, '01_amp_family_distribution.png'))
    plt.close()
    print('  [1/9] AMP family distribution saved')

# ─────────────────────────────────────────────
# Plot 2: AMP Score Distribution
# ─────────────────────────────────────────────
def plot_score_distribution(report, outdir):
    score_dist = report.get('amp_discovery', {}).get('score_distribution', {})
    if not score_dist:
        return
    scores = sorted(int(k) for k in score_dist.keys())
    counts = [score_dist[str(s)] for s in scores]
    total = sum(counts)
    colors = ['#FFC300' if s < 6 else '#2E86AB' if s < 8 else '#C73E1D' for s in scores]

    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(scores, counts, color=colors, edgecolor='white', linewidth=1.5, width=0.7)
    for bar, cnt in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + total*0.003,
                f'{cnt:,}\n({100*cnt/total:.1f}%)', ha='center', va='bottom', fontsize=9)
    ax.set_xlabel('AMP Score (out of 8 criteria)')
    ax.set_ylabel('Number of Candidates')
    ax.set_title('AMP Candidate Score Distribution')
    ax.set_xticks(scores)
    legend_patches = [
        mpatches.Patch(color='#FFC300', label='Moderate (4-5)'),
        mpatches.Patch(color='#2E86AB', label='High (6-7)'),
        mpatches.Patch(color='#C73E1D', label='Perfect (8)'),
    ]
    ax.legend(handles=legend_patches, framealpha=0.8)
    ax.set_ylim(0, max(counts) * 1.18)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, '02_amp_score_distribution.png'))
    plt.close()
    print('  [2/9] AMP score distribution saved')

# ─────────────────────────────────────────────
# Plot 3: Chemical Space (t-SNE + PCA)
# ─────────────────────────────────────────────
def plot_chemical_space(clusters_df, outdir):
    if clusters_df.empty or 'tsne1' not in clusters_df.columns:
        return
    # Sample up to 5000 for speed
    df = clusters_df.sample(min(5000, len(clusters_df)), random_state=42)
    n_clust = df['cluster'].nunique()
    colors = plt.cm.tab10(np.linspace(0, 1, n_clust))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle('Chemical Space Analysis', fontsize=15, fontweight='bold')

    # t-SNE coloured by cluster
    for ci, color in zip(sorted(df['cluster'].unique()), colors):
        mask = df['cluster'] == ci
        ax1.scatter(df.loc[mask, 'tsne1'], df.loc[mask, 'tsne2'],
                    c=[color], s=5, alpha=0.5, label=f'Cluster {ci} (n={mask.sum():,})')
    ax1.set_xlabel('t-SNE 1')
    ax1.set_ylabel('t-SNE 2')
    ax1.set_title(f't-SNE Clustering (n={len(df):,} sampled)')
    ax1.legend(markerscale=3, fontsize=8, framealpha=0.7)

    # PCA coloured by AMP score
    if 'amp_score' in df.columns:
        sc = ax2.scatter(df['pc1'], df['pc2'], c=df['amp_score'],
                         cmap='RdYlGn', s=5, alpha=0.5, vmin=4, vmax=8)
        cbar = plt.colorbar(sc, ax=ax2)
        cbar.set_label('AMP Score')
    else:
        ax2.scatter(df['pc1'], df['pc2'], s=5, alpha=0.5, color='#2E86AB')
    ax2.set_xlabel('PC1')
    ax2.set_ylabel('PC2')
    ax2.set_title('PCA Space (coloured by AMP score)')

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, '03_chemical_space.png'))
    plt.close()
    print('  [3/9] Chemical space saved')

# ─────────────────────────────────────────────
# Plot 4: Mechanism of Action
# ─────────────────────────────────────────────
def plot_moa(moa_summary, outdir):
    moa_dist = moa_summary.get('mechanisms', {})
    if not moa_dist:
        return
    labels = [k.replace('_', ' ').title() for k in moa_dist.keys()]
    values = list(moa_dist.values())
    total = sum(values)
    colors = PALETTE[:len(labels)]

    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.barh(labels[::-1], values[::-1], color=colors[::-1], edgecolor='white')
    for bar, val in zip(bars, values[::-1]):
        ax.text(bar.get_width() + total*0.003, bar.get_y() + bar.get_height()/2,
                f'{val:,} ({100*val/total:.1f}%)', va='center', fontsize=9)
    ax.set_xlabel('Number of AMPs')
    ax.set_title('Predicted Mechanisms of Action')
    ax.set_xlim(0, max(values) * 1.2)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, '04_mechanism_of_action.png'))
    plt.close()
    print('  [4/9] MOA distribution saved')

# ─────────────────────────────────────────────
# Plot 5: Pathogen Susceptibility Heatmap
# ─────────────────────────────────────────────
def plot_pathogen_spectrum(spectrum_data, outdir):
    if not spectrum_data:
        return
    pathogens = list(spectrum_data.keys())
    rates = [spectrum_data[p]['rate'] for p in pathogens]
    active = [spectrum_data[p]['active'] for p in pathogens]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle('Pathogen Susceptibility Spectrum', fontsize=15, fontweight='bold')

    short_names = [p.replace('_', ' ') for p in pathogens]
    colors = ['#C73E1D' if r >= 0.5 else '#F18F01' if r >= 0.3 else '#2E86AB' for r in rates]
    bars = ax1.barh(short_names[::-1], rates[::-1], color=colors[::-1], edgecolor='white')
    ax1.axvline(0.3, color='gray', linestyle='--', linewidth=1, alpha=0.7, label='30% threshold')
    ax1.axvline(0.5, color='red', linestyle='--', linewidth=1, alpha=0.7, label='50% threshold')
    ax1.set_xlabel('Activity Rate (fraction of AMPs active)')
    ax1.set_title('Activity Rate per Pathogen')
    ax1.set_xlim(0, 1.1)
    for bar, rate in zip(bars, rates[::-1]):
        ax1.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height()/2,
                 f'{rate:.1%}', va='center', fontsize=9)
    ax1.legend(fontsize=8)

    # Heatmap-style grid
    matrix = np.array(rates).reshape(1, -1)
    cmap = LinearSegmentedColormap.from_list('activity', ['#EEF2FF', '#2E86AB', '#C73E1D'])
    im = ax2.imshow(matrix, cmap=cmap, aspect='auto', vmin=0, vmax=1)
    ax2.set_xticks(range(len(pathogens)))
    ax2.set_xticklabels([p.replace('_', '\n') for p in pathogens], rotation=45, ha='right', fontsize=8)
    ax2.set_yticks([])
    plt.colorbar(im, ax=ax2, label='Activity Rate', orientation='horizontal', pad=0.3)
    for j, (r, a) in enumerate(zip(rates, active)):
        ax2.text(j, 0, f'{r:.1%}\n({a:,})', ha='center', va='center',
                 color='white' if r > 0.4 else 'black', fontsize=7.5, fontweight='bold')
    ax2.set_title('Activity Matrix')

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, '05_pathogen_spectrum.png'))
    plt.close()
    print('  [5/9] Pathogen spectrum saved')

# ─────────────────────────────────────────────
# Plot 6: CAZyme Class Distribution
# ─────────────────────────────────────────────
def plot_cazymes(cazyme_summary, outdir):
    classes = cazyme_summary.get('classes', {})
    if not classes:
        return
    labels = list(classes.keys())
    values = list(classes.values())
    colors = PALETTE[:len(labels)]
    full_names = {
        'GH': 'Glycoside\nHydrolases',
        'GT': 'GlycosylTransferases',
        'PL': 'Polysaccharide\nLyases',
        'CE': 'Carbohydrate\nEsterases',
        'AA': 'Auxiliary\nActivities',
        'CBM': 'Carbohydrate\nBinding Modules'
    }
    display = [full_names.get(l, l) for l in labels]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle('CAZyme Annotation', fontsize=15, fontweight='bold')

    wedges, _, autotexts = ax1.pie(
        values, colors=colors,
        autopct=lambda p: f'{p:.1f}%' if p > 5 else '',
        startangle=90, pctdistance=0.78,
        wedgeprops=dict(linewidth=1.5, edgecolor='white')
    )
    ax1.legend(wedges, [f'{d} ({v})' for d, v in zip(display, values)],
               loc='center left', bbox_to_anchor=(-0.4, 0.5), fontsize=9)
    ax1.set_title(f'Total: {sum(values):,} CAZymes')

    bars = ax2.bar(display, values, color=colors, edgecolor='white', linewidth=1)
    for bar, val in zip(bars, values):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(values)*0.01,
                 str(val), ha='center', fontsize=10, fontweight='bold')
    ax2.set_ylabel('Count')
    ax2.set_title('CAZyme Class Counts')
    ax2.set_ylim(0, max(values) * 1.15)

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, '06_cazyme_distribution.png'))
    plt.close()
    print('  [6/9] CAZyme distribution saved')

# ─────────────────────────────────────────────
# Plot 7: EPS Operon Map
# ─────────────────────────────────────────────
def plot_eps_operon(eps_df, outdir):
    if eps_df.empty:
        return
    genes = eps_df['gene'].tolist()
    detected = eps_df['detected'].tolist()
    essential = eps_df['essential'].tolist()
    n = len(genes)

    fig, ax = plt.subplots(figsize=(14, 4))
    fig.suptitle('EPS Operon Reconstruction (epsA–epsO)', fontsize=14, fontweight='bold')

    for i, (gene, det, ess) in enumerate(zip(genes, detected, essential)):
        color = '#2E86AB' if det == 'YES' and ess == 'YES' else \
                '#44BBA4' if det == 'YES' else \
                '#E0E0E0'
        edge = '#C73E1D' if ess == 'YES' else '#888888'
        # Arrow shape
        arrow = mpatches.FancyArrow(i * 1.2, 0.5, 0.9, 0,
                                     width=0.35, head_width=0.35, head_length=0.15,
                                     fc=color, ec=edge, linewidth=1.5, length_includes_head=True)
        ax.add_patch(arrow)
        ax.text(i * 1.2 + 0.45, 0.5, gene, ha='center', va='center',
                fontsize=8, fontweight='bold',
                color='white' if det == 'YES' else '#666666')
        if ess == 'YES':
            ax.text(i * 1.2 + 0.45, 0.05, '★', ha='center', va='bottom',
                    fontsize=7, color='#C73E1D')

    ax.set_xlim(-0.3, n * 1.2)
    ax.set_ylim(-0.2, 1.2)
    ax.axis('off')

    legend_patches = [
        mpatches.Patch(color='#2E86AB', label='Detected (essential)'),
        mpatches.Patch(color='#44BBA4', label='Detected (non-essential)'),
        mpatches.Patch(color='#E0E0E0', label='Not detected'),
        mpatches.Patch(color='white', label='★ = Essential gene', linewidth=0),
    ]
    ax.legend(handles=legend_patches, loc='upper right', fontsize=9, framealpha=0.9)

    detected_count = sum(1 for d in detected if d == 'YES')
    ax.set_title(f'{detected_count}/{n} genes detected', pad=5)

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, '07_eps_operon_map.png'))
    plt.close()
    print('  [7/9] EPS operon map saved')

# ─────────────────────────────────────────────
# Plot 8: PGPR Trait Radar Chart
# ─────────────────────────────────────────────
def plot_pgpr_radar(pgpr_df, outdir):
    if pgpr_df.empty:
        return
    traits = pgpr_df['trait'].tolist()
    genes_found = pgpr_df['genes_found'].tolist()
    genes_total = pgpr_df['genes_total'].tolist()
    positive = pgpr_df['positive'].tolist()
    scores = [f / max(t, 1) for f, t in zip(genes_found, genes_total)]
    n = len(traits)

    angles = np.linspace(0, 2 * np.pi, n, endpoint=False).tolist()
    angles += angles[:1]
    scores_plot = scores + scores[:1]

    short_names = [t.replace('_', '\n').title() for t in traits]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6),
                                    subplot_kw=dict(polar=True) if False else {})
    fig.suptitle('PGPR Trait Identification', fontsize=15, fontweight='bold')

    # Radar
    ax1 = fig.add_subplot(121, polar=True)
    ax1.set_theta_offset(np.pi / 2)
    ax1.set_theta_direction(-1)
    ax1.set_xticks(angles[:-1])
    ax1.set_xticklabels(short_names, size=8)
    ax1.set_ylim(0, 1)
    ax1.set_yticks([0.25, 0.5, 0.75, 1.0])
    ax1.set_yticklabels(['25%', '50%', '75%', '100%'], size=7)
    ax1.plot(angles, scores_plot, 'o-', linewidth=2, color='#2E86AB')
    ax1.fill(angles, scores_plot, alpha=0.25, color='#2E86AB')
    ax1.set_title('Gene Detection Rate per Trait', pad=15)

    # Bar chart
    ax2 = fig.add_subplot(122)
    colors = ['#2E86AB' if p == 'YES' else '#E0E0E0' for p in positive]
    bars = ax2.bar(short_names, scores, color=colors, edgecolor='white', linewidth=1.5)
    ax2.set_ylim(0, 1.2)
    ax2.set_ylabel('Gene Detection Rate')
    ax2.set_title('PGPR Trait Scores')
    ax2.axhline(0.5, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    for bar, g, t, p in zip(bars, genes_found, genes_total, positive):
        label = f'{g}/{t}\n{"✓" if p == "YES" else "✗"}'
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.03,
                 label, ha='center', va='bottom', fontsize=9,
                 color='#2E86AB' if p == 'YES' else '#888888', fontweight='bold')
    ax2.tick_params(axis='x', labelsize=8)

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, '08_pgpr_traits.png'))
    plt.close()
    print('  [8/9] PGPR radar chart saved')

# ─────────────────────────────────────────────
# Plot 9: Resistance Risk + Physicochemical Summary Dashboard
# ─────────────────────────────────────────────
def plot_summary_dashboard(report, outdir):
    resist = report.get('resistance_risk', {})
    physchem = report.get('physicochemical_summary', {})
    assembly = report.get('assembly', {})
    amp_disc = report.get('amp_discovery', {})

    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.4)
    fig.suptitle('AMP-MineFlow Pipeline Summary Dashboard', fontsize=16, fontweight='bold')

    # Panel 1: Resistance risk donut
    ax1 = fig.add_subplot(gs[0, 0])
    if resist:
        risk_order = ['very_low', 'low', 'moderate', 'high']
        risk_labels = ['Very Low', 'Low', 'Moderate', 'High']
        risk_colors = ['#44BBA4', '#2E86AB', '#F18F01', '#C73E1D']
        risk_vals = [resist.get(r, 0) for r in risk_order]
        total_r = sum(risk_vals)
        wedges, _ = ax1.pie(risk_vals, colors=risk_colors, startangle=90,
                             wedgeprops=dict(width=0.5, linewidth=2, edgecolor='white'))
        ax1.text(0, 0, f'{total_r:,}\nAMPs', ha='center', va='center',
                 fontsize=10, fontweight='bold')
        ax1.legend(wedges, [f'{l}: {v:,}' for l, v, in zip(risk_labels, risk_vals)],
                   loc='lower center', bbox_to_anchor=(0.5, -0.3), fontsize=8, ncol=2)
    ax1.set_title('Resistance Risk Profile')

    # Panel 2: Assembly stats
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.axis('off')
    asm_data = [
        ('Contigs', f"{assembly.get('num_contigs', 'NA')}"),
        ('Genome Size', f"{assembly.get('total_length_bp', 0):,} bp"),
        ('N50', f"{assembly.get('n50_bp', 0):,} bp"),
        ('GC Content', f"{assembly.get('gc_percent', 'NA')}%"),
        ('Mean Coverage', f"{assembly.get('mean_coverage', 'NA')}×"),
        ('AMP Candidates', f"{amp_disc.get('total_candidates', 0):,}"),
        ('AMP Families', f"{len(amp_disc.get('family_distribution', {}))}"),
        ('AMP Potential', amp_disc.get('amp_potential', 'NA')),
    ]
    table = ax2.table(cellText=[[k, v] for k, v in asm_data],
                       colLabels=['Metric', 'Value'],
                       cellLoc='left', loc='center',
                       bbox=[0, -0.05, 1, 1.1])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    for (r, c), cell in table.get_celld().items():
        if r == 0:
            cell.set_facecolor('#2E86AB')
            cell.set_text_props(color='white', fontweight='bold')
        elif r % 2 == 0:
            cell.set_facecolor('#F0F4F8')
        cell.set_edgecolor('#CCCCCC')
    ax2.set_title('Assembly & Discovery Summary', pad=20)

    # Panel 3: Physicochemical means
    ax3 = fig.add_subplot(gs[0, 2])
    if physchem:
        feat_labels = {
            'net_charge_pH7': 'Net Charge\n(pH 7)',
            'hydrophobic_ratio': 'Hydrophobic\nRatio',
            'amphipathic_moment': 'Amphipathic\nMoment',
        }
        x = [feat_labels[k] for k in feat_labels if k in physchem]
        means = [physchem[k]['mean'] for k in feat_labels if k in physchem]
        stds = [physchem[k]['std'] for k in feat_labels if k in physchem]
        colors_p = ['#2E86AB', '#44BBA4', '#F18F01'][:len(x)]
        bars3 = ax3.bar(x, means, yerr=stds, color=colors_p, edgecolor='white',
                        capsize=5, error_kw={'linewidth': 1.5})
        for bar, m in zip(bars3, means):
            ax3.text(bar.get_x() + bar.get_width()/2,
                     bar.get_height() + max(abs(v) for v in means) * 0.05,
                     f'{m:.3f}', ha='center', fontsize=9)
        ax3.set_ylabel('Mean Value')
        ax3.set_title('Key Physicochemical Properties')
        ax3.axhline(0, color='black', linewidth=0.8)

    # Panel 4: MOA donut
    ax4 = fig.add_subplot(gs[1, 0])
    moa_dist = report.get('mechanism_of_action', {}).get('distribution', {})
    if moa_dist:
        moa_labels_short = [k.replace('membrane_disruption_', 'MD-').replace('_', ' ').title()
                             for k in moa_dist.keys()]
        moa_vals = list(moa_dist.values())
        moa_colors = PALETTE[:len(moa_vals)]
        ax4.pie(moa_vals, colors=moa_colors, startangle=90,
                wedgeprops=dict(width=0.55, linewidth=1.5, edgecolor='white'))
        ax4.legend([f'{l}: {v:,}' for l, v in zip(moa_labels_short, moa_vals)],
                   loc='lower center', bbox_to_anchor=(0.5, -0.45),
                   fontsize=7, ncol=1)
    ax4.set_title('MOA Distribution')

    # Panel 5: Pathogen activity bar
    ax5 = fig.add_subplot(gs[1, 1:])
    spec_data = report.get('pathogen_spectrum', {})
    if spec_data:
        paths = [p.replace('_', ' ') for p in spec_data.keys()]
        rates = [spec_data[p]['rate'] for p in spec_data.keys()]
        bar_colors = ['#C73E1D' if r >= 0.5 else '#F18F01' if r >= 0.3 else '#2E86AB'
                      for r in rates]
        ax5.bar(paths, rates, color=bar_colors, edgecolor='white', linewidth=1)
        ax5.axhline(0.5, color='red', linestyle='--', linewidth=1, alpha=0.6, label='50%')
        ax5.axhline(0.3, color='gray', linestyle='--', linewidth=1, alpha=0.6, label='30%')
        ax5.set_ylim(0, 1.1)
        ax5.set_ylabel('Activity Rate')
        ax5.set_title('Pathogen Activity Rates')
        ax5.tick_params(axis='x', rotation=35, labelsize=8)
        ax5.legend(fontsize=8)

    plt.savefig(os.path.join(outdir, '09_pipeline_summary_dashboard.png'))
    plt.close()
    print('  [9/9] Summary dashboard saved')


def main():
    parser = argparse.ArgumentParser(description='Generate AMP-MineFlow plots')
    parser.add_argument('--results', default='results', help='Results directory')
    parser.add_argument('--outdir', default=None, help='Output directory for plots')
    args = parser.parse_args()

    results = args.results
    outdir = args.outdir or os.path.join(results, '14_report', 'plots')
    os.makedirs(outdir, exist_ok=True)

    print(f'Generating plots from: {results}')
    print(f'Saving plots to: {outdir}')

    # Load data
    report = load_json(os.path.join(results, '14_report', 'pipeline_summary.json'))
    fam_summary = load_json(os.path.join(results, '04_classification', 'family_summary.json'))
    cazyme_summary = load_json(os.path.join(results, '10_cazyme_annotation', 'cazyme_summary.json'))
    moa_summary = load_json(os.path.join(results, '07_moa_prediction', 'moa_summary.json'))
    clusters_df = load_tsv(os.path.join(results, '06_chemical_space', 'chemical_space_clusters.tsv'))
    eps_df = load_tsv(os.path.join(results, '11_eps_pathway', 'eps_operon_reconstruction.tsv'))
    pgpr_df = load_tsv(os.path.join(results, '12_pgpr_screening', 'pgpr_traits.tsv'))

    plot_amp_families(fam_summary, outdir)
    plot_score_distribution(report, outdir)
    plot_chemical_space(clusters_df, outdir)
    plot_moa(moa_summary, outdir)
    plot_pathogen_spectrum(report.get('pathogen_spectrum', {}), outdir)
    plot_cazymes(cazyme_summary, outdir)
    plot_eps_operon(eps_df, outdir)
    plot_pgpr_radar(pgpr_df, outdir)
    plot_summary_dashboard(report, outdir)

    print(f'\nAll 9 plots saved to: {outdir}')

if __name__ == '__main__':
    main()
