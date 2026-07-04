import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import re
import argparse
import os
from matplotlib.lines import Line2D


parser = argparse.ArgumentParser(description="Generate Manhattan plot from HyPhy results.")
parser.add_argument("--stats", required=True, help="Path to the HYPHY_selection_results.tsv file")
parser.add_argument("--gff", required=True, help="Path to the sorted GFF file")
parser.add_argument("--gea", default="no", help="Path to the GEA file, or 'no'")
args = parser.parse_args()

HAS_GEA = args.gea.lower() != "no"

if HAS_GEA:
    usecols_list = [0, 2, 5, 11, 12]
    names_list = ['gene', 'signif', 'proportion', 'pval', 'gea']
else:
    usecols_list = [0, 2, 5, 11] 
    names_list = ['gene', 'signif', 'proportion', 'pval']


df_stats = pd.read_csv(args.stats, 
                       sep=r'\s+',
                       header=0,
                       usecols=usecols_list,
                       names=names_list,
                       na_values=['NA', ''],
                       keep_default_na=True)

if not HAS_GEA:
    df_stats['gea'] = 0

df_stats['proportion'] = df_stats['proportion'].astype(str).str.replace('%', '').astype(float) / 100.0

# GFF
gff_cols = ['contig', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
df_gff = pd.read_csv(args.gff, 
                     sep='\t', comment='#', header=None, names=gff_cols)

output_dir = os.path.dirname(args.stats)


########################################

# calculate the log p-val, max at 25 
raw_log_pval = -np.log10(df_stats['pval'] + 1e-300)
df_stats['plot_pval'] = np.clip(raw_log_pval, a_min=None, a_max=25)


df_gff = df_gff[df_gff['type'] == 'CDS'].copy()


def extract_gene_id(attr_string):
    match = re.search(r'Parent=([^;]+)', attr_string)
    if match:
        raw_id = match.group(1)
        return re.sub(r'^(gene-|transcript:|mRNA-)', '', raw_id)
    return None

df_gff['gene'] = df_gff['attributes'].apply(extract_gene_id)

df_gff = df_gff.groupby(['contig', 'gene'], as_index=False).agg({
    'start': 'min',
    'end': 'max'
})


# merge & coords
df = pd.merge(df_stats, df_gff, on='gene', how='inner')

df['contig_num'] = df['contig'].str.split('_').str[-1].astype(int)


df['midpoint'] = (df['start'] + df['end']) / 2
df = df.sort_values(['contig_num', 'start'])

contig_lengths = df.groupby('contig_num')['midpoint'].max()
contig_offsets = contig_lengths.cumsum().shift(1).fillna(0)

df['abs_coord'] = df['midpoint'] + df['contig_num'].map(contig_offsets)
contig_centers = df.groupby('contig_num')['abs_coord'].mean()


# thresholds
alpha = 0.05
n_tests = len(df)

# Bonferroni
bonferroni_threshold = alpha / n_tests
log_bonferroni = -np.log10(bonferroni_threshold)

# FDR
pvals_sorted = df['pval'].sort_values()
ranks = np.arange(1, n_tests + 1)
bh_critical_values = (ranks / n_tests) * alpha

passing_pvals = pvals_sorted[pvals_sorted <= bh_critical_values]

if not passing_pvals.empty:
    fdr_threshold = passing_pvals.max()
    log_fdr = -np.log10(fdr_threshold)
else:
    fdr_threshold = None
    log_fdr = None



# plot config
size_circ = 5   
size_tri = 12    
size_tri_sig = 120
cmap = mcolors.LinearSegmentedColormap.from_list("blue_red", ["blue", "red"])

plot_configs = [
    {
        'file_suffix': 'log_pval_0.05',
        'y_column': 'plot_pval', 
        'y_label': '-log10(p-val)', 
        'title': 'Manhattan plot of genes under selection (-log10, p <= 0.05)',
        'pval_max': 0.05, 
        'bonf_line': log_bonferroni,
        'fdr_line': log_fdr
    },
]


# plot & save
for config in plot_configs:
    
    df_plot = df[df['pval'] <= config['pval_max']].copy()
    
    if df_plot.empty:
        print(f"No data to plot for threshold p <= {config['pval_max']}. Skipping.")
        continue

    df_nonsig = df_plot[df_plot['signif'] == 0]
    df_sig = df_plot[df_plot['signif'] == 1]

    fig, ax = plt.subplots(figsize=(12, 6))

    fig.patch.set_facecolor('white') 
    ax.set_facecolor('white')
    
    y_col = config['y_column']
    
    nonsig_no_gea = df_nonsig[df_nonsig['gea'] == 0]
    nonsig_gea = df_nonsig[df_nonsig['gea'] == 1]

    sig_no_gea = df_sig[df_sig['gea'] == 0]
    sig_gea = df_sig[df_sig['gea'] == 1]

    ax.scatter(nonsig_no_gea['abs_coord'], nonsig_no_gea[y_col], c='grey', marker='o', s=size_circ, alpha=0.4, edgecolors='none')
    ax.scatter(nonsig_gea['abs_coord'], nonsig_gea[y_col], c='grey', marker='^', s=size_tri, alpha=0.4, edgecolors='none')

    sc_circ = ax.scatter(sig_no_gea['abs_coord'], sig_no_gea[y_col], c=sig_no_gea['proportion'], cmap=cmap, vmin=0, vmax=1, marker='o', s=size_circ, alpha=0.9, edgecolors='none')
    sc_tri = ax.scatter(sig_gea['abs_coord'], sig_gea[y_col], c=sig_gea['proportion'], cmap=cmap, vmin=0, vmax=1, marker='^', s=size_tri_sig, alpha=0.9, edgecolors='none')

    ax.axhline(y=config['bonf_line'], color='black', linestyle='-', linewidth=1.5, 
               label=f"Bonferroni (p={bonferroni_threshold:.2e})")

    if config['fdr_line'] is not None:
        ax.axhline(y=config['fdr_line'], color='black', linestyle='--', linewidth=1.5, 
                   label=f"FDR Cutoff (p={fdr_threshold:.2e})")

    # contig shades
    for i in range(len(contig_lengths)):
        if i % 2 == 0:
            start_pos = contig_offsets.iloc[i]
            end_pos = start_pos + contig_lengths.iloc[i]
            ax.axvspan(start_pos, end_pos, color='grey', alpha=0.15, zorder=0)
    
    # x-ticks spacing logic
    tick_positions = []
    tick_labels = []
    min_dist = (df_plot['abs_coord'].max() - df_plot['abs_coord'].min()) * 0.03
    last_pos = -min_dist 

    for contig_num, center_pos in contig_centers.items():
        if (center_pos - last_pos) >= min_dist:
            tick_positions.append(center_pos)
            tick_labels.append(str(contig_num))
            last_pos = center_pos

    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, fontsize=7, rotation=45, ha='right')

    # X-axis padding
    x_min, x_max = df_plot['abs_coord'].min(), df_plot['abs_coord'].max()
    padding = (x_max - x_min) * 0.02
    ax.set_xlim(x_min - padding, x_max + padding)
    
    # Y-axis limit
    ax.set_ylim(bottom=-1)


    # adjusted legend
    legend_elements = []
    if HAS_GEA:
        legend_elements.append(Line2D([0], [0], marker='o', color='white', markerfacecolor='black', markersize=5, label='Non-GEA'))
        legend_elements.append(Line2D([0], [0], marker='^', color='white', markerfacecolor='black', markersize=7, label='GEA'))
    else:
        legend_elements.append(Line2D([0], [0], marker='o', color='white', markerfacecolor='black', markersize=5, label='Genes'))

    legend_elements.append(Line2D([0], [0], color='black', linestyle='-', linewidth=1.5, label='Bonferroni'))
    
    if config['fdr_line'] is not None:
      legend_elements.append(Line2D([0], [0], color='black', linestyle='--', linewidth=1.5, label='FDR'))
    
    legend = ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(0.99, 0.85), facecolor='white', edgecolor='black')
    for text in legend.get_texts():
        text.set_color("black")
        
    # colorbar
    fig.subplots_adjust(right=0.82) 
    cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7]) 
    cbar = fig.colorbar(sc_circ, cax=cbar_ax)
    cbar.ax.tick_params(colors='black', labelcolor='black', labelsize=9)
    cbar.set_label('Detected proportion of gene', color='black', fontsize=10, labelpad=15)

    # save
    plt.savefig(f"{output_dir}/manhattan_plot_positive_selection_genes_{config['file_suffix']}.png", dpi=300, bbox_inches='tight', facecolor=fig.get_facecolor())
    print(f"Saved plot : {output_dir}")
    
    