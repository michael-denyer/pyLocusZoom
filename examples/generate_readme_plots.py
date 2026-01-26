#!/usr/bin/env python3
"""Generate example plots for README documentation.

Note: Backend integration (plotly/bokeh) is not yet fully implemented in the main
plot methods. This script generates matplotlib plots only.
"""

import matplotlib
matplotlib.use("Agg")  # Non-interactive backend

import numpy as np
import pandas as pd

from pylocuszoom import LocusZoomPlotter

# Generate synthetic GWAS data
np.random.seed(42)
n_snps = 500
positions = np.sort(np.random.randint(1_000_000, 2_000_000, n_snps))

# Create a peak around position 1,500,000
p_values = np.ones(n_snps) * 0.5
peak_center = 1_500_000
for i, pos in enumerate(positions):
    dist = abs(pos - peak_center)
    if dist < 100_000:
        p_values[i] = 10 ** -(8 * np.exp(-dist / 30_000))
    else:
        p_values[i] = np.random.uniform(0.01, 1)

gwas_df = pd.DataFrame({
    "ps": positions,
    "p_wald": p_values,
    "rs": [f"rs{i}" for i in range(n_snps)],
})

# Create gene annotations
genes_df = pd.DataFrame({
    "chr": ["1", "1", "1", "1"],
    "start": [1_100_000, 1_400_000, 1_550_000, 1_800_000],
    "end": [1_200_000, 1_520_000, 1_650_000, 1_900_000],
    "gene_name": ["GENE1", "GENE2", "GENE3", "GENE4"],
    "strand": ["+", "-", "+", "-"],
})

# Create exon annotations
exons_df = pd.DataFrame({
    "chr": ["1", "1", "1", "1", "1", "1"],
    "start": [1_100_000, 1_150_000, 1_400_000, 1_450_000, 1_550_000, 1_600_000],
    "end": [1_120_000, 1_170_000, 1_420_000, 1_470_000, 1_580_000, 1_630_000],
    "gene_name": ["GENE1", "GENE1", "GENE2", "GENE2", "GENE3", "GENE3"],
})

print("Generating example plots...")

# 1. Basic matplotlib plot
print("1. Basic regional plot...")
plotter = LocusZoomPlotter(species="dog", log_level=None)
fig = plotter.plot(
    gwas_df,
    chrom=1,
    start=1_000_000,
    end=2_000_000,
    lead_pos=1_500_000,
    genes_df=genes_df,
    exons_df=exons_df,
    show_recombination=False,
    snp_labels=True,
    label_top_n=3,
)
fig.savefig("examples/regional_plot.png", dpi=150, bbox_inches="tight")
print("   Saved: examples/regional_plot.png")

# 2. Stacked plot
print("2. Stacked plot...")
gwas_df2 = gwas_df.copy()
gwas_df2["p_wald"] = np.ones(n_snps) * 0.5
peak_center2 = 1_700_000
for i, pos in enumerate(positions):
    dist = abs(pos - peak_center2)
    if dist < 80_000:
        gwas_df2.loc[i, "p_wald"] = 10 ** -(6 * np.exp(-dist / 25_000))
    else:
        gwas_df2.loc[i, "p_wald"] = np.random.uniform(0.05, 1)

fig = plotter.plot_stacked(
    [gwas_df, gwas_df2],
    chrom=1,
    start=1_000_000,
    end=2_000_000,
    panel_labels=["Phenotype A", "Phenotype B"],
    genes_df=genes_df,
    show_recombination=False,
    label_top_n=2,
)
fig.savefig("examples/stacked_plot.png", dpi=150, bbox_inches="tight")
print("   Saved: examples/stacked_plot.png")

# 3. eQTL overlay
print("3. eQTL overlay plot...")
eqtl_df = pd.DataFrame({
    "pos": [1_480_000, 1_500_000, 1_520_000, 1_550_000, 1_600_000],
    "p_value": [1e-5, 1e-7, 1e-6, 1e-4, 0.01],
    "gene": ["GENE2", "GENE2", "GENE2", "GENE3", "GENE3"],
})

fig = plotter.plot_stacked(
    [gwas_df],
    chrom=1,
    start=1_000_000,
    end=2_000_000,
    eqtl_df=eqtl_df,
    eqtl_gene="GENE2",
    genes_df=genes_df,
    show_recombination=False,
)
fig.savefig("examples/eqtl_overlay.png", dpi=150, bbox_inches="tight")
print("   Saved: examples/eqtl_overlay.png")

print("\nAll plots generated successfully!")
