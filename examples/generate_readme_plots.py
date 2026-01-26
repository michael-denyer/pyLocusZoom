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

# Ensure lead SNP exists at exact position
peak_center = 1_500_000
positions[250] = peak_center  # Place lead SNP in middle of array

# Create a peak around position 1,500,000
p_values = np.ones(n_snps) * 0.5
for i, pos in enumerate(positions):
    dist = abs(pos - peak_center)
    if dist < 100_000:
        p_values[i] = 10 ** -(8 * np.exp(-dist / 30_000))
    else:
        p_values[i] = np.random.uniform(0.01, 1)

# Generate synthetic LD values (R² with lead SNP at peak_center)
# LD decays with distance from lead SNP
ld_values = []
for pos in positions:
    dist = abs(pos - peak_center)
    if dist == 0:
        r2 = 1.0
    elif dist < 50_000:
        # High LD close to lead
        r2 = max(0, 0.9 * np.exp(-dist / 20_000) + np.random.uniform(-0.1, 0.1))
    elif dist < 150_000:
        # Moderate LD at medium distance
        r2 = max(0, 0.5 * np.exp(-dist / 50_000) + np.random.uniform(-0.1, 0.1))
    else:
        # Low/no LD far from lead
        r2 = max(0, np.random.uniform(0, 0.2))
    ld_values.append(min(1.0, r2))

gwas_df = pd.DataFrame(
    {
        "ps": positions,
        "p_wald": p_values,
        "rs": [f"rs{i}" for i in range(n_snps)],
        "ld_r2": ld_values,  # Pre-computed LD column
    }
)

# Create gene annotations
genes_df = pd.DataFrame(
    {
        "chr": ["1", "1", "1", "1"],
        "start": [1_100_000, 1_400_000, 1_550_000, 1_800_000],
        "end": [1_200_000, 1_520_000, 1_650_000, 1_900_000],
        "gene_name": ["GENE1", "GENE2", "GENE3", "GENE4"],
        "strand": ["+", "-", "+", "-"],
    }
)

# Create exon annotations
exons_df = pd.DataFrame(
    {
        "chr": ["1", "1", "1", "1", "1", "1"],
        "start": [1_100_000, 1_150_000, 1_400_000, 1_450_000, 1_550_000, 1_600_000],
        "end": [1_120_000, 1_170_000, 1_420_000, 1_470_000, 1_580_000, 1_630_000],
        "gene_name": ["GENE1", "GENE1", "GENE2", "GENE2", "GENE3", "GENE3"],
    }
)

print("Generating example plots...")

# 1. Basic matplotlib plot with LD coloring
print("1. Basic regional plot with LD coloring...")
plotter = LocusZoomPlotter(species="canine", log_level=None)
fig = plotter.plot(
    gwas_df,
    chrom=1,
    start=1_000_000,
    end=2_000_000,
    lead_pos=1_500_000,
    ld_col="ld_r2",  # Use pre-computed LD values for coloring
    genes_df=genes_df,
    exons_df=exons_df,
    show_recombination=False,
    snp_labels=True,
    label_top_n=1,
)
fig.savefig("examples/regional_plot.png", dpi=150, bbox_inches="tight")
print("   Saved: examples/regional_plot.png")

# 2. Stacked plot with LD coloring
print("2. Stacked plot with LD coloring...")
gwas_df2 = gwas_df.copy()
gwas_df2["p_wald"] = np.ones(n_snps) * 0.5
peak_center2 = 1_700_000
# Ensure lead SNP exists at exact position for second panel
positions[350] = peak_center2
gwas_df2.loc[350, "ps"] = peak_center2
for i, pos in enumerate(positions):
    dist = abs(pos - peak_center2)
    if dist < 80_000:
        gwas_df2.loc[i, "p_wald"] = 10 ** -(6 * np.exp(-dist / 25_000))
    else:
        gwas_df2.loc[i, "p_wald"] = np.random.uniform(0.05, 1)

# Generate LD for second GWAS (different lead SNP)
ld_values2 = []
for pos in positions:
    dist = abs(pos - peak_center2)
    if dist == 0:
        r2 = 1.0
    elif dist < 50_000:
        r2 = max(0, 0.9 * np.exp(-dist / 20_000) + np.random.uniform(-0.1, 0.1))
    elif dist < 150_000:
        r2 = max(0, 0.5 * np.exp(-dist / 50_000) + np.random.uniform(-0.1, 0.1))
    else:
        r2 = max(0, np.random.uniform(0, 0.2))
    ld_values2.append(min(1.0, r2))
gwas_df2["ld_r2"] = ld_values2

fig = plotter.plot_stacked(
    [gwas_df, gwas_df2],
    chrom=1,
    start=1_000_000,
    end=2_000_000,
    lead_positions=[1_500_000, 1_700_000],  # Lead SNPs for each panel
    ld_col="ld_r2",  # Use pre-computed LD values for coloring
    panel_labels=["Phenotype A", "Phenotype B"],
    genes_df=genes_df,
    show_recombination=False,
    label_top_n=1,
)
fig.savefig("examples/stacked_plot.png", dpi=150, bbox_inches="tight")
print("   Saved: examples/stacked_plot.png")

# 3. eQTL overlay with effect sizes
print("3. eQTL overlay plot...")
eqtl_df = pd.DataFrame(
    {
        "pos": [1_420_000, 1_450_000, 1_480_000, 1_500_000, 1_520_000, 1_550_000, 1_600_000, 1_650_000],
        "p_value": [1e-4, 1e-5, 1e-6, 1e-7, 1e-6, 1e-4, 1e-3, 0.01],
        "gene": ["GENE2", "GENE2", "GENE2", "GENE2", "GENE2", "GENE2", "GENE3", "GENE3"],
        "effect_size": [0.35, 0.28, 0.22, 0.15, -0.18, -0.25, -0.32, 0.12],
    }
)

fig = plotter.plot_stacked(
    [gwas_df],
    chrom=1,
    start=1_000_000,
    end=2_000_000,
    lead_positions=[1_500_000],  # Lead SNP for LD coloring
    ld_col="ld_r2",  # Use pre-computed LD values for coloring
    eqtl_df=eqtl_df,
    eqtl_gene="GENE2",
    genes_df=genes_df,
    show_recombination=False,
    label_top_n=1,
)
fig.savefig("examples/eqtl_overlay.png", dpi=150, bbox_inches="tight")
print("   Saved: examples/eqtl_overlay.png")

# 4. Fine-mapping/SuSiE plot
print("4. Fine-mapping/SuSiE plot with credible sets...")

# Generate synthetic fine-mapping data
# Create PIP values that peak around the lead SNP
finemapping_positions = positions.copy()
pip_values = []
cs_assignments = []

for i, pos in enumerate(finemapping_positions):
    dist = abs(pos - peak_center)
    if dist < 20_000:
        # High PIP near causal variant
        pip = max(0, 0.95 * np.exp(-dist / 8_000) + np.random.uniform(-0.05, 0.05))
        cs = 1 if pip > 0.1 else 0
    elif dist < 80_000:
        # Moderate PIP in LD region
        pip = max(0, 0.3 * np.exp(-dist / 30_000) + np.random.uniform(-0.02, 0.02))
        cs = 1 if pip > 0.05 else 0
    else:
        # Low PIP elsewhere
        pip = max(0, np.random.uniform(0, 0.02))
        cs = 0
    pip_values.append(min(1.0, pip))
    cs_assignments.append(cs)

# Add a second credible set near a different peak
for i, pos in enumerate(finemapping_positions):
    dist = abs(pos - 1_300_000)  # Second signal
    if dist < 30_000:
        pip_values[i] = max(pip_values[i], 0.7 * np.exp(-dist / 12_000))
        if pip_values[i] > 0.05:
            cs_assignments[i] = 2  # Second credible set

finemapping_df = pd.DataFrame(
    {
        "pos": finemapping_positions,
        "pip": pip_values,
        "cs": cs_assignments,
        "rs": [f"rs{i}" for i in range(n_snps)],
    }
)

fig = plotter.plot_stacked(
    [gwas_df],
    chrom=1,
    start=1_000_000,
    end=2_000_000,
    lead_positions=[1_500_000],
    ld_col="ld_r2",
    finemapping_df=finemapping_df,
    finemapping_cs_col="cs",
    genes_df=genes_df,
    show_recombination=False,
    label_top_n=1,
)
fig.savefig("examples/finemapping_plot.png", dpi=150, bbox_inches="tight")
print("   Saved: examples/finemapping_plot.png")

print("\nAll plots generated successfully!")
