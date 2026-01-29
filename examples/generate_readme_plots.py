#!/usr/bin/env python3
"""Generate example plots for README documentation.

Generates:
- Static PNG plots (matplotlib) for README display
- Interactive HTML plots (plotly/bokeh) for exploration
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

# Create gene annotations - realistic overlapping genes for multi-track display
genes_df = pd.DataFrame(
    {
        "chr": ["1", "1", "1", "1", "1", "1"],
        "start": [
            1_050_000,  # ABCB1 - long gene
            1_250_000,  # ADAM7 - overlaps with ABCB1
            1_400_000,  # SLC25A - near peak
            1_520_000,  # PDGFA - overlaps SLC25A, near peak
            1_700_000,  # TP53 - downstream
            1_850_000,  # WRAP73 - far downstream
        ],
        "end": [
            1_280_000,  # ABCB1
            1_420_000,  # ADAM7
            1_560_000,  # SLC25A
            1_680_000,  # PDGFA
            1_820_000,  # TP53
            1_960_000,  # WRAP73
        ],
        "gene_name": ["ABCB1", "ADAM7", "SLC25A", "PDGFA", "TP53", "WRAP73"],
        "strand": ["+", "-", "+", "-", "+", "-"],
    }
)

# Create exon annotations - multiple exons per gene
exons_df = pd.DataFrame(
    {
        "chr": ["1"] * 22,
        "start": [
            # ABCB1 - 5 exons
            1_055_000,
            1_100_000,
            1_160_000,
            1_220_000,
            1_260_000,
            # ADAM7 - 4 exons
            1_255_000,
            1_300_000,
            1_360_000,
            1_400_000,
            # SLC25A - 4 exons
            1_405_000,
            1_450_000,
            1_500_000,
            1_540_000,
            # PDGFA - 3 exons
            1_525_000,
            1_590_000,
            1_650_000,
            # TP53 - 4 exons
            1_705_000,
            1_740_000,
            1_780_000,
            1_800_000,
            # WRAP73 - 2 exons
            1_855_000,
            1_920_000,
        ],
        "end": [
            # ABCB1 exons
            1_075_000,
            1_125_000,
            1_185_000,
            1_245_000,
            1_280_000,
            # ADAM7 exons
            1_280_000,
            1_330_000,
            1_390_000,
            1_420_000,
            # SLC25A exons
            1_430_000,
            1_475_000,
            1_525_000,
            1_560_000,
            # PDGFA exons
            1_555_000,
            1_620_000,
            1_680_000,
            # TP53 exons
            1_725_000,
            1_765_000,
            1_800_000,
            1_820_000,
            # WRAP73 exons
            1_885_000,
            1_955_000,
        ],
        "gene_name": [
            "ABCB1",
            "ABCB1",
            "ABCB1",
            "ABCB1",
            "ABCB1",
            "ADAM7",
            "ADAM7",
            "ADAM7",
            "ADAM7",
            "SLC25A",
            "SLC25A",
            "SLC25A",
            "SLC25A",
            "PDGFA",
            "PDGFA",
            "PDGFA",
            "TP53",
            "TP53",
            "TP53",
            "TP53",
            "WRAP73",
            "WRAP73",
        ],
    }
)

print("Generating example plots...")

plotter = LocusZoomPlotter(species="canine", log_level=None)

# 1. Regional plot with recombination overlay (full LocusZoom style)
# Use a region that has recombination data (canine maps start at ~4Mb on chr1)
print("1. Regional plot with recombination overlay...")
recomb_positions = np.sort(np.random.randint(12_000_000, 14_000_000, n_snps))
recomb_peak_center = 13_000_000
recomb_positions[250] = recomb_peak_center

recomb_p_values = np.ones(n_snps) * 0.5
for i, pos in enumerate(recomb_positions):
    dist = abs(pos - recomb_peak_center)
    if dist < 100_000:
        recomb_p_values[i] = 10 ** -(8 * np.exp(-dist / 30_000))
    else:
        recomb_p_values[i] = np.random.uniform(0.01, 1)

recomb_ld_values = []
for pos in recomb_positions:
    dist = abs(pos - recomb_peak_center)
    if dist == 0:
        r2 = 1.0
    elif dist < 50_000:
        r2 = max(0, 0.9 * np.exp(-dist / 20_000) + np.random.uniform(-0.1, 0.1))
    elif dist < 150_000:
        r2 = max(0, 0.5 * np.exp(-dist / 50_000) + np.random.uniform(-0.1, 0.1))
    else:
        r2 = max(0, np.random.uniform(0, 0.2))
    recomb_ld_values.append(min(1.0, r2))

recomb_gwas_df = pd.DataFrame(
    {
        "ps": recomb_positions,
        "p_wald": recomb_p_values,
        "rs": [f"rs{i}" for i in range(n_snps)],
        "ld_r2": recomb_ld_values,
    }
)

# More realistic genes - overlapping to create multiple tracks
recomb_genes_df = pd.DataFrame(
    {
        "chr": ["1", "1", "1", "1", "1", "1"],
        "start": [
            12_050_000,  # BRCA1 - long gene spanning much of region
            12_400_000,  # NBR2 - overlaps with BRCA1 (different strand)
            12_700_000,  # RND2 - short gene
            12_900_000,  # CNTNAP1 - overlaps peak region
            13_200_000,  # EZH1 - near peak
            13_550_000,  # NAGLU - downstream
        ],
        "end": [
            12_550_000,  # BRCA1
            12_700_000,  # NBR2
            12_850_000,  # RND2
            13_250_000,  # CNTNAP1
            13_500_000,  # EZH1
            13_850_000,  # NAGLU
        ],
        "gene_name": ["BRCA1", "NBR2", "RND2", "CNTNAP1", "EZH1", "NAGLU"],
        "strand": ["+", "-", "+", "-", "+", "-"],
    }
)

# Realistic exon structure - exons are small coding regions within genes
recomb_exons_df = pd.DataFrame(
    {
        "chr": ["1"] * 26,
        "start": [
            # BRCA1 - 6 exons (tumor suppressor with many exons)
            12_050_000,
            12_120_000,
            12_200_000,
            12_300_000,
            12_420_000,
            12_510_000,
            # NBR2 - 4 exons
            12_400_000,
            12_500_000,
            12_590_000,
            12_670_000,
            # RND2 - 3 exons (small gene)
            12_705_000,
            12_770_000,
            12_830_000,
            # CNTNAP1 - 5 exons
            12_905_000,
            13_000_000,
            13_080_000,
            13_160_000,
            13_220_000,
            # EZH1 - 4 exons
            13_205_000,
            13_300_000,
            13_400_000,
            13_475_000,
            # NAGLU - 4 exons
            13_560_000,
            13_670_000,
            13_750_000,
            13_820_000,
        ],
        "end": [
            # BRCA1 exons - small boxes (10-25kb each)
            12_065_000,
            12_135_000,
            12_220_000,
            12_320_000,
            12_440_000,
            12_535_000,
            # NBR2 exons
            12_420_000,
            12_520_000,
            12_610_000,
            12_695_000,
            # RND2 exons
            12_720_000,
            12_790_000,
            12_850_000,
            # CNTNAP1 exons
            12_925_000,
            13_020_000,
            13_100_000,
            13_180_000,
            13_245_000,
            # EZH1 exons
            13_225_000,
            13_325_000,
            13_425_000,
            13_500_000,
            # NAGLU exons
            13_585_000,
            13_695_000,
            13_775_000,
            13_850_000,
        ],
        "gene_name": [
            "BRCA1",
            "BRCA1",
            "BRCA1",
            "BRCA1",
            "BRCA1",
            "BRCA1",
            "NBR2",
            "NBR2",
            "NBR2",
            "NBR2",
            "RND2",
            "RND2",
            "RND2",
            "CNTNAP1",
            "CNTNAP1",
            "CNTNAP1",
            "CNTNAP1",
            "CNTNAP1",
            "EZH1",
            "EZH1",
            "EZH1",
            "EZH1",
            "NAGLU",
            "NAGLU",
            "NAGLU",
            "NAGLU",
        ],
    }
)

fig = plotter.plot(
    recomb_gwas_df,
    chrom=1,
    start=12_000_000,
    end=14_000_000,
    lead_pos=13_000_000,
    ld_col="ld_r2",
    genes_df=recomb_genes_df,
    exons_df=recomb_exons_df,  # Show exon structure
    show_recombination=True,  # Enable recombination rate overlay
    snp_labels=True,
    label_top_n=1,
)
fig.savefig("examples/regional_plot_with_recomb.png", dpi=150, bbox_inches="tight")
print("   Saved: examples/regional_plot_with_recomb.png")

# 2. Basic matplotlib plot with LD coloring (no recombination)
print("2. Basic regional plot with LD coloring...")
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

# 3. Stacked plot with LD coloring
print("3. Stacked plot with LD coloring...")
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
    exons_df=exons_df,
    show_recombination=False,
    label_top_n=1,
)
fig.savefig("examples/stacked_plot.png", dpi=150, bbox_inches="tight")
print("   Saved: examples/stacked_plot.png")

# 4. eQTL overlay with effect sizes
# Realistic eQTL data for SLC25A gene (near the GWAS peak)
print("4. eQTL overlay plot...")
eqtl_df = pd.DataFrame(
    {
        "pos": [
            1_410_000,  # Upstream of SLC25A
            1_435_000,
            1_460_000,
            1_485_000,
            1_500_000,  # Near peak - strongest eQTL
            1_515_000,
            1_530_000,
            1_545_000,
            1_560_000,
            1_590_000,  # PDGFA region
            1_620_000,
            1_660_000,
        ],
        "p_value": [
            5e-4,
            8e-6,
            2e-8,
            5e-12,
            1e-15,  # Strong eQTL signal at peak
            3e-10,
            1e-6,
            5e-4,
            0.01,  # Declining signal
            1e-3,
            0.05,
            0.2,  # Weak signal in PDGFA
        ],
        "gene": [
            "SLC25A",
            "SLC25A",
            "SLC25A",
            "SLC25A",
            "SLC25A",
            "SLC25A",
            "SLC25A",
            "SLC25A",
            "SLC25A",
            "PDGFA",
            "PDGFA",
            "PDGFA",
        ],
        "effect_size": [
            0.18,
            0.25,
            0.38,
            0.52,
            0.65,  # Increasing effect toward peak
            0.48,
            0.32,
            0.15,
            0.08,  # Declining effect
            -0.22,
            -0.15,
            -0.08,  # Opposite direction for PDGFA
        ],
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
    eqtl_gene="SLC25A",
    genes_df=genes_df,
    exons_df=exons_df,
    show_recombination=False,
    label_top_n=1,
)
fig.savefig("examples/eqtl_overlay.png", dpi=150, bbox_inches="tight")
print("   Saved: examples/eqtl_overlay.png")

# 5. Fine-mapping/SuSiE plot
print("5. Fine-mapping/SuSiE plot with credible sets...")

# Generate realistic fine-mapping data mimicking SuSiE output
# Two independent signals: primary at peak_center (1.5Mb), secondary at 1.3Mb
finemapping_positions = positions.copy()
pip_values = np.zeros(n_snps)
cs_assignments = np.zeros(n_snps, dtype=int)

# Primary signal (CS1) - strong causal variant at peak_center
# SuSiE typically identifies a few high-PIP variants in tight LD
primary_causal_idx = 250  # Lead SNP position
pip_values[primary_causal_idx] = 0.89  # High confidence causal variant

# Add 3-4 variants in the 95% credible set for CS1
# These are in LD with the causal variant
cs1_members = []
for i, pos in enumerate(finemapping_positions):
    dist = abs(pos - peak_center)
    if i == primary_causal_idx:
        cs_assignments[i] = 1
        cs1_members.append(i)
    elif dist < 15_000 and len(cs1_members) < 5:
        # Nearby variants with moderate PIP (LD-driven)
        pip_values[i] = max(
            0.02, 0.15 * np.exp(-dist / 5_000) + np.random.uniform(0, 0.03)
        )
        if pip_values[i] > 0.01:
            cs_assignments[i] = 1
            cs1_members.append(i)
    elif dist < 50_000:
        # Background noise in LD region
        pip_values[i] = max(0, np.random.uniform(0, 0.01))

# Secondary signal (CS2) - independent signal at 1.3Mb (in ADAM7 gene region)
secondary_center = 1_300_000
# Find closest variant to secondary center
secondary_idx = np.argmin(np.abs(finemapping_positions - secondary_center))
pip_values[secondary_idx] = 0.72  # Moderately confident causal variant
cs_assignments[secondary_idx] = 2

# Add CS2 members
cs2_members = [secondary_idx]
for i, pos in enumerate(finemapping_positions):
    if i == secondary_idx:
        continue
    dist = abs(pos - secondary_center)
    if dist < 20_000 and len(cs2_members) < 4:
        pip = max(0.02, 0.12 * np.exp(-dist / 8_000) + np.random.uniform(0, 0.02))
        if pip > pip_values[i]:  # Don't overwrite CS1 members
            pip_values[i] = pip
            if cs_assignments[i] == 0:  # Don't reassign CS1 members
                cs_assignments[i] = 2
                cs2_members.append(i)
    elif dist < 60_000 and cs_assignments[i] == 0:
        # Background in CS2 region
        pip_values[i] = max(pip_values[i], np.random.uniform(0, 0.008))

# Add very low background PIP elsewhere (typical SuSiE output)
for i in range(n_snps):
    if pip_values[i] == 0:
        pip_values[i] = np.random.uniform(0, 0.002)

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
    exons_df=exons_df,
    show_recombination=False,
    label_top_n=1,
)
fig.savefig("examples/finemapping_plot.png", dpi=150, bbox_inches="tight")
print("   Saved: examples/finemapping_plot.png")

# 6. Interactive Plotly regional plot with recombination
print("6. Interactive Plotly regional plot with recombination...")
plotly_plotter = LocusZoomPlotter(species="canine", backend="plotly", log_level=None)
fig = plotly_plotter.plot(
    recomb_gwas_df,
    chrom=1,
    start=12_000_000,
    end=14_000_000,
    lead_pos=13_000_000,
    ld_col="ld_r2",
    genes_df=recomb_genes_df,
    exons_df=recomb_exons_df,
    show_recombination=True,
)
fig.write_html("examples/regional_recomb_plotly.html")
print("   Saved: examples/regional_recomb_plotly.html")

# 8. Interactive Plotly eQTL plot
print("8. Interactive Plotly eQTL plot...")
fig = plotly_plotter.plot_stacked(
    [gwas_df],
    chrom=1,
    start=1_000_000,
    end=2_000_000,
    lead_positions=[1_500_000],
    ld_col="ld_r2",
    eqtl_df=eqtl_df,
    eqtl_gene="SLC25A",
    genes_df=genes_df,
    exons_df=exons_df,
    show_recombination=False,
)
fig.write_html("examples/eqtl_plotly.html")
print("   Saved: examples/eqtl_plotly.html")

# 10. Interactive Plotly fine-mapping plot
print("10. Interactive Plotly fine-mapping plot...")
fig = plotly_plotter.plot_stacked(
    [gwas_df],
    chrom=1,
    start=1_000_000,
    end=2_000_000,
    lead_positions=[1_500_000],
    ld_col="ld_r2",
    finemapping_df=finemapping_df,
    finemapping_cs_col="cs",
    genes_df=genes_df,
    exons_df=exons_df,
    show_recombination=False,
)
fig.write_html("examples/finemapping_plotly.html")
print("   Saved: examples/finemapping_plotly.html")

# 11. Interactive Bokeh regional plot with recombination
print("11. Interactive Bokeh regional plot with recombination...")
from bokeh.io import output_file, save

bokeh_plotter = LocusZoomPlotter(species="canine", backend="bokeh", log_level=None)
fig = bokeh_plotter.plot(
    recomb_gwas_df,
    chrom=1,
    start=12_000_000,
    end=14_000_000,
    lead_pos=13_000_000,
    ld_col="ld_r2",
    genes_df=recomb_genes_df,
    exons_df=recomb_exons_df,
    show_recombination=True,
)
output_file("examples/regional_recomb_bokeh.html")
save(fig)
print("   Saved: examples/regional_recomb_bokeh.html")

# 13. Interactive Bokeh eQTL plot
print("13. Interactive Bokeh eQTL plot...")
fig = bokeh_plotter.plot_stacked(
    [gwas_df],
    chrom=1,
    start=1_000_000,
    end=2_000_000,
    lead_positions=[1_500_000],
    ld_col="ld_r2",
    eqtl_df=eqtl_df,
    eqtl_gene="SLC25A",
    genes_df=genes_df,
    exons_df=exons_df,
    show_recombination=False,
)
output_file("examples/eqtl_bokeh.html")
save(fig)
print("   Saved: examples/eqtl_bokeh.html")

# 15. Interactive Bokeh fine-mapping plot
print("15. Interactive Bokeh fine-mapping plot...")
fig = bokeh_plotter.plot_stacked(
    [gwas_df],
    chrom=1,
    start=1_000_000,
    end=2_000_000,
    lead_positions=[1_500_000],
    ld_col="ld_r2",
    finemapping_df=finemapping_df,
    finemapping_cs_col="cs",
    genes_df=genes_df,
    exons_df=exons_df,
    show_recombination=False,
)
output_file("examples/finemapping_bokeh.html")
save(fig)
print("   Saved: examples/finemapping_bokeh.html")

# 16. PheWAS plot
# Realistic PheWAS data showing pleiotropic effects across multiple trait categories
print("16. PheWAS plot...")
phewas_df = pd.DataFrame(
    {
        "phenotype": [
            # Anthropometric (5 traits)
            "Height",
            "BMI",
            "Weight",
            "Waist circumference",
            "Hip circumference",
            # Metabolic (6 traits)
            "Type 2 Diabetes",
            "Fasting Glucose",
            "HbA1c",
            "HOMA-IR",
            "Fasting Insulin",
            "Proinsulin",
            # Cardiovascular (5 traits)
            "Coronary Artery Disease",
            "Myocardial Infarction",
            "Stroke",
            "Atrial Fibrillation",
            "Heart Rate",
            # Lipids (4 traits)
            "HDL Cholesterol",
            "LDL Cholesterol",
            "Total Cholesterol",
            "Triglycerides",
            # Blood pressure (3 traits)
            "Systolic BP",
            "Diastolic BP",
            "Pulse Pressure",
            # Other (4 traits)
            "eGFR",
            "Urate",
            "CRP",
            "Vitamin D",
        ],
        "p_value": [
            # Anthropometric - moderate association with height
            1e-12,
            0.15,
            1e-6,
            0.08,
            0.25,
            # Metabolic - strong T2D signal (primary association)
            5e-45,
            2e-28,
            8e-22,
            1e-18,
            5e-15,
            3e-10,
            # Cardiovascular - secondary CAD signal
            2e-8,
            5e-6,
            0.02,
            0.35,
            0.18,
            # Lipids - modest HDL association
            1e-10,
            0.08,
            0.15,
            5e-4,
            # Blood pressure - weak signal
            0.005,
            0.02,
            0.08,
            # Other - no strong signals
            0.12,
            0.45,
            0.03,
            0.68,
        ],
        "category": [
            "Anthropometric",
            "Anthropometric",
            "Anthropometric",
            "Anthropometric",
            "Anthropometric",
            "Metabolic",
            "Metabolic",
            "Metabolic",
            "Metabolic",
            "Metabolic",
            "Metabolic",
            "Cardiovascular",
            "Cardiovascular",
            "Cardiovascular",
            "Cardiovascular",
            "Cardiovascular",
            "Lipids",
            "Lipids",
            "Lipids",
            "Lipids",
            "Blood Pressure",
            "Blood Pressure",
            "Blood Pressure",
            "Other",
            "Other",
            "Other",
            "Other",
        ],
    }
)
fig = plotter.plot_phewas(phewas_df, variant_id="rs7903146")
fig.savefig("examples/phewas_plot.png", dpi=150, bbox_inches="tight")
print("   Saved: examples/phewas_plot.png")

# 17. Forest plot (odds ratios with null at 1.0)
# Realistic meta-analysis forest plot for T2D association
print("17. Forest plot...")
forest_df = pd.DataFrame(
    {
        "study": [
            "DIAGRAM (EUR)",
            "AGEN (EAS)",
            "SIGMA (AMR)",
            "UK Biobank",
            "FinnGen",
            "PAGE (Multi-ethnic)",
            "Combined Meta-analysis",
        ],
        "effect": [1.35, 1.28, 1.42, 1.31, 1.38, 1.25, 1.33],
        "ci_lower": [1.28, 1.18, 1.22, 1.26, 1.25, 1.12, 1.29],
        "ci_upper": [1.43, 1.39, 1.65, 1.36, 1.52, 1.40, 1.37],
        "weight": [22, 18, 8, 28, 12, 12, 100],
    }
)
fig = plotter.plot_forest(
    forest_df,
    variant_id="rs7903146 (TCF7L2)",
    weight_col="weight",
    null_value=1.0,
    effect_label="Odds Ratio (T2D)",
)
fig.savefig("examples/forest_plot.png", dpi=150, bbox_inches="tight")
print("   Saved: examples/forest_plot.png")

# Manhattan plot - genome-wide view
print("18. Manhattan plot...")
np.random.seed(42)
manhattan_data = []
for chrom in range(1, 23):
    n_variants = np.random.randint(300, 600)
    chrom_positions = np.sort(np.random.randint(1e6, 2e8, n_variants))
    chrom_pvalues = np.random.uniform(0, 1, n_variants)
    # Add some significant hits
    if chrom in [6, 11, 17]:  # Chromosomes with GWAS hits
        n_hits = np.random.randint(5, 15)
        hit_indices = np.random.choice(n_variants, n_hits, replace=False)
        chrom_pvalues[hit_indices] = 10 ** np.random.uniform(-12, -8, n_hits)
    for i in range(n_variants):
        manhattan_data.append(
            {"chrom": str(chrom), "pos": chrom_positions[i], "p": chrom_pvalues[i]}
        )

manhattan_df = pd.DataFrame(manhattan_data)

manhattan_plotter = LocusZoomPlotter(species="human", log_level=None)
fig = manhattan_plotter.plot_manhattan(
    manhattan_df,
    significance_threshold=5e-8,
    figsize=(14, 4),
    title="Genome-wide Association Study",
)
fig.savefig("examples/manhattan_plot.png", dpi=150, bbox_inches="tight")
print("   Saved: examples/manhattan_plot.png")

# QQ plot
print("19. QQ plot...")
np.random.seed(42)
# Create inflated p-value distribution (lambda ~ 1.1)
n_pvalues = 5000
qq_pvalues = np.random.uniform(0, 1, n_pvalues)
# Add some true associations to create deviation at tail
n_true = 50
qq_pvalues[:n_true] = 10 ** np.random.uniform(-10, -5, n_true)

qq_df = pd.DataFrame({"p": qq_pvalues})

qq_plotter = LocusZoomPlotter(log_level=None)
fig = qq_plotter.plot_qq(
    qq_df,
    show_confidence_band=True,
    show_lambda=True,
    figsize=(5, 5),
)
fig.savefig("examples/qq_plot.png", dpi=150, bbox_inches="tight")
print("   Saved: examples/qq_plot.png")

# Interactive Plotly Manhattan plot
print("20. Interactive Plotly Manhattan plot...")
manhattan_plotter_plotly = LocusZoomPlotter(
    species="human", backend="plotly", log_level=None
)
fig = manhattan_plotter_plotly.plot_manhattan(
    manhattan_df,
    significance_threshold=5e-8,
    figsize=(14, 4),
    title="Genome-wide Association Study",
)
fig.write_html("examples/manhattan_plotly.html")
print("   Saved: examples/manhattan_plotly.html")

# Interactive Bokeh Manhattan plot
print("21. Interactive Bokeh Manhattan plot...")
from bokeh.io import save
from bokeh.resources import CDN

manhattan_plotter_bokeh = LocusZoomPlotter(
    species="human", backend="bokeh", log_level=None
)
fig = manhattan_plotter_bokeh.plot_manhattan(
    manhattan_df,
    significance_threshold=5e-8,
    figsize=(14, 4),
    title="Genome-wide Association Study",
)
save(
    fig, filename="examples/manhattan_bokeh.html", resources=CDN, title="Manhattan Plot"
)
print("   Saved: examples/manhattan_bokeh.html")

# Interactive Plotly QQ plot
print("22. Interactive Plotly QQ plot...")
qq_plotter_plotly = LocusZoomPlotter(backend="plotly", log_level=None)
fig = qq_plotter_plotly.plot_qq(
    qq_df,
    show_confidence_band=True,
    show_lambda=True,
    figsize=(5, 5),
)
fig.write_html("examples/qq_plotly.html")
print("   Saved: examples/qq_plotly.html")

# Interactive Bokeh QQ plot
print("23. Interactive Bokeh QQ plot...")
qq_plotter_bokeh = LocusZoomPlotter(backend="bokeh", log_level=None)
fig = qq_plotter_bokeh.plot_qq(
    qq_df,
    show_confidence_band=True,
    show_lambda=True,
    figsize=(5, 5),
)
save(fig, filename="examples/qq_bokeh.html", resources=CDN, title="QQ Plot")
print("   Saved: examples/qq_bokeh.html")

# Stacked Manhattan plot - multiple GWAS comparison
print("24. Stacked Manhattan plot...")
# Create 3 different GWAS datasets with different signal patterns
np.random.seed(123)
stacked_gwas_dfs = []
panel_names = ["Discovery Cohort", "Replication Cohort", "Meta-analysis"]

for i, name in enumerate(panel_names):
    gwas_data = []
    for chrom in range(1, 23):
        n_variants = np.random.randint(200, 400)
        chrom_positions = np.sort(np.random.randint(1e6, 2e8, n_variants))
        chrom_pvalues = np.random.uniform(0, 1, n_variants)
        # Add hits at different locations per cohort
        if chrom == 6:
            n_hits = np.random.randint(3, 8)
            hit_indices = np.random.choice(n_variants, n_hits, replace=False)
            chrom_pvalues[hit_indices] = 10 ** np.random.uniform(-15, -8, n_hits)
        if chrom == 11 and i >= 1:  # Replicated signal
            n_hits = np.random.randint(2, 5)
            hit_indices = np.random.choice(n_variants, n_hits, replace=False)
            chrom_pvalues[hit_indices] = 10 ** np.random.uniform(-12, -8, n_hits)
        if chrom == 17 and i == 2:  # Meta-analysis only signal
            n_hits = np.random.randint(2, 4)
            hit_indices = np.random.choice(n_variants, n_hits, replace=False)
            chrom_pvalues[hit_indices] = 10 ** np.random.uniform(-10, -8, n_hits)
        for j in range(n_variants):
            gwas_data.append(
                {"chrom": str(chrom), "pos": chrom_positions[j], "p": chrom_pvalues[j]}
            )
    stacked_gwas_dfs.append(pd.DataFrame(gwas_data))

fig = manhattan_plotter.plot_manhattan_stacked(
    stacked_gwas_dfs,
    panel_labels=panel_names,
    significance_threshold=5e-8,
    figsize=(14, 9),
    title="Multi-cohort GWAS Comparison",
)
fig.savefig("examples/manhattan_stacked.png", dpi=150, bbox_inches="tight")
print("   Saved: examples/manhattan_stacked.png")

# Interactive Plotly stacked Manhattan plot
print("25. Interactive Plotly stacked Manhattan plot...")
fig = manhattan_plotter_plotly.plot_manhattan_stacked(
    stacked_gwas_dfs,
    panel_labels=panel_names,
    significance_threshold=5e-8,
    figsize=(14, 9),
    title="Multi-cohort GWAS Comparison",
)
fig.write_html("examples/manhattan_stacked_plotly.html")
print("   Saved: examples/manhattan_stacked_plotly.html")

# Interactive Bokeh stacked Manhattan plot
print("26. Interactive Bokeh stacked Manhattan plot...")
fig = manhattan_plotter_bokeh.plot_manhattan_stacked(
    stacked_gwas_dfs,
    panel_labels=panel_names,
    significance_threshold=5e-8,
    figsize=(14, 9),
    title="Multi-cohort GWAS Comparison",
)
save(
    fig,
    filename="examples/manhattan_stacked_bokeh.html",
    resources=CDN,
    title="Stacked Manhattan Plot",
)
print("   Saved: examples/manhattan_stacked_bokeh.html")

# Side-by-side Manhattan + QQ plot
print("27. Side-by-side Manhattan + QQ plot...")
fig = manhattan_plotter.plot_manhattan_qq(
    manhattan_df,
    significance_threshold=5e-8,
    show_confidence_band=True,
    show_lambda=True,
    figsize=(16, 5),
    title="GWAS Summary",
)
fig.savefig("examples/manhattan_qq_sidebyside.png", dpi=150, bbox_inches="tight")
print("   Saved: examples/manhattan_qq_sidebyside.png")

# Interactive Plotly side-by-side Manhattan + QQ plot
print("28. Interactive Plotly side-by-side Manhattan + QQ plot...")
fig = manhattan_plotter_plotly.plot_manhattan_qq(
    manhattan_df,
    significance_threshold=5e-8,
    show_confidence_band=True,
    show_lambda=True,
    figsize=(16, 5),
    title="GWAS Summary",
)
fig.write_html("examples/manhattan_qq_plotly.html")
print("   Saved: examples/manhattan_qq_plotly.html")

# Interactive Bokeh side-by-side Manhattan + QQ plot
print("29. Interactive Bokeh side-by-side Manhattan + QQ plot...")
fig = manhattan_plotter_bokeh.plot_manhattan_qq(
    manhattan_df,
    significance_threshold=5e-8,
    show_confidence_band=True,
    show_lambda=True,
    figsize=(16, 5),
    title="GWAS Summary",
)
save(
    fig,
    filename="examples/manhattan_qq_bokeh.html",
    resources=CDN,
    title="Manhattan + QQ Plot",
)
print("   Saved: examples/manhattan_qq_bokeh.html")

print("\nAll plots generated successfully!")
print("\nInteractive HTML files can be opened in a browser to test hover tooltips.")
