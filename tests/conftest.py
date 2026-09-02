"""Pytest configuration for pyLocusZoom tests."""

import os

from hypothesis import Verbosity, settings

# Hypothesis profiles for different environments
# CI: more examples for thorough testing
settings.register_profile("ci", max_examples=100, deadline=None)
# Dev: fewer examples for faster feedback
settings.register_profile("dev", max_examples=20, deadline=None)
# Debug: verbose output for troubleshooting
settings.register_profile(
    "debug", max_examples=10, verbosity=Verbosity.verbose, deadline=None
)

settings.load_profile(os.getenv("HYPOTHESIS_PROFILE", "dev"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402


def _figure_types() -> dict[str, type]:
    """Map each built-in backend to the figure type its plotters return."""
    import bokeh.models
    import matplotlib.figure
    import plotly.graph_objects

    return {
        "matplotlib": matplotlib.figure.Figure,
        "plotly": plotly.graph_objects.Figure,
        "bokeh": bokeh.models.LayoutDOM,
    }


FIGURE_TYPES = _figure_types()


@pytest.fixture(autouse=True)
def close_matplotlib_figures():
    """Close every pyplot figure a test leaves open.

    Plotters build figures through ``plt.subplots``, which pyplot retains until
    something closes them. Without this the suite warns past 20 open figures and
    each xdist worker holds every figure it ever drew.
    """
    yield

    import matplotlib.pyplot as plt

    plt.close("all")


@pytest.fixture
def plink_assoc_file(tmp_path):
    """Three-SNP PLINK .assoc file on disk."""
    content = """CHR SNP BP A1 TEST NMISS BETA STAT P
1 rs123 1000000 A ADD 1000 0.5 2.5 0.01
1 rs456 1001000 G ADD 1000 0.3 1.5 0.1
1 rs789 1002000 T ADD 1000 -0.2 -1.0 1e-8
"""
    filepath = tmp_path / "test.assoc"
    filepath.write_text(content)
    return filepath


@pytest.fixture
def regenie_file(tmp_path):
    """Three-SNP REGENIE file on disk, p-values as LOG10P."""
    content = """CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA
1 1000000 rs123 A G 0.3 1000 ADD 0.5 0.2 6.25 2.0 NA
1 1001000 rs456 C T 0.2 1000 ADD 0.3 0.15 4.0 1.5 NA
1 1002000 rs789 G A 0.4 1000 ADD -0.2 0.1 4.0 8.0 NA
"""
    filepath = tmp_path / "test.regenie"
    filepath.write_text(content)
    return filepath


@pytest.fixture
def susie_file(tmp_path):
    """Four-variant SuSiE results file on disk, two credible sets."""
    content = """pos\tpip\tcs\tsnp
1000000\t0.85\t1\trs123
1001000\t0.12\t1\trs456
1002000\t0.02\t0\trs789
1003000\t0.45\t2\trs101
"""
    filepath = tmp_path / "susie.tsv"
    filepath.write_text(content)
    return filepath


@pytest.fixture
def bed_file(tmp_path):
    """Three-gene BED file on disk."""
    content = """chr1\t1000000\t1020000\tGENE1
chr1\t1050000\t1080000\tGENE2
chr1\t1100000\t1150000\tGENE3
"""
    filepath = tmp_path / "genes.bed"
    filepath.write_text(content)
    return filepath


@pytest.fixture
def sample_genes_df():
    """Sample gene annotations DataFrame for testing."""
    return pd.DataFrame(
        {
            "gene_name": ["GENE_A", "GENE_B", "GENE_C"],
            "chr": [1, 1, 1],
            "start": [1100000, 1400000, 1700000],
            "end": [1150000, 1500000, 1800000],
            "strand": ["+", "-", "+"],
        }
    )


@pytest.fixture
def sample_recomb_df():
    """Sample recombination rate DataFrame for testing."""
    return pd.DataFrame(
        {
            "pos": [1000000, 1200000, 1400000, 1600000, 1800000, 2000000],
            "rate": [0.5, 1.2, 2.5, 1.8, 0.8, 0.3],
        }
    )


@pytest.fixture
def regional_gwas_df():
    """100-SNP single-chromosome region in the rs/chr/ps/p_wald schema."""
    rng = np.random.default_rng(42)
    n_snps = 100
    positions = np.sort(rng.integers(1000000, 2000000, n_snps))

    return pd.DataFrame(
        {
            "rs": [f"rs{i}" for i in range(n_snps)],
            "chr": [1] * n_snps,
            "ps": positions,
            "p_wald": rng.uniform(1e-10, 1, n_snps),
        }
    )


@pytest.fixture
def tiny_regional_gwas_df():
    """Three-SNP region in the rs/ps/p_wald schema, no chr column."""
    return pd.DataFrame(
        {
            "rs": ["rs1", "rs2", "rs3"],
            "ps": [1100000, 1500000, 1900000],
            "p_wald": [1e-8, 1e-5, 1e-3],
        }
    )


@pytest.fixture
def sample_eqtl_df():
    """eQTL associations with signed effect sizes for one gene."""
    return pd.DataFrame(
        {
            "pos": [1200000, 1400000, 1600000, 1800000],
            "p_value": [1e-8, 1e-6, 1e-4, 1e-5],
            "effect_size": [0.5, -0.3, 0.8, -0.2],
            "gene": ["GENE_A", "GENE_A", "GENE_A", "GENE_A"],
        }
    )


@pytest.fixture
def sample_finemapping_df():
    """Fine-mapping PIPs spanning two credible sets plus non-CS variants."""
    return pd.DataFrame(
        {
            "pos": [1200000, 1300000, 1400000, 1500000, 1600000],
            "pip": [0.85, 0.10, 0.03, 0.45, 0.30],
            "cs": [1, 1, 0, 2, 2],
        }
    )


@pytest.fixture
def small_regional_gwas_df():
    """Five-SNP region in the rs/ps/p_wald schema, no chr column."""
    return pd.DataFrame(
        {
            "rs": ["rs1", "rs2", "rs3", "rs4", "rs5"],
            "ps": [1100000, 1300000, 1500000, 1700000, 1900000],
            "p_wald": [1e-8, 1e-6, 1e-5, 1e-4, 0.01],
        }
    )


@pytest.fixture
def manhattan_gwas_df():
    """Multi-chromosome GWAS in the chrom/pos/p schema."""
    rng = np.random.default_rng(42)
    return pd.DataFrame(
        {
            "chrom": np.repeat([1, 2, 3], [40, 30, 30]),
            "pos": np.concatenate(
                [
                    np.sort(rng.integers(int(1e6), int(1e8), 40)),
                    np.sort(rng.integers(int(1e6), int(1e8), 30)),
                    np.sort(rng.integers(int(1e6), int(1e8), 30)),
                ]
            ),
            "p": rng.uniform(1e-10, 1, 100),
        }
    )


@pytest.fixture
def regional_plotter():
    """LocusZoomPlotter with logging left at the caller's level."""
    from pylocuszoom import LocusZoomPlotter

    return LocusZoomPlotter(species=None, log_level=None)


@pytest.fixture
def canine_plotter():
    """LocusZoomPlotter on the canine default build."""
    from pylocuszoom import LocusZoomPlotter

    return LocusZoomPlotter(species="canine")


@pytest.fixture
def manhattan_plotter():
    """ManhattanPlotter on the human genome build."""
    from pylocuszoom import ManhattanPlotter

    return ManhattanPlotter(species="human")


@pytest.fixture
def labelled_gwas_df():
    """Five SNPs with neglog10p precomputed, the shape labels.py consumes."""
    return pd.DataFrame(
        {
            "rs": ["rs1", "rs2", "rs3", "rs4", "rs5"],
            "ps": [1100000, 1200000, 1300000, 1400000, 1500000],
            "neglog10p": [8, 5, 3, 6, 9],
        }
    )


@pytest.fixture
def manhattan_rs_gwas_df():
    """Fifty single-chromosome variants in the rs/chrom/pos/p schema."""
    rng = np.random.default_rng(42)
    n_snps = 50
    positions = np.sort(rng.integers(1_000_000, 2_000_000, n_snps))
    return pd.DataFrame(
        {
            "rs": [f"rs{i}" for i in range(n_snps)],
            "chrom": [1] * n_snps,
            "pos": positions,
            "p": rng.uniform(1e-10, 1, n_snps),
        }
    )


@pytest.fixture
def phewas_with_effects_df():
    """PheWAS associations carrying signed effect sizes."""
    return pd.DataFrame(
        {
            "phenotype": ["Type 2 Diabetes", "BMI", "Height", "Blood Pressure"],
            "category": ["Metabolic", "Metabolic", "Anthropometric", "Cardiovascular"],
            "p": [1e-8, 1e-4, 0.05, 1e-6],
            "effect": [0.3, -0.2, 0.1, 0.25],
        }
    )


@pytest.fixture
def sample_forest_df():
    """Per-study effects with confidence intervals, plus a meta-analysis row."""
    return pd.DataFrame(
        {
            "study": ["Study A", "Study B", "Study C", "Meta-analysis"],
            "effect": [0.25, 0.30, 0.20, 0.24],
            "ci_lower": [0.10, 0.15, 0.05, 0.18],
            "ci_upper": [0.40, 0.45, 0.35, 0.30],
        }
    )


@pytest.fixture
def stats_plotter():
    """StatsPlotter with default settings."""
    from pylocuszoom import StatsPlotter

    return StatsPlotter()


@pytest.fixture
def fake_plink(tmp_path):
    """Drive LD tests through the real PLINK boundary: subprocess plus its files.

    Returns the bfile prefix and a context-manager factory. Patching
    ``subprocess.run`` rather than ``calculate_ld`` keeps command construction,
    output parsing and the R2 merge inside the test, so a change to any of them
    fails here instead of passing against a mock.
    """
    import pathlib
    import subprocess
    from unittest.mock import patch

    bfile = tmp_path / "genotypes"
    for suffix in (".bed", ".bim", ".fam"):
        bfile.with_suffix(suffix).write_bytes(b"")

    def plink_writes(ld_body, returncode=0, stderr=""):
        def fake_run(cmd, **kwargs):
            out_prefix = cmd[cmd.index("--out") + 1]
            if ld_body is not None:
                pathlib.Path(f"{out_prefix}.ld").write_text(ld_body)
            return subprocess.CompletedProcess(cmd, returncode, "", stderr)

        return patch("subprocess.run", side_effect=fake_run)

    return str(bfile), plink_writes


@pytest.fixture
def warning_records():
    """Collect pylocuszoom WARNING messages; loguru does not feed caplog.

    The wrapper in ``pylocuszoom.logging`` drops warnings entirely while
    disabled, so the sink alone would capture nothing after any test that
    called ``disable_logging()``.
    """
    from loguru import logger as loguru_logger

    from pylocuszoom.logging import logger

    records: list[str] = []
    handler_id = loguru_logger.add(
        records.append,
        level="WARNING",
        format="{message}",
        filter=lambda record: record["name"].startswith("pylocuszoom"),
    )
    was_enabled = logger._enabled
    logger._enabled = True
    try:
        yield records
    finally:
        logger._enabled = was_enabled
        loguru_logger.remove(handler_id)
