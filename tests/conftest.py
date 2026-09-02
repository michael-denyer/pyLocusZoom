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


@pytest.fixture
def sample_gwas_df():
    """Sample GWAS results DataFrame for testing."""
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
def sample_exons_df():
    """Sample exon annotations DataFrame for testing."""
    return pd.DataFrame(
        {
            "gene_name": ["GENE_A", "GENE_A", "GENE_B", "GENE_B", "GENE_C"],
            "chr": [1, 1, 1, 1, 1],
            "start": [1100000, 1120000, 1400000, 1450000, 1700000],
            "end": [1110000, 1130000, 1420000, 1470000, 1750000],
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
def stats_plotter():
    """StatsPlotter with default settings."""
    from pylocuszoom import StatsPlotter

    return StatsPlotter()


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
