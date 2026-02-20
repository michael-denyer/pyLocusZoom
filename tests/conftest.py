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
    np.random.seed(42)
    n_snps = 100
    positions = np.sort(np.random.randint(1000000, 2000000, n_snps))

    return pd.DataFrame(
        {
            "rs": [f"rs{i}" for i in range(n_snps)],
            "chr": [1] * n_snps,
            "ps": positions,
            "p_wald": np.random.uniform(1e-10, 1, n_snps),
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
