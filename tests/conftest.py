"""Pytest configuration for snp-scope-plot tests."""

import numpy as np
import pandas as pd
import pytest


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
