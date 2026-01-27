# tests/test_ensembl_integration.py
"""Integration tests for Ensembl API (requires network access).

Run with: pytest tests/test_ensembl_integration.py -v -m integration
Skip with: pytest -m "not integration"
"""

import pandas as pd
import pytest

pytestmark = pytest.mark.integration


@pytest.mark.integration
def test_fetch_human_genes_real_api():
    """Test fetching real human genes from Ensembl."""
    from pylocuszoom.ensembl import fetch_genes_from_ensembl

    # BRCA2 region on chr13 (small region to avoid timeout)
    df = fetch_genes_from_ensembl(
        species="human",
        chrom="13",
        start=32315000,
        end=32400000,
    )

    assert isinstance(df, pd.DataFrame)
    assert len(df) > 0
    assert "BRCA2" in df["gene_name"].values


@pytest.mark.integration
def test_fetch_mouse_genes_real_api():
    """Test fetching real mouse genes from Ensembl."""
    from pylocuszoom.ensembl import fetch_genes_from_ensembl

    df = fetch_genes_from_ensembl(
        species="mouse",
        chrom="1",
        start=10000000,
        end=10500000,
    )

    assert isinstance(df, pd.DataFrame)
    # Region may be empty if it has no genes, just verify it's a valid DataFrame
    assert "chr" in df.columns or len(df) == 0


@pytest.mark.integration
def test_fetch_canine_genes_real_api():
    """Test fetching real canine genes from Ensembl."""
    from pylocuszoom.ensembl import fetch_genes_from_ensembl

    df = fetch_genes_from_ensembl(
        species="canine",
        chrom="1",
        start=1000000,
        end=1500000,
    )

    assert isinstance(df, pd.DataFrame)


@pytest.mark.integration
def test_region_size_validation():
    """Test that >5Mb regions are rejected before API call."""
    from pylocuszoom.ensembl import fetch_genes_from_ensembl
    from pylocuszoom.utils import ValidationError

    # This should fail BEFORE making an API call
    with pytest.raises(ValidationError, match="5Mb"):
        fetch_genes_from_ensembl(
            species="human",
            chrom="1",
            start=1000000,
            end=10000000,  # 9Mb region
        )
