# tests/test_ensembl.py
"""Tests for Ensembl REST API integration."""

from unittest.mock import Mock, patch

import pandas as pd
import pytest


def test_get_ensembl_species_name_canine():
    """Test species name mapping for canine."""
    from pylocuszoom.ensembl import get_ensembl_species_name

    assert get_ensembl_species_name("canine") == "canis_lupus_familiaris"
    assert get_ensembl_species_name("dog") == "canis_lupus_familiaris"


def test_get_ensembl_species_name_human():
    """Test species name mapping for human."""
    from pylocuszoom.ensembl import get_ensembl_species_name

    assert get_ensembl_species_name("human") == "homo_sapiens"
    assert get_ensembl_species_name("homo_sapiens") == "homo_sapiens"


def test_get_ensembl_species_name_unknown():
    """Test unknown species returns input unchanged."""
    from pylocuszoom.ensembl import get_ensembl_species_name

    assert get_ensembl_species_name("my_custom_species") == "my_custom_species"


def test_fetch_genes_from_ensembl_success():
    """Test fetching genes from Ensembl API with mocked response."""
    from pylocuszoom.ensembl import fetch_genes_from_ensembl

    mock_response = Mock()
    mock_response.ok = True
    mock_response.json.return_value = [
        {
            "id": "ENSG00000139618",
            "external_name": "BRCA2",
            "seq_region_name": "13",
            "start": 32315474,
            "end": 32400266,
            "strand": 1,
            "biotype": "protein_coding",
            "feature_type": "gene",
        },
        {
            "id": "ENSG00000012048",
            "external_name": "BRCA1",
            "seq_region_name": "17",
            "start": 43044295,
            "end": 43170245,
            "strand": -1,
            "biotype": "protein_coding",
            "feature_type": "gene",
        },
    ]

    with patch("pylocuszoom.ensembl.requests.get", return_value=mock_response):
        df = fetch_genes_from_ensembl("human", chrom="13", start=32000000, end=33000000)

    assert isinstance(df, pd.DataFrame)
    assert len(df) == 2
    assert "chr" in df.columns
    assert "start" in df.columns
    assert "end" in df.columns
    assert "gene_name" in df.columns
    assert "strand" in df.columns
    # Sort by start position for deterministic ordering
    df_sorted = df.sort_values("start")
    assert df_sorted["gene_name"].tolist() == ["BRCA2", "BRCA1"]


def test_fetch_genes_from_ensembl_api_error_warns():
    """Test handling of API errors - should warn and return empty DataFrame."""
    from pylocuszoom.ensembl import fetch_genes_from_ensembl

    mock_response = Mock()
    mock_response.ok = False
    mock_response.status_code = 503
    mock_response.text = "Service Unavailable"

    with patch("pylocuszoom.ensembl.requests.get", return_value=mock_response):
        df = fetch_genes_from_ensembl("human", chrom="13", start=32000000, end=33000000)

    # Should return empty DataFrame on error
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 0


def test_fetch_genes_region_too_large():
    """Test that regions > 5Mb raise ValidationError."""
    from pylocuszoom.ensembl import fetch_genes_from_ensembl
    from pylocuszoom.utils import ValidationError

    with pytest.raises(ValidationError, match="5Mb"):
        fetch_genes_from_ensembl("human", chrom="1", start=1000000, end=10000000)


def test_fetch_genes_retry_on_429():
    """Test that 429 rate limit responses trigger retry."""
    from pylocuszoom.ensembl import fetch_genes_from_ensembl

    # First call returns 429, second returns success
    mock_429 = Mock()
    mock_429.ok = False
    mock_429.status_code = 429
    mock_429.text = "Rate limited"

    mock_success = Mock()
    mock_success.ok = True
    mock_success.json.return_value = [
        {
            "id": "ENSG00000139618",
            "external_name": "BRCA2",
            "seq_region_name": "13",
            "start": 32315474,
            "end": 32400266,
            "strand": 1,
            "biotype": "protein_coding",
            "feature_type": "gene",
        },
    ]

    with patch("pylocuszoom.ensembl.requests.get", side_effect=[mock_429, mock_success]):
        with patch("pylocuszoom.ensembl.time.sleep"):  # Skip actual sleep
            df = fetch_genes_from_ensembl(
                "human", chrom="13", start=32000000, end=33000000
            )

    assert len(df) == 1
    assert df["gene_name"].iloc[0] == "BRCA2"


def test_fetch_exons_from_ensembl_success():
    """Test fetching exons from Ensembl API with mocked response."""
    from pylocuszoom.ensembl import fetch_exons_from_ensembl

    mock_response = Mock()
    mock_response.ok = True
    mock_response.json.return_value = [
        {
            "id": "ENSE00003659301",
            "Parent": "ENST00000380152",
            "seq_region_name": "13",
            "start": 32315474,
            "end": 32315667,
            "strand": 1,
            "feature_type": "exon",
        },
        {
            "id": "ENSE00003527960",
            "Parent": "ENST00000380152",
            "seq_region_name": "13",
            "start": 32316422,
            "end": 32316527,
            "strand": 1,
            "feature_type": "exon",
        },
    ]

    with patch("pylocuszoom.ensembl.requests.get", return_value=mock_response):
        df = fetch_exons_from_ensembl("human", chrom="13", start=32000000, end=33000000)

    assert isinstance(df, pd.DataFrame)
    assert len(df) == 2
    assert "chr" in df.columns
    assert "start" in df.columns
    assert "end" in df.columns
    assert "exon_id" in df.columns


def test_fetch_exons_region_too_large():
    """Test that regions > 5Mb raise ValidationError."""
    from pylocuszoom.ensembl import fetch_exons_from_ensembl
    from pylocuszoom.utils import ValidationError

    with pytest.raises(ValidationError, match="5Mb"):
        fetch_exons_from_ensembl("human", chrom="1", start=1000000, end=10000000)


# --- Caching tests ---
from pathlib import Path
import tempfile


def test_get_ensembl_cache_dir():
    """Test cache directory follows snp-scope-plot convention."""
    from pylocuszoom.ensembl import get_ensembl_cache_dir

    cache_dir = get_ensembl_cache_dir()
    assert isinstance(cache_dir, Path)
    assert "snp-scope-plot" in str(cache_dir)
    assert "ensembl" in str(cache_dir)


def test_get_cached_genes_miss():
    """Test cache miss returns None."""
    from pylocuszoom.ensembl import get_cached_genes

    with tempfile.TemporaryDirectory() as tmpdir:
        result = get_cached_genes(
            cache_dir=Path(tmpdir),
            species="human",
            chrom="13",
            start=32000000,
            end=33000000,
        )
        assert result is None


def test_save_and_load_cached_genes():
    """Test saving and loading cached genes using CSV."""
    from pylocuszoom.ensembl import save_cached_genes, get_cached_genes

    df = pd.DataFrame({
        "chr": ["13", "13"],
        "start": [32315474, 32400000],
        "end": [32400266, 32500000],
        "gene_name": ["BRCA2", "TEST"],
        "strand": ["+", "-"],
    })

    with tempfile.TemporaryDirectory() as tmpdir:
        cache_dir = Path(tmpdir)

        save_cached_genes(
            df,
            cache_dir=cache_dir,
            species="human",
            chrom="13",
            start=32000000,
            end=33000000,
        )

        # Verify CSV file created (not parquet)
        csv_files = list(cache_dir.glob("**/*.csv"))
        assert len(csv_files) == 1

        loaded = get_cached_genes(
            cache_dir=cache_dir,
            species="human",
            chrom="13",
            start=32000000,
            end=33000000,
        )

        assert loaded is not None
        assert len(loaded) == 2
        # Sort for deterministic comparison
        loaded_sorted = loaded.sort_values("start")
        assert loaded_sorted["gene_name"].tolist() == ["BRCA2", "TEST"]
