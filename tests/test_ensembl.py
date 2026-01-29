# tests/test_ensembl.py
"""Tests for Ensembl REST API integration."""

import tempfile
from pathlib import Path
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

    with patch(
        "pylocuszoom.ensembl.requests.get", side_effect=[mock_429, mock_success]
    ):
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
    from pylocuszoom.ensembl import get_cached_genes, save_cached_genes

    df = pd.DataFrame(
        {
            "chr": ["13", "13"],
            "start": [32315474, 32400000],
            "end": [32400266, 32500000],
            "gene_name": ["BRCA2", "TEST"],
            "strand": ["+", "-"],
        }
    )

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


# --- get_genes_for_region tests ---


def test_get_genes_for_region_uses_cache():
    """Test that get_genes_for_region uses cache when available."""
    from pylocuszoom.ensembl import get_genes_for_region, save_cached_genes

    # Pre-populate cache
    cached_df = pd.DataFrame(
        {
            "chr": ["13"],
            "start": [32315474],
            "end": [32400266],
            "gene_name": ["CACHED_GENE"],
            "strand": ["+"],
        }
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        cache_dir = Path(tmpdir)
        save_cached_genes(cached_df, cache_dir, "human", "13", 32000000, 33000000)

        # Mock the API call - should NOT be called due to cache hit
        with patch("pylocuszoom.ensembl.requests.get") as mock_get:
            result = get_genes_for_region(
                species="human",
                chrom="13",
                start=32000000,
                end=33000000,
                cache_dir=cache_dir,
            )

            # API should not have been called
            mock_get.assert_not_called()

        assert len(result) == 1
        assert result["gene_name"].iloc[0] == "CACHED_GENE"


def test_get_genes_for_region_fetches_and_caches():
    """Test that get_genes_for_region fetches from API and caches result."""
    from pylocuszoom.ensembl import get_genes_for_region

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
    ]

    with tempfile.TemporaryDirectory() as tmpdir:
        cache_dir = Path(tmpdir)

        with patch("pylocuszoom.ensembl.requests.get", return_value=mock_response):
            result = get_genes_for_region(
                species="human",
                chrom="13",
                start=32000000,
                end=33000000,
                cache_dir=cache_dir,
            )

        assert len(result) == 1
        assert result["gene_name"].iloc[0] == "BRCA2"

        # Verify cache file was created (CSV, not parquet)
        csv_files = list(cache_dir.glob("**/*.csv"))
        assert len(csv_files) == 1


def test_get_genes_for_region_include_exons():
    """Test fetching genes with exons included."""
    from pylocuszoom.ensembl import get_genes_for_region

    mock_genes = Mock()
    mock_genes.ok = True
    mock_genes.json.return_value = [
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

    mock_exons = Mock()
    mock_exons.ok = True
    mock_exons.json.return_value = [
        {
            "id": "ENSE00003659301",
            "Parent": "ENST00000380152",
            "seq_region_name": "13",
            "start": 32315474,
            "end": 32315667,
            "strand": 1,
            "feature_type": "exon",
        },
    ]

    with tempfile.TemporaryDirectory() as tmpdir:
        cache_dir = Path(tmpdir)

        with patch(
            "pylocuszoom.ensembl.requests.get", side_effect=[mock_genes, mock_exons]
        ):
            genes_df, exons_df = get_genes_for_region(
                species="human",
                chrom="13",
                start=32000000,
                end=33000000,
                cache_dir=cache_dir,
                include_exons=True,
            )

        assert len(genes_df) == 1
        assert len(exons_df) == 1


# --- clear_ensembl_cache tests ---


def test_clear_ensembl_cache():
    """Test clearing the Ensembl cache."""
    from pylocuszoom.ensembl import clear_ensembl_cache, save_cached_genes

    with tempfile.TemporaryDirectory() as tmpdir:
        cache_dir = Path(tmpdir)

        # Create some cache files
        df = pd.DataFrame(
            {
                "chr": ["1"],
                "start": [100],
                "end": [200],
                "gene_name": ["X"],
                "strand": ["+"],
            }
        )
        save_cached_genes(df, cache_dir, "human", "1", 100, 200)
        save_cached_genes(df, cache_dir, "mouse", "1", 100, 200)

        # Should have 2 CSV files in species subdirs
        csv_files = list(cache_dir.glob("**/*.csv"))
        assert len(csv_files) == 2

        # Clear cache
        deleted = clear_ensembl_cache(cache_dir)

        assert deleted == 2
        assert len(list(cache_dir.glob("**/*.csv"))) == 0


def test_clear_ensembl_cache_species_specific():
    """Test clearing cache for specific species only."""
    from pylocuszoom.ensembl import clear_ensembl_cache, save_cached_genes

    with tempfile.TemporaryDirectory() as tmpdir:
        cache_dir = Path(tmpdir)

        df = pd.DataFrame(
            {
                "chr": ["1"],
                "start": [100],
                "end": [200],
                "gene_name": ["X"],
                "strand": ["+"],
            }
        )
        save_cached_genes(df, cache_dir, "human", "1", 100, 200)
        save_cached_genes(df, cache_dir, "mouse", "1", 100, 200)

        # Clear only human cache
        deleted = clear_ensembl_cache(cache_dir, species="human")

        assert deleted == 1
        # Mouse cache should still exist
        assert len(list(cache_dir.glob("**/*.csv"))) == 1


# --- Consolidation tests ---


class TestNormalizeChromConsolidation:
    """Verify ensembl uses shared normalize_chrom from utils."""

    def test_ensembl_uses_utils_normalize_chrom(self):
        """Confirm _normalize_chrom was removed from ensembl module."""
        import pylocuszoom.ensembl as ensembl_module

        # After consolidation, _normalize_chrom should not exist in ensembl
        assert not hasattr(ensembl_module, "_normalize_chrom"), (
            "_normalize_chrom should be removed from ensembl.py - "
            "use normalize_chrom from utils instead"
        )


# --- Export tests ---


def test_ensembl_functions_exported():
    """Test that Ensembl functions are exported from main package."""
    from pylocuszoom import (
        clear_ensembl_cache,
        fetch_genes_from_ensembl,
        get_genes_for_region,
    )

    assert callable(get_genes_for_region)
    assert callable(fetch_genes_from_ensembl)
    assert callable(clear_ensembl_cache)
