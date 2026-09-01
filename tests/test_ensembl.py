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
    """Test cache directory follows pylocuszoom convention."""
    from pylocuszoom.ensembl import get_ensembl_cache_dir

    cache_dir = get_ensembl_cache_dir()
    assert isinstance(cache_dir, Path)
    assert "pylocuszoom" in str(cache_dir)
    assert "ensembl" in str(cache_dir)


def test_ensembl_cache_dir_honors_xdg_on_macos(tmp_path, monkeypatch):
    """macOS honors XDG_CACHE_HOME (previously hardcoded ~/.cache)."""
    from pylocuszoom.ensembl import get_ensembl_cache_dir

    monkeypatch.setattr("sys.platform", "darwin")
    monkeypatch.setattr("os.path.exists", lambda _: False)
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path))

    assert get_ensembl_cache_dir() == tmp_path / "pylocuszoom" / "ensembl"


def test_ensembl_cache_dir_routes_to_dbfs_on_databricks(monkeypatch):
    """Databricks routes the Ensembl cache under /dbfs (previously ~/.cache)."""
    from pylocuszoom.ensembl import get_ensembl_cache_dir

    monkeypatch.setattr("sys.platform", "linux")
    monkeypatch.setattr("os.path.exists", lambda p: p == "/dbfs")
    monkeypatch.setattr(Path, "mkdir", lambda self, *a, **k: None)

    assert get_ensembl_cache_dir() == Path("/dbfs/FileStore/reference_data/ensembl")


def test_cache_resolvers_share_base(tmp_path, monkeypatch):
    """Recombination and Ensembl caches resolve under one shared root."""
    from pylocuszoom.ensembl import get_ensembl_cache_dir
    from pylocuszoom.recombination import get_default_data_dir

    monkeypatch.setattr("sys.platform", "darwin")
    monkeypatch.setattr("os.path.exists", lambda _: False)
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path))
    monkeypatch.setattr(Path, "mkdir", lambda self, *a, **k: None)

    assert get_default_data_dir() == tmp_path / "pylocuszoom" / "recombination_maps"
    assert get_ensembl_cache_dir() == tmp_path / "pylocuszoom" / "ensembl"
    assert get_ensembl_cache_dir().parent == get_default_data_dir().parent


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


class TestPathTraversalProtection:
    """Tests that malicious species names cannot escape the cache directory."""

    def test_path_traversal_species_raises_validation_error(self, tmp_path):
        """Species name with '../' should raise ValidationError."""
        from pylocuszoom.ensembl import _safe_species_dir
        from pylocuszoom.utils import ValidationError

        with pytest.raises(ValidationError, match="Invalid species name"):
            _safe_species_dir(tmp_path, "../../etc")

    def test_absolute_path_species_raises_validation_error(self, tmp_path):
        """Absolute path as species should raise ValidationError."""
        from pylocuszoom.ensembl import _safe_species_dir
        from pylocuszoom.utils import ValidationError

        with pytest.raises(ValidationError, match="Invalid species name"):
            _safe_species_dir(tmp_path, "/etc/passwd")

    def test_safe_species_stays_within_cache(self, tmp_path):
        """Normal species name should resolve within cache_dir."""
        from pylocuszoom.ensembl import _safe_species_dir

        result = _safe_species_dir(tmp_path, "canis_lupus_familiaris")
        assert result.is_relative_to(tmp_path)

    def test_get_cached_genes_rejects_traversal(self, tmp_path):
        """get_cached_genes should reject path traversal species."""
        from pylocuszoom.ensembl import get_cached_genes
        from pylocuszoom.utils import ValidationError

        with pytest.raises(ValidationError, match="Invalid species name"):
            get_cached_genes(tmp_path, "../../etc/passwd", "1", 1, 100)

    def test_save_cached_genes_rejects_traversal(self, tmp_path):
        """save_cached_genes should reject path traversal species."""
        import pandas as pd

        from pylocuszoom.ensembl import save_cached_genes
        from pylocuszoom.utils import ValidationError

        df = pd.DataFrame({"gene": ["BRCA1"]})
        with pytest.raises(ValidationError, match="Invalid species name"):
            save_cached_genes(df, tmp_path, "../../etc/passwd", "1", 1, 100)

    def test_clear_ensembl_cache_rejects_traversal(self, tmp_path):
        """clear_ensembl_cache should reject path traversal species."""
        from pylocuszoom.ensembl import clear_ensembl_cache
        from pylocuszoom.utils import ValidationError

        with pytest.raises(ValidationError, match="Invalid species name"):
            clear_ensembl_cache(tmp_path, species="../../etc/passwd")


class TestEmptyResultCaching:
    """A cached empty region must reload, and a failed fetch must not be cached.

    Both bugs surfaced as ``pandas.errors.EmptyDataError`` on the second call:
    an empty DataFrame serialises to a one-byte file that cannot be parsed back.
    """

    @staticmethod
    def _ok_response(payload):
        response = Mock()
        response.ok = True
        response.status_code = 200
        response.json.return_value = payload
        return response

    @staticmethod
    def _gene_payload():
        return [
            {
                "feature_type": "gene",
                "seq_region_name": "1",
                "start": 1000,
                "end": 2000,
                "external_name": "BRCA2",
                "strand": 1,
                "id": "ENSG00000139618",
                "biotype": "protein_coding",
            }
        ]

    def test_gene_sparse_region_reloads_from_cache(self, tmp_path):
        """A region with no genes caches and reloads without raising."""
        from pylocuszoom.ensembl import get_genes_for_region

        with patch(
            "pylocuszoom.ensembl.requests.get", return_value=self._ok_response([])
        ) as mock_get:
            first = get_genes_for_region("human", "1", 1_000_000, 1_100_000, tmp_path)
            assert mock_get.call_count == 1

            second = get_genes_for_region("human", "1", 1_000_000, 1_100_000, tmp_path)

        assert first.empty
        assert second.empty
        assert mock_get.call_count == 1, "second call must be served from cache"
        assert list(second.columns) == list(first.columns)
        assert "gene_name" in second.columns

    def test_failed_fetch_is_not_cached(self, tmp_path):
        """An API outage must not poison the cache for the region."""
        import requests

        from pylocuszoom.ensembl import get_genes_for_region

        with (
            patch("pylocuszoom.ensembl.time.sleep"),
            patch(
                "pylocuszoom.ensembl.requests.get",
                side_effect=requests.exceptions.ConnectionError("network down"),
            ),
        ):
            during_outage = get_genes_for_region(
                "human", "1", 1_000_000, 1_100_000, tmp_path
            )
        assert during_outage.empty

        assert not list(tmp_path.rglob("genes_*.csv")), (
            "a failed fetch must leave no cache file"
        )

        with patch(
            "pylocuszoom.ensembl.requests.get",
            return_value=self._ok_response(self._gene_payload()),
        ) as mock_get:
            after_recovery = get_genes_for_region(
                "human", "1", 1_000_000, 1_100_000, tmp_path
            )

        assert mock_get.called, "recovery must retry the API, not reuse a failed result"
        assert len(after_recovery) == 1
        assert after_recovery["gene_name"].iloc[0] == "BRCA2"

    def test_zero_byte_cache_file_is_ignored(self, tmp_path):
        """A cache file poisoned by an older release is treated as a miss."""
        from pylocuszoom.ensembl import _cache_key, get_cached_genes

        species_dir = tmp_path / "homo_sapiens"
        species_dir.mkdir()
        key = _cache_key("homo_sapiens", "1", 1_000_000, 1_100_000)
        (species_dir / f"genes_{key}.csv").write_text("\n")

        assert get_cached_genes(tmp_path, "human", "1", 1_000_000, 1_100_000) is None


class TestAssemblyMismatch:
    """Ensembl serves one assembly per species and ignores coord_system_version.

    Asking for canine genes while working in CanFam3.1 returns ROS_Cfam_1.0
    coordinates with no error, shifting ATP9B on chr1 from 1,136,865 to
    938,796. The mismatch has to be loud, and the cache must not mix builds.
    """

    @staticmethod
    def _ok_response(payload):
        response = Mock()
        response.ok = True
        response.status_code = 200
        response.json.return_value = payload
        return response

    @staticmethod
    def _dog_payload():
        return [
            {
                "feature_type": "gene",
                "seq_region_name": "1",
                "start": 938796,
                "end": 1175952,
                "external_name": "ATP9B",
                "strand": -1,
                "id": "ENSCAFG00845000134",
                "biotype": "protein_coding",
                "assembly_name": "ROS_Cfam_1.0",
            }
        ]

    def test_assembly_token_folds_synonyms(self):
        """Equivalent spellings of one assembly compare equal."""
        from pylocuszoom.ensembl import _assembly_token

        assert _assembly_token("CanFam4.0") == _assembly_token("UU_Cfam_GSD_1.0")
        assert _assembly_token("CanFam3.1") == _assembly_token("canfam3")
        assert _assembly_token("hg38") == _assembly_token("GRCh38")
        assert _assembly_token("CanFam3.1") != _assembly_token("ROS_Cfam_1.0")

    def test_cache_key_separates_builds(self):
        """The same region under two builds must not share a cache entry."""
        from pylocuszoom.ensembl import _cache_key

        canfam3 = _cache_key("canis_lupus_familiaris", "1", 1, 100, "canfam3.1")
        canfam4 = _cache_key("canis_lupus_familiaris", "1", 1, 100, "canfam4")
        assert canfam3 != canfam4

    def test_gene_record_carries_assembly(self):
        """Every gene row records the assembly Ensembl served it on."""
        from pylocuszoom.ensembl import fetch_genes_from_ensembl

        with patch(
            "pylocuszoom.ensembl.requests.get",
            return_value=self._ok_response(self._dog_payload()),
        ):
            genes = fetch_genes_from_ensembl("canine", "1", 900_000, 1_200_000)

        assert genes["assembly"].tolist() == ["ROS_Cfam_1.0"]

    def test_exon_record_carries_assembly(self):
        """Exon rows record their assembly too, not just genes."""
        from pylocuszoom.ensembl import fetch_exons_from_ensembl

        exon = dict(self._dog_payload()[0], feature_type="exon", id="ENSCAFE00000001")
        with patch(
            "pylocuszoom.ensembl.requests.get", return_value=self._ok_response([exon])
        ):
            with pytest.warns(UserWarning, match="ROS_Cfam_1.0"):
                exons = fetch_exons_from_ensembl(
                    "canine", "1", 900_000, 1_200_000, genome_build="canfam3.1"
                )

        assert exons["assembly"].tolist() == ["ROS_Cfam_1.0"]

    def test_fetch_warns_when_assembly_differs(self):
        """A CanFam3.1 caller is told the genes came back on ROS_Cfam_1.0."""
        from pylocuszoom.ensembl import fetch_genes_from_ensembl

        with patch(
            "pylocuszoom.ensembl.requests.get",
            return_value=self._ok_response(self._dog_payload()),
        ):
            with pytest.warns(UserWarning, match="ROS_Cfam_1.0"):
                fetch_genes_from_ensembl(
                    "canine", "1", 900_000, 1_200_000, genome_build="canfam3.1"
                )

    def test_fetch_silent_when_assembly_matches(self, recwarn):
        """No warning when the caller already works in Ensembl's assembly."""
        from pylocuszoom.ensembl import fetch_genes_from_ensembl

        with patch(
            "pylocuszoom.ensembl.requests.get",
            return_value=self._ok_response(self._dog_payload()),
        ):
            fetch_genes_from_ensembl(
                "canine", "1", 900_000, 1_200_000, genome_build="ROS_Cfam_1.0"
            )

        assert [w for w in recwarn.list if w.category is UserWarning] == []

    def test_cache_hit_still_warns(self, tmp_path):
        """Reloading from cache in a fresh session repeats the warning."""
        from pylocuszoom.ensembl import get_genes_for_region

        with patch(
            "pylocuszoom.ensembl.requests.get",
            return_value=self._ok_response(self._dog_payload()),
        ) as mock_get:
            with pytest.warns(UserWarning, match="ROS_Cfam_1.0"):
                get_genes_for_region(
                    "canine",
                    "1",
                    900_000,
                    1_200_000,
                    tmp_path,
                    genome_build="canfam3.1",
                )
            with pytest.warns(UserWarning, match="ROS_Cfam_1.0"):
                cached = get_genes_for_region(
                    "canine",
                    "1",
                    900_000,
                    1_200_000,
                    tmp_path,
                    genome_build="canfam3.1",
                )

        assert mock_get.call_count == 1, "second call must be served from cache"
        assert cached["assembly"].tolist() == ["ROS_Cfam_1.0"]

    def test_plotter_passes_its_build_through(self):
        """auto_genes must hand the plotter's build to the Ensembl fetch."""
        from pylocuszoom import LocusZoomPlotter

        plotter = LocusZoomPlotter(species="canine", auto_genes=True)
        gwas = pd.DataFrame(
            {"ps": [1_000_000, 1_050_000], "p_wald": [0.5, 1e-6], "rs": ["a", "b"]}
        )

        with patch(
            "pylocuszoom.plotter.get_genes_for_region",
            return_value=pd.DataFrame(
                columns=["chr", "start", "end", "gene_name", "assembly"]
            ),
        ) as mock_fetch:
            plotter.plot(gwas, chrom=1, start=900_000, end=1_200_000)

        assert mock_fetch.call_args.kwargs["genome_build"] == plotter.genome_build
