# tests/test_ensembl.py
"""Tests for Ensembl REST API integration."""

import tempfile
from pathlib import Path
from unittest.mock import Mock, patch

import pandas as pd
import pytest

from pylocuszoom._gene_source import EXON_COLUMNS
from tests.reference_mocks import (
    gene_transcript_exon_payload,
    ok_response,
    ros_cfam_gene_payload,
)


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
    from pylocuszoom.ensembl import fetch_overlap_frames

    mock_response = ok_response(
        [
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
    )

    with patch("pylocuszoom._http.requests.get", return_value=mock_response):
        df, _ = fetch_overlap_frames("human", chrom="13", start=32000000, end=33000000)

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


def test_fetch_genes_from_ensembl_api_error_raises():
    """A service failure raises rather than passing as a gene-free region."""
    from pylocuszoom.ensembl import fetch_overlap_frames
    from pylocuszoom.exceptions import EnsemblAPIError

    mock_response = Mock()
    mock_response.ok = False
    mock_response.status_code = 503
    mock_response.text = "Service Unavailable"

    with (
        patch("pylocuszoom._http.time.sleep"),
        patch("pylocuszoom._http.requests.get", return_value=mock_response),
        pytest.raises(EnsemblAPIError),
    ):
        fetch_overlap_frames("human", chrom="13", start=32000000, end=33000000)


def test_fetch_genes_region_too_large():
    """Test that regions > 5Mb raise ValidationError."""
    from pylocuszoom.ensembl import fetch_overlap_frames
    from pylocuszoom.utils import ValidationError

    with pytest.raises(ValidationError, match="5Mb"):
        fetch_overlap_frames("human", chrom="1", start=1000000, end=10000000)


def test_fetch_genes_retry_on_429():
    """Test that 429 rate limit responses trigger retry."""
    from pylocuszoom.ensembl import fetch_overlap_frames

    # First call returns 429, second returns success
    mock_429 = Mock()
    mock_429.ok = False
    mock_429.status_code = 429
    mock_429.text = "Rate limited"

    mock_success = ok_response(
        [
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
    )

    with patch("pylocuszoom._http.requests.get", side_effect=[mock_429, mock_success]):
        with patch("pylocuszoom._http.time.sleep"):  # Skip actual sleep
            df, _ = fetch_overlap_frames(
                "human", chrom="13", start=32000000, end=33000000
            )

    assert len(df) == 1
    assert df["gene_name"].iloc[0] == "BRCA2"


def test_fetch_exons_from_ensembl_success():
    """Exons carry the symbol of the gene their transcript belongs to."""
    from pylocuszoom.ensembl import fetch_overlap_frames

    mock_response = ok_response(gene_transcript_exon_payload())

    with patch("pylocuszoom._http.requests.get", return_value=mock_response):
        _, df = fetch_overlap_frames("human", chrom="1", start=1000000, end=2000000)

    assert isinstance(df, pd.DataFrame)
    assert len(df) == 3
    assert df["gene_name"].tolist() == ["NFATC1"] * 3
    assert df["transcript_id"].tolist() == ["ENST01"] * 3
    assert df["exon_id"].tolist() == ["ENSE01", "ENSE02", "ENSE03"]


def test_exons_without_their_gene_are_dropped():
    """An exon whose gene the response omits cannot be drawn, so it is dropped."""
    from pylocuszoom.ensembl import fetch_overlap_frames

    orphaned = [
        feature
        for feature in gene_transcript_exon_payload()
        if feature["feature_type"] != "gene"
    ]
    mock_response = ok_response(orphaned)

    with patch("pylocuszoom._http.requests.get", return_value=mock_response):
        _, df = fetch_overlap_frames("human", chrom="13", start=32000000, end=33000000)

    assert df.empty
    assert list(df.columns) == list(EXON_COLUMNS)


# --- Caching tests ---


def test_ensembl_cache_dir():
    """Test cache directory follows pylocuszoom convention."""
    from pylocuszoom._gene_cache import cache_root

    cache_dir = cache_root("ensembl")
    assert isinstance(cache_dir, Path)
    assert "pylocuszoom" in str(cache_dir)
    assert "ensembl" in str(cache_dir)


def test_ensembl_cache_dir_honors_xdg_on_macos(tmp_path, monkeypatch):
    """macOS honors XDG_CACHE_HOME (previously hardcoded ~/.cache)."""
    from pylocuszoom._gene_cache import cache_root

    monkeypatch.setattr("sys.platform", "darwin")
    monkeypatch.setattr("os.path.exists", lambda _: False)
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path))

    assert cache_root("ensembl") == tmp_path / "pylocuszoom" / "ensembl"


def test_ensembl_cache_dir_routes_to_dbfs_on_databricks(monkeypatch):
    """Databricks routes the Ensembl cache under /dbfs (previously ~/.cache)."""
    from pylocuszoom._gene_cache import cache_root

    monkeypatch.setattr("sys.platform", "linux")
    monkeypatch.setattr("os.path.exists", lambda p: p == "/dbfs")
    monkeypatch.setattr(Path, "mkdir", lambda self, *a, **k: None)

    assert cache_root("ensembl") == Path("/dbfs/FileStore/reference_data/ensembl")


def test_cache_resolvers_share_base(tmp_path, monkeypatch):
    """Recombination and Ensembl caches resolve under one shared root."""
    from pylocuszoom._gene_cache import cache_root
    from pylocuszoom.recombination import get_default_data_dir

    monkeypatch.setattr("sys.platform", "darwin")
    monkeypatch.setattr("os.path.exists", lambda _: False)
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path))
    monkeypatch.setattr(Path, "mkdir", lambda self, *a, **k: None)

    assert get_default_data_dir() == tmp_path / "pylocuszoom" / "recombination_maps"
    assert cache_root("ensembl") == tmp_path / "pylocuszoom" / "ensembl"
    assert cache_root("ensembl").parent == get_default_data_dir().parent


def test_load_annotations_miss():
    """Test cache miss returns None."""
    from pylocuszoom._gene_cache import load_annotations

    with tempfile.TemporaryDirectory() as tmpdir:
        result = load_annotations(
            Path(tmpdir),
            "homo_sapiens",
            "13",
            32000000,
            33000000,
        )
        assert result is None


def test_save_and_load_cached_annotations():
    """Test saving and loading a cached entry using CSV."""
    from pylocuszoom._gene_cache import load_annotations, save_annotations
    from pylocuszoom._gene_source import GeneAnnotations

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

        save_annotations(
            GeneAnnotations(df, pd.DataFrame(columns=list(EXON_COLUMNS))),
            cache_dir,
            "homo_sapiens",
            "13",
            32000000,
            33000000,
        )

        # Verify CSV files created (not parquet), one per frame
        csv_files = list(cache_dir.glob("**/*.csv"))
        assert len(csv_files) == 2

        loaded = load_annotations(
            cache_dir,
            "homo_sapiens",
            "13",
            32000000,
            33000000,
        )

        assert loaded is not None
        assert len(loaded.genes) == 2
        assert loaded.exons.empty
        # Sort for deterministic comparison
        loaded_sorted = loaded.genes.sort_values("start")
        assert loaded_sorted["gene_name"].tolist() == ["BRCA2", "TEST"]


# --- get_genes_for_build tests ---


def test_get_genes_for_build_uses_cache():
    """Test that get_genes_for_build uses cache when available."""
    from pylocuszoom._gene_cache import save_annotations
    from pylocuszoom._gene_source import GeneAnnotations
    from pylocuszoom.ensembl import ensembl_source
    from pylocuszoom.reference_genes import get_genes_for_build

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
        save_annotations(
            GeneAnnotations(cached_df, pd.DataFrame(columns=list(EXON_COLUMNS))),
            cache_dir,
            "homo_sapiens",
            "13",
            32000000,
            33000000,
        )

        # Mock the API call - should NOT be called due to cache hit
        with patch("pylocuszoom._http.requests.get") as mock_get:
            result = get_genes_for_build(
                ensembl_source("human"),
                "13",
                32000000,
                33000000,
                cache_dir=cache_dir,
            )

            # API should not have been called
            mock_get.assert_not_called()

        assert len(result.genes) == 1
        assert result.genes["gene_name"].iloc[0] == "CACHED_GENE"


def test_get_genes_for_build_fetches_and_caches():
    """Test that get_genes_for_build fetches from API and caches result."""
    from pylocuszoom.ensembl import ensembl_source
    from pylocuszoom.reference_genes import get_genes_for_build

    mock_response = ok_response(
        [
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
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        cache_dir = Path(tmpdir)

        with patch("pylocuszoom._http.requests.get", return_value=mock_response):
            result = get_genes_for_build(
                ensembl_source("human"),
                "13",
                32000000,
                33000000,
                cache_dir=cache_dir,
            )

        assert len(result.genes) == 1
        assert result.genes["gene_name"].iloc[0] == "BRCA2"

        # Verify cache files were created (CSV, not parquet)
        csv_files = list(cache_dir.glob("**/*.csv"))
        assert len(csv_files) == 2


def test_get_genes_for_build_returns_exons():
    """One request serves both frames, and both are cached."""
    from pylocuszoom.ensembl import ensembl_source
    from pylocuszoom.reference_genes import get_genes_for_build

    mock_response = ok_response(gene_transcript_exon_payload())

    with tempfile.TemporaryDirectory() as tmpdir:
        cache_dir = Path(tmpdir)

        with patch("pylocuszoom._http.requests.get", return_value=mock_response):
            genes_df, exons_df = get_genes_for_build(
                ensembl_source("human"),
                "1",
                1000000,
                2000000,
                cache_dir=cache_dir,
            )

        assert len(genes_df) == 1
        assert len(exons_df) == 3
        assert exons_df["gene_name"].unique().tolist() == ["NFATC1"]

        with patch("pylocuszoom._http.requests.get") as mock_get:
            cached = get_genes_for_build(
                ensembl_source("human"),
                "1",
                1000000,
                2000000,
                cache_dir=cache_dir,
            )

        assert not mock_get.called, "a cached entry must not refetch its exons"
        assert len(cached.exons) == 3


# --- clear_gene_cache tests ---


def test_clear_gene_cache():
    """Test clearing the Ensembl cache."""
    from pylocuszoom._gene_cache import save_annotations
    from pylocuszoom._gene_source import GeneAnnotations
    from pylocuszoom.reference_genes import clear_gene_cache

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
        entry = GeneAnnotations(df, pd.DataFrame(columns=list(EXON_COLUMNS)))
        save_annotations(entry, cache_dir, "homo_sapiens", "1", 100, 200)
        save_annotations(entry, cache_dir, "mus_musculus", "1", 100, 200)

        # Two entries of two frames each
        csv_files = list(cache_dir.glob("**/*.csv"))
        assert len(csv_files) == 4

        # Clear cache
        deleted = clear_gene_cache("ensembl", cache_dir)

        assert deleted == 4
        assert len(list(cache_dir.glob("**/*.csv"))) == 0


def test_clear_gene_cache_species_specific():
    """Test clearing cache for specific species only."""
    from pylocuszoom._gene_cache import save_annotations
    from pylocuszoom._gene_source import GeneAnnotations
    from pylocuszoom.reference_genes import clear_gene_cache

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
        entry = GeneAnnotations(df, pd.DataFrame(columns=list(EXON_COLUMNS)))
        save_annotations(entry, cache_dir, "homo_sapiens", "1", 100, 200)
        save_annotations(entry, cache_dir, "mus_musculus", "1", 100, 200)

        # Clear only human cache
        deleted = clear_gene_cache("ensembl", cache_dir, "homo_sapiens")

        assert deleted == 2
        # Mouse cache should still exist
        assert len(list(cache_dir.glob("**/*.csv"))) == 2


class TestPathTraversalProtection:
    """Tests that malicious species names cannot escape the cache directory."""

    def test_path_traversal_species_raises_validation_error(self, tmp_path):
        """Species name with '../' should raise ValidationError."""
        from pylocuszoom._gene_cache import safe_species_dir
        from pylocuszoom.utils import ValidationError

        with pytest.raises(ValidationError, match="Invalid species name"):
            safe_species_dir(tmp_path, "../../etc")

    def test_absolute_path_species_raises_validation_error(self, tmp_path):
        """Absolute path as species should raise ValidationError."""
        from pylocuszoom._gene_cache import safe_species_dir
        from pylocuszoom.utils import ValidationError

        with pytest.raises(ValidationError, match="Invalid species name"):
            safe_species_dir(tmp_path, "/etc/passwd")

    def test_safe_species_stays_within_cache(self, tmp_path):
        """Normal species name should resolve within cache_dir."""
        from pylocuszoom._gene_cache import safe_species_dir

        result = safe_species_dir(tmp_path, "canis_lupus_familiaris")
        assert result.is_relative_to(tmp_path)

    def test_load_annotations_rejects_traversal(self, tmp_path):
        """load_annotations should reject path traversal species."""
        from pylocuszoom._gene_cache import load_annotations
        from pylocuszoom.utils import ValidationError

        with pytest.raises(ValidationError, match="Invalid species name"):
            load_annotations(tmp_path, "../../etc/passwd", "1", 1, 100)

    def test_save_annotations_rejects_traversal(self, tmp_path):
        """save_annotations should reject path traversal species."""
        import pandas as pd

        from pylocuszoom._gene_cache import save_annotations
        from pylocuszoom._gene_source import GeneAnnotations
        from pylocuszoom.utils import ValidationError

        annotations = GeneAnnotations(pd.DataFrame({"gene": ["BRCA1"]}), pd.DataFrame())
        with pytest.raises(ValidationError, match="Invalid species name"):
            save_annotations(annotations, tmp_path, "../../etc/passwd", "1", 1, 100)

    def test_clear_gene_cache_rejects_traversal(self, tmp_path):
        """clear_gene_cache should reject path traversal species."""
        from pylocuszoom.reference_genes import clear_gene_cache
        from pylocuszoom.utils import ValidationError

        with pytest.raises(ValidationError, match="Invalid species name"):
            clear_gene_cache("ensembl", tmp_path, "../../etc/passwd")


class TestEmptyResultCaching:
    """A cached empty region must reload, and a failed fetch must not be cached.

    Both bugs surfaced as ``pandas.errors.EmptyDataError`` on the second call:
    an empty DataFrame serialises to a one-byte file that cannot be parsed back.
    """

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
        from pylocuszoom.ensembl import ensembl_source
        from pylocuszoom.reference_genes import get_genes_for_build

        with patch(
            "pylocuszoom._http.requests.get", return_value=ok_response([])
        ) as mock_get:
            first = get_genes_for_build(
                ensembl_source("human"), "1", 1_000_000, 1_100_000, cache_dir=tmp_path
            )
            assert mock_get.call_count == 1

            second = get_genes_for_build(
                ensembl_source("human"), "1", 1_000_000, 1_100_000, cache_dir=tmp_path
            )

        assert first.genes.empty
        assert second.genes.empty
        assert mock_get.call_count == 1, "second call must be served from cache"
        assert list(second.genes.columns) == list(first.genes.columns)
        assert "gene_name" in second.genes.columns

    def test_failed_fetch_is_not_cached(self, tmp_path):
        """An API outage must not poison the cache for the region."""
        import requests

        from pylocuszoom.ensembl import ensembl_source
        from pylocuszoom.exceptions import EnsemblAPIError
        from pylocuszoom.reference_genes import get_genes_for_build

        with (
            patch("pylocuszoom._http.time.sleep"),
            patch(
                "pylocuszoom._http.requests.get",
                side_effect=requests.exceptions.ConnectionError("network down"),
            ),
        ):
            with pytest.raises(EnsemblAPIError):
                get_genes_for_build(
                    ensembl_source("human"),
                    "1",
                    1_000_000,
                    1_100_000,
                    cache_dir=tmp_path,
                )

        assert not list(tmp_path.rglob("genes_*.csv")), (
            "a failed fetch must leave no cache file"
        )

        with patch(
            "pylocuszoom._http.requests.get",
            return_value=ok_response(self._gene_payload()),
        ) as mock_get:
            after_recovery = get_genes_for_build(
                ensembl_source("human"), "1", 1_000_000, 1_100_000, cache_dir=tmp_path
            )

        assert mock_get.called, "recovery must retry the API, not reuse a failed result"
        assert len(after_recovery.genes) == 1
        assert after_recovery.genes["gene_name"].iloc[0] == "BRCA2"

    def test_zero_byte_cache_file_is_ignored(self, tmp_path):
        """A cache file poisoned by an older release is treated as a miss."""
        from pylocuszoom._gene_cache import cache_key, load_annotations

        species_dir = tmp_path / "homo_sapiens"
        species_dir.mkdir()
        key = cache_key("homo_sapiens", "1", 1_000_000, 1_100_000)
        (species_dir / f"genes_{key}.csv").write_text("\n")
        (species_dir / f"exons_{key}.csv").write_text("\n")

        assert (
            load_annotations(tmp_path, "homo_sapiens", "1", 1_000_000, 1_100_000)
            is None
        )


class TestAssemblyMismatch:
    """Ensembl serves one assembly per species and ignores coord_system_version.

    Asking for canine genes while working in CanFam3.1 returns ROS_Cfam_1.0
    coordinates with no error, shifting ATP9B on chr1 from 1,136,865 to
    938,796. The mismatch has to be loud, and the cache must not mix builds.
    """

    def test_assembly_token_folds_synonyms(self):
        """Equivalent spellings of one assembly compare equal."""
        from pylocuszoom.utils import assembly_token

        assert assembly_token("CanFam4.0") == assembly_token("UU_Cfam_GSD_1.0")
        assert assembly_token("CanFam3.1") == assembly_token("canfam3")
        assert assembly_token("hg38") == assembly_token("GRCh38")
        assert assembly_token("CanFam3.1") != assembly_token("ROS_Cfam_1.0")

    def test_two_builds_do_not_share_a_cache_entry(self, tmp_path):
        """The same region under two builds must not share a cache entry."""
        from pylocuszoom.ensembl import ensembl_source
        from pylocuszoom.reference_genes import get_genes_for_build

        with patch(
            "pylocuszoom._http.requests.get",
            return_value=ok_response(ros_cfam_gene_payload()),
        ) as mock_get:
            for build in ("canfam3.1", "canfam4"):
                with pytest.warns(UserWarning, match="ROS_Cfam_1.0"):
                    get_genes_for_build(
                        ensembl_source("canine", build),
                        "1",
                        1,
                        100,
                        cache_dir=tmp_path,
                    )

        assert mock_get.call_count == 2, "each build must fetch its own entry"

    def test_gene_record_carries_assembly(self):
        """Every gene row records the assembly Ensembl served it on."""
        from pylocuszoom.ensembl import fetch_overlap_frames

        with patch(
            "pylocuszoom._http.requests.get",
            return_value=ok_response(ros_cfam_gene_payload()),
        ):
            genes, _ = fetch_overlap_frames("canine", "1", 900_000, 1_200_000)

        assert genes["assembly"].tolist() == ["ROS_Cfam_1.0"]

    def test_exon_record_carries_assembly(self):
        """Exon rows record their assembly too, not just genes."""
        from pylocuszoom.ensembl import fetch_overlap_frames

        gene = ros_cfam_gene_payload()[0]
        transcript = dict(
            gene,
            feature_type="transcript",
            id="ENSCAFT00000001",
            Parent=gene["id"],
        )
        exon = dict(
            gene,
            feature_type="exon",
            id="ENSCAFE00000001",
            Parent=transcript["id"],
        )
        payload = [gene, transcript, exon]
        with patch("pylocuszoom._http.requests.get", return_value=ok_response(payload)):
            with pytest.warns(UserWarning, match="ROS_Cfam_1.0"):
                _, exons = fetch_overlap_frames(
                    "canine", "1", 900_000, 1_200_000, genome_build="canfam3.1"
                )

        assert exons["assembly"].tolist() == ["ROS_Cfam_1.0"]

    def test_fetch_warns_when_assembly_differs(self):
        """A CanFam3.1 caller is told the genes came back on ROS_Cfam_1.0."""
        from pylocuszoom.ensembl import fetch_overlap_frames

        with patch(
            "pylocuszoom._http.requests.get",
            return_value=ok_response(ros_cfam_gene_payload()),
        ):
            with pytest.warns(UserWarning, match="ROS_Cfam_1.0"):
                fetch_overlap_frames(
                    "canine", "1", 900_000, 1_200_000, genome_build="canfam3.1"
                )

    def test_fetch_silent_when_assembly_matches(self, recwarn):
        """No warning when the caller already works in Ensembl's assembly."""
        from pylocuszoom.ensembl import fetch_overlap_frames

        with patch(
            "pylocuszoom._http.requests.get",
            return_value=ok_response(ros_cfam_gene_payload()),
        ):
            fetch_overlap_frames(
                "canine", "1", 900_000, 1_200_000, genome_build="ROS_Cfam_1.0"
            )

        assert [w for w in recwarn.list if w.category is UserWarning] == []

    def test_cache_hit_still_warns(self, tmp_path):
        """Reloading from cache in a fresh session repeats the warning."""
        from pylocuszoom.ensembl import ensembl_source
        from pylocuszoom.reference_genes import get_genes_for_build

        with patch(
            "pylocuszoom._http.requests.get",
            return_value=ok_response(ros_cfam_gene_payload()),
        ) as mock_get:
            source = ensembl_source("canine", "canfam3.1")
            with pytest.warns(UserWarning, match="ROS_Cfam_1.0"):
                get_genes_for_build(source, "1", 900_000, 1_200_000, cache_dir=tmp_path)
            with pytest.warns(UserWarning, match="ROS_Cfam_1.0"):
                cached = get_genes_for_build(
                    source, "1", 900_000, 1_200_000, cache_dir=tmp_path
                )

        assert mock_get.call_count == 1, "second call must be served from cache"
        assert cached.genes["assembly"].tolist() == ["ROS_Cfam_1.0"]
