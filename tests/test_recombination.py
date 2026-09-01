"""Tests for recombination rate overlay module."""

import io
import tarfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest
import requests

from pylocuszoom._liftover import InMemoryLifter, liftover_positions
from pylocuszoom.exceptions import DataDownloadError
from pylocuszoom.recombination import (
    RECOMB_COLOR,
    _download_with_progress,
    _normalize_build,
    _publish_map_generation,
    download_canine_recombination_maps,
    download_liftover_chain,
    ensure_recomb_header,
    ensure_recomb_maps,
    get_default_data_dir,
    get_recombination_rate_for_region,
    liftover_recombination_map,
    load_recombination_map,
)


class TestEnsureRecombHeader:
    """Pure header-detection tests (no tarballs, no download)."""

    def test_prepends_header_when_first_token_numeric(self):
        content = "1\t1000\t0.5\t0.1\n1\t2000\t0.6\t0.2\n"
        result = ensure_recomb_header(content, "chr1_recomb.tsv")
        assert result == "chr\tpos\trate\tcM\n" + content

    def test_keeps_content_when_known_header_present(self):
        content = "chr\tpos\trate\tcM\n1\t1000\t0.5\t0.1\n"
        assert ensure_recomb_header(content, "chr1_recomb.tsv") == content

    def test_accepts_hash_prefixed_and_alternate_header_names(self):
        for header in ("#chrom", "position", "BP", "chromosome"):
            content = f"{header}\tx\ty\tz\n1\t2\t3\t4\n"
            assert ensure_recomb_header(content, "f.tsv") == content

    def test_rejects_unrecognised_first_token(self):
        content = "<html><body>404 Not Found</body></html>\n"
        with pytest.raises(DataDownloadError, match="refusing to treat as header"):
            ensure_recomb_header(content, "chr1_recomb.tsv")


class TestGetDefaultDataDir:
    """Tests for get_default_data_dir function."""

    def test_returns_path_object(self):
        """Should return a Path object."""
        result = get_default_data_dir()
        assert isinstance(result, Path)

    def test_path_ends_with_recombination_maps(self):
        """Path should end with recombination_maps directory."""
        result = get_default_data_dir()
        assert result.name == "recombination_maps"

    @patch.dict("os.environ", {"DATABRICKS_RUNTIME_VERSION": "12.2"})
    @patch("os.path.exists")
    def test_uses_dbfs_on_databricks(self, mock_exists):
        """Should use /dbfs path when on Databricks."""
        mock_exists.return_value = True
        result = get_default_data_dir()
        assert "/dbfs" in str(result)


class TestRecombColor:
    """Tests for recombination color constant."""

    def test_uses_correct_color(self):
        """Should use the defined recombination color."""
        assert RECOMB_COLOR == "#7FCDFF"  # Light blue


class TestLoadRecombinationMap:
    """Tests for load_recombination_map function."""

    def test_raises_for_missing_file(self, tmp_path):
        """Should raise FileNotFoundError when map file doesn't exist."""
        with pytest.raises(FileNotFoundError, match="Recombination map not found"):
            load_recombination_map(chrom=1, data_dir=str(tmp_path))

    def test_loads_valid_map_file(self, tmp_path):
        """Should load and parse valid recombination map file."""
        # Create test file
        map_content = "chr\tpos\trate\tcM\n1\t1000\t0.5\t0.001\n1\t5000\t1.2\t0.005\n"
        map_file = tmp_path / "chr1_recomb.tsv"
        map_file.write_text(map_content)

        result = load_recombination_map(chrom=1, data_dir=str(tmp_path))

        assert len(result) == 2
        assert "pos" in result.columns
        assert "rate" in result.columns
        assert result["pos"].iloc[0] == 1000
        assert result["rate"].iloc[0] == 0.5

    def test_handles_chr_prefix_in_argument(self, tmp_path):
        """Should handle 'chr' prefix in chromosome argument."""
        map_content = "chr\tpos\trate\tcM\n1\t1000\t0.5\t0.001\n"
        map_file = tmp_path / "chr1_recomb.tsv"
        map_file.write_text(map_content)

        # Should work with "chr1" argument
        result = load_recombination_map(chrom="chr1", data_dir=str(tmp_path))
        assert len(result) == 1

    def test_non_numeric_values_produce_warning(self, tmp_path):
        """Non-numeric values in pos/rate should produce a warning."""
        map_content = "chr\tpos\trate\tcM\n1\t1000\t0.5\t0.001\n1\tBAD\t1.2\t0.005\n"
        map_file = tmp_path / "chr1_recomb.tsv"
        map_file.write_text(map_content)

        log_capture = io.StringIO()
        from pylocuszoom.logging import logger as plz_logger

        plz_logger.enable("WARNING", sink=log_capture)
        try:
            result = load_recombination_map(chrom=1, data_dir=str(tmp_path))
        finally:
            plz_logger.enable("INFO")

        # Only valid row should remain
        assert len(result) == 1
        assert result["pos"].iloc[0] == 1000

        # Warning should mention non-numeric values
        log_output = log_capture.getvalue()
        assert "non-numeric values" in log_output
        assert "chr1" in log_output


class TestGetRecombinationRateForRegion:
    """Tests for get_recombination_rate_for_region function."""

    def test_filters_to_region(self, tmp_path):
        """Should return only data within specified region."""
        # Create test file with data spanning 1000-10000
        map_content = (
            "chr\tpos\trate\tcM\n"
            "1\t1000\t0.5\t0.001\n"
            "1\t3000\t1.2\t0.003\n"
            "1\t5000\t2.0\t0.005\n"
            "1\t7000\t1.5\t0.007\n"
            "1\t10000\t0.8\t0.010\n"
        )
        map_file = tmp_path / "chr1_recomb.tsv"
        map_file.write_text(map_content)

        result = get_recombination_rate_for_region(
            chrom=1, start=2000, end=6000, data_dir=str(tmp_path)
        )

        # Should only include positions 3000 and 5000
        assert len(result) == 2
        assert 3000 in result["pos"].values
        assert 5000 in result["pos"].values
        assert 1000 not in result["pos"].values

    def test_returns_only_pos_and_rate_columns(self, tmp_path):
        """Should return DataFrame with only pos and rate columns."""
        map_content = "chr\tpos\trate\tcM\n1\t1000\t0.5\t0.001\n"
        map_file = tmp_path / "chr1_recomb.tsv"
        map_file.write_text(map_content)

        result = get_recombination_rate_for_region(
            chrom=1, start=0, end=2000, data_dir=str(tmp_path)
        )

        assert list(result.columns) == ["pos", "rate"]


class TestNormalizeBuild:
    """Tests for _normalize_build function."""

    def test_none_returns_none(self):
        """None input returns None."""
        assert _normalize_build(None) is None

    def test_canfam4_variations(self):
        """Various CanFam4 names normalize correctly."""
        assert _normalize_build("canfam4") == "canfam4"
        assert _normalize_build("CanFam4.0") == "canfam4"
        assert _normalize_build("UU_Cfam_GSD_1.0") == "canfam4"

    def test_canfam3_variations(self):
        """Various CanFam3 names normalize correctly."""
        assert _normalize_build("canfam3") == "canfam3"
        assert _normalize_build("CanFam3.1") == "canfam3"

    def test_unknown_build_lowercase(self):
        """Unknown build returns lowercase."""
        assert _normalize_build("hg38") == "hg38"

    def test_unknown_build_does_not_raise(self):
        """Unknown build returns lowercase without raising."""
        result = _normalize_build("FelCat9")
        assert result == "felcat9"
        assert isinstance(result, str)


class TestDownloadLiftoverChain:
    """Tests for download_liftover_chain function."""

    def test_returns_existing_file(self, tmp_path, monkeypatch):
        """Returns existing chain file without re-downloading."""
        # Create mock chain file
        monkeypatch.setattr(
            "pylocuszoom.recombination.get_default_data_dir", lambda: tmp_path
        )
        chain_file = tmp_path / "canFam3ToCanFam4.over.chain.gz"
        chain_file.write_bytes(b"mock chain data")

        result = download_liftover_chain(force=False)
        assert result == chain_file

    @patch("pylocuszoom.recombination._download_with_progress")
    def test_downloads_when_missing(self, mock_download, tmp_path, monkeypatch):
        """Downloads chain file when not present."""
        monkeypatch.setattr(
            "pylocuszoom.recombination.get_default_data_dir", lambda: tmp_path
        )

        # Mock the download to create the file
        def create_file(url, dest, desc):
            dest.write_bytes(b"mock chain data")

        mock_download.side_effect = create_file

        result = download_liftover_chain(force=False)
        mock_download.assert_called_once()
        assert result.exists()

    @patch("pylocuszoom.recombination._download_with_progress")
    def test_force_redownload(self, mock_download, tmp_path, monkeypatch):
        """Force=True re-downloads even if file exists."""
        monkeypatch.setattr(
            "pylocuszoom.recombination.get_default_data_dir", lambda: tmp_path
        )

        # Create existing file
        chain_file = tmp_path / "canFam3ToCanFam4.over.chain.gz"
        chain_file.write_bytes(b"old data")

        def create_file(url, dest, desc):
            dest.write_bytes(b"new data")

        mock_download.side_effect = create_file

        download_liftover_chain(force=True)

        # Observable behaviour: the existing file's contents were
        # overwritten with the freshly-fetched bytes.
        assert chain_file.read_bytes() == b"new data", (
            "force=True must replace existing chain file contents"
        )


class TestLiftoverPositions:
    """Pure liftover math via an in-memory lifter (no pyliftover, no network)."""

    def test_lifts_mapped_positions_with_chrom_arg(self):
        df = pd.DataFrame({"pos": [1000, 2000], "rate": [0.5, 0.6]})
        lifter = InMemoryLifter({("chr1", 1000): 1100, ("chr1", 2000): 2200})
        result = liftover_positions(df, lifter, chrom=1)
        assert list(result["pos"]) == [1100, 2200]
        assert list(result["rate"]) == [0.5, 0.6]

    def test_drops_unmapped_positions(self):
        df = pd.DataFrame({"pos": [1000, 2000], "rate": [0.5, 0.6]})
        lifter = InMemoryLifter({("chr1", 1000): 1100})  # 2000 fails to map
        result = liftover_positions(df, lifter, chrom=1)
        assert list(result["pos"]) == [1100]
        assert list(result["rate"]) == [0.5]

    def test_uses_chr_column_when_present(self):
        df = pd.DataFrame({"chr": [1, 2], "pos": [1000, 3000], "rate": [0.1, 0.2]})
        lifter = InMemoryLifter({("chr1", 1000): 1100, ("chr2", 3000): 3300})
        result = liftover_positions(df, lifter)
        assert set(result["pos"]) == {1100, 3300}

    def test_result_sorted_by_lifted_position(self):
        df = pd.DataFrame({"pos": [1000, 2000], "rate": [0.5, 0.6]})
        lifter = InMemoryLifter({("chr1", 1000): 5000, ("chr1", 2000): 1000})
        result = liftover_positions(df, lifter, chrom=1)
        assert list(result["pos"]) == [1000, 5000]

    def test_takes_first_of_multiple_mappings(self):
        class MultiLifter:
            def convert(self, chrom, pos):
                return [(chrom, 1111, "+", 0), (chrom, 9999, "+", 0)]

        df = pd.DataFrame({"pos": [1000], "rate": [0.5]})
        result = liftover_positions(df, MultiLifter(), chrom=1)
        assert list(result["pos"]) == [1111]

    def test_requires_chr_column_or_chrom(self):
        df = pd.DataFrame({"pos": [1000], "rate": [0.5]})
        with pytest.raises(ValueError, match="Either 'chr' column or chrom"):
            liftover_positions(df, InMemoryLifter({}))


class TestLiftoverRecombinationMap:
    """Tests for liftover_recombination_map function."""

    def test_raises_import_error_when_pyliftover_missing(self):
        """Should raise ImportError with install instructions when pyliftover unavailable."""
        import sys

        # Temporarily hide pyliftover from imports
        original = sys.modules.get("pyliftover")
        sys.modules["pyliftover"] = None  # type: ignore[assignment]
        # Force re-import of the function's local import
        try:
            # The import happens inside liftover_recombination_map, so we need
            # to call it. Create minimal valid input.
            df = pd.DataFrame({"pos": [1000000], "rate": [0.5]})
            with pytest.raises(ImportError, match="pip install pyliftover"):
                liftover_recombination_map(df, chrom=1)
        finally:
            if original is not None:
                sys.modules["pyliftover"] = original
            else:
                sys.modules.pop("pyliftover", None)

    @patch("pylocuszoom.recombination.download_liftover_chain")
    @patch("pyliftover.LiftOver")
    def test_lifts_positions(self, mock_liftover_class, mock_download, tmp_path):
        """Successfully lifts over positions."""
        # Mock chain download
        chain_file = tmp_path / "chain.gz"
        chain_file.touch()
        mock_download.return_value = chain_file

        # Mock LiftOver
        mock_lo = MagicMock()
        mock_lo.convert_coordinate.side_effect = [
            [("chr1", 1000100, "+", 1)],  # First position maps
            [("chr1", 1500100, "+", 1)],  # Second position maps
        ]
        mock_liftover_class.return_value = mock_lo

        df = pd.DataFrame(
            {
                "pos": [1000000, 1500000],
                "rate": [0.5, 1.0],
            }
        )

        result = liftover_recombination_map(df, chrom=1)

        assert len(result) == 2
        assert result["pos"].iloc[0] == 1000100
        assert result["pos"].iloc[1] == 1500100

    @patch("pylocuszoom.recombination.download_liftover_chain")
    @patch("pyliftover.LiftOver")
    def test_drops_unmapped_positions(
        self, mock_liftover_class, mock_download, tmp_path
    ):
        """Positions that fail to map are dropped."""
        chain_file = tmp_path / "chain.gz"
        chain_file.touch()
        mock_download.return_value = chain_file

        mock_lo = MagicMock()
        mock_lo.convert_coordinate.side_effect = [
            [("chr1", 1000100, "+", 1)],  # Maps
            [],  # Fails to map
            [("chr1", 2000100, "+", 1)],  # Maps
        ]
        mock_liftover_class.return_value = mock_lo

        df = pd.DataFrame(
            {
                "pos": [1000000, 1500000, 2000000],
                "rate": [0.5, 1.0, 1.5],
            }
        )

        result = liftover_recombination_map(df, chrom=1)

        assert len(result) == 2
        assert 1500100 not in result["pos"].values

    @patch("pylocuszoom.recombination.download_liftover_chain")
    @patch("pyliftover.LiftOver")
    def test_uses_chr_column_if_present(
        self, mock_liftover_class, mock_download, tmp_path
    ):
        """Uses chr column from DataFrame if present."""
        chain_file = tmp_path / "chain.gz"
        chain_file.touch()
        mock_download.return_value = chain_file

        mock_lo = MagicMock()
        mock_lo.convert_coordinate.return_value = [("chr1", 1000100, "+", 1)]
        mock_liftover_class.return_value = mock_lo

        df = pd.DataFrame(
            {
                "chr": [1],
                "pos": [1000000],
                "rate": [0.5],
            }
        )

        liftover_recombination_map(df)  # No chrom argument needed

        # Should have used chr column
        mock_lo.convert_coordinate.assert_called()

    @patch("pylocuszoom.recombination.download_liftover_chain")
    @patch("pyliftover.LiftOver")
    def test_requires_chr_or_chrom_param(
        self, mock_liftover_class, mock_download, tmp_path
    ):
        """Raises ValueError if neither chr column nor chrom param."""
        chain_file = tmp_path / "chain.gz"
        chain_file.touch()
        mock_download.return_value = chain_file
        mock_liftover_class.return_value = MagicMock()

        df = pd.DataFrame(
            {
                "pos": [1000000],
                "rate": [0.5],
            }
        )

        with pytest.raises(ValueError, match="chr"):
            liftover_recombination_map(df)


class TestDownloadCanineRecombinationMaps:
    """Tests for download_canine_recombination_maps function."""

    def test_returns_existing_complete_data(self, tmp_path, monkeypatch):
        """Returns existing directory if all files present."""
        monkeypatch.setattr(
            "pylocuszoom.recombination.get_default_data_dir", lambda: tmp_path
        )

        # Create the complete 38-autosome map set.
        for i in range(1, 39):
            (tmp_path / f"chr{i}_recomb.tsv").touch()

        result = download_canine_recombination_maps(force=False)
        assert result == tmp_path

    @patch("pylocuszoom.recombination._download_with_progress")
    def test_rejects_wrong_39_file_manifest(self, mock_download, tmp_path):
        """A count of 39 files is not proof that the canine set is complete."""
        for i in range(1, 40):
            (tmp_path / f"chr{i}_recomb.tsv").touch()

        mock_download.side_effect = DataDownloadError("download attempted")

        with pytest.raises(DataDownloadError, match="download attempted"):
            download_canine_recombination_maps(output_dir=str(tmp_path), force=False)


class TestPublishMapGeneration:
    """Atomic publication tests for the recombination map set."""

    @staticmethod
    def _write_maps(path: Path, content: str, *, complete: bool = True) -> None:
        path.mkdir(parents=True, exist_ok=True)
        stop = 39 if complete else 38
        for chrom in range(1, stop):
            (path / f"chr{chrom}_recomb.tsv").write_text(content)

    def test_switches_complete_generations_and_preserves_ancillary_files(
        self, tmp_path
    ):
        output = tmp_path / "maps"
        self._write_maps(output, "old")
        (output / "canFam3ToCanFam4.over.chain.gz").write_text("chain")

        first_staging = tmp_path / "first-staging"
        self._write_maps(first_staging, "first")
        _publish_map_generation(first_staging, output)

        first_generation = output.resolve()
        assert output.is_symlink()
        assert (output / "chr1_recomb.tsv").read_text() == "first"
        assert (output / "canFam3ToCanFam4.over.chain.gz").read_text() == "chain"

        second_staging = tmp_path / "second-staging"
        self._write_maps(second_staging, "second")
        _publish_map_generation(second_staging, output)

        assert output.resolve() != first_generation
        assert (output / "chr1_recomb.tsv").read_text() == "second"
        assert not first_generation.exists()

    def test_incomplete_generation_leaves_active_maps_unchanged(self, tmp_path):
        output = tmp_path / "maps"
        self._write_maps(output, "old")
        staging = tmp_path / "incomplete-staging"
        self._write_maps(staging, "new", complete=False)

        with pytest.raises(DataDownloadError, match="complete canine map set"):
            _publish_map_generation(staging, output)

        assert not output.is_symlink()
        assert (output / "chr1_recomb.tsv").read_text() == "old"
        assert (output / "chr38_recomb.tsv").read_text() == "old"


class TestEnsureRecombMaps:
    """Tests for ensure_recomb_maps function."""

    @patch("pylocuszoom.recombination.download_canine_recombination_maps")
    @patch("pylocuszoom.recombination.get_default_data_dir")
    def test_ensure_recomb_maps_downloads_if_missing(
        self, mock_get_dir, mock_download, tmp_path
    ):
        """Test that ensure_recomb_maps triggers download when maps missing."""
        mock_get_dir.return_value = tmp_path / "recomb_data"
        mock_download.return_value = tmp_path / "recomb_data"

        result = ensure_recomb_maps(species="canine")

        mock_download.assert_called_once()
        assert result == tmp_path / "recomb_data"

    @patch("pylocuszoom.recombination.download_canine_recombination_maps")
    @patch("pylocuszoom.recombination.get_default_data_dir")
    def test_ensure_recomb_maps_skips_download_if_exists(
        self, mock_get_dir, mock_download, tmp_path
    ):
        """Test that ensure_recomb_maps skips download when maps exist."""
        data_dir = tmp_path / "recomb_data"
        data_dir.mkdir()
        # Create the complete 38-autosome map set.
        for i in range(1, 39):
            (data_dir / f"chr{i}_recomb.tsv").touch()

        mock_get_dir.return_value = data_dir

        result = ensure_recomb_maps(species="canine")

        mock_download.assert_not_called()
        assert result == data_dir

    def test_ensure_recomb_maps_non_canine_returns_none(self):
        """Test that ensure_recomb_maps returns None for non-canine species."""
        result = ensure_recomb_maps(species="human")
        assert result is None

    @patch("pylocuszoom.recombination.download_canine_recombination_maps")
    @patch("pylocuszoom.recombination.get_default_data_dir")
    def test_ensure_recomb_maps_warns_and_returns_none_on_download_error(
        self, mock_get_dir, mock_download, tmp_path
    ):
        mock_get_dir.return_value = tmp_path / "recomb_data"
        mock_download.side_effect = DataDownloadError("Network error")

        with pytest.warns(UserWarning, match="recombination maps.*Network error"):
            result = ensure_recomb_maps(species="canine")

        assert result is None

    @patch("pylocuszoom.recombination.download_canine_recombination_maps")
    @patch("pylocuszoom.recombination.get_default_data_dir")
    def test_ensure_recomb_maps_warns_and_returns_none_on_io_error(
        self, mock_get_dir, mock_download, tmp_path
    ):
        mock_get_dir.return_value = tmp_path / "recomb_data"
        mock_download.side_effect = OSError("Disk full")

        with pytest.warns(UserWarning, match="recombination maps.*Disk full"):
            result = ensure_recomb_maps(species="canine")

        assert result is None

    @patch("pylocuszoom.recombination._download_with_progress")
    def test_ensure_recomb_maps_warns_and_returns_none_on_corrupt_archive(
        self, mock_download, tmp_path
    ):
        def write_garbage(url, dest_path, desc=None):
            dest_path.write_bytes(b"not a gzip archive")

        mock_download.side_effect = write_garbage

        with pytest.warns(UserWarning, match="recombination maps"):
            result = ensure_recomb_maps(species="canine", data_dir=str(tmp_path / "x"))

        assert result is None


class TestDownloadWithProgress:
    """A download must never leave a truncated file at the destination."""

    @staticmethod
    def _streaming_response(chunks):
        response = MagicMock()
        response.headers = {"content-length": str(sum(len(c) for c in chunks))}
        response.iter_content.return_value = iter(chunks)
        return response

    def test_complete_download_lands_at_dest(self, tmp_path):
        dest = tmp_path / "file.gz"
        response = self._streaming_response([b"abcd", b"efgh"])

        with patch("pylocuszoom.recombination.requests.get", return_value=response):
            _download_with_progress("https://example.invalid/f", dest)

        assert dest.read_bytes() == b"abcdefgh"
        assert [p.name for p in tmp_path.iterdir()] == ["file.gz"]

    def test_interrupted_download_leaves_nothing_behind(self, tmp_path):
        dest = tmp_path / "file.gz"

        def chunks():
            yield b"abcd"
            raise requests.ConnectionError("connection reset")

        response = self._streaming_response([])
        response.iter_content.return_value = chunks()

        with (
            patch("pylocuszoom.recombination.requests.get", return_value=response),
            pytest.raises(DataDownloadError, match="example.invalid/f") as exc_info,
        ):
            _download_with_progress("https://example.invalid/f", dest)

        assert isinstance(exc_info.value.__cause__, requests.ConnectionError)
        assert not dest.exists()
        assert list(tmp_path.iterdir()) == []

    def test_http_error_raises_download_error(self, tmp_path):
        dest = tmp_path / "file.gz"
        response = self._streaming_response([])
        original = requests.HTTPError("404 Client Error")
        response.raise_for_status.side_effect = original

        with (
            patch("pylocuszoom.recombination.requests.get", return_value=response),
            pytest.raises(DataDownloadError, match="example.invalid/f") as exc_info,
        ):
            _download_with_progress("https://example.invalid/f", dest)

        assert exc_info.value.__cause__ is original
        assert list(tmp_path.iterdir()) == []


class TestDownloadLiftoverChainDownloadError:
    """download_liftover_chain surfaces failures as DataDownloadError."""

    @patch("pylocuszoom.recombination._download_with_progress")
    def test_download_error_propagates(self, mock_download, tmp_path, monkeypatch):
        monkeypatch.setattr(
            "pylocuszoom.recombination.get_default_data_dir", lambda: tmp_path
        )
        mock_download.side_effect = DataDownloadError("404 Client Error")

        with pytest.raises(DataDownloadError, match="404 Client Error"):
            download_liftover_chain(force=True)


class TestTarTraversalOSError:
    """Tests for OSError handling in tar path traversal check."""

    @patch("pylocuszoom.recombination._download_with_progress")
    def test_oserror_in_tar_member_path_is_caught(self, mock_download, tmp_path):
        """OSError during Path.resolve() in tar extraction is caught and skipped."""
        # Create a real tar.gz in memory with a normal file
        tar_buffer = io.BytesIO()
        with tarfile.open(fileobj=tar_buffer, mode="w:gz") as tar:
            # Add a dummy file
            info = tarfile.TarInfo(name="test_file.txt")
            info.size = 4
            tar.addfile(info, io.BytesIO(b"test"))
        tar_buffer.seek(0)

        # Mock download to write the tar to disk
        def write_tar(url, dest, desc):
            dest.write_bytes(tar_buffer.getvalue())

        mock_download.side_effect = write_tar

        # Patch Path.resolve to raise OSError for the member path
        original_resolve = Path.resolve

        def patched_resolve(self_path, *args, **kwargs):
            if "test_file.txt" in str(self_path):
                raise OSError("Invalid path")
            return original_resolve(self_path, *args, **kwargs)

        log_capture = io.StringIO()
        from pylocuszoom.logging import logger as plz_logger

        plz_logger.enable("WARNING", sink=log_capture)
        try:
            with patch.object(Path, "resolve", patched_resolve):
                with pytest.raises(
                    DataDownloadError, match="Could not find chromosome"
                ):
                    download_canine_recombination_maps(
                        output_dir=str(tmp_path / "output"), force=True
                    )
        finally:
            plz_logger.enable("INFO")

        assert "Skipping unsafe path in archive" in log_capture.getvalue()


class TestDownloadCanineRecombHeaderDetection:
    """Regression: download_canine_recombination_maps must reject unknown
    non-numeric first tokens (corrupted mirror / HTML error body) and
    accept all plausible header variants (case-insensitive, optional '#').
    """

    @staticmethod
    def _make_tarball(tar_path: Path, filename: str, content: str) -> None:
        """Create a minimal .tar.gz containing one chromosome map file."""
        with tarfile.open(tar_path, "w:gz") as tar:
            data = content.encode("utf-8")
            info = tarfile.TarInfo(name=filename)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))

    def _fake_download(self, filename: str, content: str):
        """Return a _download_with_progress mock that writes our fake tarball."""

        def side_effect(url, dest_path, desc=None):
            with tarfile.open(dest_path, "w:gz") as tar:
                for chrom in map(str, range(1, 39)):
                    map_filename = f"chr{chrom}.txt"
                    map_content = (
                        content
                        if map_filename == filename
                        else (f"chr\tpos\trate\tcM\n{chrom}\t1000\t0.5\t0.1\n")
                    )
                    data = map_content.encode("utf-8")
                    info = tarfile.TarInfo(name=map_filename)
                    info.size = len(data)
                    tar.addfile(info, io.BytesIO(data))

        return side_effect

    @patch("pylocuszoom.recombination._download_with_progress")
    def test_raises_on_html_corrupted_body(self, mock_download, tmp_path):
        """An HTML error body masquerading as a map is a download failure."""
        html_content = "<html><body>502 Bad Gateway</body></html>\n"
        mock_download.side_effect = self._fake_download("chr1.txt", html_content)

        with pytest.raises(DataDownloadError, match="Unrecognised first token"):
            download_canine_recombination_maps(tmp_path / "out")

    @patch("pylocuszoom.recombination._download_with_progress")
    def test_plot_warns_and_renders_without_overlay_on_corrupt_archive(
        self, mock_download, tmp_path
    ):
        """A corrupted mirror must be loud, but must not stop the plot."""
        from pylocuszoom import LocusZoomPlotter

        html_content = "<html><body>502 Bad Gateway</body></html>\n"
        mock_download.side_effect = self._fake_download("chr1.txt", html_content)
        plotter = LocusZoomPlotter(
            species="canine", recomb_data_dir=str(tmp_path / "out"), log_level=None
        )
        gwas = pd.DataFrame(
            {"ps": [1_000_000, 1_050_000], "p_wald": [0.5, 1e-6], "rs": ["a", "b"]}
        )

        with pytest.warns(UserWarning, match="Unrecognised first token"):
            fig = plotter.plot(gwas, chrom=1, start=900_000, end=1_200_000)

        assert fig is not None
        assert not (tmp_path / "out").exists()

    @patch("pylocuszoom.recombination._download_with_progress")
    def test_accepts_lowercase_chr_header(self, mock_download, tmp_path):
        """Pre-existing canonical header form must still be accepted."""
        content = "chr\tpos\trate\tcM\n1\t1000\t0.5\t0.1\n"
        mock_download.side_effect = self._fake_download("chr1.txt", content)

        result = download_canine_recombination_maps(tmp_path / "out")
        assert (result / "chr1_recomb.tsv").exists()

    @patch("pylocuszoom.recombination._download_with_progress")
    def test_accepts_capitalised_chromosome_header(self, mock_download, tmp_path):
        """'Chromosome' is used by several mirrors; must pass."""
        content = "Chromosome\tPosition\tRate\tcM\n1\t1000\t0.5\t0.1\n"
        mock_download.side_effect = self._fake_download("chr1.txt", content)

        result = download_canine_recombination_maps(tmp_path / "out")
        assert (result / "chr1_recomb.tsv").exists()

    @patch("pylocuszoom.recombination._download_with_progress")
    def test_accepts_hash_prefixed_header(self, mock_download, tmp_path):
        """Some maps use '#chr' as a commented header; must pass."""
        content = "#chr\tpos\trate\tcM\n1\t1000\t0.5\t0.1\n"
        mock_download.side_effect = self._fake_download("chr1.txt", content)

        result = download_canine_recombination_maps(tmp_path / "out")
        assert (result / "chr1_recomb.tsv").exists()

    @patch("pylocuszoom.recombination._download_with_progress")
    def test_accepts_numeric_first_token_prepends_header(self, mock_download, tmp_path):
        """Numeric first token means no header; one is prepended."""
        content = "1\t1000\t0.5\t0.1\n1\t2000\t0.6\t0.2\n"
        mock_download.side_effect = self._fake_download("chr1.txt", content)

        result = download_canine_recombination_maps(tmp_path / "out")
        out_file = result / "chr1_recomb.tsv"
        assert out_file.exists()
        assert out_file.read_text().startswith("chr\tpos\trate\tcM\n")
