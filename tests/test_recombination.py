"""Tests for recombination rate overlay module."""

import io
import tarfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest
import requests

from pylocuszoom.recombination import (
    RECOMB_COLOR,
    _normalize_build,
    download_canine_recombination_maps,
    download_liftover_chain,
    ensure_recomb_maps,
    get_default_data_dir,
    get_recombination_rate_for_region,
    liftover_recombination_map,
    load_recombination_map,
)


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
        mock_download.assert_called_once()


class TestLiftoverRecombinationMap:
    """Tests for liftover_recombination_map function."""

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

        # Create 39 mock files (38 autosomes + X)
        for i in range(1, 39):
            (tmp_path / f"chr{i}_recomb.tsv").touch()
        (tmp_path / "chrX_recomb.tsv").touch()

        result = download_canine_recombination_maps(force=False)
        assert result == tmp_path


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
        # Create 39 fake chromosome files
        for i in range(1, 40):
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
    def test_ensure_recomb_maps_handles_network_error(
        self, mock_get_dir, mock_download, tmp_path
    ):
        """Test that ensure_recomb_maps returns None on network errors."""
        mock_get_dir.return_value = tmp_path / "recomb_data"
        mock_download.side_effect = requests.RequestException("Network error")

        result = ensure_recomb_maps(species="canine")
        assert result is None

    @patch("pylocuszoom.recombination.download_canine_recombination_maps")
    @patch("pylocuszoom.recombination.get_default_data_dir")
    def test_ensure_recomb_maps_handles_io_error(
        self, mock_get_dir, mock_download, tmp_path
    ):
        """Test that ensure_recomb_maps returns None on IO errors."""
        mock_get_dir.return_value = tmp_path / "recomb_data"
        mock_download.side_effect = IOError("Disk full")

        result = ensure_recomb_maps(species="canine")
        assert result is None

    @patch("pylocuszoom.recombination.download_canine_recombination_maps")
    @patch("pylocuszoom.recombination.get_default_data_dir")
    def test_ensure_recomb_maps_handles_tar_error(
        self, mock_get_dir, mock_download, tmp_path
    ):
        """Test that ensure_recomb_maps returns None on corrupt archive."""
        mock_get_dir.return_value = tmp_path / "recomb_data"
        mock_download.side_effect = tarfile.TarError("Corrupt archive")

        result = ensure_recomb_maps(species="canine")
        assert result is None

    @pytest.mark.parametrize(
        "exception,expected_msg",
        [
            (requests.RequestException("Connection refused"), "recombination maps"),
            (tarfile.TarError("Bad header"), "recombination maps"),
            (OSError("Permission denied"), "recombination maps"),
        ],
    )
    @patch("pylocuszoom.recombination.download_canine_recombination_maps")
    @patch("pylocuszoom.recombination.get_default_data_dir")
    def test_ensure_recomb_maps_error_message(
        self, mock_get_dir, mock_download, tmp_path, exception, expected_msg
    ):
        """Download errors produce log message with context."""
        mock_get_dir.return_value = tmp_path / "recomb_data"
        mock_download.side_effect = exception

        log_capture = io.StringIO()
        from pylocuszoom.logging import logger as plz_logger

        plz_logger.enable("ERROR", sink=log_capture)
        try:
            ensure_recomb_maps(species="canine")
        finally:
            plz_logger.enable("INFO")

        assert expected_msg in log_capture.getvalue()


class TestDownloadLiftoverChainHTTPError:
    """Tests for HTTPError handling in download_liftover_chain."""

    @patch("pylocuszoom.recombination._download_with_progress")
    def test_http_error_includes_url(self, mock_download, tmp_path, monkeypatch):
        """HTTPError should include the URL in the error message."""
        monkeypatch.setattr(
            "pylocuszoom.recombination.get_default_data_dir", lambda: tmp_path
        )
        mock_download.side_effect = requests.HTTPError("404 Client Error")

        with pytest.raises(requests.HTTPError, match="liftover chain file"):
            download_liftover_chain(force=True)

    @patch("pylocuszoom.recombination._download_with_progress")
    def test_http_error_preserves_original(self, mock_download, tmp_path, monkeypatch):
        """HTTPError should chain the original exception."""
        monkeypatch.setattr(
            "pylocuszoom.recombination.get_default_data_dir", lambda: tmp_path
        )
        original = requests.HTTPError("404 Client Error")
        mock_download.side_effect = original

        with pytest.raises(requests.HTTPError) as exc_info:
            download_liftover_chain(force=True)

        assert exc_info.value.__cause__ is original


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
                with pytest.raises(RuntimeError, match="Could not find chromosome"):
                    download_canine_recombination_maps(
                        output_dir=str(tmp_path / "output"), force=True
                    )
        finally:
            plz_logger.enable("INFO")

        assert "Skipping unsafe path in archive" in log_capture.getvalue()
