"""Recombination maps: loading, region lookup, liftover and overlay status."""

import io
from pathlib import Path
from unittest.mock import MagicMock, Mock, patch

import pandas as pd
import pytest

from pylocuszoom._liftover import InMemoryLifter
from pylocuszoom.colors import RECOMB_COLOR
from pylocuszoom.exceptions import DataDownloadError
from pylocuszoom.recombination import (
    RecombStatus,
    download_liftover_chain,
    get_default_data_dir,
    get_recombination_rate_for_region,
    liftover_recombination_map,
    load_recombination_map,
    recomb_for_region,
)
from tests.conftest import write_canine_map_set


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

    def test_missing_map_names_the_function_that_downloads_it(self, tmp_path):
        """The remedy in the message has to be a function that exists."""
        with pytest.raises(FileNotFoundError, match=r"ensure_recomb_maps\(species="):
            load_recombination_map(1, species="canine", data_dir=str(tmp_path))

    def test_missing_map_for_a_species_without_maps_says_so(self, tmp_path):
        """A species with no built-in source is told to supply its own maps."""
        with pytest.raises(FileNotFoundError, match="no built-in recombination maps"):
            load_recombination_map(1, species="feline", data_dir=str(tmp_path))

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


class TestDownloadLiftoverChain:
    """Tests for download_liftover_chain function."""

    def test_returns_existing_file(self, tmp_path, monkeypatch):
        """Returns existing chain file without re-downloading."""
        # Create mock chain file
        monkeypatch.setattr("pylocuszoom.recombination.get_chain_dir", lambda: tmp_path)
        chain_file = tmp_path / "canFam3ToCanFam4.over.chain.gz"
        chain_file.write_bytes(b"mock chain data")

        result = download_liftover_chain(force=False)
        assert result == chain_file

    @patch("pylocuszoom.recombination.download_file")
    def test_downloads_when_missing(self, mock_download, tmp_path, monkeypatch):
        """Downloads chain file when not present."""
        monkeypatch.setattr("pylocuszoom.recombination.get_chain_dir", lambda: tmp_path)

        # Mock the download to create the file
        def create_file(url, dest, desc):
            dest.write_bytes(b"mock chain data")

        mock_download.side_effect = create_file

        result = download_liftover_chain(force=False)
        mock_download.assert_called_once()
        assert result.exists()

    @patch("pylocuszoom.recombination.download_file")
    def test_force_redownload(self, mock_download, tmp_path, monkeypatch):
        """Force=True re-downloads even if file exists."""
        monkeypatch.setattr("pylocuszoom.recombination.get_chain_dir", lambda: tmp_path)

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


class TestLiftoverRecombinationMap:
    """Tests for liftover_recombination_map function."""

    def test_raises_typed_error_when_pyliftover_missing(self):
        """A missing extra is OptionalDependencyMissing, not a bare ImportError."""
        import sys

        from pylocuszoom.exceptions import OptionalDependencyMissing

        # Temporarily hide pyliftover from imports
        original = sys.modules.get("pyliftover")
        sys.modules["pyliftover"] = None  # type: ignore[assignment]
        # Force re-import of the function's local import
        try:
            # The import happens inside liftover_recombination_map, so we need
            # to call it. Create minimal valid input.
            df = pd.DataFrame({"pos": [1000000], "rate": [0.5]})
            with pytest.raises(
                OptionalDependencyMissing, match="pip install pyliftover"
            ):
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

        # pyliftover hits are 0-based; the lifted 1-based position adds one.
        assert len(result) == 2
        assert result["pos"].iloc[0] == 1000101
        assert result["pos"].iloc[1] == 1500101

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


class TestRecombForRegion:
    """One status per way the overlay can be unavailable, and no warnings."""

    @staticmethod
    def _failing_download(exc):
        return patch(
            "pylocuszoom.recombination.download_recombination_maps",
            Mock(side_effect=exc),
        )

    @patch("pylocuszoom.recombination.get_default_data_dir")
    def test_a_download_failure_is_a_status_with_the_cause_in_it(
        self, mock_get_dir, tmp_path
    ):
        mock_get_dir.return_value = tmp_path / "recomb_data"

        with self._failing_download(DataDownloadError("Network error")):
            result = recomb_for_region(1, 1_000_000, 2_000_000, species="canine")

        assert result.status is RecombStatus.DOWNLOAD_FAILED
        assert "Network error" in result.detail
        assert result.frame is None

    @patch("pylocuszoom.recombination.get_default_data_dir")
    def test_an_os_error_is_also_a_download_failure(self, mock_get_dir, tmp_path):
        mock_get_dir.return_value = tmp_path / "recomb_data"

        with self._failing_download(OSError("Disk full")):
            result = recomb_for_region(1, 1_000_000, 2_000_000, species="canine")

        assert result.status is RecombStatus.DOWNLOAD_FAILED
        assert "Disk full" in result.detail

    def test_a_species_with_no_built_in_maps_says_so_instead_of_going_quiet(self):
        """This case used to reach the user as a debug log and nothing else."""
        result = recomb_for_region(1, 1_000_000, 2_000_000, species="human")

        assert result.status is RecombStatus.NO_MAPS_FOR_SPECIES
        assert "human" in result.detail

    def test_a_chromosome_the_map_set_does_not_cover_is_its_own_status(self, tmp_path):
        with patch(
            "pylocuszoom.recombination.get_default_data_dir", return_value=tmp_path
        ):
            result = recomb_for_region(
                41, 1_000_000, 2_000_000, species="canine", data_dir=str(tmp_path)
            )

        assert result.status is RecombStatus.NO_MAP_FOR_CHROMOSOME

    def test_a_region_that_resolves_carries_its_frame(self, tmp_path):
        (tmp_path / "chr1_recomb.tsv").write_text(
            "chr\tpos\trate\tcM\n1\t1500000\t0.5\t0.1\n"
        )

        with patch(
            "pylocuszoom.recombination.get_default_data_dir", return_value=tmp_path
        ):
            result = recomb_for_region(
                1, 1_000_000, 2_000_000, species="canine", data_dir=str(tmp_path)
            )

        assert result.status is RecombStatus.OK
        assert result.detail == ""
        assert list(result.frame["pos"]) == [1500000]

    @patch("pylocuszoom.recombination.get_default_data_dir")
    def test_no_outcome_warns(self, mock_get_dir, tmp_path, recwarn):
        """The caller renders one policy; this layer only reports."""
        mock_get_dir.return_value = tmp_path / "recomb_data"

        with self._failing_download(DataDownloadError("Network error")):
            recomb_for_region(1, 1_000_000, 2_000_000, species="canine")

        assert list(recwarn) == []


class TestDownloadLiftoverChainDownloadError:
    """download_liftover_chain surfaces failures as DataDownloadError."""

    @patch("pylocuszoom.recombination.download_file")
    def test_download_error_propagates(self, mock_download, tmp_path, monkeypatch):
        monkeypatch.setattr(
            "pylocuszoom.recombination.get_default_data_dir", lambda: tmp_path
        )
        mock_download.side_effect = DataDownloadError("404 Client Error")

        with pytest.raises(DataDownloadError, match="404 Client Error"):
            download_liftover_chain(force=True)


@pytest.mark.parametrize("species", [None, "canine", "human"])
def test_custom_maps_are_read_only_and_need_no_complete_bundle(
    tmp_path, monkeypatch, species
):
    (tmp_path / "chr1_recomb.tsv").write_text("chr\tpos\trate\tcM\n1\t150\t42\t0.1\n")
    (tmp_path / "notes.txt").write_text("caller data")

    def unexpected_download(*args, **kwargs):
        raise AssertionError("custom maps must not download")

    monkeypatch.setattr("pylocuszoom.recombination.download_file", unexpected_download)
    before = {path.name: path.read_bytes() for path in tmp_path.iterdir()}
    result = recomb_for_region(1, 100, 200, species=species, data_dir=str(tmp_path))
    assert result.status is RecombStatus.OK
    assert result.frame["rate"].tolist() == [42]
    assert {path.name: path.read_bytes() for path in tmp_path.iterdir()} == before


def test_custom_maps_are_already_in_the_requested_build(tmp_path, monkeypatch):
    (tmp_path / "chr1_recomb.tsv").write_text("chr\tpos\trate\tcM\n1\t150\t42\t0.1\n")

    def unexpected_liftover(*args, **kwargs):
        raise AssertionError("custom map coordinates must not be reinterpreted")

    monkeypatch.setattr(
        "pylocuszoom.recombination.ensure_recomb_maps", lambda **kwargs: tmp_path
    )
    monkeypatch.setattr(
        "pylocuszoom.recombination.liftover_recombination_map", unexpected_liftover
    )
    result = recomb_for_region(
        1, 100, 200, data_dir=str(tmp_path), genome_build="canfam4"
    )
    assert result.frame["pos"].tolist() == [150]


def test_managed_maps_do_not_silently_use_an_unsupported_build(tmp_path, monkeypatch):
    write_canine_map_set(tmp_path, "chr\tpos\trate\tcM\n1\t150\t1\t0.1\n")
    monkeypatch.setattr(
        "pylocuszoom.recombination.get_default_data_dir", lambda: tmp_path
    )
    result = recomb_for_region(1, 100, 200, genome_build="unknown-build")
    assert result.status.value == "build_unavailable"
    assert result.frame is None
    assert "unknown-build" in result.detail


def test_managed_maps_use_their_known_liftover_chain(tmp_path, monkeypatch):
    write_canine_map_set(tmp_path, "chr\tpos\trate\tcM\n1\t150\t1\t0.1\n")
    monkeypatch.setattr(
        "pylocuszoom.recombination.get_default_data_dir", lambda: tmp_path
    )

    def lift(frame, from_build, to_build, chrom):
        assert from_build == "canfam3"
        assert to_build == "canfam4"
        return frame.assign(pos=frame["pos"] + 100)

    monkeypatch.setattr("pylocuszoom.recombination.liftover_recombination_map", lift)
    result = recomb_for_region(1, 200, 300, genome_build="canfam4")
    assert result.status is RecombStatus.OK
    assert result.frame["pos"].tolist() == [250]


class TestLiftoverChainFailures:
    """A chain that cannot be fetched or read skips the overlay like a map does."""

    CHAIN = "chain 1000 chr1 5000 + 0 1000 chr1 5000 + 100 1100 1\n1000\n\n"

    @pytest.fixture
    def chain_dir(self, tmp_path, monkeypatch):
        """Install managed canfam3 maps and return the (not yet created) chain dir."""
        maps = tmp_path / "recombination_maps"
        write_canine_map_set(maps, "chr\tpos\trate\tcM\n1\t150\t1.0\t0.1\n")
        chain_dir = tmp_path / "liftover"
        monkeypatch.setattr(
            "pylocuszoom.recombination.get_default_data_dir", lambda: maps
        )
        monkeypatch.setattr(
            "pylocuszoom.recombination.get_chain_dir", lambda: chain_dir
        )
        return chain_dir

    def test_a_chain_download_failure_is_a_download_failure(
        self, chain_dir, monkeypatch
    ):
        monkeypatch.setattr(
            "pylocuszoom.recombination.download_file",
            Mock(side_effect=DataDownloadError("simulated chain 404")),
        )

        result = recomb_for_region(1, 1, 5000, species="canine", genome_build="canfam4")

        assert result.status is RecombStatus.DOWNLOAD_FAILED
        assert "simulated chain 404" in result.detail

    def test_a_corrupt_cached_chain_is_downloaded_again(self, chain_dir, monkeypatch):
        import gzip

        chain_dir.mkdir()
        (chain_dir / "canFam3ToCanFam4.over.chain.gz").write_bytes(b"not a chain")

        def download(url, dest, desc=None):
            Path(dest).write_bytes(gzip.compress(self.CHAIN.encode()))

        fetch = Mock(side_effect=download)
        monkeypatch.setattr("pylocuszoom.recombination.download_file", fetch)

        result = recomb_for_region(1, 1, 5000, species="canine", genome_build="canfam4")

        assert fetch.call_count == 1
        assert result.status is RecombStatus.OK
        assert result.frame["pos"].tolist() == [250]

    def test_a_chain_that_stays_unreadable_is_a_download_failure(
        self, chain_dir, monkeypatch
    ):
        def download(url, dest, desc=None):
            Path(dest).write_bytes(b"<html>502</html>")

        monkeypatch.setattr("pylocuszoom.recombination.download_file", download)

        result = recomb_for_region(1, 1, 5000, species="canine", genome_build="canfam4")

        assert result.status is RecombStatus.DOWNLOAD_FAILED
        assert "unreadable" in result.detail


def test_a_lift_that_drops_the_whole_map_is_not_ok(tmp_path, recwarn):
    """An empty overlay after liftover must say why, not report OK."""
    (tmp_path / "chr1_recomb.tsv").write_text(
        "chr\tpos\trate\tcM\n1\t150\t42\t0.1\n1\t180\t40\t0.2\n"
    )

    result = recomb_for_region(
        1,
        100,
        200,
        species="canine",
        data_dir=str(tmp_path),
        lifter=InMemoryLifter({("chr2", 0): 0}),
    )

    assert result.status is not RecombStatus.OK
    assert "chr1" in result.detail
    assert list(recwarn) == []
