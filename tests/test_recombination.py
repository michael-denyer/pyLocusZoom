"""Tests for recombination rate overlay module."""

import io
import tarfile
from pathlib import Path
from unittest.mock import MagicMock, Mock, patch

import pandas as pd
import pytest

from pylocuszoom._liftover import InMemoryLifter
from pylocuszoom.colors import RECOMB_COLOR
from pylocuszoom.exceptions import DataDownloadError, ValidationError
from pylocuszoom.recombination import (
    CANINE_SOURCE,
    RecombStatus,
    _publish_map_generation,
    _stage_archive,
    download_canine_recombination_maps,
    download_liftover_chain,
    ensure_recomb_header,
    ensure_recomb_maps,
    get_default_data_dir,
    get_recombination_rate_for_region,
    liftover_recombination_map,
    load_recombination_map,
    recomb_for_region,
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

    @patch("pylocuszoom.recombination.download_file")
    def test_rejects_wrong_39_file_manifest(self, mock_download, tmp_path, monkeypatch):
        """A count of 39 files is not proof that the canine set is complete."""
        monkeypatch.setattr(
            "pylocuszoom.recombination.get_default_data_dir", lambda: tmp_path
        )
        for i in range(1, 40):
            (tmp_path / f"chr{i}_recomb.tsv").touch()

        mock_download.side_effect = DataDownloadError("download attempted")

        with pytest.raises(DataDownloadError, match="download attempted"):
            download_canine_recombination_maps(force=False)

    @staticmethod
    def _fake_archive(url, dest, desc=None):
        TestStageArchive._tar(
            dest,
            {
                f"chr{i}.txt": f"chr\tpos\trate\tcM\n{i}\t1\t1\t0\n"
                for i in range(1, 39)
            },
        )

    def test_caller_files_in_output_dir_survive(self, tmp_path, monkeypatch):
        """output_dir is the caller's; the library never moves or deletes it."""
        monkeypatch.setattr(
            "pylocuszoom.recombination.download_file", self._fake_archive
        )
        caller_dir = tmp_path / "my_project_data"
        caller_dir.mkdir()
        (caller_dir / "genotypes.bed").write_text("precious")

        with pytest.raises(ValidationError, match="my_project_data"):
            download_canine_recombination_maps(output_dir=str(caller_dir))

        assert (caller_dir / "genotypes.bed").read_text() == "precious"
        assert sorted(p.name for p in tmp_path.iterdir()) == ["my_project_data"]

    def test_a_custom_map_beside_a_complete_set_survives(self, tmp_path, monkeypatch):
        """An extra caller map makes the set inexact; it must not be wiped."""
        monkeypatch.setattr(
            "pylocuszoom.recombination.download_file", self._fake_archive
        )
        for i in range(1, 39):
            (tmp_path / f"chr{i}_recomb.tsv").write_text("old")
        (tmp_path / "chrX_recomb.tsv").write_text("custom")

        with pytest.raises(ValidationError):
            download_canine_recombination_maps(output_dir=str(tmp_path))

        assert (tmp_path / "chrX_recomb.tsv").read_text() == "custom"

    def test_new_output_dir_receives_the_maps(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            "pylocuszoom.recombination.download_file", self._fake_archive
        )
        output = tmp_path / "maps"

        download_canine_recombination_maps(output_dir=str(output))

        assert {p.name for p in output.iterdir()} == CANINE_SOURCE.filenames

    def test_force_refreshes_an_output_dir_holding_only_maps(
        self, tmp_path, monkeypatch
    ):
        monkeypatch.setattr(
            "pylocuszoom.recombination.download_file", self._fake_archive
        )
        for i in range(1, 39):
            (tmp_path / f"chr{i}_recomb.tsv").write_text("old")

        download_canine_recombination_maps(output_dir=str(tmp_path), force=True)

        assert (tmp_path / "chr1_recomb.tsv").read_text().startswith("chr\tpos")


class TestStageArchive:
    """Archive members become canonical map files, never extracted paths."""

    @staticmethod
    def _tar(path: Path, members: dict[str, str]) -> None:
        with tarfile.open(path, "w:gz") as tar:
            for name, content in members.items():
                data = content.encode()
                info = tarfile.TarInfo(name=name)
                info.size = len(data)
                tar.addfile(info, io.BytesIO(data))

    def test_stages_only_maps_with_canonical_names_and_headers(self, tmp_path):
        archive = tmp_path / "maps.tar.gz"
        self._tar(
            archive,
            {
                "nested/chr1.txt": "1\t1000\t0.5\t0.1\n",
                "chr7_average_canFam3.1.txt": "chr\tpos\trate\tcM\n7\t1000\t0.5\t0.1\n",
                "README.txt": "ignored",
            },
        )
        staging = tmp_path / "out"
        staging.mkdir()
        _stage_archive(archive, CANINE_SOURCE, staging)
        assert {p.name for p in staging.iterdir()} == {
            "chr1_recomb.tsv",
            "chr7_recomb.tsv",
        }
        assert (
            staging / "chr1_recomb.tsv"
        ).read_text() == "chr\tpos\trate\tcM\n1\t1000\t0.5\t0.1\n"
        assert (
            staging / "chr7_recomb.tsv"
        ).read_text() == "chr\tpos\trate\tcM\n7\t1000\t0.5\t0.1\n"

    def test_rejects_traversal_without_writing_outside_staging(self, tmp_path):
        archive = tmp_path / "maps.tar.gz"
        self._tar(archive, {"../chr1.txt": "1\t1\t1\t1\n"})
        staging = tmp_path / "out"
        staging.mkdir()
        with pytest.raises(DataDownloadError, match="unsafe member"):
            _stage_archive(archive, CANINE_SOURCE, staging)
        assert not (tmp_path / "chr1.txt").exists()
        assert list(staging.iterdir()) == []

    def test_a_file_that_is_not_a_tarball_names_its_source(self, tmp_path):
        archive = tmp_path / "maps.tar.gz"
        archive.write_bytes(b"not a gzip archive")
        with pytest.raises(DataDownloadError, match="dog_genetic_maps"):
            _stage_archive(archive, CANINE_SOURCE, tmp_path)

    def test_a_filename_without_a_chromosome_is_an_error(self, tmp_path):
        archive = tmp_path / "maps.tar.gz"
        self._tar(archive, {"chr_notes.txt": "chr\tpos\trate\tcM\n1\t1\t1\t1\n"})
        staging = tmp_path / "out"
        staging.mkdir()
        with pytest.raises(DataDownloadError, match="does not name a chromosome"):
            _stage_archive(archive, CANINE_SOURCE, staging)
        assert list(staging.iterdir()) == []


class TestPublishMapGeneration:
    """Atomic publication tests for the recombination map set."""

    @staticmethod
    def _write_maps(path: Path, content: str, *, complete: bool = True) -> None:
        path.mkdir(parents=True, exist_ok=True)
        stop = 39 if complete else 38
        for chrom in range(1, stop):
            (path / f"chr{chrom}_recomb.tsv").write_text(content)

    def test_switches_complete_generations(self, tmp_path):
        output = tmp_path / "maps"
        self._write_maps(output, "old")

        first_staging = tmp_path / "first-staging"
        self._write_maps(first_staging, "first")
        _publish_map_generation(first_staging, output, CANINE_SOURCE)

        assert not output.is_symlink(), "the published set is a plain directory"
        assert (output / "chr1_recomb.tsv").read_text() == "first"

        second_staging = tmp_path / "second-staging"
        self._write_maps(second_staging, "second")
        _publish_map_generation(second_staging, output, CANINE_SOURCE)

        assert (output / "chr1_recomb.tsv").read_text() == "second"
        assert list(tmp_path.glob(".maps.previous-*")) == [], (
            "the replaced set must not be left behind"
        )

    def test_replaces_a_legacy_symlinked_generation(self, tmp_path):
        """Releases before this one published the maps behind a symlink."""
        generation = tmp_path / ".maps.generation-old"
        self._write_maps(generation, "old")
        output = tmp_path / "maps"
        output.symlink_to(generation.name, target_is_directory=True)

        staging = tmp_path / "staging"
        self._write_maps(staging, "new")
        _publish_map_generation(staging, output, CANINE_SOURCE)

        assert not output.is_symlink()
        assert (output / "chr1_recomb.tsv").read_text() == "new"
        assert (generation / "chr1_recomb.tsv").read_text() == "old", (
            "a replaced symlink does not confer ownership of its target"
        )

    def test_incomplete_generation_leaves_active_maps_unchanged(self, tmp_path):
        output = tmp_path / "maps"
        self._write_maps(output, "old")
        staging = tmp_path / "incomplete-staging"
        self._write_maps(staging, "new", complete=False)

        with pytest.raises(DataDownloadError, match="complete canine map set"):
            _publish_map_generation(staging, output, CANINE_SOURCE)

        assert (output / "chr1_recomb.tsv").read_text() == "old"
        assert (output / "chr38_recomb.tsv").read_text() == "old"


class TestEnsureRecombMaps:
    """Tests for ensure_recomb_maps function."""

    @staticmethod
    def _patched_download(mock_download):
        """Stand in for the network-facing download step."""
        return patch(
            "pylocuszoom.recombination.download_recombination_maps", mock_download
        )

    @patch("pylocuszoom.recombination.get_default_data_dir")
    def test_ensure_recomb_maps_downloads_if_missing(self, mock_get_dir, tmp_path):
        """Test that ensure_recomb_maps triggers download when maps missing."""
        mock_get_dir.return_value = tmp_path / "recomb_data"
        mock_download = Mock(return_value=tmp_path / "recomb_data")

        with self._patched_download(mock_download):
            result = ensure_recomb_maps(species="canine")

        mock_download.assert_called_once()
        assert result == tmp_path / "recomb_data"

    @patch("pylocuszoom.recombination.get_default_data_dir")
    def test_ensure_recomb_maps_skips_download_if_exists(self, mock_get_dir, tmp_path):
        """Test that ensure_recomb_maps skips download when maps exist."""
        data_dir = tmp_path / "recomb_data"
        data_dir.mkdir()
        # Create the complete 38-autosome map set.
        for i in range(1, 39):
            (data_dir / f"chr{i}_recomb.tsv").touch()

        mock_get_dir.return_value = data_dir
        mock_download = Mock()

        with self._patched_download(mock_download):
            result = ensure_recomb_maps(species="canine")

        mock_download.assert_not_called()
        assert result == data_dir

    def test_ensure_recomb_maps_non_canine_returns_none(self):
        """Test that ensure_recomb_maps returns None for non-canine species."""
        result = ensure_recomb_maps(species="human")
        assert result is None

    @patch("pylocuszoom.recombination.get_default_data_dir")
    def test_ensure_recomb_maps_propagates_a_download_error(
        self, mock_get_dir, tmp_path
    ):
        """The data layer raises; recomb_for_region is what degrades."""
        mock_get_dir.return_value = tmp_path / "recomb_data"
        mock_download = Mock(side_effect=DataDownloadError("Network error"))

        with self._patched_download(mock_download):
            with pytest.raises(DataDownloadError, match="Network error"):
                ensure_recomb_maps(species="canine")

    @patch("pylocuszoom.recombination.get_default_data_dir")
    def test_ensure_recomb_maps_propagates_an_io_error(self, mock_get_dir, tmp_path):
        mock_get_dir.return_value = tmp_path / "recomb_data"
        mock_download = Mock(side_effect=OSError("Disk full"))

        with self._patched_download(mock_download):
            with pytest.raises(OSError, match="Disk full"):
                ensure_recomb_maps(species="canine")


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


class TestEnsureRecombMapsCorruptArchive:
    """A corrupt archive surfaces as a status, not a crash, through the plotter."""

    @patch("pylocuszoom.recombination.download_file")
    def test_a_corrupt_archive_is_a_download_failure(
        self, mock_download, tmp_path, monkeypatch
    ):
        monkeypatch.setattr(
            "pylocuszoom.recombination.get_default_data_dir", lambda: tmp_path / "x"
        )

        def write_garbage(url, dest_path, desc=None):
            dest_path.write_bytes(b"not a gzip archive")

        mock_download.side_effect = write_garbage

        result = recomb_for_region(1, 1_000_000, 2_000_000, species="canine")

        assert result.status is RecombStatus.DOWNLOAD_FAILED


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


class TestArchiveWithoutMaps:
    def test_non_map_members_are_ignored(self, tmp_path, monkeypatch):
        def download(url, dest, desc):
            TestStageArchive._tar(dest, {"README.txt": "notes"})

        monkeypatch.setattr("pylocuszoom.recombination.download_file", download)
        with pytest.raises(DataDownloadError, match="Could not find chromosome"):
            download_canine_recombination_maps(tmp_path / "output")
        assert not (tmp_path / "output").exists()


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
        """Return a download_file mock that writes our fake tarball."""

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

    @patch("pylocuszoom.recombination.download_file")
    def test_raises_on_html_corrupted_body(self, mock_download, tmp_path):
        """An HTML error body masquerading as a map is a download failure."""
        html_content = "<html><body>502 Bad Gateway</body></html>\n"
        mock_download.side_effect = self._fake_download("chr1.txt", html_content)

        with pytest.raises(DataDownloadError, match="Unrecognised first token"):
            download_canine_recombination_maps(tmp_path / "out")

    @patch("pylocuszoom.recombination.download_file")
    def test_plot_warns_and_renders_without_overlay_on_corrupt_archive(
        self, mock_download, tmp_path, monkeypatch
    ):
        """A corrupted mirror must be loud, but must not stop the plot."""
        from pylocuszoom import LocusZoomPlotter

        html_content = "<html><body>502 Bad Gateway</body></html>\n"
        mock_download.side_effect = self._fake_download("chr1.txt", html_content)
        monkeypatch.setattr(
            "pylocuszoom.recombination.get_default_data_dir", lambda: tmp_path / "out"
        )
        plotter = LocusZoomPlotter(species="canine", log_level=None)
        gwas = pd.DataFrame(
            {"pos": [1_000_000, 1_050_000], "p_value": [0.5, 1e-6], "rs": ["a", "b"]}
        )

        with pytest.warns(UserWarning, match="Unrecognised first token"):
            fig = plotter.plot(gwas, chrom=1, start=900_000, end=1_200_000)

        assert fig is not None
        assert not (tmp_path / "out").exists()

    @patch("pylocuszoom.recombination.download_file")
    def test_accepts_lowercase_chr_header(self, mock_download, tmp_path):
        """Pre-existing canonical header form must still be accepted."""
        content = "chr\tpos\trate\tcM\n1\t1000\t0.5\t0.1\n"
        mock_download.side_effect = self._fake_download("chr1.txt", content)

        result = download_canine_recombination_maps(tmp_path / "out")
        assert (result / "chr1_recomb.tsv").exists()

    @patch("pylocuszoom.recombination.download_file")
    def test_accepts_capitalised_chromosome_header(self, mock_download, tmp_path):
        """'Chromosome' is used by several mirrors; must pass."""
        content = "Chromosome\tPosition\tRate\tcM\n1\t1000\t0.5\t0.1\n"
        mock_download.side_effect = self._fake_download("chr1.txt", content)

        result = download_canine_recombination_maps(tmp_path / "out")
        assert (result / "chr1_recomb.tsv").exists()

    @patch("pylocuszoom.recombination.download_file")
    def test_accepts_hash_prefixed_header(self, mock_download, tmp_path):
        """Some maps use '#chr' as a commented header; must pass."""
        content = "#chr\tpos\trate\tcM\n1\t1000\t0.5\t0.1\n"
        mock_download.side_effect = self._fake_download("chr1.txt", content)

        result = download_canine_recombination_maps(tmp_path / "out")
        assert (result / "chr1_recomb.tsv").exists()

    @patch("pylocuszoom.recombination.download_file")
    def test_accepts_numeric_first_token_prepends_header(self, mock_download, tmp_path):
        """Numeric first token means no header; one is prepended."""
        content = "1\t1000\t0.5\t0.1\n1\t2000\t0.6\t0.2\n"
        mock_download.side_effect = self._fake_download("chr1.txt", content)

        result = download_canine_recombination_maps(tmp_path / "out")
        out_file = result / "chr1_recomb.tsv"
        assert out_file.exists()
        assert out_file.read_text().startswith("chr\tpos\trate\tcM\n")


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


@pytest.mark.parametrize(
    "member_type", [tarfile.SYMTYPE, tarfile.LNKTYPE, tarfile.FIFOTYPE]
)
def test_archive_rejects_links_and_special_files_before_publication(
    tmp_path, monkeypatch, member_type
):
    output = tmp_path / "maps"
    escaped = tmp_path / "escaped"
    escaped.mkdir()

    def download(url, dest, desc):
        with tarfile.open(dest, "w:gz") as archive:
            link = tarfile.TarInfo("bridge")
            link.type = member_type
            link.linkname = str(escaped)
            archive.addfile(link)
            for name, body in [("bridge/witness.txt", b"escaped")] + [
                (f"chr{chrom}.txt", f"{chrom}\t150\t1\t0.1\n".encode())
                for chrom in range(1, 39)
            ]:
                member = tarfile.TarInfo(name)
                member.size = len(body)
                archive.addfile(member, io.BytesIO(body))

    monkeypatch.setattr("pylocuszoom.recombination.download_file", download)
    with pytest.raises(DataDownloadError, match="regular file|unsafe member"):
        download_canine_recombination_maps(output)
    assert list(escaped.iterdir()) == []
    assert not output.exists()


def test_archive_rejects_duplicate_chromosome_maps(tmp_path, monkeypatch):
    def download(url, dest, desc):
        with tarfile.open(dest, "w:gz") as archive:
            for chrom in list(range(1, 39)) + [1]:
                member = tarfile.TarInfo(f"chr{chrom}.txt")
                data = f"{chrom}\t150\t1\t0.1\n".encode()
                member.size = len(data)
                archive.addfile(member, io.BytesIO(data))

    monkeypatch.setattr("pylocuszoom.recombination.download_file", download)
    with pytest.raises(DataDownloadError, match="Duplicate"):
        download_canine_recombination_maps(tmp_path / "maps")
    assert not (tmp_path / "maps").exists()


def test_managed_maps_do_not_silently_use_an_unsupported_build(tmp_path, monkeypatch):
    TestPublishMapGeneration._write_maps(
        tmp_path, "chr\tpos\trate\tcM\n1\t150\t1\t0.1\n"
    )
    monkeypatch.setattr(
        "pylocuszoom.recombination.get_default_data_dir", lambda: tmp_path
    )
    result = recomb_for_region(1, 100, 200, genome_build="unknown-build")
    assert result.status.value == "build_unavailable"
    assert result.frame is None
    assert "unknown-build" in result.detail


def test_managed_maps_use_their_known_liftover_chain(tmp_path, monkeypatch):
    TestPublishMapGeneration._write_maps(
        tmp_path, "chr\tpos\trate\tcM\n1\t150\t1\t0.1\n"
    )
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
        TestPublishMapGeneration._write_maps(
            maps, "chr\tpos\trate\tcM\n1\t150\t1.0\t0.1\n"
        )
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
