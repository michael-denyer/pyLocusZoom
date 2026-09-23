"""Recombination maps: loading, region lookup, liftover and why a lookup fails."""

from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from pylocuszoom._liftover import InMemoryLifter
from pylocuszoom.colors import RECOMB_COLOR
from pylocuszoom.exceptions import (
    DataDownloadError,
    PyLocusZoomError,
    RecombinationMapNotFound,
    ValidationError,
)
from pylocuszoom.recombination import (
    get_default_data_dir,
    get_recombination_rate_for_region,
    load_recombination_map,
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

    def test_non_numeric_values_produce_warning(self, tmp_path, warning_records):
        """Non-numeric values in pos/rate should produce a warning."""
        map_content = "chr\tpos\trate\tcM\n1\t1000\t0.5\t0.001\n1\tBAD\t1.2\t0.005\n"
        (tmp_path / "chr1_recomb.tsv").write_text(map_content)

        result = load_recombination_map(chrom=1, data_dir=str(tmp_path))

        # Only valid row should remain
        assert len(result) == 1
        assert result["pos"].iloc[0] == 1000
        assert any(
            "non-numeric values" in record and "chr1" in record
            for record in warning_records
        )


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


class TestRegionLookupFailures:
    """Every reason there is no frame is a typed error, and none of them warns."""

    @staticmethod
    def _failing_download(exc):
        return patch(
            "pylocuszoom.recombination.download_file",
            Mock(side_effect=exc),
        )

    def test_a_download_failure_carries_its_cause(self, cache_home):
        with (
            self._failing_download(DataDownloadError("Network error")),
            pytest.raises(DataDownloadError, match="Network error"),
        ):
            get_recombination_rate_for_region(1, 1_000_000, 2_000_000, species="canine")

    def test_an_unwritable_cache_is_a_download_error_not_an_os_error(
        self, tmp_path, monkeypatch
    ):
        blocker = tmp_path / "not_a_directory"
        blocker.write_text("")
        monkeypatch.setenv("XDG_CACHE_HOME", str(blocker / "cache"))

        with (
            self._failing_download(AssertionError("must not download")),
            pytest.raises(DataDownloadError, match="Could not write"),
        ):
            get_recombination_rate_for_region(1, 1_000_000, 2_000_000, species="canine")

    def test_a_species_with_no_built_in_maps_says_so(self):
        with pytest.raises(RecombinationMapNotFound, match="'human'"):
            get_recombination_rate_for_region(1, 1_000_000, 2_000_000, species="human")

    def test_a_chromosome_the_map_set_does_not_cover_is_its_own_error(self, tmp_path):
        with pytest.raises(RecombinationMapNotFound, match="chr41"):
            get_recombination_rate_for_region(
                41, 1_000_000, 2_000_000, species="canine", data_dir=str(tmp_path)
            )

    def test_a_missing_map_is_still_a_file_not_found_error(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_recombination_map(1, data_dir=str(tmp_path))

    def test_no_failure_warns(self, cache_home, recwarn):
        """The caller renders one skip policy; this layer only raises."""
        with (
            self._failing_download(DataDownloadError("Network error")),
            pytest.raises(PyLocusZoomError),
        ):
            get_recombination_rate_for_region(1, 1_000_000, 2_000_000, species="canine")

        assert list(recwarn) == []


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
    frame = get_recombination_rate_for_region(
        1, 100, 200, species=species, data_dir=str(tmp_path)
    )
    assert frame["rate"].tolist() == [42]
    assert {path.name: path.read_bytes() for path in tmp_path.iterdir()} == before


def test_custom_maps_are_already_in_the_requested_build(tmp_path, monkeypatch):
    (tmp_path / "chr1_recomb.tsv").write_text("chr\tpos\trate\tcM\n1\t150\t42\t0.1\n")

    def unexpected_liftover(*args, **kwargs):
        raise AssertionError("custom map coordinates must not be reinterpreted")

    monkeypatch.setattr("pylocuszoom.recombination.chain_lifter", unexpected_liftover)
    frame = get_recombination_rate_for_region(
        1, 100, 200, data_dir=str(tmp_path), genome_build="canfam4"
    )
    assert frame["pos"].tolist() == [150]


def test_managed_maps_do_not_silently_use_an_unsupported_build(cache_home):
    write_canine_map_set(
        cache_home / "recombination_maps", "chr\tpos\trate\tcM\n1\t150\t1\t0.1\n"
    )
    with pytest.raises(ValidationError, match="unknown-build"):
        get_recombination_rate_for_region(1, 100, 200, genome_build="unknown-build")


def test_managed_maps_use_their_known_liftover_chain(cache_home, monkeypatch):
    write_canine_map_set(
        cache_home / "recombination_maps", "chr\tpos\trate\tcM\n1\t150\t1\t0.1\n"
    )

    def registered_chain(source, target):
        assert (source.key, target.key) == ("canfam3", "canfam4")
        return InMemoryLifter({("chr1", 149): 249})

    monkeypatch.setattr("pylocuszoom.recombination.chain_lifter", registered_chain)
    frame = get_recombination_rate_for_region(1, 200, 300, genome_build="canfam4")
    assert frame["pos"].tolist() == [250]


class TestLiftoverChainFailures:
    """A chain that cannot be fetched or read fails the lookup like a map does."""

    CHAIN = "chain 1000 chr1 5000 + 0 1000 chr1 5000 + 100 1100 1\n1000\n\n"

    @pytest.fixture
    def chain_dir(self, cache_home):
        """Install managed canfam3 maps and return the (not yet created) chain dir."""
        write_canine_map_set(
            cache_home / "recombination_maps",
            "chr\tpos\trate\tcM\n1\t150\t1.0\t0.1\n",
        )
        return cache_home / "liftover"

    def test_a_chain_download_failure_is_a_download_error(self, chain_dir, monkeypatch):
        monkeypatch.setattr(
            "pylocuszoom._liftover.download_file",
            Mock(side_effect=DataDownloadError("simulated chain 404")),
        )

        with pytest.raises(DataDownloadError, match="simulated chain 404"):
            get_recombination_rate_for_region(
                1, 1, 5000, species="canine", genome_build="canfam4"
            )

    def test_a_corrupt_cached_chain_is_downloaded_again(self, chain_dir, monkeypatch):
        import gzip

        chain_dir.mkdir()
        (chain_dir / "canFam3ToCanFam4.over.chain.gz").write_bytes(b"not a chain")

        def download(url, dest, desc=None):
            Path(dest).write_bytes(gzip.compress(self.CHAIN.encode()))

        fetch = Mock(side_effect=download)
        monkeypatch.setattr("pylocuszoom._liftover.download_file", fetch)

        frame = get_recombination_rate_for_region(
            1, 1, 5000, species="canine", genome_build="canfam4"
        )

        assert fetch.call_count == 1
        assert frame["pos"].tolist() == [250]

    def test_a_chain_that_stays_unreadable_is_a_download_error(
        self, chain_dir, monkeypatch
    ):
        def download(url, dest, desc=None):
            Path(dest).write_bytes(b"<html>502</html>")

        monkeypatch.setattr("pylocuszoom._liftover.download_file", download)

        with pytest.raises(DataDownloadError, match="unreadable"):
            get_recombination_rate_for_region(
                1, 1, 5000, species="canine", genome_build="canfam4"
            )


def test_a_lift_that_drops_the_whole_map_says_why(tmp_path, recwarn):
    """An empty overlay after liftover must say why, not return an empty frame."""
    (tmp_path / "chr1_recomb.tsv").write_text(
        "chr\tpos\trate\tcM\n1\t150\t42\t0.1\n1\t180\t40\t0.2\n"
    )

    with pytest.raises(ValidationError, match="chr1 is unknown to the liftover chain"):
        get_recombination_rate_for_region(
            1,
            100,
            200,
            species="canine",
            data_dir=str(tmp_path),
            lifter=InMemoryLifter({("chr2", 0): 0}),
        )

    assert list(recwarn) == []


def test_a_lifted_map_is_sorted_by_its_new_positions(tmp_path):
    """The overlay line is drawn left to right, whatever order the chain gives."""
    (tmp_path / "chr1_recomb.tsv").write_text(
        "chr\tpos\trate\tcM\n1\t1000\t0.5\t0.1\n1\t2000\t0.6\t0.2\n"
    )
    lifter = InMemoryLifter({("chr1", 999): 4_999, ("chr1", 1_999): 999})

    frame = get_recombination_rate_for_region(
        1, 1, 10_000, species="canine", data_dir=str(tmp_path), lifter=lifter
    )

    assert frame["pos"].tolist() == [1_000, 5_000]
    assert frame["rate"].tolist() == [0.6, 0.5]
