"""Fetching, unpacking and publishing the managed recombination map set."""

import io
import tarfile
from pathlib import Path
from unittest.mock import Mock, patch

import pandas as pd
import pytest

from pylocuszoom.exceptions import DataDownloadError, ValidationError
from pylocuszoom.recombination import (
    CANINE_SOURCE,
    _publish_map_generation,
    _stage_archive,
    download_canine_recombination_maps,
    ensure_recomb_header,
    ensure_recomb_maps,
    get_recombination_rate_for_region,
)
from tests.conftest import write_canine_map_set


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

    def test_switches_complete_generations(self, tmp_path):
        output = tmp_path / "maps"
        write_canine_map_set(output, "old")

        first_staging = tmp_path / "first-staging"
        write_canine_map_set(first_staging, "first")
        _publish_map_generation(first_staging, output, CANINE_SOURCE)

        assert not output.is_symlink(), "the published set is a plain directory"
        assert (output / "chr1_recomb.tsv").read_text() == "first"

        second_staging = tmp_path / "second-staging"
        write_canine_map_set(second_staging, "second")
        _publish_map_generation(second_staging, output, CANINE_SOURCE)

        assert (output / "chr1_recomb.tsv").read_text() == "second"
        assert list(tmp_path.glob(".maps.previous-*")) == [], (
            "the replaced set must not be left behind"
        )

    def test_replaces_a_legacy_symlinked_generation(self, tmp_path):
        """Releases before this one published the maps behind a symlink."""
        generation = tmp_path / ".maps.generation-old"
        write_canine_map_set(generation, "old")
        output = tmp_path / "maps"
        output.symlink_to(generation.name, target_is_directory=True)

        staging = tmp_path / "staging"
        write_canine_map_set(staging, "new")
        _publish_map_generation(staging, output, CANINE_SOURCE)

        assert not output.is_symlink()
        assert (output / "chr1_recomb.tsv").read_text() == "new"
        assert (generation / "chr1_recomb.tsv").read_text() == "old", (
            "a replaced symlink does not confer ownership of its target"
        )

    def test_a_stray_map_in_the_target_does_not_make_the_set_incomplete(self, tmp_path):
        """Otherwise every later ensure_recomb_maps would download again."""
        output = tmp_path / "maps"
        write_canine_map_set(output, "old")
        (output / "chr39_recomb.tsv").write_text("stray")
        staging = tmp_path / "staging"
        write_canine_map_set(staging, "new")

        _publish_map_generation(staging, output, CANINE_SOURCE)

        assert {p.name for p in output.iterdir()} == CANINE_SOURCE.filenames

    def test_incomplete_generation_leaves_active_maps_unchanged(self, tmp_path):
        output = tmp_path / "maps"
        write_canine_map_set(output, "old")
        staging = tmp_path / "incomplete-staging"
        write_canine_map_set(staging, "new", complete=False)

        with pytest.raises(DataDownloadError, match="complete canine map set"):
            _publish_map_generation(staging, output, CANINE_SOURCE)

        assert (output / "chr1_recomb.tsv").read_text() == "old"
        assert (output / "chr38_recomb.tsv").read_text() == "old"


class TestConcurrentPublication:
    """Two writers publishing at once converge, and readers never see a gap.

    Map directories are shared: pytest-xdist workers share XDG_CACHE_HOME and
    every Databricks notebook shares /dbfs/FileStore/reference_data.
    """

    @pytest.fixture
    def output(self, tmp_path):
        return tmp_path / "recombination_maps"

    @staticmethod
    def _staged(tmp_path, tag):
        staging = tmp_path / f"staging_{tag}"
        write_canine_map_set(staging, f"chr\tpos\trate\tcM\n1\t1\t1\t{tag}\n")
        return staging

    def _interleave(self, monkeypatch, output, name, other_writer):
        """Run other_writer at writer A's first os.<name>, checking every step.

        Records whether the published directory existed at every rename or
        replace either writer made.
        """
        import os

        real = getattr(os, name)
        seen = []
        state = {"fired": False}

        def interleaving(src, dst, *args, **kwargs):
            seen.append(output.is_dir())
            if not state["fired"]:
                state["fired"] = True
                other_writer()
            return real(src, dst, *args, **kwargs)

        monkeypatch.setattr(os, name, interleaving)
        return seen

    def test_a_writer_that_loses_the_first_install_still_succeeds(
        self, tmp_path, output, monkeypatch
    ):
        staged_a, staged_b = self._staged(tmp_path, "A"), self._staged(tmp_path, "B")
        self._interleave(
            monkeypatch,
            output,
            "rename",
            lambda: _publish_map_generation(staged_b, output, CANINE_SOURCE),
        )

        _publish_map_generation(staged_a, output, CANINE_SOURCE)

        assert {p.name for p in output.iterdir()} == CANINE_SOURCE.filenames
        assert sorted(p.name for p in tmp_path.iterdir()) == [
            "recombination_maps",
            "staging_A",
        ], "nothing moved aside, and the loser's staging is its caller's to drop"

    def test_refreshes_racing_on_a_published_set_never_open_a_gap(
        self, tmp_path, output, monkeypatch
    ):
        _publish_map_generation(self._staged(tmp_path, "old"), output, CANINE_SOURCE)
        staged_a, staged_b = self._staged(tmp_path, "A"), self._staged(tmp_path, "B")
        seen = self._interleave(
            monkeypatch,
            output,
            "replace",
            lambda: _publish_map_generation(staged_b, output, CANINE_SOURCE),
        )

        _publish_map_generation(staged_a, output, CANINE_SOURCE)

        assert seen and all(seen), "a reader found the maps directory missing"
        assert {p.name for p in output.iterdir()} == CANINE_SOURCE.filenames
        assert not list(tmp_path.glob(".*previous*"))
        assert {
            (output / name).read_text().split()[-1] for name in CANINE_SOURCE.filenames
        } <= {"A", "B"}


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
        """The data layer raises; the plotter is what degrades."""
        mock_get_dir.return_value = tmp_path / "recomb_data"
        mock_download = Mock(side_effect=DataDownloadError("Network error"))

        with self._patched_download(mock_download):
            with pytest.raises(DataDownloadError, match="Network error"):
                ensure_recomb_maps(species="canine")

    def test_ensure_recomb_maps_reports_an_unwritable_cache_as_a_download_error(
        self, tmp_path, monkeypatch
    ):
        blocker = tmp_path / "not_a_directory"
        blocker.write_text("")
        monkeypatch.setenv("XDG_CACHE_HOME", str(blocker / "cache"))

        with pytest.raises(DataDownloadError, match="Could not write"):
            ensure_recomb_maps(species="canine")


class TestEnsureRecombMapsCorruptArchive:
    """A corrupt archive is a typed download error, not a crash."""

    @patch("pylocuszoom.recombination.download_file")
    def test_a_corrupt_archive_is_a_download_error(
        self, mock_download, tmp_path, monkeypatch
    ):
        monkeypatch.setattr(
            "pylocuszoom.recombination.get_default_data_dir", lambda: tmp_path / "x"
        )

        def write_garbage(url, dest_path, desc=None):
            dest_path.write_bytes(b"not a gzip archive")

        mock_download.side_effect = write_garbage

        with pytest.raises(DataDownloadError, match="not a valid tar.gz"):
            get_recombination_rate_for_region(1, 1_000_000, 2_000_000, species="canine")


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
        plotter = LocusZoomPlotter(species="canine")
        gwas = pd.DataFrame(
            {
                "chr": 1,
                "pos": [1_000_000, 1_050_000],
                "p_value": [0.5, 1e-6],
                "rs": ["a", "b"],
            }
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
