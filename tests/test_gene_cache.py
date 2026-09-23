"""Atomic publication of complete gene and exon cache entries."""

import pandas as pd
import pytest

from pylocuszoom._gene_cache import load_annotations, save_annotations
from pylocuszoom._gene_source import GeneAnnotations


def test_failed_exon_write_keeps_previous_complete_entry(tmp_path, monkeypatch):
    old = GeneAnnotations(
        pd.DataFrame({"name": ["old"]}), pd.DataFrame({"name": ["old"]})
    )
    new = GeneAnnotations(
        pd.DataFrame({"name": ["new"]}), pd.DataFrame({"name": ["new"]})
    )
    save_annotations(old, tmp_path, "human", "1", 100, 200)
    original = pd.DataFrame.to_csv

    def fail_new_exons(frame, *args, **kwargs):
        if frame is new.exons:
            raise OSError("interrupted exon write")
        return original(frame, *args, **kwargs)

    monkeypatch.setattr(pd.DataFrame, "to_csv", fail_new_exons)
    save_annotations(new, tmp_path, "human", "1", 100, 200)
    loaded = load_annotations(tmp_path, "human", "1", 100, 200)
    assert loaded is not None
    assert loaded.genes["name"].tolist() == ["old"]
    assert loaded.exons["name"].tolist() == ["old"]
    assert not list(tmp_path.rglob("*.part"))


def test_readers_see_complete_entry_while_another_writer_is_staging(
    tmp_path, monkeypatch
):
    from concurrent.futures import ThreadPoolExecutor
    from threading import Event

    old = GeneAnnotations(
        pd.DataFrame({"name": ["old"]}), pd.DataFrame({"name": ["old"]})
    )
    new = GeneAnnotations(
        pd.DataFrame({"name": ["new"]}), pd.DataFrame({"name": ["new"]})
    )
    save_annotations(old, tmp_path, "human", "1", 100, 200)
    staging, resume = Event(), Event()
    original = pd.DataFrame.to_csv

    def pause_before_exons(frame, *args, **kwargs):
        if frame is new.exons:
            staging.set()
            assert resume.wait(5)
        return original(frame, *args, **kwargs)

    monkeypatch.setattr(pd.DataFrame, "to_csv", pause_before_exons)
    with ThreadPoolExecutor(max_workers=1) as executor:
        writer = executor.submit(
            save_annotations, new, tmp_path, "human", "1", 100, 200
        )
        try:
            assert staging.wait(5)
            during_write = load_annotations(tmp_path, "human", "1", 100, 200)
        finally:
            resume.set()
        writer.result(timeout=10)
    assert during_write.genes["name"].tolist() == ["old"]
    assert during_write.exons["name"].tolist() == ["old"]
    published = load_annotations(tmp_path, "human", "1", 100, 200)
    assert published.genes["name"].tolist() == ["new"]
    assert published.exons["name"].tolist() == ["new"]


def test_corrupt_or_incomplete_archives_are_cache_misses(tmp_path):
    from zipfile import ZipFile

    entry = GeneAnnotations(
        pd.DataFrame({"name": ["gene"]}), pd.DataFrame({"name": ["exon"]})
    )
    save_annotations(entry, tmp_path, "human", "1", 100, 200)
    archive_path = next(tmp_path.rglob("*.zip"))
    archive_path.write_bytes(b"truncated archive")
    assert load_annotations(tmp_path, "human", "1", 100, 200) is None
    with ZipFile(archive_path, "w") as archive:
        archive.writestr("genes.csv", "name\ngene\n")
    assert load_annotations(tmp_path, "human", "1", 100, 200) is None


def test_legacy_pairs_are_not_reused_but_can_be_cleared(tmp_path):
    from pylocuszoom._gene_cache import cache_key, clear_cache

    species = tmp_path / "human"
    species.mkdir()
    key = cache_key("human", "1", 100, 200)
    (species / f"genes_{key}.csv").write_text("name\nnew\n")
    (species / f"exons_{key}.csv").write_text("name\nold\n")
    assert load_annotations(tmp_path, "human", "1", 100, 200) is None
    assert clear_cache(tmp_path) == 2
    assert list(species.iterdir()) == []


def test_encrypted_cache_member_is_a_miss(tmp_path):
    entry = GeneAnnotations(
        pd.DataFrame({"name": ["gene"]}), pd.DataFrame({"name": ["exon"]})
    )
    save_annotations(entry, tmp_path, "human", "1", 100, 200)
    archive_path = next(tmp_path.rglob("*.zip"))
    data = bytearray(archive_path.read_bytes())
    # General-purpose bit 0 marks a ZIP member encrypted. Set it on both
    # copies of the first member's header, as a corrupt cache could contain.
    for signature, offset in ((b"PK\x03\x04", 6), (b"PK\x01\x02", 8)):
        data[data.index(signature) + offset] |= 1
    archive_path.write_bytes(data)
    assert load_annotations(tmp_path, "human", "1", 100, 200) is None


def test_unsupported_cache_compression_is_a_miss(tmp_path):
    entry = GeneAnnotations(
        pd.DataFrame({"name": ["gene"]}), pd.DataFrame({"name": ["exon"]})
    )
    save_annotations(entry, tmp_path, "human", "1", 100, 200)
    archive_path = next(tmp_path.rglob("*.zip"))
    data = bytearray(archive_path.read_bytes())
    for signature, offset in ((b"PK\x03\x04", 8), (b"PK\x01\x02", 10)):
        index = data.index(signature) + offset
        data[index : index + 2] = (99).to_bytes(2, "little")
    archive_path.write_bytes(data)
    assert load_annotations(tmp_path, "human", "1", 100, 200) is None


def test_truncated_cache_member_is_a_miss(tmp_path):
    import struct
    from zipfile import ZIP_DEFLATED, ZipFile

    from pylocuszoom._gene_cache import _entry_file

    archive_path = _entry_file(tmp_path, "human", "1", 100, 200, "")
    archive_path.parent.mkdir()
    with ZipFile(archive_path, "w", compression=ZIP_DEFLATED) as archive:
        archive.writestr("genes.csv", "name\ngene\n")
        archive.writestr("exons.csv", "name\nexon\n")
    data = bytearray(archive_path.read_bytes())
    # Preserve the directory while making the compressed stream end before
    # it can deliver the member's declared uncompressed length.
    local = data.index(b"PK\x03\x04")
    central = data.index(b"PK\x01\x02")
    struct.pack_into("<I", data, local + 18, 1)
    struct.pack_into("<I", data, central + 20, 1)
    archive_path.write_bytes(data)
    assert load_annotations(tmp_path, "human", "1", 100, 200) is None


def test_cache_cleanup_failure_does_not_abort_annotation_save(tmp_path, monkeypatch):
    from pathlib import Path

    old = GeneAnnotations(
        pd.DataFrame({"name": ["old"]}), pd.DataFrame({"name": ["old"]})
    )
    new = GeneAnnotations(
        pd.DataFrame({"name": ["new"]}), pd.DataFrame({"name": ["new"]})
    )
    save_annotations(old, tmp_path, "human", "1", 100, 200)
    original_to_csv = pd.DataFrame.to_csv
    original_unlink = Path.unlink

    def fail_new_write(frame, *args, **kwargs):
        if frame is new.genes:
            raise OSError("interrupted cache write")
        return original_to_csv(frame, *args, **kwargs)

    def fail_cleanup(path, *args, **kwargs):
        if path.name.endswith(".part"):
            raise PermissionError("cache directory became read-only")
        return original_unlink(path, *args, **kwargs)

    with monkeypatch.context() as patch:
        patch.setattr(pd.DataFrame, "to_csv", fail_new_write)
        patch.setattr(Path, "unlink", fail_cleanup)
        save_annotations(new, tmp_path, "human", "1", 100, 200)
    loaded = load_annotations(tmp_path, "human", "1", 100, 200)
    assert loaded.genes["name"].tolist() == ["old"]
    assert loaded.exons["name"].tolist() == ["old"]


def test_invalid_deflate_stream_is_a_miss(tmp_path):
    entry = GeneAnnotations(
        pd.DataFrame({"name": ["gene"]}), pd.DataFrame({"name": ["exon"]})
    )
    save_annotations(entry, tmp_path, "human", "1", 100, 200)
    archive_path = next(tmp_path.rglob("*.zip"))
    data = bytearray(archive_path.read_bytes())
    # Stored CSV bytes cannot be decoded as a deflate stream.
    for signature, offset in ((b"PK\x03\x04", 8), (b"PK\x01\x02", 10)):
        index = data.index(signature) + offset
        data[index : index + 2] = (8).to_bytes(2, "little")
    archive_path.write_bytes(data)
    assert load_annotations(tmp_path, "human", "1", 100, 200) is None


@pytest.mark.parametrize("use_cache", [False, True])
def test_an_unwritable_cache_base_still_serves_genes(tmp_path, monkeypatch, use_cache):
    """Building the cache path must not create it; only a write may fail, quietly."""
    from pylocuszoom._gene_source import GeneSource, empty_annotations
    from pylocuszoom.reference_genes import get_genes_for_build

    blocker = tmp_path / "not_a_directory"
    blocker.write_text("")
    monkeypatch.setenv("XDG_CACHE_HOME", str(blocker / "cache"))
    source = GeneSource(
        name="ensembl",
        cache_species="x",
        build_token="",
        fetch=lambda chrom, start, end: empty_annotations(),
    )

    annotations = get_genes_for_build(source, 1, 1, 100, use_cache=use_cache)

    assert annotations.genes.empty
