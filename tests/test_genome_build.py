"""Tests for the genome-build table and the names it folds together."""

import pytest

from pylocuszoom.genome_build import (
    GENOME_BUILDS,
    assembly_token,
    resolve_build,
    ucsc_chrom,
)
from pylocuszoom.species import SPECIES


def test_assembly_token_folds_synonyms():
    """Equivalent spellings of one assembly compare equal."""
    assert assembly_token("CanFam4.0") == assembly_token("UU_Cfam_GSD_1.0")
    assert assembly_token("CanFam3.1") == assembly_token("canfam3")
    assert assembly_token("hg38") == assembly_token("GRCh38")
    assert assembly_token("CanFam3.1") != assembly_token("ROS_Cfam_1.0")


def test_assembly_token_keeps_an_unknown_name_comparable():
    assert assembly_token("Sscrofa_11.1") == assembly_token("sscrofa11.1")


@pytest.mark.parametrize(
    "name,key",
    [
        ("canFam3", "canfam3"),
        ("CanFam3.1", "canfam3"),
        ("UU_Cfam_GSD_1.0", "canfam4"),
        ("felCat9", "felcat9"),
        ("hg19", "grch37"),
    ],
)
def test_resolve_build_accepts_every_spelling(name, key):
    assert resolve_build(name) is GENOME_BUILDS[key]


def test_resolve_build_is_idempotent_and_none_for_unknown():
    build = GENOME_BUILDS["canfam4"]

    assert resolve_build(build) is build
    assert resolve_build(None) is None
    assert resolve_build("Sscrofa_11.1") is None


@pytest.mark.parametrize("species", [s for s in SPECIES.values() if s.default_build])
def test_every_default_build_is_a_build_of_its_species(species):
    assert resolve_build(species.default_build).species == species.key


def test_every_chain_targets_a_known_build():
    for build in GENOME_BUILDS.values():
        for target, _ in build.liftover_chains:
            assert GENOME_BUILDS[target].species == build.species


def test_chain_url_is_none_between_unchained_builds():
    canfam3, canfam4 = GENOME_BUILDS["canfam3"], GENOME_BUILDS["canfam4"]

    assert canfam3.chain_url(canfam4).endswith("canFam3ToCanFam4.over.chain.gz")
    assert canfam4.chain_url(canfam3) is None


@pytest.mark.parametrize(
    "chrom,build,expected",
    [
        (39, "canfam3", "chrX"),
        ("chr41", "canfam4", "chrX"),
        (19, "felcat9", "chrX"),
        (39, None, "chr39"),
        (12, "canfam3", "chr12"),
    ],
)
def test_ucsc_chrom_applies_the_build_renames(chrom, build, expected):
    assert ucsc_chrom(chrom, resolve_build(build)) == expected
