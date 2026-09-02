"""The species record is the one place a species is interpreted.

The regression these tests pin: ``species="dog"`` used to resolve correctly for
Ensembl and then degrade silently in three other subsystems, because each kept
its own table and they disagreed. Every assertion below reads one subsystem
through the alias, so a table that drifts back out of the record fails here.
"""

import pytest

from pylocuszoom import LocusZoomPlotter, Species, resolve_species
from pylocuszoom.exceptions import ValidationError
from pylocuszoom.ld import build_ld_command
from pylocuszoom.manhattan import get_chromosome_order
from pylocuszoom.recombination import RECOMB_SOURCES


def _plink_flags(plotter):
    cmd = build_ld_command(
        plink_path="plink",
        bfile_path="data",
        lead_snp="rs1",
        output_path="out",
        species=plotter.species,
    )
    return cmd[1 : cmd.index("--bfile")]


class TestDogAliasReachesEverySubsystem:
    """The three degradations the review measured, one test each."""

    def test_dog_gives_plink_the_canine_chromosome_set(self):
        plotter = LocusZoomPlotter(species="dog", log_level=None)
        assert _plink_flags(plotter) == ["--dog"]

    def test_dog_resolves_the_canine_recombination_source(self):
        plotter = LocusZoomPlotter(species="dog", log_level=None)
        assert RECOMB_SOURCES[plotter.species.key].native_build == "canfam3"

    def test_dog_defaults_to_the_canfam3_build(self):
        plotter = LocusZoomPlotter(species="dog", log_level=None)
        assert plotter.genome_build == "canfam3.1"

    def test_case_is_folded(self):
        plotter = LocusZoomPlotter(species="Canine", log_level=None)
        assert _plink_flags(plotter) == ["--dog"]
        assert plotter.genome_build == "canfam3.1"

    def test_cat_reaches_the_feline_record(self):
        plotter = LocusZoomPlotter(species="cat", log_level=None)
        assert _plink_flags(plotter) == ["--chr-set", "18"]
        assert plotter.genome_build == "felCat9"
        assert plotter.species.ensembl_name == "felis_catus"

    def test_unknown_species_is_an_ensembl_only_record(self):
        plotter = LocusZoomPlotter(species="Sus_scrofa", log_level=None)
        assert plotter.species == Species(key="sus_scrofa", ensembl_name="sus_scrofa")
        assert _plink_flags(plotter) == []
        assert plotter.genome_build is None

    def test_empty_name_is_rejected(self):
        with pytest.raises(ValidationError, match="non-empty"):
            resolve_species("")


class TestResolveSpecies:
    def test_none_stays_none(self):
        assert resolve_species(None) is None

    def test_a_record_passes_through_unchanged(self):
        canine = resolve_species("canine")
        assert resolve_species(canine) is canine

    def test_every_alias_reaches_a_record_in_the_registry(self):
        for name in ("dog", "canis_familiaris", "cat", "human", "mouse", "rat"):
            assert isinstance(resolve_species(name), Species)

    def test_the_record_is_frozen(self):
        with pytest.raises(AttributeError):
            resolve_species("canine").key = "feline"


class TestChromosomeOrderComesFromTheRecord:
    def test_canine_order_spans_38_autosomes(self):
        order = get_chromosome_order(species="dog")
        assert order[:2] == ["1", "2"]
        assert order[-4:] == ["38", "X", "Y", "MT"]

    def test_a_species_without_a_built_in_order_says_so(self):
        with pytest.raises(ValidationError, match="No built-in chromosome order"):
            get_chromosome_order(species="mouse")

    def test_the_returned_order_is_not_the_record_s_own_storage(self):
        order = get_chromosome_order(species="canine")
        order.append("bogus")
        assert "bogus" not in get_chromosome_order(species="canine")
