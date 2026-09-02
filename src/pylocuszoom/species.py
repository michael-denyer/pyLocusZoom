"""One record per species, resolved once at the API boundary.

``species`` used to be a bare string that five modules interpreted with five
private tables. ``ld`` compared it literally against ``"canine"``,
``recombination`` keyed on ``"canine"``, ``ensembl`` also accepted ``"dog"``,
``manhattan`` accepted both under a second table of the same name mapping the
opposite direction, and ``plotter`` carried an inline default-build dict. The
tables disagreed, so ``species="dog"`` produced correct Ensembl genes, a PLINK
run without ``--dog``, no recombination overlay, and no assembly-mismatch
warning, with no error at any layer.

Everything one species knows now lives in one frozen record. Callers resolve a
name once, at the boundary, and thread the record. ``resolve_species`` is
idempotent, so a public function can call it on an argument that may already be
a record without knowing which it got.

A species the library reaches only through Ensembl (Ensembl serves far more
than this table names) still works: ``resolve_species`` builds an Ensembl-only
record for a name the table does not carry, with no PLINK flags, no default
build and no chromosome order. The entry points that need one of those say so
when they find it empty, and the recombination path reports the missing maps
as a warning, so an unknown name is audible without being fatal.
"""

from dataclasses import dataclass

from .exceptions import ValidationError

CANINE_CHROMOSOMES = tuple(str(i) for i in range(1, 39)) + ("X", "Y", "MT")
FELINE_CHROMOSOMES = (
    "A1",
    "A2",
    "A3",
    "B1",
    "B2",
    "B3",
    "B4",
    "C1",
    "C2",
    "D1",
    "D2",
    "D3",
    "D4",
    "E1",
    "E2",
    "E3",
    "X",
    "Y",
    "MT",
)
HUMAN_CHROMOSOMES = tuple(str(i) for i in range(1, 23)) + ("X", "Y", "MT")


@dataclass(frozen=True)
class Species:
    """Everything pyLocusZoom knows about one species.

    Attributes:
        key: Canonical name, and the key every table in the package uses.
        ensembl_name: Ensembl's own species name, for gene annotation.
        aliases: Other names callers may pass for this species.
        plink_flags: Flags PLINK needs to read this species' chromosome set.
            Empty means PLINK's default (human) set is correct.
        default_build: Genome build assumed when the caller names none, or
            None when the species has no one obvious reference.
        chromosomes: Display order for whole-genome plots. Empty means the
            package has no built-in order and the caller must supply one.
    """

    key: str
    ensembl_name: str
    aliases: tuple[str, ...] = ()
    plink_flags: tuple[str, ...] = ()
    default_build: str | None = None
    chromosomes: tuple[str, ...] = ()


SPECIES: dict[str, Species] = {
    record.key: record
    for record in (
        Species(
            key="canine",
            ensembl_name="canis_lupus_familiaris",
            aliases=("dog", "canis_familiaris"),
            plink_flags=("--dog",),
            default_build="canfam3.1",
            chromosomes=CANINE_CHROMOSOMES,
        ),
        Species(
            key="feline",
            ensembl_name="felis_catus",
            aliases=("cat",),
            plink_flags=("--chr-set", "18"),
            default_build="felCat9",
            chromosomes=FELINE_CHROMOSOMES,
        ),
        Species(
            key="human",
            ensembl_name="homo_sapiens",
            chromosomes=HUMAN_CHROMOSOMES,
        ),
        Species(key="mouse", ensembl_name="mus_musculus"),
        Species(key="rat", ensembl_name="rattus_norvegicus"),
    )
}

_BY_NAME: dict[str, Species] = {
    name: record
    for record in SPECIES.values()
    for name in (record.key, *record.aliases)
}


def resolve_species(species: str | Species | None) -> Species | None:
    """Resolve a species name to its record, folding case and aliases.

    Idempotent: a record passes through unchanged, so a public function can
    call this on an argument a caller may already have resolved.

    Args:
        species: Canonical name, alias, an already-resolved record, or None
            for no species-specific behaviour.

    A name the table does not carry becomes an Ensembl-only record: the gene
    track works for any Ensembl species, and the README promises that.

    Returns:
        The record, or None if None was passed.

    Raises:
        ValidationError: If the name is empty.
    """
    if species is None or isinstance(species, Species):
        return species
    name = species.lower()
    if not name:
        raise ValidationError("species must be a non-empty name or None")
    return _BY_NAME.get(name, Species(key=name, ensembl_name=name))


def ensembl_species_name(species: str | Species) -> str:
    """Return the Ensembl species name for a record or a caller's name."""
    return resolve_species(species).ensembl_name
