"""One record per genome build, the way ``species.py`` holds one per species.

Build facts used to live in six string tables across four modules: the
synonym folding in ``utils``, the UCSC routing table in ``reference_genes``,
the assembly names in ``ucsc``, the recombination maps' native build and chain
keys, and a free-standing chain URL. "canFam4 is UU_Cfam_GSD_1.0" was written
down twice, and adding a chain meant editing three of them. Each fact now lives
in one row here.

A build this table does not carry is still usable: Ensembl serves many more
assemblies than the package models, so ``resolve_build`` returns None for an
unknown name and ``assembly_token`` still folds its spelling for comparison.
"""

import re
from dataclasses import dataclass
from typing import Optional, Union

# PLINK --dog codes X as 39 and the pseudoautosomal XY as 41, and
# --chr-set 18 codes them 19 and 21. UCSC assemblies, and so UCSC chains and
# tracks, name all of them chrX.
_CANINE_UCSC_RENAMES = (("39", "X"), ("41", "X"), ("XY", "X"))
_FELINE_UCSC_RENAMES = (("19", "X"), ("21", "X"), ("XY", "X"))


@dataclass(frozen=True)
class GenomeBuild:
    """Everything pyLocusZoom knows about one genome build.

    Attributes:
        key: Canonical token, the value ``assembly_token`` folds every
            spelling of this build to.
        species: ``Species.key`` of the species the build belongs to.
        assembly_name: The assembly's own name, recorded on every gene row.
        aliases: Other spellings callers and services use for this build.
        ucsc_genome: UCSC genome that serves gene annotations for this build.
            Set only for builds Ensembl no longer serves; every other build's
            genes come from Ensembl.
        liftover_chains: ``(target key, chain URL)`` pairs for the builds
            this one can be lifted to.
        ucsc_chrom_renames: ``(code, UCSC name)`` pairs for chromosome codes
            UCSC spells differently, such as PLINK's numeric X codes. The UCSC
            name has no ``chr`` prefix.
    """

    key: str
    species: str
    assembly_name: str
    aliases: tuple[str, ...] = ()
    ucsc_genome: Optional[str] = None
    liftover_chains: tuple[tuple[str, str], ...] = ()
    ucsc_chrom_renames: tuple[tuple[str, str], ...] = ()

    def chain_url(self, target: "GenomeBuild") -> Optional[str]:
        """Return the chain URL lifting this build to ``target``, or None."""
        return dict(self.liftover_chains).get(target.key)


GENOME_BUILDS: dict[str, GenomeBuild] = {
    build.key: build
    for build in (
        GenomeBuild(
            key="canfam3",
            species="canine",
            assembly_name="CanFam3.1",
            ucsc_genome="canFam3",
            liftover_chains=(
                (
                    "canfam4",
                    "https://hgdownload.soe.ucsc.edu/gbdb/canFam3/liftOver/"
                    "canFam3ToCanFam4.over.chain.gz",
                ),
            ),
            ucsc_chrom_renames=_CANINE_UCSC_RENAMES,
        ),
        GenomeBuild(
            key="canfam4",
            species="canine",
            assembly_name="UU_Cfam_GSD_1.0",
            aliases=("CanFam4.0",),
            ucsc_genome="canFam4",
            ucsc_chrom_renames=_CANINE_UCSC_RENAMES,
        ),
        GenomeBuild(
            key="roscfam1",
            species="canine",
            assembly_name="ROS_Cfam_1.0",
            ucsc_chrom_renames=_CANINE_UCSC_RENAMES,
        ),
        GenomeBuild(
            key="felcat9",
            species="feline",
            assembly_name="Felis_catus_9.0",
            aliases=("felCat9.0",),
            ucsc_genome="felCat9",
            ucsc_chrom_renames=_FELINE_UCSC_RENAMES,
        ),
        GenomeBuild(
            key="fca126",
            species="feline",
            assembly_name="F.catus_Fca126_mat1.0",
            ucsc_chrom_renames=_FELINE_UCSC_RENAMES,
        ),
        GenomeBuild(
            key="grch37",
            species="human",
            assembly_name="GRCh37",
            aliases=("hg19", "GRCh37.p13"),
        ),
        GenomeBuild(
            key="grch38",
            species="human",
            assembly_name="GRCh38",
            aliases=("hg38", "GRCh38.p14"),
        ),
        GenomeBuild(
            key="grcm39", species="mouse", assembly_name="GRCm39", aliases=("mm39",)
        ),
    )
}


def _fold(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", name.lower())


_BY_TOKEN: dict[str, GenomeBuild] = {
    _fold(name): build
    for build in GENOME_BUILDS.values()
    for name in (
        build.key,
        build.assembly_name,
        build.ucsc_genome or "",
        *build.aliases,
    )
    if name
}


def assembly_token(name: str) -> str:
    """Reduce an assembly or build name to a comparable token.

    Strips punctuation and case, then folds every spelling of a known build to
    its key, so ``"CanFam4.0"``, ``"canfam4"`` and ``"UU_Cfam_GSD_1.0"`` all
    compare equal. An unknown name keeps its folded spelling.
    """
    token = _fold(name)
    build = _BY_TOKEN.get(token)
    return build.key if build is not None else token


def resolve_build(build: Union[str, GenomeBuild, None]) -> Optional[GenomeBuild]:
    """Resolve a build name to its record, folding case, punctuation and aliases.

    Idempotent: a record passes through unchanged.

    Args:
        build: Build name, alias, UCSC genome name, record, or None.

    Returns:
        The record, or None for None or a build this table does not carry.
    """
    if build is None or isinstance(build, GenomeBuild):
        return build
    return _BY_TOKEN.get(_fold(build))


def ucsc_chrom(chrom: Union[int, str], build: Optional[GenomeBuild] = None) -> str:
    """Return the UCSC name for a chromosome code, applying the build's renames.

    Args:
        chrom: Chromosome as a frame spells it (12, "12", "chr12", "X", 39).
        build: Build whose renames apply, or None for none.
    """
    name = str(chrom)
    if name.lower().startswith("chr"):
        name = name[3:]
    renames = dict(build.ucsc_chrom_renames) if build is not None else {}
    return f"chr{renames.get(name, name)}"
